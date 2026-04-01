#!/usr/bin/env python3
"""
VASP Workflow Manager - Batch submission with dynamic monitoring
No MongoDB/FireWorks - uses local JSON database for job tracking

Features:
- Sequential VASP workflow: Relax → SC → (PARCHG for semiconductors) → ELF
- PARCHG calculations for semiconductors
- Dynamic job submission with concurrency control

Note: Pre-screening should be done separately with prescreen.py before running this workflow.
"""

import os
import sys
import json
import time
import argparse
import zipfile
import warnings
import subprocess
import shutil
import numpy as np
from pathlib import Path
from io import StringIO
from datetime import datetime

from pymatgen.core import Structure, SETTINGS as PMG_SETTINGS
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.ase import AseAtomsAdaptor

try:
    from pyxtal import pyxtal
    PYXTAL_AVAILABLE = True
except ImportError:
    PYXTAL_AVAILABLE = False
    print("WARNING: PyXtal not available. Structures will not be symmetrized.")

warnings.filterwarnings('ignore', category=UserWarning, message='.*POTCAR data with symbol.*')
warnings.filterwarnings('ignore', message='Using UFloat objects with std_dev==0')

DEFAULT_PMG_VASP_PSP_DIR_CANDIDATES = [
    Path('/projects/mmi/Ridwan/potcarFiles/VASP6.4/potpaw_PBE'),
    Path.home() / 'apps' / 'PBE64',
]


def _slurm_env_value(name, default=""):
    """Read optional SLURM-related env vars safely for script templating."""
    value = os.environ.get(name, default)
    if value is None:
        return default
    # Collapse whitespace/newlines to avoid malformed #SBATCH lines.
    return " ".join(str(value).split())


def resolve_pmg_vasp_psp_dir():
    """
    Resolve POTCAR root directory for pymatgen.

    Preference order:
    1) Existing PMG_VASP_PSP_DIR environment variable
    2) Known cluster-local defaults in DEFAULT_PMG_VASP_PSP_DIR_CANDIDATES
    """
    configured = os.environ.get('PMG_VASP_PSP_DIR')
    if configured:
        return str(Path(configured).expanduser())

    for candidate in DEFAULT_PMG_VASP_PSP_DIR_CANDIDATES:
        potcar_family = candidate / 'POT_GGA_PAW_PBE'
        if potcar_family.exists():
            return str(candidate)

    return None


def check_electronic_convergence_outcar(outcar_path):
    """
    Check electronic convergence from OUTCAR file.
    
    For timed-out jobs, vasprun.xml is incomplete/corrupted. This function checks
    OUTCAR for "aborting loop because EDIFF is reached" marker, which indicates
    electronic SCF converged in at least one ionic step.
    
    Returns:
        bool: True if electronic convergence was achieved
    """
    if not outcar_path.exists():
        return False
    
    try:
        with open(outcar_path, 'r') as f:
            content = f.read()
        return 'aborting loop because EDIFF is reached' in content
    except Exception:
        return False


def parse_band_gap_from_vasprun(vasprun_path):
    """
    Parse band gap from vasprun.xml using pymatgen's vasprun parser.
    
    Handles both semiconductors and metals with ISMEAR=0 where efermi is None.
    
    Returns:
        tuple: (band_gap, is_semiconductor)
    """
    try:
        vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=True)
        
        if not vr.converged_electronic:
            print(f"    Warning: Electronic SCF not converged")
            return None, False
        
        # For systems with ISMEAR=0, efermi is None
        efermi = vr.efermi
        if efermi is None:
            eigenvalues = vr.eigenvalues
            spin = Spin.up if Spin.up in eigenvalues else list(eigenvalues.keys())[0]
            eigs = eigenvalues[spin]
            
            energies = eigs[:, :, 0]
            occupations = eigs[:, :, 1]
            
            flat_energies = energies.flatten()
            flat_occupations = occupations.flatten()
            sort_idx = np.argsort(flat_energies)
            sorted_energies = flat_energies[sort_idx]
            sorted_occupations = flat_occupations[sort_idx]
            
            # Find Fermi level: energy where occupation drops below threshold
            # For semiconductors: clear jump from 1.0 to 0.0
            # For metals: gradual decrease through partial occupations
            occ_threshold = 0.5
            fermi_idx = np.where(sorted_occupations < occ_threshold)[0]
            
            if len(fermi_idx) == 0:
                efermi = sorted_energies[-1]
            elif fermi_idx[0] == 0:
                efermi = sorted_energies[0]
            else:
                last_occ_idx = fermi_idx[0] - 1
                first_unocc_idx = fermi_idx[0]
                
                e_below = sorted_energies[last_occ_idx]
                e_above = sorted_energies[first_unocc_idx]
                occ_below = sorted_occupations[last_occ_idx]
                occ_above = sorted_occupations[first_unocc_idx]
                
                if occ_below > occ_above:
                    efermi = (e_below + e_above) / 2
                else:
                    efermi = e_below
        
        band_structure = vr.get_band_structure(line_mode=False, efermi=efermi)
        
        if band_structure is None:
            print(f"    Warning: Band structure not available")
            print(f"      VASP finished but didn't calculate eigenvalues properly")
            return None, False
        
        band_gap_dict = band_structure.get_band_gap()
        if band_gap_dict is None or 'energy' not in band_gap_dict:
            print(f"    Warning: Band gap data not available in band structure")
            return None, False
        
        band_gap = band_gap_dict['energy']
        is_semiconductor = not band_structure.is_metal()
        
        # Display efermi info
        if vr.efermi is not None:
            print(f"    e_fermi = {vr.efermi:.4f} eV, band_gap = {band_gap:.4f} eV, is_metal = {band_structure.is_metal()}")
        else:
            print(f"    e_fermi = {efermi:.4f} eV (calculated), band_gap = {band_gap:.4f} eV, is_metal = {band_structure.is_metal()}")
        return band_gap, is_semiconductor
        
    except Exception as e:
        error_msg = str(e).split('\n')[0]
        print(f"    Warning: Could not parse vasprun.xml: {error_msg}")
        print(f"      File: {vasprun_path}")
        return None, False


def get_parchg_incar_settings_band(iband, grid):
    """Get INCAR settings for band-specific PARCHG."""
    return {
        'IBAND': iband,
        'LPARD': True,
        'LAECHG': False,
        'LVHAR': False,
        'LCHARG': False,
        'NGXF': grid[0],
        'NGYF': grid[1],
        'NGZF': grid[2],
        'LWAVE': False,
        'ISTART': 1,
        'ICHARG': 11
    }


def get_parchg_incar_settings_energy(eint, grid):
    """Get INCAR settings for energy-window PARCHG."""
    return {
        'EINT': eint,
        'NBMOD': -3,
        'LPARD': True,
        'LAECHG': False,
        'LVHAR': False,
        'LCHARG': False,
        'NGXF': grid[0],
        'NGYF': grid[1],
        'NGZF': grid[2],
        'LWAVE': False,
        'ISTART': 1,
        'ICHARG': 11
    }


class WorkflowDatabase:
    """Simple JSON-based database for tracking job states."""
    
    def __init__(self, db_path):
        self.db_path = Path(db_path)
        self.data = {'structures': {}, 'config': {}}
        self.load()
    
    def load(self):
        """Load database from JSON file."""
        if self.db_path.exists():
            with open(self.db_path, 'r') as f:
                self.data = json.load(f)
    
    def save(self):
        """Save database to JSON file."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = self.db_path.with_suffix('.tmp')
        with open(tmp_path, 'w') as f:
            json.dump(self.data, f, indent=2)
        tmp_path.replace(self.db_path)
    
    def add_structure(self, struct_id, comp_name, struct_idx, base_dir, chemsys=None):
        """Add a new structure to track."""
        self.data['structures'][struct_id] = {
            'composition': comp_name,
            'chemsys': chemsys,
            'structure_idx': struct_idx,
            'state': 'PENDING',
            'relax_job_id': None,
            'spe_job_id': None,
            'band_gap': None,
            'is_semiconductor': None,
            'relax_dir': str(base_dir / struct_id / 'Relax'),
            'spe_dir': str(base_dir / struct_id / 'SPE'),
            'last_updated': datetime.now().isoformat(),
            'error': None
        }
        self.save()
    
    def update_state(self, struct_id, state, **kwargs):
        """Update structure state and additional fields."""
        if struct_id in self.data['structures']:
            self.data['structures'][struct_id]['state'] = state
            self.data['structures'][struct_id]['last_updated'] = datetime.now().isoformat()
            for key, value in kwargs.items():
                self.data['structures'][struct_id][key] = value
            self.save()
    
    def get_structure(self, struct_id):
        """Get structure data."""
        return self.data['structures'].get(struct_id)
    
    def get_by_state(self, state):
        """Get all structures in a specific state."""
        return [sid for sid, sdata in self.data['structures'].items() 
                if sdata['state'] == state]
    
    def get_running_count(self):
        """Count structures currently running (any stage)."""
        running_states = ['RELAX_RUNNING', 'SC_RUNNING', 'ELF_RUNNING', 'PARCHG_RUNNING']
        return sum(1 for s in self.data['structures'].values() 
                   if s['state'] in running_states)
    
    def get_stats(self):
        """Get overall statistics."""
        states = {}
        for s in self.data['structures'].values():
            state = s['state']
            states[state] = states.get(state, 0) + 1
        return {
            'total': len(self.data['structures']),
            'states': states,
            'running': self.get_running_count()
        }


class VASPWorkflowManager:
    """Manages VASP job submission and monitoring with batch control."""
    
    def __init__(self, db_path, max_concurrent=10, check_interval=60):
        self.db = WorkflowDatabase(db_path)
        self.max_concurrent = max_concurrent
        self.check_interval = check_interval
        self.pmg_vasp_psp_dir = resolve_pmg_vasp_psp_dir()
        if self.pmg_vasp_psp_dir:
            os.environ['PMG_VASP_PSP_DIR'] = self.pmg_vasp_psp_dir
            PMG_SETTINGS['PMG_VASP_PSP_DIR'] = self.pmg_vasp_psp_dir

    def ensure_pmg_vasp_psp_dir(self):
        """Ensure PMG_VASP_PSP_DIR is available for pymatgen POTCAR generation."""
        if self.pmg_vasp_psp_dir:
            os.environ['PMG_VASP_PSP_DIR'] = self.pmg_vasp_psp_dir
            PMG_SETTINGS['PMG_VASP_PSP_DIR'] = self.pmg_vasp_psp_dir
            return self.pmg_vasp_psp_dir

        self.pmg_vasp_psp_dir = resolve_pmg_vasp_psp_dir()
        if self.pmg_vasp_psp_dir:
            os.environ['PMG_VASP_PSP_DIR'] = self.pmg_vasp_psp_dir
            PMG_SETTINGS['PMG_VASP_PSP_DIR'] = self.pmg_vasp_psp_dir
            return self.pmg_vasp_psp_dir

        candidate_list = ', '.join(str(p) for p in DEFAULT_PMG_VASP_PSP_DIR_CANDIDATES)
        raise RuntimeError(
            "PMG_VASP_PSP_DIR is not set and no default POTCAR directory was found. "
            f"Set PMG_VASP_PSP_DIR or place POTCAR files under one of: {candidate_list}"
        )
    
    def read_structures_from_zip(self, zip_path, max_structures=None):
        """Read CIF structures from zip file."""
        structures = []
        with zipfile.ZipFile(zip_path, 'r') as zf:
            cif_files = sorted([f for f in zf.namelist() if f.endswith('.cif')])
            if max_structures:
                cif_files = cif_files[:max_structures]
            
            for cif_file in cif_files:
                try:
                    with zf.open(cif_file) as f:
                        cif_content = f.read().decode('utf-8')
                        parser = CifParser(StringIO(cif_content))
                        structure = parser.parse_structures(primitive=True)[0]
                        structures.append(structure)
                except Exception as e:
                    print(f"  Warning: Could not parse {cif_file}: {e}")
        return structures
    
    def create_vasp_inputs(self, structure, job_dir, job_type='relax', parchg_settings=None):
        """
        Create VASP input files using pymatgen.
        
        Note: 
        - For relax jobs, structure is symmetrized using PyXtal before VASP input generation
        - For SC, ELF, and PARCHG jobs, POSCAR is deleted after generation
        - The relaxed structure (CONTCAR from Relax) will be copied by the SLURM script
        """
        job_dir = Path(job_dir)
        job_dir.mkdir(parents=True, exist_ok=True)
        self.ensure_pmg_vasp_psp_dir()
        
        if job_type == 'relax':
            # Symmetrize structure using PyXtal with progressive tolerance
            if PYXTAL_AVAILABLE:
                tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
                symmetrized = False
                for tol in tolerances:
                    try:
                        adaptor = AseAtomsAdaptor()
                        xtal = pyxtal()
                        xtal.from_seed(structure, tol=tol)
                        if not xtal.valid:
                            continue
                        if len(xtal.check_short_distances(r=0.5)) > 0:
                            continue
                        atoms = xtal.to_ase()
                        structure = adaptor.get_structure(atoms)
                        symmetrized = True
                        break
                    except Exception:
                        continue
                
                if not symmetrized:
                    print(f"    Warning: Could not symmetrize structure with tolerances {tolerances}")
                    print(f"    Proceeding with original structure...")
            
            vis = MPRelaxSet(structure, 
                user_incar_settings={
                    'PREC': 'Accurate',
                    'ALGO': 'Fast',
                    'ADDGRID': True,
                    'EDIFF': 1e-6,
                    'EDIFFG': -0.01,
                    'IBRION': 2,
                    'ISIF': 3,
                    'NELM': 60,
                    'NSW': 30,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ISPIN': 1,
                    'NPAR': 4,
                    'POTIM': 0.3,
                    'LWAVE': False,
                    'LCHARG': False,
                    'LAECHG': False,
                    'LASPH': True,
                    'LORBIT': 0,
                },
                user_kpoints_settings={'reciprocal_density': 64}
            )
        elif job_type == 'sc':
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    'PREC': 'Accurate',
                    'ALGO': 'Normal',
                    'ADDGRID': True,
                    'EDIFF': 1e-6,
                    'IBRION': -1,
                    'NSW': 0,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ISPIN': 1,
                    'ISYM': 0,
                    'NPAR': 4,
                    'LWAVE': True,
                    'LCHARG': True,
                    'LAECHG': False,
                    'LASPH': True,
                    'LORBIT': 0,
                },
                user_kpoints_settings={'reciprocal_density': 64}
            )
        elif job_type == 'elf':
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    'PREC': 'Accurate',
                    'ALGO': 'Normal',
                    'ADDGRID': True,
                    'EDIFF': 1e-6,
                    'IBRION': -1,
                    'NSW': 0,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ISPIN': 1,
                    'ISYM': 0,
                    'NPAR': 4,
                    'LELF': True,
                    'ISTART': 1,
                    'ICHARG': 11,
                    'NEDOS': 1000,
                    'LWAVE': False,
                    'LCHARG': False,
                    'LAECHG': False,
                    'LASPH': True,
                    'LORBIT': 10,
                },
                user_kpoints_settings={'reciprocal_density': 64}
            )
        elif job_type == 'parchg':
            if not parchg_settings:
                raise ValueError("parchg_settings required for PARCHG job")
            
            base_settings = {
                'PREC': 'Accurate',
                'ALGO': 'Normal',
                'EDIFF': 1e-6,
                'IBRION': -1,
                'NSW': 0,
                'ISMEAR': 1,
                'SIGMA': 0.2,
                'ISPIN': 1,
                'ISYM': 0,
                'NPAR': 4,
            }
            base_settings.update(parchg_settings)
            vis = MPStaticSet(structure, 
                user_incar_settings=base_settings,
                user_kpoints_settings={'reciprocal_density': 64}
            )
        else:
            raise ValueError(f"Unknown job_type: {job_type}")
        
        vis.write_input(job_dir)
        
        # For SC, ELF, and PARCHG, delete the generated POSCAR
        # We will copy CONTCAR from Relax in the SLURM script
        if job_type in ['sc', 'elf', 'parchg']:
            poscar_path = Path(job_dir) / 'POSCAR'
            if poscar_path.exists():
                poscar_path.unlink()
        
        return job_dir
    
    def create_slurm_script(self, job_dir, job_name, job_type='relax', prev_dir=None, relax_dir=None, parchg_label=None):
        """Create SLURM submission script."""
        job_dir = Path(job_dir).resolve()
        script_path = job_dir / 'job.sh'
        psp_dir = self.ensure_pmg_vasp_psp_dir()
        slurm_partitions = _slurm_env_value('VASP_JOB_PARTITIONS', 'Apus,Orion,Nebula')
        slurm_constraint = _slurm_env_value('VASP_JOB_CONSTRAINT', '')
        slurm_exclude = _slurm_env_value('VASP_JOB_EXCLUDE_NODES', '')
        constraint_line = f"#SBATCH --constraint={slurm_constraint}\n" if slurm_constraint else ""
        exclude_line = f"#SBATCH --exclude={slurm_exclude}\n" if slurm_exclude else ""
        
        # Convert all paths to absolute paths
        if prev_dir:
            prev_dir = Path(prev_dir).resolve()
        if relax_dir:
            relax_dir = Path(relax_dir).resolve()
        
        script = f"""#!/bin/bash
#SBATCH --job-name={job_name}_{job_type}
#SBATCH --partition={slurm_partitions}
{constraint_line}{exclude_line}#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output={job_dir}/vasp_%j.out
#SBATCH --error={job_dir}/vasp_%j.err

# Load modules
module purge
module load intel/mkl/2024.0 intel/2024 intel-mpi/2021.11
ulimit -s unlimited

# Set environment
export OMP_NUM_THREADS=1
export PMG_VASP_PSP_DIR={psp_dir}

# Intel MPI settings for SLURM
if [ -e /opt/slurm/lib/libpmi.so ]; then
  export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi.so
else
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so.0
fi
export I_MPI_FABRICS=shm:ofi

# VASP executable (use srun for SLURM-native MPI launching)
VASP_CMD="srun --mpi=pmi2 /projects/mmi/Ridwan/potcarFiles/VASP6.4/vasp.6.4.3/bin/vasp_std"


# Change to job directory
cd {job_dir}

"""
        
        # SPE job type: unified SC-PARCHG-ELF sequential calculation
        if job_type == 'spe':
            script += f"""
# ========================================
# SPE Job: SC -> PARCHG (5x) -> ELF
# ========================================

# Copy CONTCAR from relaxation as POSCAR
echo "Checking for relaxed structure..."
if [ ! -f "{relax_dir}/CONTCAR" ] || [ ! -s "{relax_dir}/CONTCAR" ]; then
    echo "ERROR: CONTCAR not found or empty in {relax_dir}"
    echo "VASP calculation failed at SC" > VASP_FAILED
    exit 1
fi

cp "{relax_dir}/CONTCAR" POSCAR
echo "Copied CONTCAR from {relax_dir} as POSCAR (size: $(du -h POSCAR | cut -f1))"

# Verify copy succeeded
if [ ! -f "POSCAR" ] || [ ! -s "POSCAR" ]; then
    echo "ERROR: Failed to copy CONTCAR as POSCAR"
    echo "VASP calculation failed at SC" > VASP_FAILED
    exit 1
fi
echo "Successfully verified POSCAR is present"

# ========================================
# Stage 1: SC Calculation
# ========================================
echo ""
echo "========================================  "
echo "Stage 1: SC Calculation"
echo "========================================"
cp INCAR-SC INCAR

$VASP_CMD

EXIT_CODE=$?
echo "SC calculation exit code: $EXIT_CODE"

if [ $EXIT_CODE -ne 0 ]; then
    echo "VASP calculation failed at SC" > VASP_FAILED
    echo "ERROR: SC calculation failed with exit code $EXIT_CODE"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

# Verify SC outputs
if [ ! -f "CHGCAR" ] || [ ! -s "CHGCAR" ]; then
    echo "VASP calculation failed at SC" > VASP_FAILED
    echo "ERROR: CHGCAR missing or empty"
    [ ! -f "CHGCAR" ] && echo "  CHGCAR does not exist"
    [ -f "CHGCAR" ] && [ ! -s "CHGCAR" ] && echo "  CHGCAR is empty"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

if [ ! -f "WAVECAR" ] || [ ! -s "WAVECAR" ]; then
    echo "VASP calculation failed at SC" > VASP_FAILED
    echo "ERROR: WAVECAR missing or empty"
    [ ! -f "WAVECAR" ] && echo "  WAVECAR does not exist"
    [ -f "WAVECAR" ] && [ ! -s "WAVECAR" ] && echo "  WAVECAR is empty"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

if [ ! -f "vasprun.xml" ] || [ ! -s "vasprun.xml" ]; then
    echo "VASP calculation failed at SC" > VASP_FAILED
    echo "ERROR: vasprun.xml missing or empty"
    [ ! -f "vasprun.xml" ] && echo "  vasprun.xml does not exist"
    [ -f "vasprun.xml" ] && [ ! -s "vasprun.xml" ] && echo "  vasprun.xml is empty"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

echo "Verified CHGCAR ($(du -h CHGCAR | cut -f1)), WAVECAR ($(du -h WAVECAR | cut -f1)), and vasprun.xml exist"

# Save SC outputs with -SC suffix
mv OSZICAR OSZICAR-SC
cp vasprun.xml vasprun.xml-SC
cp CHGCAR CHGCAR-SC
cp WAVECAR WAVECAR-SC
touch SC_DONE
echo "SC calculation completed successfully"

# Cleanup SC calculation outputs for next stage
rm -f OSZICAR vasprun.xml OUTCAR 2>/dev/null

# ========================================
# Stage 2: Generate PARCHG INCARs
# ========================================
echo ""
echo "========================================"
echo "Stage 2: Generating PARCHG INCARs"
echo "========================================"

python3 generate_parchg_incars.py vasprun.xml-SC 64

if [ $? -ne 0 ]; then
    echo "VASP calculation failed at PARCHG INCAR generation" > VASP_FAILED
    echo "ERROR: Failed to generate PARCHG INCARs"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

# ========================================
# Stage 3: PARCHG Calculations
# ========================================
echo ""
echo "========================================"
echo "Stage 3: PARCHG Calculations"
echo "========================================"

for label in e0025 e05 e10 band0 band1; do
    echo ""
    echo "Running PARCHG-$label..."
    cp INCAR-PARCHG-$label INCAR
    cp CHGCAR-SC CHGCAR
    cp WAVECAR-SC WAVECAR
    
    $VASP_CMD
    
    EXIT_CODE=$?
    echo "PARCHG-$label exit code: $EXIT_CODE"
    
    if [ $EXIT_CODE -ne 0 ]; then
        echo "VASP calculation failed at PARCHG-$label" > VASP_FAILED
        echo "ERROR: PARCHG-$label failed with exit code $EXIT_CODE"
        # Clean up large intermediate files to save disk space
        rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
        exit 1
    fi
    
    if [ ! -f "PARCHG" ] || [ ! -s "PARCHG" ]; then
        echo "VASP calculation failed at PARCHG-$label" > VASP_FAILED
        echo "ERROR: PARCHG file missing or empty for $label"
        # Clean up large intermediate files to save disk space
        rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
        exit 1
    fi
    
    mv PARCHG PARCHG-$label
    echo "PARCHG-$label completed successfully"
    
    # Cleanup before next PARCHG run
    rm -f CHG WFULL TMPCAR OSZICAR vasprun.xml OUTCAR 2>/dev/null
done

touch PARCHG_DONE
echo "All PARCHG calculations completed successfully"

# Cleanup before ELF calculation
rm -f CHG WFULL TMPCAR OSZICAR vasprun.xml OUTCAR 2>/dev/null

# ========================================
# Stage 4: ELF Calculation
# ========================================
echo ""
echo "========================================"
echo "Stage 4: ELF Calculation"
echo "========================================"
cp INCAR-ELF INCAR
cp CHGCAR-SC CHGCAR
cp WAVECAR-SC WAVECAR

$VASP_CMD

EXIT_CODE=$?
echo "ELF calculation exit code: $EXIT_CODE"

if [ $EXIT_CODE -ne 0 ]; then
    echo "VASP calculation failed at ELF" > VASP_FAILED
    echo "ERROR: ELF calculation failed with exit code $EXIT_CODE"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

if [ ! -f "ELFCAR" ] || [ ! -s "ELFCAR" ]; then
    echo "VASP calculation failed at ELF" > VASP_FAILED
    echo "ERROR: ELFCAR missing or empty"
    [ ! -f "ELFCAR" ] && echo "  ELFCAR does not exist"
    [ -f "ELFCAR" ] && [ ! -s "ELFCAR" ] && echo "  ELFCAR is empty"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

echo "Verified ELFCAR exists ($(du -h ELFCAR | cut -f1))"
echo "ELF calculation completed successfully"
touch ELF_DONE

# ========================================
# Stage 5: Compress PARCHG files
# ========================================
echo ""
echo "========================================"
echo "Stage 5: Compress PARCHG files"
echo "========================================"

# Validate and compress PARCHG files to save disk space
# Only compress if we have all 5 PARCHG files (band0, band1, e0025, e05, e10)
echo "Validating PARCHG files..."

PARCHG_FILES="PARCHG-band0 PARCHG-band1 PARCHG-e0025 PARCHG-e05 PARCHG-e10"
ALL_VALID=true

# Check if all 5 expected PARCHG files exist and are non-empty
for pf in $PARCHG_FILES; do
    if [ ! -f "$pf" ] || [ ! -s "$pf" ]; then
        echo "  Warning: $pf is missing or empty"
        ALL_VALID=false
    fi
done

if [ "$ALL_VALID" = true ]; then
    if [ ! -f "PARCHG.tar.gz" ]; then
        echo "  All 5 PARCHG files validated (non-empty)"
        if tar czf PARCHG.tar.gz $PARCHG_FILES 2>/dev/null; then
            echo "  Created PARCHG.tar.gz"
            
            # Remove original PARCHG-* files to save space
            rm -f PARCHG-*
            echo "  Removed original PARCHG-* files (archive preserved)"
            echo "  Disk space saved: PARCHG files compressed into single archive"
        else
            echo "  ERROR: Failed to create PARCHG.tar.gz"
            echo "  Keeping original PARCHG files"
        fi
    else
        echo "  PARCHG.tar.gz already exists"
    fi
else
    echo "  ERROR: Not all PARCHG files are valid - skipping compression"
    echo "  This indicates PARCHG calculations may have failed"
    echo "  Keeping original files for debugging"
fi

# Mark entire SPE workflow as complete
touch VASP_DONE

# ========================================
# Stage 6: Cleanup
# ========================================
echo ""
echo "========================================"
echo "Stage 6: Cleanup"
echo "========================================"

# Clean up SPE directory
rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null
rm -f CHGCAR-SC WAVECAR-SC vasprun.xml 2>/dev/null
echo "Cleaned up SPE/ directory"

# Clean up Relax directory
if [ -d "{relax_dir}" ]; then
    rm -f "{relax_dir}/CHGCAR" "{relax_dir}/CHG" "{relax_dir}/WAVECAR" "{relax_dir}/WFULL" "{relax_dir}/AECCAR"* "{relax_dir}/TMPCAR" 2>/dev/null
    echo "Cleaned up Relax/ directory"
fi

echo ""
echo "========================================"
echo "SPE Job Completed Successfully"
echo "========================================"
"""
            with open(script_path, 'w') as f:
                f.write(script)
            os.chmod(script_path, 0o755)
            return script_path
        
        # Relax job type - only job type besides 'spe' that's still used
        if job_type == 'relax':
            script += f"""
# Run VASP
echo "Starting VASP calculation: {job_type}"
echo "Working directory: $(pwd)"
echo "VASP command: $VASP_CMD"
echo "Start time: $(date)"

$VASP_CMD

EXIT_CODE=$?

echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"

# Check if successful
if [ $EXIT_CODE -eq 0 ]; then
    # Verify critical files for Relax calculation
    if [ -f "CONTCAR" ] && [ -s "CONTCAR" ]; then
        echo "VASP calculation completed successfully"
        echo "Verified CONTCAR exists"
        
        # Clean up large unnecessary files to save disk space
        rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null
        
        touch VASP_DONE
    else
        echo "VASP calculation failed: CONTCAR missing/empty"
        # Clean up large intermediate files to save disk space
        rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
        touch VASP_FAILED
    fi
else
    echo "VASP calculation failed with exit code $EXIT_CODE"
    # Clean up large intermediate files to save disk space
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    touch VASP_FAILED
fi
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        
        os.chmod(script_path, 0o755)
        return script_path
    
    def submit_job(self, script_path):
        """Submit a SLURM job and return job ID."""
        result = subprocess.run(
            ['sbatch', str(script_path)],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            return job_id
        else:
            raise RuntimeError(f"sbatch failed: {result.stderr}")
    
    def check_job_status(self, job_id):
        """Check SLURM job status. Returns: RUNNING, COMPLETED, FAILED, or NOTFOUND."""
        result = subprocess.run(
            ['squeue', '-j', job_id, '-h', '-o', '%T'],
            capture_output=True,
            text=True
        )
        
        if result.stdout.strip():
            slurm_state = result.stdout.strip()
            if slurm_state in ['RUNNING', 'PENDING', 'CONFIGURING']:
                return 'RUNNING'
            else:
                return 'RUNNING'  # Other states still in queue
        else:
            return 'NOTFOUND'
    
    def check_local_status(self, job_dir):
        """Check local directory for completion markers."""
        job_dir = Path(job_dir)
        if (job_dir / 'VASP_DONE').exists():
            return 'DONE'
        elif (job_dir / 'VASP_FAILED').exists():
            return 'FAILED'
        else:
            return 'UNKNOWN'
    
    def submit_relax(self, struct_id, structure):
        """Submit relaxation job for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        relax_dir = Path(sdata['relax_dir'])
        job_name = struct_id
        
        print(f"  Submitting Relax: {struct_id}")
        
        try:
            self.create_vasp_inputs(structure, relax_dir, 'relax')
            script = self.create_slurm_script(relax_dir, job_name, 'relax')
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'RELAX_RUNNING', relax_job_id=job_id)
            print(f"    Relax job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'RELAX_FAILED', error=str(e))
            return False
    
    def submit_spe(self, struct_id, structure):
        """
        Submit unified SPE (SC-PARCHG-ELF) job for a structure.
        
        The SPE job runs sequentially in one SLURM script:
        1. SC calculation
        2. PARCHG calculations (5 energy windows)
        3. ELF calculation
        
        Stage markers (SC_DONE, PARCHG_DONE, VASP_DONE) track progress.
        """
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        spe_dir = Path(sdata['spe_dir'])
        relax_dir = Path(sdata['relax_dir'])
        job_name = struct_id
        
        print(f"  Submitting SPE: {struct_id}")
        
        try:
            spe_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate INCAR-SC using sc job type
            sc_temp_dir = spe_dir / 'temp_sc'
            sc_temp_dir.mkdir(exist_ok=True)
            self.create_vasp_inputs(structure, sc_temp_dir, 'sc')
            shutil.move(sc_temp_dir / 'INCAR', spe_dir / 'INCAR-SC')
            shutil.rmtree(sc_temp_dir)
            
            # Generate INCAR-ELF using elf job type
            elf_temp_dir = spe_dir / 'temp_elf'
            elf_temp_dir.mkdir(exist_ok=True)
            self.create_vasp_inputs(structure, elf_temp_dir, 'elf')
            shutil.move(elf_temp_dir / 'INCAR', spe_dir / 'INCAR-ELF')
            shutil.rmtree(elf_temp_dir)
            
            # Generate POSCAR, POTCAR, KPOINTS using sc job type
            # (These will be deleted, CONTCAR will be copied in the script)
            self.create_vasp_inputs(structure, spe_dir, 'sc')
            
            # Copy PARCHG INCAR generator helper script
            helper_script = Path(__file__).parent / 'generate_parchg_incars.py'
            if helper_script.exists():
                shutil.copy2(helper_script, spe_dir / 'generate_parchg_incars.py')
            else:
                raise FileNotFoundError("generate_parchg_incars.py not found in workflow directory")
            
            # Create unified SPE SLURM script
            script = self.create_slurm_script(spe_dir, job_name, 'spe', relax_dir=relax_dir)
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'SC_RUNNING', spe_job_id=job_id)
            print(f"    SPE job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'SC_FAILED', error=str(e))
            return False
    
    def update_structure_status(self, struct_id):
        """Check and update status of a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return
        
        state = sdata['state']
        
        if state == 'RELAX_RUNNING':
            job_status = self.check_job_status(sdata['relax_job_id'])
            if job_status == 'NOTFOUND':
                local_status = self.check_local_status(sdata['relax_dir'])
                if local_status == 'DONE':
                    # Check convergence before marking as done
                    vasprun_path = Path(sdata['relax_dir']) / 'vasprun.xml'
                    if vasprun_path.exists():
                        try:
                            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
                            if not vr.converged_electronic:
                                self.db.update_state(struct_id, 'RELAX_FAILED', 
                                                   error='Electronic SCF not converged')
                                print(f"  {struct_id}: Relax FAILED (electronic not converged)")
                            else:
                                self.db.update_state(struct_id, 'RELAX_DONE')
                                print(f"  {struct_id}: Relax completed (electronic converged)")
                        except Exception as e:
                            self.db.update_state(struct_id, 'RELAX_FAILED', 
                                               error=f'Could not check convergence: {e}')
                            print(f"  {struct_id}: Relax FAILED (convergence check error)")
                    else:
                        self.db.update_state(struct_id, 'RELAX_FAILED', 
                                           error='vasprun.xml not found')
                        print(f"  {struct_id}: Relax FAILED (vasprun.xml missing)")
                elif local_status == 'FAILED':
                    self.db.update_state(struct_id, 'RELAX_FAILED')
                    print(f"  {struct_id}: Relax failed")
                else:
                    # Job not in queue and no completion marker - check if timed out
                    relax_dir = Path(sdata['relax_dir'])
                    err_files = list(relax_dir.glob('vasp_*.err'))
                    is_timeout = False
                    
                    if err_files:
                        err_file = max(err_files, key=lambda p: p.stat().st_mtime)
                        try:
                            with open(err_file, 'r') as f:
                                if 'DUE TO TIME LIMIT' in f.read():
                                    is_timeout = True
                        except Exception:
                            pass
                    
                    if is_timeout:
                        # Check if electronic converged (using OUTCAR) and CONTCAR exists
                        outcar_path = relax_dir / 'OUTCAR'
                        contcar_path = relax_dir / 'CONTCAR'
                        
                        if not outcar_path.exists():
                            self.db.update_state(struct_id, 'RELAX_FAILED',
                                               error='Job timed out, OUTCAR not found')
                            print(f"  {struct_id}: Relax FAILED (timeout, OUTCAR missing)")
                        elif not contcar_path.exists() or contcar_path.stat().st_size == 0:
                            self.db.update_state(struct_id, 'RELAX_FAILED',
                                               error='Job timed out, CONTCAR missing/empty')
                            print(f"  {struct_id}: Relax FAILED (timeout, CONTCAR missing)")
                        elif check_electronic_convergence_outcar(outcar_path):
                            self.db.update_state(struct_id, 'RELAX_TMOUT',
                                               error='Relaxation timed out but electronic converged')
                            print(f"  {struct_id}: Relax TMOUT (electronic converged, proceeding)")
                        else:
                            self.db.update_state(struct_id, 'RELAX_FAILED',
                                               error='Job timed out, electronic not converged')
                            print(f"  {struct_id}: Relax FAILED (timeout, electronic not converged)")
                    else:
                        self.db.update_state(struct_id, 'RELAX_FAILED', 
                                           error='Job terminated without completion marker (crash)')
                        print(f"  {struct_id}: Relax FAILED (crash)")
        
        elif state == 'SC_RUNNING':
            # SPE workflow: check stage markers
            spe_dir = Path(sdata['spe_dir'])
            job_status = self.check_job_status(sdata['spe_job_id'])
            
            # Check for stage completion marker
            if (spe_dir / 'SC_DONE').exists():
                # SC stage completed, parse band gap before transition to PARCHG_RUNNING
                vasprun_sc = spe_dir / 'vasprun.xml-SC'
                if vasprun_sc.exists():
                    band_gap, is_semiconductor = parse_band_gap_from_vasprun(vasprun_sc)
                    self.db.update_state(struct_id, 'PARCHG_RUNNING', 
                                       band_gap=band_gap, 
                                       is_semiconductor=is_semiconductor)
                else:
                    self.db.update_state(struct_id, 'PARCHG_RUNNING')
                print(f"  {struct_id}: SC completed, PARCHG stage starting")
            elif (spe_dir / 'VASP_FAILED').exists():
                error_msg = (spe_dir / 'VASP_FAILED').read_text().strip()
                if 'at SC' in error_msg:
                    self.db.update_state(struct_id, 'SC_FAILED', error=error_msg)
                    print(f"  {struct_id}: SC FAILED - {error_msg}")
            elif job_status == 'NOTFOUND':
                # Job finished but no markers - likely timed out or crashed
                err_files = list(spe_dir.glob('vasp_*.err'))
                is_timeout = False
                if err_files:
                    err_file = max(err_files, key=lambda p: p.stat().st_mtime)
                    try:
                        with open(err_file, 'r') as f:
                            if 'DUE TO TIME LIMIT' in f.read():
                                is_timeout = True
                    except Exception:
                        pass
                
                if is_timeout:
                    self.db.update_state(struct_id, 'SC_FAILED',
                                       error='SPE job timed out during SC stage')
                    print(f"  {struct_id}: SC FAILED (timeout)")
                else:
                    self.db.update_state(struct_id, 'SC_FAILED',
                                       error='SPE job terminated without completion marker (crash)')
                    print(f"  {struct_id}: SC FAILED (crash)")
        
        elif state == 'PARCHG_RUNNING':
            # SPE workflow: PARCHG stage running, check for completion
            spe_dir = Path(sdata['spe_dir'])
            job_status = self.check_job_status(sdata['spe_job_id'])
            
            if (spe_dir / 'PARCHG_DONE').exists():
                # PARCHG stage completed, transition to ELF_RUNNING (same job continues)
                self.db.update_state(struct_id, 'ELF_RUNNING')
                print(f"  {struct_id}: PARCHG completed, ELF stage starting")
            elif (spe_dir / 'VASP_FAILED').exists():
                error_msg = (spe_dir / 'VASP_FAILED').read_text().strip()
                if 'at PARCHG' in error_msg:
                    self.db.update_state(struct_id, 'PARCHG_FAILED', error=error_msg)
                    print(f"  {struct_id}: PARCHG FAILED - {error_msg}")
            elif job_status == 'NOTFOUND':
                # Job finished but no markers - likely timed out or crashed
                err_files = list(spe_dir.glob('vasp_*.err'))
                is_timeout = False
                if err_files:
                    err_file = max(err_files, key=lambda p: p.stat().st_mtime)
                    try:
                        with open(err_file, 'r') as f:
                            if 'DUE TO TIME LIMIT' in f.read():
                                is_timeout = True
                    except Exception:
                        pass
                
                if is_timeout:
                    self.db.update_state(struct_id, 'PARCHG_FAILED',
                                       error='SPE job timed out during PARCHG stage')
                    print(f"  {struct_id}: PARCHG FAILED (timeout)")
                else:
                    self.db.update_state(struct_id, 'PARCHG_FAILED',
                                       error='SPE job terminated without completion marker (crash)')
                    print(f"  {struct_id}: PARCHG FAILED (crash)")
        
        elif state == 'ELF_RUNNING':
            # SPE workflow: ELF stage running, check for completion
            spe_dir = Path(sdata['spe_dir'])
            job_status = self.check_job_status(sdata['spe_job_id'])
            
            if (spe_dir / 'VASP_DONE').exists():
                self.db.update_state(struct_id, 'ELF_DONE')
                print(f"  {struct_id}: ELF completed")
            elif (spe_dir / 'VASP_FAILED').exists():
                error_msg = (spe_dir / 'VASP_FAILED').read_text().strip()
                if 'at ELF' in error_msg:
                    self.db.update_state(struct_id, 'ELF_FAILED', error=error_msg)
                    print(f"  {struct_id}: ELF FAILED - {error_msg}")
            elif job_status == 'NOTFOUND':
                # Job finished but no markers - likely timed out or crashed
                err_files = list(spe_dir.glob('vasp_*.err'))
                is_timeout = False
                if err_files:
                    err_file = max(err_files, key=lambda p: p.stat().st_mtime)
                    try:
                        with open(err_file, 'r') as f:
                            if 'DUE TO TIME LIMIT' in f.read():
                                is_timeout = True
                    except Exception:
                        pass
                
                if is_timeout:
                    self.db.update_state(struct_id, 'ELF_FAILED',
                                       error='SPE job timed out during ELF stage')
                    print(f"  {struct_id}: ELF FAILED (timeout)")
                else:
                    self.db.update_state(struct_id, 'ELF_FAILED',
                                       error='SPE job terminated without completion marker (crash)')
                    print(f"  {struct_id}: ELF FAILED (crash)")
    
    def monitor_and_submit(self, structures_dict):
        """Main monitoring loop that checks status and submits new jobs."""
        print("\n" + "="*70)
        print("Starting workflow monitoring loop...")
        print(f"Max concurrent structures: {self.max_concurrent}")
        print(f"Check interval: {self.check_interval}s")
        print("="*70 + "\n")
        sys.stdout.flush()
        
        while True:
            print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking job status...")
            sys.stdout.flush()
            
            # Update status of all running jobs
            for struct_id in list(self.db.data['structures'].keys()):
                self.update_structure_status(struct_id)
            
            # Check if we can submit new jobs
            running_count = self.db.get_running_count()
            print(f"Currently running: {running_count}/{self.max_concurrent}")
            
            # Try to submit next stage for completed jobs
            for struct_id in list(self.db.data['structures'].keys()):
                if running_count >= self.max_concurrent:
                    break
                
                sdata = self.db.get_structure(struct_id)
                state = sdata['state']
                structure = structures_dict.get(struct_id)
                
                if not structure:
                    continue
                
                if state == 'PENDING':
                    if self.submit_relax(struct_id, structure):
                        running_count += 1
                
                elif state == 'RELAX_DONE' or state == 'RELAX_TMOUT':
                    # Submit unified SPE job (SC->PARCHG->ELF)
                    if self.submit_spe(struct_id, structure):
                        running_count += 1
            
            # Print statistics
            stats = self.db.get_stats()
            print("\nStatistics:")
            for state, count in sorted(stats['states'].items()):
                print(f"  {state}: {count}")
            sys.stdout.flush()
            
            # Check if all done
            pending_count = len(self.db.get_by_state('PENDING'))
            if running_count == 0 and pending_count == 0:
                completed = len(self.db.get_by_state('ELF_DONE'))
                total = stats['total']
                failed_count = (len(self.db.get_by_state('ELF_FAILED')) + 
                               len(self.db.get_by_state('PARCHG_FAILED')) + 
                               len(self.db.get_by_state('RELAX_FAILED')) + 
                               len(self.db.get_by_state('SC_FAILED')))
                if completed + failed_count >= total:
                    print("\n" + "="*70)
                    print("All workflows completed!")
                    print("="*70)
                    sys.stdout.flush()
                    break
            
            print(f"\nSleeping for {self.check_interval}s...")
            sys.stdout.flush()
            time.sleep(self.check_interval)
    
    def initialize_structures(self, results_dir, output_dir,
                             max_compositions=None, max_structures=5):
        """Scan CIF directory and initialize database."""
        results_dir = Path(results_dir)
        output_dir = Path(output_dir)
        
        print("="*70)
        print("Initializing VASP Workflow")
        print("="*70)
        print(f"Results directory: {results_dir}")
        print(f"Output directory: {output_dir}")
        print(f"Max concurrent: {self.max_concurrent}")
        print(f"Max compositions: {max_compositions or 'all'}")
        print(f"Max structures: {max_structures}")

        print("="*70 + "\n")
        
        structures_dict = {}

        cif_files = sorted(results_dir.glob("*.cif"))
        if not cif_files:
            print("  No CIF files found in results directory")
            return structures_dict

        # Group CIFs by composition (prefix before "_s")
        comp_map = {}
        for cif_path in cif_files:
            name = cif_path.stem
            if "_s" not in name:
                print(f"  Skipping {cif_path.name} (no _s index)")
                continue
            comp_name = name.split("_s")[0]
            comp_map.setdefault(comp_name, []).append(cif_path)

        comp_names = sorted(comp_map.keys())
        if max_compositions:
            comp_names = comp_names[:max_compositions]

        for comp_name in comp_names:
            files = sorted(comp_map[comp_name])
            if max_structures:
                files = files[:max_structures]
            print(f"Scanning {comp_name} ({len(files)} files)...")

            added_count = 0
            for cif_path in files:
                struct_id = cif_path.stem  # e.g., Gd2Be2Fe2H13_s036
                try:
                    structure = Structure.from_file(str(cif_path))
                except Exception as e:
                    print(f"  Warning: Could not parse {cif_path.name}: {e}")
                    continue

                structures_dict[struct_id] = structure

                if struct_id not in self.db.data['structures']:
                    elements = sorted([str(el) for el in structure.composition.elements])
                    chemsys = '-'.join(elements)

                    self.db.add_structure(
                        struct_id, comp_name, None,
                        output_dir / comp_name,
                        chemsys=chemsys
                    )
                    added_count += 1

            if added_count > 0:
                print(f"  Added {added_count} structures")
        
        # Load structures from database that aren't in structures_dict yet
        # This handles resume scenarios where structures exist in DB but weren't loaded from ZIP
        print("\nChecking database for additional structures...")
        loaded_from_contcar = 0
        skipped_count = 0
        
        for struct_id, sdata in self.db.data['structures'].items():
            if struct_id in structures_dict:
                continue  # Already loaded from ZIP
            
            # Try to load from Relax/CONTCAR for structures that have been processed
            if sdata['state'] not in ['PENDING', 'RELAX_RUNNING']:
                relax_dir = Path(sdata['relax_dir'])
                contcar_path = relax_dir / 'CONTCAR'
                
                if contcar_path.exists():
                    try:
                        structure = Structure.from_file(str(contcar_path))
                        structures_dict[struct_id] = structure
                        loaded_from_contcar += 1
                    except Exception as e:
                        print(f"  Warning: Could not load {struct_id} from CONTCAR: {e}")
                        skipped_count += 1
                else:
                    skipped_count += 1
        
        if loaded_from_contcar > 0:
            print(f"  Loaded {loaded_from_contcar} structures from CONTCAR files (resume)")
        if skipped_count > 0:
            print(f"  Skipped {skipped_count} structures (no CONTCAR available)")
        
        self.db.data['config'] = {
            'max_concurrent': self.max_concurrent,
            'results_dir': str(results_dir),
            'output_dir': str(output_dir),
            'max_structures': max_structures
        }
        self.db.save()
        
        print(f"\nTotal structures ready for workflow: {len(structures_dict)}")
        
        # Report structures in database but not in structures_dict
        missing_count = len(self.db.data['structures']) - len(structures_dict)
        if missing_count > 0:
            print(f"  Note: {missing_count} structures in database but not in structures_dict")
            print(f"        (These will be skipped during monitoring)")
        
        return structures_dict


def main():
    parser = argparse.ArgumentParser(
        description="VASP Workflow Manager - Batch submission with monitoring"
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help="Directory containing CIF files (e.g., cif_ehull_1e-3_1e-2)"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='/scratch/$USER/VASP_JOBS',
        help="Output directory for VASP jobs"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="JSON database file path"
    )
    parser.add_argument(
        '--max-concurrent',
        type=int,
        default=10,
        help="Max concurrent structures running"
    )
    parser.add_argument(
        '--max-compositions',
        type=int,
        default=None,
        help="Max compositions to process"
    )
    parser.add_argument(
        '--max-structures',
        type=int,
        default=5,
        help="Max structures per composition"
    )
    parser.add_argument(
        '--check-interval',
        type=int,
        default=60,
        help="Status check interval in seconds"
    )
    parser.add_argument(
        '--init-only',
        action='store_true',
        help="Only initialize database, don't start monitoring"
    )
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = output_dir / args.db
    
    # Create workflow manager
    manager = VASPWorkflowManager(
        db_path=db_path,
        max_concurrent=args.max_concurrent,
        check_interval=args.check_interval
    )
    
    # Initialize structures
    structures_dict = manager.initialize_structures(
        results_dir=results_dir,
        output_dir=output_dir,
        max_compositions=args.max_compositions,
        max_structures=args.max_structures
    )
    
    if args.init_only:
        print("\n" + "="*70)
        print("Initialization complete!")
        print(f"Database: {db_path}")
        print("="*70)
        return
    
    # Start monitoring and submission loop
    try:
        manager.monitor_and_submit(structures_dict)
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user.")
        print(f"Database saved to: {db_path}")
        print("Resume with same command to continue.")


if __name__ == '__main__':
    main()
