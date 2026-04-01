#!/usr/bin/env python3
"""
Electronic Structure Workflow Manager - SCF/PARCHG/ELF/BAND/PDOS calculations

This workflow runs sequential electronic structure calculations on refined electride structures:
1. SCF (Self-Consistent Field) - conventional cell
2. PARCHG (Partial Charge Density) - 5 energy windows, conventional cell
3. ELF (Electron Localization Function) - conventional cell
4. BAND (Band Structure) - primitive cell
5. PDOS (Projected Density of States) - primitive cell

Input:
  - Structure IDs from user input or all RELAX_DONE from REFINE_VASP_JOBS/workflow.json
  - Refined CONTCARs from REFINE_VASP_JOBS/{composition}/{struct_id}/Relax/CONTCAR

Output:
  - Electronic structure data in ELECTRONIC_JOBS/{composition}/{struct_id}/SPE/
  - Band structure in ELECTRONIC_JOBS/{composition}/{struct_id}/BAND/
  - PDOS in ELECTRONIC_JOBS/{composition}/{struct_id}/DOS/

Note: Uses conventional cell for SCF/PARCHG/ELF (visualization) and primitive cell for BAND/PDOS (efficiency)
"""

import os
import sys
import json
import time
import argparse
import warnings
import subprocess
import shutil
import numpy as np
from pathlib import Path
from datetime import datetime

from pymatgen.core import Structure, SETTINGS as PMG_SETTINGS
from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin

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

DEFAULT_VASP_CMD = "srun --mpi=pmi2 /projects/mmi/Ridwan/potcarFiles/VASP6.4/vasp.6.4.3/bin/vasp_std"
DEFAULT_SLURM_EXCLUDE_NODES = "str-c88,str-c89,str-c90,str-c91,str-c92,str-c93,str-c94,str-c95,str-c96,str-c97"


def _slurm_env_value(name, default=""):
    """Read optional SLURM/environment values safely for script templating."""
    value = os.environ.get(name, default)
    if value is None:
        return default
    # Collapse whitespace/newlines to avoid malformed script lines.
    return " ".join(str(value).split())


def resolve_pmg_vasp_psp_dir():
    """
    Resolve POTCAR root directory for pymatgen.

    Preference order:
    1) Existing PMG_VASP_PSP_DIR environment variable
    2) Known cluster-local defaults
    """
    configured = os.environ.get('PMG_VASP_PSP_DIR')
    if configured:
        return str(Path(configured).expanduser())

    for candidate in DEFAULT_PMG_VASP_PSP_DIR_CANDIDATES:
        if (candidate / 'POT_GGA_PAW_PBE').exists():
            return str(candidate)

    return None


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


def load_structure_from_contcar(contcar_path):
    """Load structure from CONTCAR and symmetrize with PyXtal."""
    if not contcar_path.exists():
        return None, None
    
    try:
        structure = Structure.from_file(str(contcar_path))
    except Exception as e:
        print(f"  Warning: Could not load CONTCAR: {e}")
        return None, None
    
    # Symmetrize with PyXtal
    sym_structure = structure
    if PYXTAL_AVAILABLE:
        tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
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
                sym_structure = adaptor.get_structure(atoms)
                break
            except Exception:
                continue
    
    # Get conventional and primitive cells
    try:
        sga = SpacegroupAnalyzer(sym_structure)
        conv_structure = sga.get_conventional_standard_structure()
        prim_structure = sga.get_primitive_standard_structure()
        return conv_structure, prim_structure
    except Exception as e:
        print(f"  Warning: Could not get conventional/primitive cells: {e}")
        return sym_structure, sym_structure


class WorkflowDatabase:
    """JSON-based database for tracking electronic structure job states."""
    
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
    
    def add_structure(self, struct_id, comp_name, base_dir):
        """Add a new structure to track."""
        self.data['structures'][struct_id] = {
            'composition': comp_name,
            'state': 'PENDING',
            'spe_job_id': None,
            'band_job_id': None,
            'dos_job_id': None,
            'band_gap': None,
            'is_semiconductor': None,
            'spe_dir': str(base_dir / struct_id / 'SPE'),
            'band_dir': str(base_dir / struct_id / 'BAND'),
            'dos_dir': str(base_dir / struct_id / 'DOS'),
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
        running_states = ['SC_RUNNING', 'PARCHG_RUNNING', 'ELF_RUNNING', 
                         'BAND_RUNNING', 'PDOS_RUNNING']
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


class ElectronicWorkflowManager:
    """Manages electronic structure calculations (SCF/PARCHG/ELF/BAND/PDOS)."""
    
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
    
    def create_vasp_inputs(self, structure, job_dir, job_type='sc'):
        """Create VASP input files for electronic structure calculations."""
        job_dir = Path(job_dir)
        job_dir.mkdir(parents=True, exist_ok=True)

        common_static_incar = {
            'PREC': 'Accurate',
            'ALGO': 'Normal',
            'ADDGRID': True,
            'EDIFF': 1e-7,
            'IBRION': -1,
            'NSW': 0,
            'ISPIN': 1,
            'LASPH': True,
            'LREAL': False,
            'NPAR': 4,
        }
        
        if job_type == 'sc':
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    **common_static_incar,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ISYM': 0,
                    'LWAVE': True,
                    'LCHARG': True,
                    'LAECHG': False,
                    'LORBIT': 0,
                },
                user_kpoints_settings={'reciprocal_density': 250}
            )
        elif job_type == 'elf':
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    **common_static_incar,
                    'ISMEAR': -5,
                    'ISYM': 0,
                    'LELF': True,
                    'ISTART': 1,
                    'ICHARG': 11,
                    'NEDOS': 1000,
                    'LWAVE': False,
                    'LCHARG': False,
                    'LAECHG': False,
                    'LORBIT': 10,
                },
                user_kpoints_settings={'reciprocal_density': 250}
            )
        elif job_type == 'band':
            # Band structure with line-mode k-path
            # divisions=40: number of k-points between each pair of high-symmetry points
            kpath = HighSymmKpath(structure)
            kpoints = Kpoints.automatic_linemode(divisions=40, ibz=kpath)
            
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    **common_static_incar,
                    'ISMEAR': 0,
                    'SIGMA': 0.05,
                    'ICHARG': 2,
                    'LORBIT': 11,
                    'LWAVE': False,
                    'LCHARG': False,
                }
            )
            self.ensure_pmg_vasp_psp_dir()
            vis.write_input(job_dir)
            # Overwrite KPOINTS with explicit line-mode k-path
            kpoints.write_file(str(job_dir / 'KPOINTS'))
            # Keep POSCAR (primitive cell for band structure)
            return job_dir
        elif job_type == 'dos':
            # PDOS calculation with primitive cell
            vis = MPStaticSet(structure, 
                user_incar_settings={
                    **common_static_incar,
                    'ISMEAR': -5,
                    'ICHARG': 2,
                    'LORBIT': 11,
                    'NEDOS': 3000,
                    'LWAVE': False,
                    'LCHARG': False,
                },
                user_kpoints_settings={'reciprocal_density': 250}
            )
        else:
            raise ValueError(f"Unknown job_type: {job_type}")
        
        self.ensure_pmg_vasp_psp_dir()
        vis.write_input(job_dir)
        
        
        return job_dir
    
    def create_spe_script(self, struct_id, spe_dir):
        """Create SLURM job script for SPE workflow (SCF→PARCHG→ELF)."""
        spe_dir = Path(spe_dir).resolve()
        script_path = spe_dir / 'job.sh'
        psp_dir = self.ensure_pmg_vasp_psp_dir()
        slurm_partitions = _slurm_env_value('VASP_JOB_PARTITIONS', 'Apus,Orion,Nebula')
        slurm_constraint = _slurm_env_value('VASP_JOB_CONSTRAINT', '')
        slurm_exclude = _slurm_env_value('VASP_JOB_EXCLUDE_NODES', DEFAULT_SLURM_EXCLUDE_NODES)
        vasp_cmd = _slurm_env_value('VASP_CMD', DEFAULT_VASP_CMD)
        constraint_line = f"#SBATCH --constraint={slurm_constraint}\n" if slurm_constraint else ""
        exclude_line = f"#SBATCH --exclude={slurm_exclude}\n" if slurm_exclude else ""
        
        script = f"""#!/bin/bash
#SBATCH --job-name={struct_id}_spe
#SBATCH --partition={slurm_partitions}
{constraint_line}{exclude_line}#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output={spe_dir}/vasp_%j.out
#SBATCH --error={spe_dir}/vasp_%j.err

# Load modules
module purge
module load intel/mkl/2024.0 intel/2024 intel-mpi/2021.11
ulimit -s unlimited

# Set environment
export OMP_NUM_THREADS=1
export PMG_VASP_PSP_DIR={psp_dir}

# Intel MPI settings
if [ -e /opt/slurm/lib/libpmi.so ]; then
  export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi.so
else
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so.0
fi
export I_MPI_FABRICS=shm:ofi

# VASP executable (use srun for SLURM-native MPI launching)
VASP_CMD="{vasp_cmd}"

echo "========================================"
echo "SPE Workflow: {struct_id}"
echo "========================================"
echo "Start time: $(date)"
echo ""

# ========================================
# Stage 1: SC Calculation
# ========================================
echo ""
echo "========================================"
echo "Stage 1: SC Calculation"
echo "========================================"
cd {spe_dir}

if [ ! -f "POSCAR" ] || [ ! -s "POSCAR" ]; then
    echo "ERROR: POSCAR not found in SPE directory"
    touch SC_FAILED
    exit 1
fi

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

python3 generate_parchg_incars.py vasprun.xml-SC 250

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

echo ""
echo "========================================"
echo "SPE Job Completed Successfully"
echo "========================================"
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        return script_path
    
    def create_band_script(self, struct_id, band_dir):
        """Create SLURM job script for band structure calculation."""
        band_dir = Path(band_dir).resolve()
        script_path = band_dir / 'job.sh'
        psp_dir = self.ensure_pmg_vasp_psp_dir()
        slurm_partitions = _slurm_env_value('VASP_JOB_PARTITIONS', 'Apus,Orion,Nebula')
        slurm_constraint = _slurm_env_value('VASP_JOB_CONSTRAINT', '')
        slurm_exclude = _slurm_env_value('VASP_JOB_EXCLUDE_NODES', DEFAULT_SLURM_EXCLUDE_NODES)
        vasp_cmd = _slurm_env_value('VASP_CMD', DEFAULT_VASP_CMD)
        constraint_line = f"#SBATCH --constraint={slurm_constraint}\n" if slurm_constraint else ""
        exclude_line = f"#SBATCH --exclude={slurm_exclude}\n" if slurm_exclude else ""
        
        script = f"""#!/bin/bash
#SBATCH --job-name={struct_id}_band
#SBATCH --partition={slurm_partitions}
{constraint_line}{exclude_line}#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output={band_dir}/vasp_%j.out
#SBATCH --error={band_dir}/vasp_%j.err

# Load modules
module purge
module load intel/mkl/2024.0 intel/2024 intel-mpi/2021.11
ulimit -s unlimited

# Set environment
export OMP_NUM_THREADS=1
export PMG_VASP_PSP_DIR={psp_dir}

# Intel MPI settings
if [ -e /opt/slurm/lib/libpmi.so ]; then
  export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi.so
else
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so.0
fi
export I_MPI_FABRICS=shm:ofi

# VASP executable (use srun for SLURM-native MPI launching)
VASP_CMD="{vasp_cmd}"

echo "========================================"
echo "Band Structure: {struct_id}"
echo "========================================"
echo "Start time: $(date)"
echo ""

cd {band_dir}

if [ ! -f "POSCAR" ] || [ ! -s "POSCAR" ]; then
    echo "ERROR: POSCAR not found"
    touch BAND_FAILED
    exit 1
fi

echo "Running band structure calculation (primitive cell)"

$VASP_CMD

EXIT_CODE=$?
echo "Band structure exit code: $EXIT_CODE"

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: Band structure failed with exit code $EXIT_CODE"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    touch BAND_FAILED
    exit 1
fi

if [ ! -f "vasprun.xml" ] || [ ! -s "vasprun.xml" ]; then
    echo "ERROR: vasprun.xml missing or empty"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    touch BAND_FAILED
    exit 1
fi

echo "Band structure completed successfully"
touch BAND_DONE

# Cleanup
rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null

echo ""
echo "========================================"
echo "Band Structure Complete"
echo "========================================"
echo "End time: $(date)"
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        return script_path
    
    def create_dos_script(self, struct_id, dos_dir):
        """Create SLURM job script for PDOS calculation."""
        dos_dir = Path(dos_dir).resolve()
        script_path = dos_dir / 'job.sh'
        psp_dir = self.ensure_pmg_vasp_psp_dir()
        slurm_partitions = _slurm_env_value('VASP_JOB_PARTITIONS', 'Apus,Orion,Nebula')
        slurm_constraint = _slurm_env_value('VASP_JOB_CONSTRAINT', '')
        slurm_exclude = _slurm_env_value('VASP_JOB_EXCLUDE_NODES', DEFAULT_SLURM_EXCLUDE_NODES)
        vasp_cmd = _slurm_env_value('VASP_CMD', DEFAULT_VASP_CMD)
        constraint_line = f"#SBATCH --constraint={slurm_constraint}\n" if slurm_constraint else ""
        exclude_line = f"#SBATCH --exclude={slurm_exclude}\n" if slurm_exclude else ""
        
        script = f"""#!/bin/bash
#SBATCH --job-name={struct_id}_dos
#SBATCH --partition={slurm_partitions}
{constraint_line}{exclude_line}#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output={dos_dir}/vasp_%j.out
#SBATCH --error={dos_dir}/vasp_%j.err

# Load modules
module purge
module load intel/mkl/2024.0 intel/2024 intel-mpi/2021.11
ulimit -s unlimited

# Set environment
export OMP_NUM_THREADS=1
export PMG_VASP_PSP_DIR={psp_dir}

# Intel MPI settings
if [ -e /opt/slurm/lib/libpmi.so ]; then
  export I_MPI_PMI_LIBRARY=/opt/slurm/lib/libpmi.so
else
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so.0
fi
export I_MPI_FABRICS=shm:ofi

# VASP executable (use srun for SLURM-native MPI launching)
VASP_CMD="{vasp_cmd}"

echo "========================================"
echo "PDOS: {struct_id}"
echo "========================================"
echo "Start time: $(date)"
echo ""

cd {dos_dir}

if [ ! -f "POSCAR" ] || [ ! -s "POSCAR" ]; then
    echo "ERROR: POSCAR not found"
    touch PDOS_FAILED
    exit 1
fi

echo "Running PDOS calculation (primitive cell)"

$VASP_CMD

EXIT_CODE=$?
echo "PDOS exit code: $EXIT_CODE"

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: PDOS failed with exit code $EXIT_CODE"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    touch PDOS_FAILED
    exit 1
fi

if [ ! -f "vasprun.xml" ] || [ ! -s "vasprun.xml" ]; then
    echo "ERROR: vasprun.xml missing or empty"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    touch PDOS_FAILED
    exit 1
fi

if [ ! -f "DOSCAR" ] || [ ! -s "DOSCAR" ]; then
    echo "ERROR: DOSCAR missing or empty"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    touch PDOS_FAILED
    exit 1
fi

echo "PDOS completed successfully"
touch PDOS_DONE

# Cleanup
rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null

echo ""
echo "========================================"
echo "PDOS Complete"
echo "========================================"
echo "End time: $(date)"
"""
        
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        return script_path
    
    def submit_job(self, script_path):
        """Submit SLURM job."""
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
        """Check SLURM job status."""
        result = subprocess.run(
            ['squeue', '-j', job_id, '-h', '-o', '%T'],
            capture_output=True,
            text=True
        )
        return 'RUNNING' if result.stdout.strip() else 'NOTFOUND'
    
    def check_spe_status(self, spe_dir):
        """Check SPE job status (SC/PARCHG/ELF stages)."""
        spe_dir = Path(spe_dir)
        
        # Check for VASP_FAILED marker (written at any stage failure)
        if (spe_dir / 'VASP_FAILED').exists():
            # Read the failure message to determine which stage failed
            try:
                with open(spe_dir / 'VASP_FAILED', 'r') as f:
                    failure_msg = f.read().strip()
                    if 'at SC' in failure_msg or 'at PARCHG INCAR generation' in failure_msg:
                        return 'SC_FAILED'
                    elif 'at PARCHG' in failure_msg:
                        return 'PARCHG_FAILED'
                    elif 'at ELF' in failure_msg:
                        return 'ELF_FAILED'
                    else:
                        return 'SC_FAILED'  # Default to SC_FAILED if unclear
            except Exception:
                return 'SC_FAILED'
        
        # Check progress through stages
        if (spe_dir / 'SC_DONE').exists():
            if (spe_dir / 'PARCHG_DONE').exists():
                if (spe_dir / 'ELF_DONE').exists():
                    return 'ELF_DONE'  # All SPE stages complete
                else:
                    return 'ELF_RUNNING'
            else:
                return 'PARCHG_RUNNING'
        else:
            return 'SC_RUNNING'
    
    def check_band_status(self, band_dir):
        """Check BAND job status."""
        band_dir = Path(band_dir)
        if (band_dir / 'BAND_DONE').exists():
            return 'DONE'
        elif (band_dir / 'BAND_FAILED').exists():
            return 'FAILED'
        else:
            return 'UNKNOWN'
    
    def check_dos_status(self, dos_dir):
        """Check DOS job status."""
        dos_dir = Path(dos_dir)
        if (dos_dir / 'PDOS_DONE').exists():
            return 'DONE'
        elif (dos_dir / 'PDOS_FAILED').exists():
            return 'FAILED'
        else:
            return 'UNKNOWN'
    
    def submit_spe_job(self, struct_id, conv_structure):
        """Submit SPE job (SC→PARCHG→ELF) for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        spe_dir = Path(sdata['spe_dir'])
        
        print(f"  Submitting SPE job: {struct_id}")
        
        try:
            # Create SPE inputs (conventional cell)
            self.create_vasp_inputs(conv_structure, spe_dir, 'sc')
            shutil.copy2(spe_dir / 'INCAR', spe_dir / 'INCAR-SC')
            
            elf_temp = spe_dir / 'temp_elf'
            elf_temp.mkdir(exist_ok=True)
            self.create_vasp_inputs(conv_structure, elf_temp, 'elf')
            shutil.move(elf_temp / 'INCAR', spe_dir / 'INCAR-ELF')
            shutil.rmtree(elf_temp)
            
            # Copy PARCHG helper script
            helper_script = Path(__file__).parent / 'generate_parchg_incars.py'
            if helper_script.exists():
                shutil.copy2(helper_script, spe_dir / 'generate_parchg_incars.py')
            else:
                raise FileNotFoundError("generate_parchg_incars.py not found in workflow directory")
            
            # Create SLURM script
            script = self.create_spe_script(struct_id, spe_dir)
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'SC_RUNNING', spe_job_id=job_id)
            print(f"    SPE job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'SC_FAILED', error=str(e))
            return False
    
    def submit_band_job(self, struct_id, prim_structure):
        """Submit BAND job for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        band_dir = Path(sdata['band_dir'])
        
        print(f"  Submitting BAND job: {struct_id}")
        
        try:
            # Create BAND inputs (primitive cell)
            self.create_vasp_inputs(prim_structure, band_dir, 'band')
            
            # Create SLURM script
            script = self.create_band_script(struct_id, band_dir)
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'BAND_RUNNING', band_job_id=job_id)
            print(f"    BAND job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'BAND_FAILED', error=str(e))
            return False
    
    def submit_dos_job(self, struct_id, prim_structure):
        """Submit DOS job for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        dos_dir = Path(sdata['dos_dir'])
        
        print(f"  Submitting DOS job: {struct_id}")
        
        try:
            # Create DOS inputs (primitive cell)
            self.create_vasp_inputs(prim_structure, dos_dir, 'dos')
            
            # Create SLURM script
            script = self.create_dos_script(struct_id, dos_dir)
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'PDOS_RUNNING', dos_job_id=job_id)
            print(f"    DOS job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'PDOS_FAILED', error=str(e))
            return False
    
    def update_structure_status(self, struct_id):
        """Update status of a structure with separate SC/PARCHG/ELF/BAND/DOS stages."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return
        
        state = sdata['state']
        spe_dir = Path(sdata['spe_dir'])
        band_dir = Path(sdata['band_dir'])
        dos_dir = Path(sdata['dos_dir'])
        
        # Check SPE job stages (SC → PARCHG → ELF)
        if state in ['SC_RUNNING', 'PARCHG_RUNNING', 'ELF_RUNNING']:
            # Check local status first (stage markers)
            local_status = self.check_spe_status(spe_dir)
            
            # Handle failures immediately
            if local_status == 'SC_FAILED':
                self.db.update_state(struct_id, 'SC_FAILED', error='SC calculation failed')
                print(f"  {struct_id}: SC FAILED")
            elif local_status == 'PARCHG_FAILED':
                self.db.update_state(struct_id, 'PARCHG_FAILED', error='PARCHG calculation failed')
                print(f"  {struct_id}: PARCHG FAILED")
            elif local_status == 'ELF_FAILED':
                self.db.update_state(struct_id, 'ELF_FAILED', error='ELF calculation failed')
                print(f"  {struct_id}: ELF FAILED")
            # Handle stage transitions (only if different from current state)
            elif local_status == 'PARCHG_RUNNING' and state != 'PARCHG_RUNNING':
                # Parse band gap from SC vasprun.xml
                vasprun_sc = spe_dir / 'vasprun.xml-SC'
                if vasprun_sc.exists():
                    band_gap, is_semiconductor = parse_band_gap_from_vasprun(vasprun_sc)
                    self.db.update_state(struct_id, 'PARCHG_RUNNING', 
                                       band_gap=band_gap, 
                                       is_semiconductor=is_semiconductor)
                else:
                    self.db.update_state(struct_id, 'PARCHG_RUNNING')
                print(f"  {struct_id}: SC completed, PARCHG running")
            elif local_status == 'ELF_RUNNING' and state != 'ELF_RUNNING':
                self.db.update_state(struct_id, 'ELF_RUNNING')
                print(f"  {struct_id}: PARCHG completed, ELF running")
            elif local_status == 'ELF_DONE':
                self.db.update_state(struct_id, 'ELF_DONE')
                print(f"  {struct_id}: ELF completed (SPE workflow done)")
            # Check if job terminated without leaving markers
            elif local_status == 'SC_RUNNING':
                job_status = self.check_job_status(sdata['spe_job_id'])
                if job_status == 'NOTFOUND':
                    # Job finished but no markers - crashed
                    self.db.update_state(struct_id, 'SC_FAILED', error='Job crashed without stage markers')
                    print(f"  {struct_id}: SC FAILED (crash)")
        
        # Check BAND job
        elif state == 'BAND_RUNNING':
            job_status = self.check_job_status(sdata['band_job_id'])
            if job_status == 'NOTFOUND':
                local_status = self.check_band_status(band_dir)
                if local_status == 'DONE':
                    self.db.update_state(struct_id, 'BAND_DONE')
                    print(f"  {struct_id}: BAND completed")
                elif local_status == 'FAILED':
                    self.db.update_state(struct_id, 'BAND_FAILED', error='Band structure failed')
                    print(f"  {struct_id}: BAND FAILED")
                else:
                    self.db.update_state(struct_id, 'BAND_FAILED', error='Job crashed')
                    print(f"  {struct_id}: BAND FAILED (crash)")
        
        # Check DOS job
        elif state == 'PDOS_RUNNING':
            job_status = self.check_job_status(sdata['dos_job_id'])
            if job_status == 'NOTFOUND':
                local_status = self.check_dos_status(dos_dir)
                if local_status == 'DONE':
                    self.db.update_state(struct_id, 'PDOS_DONE')
                    print(f"  {struct_id}: PDOS completed")
                elif local_status == 'FAILED':
                    self.db.update_state(struct_id, 'PDOS_FAILED', error='PDOS failed')
                    print(f"  {struct_id}: PDOS FAILED")
                else:
                    self.db.update_state(struct_id, 'PDOS_FAILED', error='Job crashed')
                    print(f"  {struct_id}: PDOS FAILED (crash)")
    
    def monitor_and_submit(self, structures_dict):
        """Main monitoring loop with three separate job submissions."""
        print("\n" + "="*70)
        print("Starting electronic workflow monitoring...")
        print(f"Max concurrent: {self.max_concurrent}")
        print(f"Check interval: {self.check_interval}s")
        print("Workflow: SPE → BAND+DOS (parallel)")
        print("="*70 + "\n")
        sys.stdout.flush()
        
        while True:
            print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking status...")
            sys.stdout.flush()
            
            # Update status of all structures
            for struct_id in list(self.db.data['structures'].keys()):
                self.update_structure_status(struct_id)
            
            running_count = self.db.get_running_count()
            print(f"Currently running: {running_count}/{self.max_concurrent}")
            
            # Submit new jobs (respecting max_concurrent limit)
            for struct_id in list(self.db.data['structures'].keys()):
                if running_count >= self.max_concurrent:
                    break
                
                sdata = self.db.get_structure(struct_id)
                state = sdata['state']
                conv_struct, prim_struct = structures_dict.get(struct_id, (None, None))
                
                if not conv_struct or not prim_struct:
                    continue
                
                # Submit SPE job for PENDING structures
                if state == 'PENDING':
                    if self.submit_spe_job(struct_id, conv_struct):
                        running_count += 1
                
                # Submit BAND job for ELF_DONE structures (SPE complete)
                elif state == 'ELF_DONE':
                    if self.submit_band_job(struct_id, prim_struct):
                        running_count += 1
                
                # Submit DOS job for BAND_DONE structures (after BAND submitted)
                elif state == 'BAND_DONE':
                    if self.submit_dos_job(struct_id, prim_struct):
                        running_count += 1
            
            # Check for completed workflows
            for struct_id in list(self.db.data['structures'].keys()):
                sdata = self.db.get_structure(struct_id)
                band_dir = Path(sdata['band_dir'])
                dos_dir = Path(sdata['dos_dir'])
                
                # Mark as complete if both BAND and DOS are done
                if sdata['state'] == 'PDOS_DONE':
                    if (band_dir / 'BAND_DONE').exists():
                        self.db.update_state(struct_id, 'COMPLETE')
                        print(f"  {struct_id}: All stages complete (SPE+BAND+DOS)")
            
            stats = self.db.get_stats()
            print("\nStatistics:")
            for state, count in sorted(stats['states'].items()):
                print(f"  {state}: {count}")
            sys.stdout.flush()
            
            # Check if all done
            pending = len(self.db.get_by_state('PENDING'))
            elf_done = len(self.db.get_by_state('ELF_DONE'))
            band_done = len(self.db.get_by_state('BAND_DONE'))
            
            if running_count == 0 and pending == 0 and elf_done == 0 and band_done == 0:
                completed = len(self.db.get_by_state('COMPLETE'))
                failed = (len(self.db.get_by_state('SC_FAILED')) + 
                         len(self.db.get_by_state('PARCHG_FAILED')) + 
                         len(self.db.get_by_state('ELF_FAILED')) + 
                         len(self.db.get_by_state('BAND_FAILED')) + 
                         len(self.db.get_by_state('PDOS_FAILED')))
                total = stats['total']
                if completed + failed >= total:
                    print("\n" + "="*70)
                    print("All workflows completed!")
                    print(f"Completed: {completed}/{total}")
                    print(f"Failed: {failed}/{total}")
                    print("="*70)
                    sys.stdout.flush()
                    break
            
            print(f"\nSleeping for {self.check_interval}s...")
            sys.stdout.flush()
            time.sleep(self.check_interval)
    
    def initialize_structures(self, refine_jobs_dir, output_dir, structure_ids=None):
        """Initialize structures from refined workflow."""
        refine_jobs_dir = Path(refine_jobs_dir)
        output_dir = Path(output_dir)
        
        print("="*70)
        print("Electronic Structure Workflow Initialization")
        print("="*70)
        print(f"Refined jobs dir: {refine_jobs_dir}")
        print(f"Output directory: {output_dir}")
        print(f"Max concurrent: {self.max_concurrent}")
        print("="*70 + "\n")
        
        # Load refined workflow database
        refine_db_path = refine_jobs_dir / 'workflow.json'
        if not refine_db_path.exists():
            print(f"ERROR: Refined workflow database not found: {refine_db_path}")
            return {}
        
        with open(refine_db_path, 'r') as f:
            refine_db = json.load(f)
        
        # Filter structures
        if structure_ids:
            print(f"User-specified structure IDs: {len(structure_ids)}")
            target_structures = {sid: sdata for sid, sdata in refine_db['structures'].items() 
                               if sid in structure_ids}
        else:
            print("No structure IDs specified - processing all RELAX_DONE structures")
            target_structures = {sid: sdata for sid, sdata in refine_db['structures'].items() 
                               if sdata['state'] == 'RELAX_DONE'}
        
        print(f"Found {len(target_structures)} structures to process\n")
        
        # Load structures
        structures_dict = {}
        loaded = 0
        failed = 0
        
        for struct_id, sdata in target_structures.items():
            relax_dir = Path(sdata['relax_dir'])
            contcar_path = relax_dir / 'CONTCAR'
            
            if not contcar_path.exists():
                print(f"  Warning: CONTCAR not found for {struct_id}")
                failed += 1
                continue
            
            conv_struct, prim_struct = load_structure_from_contcar(contcar_path)
            if conv_struct is None or prim_struct is None:
                print(f"  Warning: Could not load structure for {struct_id}")
                failed += 1
                continue
            
            structures_dict[struct_id] = (conv_struct, prim_struct)
            
            if struct_id not in self.db.data['structures']:
                comp_name = sdata['composition']
                self.db.add_structure(struct_id, comp_name, output_dir / comp_name)
            
            loaded += 1
            if loaded % 50 == 0:
                print(f"  Loaded {loaded} structures...")
        
        print(f"\nSuccessfully loaded: {loaded}")
        if failed > 0:
            print(f"Failed to load: {failed}")
        
        self.db.data['config'] = {
            'max_concurrent': self.max_concurrent,
            'refine_jobs_dir': str(refine_jobs_dir),
            'output_dir': str(output_dir)
        }
        self.db.save()
        
        return structures_dict


def main():
    parser = argparse.ArgumentParser(
        description="Electronic Structure Workflow Manager - SCF/PARCHG/ELF/BAND/PDOS"
    )
    parser.add_argument(
        '--refine-jobs',
        type=str,
        required=True,
        help="Path to REFINE_VASP_JOBS directory"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./ELECTRONIC_JOBS',
        help="Output directory for electronic structure jobs (default: ./ELECTRONIC_JOBS)"
    )
    parser.add_argument(
        '--structure-ids',
        type=str,
        nargs='+',
        default=None,
        help="Specific structure IDs to process (default: all RELAX_DONE)"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="JSON database file (default: workflow.json)"
    )
    parser.add_argument(
        '--max-concurrent',
        type=int,
        default=10,
        help="Max concurrent structures (default: 10)"
    )
    parser.add_argument(
        '--check-interval',
        type=int,
        default=60,
        help="Status check interval in seconds (default: 60)"
    )
    parser.add_argument(
        '--init-only',
        action='store_true',
        help="Only initialize database, don't start monitoring"
    )
    
    args = parser.parse_args()
    
    refine_jobs_dir = Path(args.refine_jobs).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = output_dir / args.db
    
    if not refine_jobs_dir.exists():
        print(f"ERROR: Refined jobs directory not found: {refine_jobs_dir}")
        sys.exit(1)
    
    manager = ElectronicWorkflowManager(
        db_path=db_path,
        max_concurrent=args.max_concurrent,
        check_interval=args.check_interval
    )
    
    structures_dict = manager.initialize_structures(
        refine_jobs_dir=refine_jobs_dir,
        output_dir=output_dir,
        structure_ids=args.structure_ids
    )
    
    if not structures_dict:
        print("\nERROR: No structures to process!")
        sys.exit(1)
    
    if args.init_only:
        print("\n" + "="*70)
        print("Initialization complete!")
        print(f"Database: {db_path}")
        print("="*70)
        return
    
    try:
        manager.monitor_and_submit(structures_dict)
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user.")
        print(f"Database saved to: {db_path}")
        print("Resume with same command to continue.")


if __name__ == '__main__':
    main()
