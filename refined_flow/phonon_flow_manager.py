#!/usr/bin/env python3
"""
Phonon Calculation Workflow Manager - Phonopy Finite Displacement Method

This workflow runs phonon calculations using phonopy's finite displacement approach:
1. Generate supercell with phonopy
2. Create displaced structures (POSCAR-001, POSCAR-002, ...)
3. Run VASP static calculations on each displacement
4. Monitor all displacement calculations
5. Collect results for phonopy post-processing

Input:
  - Structure IDs (must be RELAX_DONE in REFINE_VASP_JOBS/workflow.json)
  - Refined CONTCARs from REFINE_VASP_JOBS/{composition}/{struct_id}/Relax/CONTCAR

Output:
  - Phonon calculations in PHONON_JOBS/{composition}/{struct_id}/PHON/
  - Displacement subdirectories: {001, 002, 003, ...}
  - Each with VASP static calculation results

Usage:
    python3 refined_flow/phonon_flow_manager.py \
        --refine-jobs REFINE_VASP-out-Boron \
        --output-dir PHONON_Boron_check2 \
        --structure-ids Na2Cs1B2H8_s007 \
        --max-concurrent 0
"""

import os
import sys
import json
import time
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

import numpy as np
from phonopy import Phonopy
from phonopy.physical_units import get_physical_units
VaspToTHz = get_physical_units().DefaultToTHz
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core import Structure, SETTINGS as PMG_SETTINGS
from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun

DEFAULT_PMG_VASP_PSP_DIR_CANDIDATES = [
    Path('/projects/mmi/Ridwan/potcarFiles/VASP6.4/potpaw_PBE'),
    Path.home() / 'apps' / 'PBE64',
]

DEFAULT_VASP_CMD = "srun --mpi=pmi2 /projects/mmi/Ridwan/potcarFiles/VASP6.4/vasp.6.4.3/bin/vasp_std"
DEFAULT_SLURM_EXCLUDE_NODES = "str-c88,str-c89,str-c90,str-c91,str-c92,str-c93,str-c94,str-c95,str-c96,str-c97"
DEFAULT_VASP_JOB_TIME = "12:00:00"
MAX_PHONON_SUPERCELL_ATOMS = 300


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


def load_structure_from_contcar(contcar_path):
    """Load structure from CONTCAR."""
    if not contcar_path.exists():
        return None
    
    try:
        structure = Structure.from_file(str(contcar_path))
        return structure
    except Exception as e:
        print(f"  Warning: Could not load CONTCAR: {e}")
        return None


def calculate_supercell_dim_from_kpoints(structure, kpoint_density=250):
    """
    Calculate supercell dimensions based on k-point mesh.
    
    Similar to the approach in generate_supercell.sh that uses vaspkit
    to get k-mesh from k_spacing and uses those as supercell dimensions.
    
    Args:
        structure: Pymatgen Structure
        kpoint_density: Reciprocal density (k-points per Å^-3)
    
    Returns:
        Supercell matrix as 3x3 diagonal array
    """
    from pymatgen.io.vasp.inputs import Kpoints
    
    # Generate k-point mesh based on reciprocal density
    kpoints = Kpoints.automatic_density(structure, kpoint_density)
    kpts = kpoints.kpts[0]  # Get the [k1, k2, k3] mesh
    
    # Use k-point mesh as supercell dimensions
    # This ensures phonon q-point mesh matches electronic k-point mesh
    supercell_matrix = np.diag(kpts)
    
    return supercell_matrix


def pymatgen_to_phonopy_atoms(structure):
    """Convert Pymatgen Structure to PhonopyAtoms object."""
    lattice = structure.lattice.matrix
    positions = structure.frac_coords
    numbers = [site.specie.Z for site in structure]
    symbols = [str(site.specie.symbol) for site in structure]
    
    return PhonopyAtoms(
        cell=lattice,
        scaled_positions=positions,
        numbers=numbers,
        symbols=symbols
    )


def estimate_supercell_atom_count(structure, supercell_matrix):
    """Estimate atom count in the phonon supercell before generating displacements."""
    multiplier = int(round(abs(np.linalg.det(np.array(supercell_matrix, dtype=float)))))
    return len(structure) * multiplier


def generate_supercell_with_phonopy(structure, phon_dir, supercell_matrix=None, auto_dim=False):
    """
    Generate phonopy supercell and displaced structures using phonopy Python API.
    
    Args:
        structure: Pymatgen Structure
        phon_dir: Directory for phonopy files
        supercell_matrix: Supercell matrix (3x3 array or list)
                         If None and auto_dim=True, calculate from k-points
        auto_dim: If True, automatically calculate supercell dimensions
    
    Returns:
        Tuple of (n_displacements, supercell_matrix_used)
    """
    phon_dir = Path(phon_dir)
    phon_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine supercell matrix
    if supercell_matrix is None:
        if auto_dim:
            # Auto-calculate based on k-point mesh
            supercell_matrix = calculate_supercell_dim_from_kpoints(structure, kpoint_density=200)
            print(f"    Auto-calculated supercell dimensions: {np.diag(supercell_matrix).astype(int)}") 
            
        else:
            # Default 2x2x2
            supercell_matrix = np.diag([2, 2, 2])
    else:
        supercell_matrix = np.array(supercell_matrix)

    estimated_atoms = estimate_supercell_atom_count(structure, supercell_matrix)
    print(f"    Estimated phonon supercell size: {estimated_atoms} atoms")
    if estimated_atoms > MAX_PHONON_SUPERCELL_ATOMS:
        raise RuntimeError(
            f"Refusing to initialize phonon calculation: estimated supercell has "
            f"{estimated_atoms} atoms (> {MAX_PHONON_SUPERCELL_ATOMS})"
        )
    
    # Write original POSCAR
    poscar_path = phon_dir / 'POSCAR'
    poscar = Poscar(structure)
    poscar.write_file(str(poscar_path))
    
    # Convert to phonopy atoms format
    unitcell = pymatgen_to_phonopy_atoms(structure)
    
    # Create Phonopy instance
    # Use primitive_matrix='auto' to let phonopy find primitive cell
    phonon = Phonopy(
        unitcell,
        supercell_matrix=supercell_matrix,
        primitive_matrix='auto',
        factor=VaspToTHz  # VASP unit conversion factor (eV/Å²/amu -> THz)
    )
    
    # Generate displacements (default 0.01)
    phonon.generate_displacements(distance=0.01)
    
    # Save phonopy displacement info
    phonon.save(filename=str(phon_dir / 'phonopy_disp.yaml'))
    
    # Get supercells with displacements
    supercells = phonon.supercells_with_displacements
    n_displacements = len(supercells)
    
    if n_displacements == 0:
        raise RuntimeError("phonopy did not generate any displaced structures")
    
    # Write displaced structures as POSCAR-001, POSCAR-002, ...
    for i, supercell in enumerate(supercells, start=1):
        disp_poscar_path = phon_dir / f'POSCAR-{i:03d}'
        
        # Convert phonopy atoms to pymatgen Structure
        lattice = supercell.cell
        positions = supercell.scaled_positions
        symbols = supercell.symbols
        
        disp_structure = Structure(
            lattice=lattice,
            species=symbols,
            coords=positions,
            coords_are_cartesian=False
        )
        
        # Write POSCAR
        disp_poscar = Poscar(disp_structure)
        disp_poscar.write_file(str(disp_poscar_path))
    
    return n_displacements, supercell_matrix


class WorkflowDatabase:
    """JSON-based database for tracking phonon job states."""
    
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
    
    def add_structure(self, struct_id, comp_name, phon_dir, n_displacements):
        """Add a new structure to track."""
        self.data['structures'][struct_id] = {
            'composition': comp_name,
            'state': 'PENDING',
            'phon_dir': str(phon_dir),
            'n_displacements': n_displacements,
            'displacement_jobs': {},
            'completed_displacements': 0,
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
        """Count displacement jobs currently running."""
        count = 0
        for sdata in self.data['structures'].values():
            if sdata['state'] == 'PHON_RUNNING':
                for disp_id, job_info in sdata.get('displacement_jobs', {}).items():
                    if job_info.get('status') == 'RUNNING':
                        count += 1
        return count


class PhononWorkflowManager:
    """Manages phonon calculations using phonopy finite displacement method."""
    
    def __init__(self, db_path, max_concurrent=20, check_interval=60, supercell_matrix=None, auto_dim=False):
        self.db = WorkflowDatabase(db_path)
        self.max_concurrent = max_concurrent
        self.check_interval = check_interval
        self.supercell_matrix = supercell_matrix
        self.auto_dim = auto_dim
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
    
    def create_vasp_inputs(self, structure, disp_dir):
        """Create VASP input files for displacement calculation."""
        disp_dir = Path(disp_dir)
        disp_dir.mkdir(parents=True, exist_ok=True)
        
        # High-precision static calculation settings aligned with refined relaxation
        vis = MPStaticSet(
    structure,
    user_incar_settings={
                'PREC': 'Accurate',
                'ALGO': 'Normal',
                'ADDGRID': True,

                'EDIFF': 1e-8,
                'NELM': 100,

                'IBRION': -1,
                'NSW': 0,

                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ISPIN': 1,

                'LASPH': True,
                'LREAL': False,
                'ISYM': 0,
                'NPAR': 4,

                'LWAVE': False,
                'LCHARG': False,
                'LAECHG': False,
                'LORBIT': 0,
            },
            user_kpoints_settings={'reciprocal_density': 250}
        )
        self.ensure_pmg_vasp_psp_dir()
        vis.write_input(disp_dir)
    
    def create_displacement_job_script(self, struct_id, disp_id, disp_dir):
        """Create SLURM job script for a single displacement calculation."""
        disp_dir = Path(disp_dir).resolve()
        script_path = disp_dir / 'job.sh'
        psp_dir = self.ensure_pmg_vasp_psp_dir()
        slurm_partitions = _slurm_env_value('VASP_JOB_PARTITIONS', 'Apus,Orion,Nebula')
        slurm_constraint = _slurm_env_value('VASP_JOB_CONSTRAINT', '')
        slurm_exclude = _slurm_env_value('VASP_JOB_EXCLUDE_NODES', DEFAULT_SLURM_EXCLUDE_NODES)
        slurm_time = _slurm_env_value('VASP_JOB_TIME', DEFAULT_VASP_JOB_TIME)
        vasp_cmd = _slurm_env_value('VASP_CMD', DEFAULT_VASP_CMD)
        constraint_line = f"#SBATCH --constraint={slurm_constraint}\n" if slurm_constraint else ""
        exclude_line = f"#SBATCH --exclude={slurm_exclude}\n" if slurm_exclude else ""
        
        script = f"""#!/bin/bash
#SBATCH --job-name={struct_id}_d{disp_id}
#SBATCH --partition={slurm_partitions}
{constraint_line}{exclude_line}#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time={slurm_time}
#SBATCH --output={disp_dir}/vasp_%j.out
#SBATCH --error={disp_dir}/vasp_%j.err

# Load modules
module purge
module load intel/mkl/2024.0 intel/2024 intel-mpi/2021.11
ulimit -s unlimited

# Set environment
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
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
echo "Phonon Displacement: {struct_id} - {disp_id}"
echo "========================================"
echo "Start time: $(date)"
echo ""

cd {disp_dir}

if [ ! -f "POSCAR" ] || [ ! -s "POSCAR" ]; then
    echo "ERROR: POSCAR not found"
    echo "VASP calculation failed" > VASP_FAILED
    exit 1
fi

echo "Running VASP static calculation for displacement {disp_id}..."

$VASP_CMD

EXIT_CODE=$?
echo "VASP exit code: $EXIT_CODE"

if [ $EXIT_CODE -ne 0 ]; then
    echo "VASP calculation failed" > VASP_FAILED
    echo "ERROR: VASP failed with exit code $EXIT_CODE"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

if [ ! -f "vasprun.xml" ] || [ ! -s "vasprun.xml" ]; then
    echo "VASP calculation failed" > VASP_FAILED
    echo "ERROR: vasprun.xml missing or empty"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

# Verify vasprun.xml is valid
if ! grep -q "</modeling>" vasprun.xml; then
    echo "VASP calculation failed" > VASP_FAILED
    echo "ERROR: vasprun.xml is incomplete"
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR 2>/dev/null
    exit 1
fi

# Check if calculation completed (convergence will be verified by Python)
echo "Displacement {disp_id} completed successfully"
touch VASP_DONE

# Cleanup
rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null

echo ""
echo "========================================"
echo "Displacement Complete"
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
    
    def check_displacement_status(self, disp_dir):
        """Check displacement calculation status with electronic convergence verification."""
        disp_dir = Path(disp_dir)

        if (disp_dir / 'VASP_FAILED').exists():
            return 'FAILED'

        vasprun_path = disp_dir / 'vasprun.xml'
        has_done_marker = (disp_dir / 'VASP_DONE').exists()

        # Primary path: explicit completion marker from job script.
        # Fallback path: valid, electronically converged vasprun.xml (marker may be missing
        # if a wrapper script got interrupted after VASP actually completed).
        if has_done_marker or (vasprun_path.exists() and vasprun_path.stat().st_size > 0):
            if not (vasprun_path.exists() and vasprun_path.stat().st_size > 0):
                return 'FAILED'

            try:
                vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
                if not vr.converged_electronic:
                    return 'FAILED'
                return 'DONE'
            except Exception:
                return 'FAILED'

        else:
            return 'UNKNOWN'
    
    def initialize_phonon_calculation(self, struct_id, structure):
        """Initialize phonon calculation for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        phon_dir = Path(sdata['phon_dir'])
        
        print(f"  Initializing phonon calculation: {struct_id}")
        
        try:
            # Generate supercell and displaced structures with phonopy Python API
            n_displacements, supercell_matrix_used = generate_supercell_with_phonopy(
                structure, phon_dir, 
                supercell_matrix=self.supercell_matrix,
                auto_dim=self.auto_dim
            )
            
            print(f"    Generated {n_displacements} displaced structures")
            print(f"    Supercell matrix used: {np.diag(supercell_matrix_used).astype(int)}")
            
            # Update database
            self.db.update_state(struct_id, 'INITIALIZED', 
                               n_displacements=n_displacements,
                               supercell_matrix=supercell_matrix_used.tolist())
            
            # Create VASP inputs for each displacement
            for i in range(1, n_displacements + 1):
                disp_id = f"{i:03d}"
                disp_poscar = phon_dir / f'POSCAR-{disp_id}'
                disp_dir = phon_dir / disp_id
                
                if not disp_poscar.exists():
                    print(f"    Warning: {disp_poscar} not found")
                    continue
                
                # Load displaced structure
                disp_structure = Structure.from_file(str(disp_poscar))
                
                # Create VASP inputs
                self.create_vasp_inputs(disp_structure, disp_dir)
                
                # Create job script
                self.create_displacement_job_script(struct_id, disp_id, disp_dir)
            
            print(f"    Created VASP inputs for all displacements")
            return True
            
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'INIT_FAILED', error=str(e))
            return False
    
    def submit_displacement_jobs(self, struct_id, max_to_submit):
        """Submit displacement jobs for a structure."""
        sdata = self.db.get_structure(struct_id)
        if not sdata or sdata['state'] not in ['INITIALIZED', 'PHON_RUNNING']:
            return 0
        
        phon_dir = Path(sdata['phon_dir'])
        n_displacements = sdata['n_displacements']
        displacement_jobs = sdata.get('displacement_jobs', {})
        
        submitted = 0
        
        for i in range(1, n_displacements + 1):
            if submitted >= max_to_submit:
                break
            
            disp_id = f"{i:03d}"
            
            # Skip if already submitted
            if disp_id in displacement_jobs:
                continue
            
            disp_dir = phon_dir / disp_id
            job_script = disp_dir / 'job.sh'
            
            if not job_script.exists():
                continue
            
            try:
                job_id = self.submit_job(job_script)
                displacement_jobs[disp_id] = {
                    'job_id': job_id,
                    'status': 'RUNNING',
                    'submitted_at': datetime.now().isoformat()
                }
                submitted += 1
                
                if submitted == 1:
                    print(f"  {struct_id}: Submitting displacement jobs...")
                
                if submitted % 10 == 0:
                    print(f"    Submitted {submitted} jobs...")
                
            except Exception as e:
                print(f"    Error submitting {disp_id}: {e}")
        
        if submitted > 0:
            self.db.update_state(struct_id, 'PHON_RUNNING', 
                               displacement_jobs=displacement_jobs)
            print(f"    Total submitted: {submitted} displacement jobs")
        
        return submitted
    
    def update_structure_status(self, struct_id):
        """Update status of a structure's displacement calculations."""
        sdata = self.db.get_structure(struct_id)
        if not sdata or sdata['state'] != 'PHON_RUNNING':
            return
        
        phon_dir = Path(sdata['phon_dir'])
        n_displacements = sdata['n_displacements']
        displacement_jobs = sdata.get('displacement_jobs', {})
        
        completed = 0
        failed = 0
        running = 0
        failed_displacements = []
        
        for i in range(1, n_displacements + 1):
            disp_id = f"{i:03d}"
            disp_dir = phon_dir / disp_id
            
            if disp_id not in displacement_jobs:
                continue
            
            job_info = displacement_jobs[disp_id]
            
            # Check job status
            if job_info['status'] == 'RUNNING':
                job_status = self.check_job_status(job_info['job_id'])
                
                if job_status == 'NOTFOUND':
                    # Job finished, check results
                    local_status = self.check_displacement_status(disp_dir)
                    
                    if local_status == 'DONE':
                        job_info['status'] = 'DONE'
                        completed += 1
                    elif local_status == 'FAILED':
                        job_info['status'] = 'FAILED'
                        job_info['error'] = 'Electronic SCF not converged or vasprun.xml invalid'
                        failed += 1
                        failed_displacements.append(disp_id)
                    else:
                        job_info['status'] = 'FAILED'
                        job_info['error'] = 'Job terminated without completion marker'
                        failed += 1
                        failed_displacements.append(disp_id)
                else:
                    running += 1
            elif job_info['status'] == 'DONE':
                completed += 1
            elif job_info['status'] == 'FAILED':
                failed += 1
                if disp_id not in failed_displacements:
                    failed_displacements.append(disp_id)
        
        # Update database
        self.db.update_state(struct_id, 'PHON_RUNNING',
                           displacement_jobs=displacement_jobs,
                           completed_displacements=completed)
        
        # Check if all done
        if completed == n_displacements:
            self.db.update_state(struct_id, 'PHON_DONE')
            print(f"  {struct_id}: All {n_displacements} displacements completed (all electronically converged)")
        elif failed > 0 and completed + failed == n_displacements:
            error_msg = f'{failed}/{n_displacements} displacements failed: {", ".join(failed_displacements)}'
            self.db.update_state(struct_id, 'PHON_FAILED', error=error_msg)
            print(f"  {struct_id}: FAILED ({failed}/{n_displacements} displacements failed)")
            print(f"    Failed: {', '.join(failed_displacements)}")
    
    def monitor_and_submit(self, structures_dict):
        """Main monitoring loop."""
        print("\n" + "="*70)
        print("Starting phonon workflow monitoring...")
        print(f"Max concurrent displacement jobs: {self.max_concurrent}")
        print(f"Check interval: {self.check_interval}s")
        if self.auto_dim:
            print(f"Supercell dimensions: AUTO (based on k-point mesh)")
        elif self.supercell_matrix is not None:
            print(f"Supercell matrix: {np.diag(self.supercell_matrix).astype(int)}")
        else:
            print(f"Supercell dimensions: [2, 2, 2] (default)")
        print("="*70 + "\n")
        sys.stdout.flush()
        
        while True:
            print(f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Checking status...")
            sys.stdout.flush()
            
            # Update status of running structures
            for struct_id in list(self.db.data['structures'].keys()):
                sdata = self.db.get_structure(struct_id)
                if sdata['state'] == 'PHON_RUNNING':
                    self.update_structure_status(struct_id)
            
            running_count = self.db.get_running_count()
            print(f"Currently running displacement jobs: {running_count}/{self.max_concurrent}")
            
            # Initialize pending structures
            for struct_id in list(self.db.data['structures'].keys()):
                sdata = self.db.get_structure(struct_id)
                if sdata['state'] == 'PENDING':
                    structure = structures_dict.get(struct_id)
                    if structure:
                        self.initialize_phonon_calculation(struct_id, structure)
            
            # Submit displacement jobs
            for struct_id in list(self.db.data['structures'].keys()):
                if running_count >= self.max_concurrent:
                    break
                
                sdata = self.db.get_structure(struct_id)
                if sdata['state'] in ['INITIALIZED', 'PHON_RUNNING']:
                    max_to_submit = self.max_concurrent - running_count
                    submitted = self.submit_displacement_jobs(struct_id, max_to_submit)
                    running_count += submitted
            
            # Print statistics
            stats = {}
            for sdata in self.db.data['structures'].values():
                state = sdata['state']
                stats[state] = stats.get(state, 0) + 1
            
            print("\nStatistics:")
            for state, count in sorted(stats.items()):
                print(f"  {state}: {count}")
            sys.stdout.flush()
            
            # Check if all done
            pending = len(self.db.get_by_state('PENDING'))
            initialized = len(self.db.get_by_state('INITIALIZED'))
            running = len(self.db.get_by_state('PHON_RUNNING'))
            
            if running_count == 0 and pending == 0 and initialized == 0 and running == 0:
                completed = len(self.db.get_by_state('PHON_DONE'))
                failed = (len(self.db.get_by_state('INIT_FAILED')) + 
                         len(self.db.get_by_state('PHON_FAILED')))
                total = len(self.db.data['structures'])
                
                if completed + failed >= total:
                    print("\n" + "="*70)
                    print("All phonon calculations completed!")
                    print(f"Completed: {completed}/{total}")
                    print(f"Failed: {failed}/{total}")
                    print("="*70)
                    sys.stdout.flush()
                    break
            
            print(f"\nSleeping for {self.check_interval}s...")
            sys.stdout.flush()
            time.sleep(self.check_interval)
    
    def initialize_structures(self, refine_jobs_dir, output_dir, structure_ids):
        """Initialize structures from refined workflow."""
        refine_jobs_dir = Path(refine_jobs_dir)
        output_dir = Path(output_dir)
        
        print("="*70)
        print("Phonon Workflow Initialization")
        print("="*70)
        print(f"Refined jobs dir: {refine_jobs_dir}")
        print(f"Output directory: {output_dir}")
        print(f"Max concurrent displacement jobs: {self.max_concurrent}")
        if self.auto_dim:
            print(f"Supercell dimensions: AUTO (based on k-point mesh)")
        elif self.supercell_matrix is not None:
            print(f"Supercell matrix: {np.diag(self.supercell_matrix).astype(int)}")
        else:
            print(f"Supercell dimensions: [2, 2, 2] (default)")
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
            print(f"ERROR: No structure IDs specified. Please provide --structure-ids")
            return {}
        
        print(f"Found {len(target_structures)} structures to process\n")
        
        # Load structures
        structures_dict = {}
        loaded = 0
        failed = 0
        
        for struct_id, sdata in target_structures.items():
            # Check if relaxation is done
            if sdata['state'] != 'RELAX_DONE':
                print(f"  {struct_id}: Skipping (state: {sdata['state']})")
                failed += 1
                continue
            
            relax_dir = Path(sdata['relax_dir'])
            contcar_path = relax_dir / 'CONTCAR'
            
            if not contcar_path.exists():
                print(f"  Warning: CONTCAR not found for {struct_id}")
                failed += 1
                continue
            
            structure = load_structure_from_contcar(contcar_path)
            if structure is None:
                print(f"  Warning: Could not load structure for {struct_id}")
                failed += 1
                continue
            
            structures_dict[struct_id] = structure
            
            if struct_id not in self.db.data['structures']:
                comp_name = sdata['composition']
                phon_dir = output_dir / comp_name / struct_id / 'PHON'
                self.db.add_structure(struct_id, comp_name, phon_dir, 0)
            
            loaded += 1
        
        print(f"\nSuccessfully loaded: {loaded}")
        if failed > 0:
            print(f"Failed to load or not RELAX_DONE: {failed}")
        
        self.db.data['config'] = {
            'max_concurrent': self.max_concurrent,
            'refine_jobs_dir': str(refine_jobs_dir),
            'output_dir': str(output_dir),
            'supercell_matrix': self.supercell_matrix.tolist() if self.supercell_matrix is not None else None,
            'auto_dim': self.auto_dim
        }
        self.db.save()
        
        return structures_dict


def main():
    parser = argparse.ArgumentParser(
        description="Phonon Workflow Manager - Phonopy Finite Displacement Method"
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
        default='./PHONON_JOBS',
        help="Output directory for phonon jobs (default: ./PHONON_JOBS)"
    )
    parser.add_argument(
        '--structure-ids',
        type=str,
        nargs='+',
        required=True,
        help="Structure IDs to process (must be RELAX_DONE)"
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
        default=20,
        help="Max concurrent displacement jobs (default: 20)"
    )
    parser.add_argument(
        '--check-interval',
        type=int,
        default=60,
        help="Status check interval in seconds (default: 60)"
    )
    parser.add_argument(
        '--supercell-dim',
        type=int,
        nargs=3,
        default=None,
        metavar=('N1', 'N2', 'N3'),
        help="Manually specify supercell dimensions (e.g., --supercell-dim 2 2 2). "
             "If not specified, auto-calculates from k-point mesh (default behavior)"
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
    
    # Prepare supercell matrix
    if args.supercell_dim is not None:
        # User explicitly specified dimensions - use them
        supercell_matrix = np.diag(args.supercell_dim)
        auto_dim = False
    else:
        # Default: auto-calculate from k-point mesh
        supercell_matrix = None
        auto_dim = True
    
    manager = PhononWorkflowManager(
        db_path=db_path,
        max_concurrent=args.max_concurrent,
        check_interval=args.check_interval,
        supercell_matrix=supercell_matrix,
        auto_dim=auto_dim
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
