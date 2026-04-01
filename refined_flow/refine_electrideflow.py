#!/usr/bin/env python3
"""
Refined Electride Workflow Manager - High-precision 3-step relaxation for electride candidates

This workflow runs on structures identified as electrides after filtering by filter_comb_db.py.
Uses finer k-point grid and consecutive relaxation steps for more accurate structures.

Key differences from workflow_manager.py:
- Reads electride candidates from electride_candidates.db (PyXtal/ASE database)
- Uses reciprocal_density=250 (vs 64) for all jobs
- Uses primitive cells from PyXtal symmetrization (faster VASP calculations)
- RELAX: 3 consecutive steps with progressive parameters (NELM=60, NSW=100 each)
  * Step 1: ISIF=2, EDIFFG=-0.02, IBRION=2, POTIM=0.3, EDIFF=1e-7
  * Step 2: ISIF=3, EDIFFG=-0.01, IBRION=2, POTIM=0.2, EDIFF=1e-8
  * Step 3: ISIF=3, EDIFFG=-0.005, IBRION=2, POTIM=0.1, EDIFF=1e-8
  * All steps: PREC=Accurate, LREAL=.FALSE., LASPH=.TRUE., NPAR=4
- SPE workflow (SC/PARCHG/ELF) is skipped - only refined relaxation
- SLURM time limit: 4 hours (vs 20 minutes)
- Output to REFINE_VASP_JOBS folder

Input:
  - electride_candidates.db or .csv: Filtered electride candidates (structure IDs + compositions)
  - VASP_JOBS/: Original workflow directory with already-relaxed and symmetrized POSCARs
    Directory structure: VASP_JOBS/{composition}/{struct_id}/Relax/POSCAR

Output: High-quality refined relaxed structures (CONTCAR) for further analysis

Note: Structures are loaded from original VASP_JOBS POSCARs (not re-generated from CIFs) 
to ensure consistency with the original workflow's symmetrization.
"""

import os
import sys
import json
import time
import argparse
import warnings
import subprocess
import shutil
import pandas as pd
from pathlib import Path
from datetime import datetime

from pymatgen.core import Structure, SETTINGS as PMG_SETTINGS
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.ase import AseAtomsAdaptor

try:
    from pyxtal import pyxtal
    from pyxtal.db import database_topology
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


def load_electride_candidates(input_path, vasp_jobs_dir):
    """
    Load electride candidates from CSV or database and extract structures from original VASP_JOBS.
    
    This function:
    1. Reads structure IDs and compositions from electride_candidates.csv or .db
    2. Filters for high-symmetry structures (space group > 15)
    3. Loads relaxed CONTCARs from original workflow
    4. Applies PyXtal symmetrization
    
    Note: We load CONTCARs from VASP_JOBS/{comp}/{struct_id}/Relax/CONTCAR and apply
    the same PyXtal symmetrization as the original workflow. This uses the relaxed
    structures as starting point for high-precision refinement.
    
    Filtering: Only processes orthorhombic and higher symmetry systems (space group > 15),
    excluding triclinic (1-2) and monoclinic (3-15) space groups.
    
    Args:
        input_path: Path to electride_candidates.csv or electride_candidates.db
        vasp_jobs_dir: Path to original VASP_JOBS directory
    
    Returns:
        dict: {structure_id: structure} mapping with PyXtal-symmetrized structures
    """
    vasp_jobs_dir = Path(vasp_jobs_dir)
    
    if not vasp_jobs_dir.exists():
        print(f"ERROR: VASP_JOBS directory not found: {vasp_jobs_dir}")
        return {}
    
    input_path = Path(input_path)
    candidates = []  # List of (structure_id, composition) tuples
    
    # Load structure IDs and compositions from CSV or database
    if input_path.suffix == '.csv':
        print(f"Loading structure IDs from CSV: {input_path}")
        try:
            df = pd.read_csv(input_path)
            
            # Filter for high-symmetry structures (space group > 15)
            if 'spacegroup' in df.columns:
                original_count = len(df)
                df = df[df['spacegroup'] > 15]
                filtered_count = original_count - len(df)
                print(f"  Filtered {filtered_count} low-symmetry structures (space group <= 15)")
                print(f"  Focusing on {len(df)} high-symmetry structures (space group > 15)")
            
            for _, row in df.iterrows():
                struct_id = row['formula']
                composition = row.get('composition', struct_id.rsplit('_s', 1)[0])
                candidates.append((struct_id, composition))
            print(f"  Found {len(candidates)} structure IDs in CSV")
        except Exception as e:
            print(f"ERROR: Could not read CSV: {e}")
            return {}
    
    elif input_path.suffix == '.db':
        print(f"Loading structure IDs from database: {input_path}")
        if not PYXTAL_AVAILABLE:
            print("ERROR: PyXtal not available. Cannot load database.")
            return {}
        
        try:
            db = database_topology(str(input_path))
            total_count = 0
            filtered_count = 0
            
            for row in db.db.select():
                try:
                    total_count += 1
                    struct_id = row.structure_id
                    composition = getattr(row, 'composition', struct_id.rsplit('_s', 1)[0])
                    
                    # Filter for high-symmetry structures (space group > 15)
                    space_group_num = getattr(row, 'space_group_number', None)
                    if space_group_num is not None and space_group_num <= 15:
                        filtered_count += 1
                        continue
                    
                    candidates.append((struct_id, composition))
                except AttributeError:
                    continue
            
            if filtered_count > 0:
                print(f"  Filtered {filtered_count} low-symmetry structures (space group <= 15)")
                print(f"  Focusing on {len(candidates)} high-symmetry structures (space group > 15)")
            print(f"  Found {len(candidates)} structure IDs in database")
        except Exception as e:
            print(f"ERROR: Could not read database: {e}")
            return {}
    else:
        print(f"ERROR: Unknown file type: {input_path}")
        print("  Supported: .csv or .db")
        return {}
    
    # Load structures from original VASP_JOBS directory
    print(f"\nLoading structures from: {vasp_jobs_dir}")
    print("  (Using relaxed CONTCARs from original workflow with PyXtal symmetrization)\n")
    
    structures_dict = {}
    loaded_count = 0
    failed_count = 0
    symmetrized_count = 0
    
    for struct_id, composition in candidates:
        # VASP_JOBS directory structure: VASP_JOBS/{comp}/{struct_id}/Relax/CONTCAR
        contcar_path = vasp_jobs_dir / composition / struct_id / 'Relax' / 'CONTCAR'
        
        if contcar_path.exists():
            try:
                structure = Structure.from_file(str(contcar_path))
                
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
                            symmetrized_count += 1
                            break
                        except Exception:
                            continue
                    
                    if not symmetrized:
                        pass
                
                structures_dict[struct_id] = structure
                loaded_count += 1
                
                if loaded_count % 100 == 0:
                    print(f"  Loaded {loaded_count} structures...")
            except Exception as e:
                print(f"  Warning: Could not load {struct_id} from CONTCAR: {e}")
                failed_count += 1
        else:
            print(f"  Warning: CONTCAR not found for {struct_id}: {contcar_path}")
            failed_count += 1
    
    print(f"\nSuccessfully loaded: {loaded_count}")
    if PYXTAL_AVAILABLE:
        print(f"  PyXtal symmetrized: {symmetrized_count}")
        print(f"  Used original structure: {loaded_count - symmetrized_count}")
    if failed_count > 0:
        print(f"Failed to load: {failed_count}")
    
    return structures_dict


class WorkflowDatabase:
    """
    Simple JSON-based database for tracking job states.
    
    Job Status Design for 3-Step Relaxation:
    ----------------------------------------
    The refined workflow uses a single SLURM job for all 3 relaxation steps,
    so we use single-level status tracking:
    
    - PENDING: Not yet submitted
    - RELAX_RUNNING: 3-step job is running
    - RELAX_DONE: All 3 steps completed successfully (final state)
    - RELAX_TMOUT: Step 3 timed out but produced CONTCAR (usable, final state)
    - RELAX_FAILED: Job failed at any step
    
    Failure debugging: POSCAR-1, POSCAR-2, POSCAR-3 files in job directory
    indicate which step was attempted, allowing identification of failure point.
    """
    
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
            'relax_dir': str(base_dir / struct_id / 'Relax'),
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
        """Count structures currently running."""
        running_states = ['RELAX_RUNNING']
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
    """Manages refined VASP job submission and monitoring for electride candidates."""
    
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
    
    def read_structure_from_contcar(self, contcar_path):
        """Read structure from CONTCAR file."""
        try:
            return Structure.from_file(str(contcar_path))
        except Exception as e:
            print(f"  Warning: Could not read CONTCAR: {e}")
            return None
    
    def create_vasp_inputs(self, structure, job_dir, job_type='relax', step=1):
        """
        Create VASP input files for refined 3-step relaxation.
        
        Progressive relaxation parameters:
        - Step 1: ISIF=2, EDIFFG=-0.02, IBRION=2, POTIM=0.3, EDIFF=1e-7
        - Step 2: ISIF=3, EDIFFG=-0.01, IBRION=2, POTIM=0.2, EDIFF=1e-8
        - Step 3: ISIF=3, EDIFFG=-0.005, IBRION=2, POTIM=0.1, EDIFF=1e-8
        - All steps: PREC=Accurate, LREAL=.FALSE., LASPH=.TRUE., NPAR=4
        
        Key differences from workflow_manager.py:
        - Higher k-point density (250 vs 64)
        - More electronic/ionic steps (NELM=60, NSW=100)
        - Progressive tightening of convergence criteria
        
        Note: Structure is already PyXtal-symmetrized in load_electride_candidates(),
        so we use it directly here.
        """
        job_dir = Path(job_dir)
        job_dir.mkdir(parents=True, exist_ok=True)
        
        if job_type != 'relax':
            raise ValueError(f"Only 'relax' job_type supported in refined workflow, got: {job_type}")
        
        # Progressive relaxation parameters for each step
        if step == 1:
            ediffg = -0.02
            isif = 2
            ibrion = 2
            potim = 0.3
            ediff = 1e-7
        elif step == 2:
            ediffg = -0.01
            isif = 3
            ibrion = 2
            potim = 0.2
            ediff = 1e-8
        elif step == 3:
            ediffg = -0.005
            isif = 3
            ibrion = 2
            potim = 0.1
            ediff = 1e-8
        else:
            raise ValueError(f"Invalid step number: {step}. Must be 1, 2, or 3.")
        
        vis = MPRelaxSet(structure, 
            user_incar_settings={
                'PREC': 'Accurate',
                'ALGO': 'normal',
                'ADDGRID': True,
                'EDIFF': ediff,
                'EDIFFG': ediffg,
                'IBRION': ibrion,
                'ISIF': isif,
                'NELM': 60,
                'NSW': 100,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ISPIN': 1,
                'POTIM': potim,
                'LWAVE': False,
                'LCHARG': False,
                'LAECHG': False,
                'LASPH': True,
                'LREAL': False,
                'NPAR': 4,
                'LORBIT': 0,
            },
            user_kpoints_settings={'reciprocal_density': 250}
        )
        
        self.ensure_pmg_vasp_psp_dir()
        vis.write_input(job_dir)
        return job_dir
    
    def create_slurm_script(self, job_dir, job_name, job_type='relax'):
        """
        Create SLURM submission script for 3-step progressive relaxation.
        
        Progressive Relaxation Strategy:
        --------------------------------
        - Step 1: ISIF=2, EDIFFG=-0.02, IBRION=2, POTIM=0.3 (quick ionic, CG)
        - Step 2: ISIF=3, EDIFFG=-0.01, IBRION=2, POTIM=0.2 (full relax)
        - Step 3: ISIF=3, EDIFFG=-0.005, IBRION=2, POTIM=0.1 (tight)
        - Step 3 adaptive retries: same INCAR-3 settings, up to 3 reruns
        
        Timeout/Failure Handling:
        -------------------------
        1. Before each step: Save POSCAR-{1,2,3} for debugging
        2. Steps 1-2: Timeout OK if CONTCAR exists + electronic converged (continue)
        3. Step 3 requires "reached required accuracy" at least once
        4. After first accuracy hit, rerun and require ionic steps <= 3 to finish
        5. Step 3: Timeout → RELAX_TMOUT marker (CONTCAR usable)
        6. Any VASP failure → VASP_FAILED + RELAX_FAIL_REASON.txt
        7. Exit codes 140, 143 indicate SLURM timeout
        
        Markers created:
        - VASP_DONE: All 3 steps completed successfully
        - RELAX_TMOUT: Step 3 timed out with valid CONTCAR
        - VASP_FAILED: VASP error at any step
        - RELAX_FAIL_REASON.txt: explicit step/reason for monitor reporting
        """
        job_dir = Path(job_dir).resolve()
        script_path = job_dir / 'job.sh'
        psp_dir = self.ensure_pmg_vasp_psp_dir()
        
        script = f"""#!/bin/bash
#SBATCH --job-name={job_name}_{job_type}_refine
#SBATCH --partition=Apus,Orion,Nebula
#SBATCH --exclude=str-c88,str-c89,str-c90,str-c91,str-c92,str-c93,str-c94,str-c95,str-c96,str-c97
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --output={job_dir}/vasp_%j.out
#SBATCH --error={job_dir}/vasp_%j.err

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
        
        # Relax job type - only job type used in refined workflow
        if job_type == 'relax':
            script += f"""
# Run 3 consecutive VASP relaxation steps with timeout/failure handling
echo "Starting 3-step VASP relaxation"
echo "Working directory: $(pwd)"
echo "VASP command: $VASP_CMD"
echo "Start time: $(date)"
FAIL_REASON_FILE="RELAX_FAIL_REASON.txt"
rm -f VASP_DONE VASP_FAILED RELAX_TMOUT "$FAIL_REASON_FILE" 2>/dev/null

cleanup_intermediate_files() {{
    rm -f WAVECAR CHGCAR CHG WFULL TMPCAR AECCAR* 2>/dev/null
}}

cleanup_failed_files() {{
    rm -f CHGCAR CHG WAVECAR vasprun.xml WFULL AECCAR* TMPCAR OSZICAR OUTCAR 2>/dev/null
}}

record_failure() {{
    local step="$1"
    local reason="$2"
    echo "step=$step" > "$FAIL_REASON_FILE"
    echo "reason=$reason" >> "$FAIL_REASON_FILE"
    touch VASP_FAILED
}}

fail_job() {{
    local step="$1"
    local reason="$2"
    record_failure "$step" "$reason"
    exit 1
}}

# Relaxation step 1
echo ""
echo "========================================"
echo "Relaxation Step 1/3"
echo "  ISIF=2, EDIFFG=-0.02, IBRION=2, POTIM=0.3"
echo "========================================"

# Save initial POSCAR for debugging
cp POSCAR POSCAR-1

# Use INCAR-1 for step 1
cp INCAR-1 INCAR

$VASP_CMD

EXIT_CODE=$?
echo "Step 1 exit code: $EXIT_CODE"

# Check if SLURM timeout occurred
if [ $EXIT_CODE -eq 140 ] || [ $EXIT_CODE -eq 143 ]; then
    echo "WARNING: Step 1 timed out (exit code $EXIT_CODE)"
    if [ -f "CONTCAR" ] && [ -s "CONTCAR" ]; then
        echo "CONTCAR exists, proceeding to step 2 with partial relaxation"
    else
        echo "ERROR: Step 1 timed out without producing CONTCAR"
        cleanup_failed_files
        fail_job "step1" "timeout_without_contcar"
    fi
elif [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: Relaxation step 1 failed with exit code $EXIT_CODE"
    cleanup_failed_files
    fail_job "step1" "vasp_exit_code_$EXIT_CODE"
fi

# Verify CONTCAR from step 1
if [ ! -f "CONTCAR" ] || [ ! -s "CONTCAR" ]; then
    echo "ERROR: CONTCAR missing/empty after step 1"
    cleanup_failed_files
    fail_job "step1" "missing_contcar"
fi

# Check electronic convergence from OUTCAR
echo ""
echo "Checking electronic convergence..."
if [ ! -f "OUTCAR" ] || [ ! -s "OUTCAR" ]; then
    echo "ERROR: OUTCAR not found - cannot verify convergence"
    cleanup_failed_files
    fail_job "step1" "missing_outcar"
fi

if grep -q "aborting loop because EDIFF is reached" OUTCAR; then
    echo "  Electronic SCF converged in step 1"
else
    echo "ERROR: Electronic SCF did not converge in step 1"
    cleanup_failed_files
    fail_job "step1" "electronic_scf_not_converged"
fi

echo "Step 1 completed (exit code: $EXIT_CODE)"
echo "Copying CONTCAR -> POSCAR for step 2"
cp CONTCAR POSCAR

# Clean large intermediate files between steps
cleanup_intermediate_files

# Relaxation step 2
echo ""
echo "========================================"
echo "Relaxation Step 2/3"
echo "  ISIF=3, EDIFFG=-0.01, IBRION=2, POTIM=0.2"
echo "========================================"

# Save POSCAR for debugging (should be identical to step 1 CONTCAR)
cp POSCAR POSCAR-2

# Use INCAR-2 for step 2
cp INCAR-2 INCAR

$VASP_CMD

EXIT_CODE=$?
echo "Step 2 exit code: $EXIT_CODE"

# Check if SLURM timeout occurred
if [ $EXIT_CODE -eq 140 ] || [ $EXIT_CODE -eq 143 ]; then
    echo "WARNING: Step 2 timed out (exit code $EXIT_CODE)"
    if [ -f "CONTCAR" ] && [ -s "CONTCAR" ]; then
        echo "CONTCAR exists, proceeding to step 3 with partial relaxation"
    else
        echo "ERROR: Step 2 timed out without producing CONTCAR"
        cleanup_failed_files
        fail_job "step2" "timeout_without_contcar"
    fi
elif [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: Relaxation step 2 failed with exit code $EXIT_CODE"
    cleanup_failed_files
    fail_job "step2" "vasp_exit_code_$EXIT_CODE"
fi

# Verify CONTCAR from step 2
if [ ! -f "CONTCAR" ] || [ ! -s "CONTCAR" ]; then
    echo "ERROR: CONTCAR missing/empty after step 2"
    cleanup_failed_files
    fail_job "step2" "missing_contcar"
fi

# Check electronic convergence from OUTCAR
echo ""
echo "Checking electronic convergence..."
if [ ! -f "OUTCAR" ] || [ ! -s "OUTCAR" ]; then
    echo "ERROR: OUTCAR not found - cannot verify convergence"
    cleanup_failed_files
    fail_job "step2" "missing_outcar"
fi

if grep -q "aborting loop because EDIFF is reached" OUTCAR; then
    echo "  Electronic SCF converged in step 2"
else
    echo "ERROR: Electronic SCF did not converge in step 2"
    cleanup_failed_files
    fail_job "step2" "electronic_scf_not_converged"
fi

echo "Step 2 completed (exit code: $EXIT_CODE)"
echo "Copying CONTCAR -> POSCAR for step 3"
cp CONTCAR POSCAR

# Clean large intermediate files between steps
cleanup_intermediate_files

# Relaxation step 3 (adaptive reruns)
echo ""
echo "========================================"
echo "Relaxation Step 3/3 (Adaptive)"
echo "  ISIF=3, EDIFFG=-0.005, IBRION=2, POTIM=0.1"
echo "  Stop criteria:"
echo "    (1) OUTCAR contains 'reached required accuracy' at least once"
echo "    (2) In a follow-up run, ionic relax steps <= 3"
echo "  Maximum reruns with same settings: 3"
echo "========================================"

# Save POSCAR for debugging (should be identical to step 2 CONTCAR)
cp POSCAR POSCAR-3

# Use INCAR-3 for step 3 (final high-precision; reused for reruns)
cp INCAR-3 INCAR

STEP3_MAX_RERUNS=3
STEP3_TOTAL_ATTEMPTS=$((STEP3_MAX_RERUNS + 1))
STEP3_ATTEMPT=1
STEP3_REQUIRED_ACCURACY_HIT=0
STEP3_IONIC_STEP_THRESHOLD=3

get_ionic_steps() {{
    if [ ! -f "OSZICAR" ] || [ ! -s "OSZICAR" ]; then
        echo 0
        return
    fi
    awk '/^[[:space:]]*[0-9]+[[:space:]]+F=/{{if($1+0>max) max=$1+0}} END{{print max+0}}' OSZICAR
}}

while [ $STEP3_ATTEMPT -le $STEP3_TOTAL_ATTEMPTS ]; do
    echo ""
    echo "Running Step 3 attempt $STEP3_ATTEMPT/$STEP3_TOTAL_ATTEMPTS"
    $VASP_CMD

    EXIT_CODE=$?
    echo "Step 3 attempt $STEP3_ATTEMPT exit code: $EXIT_CODE"

    # Check if SLURM timeout occurred in final step
    if [ $EXIT_CODE -eq 140 ] || [ $EXIT_CODE -eq 143 ]; then
        echo "WARNING: Step 3 timed out (exit code $EXIT_CODE)"
        if [ -f "CONTCAR" ] && [ -s "CONTCAR" ]; then
            echo "CONTCAR exists - marking as RELAX_TMOUT"
            echo "Partial relaxation completed, may proceed to analysis"
            echo "Cleaning up large intermediate files..."
            rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR OSZICAR OUTCAR 2>/dev/null
            touch RELAX_TMOUT
            touch VASP_DONE
            exit 0
        else
            echo "ERROR: Step 3 timed out without producing CONTCAR"
            cleanup_failed_files
            fail_job "step3" "timeout_without_contcar"
        fi
    elif [ $EXIT_CODE -ne 0 ]; then
        echo "ERROR: Relaxation step 3 failed with exit code $EXIT_CODE"
        cleanup_failed_files
        fail_job "step3" "vasp_exit_code_$EXIT_CODE"
    fi

    if [ ! -f "CONTCAR" ] || [ ! -s "CONTCAR" ]; then
        echo "ERROR: CONTCAR missing/empty after step 3 attempt $STEP3_ATTEMPT"
        cleanup_failed_files
        fail_job "step3" "missing_contcar"
    fi

    if [ ! -f "OUTCAR" ] || [ ! -s "OUTCAR" ]; then
        echo "ERROR: OUTCAR not found after step 3 attempt $STEP3_ATTEMPT"
        cleanup_failed_files
        fail_job "step3" "missing_outcar"
    fi

    STEP3_IONIC_STEPS=$(get_ionic_steps)
    if [ "$STEP3_IONIC_STEPS" -le 0 ]; then
        echo "ERROR: Could not determine ionic relaxation steps from OSZICAR"
        cleanup_failed_files
        fail_job "step3" "ionic_step_count_unavailable"
    fi

    HAS_REQUIRED_ACCURACY=0
    if grep -q "reached required accuracy" OUTCAR; then
        HAS_REQUIRED_ACCURACY=1
    fi

    if [ $STEP3_REQUIRED_ACCURACY_HIT -eq 0 ]; then
        if [ $HAS_REQUIRED_ACCURACY -eq 1 ]; then
            STEP3_REQUIRED_ACCURACY_HIT=1
            echo "  First required-accuracy hit detected."
            echo "  Forcing one more rerun to apply ionic-step criterion (<= $STEP3_IONIC_STEP_THRESHOLD)."
        else
            echo "  Required accuracy not reached yet (ionic steps: $STEP3_IONIC_STEPS)."
        fi
    else
        if [ "$STEP3_IONIC_STEPS" -le "$STEP3_IONIC_STEP_THRESHOLD" ]; then
            echo "  Post-accuracy ionic-step criterion met: $STEP3_IONIC_STEPS <= $STEP3_IONIC_STEP_THRESHOLD"
            break
        fi
        echo "  Post-accuracy ionic-step criterion not met: $STEP3_IONIC_STEPS > $STEP3_IONIC_STEP_THRESHOLD"
    fi

    if [ $STEP3_ATTEMPT -ge $STEP3_TOTAL_ATTEMPTS ]; then
        cleanup_failed_files
        if [ $STEP3_REQUIRED_ACCURACY_HIT -eq 0 ]; then
            echo "ERROR: Step 3 never reached required accuracy after $STEP3_TOTAL_ATTEMPTS attempts"
            fail_job "step3" "required_accuracy_not_reached_after_$STEP3_TOTAL_ATTEMPTS_attempts"
        else
            echo "ERROR: Step 3 reached required accuracy, but ionic steps stayed > $STEP3_IONIC_STEP_THRESHOLD after $STEP3_TOTAL_ATTEMPTS attempts"
            fail_job "step3" "post_accuracy_ionic_steps_gt_$STEP3_IONIC_STEP_THRESHOLD_after_$STEP3_TOTAL_ATTEMPTS_attempts"
        fi
    fi

    echo "Step 3 rerun required..."
    echo "Copying CONTCAR -> POSCAR before rerun"
    cp CONTCAR POSCAR
    cleanup_intermediate_files
    STEP3_ATTEMPT=$((STEP3_ATTEMPT + 1))
done

echo ""
echo "========================================"
echo "All 3 Relaxation Steps Completed Successfully"
echo "========================================"
echo "End time: $(date)"
echo "Final structure saved in CONTCAR"

# Clean up large unnecessary files to save disk space
echo ""
echo "Cleaning up large intermediate files..."
rm -f CHGCAR CHG WAVECAR WFULL AECCAR* TMPCAR 2>/dev/null
echo "Cleanup complete - kept: POSCAR, CONTCAR, INCAR, KPOINTS, POTCAR, OUTCAR, OSZICAR, vasprun.xml, POSCAR-*"

touch VASP_DONE
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
        """Check SLURM job status."""
        result = subprocess.run(
            ['squeue', '-j', job_id, '-h', '-o', '%T'],
            capture_output=True,
            text=True
        )
        
        if result.stdout.strip():
            return 'RUNNING'
        else:
            return 'NOTFOUND'
    
    def check_local_status(self, job_dir):
        """Check local directory for completion markers."""
        job_dir = Path(job_dir)
        if (job_dir / 'VASP_DONE').exists():
            # Check if it's a timeout completion
            if (job_dir / 'RELAX_TMOUT').exists():
                return 'TIMEOUT'
            return 'DONE'
        elif (job_dir / 'VASP_FAILED').exists():
            return 'FAILED'
        else:
            return 'UNKNOWN'

    def read_failure_context(self, job_dir):
        """Read failure context written by the SLURM script."""
        reason_file = Path(job_dir) / 'RELAX_FAIL_REASON.txt'
        if not reason_file.exists():
            return None, None

        step = None
        reason = None
        try:
            with open(reason_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('step='):
                        step = line.split('=', 1)[1].strip()
                    elif line.startswith('reason='):
                        reason = line.split('=', 1)[1].strip()
        except Exception:
            return None, None

        return step, reason
    
    def submit_relax(self, struct_id, structure):
        """
        Submit relaxation job for a structure.
        
        Generates 3 INCAR files (INCAR-1, INCAR-2, INCAR-3) with progressive parameters:
        - INCAR-1: ISIF=2, EDIFFG=-0.02, IBRION=2, POTIM=0.3, EDIFF=1e-7
        - INCAR-2: ISIF=3, EDIFFG=-0.01, IBRION=2, POTIM=0.2, EDIFF=1e-8
        - INCAR-3: ISIF=3, EDIFFG=-0.005, IBRION=2, POTIM=0.1, EDIFF=1e-8
        - All steps: PREC=Accurate, LREAL=.FALSE., LASPH=.TRUE., NPAR=4
        """
        sdata = self.db.get_structure(struct_id)
        if not sdata:
            return False
        
        relax_dir = Path(sdata['relax_dir'])
        job_name = struct_id
        
        print(f"  Submitting Relax: {struct_id}")
        
        try:
            for step_num in [1, 2, 3]:
                self.create_vasp_inputs(structure, relax_dir, 'relax', step=step_num)
                shutil.copy2(relax_dir / 'INCAR', relax_dir / f'INCAR-{step_num}')
            
            script = self.create_slurm_script(relax_dir, job_name, 'relax')
            job_id = self.submit_job(script)
            
            self.db.update_state(struct_id, 'RELAX_RUNNING', relax_job_id=job_id)
            print(f"    Relax job ID: {job_id}")
            return True
        except Exception as e:
            print(f"    Error: {e}")
            self.db.update_state(struct_id, 'RELAX_FAILED', error=str(e))
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
                relax_dir = Path(sdata['relax_dir'])
                
                if local_status == 'DONE':
                    # Successfully completed all 3 steps
                    vasprun_path = relax_dir / 'vasprun.xml'
                    if vasprun_path.exists():
                        try:
                            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
                            if not vr.converged_electronic:
                                self.db.update_state(struct_id, 'RELAX_FAILED', 
                                                   error='Electronic SCF not converged in final step')
                                print(f"  {struct_id}: Relax FAILED (electronic not converged)")
                            else:
                                self.db.update_state(struct_id, 'RELAX_DONE')
                                print(f"  {struct_id}: 3-step relaxation completed")
                        except Exception as e:
                            self.db.update_state(struct_id, 'RELAX_FAILED', 
                                               error=f'Convergence check error: {e}')
                    else:
                        self.db.update_state(struct_id, 'RELAX_FAILED', 
                                           error='vasprun.xml not found')
                
                elif local_status == 'TIMEOUT':
                    # Step 3 timed out but CONTCAR exists (marked by SLURM script)
                    self.db.update_state(struct_id, 'RELAX_TMOUT',
                                       error='Step 3 timed out but produced CONTCAR')
                    print(f"  {struct_id}: Relax TMOUT (step 3 timeout, CONTCAR available)")
                
                elif local_status == 'FAILED':
                    # VASP_FAILED marker exists - load explicit context from script.
                    step, reason = self.read_failure_context(relax_dir)
                    if step and reason:
                        error_msg = f'VASP failed at {step}: {reason}'
                        print(f"  {struct_id}: Relax FAILED ({step}: {reason})")
                    elif step:
                        error_msg = f'VASP failed at {step}'
                        print(f"  {struct_id}: Relax FAILED ({step})")
                    else:
                        error_msg = 'VASP failed (no failure context found)'
                        print(f"  {struct_id}: Relax FAILED (no failure context)")
                    self.db.update_state(struct_id, 'RELAX_FAILED', error=error_msg)
                
                else:
                    # Job not in queue and no completion marker - crashed
                    self.db.update_state(struct_id, 'RELAX_FAILED', 
                                       error='Job terminated without completion marker (crash)')
                    print(f"  {struct_id}: Relax FAILED (crash)")
    
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
            
            for struct_id in list(self.db.data['structures'].keys()):
                self.update_structure_status(struct_id)
            
            running_count = self.db.get_running_count()
            print(f"Currently running: {running_count}/{self.max_concurrent}")
            
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
            
            stats = self.db.get_stats()
            print("\nStatistics:")
            for state, count in sorted(stats['states'].items()):
                print(f"  {state}: {count}")
            sys.stdout.flush()
            
            pending_count = len(self.db.get_by_state('PENDING'))
            if running_count == 0 and pending_count == 0:
                completed = len(self.db.get_by_state('RELAX_DONE'))
                tmout = len(self.db.get_by_state('RELAX_TMOUT'))
                total = stats['total']
                failed_count = len(self.db.get_by_state('RELAX_FAILED'))
                if completed + tmout + failed_count >= total:
                    print("\n" + "="*70)
                    print("All workflows completed!")
                    print(f"Successfully refined: {completed}/{total}")
                    print(f"Timed out (usable): {tmout}/{total}")
                    print(f"Failed: {failed_count}/{total}")
                    print("="*70)
                    sys.stdout.flush()
                    break
            
            print(f"\nSleeping for {self.check_interval}s...")
            sys.stdout.flush()
            time.sleep(self.check_interval)
    
    def initialize_structures(self, input_path, vasp_jobs_dir, output_dir, max_structures=0):
        """
        Initialize database with electride candidates.
        
        Reads structure IDs from electride_candidates.csv or .db
        Loads already-relaxed and symmetrized POSCARs from original VASP_JOBS directory
        Maintains composition-based folder structure like workflow_manager.py
        
        Args:
            input_path: Path to electride_candidates.csv or electride_candidates.db
            vasp_jobs_dir: Path to original VASP_JOBS directory (required)
            output_dir: Output directory for refined calculations (REFINE_VASP_JOBS)
            max_structures: Max total electrides to process (0 or None = all)
        """
        output_dir = Path(output_dir)
        
        print("="*70)
        print("Initializing Refined Electride Workflow")
        print("="*70)
        print(f"Input file: {input_path}")
        print(f"VASP_JOBS dir: {vasp_jobs_dir}")
        print(f"Output directory: {output_dir}")
        print(f"Max concurrent: {self.max_concurrent}")
        if max_structures and max_structures > 0:
            print(f"Max structures: {max_structures}")
        else:
            print(f"Max structures: all")
        print("="*70 + "\n")
        
        # Load electride candidates (structure IDs from input, structures from VASP_JOBS)
        structures_dict = load_electride_candidates(input_path, vasp_jobs_dir)
        
        if not structures_dict:
            print("ERROR: No electride candidates loaded!")
            return {}
        
        # Sort by structure ID and limit if max_structures is specified
        sorted_struct_ids = sorted(structures_dict.keys())
        if max_structures and max_structures > 0:
            sorted_struct_ids = sorted_struct_ids[:max_structures]
            structures_dict = {sid: structures_dict[sid] for sid in sorted_struct_ids}
            print(f"\nLimiting to first {len(sorted_struct_ids)} electrides (max_structures={max_structures})\n")
        
        added_count = 0
        
        # Group by composition for better organization
        by_composition = {}
        for struct_id in structures_dict.keys():
            comp_name = struct_id.rsplit('_s', 1)[0]
            if comp_name not in by_composition:
                by_composition[comp_name] = []
            by_composition[comp_name].append(struct_id)
        
        print(f"Processing {len(structures_dict)} electrides from {len(by_composition)} compositions\n")
        
        for comp_name in sorted(by_composition.keys()):
            struct_ids = by_composition[comp_name]
            print(f"Composition {comp_name}: {len(struct_ids)} structures")
            
            for struct_id in struct_ids:
                struct_idx = int(struct_id.rsplit('_s', 1)[1])
                structure = structures_dict[struct_id]
                
                if struct_id not in self.db.data['structures']:
                    elements = sorted([str(el) for el in structure.composition.elements])
                    chemsys = '-'.join(elements)
                    
                    # Maintain composition-based folder structure: output_dir/comp_name/struct_id/
                    self.db.add_structure(
                        struct_id, comp_name, struct_idx,
                        output_dir / comp_name,
                        chemsys=chemsys
                    )
                    added_count += 1
        
        self.db.data['config'] = {
            'max_concurrent': self.max_concurrent,
            'input_path': str(input_path),
            'vasp_jobs_dir': str(vasp_jobs_dir),
            'output_dir': str(output_dir),
            'max_structures': max_structures
        }
        self.db.save()
        
        print(f"\nStructures loaded from VASP_JOBS: {len(structures_dict)}")
        print(f"New structures added to workflow DB: {added_count}")
        print(f"\nFolder structure: {output_dir}/{{composition}}/{{struct_id}}/Relax")
        print(f"Output: High-quality CONTCAR from 3-step relaxation")
        
        return structures_dict


def main():
    parser = argparse.ArgumentParser(
        description="Refined Electride Workflow Manager - High-precision VASP calculations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Using database file
  python3 refine_electrideflow.py \\
      --input electride_candidates.db \\
      --vasp-jobs-dir VASP_JOBS/
  
  # Using CSV file
  python3 refine_electrideflow.py \\
      --input electride_candidates.csv \\
      --vasp-jobs-dir VASP_JOBS/
  
  # Limit to first 100 structures
  python3 refine_electrideflow.py \\
      --input electride_candidates.db \\
      --vasp-jobs-dir VASP_JOBS/ \\
      --max-structures 100

Note on structure source:
  - Structure IDs are read from the input file (CSV or database)
  - Actual structures (POSCARs) are loaded from original VASP_JOBS directory
  - This ensures we use the same symmetrized structures that the original
    workflow successfully created, avoiding any symmetrization discrepancies

Filtering:
  - Only processes high-symmetry structures (space group > 15)
  - Excludes triclinic (space groups 1-2) and monoclinic (3-15)
  - Focuses computational resources on orthorhombic and higher symmetry systems

Prerequisites:
  1. Run workflow_manager.py to complete initial VASP calculations
  2. Run analyze.py to identify electride candidates
  3. Run filter_comb_db.py to create electride_candidates.db/csv
"""
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help="Path to electride_candidates.csv or electride_candidates.db"
    )
    parser.add_argument(
        '--vasp-jobs-dir',
        type=str,
        required=True,
        help="Path to original VASP_JOBS directory (required)"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./REFINE_VASP_JOBS',
        help="Output directory for refined VASP jobs (default: ./REFINE_VASP_JOBS)"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="JSON database file path (default: workflow.json)"
    )
    parser.add_argument(
        '--max-concurrent',
        type=int,
        default=10,
        help="Max concurrent structures running (default: 10)"
    )
    parser.add_argument(
        '--check-interval',
        type=int,
        default=60,
        help="Status check interval in seconds (default: 60)"
    )
    parser.add_argument(
        '--max-structures',
        type=int,
        default=0,
        help="Max total electrides to process (default: 0 = all)"
    )
    parser.add_argument(
        '--init-only',
        action='store_true',
        help="Only initialize database, don't start monitoring"
    )
    
    args = parser.parse_args()
    
    input_path = Path(args.input).expanduser()
    vasp_jobs_dir = Path(args.vasp_jobs_dir).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = output_dir / args.db
    
    # Check if input file exists
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}")
        print("\nPlease run filter_comb_db.py first to create electride_candidates.db or .csv")
        sys.exit(1)
    
    # Check if VASP_JOBS directory exists
    if not vasp_jobs_dir.exists():
        print(f"ERROR: VASP_JOBS directory not found: {vasp_jobs_dir}")
        print("\nPlease provide valid path to original VASP_JOBS directory")
        sys.exit(1)
    
    # Create workflow manager
    manager = VASPWorkflowManager(
        db_path=db_path,
        max_concurrent=args.max_concurrent,
        check_interval=args.check_interval
    )
    
    # Initialize structures
    structures_dict = manager.initialize_structures(
        input_path=input_path,
        vasp_jobs_dir=vasp_jobs_dir,
        output_dir=output_dir,
        max_structures=args.max_structures
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
    
    # Start monitoring and submission loop
    try:
        manager.monitor_and_submit(structures_dict)
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user.")
        print(f"Database saved to: {db_path}")
        print("Resume with same command to continue.")


if __name__ == '__main__':
    main()
