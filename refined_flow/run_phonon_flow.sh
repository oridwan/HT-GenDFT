#!/bin/bash
# Phonon Workflow - Submits workflow manager as SLURM job
# Usage: bash run_phonon_flow.sh [options]
# Cs2K1B2H8_s001
''' 
STRUCT_IDS=$(python3 - <<'PY'
import json
p="/scratch/oridwan/SuperConductorFlow/REFINE_VASP-out-Boron_2/workflow.json"
d=json.load(open(p))["structures"]
print(" ".join(sid for sid,s in d.items() if s.get("state")=="RELAX_DONE"))
PY
)

bash /scratch/oridwan/SuperConductorFlow/refined_flow/run_phonon_flow.sh \
  --refine-jobs /scratch/oridwan/SuperConductorFlow/REFINE_VASP-out-Boron_3 \
  --output-dir /scratch/oridwan/SuperConductorFlow/PHONON_Boron_check3 \
  --structure-ids $STRUCT_IDS


'''
set -e

# Resolve script directory so this launcher works from any current directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
REFINE_JOBS="./REFINE_VASP-out-Boron_2"
OUTPUT_DIR="./PHONON_Boron_check2"
MAX_CONCURRENT=100
CHECK_INTERVAL=60
SUPERCELL_DIM=""
DB_NAME="workflow.json"
CONDA_ENV="mattersim"
STRUCTURE_IDS=""
VASP_JOB_TIME="12:00:00"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "========================================================================"
echo "Phonon Workflow - Phonopy Finite Displacement Method"
echo "========================================================================"
echo ""
echo "Workflow:"
echo "  1. Generate supercell with phonopy"
echo "  2. Create displaced structures (finite displacements)"
echo "  3. Run VASP static calculations on each displacement"
echo "  4. Monitor all displacement calculations"
echo ""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --refine-jobs)
            REFINE_JOBS="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --structure-ids)
            shift
            while [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-- ]]; do
                STRUCTURE_IDS="$STRUCTURE_IDS $1"
                shift
            done
            ;;
        --max-concurrent)
            MAX_CONCURRENT="$2"
            shift 2
            ;;
        --check-interval)
            CHECK_INTERVAL="$2"
            shift 2
            ;;
        --supercell-dim)
            SUPERCELL_DIM="$2 $3 $4"
            shift 4
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --vasp-job-time)
            VASP_JOB_TIME="$2"
            shift 2
            ;;
        --help)
            echo "Usage: bash run_phonon_flow.sh [options]"
            echo ""
            echo "This script submits the phonon workflow manager as a SLURM job."
            echo "It runs phonon calculations using phonopy finite displacement method."
            echo ""
            echo "Options:"
            echo "  --refine-jobs DIR          Path to refined jobs directory (default: ./REFINE_VASP-out-Boron)"
            echo "  --output-dir DIR           Output directory for phonon jobs (default: ./PHONON_JOBS_Boron)"
            echo "  --structure-ids ID1 ID2... Structure IDs to process (REQUIRED, must be RELAX_DONE)"
            echo "  --max-concurrent N         Max concurrent displacement jobs (default: 100)"
            echo "  --check-interval SECONDS   Status check interval (default: 60)"
            echo "  --supercell-dim N1 N2 N3   Manual supercell dimensions (default: auto from k-mesh)"
            echo "  --conda-env NAME           Conda environment name (default: mattersim)"
            echo "  --vasp-job-time HH:MM:SS   Per-displacement VASP walltime (default: 12:00:00)"
            echo ""
            echo "Examples:"
            echo "  # Process structures with auto supercell (default behavior)"
            echo "  bash run_phonon_flow.sh --structure-ids Ba2N_s001 Ca5P3_s021"
            echo ""
            echo "  # Manual supercell dimensions (overrides auto)"
            echo "  bash run_phonon_flow.sh --structure-ids Ba2N_s001 --supercell-dim 3 3 3"
            echo ""
            echo "  # Custom concurrency"
            echo "  bash run_phonon_flow.sh --structure-ids Ba2N_s001 --max-concurrent 10"
            echo ""
            echo "Prerequisites:"
            echo "  1. Complete refined relaxation workflow (run_refineflow.sh)"
            echo "  2. Structures must be in RELAX_DONE state"
            echo "  3. phonopy must be installed in conda environment"
            echo ""
            echo "Monitoring:"
            echo "  Check status:    python workflow_status.py ./PHONON_JOBS/workflow.json"
            echo "  View log:        tail -f phonon_workflow_<JOBID>.out"
            echo "  Check queue:     squeue -u \$USER"
            echo "  Cancel:          scancel <JOBID>"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if structure IDs provided
if [ -z "$STRUCTURE_IDS" ]; then
    echo -e "${RED}Error: No structure IDs specified${NC}"
    echo ""
    echo "Please provide structure IDs with --structure-ids option"
    echo "Example: bash run_phonon_flow.sh --structure-ids Ba2N_s001 Ca5P3_s021"
    echo ""
    echo "Use --help for more information"
    exit 1
fi

# Expand paths
REFINE_JOBS=$(eval echo "$REFINE_JOBS")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check if refined jobs directory exists
if [ ! -d "$REFINE_JOBS" ]; then
    echo -e "${RED}Error: Refined jobs directory not found: $REFINE_JOBS${NC}"
    echo ""
    echo "Please provide valid path to REFINE_VASP_JOBS directory"
    echo "Use: --refine-jobs /path/to/REFINE_VASP_JOBS"
    exit 1
fi

# Check if refined workflow database exists
REFINE_DB="$REFINE_JOBS/workflow.json"
if [ ! -f "$REFINE_DB" ]; then
    echo -e "${RED}Error: Refined workflow database not found: $REFINE_DB${NC}"
    echo ""
    echo "Please complete the refined relaxation workflow first (run_refineflow.sh)"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if database exists
DB_PATH="$OUTPUT_DIR/$DB_NAME"
if [ -f "$DB_PATH" ]; then
    echo -e "${YELLOW}Database already exists: $DB_PATH${NC}"
    echo "The workflow will resume from previous state."
    echo ""
fi

# Print configuration
echo -e "${GREEN}Configuration:${NC}"
echo "  Refined jobs: $REFINE_JOBS"
echo "  Output dir: $OUTPUT_DIR"
echo "  Max concurrent: $MAX_CONCURRENT"
echo "  Structure IDs:$STRUCTURE_IDS"
if [ -n "$SUPERCELL_DIM" ]; then
    echo "  Supercell dimensions: $SUPERCELL_DIM (manual)"
else
    echo "  Supercell: AUTO (from k-point mesh)"
fi
echo "  Check interval: ${CHECK_INTERVAL}s"
echo "  Database: $DB_PATH"
echo "  Conda env: $CONDA_ENV"
echo "  Displacement walltime: $VASP_JOB_TIME"
echo ""

# Export variables for SLURM script
export REFINE_JOBS
export OUTPUT_DIR
export MAX_CONCURRENT
export CHECK_INTERVAL
export SUPERCELL_DIM
export DB_NAME
export CONDA_ENV
export STRUCTURE_IDS
export VASP_JOB_TIME

# Submit the workflow manager job
echo "Submitting phonon workflow manager as SLURM job..."
SUBMIT_SCRIPT="$SCRIPT_DIR/submit_phonon_flow.sh"
if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo -e "${RED}Error: Submit script not found: $SUBMIT_SCRIPT${NC}"
    exit 1
fi
JOBID=$(sbatch "$SUBMIT_SCRIPT" | awk '{print $NF}')

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Phonon workflow manager submitted successfully!${NC}"
    echo ""
    echo "Job ID: $JOBID"
    echo ""
    echo "Monitoring commands:"
    echo "  View log:        tail -f phonon_workflow_${JOBID}.out"
    echo "  Check status:    squeue -j $JOBID"
    echo "  Check workflow:  python workflow_status.py $OUTPUT_DIR/workflow.json"
    echo "  Cancel job:      scancel $JOBID"
    echo ""
    echo "The workflow will:"
    if [ -n "$SUPERCELL_DIM" ]; then
        echo "  - Generate supercells with phonopy (dim: $SUPERCELL_DIM)"
    else
        echo "  - Generate supercells with phonopy (AUTO dimensions from k-mesh)"
    fi
    echo "  - Create displaced structures"
    echo "  - Submit VASP calculations for each displacement"
    echo "  - Monitor completion (max $MAX_CONCURRENT concurrent jobs)"
    echo "  - Save results to $OUTPUT_DIR"
    echo ""
else
    echo -e "${RED}Error: Failed to submit job${NC}"
    exit 1
fi
