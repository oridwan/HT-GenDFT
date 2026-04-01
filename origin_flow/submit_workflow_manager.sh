#!/bin/bash
#SBATCH --job-name=vaspflow_manager
#SBATCH --partition=Apus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=10-00:00:00
#SBATCH --output=workflow_manager_%j.out
#SBATCH --error=workflow_manager_%j.err
#SBATCH --chdir=/scratch/oridwan/SuperConductorFlow
#SBATCH --exclude=str-c88,str-c89,str-c90,str-c91,str-c92,str-c93,str-c94,str-c95,str-c96,str-c97

# VASPflow Workflow Manager - VASP calculations only

# Pre-screening should be done separately with prescreen.py

set -e

# Fixed project root for simple, predictable execution under SLURM
REPO_ROOT="/scratch/oridwan/SuperConductorFlow"
cd "$REPO_ROOT"

# Default parameters (can be overridden by environment variables)
RESULTS_DIR=${RESULTS_DIR:-"/scratch/oridwan/SuperConductorFlow/prescreen_new/filter/cif_1e_1_to_5e_2/"}
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/oridwan/SuperConductorFlow/VASP-out-Boron"}
MAX_CONCURRENT=${MAX_CONCURRENT:-200}
MAX_COMPOSITIONS=${MAX_COMPOSITIONS:-""}
MAX_STRUCTURES=${MAX_STRUCTURES:-100}
CHECK_INTERVAL=${CHECK_INTERVAL:-60}
CONDA_ENV=${CONDA_ENV:-"base"}
DB_NAME=${DB_NAME:-"workflow.json"}
VASP_PP_PATH=${VASP_PP_PATH:-"/projects/mmi/Ridwan/potcarFiles/VASP6.4/potpaw_PBE"}
VASP_JOB_CONSTRAINT=${VASP_JOB_CONSTRAINT:-""}
VASP_JOB_EXCLUDE_NODES=${VASP_JOB_EXCLUDE_NODES:-"str-c88,str-c89,str-c90,str-c91,str-c92,str-c93,str-c94,str-c95,str-c96,str-c97"}

echo "========================================================================"
echo "VASPflow Workflow Manager (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "Repo root: $REPO_ROOT"
echo ""

# Load environment
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
    source ~/.bashrc
fi

conda activate "$CONDA_ENV"

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment '$CONDA_ENV'"
    exit 1
fi

echo "Conda environment: $CONDA_ENV"
echo "Python: $(which python3)"
echo ""

# Export pseudopotential paths used by pymatgen/VASP
export VASP_PP_PATH="$VASP_PP_PATH"
export PMG_VASP_PSP_DIR="$VASP_PP_PATH"
export VASP_JOB_CONSTRAINT="$VASP_JOB_CONSTRAINT"
export VASP_JOB_EXCLUDE_NODES="$VASP_JOB_EXCLUDE_NODES"
echo "VASP_PP_PATH: $VASP_PP_PATH"
echo "PMG_VASP_PSP_DIR: $PMG_VASP_PSP_DIR"
echo "VASP_JOB_CONSTRAINT: ${VASP_JOB_CONSTRAINT:-<none>}"
echo "VASP_JOB_EXCLUDE_NODES: ${VASP_JOB_EXCLUDE_NODES:-<none>}"
echo ""

# Check pymatgen
if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "Error: pymatgen not found"
    exit 1
fi

# Expand paths
RESULTS_DIR=$(eval echo "$RESULTS_DIR")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check results directory
if [ ! -d "$RESULTS_DIR" ]; then
    echo "Error: Results directory not found: $RESULTS_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Resolve workflow manager path from repository root
MANAGER_PY="$REPO_ROOT/origin_flow/workflow_manager.py"
if [ ! -f "$MANAGER_PY" ]; then
    echo "Error: workflow_manager.py not found: $MANAGER_PY"
    exit 1
fi

# Build command
CMD=(
    python3 "$MANAGER_PY"
    --results-dir "$RESULTS_DIR"
    --output-dir "$OUTPUT_DIR"
    --max-concurrent "$MAX_CONCURRENT"
    --max-structures "$MAX_STRUCTURES"
    --check-interval "$CHECK_INTERVAL"
    --db "$DB_NAME"
)

if [ -n "$MAX_COMPOSITIONS" ]; then
    CMD+=(--max-compositions "$MAX_COMPOSITIONS")
fi

# Print configuration
echo "Configuration:"
echo "  Results dir: $RESULTS_DIR"
echo "  Output dir: $OUTPUT_DIR"
echo "  Max concurrent: $MAX_CONCURRENT"
echo "  Max structures: $MAX_STRUCTURES"
echo "  Max compositions: ${MAX_COMPOSITIONS:-all}"
echo "  Check interval: ${CHECK_INTERVAL}s"
echo "  Database: $OUTPUT_DIR/$DB_NAME"
echo ""

# Check if resuming
if [ -f "$OUTPUT_DIR/$DB_NAME" ]; then
    echo "Database exists - resuming from previous state"
    echo ""
fi

echo "========================================================================"
echo "Starting workflow manager..."
echo "========================================================================"
echo ""
echo "Command: ${CMD[*]}"
echo ""

# Run workflow manager
"${CMD[@]}"

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Workflow manager finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

exit $EXIT_CODE
