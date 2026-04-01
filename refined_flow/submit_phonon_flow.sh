#!/bin/bash
#SBATCH --job-name=phonon_workflow
#SBATCH --partition=Apus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=5-00:00:00
#SBATCH --output=phonon_workflow_%j.out
#SBATCH --error=phonon_workflow_%j.err

# Phonon Workflow Manager
# Runs phonopy finite displacement phonon calculations on refined electride structures

set -e

# Resolve workflow manager path (avoid SLURM spool script location)
PROJECT_ROOT="${PROJECT_ROOT:-/scratch/oridwan/SuperConductorFlow}"
FLOW_DIR="${FLOW_DIR:-$PROJECT_ROOT/refined_flow}"

# Default parameters
REFINE_JOBS=${REFINE_JOBS:-"/scratch/oridwan/SuperConductorFlow/REFINE_VASP-out-Boron_3"}
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/oridwan/SuperConductorFlow/PHONON_Boron_check-3"}
MAX_CONCURRENT=${MAX_CONCURRENT:-100}
CHECK_INTERVAL=${CHECK_INTERVAL:-60}
SUPERCELL_DIM=${SUPERCELL_DIM:-""}
CONDA_ENV=${CONDA_ENV:-"mattersim"}
DB_NAME=${DB_NAME:-"workflow.json"}
STRUCTURE_IDS=${STRUCTURE_IDS:-""}
VASP_JOB_TIME=${VASP_JOB_TIME:-"12:00:00"}
export VASP_JOB_TIME

echo "========================================================================"
echo "Phonon Workflow Manager (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo ""
echo "Workflow:"
if [ -n "$SUPERCELL_DIM" ]; then
    echo "  1. Generate supercell with phonopy (dim: $SUPERCELL_DIM)"
else
    echo "  1. Generate supercell with phonopy (AUTO from k-mesh)"
fi
echo "  2. Create displaced structures"
echo "  3. Submit VASP static calculations for each displacement"
echo "  4. Monitor and manage completion"
echo ""

# Load environment
source ~/.bashrc
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment '$CONDA_ENV'"
    exit 1
fi

echo "Conda environment: $CONDA_ENV"
echo "Python: $(which python3)"
echo ""

# Check dependencies
if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "Error: pymatgen not found"
    exit 1
fi

if ! python3 -c "import phonopy" 2>/dev/null; then
    echo "Error: phonopy not found. Please install: conda install -c conda-forge phonopy"
    exit 1
fi

echo "phonopy version: $(python3 -c 'import phonopy; print(phonopy.__version__)')"
echo ""

# Expand paths
REFINE_JOBS=$(eval echo "$REFINE_JOBS")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check refined jobs directory
if [ ! -d "$REFINE_JOBS" ]; then
    echo "Error: Refined jobs directory not found: $REFINE_JOBS"
    exit 1
fi

# Check refined workflow database
REFINE_DB="$REFINE_JOBS/workflow.json"
if [ ! -f "$REFINE_DB" ]; then
    echo "Error: Refined workflow database not found: $REFINE_DB"
    echo "Please complete the refined relaxation workflow first"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if structure IDs provided
if [ -z "$STRUCTURE_IDS" ]; then
    echo "Error: No structure IDs specified"
    echo "Please provide --structure-ids option"
    exit 1
fi

# Build command
if [ ! -f "$FLOW_DIR/phonon_flow_manager.py" ]; then
    echo "Error: phonon_flow_manager.py not found at: $FLOW_DIR/phonon_flow_manager.py"
    echo "Set FLOW_DIR or PROJECT_ROOT to the correct absolute path"
    exit 1
fi

CMD="python3 $FLOW_DIR/phonon_flow_manager.py"
CMD="$CMD --refine-jobs $REFINE_JOBS"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --max-concurrent $MAX_CONCURRENT"
CMD="$CMD --check-interval $CHECK_INTERVAL"
CMD="$CMD --db $DB_NAME"

# Add supercell dimensions if manually specified
if [ -n "$SUPERCELL_DIM" ]; then
    CMD="$CMD --supercell-dim $SUPERCELL_DIM"
fi

CMD="$CMD --structure-ids $STRUCTURE_IDS"

# Print configuration
echo "Configuration:"
echo "  Refined jobs: $REFINE_JOBS"
echo "  Output dir: $OUTPUT_DIR"
echo "  Max concurrent: $MAX_CONCURRENT"
echo "  Structure IDs:$STRUCTURE_IDS"
echo "  Displacement walltime: $VASP_JOB_TIME"
if [ -n "$SUPERCELL_DIM" ]; then
    echo "  Supercell dim: $SUPERCELL_DIM (manual)"
else
    echo "  Supercell: AUTO (from k-point mesh)"
fi
echo "  Check interval: ${CHECK_INTERVAL}s"
echo "  Database: $OUTPUT_DIR/$DB_NAME"
echo ""

# Check if resuming
if [ -f "$OUTPUT_DIR/$DB_NAME" ]; then
    echo "Database exists - resuming from previous state"
    echo ""
fi

echo "========================================================================"
echo "Starting phonon workflow manager..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run workflow manager
$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Phonon workflow manager finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

exit $EXIT_CODE
