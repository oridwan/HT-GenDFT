#!/bin/bash
#SBATCH --job-name=electronic_workflow
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=electronic_workflow_%j.out
#SBATCH --error=electronic_workflow_%j.err

# Electronic Structure Workflow Manager
# Runs SCF/PARCHG/ELF/BAND/PDOS calculations on refined electride structures

set -e

# Resolve robust paths under SLURM (script may execute from spool copy)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
if [ -f "$SUBMIT_DIR/refined_flow/electronic_flow_manager.py" ]; then
    FLOW_DIR="$SUBMIT_DIR/refined_flow"
elif [ -f "$SUBMIT_DIR/electronic_flow_manager.py" ]; then
    FLOW_DIR="$SUBMIT_DIR"
elif [ -f "$SCRIPT_DIR/electronic_flow_manager.py" ]; then
    FLOW_DIR="$SCRIPT_DIR"
else
    FLOW_DIR="$SCRIPT_DIR"
fi

# Default parameters (can be overridden by environment variables)
REFINE_JOBS=${REFINE_JOBS:-"/scratch/oridwan/SuperConductorFlow/REFINE_VASP-out-Boron_2"}
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/oridwan/SuperConductorFlow/ELECTRONIC_JOBS-Boron"}
MAX_CONCURRENT=${MAX_CONCURRENT:-100}
CHECK_INTERVAL=${CHECK_INTERVAL:-60}
CONDA_ENV=${CONDA_ENV:-"mattersim"}
DB_NAME=${DB_NAME:-"workflow.json"}
STRUCTURE_IDS=${STRUCTURE_IDS:-""}

echo "========================================================================"
echo "Electronic Structure Workflow Manager (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo ""
echo "Workflow stages per structure:"
echo "  1. SCF (Self-Consistent Field) - conventional cell"
echo "  2. PARCHG (5 energy windows) - conventional cell"
echo "  3. ELF (Electron Localization) - conventional cell"
echo "  4. BAND (Band Structure) - primitive cell"
echo "  5. PDOS (Projected DOS) - primitive cell"
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

if ! python3 -c "from pymatgen.symmetry.bandstructure import HighSymmKpath" 2>/dev/null; then
    echo "Error: pymatgen.symmetry.bandstructure not available"
    exit 1
fi

# Expand paths
REFINE_JOBS=$(eval echo "$REFINE_JOBS")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check refined jobs directory
if [ ! -d "$REFINE_JOBS" ]; then
    echo "Error: Refined jobs directory not found: $REFINE_JOBS"
    echo "Please provide valid path to REFINE_VASP_JOBS directory"
    exit 1
fi

# Check refined workflow database
REFINE_DB="$REFINE_JOBS/workflow.json"
if [ ! -f "$REFINE_DB" ]; then
    echo "Error: Refined workflow database not found: $REFINE_DB"
    echo "Please complete the refined relaxation workflow first"
    exit 1
fi

# Count RELAX_DONE structures
RELAX_DONE_COUNT=$(grep -c '"state": "RELAX_DONE"' "$REFINE_DB" || true)
echo "Structures in RELAX_DONE state: $RELAX_DONE_COUNT"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check workflow manager script path
FLOW_SCRIPT="$FLOW_DIR/electronic_flow_manager.py"
if [ ! -f "$FLOW_SCRIPT" ]; then
    echo "Error: electronic workflow manager script not found: $FLOW_SCRIPT"
    echo "SLURM submit dir: $SUBMIT_DIR"
    echo "Script dir: $SCRIPT_DIR"
    exit 1
fi

# Build command
CMD="python3 $FLOW_SCRIPT"
CMD="$CMD --refine-jobs $REFINE_JOBS"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --max-concurrent $MAX_CONCURRENT"
CMD="$CMD --check-interval $CHECK_INTERVAL"
CMD="$CMD --db $DB_NAME"

# Add structure IDs if specified
if [ -n "$STRUCTURE_IDS" ]; then
    CMD="$CMD --structure-ids$STRUCTURE_IDS"
fi

# Print configuration
echo "Configuration:"
echo "  Refined jobs: $REFINE_JOBS"
echo "  Output dir: $OUTPUT_DIR"
echo "  Max concurrent: $MAX_CONCURRENT"
if [ -n "$STRUCTURE_IDS" ]; then
    echo "  Structure IDs:$STRUCTURE_IDS"
else
    echo "  Structure IDs: all RELAX_DONE structures"
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
echo "Starting electronic workflow manager..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run workflow manager
$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Electronic workflow manager finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

exit $EXIT_CODE
