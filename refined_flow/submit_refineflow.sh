#!/bin/bash
#SBATCH --job-name=refine_electride
#SBATCH --partition=Apus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=10-00:00:00
#SBATCH --output=refine_workflow_%j.out
#SBATCH --error=refine_workflow_%j.err

# Refined Electride Workflow Manager - High-precision VASP calculations
# Runs on electride candidates identified by analyze.py
set -e

# Resolve robust paths under SLURM (script may execute from spool copy)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
if [ -f "$SUBMIT_DIR/refined_flow/refine_electrideflow.py" ]; then
    FLOW_DIR="$SUBMIT_DIR/refined_flow"
elif [ -f "$SUBMIT_DIR/refine_electrideflow.py" ]; then
    FLOW_DIR="$SUBMIT_DIR"
elif [ -f "$SCRIPT_DIR/refine_electrideflow.py" ]; then
    FLOW_DIR="$SCRIPT_DIR"
else
    FLOW_DIR="$SCRIPT_DIR"
fi

# Default parameters (can be overridden by environment variables)
INPUT_FILE=${INPUT_FILE:-"/scratch/oridwan/SuperConductorFlow/VASP-out-Boron/electride_analysis_candidates-Strict.csv"}
VASP_JOBS_DIR=${VASP_JOBS_DIR:-"/scratch/oridwan/SuperConductorFlow/VASP-out-Boron"}
OUTPUT_DIR=${OUTPUT_DIR:-"/scratch/oridwan/SuperConductorFlow/REFINE_VASP-out-Boron_2"}
MAX_CONCURRENT=${MAX_CONCURRENT:-100}
MAX_STRUCTURES=${MAX_STRUCTURES:-0}
CHECK_INTERVAL=${CHECK_INTERVAL:-60}
CONDA_ENV=${CONDA_ENV:-"mattersim"}
DB_NAME=${DB_NAME:-"workflow.json"}

echo "========================================================================"
echo "Refined Electride Workflow Manager (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo ""
echo "REFINED Settings:"
echo "  - K-point density: 250 (vs 64 in original)"
echo "  - RELAX NELM: 60 (more electronic steps)"
echo "  - RELAX NSW: 100 (more ionic steps, vs 30 in original)"
echo "  - 3-step relaxation with adaptive Step-3 reruns (max 3 reruns)"
echo "  - INCAR policy: IBRION=2 for all steps, NPAR=4"
echo "  - Step-3 stop criteria: reached required accuracy + follow-up ionic steps <= 3"
echo "  - SLURM time per structure job: 12 hours (vs 20 minutes in original)"
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

# Check pymatgen
if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "Error: pymatgen not found"
    exit 1
fi

# Expand paths
INPUT_FILE=$(eval echo "$INPUT_FILE")
VASP_JOBS_DIR=$(eval echo "$VASP_JOBS_DIR")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check input file
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    echo "Please run filter_comb_db.py first to create electride_candidates.db or .csv"
    exit 1
fi

# Count electride candidates
if [[ "$INPUT_FILE" == *.csv ]]; then
    ELECTRIDE_COUNT=$(($(wc -l < "$INPUT_FILE") - 1))
    echo "Electride candidates found: $ELECTRIDE_COUNT"
elif [[ "$INPUT_FILE" == *.db ]] && command -v ase &> /dev/null; then
    ELECTRIDE_COUNT=$(ase db "$INPUT_FILE" -c | tail -1 | awk '{print $1}')
    echo "Electride candidates found: $ELECTRIDE_COUNT"
else
    echo "Warning: Cannot count structures in $INPUT_FILE"
fi
echo ""

# Check VASP_JOBS directory (required)
if [ ! -d "$VASP_JOBS_DIR" ]; then
    echo "Error: VASP_JOBS directory not found: $VASP_JOBS_DIR"
    echo "Please provide valid path to original VASP_JOBS directory"
    exit 1
fi

# Check refined workflow manager script path
FLOW_SCRIPT="$FLOW_DIR/refine_electrideflow.py"
if [ ! -f "$FLOW_SCRIPT" ]; then
    echo "Error: refine workflow manager script not found: $FLOW_SCRIPT"
    echo "SLURM submit dir: $SUBMIT_DIR"
    echo "Script dir: $SCRIPT_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build command
CMD="python3 $FLOW_SCRIPT"
CMD="$CMD --input $INPUT_FILE"
CMD="$CMD --vasp-jobs-dir $VASP_JOBS_DIR"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --max-concurrent $MAX_CONCURRENT"
CMD="$CMD --max-structures $MAX_STRUCTURES"
CMD="$CMD --check-interval $CHECK_INTERVAL"
CMD="$CMD --db $DB_NAME"

# Print configuration
echo "Configuration:"
echo "  Input file: $INPUT_FILE"
echo "  VASP_JOBS dir: $VASP_JOBS_DIR"
echo "  Output dir: $OUTPUT_DIR"
echo "  Max concurrent: $MAX_CONCURRENT"
if [ "$MAX_STRUCTURES" -gt 0 ] 2>/dev/null; then
    echo "  Max structures: $MAX_STRUCTURES"
else
    echo "  Max structures: all"
fi
echo "  Check interval: ${CHECK_INTERVAL}s"
echo "  Database: $OUTPUT_DIR/$DB_NAME"
echo ""
echo "Note: Structures (CONTCARs) will be loaded from original VASP_JOBS directory"
echo ""

# Check if resuming
if [ -f "$OUTPUT_DIR/$DB_NAME" ]; then
    echo "Database exists - resuming from previous state"
    echo ""
fi

echo "========================================================================"
echo "Starting refined workflow manager..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run workflow manager
$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Refined workflow manager finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

exit $EXIT_CODE
