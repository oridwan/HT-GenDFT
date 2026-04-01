#!/bin/bash
# Refined Electride Workflow - Submits workflow manager as SLURM job
# Usage: bash run_refineflow.sh [options]

set -e

# Default values
INPUT_FILE="./electride_candidates.db"
VASP_JOBS_DIR="./VASP_JOBS"
OUTPUT_DIR="./REFINE_VASP_JOBS"
MAX_CONCURRENT=10
MAX_STRUCTURES=0
CHECK_INTERVAL=60
DB_NAME="workflow.json"
CONDA_ENV="vaspflow"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "========================================================================"
echo "Refined Electride Workflow - High-Precision VASP Calculations"
echo "========================================================================"
echo ""
echo "Settings:"
echo "  - K-point density: 250 (vs 64 in original)"
echo "  - RELAX NELM: 60 (more electronic steps)"
echo "  - RELAX NSW: 100 (more ionic steps, vs 30 in original)"
echo "  - 3-step consecutive relaxation workflow"
echo "  - SLURM time: 4 hours (vs 20 minutes in original)"
echo ""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_FILE="$2"
            shift 2
            ;;
        --vasp-jobs-dir)
            VASP_JOBS_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --max-concurrent)
            MAX_CONCURRENT="$2"
            shift 2
            ;;
        --max-structures)
            MAX_STRUCTURES="$2"
            shift 2
            ;;
        --check-interval)
            CHECK_INTERVAL="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --help)
            echo "Usage: bash run_refineflow.sh [options]"
            echo ""
            echo "This script submits the refined electride workflow manager as a SLURM job."
            echo "It runs high-precision VASP calculations on filtered electride candidates."
            echo ""
            echo "Options:"
            echo "  --input FILE               Path to electride_candidates.csv or .db (default: ./electride_candidates.db)"
            echo "  --vasp-jobs-dir DIR        Path to original VASP_JOBS directory (default: ./VASP_JOBS)"
            echo "  --output-dir DIR           Output directory for refined jobs (default: ./REFINE_VASP_JOBS)"
            echo "  --max-concurrent N         Max concurrent structures (default: 10)"
            echo "  --max-structures N         Max total electrides to process (default: 0 = all)"
            echo "  --check-interval SECONDS   Status check interval (default: 60)"
            echo "  --conda-env NAME           Conda environment name (default: vaspflow)"
            echo ""
            echo "Example:"
            echo "  bash run_refineflow.sh --max-concurrent 20"
            echo "  bash run_refineflow.sh --input electride_candidates.db"
            echo "  bash run_refineflow.sh --input electride_candidates.csv --vasp-jobs-dir /path/to/VASP_JOBS/"
            echo ""
            echo "Note:"
            echo "  - Structure IDs are read from the input file (CSV or database)"
            echo "  - Actual structures (POSCARs) are loaded from original VASP_JOBS directory"
            echo "  - This ensures we use the same symmetrized structures as the original workflow"
            echo ""
            echo "Prerequisites:"
            echo "  1. Run workflow_manager.py to complete initial VASP calculations"
            echo "  2. Run analyze.py to identify electride candidates"
            echo "  3. Run filter_comb_db.py to create electride_candidates.db/csv"
            echo ""
            echo "Monitoring:"
            echo "  Check status:    python workflow_status.py ./REFINE_VASP_JOBS/workflow.json"
            echo "  View log:        tail -f refine_workflow_<JOBID>.out"
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

# Expand paths
INPUT_FILE=$(eval echo "$INPUT_FILE")
VASP_JOBS_DIR=$(eval echo "$VASP_JOBS_DIR")
OUTPUT_DIR=$(eval echo "$OUTPUT_DIR")

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Error: Input file not found: $INPUT_FILE${NC}"
    echo ""
    echo "Please run filter_comb_db.py first to create electride_candidates.db or .csv"
    exit 1
fi

# Count electride candidates
if [[ "$INPUT_FILE" == *.csv ]]; then
    ELECTRIDE_COUNT=$(($(wc -l < "$INPUT_FILE") - 1))
    echo -e "${GREEN}Found $ELECTRIDE_COUNT electride candidates in CSV${NC}"
elif [[ "$INPUT_FILE" == *.db ]] && command -v ase &> /dev/null; then
    ELECTRIDE_COUNT=$(ase db "$INPUT_FILE" -c | tail -1 | awk '{print $1}')
    echo -e "${GREEN}Found $ELECTRIDE_COUNT electride candidates in database${NC}"
else
    echo -e "${YELLOW}Warning: Cannot count structures in $INPUT_FILE${NC}"
fi
echo ""

# Check if VASP_JOBS directory exists (required)
if [ ! -d "$VASP_JOBS_DIR" ]; then
    echo -e "${RED}Error: VASP_JOBS directory not found: $VASP_JOBS_DIR${NC}"
    echo ""
    echo "Please provide valid path to original VASP_JOBS directory"
    echo "Use: --vasp-jobs-dir /path/to/VASP_JOBS/"
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
echo "  Database: $DB_PATH"
echo "  Conda env: $CONDA_ENV"
echo ""
echo -e "${YELLOW}Note: Structures (CONTCARs) will be loaded from original VASP_JOBS directory${NC}"
echo ""

# Export variables for SLURM script
export INPUT_FILE
export VASP_JOBS_DIR
export OUTPUT_DIR
export MAX_CONCURRENT
export MAX_STRUCTURES
export CHECK_INTERVAL
export DB_NAME
export CONDA_ENV

# Submit the workflow manager job
echo "Submitting refined workflow manager as SLURM job..."
JOBID=$(sbatch submit_refineflow.sh | awk '{print $NF}')

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Refined workflow manager submitted successfully!${NC}"
    echo ""
    echo "Job ID: $JOBID"
    echo ""
    echo "Monitoring commands:"
    echo "  View log:        tail -f refine_workflow_${JOBID}.out"
    echo "  Check status:    squeue -j $JOBID"
    echo "  Check workflow:  python workflow_status.py $OUTPUT_DIR/workflow.json"
    echo "  Cancel job:      scancel $JOBID"
    echo ""
    echo "The refined workflow will:"
    echo "  - Run on identified electride candidates only"
    echo "  - Use higher k-point density (250)"
    echo "  - Use more electronic/ionic steps"
    echo "  - Save results to $OUTPUT_DIR"
    echo ""
else
    echo -e "${RED}Error: Failed to submit job${NC}"
    exit 1
fi

