#!/bin/bash
# Electronic Structure Workflow - Submits workflow manager as SLURM job
# Usage: bash run_electronic_flow.sh [options]

set -e

# Default values
REFINE_JOBS="./REFINE_VASP_JOBS"
OUTPUT_DIR="./ELECTRONIC_JOBS"
MAX_CONCURRENT=10
CHECK_INTERVAL=60
DB_NAME="workflow.json"
CONDA_ENV="vaspflow"
STRUCTURE_IDS=""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "========================================================================"
echo "Electronic Structure Workflow - SCF/PARCHG/ELF/BAND/PDOS"
echo "========================================================================"
echo ""
echo "Sequential workflow per structure:"
echo "  1. SCF (Self-Consistent Field) - conventional cell"
echo "  2. PARCHG (Partial Charge Density) - 5 windows, conventional cell"
echo "  3. ELF (Electron Localization Function) - conventional cell"
echo "  4. BAND (Band Structure) - primitive cell"
echo "  5. PDOS (Projected Density of States) - primitive cell"
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
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --help)
            echo "Usage: bash run_electronic_flow.sh [options]"
            echo ""
            echo "This script submits the electronic structure workflow manager as a SLURM job."
            echo "It runs SCF/PARCHG/ELF/BAND/PDOS calculations on refined electride structures."
            echo ""
            echo "Options:"
            echo "  --refine-jobs DIR          Path to REFINE_VASP_JOBS directory (default: ./REFINE_VASP_JOBS)"
            echo "  --output-dir DIR           Output directory for electronic jobs (default: ./ELECTRONIC_JOBS)"
            echo "  --structure-ids ID1 ID2... Specific structure IDs to process (default: all RELAX_DONE)"
            echo "  --max-concurrent N         Max concurrent structures (default: 10)"
            echo "  --check-interval SECONDS   Status check interval (default: 60)"
            echo "  --conda-env NAME           Conda environment name (default: vaspflow)"
            echo ""
            echo "Examples:"
            echo "  # Process all RELAX_DONE structures"
            echo "  bash run_electronic_flow.sh"
            echo ""
            echo "  # Process specific structures"
            echo "  bash run_electronic_flow.sh --structure-ids Ba2N_s001 Ca5P3_s002 Li3N_s003"
            echo ""
            echo "  # Custom settings"
            echo "  bash run_electronic_flow.sh --max-concurrent 5 --refine-jobs /path/to/REFINE_VASP_JOBS"
            echo ""
            echo "Prerequisites:"
            echo "  1. Complete refined relaxation workflow (run_refineflow.sh)"
            echo "  2. Structures must be in RELAX_DONE state"
            echo ""
            echo "Monitoring:"
            echo "  Check status:    python workflow_status.py ./ELECTRONIC_JOBS/workflow.json"
            echo "  View log:        tail -f electronic_workflow_<JOBID>.out"
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

# Count RELAX_DONE structures
RELAX_DONE_COUNT=$(grep -c '"state": "RELAX_DONE"' "$REFINE_DB" || true)
echo -e "${GREEN}Found $RELAX_DONE_COUNT structures in RELAX_DONE state${NC}"
echo ""

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
if [ -n "$STRUCTURE_IDS" ]; then
    echo "  Structure IDs:$STRUCTURE_IDS"
else
    echo "  Structure IDs: all RELAX_DONE structures"
fi
echo "  Check interval: ${CHECK_INTERVAL}s"
echo "  Database: $DB_PATH"
echo "  Conda env: $CONDA_ENV"
echo ""

# Export variables for SLURM script
export REFINE_JOBS
export OUTPUT_DIR
export MAX_CONCURRENT
export CHECK_INTERVAL
export DB_NAME
export CONDA_ENV
export STRUCTURE_IDS

# Submit the workflow manager job
echo "Submitting electronic workflow manager as SLURM job..."
JOBID=$(sbatch submit_electronic_flow.sh | awk '{print $NF}')

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Electronic workflow manager submitted successfully!${NC}"
    echo ""
    echo "Job ID: $JOBID"
    echo ""
    echo "Monitoring commands:"
    echo "  View log:        tail -f electronic_workflow_${JOBID}.out"
    echo "  Check status:    squeue -j $JOBID"
    echo "  Check workflow:  python workflow_status.py $OUTPUT_DIR/workflow.json"
    echo "  Cancel job:      scancel $JOBID"
    echo ""
    echo "The workflow will:"
    echo "  - Run SCF/PARCHG/ELF on conventional cell (visualization)"
    echo "  - Run BAND/PDOS on primitive cell (computational efficiency)"
    echo "  - Save results to $OUTPUT_DIR"
    echo ""
else
    echo -e "${RED}Error: Failed to submit job${NC}"
    exit 1
fi
