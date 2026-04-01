#!/bin/bash
# VASPflow Pre-screening - Submits MatterSim pre-screening as SLURM job
# Usage: bash run_prescreen.sh [options]

set -e

# Default values
RESULTS_DIR="../results/quaternary_boron_hydride/"
OUTPUT_DIR="./VASP_JOBS_Boron/"
MAX_COMPOSITIONS=""
MAX_STRUCTURES=0
CONDA_ENV="mattersim"
MP_API_KEY="${MP_API_KEY:-}"
HULL_THRESHOLD=0.1
DEVICE="cuda"
BATCH_SIZE=32
MAX_ATOMS_GPU=2048
PURE_PBE="--pure-pbe"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo "========================================================================"
echo "VASPflow - MatterSim Pre-screening"
echo "========================================================================"
echo ""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --results-dir)
            RESULTS_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --max-compositions)
            MAX_COMPOSITIONS="$2"
            shift 2
            ;;
        --max-structures)
            MAX_STRUCTURES="$2"
            shift 2
            ;;
        --mp-api-key)
            MP_API_KEY="$2"
            shift 2
            ;;
        --hull-threshold)
            HULL_THRESHOLD="$2"
            shift 2
            ;;
        --device)
            DEVICE="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --max-atoms-gpu)
            MAX_ATOMS_GPU="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --pure-pbe)
            PURE_PBE="--pure-pbe"
            shift
            ;;
        --help)
            echo "Usage: bash run_prescreen.sh [options]"
            echo ""
            echo "This script submits pre-screening as a SLURM job."
            echo "Pre-screening uses MatterSim to filter structures by thermodynamic stability."
            echo ""
            echo "Options:"
            echo "  --results-dir DIR          MatterGen results directory"
            echo "  --output-dir DIR           Output directory for results"
            echo "  --max-compositions N       Max compositions to process (default: all)"
            echo "  --max-structures N         Max structures per composition (default: 5, use 0 for all)"
            echo "  --mp-api-key KEY           Materials Project API key"
            echo "  --hull-threshold FLOAT     Energy above hull threshold in eV/atom (default: 0.1)"
            echo "  --device DEVICE            MatterSim device: cpu or cuda (default: cuda)"
            echo "                             Note: GPU is HIGHLY RECOMMENDED (10-50x faster than CPU)"
            echo "  --batch-size N             Batch size for GPU parallel relaxation (default: 32)"
            echo "  --max-atoms-gpu N          Max total atoms on GPU simultaneously (default: 2048)"
            echo "                             Adjust for GPU memory: 2048 (V100 16GB), 4096 (A100 40GB),"
            echo "                             8192 (A100 80GB / H100)"
            echo "  --conda-env NAME           Conda environment name (default: mattersim)"
            echo "  --pure-pbe                 Filter MP entries to pure GGA-PBE only (exclude PBE+U)"
            echo "                             Default: accept both PBE and PBE+U (recommended)"
            echo ""
            echo "Example:"
            echo "  bash run_prescreen.sh --max-structures 10 --hull-threshold 0.05"
            echo "  bash run_prescreen.sh --device cpu --mp-api-key YOUR_KEY  # Use CPU if no GPU"
            echo ""
            echo "Monitoring:"
            echo "  Check status:    squeue -u \$USER | grep prescreen"
            echo "  View log:        tail -f prescreen_<JOBID>.out"
            echo "  Check results:   cat ./VASP_JOBS/prescreening_stability.json | jq '.summary'"
            echo ""
            exit 0
            ;;
        *)
            echo -e "${RED}Error: Unknown option $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate inputs
if [ ! -d "$RESULTS_DIR" ]; then
    echo -e "${RED}Error: Results directory not found: $RESULTS_DIR${NC}"
    exit 1
fi

if [ -z "$MP_API_KEY" ]; then
    echo -e "${YELLOW}Warning: MP_API_KEY not set${NC}"
    echo "  Pre-screening requires Materials Project API key"
    echo "  Get your free key at: https://next-gen.materialsproject.org/api"
    echo "  Set it with: export MP_API_KEY=your_key"
    echo "  Or pass it with: --mp-api-key YOUR_KEY"
    echo ""
    read -p "Continue without API key? (y/N) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Print configuration
echo "Configuration:"
echo "  Results directory:    $RESULTS_DIR"
echo "  Output directory:     $OUTPUT_DIR"
echo "  Max compositions:     ${MAX_COMPOSITIONS:-all}"
echo "  Max structures:       ${MAX_STRUCTURES} (0 = all)"
echo "  Hull threshold:       ${HULL_THRESHOLD} eV/atom"
echo "  Device:               $DEVICE"
echo "  Batch size:           $BATCH_SIZE"
echo "  Max atoms on GPU:     $MAX_ATOMS_GPU"
echo "  Conda environment:    $CONDA_ENV"
echo "  MP API key:           ${MP_API_KEY:+[SET]}${MP_API_KEY:-[NOT SET]}"
echo "  Functional filter:    ${PURE_PBE:+Pure PBE only}${PURE_PBE:-Mixed PBE/PBE+U (default)}"
echo ""

# Check if database already exists
if [ -f "$OUTPUT_DIR/prescreening_stability.json" ]; then
    echo -e "${YELLOW}Warning: Pre-screening results already exist:${NC}"
    echo "  $OUTPUT_DIR/prescreening_stability.json"
    echo ""
    read -p "Overwrite existing results? (y/N) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborting."
        exit 0
    fi
fi

# Submit job
echo "Submitting pre-screening job to SLURM..."
echo ""

export RESULTS_DIR
export OUTPUT_DIR
export MAX_COMPOSITIONS
export MAX_STRUCTURES
export CONDA_ENV
export MP_API_KEY
export HULL_THRESHOLD
export DEVICE
export BATCH_SIZE
export MAX_ATOMS_GPU
export PURE_PBE

JOB_ID=$(sbatch submit_prescreen.sh | awk '{print $NF}')

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Pre-screening job submitted successfully!${NC}"
    echo ""
    echo "Job ID: $JOB_ID"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u $USER | grep prescreen"
    echo "  tail -f prescreen_${JOB_ID}.out"
    echo ""
    echo "When complete, check results:"
    echo "  cat $OUTPUT_DIR/prescreening_stability.json | jq '.summary'"
    echo ""
    echo "Then proceed with VASP workflow:"
    echo "  bash run_workflow.sh --max-concurrent 20"
    echo ""
else
    echo -e "${RED}Error: Failed to submit job${NC}"
    exit 1
fi

