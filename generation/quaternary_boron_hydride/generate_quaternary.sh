#!/bin/bash
# FOR NEW SYSTEM: Change only these 3 items:
# 1. Line 13: --job-name (e.g., gen_ternary_oxide_csp)
# 2. Line 33: COMPOSITIONS_FILE (your new compositions JSON file)
# 3. Line 34: OUTPUT_DIR (your results directory)
# 4. START_INDEX ,MAX_COMPOSITIONS 
# 5. working directory (line 35)
# Then run: sbatch generate_quaternary.sh

#SBATCH --job-name=gen_qua_hyd_csp
#SBATCH --partition=Hydrus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/gen_qua_hyd_csp_%j.out
#SBATCH --error=logs/gen_qua_hyd_csp_%j.err

# Load modules
module purge
module load cuda/11.8

PROJECT_ROOT=/projects/mmi/Ridwan/Github/SuperConductorFlow

# Activate conda environment (non-interactive shell safe)
if [ -f "$HOME/.bashrc" ]; then
    # shellcheck disable=SC1090
    source "$HOME/.bashrc"
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found in PATH"
    exit 1
fi
eval "$(conda shell.bash hook)"
conda activate mattersim
if [ ! -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
    echo "ERROR: Python venv not found: $PROJECT_ROOT/.venv/bin/activate"
    exit 1
fi
source "$PROJECT_ROOT/.venv/bin/activate"

# Navigate to working directory
cd "$PROJECT_ROOT/generation/quaternary_boron_hydride"

# Create logs directory if it doesn't exist
mkdir -p logs

# Memory optimization
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# Generation parameters
CSP_MODEL="$PROJECT_ROOT/generation/MatterGen_checkpoints/18-55-08"
COMPOSITIONS_FILE="quaternary_boron_hydride_compositions.json"
OUTPUT_DIR="$PROJECT_ROOT/results/quaternary_boron_hydride"
STRUCTURES_PER_ATOM=2.0
MAX_BATCH_SIZE=100  # Maximum structures per GPU batch (splits larger requests to avoid OOM)
TIMEOUT_PER_BATCH=1800  # Base timeout in seconds per batch (30 minutes) - auto-scaled for multiple batches
MAX_COMPOSITIONS=-1  # TEST: Changed from -1 (all) to 10 (just first compositions) 6388
START_INDEX=0
BASE_START_INDEX="${BASE_START_INDEX:-$START_INDEX}"
GLOBAL_END_INDEX="${GLOBAL_END_INDEX:-0}"

# If running as a Slurm array, compute task-local bounds from exported base offset
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    if [ -z "$COMPOSITIONS_PER_JOB" ]; then
        echo "ERROR: COMPOSITIONS_PER_JOB is not set for array task $SLURM_ARRAY_TASK_ID"
        exit 1
    fi
    MAX_COMPOSITIONS=$COMPOSITIONS_PER_JOB
    START_INDEX=$((BASE_START_INDEX + SLURM_ARRAY_TASK_ID * MAX_COMPOSITIONS))

    if [ "$GLOBAL_END_INDEX" -gt 0 ]; then
        REMAINING=$((GLOBAL_END_INDEX - START_INDEX))
        if [ "$REMAINING" -le 0 ]; then
            echo "No compositions assigned to array task $SLURM_ARRAY_TASK_ID (start=$START_INDEX, end=$GLOBAL_END_INDEX). Exiting."
            exit 0
        fi
        if [ "$REMAINING" -lt "$MAX_COMPOSITIONS" ]; then
            MAX_COMPOSITIONS=$REMAINING
        fi
    fi
fi

echo "=========================================="
echo "Quaternary Hydride Generation (CSP Mode)"
echo "=========================================="
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Array task ID: $SLURM_ARRAY_TASK_ID"
    echo "Base start index: $BASE_START_INDEX"
    if [ "$GLOBAL_END_INDEX" -gt 0 ]; then
        echo "Global end index: $GLOBAL_END_INDEX"
    fi
fi
echo "CSP Model: $CSP_MODEL"
echo "Compositions: $COMPOSITIONS_FILE"
echo "Output: $OUTPUT_DIR"
echo "Structures per atom: $STRUCTURES_PER_ATOM"
echo "Max batch size: $MAX_BATCH_SIZE (auto-split for large requests)"
echo "Max compositions: $MAX_COMPOSITIONS"
echo "Start index: $START_INDEX"
echo "=========================================="

# Check if CSP model exists
if [ ! -d "$CSP_MODEL" ]; then
    echo "ERROR: CSP model directory not found: $CSP_MODEL"
    echo ""
    echo "Please update CSP_MODEL path in this script with your fine-tuned CSP checkpoint"
    echo "Example: outputs/singlerun/2025-10-16/18-55-08"
    exit 1
fi

# Check if compositions file exists
if [ ! -f "$COMPOSITIONS_FILE" ]; then
    echo "ERROR: Compositions file not found: $COMPOSITIONS_FILE"
    echo ""
    echo "Please run search_quaternary_electrides.py first to generate compositions"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run generation
# The script will skip compositions that already have generated structures
# This allows you to resubmit the job if it times out
python generate_structures_batch.py \
    --compositions "$COMPOSITIONS_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --model "$CSP_MODEL" \
    --structures-per-atom $STRUCTURES_PER_ATOM \
    --max-batch-size $MAX_BATCH_SIZE \
    --timeout $TIMEOUT_PER_BATCH \
    --max-compositions $MAX_COMPOSITIONS \
    --start-index $START_INDEX

echo ""
echo " Structure generation completed"

# ============================================================
# STEP 2: Create Summary
# ============================================================
echo ""
echo "STEP 2: Creating summary..."
echo "------------------------------------------------------------"

# Call the Python script that creates summary statistics
python summarize_generation.py \
    --output-dir "$OUTPUT_DIR"

SUMMARY_STATUS=$?

if [ $SUMMARY_STATUS -ne 0 ]; then
    echo ""
    echo "WARNING: Summary generation failed with exit code $SUMMARY_STATUS"
fi

echo ""
echo " Summary created"

# ============================================================
# COMPLETION
# ============================================================
echo ""
echo "=========================================="
echo "Quaternary hydride generation completed!"
echo "=========================================="
echo ""
echo "Results location: $OUTPUT_DIR/"
echo ""
echo "Generated files:"
echo "  • Structure directories: ${OUTPUT_DIR}/*_structures/"
echo "  • Generation statistics: ${OUTPUT_DIR}/generation_statistics.json"
echo "  • Generation summary: ${OUTPUT_DIR}/generation_summary.json"
if [ -f "${OUTPUT_DIR}/failed_compositions.txt" ]; then
    echo "  • Failed compositions: ${OUTPUT_DIR}/failed_compositions.txt"
fi
echo ""
echo "Next steps:"
echo "  1. Review generation_summary.json"
echo "  2. Evaluate structures: sbatch ../mattergen_evaluate.sh"
echo "  3. Filter by stability and electride probability"
echo ""
echo "=========================================="
