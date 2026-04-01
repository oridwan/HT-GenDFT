#!/bin/bash
# Submit multiple parallel prescreening jobs across multiple partitions
# Each job processes a different subset of compositions and partitions are assigned round-robin.
# change RESULTS_DIR, OUTPUT_DIR

set -euo pipefail

# Default configuration
COMPOSITIONS_PER_JOB=500
RESULTS_DIR="../results/quaternary_boron_hydride"
OUTPUT_DIR="./VASP_JOBS_Boron"
BATCH_SIZE=32
HULL_THRESHOLD=0.1
DEVICE="cuda"
MAX_STRUCTURES=0
PURE_PBE=""
PARTITIONS="GPU,Hydrus"
CPUS_PER_TASK=2
MEMORY="32G"
TIME_LIMIT="10-00:00:00"
GPUS_PER_JOB=1
CONDA_ENV="mattersim"
MP_API_KEY="${MP_API_KEY:-}"

# Parse command-line arguments
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
        --batch-size)
            BATCH_SIZE="$2"
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
        --compositions-per-job)
            COMPOSITIONS_PER_JOB="$2"
            shift 2
            ;;
        --max-structures)
            MAX_STRUCTURES="$2"
            shift 2
            ;;
        --partitions)
            PARTITIONS="$2"
            shift 2
            ;;
        --cpus-per-task)
            CPUS_PER_TASK="$2"
            shift 2
            ;;
        --mem)
            MEMORY="$2"
            shift 2
            ;;
        --time)
            TIME_LIMIT="$2"
            shift 2
            ;;
        --gpus)
            GPUS_PER_JOB="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --mp-api-key)
            MP_API_KEY="$2"
            shift 2
            ;;
        --pure-pbe)
            PURE_PBE="--pure-pbe"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Submit multiple parallel prescreening jobs across SLURM partitions"
            echo ""
            echo "Options:"
            echo "  --results-dir DIR          Path to MatterGen results (default: ../results/quaternary_hydride)"
            echo "  --output-dir DIR           Path to output directory (default: ./VASP_JOBS)"
            echo "  --batch-size N             Structures per batch in Python (default: 32)"
            echo "  --hull-threshold FLOAT     Energy above hull threshold in eV/atom (default: 0.1)"
            echo "  --device DEVICE            Device: cpu or cuda (default: cuda)"
            echo "  --compositions-per-job N   Compositions per SLURM job (default: 300)"
            echo "  --max-structures N         Max structures per composition (default: 0 = all)"
            echo "  --partitions LIST          Comma-separated partitions (default: GPU,Hydrus)"
            echo "  --cpus-per-task N          CPUs per SLURM job (default: 8)"
            echo "  --mem MEM                  Memory per SLURM job (default: 64G)"
            echo "  --time LIMIT               Time per SLURM job (default: 3-00:00:00)"
            echo "  --gpus N                   GPUs per job, generic request (default: 1)"
            echo "  --conda-env NAME           Conda environment to activate (default: mattersim)"
            echo "  --mp-api-key KEY           Materials Project API key (or use MP_API_KEY env var)"
            echo "  --pure-pbe                 Filter MP entries to pure GGA-PBE only (exclude PBE+U)"
            echo "  -h, --help                 Show this help message"
            echo ""
            echo "Example:"
            echo "  $0 --results-dir ../results/quaternary_hydride --output-dir ./VASP_JOBS --partitions GPU,Hydrus --compositions-per-job 250"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

if [ ! -d "$RESULTS_DIR" ]; then
    echo "ERROR: Results directory not found: $RESULTS_DIR"
    exit 1
fi

if [ -z "$MP_API_KEY" ]; then
    echo "ERROR: MP_API_KEY is required"
    echo "Set env var: export MP_API_KEY=..."
    echo "Or pass: --mp-api-key ..."
    exit 1
fi

IFS=',' read -r -a PARTITION_ARR <<< "$PARTITIONS"
if [ ${#PARTITION_ARR[@]} -eq 0 ]; then
    echo "ERROR: No partitions provided"
    exit 1
fi

# Count all composition directories and those with CIF zip
TOTAL_COMPOSITIONS=$(find "$RESULTS_DIR" -maxdepth 1 -type d -name "*_structures" | wc -l)
WITH_ZIP=$(find "$RESULTS_DIR" -maxdepth 2 -type f -path "*/generated_crystals_cif.zip" | wc -l)

if [ "$TOTAL_COMPOSITIONS" -eq 0 ]; then
    echo "ERROR: No *_structures directories found in $RESULTS_DIR"
    exit 1
fi

NUM_JOBS=$(( (TOTAL_COMPOSITIONS + COMPOSITIONS_PER_JOB - 1) / COMPOSITIONS_PER_JOB ))
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "Parallel Prescreening Job Submission"
echo "========================================"
echo "Results directory: $RESULTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Total composition dirs: $TOTAL_COMPOSITIONS"
echo "Dirs with CIF zip: $WITH_ZIP"
echo "Compositions per job: $COMPOSITIONS_PER_JOB"
echo "Number of jobs: $NUM_JOBS"
echo "Partitions: $PARTITIONS"
echo "CPUs per task: $CPUS_PER_TASK"
echo "Memory: $MEMORY"
echo "Time limit: $TIME_LIMIT"
echo "GPUs per job: $GPUS_PER_JOB"
echo "Device: $DEVICE"
echo "Batch size: $BATCH_SIZE"
echo "Functional filter: ${PURE_PBE:+Pure PBE only}${PURE_PBE:-Mixed PBE/PBE+U}"
echo "========================================"
echo ""

for i in $(seq 0 $((NUM_JOBS - 1))); do
    START_IDX=$((i * COMPOSITIONS_PER_JOB))

    p_idx=$((i % ${#PARTITION_ARR[@]}))
    PARTITION="${PARTITION_ARR[$p_idx]}"

    BATCH_OUTPUT_FILE="${OUTPUT_DIR}/prescreening_stability_batch${i}.json"
    BATCH_CHECKPOINT_FILE="${OUTPUT_DIR}/prescreening_checkpoint_batch${i}.json"

    # Skip fully completed batches (output exists and no active checkpoint).
    # If checkpoint exists, submit to allow resume.
    if [ -s "$BATCH_OUTPUT_FILE" ] && [ ! -f "$BATCH_CHECKPOINT_FILE" ]; then
        echo "Skipping batch $i on $PARTITION: found completed output $BATCH_OUTPUT_FILE"
        continue
    fi

    TEMP_SCRIPT="submit_prescreen_batch${i}.sh"

    cat > "$TEMP_SCRIPT" <<EOF_HEADER
#!/bin/bash
#SBATCH --job-name=ps_b${i}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --gres=gpu:${GPUS_PER_JOB}
#SBATCH --mem=${MEMORY}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=prescreen_batch${i}_%j.out
#SBATCH --error=prescreen_batch${i}_%j.err

set -euo pipefail

# Some cluster bashrc files reference vars before defining them.
# Temporarily disable nounset while sourcing shell init files.
set +u
source ~/.bashrc
set -u
conda activate ${CONDA_ENV}

if [ \$? -ne 0 ]; then
    echo "ERROR: Failed to activate conda env: ${CONDA_ENV}"
    exit 1
fi

if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader || true
fi

export OMP_NUM_THREADS=${CPUS_PER_TASK}
export MKL_NUM_THREADS=${CPUS_PER_TASK}
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export MP_API_KEY="${MP_API_KEY}"

START_INDEX=${START_IDX}
MAX_COMPS=${COMPOSITIONS_PER_JOB}
BATCH_ID=${i}
RESULTS_DIR="${RESULTS_DIR}"
OUTPUT_DIR="${OUTPUT_DIR}"
BATCH_SIZE=${BATCH_SIZE}
HULL_THRESHOLD=${HULL_THRESHOLD}
MAX_STRUCTURES=${MAX_STRUCTURES}
DEVICE="${DEVICE}"
PURE_PBE="${PURE_PBE}"

echo "========================================"
echo "Prescreen batch \$BATCH_ID on partition ${PARTITION}"
echo "Composition range: \$START_INDEX - \$((START_INDEX + MAX_COMPS - 1))"
echo "========================================"

python3 prescreen.py \
    --results-dir "\$RESULTS_DIR" \
    --output-dir "\$OUTPUT_DIR" \
    --hull-threshold "\$HULL_THRESHOLD" \
    --device "\$DEVICE" \
    --batch-size "\$BATCH_SIZE" \
    --start-composition "\$START_INDEX" \
    --max-compositions "\$MAX_COMPS" \
    --max-structures "\$MAX_STRUCTURES" \
    --batch-id "\$BATCH_ID" \
    \$PURE_PBE
EOF_HEADER

    JOB_ID=$(sbatch "$TEMP_SCRIPT" | awk '{print $NF}')
    echo "Submitted batch $i to partition $PARTITION (compositions $START_IDX-$((START_IDX + COMPOSITIONS_PER_JOB - 1))): Job ID $JOB_ID"

    rm -f "$TEMP_SCRIPT"
    sleep 300

done

echo ""
echo "========================================"
echo "All jobs submitted"
echo "Monitor: squeue -u \$USER | grep ps_b"
echo "Progress files: ls ${OUTPUT_DIR}/prescreening_checkpoint_batch*.json"
echo "Result files: ls ${OUTPUT_DIR}/prescreening_stability_batch*.json"
echo "DB files: ls ${OUTPUT_DIR}/prescreening_structures_batch*.db"
echo "Merge when done: python3 merge_prescreen_batches.py --output-dir ${OUTPUT_DIR}"
echo "========================================"
