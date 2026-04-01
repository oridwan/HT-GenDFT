#!/bin/bash
# Submit multiple parallel prescreening jobs across multiple GPU nodes
# Each job processes a different subset of compositions

# Default configuration
COMPOSITIONS_PER_JOB=100
RESULTS_DIR="../results/quaternary_boron_hydride"
OUTPUT_DIR="./VASP_JOBS_Boron"
BATCH_SIZE=32
HULL_THRESHOLD=0.1
DEVICE="cuda"
MAX_STRUCTURES=0
MAX_ATOMS_GPU=2048
PURE_PBE="--pure-pbe"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRESCREEN_SCRIPT="$SCRIPT_DIR/prescreen.py"

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
        --max-atoms-gpu)
            MAX_ATOMS_GPU="$2"
            shift 2
            ;;
        --pure-pbe)
            PURE_PBE="--pure-pbe"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Submit multiple parallel prescreening jobs across GPU nodes"
            echo ""
            echo "Options:"
            echo "  --results-dir DIR          Path to MatterGen results (default: ./mattergen_results/ternary_csp_electrides)"
            echo "  --output-dir DIR           Path to output directory (default: ./VASP_JOBS)"
            echo "  --batch-size N             Batch size for GPU parallel relaxation (default: 32)"
            echo "  --max-atoms-gpu N          Max total atoms on GPU simultaneously (default: 2048)"
            echo "                             Adjust for GPU: 2048 (V100 16GB), 4096 (A100 40GB), 8192 (A100 80GB/H100)"
            echo "  --hull-threshold FLOAT     Energy above hull threshold in eV/atom (default: 0.1)"
            echo "  --device DEVICE            Device: cpu or cuda (default: cuda)"
            echo "  --compositions-per-job N   Compositions per parallel job (default: 400)"
            echo "  --max-structures N         Max structures per composition (default: 5, 0=all)"
            echo "  --pure-pbe                 Filter MP entries to pure GGA-PBE only (exclude PBE+U) [default]"
            echo "  -h, --help                 Show this help message"
            echo ""
            echo "Example:"
            echo "  $0 --results-dir ./results --batch-size 64 --max-structures 10 --compositions-per-job 200"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Auto-detect total compositions from results directory
if [ ! -d "$RESULTS_DIR" ]; then
    echo "ERROR: Results directory not found: $RESULTS_DIR"
    echo "Please set RESULTS_DIR to your MatterGen results directory"
    exit 1
fi

if [ ! -f "$PRESCREEN_SCRIPT" ]; then
    echo "ERROR: prescreen.py not found: $PRESCREEN_SCRIPT"
    exit 1
fi

# Count composition directories
TOTAL_COMPOSITIONS=$(find "$RESULTS_DIR" -maxdepth 1 -type d -name "*_structures" | wc -l)

if [ $TOTAL_COMPOSITIONS -eq 0 ]; then
    echo "ERROR: No *_structures directories found in $RESULTS_DIR"
    exit 1
fi

echo "Found $TOTAL_COMPOSITIONS compositions in $RESULTS_DIR"

# Calculate number of jobs needed
NUM_JOBS=$(( (TOTAL_COMPOSITIONS + COMPOSITIONS_PER_JOB - 1) / COMPOSITIONS_PER_JOB ))

echo "========================================"
echo "Parallel Prescreening Job Submission"
echo "========================================"
echo "Results directory: $RESULTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Total compositions: $TOTAL_COMPOSITIONS"
echo "Compositions per job: $COMPOSITIONS_PER_JOB"
echo "Number of parallel jobs: $NUM_JOBS"
echo "Batch size: $BATCH_SIZE"
echo "Max atoms on GPU: $MAX_ATOMS_GPU"
echo "Hull threshold: $HULL_THRESHOLD eV/atom"
echo "Max structures per composition: $MAX_STRUCTURES"
echo "Device: $DEVICE"
echo "Functional filter: ${PURE_PBE:+Pure PBE only}${PURE_PBE:-Mixed PBE/PBE+U}"
echo "========================================"
echo ""

# Create temporary SLURM scripts for each batch
for i in $(seq 0 $((NUM_JOBS - 1))); do
    START_IDX=$((i * COMPOSITIONS_PER_JOB))
    
    TEMP_SCRIPT="submit_prescreen_batch${i}.sh"
    
    # Create SLURM script for this batch with batch-specific output files
    cat > $TEMP_SCRIPT <<EOF_HEADER
#!/bin/bash
#SBATCH --job-name=prescreen_b${i}
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=prescreen_batch${i}_%j.out
#SBATCH --error=prescreen_batch${i}_%j.err

# Load conda
eval "\$(conda shell.bash hook)"
conda activate mattersim

if [ \$? -ne 0 ]; then
    echo "ERROR: Failed to activate mattersim conda environment"
    exit 1
fi

# Load CUDA module for shared libraries (required for PyTorch)
module load cuda/11.8

# Print GPU info
nvidia-smi

# Set CUDA memory configuration for better GPU memory management
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

# Get batch parameters
EOF_HEADER
    
    # Add batch-specific parameters
    cat >> $TEMP_SCRIPT <<EOF_PARAMS
START_INDEX=$START_IDX
MAX_COMPS=$COMPOSITIONS_PER_JOB
BATCH_ID=$i
RESULTS_DIR="$RESULTS_DIR"
OUTPUT_DIR="$OUTPUT_DIR"
BATCH_SIZE=$BATCH_SIZE
MAX_ATOMS_GPU=$MAX_ATOMS_GPU
HULL_THRESHOLD=$HULL_THRESHOLD
MAX_STRUCTURES=$MAX_STRUCTURES
DEVICE="$DEVICE"
PURE_PBE="$PURE_PBE"
PRESCREEN_SCRIPT="$PRESCREEN_SCRIPT"

EOF_PARAMS
    
    # Add execution command
    cat >> $TEMP_SCRIPT << 'EOF_FOOTER'
echo "========================================"
echo "Prescreening Batch Job"
echo "========================================"
echo "Batch ID: $BATCH_ID"
echo "Start composition index: $START_INDEX"
echo "Max compositions: $MAX_COMPS"
echo "Composition range: $START_INDEX to $((START_INDEX + MAX_COMPS - 1))"
echo "Results directory: $RESULTS_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Batch size: $BATCH_SIZE"
echo "Max atoms on GPU: $MAX_ATOMS_GPU"
echo "Hull threshold: $HULL_THRESHOLD eV/atom"
echo "Max structures per composition: $MAX_STRUCTURES"
echo "Device: $DEVICE"
echo "Script: $PRESCREEN_SCRIPT"
echo "========================================"

# Run prescreening
python3 "$PRESCREEN_SCRIPT" \
    --results-dir "$RESULTS_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --hull-threshold $HULL_THRESHOLD \
    --device $DEVICE \
    --batch-size $BATCH_SIZE \
    --max-atoms-gpu $MAX_ATOMS_GPU \
    --start-composition $START_INDEX \
    --max-compositions $MAX_COMPS \
    --max-structures $MAX_STRUCTURES \
    --batch-id $BATCH_ID \
    $PURE_PBE

EXIT_CODE=$?

echo ""
echo "========================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Batch completed successfully"
else
    echo "Batch failed with exit code $EXIT_CODE"
fi
echo "========================================"

exit $EXIT_CODE
EOF_FOOTER
    
    # Submit the job
    JOB_ID=$(sbatch $TEMP_SCRIPT | awk '{print $NF}')
    
    echo "Submitted batch $i (compositions $START_IDX-$((START_IDX + COMPOSITIONS_PER_JOB - 1))): Job ID $JOB_ID"
    
    # Clean up temporary script
    rm $TEMP_SCRIPT
    
    # Small delay to avoid overwhelming the scheduler
    sleep 0.5
done

echo ""
echo "========================================"
echo "All jobs submitted!"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check progress: ls $OUTPUT_DIR/prescreening_checkpoint_batch*.json"
echo "View results: ls $OUTPUT_DIR/prescreening_stability_batch*.json"
echo "View databases: ls $OUTPUT_DIR/prescreening_structures_batch*.db"
echo "View shared MP cache: ls $OUTPUT_DIR/mp_mattersim.json"
echo ""
echo "After all jobs complete, merge results with:"
echo "  python3 merge_prescreen_batches.py --output-dir $OUTPUT_DIR"
echo "  (merges JSON and database; MP cache is already shared)"
echo "========================================"
