#!/bin/bash
#SBATCH --job-name=prescreen
#SBATCH --partition=Hydrus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=10-00:00:00
#SBATCH --output=prescreen_%j.out
#SBATCH --error=prescreen_%j.err

# VASPflow Pre-screening with MatterSim
# Fast thermodynamic stability screening before expensive VASP calculations

set -e

# Default parameters
RESULTS_DIR=${RESULTS_DIR:-"./mattergen_results/ternary_csp_electrides"}
OUTPUT_DIR=${OUTPUT_DIR:-"./VASP_JOBS"}
MAX_COMPOSITIONS=${MAX_COMPOSITIONS:-""}
MAX_STRUCTURES=${MAX_STRUCTURES:-0}
CONDA_ENV=${CONDA_ENV:-"mattersim"}
MP_API_KEY=${MP_API_KEY:-""}
HULL_THRESHOLD=${HULL_THRESHOLD:-0.1}
DEVICE=${DEVICE:-"cuda"}
BATCH_SIZE=${BATCH_SIZE:-32}
PURE_PBE=${PURE_PBE:-""}

echo "========================================================================"
echo "VASPflow Pre-screening (MatterSim)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "GPUs: ${SLURM_GPUS:-0} (${CUDA_VISIBLE_DEVICES:-none})"
echo "Memory: 128 GB"
echo "Start time: $(date)"
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

# Check GPU availability and auto-detect device
if command -v nvidia-smi &> /dev/null; then
    echo "GPU Information:"
    echo "========================================================================"
    nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader
    
    # Get GPU memory info
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
    GPU_MEM_TOTAL=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1)
    GPU_MEM_FREE=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits | head -1)
    
    echo ""
    echo "Allocated GPU: $GPU_NAME"
    echo "Total VRAM: $GPU_MEM_TOTAL MiB (~$((GPU_MEM_TOTAL/1024)) GB)"
    echo "Free VRAM: $GPU_MEM_FREE MiB (~$((GPU_MEM_FREE/1024)) GB)"
    echo "========================================================================"
    echo ""
    
    # Auto-detect: use cuda if GPU available, otherwise fallback to cpu
    if [ "$DEVICE" = "cuda" ]; then
        if ! python3 -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
            echo "Warning: CUDA requested but not available in PyTorch, falling back to CPU"
            DEVICE="cpu"
        else
            echo "  Using CUDA device for MatterSim"
            echo "  PyTorch detects GPU with full VRAM access"
        fi
    fi
else
    echo "No GPU detected (nvidia-smi not found)"
    if [ "$DEVICE" = "cuda" ]; then
        echo "Warning: CUDA requested but no GPU available, falling back to CPU"
        DEVICE="cpu"
    fi
fi
echo ""

# Set number of threads for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# PyTorch CUDA memory allocator settings (prevents fragmentation and OOM)
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "Computation device: $DEVICE"
echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo "  PYTORCH_CUDA_ALLOC_CONF: $PYTORCH_CUDA_ALLOC_CONF"
echo ""

# Check dependencies
if ! python3 -c "import mattersim" 2>/dev/null; then
    echo "Error: MatterSim not found"
    exit 1
fi

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

# Build command
CMD="python3 prescreen.py"
CMD="$CMD --results-dir $RESULTS_DIR"
CMD="$CMD --output-dir $OUTPUT_DIR"
CMD="$CMD --max-structures $MAX_STRUCTURES"
CMD="$CMD --hull-threshold $HULL_THRESHOLD"
CMD="$CMD --device $DEVICE"
CMD="$CMD --batch-size $BATCH_SIZE"

if [ -n "$MAX_COMPOSITIONS" ]; then
    CMD="$CMD --max-compositions $MAX_COMPOSITIONS"
fi

if [ -n "$MP_API_KEY" ]; then
    CMD="$CMD --mp-api-key $MP_API_KEY"
fi

if [ -n "$PURE_PBE" ]; then
    CMD="$CMD $PURE_PBE"
fi

# Print configuration
echo "Configuration:"
echo "  Results dir: $RESULTS_DIR"
echo "  Output dir: $OUTPUT_DIR"
echo "  Max structures: $MAX_STRUCTURES"
echo "  Max compositions: ${MAX_COMPOSITIONS:-all}"
echo "  Hull threshold: ${HULL_THRESHOLD} eV/atom"
echo "  Device: $DEVICE"
echo "  Batch size: $BATCH_SIZE"
if [ -n "$PURE_PBE" ]; then
    echo "  Functional filtering: Pure GGA-PBE only"
else
    echo "  Functional filtering: Mixed PBE/PBE+U"
fi

# Check MP API key
if [ -z "$MP_API_KEY" ]; then
    echo "ERROR: MP_API_KEY environment variable not set"
    echo "Set it in ~/.bashrc: export MP_API_KEY=your_32_character_key"
    exit 1
fi
echo "  MP API key: ${MP_API_KEY:0:8}..." # Show first 8 chars only
echo ""

echo "========================================================================"
echo "Starting pre-screening..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run pre-screening
$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Pre-screening finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "Next step: Run workflow manager"
    echo "  sbatch submit_workflow_manager.sh"
fi

exit $EXIT_CODE

