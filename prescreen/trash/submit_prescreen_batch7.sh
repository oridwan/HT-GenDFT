#!/bin/bash
#SBATCH --job-name=ps_b7
#SBATCH --partition=Hydrus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=prescreen_batch7_%j.out
#SBATCH --error=prescreen_batch7_%j.err

set -euo pipefail

echo "Job start: $(date)"
echo "Host: $(hostname)"
echo "PWD: $(pwd)"

# Robust conda activation: do not let a non-interactive .bashrc abort the job silently.
set +e
set +u
source ~/.bashrc
BASHRC_RC=$?
set -u
set -e
if [ "$BASHRC_RC" -ne 0 ]; then
    echo "WARNING: source ~/.bashrc returned non-zero ($BASHRC_RC); continuing with explicit conda init"
fi

if ! command -v conda >/dev/null 2>&1; then
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    fi
fi

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda command not found after initialization"
    exit 1
fi

conda activate mattersim || {
    echo "ERROR: Failed to activate conda env: mattersim"
    exit 1
}

if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader || true
fi

export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export MP_API_KEY="Yw15APriCn1jSHZE9pzSVQUMtrMwJSew"

START_INDEX=3500
MAX_COMPS=500
BATCH_ID=7
RESULTS_DIR="../results/quaternary_hydride"
OUTPUT_DIR="./VASP_JOBS"
BATCH_SIZE=32
HULL_THRESHOLD=0.1
MAX_STRUCTURES=0
DEVICE="cuda"
PURE_PBE=""

echo "========================================"
echo "Prescreen batch $BATCH_ID on partition Hydrus"
echo "Composition range: $START_INDEX - $((START_INDEX + MAX_COMPS - 1))"
echo "========================================"

python3 -u prescreen.py     --results-dir "$RESULTS_DIR"     --output-dir "$OUTPUT_DIR"     --hull-threshold "$HULL_THRESHOLD"     --device "$DEVICE"     --batch-size "$BATCH_SIZE"     --start-composition "$START_INDEX"     --max-compositions "$MAX_COMPS"     --max-structures "$MAX_STRUCTURES"     --batch-id "$BATCH_ID"     $PURE_PBE
