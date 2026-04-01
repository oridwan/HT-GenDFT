#!/bin/bash
#SBATCH --job-name=ms_hull_origin
#SBATCH --partition=GPU,Hydrus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=mattersim_e_hull_origin_%j.out
#SBATCH --error=mattersim_e_hull_origin_%j.err

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

VASP_JOBS_DIR="${VASP_JOBS_DIR:-$REPO_ROOT/VASP-out-Boron}"
CANDIDATE_CSV="${CANDIDATE_CSV:-$VASP_JOBS_DIR/electride_analysis_candidates-Strict.csv}"
OUTPUT_JSON="${OUTPUT_JSON:-$VASP_JOBS_DIR/mattersim_stability_candidates.json}"
CONDA_ENV="${CONDA_ENV:-mattersim}"
DEVICE="${DEVICE:-cuda}"
PURE_PBE="${PURE_PBE:-1}"

source ~/.bashrc
conda activate "$CONDA_ENV"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export PYTHONUNBUFFERED=1
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/$USER/matplotlib_$SLURM_JOB_ID}"
mkdir -p "$MPLCONFIGDIR"

ARGS=(
    --vasp-jobs "$VASP_JOBS_DIR"
    --candidate-csv "$CANDIDATE_CSV"
    --output "$OUTPUT_JSON"
    --mp-api-key "$MP_API_KEY"
    --device "$DEVICE"
)

if [ -n "$PURE_PBE" ]; then
    ARGS+=(--pure-pbe)
fi

python3 "$SCRIPT_DIR/compute_mattersim_e_hull_origin.py" "${ARGS[@]}"
