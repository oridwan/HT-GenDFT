#!/bin/bash
#SBATCH --job-name=mattersim_e_hull
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=80G
#SBATCH --time=96:00:00
#SBATCH --output=mattersim_e_hull_%j.out
#SBATCH --error=mattersim_e_hull_%j.err

# Usage:
#   export MP_API_KEY=your_32_character_key
#   sbatch refined_flow/submit_mattersim_e_hull.sh
#
# Optional:
#   REFINE_JOBS_DIR=/scratch/oridwan/SuperConductorFlow/REFINE_VASP-out-Boron_2 sbatch refined_flow/submit_mattersim_e_hull.sh
#   PURE_PBE= sbatch refined_flow/submit_mattersim_e_hull.sh
#
# Help:
#   bash refined_flow/submit_mattersim_e_hull.sh --help

set -e

if [ "${1:-}" = "--help" ]; then
    cat <<'EOF'
Usage:
  export MP_API_KEY=your_32_character_key
  sbatch refined_flow/submit_mattersim_e_hull.sh

Optional:
  REFINE_JOBS_DIR=/scratch/oridwan/SuperConductorFlow/REFINE_VASP-out-Boron_2 sbatch refined_flow/submit_mattersim_e_hull.sh
  DEVICE=cpu sbatch refined_flow/submit_mattersim_e_hull.sh
  PURE_PBE= sbatch refined_flow/submit_mattersim_e_hull.sh
EOF
    exit 0
fi

SCRIPT_DIR="${FLOW_SCRIPT_DIR:-}"
if [ -z "$SCRIPT_DIR" ] || [ ! -f "$SCRIPT_DIR/compute_mattersim_e_hull.py" ]; then
    if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/refined_flow/compute_mattersim_e_hull.py" ]; then
        SCRIPT_DIR="$SLURM_SUBMIT_DIR/refined_flow"
    else
        SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    fi
fi

REFINE_JOBS_DIR="${REFINE_JOBS_DIR:-./REFINE_VASP-out-Boron_2}"
DEVICE="${DEVICE:-cuda}"
CONDA_ENV="${CONDA_ENV:-mattersim}"
PURE_PBE="${PURE_PBE:-1}"

source ~/.bashrc
conda activate "$CONDA_ENV"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"
export PYTHONUNBUFFERED=1
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/$USER/matplotlib_$SLURM_JOB_ID}"
mkdir -p "$MPLCONFIGDIR"

ARGS=(
    --refine-jobs "$REFINE_JOBS_DIR"
    --mp-api-key "$MP_API_KEY"
    --device "$DEVICE"
)

if [ -n "$PURE_PBE" ]; then
    ARGS+=(--pure-pbe)
fi

python3 "$SCRIPT_DIR/compute_mattersim_e_hull.py" "${ARGS[@]}"
