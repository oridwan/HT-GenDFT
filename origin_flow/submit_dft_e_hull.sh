#!/bin/bash
#SBATCH --job-name=dft_e_hull
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=dft_e_hull_%j.out
#SBATCH --error=dft_e_hull_%j.err

# DFT Energy Above Hull Calculation
# Computes DFT-level thermodynamic stability using VASP energies + MP competing phases

echo "======================================================================"
echo "DFT Energy Above Hull Calculation"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo "======================================================================"

# Activate conda environment
source ~/.bashrc
conda activate mattersim

if [ $? -ne 0 ]; then
    echo "ERROR: Could not activate mattersim environment"
    exit 1
fi

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python: $(which python3)"

# Set number of threads for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo ""

# Check MP API key
if [ -z "$MP_API_KEY" ]; then
    echo "ERROR: MP_API_KEY environment variable not set"
    echo "Set it in ~/.bashrc: export MP_API_KEY=your_32_character_key"
    exit 1
fi

echo "MP_API_KEY: ${MP_API_KEY:0:8}..." # Show first 8 chars only
echo "======================================================================"
echo ""

# Default parameters (can be overridden by command line)
VASP_JOBS="${VASP_JOBS:-./VASP_JOBS}"
OUTPUT="${OUTPUT:-dft_stability_results.json}"
PRESCREEN_RESULTS="${PRESCREEN_RESULTS:-./VASP_JOBS/prescreening_stability.json}"

# Run the computation
echo "Running DFT hull calculation..."
echo "  VASP jobs: $VASP_JOBS"
echo "  Output: $OUTPUT"
echo "  Pre-screening filter: $PRESCREEN_RESULTS"
if [ -n "$PURE_PBE" ]; then
    echo "  Functional filtering: Pure GGA-PBE only"
else
    echo "  Functional filtering: Mixed PBE/PBE+U"
fi
if [ -n "$HULL_THRESHOLD" ]; then
    echo "  Hull threshold: $HULL_THRESHOLD eV/atom"
else
    echo "  Hull threshold: 0.1 eV/atom (default)"
fi
if [ -n "$OUTLIER_THRESHOLD" ]; then
    echo "  Outlier threshold: $OUTLIER_THRESHOLD eV/atom"
else
    echo "  Outlier threshold: 0.5 eV/atom (default)"
fi
echo ""

# Build command with optional arguments
CMD="python3 compute_dft_e_hull.py --vasp-jobs \"$VASP_JOBS\" --output \"$OUTPUT\" --prescreen-results \"$PRESCREEN_RESULTS\""

if [ -n "$PURE_PBE" ]; then
    CMD="$CMD $PURE_PBE"
fi

if [ -n "$HULL_THRESHOLD" ]; then
    CMD="$CMD --hull-threshold $HULL_THRESHOLD"
fi

if [ -n "$OUTLIER_THRESHOLD" ]; then
    CMD="$CMD --outlier-threshold $OUTLIER_THRESHOLD"
fi

eval $CMD

EXIT_CODE=$?

echo ""
echo "======================================================================"
echo "Job completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE

