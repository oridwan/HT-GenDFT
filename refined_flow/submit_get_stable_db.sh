#!/bin/bash
#SBATCH --job-name=get_stable_db
#SBATCH --partition=Orion,Apus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=get_stable_db_%j.out
#SBATCH --error=get_stable_db_%j.err

# Extract stable electride candidates confirmed by both MatterSim and VASP-DFT
# Reads pre-computed results from JSON files in REFINE_VASP_JOBS directory
# No MatterSim or heavy computation required - just reads JSONs and CONTCARs

echo "========================================================================"
echo "Get Stable Electride Database"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "========================================================================"
echo ""

# Load conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate ${CONDA_ENV:-vaspflow}

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python: $(which python3)"
echo ""

# Set threading for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# Configuration from environment variables
echo "Configuration:"
echo "  Script directory: $SCRIPT_DIR"
echo "  Working directory: $CURRENT_DIR"
echo "  Refined VASP jobs: $REFINE_JOBS"
echo "  Electride CSV: $ELECTRIDE_CSV"
echo "  Threshold: $THRESHOLD eV/atom"
echo "========================================================================"
echo ""

# Change to working directory
cd "$CURRENT_DIR"

# Run the script
echo "Running get_stable_ele_db.py..."
echo ""

python3 "$SCRIPT_DIR/get_stable_ele_db.py" \
    --refine-jobs "$REFINE_JOBS" \
    --electride-csv "$ELECTRIDE_CSV" \
    --threshold "$THRESHOLD"

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Job Completed"
echo "========================================================================"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Status: SUCCESS"
    echo ""
    echo "Output files:"
    echo "  $REFINE_JOBS/stable_electrides.csv"
    echo "  $REFINE_JOBS/stable_electrides.db"
else
    echo "Status: FAILED"
fi

echo "========================================================================"

exit $EXIT_CODE
