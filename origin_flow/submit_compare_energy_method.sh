#!/bin/bash
#SBATCH --job-name=compare_energy
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-02:00:00
#SBATCH --output=compare_energy_method_%j.out
#SBATCH --error=compare_energy_method_%j.err

# Energy Method Comparison Script
# Compares MatterSim vs DFT energies to diagnose non-linear hull behavior

echo "======================================================================"
echo "Energy Method Comparison"
echo "======================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo "======================================================================"

# Activate conda environment
source ~/.bashrc
CONDA_ENV="${CONDA_ENV:-vaspflow}"
conda activate "$CONDA_ENV"

if [ $? -ne 0 ]; then
    echo "ERROR: Could not activate $CONDA_ENV environment"
    exit 1
fi

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Python: $(which python3)"

# Set number of threads for CPU parallelization (for numpy/matplotlib)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo ""

# Get parameters from environment or use defaults
VASP_JOBS="${VASP_JOBS:-./VASP_JOBS}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
OUTLIER_THRESHOLD="${OUTLIER_THRESHOLD:-}"

# Check if VASP_JOBS directory exists
if [ ! -d "$VASP_JOBS" ]; then
    echo "ERROR: VASP jobs directory not found: $VASP_JOBS"
    exit 1
fi

echo "======================================================================"
echo "Job Parameters"
echo "======================================================================"
echo "VASP jobs directory: $VASP_JOBS"
if [ -n "$OUTPUT_DIR" ]; then
    echo "Output directory: $OUTPUT_DIR"
else
    echo "Output directory: $VASP_JOBS (default)"
fi
if [ -n "$OUTLIER_THRESHOLD" ]; then
    echo "Outlier threshold: $OUTLIER_THRESHOLD eV/atom"
else
    echo "Outlier threshold: 0.5 eV/atom (default)"
fi
echo ""

# Check input files
echo "Checking input files..."
FILES_CHECKED=0
FILES_FOUND=0

check_file() {
    local file=$1
    local desc=$2
    FILES_CHECKED=$((FILES_CHECKED + 1))
    if [ -f "$file" ]; then
        echo "    $desc: $file"
        FILES_FOUND=$((FILES_FOUND + 1))
    else
        echo "    $desc: $file (missing)"
    fi
}

check_file "$VASP_JOBS/mp_mattersim.json" "MatterSim MP cache"
check_file "$VASP_JOBS/mp_vaspdft.json" "DFT MP cache"
check_file "$VASP_JOBS/prescreening_stability.json" "Prescreening results"
check_file "$VASP_JOBS/workflow.json" "Workflow database"

echo ""
if [ $FILES_FOUND -eq 0 ]; then
    echo "ERROR: No input files found - cannot proceed"
    exit 1
elif [ $FILES_FOUND -lt $FILES_CHECKED ]; then
    echo "WARNING: $((FILES_CHECKED - FILES_FOUND)) file(s) missing - results may be partial"
fi
echo "======================================================================"
echo ""

# Build command
CMD="python3 compare_energy_methods.py --vasp-jobs \"$VASP_JOBS\""
if [ -n "$OUTPUT_DIR" ]; then
    CMD="$CMD --output-dir \"$OUTPUT_DIR\""
fi
if [ -n "$OUTLIER_THRESHOLD" ]; then
    CMD="$CMD --outlier-threshold $OUTLIER_THRESHOLD"
fi

# Run the comparison
echo "Running energy method comparison..."
echo "Command: $CMD"
echo ""

eval $CMD

EXIT_CODE=$?

echo ""
echo "======================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "Job completed successfully"
    echo ""
    echo "Generated files:"
    if [ -n "$OUTPUT_DIR" ]; then
        OUTPUT_LOC="$OUTPUT_DIR"
    else
        OUTPUT_LOC="$VASP_JOBS"
    fi
    
    echo "  Plots:"
    [ -f "$OUTPUT_LOC/mp_phases_comparison_scatter.png" ] && echo "    - mp_phases_comparison_scatter.png"
    [ -f "$OUTPUT_LOC/mp_phases_comparison_residuals.png" ] && echo "    - mp_phases_comparison_residuals.png"
    [ -f "$OUTPUT_LOC/generated_structures_comparison_scatter.png" ] && echo "    - generated_structures_comparison_scatter.png"
    [ -f "$OUTPUT_LOC/generated_structures_comparison_residuals.png" ] && echo "    - generated_structures_comparison_residuals.png"
    
    echo "  Statistics:"
    [ -f "$OUTPUT_LOC/energy_method_comparison.json" ] && echo "    - energy_method_comparison.json"
    
    echo ""
    echo "To view plots:"
    echo "  cd $OUTPUT_LOC"
    echo "  open *.png"
else
    echo "Job failed with exit code: $EXIT_CODE"
fi
echo "End time: $(date)"
echo "======================================================================"

exit $EXIT_CODE

