#!/bin/bash

# Wrapper script to submit compare_energy_methods.py job to SLURM
# This script compares energies from different calculation methods to diagnose non-linear hull behavior

# Parse arguments
VASP_JOBS="./VASP_JOBS"
OUTPUT_DIR=""
OUTLIER_THRESHOLD=""
CONDA_ENV="vaspflow"

print_usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --vasp-jobs DIR         VASP jobs directory (default: ./VASP_JOBS)"
    echo "  --output-dir DIR        Output directory for plots (default: same as vasp-jobs)"
    echo "  --outlier-threshold VAL E_hull outlier threshold in eV/atom (default: 0.5)"
    echo "  --conda-env ENV         Conda environment name (default: vaspflow)"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --vasp-jobs ./VASP_JOBS --output-dir ./analysis_plots"
    echo "  $0 --outlier-threshold 0.5"
    echo ""
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --vasp-jobs)
            VASP_JOBS="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --outlier-threshold)
            OUTLIER_THRESHOLD="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -h|--help)
            print_usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Check if VASP_JOBS directory exists
if [ ! -d "$VASP_JOBS" ]; then
    echo "ERROR: VASP jobs directory not found: $VASP_JOBS"
    exit 1
fi

# Check required files
REQUIRED_FILES=(
    "$VASP_JOBS/mp_mattersim.json"
    "$VASP_JOBS/prescreening_stability.json"
    "$VASP_JOBS/workflow.json"
)

MISSING_FILES=0
for FILE in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$FILE" ]; then
        echo "WARNING: Required file not found: $FILE"
        MISSING_FILES=$((MISSING_FILES + 1))
    fi
done

if [ $MISSING_FILES -gt 0 ]; then
    echo "WARNING: $MISSING_FILES required file(s) missing"
    echo "The job will still be submitted but may produce partial results"
    echo ""
fi

# Build sbatch export variables
EXPORT_VARS="VASP_JOBS=$VASP_JOBS,CONDA_ENV=$CONDA_ENV"
if [ -n "$OUTPUT_DIR" ]; then
    EXPORT_VARS="$EXPORT_VARS,OUTPUT_DIR=$OUTPUT_DIR"
fi
if [ -n "$OUTLIER_THRESHOLD" ]; then
    EXPORT_VARS="$EXPORT_VARS,OUTLIER_THRESHOLD=$OUTLIER_THRESHOLD"
fi

echo "======================================================================"
echo "Submitting Energy Method Comparison Job"
echo "======================================================================"
echo "VASP jobs directory: $VASP_JOBS"
echo "Output directory: ${OUTPUT_DIR:-$VASP_JOBS (default)}"
if [ -n "$OUTLIER_THRESHOLD" ]; then
    echo "Outlier threshold: $OUTLIER_THRESHOLD eV/atom"
else
    echo "Outlier threshold: 0.2 eV/atom (default)"
fi
echo "Conda environment: $CONDA_ENV"
echo "======================================================================"
echo ""

# Submit job
JOB_ID=$(sbatch --export=ALL,$EXPORT_VARS submit_compare_energy_method.sh | awk '{print $NF}')

if [ -z "$JOB_ID" ]; then
    echo "ERROR: Job submission failed"
    exit 1
fi

echo "Job submitted successfully!"
echo "  Job ID: $JOB_ID"
echo ""
echo "Monitor job status:"
echo "  squeue -j $JOB_ID"
echo ""
echo "View output:"
echo "  tail -f compare_energy_method_${JOB_ID}.out"
echo ""
echo "Cancel job:"
echo "  scancel $JOB_ID"
echo "======================================================================"

