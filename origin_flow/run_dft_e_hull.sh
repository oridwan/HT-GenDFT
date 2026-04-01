#!/bin/bash
# Wrapper script to submit DFT energy_above_hull calculation as SLURM job

set -e

# Default values
VASP_JOBS="./VASP_JOBS"
OUTPUT="dft_stability_results.json"
PRESCREEN_RESULTS="./VASP_JOBS/prescreening_stability.json"
PURE_PBE=""
HULL_THRESHOLD=""
OUTLIER_THRESHOLD=""
SUBMIT_SCRIPT="submit_dft_e_hull.sh"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vasp-jobs)
            VASP_JOBS="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --prescreen-results)
            PRESCREEN_RESULTS="$2"
            shift 2
            ;;
        --pure-pbe)
            PURE_PBE="--pure-pbe"
            shift 1
            ;;
        --hull-threshold)
            HULL_THRESHOLD="$2"
            shift 2
            ;;
        --outlier-threshold)
            OUTLIER_THRESHOLD="$2"
            shift 2
            ;;
        -h|--help)
            cat << EOF
Usage: $0 [OPTIONS]

Compute DFT energy_above_hull for VASP-relaxed structures.

Options:
    --vasp-jobs PATH           VASP jobs directory (default: ./VASP_JOBS)
    --output FILE              Output JSON file (default: dft_stability_results.json)
    --prescreen-results FILE   Pre-screening filter (default: ./VASP_JOBS/prescreening_stability.json)
    --pure-pbe                 Filter MP entries to pure GGA-PBE only (exclude PBE+U)
    --hull-threshold VALUE     E_hull threshold for stability analysis in eV/atom (default: 0.1)
    --outlier-threshold VALUE  E_hull outlier threshold for plot filtering in eV/atom (default: 0.5)
    -h, --help                 Show this help message

Example:
    # Basic usage (mixed PBE/PBE+U for best accuracy)
    bash run_dft_e_hull.sh

    # Pure PBE filtering (strict functional consistency)
    bash run_dft_e_hull.sh --pure-pbe

    # Custom thresholds
    bash run_dft_e_hull.sh --hull-threshold 0.15 --outlier-threshold 1.0

    # Custom paths
    bash run_dft_e_hull.sh \\
        --vasp-jobs ./VASP_JOBS \\
        --output my_dft_results.json \\
        --prescreen-results ./VASP_JOBS/prescreening_stability.json

Notes:
    - Requires mattersim conda environment
    - Requires MP_API_KEY environment variable (32 characters)
    - Only processes structures that completed relaxation (RELAX_DONE or later)
    - If --prescreen-results provided, only processes structures that passed pre-screening
    - Uses VASP energies from vasprun.xml + MP DFT energies for competing phases
    - Default: Mixed PBE/PBE+U (MP recommended methodology for accuracy)
    - Use --pure-pbe only if your VASP uses pure PBE without +U corrections

EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if submit script exists
if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "ERROR: Submit script not found: $SUBMIT_SCRIPT"
    exit 1
fi

# Export variables for SLURM script
export VASP_JOBS
export OUTPUT
export PRESCREEN_RESULTS
export PURE_PBE
export HULL_THRESHOLD
export OUTLIER_THRESHOLD

echo "======================================================================"
echo "Submitting DFT Energy Above Hull Calculation"
echo "======================================================================"
echo "VASP jobs: $VASP_JOBS"
echo "Output: $OUTPUT"
echo "Pre-screening filter: $PRESCREEN_RESULTS"
if [ -n "$PURE_PBE" ]; then
    echo "Functional filtering: Pure GGA-PBE only (PBE+U excluded)"
else
    echo "Functional filtering: Mixed PBE/PBE+U (recommended)"
fi
if [ -n "$HULL_THRESHOLD" ]; then
    echo "Hull threshold: $HULL_THRESHOLD eV/atom"
else
    echo "Hull threshold: 0.1 eV/atom (default)"
fi
if [ -n "$OUTLIER_THRESHOLD" ]; then
    echo "Outlier threshold: $OUTLIER_THRESHOLD eV/atom"
else
    echo "Outlier threshold: 0.5 eV/atom (default)"
fi
echo "======================================================================"

# Submit job
JOB_ID=$(sbatch "$SUBMIT_SCRIPT" | awk '{print $4}')

if [ -n "$JOB_ID" ]; then
    echo "Job submitted successfully!"
    echo "Job ID: $JOB_ID"
    echo ""
    echo "Monitor with:"
    echo "  squeue -j $JOB_ID"
    echo "  tail -f dft_e_hull_${JOB_ID}.out"
    echo ""
    echo "When complete, results will be in: $VASP_JOBS/$OUTPUT"
else
    echo "ERROR: Job submission failed"
    exit 1
fi

