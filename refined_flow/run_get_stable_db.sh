#!/bin/bash
#
# Wrapper script to submit get_stable_ele_db.py via SLURM
#
# Usage:
#   ./run_get_stable_db.sh [OPTIONS]
#
# Options:
#   --refine-jobs DIR       Refined VASP jobs directory (default: ./REFINE_VASP_JOBS)
#   --electride-csv PATH    Electride candidates CSV (default: ./electride_candidates.csv)
#   --threshold VALUE       Energy above hull threshold in eV/atom (default: 0.005)
#   --conda-env ENV         Conda environment name (default: vaspflow)
#
# Example:
#   ./run_get_stable_db.sh --refine-jobs REFINE_VASP_JOBS/ --threshold 0.01
#

set -e

# Default values
REFINE_JOBS="./REFINE_VASP_JOBS"
ELECTRIDE_CSV="./electride_candidates.csv"
THRESHOLD="0.005"
CONDA_ENV="vaspflow"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --refine-jobs)
            REFINE_JOBS="$2"
            shift 2
            ;;
        --electride-csv)
            ELECTRIDE_CSV="$2"
            shift 2
            ;;
        --threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --refine-jobs DIR       Refined VASP jobs directory (default: ./REFINE_VASP_JOBS)"
            echo "  --electride-csv PATH    Electride candidates CSV (default: ./electride_candidates.csv)"
            echo "  --threshold VALUE       Energy above hull threshold in eV/atom (default: 0.005)"
            echo "  --conda-env ENV         Conda environment name (default: vaspflow)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Get script directory (where this script and submit script are located)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get current working directory (where command is run from)
CURRENT_DIR="$(pwd)"

# Convert to absolute paths
if [[ ! "$REFINE_JOBS" = /* ]]; then
    REFINE_JOBS="$CURRENT_DIR/$REFINE_JOBS"
fi

if [[ ! "$ELECTRIDE_CSV" = /* ]]; then
    ELECTRIDE_CSV="$CURRENT_DIR/$ELECTRIDE_CSV"
fi

# Check if refine jobs directory exists
if [[ ! -d "$REFINE_JOBS" ]]; then
    echo "ERROR: Refined VASP jobs directory not found: $REFINE_JOBS"
    exit 1
fi

# Export environment variables for SLURM script
export SCRIPT_DIR
export CURRENT_DIR
export REFINE_JOBS
export ELECTRIDE_CSV
export THRESHOLD
export CONDA_ENV

echo "========================================"
echo "Submitting get_stable_ele_db.py"
echo "========================================"
echo "Script directory: $SCRIPT_DIR"
echo "Working directory: $CURRENT_DIR"
echo "Refined VASP jobs: $REFINE_JOBS"
echo "Electride CSV: $ELECTRIDE_CSV"
echo "Threshold: $THRESHOLD eV/atom"
echo "Conda environment: $CONDA_ENV"
echo "========================================"
echo ""

# Submit SLURM job
sbatch "$SCRIPT_DIR/submit_get_stable_db.sh"

