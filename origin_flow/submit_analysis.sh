#!/bin/bash
# Convenience script to submit electride analysis job
# Uses HT-electride methodology

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

# Default values
VASP_JOBS_DIR="${VASP_JOBS_DIR:-}"
BADER_EXE="${BADER_EXE:-}"
OUTPUT_CSV="${OUTPUT_CSV:-}"
ELF_THRESHOLD="${ELF_THRESHOLD:-0.6}"
CONDA_ENV="${CONDA_ENV:-mattersim}"

# Help message
show_help() {
    cat << EOF
Submit electride analysis job to SLURM

Usage: bash submit_analysis.sh [OPTIONS]

Options:
  --vasp-jobs DIR        Path to VASP jobs directory (default: <repo>/VASP-out-Boron)
  --bader-exe PATH       Path to bader executable (default: ~/miniconda3/envs/mattersim/bin/bader)
  --output FILE          Output CSV file (default: <vasp-jobs>/electride_analysis.csv)
  --threshold VALUE      ELF threshold (default: 0.6)
  --conda-env NAME       Conda environment name (default: mattersim)
  -h, --help             Show this help message

Examples:
  # Submit with default settings
  bash submit_analysis.sh

  # Use custom VASP jobs directory
  bash submit_analysis.sh --vasp-jobs /path/to/VASP_JOBS

  # Adjust ELF threshold
  bash submit_analysis.sh --threshold 0.7

  # Custom output file
  bash submit_analysis.sh --output my_electrides.csv

Environment Variables:
  All options can also be set via environment variables:
    VASP_JOBS_DIR, BADER_EXE, OUTPUT_CSV, ELF_THRESHOLD, CONDA_ENV

Method:
  This script uses Electride.py in batch mode (-p option) to analyze all
  completed ELF calculations. Results are output as a CSV table showing
  interstitial electron volumes for different energy windows (e0025, e05, e1)
  and bands (band0, band1).

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vasp-jobs)
            VASP_JOBS_DIR="$2"
            shift 2
            ;;
        --bader-exe)
            BADER_EXE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_CSV="$2"
            shift 2
            ;;
        --threshold)
            ELF_THRESHOLD="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

if [ -z "$VASP_JOBS_DIR" ]; then
    VASP_JOBS_DIR="$REPO_ROOT/VASP-out-Boron"
fi
if [ -z "$BADER_EXE" ]; then
    BADER_EXE="$HOME/miniconda3/envs/mattersim/bin/bader"
fi
if [ -z "$OUTPUT_CSV" ]; then
    OUTPUT_CSV="$VASP_JOBS_DIR/electride_analysis.csv"
fi

# Export variables for SLURM script
export VASP_JOBS_DIR
export BADER_EXE
export OUTPUT_CSV
export ELF_THRESHOLD
export CONDA_ENV
export REPO_ROOT
export SCRIPT_DIR

# Check if VASP_JOBS directory exists
if [ ! -d "$VASP_JOBS_DIR" ]; then
    echo "ERROR: VASP jobs directory not found: $VASP_JOBS_DIR"
    exit 1
fi

# Check if workflow database exists
WORKFLOW_DB="$VASP_JOBS_DIR/workflow.json"
if [ ! -f "$WORKFLOW_DB" ]; then
    echo "ERROR: Workflow database not found: $WORKFLOW_DB"
    exit 1
fi

# Count completed ELF calculations from workflow database
ELF_COUNT=$(python3 << EOF
import json
with open("$WORKFLOW_DB", 'r') as f:
    data = json.load(f)
count = sum(1 for s in data['structures'].values() if s['state'] == 'ELF_DONE')
print(count)
EOF
)

# Print configuration
echo "========================================================================"
echo "Submitting Electride Analysis Job"
echo "========================================================================"
echo "Configuration:"
echo "  VASP jobs directory: $VASP_JOBS_DIR"
echo "  Completed ELF calculations: $ELF_COUNT"
echo "  Bader executable: $BADER_EXE"
echo "  Output CSV: $OUTPUT_CSV"
echo "  ELF threshold: $ELF_THRESHOLD"
echo "  Conda environment: $CONDA_ENV"
echo "========================================================================"
echo ""

if [ $ELF_COUNT -eq 0 ]; then
    echo "WARNING: No completed ELF calculations found in $VASP_JOBS_DIR"
    echo "         Make sure the workflow has finished ELF stage"
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Submit job
JOBID=$(sbatch "$SCRIPT_DIR/analyze.sh" | awk '{print $NF}')

if [ $? -eq 0 ]; then
    echo "Job submitted successfully!"
    echo "Job ID: $JOBID"
    echo ""
    echo "Monitor the job with:"
    echo "  squeue -j $JOBID"
    echo ""
    echo "Check output:"
    echo "  tail -f electride_analysis_${JOBID}.out"
    echo ""
    echo "View results when complete:"
    echo "  cat $OUTPUT_CSV"
    echo "  cat electride_analysis_detailed.log  # detailed per-structure output"
else
    echo "ERROR: Failed to submit job"
    exit 1
fi
