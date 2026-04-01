#!/bin/bash
#SBATCH --job-name=analyze_electrides
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --output=electride_analysis_%j.out
#SBATCH --error=electride_analysis_%j.err

# SLURM script to analyze completed ELF calculations for electride candidates

# Resolve paths robustly under sbatch, which runs a spooled script copy.
REPO_ROOT=${REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}
SCRIPT_DIR=${SCRIPT_DIR:-}
if [ -n "$REPO_ROOT" ] && [ -z "$SCRIPT_DIR" ]; then
    SCRIPT_DIR="$REPO_ROOT/origin_flow"
fi
if [ -z "$REPO_ROOT" ] || [ ! -f "$SCRIPT_DIR/analyze.py" ]; then
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
    REPO_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
fi

# Configuration (can be overridden by environment variables)
VASP_JOBS_DIR=${VASP_JOBS_DIR:-"$REPO_ROOT/VASP-out-Boron"}
BADER_EXE=${BADER_EXE:-"$HOME/miniconda3/envs/mattersim/bin/bader"}
OUTPUT_CSV=${OUTPUT_CSV:-"$VASP_JOBS_DIR/electride_analysis.csv"}
PYXTAL_DB=${PYXTAL_DB:-"$VASP_JOBS_DIR/electride_data.db"}
PRESCREENING=${PRESCREENING:-"$REPO_ROOT/prescreen_new/VASP_JOBS_Boron/prescreening_stability.json"}
ELF_THRESHOLD=${ELF_THRESHOLD:-0.6}
CONDA_ENV=${CONDA_ENV:-"mattersim"}

echo "========================================================================"
echo "Electride Analysis Job"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""
echo "Configuration:"
echo "  VASP jobs directory: $VASP_JOBS_DIR"
echo "  Bader executable: $BADER_EXE"
echo "  Output CSV: $OUTPUT_CSV"
echo "  PyXtal database: $PYXTAL_DB"
echo "  Prescreening results: $PRESCREENING"
echo "  ELF threshold: $ELF_THRESHOLD"
echo "  Conda environment: $CONDA_ENV"
echo ""
echo "Features:"
echo "  - Incremental analysis (skips already-analyzed structures)"
echo "  - Electride criteria: (e0025 > 0 OR e05 > 0) AND (e10 > 0 OR band0 > 0)"
echo "  - Uses MatterSim e_above_hull from prescreening"
echo "  - Adds spacegroup from CONTCAR"
echo "  - Saves to PyXtal database"
echo "  - Parallel processing with ${SLURM_CPUS_PER_TASK} workers"
echo "========================================================================"
echo ""

# Activate conda environment
source ~/.bashrc
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate conda environment: $CONDA_ENV"
    exit 1
fi

echo "Conda environment activated: $(which python3)"
echo ""

# Keep relative logs in the repo instead of the caller's working directory.
cd "$REPO_ROOT" || exit 1

# Set number of threads for CPU parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}

echo "Parallelization settings:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS: $MKL_NUM_THREADS"
echo ""

# Check if VASP_JOBS directory exists
if [ ! -d "$VASP_JOBS_DIR" ]; then
    echo "ERROR: VASP jobs directory not found: $VASP_JOBS_DIR"
    exit 1
fi

# Check if bader executable exists
if [ ! -x "$BADER_EXE" ]; then
    echo "WARNING: Bader executable not found or not executable: $BADER_EXE"
    echo "         Will try 'bader' from PATH"
    BADER_CMD="bader"
else
    BADER_CMD="$BADER_EXE"
fi

# Check if workflow database exists
WORKFLOW_DB="${VASP_JOBS_DIR}/workflow.json"
if [ ! -f "$WORKFLOW_DB" ]; then
    echo "ERROR: Workflow database not found: $WORKFLOW_DB"
    exit 1
fi

# Run analysis script
echo "Running electride analysis..."
echo ""

python3 "$SCRIPT_DIR/analyze.py" \
    --db "$WORKFLOW_DB" \
    --bader-exe "$BADER_CMD" \
    --threshold "$ELF_THRESHOLD" \
    --output "$OUTPUT_CSV" \
    --pyxtal-db "$PYXTAL_DB" \
    --prescreening "$PRESCREENING" \
    --workers "${SLURM_CPUS_PER_TASK:-32}" \
    2>&1 | tee electride_analysis_detailed.log

EXIT_CODE=$?

# Report output files
if [ $EXIT_CODE -eq 0 ] && [ -f "$OUTPUT_CSV" ]; then
    echo ""
    echo "Final CSV output: ${OUTPUT_CSV}"
else
    echo "WARNING: No output CSV generated"
fi

if [ $EXIT_CODE -eq 0 ] && [ -f "$PYXTAL_DB" ]; then
    echo "PyXtal database: ${PYXTAL_DB}"
fi

echo ""
echo "========================================================================"
echo "Job completed"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "========================================================================"

exit $EXIT_CODE
