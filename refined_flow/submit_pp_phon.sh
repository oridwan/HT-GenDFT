#!/bin/bash
#SBATCH --job-name=phonon_pp
#SBATCH --partition=Apus,Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=phonon_postproc_%j.out
#SBATCH --error=phonon_postproc_%j.err

# Phonon Post-Processing SLURM Script

set -e

# Default parameters (can be overridden by environment variables)
PHONON_JOBS=${PHONON_JOBS:-"./PHONON_JOBS"}
STRUCTURE_IDS=${STRUCTURE_IDS:-""}
OUTPUT_SUMMARY=${OUTPUT_SUMMARY:-"phonon_postproc_summary.json"}
CONDA_ENV=${CONDA_ENV:-"vaspflow"}

echo "========================================================================"
echo "Phonon Post-Processing (SLURM Job)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo ""

# Load environment
source ~/.bashrc
conda activate $CONDA_ENV

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment '$CONDA_ENV'"
    exit 1
fi

echo "Conda environment: $CONDA_ENV"
echo "Python: $(which python3)"
echo ""

# Set thread limits for numpy/scipy parallel operations
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Parallel threads set to: $SLURM_CPUS_PER_TASK"
echo ""

# Check dependencies
echo "Checking dependencies..."
if ! python3 -c "import pymatgen" 2>/dev/null; then
    echo "Error: pymatgen not found"
    exit 1
fi

if ! python3 -c "import phonopy" 2>/dev/null; then
    echo "Error: phonopy not found"
    exit 1
fi

echo "pymatgen version: $(python3 -c 'from importlib.metadata import version; print(version("pymatgen"))' 2>/dev/null || echo "installed")"
echo "phonopy version: $(python3 -c 'import phonopy; print(phonopy.__version__)')"
echo "numpy version: $(python3 -c 'import numpy; print(numpy.__version__)')"
echo ""

# Expand paths
PHONON_JOBS=$(eval echo "$PHONON_JOBS")

# Check phonon jobs directory
if [ ! -d "$PHONON_JOBS" ]; then
    echo "Error: PHONON_JOBS directory not found: $PHONON_JOBS"
    exit 1
fi

# Check workflow database
WORKFLOW_DB="$PHONON_JOBS/workflow.json"
if [ ! -f "$WORKFLOW_DB" ]; then
    echo "Error: Workflow database not found: $WORKFLOW_DB"
    echo "Please complete the phonon workflow first"
    exit 1
fi

# Count PHON_DONE structures
if command -v jq &> /dev/null; then
    PHON_DONE_COUNT=$(jq '[.structures | to_entries[] | select(.value.state == "PHON_DONE")] | length' "$WORKFLOW_DB")
    echo "PHON_DONE structures found: $PHON_DONE_COUNT"
    
    if [ "$PHON_DONE_COUNT" -eq 0 ]; then
        echo "Warning: No PHON_DONE structures found in database"
        echo "Please complete phonon calculations first (run_phonon_flow.sh)"
        exit 0
    fi
    echo ""
fi

# Build command
CMD="python3 postproc_phonon.py"
CMD="$CMD --phonon-jobs $PHONON_JOBS"
CMD="$CMD --output-summary $OUTPUT_SUMMARY"

# Add structure IDs if specified
if [ -n "$STRUCTURE_IDS" ]; then
    CMD="$CMD --structure-ids$STRUCTURE_IDS"
    echo "Processing specific structures:$STRUCTURE_IDS"
else
    echo "Processing all PHON_DONE structures"
fi

echo ""

# Print configuration
echo "Configuration:"
echo "  PHONON_JOBS: $PHONON_JOBS"
echo "  Output summary: $OUTPUT_SUMMARY"
echo "  Database: $WORKFLOW_DB"
echo ""

echo "========================================================================"
echo "Starting phonon post-processing..."
echo "========================================================================"
echo ""
echo "Command: $CMD"
echo ""

# Run post-processing
$CMD

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Phonon post-processing finished"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "Output files created in each PHON directory:"
    echo "  - phonopy_params.yaml (force constants)"
    echo "  - phonon_band_dos.png (combined band structure + DOS plot)"
    echo "  - phonon_band_dos.pdf (vector format)"
    echo "  - phonon_band.dat (band structure data)"
    echo "  - band_kpath.dat (k-path metadata: segments + labels)"
    echo "  - phonon_dos.dat (DOS data: total + element-projected)"
    echo "  - thermal.dat (thermal properties)"
    echo ""
    echo "Summary saved to: $PHONON_JOBS/$OUTPUT_SUMMARY"
fi

exit $EXIT_CODE
