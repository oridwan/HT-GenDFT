#!/bin/bash
#
# Wrapper script to submit MatterSim energy above hull calculation for refined structures
#
# Usage:
#   bash run_mattersim_e_hull.sh [OPTIONS]
#
# Options:
#   --refine-jobs DIR       Refined VASP jobs directory (default: ./REFINE_VASP_JOBS)
#   --device DEVICE         Device for MatterSim: cpu or cuda (default: cuda)
#   --conda-env ENV         Conda environment name (default: mattersim)
#   --pure-pbe              Filter MP to pure GGA-PBE only (exclude PBE+U)
#   --help                  Show this help message
#
# Environment variables:
#   MP_API_KEY              Materials Project API key (required)
#
# Examples:
#   # Basic usage (GPU by default, mixed PBE/PBE+U)
#   export MP_API_KEY=your_32_character_key
#   bash run_mattersim_e_hull.sh
#
#   # With pure GGA-PBE filtering (to match pure PBE DFT)
#   bash run_mattersim_e_hull.sh --pure-pbe
#
#   # Custom path with CPU mode (not recommended, slower)
#   bash run_mattersim_e_hull.sh --refine-jobs /scratch/$USER/REFINE_VASP_JOBS --device cpu
#

set -e

# Resolve wrapper script directory so sbatch works from any CWD
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default configuration
REFINE_JOBS_DIR="./REFINE_VASP-out-Boron_2"
DEVICE="cuda"
CONDA_ENV="mattersim"
PURE_PBE="--pure-pbe"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --refine-jobs)
            REFINE_JOBS_DIR="$2"
            shift 2
            ;;
        --device)
            DEVICE="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --pure-pbe)
            PURE_PBE="--pure-pbe"
            shift
            ;;
        --help)
            head -n 30 "$0" | tail -n 24
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Verify MP API key
if [ -z "$MP_API_KEY" ]; then
    echo "ERROR: MP_API_KEY environment variable not set"
    echo "Set it with: export MP_API_KEY=your_32_character_key"
    echo "Get key from: https://next-gen.materialsproject.org/api"
    exit 1
fi

# Validate device
if [ "$DEVICE" != "cpu" ] && [ "$DEVICE" != "cuda" ]; then
    echo "ERROR: Device must be 'cpu' or 'cuda', got: $DEVICE"
    exit 1
fi

# Convert to absolute path
REFINE_JOBS_DIR=$(cd "$(dirname "$REFINE_JOBS_DIR")" 2>/dev/null && pwd)/$(basename "$REFINE_JOBS_DIR") || REFINE_JOBS_DIR=$(realpath "$REFINE_JOBS_DIR" 2>/dev/null) || REFINE_JOBS_DIR="$REFINE_JOBS_DIR"

echo "========================================================================"
echo "MatterSim Energy Above Hull - Refined Structures"
echo "========================================================================"
echo "Configuration:"
echo "  Refined VASP jobs: $REFINE_JOBS_DIR"
echo "  Device: $DEVICE"
echo "  Conda environment: $CONDA_ENV"
if [ -n "$PURE_PBE" ]; then
    echo "  Functional filter: Pure GGA-PBE only (PBE+U excluded)"
else
    echo "  Functional filter: Mixed PBE/PBE+U (recommended)"
fi
echo "  MP API key: ${MP_API_KEY:0:8}... (${#MP_API_KEY} chars)"
echo "========================================================================"
echo ""

# Validate required files
if [ ! -d "$REFINE_JOBS_DIR" ]; then
    echo "ERROR: Refined VASP jobs directory not found: $REFINE_JOBS_DIR"
    exit 1
fi

if [ ! -f "$REFINE_JOBS_DIR/workflow.json" ]; then
    echo "ERROR: Workflow database not found: $REFINE_JOBS_DIR/workflow.json"
    echo "Run refined relaxation workflow first"
    exit 1
fi

# Count RELAX_DONE structures
echo "Checking for completed structures..."
RELAX_DONE_COUNT=$(grep -o '"state": "RELAX_DONE"' "$REFINE_JOBS_DIR/workflow.json" | wc -l)

if [ "$RELAX_DONE_COUNT" -eq 0 ]; then
    echo "ERROR: No RELAX_DONE structures found in workflow.json"
    echo "Complete refined relaxation workflow first"
    exit 1
fi

echo "  Found $RELAX_DONE_COUNT RELAX_DONE structures"
echo ""

# Export configuration for SLURM script
export REFINE_JOBS_DIR
export DEVICE
export CONDA_ENV
export PURE_PBE
export MP_API_KEY
export FLOW_SCRIPT_DIR="$SCRIPT_DIR"

# Submit job
echo "Submitting SLURM job..."
JOBID=$(sbatch --export=ALL,FLOW_SCRIPT_DIR="$SCRIPT_DIR" "$SCRIPT_DIR/submit_mattersim_e_hull.sh" | awk '{print $NF}')

if [ -z "$JOBID" ]; then
    echo "ERROR: Failed to submit job"
    exit 1
fi

echo ""
echo "========================================================================"
echo "Job Submitted Successfully"
echo "========================================================================"
echo "  Job ID: $JOBID"
echo "  Structures to process: $RELAX_DONE_COUNT"
echo ""
echo "Monitoring commands:"
echo "  View log:        tail -f mattersim_e_hull_${JOBID}.out"
echo "  Check status:    squeue -j $JOBID"
echo "  Cancel job:      scancel $JOBID"
echo ""
echo "Output files (after completion):"
echo "  $REFINE_JOBS_DIR/mattersim_stability_results.json (MatterSim energies & hulls)"
echo "  $REFINE_JOBS_DIR/mp_mattersim.json (MP phases relaxed with MatterSim)"
echo "========================================================================"
