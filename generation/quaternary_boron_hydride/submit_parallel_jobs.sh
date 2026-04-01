#!/bin/bash
# Submit parallel generation jobs with optional START_INDEX and MAX_COMPOSITIONS
# Usage:
#   bash submit_parallel_jobs.sh                       # Process all compositions
#   bash submit_parallel_jobs.sh 0 1300                # Process compositions [0, 2700)
#   bash submit_parallel_jobs.sh 2700 4200             # Process compositions [2700, 4200)
#   bash submit_parallel_jobs.sh 2700 4130             # Process compositions [2700, end)

# Configuration
COMPOSITIONS_PER_JOB=50
COMPOSITIONS_FILE="quaternary_boron_hydride_compositions.json"

# Command-line arguments
START_INDEX=${1:-0}
END_INDEX=${2:-0}  # 0 means process all remaining

# Auto-detect total compositions from JSON file
if [ ! -f "$COMPOSITIONS_FILE" ]; then
    echo "ERROR: Compositions file not found: $COMPOSITIONS_FILE"
    echo "Please run search_quaternary_electrides.py first"
    exit 1
fi

# Count compositions in JSON file
if command -v python3 &> /dev/null; then
    TOTAL_COMPOSITIONS=$(python3 -c "import json; print(len(json.load(open('$COMPOSITIONS_FILE'))))")
elif command -v jq &> /dev/null; then
    TOTAL_COMPOSITIONS=$(jq 'length' "$COMPOSITIONS_FILE")
else
    echo "ERROR: Neither python3 nor jq found. Please install one of them."
    exit 1
fi

echo "Found $TOTAL_COMPOSITIONS total compositions in $COMPOSITIONS_FILE"
echo "Starting from index: $START_INDEX"

# Normalize END_INDEX
if [ "$END_INDEX" -eq 0 ] || [ "$END_INDEX" -gt "$TOTAL_COMPOSITIONS" ]; then
    END_INDEX=$TOTAL_COMPOSITIONS
fi

if [ "$END_INDEX" -le "$START_INDEX" ]; then
    echo "ERROR: END_INDEX must be greater than START_INDEX"
    echo "Got START_INDEX=$START_INDEX, END_INDEX=$END_INDEX"
    exit 1
fi

# Calculate compositions to process
COMPOSITIONS_TO_PROCESS=$((END_INDEX - START_INDEX))

echo "Processing: $COMPOSITIONS_TO_PROCESS compositions ([${START_INDEX}, ${END_INDEX}))"

# Calculate number of jobs needed
NUM_JOBS=$(( (COMPOSITIONS_TO_PROCESS + COMPOSITIONS_PER_JOB - 1) / COMPOSITIONS_PER_JOB ))

echo "Submitting $NUM_JOBS parallel generation jobs"
echo "Each job will process $COMPOSITIONS_PER_JOB compositions"
echo "========================================="

ARRAY_SPEC="0-$((NUM_JOBS - 1))"
JOB_ID=$(sbatch --array=$ARRAY_SPEC \
    --export=ALL,COMPOSITIONS_PER_JOB=$COMPOSITIONS_PER_JOB,BASE_START_INDEX=$START_INDEX,GLOBAL_END_INDEX=$END_INDEX \
    generate_quaternary.sh | awk '{print $NF}')

echo "Submitted array job $JOB_ID with array spec: $ARRAY_SPEC"

echo "=========================================="
echo "All jobs submitted!"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check progress: ls -d ../results/quaternary_boron_hydride/*_structures | wc -l"
echo "=========================================="
