#!/bin/bash
# Submit multiple parallel generation jobs for quaternary hydrides
# Each job processes a different subset of compositions

# Configuration
COMPOSITIONS_PER_JOB=400
COMPOSITIONS_FILE="quaternary_hydride_compositions.json"
MAX_CONCURRENT_JOBS=12

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

echo "Found $TOTAL_COMPOSITIONS compositions in $COMPOSITIONS_FILE"

# Calculate number of jobs needed
NUM_JOBS=$(( (TOTAL_COMPOSITIONS + COMPOSITIONS_PER_JOB - 1) / COMPOSITIONS_PER_JOB ))

echo "Submitting $NUM_JOBS parallel generation jobs"
echo "Each job will process $COMPOSITIONS_PER_JOB compositions"
echo "Max concurrent jobs: $MAX_CONCURRENT_JOBS"
echo "=========================================="

ARRAY_SPEC="0-$((NUM_JOBS - 1))%${MAX_CONCURRENT_JOBS}"
JOB_ID=$(sbatch --array=$ARRAY_SPEC \
    --export=ALL,COMPOSITIONS_PER_JOB=$COMPOSITIONS_PER_JOB \
    generate_quaternary.sh | awk '{print $NF}')

echo "Submitted array job $JOB_ID with array spec: $ARRAY_SPEC"

echo "=========================================="
echo "All jobs submitted!"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check progress: ls -d ../results/quaternary_hydride/*_structures | wc -l"
echo "=========================================="

