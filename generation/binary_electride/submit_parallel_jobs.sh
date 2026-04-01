#!/bin/bash
# Submit multiple parallel generation jobs for binary electrides
# Each job processes a different subset of compositions

# Configuration
COMPOSITIONS_PER_JOB=400
COMPOSITIONS_FILE="binary_electride_compositions.json"

# Auto-detect total compositions from JSON file
if [ ! -f "$COMPOSITIONS_FILE" ]; then
    echo "ERROR: Compositions file not found: $COMPOSITIONS_FILE"
    echo "Please run search_binary_electrides.py first"
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
echo "=========================================="

for i in $(seq 0 $((NUM_JOBS - 1))); do
    START_IDX=$((i * COMPOSITIONS_PER_JOB))
    
    TEMP_SCRIPT="generate_binary_csp_batch${i}.sh"
    
    sed "s/^START_INDEX=.*/START_INDEX=$START_IDX/" generate_binary_csp.sh | \
    sed "s/^MAX_COMPOSITIONS=.*/MAX_COMPOSITIONS=$COMPOSITIONS_PER_JOB/" | \
    sed "s/gen_bin_ele_csp_%j/gen_bin_ele_csp_batch${i}_%j/" > $TEMP_SCRIPT
    
    JOB_ID=$(sbatch $TEMP_SCRIPT | awk '{print $NF}')
    
    echo "Submitted batch $i (compositions $START_IDX-$((START_IDX + COMPOSITIONS_PER_JOB - 1))): Job ID $JOB_ID"
    
    rm $TEMP_SCRIPT
    
    sleep 1
done

echo "=========================================="
echo "All jobs submitted!"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check progress: ls -d ../results/binary_csp_electrides/*_structures | wc -l"
echo "=========================================="

