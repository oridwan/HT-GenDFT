#!/bin/bash
# Quick generation status checks - bash version
# Provides quick snapshots of generation progress

cd "$(dirname "$0")"
RESULTS_DIR="../../results/quaternary_boron_hydride"

print_section() {
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "════════════════════════════════════════════════════════════════"
}

print_section "QUICK STATUS SUMMARY"

# Count folders
total_folders=$(find "$RESULTS_DIR" -maxdepth 1 -type d -name '*_structures' 2>/dev/null | wc -l)
echo "Total generated folders:        $total_folders"

# Count complete (has extxyz)
complete=$(find "$RESULTS_DIR" -maxdepth 2 -name 'generated_crystals.extxyz' 2>/dev/null | wc -l)
echo "Complete (with extxyz):         $complete"

# Count with zip
with_zip=$(find "$RESULTS_DIR" -maxdepth 2 -name 'generated_crystals_cif.zip' 2>/dev/null | wc -l)
echo "With CIF ZIP:                   $with_zip"

# Count only partial (< 5 files in folder)
partial=0
for dir in "$RESULTS_DIR"/*_structures; do
    if [ -d "$dir" ]; then
        count=$(find "$dir" -maxdepth 1 -type f | wc -l)
        if [ "$count" -lt 5 ] && [ "$count" -gt 0 ]; then
            ((partial++))
        fi
    fi
done
echo "Partial/timeout folders:        $partial"

# Count empty (directories with 0 files)
empty=$(find "$RESULTS_DIR" -maxdepth 1 -type d -name '*_structures' 2>/dev/null | while read dir; do
  [ -z "$(find "$dir" -maxdepth 1 -type f 2>/dev/null)" ] && echo "$dir"
done | wc -l)
echo "Empty/failed folders:           $empty"

print_section "STORAGE STATISTICS"

total_size=$(du -sh "$RESULTS_DIR" 2>/dev/null | awk '{print $1}')
echo "Total size:                     $total_size"

avg_per_folder=$(find "$RESULTS_DIR" -maxdepth 1 -type d -name '*_structures' 2>/dev/null | \
    xargs -I {} du -s {} 2>/dev/null | awk -v n="$total_folders" '{sum+=$1} END {if(n>0) printf "%.1f MB\n", sum/n/1024; else print "N/A"}')
echo "Average per folder:             $avg_per_folder"

print_section "FILE STATISTICS"

extxyz_count=$(find "$RESULTS_DIR" -maxdepth 2 -name 'generated_crystals.extxyz' 2>/dev/null | wc -l)
echo "Files with extxyz:              $extxyz_count"

zip_count=$(find "$RESULTS_DIR" -maxdepth 2 -name 'generated_crystals_cif.zip' 2>/dev/null | wc -l)
echo "Files with CIF ZIP:             $zip_count"

both=$(comm -12 <(find "$RESULTS_DIR" -maxdepth 2 -name 'generated_crystals.extxyz' 2>/dev/null | sed 's|/generated_crystals.*||' | sort) \
                <(find "$RESULTS_DIR" -maxdepth 2 -name 'generated_crystals_cif.zip' 2>/dev/null | sed 's|/generated_crystals.*||' | sort) | wc -l)
echo "With both (extxyz + ZIP):       $both"

print_section "FAILURE ANALYSIS"

if [ -f "$RESULTS_DIR/failed_compositions.txt" ]; then
    failed=$(wc -l < "$RESULTS_DIR/failed_compositions.txt")
    echo "Failed compositions recorded:   $failed"
    echo "First 5 failed:"
    head -5 "$RESULTS_DIR/failed_compositions.txt" | sed 's/^/  - /'
    if [ "$failed" -gt 5 ]; then
        echo "  ... and $((failed - 5)) more"
    fi
else
    echo "No failed_compositions.txt found"
fi

print_section "MISSING COMPOSITIONS"

# Use Python to find missing ranges
python3 -c "
import json
from pathlib import Path

comps_file = Path('quaternary_boron_hydride_compositions.json')
results = Path('$RESULTS_DIR')

if comps_file.exists():
    with open(comps_file) as f:
        comps = json.load(f)
    
    total = len(comps)
    generated = len(list(results.glob('*_structures')))
    missing = total - generated
    pct = generated / total * 100 if total > 0 else 0
    
    print(f'Total compositions:            {total}')
    print(f'Generated (folder exists):     {generated} ({pct:.1f}%)')
    print(f'Missing (never attempted):     {missing} ({100-pct:.1f}%)')
" || echo "Could not analyze missing compositions"

print_section "RECOMMENDATIONS"

if [ "$total_folders" -lt 6000 ]; then
    echo "→ Run: ./submit_parallel_jobs.sh"
    echo "  to generate missing compositions"
fi

if [ "$partial" -gt 100 ]; then
    echo "→ Check partial folders for timeouts"
    echo "  May need longer timeout or smaller batch size"
fi

echo ""
echo "For detailed report:"
echo "  python3 generation_report.py"
echo ""
echo "For detailed coverage analysis:"
echo "  python3 check_coverage.py --show-ranges"
echo ""
