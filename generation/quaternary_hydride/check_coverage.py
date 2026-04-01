#!/usr/bin/env python3
"""
Check which compositions from the JSON file have been generated.
Reports missing compositions and statistics.
"""
import json
import os
from pathlib import Path
import argparse

def main():
    parser = argparse.ArgumentParser(description='Check structure generation coverage')
    parser.add_argument('--compositions', default='quaternary_hydride_compositions.json',
                        help='Path to compositions JSON file')
    parser.add_argument('--output-dir', default='../../results/quaternary_hydride',
                        help='Path to results directory')
    parser.add_argument('--show-missing', action='store_true',
                        help='Print list of missing compositions')
    parser.add_argument('--show-ranges', action='store_true',
                        help='Show missing composition index ranges')
    args = parser.parse_args()
    
    # Load compositions
    with open(args.compositions, 'r') as f:
        compositions = json.load(f)
    
    total = len(compositions)
    print(f"Total compositions in file: {total}")
    print(f"Index range: 0 to {total - 1}")
    print()
    
    # Check which have been generated
    output_dir = Path(args.output_dir)
    generated = set()
    missing = []
    
    for idx, comp_data in enumerate(compositions):
        formula = comp_data['formula']
        struct_dir = output_dir / f"{formula}_structures"
        
        if struct_dir.exists():
            generated.add(idx)
        else:
            missing.append((idx, formula))
    
    # Statistics
    num_generated = len(generated)
    num_missing = len(missing)
    coverage_pct = (num_generated / total) * 100
    
    print("=" * 60)
    print("COVERAGE SUMMARY")
    print("=" * 60)
    print(f"Total compositions:     {total:6d}")
    print(f"Generated:              {num_generated:6d} ({coverage_pct:5.2f}%)")
    print(f"Missing:                {num_missing:6d} ({100-coverage_pct:5.2f}%)")
    print("=" * 60)
    print()
    
    if num_missing == 0:
        print("✓ All compositions have been generated!")
        return 0
    
    # Show missing ranges
    if args.show_ranges or num_missing <= 50:
        print("Missing composition index ranges:")
        ranges = []
        start = missing[0][0]
        prev = missing[0][0]
        
        for idx, formula in missing[1:]:
            if idx == prev + 1:
                prev = idx
            else:
                ranges.append((start, prev))
                start = idx
                prev = idx
        ranges.append((start, prev))
        
        for start, end in ranges:
            if start == end:
                print(f"  Index {start}: {compositions[start]['formula']}")
            else:
                count = end - start + 1
                print(f"  Indices {start}-{end} ({count} compositions)")
        print()
    
    # Show individual missing compositions
    if args.show_missing:
        print(f"Missing compositions (showing first 100):")
        for idx, formula in missing[:100]:
            print(f"  [{idx:4d}] {formula}")
        if num_missing > 100:
            print(f"  ... and {num_missing - 100} more")
        print()
    
    # Provide recommendation
    print("To generate missing compositions:")
    if num_missing > 0:
        first_missing = missing[0][0]
        last_missing = missing[-1][0]
        
        if len(ranges) == 1 and ranges[0][0] == ranges[0][1]:
            # Single missing composition
            print(f"  python generate_structures_batch.py \\")
            print(f"    --compositions quaternary_hydride_compositions.json \\")
            print(f"    --output-dir ../../results/quaternary_hydride \\")
            print(f"    --model ../../MatterGen_checkpoints/18-55-08 \\")
            print(f"    --start-index {first_missing} \\")
            print(f"    --max-compositions 1")
        elif len(ranges) == 1:
            # Single contiguous range
            count = ranges[0][1] - ranges[0][0] + 1
            print(f"  python generate_structures_batch.py \\")
            print(f"    --compositions quaternary_hydride_compositions.json \\")
            print(f"    --output-dir ../../results/quaternary_hydride \\")
            print(f"    --model ../../MatterGen_checkpoints/18-55-08 \\")
            print(f"    --start-index {first_missing} \\")
            print(f"    --max-compositions {count}")
        else:
            # Multiple ranges - suggest reprocessing all or using array job
            print(f"  # Option 1: Reprocess all (skips existing)")
            print(f"  python generate_structures_batch.py \\")
            print(f"    --compositions quaternary_hydride_compositions.json \\")
            print(f"    --output-dir ../../results/quaternary_hydride \\")
            print(f"    --model ../../MatterGen_checkpoints/18-55-08 \\")
            print(f"    --start-index 0 \\")
            print(f"    --max-compositions {total}")
            print()
            print(f"  # Option 2: Process each range separately")
            for start, end in ranges[:5]:
                count = end - start + 1
                print(f"  # Range {start}-{end}:")
                print(f"  python generate_structures_batch.py ... --start-index {start} --max-compositions {count}")
    
    return 0

if __name__ == '__main__':
    exit(main())
