#!/usr/bin/env python3
"""
Comprehensive generation status report for quaternary hydrides.
Analyzes all compositions, directories, files, and provides detailed statistics.
"""
import json
import os
from pathlib import Path
from collections import defaultdict
import argparse
from datetime import datetime

def get_dir_size(path):
    """Get total size of directory in bytes."""
    total = 0
    try:
        for entry in path.rglob('*'):
            if entry.is_file():
                total += entry.stat().st_size
    except:
        pass
    return total

def format_size(bytes_size):
    """Format bytes to human readable size."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_size < 1024:
            return f"{bytes_size:.1f}{unit}"
        bytes_size /= 1024
    return f"{bytes_size:.1f}TB"

def categorize_folder(folder_path):
    """Categorize folder by its contents and status."""
    files = list(folder_path.glob('*'))
    file_names = {f.name for f in files}
    
    has_extxyz = 'generated_crystals.extxyz' in file_names
    has_zip = 'generated_crystals_cif.zip' in file_names
    has_structures = has_extxyz or has_zip
    file_count = len(files)
    
    if has_structures and file_count >= 5:
        return 'complete', f"{file_count} files"
    elif file_count > 0 and file_count < 5:
        return 'partial', f"{file_count} files (TIMEOUT/INCOMPLETE)"
    elif file_count == 0:
        return 'empty', "0 files (FAILED)"
    else:
        return 'unknown', f"{file_count} files"

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive generation report')
    parser.add_argument('--compositions', default='quaternary_hydride_compositions.json',
                        help='Path to compositions JSON file')
    parser.add_argument('--output-dir', default='../../results/quaternary_hydride',
                        help='Path to results directory')
    parser.add_argument('--failed-file', default='../../results/quaternary_hydride/failed_compositions.txt',
                        help='Path to failed compositions file')
    parser.add_argument('--save-report', 
                        help='Save report to file (default: print to stdout)')
    args = parser.parse_args()
    
    # Load compositions
    comp_path = Path(args.compositions)
    if not comp_path.exists():
        print(f"ERROR: Compositions file not found: {comp_path}")
        return 1
    
    with open(comp_path, 'r') as f:
        compositions = json.load(f)
    
    total_comps = len(compositions)
    
    # Check results directory
    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        print(f"ERROR: Output directory not found: {output_dir}")
        return 1
    
    # Load failed compositions if available
    failed_comps = set()
    failed_file = Path(args.failed_file)
    if failed_file.exists():
        with open(failed_file, 'r') as f:
            failed_comps = {line.strip() for line in f if line.strip()}
    
    # Analyze all compositions
    generated = {}
    missing = []
    failed = []
    status_counts = defaultdict(int)
    folder_stats = defaultdict(int)
    file_type_counts = {'extxyz': 0, 'zip': 0, 'both': 0}
    
    for idx, comp_data in enumerate(compositions):
        formula = comp_data['formula']
        struct_dir = output_dir / f"{formula}_structures"
        
        if formula in failed_comps:
            failed.append((idx, formula))
        elif struct_dir.exists():
            cat, desc = categorize_folder(struct_dir)
            status_counts[cat] += 1
            folder_stats[formula] = {
                'index': idx,
                'path': struct_dir,
                'status': cat,
                'description': desc,
                'size': get_dir_size(struct_dir),
                'files': len(list(struct_dir.glob('*')))
            }
            
            # Count file types
            has_extxyz = (struct_dir / 'generated_crystals.extxyz').exists()
            has_zip = (struct_dir / 'generated_crystals_cif.zip').exists()
            if has_extxyz and has_zip:
                file_type_counts['both'] += 1
            elif has_extxyz:
                file_type_counts['extxyz'] += 1
            elif has_zip:
                file_type_counts['zip'] += 1
            
            generated[formula] = (idx, struct_dir)
        else:
            missing.append((idx, formula))
    
    # Calculate statistics
    num_generated = len(generated)
    num_missing = len(missing)
    num_failed = len(failed)
    num_complete = status_counts['complete']
    num_partial = status_counts['partial']
    num_empty = status_counts['empty']
    
    total_size = sum(s['size'] for s in folder_stats.values())
    complete_size = sum(s['size'] for s in folder_stats.values() if s['status'] == 'complete')
    
    coverage_pct = (num_generated / total_comps) * 100
    success_rate = (num_complete / num_generated * 100) if num_generated > 0 else 0
    
    # Build report
    lines = []
    lines.append("=" * 80)
    lines.append("QUATERNARY HYDRIDE GENERATION - COMPREHENSIVE STATUS REPORT")
    lines.append("=" * 80)
    lines.append(f"Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    
    # Overview
    lines.append("OVERVIEW")
    lines.append("-" * 80)
    lines.append(f"Total compositions in file:     {total_comps:6d}")
    lines.append(f"Generated folders:              {num_generated:6d} ({coverage_pct:5.2f}%)")
    lines.append(f"Missing (never attempted):      {num_missing:6d} ({100-coverage_pct:5.2f}%)")
    lines.append(f"Failed (recorded failures):     {num_failed:6d}")
    lines.append(f"Total processed (gen+failed):   {num_generated + num_failed:6d}")
    lines.append("")
    
    # Breakdown of generated folders
    lines.append("BREAKDOWN OF GENERATED FOLDERS")
    lines.append("-" * 80)
    lines.append(f"Complete (has structures):      {num_complete:6d} ({num_complete/num_generated*100:5.2f}%) - Success")
    lines.append(f"Partial (< 5 files):            {num_partial:6d} ({num_partial/num_generated*100:5.2f}%) - Timeout/Incomplete")
    lines.append(f"Empty (0 files):                {num_empty:6d} ({num_empty/num_generated*100:5.2f}%) - Failed Generation")
    lines.append("")
    
    # File format statistics
    if num_complete > 0:
        lines.append("OUTPUT FILE TYPES (Complete folders)")
        lines.append("-" * 80)
        lines.append(f"With both CIF ZIP + XYZ:        {file_type_counts['both']:6d} ({file_type_counts['both']/num_complete*100:5.2f}%)")
        lines.append(f"XYZ only (no CIF ZIP):          {file_type_counts['extxyz']:6d} ({file_type_counts['extxyz']/num_complete*100:5.2f}%)")
        lines.append(f"CIF ZIP only (no XYZ):          {file_type_counts['zip']:6d} ({file_type_counts['zip']/num_complete*100:5.2f}%)")
        lines.append("")
    
    # Storage statistics
    lines.append("STORAGE STATISTICS")
    lines.append("-" * 80)
    lines.append(f"Total size of all folders:      {format_size(total_size):>10s}")
    lines.append(f"Average per generated folder:   {format_size(total_size/num_generated if num_generated > 0 else 0):>10s}")
    lines.append(f"Complete folders total:         {format_size(complete_size):>10s}")
    lines.append("")
    
    # Missing composition ranges
    if num_missing > 0:
        lines.append("MISSING COMPOSITION RANGES")
        lines.append("-" * 80)
        ranges = []
        if missing:
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
            
            if len(ranges) > 20:
                lines.append(f"Showing first 20 of {len(ranges)} ranges:")
            
            for i, (start, end) in enumerate(ranges[:20], 1):
                count = end - start + 1
                if start == end:
                    lines.append(f"  {i:2d}. Index {start:6d}: {compositions[start]['formula']}")
                else:
                    lines.append(f"  {i:2d}. Indices {start:6d}-{end:6d} ({count:4d} compositions)")
            
            if len(ranges) > 20:
                lines.append(f"  ... and {len(ranges) - 20} more ranges")
        lines.append("")
    
    # Failed compositions
    if num_failed > 0:
        lines.append("FAILED COMPOSITIONS (recorded failures)")
        lines.append("-" * 80)
        if num_failed <= 20:
            for idx, formula in sorted(failed)[:20]:
                lines.append(f"  [{idx:4d}] {formula}")
        else:
            for idx, formula in sorted(failed)[:10]:
                lines.append(f"  [{idx:4d}] {formula}")
            lines.append(f"  ... and {num_failed - 10} more")
        lines.append("")
    
    # Problematic folders (partial/empty)
    problematic = [(f, info) for f, info in folder_stats.items() 
                   if info['status'] in ['partial', 'empty']]
    if problematic:
        lines.append("PROBLEMATIC FOLDERS (Partial or Empty)")
        lines.append("-" * 80)
        if len(problematic) <= 15:
            for formula, info in sorted(problematic, key=lambda x: x[1]['index'])[:15]:
                lines.append(f"  [{info['index']:4d}] {formula:20s} - {info['status']:8s} ({info['description']})")
        else:
            for formula, info in sorted(problematic, key=lambda x: x[1]['index'])[:10]:
                lines.append(f"  [{info['index']:4d}] {formula:20s} - {info['status']:8s} ({info['description']})")
            lines.append(f"  ... and {len(problematic) - 10} more")
        lines.append("")
    
    # Recommendations
    lines.append("RECOMMENDATIONS")
    lines.append("-" * 80)
    if num_missing > 0:
        lines.append(f"1. Generate missing {num_missing} compositions:")
        lines.append(f"   ./submit_parallel_jobs.sh")
        lines.append(f"   Or manually generate ranges (see MISSING COMPOSITION RANGES above)")
        lines.append("")
    
    if num_partial > 0:
        lines.append(f"2. Check {num_partial} partial/timeout folders for debugging:")
        lines.append(f"   python3 -c \"import json; data=json.load(open('generation_report.json'))\"")
        lines.append(f"   May need to retry with longer timeout or smaller batch size")
        lines.append("")
    
    if num_empty > 0:
        lines.append(f"3. Investigate {num_empty} empty failed folders:")
        lines.append(f"   Check failed_compositions_detailed.txt for error details")
        lines.append(f"   May need model checkpoint update or composition validation")
        lines.append("")
    
    if num_generated == total_comps:
        lines.append("✓ ALL COMPOSITIONS PROCESSED!")
        if num_complete == total_comps:
            lines.append("✓ ALL COMPOSITIONS SUCCESSFULLY GENERATED!")
        else:
            lines.append(f"  Note: {total_comps - num_complete} folders need attention")
    else:
        pct_done = num_generated / total_comps * 100
        remaining = total_comps - num_generated
        lines.append(f"Progress: {pct_done:.1f}% complete ({remaining} compositions remaining)")
    
    lines.append("")
    lines.append("=" * 80)
    
    # Print and optionally save
    report_text = "\n".join(lines)
    print(report_text)
    
    if args.save_report:
        with open(args.save_report, 'w') as f:
            f.write(report_text)
        print(f"\nReport saved to: {args.save_report}")
    
    # Also save machine-readable JSON summary
    json_summary = {
        'timestamp': datetime.now().isoformat(),
        'total_compositions': total_comps,
        'generated': num_generated,
        'missing': num_missing,
        'failed': num_failed,
        'coverage_percent': coverage_pct,
        'status_breakdown': {
            'complete': num_complete,
            'partial': num_partial,
            'empty': num_empty
        },
        'file_types': file_type_counts,
        'storage': {
            'total_bytes': total_size,
            'total_formatted': format_size(total_size),
            'complete_bytes': complete_size,
            'complete_formatted': format_size(complete_size)
        },
        'missing_ranges': [(start, end) for start, end in ranges] if 'ranges' in locals() else [],
        'failed_compositions': [f for _, f in sorted(failed)] if failed else [],
        'problematic_folders': {f: {
            'index': info['index'],
            'status': info['status'],
            'files': info['files'],
            'size': info['size']
        } for f, info in problematic} if problematic else {}
    }
    
    json_file = output_dir / 'generation_report.json'
    with open(json_file, 'w') as f:
        json.dump(json_summary, f, indent=2)
    print(f"JSON summary saved to: {json_file}")
    
    return 0

if __name__ == '__main__':
    exit(main())
