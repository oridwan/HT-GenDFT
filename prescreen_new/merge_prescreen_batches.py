#!/usr/bin/env python3
"""
Merge parallel prescreening batch results into a single file.

Usage:
python3 merge_prescreen_batches.py \
  --output-dir VASP_JOBS_Boron \
  --skip-database \
  --keep-batches \
  --output-file prescreening_stability.json
"""

import argparse
import json
import shutil
import sqlite3
from pathlib import Path
from typing import Dict, List, Any
import sys
import re


def _batch_sort_key(path: Path):
    """Sort batch files numerically when names include 'batch<N>'."""
    match = re.search(r"batch(\d+)", path.name)
    if match:
        return (0, int(match.group(1)))
    return (1, path.name)


def load_batch_results(output_dir: Path) -> List[Dict[str, Any]]:
    """
    Load all prescreening_stability_batch*.json files.
    
    Args:
        output_dir: Directory containing batch files
        
    Returns:
        List of batch result dictionaries
    """
    batch_files = sorted(output_dir.glob('prescreening_stability_batch*.json'), key=_batch_sort_key)
    
    if not batch_files:
        print(f"ERROR: No batch files found in {output_dir}")
        print("Expected files matching: prescreening_stability_batch*.json")
        sys.exit(1)
    
    print(f"Found {len(batch_files)} batch files:")
    for bf in batch_files:
        print(f"  - {bf.name}")
    print()
    
    batch_results = []
    for bf in batch_files:
        try:
            with open(bf, 'r') as f:
                data = json.load(f)
                if not isinstance(data, dict):
                    print(f"WARNING: Skipping {bf.name} (top-level JSON is not an object)")
                    continue
                if 'results' not in data or not isinstance(data.get('results'), list):
                    print(f"WARNING: Skipping {bf.name} (missing/invalid 'results' list)")
                    continue
                data['_source_file'] = bf.name
                batch_results.append(data)
        except Exception as e:
            print(f"WARNING: Failed to load {bf.name}: {e}")
            continue
    
    return batch_results


def merge_batch_results(batch_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Merge batch results into a single result dictionary.
    
    Args:
        batch_results: List of batch result dictionaries
        
    Returns:
        Merged result dictionary
    """
    if not batch_results:
        print("ERROR: No valid batch results to merge")
        sys.exit(1)
    
    first_summary = batch_results[0].get('summary', {})
    if not isinstance(first_summary, dict):
        first_summary = {}

    merged = {
        'summary': first_summary.copy(),
        'results': []
    }
    
    # Track unique structures by structure_id (+ composition safeguard)
    seen_structures = set()
    duplicate_count = 0
    total_result_records = 0
    source_files = []

    for batch_data in batch_results:
        source_name = batch_data.get('_source_file')
        if source_name:
            source_files.append(source_name)
        for result in batch_data.get('results', []):
            total_result_records += 1
            composition = result.get('composition')
            struct_id = result.get('structure_id')

            # Skip malformed rows that can't be identified
            if struct_id is None and composition is None:
                continue

            key = (composition, struct_id)
            
            if key in seen_structures:
                duplicate_count += 1
                continue
            
            seen_structures.add(key)
            merged['results'].append(result)
    
    # Sort by composition, then structure_id
    merged['results'].sort(key=lambda x: (
        str(x.get('composition', '')),
        str(x.get('structure_id', ''))
    ))
    
    # Aggregate summary statistics from all batches
    def _sum_summary_int(key: str) -> int:
        total = 0
        for batch in batch_results:
            summary = batch.get('summary', {})
            if not isinstance(summary, dict):
                continue
            value = summary.get(key, 0)
            try:
                total += int(value)
            except Exception:
                continue
        return total

    total_loaded = _sum_summary_int('total_structures_loaded')
    total_duplicates = _sum_summary_int('duplicate_structures_removed')
    total_unique = _sum_summary_int('unique_structures_processed')
    
    # Count structures that passed prescreening
    passed_count = sum(1 for r in merged['results'] if r.get('passed_prescreening', False))
    failed_count = len(merged['results']) - passed_count
    
    # Count duplicates in merged results (structures with 'duplicate_of' field)
    duplicate_records_count = sum(1 for r in merged['results'] if 'duplicate_of' in r)
    
    # Update summary
    merged['summary']['total_batches'] = len(batch_results)
    merged['summary']['source_batch_files'] = source_files
    merged['summary']['total_structures_loaded'] = total_loaded
    merged['summary']['duplicate_structures_removed'] = total_duplicates
    merged['summary']['unique_structures_processed'] = total_unique
    merged['summary']['total_structures'] = len(merged['results'])  # Keep for compatibility
    merged['summary']['passed_prescreening'] = passed_count
    merged['summary']['failed_prescreening'] = failed_count
    merged['summary']['total_result_records_before_dedup'] = total_result_records
    merged['summary']['cross_batch_duplicates_skipped'] = duplicate_count
    
    print(f"Merged {len(batch_results)} batches:")
    print(f"  Total structures loaded: {total_loaded}")
    print(f"  Duplicates removed: {total_duplicates} ({duplicate_records_count} duplicate records in merged results)")
    print(f"  Unique structures processed: {total_unique}")
    print(f"  Total result records: {len(merged['results'])} (from {total_result_records} input rows)")
    print(f"  Cross-batch duplicates skipped: {duplicate_count}")
    print(f"  Passed prescreening: {passed_count}")
    print(f"  Failed prescreening: {failed_count}")
    
    return merged


def check_mp_cache(output_dir: Path) -> bool:
    """
    Check if shared MP cache exists (no merging needed - batches use shared file).
    
    Args:
        output_dir: Directory containing MP cache
        
    Returns:
        bool: True if cache exists
    """
    mp_cache = output_dir / 'mp_mattersim.json'
    
    if not mp_cache.exists():
        print("\nNo MP cache file found (mp_mattersim.json)")
        print("  Note: Parallel batches share a single MP cache with file locking")
        return False
    
    print(f"\n{'='*70}")
    print("Checking Shared MP MatterSim Cache")
    print(f"{'='*70}")
    
    try:
        with open(mp_cache, 'r') as f:
            cache_data = json.load(f)
        
        print(f"MP cache: {mp_cache}")
        print(f"  Total entries: {len(cache_data)}")
        
        # Show breakdown by chemical system
        chemsys_counts = {}
        for entry in cache_data:
            chemsys = entry.get('chemsys', 'unknown')
            chemsys_counts[chemsys] = chemsys_counts.get(chemsys, 0) + 1
        
        print(f"  Breakdown by chemical system:")
        for chemsys in sorted(chemsys_counts.keys())[:10]:  # Show first 10
            print(f"    {chemsys}: {chemsys_counts[chemsys]} phases")
        if len(chemsys_counts) > 10:
            print(f"    ... and {len(chemsys_counts) - 10} more systems")
        
        return True
    except Exception as e:
        print(f"  ERROR: Could not read MP cache: {e}")
        return False


def merge_pyxtal_databases(output_dir: Path) -> bool:
    """
    Merge PyXtal database files from parallel batches using direct SQL.
    
    PyXtal wraps ASE database, which is SQLite underneath. For bulk operations,
    we use SQL directly for speed, avoiding expensive per-structure operations.
    
    Args:
        output_dir: Directory containing batch database files
        
    Returns:
        bool: True if successful, False otherwise
    """
    batch_dbs = sorted(output_dir.glob('prescreening_structures_batch*.db'))
    
    if not batch_dbs:
        print("\nNo batch database files found to merge")
        return False
    
    print(f"\n{'='*70}")
    print("Merging PyXtal Database Files")
    print(f"{'='*70}")
    print(f"Found {len(batch_dbs)} database files:")
    for db in batch_dbs:
        print(f"  - {db.name}")
    print()
    
    merged_db_path = output_dir / 'prescreening_structures.db'
    
    # Remove existing merged database if it exists
    if merged_db_path.exists():
        merged_db_path.unlink()
        print(f"Removed existing merged database\n")
    
    try:
        # Copy first batch as base (fast file copy)
        shutil.copy2(batch_dbs[0], merged_db_path)
        print(f"Initialized merged database from {batch_dbs[0].name}")
        
        # Get count from first batch
        conn = sqlite3.connect(str(merged_db_path))
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM systems")
        initial_count = cursor.fetchone()[0]
        print(f"  Initial structures: {initial_count}\n")
        
        # Merge remaining databases using SQL (fast bulk operations)
        if len(batch_dbs) > 1:
            for i, batch_db in enumerate(batch_dbs[1:], 1):
                print(f"Merging {batch_db.name}...")
                
                try:
                    # Use unique alias for each batch to avoid conflicts
                    batch_alias = f"batch_{i}_{batch_db.stem}"
                    
                    # Attach the batch database (use string formatting, ATTACH doesn't support parameterized queries)
                    cursor.execute(f"ATTACH DATABASE '{batch_db}' AS {batch_alias}")
                    
                    # Get table names from batch database
                    cursor.execute(f"SELECT name FROM {batch_alias}.sqlite_master WHERE type='table'")
                    tables = cursor.fetchall()
                    
                    batch_structures = 0
                    for (table_name,) in tables:
                        # Copy data using INSERT OR IGNORE (skips duplicates automatically)
                        try:
                            cursor.execute(f"INSERT OR IGNORE INTO {table_name} SELECT * FROM {batch_alias}.{table_name}")
                            added = cursor.rowcount
                            if table_name == 'systems':
                                batch_structures = added
                            print(f"  Added {added} entries to table '{table_name}'")
                        except sqlite3.Error as e:
                            print(f"  Warning: Could not merge table '{table_name}': {e}")
                    
                    # Commit before detaching (critical to release locks)
                    conn.commit()
                    
                    # Now detach the batch database
                    cursor.execute(f"DETACH DATABASE {batch_alias}")
                    
                    print(f"  → {batch_structures} new structures from {batch_db.name}\n")
                    
                except sqlite3.Error as e:
                    print(f"  ERROR: {e}")
                    print(f"  This database may be corrupted (check for batch job crash)\n")
                    # Try to commit and detach if attached
                    try:
                        conn.commit()
                        if 'batch_alias' in locals():
                            cursor.execute(f"DETACH DATABASE {batch_alias}")
                    except:
                        pass
                    continue
        
        # Get final count
        cursor.execute("SELECT COUNT(*) FROM systems")
        final_count = cursor.fetchone()[0]
        conn.close()
        
        print(f"Database merge complete!")
        print(f"  Total structures in merged database: {final_count}")
        print(f"  Merged database: {merged_db_path}")
        return True
        
    except Exception as e:
        print(f"\nERROR merging databases: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Merge parallel prescreening batch results"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./VASP_JOBS',
        help="Directory containing batch result files"
    )
    parser.add_argument(
        '--keep-batches',
        action='store_true',
        help="Keep batch files after merging (default: delete them)"
    )
    parser.add_argument(
        '--output-file',
        type=str,
        default=None,
        help="Output filename (default: prescreening_stability.json)"
    )
    parser.add_argument(
        '--skip-database',
        action='store_true',
        help="Skip merging PyXtal database files (only merge JSON)"
    )
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir).expanduser()
    
    if not output_dir.exists():
        print(f"ERROR: Output directory not found: {output_dir}")
        sys.exit(1)
    
    print("="*70)
    print("Merging Prescreening Batch Results")
    print("="*70)
    print(f"Output directory: {output_dir}\n")
    
    # Load all batch results
    batch_results = load_batch_results(output_dir)
    
    # Merge results
    merged_data = merge_batch_results(batch_results)
    
    # Write merged file
    output_filename = args.output_file or 'prescreening_stability.json'
    output_file = output_dir / output_filename
    
    with open(output_file, 'w') as f:
        json.dump(merged_data, f, indent=2)
    
    print(f"\nMerged results saved to: {output_file}")
    
    # Check shared MP cache (no merging needed - file is shared across batches)
    print()
    check_mp_cache(output_dir)
    
    # Merge PyXtal databases unless skipped
    db_success = False
    if not args.skip_database:
        db_success = merge_pyxtal_databases(output_dir)
        if not db_success:
            print("\nWARNING: Database merge failed or was skipped")
            print("  Batch database files will be kept for manual inspection")
    else:
        print("\nSkipping database merge (--skip-database)")
        db_success = True  # Consider it successful if skipped
    
    # Cleanup batch files ONLY if merge was successful
    if not args.keep_batches and db_success:
        print("\nCleaning up batch files...")
        batch_files = list(output_dir.glob('prescreening_stability_batch*.json'))
        checkpoint_files = list(output_dir.glob('prescreening_checkpoint_batch*.json'))
        db_files = list(output_dir.glob('prescreening_structures_batch*.db'))
        mp_lock_files = list(output_dir.glob('mp_mattersim.lock'))
        chemsys_lock_files = list(output_dir.glob('mp_cache_*.lock'))
        
        all_files = batch_files + checkpoint_files + db_files + mp_lock_files + chemsys_lock_files
        for bf in all_files:
            try:
                bf.unlink()
                print(f"  Deleted: {bf.name}")
            except Exception as e:
                print(f"  WARNING: Failed to delete {bf.name}: {e}")
        
        print(f"Deleted {len(all_files)} batch and checkpoint files")
    elif not args.keep_batches and not db_success:
        print("\nBatch files KEPT due to merge failure")
        print("  Fix the database issues, then re-run merge with --keep-batches")
    else:
        print("\nBatch files kept (--keep-batches flag set)")
    
    print()
    print("="*70)
    print("Merge complete!")
    print("="*70)


if __name__ == '__main__':
    main()
