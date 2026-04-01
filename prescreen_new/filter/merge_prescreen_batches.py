#!/usr/bin/env python3
"""
Merge parallel prescreening batch results into a single file.

Usage:
    python3 merge_prescreen_batches.py --output-dir ./VASP_JOBS [--skip-database]
"""

import argparse
import json
import shutil
import sqlite3
from pathlib import Path
from typing import Dict, List, Any
import sys


def load_batch_results(output_dir: Path) -> List[Dict[str, Any]]:
    """
    Load all prescreening_stability_batch*.json files.
    
    Args:
        output_dir: Directory containing batch files
        
    Returns:
        List of batch result dictionaries
    """
    batch_files = sorted(output_dir.glob('prescreening_stability_batch*.json'))
    
    if not batch_files:
        print(f"WARNING: No batch JSON files found in {output_dir}")
        print("Expected files matching: prescreening_stability_batch*.json")
        print("Proceeding without JSON merge.")
        return []
    
    print(f"Found {len(batch_files)} batch files:")
    for bf in batch_files:
        print(f"  - {bf.name}")
    print()
    
    batch_results = []
    for bf in batch_files:
        try:
            with open(bf, 'r') as f:
                data = json.load(f)
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
    
    # Use first batch as template
    merged = {
        'summary': batch_results[0]['summary'].copy(),
        'results': []
    }
    
    # Track unique structures by composition + structure_id
    seen_structures = set()
    duplicate_count = 0
    
    for batch_data in batch_results:
        for result in batch_data.get('results', []):
            composition = result['composition']
            struct_id = result['structure_id']
            key = (composition, struct_id)
            
            if key in seen_structures:
                duplicate_count += 1
                continue
            
            seen_structures.add(key)
            merged['results'].append(result)
    
    # Sort by composition, then structure_id
    merged['results'].sort(key=lambda x: (x['composition'], x['structure_id']))
    
    # Update summary
    merged['summary']['total_batches'] = len(batch_results)
    merged['summary']['total_structures'] = len(merged['results'])
    
    print(f"Merged {len(batch_results)} batches:")
    print(f"  Total unique structures: {len(merged['results'])}")
    print(f"  Duplicates skipped: {duplicate_count}")
    
    # Count structures that passed prescreening
    passed_count = sum(1 for r in merged['results'] if r.get('passed_prescreening', False))
    failed_count = len(merged['results']) - passed_count
    
    merged['summary']['passed_prescreening'] = passed_count
    merged['summary']['failed_prescreening'] = failed_count
    
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

        # Helper mapping table used to remap source IDs -> merged IDs.
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS _id_map_work (old_id INTEGER PRIMARY KEY, new_id INTEGER NOT NULL)"
        )

        # Build column list dynamically so merge remains robust to schema drift.
        cursor.execute("PRAGMA table_info(systems)")
        systems_columns = [row[1] for row in cursor.fetchall() if row[1] != "id"]
        systems_col_sql = ", ".join(systems_columns)
        systems_select_sql = ", ".join([f"s.{col}" for col in systems_columns])

        # Merge remaining databases using ID remapping.
        if len(batch_dbs) > 1:
            for i, batch_db in enumerate(batch_dbs[1:], 1):
                print(f"Merging {batch_db.name}...")
                batch_alias = f"batch_{i}"
                attached = False

                try:
                    # ATTACH doesn't support parameterized DB name.
                    cursor.execute(f"ATTACH DATABASE '{batch_db}' AS {batch_alias}")
                    attached = True

                    # Track current max ID so we can map only newly inserted systems.
                    cursor.execute("SELECT COALESCE(MAX(id), 0) FROM systems")
                    max_id_before = cursor.fetchone()[0]

                    # Insert systems rows without source IDs to avoid PK collisions.
                    cursor.execute(
                        f"""
                        INSERT OR IGNORE INTO systems ({systems_col_sql})
                        SELECT {systems_select_sql}
                        FROM {batch_alias}.systems AS s
                        """
                    )
                    cursor.execute("SELECT changes()")
                    batch_structures = cursor.fetchone()[0]

                    # Build old_id -> new_id mapping for newly inserted rows only.
                    # Avoiding a cross-database join here is more robust for large batches.
                    cursor.execute("SELECT id, unique_id FROM systems WHERE id > ?", (max_id_before,))
                    new_id_by_unique_id = {uid: new_id for new_id, uid in cursor.fetchall()}

                    cursor.execute(f"SELECT id, unique_id FROM {batch_alias}.systems")
                    id_map_rows = []
                    for old_id, unique_id in cursor.fetchall():
                        new_id = new_id_by_unique_id.get(unique_id)
                        if new_id is not None:
                            id_map_rows.append((old_id, new_id))

                    cursor.execute("DELETE FROM _id_map_work")
                    if id_map_rows:
                        cursor.executemany(
                            "INSERT INTO _id_map_work (old_id, new_id) VALUES (?, ?)",
                            id_map_rows,
                        )

                    # Copy dependent tables with remapped IDs.
                    cursor.execute(
                        f"""
                        INSERT INTO species (Z, n, id)
                        SELECT t.Z, t.n, m.new_id
                        FROM {batch_alias}.species AS t
                        JOIN _id_map_work AS m ON t.id = m.old_id
                        """
                    )
                    cursor.execute("SELECT changes()")
                    added_species = cursor.fetchone()[0]

                    cursor.execute(
                        f"""
                        INSERT INTO keys ("key", id)
                        SELECT t."key", m.new_id
                        FROM {batch_alias}.keys AS t
                        JOIN _id_map_work AS m ON t.id = m.old_id
                        """
                    )
                    cursor.execute("SELECT changes()")
                    added_keys = cursor.fetchone()[0]

                    cursor.execute(
                        f"""
                        INSERT INTO text_key_values ("key", value, id)
                        SELECT t."key", t.value, m.new_id
                        FROM {batch_alias}.text_key_values AS t
                        JOIN _id_map_work AS m ON t.id = m.old_id
                        """
                    )
                    cursor.execute("SELECT changes()")
                    added_text_kv = cursor.fetchone()[0]

                    cursor.execute(
                        f"""
                        INSERT INTO number_key_values ("key", value, id)
                        SELECT t."key", t.value, m.new_id
                        FROM {batch_alias}.number_key_values AS t
                        JOIN _id_map_work AS m ON t.id = m.old_id
                        """
                    )
                    cursor.execute("SELECT changes()")
                    added_number_kv = cursor.fetchone()[0]

                    conn.commit()
                    cursor.execute(f"DETACH DATABASE {batch_alias}")
                    attached = False

                    print(f"  Added {batch_structures} entries to table 'systems'")
                    print(f"  Added {added_species} entries to table 'species'")
                    print(f"  Added {added_keys} entries to table 'keys'")
                    print(f"  Added {added_text_kv} entries to table 'text_key_values'")
                    print(f"  Added {added_number_kv} entries to table 'number_key_values'")
                    print(f"  → {batch_structures} new structures from {batch_db.name}\n")

                except sqlite3.Error as e:
                    print(f"  ERROR: {e}")
                    print(f"  This database may be corrupted (check for batch job crash)\n")
                    try:
                        conn.rollback()
                        if attached:
                            cursor.execute(f"DETACH DATABASE {batch_alias}")
                    except Exception:
                        pass
                    continue
        
        cursor.execute("DROP TABLE IF EXISTS _id_map_work")
        conn.commit()

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
    
    # Load and merge batch JSON results (if present)
    batch_results = load_batch_results(output_dir)
    if batch_results:
        merged_data = merge_batch_results(batch_results)

        # Write merged file
        output_filename = args.output_file or 'prescreening_stability.json'
        output_file = output_dir / output_filename

        with open(output_file, 'w') as f:
            json.dump(merged_data, f, indent=2)

        print(f"\nMerged results saved to: {output_file}")
    else:
        print("\nSkipping JSON merge/write (no batch JSON files found).")
    
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
    
    if db_success:
        print("\nBatch and database files kept (no cleanup performed)")
    else:
        print("\nBatch and database files kept due to merge failure")
    
    print()
    print("="*70)
    print("Merge complete!")
    print("="*70)


if __name__ == '__main__':
    main()
