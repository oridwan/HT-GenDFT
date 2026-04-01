#!/usr/bin/env python3
"""
Extract electride structures from PyXtal database to CIF files.

Usage:
    python3 extract_electride_struct.py --db electride_data.db --output-dir electride_structures
    
Features:
- Extracts only structures marked as electrides (is_electride=True)
- Exports each structure as structure_id.cif
- Creates output directory if it doesn't exist
- Reports statistics on extracted structures
"""

import os
import sys
import argparse
import sqlite3
import json
from pathlib import Path

try:
    from pyxtal.db import database_topology
    PYXTAL_AVAILABLE = True
except ImportError:
    print("ERROR: PyXtal not available. Please install: conda install -c conda-forge pyxtal")
    sys.exit(1)


def extract_electride_structures(db_path, output_dir, electrides_only=True):
    """
    Extract structures from PyXtal database to CIF files.
    
    Args:
        db_path: Path to PyXtal database file
        output_dir: Output directory for CIF files
        electrides_only: If True, only extract structures with is_electride=True
    
    Returns:
        tuple: (total_count, extracted_count, skipped_count)
    """
    if not os.path.exists(db_path):
        print(f"ERROR: Database not found: {db_path}")
        return 0, 0, 0
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    print(f"Loaded database: {db_path}")
    print(f"Output directory: {output_path}")
    print("")
    
    # Query database using sqlite3 directly (PyXtal database uses SQLite)
    try:
        # Connect to SQLite database directly
        conn = sqlite3.connect(str(db_path))
        cursor = conn.cursor()
        
        # Query the systems table (PyXtal/ASE database structure)
        if electrides_only:
            # Query only electrides - check key_value_pairs field for is_electride flag
            query = """
                SELECT id, unique_id, key_value_pairs FROM systems 
                WHERE key_value_pairs LIKE '%"is_electride": true%'
                ORDER BY unique_id
            """
            cursor.execute(query)
            print("Extracting electrides only (is_electride=True)...")
        else:
            # Query all structures
            query = "SELECT id, unique_id, key_value_pairs FROM systems ORDER BY unique_id"
            cursor.execute(query)
            print("Extracting all structures...")
        
        rows = cursor.fetchall()
        total_count = len(rows)
        
        if total_count == 0:
            cursor.close()
            conn.close()
            print("No structures found in database")
            return 0, 0, 0
        
        print(f"Found {total_count} structures to extract")
        print("")
        
        # Load PyXtal database for structure retrieval
        db = database_topology(db_path)
        
        extracted_count = 0
        skipped_count = 0
        
        for row in rows:
            db_id, unique_id, key_value_pairs = row
            
            try:
                # Parse key_value_pairs to get structure_id
                kvp_dict = json.loads(key_value_pairs) if key_value_pairs else {}
                structure_id = kvp_dict.get('structure_id', unique_id)
                
                # Get structure from database by ID
                xtal = db.get_pyxtal(db_id)
                
                if xtal is None:
                    print(f"  Warning: Could not retrieve structure {structure_id} (ID={db_id})")
                    skipped_count += 1
                    continue
                
                # Export to CIF
                cif_path = output_path / f"{structure_id}.cif"
                xtal.to_file(str(cif_path), fmt='cif')
                
                extracted_count += 1
                
                # Print progress every 100 structures
                if extracted_count % 100 == 0:
                    print(f"  Extracted {extracted_count}/{total_count} structures...")
                
            except Exception as e:
                print(f"  Error extracting {structure_id}: {e}")
                skipped_count += 1
                continue
        
        cursor.close()
        conn.close()
        
        return total_count, extracted_count, skipped_count
        
    except Exception as e:
        print(f"ERROR: Database query failed: {e}")
        import traceback
        traceback.print_exc()
        return 0, 0, 0


def main():
    parser = argparse.ArgumentParser(
        description="Extract electride structures from PyXtal database to CIF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract electrides only (default)
  python3 extract_electride_struct.py --db electride_data.db --output-dir electride_structures
  
  # Extract all structures (including non-electrides)
  python3 extract_electride_struct.py --db electride_data.db --output-dir all_structures --all
  
  # Custom output directory
  python3 extract_electride_struct.py --db ./VASP_JOBS/electride_data.db --output-dir ./CIF_files
"""
    )
    
    parser.add_argument(
        '--db',
        type=str,
        required=True,
        help="Path to PyXtal database file (e.g., electride_data.db)"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./electride_structures',
        help="Output directory for CIF files (default: ./electride_structures)"
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help="Extract all structures (default: electrides only)"
    )
    
    args = parser.parse_args()
    
    print("="*70)
    print("Extract Electride Structures from PyXtal Database")
    print("="*70)
    print(f"Database: {args.db}")
    print(f"Output directory: {args.output_dir}")
    print(f"Mode: {'All structures' if args.all else 'Electrides only'}")
    print("="*70)
    print("")
    
    # Extract structures
    electrides_only = not args.all
    total, extracted, skipped = extract_electride_structures(
        args.db,
        args.output_dir,
        electrides_only=electrides_only
    )
    
    # Print summary
    print("")
    print("="*70)
    print("Extraction Summary")
    print("="*70)
    print(f"Total structures in database: {total}")
    print(f"Successfully extracted: {extracted}")
    print(f"Skipped/Failed: {skipped}")
    print("")
    print(f"CIF files saved to: {args.output_dir}")
    print("="*70)
    
    if extracted == 0 and total > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()

