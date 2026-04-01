#!/usr/bin/env python3
"""
Filter electride database with stricter interstitial volume criteria.

Creates electride_candidates.db with:
- Only structures meeting volume thresholds
- Symmetrized structures (PyXtal)
- Interstitial volume data (e0025, e05, e10, band0, band1)
- Thermodynamic stability (e_above_hull)
- Duplicate removal
- Sorted by space group (descending - high symmetry first)

Also generates electride_candidates.csv for quick viewing.

Usage:
    python3 filter_comb_db.py --input electride_data.db --csv electride_analysis.csv
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

try:
    from pyxtal.db import database_topology
    from pyxtal import pyxtal
    from pyxtal.util import new_struc_wo_energy
except ImportError as e:
    print(f"ERROR: Required package not available: {e}")
    print("Install with: conda install -c conda-forge pyxtal")
    sys.exit(1)


def filter_and_create_db(input_db_path, csv_path, output_db_path, 
                         min_energy_volume=20.0, min_band_volume=20.0):
    """
    Filter electride database with stricter volume criteria.
    
    Args:
        input_db_path: Input PyXtal database (electride_data.db)
        csv_path: CSV with interstitial volumes (electride_analysis.csv)
        output_db_path: Output PyXtal database (electride_candidates.db)
        min_energy_volume: Minimum for max(e0025, e05, e10)
        min_band_volume: Minimum for max(band0, band1)
    
    Returns:
        tuple: (total_count, added_count, duplicate_count, filtered_count, failed_count)
    """
    print("="*70)
    print("Filtering Electride Database with Stricter Criteria")
    print("="*70)
    print(f"Input database: {input_db_path}")
    print(f"Input CSV: {csv_path}")
    print(f"Output database: {output_db_path}")
    print(f"Criteria:")
    print(f"  max(e0025, e05, e10) >= {min_energy_volume} Å³")
    print(f"  max(band0, band1) >= {min_band_volume} Å³")
    print("="*70)
    print("")
    
    # Load CSV with interstitial volumes
    print("Loading CSV data...")
    df = pd.read_csv(csv_path)
    print(f"  Loaded {len(df)} structures from CSV")
    
    # Load input database
    print("Loading input database...")
    db = database_topology(input_db_path)
    total_count = db.db.count()
    print(f"  Database contains {total_count} structures")
    print("")
    
    # Create output database (delete if exists)
    output_path = Path(output_db_path)
    if output_path.exists():
        print(f"Removing existing output database: {output_db_path}")
        output_path.unlink()
    
    db_out = database_topology(output_db_path)
    
    # Track duplicates
    xtals = []
    
    structures_to_add = []
    
    # Counters
    duplicate_count = 0
    failed_count = 0
    filtered_count = 0
    
    # Tolerances for PyXtal symmetrization (coarse to fine)
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    
    print("Processing structures...")
    
    for row in db.db.select():
        # Get structure_id directly (works for both PyXtal and ASE structures)
        try:
            struct_id = row.structure_id
        except AttributeError:
            print(f"  Warning: Skipping row {row.id} (no structure_id attribute)")
            failed_count += 1
            continue
        
        # Find corresponding CSV row
        df_row = df[df['formula'] == struct_id]
        
        if len(df_row) == 0:
            print(f"  Warning: {struct_id} not found in CSV")
            failed_count += 1
            continue
        
        csv_row = df_row.iloc[0]
        
        # Apply stricter filtering
        energy_volume = max(csv_row['e0025'], csv_row['e05'], csv_row['e10'])
        band_volume = max(csv_row['band0'], csv_row['band1'])
        
        if energy_volume < min_energy_volume:
            filtered_count += 1
            continue
        
        if band_volume < min_band_volume:
            filtered_count += 1
            continue
        
        # Get atoms object
        try:
            atoms = row.toatoms()
        except Exception as e:
            print(f"  Warning: Skipping {struct_id} (could not get atoms: {e})")
            failed_count += 1
            continue
        
        # Try to symmetrize with progressive tolerances
        xtal = pyxtal()
        symmetrized = False
        
        for tol in tolerances:
            try:
                xtal.from_seed(atoms, tol=tol)
                
                if not xtal.valid:
                    continue
                
                # Check for duplicates
                if new_struc_wo_energy(xtal, xtals):
                    # Store structure for later sorting
                    dicts = {
                        'e0025': float(csv_row['e0025']),
                        'e05': float(csv_row['e05']),
                        'e10': float(csv_row['e10']),
                        'band0': float(csv_row['band0']),
                        'band1': float(csv_row['band1']),
                        'e_above_hull': float(csv_row['e_above_hull']),
                        'structure_id': struct_id,
                        'composition': csv_row['composition'],
                        'is_electride': True,
                    }
                    
                    space_group_num = xtal.group.number
                    structures_to_add.append({
                        'xtal': xtal,
                        'dicts': dicts,
                        'space_group': space_group_num,
                        'is_pyxtal': True
                    })
                    
                    xtals.append(xtal)
                    
                    if len(structures_to_add) % 100 == 0:
                        print(f"  Collected {len(structures_to_add)} structures...")
                    
                    symmetrized = True
                    break
                else:
                    duplicate_count += 1
                    symmetrized = True
                    break
                    
            except Exception as e:
                continue
        
        if not symmetrized:
            # Save structure without symmetrization (ASE atoms directly)
            try:
                # All structures in filtered DB pass strict criteria, so is_electride=True
                # Check for duplicates using composition and all 5 interstitial volumes
                composition = csv_row['composition']
                e0025 = float(csv_row['e0025'])
                e05 = float(csv_row['e05'])
                e10 = float(csv_row['e10'])
                band0 = float(csv_row['band0'])
                band1 = float(csv_row['band1'])
                e_hull = float(csv_row['e_above_hull'])
                
                is_duplicate = False
                for prev_row in db_out.db.select():
                    try:
                        prev_comp = getattr(prev_row, 'composition', None)
                        if prev_comp is None or prev_comp != composition:
                            continue
                        
                        # Check if all 5 interstitial volumes are close (within 0.5 Å³)
                        prev_e0025 = float(getattr(prev_row, 'e0025', -999))
                        prev_e05 = float(getattr(prev_row, 'e05', -999))
                        prev_e10 = float(getattr(prev_row, 'e10', -999))
                        prev_band0 = float(getattr(prev_row, 'band0', -999))
                        prev_band1 = float(getattr(prev_row, 'band1', -999))
                        
                        if (abs(e0025 - prev_e0025) < 0.5 and
                            abs(e05 - prev_e05) < 0.5 and
                            abs(e10 - prev_e10) < 0.5 and
                            abs(band0 - prev_band0) < 0.5 and
                            abs(band1 - prev_band1) < 0.5):
                            is_duplicate = True
                            break
                    except (AttributeError, ValueError):
                        continue
                
                if not is_duplicate:
                    # Collect structure for later sorting
                    structures_to_add.append({
                        'atoms': atoms,
                        'struct_id': struct_id,
                        'composition': composition,
                        'e0025': e0025,
                        'e05': e05,
                        'e10': e10,
                        'band0': band0,
                        'band1': band1,
                        'e_above_hull': e_hull,
                        'is_electride': True,
                        'space_group': 1,
                        'is_pyxtal': False
                    })
                    
                    if len(structures_to_add) % 100 == 0:
                        print(f"  Collected {len(structures_to_add)} structures...")
                    
                    print(f"  Collected without symmetrization: {struct_id}")
                else:
                    duplicate_count += 1
            except Exception as e:
                print(f"  Warning: Failed to save {struct_id}: {e}")
                failed_count += 1
    
    print(f"  Collected {len(structures_to_add)} structures (complete)")
    print("")
    
    # Sort structures by space group (descending - high symmetry first)
    print("Sorting structures by space group (descending)...")
    structures_to_add.sort(key=lambda x: x['space_group'], reverse=True)
    print(f"  Sorted {len(structures_to_add)} structures")
    print("")
    
    # Write sorted structures to database
    print("Writing sorted structures to database...")
    added_count = 0
    for struct_data in structures_to_add:
        try:
            if struct_data['is_pyxtal']:
                # PyXtal-symmetrized structure
                db_out.add_xtal(struct_data['xtal'], struct_data['dicts'])
            else:
                # ASE-only structure (P1)
                db_out.db.write(
                    struct_data['atoms'],
                    structure_id=struct_data['struct_id'],
                    composition=struct_data['composition'],
                    e0025=struct_data['e0025'],
                    e05=struct_data['e05'],
                    e10=struct_data['e10'],
                    band0=struct_data['band0'],
                    band1=struct_data['band1'],
                    e_above_hull=struct_data['e_above_hull'],
                    is_electride=struct_data['is_electride'],
                    symmetrized=False,
                    space_group_number=1
                )
            added_count += 1
            
            if added_count % 100 == 0:
                print(f"  Written {added_count}/{len(structures_to_add)} structures...")
        except Exception as e:
            print(f"  Error writing structure: {e}")
            failed_count += 1
    
    print(f"  Written {added_count} structures to database")
    print("")
    
    return total_count, added_count, duplicate_count, filtered_count, failed_count


def export_to_csv(db_path, output_csv):
    """
    Export database to CSV file.
    
    Args:
        db_path: Path to PyXtal database
        output_csv: Output CSV file path
    """
    print(f"Exporting to CSV: {output_csv}")
    
    db = database_topology(db_path)
    
    results = []
    failed_retrieval = 0
    
    for row in db.db.select():
        try:
            # Get structure_id directly
            struct_id = row.structure_id
            
            # Get atoms object
            try:
                atoms = row.toatoms()
            except:
                print(f"  Warning: Skipping {struct_id} (could not get atoms)")
                failed_retrieval += 1
                continue
            
            # Get metadata directly from row attributes
            composition = getattr(row, 'composition', 'unknown')
            e0025 = float(getattr(row, 'e0025', 0))
            e05 = float(getattr(row, 'e05', 0))
            e10 = float(getattr(row, 'e10', 0))
            band0 = float(getattr(row, 'band0', 0))
            band1 = float(getattr(row, 'band1', 0))
            e_above_hull = float(getattr(row, 'e_above_hull', float('inf')))
            
            # Get space group - check if structure was symmetrized
            was_symmetrized = getattr(row, 'symmetrized', True)
            
            if was_symmetrized:
                # Try to get PyXtal structure for space group
                xtal = None
                for tol in [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]:
                    try:
                        xtal_temp = pyxtal()
                        xtal_temp.from_seed(atoms, tol=tol)
                        if xtal_temp.valid:
                            xtal = xtal_temp
                            break
                    except:
                        continue
                
                if xtal is None:
                    spacegroup = getattr(row, 'space_group_number', 1)
                else:
                    spacegroup = xtal.group.number
            else:
                # Structure was saved without symmetrization
                spacegroup = getattr(row, 'space_group_number', 1)
            
            result = {
                'formula': struct_id,
                'composition': composition,
                'e0025': e0025,
                'e05': e05,
                'e10': e10,
                'band0': band0,
                'band1': band1,
                'spacegroup': spacegroup,
                'e_above_hull': e_above_hull,
            }
            
            results.append(result)
            
        except Exception as e:
            print(f"  Error exporting row: {e}")
            failed_retrieval += 1
            continue
    
    # Create DataFrame and sort by space group (descending)
    df = pd.DataFrame(results)
    df = df.sort_values('spacegroup', ascending=False)
    
    # Write CSV
    df.to_csv(output_csv, index=False)
    
    print(f"  Exported {len(df)} structures to CSV")
    if failed_retrieval > 0:
        print(f"  Failed to retrieve: {failed_retrieval}")
    print("")


def main():
    parser = argparse.ArgumentParser(
        description="Filter electride database with stricter interstitial volume criteria",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default filtering (20 Å³ thresholds)
  python3 filter_comb_db.py --input electride_data.db --csv electride_analysis.csv
  
  # Custom thresholds
  python3 filter_comb_db.py --input electride_data.db --csv electride_analysis.csv \\
      --min-energy 25 --min-band 25
  
  # Custom output names
  python3 filter_comb_db.py --input electride_data.db --csv electride_analysis.csv \\
      --output electride_candidates.db --output-csv electride_candidates.csv

Filtering Criteria:
  - max(e0025, e05, e10) >= min-energy (default: 20 Å³)
  - max(band0, band1) >= min-band (default: 20 Å³)
  - Duplicate structures removed
  - Symmetrized with PyXtal (adaptive tolerance)

Output:
  - electride_candidates.db: PyXtal database with filtered structures
  - electride_candidates.csv: CSV table for quick viewing
  
The database contains:
  - Symmetrized crystal structures
  - Space group information
  - Interstitial volumes (e0025, e05, e10, band0, band1)
  - Thermodynamic stability (e_above_hull)
  - Sorted by space group number (descending - high symmetry first)
  
Query examples:
  # View all candidates
  ase db electride_candidates.db
  
  # Sort by space group (high symmetry first)
  ase db electride_candidates.db formula,space_group_number,e_above_hull -s space_group_number-
  
  # Filter by stability
  ase db electride_candidates.db e_above_hull<0.05
"""
    )
    
    parser.add_argument(
        '--input',
        type=str,
        default='electride_data.db',
        help="Input PyXtal database (default: electride_data.db)"
    )
    parser.add_argument(
        '--csv',
        type=str,
        default='electride_analysis.csv',
        help="Input CSV with interstitial volumes (default: electride_analysis.csv)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='electride_candidates.db',
        help="Output PyXtal database (default: electride_candidates.db)"
    )
    parser.add_argument(
        '--output-csv',
        type=str,
        default='electride_candidates.csv',
        help="Output CSV file (default: electride_candidates.csv)"
    )
    parser.add_argument(
        '--min-energy',
        type=float,
        default=20.0,
        help="Minimum for max(e0025, e05, e10) in Å³ (default: 20)"
    )
    parser.add_argument(
        '--min-band',
        type=float,
        default=20.0,
        help="Minimum for max(band0, band1) in Å³ (default: 20)"
    )
    
    args = parser.parse_args()
    
    # Check input files exist
    if not Path(args.input).exists():
        print(f"ERROR: Input database not found: {args.input}")
        sys.exit(1)
    
    if not Path(args.csv).exists():
        print(f"ERROR: Input CSV not found: {args.csv}")
        sys.exit(1)
    
    # Filter and create database
    total, added, duplicates, filtered, failed = filter_and_create_db(
        args.input,
        args.csv,
        args.output,
        min_energy_volume=args.min_energy,
        min_band_volume=args.min_band
    )
    
    # Export to CSV
    if added > 0:
        export_to_csv(args.output, args.output_csv)
    
    # Print summary
    print("="*70)
    print("Filtering Summary")
    print("="*70)
    print(f"Input structures: {total}")
    print(f"Filtered out (below thresholds): {filtered}")
    print(f"Duplicates removed: {duplicates}")
    print(f"Added to output database: {added}")
    if failed > 0:
        print(f"Failed to process: {failed}")
    print("")
    print(f"Output database: {args.output} ({added} structures)")
    print(f"Output CSV: {args.output_csv}")
    print("="*70)
    print("")
    
    if added > 0:
        print("Query examples:")
        print(f"  # View all candidates")
        print(f"  ase db {args.output}")
        print(f"")
        print(f"  # Sort by space group")
        print(f"  ase db {args.output} formula,space_group_number,e_above_hull -s space_group_number-")
        print(f"")
        print(f"  # Filter by stability (< 50 meV/atom)")
        print(f"  ase db {args.output} 'e_above_hull<0.05'")
        print("")


if __name__ == '__main__':
    main()

