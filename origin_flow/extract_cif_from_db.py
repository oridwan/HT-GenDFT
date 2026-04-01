#!/usr/bin/env python3
"""
Extract specific structures from PyXtal/ASE database as CIF files.

Usage:
    python3 extract_cif_from_db.py --db electride_data.db Ba8S7_s024 Ba8S7_s029
    python3 extract_cif_from_db.py --db electride_candidates.db --output-dir cifs/ Y10N9_s001
"""

import sys
import argparse
from pathlib import Path
from ase.io import write as ase_write
from ase.geometry import get_distances

try:
    from pyxtal.db import database_topology
except ImportError as e:
    print(f"ERROR: Required package not available: {e}")
    print("Install with: conda install -c conda-forge pyxtal")
    sys.exit(1)


def extract_structure_info(row):
    """Extract structure information from database row."""
    struct_id = row.structure_id
    
    try:
        atoms = row.toatoms()
    except Exception as e:
        return None, f"Could not get atoms: {e}"
    
    # Get metadata
    composition = getattr(row, 'composition', 'unknown')
    spacegroup = getattr(row, 'space_group_number', None)
    e_above_hull = getattr(row, 'e_above_hull', None)
    symmetrized = getattr(row, 'symmetrized', None)
    
    # Get interstitial volumes if available
    e0025 = getattr(row, 'e0025', None)
    e05 = getattr(row, 'e05', None)
    e10 = getattr(row, 'e10', None)
    band0 = getattr(row, 'band0', None)
    band1 = getattr(row, 'band1', None)
    
    # Calculate minimum interatomic distance
    min_dist = None
    closest_pair = None
    if len(atoms) > 1:
        try:
            # Get all pairwise distances
            positions = atoms.get_positions()
            cell = atoms.get_cell()
            pbc = atoms.get_pbc()
            
            distances = []
            pairs = []
            for i in range(len(atoms)):
                for j in range(i+1, len(atoms)):
                    # Get distance considering PBC
                    d_vec = positions[j] - positions[i]
                    if any(pbc):
                        # Apply minimum image convention
                        for k in range(3):
                            if pbc[k]:
                                d_vec[k] -= cell[k, k] * round(d_vec[k] / cell[k, k])
                    d = (d_vec**2).sum()**0.5
                    distances.append(d)
                    pairs.append((i, j, atoms.symbols[i], atoms.symbols[j]))
            
            if distances:
                min_idx = distances.index(min(distances))
                min_dist = distances[min_idx]
                i, j, sym_i, sym_j = pairs[min_idx]
                closest_pair = f"{sym_i}{i} - {sym_j}{j}"
        except:
            pass
    
    info = {
        'structure_id': struct_id,
        'composition': composition,
        'natoms': len(atoms),
        'spacegroup': spacegroup,
        'symmetrized': symmetrized,
        'e_above_hull': e_above_hull,
        'min_distance': min_dist,
        'closest_pair': closest_pair,
        'e0025': e0025,
        'e05': e05,
        'e10': e10,
        'band0': band0,
        'band1': band1,
        'atoms': atoms
    }
    
    return info, None


def main():
    parser = argparse.ArgumentParser(
        description="Extract specific structures from database as CIF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract specific structures
  python3 extract_cif_from_db.py --db electride_data.db Ba8S7_s024 Ba8S7_s029
  
  # Extract to custom directory
  python3 extract_cif_from_db.py --db electride_candidates.db --output-dir cifs/ Y10N9_s001
  
  # Extract all structures (warning: may create many files)
  python3 extract_cif_from_db.py --db electride_candidates.db --all
"""
    )
    
    parser.add_argument(
        '--db',
        type=str,
        required=True,
        help="Path to PyXtal/ASE database"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='.',
        help="Output directory for CIF files (default: current directory)"
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help="Extract all structures from database"
    )
    parser.add_argument(
        'structure_ids',
        nargs='*',
        help="Structure IDs to extract"
    )
    
    args = parser.parse_args()
    
    # Check database exists
    if not Path(args.db).exists():
        print(f"ERROR: Database not found: {args.db}")
        sys.exit(1)
    
    # Check that either --all or structure_ids are provided
    if not args.all and not args.structure_ids:
        print("ERROR: Must provide either --all or structure IDs")
        parser.print_help()
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("Extract Structures from Database")
    print("="*70)
    print(f"Database: {args.db}")
    print(f"Output directory: {output_dir}")
    print("")
    
    # Load database
    print("Loading database...")
    db = database_topology(args.db)
    total_count = db.db.count()
    print(f"  Database contains {total_count} structures")
    print("")
    
    # Determine which structures to extract
    if args.all:
        print("Extracting all structures...")
        target_ids = None  # Will process all
    else:
        target_ids = set(args.structure_ids)
        print(f"Extracting {len(target_ids)} specified structures...")
        print("")
    
    # Extract structures
    extracted_count = 0
    not_found = []
    failed = []
    
    for row in db.db.select():
        try:
            struct_id = row.structure_id
        except AttributeError:
            continue
        
        # Skip if not in target list (if targets specified)
        if target_ids is not None and struct_id not in target_ids:
            continue
        
        # Extract structure info
        info, error = extract_structure_info(row)
        
        if error:
            print(f"ERROR: {struct_id} - {error}")
            failed.append((struct_id, error))
            continue
        
        # Print structure info
        print(f"Structure: {info['structure_id']}")
        print(f"  Composition: {info['composition']}")
        print(f"  Number of atoms: {info['natoms']}")
        if info['spacegroup']:
            print(f"  Space group: {info['spacegroup']}")
        if info['symmetrized'] is not None:
            print(f"  Symmetrized: {info['symmetrized']}")
        if info['e_above_hull'] is not None:
            print(f"  E above hull: {info['e_above_hull']:.4f} eV/atom")
        if info['min_distance']:
            print(f"  Min interatomic distance: {info['min_distance']:.3f} Å")
            if info['min_distance'] < 0.5:
                print(f"    WARNING: Very short distance! Closest pair: {info['closest_pair']}")
            else:
                print(f"    Closest pair: {info['closest_pair']}")
        if info['e0025'] is not None:
            print(f"  Interstitial volumes:")
            print(f"    e0025={info['e0025']:.2f}, e05={info['e05']:.2f}, e10={info['e10']:.2f}")
            print(f"    band0={info['band0']:.2f}, band1={info['band1']:.2f}")
        
        # Write CIF file
        cif_path = output_dir / f"{info['structure_id']}.cif"
        try:
            ase_write(str(cif_path), info['atoms'], format='cif')
            print(f"  Saved: {cif_path}")
            extracted_count += 1
        except Exception as e:
            print(f"  ERROR: Could not write CIF: {e}")
            failed.append((struct_id, f"Could not write CIF: {e}"))
        
        print("")
    
    # Check for structures not found
    if target_ids is not None:
        found_ids = set()
        for row in db.db.select():
            try:
                found_ids.add(row.structure_id)
            except AttributeError:
                pass
        not_found = list(target_ids - found_ids)
    
    # Print summary
    print("="*70)
    print("Summary")
    print("="*70)
    print(f"Extracted: {extracted_count}")
    if not_found:
        print(f"Not found: {len(not_found)}")
        for struct_id in not_found:
            print(f"  - {struct_id}")
    if failed:
        print(f"Failed: {len(failed)}")
        for struct_id, error in failed:
            print(f"  - {struct_id}: {error}")
    print("="*70)


if __name__ == '__main__':
    main()

