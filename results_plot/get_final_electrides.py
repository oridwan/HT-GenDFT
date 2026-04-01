#!/usr/bin/env python3
"""
Add additional properties to stable electrides and create final output files.

This script:
1. Loads stable_electrides.csv and stable_electrides.db from results directory
2. Removes duplicate structures using PyXtal structure comparison (same_group=True, d_tol=0.05)
3. Adds full_dft_e_hull considering all new structures + MP reference phases
4. Adds N_excess (excess valence electrons based on composition)
5. Adds band_gap from electronic workflow workflow.json files
6. Adds electride column (default True, can be updated later)
7. Outputs final_electrides.csv and final_electrides.db in the same directory

Usage:
    python3 get_final_electrides.py --results-dir Bin-Ele-HT/
    python3 get_final_electrides.py --results-dir Ter-Ele-HT/
"""

import json
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
from tabulate import tabulate

from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.structure_matcher import StructureMatcher

try:
    from pyxtal import pyxtal
    from pyxtal.db import database_topology
    PYXTAL_AVAILABLE = True
except ImportError:
    PYXTAL_AVAILABLE = False
    print("WARNING: PyXtal not available. Database operations will be skipped.")


VALENCE_ELECTRONS = {
    'H': 1, 'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1,
    'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
    'B': 3, 'Al': 3, 'Ga': 3, 'In': 3, 'Tl': 3,
    'Sc': 3, 'Y': 3,
    # 'La': 3, 'Ce': 3, 'Pr': 3, 'Nd': 3, 'Pm': 3, 'Sm': 3, 'Eu': 3,
    # 'Gd': 3, 'Tb': 3, 'Dy': 3, 'Ho': 3, 'Er': 3, 'Tm': 3, 'Yb': 3, 'Lu': 3,
    'C': 4, 'Si': 4, 'Ge': 4, 'Sn': 4, 'Pb': 4,
    'N': 3, 'P': 3, 'As': 3, 'Sb': 3, 'Bi': 3,
    'O': 2, 'S': 2, 'Se': 2, 'Te': 2, 'Po': 2,
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'At': 1,
}


def detect_system_type(df: pd.DataFrame, workflow: Dict) -> str:
    """
    Detect if the system is binary or ternary by counting elements in compositions.
    
    Returns:
        'binary', 'ternary', or 'unknown'
    """
    if len(df) > 0 and 'composition' in df.columns:
        first_comp = Composition(df.iloc[0]['composition'])
        n_elements = len(first_comp.elements)
        if n_elements == 2:
            return 'binary'
        elif n_elements == 3:
            return 'ternary'
    
    if len(workflow) > 0:
        first_struct = next(iter(workflow.values()))
        if 'chemsys' in first_struct:
            n_elements = len(first_struct['chemsys'].split('-'))
            if n_elements == 2:
                return 'binary'
            elif n_elements == 3:
                return 'ternary'
    
    return 'unknown'


def load_csv_data(csv_path: Path) -> pd.DataFrame:
    """Load CSV file into DataFrame."""
    if not csv_path.exists():
        print(f"WARNING: {csv_path} not found, skipping.")
        return pd.DataFrame()
    
    return pd.read_csv(csv_path)


def load_workflow_json(json_path: Path) -> Dict:
    """Load workflow.json file."""
    if not json_path.exists():
        print(f"WARNING: {json_path} not found, band_gap will be unavailable.")
        return {}
    
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    return data.get('structures', {})


def load_mp_phases(json_path: Path) -> List[Dict]:
    """Load MP reference phases."""
    if not json_path.exists():
        print(f"WARNING: {json_path} not found.")
        return []
    
    with open(json_path, 'r') as f:
        return json.load(f)


def calculate_excess_electrons(composition: str) -> float:
    """
    Calculate excess valence electrons in a composition.
    
    - Sum valence electrons from electropositive elements (Groups I, II, III)
    - Subtract valence electrons needed by electronegative elements (Groups V, VI, VII)
    
    Returns:
        Excess electrons (should be 0 < N_excess <= 4 for electride candidates)
    """
    comp = Composition(composition)
    
    excess = 0.0
    for element, amount in comp.items():
        elem_symbol = str(element)
        if elem_symbol in VALENCE_ELECTRONS:
            val = VALENCE_ELECTRONS[elem_symbol]
            
            # Electronegative elements (N, P, O, S, F, Cl, etc.) subtract electrons
            if elem_symbol in ['N', 'P', 'As', 'Sb', 'Bi', 'O', 'S', 'Se', 'Te', 'Po', 'F', 'Cl', 'Br', 'I', 'At', 'H']:
                excess -= abs(val) * amount
            else:
                # Electropositive elements contribute electrons
                excess += val * amount
    
    return excess


def deduplicate_structures(
    db_path: Path,
    df: pd.DataFrame
) -> Tuple[List[Dict], List[int]]:
    """
    Deduplicate structures within a single database.
    
    Uses PyXtal for structure loading with progressive tolerances and Pymatgen's StructureMatcher
    for pairwise comparison (same_group constraint via space group check, d_tol=0.05 via ltol/stol).
    
    Returns:
        unique_structures: List of unique structure data dicts
        keep_indices: Indices to keep from df
    """
    if not PYXTAL_AVAILABLE:
        print("WARNING: PyXtal not available. Skipping deduplication, keeping all structures.")
        return [], list(range(len(df)))
    
    print("="*70)
    print("Running to remove duplicates . . .")
    print("="*70)
    
    unique_structures = []
    keep_indices = []
    processed_xtals = []
    
    if not db_path.exists() or len(df) == 0:
        print(f"\nDatabase not found or empty, keeping all CSV entries")
        return [], list(range(len(df)))
    
    print(f"\nLoading database: {db_path}")
    db = database_topology(str(db_path))
    
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    structure_matcher = StructureMatcher(ltol=0.05, stol=0.05, angle_tol=2)
    
    for idx, row in df.iterrows():
        struct_id = row['formula']
        
        try:
            db_row = db.db.get(structure_id=struct_id)
            atoms = db_row.toatoms()
            
            xtal = None
            for tol in tolerances:
                try:
                    xtal_attempt = pyxtal()
                    xtal_attempt.from_seed(atoms, tol=tol)
                    
                    if xtal_attempt.valid:
                        if len(xtal_attempt.check_short_distances(r=0.5)) == 0:
                            xtal = xtal_attempt
                            break
                except:
                    continue
            
            if xtal is None or not xtal.valid:
                print(f"  {struct_id}: Failed to load with all tolerances, keeping in CSV")
                keep_indices.append(idx)
                continue
            
            is_unique = True
            pmg_struct = xtal.to_pymatgen()
            current_comp = pmg_struct.composition
            
            for existing_xtal, existing_comp in processed_xtals:
                try:
                    if current_comp != existing_comp:
                        continue
                    if xtal.group.number != existing_xtal.group.number:
                        continue
                    existing_pmg_struct = existing_xtal.to_pymatgen()
                    
                    if structure_matcher.fit(pmg_struct, existing_pmg_struct):
                        is_unique = False
                        print(f"  {struct_id}: Duplicate (matches existing structure), removing")
                        break
                        
                except Exception as e:
                    print(f"  {struct_id}: Warning in duplicate check: {e}")
                    continue
            
            if is_unique:
                keep_indices.append(idx)
                processed_xtals.append((xtal, current_comp))
                unique_structures.append({
                    'structure_id': struct_id,
                    'xtal': xtal,
                    'atoms': atoms
                })
        
        except Exception as e:
            print(f"  {struct_id}: Error loading from database: {e}, keeping in CSV")
            keep_indices.append(idx)
    
    print(f"\nDeduplication summary:")
    removed = len(df) - len(keep_indices)
    print(f"  Original: {len(df)}")
    print(f"  Unique: {len(keep_indices)}")
    print(f"  Removed: {removed}")
    print("="*70 + "\n")
    
    return unique_structures, keep_indices


def calculate_full_dft_e_hull(
    df: pd.DataFrame,
    mp_phases: List[Dict]
) -> pd.Series:
    """
    Calculate full DFT energy above hull considering all new structures + MP reference phases.
    
    Optimized by grouping structures by chemical system and building separate phase diagrams.
    
    Returns:
        Series of full_dft_e_hull values (eV/atom)
    """
    print("="*70)
    print("Calculating full DFT energy above hull")
    print("="*70)
    
    # Group structures by chemical system
    structures_by_chemsys = {}
    for idx, row in df.iterrows():
        comp = Composition(row['composition'])
        elements = sorted([str(e) for e in comp.elements])
        chemsys = '-'.join(elements)
        
        if chemsys not in structures_by_chemsys:
            structures_by_chemsys[chemsys] = []
        structures_by_chemsys[chemsys].append((idx, row, comp))
    
    print(f"Found {len(structures_by_chemsys)} unique chemical systems")
    print(f"Total structures: {len(df)}")
    print()
    
    # Initialize results
    full_e_hull_values = pd.Series(index=df.index, dtype=float)
    total_stable = 0
    
    # Process each chemical system separately
    for chemsys_idx, (chemsys, struct_list) in enumerate(structures_by_chemsys.items(), 1):
        elements = chemsys.split('-')
        print(f"[{chemsys_idx}/{len(structures_by_chemsys)}] Processing {chemsys} ({len(struct_list)} structures)...")
        
        # Filter MP phases to only relevant ones for this chemsys
        relevant_mp_phases = []
        for phase in mp_phases:
            phase_comp = Composition(phase['composition'])
            phase_elements = set(str(e) for e in phase_comp.elements)
            
            # Include if all phase elements are in our chemsys
            if phase_elements.issubset(set(elements)):
                relevant_mp_phases.append(phase)
        
        print(f"  Relevant MP phases: {len(relevant_mp_phases)}")
        
        # Create entries for this chemsys
        chemsys_entries = []
        struct_id_to_entry = {}
        
        # Add MP reference phases
        for phase in relevant_mp_phases:
            comp = Composition(phase['composition'])
            total_energy = phase['energy']
            entry_id = phase.get('entry_id', phase.get('mp_id', 'unknown'))
            
            entry = ComputedEntry(
                composition=comp,
                energy=total_energy,
                entry_id=entry_id
            )
            chemsys_entries.append(entry)
        
        # Add structures from this chemsys
        for idx, row, comp in struct_list:
            energy_per_atom = row['vasp_energy_per_atom']
            n_atoms = comp.num_atoms
            total_energy = energy_per_atom * n_atoms
            struct_id = row['formula']
            
            entry = ComputedEntry(
                composition=comp,
                energy=total_energy,
                entry_id=struct_id
            )
            chemsys_entries.append(entry)
            struct_id_to_entry[struct_id] = entry
        
        # Build phase diagram for this chemsys
        print(f"  Building phase diagram ({len(chemsys_entries)} total entries)...")
        pd_chemsys = PhaseDiagram(chemsys_entries)
        
        # Calculate e_hull for each structure in this chemsys
        stable_count = 0
        for idx, row, comp in struct_list:
            struct_id = row['formula']
            entry = struct_id_to_entry[struct_id]
            
            e_hull = pd_chemsys.get_e_above_hull(entry)
            full_e_hull_values[idx] = e_hull
            
            if e_hull < 1e-6:
                stable_count += 1
        
        print(f"  Stable in {chemsys}: {stable_count}")
        total_stable += stable_count
    
    print()
    print(f"Total stable structures (full_dft_e_hull = 0.0 eV/atom): {total_stable}")
    print("="*70 + "\n")
    
    return full_e_hull_values


def add_band_gap_from_workflow(
    df: pd.DataFrame,
    workflow: Dict
) -> pd.Series:
    """
    Add band_gap column from workflow.json file.
    
    Returns:
        Series of band_gap values (eV)
    """
    print("="*70)
    print("Adding band gap from workflow.json")
    print("="*70)
    
    band_gaps = []
    missing_count = 0
    
    for idx, row in df.iterrows():
        struct_id = row['formula']
        
        if struct_id in workflow:
            band_gap = workflow[struct_id].get('band_gap', None)
            if band_gap is None:
                band_gaps.append(np.nan)
                missing_count += 1
            else:
                band_gaps.append(band_gap)
        else:
            band_gaps.append(np.nan)
            missing_count += 1
    
    print(f"Band gap data:")
    print(f"  Found: {len(df) - missing_count}")
    print(f"  Missing: {missing_count}")
    print("="*70 + "\n")
    
    return pd.Series(band_gaps, index=df.index)


def save_to_database(
    df: pd.DataFrame,
    unique_structures: List[Dict],
    output_db_path: Path
):
    """
    Save deduplicated structures to PyXtal database.
    
    Args:
        df: DataFrame with all columns
        unique_structures: List of unique structure data from deduplication
        output_db_path: Output database path
    """
    if not PYXTAL_AVAILABLE:
        print("WARNING: PyXtal not available. Skipping database creation.")
        return
    
    print("="*70)
    print(f"Creating output database: {output_db_path}")
    print("="*70)
    
    # Create new database
    db = database_topology(str(output_db_path))
    
    # Create mapping from structure_id to row data
    struct_id_to_row = {row['formula']: row for _, row in df.iterrows()}
    
    saved_count = 0
    failed_count = 0
    
    for struct_data in unique_structures:
        struct_id = struct_data['structure_id']
        xtal = struct_data['xtal']
        
        if struct_id not in struct_id_to_row:
            print(f"  {struct_id}: WARNING - not found in merged DataFrame")
            failed_count += 1
            continue
        
        row = struct_id_to_row[struct_id]
        
        kvp_data = {
            'structure_id': struct_id,
            'composition': row['composition'],
            'space_group_number': int(row['spacegroup']),
            'e0025': float(row['e0025']),
            'e05': float(row['e05']),
            'e10': float(row['e10']),
            'band0': float(row['band0']),
            'band1': float(row['band1']),
            'mattersim_energy_per_atom': float(row['mattersim_energy_per_atom']),
            'mattersim_e_hull': float(row['mattersim_e_hull']),
            'vasp_energy_per_atom': float(row['vasp_energy_per_atom']),
            'dft_e_hull': float(row['dft_e_hull']),
            'full_dft_e_hull': float(row['full_dft_e_hull']),
            'N_excess': float(row['N_excess']),
            'electride': bool(row['electride']),
            'symmetrized': True
        }
        
        # Only add band_gap if it's not NaN (for SC_FAILED cases)
        if pd.notna(row['band_gap']):
            kvp_data['band_gap'] = float(row['band_gap'])
        
        try:
            db.add_xtal(xtal, kvp=kvp_data)
            saved_count += 1
            
            if saved_count % 100 == 0:
                print(f"  Saved {saved_count} structures...")
        except Exception as e:
            print(f"  {struct_id}: Failed to save: {e}")
            failed_count += 1
    
    # Commit database
    if hasattr(db, 'db') and hasattr(db.db, 'commit'):
        db.db.commit()
    
    print(f"\nDatabase summary:")
    print(f"  Saved: {saved_count}")
    print(f"  Failed: {failed_count}")
    print("="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Add additional properties to stable electrides and create final output files"
    )
    
    parser.add_argument('--results-dir', type=str, required=True,
                        help='Results directory containing stable_electrides.csv/db, workflow.json, and mp_vaspdft.json')
    
    args = parser.parse_args()
    
    # Convert to Path objects
    results_dir = Path(args.results_dir)
    
    print("\n" + "="*70)
    print("PROCESSING STABLE ELECTRIDES")
    print("="*70 + "\n")
    
    # Validate results directory
    if not results_dir.exists():
        print(f"ERROR: Results directory not found: {results_dir}")
        sys.exit(1)
    
    # Define file paths
    csv_path = results_dir / 'stable_electrides.csv'
    db_path = results_dir / 'stable_electrides.db'
    workflow_path = results_dir / 'workflow.json'
    mp_phases_path = results_dir / 'mp_vaspdft.json'
    output_csv_path = results_dir / 'final_electrides.csv'
    output_db_path = results_dir / 'final_electrides.db'
    
    # Load data
    print(f"Results directory: {results_dir}")
    df = load_csv_data(csv_path)
    workflow = load_workflow_json(workflow_path)
    mp_phases = load_mp_phases(mp_phases_path)
    
    system_type = detect_system_type(df, workflow)
    
    print(f"  System type: {system_type}")
    print(f"  Structures: {len(df)}")
    print(f"  Workflow entries: {len(workflow)}")
    print(f"  MP phases: {len(mp_phases)}")
    print()
    
    # Deduplicate structures
    unique_structures, keep_indices = deduplicate_structures(db_path, df)
    
    # Filter DataFrame to keep only unique structures
    df = df.iloc[keep_indices].copy()
    
    print(f"Unique structures: {len(df)}")
    print()
    
    # Add N_excess column
    print("Calculating excess electrons (N_excess)...")
    df['N_excess'] = df['composition'].apply(calculate_excess_electrons)
    print(f"  N_excess range: {df['N_excess'].min():.2f} to {df['N_excess'].max():.2f}")
    print()
    
    # Add band_gap column
    df['band_gap'] = add_band_gap_from_workflow(df, workflow)
    
    # Add electride column (default True)
    df['electride'] = True
    
    # Calculate full DFT e_hull
    df['full_dft_e_hull'] = calculate_full_dft_e_hull(df, mp_phases)
    
    # Sort by full_dft_e_hull (ascending), then by spacegroup (descending)
    df = df.sort_values(
        by=['full_dft_e_hull', 'spacegroup'],
        ascending=[True, False]
    ).reset_index(drop=True)
    
    # Save CSV
    print(f"Saving CSV: {output_csv_path}")
    df.to_csv(output_csv_path, index=False)
    print()
    
    # Save database
    if PYXTAL_AVAILABLE and len(unique_structures) > 0:
        save_to_database(df, unique_structures, output_db_path)
    
    # Print summary
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"System type: {system_type}")
    print(f"Total unique structures: {len(df)}")
    print()
    print(f"Stable structures (full_dft_e_hull = 0.0 eV/atom): {(df['full_dft_e_hull'] < 1e-6).sum()}")
    print(f"Metastable (full_dft_e_hull <= 0.05 eV/atom): {(df['full_dft_e_hull'] <= 0.05).sum()}")
    print()
    print(f"N_excess statistics:")
    print(f"  Mean: {df['N_excess'].mean():.2f}")
    print(f"  Range: [{df['N_excess'].min():.2f}, {df['N_excess'].max():.2f}]")
    print()
    print(f"Band gap statistics:")
    valid_gaps = df['band_gap'].dropna()
    if len(valid_gaps) > 0:
        print(f"  Available: {len(valid_gaps)}/{len(df)}")
        print(f"  Mean: {valid_gaps.mean():.3f} eV")
        print(f"  Range: [{valid_gaps.min():.3f}, {valid_gaps.max():.3f}] eV")
    else:
        print(f"  No band gap data available")
    print()
    print(f"Output files:")
    print(f"  CSV: {output_csv_path}")
    if PYXTAL_AVAILABLE:
        print(f"  Database: {output_db_path}")
    print("="*70 + "\n")
    
    # Show top structures
    if len(df) > 0:
        print("Top 20 structures (sorted by full_dft_e_hull):")
        display_cols = ['formula', 'composition', 'spacegroup', 'full_dft_e_hull', 
                       'dft_e_hull', 'N_excess', 'band_gap', 'electride']
        available_cols = [col for col in display_cols if col in df.columns]
        print(tabulate(df[available_cols].head(20), headers='keys', 
                      tablefmt='psql', showindex=False, floatfmt='.4f'))
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
