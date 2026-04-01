#!/usr/bin/env python3
"""
Extract stable electride candidates confirmed by both MatterSim and VASP-DFT.

This script:
1. Reads hull_comparison.json from refined workflow
2. Filters structures with mattersim_e_hull < threshold AND dft_e_hull < threshold
3. Loads CONTCAR from refined 3-step relaxations
4. Extracts spacegroup using PyXtal symmetrization
5. Loads electride volume data (e0025, e05, e10, band0, band1) from electride_analysis.csv
6. Creates stable_electrides.csv and stable_electrides.db

Output format matches electride_candidates.csv/db from analyze.py with extra columns:
- e0025, e05, e10, band0, band1 (interstitial volumes from electride analysis)
- mattersim_energy_per_atom, mattersim_e_hull
- vasp_energy_per_atom, dft_e_hull

Usage:
    python3 refined_flow/get_stable_ele_db.py  --refine-jobs REFINE_VASP_JOBS/  --electride-csv ./electride_candidates.csv  --threshold 0.005
"""

import os
import sys
import json
import argparse
import pandas as pd
from pathlib import Path
from tabulate import tabulate

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

try:
    from pyxtal import pyxtal
    from pyxtal.db import database_topology
    PYXTAL_AVAILABLE = True
except ImportError:
    PYXTAL_AVAILABLE = False
    print("WARNING: PyXtal not available. Database creation will be skipped.")


def load_hull_comparison(hull_file):
    """Load hull_comparison.json and return matched structures."""
    with open(hull_file, 'r') as f:
        data = json.load(f)
    return data.get('matched_structures', [])


def load_mattersim_results(json_file):
    """
    Load mattersim_stability_results.json for mattersim_energy_per_atom.
    
    Returns:
        dict: {structure_id: mattersim_energy_per_atom}
    """
    if not os.path.exists(json_file):
        print(f"WARNING: MatterSim results file not found: {json_file}")
        return {}
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        energy_map = {}
        for result in data.get('results', []):
            struct_id = result.get('structure_id')
            energy = result.get('mattersim_energy_per_atom')
            if struct_id and energy is not None:
                energy_map[struct_id] = energy
        
        print(f"Loaded MatterSim energies for {len(energy_map)} structures")
        return energy_map
    except Exception as e:
        print(f"WARNING: Could not load MatterSim results: {e}")
        return {}


def load_dft_results(json_file):
    """
    Load dft_stability_results.json for vasp_energy_per_atom.
    
    Returns:
        dict: {structure_id: vasp_energy_per_atom}
    """
    if not os.path.exists(json_file):
        print(f"WARNING: DFT results file not found: {json_file}")
        return {}
    
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        energy_map = {}
        for result in data.get('results', []):
            struct_id = result.get('structure_id')
            energy = result.get('vasp_energy_per_atom')
            if struct_id and energy is not None:
                energy_map[struct_id] = energy
        
        print(f"Loaded VASP DFT energies for {len(energy_map)} structures")
        return energy_map
    except Exception as e:
        print(f"WARNING: Could not load DFT results: {e}")
        return {}


def load_electride_volumes(csv_path):
    """
    Load electride volume data from electride_analysis.csv.
    
    Args:
        csv_path: Path to electride_analysis.csv
    
    Returns:
        dict: {structure_id: {'e0025': float, 'e05': float, 'e10': float, 'band0': float, 'band1': float}}
    """
    if not os.path.exists(csv_path):
        print(f"WARNING: Electride analysis CSV not found: {csv_path}")
        return {}
    
    try:
        df = pd.read_csv(csv_path)
        volumes = {}
        
        for _, row in df.iterrows():
            struct_id = row['formula']
            volumes[struct_id] = {
                'e0025': float(row.get('e0025', 0)),
                'e05': float(row.get('e05', 0)),
                'e10': float(row.get('e10', 0)),
                'band0': float(row.get('band0', 0)),
                'band1': float(row.get('band1', 0))
            }
        
        print(f"Loaded electride volumes for {len(volumes)} structures")
        return volumes
    except Exception as e:
        print(f"WARNING: Could not load electride analysis CSV: {e}")
        return {}


def load_workflow_database(db_path):
    """Load workflow.json database."""
    with open(db_path, 'r') as f:
        return json.load(f)


def get_spacegroup_from_structure(structure):
    """
    Extract spacegroup from structure using PyXtal with progressive tolerances.
    
    Args:
        structure: Pymatgen Structure object
    
    Returns:
        int: Space group number (1 if symmetrization fails)
    """
    if not PYXTAL_AVAILABLE:
        return 1
    
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    
    for tol in tolerances:
        try:
            xtal = pyxtal()
            xtal.from_seed(structure, tol=tol)
            if not xtal.valid:
                continue
            if len(xtal.check_short_distances(r=0.5)) > 0:
                continue
            return xtal.group.number
        except Exception:
            continue
    
    # If all tolerances failed, return P1
    return 1


def save_to_database(struct_id, structure, composition, spacegroup, 
                     mattersim_e_hull, dft_e_hull, 
                     mattersim_energy, vasp_energy,
                     e0025, e05, e10, band0, band1, db):
    """
    Save structure to PyXtal database with adaptive tolerance.
    
    Args:
        struct_id: Structure identifier
        structure: Pymatgen Structure object
        composition: Composition string
        spacegroup: Space group number
        mattersim_e_hull: MatterSim energy above hull
        dft_e_hull: DFT energy above hull
        mattersim_energy: MatterSim energy per atom
        vasp_energy: VASP energy per atom
        e0025, e05, e10, band0, band1: Interstitial volumes from electride analysis
        db: PyXtal database_topology instance
    
    Returns:
        bool: True if saved successfully
    """
    if db is None or not PYXTAL_AVAILABLE:
        return False
    
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    
    # Common key-value pairs for database
    kvp_data = {
        'structure_id': struct_id,
        'composition': composition,
        'space_group_number': spacegroup,
        'e0025': e0025,
        'e05': e05,
        'e10': e10,
        'band0': band0,
        'band1': band1,
        'mattersim_energy_per_atom': mattersim_energy,
        'mattersim_e_hull': mattersim_e_hull,
        'vasp_energy_per_atom': vasp_energy,
        'dft_e_hull': dft_e_hull,
        'symmetrized': True
    }
    
    # Try progressive tolerances for symmetrization
    for tol in tolerances:
        try:
            xtal = pyxtal()
            xtal.from_seed(structure, tol=tol)
            if not xtal.valid:
                continue
            if len(xtal.check_short_distances(r=0.5)) > 0:
                continue
            
            db.add_xtal(xtal, kvp=kvp_data)
            return True
        except Exception:
            continue
    
    # If all tolerances failed, save without symmetrization
    try:
        adaptor = AseAtomsAdaptor()
        atoms = adaptor.get_atoms(structure)
        
        kvp_data['symmetrized'] = False
        db.db.write(atoms, **kvp_data)
        return True
    except Exception as e:
        print(f"    Warning: Could not save {struct_id} to database: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Extract stable electride candidates confirmed by MatterSim and DFT"
    )
    parser.add_argument(
        '--refine-jobs',
        type=str,
        default='./REFINE_VASP_JOBS',
        help="Refined VASP jobs directory (default: ./REFINE_VASP_JOBS)"
    )
    parser.add_argument(
        '--hull-comparison',
        type=str,
        default='hull_comparison.json',
        help="Hull comparison JSON file (default: hull_comparison.json in refine-jobs)"
    )
    parser.add_argument(
        '--workflow-db',
        type=str,
        default='workflow.json',
        help="Workflow database (default: workflow.json in refine-jobs)"
    )
    parser.add_argument(
        '--electride-csv',
        type=str,
        default='./electride_candidates.csv',
        help="Electride candidates CSV with volume data (default: ./electride_candidates.csv relative to refine-jobs)"
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.005,
        help="Energy above hull threshold in eV/atom (default: 0.005)"
    )
    parser.add_argument(
        '--output-csv',
        type=str,
        default='stable_electrides.csv',
        help="Output CSV file (default: stable_electrides.csv)"
    )
    parser.add_argument(
        '--output-db',
        type=str,
        default='stable_electrides.db',
        help="Output PyXtal database file (default: stable_electrides.db)"
    )
    
    args = parser.parse_args()
    
    # Resolve paths
    refine_jobs = Path(args.refine_jobs).expanduser()
    
    hull_file = Path(args.hull_comparison).expanduser()
    if not hull_file.is_absolute():
        hull_file = refine_jobs / args.hull_comparison
    
    workflow_db_path = Path(args.workflow_db).expanduser()
    if not workflow_db_path.is_absolute():
        workflow_db_path = refine_jobs / args.workflow_db
    
    electride_csv_path = Path(args.electride_csv).expanduser()
    if not electride_csv_path.is_absolute():
        electride_csv_path = refine_jobs / args.electride_csv
    
    output_csv = refine_jobs / args.output_csv
    output_db = refine_jobs / args.output_db
    
    # Additional result files for energy_per_atom values
    mattersim_results_file = refine_jobs / 'mattersim_stability_results.json'
    dft_results_file = refine_jobs / 'dft_stability_results.json'
    
    print("="*70)
    print("Extract Stable Electride Candidates")
    print("="*70)
    print(f"Refined VASP jobs: {refine_jobs}")
    print(f"Hull comparison: {hull_file}")
    print(f"MatterSim results: {mattersim_results_file}")
    print(f"DFT results: {dft_results_file}")
    print(f"Workflow database: {workflow_db_path}")
    print(f"Electride candidates: {electride_csv_path}")
    print(f"Stability threshold: {args.threshold} eV/atom")
    print(f"Output CSV: {output_csv}")
    print(f"Output database: {output_db}")
    print("="*70 + "\n")
    
    # Check input files
    if not hull_file.exists():
        print(f"ERROR: Hull comparison file not found: {hull_file}")
        return 1
    
    if not workflow_db_path.exists():
        print(f"ERROR: Workflow database not found: {workflow_db_path}")
        return 1
    
    # Load data
    print("Loading hull comparison data...")
    matched_structures = load_hull_comparison(hull_file)
    print(f"  Found {len(matched_structures)} matched structures\n")
    
    print("Loading energy_per_atom data from result files...")
    mattersim_energies = load_mattersim_results(mattersim_results_file)
    dft_energies = load_dft_results(dft_results_file)
    print()
    
    print("Loading workflow database...")
    workflow_db = load_workflow_database(workflow_db_path)
    print(f"  Found {len(workflow_db['structures'])} structures in workflow\n")
    
    print("Loading electride volume data...")
    electride_volumes = load_electride_volumes(electride_csv_path)
    print()
    
    # Filter stable structures
    print(f"Filtering structures with mattersim_e_hull < {args.threshold} AND dft_e_hull < {args.threshold}...")
    stable_structures = []
    
    for struct in matched_structures:
        mattersim_e_hull = struct.get('mattersim_e_hull')
        dft_e_hull = struct.get('dft_e_hull')
        
        if mattersim_e_hull is None or dft_e_hull is None:
            continue
        
        if mattersim_e_hull < args.threshold and dft_e_hull < args.threshold:
            stable_structures.append(struct)
    
    print(f"  Found {len(stable_structures)} stable electride candidates\n")
    
    if not stable_structures:
        print("No stable structures found matching criteria.")
        return 0
    
    # Initialize PyXtal database
    if PYXTAL_AVAILABLE:
        # Remove existing database to start fresh
        if output_db.exists():
            output_db.unlink()
        pyxtal_db = database_topology(str(output_db))
    else:
        pyxtal_db = None
    
    # Process each stable structure
    print("Processing stable structures...")
    print("  Loading CONTCARs and extracting spacegroups...\n")
    
    results = {
        'formula': [],
        'composition': [],
        'e0025': [],
        'e05': [],
        'e10': [],
        'band0': [],
        'band1': [],
        'spacegroup': [],
        'mattersim_energy_per_atom': [],
        'mattersim_e_hull': [],
        'vasp_energy_per_atom': [],
        'dft_e_hull': []
    }
    
    db_saved = 0
    db_failed = 0
    processed = 0
    failed = 0
    
    for struct in stable_structures:
        struct_id = struct['structure_id']
        
        # Get workflow data
        sdata = workflow_db['structures'].get(struct_id)
        if not sdata:
            print(f"  {struct_id}: WARNING - not found in workflow database")
            failed += 1
            continue
        
        # Check if relaxation completed
        state = sdata.get('state', '')
        if state not in ['RELAX_DONE', 'RELAX_TMOUT']:
            print(f"  {struct_id}: WARNING - relaxation not completed (state: {state})")
            failed += 1
            continue
        
        # Load CONTCAR
        relax_dir = Path(sdata['relax_dir'])
        contcar_path = relax_dir / 'CONTCAR'
        
        if not contcar_path.exists():
            print(f"  {struct_id}: WARNING - CONTCAR not found")
            failed += 1
            continue
        
        try:
            structure = Structure.from_file(str(contcar_path))
        except Exception as e:
            print(f"  {struct_id}: WARNING - Could not load CONTCAR: {e}")
            failed += 1
            continue
        
        # Extract spacegroup
        spacegroup = get_spacegroup_from_structure(structure)
        
        # Get composition
        composition = sdata.get('composition', struct_id.rsplit('_s', 1)[0])
        
        # Get e_hull values from hull_comparison.json
        mattersim_e_hull = struct.get('mattersim_e_hull')
        dft_e_hull = struct.get('dft_e_hull')
        
        # Get energy_per_atom from separate result files
        mattersim_energy = mattersim_energies.get(struct_id)
        vasp_energy = dft_energies.get(struct_id)
        
        # Get electride volumes
        vol_data = electride_volumes.get(struct_id, {})
        e0025 = vol_data.get('e0025', 0.0)
        e05 = vol_data.get('e05', 0.0)
        e10 = vol_data.get('e10', 0.0)
        band0 = vol_data.get('band0', 0.0)
        band1 = vol_data.get('band1', 0.0)
        
        # Add to results
        results['formula'].append(struct_id)
        results['composition'].append(composition)
        results['e0025'].append(e0025)
        results['e05'].append(e05)
        results['e10'].append(e10)
        results['band0'].append(band0)
        results['band1'].append(band1)
        results['spacegroup'].append(spacegroup)
        results['mattersim_energy_per_atom'].append(mattersim_energy)
        results['mattersim_e_hull'].append(mattersim_e_hull)
        results['vasp_energy_per_atom'].append(vasp_energy)
        results['dft_e_hull'].append(dft_e_hull)
        
        # Save to database
        if pyxtal_db is not None:
            saved = save_to_database(
                struct_id, structure, composition, spacegroup,
                mattersim_e_hull, dft_e_hull,
                mattersim_energy, vasp_energy,
                e0025, e05, e10, band0, band1, pyxtal_db
            )
            if saved:
                db_saved += 1
            else:
                db_failed += 1
        
        processed += 1
        
        if processed % 10 == 0:
            print(f"  Processed {processed}/{len(stable_structures)} structures...")
    
    # Create DataFrame and save CSV
    df = pd.DataFrame(results)
    
    # Sort by spacegroup in descending order (highest symmetry first)
    df = df.sort_values(by='spacegroup', ascending=False).reset_index(drop=True)
    
    df.to_csv(output_csv, index=False)
    
    # Commit database
    if pyxtal_db is not None and hasattr(pyxtal_db, 'db'):
        if hasattr(pyxtal_db.db, 'commit'):
            pyxtal_db.db.commit()
    
    # Print summary
    print("\n" + "="*70)
    print("Summary")
    print("="*70)
    print(f"Total matched structures: {len(matched_structures)}")
    print(f"Stable (both e_hull < {args.threshold}): {len(stable_structures)}")
    print(f"Successfully processed: {processed}")
    print(f"Failed to process: {failed}")
    if PYXTAL_AVAILABLE:
        print(f"Saved to database: {db_saved}")
        print(f"Database save failed: {db_failed}")
    print()
    print(f"Output CSV: {output_csv}")
    if PYXTAL_AVAILABLE:
        print(f"Output database: {output_db}")
    print("="*70 + "\n")
    
    # Show results table
    if len(df) > 0:
        print("Stable Electride Candidates (sorted by spacegroup, descending):")
        # Select columns for display
        display_cols = ['formula', 'composition', 'e0025', 'e05', 'e10', 
                       'band0', 'band1', 'spacegroup', 'mattersim_e_hull', 'dft_e_hull']
        print(tabulate(df[display_cols].head(30), headers='keys', tablefmt='psql', showindex=False))
        
        # Statistics
        print("\n" + "="*70)
        print("Statistics")
        print("="*70)
        print(f"Total stable candidates: {len(df)}")
        print(f"Spacegroup distribution:")
        sg_counts = df['spacegroup'].value_counts().head(10)
        for sg, count in sg_counts.items():
            print(f"  SG {sg}: {count} structures")
        
        print(f"\nEnergy statistics:")
        print(f"  MatterSim e_hull: mean={df['mattersim_e_hull'].mean():.6f}, "
              f"max={df['mattersim_e_hull'].max():.6f} eV/atom")
        print(f"  DFT e_hull: mean={df['dft_e_hull'].mean():.6f}, "
              f"max={df['dft_e_hull'].max():.6f} eV/atom")
        print("="*70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

