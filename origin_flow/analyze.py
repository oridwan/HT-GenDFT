#!/usr/bin/env python3
"""
Analyze completed ELF calculations for electride candidates.

Workflow:
1. Extract PARCHG-* files from PARCHG.tar.gz
2. Run Electride.py analysis
3. Clean up extracted PARCHG-* files
4. PARCHG.tar.gz is ALWAYS preserved for data traceability
5. Save electride candidates to PyXtal database

Features:
- Incremental analysis: skips structures already in CSV or database
- Adds e_above_hull from MatterSim prescreening results
- Adds spacegroup from PyXtal analysis of CONTCAR
- Saves to PyXtal database with adaptive tolerance

Note: Uses MatterSim e_above_hull from prescreening (no DFT hull calculation needed).
      This is faster and more robust than VASP DFT hull calculation.
      
python origin_flow/analyze.py \
  --db VASP-out/workflow.json \
  --bader-exe "$HOME/miniconda3/envs/mattersim/bin/bader" \
  --threshold 0.6 \
  --output VASP-out/electride_analysis.csv \
  --pyxtal-db VASP-out/electride_data.db \
  --prescreening prescreen/VASP_JOBS/prescreening_stability.json \
  --workers 96
  

"""

import sys
import os
import json
import tarfile
import argparse
import multiprocessing as mp
import pandas as pd
from tabulate import tabulate
from pathlib import Path
from Electride import electride
from pyxtal import pyxtal
from pyxtal.db import database_topology
from pymatgen.io.vasp.outputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor


def load_existing_results(csv_path):
    """Load existing CSV results if available."""
    if not os.path.exists(csv_path):
        return None, set()
    
    try:
        df = pd.read_csv(csv_path)
        analyzed_ids = set(df['formula'].values)
        print(f"Found existing results: {len(analyzed_ids)} structures already analyzed")
        return df, analyzed_ids
    except Exception as e:
        print(f"Warning: Could not load existing CSV: {e}")
        return None, set()


def load_existing_database(db_path):
    """
    Load existing PyXtal database if available.
    
    Note: We rely on CSV for tracking already-analyzed structures since PyXtal 
    database_topology API doesn't provide a simple get_all() method. The database 
    is just used for storing structures, not for incremental tracking.
    """
    if not os.path.exists(db_path):
        return None
    
    try:
        db = database_topology(db_path)
        print(f"Found existing database: {db_path}")
        return db
    except Exception as e:
        print(f"Warning: Could not load existing database: {e}")
        return None


def load_prescreening_stability(prescreen_json_path):
    """Load MatterSim prescreening results for e_above_hull values."""
    if not os.path.exists(prescreen_json_path):
        print(f"Warning: Prescreening stability file not found: {prescreen_json_path}")
        return {}
    
    try:
        with open(prescreen_json_path, 'r') as f:
            data = json.load(f)
        
        e_hull_map = {}
        for result in data.get('results', []):
            struct_id = result['structure_id']
            e_hull_map[struct_id] = result.get('energy_above_hull', None)
        
        print(f"Loaded MatterSim prescreening: {len(e_hull_map)} structures with e_above_hull")
        return e_hull_map
    except Exception as e:
        print(f"Warning: Could not load prescreening stability file: {e}")
        return {}


def get_spacegroup_from_contcar(relax_dir):
    """Extract spacegroup from CONTCAR using PyXtal."""
    contcar_path = os.path.join(relax_dir, 'CONTCAR')
    if not os.path.exists(contcar_path):
        return None
    
    try:
        poscar = Poscar.from_file(contcar_path)
        struct = poscar.structure
        
        # Try PyXtal symmetrization with progressive tolerances
        tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
        for tol in tolerances:
            try:
                xtal = pyxtal()
                xtal.from_seed(struct, tol=tol)
                if not xtal.valid:
                    continue
                if len(xtal.check_short_distances(r=0.5)) > 0:
                    continue
                return xtal.group.number
            except:
                continue
        
        # If all symmetrization attempts failed, use P1 (no symmetry)
        return 1
    except Exception as e:
        return None


def load_workflow_db(db_path):
    """Load workflow database and extract ELF_DONE structures."""
    with open(db_path, 'r') as f:
        data = json.load(f)
    
    elf_done = []
    metals = []
    semiconductors = []
    
    for struct_id, sdata in data['structures'].items():
        if sdata['state'] == 'ELF_DONE':
            spe_dir = sdata.get('spe_dir')
            relax_dir = sdata.get('relax_dir', '')
            elf_done.append({
                'id': struct_id,
                'dir': spe_dir,
                'relax_dir': relax_dir,
                'composition': sdata.get('composition', ''),
                'is_semi': sdata.get('is_semiconductor', False)
            })
            if sdata.get('is_semiconductor'):
                semiconductors.append(struct_id)
            else:
                metals.append(struct_id)
    
    print(f"Total ELF_DONE: {len(elf_done)}")
    print(f"  Semiconductors (with PARCHG): {len(semiconductors)}")
    print(f"  Metals (ELFCAR only): {len(metals)}")
    print("")
    
    return elf_done


def check_and_fix_parchg_overflow(parchg_file):
    """
    Check if PARCHG file has VASP overflow markers (***********) and fix them.
    
    VASP outputs *********** when a charge value overflows the fixed-width format.
    Bader crashes on these files. This function replaces overflow markers with -9999.0.
    
    Args:
        parchg_file: Path to PARCHG file
    
    Returns:
        tuple: (has_overflow, was_fixed)
    """
    parchg_path = Path(parchg_file)
    if not parchg_path.exists():
        return False, False
    
    try:
        with open(parchg_path, 'r') as f:
            content = f.read()
        
        # Check for overflow markers
        if '***********' not in content:
            return False, False
        
        # Count overflow occurrences
        overflow_count = content.count('***********')
        
        # Replace overflow markers with a safe value
        # Use -9999.0 padded to 11 chars to maintain format
        fixed_content = content.replace('***********', '-9999.0    ')
        
        # Write fixed content back
        with open(parchg_path, 'w') as f:
            f.write(fixed_content)
        
        return True, True
    except Exception as e:
        return True, False


def extract_parchg_files(elf_dir, struct_id=""):
    """
    Extract PARCHG files from PARCHG.tar.gz for analysis.
    
    Args:
        elf_dir: Path to ELF directory
        struct_id: Structure ID for logging (optional)
    
    Returns:
        bool: True if extraction was performed, False otherwise
    """
    elf_path = Path(elf_dir)
    tar_archive = elf_path / "PARCHG.tar.gz"
    prefix = f"  [{struct_id}]" if struct_id else "  "
    
    # Check if PARCHG-* FILES (not directories) already exist
    existing_parchg_files = [f for f in elf_path.glob("PARCHG-*") if f.is_file()]
    
    if existing_parchg_files:
        # Check existing files for overflow issues
        for parchg_file in existing_parchg_files:
            has_overflow, was_fixed = check_and_fix_parchg_overflow(parchg_file)
            if has_overflow and was_fixed:
                print(f"{prefix} Fixed VASP overflow in {parchg_file.name}")
            elif has_overflow:
                print(f"{prefix} Warning: Could not fix overflow in {parchg_file.name}")
        return False  # Files already available, no extraction needed
    
    # If no tar archive exists
    if not tar_archive.exists():
        print(f"{prefix} Warning: No PARCHG.tar.gz and no PARCHG-* files")
        return False
    
    # Extract tar archive
    try:
        with tarfile.open(tar_archive, 'r:gz') as tar:
            tar.extractall(path=elf_path)
        
        # Verify extraction created expected files
        extracted_files = [f for f in elf_path.glob("PARCHG-*") if f.is_file()]
        if extracted_files:
            print(f"{prefix} Extracted {len(extracted_files)} PARCHG files")
            
            # Check extracted files for overflow issues
            for parchg_file in extracted_files:
                has_overflow, was_fixed = check_and_fix_parchg_overflow(parchg_file)
                if has_overflow and was_fixed:
                    print(f"{prefix} Fixed VASP overflow in {parchg_file.name}")
                elif has_overflow:
                    print(f"{prefix} Warning: Could not fix overflow in {parchg_file.name}")
            
            return True
        else:
            print(f"{prefix} Warning: Extraction completed but no PARCHG-* files found")
            return False
    except Exception as e:
        print(f"{prefix} Warning: Could not extract {tar_archive.name}: {e}")
        return False


def cleanup_parchg_files(elf_dir, was_extracted, struct_id=""):
    """
    Clean up extracted PARCHG-* files after analysis.
    
    IMPORTANT: Original PARCHG.tar.gz is NEVER deleted - only extracted PARCHG-* files are removed.
    This preserves the compressed archive for future data traceability.
    
    Args:
        elf_dir: ELF directory path
        was_extracted: If True, delete extracted files. If False, keep them.
        struct_id: Structure ID for logging (optional)
    """
    if not was_extracted:
        return  # Don't delete if we didn't extract
    
    elf_path = Path(elf_dir)
    tar_archive = elf_path / "PARCHG.tar.gz"
    prefix = f"  [{struct_id}]" if struct_id else "  "
    
    # Safety check: Verify tar.gz exists before deleting extracted files
    if not tar_archive.exists():
        print(f"{prefix} Warning: PARCHG.tar.gz not found, keeping extracted files for safety")
        return
    
    # Find all extracted PARCHG-* files (NOT the tar.gz!)
    parchg_files = sorted(elf_path.glob("PARCHG-*"))
    
    if not parchg_files:
        return  # Nothing to clean
    
    try:
        # Delete only extracted PARCHG-* files (tar.gz is preserved)
        for parchg_file in parchg_files:
            if parchg_file.is_file() and parchg_file.suffix != '.gz':  # Extra safety
                parchg_file.unlink()
        print(f"{prefix} Cleaned up {len(parchg_files)} extracted PARCHG files (tar.gz preserved)")
    except Exception as e:
        print(f"{prefix} Warning: Could not clean up PARCHG files: {e}")


def analyze_structure(struct_info, bader_cmd, threshold, keep_extracted=False):
    """
    Analyze a single structure for electride character.
    
    Args:
        struct_info: Dict with structure info (id, dir, relax_dir, composition, is_semi)
        bader_cmd: Bader executable command
        threshold: ELF threshold
        keep_extracted: If True, keep extracted PARCHG files. If False, delete after analysis.
    
    Returns:
        dict: Result dictionary with all analysis data
    """
    struct_id = struct_info['id']
    elf_dir = struct_info['dir']
    relax_dir = struct_info.get('relax_dir', '')
    composition = struct_info.get('composition', '')
    
    # Print structure being analyzed
    print(f"\nAnalyzing: {struct_id} ({elf_dir})")
    
    result = {
        'struct_id': struct_id,
        'composition': composition,
        'relax_dir': relax_dir,
        'volumes': [0, 0, 0, 0, 0],
        'spacegroup': None,
        'error': None,
        'is_electride': False
    }
    
    # Check if ELFCAR exists
    elfcar_path = os.path.join(elf_dir, 'ELFCAR')
    if not os.path.exists(elfcar_path):
        print(f"  [{struct_id}] ERROR: ELFCAR not found")
        result['error'] = f"ELFCAR not found"
        return result
    
    # Extract PARCHG files from tar.gz archive
    was_extracted = extract_parchg_files(elf_dir, struct_id)
    
    # Get spacegroup from CONTCAR
    if relax_dir:
        result['spacegroup'] = get_spacegroup_from_contcar(relax_dir)
        if result['spacegroup']:
            print(f"  [{struct_id}] Spacegroup: {result['spacegroup']}")
    
    # Run electride analysis
    try:
        test = electride(elf_dir, cmd=bader_cmd, ELF_min=threshold)
        
        # Extract volumes for all 5 PARCHG types: e0025, e05, e10, band0, band1
        volumes = [float(test.volume[i]) for i in range(5)]
        result['volumes'] = volumes
        
        # Clean up extracted files (delete them, tar.gz remains)
        if not keep_extracted:
            cleanup_parchg_files(elf_dir, was_extracted, struct_id)
        
        # Electride criteria: ((e0025 > 10 OR e05 > 10) AND (e10 > 10 OR band0 > 10)) AND (spacegroup > 15)
        is_electride = ((volumes[0] > 10 or volumes[1] > 10) and 
                       (volumes[2] > 10 or volumes[3] > 10) and
                       (result['spacegroup'] is not None and result['spacegroup'] > 15))
        result['is_electride'] = is_electride
        
        if test.error:
            print(f"  [{struct_id}] Warning: {test.error}")
            result['error'] = test.error
        elif is_electride:
            print(f"  [{struct_id}] *** ELECTRIDE CANDIDATE *** (e0025={volumes[0]:.2f}, e05={volumes[1]:.2f}, e10={volumes[2]:.2f}, band0={volumes[3]:.2f}, band1={volumes[4]:.2f})")
        else:
            print(f"  [{struct_id}] Not electride (e0025={volumes[0]:.2f}, e05={volumes[1]:.2f}, e10={volumes[2]:.2f}, band0={volumes[3]:.2f}, band1={volumes[4]:.2f})")
            
    except Exception as e:
        print(f"  [{struct_id}] ERROR: {e}")
        result['error'] = str(e)
        # Clean up on error
        if not keep_extracted:
            cleanup_parchg_files(elf_dir, was_extracted, struct_id)
    
    return result


def analyze_structure_wrapper(args):
    """Wrapper for multiprocessing - unpacks arguments."""
    struct_info, bader_cmd, threshold, keep_extracted = args
    return analyze_structure(struct_info, bader_cmd, threshold, keep_extracted)


def save_to_database(struct_id, relax_dir, composition, e_above_hull, is_electride, db, tolerances):
    """
    Save structure to PyXtal database with adaptive tolerance.
    
    Note: Database saving is optional - the main output is the CSV file.
    Failures here don't affect the CSV analysis results.
    """
    if db is None:
        return False
    
    contcar_path = os.path.join(relax_dir, 'CONTCAR')
    if not os.path.exists(contcar_path):
        return False
    
    try:
        poscar = Poscar.from_file(contcar_path)
        struct = poscar.structure
        
        # Try progressive tolerances
        for tol in tolerances:
            try:
                xtal = pyxtal()
                xtal.from_seed(struct, tol=tol)
                if not xtal.valid:
                    continue
                if len(xtal.check_short_distances(r=0.5)) > 0:
                    continue
                
                db.add_xtal(
                    xtal,
                    kvp={
                        'structure_id': struct_id,
                        'e_above_hull': e_above_hull,
                        'composition': composition,
                        'symmetrized': True,
                        'is_electride': is_electride
                    }
                )
                return True
            except:
                continue
        
        # If all tolerances failed, save structure without symmetrization
        # Use ASE database directly to avoid PyXtal artifacts
        try:
            # Convert pymatgen Structure to ASE Atoms
            adaptor = AseAtomsAdaptor()
            atoms = adaptor.get_atoms(struct)
            
            db.db.write(
                atoms,
                structure_id=struct_id,
                e_above_hull=e_above_hull,
                composition=composition,
                symmetrized=False,
                is_electride=is_electride,
                space_group_number=1
            )
            print(f"    Saved to database without symmetrization for {struct_id}")
            return True
        except Exception as e2:
            print(f"    Warning: Could not save {struct_id} even without symmetrization: {e2}")
            return False
    except Exception as e:
        print(f"  Warning: Could not save {struct_id} to database: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Analyze completed ELF calculations for electride candidates (incremental mode supported)"
    )
    parser.add_argument(
        '--db',
        type=str,
        required=True,
        help="Path to workflow database (workflow.json)"
    )
    parser.add_argument(
        '--bader-exe',
        type=str,
        default='bader',
        help="Path to bader executable (default: bader)"
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=0.6,
        help="ELF threshold (default: 0.6)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='electride_analysis.csv',
        help="Output CSV file (default: electride_analysis.csv)"
    )
    parser.add_argument(
        '--pyxtal-db',
        type=str,
        default='electride_data.db',
        help="PyXtal database file (default: electride_data.db)"
    )
    parser.add_argument(
        '--prescreening',
        type=str,
        default='./VASP_JOBS/prescreening_stability.json',
        help="MatterSim prescreening results JSON (default: ./VASP_JOBS/prescreening_stability.json)"
    )
    parser.add_argument(
        '--keep-extracted',
        action='store_true',
        help="Keep extracted PARCHG-* files after analysis (default: delete to save space). "
             "Note: PARCHG.tar.gz is always preserved regardless of this flag."
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=None,
        help="Number of parallel workers (default: number of CPU cores). "
             "Set to 1 for sequential processing."
    )
    
    args = parser.parse_args()
    
    # Determine number of workers
    if args.workers is None:
        n_workers = mp.cpu_count()
    else:
        n_workers = max(1, args.workers)
    
    # Check if database exists
    if not os.path.exists(args.db):
        print(f"ERROR: Workflow database not found: {args.db}")
        sys.exit(1)
    
    print("="*70)
    print("Electride Analysis - Completed ELF Calculations")
    print("="*70)
    print(f"Workflow database: {args.db}")
    print(f"Bader executable: {args.bader_exe}")
    print(f"ELF threshold: {args.threshold}")
    print(f"Output CSV: {args.output}")
    print(f"PyXtal database: {args.pyxtal_db}")
    print(f"Prescreening results: {args.prescreening}")
    print(f"Parallel workers: {n_workers}")
    print()
    print("Electride Criteria: ((e0025 > 0 OR e05 > 0) AND (e10 > 0 OR band0 > 0)) AND (spacegroup > 15)")
    print()
    print("Workflow: Extract → Analyze → Cleanup → Save to DB")
    print("  - PARCHG.tar.gz is always preserved")
    print("  - Only temporary extracted files are deleted")
    print("  - Incremental mode: skips already-analyzed structures")
    print("  - Parallel processing for faster analysis")
    print("="*70)
    print("")
    
    # Load existing results (CSV for tracking, PyXtal database for storage)
    existing_df, analyzed_csv_ids = load_existing_results(args.output)
    existing_pyxtal_db = load_existing_database(args.pyxtal_db)
    
    # Use CSV for tracking already-analyzed structures
    skip_ids = analyzed_csv_ids
    
    # Load MatterSim prescreening stability results
    e_hull_map = load_prescreening_stability(args.prescreening)
    
    # Load workflow database
    elf_structures = load_workflow_db(args.db)
    
    if not elf_structures:
        print("No completed ELF calculations found.")
        sys.exit(0)
    
    # Filter out already-analyzed structures
    structures_to_analyze = [s for s in elf_structures if s['id'] not in skip_ids]
    
    print(f"Found {len(elf_structures)} ELF_DONE structures")
    print(f"  Already analyzed: {len(skip_ids)}")
    print(f"  New to analyze: {len(structures_to_analyze)}")
    print("")
    
    if len(structures_to_analyze) == 0:
        print("No new structures to analyze. All up to date!")
        if existing_df is not None:
            print(f"\nExisting results in {args.output}:")
            print(tabulate(existing_df.head(20), headers='keys', tablefmt='psql', showindex=False))
        sys.exit(0)
    
    # Initialize or open PyXtal database
    if existing_pyxtal_db:
        pyxtal_db = existing_pyxtal_db
    else:
        pyxtal_db = database_topology(args.pyxtal_db)
    
    # Analyze new structures with multiprocessing
    bader_cmd = args.bader_exe + ' '
    
    # Prepare arguments for parallel processing
    analysis_args = [
        (struct_info, bader_cmd, args.threshold, args.keep_extracted)
        for struct_info in structures_to_analyze
    ]
    
    print(f"Starting parallel analysis with {n_workers} workers...")
    print("")
    
    # Run analysis in parallel
    if n_workers > 1 and len(structures_to_analyze) > 1:
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(analyze_structure_wrapper, analysis_args)
    else:
        # Sequential processing
        results = [analyze_structure_wrapper(arg) for arg in analysis_args]
    
    # Process results and print summary
    new_results = {
        'formula': [],
        'composition': [],
        'e0025': [],
        'e05': [],
        'e10': [],
        'band0': [],
        'band1': [],
        'spacegroup': [],
        'e_above_hull': []
    }
    
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    db_saved_count = 0
    db_failed_count = 0
    electride_count = 0
    error_count = 0
    
    print("\nCollecting results...")
    for result in results:
        struct_id = result['struct_id']
        volumes = result['volumes']
        spacegroup = result['spacegroup']
        composition = result['composition']
        is_electride = result['is_electride']
        error = result['error']
        relax_dir = result['relax_dir']
        
        # Get e_above_hull from prescreening
        e_hull = e_hull_map.get(struct_id, None)
        
        # Count statistics
        if error:
            error_count += 1
        if is_electride:
            electride_count += 1
        
        # Save to results
        new_results['formula'].append(struct_id)
        new_results['composition'].append(composition)
        new_results['e0025'].append(volumes[0])
        new_results['e05'].append(volumes[1])
        new_results['e10'].append(volumes[2])
        new_results['band0'].append(volumes[3])
        new_results['band1'].append(volumes[4])
        new_results['spacegroup'].append(spacegroup if spacegroup else 0)
        new_results['e_above_hull'].append(e_hull if e_hull is not None else float('inf'))
        
        # Save to PyXtal database (optional, CSV is the main output)
        if is_electride and relax_dir and pyxtal_db is not None:
            saved = save_to_database(
                struct_id, relax_dir, composition,
                e_hull, is_electride, 
                pyxtal_db, tolerances
            )
            if saved:
                db_saved_count += 1
            else:
                db_failed_count += 1
    
    print(f"Analysis complete: {electride_count} electrides, {error_count} errors")
    
    # Commit database
    if hasattr(pyxtal_db, 'db') and hasattr(pyxtal_db.db, 'commit'):
        pyxtal_db.db.commit()
    
    # Create DataFrame for new results
    new_df = pd.DataFrame(new_results)
    
    # Merge with existing results if any
    if existing_df is not None:
        df = pd.concat([existing_df, new_df], ignore_index=True)
    else:
        df = new_df
    
    # Sort by e_above_hull (low to high)
    df = df.sort_values(by='e_above_hull', ascending=True).reset_index(drop=True)
    
    # Save updated CSV
    df.to_csv(args.output, index=False)
    
    print("")
    print("="*70)
    print("Summary Results")
    print("="*70)
    print(f"Total structures in CSV: {len(df)}")
    print(f"New structures analyzed: {len(new_df)}")
    print(f"PyXtal database: {db_saved_count} saved, {db_failed_count} failed")
    print("")
    print(f"Results saved to: {args.output}")
    print(f"Database saved to: {args.pyxtal_db}")
    print("")
    
    # Show top results (sorted by e_above_hull)
    print("Top 20 structures (sorted by e_above_hull):")
    print(tabulate(df.head(20), headers='keys', tablefmt='psql', showindex=False))
    
    # Summary statistics: ((e0025 > 0 OR e05 > 0) AND (e10 > 0 OR band0 > 0)) AND (spacegroup > 15)
    electrides = df[((df['e0025'] > 0) | (df['e05'] > 0)) & ((df['e10'] > 0) | (df['band0'] > 0)) & (df['spacegroup'] > 15)]
    print("")
    print("="*70)
    print("Statistics")
    print("="*70)
    print(f"Total structures analyzed: {len(df)}")
    print(f"Electride candidates: {len(electrides)} ({100*len(electrides)/len(df):.1f}%)")
    print(f"  (Criteria: ((e0025 > 0 OR e05 > 0) AND (e10 > 0 OR band0 > 0)) AND (spacegroup > 15))")
    
    if len(electrides) > 0:
        print("")
        print("Top electride candidates (sorted by e_above_hull):")
        print(tabulate(
            electrides.head(10),
            headers='keys',
            tablefmt='psql',
            showindex=False
        ))
    
    print("="*70)


if __name__ == '__main__':
    main()

