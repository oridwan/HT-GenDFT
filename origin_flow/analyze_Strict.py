#!/usr/bin/env python3
"""
Analyze completed ELF calculations with strict candidate filtering.

Workflow:
1. Extract PARCHG-* files from PARCHG.tar.gz
2. Run Electride.py analysis
3. Clean up extracted PARCHG-* files
4. PARCHG.tar.gz is ALWAYS preserved for data traceability
5. Save strict electride candidates to PyXtal database

Features:
- Incremental analysis: skips structures already in CSV or database
- Adds e_above_hull from MatterSim prescreening results
- Adds spacegroup from PyXtal analysis of CONTCAR
- Saves to PyXtal database with adaptive tolerance
- Writes two CSVs: all analyzed structures and strict candidates only
- Supports excluding candidates with selected elements
- Supports configurable space-group selection criteria

Note: Uses MatterSim e_above_hull from prescreening (no DFT hull calculation needed).
      This is faster and more robust than VASP DFT hull calculation.
      
python origin_flow/analyze_Strict.py \
  --db VASP-out-Boron/workflow.json \
  --bader-exe "$HOME/miniconda3/envs/mattersim/bin/bader" \
  --threshold 0.6 \
  --output VASP-out-Boron/electride_analysis.csv \
  --candidate-output VASP-out-Boron/electride_analysis_candidates-Strict.csv \
  --pyxtal-db VASP-out-Boron/electride_data.db \
  --spacegroup-criteria "spg>15" \
  --prescreening prescreen_new/VASP_JOBS_Boron/prescreening_stability.json \
  --workers 32
  
--exclude-elements Cs Rb \
"""

import sys
import os
import json
import tarfile
import argparse
import multiprocessing as mp
import re
import pandas as pd
from tabulate import tabulate
from pathlib import Path
from Electride import electride
from pyxtal import pyxtal
from pyxtal.db import database_topology
from pymatgen.io.vasp.outputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

STRICT_VOLUME_THRESHOLD = 10.0


def derive_candidate_output_path(output_path):
    """Derive candidate-only CSV path from the all-structures CSV path."""
    output_path = Path(output_path)
    return output_path.with_name(f"{output_path.stem}_candidates{output_path.suffix}")


def structure_id_to_composition(struct_id):
    """Convert a structure ID like K1Sr1B2H8_s014 to its composition string."""
    if not struct_id:
        return ""
    return struct_id.rsplit('_s', 1)[0]


def extract_elements(formula):
    """Extract unique element symbols from a composition string."""
    if not formula:
        return []
    return sorted(set(re.findall(r'([A-Z][a-z]?)', str(formula))))


def normalize_excluded_elements(raw_values):
    """Normalize excluded-element CLI values."""
    elements = []
    seen = set()
    for value in raw_values or []:
        for token in str(value).replace(',', ' ').split():
            token = token.strip()
            if not token:
                continue
            if not re.fullmatch(r'[A-Z][a-z]?', token):
                raise ValueError(f"Invalid element symbol in --exclude-elements: {token}")
            if token not in seen:
                seen.add(token)
                elements.append(token)
    return elements


def normalize_spacegroup_criteria(criteria):
    """Normalize space-group criteria strings like 'spg>15' or '>=38'."""
    text = str(criteria or "").strip()
    if not text:
        return ""
    text = text.replace(" ", "")
    text = re.sub(r'^(spacegroup|spg)', '', text, flags=re.IGNORECASE)
    return text


def matches_spacegroup_criteria(spacegroup, criteria):
    """Check whether a space group matches a simple criterion expression."""
    normalized = normalize_spacegroup_criteria(criteria)
    if not normalized:
        return True
    if spacegroup is None:
        return False

    try:
        sg = int(spacegroup)
    except Exception:
        return False

    range_match = re.fullmatch(r'(\d+)-(\d+)', normalized)
    if range_match:
        lo = int(range_match.group(1))
        hi = int(range_match.group(2))
        return lo <= sg <= hi

    comp_match = re.fullmatch(r'(>=|<=|==|!=|>|<|=)?(\d+)', normalized)
    if not comp_match:
        raise ValueError(
            f"Unsupported --spacegroup-criteria '{criteria}'. Examples: 'spg>15', '>=38', '38-62', '=166'"
        )

    op = comp_match.group(1) or '='
    value = int(comp_match.group(2))
    if op == '>':
        return sg > value
    if op == '>=':
        return sg >= value
    if op == '<':
        return sg < value
    if op == '<=':
        return sg <= value
    if op in ('=', '=='):
        return sg == value
    if op == '!=':
        return sg != value
    raise ValueError(f"Unsupported space-group operator: {op}")


def evaluate_candidate(volumes, spacegroup, composition, excluded_elements, spacegroup_criteria):
    """Evaluate strict electride-candidate filters for one structure."""
    excluded_found = [elem for elem in excluded_elements if elem in set(extract_elements(composition))]
    passes_volume_filter = (
        (volumes[0] > STRICT_VOLUME_THRESHOLD or volumes[1] > STRICT_VOLUME_THRESHOLD) and
        (volumes[2] > STRICT_VOLUME_THRESHOLD or volumes[3] > STRICT_VOLUME_THRESHOLD)
    )
    passes_spacegroup_filter = matches_spacegroup_criteria(spacegroup, spacegroup_criteria)
    passes_excluded_elements = len(excluded_found) == 0
    is_candidate = passes_volume_filter and passes_spacegroup_filter and passes_excluded_elements
    return {
        'excluded_elements_found': ",".join(excluded_found),
        'passes_volume_filter': passes_volume_filter,
        'passes_spacegroup_filter': passes_spacegroup_filter,
        'passes_excluded_elements': passes_excluded_elements,
        'is_candidate': is_candidate,
    }


def enrich_results_dataframe(df, excluded_elements, spacegroup_criteria):
    """Add strict-filter columns and return all-results and candidates-only DataFrames."""
    if df is None:
        empty_df = pd.DataFrame()
        return empty_df, empty_df

    df = df.copy()
    if 'composition' not in df.columns:
        df['composition'] = df['formula'].map(structure_id_to_composition)
    else:
        df['composition'] = df['composition'].fillna(df['formula'].map(structure_id_to_composition))

    for col in ['e0025', 'e05', 'e10', 'band0', 'band1', 'spacegroup']:
        if col not in df.columns:
            df[col] = 0
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

    if 'e_above_hull' not in df.columns:
        df['e_above_hull'] = float('inf')
    df['e_above_hull'] = pd.to_numeric(df['e_above_hull'], errors='coerce').fillna(float('inf'))

    candidate_flags = [
        evaluate_candidate(
            [row.e0025, row.e05, row.e10, row.band0, row.band1],
            row.spacegroup,
            row.composition,
            excluded_elements,
            spacegroup_criteria,
        )
        for row in df.itertuples(index=False)
    ]
    flags_df = pd.DataFrame(candidate_flags)
    for col in flags_df.columns:
        df[col] = flags_df[col].values

    df['excluded_elements_filter'] = ",".join(excluded_elements)
    df['spacegroup_criteria'] = normalize_spacegroup_criteria(spacegroup_criteria)
    df['strict_volume_threshold'] = STRICT_VOLUME_THRESHOLD

    df = df.sort_values(by='e_above_hull', ascending=True).reset_index(drop=True)
    candidates_df = df[df['is_candidate']].copy().reset_index(drop=True)
    return df, candidates_df

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


def analyze_structure(struct_info, bader_cmd, threshold, excluded_elements, spacegroup_criteria, keep_extracted=False):
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
        'is_electride': False,
        'is_candidate': False,
        'excluded_elements_found': '',
        'passes_volume_filter': False,
        'passes_spacegroup_filter': False,
        'passes_excluded_elements': True,
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
        
        candidate_info = evaluate_candidate(
            volumes,
            result['spacegroup'],
            composition or structure_id_to_composition(struct_id),
            excluded_elements,
            spacegroup_criteria,
        )
        result.update(candidate_info)
        result['is_electride'] = candidate_info['is_candidate']
        result['is_candidate'] = candidate_info['is_candidate']

        if test.error:
            print(f"  [{struct_id}] Warning: {test.error}")
            result['error'] = test.error
        elif result['is_candidate']:
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
    struct_info, bader_cmd, threshold, excluded_elements, spacegroup_criteria, keep_extracted = args
    return analyze_structure(
        struct_info,
        bader_cmd,
        threshold,
        excluded_elements,
        spacegroup_criteria,
        keep_extracted,
    )


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
        help="All-structures CSV output file (default: electride_analysis.csv)"
    )
    parser.add_argument(
        '--candidate-output',
        type=str,
        default=None,
        help="Candidate-only CSV output file (default: <output>_candidates.csv)"
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
    parser.add_argument(
        '--exclude-elements',
        type=str,
        nargs='*',
        default=[],
        help="Exclude structures containing any of these elements from the candidate CSV, e.g. --exclude-elements Cs Rb"
    )
    parser.add_argument(
        '--spacegroup-criteria',
        type=str,
        default='spg>15',
        help="Space-group criterion for candidate selection, e.g. 'spg>15', '>=38', '38-62' (default: spg>15)"
    )
    
    args = parser.parse_args()
    args.exclude_elements = normalize_excluded_elements(args.exclude_elements)
    candidate_output = args.candidate_output or str(derive_candidate_output_path(args.output))
    spacegroup_criteria = normalize_spacegroup_criteria(args.spacegroup_criteria)
    
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
    print(f"All-structures CSV: {args.output}")
    print(f"Candidate CSV: {candidate_output}")
    print(f"PyXtal database: {args.pyxtal_db}")
    print(f"Prescreening results: {args.prescreening}")
    print(f"Parallel workers: {n_workers}")
    print(f"Excluded elements: {args.exclude_elements if args.exclude_elements else 'None'}")
    print(f"Space-group criteria: {spacegroup_criteria or 'None'}")
    print()
    print(
        f"Strict candidate criteria: "
        f"((e0025 > {STRICT_VOLUME_THRESHOLD:g} OR e05 > {STRICT_VOLUME_THRESHOLD:g}) "
        f"AND (e10 > {STRICT_VOLUME_THRESHOLD:g} OR band0 > {STRICT_VOLUME_THRESHOLD:g})) "
        f"AND ({spacegroup_criteria or 'no space-group filter'}) "
        f"AND (exclude none of {args.exclude_elements if args.exclude_elements else '[]'})"
    )
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
    
    pyxtal_db = None
    results = []
    if len(structures_to_analyze) == 0:
        print("No new structures to analyze. Regenerating CSV outputs from existing results.")
        if existing_df is not None:
            print(f"\nExisting results in {args.output}:")
            print(tabulate(existing_df.head(20), headers='keys', tablefmt='psql', showindex=False))
    else:
        # Initialize or open PyXtal database
        if existing_pyxtal_db:
            pyxtal_db = existing_pyxtal_db
        else:
            pyxtal_db = database_topology(args.pyxtal_db)
        
        # Analyze new structures with multiprocessing
        bader_cmd = args.bader_exe + ' '
        
        # Prepare arguments for parallel processing
        analysis_args = [
            (struct_info, bader_cmd, args.threshold, args.exclude_elements, spacegroup_criteria, args.keep_extracted)
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
    
    print(f"Analysis complete: {electride_count} strict candidates, {error_count} errors")
    
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

    df, candidates_df = enrich_results_dataframe(df, args.exclude_elements, spacegroup_criteria)

    # Save updated CSVs
    df.to_csv(args.output, index=False)
    candidates_df.to_csv(candidate_output, index=False)
    
    print("")
    print("="*70)
    print("Summary Results")
    print("="*70)
    print(f"Total structures in all-results CSV: {len(df)}")
    print(f"Total structures in candidate CSV: {len(candidates_df)}")
    print(f"New structures analyzed: {len(new_df)}")
    print(f"PyXtal database: {db_saved_count} saved, {db_failed_count} failed")
    print("")
    print(f"All-results CSV saved to: {args.output}")
    print(f"Candidate CSV saved to: {candidate_output}")
    print(f"Database saved to: {args.pyxtal_db}")
    print("")
    
    # Show top results (sorted by e_above_hull)
    print("Top 20 structures (sorted by e_above_hull):")
    print(tabulate(df.head(20), headers='keys', tablefmt='psql', showindex=False))
    
    print("")
    print("="*70)
    print("Statistics")
    print("="*70)
    print(f"Total structures analyzed: {len(df)}")
    candidate_fraction = (100 * len(candidates_df) / len(df)) if len(df) else 0.0
    print(f"Strict electride candidates: {len(candidates_df)} ({candidate_fraction:.1f}%)")
    print(
        f"  Criteria: ((e0025 > {STRICT_VOLUME_THRESHOLD:g} OR e05 > {STRICT_VOLUME_THRESHOLD:g}) "
        f"AND (e10 > {STRICT_VOLUME_THRESHOLD:g} OR band0 > {STRICT_VOLUME_THRESHOLD:g})) "
        f"AND ({spacegroup_criteria or 'no space-group filter'})"
    )
    if args.exclude_elements:
        print(f"  Excluding elements: {', '.join(args.exclude_elements)}")
    
    if len(candidates_df) > 0:
        print("")
        print("Top strict electride candidates (sorted by e_above_hull):")
        print(tabulate(
            candidates_df.head(10),
            headers='keys',
            tablefmt='psql',
            showindex=False
        ))
    
    print("="*70)


if __name__ == '__main__':
    main()
