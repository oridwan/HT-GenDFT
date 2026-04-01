#!/usr/bin/env python3
"""
Compute MatterSim energy_above_hull for origin-flow electride candidates only.

This follows the same workflow as refined_flow/compute_mattersim_e_hull.py:
1. Queries MP API for GGA/GGA+U reference phases (using legacy pymatgen MPRester)
2. Relaxes MP phases with MatterSim (fmax=0.001, max_steps=800)
3. Loads VASP-relaxed CONTCARs for shortlisted origin electride candidates
4. Relaxes those candidate structures with MatterSim (same tight convergence)
5. Computes energy_above_hull using MatterSim-relaxed MP phases
6. Outputs candidate-only MatterSim hull results

Usage:
    python3 origin_flow/compute_mattersim_e_hull_origin.py \
        --vasp-jobs VASP-out-Boron \
        --candidate-csv VASP-out-Boron/electride_analysis_candidates-Strict.csv \
        --device cuda \
        --output mattersim_stability_candidates.json \
        --pure-pbe
"""

import os
import sys
import csv
import json
import argparse
import warnings
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Structure, Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar

from ase.optimize import FIRE
try:
    from ase.filters import UnitCellFilter
except ImportError:
    from ase.constraints import UnitCellFilter

try:
    from pyxtal import pyxtal
    PYXTAL_AVAILABLE = True
except ImportError:
    PYXTAL_AVAILABLE = False

try:
    from mattersim.forcefield import MatterSimCalculator
    MATTERSIM_AVAILABLE = True
except ImportError:
    MATTERSIM_AVAILABLE = False

try:
    from pymatgen.ext.matproj import MPRester
except ImportError:
    print("ERROR: pymatgen package with MPRester is required")
    print("Install with: pip install pymatgen")
    sys.exit(1)

import numpy as np

warnings.filterwarnings('ignore', category=UserWarning, message='.*POTCAR data with symbol.*')
warnings.filterwarnings('ignore', message='Using UFloat objects with std_dev==0')
warnings.filterwarnings('ignore', category=DeprecationWarning, module='pkg_resources')
warnings.filterwarnings('ignore', category=UserWarning, message='.*pkg_resources is deprecated.*')


def relax_structure_mattersim(pmg_struct, calculator, structure_id=None, fmax=0.001, max_steps=800):
    """
    Relax structure using MatterSim + FIRE optimizer (tighter convergence than prescreen).
    
    Args:
        pmg_struct: Pymatgen Structure object
        calculator: MatterSimCalculator instance (reused)
        structure_id: Structure ID for error reporting
        fmax: Force convergence criterion (eV/Angstrom) - tighter than prescreen
        max_steps: Maximum optimization steps - more than prescreen
    
    Returns:
        tuple: (relaxed_structure, energy_per_atom, error_message)
    
    Note:
        - Symmetrizes structure using PyXtal with progressive tolerance (same as prescreen.py)
        - Falls back to direct conversion if all tolerances fail
        - Uses UnitCellFilter to relax both cell and atomic positions
    """
    sid_prefix = f"[{structure_id}] " if structure_id else ""
    
    # Symmetrize structure using PyXtal with progressive tolerance
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    atoms = None
    
    for tol in tolerances:
        try:
            xtal = pyxtal()
            xtal.from_seed(pmg_struct, tol=tol)
            if not xtal.valid:
                continue
            if len(xtal.check_short_distances(r=0.5)) > 0:
                continue
            atoms = xtal.to_ase()
            break
        except Exception:
            continue
    
    # If all tolerances failed, try direct conversion without symmetrization
    if atoms is None:
        try:
            adaptor = AseAtomsAdaptor()
            atoms = adaptor.get_atoms(pmg_struct)
        except Exception as e:
            return None, None, f"{sid_prefix}Structure conversion failed: {e}"
    
    # Attach calculator
    atoms.calc = calculator
    
    # Check initial energy
    try:
        initial_energy = atoms.get_potential_energy()
        if not np.isfinite(initial_energy):
            return None, None, f"{sid_prefix}Initial energy not finite: {initial_energy}"
    except Exception as e:
        return None, None, f"{sid_prefix}Initial energy calculation failed: {e}"
    
    # Relax with error handling
    try:
        ecf = UnitCellFilter(atoms)
        dyn = FIRE(ecf, a=0.1, logfile=None)
        dyn.run(fmax=fmax, steps=max_steps)
        
        # Get final energy
        energy = atoms.get_potential_energy()
        
        # Check for NaN/inf
        if not np.isfinite(energy):
            try:
                del atoms, ecf, dyn
            except:
                pass
            return None, None, f"{sid_prefix}Final energy not finite: {energy}"
        
        energy_per_atom = energy / len(atoms)
        
        # Convert back to pymatgen
        adaptor = AseAtomsAdaptor()
        relaxed_structure = adaptor.get_structure(atoms)
        
    except Exception as e:
        try:
            del atoms
            if 'ecf' in locals():
                del ecf
            if 'dyn' in locals():
                del dyn
        except:
            pass
        return None, None, f"{sid_prefix}Relaxation failed: {e}"
    
    # Clean up
    del atoms, ecf, dyn
    
    return relaxed_structure, energy_per_atom, None


def load_workflow_database(db_path):
    """Load workflow.json database."""
    with open(db_path, 'r') as f:
        return json.load(f)


def load_candidate_structure_ids(csv_path):
    """Load candidate structure IDs from a CSV file."""
    csv_path = Path(csv_path)
    with csv_path.open(newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        id_col = 'formula' if 'formula' in fieldnames else ('structure_id' if 'structure_id' in fieldnames else None)
        if id_col is None:
            raise ValueError(
                f"Candidate CSV must contain 'formula' or 'structure_id': {csv_path}"
            )

        structure_ids = []
        seen = set()
        for row in reader:
            struct_id = (row.get(id_col) or '').strip()
            if not struct_id or struct_id in seen:
                continue
            seen.add(struct_id)
            structure_ids.append(struct_id)

    return structure_ids


def resolve_path_with_fallbacks(path_str, base_dir=None):
    """Resolve a possibly-relative path against a base dir."""
    path = Path(path_str).expanduser()
    if path.is_absolute():
        return path
    if base_dir is None:
        return path
    base_dir = Path(base_dir)
    direct = base_dir / path
    if direct.exists():
        return direct
    by_name = base_dir / path.name
    if by_name.exists():
        return by_name
    return direct


def save_candidate_relaxed_structure(structure, output_root, structure_id):
    """Save a MatterSim-relaxed candidate structure as CONTCAR."""
    output_root = Path(output_root)
    structure_dir = output_root / str(structure_id)
    structure_dir.mkdir(parents=True, exist_ok=True)
    contcar_path = structure_dir / "CONTCAR"
    Poscar(structure).write_file(str(contcar_path))
    return str(contcar_path.resolve())


def save_mp_relaxed_structure(structure, output_root, mp_id):
    """Save a MatterSim-relaxed MP reference structure as CIF."""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    cif_path = output_root / f"{mp_id}.cif"
    CifWriter(structure).write_file(str(cif_path))
    return str(cif_path.resolve())


def load_mp_cache(cache_file):
    """
    Load MP cache from file.
    
    Returns:
        dict: Cache data keyed by (chemsys, entry_id)
    """
    cache_file = Path(cache_file)
    cached_data = {}
    
    if cache_file.exists():
        with open(cache_file, 'r') as f:
            cache_list = json.load(f)
            for item in cache_list:
                key = (item.get('chemsys', ''), item['entry_id'])
                cached_data[key] = item
    
    return cached_data


def save_mp_cache(cache_file, entries_data):
    """
    Save MP cache entries to file.
    
    Args:
        cache_file: Path to cache file
        entries_data: List of cache entries to save
    """
    cache_file = Path(cache_file)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(cache_file, 'w') as f:
        json.dump(entries_data, f, indent=2)
    
    print(f"    Saved cache: {cache_file} ({len(entries_data)} total entries)")


def get_mp_stable_phases_mattersim(
    chemsys,
    mp_api_key,
    cached_data,
    calculator,
    mp_relaxed_dir=None,
    fmax=0.001,
    max_steps=800,
    pure_pbe=False,
):
    """
    Get MP GGA phases and relax with MatterSim using tight convergence.
    
    Uses legacy pymatgen.ext.matproj.MPRester for complete GGA entry coverage.
    Strict filtering: only entries with '-GGA' or '-GGA+U' suffix.
    Tighter convergence than prescreen (fmax=0.001, max_steps=800) to match refined VASP.
    
    Args:
        chemsys: Chemical system string (e.g., 'B-Li-N')
        mp_api_key: Materials Project API key
        cached_data: Pre-loaded cache dict keyed by (chemsys, entry_id)
        calculator: MatterSimCalculator instance (reused)
        fmax: Force convergence criterion (default: 0.001, tighter than prescreen)
        max_steps: Maximum optimization steps (default: 800, more than prescreen)
        pure_pbe: If True, filter to GGA-PBE only (exclude PBE+U)
                  If False (default), accept both PBE and PBE+U
    
    Returns:
        tuple: (List of PDEntry objects, List of new cache entries to add)
    """
            
    # Extract entries for this chemsys from cache
    elements = chemsys.split('-')
    entries = []
    cached_count = 0
    
    for key, item in cached_data.items():
        item_chemsys = key[0]
        item_elements = sorted(item_chemsys.split('-'))
        if set(item_elements).issubset(set(elements)):
            entry = PDEntry(
                composition=Composition(item['composition']),
                energy=item['energy'],
                name=item['entry_id']
            )
            entries.append(entry)
            cached_count += 1
    
    # Check if we have all required subsystems cached
    required_chemsys = set()
    for n in range(1, len(elements) + 1):
        from itertools import combinations
        for combo in combinations(elements, n):
            required_chemsys.add('-'.join(sorted(combo)))
    
    cached_chemsys = set(item[0] for item in cached_data.keys())
    missing_chemsys = required_chemsys - cached_chemsys
    
    if cached_count > 0 and not missing_chemsys:
        print(f"    Using cached data: {cached_count} GGA phases for {chemsys}")
        return entries, []
    
    # Need to query MP for missing data
    print(f"    Querying MP for GGA phases in {chemsys} (using legacy API for complete coverage)...")
    if missing_chemsys:
        print(f"      Missing subsystems: {sorted(missing_chemsys)}")
    
    try:
        # Use legacy pymatgen MPRester which returns complete GGA entries
        mpr = MPRester(mp_api_key)
        
        # Get ALL entries in this chemical system
        computed_entries = mpr.get_entries_in_chemsys(elements)
        
        print(f"    Retrieved {len(computed_entries)} entries from MP (filtering for GGA...)")
        
        # Filter for GGA only (entry_id ending with '-GGA' or '-GGA+U')
        mp_phases = []
        seen_entries = {}
        skipped_structure_retrieval = []
        
        for comp_entry in computed_entries:
            entry_id = str(comp_entry.entry_id)
            
            # Skip if already seen
            if entry_id in seen_entries:
                continue
            
            # Only accept entries ending with '-GGA' or '-GGA+U'
            is_pure_gga = entry_id.endswith('-GGA')
            is_gga_u = entry_id.endswith('-GGA+U')
            
            # Skip non-GGA entries (r2SCAN, SCAN, or no suffix)
            if not is_pure_gga and not is_gga_u:
                continue
            
            # Skip +U if pure_pbe requested
            if pure_pbe and is_gga_u:
                continue
            
            # Extract base MP ID (e.g., 'mp-540703' from 'mp-540703-GGA')
            parts = entry_id.split('-')
            if len(parts) >= 2:
                mp_id = parts[0] + '-' + parts[1]
            else:
                mp_id = entry_id
            
            # Get structure
            structure = None
            try:
                structure = mpr.get_structure_by_material_id(mp_id)
            except Exception as e:
                skipped_structure_retrieval.append(f"{mp_id} ({comp_entry.composition.reduced_formula})")
                continue
            
            if structure is None:
                continue
            
            mp_phases.append((mp_id, structure, is_gga_u))
            seen_entries[entry_id] = True
        
        if skipped_structure_retrieval:
            print(f"    WARNING: Could not retrieve structures for {len(skipped_structure_retrieval)} phases")
                    
        print(f"    Filtered to {len(mp_phases)} GGA phases with structures (strict '-GGA'/'-GGA+U' suffix)")
        
        # Verify terminal phases
        elements_found = set()
        for mp_id, structure, is_gga_u in mp_phases:
            if len(structure.composition.elements) == 1:
                elements_found.add(str(structure.composition.elements[0]))
        
        expected_elements = set(elements)
        if elements_found != expected_elements:
            missing = expected_elements - elements_found
            print(f"    WARNING: Missing terminal phases for elements: {sorted(missing)}")
        else:
            print(f"      All terminal phases present: {sorted(elements_found)}")
                    
        print(f"    Relaxing {len(mp_phases)} phases with MatterSim (fmax={fmax}, max_steps={max_steps})...")
        
        new_entries_data = []
        success_count = 0
        failed_count = 0
        
        for mp_id, structure, is_gga_u in mp_phases:
            # Determine chemsys for this phase
            doc_elements = sorted([str(el) for el in structure.composition.elements])
            doc_chemsys = '-'.join(doc_elements)
            
            # Skip if already cached
            entry_id_cached = f"mp_mattersim_{mp_id}"
            cache_key = (doc_chemsys, entry_id_cached)
            if cache_key in cached_data:
                continue
            
            # Relax with MatterSim using tight convergence
            relaxed_struct, energy_per_atom, error_msg = relax_structure_mattersim(
                structure, calculator, structure_id=mp_id, fmax=fmax, max_steps=max_steps
            )
            
            if relaxed_struct is None or energy_per_atom is None:
                print(f"      Warning: Failed to relax {mp_id}: {error_msg}")
                failed_count += 1
                continue
            
            total_energy = energy_per_atom * relaxed_struct.composition.num_atoms
            
            entry = PDEntry(
                composition=relaxed_struct.composition,
                energy=total_energy,
                name=f"mp_mattersim_{mp_id}"
            )
            entries.append(entry)
            
            structure_path = None
            try:
                if mp_relaxed_dir is not None:
                    structure_path = save_mp_relaxed_structure(relaxed_struct, mp_relaxed_dir, mp_id)
            except Exception as e:
                print(f"      Warning: Failed to save relaxed MP structure {mp_id}: {e}")

            new_entries_data.append({
                'chemsys': doc_chemsys,
                'composition': {str(el): float(amt) for el, amt in relaxed_struct.composition.items()},
                'energy': total_energy,
                'entry_id': f"mp_mattersim_{mp_id}",
                'mp_id': mp_id,
                'structure_path': structure_path,
                'structure_format': 'cif' if structure_path else None
            })
            success_count += 1
        
        print(f"    Successfully relaxed {success_count} GGA phases")
        if failed_count > 0:
            print(f"    Failed to relax {failed_count} phases")
        
        return entries, new_entries_data
                
    except Exception as e:
        print(f"    Error querying MP for {chemsys}: {e}")
        import traceback
        traceback.print_exc()
        return entries, []


def compute_energy_above_hull(structure, energy_per_atom, mp_entries):
    """
    Compute energy above hull using PDEntry.
    
    Args:
        structure: Pymatgen Structure
        energy_per_atom: Energy per atom (eV/atom)
        mp_entries: List of PDEntry objects for reference phases
    
    Returns:
        float: Energy above hull (eV/atom)
    """
    composition = structure.composition
    total_energy = energy_per_atom * composition.num_atoms
    
    entry = PDEntry(
        composition=composition,
        energy=total_energy,
        name='generated'
    )
    
    pd = PhaseDiagram(mp_entries)
    decomp, e_above_hull = pd.get_decomp_and_e_above_hull(entry, allow_negative=True)
    
    return float(e_above_hull)


def main():
    parser = argparse.ArgumentParser(
        description="Compute MatterSim energy_above_hull for origin-flow electride candidates only"
    )
    parser.add_argument(
        '--vasp-jobs',
        type=str,
        default='./VASP-out-Boron',
        help="Origin VASP jobs directory (default: ./VASP-out-Boron)"
    )
    parser.add_argument(
        '--candidate-csv',
        type=str,
        default='electride_analysis_candidates-Strict.csv',
        help="Candidate CSV used to restrict the run to shortlisted electrides only. "
             "Default: electride_analysis_candidates-Strict.csv in the VASP jobs directory."
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="Workflow database (default: workflow.json in vasp-jobs)"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key (default: from MP_API_KEY env)"
    )
    parser.add_argument(
        '--device',
        type=str,
        default='cpu',
        choices=['cpu', 'cuda'],
        help="Device for MatterSim (default: cpu)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='mattersim_stability_candidates.json',
        help="Output JSON file (default: mattersim_stability_candidates.json)"
    )
    parser.add_argument(
        '--pure-pbe',
        action='store_true',
        help="Filter MP entries to pure GGA-PBE only (exclude PBE+U). "
             "Default: accept both PBE and PBE+U for accurate phase diagrams. "
             "Use this flag to match DFT calculations using pure PBE without +U corrections."
    )
    
    args = parser.parse_args()

    if not MATTERSIM_AVAILABLE:
        print("ERROR: MatterSim not available in the current Python environment")
        print("Run this script from the mattersim conda environment")
        return 1
    
    vasp_jobs = Path(args.vasp_jobs).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = vasp_jobs / args.db
    
    output_path = vasp_jobs / args.output
    candidate_csv_path = resolve_path_with_fallbacks(args.candidate_csv, base_dir=vasp_jobs)
    mp_cache_local = vasp_jobs / "mp_mattersim_candidates.json"
    candidate_relaxed_dir = vasp_jobs / "mattersim_relaxed_candidates"
    mp_relaxed_dir = vasp_jobs / "mp_mattersim_relaxed"
    
    # Get MP API key
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: Materials Project API key required!")
        print("  Use: --mp-api-key YOUR_KEY or export MP_API_KEY=YOUR_KEY")
        return 1
    
    print("="*70)
    print("MatterSim Energy Above Hull Calculation (Origin Electride Candidates)")
    print("="*70)
    print(f"Origin VASP jobs: {vasp_jobs}")
    print(f"Workflow database: {db_path}")
    print(f"Candidate CSV: {candidate_csv_path}")
    print(f"MP MatterSim cache: {mp_cache_local}")
    print(f"Candidate relaxed structures: {candidate_relaxed_dir}")
    print(f"MP relaxed reference structures: {mp_relaxed_dir}")
    print(f"Output file: {output_path}")
    print(f"Device: {args.device}")
    print(f"Convergence: fmax=0.001 eV/Å, max_steps=800 (tighter than prescreen)")
    
    if args.pure_pbe:
        print(f"Functional filtering: Pure GGA-PBE only (PBE+U/R2SCAN/SCAN excluded)")
    else:
        print(f"Functional filtering: Mixed PBE/PBE+U (MP recommended methodology)")
    
    print("="*70 + "\n")
    
    # Load workflow database
    print("Loading workflow database...")
    if not db_path.exists():
        print(f"ERROR: Workflow database not found: {db_path}")
        return 1
    
    db = load_workflow_database(db_path)
    print(f"Found {len(db['structures'])} structures in database")

    if not candidate_csv_path.exists():
        print(f"ERROR: Candidate CSV not found: {candidate_csv_path}")
        return 1

    print(f"Loading candidate filter: {candidate_csv_path}")
    try:
        candidate_structure_ids = set(load_candidate_structure_ids(candidate_csv_path))
    except Exception as e:
        print(f"ERROR: Could not load candidate CSV: {e}")
        return 1

    db_structure_ids = set(db['structures'].keys())
    missing_candidate_ids = sorted(candidate_structure_ids - db_structure_ids)
    candidate_structure_ids &= db_structure_ids

    if not candidate_structure_ids:
        print("ERROR: No candidate structure IDs from the CSV were found in workflow.json")
        return 1

    print(f"Candidate IDs requested: {len(candidate_structure_ids) + len(missing_candidate_ids)}")
    print(f"Candidate IDs found in workflow database: {len(candidate_structure_ids)}")
    if missing_candidate_ids:
        print(f"Missing candidate IDs in workflow database: {len(missing_candidate_ids)}")
    print()
    
    # Scan for completed origin states compatible with using Relax/CONTCAR as input.
    print("Scanning for completed candidate structures...")
    completed_structures = []
    completed_states = {
        'RELAX_DONE', 'RELAX_TMOUT', 'SC_RUNNING', 'SC_DONE',
        'PARCHG_RUNNING', 'PARCHG_DONE', 'PARCHG_FAILED', 'PARCHG_SKIPPED',
        'ELF_RUNNING', 'ELF_DONE', 'ELF_FAILED'
    }
    
    for struct_id, sdata in db['structures'].items():
        if struct_id not in candidate_structure_ids:
            continue
        if sdata['state'] in completed_states:
            completed_structures.append(struct_id)
    
    skipped_not_completed = len(candidate_structure_ids) - len(completed_structures)
    print(f"Found {len(completed_structures)} completed candidate structures")
    if skipped_not_completed:
        print(f"Skipped {skipped_not_completed} candidate structures that are not in a completed state")
    print()
    
    if not completed_structures:
        print("No completed candidate structures found. Run the origin workflow first.")
        return 0
    
    # Group by chemical system
    print("Grouping structures by chemical system...")
    structures_by_chemsys = defaultdict(list)
    
    for struct_id in completed_structures:
        sdata = db['structures'][struct_id]
        chemsys = sdata.get('chemsys')
        
        if not chemsys:
            comp = Composition(sdata['composition'])
            elements = sorted([str(el) for el in comp.elements])
            chemsys = '-'.join(elements)
        
        structures_by_chemsys[chemsys].append(struct_id)
    
    required_chemsys = set(structures_by_chemsys.keys())
    unique_chemsys = sorted(required_chemsys)
    print(f"Found {len(unique_chemsys)} unique chemical systems")
    print(f"  {unique_chemsys}\n")
    
    # Create MatterSim calculator (for both MP phases and structures)
    print("Creating MatterSimCalculator...")
    try:
        import torch
        import gc
        calc = MatterSimCalculator(
            load_path="MatterSim-v1.0.0-5M.pth",
            device=args.device
        )
        print("  Calculator created successfully\n")
    except Exception as e:
        print(f"ERROR: Failed to create MatterSimCalculator: {e}")
        return 1
    
    # Load existing MP cache
    print("Loading existing MP cache...")
    cached_data = load_mp_cache(mp_cache_local)
    print(f"  Loaded {len(cached_data)} cached entries\n")
    
    # Pre-populate MP cache with all required chemical systems
    # This ensures complete GGA coverage for all subsystems
    print("="*70)
    print("Pre-populating MP Cache for Complete GGA Coverage")
    print("="*70)
    print(f"Chemical systems: {len(unique_chemsys)}")
    print(f"Strategy: Query ALL subsystems to match DFT cache completeness")
    print("="*70 + "\n")
    
    # Build list of all required subsystems (including elementals, binaries, etc.)
    all_required_chemsys = set()
    for chemsys in unique_chemsys:
        elements = chemsys.split('-')
        # Add all subsystems
        for n in range(1, len(elements) + 1):
            from itertools import combinations
            for combo in combinations(elements, n):
                all_required_chemsys.add('-'.join(sorted(combo)))
    
    print(f"Total unique subsystems needed: {len(all_required_chemsys)}")
    print(f"  (including elementals, binaries, ternaries, etc.)")
    print()
    
    # Pre-populate cache by querying all subsystems
    for chemsys in sorted(all_required_chemsys):
        # Check if this chemsys is already complete in cache
        cached_count = 0
        missing_saved_structures = 0
        for key, item in cached_data.items():
            item_chemsys = key[0]
            if item_chemsys == chemsys:
                structure_path = item.get('structure_path')
                if structure_path and Path(structure_path).exists():
                    cached_count += 1
                else:
                    missing_saved_structures += 1
        
        if cached_count > 0 and missing_saved_structures == 0:
            print(f"  {chemsys}: Using {cached_count} cached phases")
            continue
        if cached_count > 0 or missing_saved_structures > 0:
            print(
                f"  {chemsys}: Rebuilding cached phases to save missing MP structure files "
                f"(ready={cached_count}, missing_files={missing_saved_structures})"
            )
        
        print(f"  {chemsys}: Querying and relaxing MP phases...")
        sys.stdout.flush()
        
        try:
            # Force query by using empty cache dict
            _, new_cache_entries = get_mp_stable_phases_mattersim(
                chemsys,
                mp_api_key,
                {},
                calc,
                mp_relaxed_dir=mp_relaxed_dir,
                fmax=0.001,
                max_steps=800,
                pure_pbe=args.pure_pbe,
            )
            
            # Add to cache
            for new_entry in new_cache_entries:
                key = (new_entry['chemsys'], new_entry['entry_id'])
                cached_data[key] = new_entry
            
            print(f"    → Added {len(new_cache_entries)} phases to cache")
            
        except Exception as e:
            print(f"    → Error: {e}")
        
        sys.stdout.flush()
    
    print()
    print("="*70)
    print("MP cache pre-population complete")
    print(f"Total cached entries: {len(cached_data)}")
    print("="*70 + "\n")
    
    # Save updated cache after pre-population
    if cached_data:
        print("Saving pre-populated MP cache...")
        all_cache_list = list(cached_data.values())
        save_mp_cache(mp_cache_local, all_cache_list)
        print()
    
    # Now fetch MP entries for actual structure processing
    print("="*70)
    print("Loading MP GGA Phases for Structure Processing")
    print("="*70)
    print(f"Chemical systems: {len(unique_chemsys)}")
    print("="*70 + "\n")
    
    mp_entries_cache = {}
    
    for chemsys in unique_chemsys:
        print(f"Loading MP phases for {chemsys}...")
        sys.stdout.flush()
        try:
            mp_entries, _ = get_mp_stable_phases_mattersim(
                chemsys,
                mp_api_key,
                cached_data,
                calc,
                mp_relaxed_dir=mp_relaxed_dir,
                fmax=0.001,
                max_steps=800,
                pure_pbe=args.pure_pbe,
            )
            mp_entries_cache[chemsys] = mp_entries
            
            print(f"  → {len(mp_entries)} stable phases ready\n")
            sys.stdout.flush()
        except Exception as e:
            print(f"  → Error: {e}\n")
            sys.stdout.flush()
            mp_entries_cache[chemsys] = []
    
    print("="*70)
    print("MP GGA phases loading complete")
    print("="*70 + "\n")
    
    # Process structures
    print("="*70)
    print("Relaxing Origin Candidate Structures with MatterSim")
    print("="*70)
    print(f"Processing {len(completed_structures)} structures...\n")
    
    results = []
    processed = 0
    failed = 0
    
    for chemsys, struct_ids in sorted(structures_by_chemsys.items()):
        print(f"\nProcessing {chemsys} ({len(struct_ids)} structures)...")
        
        mp_entries = mp_entries_cache.get(chemsys, [])
        print(f"  Using {len(mp_entries)} MP reference phases\n")
        
        for struct_id in struct_ids:
            sdata = db['structures'][struct_id]
            relax_dir = Path(sdata['relax_dir'])
            contcar_path = relax_dir / 'CONTCAR'
            
            print(f"  {struct_id}:")
            
            # Load VASP-relaxed structure
            if not contcar_path.exists():
                print(f"    FAILED: CONTCAR not found")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'input_structure_path': str(contcar_path.resolve()),
                    'relaxed_structure_path': None,
                    'relaxed_structure_format': None,
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'selected_candidate': True,
                    'error': 'CONTCAR not found'
                })
                continue
            
            try:
                structure = Structure.from_file(str(contcar_path))
            except Exception as e:
                print(f"    FAILED: Could not load CONTCAR: {e}")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'input_structure_path': str(contcar_path.resolve()),
                    'relaxed_structure_path': None,
                    'relaxed_structure_format': None,
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'selected_candidate': True,
                    'error': f'CONTCAR parsing failed: {e}'
                })
                continue
            
            # Relax with MatterSim (tighter convergence)
            relaxed_struct, energy_per_atom, error_msg = relax_structure_mattersim(
                structure, calc, structure_id=struct_id, fmax=0.001, max_steps=800
            )
            
            if relaxed_struct is None or energy_per_atom is None:
                print(f"    FAILED: {error_msg}")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'input_structure_path': str(contcar_path.resolve()),
                    'relaxed_structure_path': None,
                    'relaxed_structure_format': None,
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'selected_candidate': True,
                    'error': error_msg
                })
                continue
            
            print(f"    MatterSim energy: {energy_per_atom:.6f} eV/atom")

            relaxed_structure_path = None
            try:
                relaxed_structure_path = save_candidate_relaxed_structure(
                    relaxed_struct, candidate_relaxed_dir, struct_id
                )
                print(f"    Saved relaxed structure: {relaxed_structure_path}")
            except Exception as e:
                print(f"    WARNING: Failed to save relaxed structure: {e}")
            
            # Compute hull
            if not mp_entries:
                e_hull = None
                print(f"    WARNING: No MP reference phases, skipping hull calculation")
            else:
                try:
                    e_hull = compute_energy_above_hull(relaxed_struct, energy_per_atom, mp_entries)
                    print(f"    E_hull: {e_hull:.6f} eV/atom")
                except Exception as e:
                    print(f"    FAILED: Hull calculation error: {e}")
                    e_hull = None
            
            is_stable = bool(e_hull < 0.001) if e_hull is not None else None
            processed += 1
            
            results.append({
                'structure_id': struct_id,
                'composition': sdata['composition'],
                'chemsys': chemsys,
                'input_structure_path': str(contcar_path.resolve()),
                'relaxed_structure_path': relaxed_structure_path,
                'relaxed_structure_format': 'poscar' if relaxed_structure_path else None,
                'mattersim_energy_per_atom': float(energy_per_atom),
                'energy_above_hull': float(e_hull) if e_hull is not None else None,
                'is_stable': is_stable,
                'selected_candidate': True,
                'error': None
            })
    
    # Clean up calculator
    del calc
    if args.device == 'cuda':
        import torch
        torch.cuda.empty_cache()
        import gc
        gc.collect()
    
    # Save results
    print("\n" + "="*70)
    print("Saving MatterSim hull results...")
    
    output_data = {
        'summary': {
            'total_structures': len(completed_structures),
            'candidate_csv': str(candidate_csv_path),
            'candidate_ids_requested': len(candidate_structure_ids) + len(missing_candidate_ids),
            'candidate_ids_found_in_workflow': len(candidate_structure_ids),
            'missing_candidate_ids_in_workflow': missing_candidate_ids,
            'successfully_processed': processed,
            'failed': failed,
            'energy_reference': 'MatterSim-v1.0.0-5M',
            'convergence': 'fmax=0.001, max_steps=800',
            'mp_reference_source': 'MP API GGA phases relaxed with MatterSim (legacy pymatgen MPRester)',
            'mp_functional_filter': 'pure_pbe' if args.pure_pbe else 'mixed_pbe_pbeU',
            'mp_cache': str(mp_cache_local),
            'candidate_relaxed_structures_dir': str(candidate_relaxed_dir.resolve()),
            'mp_relaxed_structures_dir': str(mp_relaxed_dir.resolve())
        },
        'results': results
    }
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"MatterSim results saved to: {output_path}")
    print("="*70)
    print(f"\nMatterSim Hull Summary:")
    print(f"  Total structures: {len(completed_structures)}")
    print(f"  Successfully processed: {processed}")
    print(f"  Failed: {failed}")
    print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    exit(main())
