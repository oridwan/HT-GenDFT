#!/usr/bin/env python3
"""
MatterSim Pre-screening for VASPflow

Performs fast thermodynamic stability screening using MatterSim before expensive VASP calculations.
Outputs: prescreening_stability.json with structures passing energy_above_hull threshold.

IMPORTANT: Uses Legacy MPRester for Complete GGA Coverage
- Uses pymatgen.ext.matproj.MPRester (legacy API, not mp_api.client)
- The modern mp_api.client misses many stable GGA phases needed for accurate hull calculations
- Queries ALL GGA/GGA+U entries via get_entries_in_chemsys() with strict suffix filtering
- Strict filtering: only accepts entries with '-GGA' or '-GGA+U' suffix in entry_id
- MP structures are re-relaxed with MatterSim to obtain consistent energy reference
- Optional --pure-pbe flag to exclude PBE+U (use pure GGA-PBE only)

GPU Batch Processing:
- Uses MatterSim's BatchRelaxer for true parallel GPU relaxation
- Multiple structures are relaxed simultaneously on GPU
- --batch-size controls GPU parallelism (e.g., 32 structures at once)
- --max-atoms-gpu sets memory limit (default 2048)

Duplicate Detection:
- Removes duplicate structures within each composition BEFORE relaxation
- Uses pymatgen StructureMatcher (ltol=0.2, stol=0.2, angle_tol=5)
- Significantly reduces computational cost by avoiding redundant relaxations
- Only compares structures with same composition (not across compositions)
"""

import os
import sys
import json
import argparse
import zipfile
import warnings
import fcntl
import inspect
import numpy as np
from pathlib import Path
from io import StringIO
from datetime import datetime

from pymatgen.core import Composition
from pymatgen.io.cif import CifParser
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry

try:
    from pymatgen.ext.matproj import MPRester
except ImportError:
    print("ERROR: pymatgen package with MPRester is required")
    print("Install with: pip install pymatgen")
    sys.exit(1)
from pymatgen.io.ase import AseAtomsAdaptor

from pyxtal import pyxtal
from pyxtal.db import database_topology

# MatterSim 1.2.0 expects ase.constraints.Filter, but newer ASE moved Filter to ase.filters.
# Add a small compatibility alias before importing MatterSim.
try:
    import ase.constraints as _ase_constraints
    from ase.filters import Filter as _ase_filter
    if not hasattr(_ase_constraints, "Filter"):
        _ase_constraints.Filter = _ase_filter
except Exception:
    pass

# ASE optimizer API compatibility:
# Some ASE versions require Optimizer.converged(self, gradient), while MatterSim can call
# converged() without arguments through its optimizer wrapper. Patch signature to accept
# both call styles.
try:
    from ase.optimize.optimize import Optimizer as _ase_optimizer
    _converged_sig = inspect.signature(_ase_optimizer.converged)
    _params = list(_converged_sig.parameters.values())
    if len(_params) >= 2 and _params[1].name == "gradient" and _params[1].default is inspect._empty:
        _orig_converged = _ase_optimizer.converged

        def _converged_compat(self, gradient=None):
            if gradient is None:
                try:
                    gradient = self.optimizable.get_gradient()
                except Exception:
                    gradient = None
            return _orig_converged(self, gradient)

        _ase_optimizer.converged = _converged_compat
except Exception:
    pass

try:
    from mattersim.forcefield.potential import Potential
    from mattersim.applications.batch_relax import BatchRelaxer
    MATTERSIM_AVAILABLE = True
except Exception as e:
    MATTERSIM_AVAILABLE = False
    print(f"ERROR: MatterSim not available! ({type(e).__name__}: {e})")
    print("Hint: MatterSim 1.2.0 can require older ASE API symbols.")
    sys.exit(1)

warnings.filterwarnings('ignore', category=UserWarning, message='.*POTCAR data with symbol.*')
warnings.filterwarnings('ignore', message='Using UFloat objects with std_dev==0')
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in det')


def validate_structure(pmg_struct):
    """
    Pre-validate structure to catch obviously malformed structures before MatterSim.
    Returns (is_valid, error_message).
    """
    try:
        # Check 1: Reasonable lattice parameters (not too small/large)
        lattice = pmg_struct.lattice
        abc = lattice.abc
        if any(a < 1.0 or a > 50.0 for a in abc):
            return False, f"Unrealistic lattice parameters: {abc}"
        
        # Check 2: Reasonable angles
        angles = lattice.angles
        if any(angle < 10.0 or angle > 170.0 for angle in angles):
            return False, f"Unrealistic angles: {angles}"
        
        # Check 3: Atoms too close (min distance check)
        dist_matrix = pmg_struct.distance_matrix
        min_dist = dist_matrix[dist_matrix > 0].min() if len(dist_matrix) > 1 else 1.0
        if min_dist < 0.5:  # Too close (< 0.5 Angstrom)
            return False, f"Atoms too close: min_dist={min_dist:.3f} Å"
        
        # Check 4: Try spglib symmetry analysis (catches spglib failures early)
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        try:
            sga = SpacegroupAnalyzer(pmg_struct, symprec=0.1)
            _ = sga.get_space_group_symbol()
        except Exception as e:
            # If spglib fails here, it will definitely fail later
            return False, f"Spglib symmetry detection failed: {str(e)[:100]}"
        
        return True, None
        
    except Exception as e:
        return False, f"Validation error: {str(e)[:100]}"


def relax_structure_mattersim(structures_dict, potential, fmax=0.01, max_steps=500, max_natoms_per_batch=2048):
    """
    Batch relax multiple structures using MatterSim BatchRelaxer for efficient GPU utilization.
    
    Args:
        structures_dict: List of dicts with keys: 'id', 'structure' (pymatgen), 'composition', 'chemsys'
        potential: MatterSim Potential instance (reused for entire batch)
        fmax: Force convergence criterion (eV/Angstrom)
        max_steps: Maximum optimization steps per structure
        max_natoms_per_batch: Maximum total atoms in GPU batch (controls memory usage)
    
    Returns:
        list of dicts: Results for each structure with keys:
            - 'structure_id': Structure ID
            - 'relaxed_structure': Pymatgen Structure (None if failed)
            - 'energy_per_atom': Energy per atom (None if failed)
            - 'error': Error message (None if successful)
            - 'composition': Original composition
            - 'chemsys': Chemical system
    
    Note:
        - Symmetrizes each structure using PyXtal before relaxation
        - Uses BatchRelaxer for true parallel GPU processing
        - Falls back to direct conversion if PyXtal fails
        - Has extensive error handling for validation and relaxation failures
    """
    
    results = []
    atoms_list = []
    valid_indices = []  # Track which structures successfully converted to ASE
    
    adaptor = AseAtomsAdaptor()
    tolerances = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
    
    # Convert all structures to ASE atoms with validation
    for idx, item in enumerate(structures_dict):
        struct_id = item['id']
        pmg_struct = item['structure']
        
        # Try PyXtal symmetrization with progressive tolerances
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
        
        # If all tolerances failed, try direct conversion
        if atoms is None:
            try:
                atoms = adaptor.get_atoms(pmg_struct)
            except Exception as e:
                results.append({
                    'structure_id': struct_id,
                    'relaxed_structure': None,
                    'energy_per_atom': None,
                    'error': f"Structure conversion failed (PyXtal failed at all tolerances {tolerances}, direct conversion also failed): {e}",
                    'composition': item['composition'],
                    'chemsys': item['chemsys']
                })
                continue
        
        # Store original index for mapping results back
        atoms.info['structure_index'] = idx
        atoms.info['structure_id'] = struct_id
        atoms_list.append(atoms)
        valid_indices.append(idx)
    
    # If all structures failed conversion, return early
    if not atoms_list:
        return results
    
    actual_atom_counts = [len(atoms) for atoms in atoms_list]
    max_single_structure = max(actual_atom_counts)
    total_actual_atoms = sum(actual_atom_counts)
    avg_actual_atoms = total_actual_atoms / len(atoms_list)
    
    if max_natoms_per_batch < max_single_structure:
        print(f"    Adjusting max_natoms_per_batch: {max_natoms_per_batch} -> {max_single_structure} "
              f"(largest structure has {max_single_structure} atoms after symmetrization)")
        max_natoms_per_batch = max_single_structure
    
    print(f"    Actual atoms after symmetrization: total={total_actual_atoms}, "
          f"avg={avg_actual_atoms:.0f}, max={max_single_structure}, "
          f"max_natoms_per_batch={max_natoms_per_batch}")
    
    # Create BatchRelaxer and relax all structures
    try:
        relaxer = BatchRelaxer(
            potential=potential,
            optimizer="FIRE",
            filter="ExpCellFilter",
            fmax=fmax,
            max_natoms_per_batch=max_natoms_per_batch,
            max_n_steps=max_steps
        )
        
        trajectories = relaxer.relax(atoms_list)
        
        # Extract results from trajectories
        for traj_pos, orig_idx in enumerate(valid_indices):
            item = structures_dict[orig_idx]
            struct_id = item['id']
            
            if traj_pos in trajectories:
                trajectory = trajectories[traj_pos]
                relaxed_atoms = trajectory[-1]
                
                try:
                    energy = relaxed_atoms.info.get('total_energy', None)
                    
                    # Check for NaN/inf values
                    if energy is None or not np.isfinite(energy):
                        results.append({
                            'structure_id': struct_id,
                            'relaxed_structure': None,
                            'energy_per_atom': None,
                            'error': f"Relaxation produced non-finite energy: {energy}",
                            'composition': item['composition'],
                            'chemsys': item['chemsys']
                        })
                        continue
                    
                    energy_per_atom = energy / len(relaxed_atoms)
                    
                    relaxed_structure = adaptor.get_structure(relaxed_atoms)
                    
                    results.append({
                        'structure_id': struct_id,
                        'relaxed_structure': relaxed_structure,
                        'energy_per_atom': energy_per_atom,
                        'error': None,
                        'composition': item['composition'],
                        'chemsys': item['chemsys']
                    })
                    
                except Exception as e:
                    results.append({
                        'structure_id': struct_id,
                        'relaxed_structure': None,
                        'energy_per_atom': None,
                        'error': f"Post-relaxation processing failed: {e}",
                        'composition': item['composition'],
                        'chemsys': item['chemsys']
                    })
            else:
                # Structure not in trajectories (relaxation failed internally)
                results.append({
                    'structure_id': struct_id,
                    'relaxed_structure': None,
                    'energy_per_atom': None,
                    'error': "BatchRelaxer failed to produce trajectory",
                    'composition': item['composition'],
                    'chemsys': item['chemsys']
                })
    
    except Exception as e:
        # Batch relaxation completely failed
        for idx in valid_indices:
            item = structures_dict[idx]
            results.append({
                'structure_id': item['id'],
                'relaxed_structure': None,
                'energy_per_atom': None,
                'error': f"Batch relaxation failed: {e}",
                'composition': item['composition'],
                'chemsys': item['chemsys']
            })
    
    return results


def deduplicate_structures(structures_list):
    """
    Remove duplicate structures within each composition using StructureMatcher.
    
    Args:
        structures_list: List of dicts with keys: 'id', 'composition', 'chemsys', 'structure'
    
    Returns:
        tuple: (deduplicated_list, duplicate_count, duplicate_details)
    
    Note:
        - Only compares structures with the same composition
        - Uses reasonable tolerances: ltol=0.2, stol=0.2, angle_tol=5
        - Keeps the first structure when duplicates are found
        - Important to do this BEFORE expensive MatterSim relaxations
    """
    from pymatgen.analysis.structure_matcher import StructureMatcher
    
    # Group structures by composition
    by_composition = {}
    for item in structures_list:
        comp = item['composition']
        if comp not in by_composition:
            by_composition[comp] = []
        by_composition[comp].append(item)
    
    # Initialize matcher with reasonable tolerances
    matcher = StructureMatcher(ltol=0.2, stol=0.2, angle_tol=5)
    
    deduplicated = []
    duplicate_count = 0
    duplicate_details = []
    
    # Process each composition separately
    for comp, comp_structures in by_composition.items():
        if len(comp_structures) == 1:
            deduplicated.extend(comp_structures)
            continue
        
        # Track which structures to keep
        unique_structures = []
        duplicate_of = {}
        
        for i, item_i in enumerate(comp_structures):
            struct_i = item_i['structure']
            is_duplicate = False
            
            for j, item_j in enumerate(unique_structures):
                struct_j = item_j['structure']
                
                try:
                    if matcher.fit(struct_i, struct_j):
                        is_duplicate = True
                        duplicate_of[i] = j
                        duplicate_count += 1
                        duplicate_details.append({
                            'duplicate': item_i['id'],
                            'original': item_j['id'],
                            'composition': comp
                        })
                        break
                except Exception as e:
                    continue
            
            if not is_duplicate:
                unique_structures.append(item_i)
        
        deduplicated.extend(unique_structures)
    
    return deduplicated, duplicate_count, duplicate_details


def load_mp_cache_locked(cache_file):
    """
    Load MP cache with file locking for safe concurrent access.
    
    Returns:
        dict: Cache data keyed by (chemsys, entry_id)
    """
    cache_file = Path(cache_file)
    cached_data = {}
    
    if cache_file.exists():
        with open(cache_file, 'r') as f:
            # Acquire shared lock for reading
            fcntl.flock(f.fileno(), fcntl.LOCK_SH)
            try:
                cache_list = json.load(f)
                # Convert list to dict keyed by (chemsys, mp_id) for fast lookup
                for item in cache_list:
                    key = (item.get('chemsys', ''), item['entry_id'])
                    cached_data[key] = item
            finally:
                fcntl.flock(f.fileno(), fcntl.LOCK_UN)
    
    return cached_data


def save_mp_cache_locked(cache_file, new_entries_data):
    """
    Save MP cache entries with file locking for safe concurrent writes.
    
    Args:
        cache_file: Path to cache file
        new_entries_data: List of new cache entries to add
    """
    cache_file = Path(cache_file)
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Use a lock file to coordinate access
    lock_file = cache_file.with_suffix('.lock')
    
    with open(lock_file, 'w') as lock_f:
        # Acquire exclusive lock
        fcntl.flock(lock_f.fileno(), fcntl.LOCK_EX)
        try:
            # Load existing cache
            existing_cache = {}
            if cache_file.exists():
                with open(cache_file, 'r') as f:
                    cache_list = json.load(f)
                    for item in cache_list:
                        key = (item.get('chemsys', ''), item['entry_id'])
                        existing_cache[key] = item
            
            # Add new entries (avoid duplicates)
            for new_item in new_entries_data:
                key = (new_item['chemsys'], new_item['entry_id'])
                if key not in existing_cache:
                    existing_cache[key] = new_item
            
            # Write back to file
            all_cache_data = list(existing_cache.values())
            with open(cache_file, 'w') as f:
                json.dump(all_cache_data, f, indent=2)
            
            print(f"    Updated cache: {cache_file} ({len(all_cache_data)} total entries)")
        finally:
            fcntl.flock(lock_f.fileno(), fcntl.LOCK_UN)


def get_mp_stable_phases_mattersim(chemsys, mp_api_key, cache_file, potential, pure_pbe=False, max_atoms_gpu=2048):
    """
    Get MP GGA phase structures and relax with MatterSim for consistent energy reference.
    
    CRITICAL: Uses Legacy MPRester for Complete GGA Coverage
    - Uses pymatgen.ext.matproj.MPRester (legacy API, not mp_api.client)
    - The modern mp_api.client misses many stable GGA phases needed for accurate hull calculations
    - Queries ALL GGA/GGA+U entries via get_entries_in_chemsys()
    - Strict filtering: only accepts entries with '-GGA' or '-GGA+U' suffix in entry_id
    - MP structures are re-relaxed with MatterSim to obtain MatterSim energies
    - Optional pure_pbe flag to exclude PBE+U (use pure GGA-PBE only)
    
    Key optimizations:
    - Uses single SHARED global cache file (mp_mattersim.json) across all parallel batches
    - Per-chemsys locking prevents duplicate MP queries across batches
    - Uses BatchRelaxer for efficient parallel GPU relaxation
    - Returns PDEntry objects for phase diagram analysis
    
    Args:
        chemsys: Chemical system string (e.g., 'B-Li-N')
        mp_api_key: Materials Project API key
        cache_file: Path to SHARED global cache file (mp_mattersim.json)
        potential: MatterSim Potential instance (for batch relaxation)
        pure_pbe: If True, filter to GGA-PBE only (exclude PBE+U)
                  If False (default), accept both PBE and PBE+U
        max_atoms_gpu: Maximum total atoms on GPU simultaneously (default: 2048)
    
    Returns:
        List of PDEntry objects for GGA phases in this chemical system
    """
    cache_file = Path(cache_file)
    
    # Use per-chemical-system lock to prevent duplicate queries
    # This ensures only ONE batch queries MP for this chemsys at a time
    chemsys_lock_file = cache_file.parent / f"mp_cache_{chemsys.replace('-', '')}.lock"
    
    with open(chemsys_lock_file, 'w') as chemsys_lock:
        # Acquire exclusive lock for this chemical system
        fcntl.flock(chemsys_lock.fileno(), fcntl.LOCK_EX)
        
        try:
            # Double-check pattern: read cache AFTER acquiring lock
            # Another batch may have populated it while we waited for lock
            cached_data = load_mp_cache_locked(cache_file)
    
            # Extract entries for this chemsys from cache
            elements = chemsys.split('-')
            entries = []
            cached_count = 0
            
            for key, item in cached_data.items():
                item_chemsys = key[0]
                # Include if chemsys matches OR if it's a subset (e.g., 'B', 'Li' for 'B-Li-N')
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
                return entries
    
            # Need to query MP for missing data
            print(f"    Querying MP for GGA phases in {chemsys}...")
            if missing_chemsys:
                print(f"      Missing subsystems: {sorted(missing_chemsys)}")
    
            try:
                # Use legacy pymatgen MPRester which returns complete GGA entries
                mpr = MPRester(mp_api_key)
                
                # Get ALL entries in this chemical system
                # This automatically includes all subsystems (elementals, binaries, etc.)
                computed_entries = mpr.get_entries_in_chemsys(elements)
                
                print(f"    Retrieved {len(computed_entries)} entries from MP (filtering for GGA...)")
                
                # Process entries: filter for GGA only (entry_id ending with '-GGA' or '-GGA+U')
                mp_phases = []
                seen_entries = {}  # Track by entry_id to avoid duplicates
                skipped_structure_retrieval = []
                
                for comp_entry in computed_entries:
                    entry_id = str(comp_entry.entry_id)
                    
                    # Skip if already seen
                    if entry_id in seen_entries:
                        continue
                    
                    # Only accept entries ending with '-GGA' or '-GGA+U' (strict filtering)
                    is_pure_gga = entry_id.endswith('-GGA')
                    is_gga_u = entry_id.endswith('-GGA+U')
                    
                    # Skip non-GGA entries (r2SCAN, SCAN, or no suffix)
                    if not is_pure_gga and not is_gga_u:
                        continue
                    
                    # Skip +U if pure_pbe requested
                    if pure_pbe and is_gga_u:
                        continue
                    
                    has_U = is_gga_u
                    
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
                        skipped_structure_retrieval.append(
                            f"{mp_id} ({comp_entry.composition.reduced_formula}): {str(e).split(chr(10))[0][:120]}"
                        )
                        continue
                    
                    if structure is None:
                        skipped_structure_retrieval.append(
                            f"{mp_id} ({comp_entry.composition.reduced_formula}): returned None"
                        )
                        continue
                    
                    mp_phases.append((mp_id, structure, has_U))
                    seen_entries[entry_id] = True
                
                # Check if we have all elemental (terminal) phases
                elements_found = set()
                for mp_id, structure, has_U in mp_phases:
                    if len(structure.composition.elements) == 1:  # Pure elemental phase
                        elements_found.add(str(structure.composition.elements[0]))
                
                expected_elements = set(elements)
                missing_elements = expected_elements - elements_found
                
                # FALLBACK: If missing any elemental phases, use modern API through legacy wrapper
                if missing_elements and hasattr(mpr, 'materials'):
                    print(f"    WARNING: Missing terminal phases for {sorted(missing_elements)}, trying modern API fallback...")
                    try:
                        # Query each missing element individually via modern API
                        for elem in sorted(missing_elements):
                            elem_docs = mpr.materials.summary.search(
                                elements=[elem],
                                num_elements=(1, 1),  # Only elemental phases
                                fields=["material_id", "formula_pretty", "structure", "energy_per_atom"]
                            )
                            
                            if elem_docs:
                                print(f"      Fallback: Found {len(elem_docs)} {elem} elemental phases")
                                for doc in elem_docs:
                                    mp_id = doc.material_id
                                    structure = doc.structure
                                    
                                    if mp_id not in seen_entries and structure is not None:
                                        mp_phases.append((mp_id, structure, False))  # Assume no +U for fallback
                                        seen_entries[mp_id] = True
                            else:
                                print(f"      Fallback WARNING: No {elem} elemental phases found even in modern API!")
                        
                        print(f"    Modern API fallback: Now have {len(mp_phases)} total phases")
                    except Exception as e:
                        print(f"    Modern API fallback failed: {e}")
                
                if skipped_structure_retrieval:
                    print(f"    WARNING: Could not retrieve structures for {len(skipped_structure_retrieval)} phases:")
                    for skip_msg in skipped_structure_retrieval:
                        print(f"      - {skip_msg}")
                
                print(f"    Filtered to {len(mp_phases)} GGA phases (strict '-GGA'/'-GGA+U' suffix)")
                
                # Verify we have terminal (elemental) phases
                elements_found = set()
                all_elements_found = set()
                for mp_id, structure, has_U in mp_phases:
                    for el in structure.composition.elements:
                        all_elements_found.add(str(el))
                    if len(structure.composition.elements) == 1:  # Pure elemental phase
                        elements_found.add(str(structure.composition.elements[0]))
                
                expected_elements = set(elements)
                missing_terminal = expected_elements - elements_found
                missing_any = expected_elements - all_elements_found
                
                if missing_any:
                    print(f"    WARNING: NO phases at all for elements: {sorted(missing_any)}")
                    print(f"      Phase diagram will be incomplete - hull computation will be skipped")
                    print(f"      Possible causes: structure retrieval failed, or MP has no data")
                elif missing_terminal:
                    print(f"    WARNING: Missing terminal (elemental) phases for: {sorted(missing_terminal)}")
                    print(f"      These elements appear in compounds but not as pure phases")
                else:
                    print(f"      All terminal phases present: {sorted(elements_found)}")
                
                # Prepare structures for batch relaxation (skip already cached)
                structures_to_relax = []
                for mp_id, structure, has_U in mp_phases:
                    doc_elements = sorted([str(el) for el in structure.composition.elements])
                    doc_chemsys = '-'.join(doc_elements)
                    
                    # Skip if already cached
                    entry_id_cached = f"mp_mattersim_{mp_id}"
                    cache_key = (doc_chemsys, entry_id_cached)
                    if cache_key in cached_data:
                        continue
                    
                    structures_to_relax.append({
                        'id': mp_id,
                        'structure': structure,
                        'composition': doc_chemsys,  # Store for later use
                        'chemsys': doc_chemsys
                    })
                
                if not structures_to_relax:
                    print(f"    All {len(mp_phases)} phases already cached")
                    return entries
                
                print(f"    Batch relaxing {len(structures_to_relax)} uncached phases with MatterSim...")
                
                batch_results = relax_structure_mattersim(
                    structures_to_relax,
                    potential,
                    fmax=0.01,
                    max_steps=500,
                    max_natoms_per_batch=max_atoms_gpu
                )
                
                # Process batch results
                new_entries_data = []
                success_count = 0
                failed_count = 0
                
                for result in batch_results:
                    mp_id = result['structure_id']
                    relaxed_struct = result['relaxed_structure']
                    energy_per_atom = result['energy_per_atom']
                    error_msg = result['error']
                    doc_chemsys = result['chemsys']
                    
                    if relaxed_struct is None or energy_per_atom is None:
                        print(f"      Warning: Failed to relax {mp_id}: {error_msg}")
                        failed_count += 1
                        continue
                    
                    total_energy = energy_per_atom * relaxed_struct.composition.num_atoms
                    
                    # Create PDEntry
                    entry = PDEntry(
                        composition=relaxed_struct.composition,
                        energy=total_energy,
                        name=f"mp_mattersim_{mp_id}"
                    )
                    entries.append(entry)
                    
                    # Store for caching
                    new_entries_data.append({
                        'chemsys': doc_chemsys,
                        'composition': {str(el): float(amt) for el, amt in relaxed_struct.composition.items()},
                        'energy': total_energy,
                        'entry_id': f"mp_mattersim_{mp_id}",
                        'mp_id': mp_id
                    })
                    success_count += 1
                
                print(f"    Successfully relaxed {success_count}/{len(structures_to_relax)} GGA phases")
                if failed_count > 0:
                    print(f"    Failed to relax {failed_count}/{len(structures_to_relax)} phases")
                
                # Update cache with new entries (with file locking)
                if new_entries_data:
                    save_mp_cache_locked(cache_file, new_entries_data)
                
                return entries
                
            except Exception as e:
                print(f"    Error querying MP for {chemsys}: {e}")
                import traceback
                traceback.print_exc()
                return entries
        
        finally:
            # Release the per-chemsys lock
            fcntl.flock(chemsys_lock.fileno(), fcntl.LOCK_UN)


def compute_energy_above_hull(structure, energy_per_atom, mp_entries):
    """
    Compute energy above hull for a structure using PDEntry.
    
    Args:
        structure: Pymatgen Structure
        energy_per_atom: Energy per atom (eV/atom)
        mp_entries: List of PDEntry objects for reference phases
    
    Returns:
        float: Energy above hull (eV/atom), or None if phase diagram is incomplete
    
    Raises:
        ValueError: If mp_entries is empty or phase diagram construction fails
    """
    composition = structure.composition
    total_energy = energy_per_atom * composition.num_atoms
    
    # Check that mp_entries covers all elements in the structure
    pd_elements = set()
    for mp_entry in mp_entries:
        for el in mp_entry.composition.elements:
            pd_elements.add(str(el))
    
    struct_elements = set(str(el) for el in composition.elements)
    missing_elements = struct_elements - pd_elements
    
    if missing_elements:
        print(f"    WARNING: MP reference phases missing elements {sorted(missing_elements)} "
              f"for composition {composition.reduced_formula}")
        print(f"      Phase diagram elements: {sorted(pd_elements)}")
        print(f"      Structure elements: {sorted(struct_elements)}")
        print(f"      Cannot compute energy above hull - returning None")
        return None
    
    entry = PDEntry(
        composition=composition,
        energy=total_energy,
        name='generated'
    )
    
    pd = PhaseDiagram(mp_entries)
    decomp, e_above_hull = pd.get_decomp_and_e_above_hull(entry, allow_negative=True)
    
    return float(e_above_hull)


def read_structures_from_zip(zip_path, max_structures=None):
    """Read CIF structures from zip file."""
    structures = []
    with zipfile.ZipFile(zip_path, 'r') as zf:
        cif_files = sorted([f for f in zf.namelist() if f.endswith('.cif')])
        if max_structures:
            cif_files = cif_files[:max_structures]
        
        for cif_file in cif_files:
            try:
                with zf.open(cif_file) as f:
                    cif_content = f.read().decode('utf-8')
                    parser = CifParser(StringIO(cif_content))
                    structure = parser.parse_structures(primitive=True)[0]
                    structures.append(structure)
            except Exception as e:
                print(f"  Warning: Could not parse {cif_file}: {e}")
    return structures


def main():
    parser = argparse.ArgumentParser(
        description="MatterSim Pre-screening for VASPflow"
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help="MatterGen results directory"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key"
    )
    parser.add_argument(
        '--hull-threshold',
        type=float,
        default=0.1,
        help="Energy above hull threshold (eV/atom)"
    )
    parser.add_argument(
        '--device',
        type=str,
        default='cpu',
        choices=['cpu', 'cuda'],
        help="Device for MatterSim"
    )
    parser.add_argument(
        '--start-composition',
        type=int,
        default=0,
        help="Starting composition index for parallel batch processing (default: 0)"
    )
    parser.add_argument(
        '--max-compositions',
        type=int,
        default=None,
        help="Max compositions to process (count from start-composition)"
    )
    parser.add_argument(
        '--max-structures',
        type=int,
        default=0,
        help="Max structures per composition"
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=32,
        help="Batch size for GPU parallel relaxation (structures relaxed simultaneously)"
    )
    parser.add_argument(
        '--max-atoms-gpu',
        type=int,
        default=2048,
        help="Maximum total atoms on GPU simultaneously (default: 2048 for V100). "
             "Increase for GPUs with more VRAM: 4096 for A100"
    )
    parser.add_argument(
        '--mp-cache-dir',
        type=str,
        default=None,
        help="Directory to save MP cache file (default: same as output-dir)"
    )
    parser.add_argument(
        '--batch-id',
        type=str,
        default=None,
        help="Batch ID for parallel runs (creates unique checkpoint/output files)"
    )
    parser.add_argument(
        '--pure-pbe',
        action='store_true',
        help="Filter MP entries to pure GGA-PBE only (exclude PBE+U). "
             "Default: accept both PBE and PBE+U for accurate phase diagrams. "
             "Use this flag to match DFT calculations using pure PBE without +U corrections."
    )
    
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: Materials Project API key required!")
        print("  Use: --mp-api-key YOUR_KEY or export MP_API_KEY=YOUR_KEY")
        return 1
    
    # Use SHARED global cache file (with file locking for parallel safety)
    # All parallel batches share the same cache to avoid duplicate MP API calls
    mp_cache_dir = Path(args.mp_cache_dir).expanduser() if args.mp_cache_dir else output_dir
    mp_cache_file = mp_cache_dir / "mp_mattersim.json"
    
    print("="*70)
    print("MatterSim Pre-screening for VASPflow")
    print("="*70)
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Hull threshold: {args.hull_threshold} eV/atom")
    print(f"Device: {args.device}")
    print(f"Batch size: {args.batch_size}")
    print(f"Start composition index: {args.start_composition}")
    print(f"Max compositions: {args.max_compositions or 'all'}")
    print(f"Max structures per composition: {args.max_structures}")
    print(f"MP cache file: {mp_cache_file}")
    print(f"MP API: Legacy pymatgen.ext.matproj.MPRester (complete GGA coverage)")
    
    if args.pure_pbe:
        print(f"Functional filtering: Pure GGA-PBE only (PBE+U/R2SCAN/SCAN excluded)")
    else:
        print(f"Functional filtering: Mixed PBE/PBE+U (MP recommended methodology)")
    
    print("="*70 + "\n")
    
    
    # Scan structures
    comp_dirs = sorted(results_dir.glob("*_structures"))
    
    # Apply start index for parallel batch processing
    if args.start_composition > 0:
        comp_dirs = comp_dirs[args.start_composition:]
    
    # Apply max compositions limit
    if args.max_compositions:
        comp_dirs = comp_dirs[:args.max_compositions]
    
    all_structures = []
    
    for comp_dir in comp_dirs:
        comp_name = comp_dir.name.replace("_structures", "")
        zip_path = comp_dir / "generated_crystals_cif.zip"
        
        if not zip_path.exists():
            continue
        
        print(f"Loading {comp_name}...")
        structures = read_structures_from_zip(zip_path, args.max_structures)
        
        for idx, structure in enumerate(structures, 1):
            struct_id = f"{comp_name}_s{idx:03d}"
            elements = sorted([str(el) for el in structure.composition.elements])
            chemsys = '-'.join(elements)
            
            all_structures.append({
                'id': struct_id,
                'composition': comp_name,
                'chemsys': chemsys,
                'structure': structure
            })
        
        print(f"  Added {len(structures)} structures")
    
    print(f"\nTotal structures loaded: {len(all_structures)}\n")
    
    # Deduplicate structures within each composition
    print("="*70)
    print("Deduplicating Structures (StructureMatcher)")
    print("="*70)
    print(f"Checking for duplicates within each composition...")
    print(f"Tolerances: ltol=0.2, stol=0.2, angle_tol=5")
    sys.stdout.flush()
    
    all_structures, dup_count, dup_details = deduplicate_structures(all_structures)
    
    # Store duplicate information for results tracking
    duplicate_records = []
    if dup_count > 0:
        print(f"\nRemoved {dup_count} duplicate structures:")
        # Group by composition for cleaner output
        dup_by_comp = {}
        for dup in dup_details:
            comp = dup['composition']
            if comp not in dup_by_comp:
                dup_by_comp[comp] = []
            dup_by_comp[comp].append(dup)
            
            # Create result entry for each duplicate
            duplicate_records.append({
                'structure_id': dup['duplicate'],
                'composition': dup['composition'],
                'chemsys': None,
                'mattersim_energy_per_atom': None,
                'energy_above_hull': None,
                'is_stable': None,
                'passed_prescreening': False,
                'duplicate_of': dup['original'],
                'error': f"Duplicate of {dup['original']} (removed before relaxation)"
            })
        
        for comp, dups in sorted(dup_by_comp.items()):
            print(f"  {comp}: {len(dups)} duplicates removed")
            for dup in dups[:5]:  # Show first 5
                print(f"    - {dup['duplicate']} (duplicate of {dup['original']})")
            if len(dups) > 5:
                print(f"    ... and {len(dups) - 5} more")
    else:
        print("\nNo duplicates found")
    
    print(f"\nUnique structures after deduplication: {len(all_structures)}")
    print("="*70 + "\n")
    
    # Pre-fetch MP stable phases for all chemical systems
    unique_chemsys = sorted(set(s['chemsys'] for s in all_structures))
    
    print("="*70)
    print("Pre-fetching MP GGA Phases (Legacy pymatgen MPRester)")
    print("="*70)
    print(f"Unique chemical systems: {len(unique_chemsys)}")
    print("Strategy: Query ALL GGA/GGA+U structures with strict suffix filtering")
    print("  - Uses legacy pymatgen.ext.matproj.MPRester for complete GGA coverage")
    print("  - Only accepts entries with '-GGA' or '-GGA+U' suffix")
    if args.pure_pbe:
        print("  - Filtering to pure GGA-PBE only (PBE+U excluded)")
    else:
        print("  - Accepts both PBE and PBE+U (recommended for accuracy)")
    print("  - MP structures are re-relaxed with MatterSim for consistent energies")
    print("="*70 + "\n")
    
    # Create MatterSim Potential for MP phases (one-time)
    print(f"Creating MatterSim Potential (device={args.device})...")
    try:
        import torch
        import gc
        potential_mp = Potential.from_checkpoint(
            checkpoint_path="MatterSim-v1.0.0-5M.pth",
            device=args.device
        )
        print("  Potential created successfully\n")
    except Exception as e:
        print(f"ERROR: Failed to create MatterSim Potential: {e}")
        return 1
    
    mp_entries_cache = {}
    for chemsys in unique_chemsys:
        print(f"Fetching MP stable phases for {chemsys}...")
        sys.stdout.flush()
        try:
            mp_entries = get_mp_stable_phases_mattersim(
                chemsys, mp_api_key, mp_cache_file, potential_mp, 
                pure_pbe=args.pure_pbe, max_atoms_gpu=args.max_atoms_gpu
            )
            mp_entries_cache[chemsys] = mp_entries
            print(f"  → {len(mp_entries)} stable phases ready\n")
            sys.stdout.flush()
        except Exception as e:
            print(f"  → Error: {e}\n")
            sys.stdout.flush()
            mp_entries_cache[chemsys] = []
    
    # Clean up MP potential
    del potential_mp
    if args.device == 'cuda':
        torch.cuda.empty_cache()
        gc.collect()
    
    print("="*70)
    print("MP stable phases pre-fetch complete")
    print("="*70 + "\n")
    
    # Run pre-screening with batch processing
    print("="*70)
    print("Running MatterSim Pre-screening (True GPU Batch Processing)")
    print("="*70)
    print(f"Processing {len(all_structures)} structures in batches of {args.batch_size}...")
    print(f"Using BatchRelaxer for parallel GPU relaxation within each batch\n")
    sys.stdout.flush()
    
    # Check for existing checkpoint
    batch_suffix = f"_batch{args.batch_id}" if args.batch_id else ""
    checkpoint_file = output_dir / f'prescreening_checkpoint{batch_suffix}.json'
    db_file = output_dir / f'prescreening_structures{batch_suffix}.db'
    
    if checkpoint_file.exists():
        print(f"Found checkpoint file, resuming from previous run...")
        with open(checkpoint_file, 'r') as f:
            checkpoint_data = json.load(f)
            results = checkpoint_data.get('results', [])
            passed = checkpoint_data.get('passed', 0)
            failed = checkpoint_data.get('failed', 0)
            processed_ids = {r['structure_id'] for r in results}
            currently_processing = checkpoint_data.get('currently_processing', None)
        
        print(f"Resuming: {len(results)}/{len(all_structures) + dup_count} already processed")
        
        # Detect if we crashed during a structure (segfault detection)
        if currently_processing:
            print(f"\n   SEGFAULT DETECTED: Crashed while processing '{currently_processing}'")
            
            # Find the problematic structure
            problem_struct = None
            for s in all_structures:
                if s['id'] == currently_processing:
                    problem_struct = s
                    break
            
            if problem_struct:
                results.append({
                    'structure_id': currently_processing,
                    'composition': problem_struct['composition'],
                    'chemsys': problem_struct['chemsys'],
                    'mattersim_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': False,
                    'passed_prescreening': False,
                    'error': 'Auto-skipped due to segfault (spglib crash)'
                })
                processed_ids.add(currently_processing)
                failed += 1
                print(f"     Added '{currently_processing}' to failed list\n")
            else:
                print(f"      Could not find structure in loaded data\n")
        
        # Open existing database for incremental updates
        if db_file.exists():
            print(f"Found existing database: {db_file}")
            db = database_topology(str(db_file))
        else:
            print(f"Creating new database: {db_file}")
            db = database_topology(str(db_file))
        print()
        sys.stdout.flush()
    else:
        results = []
        passed = 0
        failed = 0
        processed_ids = set()
        # Create new database for incremental updates
        db = database_topology(str(db_file))
        print(f"Created new database: {db_file}\n")
    
    # Add duplicate records to results (these are already "processed")
    if duplicate_records:
        print(f"Adding {len(duplicate_records)} duplicate records to results...\n")
        for dup_rec in duplicate_records:
            if dup_rec['structure_id'] not in processed_ids:
                results.append(dup_rec)
                processed_ids.add(dup_rec['structure_id'])
                failed += 1  # Count duplicates as "failed" (not proceeding to VASP)
    
    # Filter unprocessed structures
    structures_to_process = [s for s in all_structures if s['id'] not in processed_ids]
    
    if not structures_to_process:
        print("All structures already processed!\n")
    else:
        print(f"Processing {len(structures_to_process)} remaining structures...\n")
        
        # Process in batches
        import torch
        import gc
        from tqdm import tqdm
        
        n_batches = (len(structures_to_process) + args.batch_size - 1) // args.batch_size
        
        for batch_idx in range(n_batches):
            batch_start = batch_idx * args.batch_size
            batch_end = min((batch_idx + 1) * args.batch_size, len(structures_to_process))
            batch = structures_to_process[batch_start:batch_end]
            
            print(f"\n{'='*70}")
            print(f"Batch {batch_idx + 1}/{n_batches}: Processing {len(batch)} structures")
            print(f"{'='*70}\n")
            
            # Mark the first structure in this batch as currently processing
            # This helps detect segfaults during calculator creation
            if batch:
                first_struct_id = batch[0]['id']
                checkpoint_data_temp = {
                    'results': [{k: v for k, v in r.items() if k != 'pmg'} for r in results],
                    'passed': passed,
                    'failed': failed,
                    'currently_processing': first_struct_id,
                    'timestamp': datetime.now().isoformat()
                }
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint_data_temp, f, indent=2)
            
            # Create MatterSim Potential for this batch (with error handling)
            try:
                potential = Potential.from_checkpoint(
                    checkpoint_path="MatterSim-v1.0.0-5M.pth",
                    device=args.device
                )
            except Exception as e:
                print(f"ERROR creating Potential for batch {batch_idx + 1}: {e}")
                print(f"Skipping entire batch of {len(batch)} structures\n")
                sys.stdout.flush()
                
                # Mark all structures in this batch as failed
                for item in batch:
                    failed += 1
                    results.append({
                        'structure_id': item['id'],
                        'pmg': None,
                        'composition': item['composition'],
                        'chemsys': item['chemsys'],
                        'mattersim_energy_per_atom': None,
                        'energy_above_hull': None,
                        'is_stable': None,
                        'passed_prescreening': False,
                        'error': f'Potential creation failed: {str(e)}'
                    })
                
                # Save checkpoint and continue to next batch
                results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
                checkpoint_data = {
                    'results': results_json,
                    'passed': passed,
                    'failed': failed,
                    'currently_processing': None,
                    'timestamp': datetime.now().isoformat()
                }
                with open(checkpoint_file, 'w') as f:
                    json.dump(checkpoint_data, f, indent=2)
                
                continue
            
            # Pre-validate all structures in batch
            valid_batch = []
            invalid_count = 0
            
            print(f"Pre-validating {len(batch)} structures...")
            for item in batch:
                struct_id = item['id']
                structure = item['structure']
                
                is_valid, validation_error = validate_structure(structure)
                if not is_valid:
                    print(f"  SKIPPING {struct_id}: {validation_error}")
                    invalid_count += 1
                    failed += 1
                    results.append({
                        'structure_id': struct_id,
                        'pmg': None,
                        'composition': item['composition'],
                        'chemsys': item['chemsys'],
                        'mattersim_energy_per_atom': None,
                        'energy_above_hull': None,
                        'is_stable': None,
                        'passed_prescreening': False,
                        'error': f"Invalid structure: {validation_error}"
                    })
                else:
                    valid_batch.append(item)
            
            if invalid_count > 0:
                print(f"  Skipped {invalid_count} invalid structures")
            print(f"  Proceeding with {len(valid_batch)} valid structures\n")
            
            # Batch relax all valid structures using BatchRelaxer
            if valid_batch:
                print(f"Batch relaxing {len(valid_batch)} structures (true GPU parallel processing)...")
                sys.stdout.flush()
                
                batch_relax_results = relax_structure_mattersim(
                    valid_batch,
                    potential,
                    fmax=0.01,
                    max_steps=500,
                    max_natoms_per_batch=args.max_atoms_gpu
                )
                
                print(f"Batch relaxation complete. Processing results...\n")
                sys.stdout.flush()
                
                # Process batch relaxation results
                for relax_result in batch_relax_results:
                    struct_id = relax_result['structure_id']
                    relaxed_struct = relax_result['relaxed_structure']
                    energy_per_atom = relax_result['energy_per_atom']
                    relax_error = relax_result['error']
                    chemsys = relax_result['chemsys']
                    
                    # Check if relaxation failed
                    if relaxed_struct is None or energy_per_atom is None:
                        print(f"  ERROR on {struct_id}: {relax_error}")
                        failed += 1
                        results.append({
                            'structure_id': struct_id,
                            'pmg': None,
                            'composition': relax_result['composition'],
                            'chemsys': chemsys,
                            'mattersim_energy_per_atom': None,
                            'energy_above_hull': None,
                            'is_stable': None,
                            'passed_prescreening': False,
                            'error': relax_error
                        })
                        continue
                    
                    # Relaxation succeeded - compute hull distance
                    mp_entries = mp_entries_cache.get(chemsys, [])
                    
                    if not mp_entries:
                        print(f"  WARNING: No MP reference phases for {chemsys}, auto-passing {struct_id}")
                        e_hull = None
                        passed_prescreening = True
                        passed += 1
                    else:
                        try:
                            e_hull = compute_energy_above_hull(relaxed_struct, energy_per_atom, mp_entries)
                        except Exception as hull_err:
                            print(f"  ERROR computing hull for {struct_id}: {hull_err}")
                            e_hull = None
                        
                        if e_hull is None:
                            print(f"  WARNING: Incomplete phase diagram for {chemsys}, auto-passing {struct_id}")
                            passed_prescreening = True
                            passed += 1
                        elif e_hull < args.hull_threshold:
                            passed_prescreening = True
                            passed += 1
                        else:
                            passed_prescreening = False
                            failed += 1
                    
                    results.append({
                        'structure_id': struct_id,
                        'pmg': relaxed_struct,
                        'composition': relax_result['composition'],
                        'chemsys': chemsys,
                        'mattersim_energy_per_atom': float(energy_per_atom),
                        'energy_above_hull': float(e_hull) if e_hull is not None else None,
                        'is_stable': e_hull < 0.001 if e_hull is not None else None,
                        'passed_prescreening': passed_prescreening
                    })
            
            # Clean up potential after batch
            del potential
            if args.device == 'cuda':
                torch.cuda.empty_cache()
                gc.collect()
            
            # Save checkpoint after each batch (clear currently_processing marker)
            results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
            checkpoint_data = {
                'results': results_json,
                'passed': passed,
                'failed': failed,
                'currently_processing': None,  # Clear marker after batch completes
                'timestamp': datetime.now().isoformat()
            }
            with open(checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
            
            # Save structures to database incrementally after each batch
            batch_db_count = 0
            batch_failed_count = 0
            batch_start_idx = max(0, len(results) - len(batch))
            tolerances_db = [5e-2, 1e-2, 1e-3, 1e-4, 1e-5]
            
            for i in range(batch_start_idx, len(results)):
                res = results[i]
                
                # Find original structure for this result
                orig_struct = None
                for s in batch:
                    if s['id'] == res['structure_id']:
                        orig_struct = s['structure']
                        break
                
                if res['pmg'] is not None:
                    # Valid structure with successful relaxation - try progressive tolerances
                    added = False
                    for tol in tolerances_db:
                        try:
                            xtal = pyxtal()
                            xtal.from_seed(res['pmg'], tol=tol)
                            if not xtal.valid:
                                continue
                            if len(xtal.check_short_distances(r=0.5)) > 0:
                                continue
                            db.add_xtal(
                                xtal,
                                kvp={
                                    'structure_id': res['structure_id'],
                                    'e_above_hull': res['energy_above_hull'],
                                    'e_mattersim': res['mattersim_energy_per_atom'],
                                    'composition': res['composition'],
                                    'chemsys': res['chemsys'],
                                    'passed_prescreening': res['passed_prescreening'],
                                    'status': 'valid',
                                    'symmetrized': True
                                }
                            )
                            batch_db_count += 1
                            added = True
                            break
                        except Exception:
                            continue
                    
                    if not added:
                        # PyXtal failed on relaxed structure - try original structure with progressive tolerances
                        tqdm.write(f"  Warning: PyXtal failed on relaxed structure {res['structure_id']}, trying original structure")
                        if orig_struct is not None:
                            orig_added = False
                            for tol in tolerances_db:
                                try:
                                    xtal_orig = pyxtal()
                                    xtal_orig.from_seed(orig_struct, tol=tol)
                                    db.add_xtal(
                                        xtal_orig,
                                        kvp={
                                            'structure_id': res['structure_id'],
                                            'e_above_hull': res['energy_above_hull'],
                                            'e_mattersim': res['mattersim_energy_per_atom'],
                                            'composition': res['composition'],
                                            'chemsys': res['chemsys'],
                                            'passed_prescreening': False,  # Never proceed to DFT (use relaxed structure failed)
                                            'status': 'failed_symmetrization_after_relax',
                                            'symmetrized': True,
                                            'note': 'original_structure_saved_due_to_relaxed_pyxtal_failure'
                                        }
                                    )
                                    batch_db_count += 1
                                    orig_added = True
                                    break
                                except Exception:
                                    continue
                            
                            if not orig_added:
                                tqdm.write(f"  Error: Could not save even original structure for {res['structure_id']} (PyXtal failed at all tolerances)")
                                batch_failed_count += 1
                        else:
                            batch_failed_count += 1
                else:
                    # Invalid/failed structure - try to save original structure for book-keeping
                    if orig_struct is not None:
                        orig_added = False
                        for tol in tolerances_db:
                            try:
                                xtal_failed = pyxtal()
                                xtal_failed.from_seed(orig_struct, tol=tol)
                                db.add_xtal(
                                    xtal_failed,
                                    kvp={
                                        'structure_id': res['structure_id'],
                                        'e_above_hull': None,
                                        'e_mattersim': None,
                                        'composition': res['composition'],
                                        'chemsys': res['chemsys'],
                                        'passed_prescreening': False,  # Never proceed to DFT
                                        'status': res.get('error', 'failed_prescreening')[:200],  # Truncate long errors
                                        'symmetrized': True,
                                        'note': 'original_structure_saved_for_bookkeeping'
                                    }
                                )
                                batch_db_count += 1
                                orig_added = True
                                break
                            except Exception:
                                continue
                        
                        if not orig_added:
                            # Even original structure failed PyXtal at all tolerances - skip database entry
                            # These are still tracked in checkpoint JSON with full error details
                            batch_failed_count += 1
                    else:
                        batch_failed_count += 1
            
            # Explicit commit after each batch
            if hasattr(db, 'db') and hasattr(db.db, 'commit'):
                db.db.commit()
            
            print(f"\nBatch {batch_idx + 1} complete.")
            print(f"  Checkpoint saved: {checkpoint_file.name}")
            print(f"  Database updated: {batch_db_count} structures added")
            if batch_failed_count > 0:
                print(f"  Could not save {batch_failed_count} structures to database (tracked in checkpoint JSON)")
            print(f"  Progress: {len(results)}/{len(all_structures)} structures\n")
    
    # Save results to JSON (exclude pmg Structure objects)
    results_json = [{k: v for k, v in r.items() if k != 'pmg'} for r in results]
    output_data = {
        'summary': {
            'total_structures_loaded': len(all_structures) + dup_count,
            'duplicate_structures_removed': dup_count,
            'unique_structures_processed': len(all_structures),
            'passed_prescreening': passed,
            'failed_prescreening': failed,
            'hull_threshold': args.hull_threshold,
            'energy_reference': 'MatterSim-v1.0.0-5M',
            'mp_api': 'legacy_pymatgen_mprester_complete_gga',
            'mp_filtering': 'pure_pbe_only' if args.pure_pbe else 'mixed_pbe_pbeU',
            'mp_suffix_filter': 'strict_-GGA_suffix' if args.pure_pbe else 'strict_-GGA_or_-GGA+U_suffix'
        },
        'results': results_json
    }
    
    output_file = output_dir / f'prescreening_stability{batch_suffix}.json'
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nJSON results saved to: {output_file}")
    
    # Database was already saved incrementally during batch processing (see lines 770-808)
    # Get final database statistics
    print("\nDatabase Statistics:")
    if db_file.exists():
        try:
            # Count structures in database
            db_check = database_topology(str(db_file))
            if hasattr(db_check, 'db') and hasattr(db_check.db, 'execute'):
                cursor = db_check.db.execute("SELECT COUNT(*) FROM systems")
                db_count = cursor.fetchone()[0]
                print(f"  PyXtal database: {db_file}")
                print(f"  Structures in database: {db_count}")
            else:
                print(f"  PyXtal database: {db_file} (count unavailable)")
        except Exception as e:
            print(f"  PyXtal database: {db_file} (statistics unavailable: {e})")
    else:
        print(f"  No database file found (all structures may have failed)")
    
    # Remove checkpoint file (no longer needed)
    try:
        if checkpoint_file.exists():
            checkpoint_file.unlink()
            print(f"  Checkpoint file removed")
    except Exception as e:
        print(f"  Warning: Could not remove checkpoint file: {e}")
    
    print("\n" + "="*70)
    print("Pre-screening Complete")
    print("="*70)
    print(f"Total structures loaded: {len(all_structures) + dup_count}")
    print(f"Duplicates removed: {dup_count}")
    print(f"Unique structures processed: {len(all_structures)}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"JSON: {output_file}")
    print(f"Database: {db_file}")
    print("="*70)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
