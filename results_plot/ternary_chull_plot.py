#!/usr/bin/env python3
"""
Plot ternary convex hull diagrams including newly discovered stable electrides.

This script plots formation energy vs composition for ternary systems, showing:
- The convex hull (connected triangular facets)
- Stable phases (on the hull)
- Metastable phases (within threshold above hull)
- Newly discovered electrides highlighted

Usage:
    python3 ternary_chull_plot.py --csv Ter-Ele-HT/stable_electrides.csv \
        --dft-results Ter-Ele-HT/dft_stability_results.json \
        --mp-phases Ter-Ele-HT/mp_vaspdft.json \
        --systems K-B-O Cs-Al-S

    python3 ternary_chull_plot.py --db Ter-Ele-HT/stable_electrides.db \
        --dft-results Ter-Ele-HT/dft_stability_results.json \
        --mp-phases Ter-Ele-HT/mp_vaspdft.json \
        --systems K-B-O Cs-Al-S \
        --output-dir ternary_hull_plots \
        --e-above-hull-max 0.05
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Polygon
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram


# Valence electrons for excess electron calculation
VALENCE_ELECTRONS = {
    'H': 1, 'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1,
    'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
    'B': 3, 'Al': 3, 'Ga': 3, 'In': 3, 'Tl': 3,
    'Sc': 3, 'Y': 3,
    'C': 4, 'Si': 4, 'Ge': 4, 'Sn': 4, 'Pb': 4,
    'N': 3, 'P': 3, 'As': 3, 'Sb': 3, 'Bi': 3,
    'O': 2, 'S': 2, 'Se': 2, 'Te': 2, 'Po': 2,
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'At': 1,
}

# Element groups matching search_ternary_electrides.py
GROUP_I_II_III = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba', 'Sc', 'Y']
GROUP_III_IV = ['B', 'Al', 'Ga', 'In', 'C', 'Si', 'Ge', 'Sn', 'Pb']
GROUP_V_VI_VII = ['N', 'P', 'As', 'Sb', 'O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I', 'H']


def reorder_elements_for_boundary(elements: List[str]) -> Tuple[List[str], List[int]]:
    """
    Reorder elements to match search_ternary_electrides.py convention:
    A (Group I/II/III metals) - B (Group III/IV semi-metals) - C (Group V/VI/VII non-metals)
    
    This ensures the N_excess boundary calculation matches the search script logic.
    
    Args:
        elements: List of 3 elements in any order (usually alphabetical)
    
    Returns:
        (reordered_elements, index_mapping): 
            - reordered_elements: [A, B, C] in chemical group order
            - index_mapping: Indices to convert from chemical order back to original order
    """
    elem_groups = []
    for elem in elements:
        if elem in GROUP_I_II_III:
            elem_groups.append((0, elem))  # Group A
        elif elem in GROUP_III_IV:
            elem_groups.append((1, elem))  # Group B
        elif elem in GROUP_V_VI_VII:
            elem_groups.append((2, elem))  # Group C
        else:
            # Fallback for elements not in defined groups
            elem_groups.append((3, elem))
    
    # Sort by group
    elem_groups_sorted = sorted(elem_groups, key=lambda x: x[0])
    reordered = [elem for _, elem in elem_groups_sorted]
    
    # Create index mapping: reordered[i] corresponds to elements[index_mapping[i]]
    index_mapping = [elements.index(elem) for elem in reordered]
    
    return reordered, index_mapping

ELECTRONEGATIVE_ELEMENTS = {'N', 'P', 'As', 'Sb', 'Bi', 'O', 'S', 'Se', 'Te', 'Po', 'F', 'Cl', 'Br', 'I', 'At', 'H'}


def formula_to_latex(formula: str) -> str:
    """Convert chemical formula to LaTeX format with subscripts.
    
    Handles both simple subscripts (Ca3 -> Ca$_{3}$) and 
    parentheses subscripts (Ca3(AlP2)2 -> Ca$_{3}$(AlP$_{2}$)$_{2}$).
    Ignores subscript 1.
    
    Examples:
        Ca7Al1P5 -> Ca$_{7}$AlP$_{5}$
        K6BO4 -> K$_{6}$BO$_{4}$
        Ca3(AlP2)2 -> Ca$_{3}$(AlP$_{2}$)$_{2}$
    """
    import re
    
    # First, handle closing parenthesis followed by number: )2 -> )$_{2}$
    def replace_paren(m):
        count = m.group(1)
        if count == '1':
            return ')'
        else:
            return f')$_{{{count}}}$'
    
    result = re.sub(r'\)(\d+)', replace_paren, formula)
    
    # Then handle element followed by number: Ca2 -> Ca$_{2}$
    def replace_elem(m):
        element = m.group(1)
        count = m.group(2)
        if count == '1':
            return element
        else:
            return f'{element}$_{{{count}}}$'
    
    result = re.sub(r'([A-Z][a-z]?)(\d+)', replace_elem, result)
    
    return result


def calculate_nexcess_boundaries_ternary(elem_A: str, elem_B: str, elem_C: str, 
                                         n_excess_max: float = 2.0, max_atoms: int = 15) -> Tuple[List[Tuple], List[Tuple]]:
    """
    Calculate exact composition boundaries for N_excess region in ternary A-B-C system.
    This matches exactly the logic in search_ternary_electrides.py.
    
    Iterates through all reduced integer combinations (l, m, n) with l+m+n <= max_atoms
    and finds compositions where 0 < N_excess <= n_excess_max, then computes the 
    convex hull of these valid compositions to get the boundary polygon.
    
    For ternary system A_l B_m C_n:
    - N_excess = val_A * l + val_B * m - val_C * n
    - x_A = l / (l + m + n), x_B = m / (l + m + n)
    
    Args:
        elem_A: Symbol of element A (vertex at bottom-left)
        elem_B: Symbol of element B (vertex at top)
        elem_C: Symbol of element C (vertex at bottom-right)
        n_excess_max: Maximum excess electrons (default 2 for ternary)
        max_atoms: Maximum total atoms in formula (default 15 for ternary)
    
    Returns:
        (boundary_points, None): List of (x_A, x_B) coordinates defining boundary polygon
    """
    from math import gcd
    from functools import reduce
    
    def gcd_multiple(numbers):
        """Calculate GCD of multiple numbers."""
        return reduce(gcd, numbers)
    
    if elem_A not in VALENCE_ELECTRONS or elem_B not in VALENCE_ELECTRONS or elem_C not in VALENCE_ELECTRONS:
        return None, None
    
    val_A = VALENCE_ELECTRONS[elem_A]
    val_B = VALENCE_ELECTRONS[elem_B]
    val_C = VALENCE_ELECTRONS[elem_C]
    
    valid_points = []
    
    # Iterate through all possible stoichiometries
    for l in range(1, max_atoms):
        for m in range(1, max_atoms):
            for n in range(1, max_atoms):
                if l + m + n > max_atoms:
                    continue
                
                # Reduce to smallest integer ratio
                g = gcd_multiple([l, m, n])
                l_p, m_p, n_p = l // g, m // g, n // g
                
                # Skip if not reduced form
                if (l_p, m_p, n_p) != (l, m, n):
                    continue
                
                # Calculate excess electrons
                n_excess = val_A * l_p + val_B * m_p - val_C * n_p
                
                # Check if in valid range
                if 0 < n_excess <= n_excess_max:
                    x_A = l_p / (l_p + m_p + n_p)
                    x_B = m_p / (l_p + m_p + n_p)
                    valid_points.append((x_A, x_B))
    
    if not valid_points:
        return None, None
    
    # Compute convex hull to get boundary polygon
    from scipy.spatial import ConvexHull
    try:
        hull = ConvexHull(valid_points)
        boundary_points = [valid_points[i] for i in hull.vertices]
        # Sort points by angle for proper polygon plotting
        center = np.mean(boundary_points, axis=0)
        angles = [np.arctan2(p[1] - center[1], p[0] - center[0]) for p in boundary_points]
        sorted_indices = np.argsort(angles)
        boundary_points = [boundary_points[i] for i in sorted_indices]
        return boundary_points, None
    except:
        # If convex hull fails (e.g., all points are collinear), return all points
        return valid_points, None


def prompt_for_confirmed_electrides(
    stable_electrides: List[Dict],
    metastable_electrides: List[Dict]
) -> List[str]:
    """Prompt user to confirm actual electrides after PARCHG review."""
    print("\n" + "="*70)
    print("ELECTRIDE CONFIRMATION")
    print("="*70)
    print("\nPlease review PARCHG files to confirm actual electrides.")
    print("Candidate electrides:")
    
    all_candidates = stable_electrides + metastable_electrides
    for i, e in enumerate(all_candidates, 1):
        status = "stable" if e in stable_electrides else f"metastable (E={e.get('e_above_hull', 0):.4f} eV/atom)"
        pearson = e['metadata'].get('pearson_symbol', '')
        structure_id = e['entry_id']
        if pearson:
            label = f"{pearson}-{structure_id}"
        else:
            label = structure_id
        print(f"  {i:2d}. {label:30s} - {status}")
    
    print("\nEnter structure_ids of CONFIRMED electrides (comma-separated):")
    print("  (or press Enter to skip - all will be marked as candidates)")
    print("  Example: Cs2Al2S3_s013,Cs6Al2S5_s013")
    
    user_input = input("> ").strip()
    
    if not user_input:
        print("  Skipping confirmation. All structures marked as candidates.")
        return []
    
    confirmed_ids = [s.strip() for s in user_input.split(',')]
    print(f"  Confirmed {len(confirmed_ids)} electrides: {', '.join(confirmed_ids)}")
    return confirmed_ids


def prompt_for_label_positions(phase_data: List[Dict], phase_type: str) -> List[Tuple[float, float]]:
    """Prompt user to input label positions for each phase interactively.
    
    Args:
        phase_data: List of phase data dicts with 'cart_x', 'cart_y', 'formula', 'entry_id'
        phase_type: Description of phase type (e.g., 'MP stable', 'New stable electride')
    
    Returns:
        List of (x_offset, y_offset) tuples in points (matplotlib offset points)
    """
    print(f"\n{'='*70}")
    print(f"INTERACTIVE LABEL POSITIONING: {phase_type}")
    print(f"{'='*70}")
    print(f"Enter label position offsets (in points) for each phase.")
    print(f"Format: x_offset,y_offset (e.g., 30,40 or -20,30)")
    print(f"Default: 20,20 (upper-right of the marker)")
    print(f"{'='*70}\n")
    
    positions = []
    for i, d in enumerate(phase_data, 1):
        formula = d.get('formula', d.get('entry_id', 'Unknown'))
        entry_id = d.get('entry_id', '')
        x_pos = d.get('cart_x', 0)
        y_pos = d.get('cart_y', 0)
        
        print(f"{i:2d}. {formula:15s} (ID: {entry_id:20s}) at ({x_pos:.3f}, {y_pos:.3f})")
        user_input = input(f"    Label offset [default: 20,20]: ").strip()
        
        if not user_input:
            positions.append((20.0, 20.0))
        else:
            try:
                parts = user_input.split(',')
                x_off = float(parts[0].strip())
                y_off = float(parts[1].strip())
                positions.append((x_off, y_off))
            except (ValueError, IndexError):
                print(f"    Invalid input, using default (20, 20)")
                positions.append((20.0, 20.0))
    
    print(f"\n{'='*70}\n")
    return positions


def load_electride_candidates_csv(csv_path: Path, systems: List[str], e_hull_max: float) -> List[Dict]:
    """Load electride candidates from CSV file filtered by chemical systems."""
    import csv
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    # Normalize systems to alphabetical order
    normalized_systems = ['-'.join(sorted(s.split('-'))) for s in systems]
    
    # Filter by chemical systems and e_hull
    filtered_rows = []
    for row in rows:
        dft_e_hull = float(row.get('dft_e_hull', float('inf')))
        if dft_e_hull <= e_hull_max:
            comp = Composition(row['composition'])
            elements = sorted([str(e) for e in comp.elements])
            chemsys = '-'.join(elements)
            
            if chemsys in normalized_systems:
                # Convert numeric fields
                row['vasp_energy_per_atom'] = float(row['vasp_energy_per_atom'])
                row['dft_e_hull'] = dft_e_hull
                if 'spacegroup' in row:
                    row['spacegroup'] = int(row['spacegroup'])
                if 'formula' not in row:
                    row['formula'] = row['composition']
                filtered_rows.append(row)
    
    return filtered_rows


def load_electride_candidates_db(db_path: Path, systems: List[str], e_hull_max: float) -> List[Dict]:
    """Load electride candidates from ASE database filtered by chemical systems."""
    from ase.db import connect
    
    db = connect(str(db_path))
    
    # Normalize systems to alphabetical order
    normalized_systems = ['-'.join(sorted(s.split('-'))) for s in systems]
    
    rows = []
    for row in db.select():
        if row.get('dft_e_hull', float('inf')) <= e_hull_max:
            formula = row.formula
            comp = Composition(formula)
            elements = sorted([str(e) for e in comp.elements])
            chemsys = '-'.join(elements)
            
            if chemsys in normalized_systems:
                rows.append({
                    'formula': formula,
                    'composition': formula,
                    'vasp_energy_per_atom': row.get('vasp_energy_per_atom'),
                    'dft_e_hull': row.get('dft_e_hull'),
                    'structure_id': row.get('structure_id'),
                    'spacegroup': row.get('space_group_number')
                })
    
    return rows


def load_dft_results(json_path: Path) -> Dict:
    """Load DFT stability results."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Handle nested structure with "results" field
    if isinstance(data, dict) and 'results' in data:
        entries_list = data['results']
    else:
        entries_list = data
    
    results = {}
    for entry in entries_list:
        structure_id = entry['structure_id']
        results[structure_id] = {
            'composition': entry['composition'],
            'chemsys': entry['chemsys'],
            'energy_per_atom': entry['vasp_energy_per_atom'],
            'e_hull': entry['energy_above_hull']
        }
    
    return results


def get_pearson_symbols_from_db(db_path: Path, structure_ids: List[str]) -> Dict[str, str]:
    """Get Pearson symbols for structure IDs from ASE database."""
    from ase.db import connect
    
    if not db_path.exists():
        return {}
    
    db = connect(str(db_path))
    pearson_map = {}
    
    for sid in structure_ids:
        try:
            row = db.get(structure_id=sid)
            pearson_map[sid] = row.get('pearson_symbol', '')
        except KeyError:
            pearson_map[sid] = ''
    
    return pearson_map


def load_mp_phases(json_path: Path) -> List[Dict]:
    """Load MP reference phases."""
    with open(json_path, 'r') as f:
        return json.load(f)


def create_entries_for_system(
    system: str,
    electride_candidates: List[Dict],
    dft_results: Dict,
    mp_phases: List[Dict]
) -> Tuple[List[ComputedEntry], Dict[str, Dict], List[str]]:
    """Create ComputedEntry objects for a given chemical system."""
    
    entries = []
    entry_metadata = {}
    temp_electride_data = []  # Temporary storage for deduplication
    
    # Normalize system to alphabetical order
    system_normalized = '-'.join(sorted(system.split('-')))
    
    # Get elements in this system (alphabetically sorted)
    elements = sorted(system_normalized.split('-'))
    
    # Collect electride candidates for deduplication
    for row in electride_candidates:
        comp = Composition(row['composition'])
        comp_elements = sorted([str(e) for e in comp.elements])
        chemsys = '-'.join(comp_elements)
        
        if chemsys == system_normalized:
            energy_per_atom = row['vasp_energy_per_atom']
            n_atoms = comp.num_atoms
            total_energy = energy_per_atom * n_atoms
            
            structure_id = row.get('structure_id', row['formula'])
            spacegroup = row.get('spacegroup', None)
            
            temp_electride_data.append({
                'structure_id': structure_id,
                'composition': comp,
                'formula': row.get('formula', row['composition']),
                'energy_per_atom': energy_per_atom,
                'total_energy': total_energy,
                'spacegroup': spacegroup,
                'n_atoms': n_atoms
            })
    
    # Deduplicate electrides: For ternary plots, keep only lowest energy per composition
    # (ignore spacegroup since same composition overlaps in top-down view)
    deduplicated_electrides = []
    comp_groups = {}
    
    for data in temp_electride_data:
        comp_key = data['composition'].reduced_formula
        if comp_key not in comp_groups:
            comp_groups[comp_key] = data
        else:
            # Keep the one with lower energy per atom
            if data['energy_per_atom'] < comp_groups[comp_key]['energy_per_atom']:
                comp_groups[comp_key] = data
    
    deduplicated_electrides = list(comp_groups.values())
    
    # Add deduplicated electride candidates
    for data in deduplicated_electrides:
        entry = ComputedEntry(
            composition=data['composition'],
            energy=data['total_energy'],
            entry_id=data['structure_id']
        )
        entries.append(entry)
        
        # Store metadata
        entry_metadata[data['structure_id']] = {
            'is_electride': True,
            'spacegroup': data['spacegroup'],
            'pearson_symbol': None,
            'formula': data['formula']
        }
    
    # Add MP reference phases (including elemental and binary entries)
    for phase in mp_phases:
        phase_chemsys = phase['chemsys']
        comp_dict = phase['composition']
        comp = Composition(comp_dict)
        
        # Normalize phase chemsys
        phase_chemsys_normalized = '-'.join(sorted(phase_chemsys.split('-')))
        
        # Get phase elements
        phase_elements = sorted(list(comp_dict.keys()))
        
        # Include if it matches the system or is a subsystem (elemental, binary)
        is_system_phase = phase_chemsys_normalized == system_normalized
        is_subsystem = all(elem in elements for elem in phase_elements)
        
        if is_system_phase or is_subsystem:
            total_energy = phase['energy']
            entry_id = phase.get('entry_id', phase.get('mp_id', 'unknown'))
            
            entry = ComputedEntry(
                composition=comp,
                energy=total_energy,
                entry_id=entry_id
            )
            entries.append(entry)
            
            # Store metadata
            entry_metadata[entry_id] = {
                'is_electride': False,
                'mp_id': phase.get('mp_id', entry_id.split('-')[0] if '-' in entry_id else entry_id)
            }
    
    return entries, entry_metadata, elements


def cart2bary(coords_cart, elements):
    """Convert Cartesian to barycentric coordinates for ternary plot.
    
    Args:
        coords_cart: Array of [x, y] Cartesian coordinates
        elements: List of 3 element names [A, B, C]
    
    Returns:
        Array of [a, b, c] barycentric coordinates (fractions of A, B, C)
    """
    # Standard ternary plot: A at top, B at bottom-left, C at bottom-right
    # A = (0.5, sqrt(3)/2), B = (0, 0), C = (1, 0)
    sqrt3 = np.sqrt(3)
    
    x, y = coords_cart
    
    # Barycentric coordinates
    c = x - y / sqrt3
    b = 2 * y / sqrt3
    a = 1 - b - c
    
    return np.array([a, b, c])


def bary2cart(coords_bary):
    """Convert barycentric to Cartesian coordinates for ternary plot.
    
    Args:
        coords_bary: Array of [a, b, c] barycentric coordinates
    
    Returns:
        Array of [x, y] Cartesian coordinates
    """
    a, b, c = coords_bary
    
    sqrt3 = np.sqrt(3)
    x = 0.5 * a + c
    y = sqrt3 / 2 * a
    
    return np.array([x, y])


def plot_ternary_hull(
    system: str,
    entries: List[ComputedEntry],
    entry_metadata: Dict[str, Dict],
    elements: List[str],
    output_dir: Path,
    e_hull_max: float = 0.05,
    confirmed_ids: List[str] = None
):
    """Plot ternary convex hull with formation energy as color."""
    
    # Separate MP entries from electride entries
    mp_entries = [e for e in entries if not entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
    
    # Create two phase diagrams
    pd_mp_only = PhaseDiagram(mp_entries)
    pd_all = PhaseDiagram(entries)
    
    # Get reference energies for pure elements
    ref_energies = {}
    for elem in elements:
        min_energy = 1e8
        for entry in entries:
            if len(entry.composition.elements) == 1 and str(entry.composition.elements[0]) == elem:
                energy_per_atom = entry.energy / entry.composition.num_atoms
                if energy_per_atom < min_energy:
                    min_energy = energy_per_atom
        ref_energies[elem] = min_energy
    
    # Prepare data for plotting
    all_data = []
    for entry in entries:
        comp = entry.composition
        
        # Skip if not in this ternary system
        comp_elems = sorted([str(e) for e in comp.elements])
        if len(comp_elems) > 3:
            continue
        
        # Get barycentric coordinates (fractions of each element)
        bary_coords = np.array([comp.get_atomic_fraction(elem) for elem in elements])
        cart_coords = bary2cart(bary_coords)
        
        # Calculate formation energy
        energy_per_atom = entry.energy / comp.num_atoms
        formation_energy = energy_per_atom - sum(bary_coords[i] * ref_energies[elements[i]] 
                                                  for i in range(3))
        
        # Get energy above hull
        e_above_hull = pd_all.get_e_above_hull(entry)
        
        # Get metadata
        metadata = entry_metadata.get(entry.entry_id, {})
        is_electride = metadata.get('is_electride', False)
        
        # Check if on hull
        on_hull = entry in pd_all.stable_entries
        
        # Check if on MP-only hull
        on_mp_hull = entry in pd_mp_only.stable_entries
        
        # Get energy above MP-only hull for electrides
        if is_electride:
            try:
                decomp = pd_mp_only.get_decomposition(entry.composition)
                ref_energy = sum(amt * pd_mp_only.get_form_energy_per_atom(e) 
                                for e, amt in decomp.items())
                e_above_mp_hull = formation_energy - ref_energy
            except:
                e_above_mp_hull = 0.0
        else:
            e_above_mp_hull = 0.0
        
        all_data.append({
            'cart_x': cart_coords[0],
            'cart_y': cart_coords[1],
            'bary': bary_coords,
            'formation_energy': formation_energy,
            'e_above_hull': e_above_hull,
            'e_above_mp_hull': e_above_mp_hull,
            'on_hull': on_hull,
            'on_mp_hull': on_mp_hull,
            'is_electride': is_electride,
            'formula': comp.reduced_formula,
            'entry_id': entry.entry_id,
            'entry': entry,
            'metadata': metadata
        })
    
    # Identify electride candidates
    hull_phases = [d for d in all_data if d['on_hull']]
    metastable_phases = [d for d in all_data if not d['on_hull'] and d['e_above_hull'] <= e_hull_max]
    new_stable = [d for d in hull_phases if d['is_electride']]
    
    # For metastable electrides, keep only lowest energy per composition
    new_meta_all = [d for d in metastable_phases if d['is_electride']]
    comp_groups = {}
    for d in new_meta_all:
        comp_key = d['formula']
        if comp_key not in comp_groups:
            comp_groups[comp_key] = d
        else:
            if d['e_above_hull'] < comp_groups[comp_key]['e_above_hull']:
                comp_groups[comp_key] = d
    new_meta = list(comp_groups.values())
    
    # Use provided confirmed_ids or empty list
    if confirmed_ids is None:
        confirmed_ids = []
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Prepare colormap for energy above hull (rainbow: purple=low, red=high)
    all_e_above_hull = [d['e_above_hull'] for d in all_data 
                        if d['on_hull'] or d['e_above_hull'] <= e_hull_max]
    vmin, vmax = 0.0, max(all_e_above_hull) if all_e_above_hull else 0.1
    cmap = plt.cm.rainbow  # Purple (stable) to red (metastable)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    
    # Get the purple color from colormap (at 0.0 energy above hull)
    stable_purple = cmap(norm(0.0))
    
    # Draw ternary axes with black color
    draw_ternary_axes(ax, elements, 'black')
    
    # Reorder elements to match search script convention for boundary calculation
    reordered_elems, index_mapping = reorder_elements_for_boundary(elements)
    boundary_points, _ = calculate_nexcess_boundaries_ternary(
        reordered_elems[0], reordered_elems[1], reordered_elems[2], 
        n_excess_max=2.0, 
        max_atoms=20
    )
    if boundary_points is not None and len(boundary_points) >= 3:
        # Convert barycentric coordinates to Cartesian
        # boundary_points are in reordered (chemical group) order: (x_A, x_B) where A, B, C are metal, semi-metal, non-metal
        # Need to convert to plot order (alphabetical): elements[0], elements[1], elements[2]
        polygon_points = []
        for x_reordered_0, x_reordered_1 in boundary_points:
            x_reordered_2 = 1 - x_reordered_0 - x_reordered_1
            x_reordered = [x_reordered_0, x_reordered_1, x_reordered_2]
            
            # Convert from reordered (chemical) to alphabetical order
            x_plot = [0, 0, 0]
            for i, idx in enumerate(index_mapping):
                x_plot[idx] = x_reordered[i]
            
            bary = np.array(x_plot)
            cart = bary2cart(bary)
            polygon_points.append(cart)
        
        # Create and fill polygon
        polygon = Polygon(polygon_points, alpha=0.35, facecolor='gold', 
                        edgecolor='none', zorder=1, label='')
        ax.add_patch(polygon)
    
    # Draw full hull facets first (purple solid lines, lower zorder)
    draw_hull_facets(ax, pd_all, all_data, elements, linestyle='-', color=stable_purple, 
                    linewidth=2.5, alpha=0.8, label='Convex hull', zorder=2)
    
    # Draw MP-only hull facets on top (black dashed lines, higher zorder)
    draw_hull_facets(ax, pd_mp_only, all_data, elements, linestyle='--', color='black', 
                    linewidth=2.5, alpha=0.8, label='MP-only hull', zorder=3)
    
    # Get MP stable phases
    mp_stable = [d for d in hull_phases if not d['is_electride']]
    
    # Separate elemental and compound MP phases
    mp_elemental = [d for d in mp_stable if len(d['entry'].composition.elements) == 1]
    mp_compounds = [d for d in mp_stable if len(d['entry'].composition.elements) > 1]
    
    # Plot all MP stable phases (markers only)
    if mp_stable:
        mp_x = [d['cart_x'] for d in mp_stable]
        mp_y = [d['cart_y'] for d in mp_stable]
        ax.scatter(mp_x, mp_y, color='black', s=180, marker='o', edgecolors='none', 
                  linewidth=0, zorder=4, label='MP stable')
    
    # Only prompt for and label compound MP phases
    if mp_compounds:
        mp_label_positions = prompt_for_label_positions(mp_compounds, "MP stable compound phases")
        
        # Label compound MP phases
        for i, d in enumerate(mp_compounds):
            formula_latex = formula_to_latex(d['formula'])
            x_offset, y_offset = mp_label_positions[i]
            
            ax.annotate(formula_latex, xy=(d['cart_x'], d['cart_y']),
                       xytext=(x_offset, y_offset), textcoords='offset points',
                       fontsize=18, ha='left', color='black', fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.4', facecolor='lightgray',
                               alpha=0.8, edgecolor='black', linewidth=1.2),
                       arrowprops=dict(arrowstyle='-', color='black', lw=1.2, alpha=0.7))
    
    # Separate new structures by stable/metastable and electride/non-electride
    new_stable = [d for d in hull_phases if d['is_electride']]
    new_meta = [d for d in metastable_phases if d['is_electride']]
    
    # Plot stable new structures
    if new_stable:
        # Prompt for label positions
        new_stable_positions = prompt_for_label_positions(new_stable, "New stable structures")
        
        # Plot markers (triangle for confirmed electride, square for candidate/non-electride, purple fill, red edge)
        for d in new_stable:
            is_confirmed = d['entry_id'] in confirmed_ids
            marker = '^' if is_confirmed else 's'  # Triangle for confirmed electride, square for candidate
            facecolor = stable_purple  # Always purple fill for stable structures
            marker_size = 220 if is_confirmed else 180
            ax.scatter(d['cart_x'], d['cart_y'], color=facecolor, s=marker_size, marker=marker, 
                      edgecolors='red', linewidth=2.5, zorder=6)
        
        # Label new stable structures
        for i, d in enumerate(new_stable):
            pearson = d['metadata'].get('pearson_symbol', '')
            structure_id = d['entry_id']
            formula_part = structure_id.split('_')[0] if '_' in structure_id else d['formula']
            formula_latex = formula_to_latex(formula_part)
            if pearson:
                label_text = f"{pearson}-{formula_latex}"
            else:
                label_text = formula_part
            
            x_offset, y_offset = new_stable_positions[i]
            
            ax.annotate(label_text, xy=(d['cart_x'], d['cart_y']),
                       xytext=(x_offset, y_offset), textcoords='offset points',
                       fontsize=18, ha='left', color=stable_purple, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.4', facecolor='lavender',
                               alpha=0.9, edgecolor=stable_purple, linewidth=2),
                       arrowprops=dict(arrowstyle='-', color=stable_purple, lw=1.5, alpha=0.7))
    
    # Plot metastable new structures
    if new_meta:
        # Prompt for label positions
        new_meta_positions = prompt_for_label_positions(new_meta, "Metastable structures")
        
        # Plot with colormap fill (based on energy above hull), red edges
        for d in new_meta:
            is_confirmed = d['entry_id'] in confirmed_ids
            marker = '^' if is_confirmed else 's'  # Triangle for confirmed electride, square for candidate
            facecolor = cmap(norm(d['e_above_hull']))  # Always filled with colormap color
            marker_size = 150 if is_confirmed else 120
            ax.scatter(d['cart_x'], d['cart_y'], c=[facecolor], s=marker_size, marker=marker, 
                      edgecolors='red', linewidth=2, zorder=5)
        
        # Label metastable structures
        for i, d in enumerate(new_meta):
            pearson = d['metadata'].get('pearson_symbol', '')
            structure_id = d['entry_id']
            formula_part = structure_id.split('_')[0] if '_' in structure_id else d['formula']
            formula_latex = formula_to_latex(formula_part)
            e_hull_val = d['e_above_hull']
            if pearson:
                label_text = f"{pearson}-{formula_latex}\n({e_hull_val:.3f} eV/atom)"
            else:
                label_text = f"{formula_part}\n({e_hull_val:.3f} eV/atom)"
            
            x_offset, y_offset = new_meta_positions[i]
            
            ax.annotate(label_text, xy=(d['cart_x'], d['cart_y']),
                       xytext=(x_offset, y_offset), textcoords='offset points',
                       fontsize=18, ha='left', color='darkorange',
                       bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                               alpha=0.8, edgecolor='darkorange', linewidth=1.2),
                       arrowprops=dict(arrowstyle='-', color='darkorange', lw=1.2, alpha=0.7))
    
    # Add legend entries (only if points exist)
    if new_stable or new_meta:
        # Count confirmed electrides and non-electrides
        all_structures = new_stable + new_meta
        confirmed_count = sum(1 for d in all_structures if d['entry_id'] in confirmed_ids)
        candidate_count = sum(1 for d in all_structures if d['entry_id'] not in confirmed_ids)
        
        if confirmed_count > 0:
            ax.scatter([], [], facecolors='none', s=220, marker='^', 
                      edgecolors='red', linewidth=2.5, label='Electride')
        if candidate_count > 0:
            ax.scatter([], [], facecolors='none', s=180, marker='s', 
                      edgecolors='red', linewidth=2, label='Not electride')
    
    # Configure plot
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Set title (only show energy range if there are metastable structures)
    if new_meta:
        title = f'{system} ($E_{{\\mathrm{{hull}}}} \\leq {e_hull_max:.2f}$ eV/atom)'
    else:
        title = f'{system}'
    ax.set_title(title, fontsize=24, fontweight='bold', pad=25)
    
    ax.legend(loc='upper left', fontsize=18, framealpha=0.0, edgecolor='none')
    
    # Save figure
    output_path = output_dir / f'{system.replace("-", "")}_ternary_hull.png'
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    # Also save as PDF
    output_pdf = output_dir / f'{system.replace("-", "")}_ternary_hull.pdf'
    fig.savefig(output_pdf, bbox_inches='tight')
    print(f"Saved: {output_pdf}")
    
    plt.close(fig)
    
    # Print summary
    electride_data = [d for d in all_data if d['is_electride']]
    electrides_below_mp_hull = [d for d in electride_data if d['e_above_mp_hull'] < -1e-6]
    
    print(f"\n{system} System Summary:")
    print(f"  Total entries: {len(entries)}")
    print(f"  MP phases: {len(mp_entries)}")
    print(f"  Stable phases (on hull): {len(hull_phases)}")
    print(f"  Electride candidates: {len(electride_data)}")
    print(f"  New stable electrides: {len(new_stable)}")
    print(f"  Electrides below original MP hull: {len(electrides_below_mp_hull)}")
    print(f"  New metastable electrides (E_hull ≤ {e_hull_max}): {len(new_meta)}")
    print(f"  Confirmed electrides: {len(confirmed_ids)}")
    
    if new_stable:
        print(f"\n  Stable electrides:")
        for d in new_stable:
            pearson = d['metadata'].get('pearson_symbol', 'N/A')
            confirmed = " [CONFIRMED]" if d['entry_id'] in confirmed_ids else ""
            print(f"    - {pearson}-{d['formula']} ({d['entry_id']}): "
                  f"E_form = {d['formation_energy']:.4f} eV/atom, "
                  f"E_MP-hull = {d['e_above_mp_hull']:.4f} eV/atom{confirmed}")
    
    if new_meta:
        print(f"\n  Metastable electrides (top 10 by stability):")
        sorted_meta = sorted(new_meta, key=lambda d: d['e_above_hull'])[:10]
        for d in sorted_meta:
            pearson = d['metadata'].get('pearson_symbol', 'N/A')
            confirmed = " [CONFIRMED]" if d['entry_id'] in confirmed_ids else ""
            print(f"    - {pearson}-{d['formula']} ({d['entry_id']}): "
                  f"E_form = {d['formation_energy']:.4f} eV/atom, "
                  f"E_hull = {d['e_above_hull']:.4f} eV/atom{confirmed}")


def draw_hull_facets(ax, phase_diagram, all_data, elements, label=None, zorder=1, **kwargs):
    """Draw the convex hull facets for a ternary phase diagram.
    
    Args:
        ax: Matplotlib axes
        phase_diagram: PhaseDiagram object
        all_data: List of all phase data dicts
        elements: List of element names
        label: Label for legend
        zorder: Drawing order (higher values are drawn on top)
        **kwargs: Plotting kwargs (linestyle, color, linewidth, etc.)
    """
    # Create a mapping from entry to data
    entry_to_data = {d['entry']: d for d in all_data}
    
    # Get stable entries
    stable_entries = phase_diagram.stable_entries
    
    # Get all stable entry coordinates
    stable_coords = []
    stable_entries_list = list(stable_entries)
    
    for entry in stable_entries_list:
        if entry in entry_to_data:
            d = entry_to_data[entry]
            stable_coords.append([d['cart_x'], d['cart_y']])
    
    if len(stable_coords) < 3:
        return  # Not enough points for facets
    
    stable_coords = np.array(stable_coords)
    
    # Use Delaunay triangulation to find facets
    from scipy.spatial import Delaunay
    if len(stable_coords) >= 3:
        try:
            tri = Delaunay(stable_coords)
            
            # Track which edges we've drawn to avoid duplicates
            drawn_edges = set()
            
            # Draw each simplex edge
            first_line = True
            for simplex in tri.simplices:
                # Draw the three edges of the triangle
                for i in range(3):
                    p1_idx = simplex[i]
                    p2_idx = simplex[(i+1)%3]
                    edge = tuple(sorted([p1_idx, p2_idx]))
                    
                    if edge not in drawn_edges:
                        drawn_edges.add(edge)
                        p1 = stable_coords[p1_idx]
                        p2 = stable_coords[p2_idx]
                        
                        # Only add label to first line
                        if first_line and label:
                            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], label=label, zorder=zorder, **kwargs)
                            first_line = False
                        else:
                            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], zorder=zorder, **kwargs)
        except:
            pass  # Skip if triangulation fails


def draw_ternary_axes(ax, elements, color='purple'):
    """Draw the triangular axes for a ternary plot."""
    # Triangle vertices: A at top, B at bottom-left, C at bottom-right
    sqrt3 = np.sqrt(3)
    vertices = np.array([
        [0.5, sqrt3/2],  # A (top)
        [0, 0],          # B (bottom-left)
        [1, 0]           # C (bottom-right)
    ])
    
    # Draw triangle
    triangle = Polygon(vertices, fill=False, edgecolor=color, linewidth=2.5)
    ax.add_patch(triangle)
    
    # Label vertices with element names (larger font)
    offset = 0.06
    ax.text(vertices[0,0], vertices[0,1] + offset, elements[0], 
           ha='center', va='bottom', fontsize=22, fontweight='bold')
    ax.text(vertices[1,0] - offset, vertices[1,1], elements[1], 
           ha='right', va='center', fontsize=22, fontweight='bold')
    ax.text(vertices[2,0] + offset, vertices[2,1], elements[2], 
           ha='left', va='center', fontsize=22, fontweight='bold')
    
    # Set limits with padding
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, sqrt3/2 + 0.1)


def main():
    parser = argparse.ArgumentParser(
        description='Plot ternary convex hull diagrams with newly discovered electrides',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Input sources
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--csv', type=Path,
                           help='Path to stable_electrides.csv')
    input_group.add_argument('--db', type=Path,
                           help='Path to stable_electrides.db')
    
    # Required files
    parser.add_argument('--dft-results', type=Path, required=True,
                       help='Path to dft_stability_results.json')
    parser.add_argument('--mp-phases', type=Path, required=True,
                       help='Path to mp_vaspdft.json')
    
    # Systems to plot
    parser.add_argument('--systems', nargs='+', required=True,
                       help='Chemical systems to plot (e.g., K-B-O Cs-Al-S)')
    
    # Output
    parser.add_argument('--output-dir', type=Path, default=Path('.'),
                       help='Output directory for plots (default: current directory)')
    
    # Energy above hull threshold
    parser.add_argument('--e-above-hull-max', type=float, default=0.05,
                       help='Maximum energy above hull to display (eV/atom, default: 0.05)')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load electride candidates
    print("Loading electride candidates...")
    if args.csv:
        electride_candidates = load_electride_candidates_csv(args.csv, args.systems, args.e_above_hull_max)
        print(f"  Loaded {len(electride_candidates)} electride candidates from CSV")
        # Get Pearson symbols from database if it exists
        db_path = args.csv.parent / 'stable_electrides.db'
        if db_path.exists():
            structure_ids = [row['structure_id'] for row in electride_candidates]
            pearson_map = get_pearson_symbols_from_db(db_path, structure_ids)
            print(f"  Loaded Pearson symbols for {len(pearson_map)} structures from database")
        else:
            pearson_map = {}
    else:
        electride_candidates = load_electride_candidates_db(args.db, args.systems, args.e_above_hull_max)
        print(f"  Loaded {len(electride_candidates)} electride candidates from database")
        # Get Pearson symbols from database
        structure_ids = [row['structure_id'] for row in electride_candidates]
        pearson_map = get_pearson_symbols_from_db(args.db, structure_ids)
        print(f"  Loaded Pearson symbols for {len(pearson_map)} structures")
    
    # Load DFT results
    print("\nLoading DFT results...")
    dft_results = load_dft_results(args.dft_results)
    print(f"  Loaded {len(dft_results)} DFT results")
    
    # Load MP phases
    print("\nLoading MP phases...")
    mp_phases = load_mp_phases(args.mp_phases)
    print(f"  Loaded {len(mp_phases)} MP phases")
    
    # Plot each system
    print("\n" + "="*70)
    for system in args.systems:
        print(f"\nProcessing {system} system...")
        print("-"*70)
        
        # Verify it's a ternary system
        if len(system.split('-')) != 3:
            print(f"  Error: {system} is not a ternary system (must have 3 elements)")
            continue
        
        # Create entries
        entries, entry_metadata, elements = create_entries_for_system(
            system, electride_candidates, dft_results, mp_phases
        )
        
        if not entries:
            print(f"  Warning: No entries found for {system} system")
            continue
        
        # Add Pearson symbols to metadata
        for entry_id, metadata in entry_metadata.items():
            if metadata.get('is_electride', False) and entry_id in pearson_map:
                metadata['pearson_symbol'] = pearson_map[entry_id]
        
        electride_count = sum(1 for m in entry_metadata.values() if m.get('is_electride', False))
        if electride_count == 0:
            print(f"  Warning: No electride candidates found for {system} system")
        
        # Create phase diagram to identify stable and metastable candidates
        pd_all = PhaseDiagram(entries)
        
        # Collect all electride candidates for confirmation prompt
        new_stable = []
        new_meta_all = []
        for entry in entries:
            metadata = entry_metadata.get(entry.entry_id, {})
            if metadata.get('is_electride', False):
                e_above_hull = pd_all.get_e_above_hull(entry)
                on_hull = entry in pd_all.stable_entries
                entry_dict = {
                    'entry_id': entry.entry_id,
                    'e_above_hull': e_above_hull,
                    'metadata': metadata,
                    'formula': entry.composition.reduced_formula
                }
                if on_hull:
                    new_stable.append(entry_dict)
                elif e_above_hull <= args.e_above_hull_max:
                    new_meta_all.append(entry_dict)
        
        # For metastable, keep only lowest energy per composition
        comp_groups = {}
        for d in new_meta_all:
            comp_key = d['formula']
            if comp_key not in comp_groups:
                comp_groups[comp_key] = d
            else:
                if d['e_above_hull'] < comp_groups[comp_key]['e_above_hull']:
                    comp_groups[comp_key] = d
        new_meta = list(comp_groups.values())
        
        # Prompt user for confirmed electrides (only once)
        confirmed_ids = prompt_for_confirmed_electrides(new_stable, new_meta)
        
        # Plot
        plot_ternary_hull(system, entries, entry_metadata, elements, args.output_dir, 
                         args.e_above_hull_max, confirmed_ids=confirmed_ids)
    
    print("\n" + "="*70)
    print("All plots completed!")
    print(f"Output directory: {args.output_dir.resolve()}")


if __name__ == '__main__':
    main()
