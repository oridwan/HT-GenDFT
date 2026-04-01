#!/usr/bin/env python3
"""
Plot binary convex hull diagrams including newly discovered stable electrides.

This script plots formation energy vs composition for binary systems, showing:
- The convex hull (connected line)
- Stable phases (on the hull)
- Metastable phases (within threshold above hull)
- Newly discovered electrides highlighted

Usage:
    python3 binary_chull_plot.py --csv Bin-Ele-HT/stable_electrides.csv \
        --dft-results Bin-Ele-HT/dft_stability_results.json \
        --mp-phases Bin-Ele-HT/mp_vaspdft.json \
        --systems Ca-P Y-N

    python3 binary_chull_plot.py --db Bin-Ele-HT/stable_electrides.db \
        --dft-results Bin-Ele-HT/dft_stability_results.json \
        --mp-phases Bin-Ele-HT/mp_vaspdft.json \
        --systems Ca-P Y-N \
        --output-dir hull_plots \
        --e-above-hull-max 0.05
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

ELECTRONEGATIVE_ELEMENTS = {'N', 'P', 'As', 'Sb', 'Bi', 'O', 'S', 'Se', 'Te', 'Po', 'F', 'Cl', 'Br', 'I', 'At', 'H'}


def formula_to_latex(formula: str) -> str:
    """Convert chemical formula to LaTeX format with subscripts.
    
    Handles both simple subscripts (Ca3 -> Ca$_{3}$) and 
    parentheses subscripts (Ca3(AlP2)2 -> Ca$_{3}$(AlP$_{2}$)$_{2}$).
    Ignores subscript 1.
    
    Examples:
        Ca5P3 -> Ca$_{5}$P$_{3}$
        Ca3P1 -> Ca$_{3}$P
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


def calculate_nexcess_boundaries_binary(elem_A: str, elem_C: str, n_excess_max: float = 4.0, max_atoms: int = 20) -> Tuple[float, float]:
    """
    Calculate exact composition boundaries for N_excess region in binary A-C system.
    This matches exactly the logic in search_binary_electrides.py.
    
    Iterates through all reduced integer combinations (l, n) with l+n <= max_atoms
    and finds compositions where 0 < N_excess <= n_excess_max.
    
    For binary system A_l C_n:
    - N_excess = val_A * l - val_C * n
    - x_C = n / (l + n) (fraction of element C)
    
    Args:
        elem_A: Symbol of element A (electropositive, usually metal)
        elem_C: Symbol of element C (electronegative, usually non-metal)
        n_excess_max: Maximum excess electrons (default 4 for binary)
        max_atoms: Maximum total atoms in formula (default 20 for binary)
    
    Returns:
        (x_C_min, x_C_max): Composition fraction boundaries for element C
    """
    from math import gcd
    
    if elem_A not in VALENCE_ELECTRONS or elem_C not in VALENCE_ELECTRONS:
        return None, None
    
    val_A = VALENCE_ELECTRONS[elem_A]
    val_C = VALENCE_ELECTRONS[elem_C]
    
    valid_x_C = []
    
    # Iterate through all possible stoichiometries
    for l in range(1, max_atoms):
        for n in range(1, max_atoms):
            if l + n > max_atoms:
                continue
            
            # Reduce to smallest integer ratio
            g = gcd(l, n)
            l_p, n_p = l // g, n // g
            
            # Skip if not reduced form
            if (l_p, n_p) != (l, n):
                continue
            
            # Calculate excess electrons
            n_excess = val_A * l_p - val_C * n_p
            
            # Check if in valid range
            if 0 < n_excess <= n_excess_max:
                x_C = n_p / (l_p + n_p)
                valid_x_C.append(x_C)
    
    if not valid_x_C:
        return None, None
    
    # Return the envelope (min and max x_C values)
    return min(valid_x_C), max(valid_x_C)


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
    print("  Example: Ca5P3_s021,Ca2P_s055")
    
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
        phase_data: List of phase data dicts with 'x', 'y', 'formula', 'entry_id'
        phase_type: Description of phase type (e.g., 'MP stable', 'New stable electride')
    
    Returns:
        List of (x_offset, y_offset) tuples in points (matplotlib offset points)
    """
    print(f"\n{'='*70}")
    print(f"INTERACTIVE LABEL POSITIONING: {phase_type}")
    print(f"{'='*70}")
    print(f"Enter label position offsets (in points) for each phase.")
    print(f"Format: x_offset,y_offset (e.g., 30,40 or -20,30)")
    print(f"Default: 0,30 (above the marker)")
    print(f"{'='*70}\n")
    
    positions = []
    for i, d in enumerate(phase_data, 1):
        formula = d.get('formula', d.get('entry_id', 'Unknown'))
        entry_id = d.get('entry_id', '')
        x_pos = d.get('x', d.get('cart_x', 0))
        y_pos = d.get('y', d.get('cart_y', 0))
        
        print(f"{i:2d}. {formula:15s} (ID: {entry_id:20s}) at ({x_pos:.3f}, {y_pos:.3f})")
        user_input = input(f"    Label offset [default: 0,30]: ").strip()
        
        if not user_input:
            positions.append((0.0, 30.0))
        else:
            try:
                parts = user_input.split(',')
                x_off = float(parts[0].strip())
                y_off = float(parts[1].strip())
                positions.append((x_off, y_off))
            except (ValueError, IndexError):
                print(f"    Invalid input, using default (0, 30)")
                positions.append((0.0, 30.0))
    
    print(f"\n{'='*70}\n")
    return positions


def load_electride_candidates_csv(csv_path: Path, systems: List[str], e_hull_max: float) -> pd.DataFrame:
    """Load electride candidates from CSV file filtered by chemical systems.
    
    Note: We load candidates with original dft_e_hull <= e_hull_max (relative to MP phases).
    The actual e_above_hull will be recalculated with all entries.
    """
    df = pd.read_csv(csv_path)
    
    # Filter by e_above_hull threshold (candidates that were close to hull relative to MP)
    df = df[df['dft_e_hull'] <= e_hull_max].copy()
    
    # Normalize systems to alphabetical order
    normalized_systems = ['-'.join(sorted(s.split('-'))) for s in systems]
    
    # Filter by chemical systems
    filtered_rows = []
    for _, row in df.iterrows():
        comp = Composition(row['composition'])
        elements = sorted([str(e) for e in comp.elements])
        chemsys = '-'.join(elements)
        
        if chemsys in normalized_systems:
            # Ensure we have structure_id and spacegroup
            if 'formula' not in row:
                row['formula'] = row['composition']
            filtered_rows.append(row)
    
    return pd.DataFrame(filtered_rows)


def load_electride_candidates_db(db_path: Path, systems: List[str], e_hull_max: float) -> pd.DataFrame:
    """Load electride candidates from ASE database filtered by chemical systems.
    
    Note: We load candidates with original dft_e_hull <= e_hull_max (relative to MP phases).
    The actual e_above_hull will be recalculated with all entries.
    """
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
    
    return pd.DataFrame(rows)


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
    electride_candidates: pd.DataFrame,
    dft_results: Dict,
    mp_phases: List[Dict]
) -> Tuple[List[ComputedEntry], Dict[str, Dict], List[str]]:
    """Create ComputedEntry objects for a given chemical system.
    
    Returns:
        entries: All entries (MP + electride candidates)
        entry_metadata: Dict mapping entry_id to {pearson_symbol, spacegroup, mp_id, etc.}
        elements: List of elements in the system [elem1, elem2]
    """
    
    entries = []
    entry_metadata = {}
    temp_electride_data = []  # Temporary storage for deduplication
    
    # Normalize system to alphabetical order
    system_normalized = '-'.join(sorted(system.split('-')))
    
    # Get elements in this system (alphabetically sorted)
    elements = sorted(system_normalized.split('-'))
    
    # Collect electride candidates for deduplication
    for _, row in electride_candidates.iterrows():
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
    
    # Deduplicate electrides: same composition, same spacegroup, E_form within 0.05 eV/atom
    deduplicated_electrides = []
    processed = set()
    
    for i, data1 in enumerate(temp_electride_data):
        if i in processed:
            continue
        
        duplicates = [data1]
        for j, data2 in enumerate(temp_electride_data):
            if j != i and j not in processed:
                # Same composition (reduced formula)
                same_comp = data1['composition'].reduced_formula == data2['composition'].reduced_formula
                # Same spacegroup (if available)
                same_sg = (data1['spacegroup'] is not None and 
                          data2['spacegroup'] is not None and 
                          data1['spacegroup'] == data2['spacegroup'])
                # Energy difference within 0.05 eV/atom
                e_diff = abs(data1['energy_per_atom'] - data2['energy_per_atom'])
                
                if same_comp and same_sg and e_diff <= 0.05:
                    duplicates.append(data2)
                    processed.add(j)
        
        # Keep the one with lowest energy
        best = min(duplicates, key=lambda d: d['energy_per_atom'])
        deduplicated_electrides.append(best)
        processed.add(i)
    
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
            'pearson_symbol': None,  # Will try to get from database
            'formula': data['formula']
        }
    
    # Add MP reference phases (including elemental entries)
    for phase in mp_phases:
        phase_chemsys = phase['chemsys']
        comp_dict = phase['composition']
        comp = Composition(comp_dict)
        
        # Normalize phase chemsys
        phase_chemsys_normalized = '-'.join(sorted(phase_chemsys.split('-')))
        
        # Include if it matches the system or is an elemental entry for one of the system elements
        phase_elements = list(comp_dict.keys())
        is_system_phase = phase_chemsys_normalized == system_normalized
        is_elemental = len(phase_elements) == 1 and phase_elements[0] in elements
        
        if is_system_phase or is_elemental:
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


def select_ground_state_phases(
    entries: List[ComputedEntry],
    entry_metadata: Dict[str, Dict],
    elements: List[str],
    e_hull_max: float
) -> Tuple[ComputedEntry, ComputedEntry, float, float]:
    """Select two ground state phases that bracket the electride composition range.
    
    Returns:
        ground_state_1: Lower composition ground state (elemental or compound)
        ground_state_2: Higher composition ground state (elemental or compound)
        x1: Composition of ground_state_1 (fraction of elem2)
        x2: Composition of ground_state_2 (fraction of elem2)
    """
    # Get electride entries
    electride_entries = [e for e in entries if entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
    
    if not electride_entries:
        # No electrides, use elemental phases
        elem1_entry = None
        elem2_entry = None
        for entry in entries:
            if len(entry.composition.elements) == 1:
                if str(entry.composition.elements[0]) == elements[0]:
                    elem1_entry = entry
                elif str(entry.composition.elements[0]) == elements[1]:
                    elem2_entry = entry
        return elem1_entry, elem2_entry, 0.0, 1.0
    
    # Calculate composition (x = fraction of elem2) for all electrides
    electride_compositions = []
    for entry in electride_entries:
        comp = entry.composition
        x = comp.get_atomic_fraction(elements[1])
        electride_compositions.append(x)
    
    # Find composition range containing electrides
    x_min = min(electride_compositions)
    x_max = max(electride_compositions)
    
    # Get stable MP phases
    mp_entries = [e for e in entries if not entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
    pd_mp = PhaseDiagram(mp_entries)
    stable_mp_entries = list(pd_mp.stable_entries)
    
    # Find bracketing phases
    ground_state_1 = None
    ground_state_2 = None
    x1 = None
    x2 = None
    
    # Sort stable MP entries by composition
    stable_with_comp = []
    for entry in stable_mp_entries:
        comp = entry.composition
        comp_elems = sorted([str(e) for e in comp.elements])
        if len(comp_elems) <= 2 and all(e in elements for e in comp_elems):
            x = comp.get_atomic_fraction(elements[1])
            stable_with_comp.append((x, entry))
    
    stable_with_comp.sort(key=lambda t: t[0])
    
    # Find ground state 1: closest stable phase with x <= x_min
    for x, entry in stable_with_comp:
        if x <= x_min + 0.01:  # Small tolerance
            ground_state_1 = entry
            x1 = x
    
    # Find ground state 2: closest stable phase with x >= x_max
    for x, entry in reversed(stable_with_comp):
        if x >= x_max - 0.01:  # Small tolerance
            ground_state_2 = entry
            x2 = x
    
    # Fallback to elemental phases if not found
    if ground_state_1 is None:
        for x, entry in stable_with_comp:
            if len(entry.composition.elements) == 1 and str(entry.composition.elements[0]) == elements[0]:
                ground_state_1 = entry
                x1 = x
                break
    
    if ground_state_2 is None:
        for x, entry in reversed(stable_with_comp):
            if len(entry.composition.elements) == 1 and str(entry.composition.elements[0]) == elements[1]:
                ground_state_2 = entry
                x2 = x
                break
    
    return ground_state_1, ground_state_2, x1, x2


def is_metal(element_symbol: str) -> bool:
    """Check if an element is a metal."""
    from pymatgen.core import Element
    try:
        elem = Element(element_symbol)
        return elem.is_metal or elem.is_alkali or elem.is_alkaline or elem.is_transition_metal or elem.is_post_transition_metal or elem.is_rare_earth_metal or elem.is_actinoid
    except:
        return False


def plot_binary_hull_combined(
    system: str,
    entries: List[ComputedEntry],
    entry_metadata: Dict[str, Dict],
    elements: List[str],
    output_dir: Path,
    e_hull_max: float = 0.05,
    confirmed_ids: List[str] = None
):
    """Plot combined binary convex hull with zoomed and full views.
    
    Args:
        system: Chemical system (e.g., 'Ca-P')
        entries: All ComputedEntry objects
        entry_metadata: Metadata dict for entries (pearson_symbol, mp_id, etc.)
        elements: [elem1, elem2] sorted alphabetically
        output_dir: Output directory
        e_hull_max: Maximum energy above hull to display (eV/atom)
        confirmed_ids: List of confirmed electride structure IDs
    """
    
    # Ensure metal is on the left (x=0)
    if is_metal(elements[1]) and not is_metal(elements[0]):
        elements = [elements[1], elements[0]]
    
    if confirmed_ids is None:
        confirmed_ids = []
    
    # Prepare phase data for both full and zoomed plots
    # We need to identify which phases will be plotted in each view
    full_mp_stable, _, _ = _prepare_phase_data(
        entries, entry_metadata, elements, e_hull_max, plot_full_hull=True
    )
    zoom_mp_stable, zoom_new_stable, zoom_new_meta = _prepare_phase_data(
        entries, entry_metadata, elements, e_hull_max, plot_full_hull=False
    )
    
    # Prompt for all label positions upfront in the correct order
    print("\n" + "="*70)
    print("LABEL POSITION SETUP")
    print("="*70)
    print("You will be prompted for label positions in the following order:")
    print("  1. MP stable phases (full hull plot)")
    print("  2. MP stable phases (zoomed-in hull plot)")
    print("  3. New stable structures (zoomed-in hull plot)")
    print("  4. Metastable structures (zoomed-in hull plot)")
    print("")
    print("Note: Full hull plot only shows MP stable phase labels.")
    print("      New stable/metastable structures are labeled only in zoomed plot.")
    print("="*70)
    
    # 1. Full hull - MP stable
    full_mp_positions = []
    if full_mp_stable:
        full_mp_positions = prompt_for_label_positions(full_mp_stable, "FULL HULL - MP stable phases")
    
    # Full hull does not show labels for new stable/metastable structures
    full_new_stable_positions = []
    full_new_meta_positions = []
    
    # 2. Zoomed hull - MP stable
    zoom_mp_positions = []
    if zoom_mp_stable:
        zoom_mp_positions = prompt_for_label_positions(zoom_mp_stable, "ZOOMED HULL - MP stable phases")
    
    # 3. Zoomed hull - New stable
    zoom_new_stable_positions = []
    if zoom_new_stable:
        zoom_new_stable_positions = prompt_for_label_positions(zoom_new_stable, "ZOOMED HULL - New stable structures")
    
    # 4. Zoomed hull - Metastable
    zoom_new_meta_positions = []
    if zoom_new_meta:
        zoom_new_meta_positions = prompt_for_label_positions(zoom_new_meta, "ZOOMED HULL - Metastable structures")
    
    # Create figure with 2 subplots (top: full, bottom: zoomed)
    fig, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(12, 16), 
                                            gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.15, 'top': 0.96})
    
    # Plot full hull on top axis
    _plot_single_hull(ax_full, entries, entry_metadata, elements, e_hull_max, 
                     confirmed_ids, plot_full_hull=True, y_max=0.1,
                     mp_label_positions=full_mp_positions,
                     new_stable_positions=full_new_stable_positions,
                     new_meta_positions=full_new_meta_positions)
    
    # Plot zoomed hull on bottom axis  
    _plot_single_hull(ax_zoom, entries, entry_metadata, elements, e_hull_max,
                     confirmed_ids, plot_full_hull=False, y_max=0.005,
                     mp_label_positions=zoom_mp_positions,
                     new_stable_positions=zoom_new_stable_positions,
                     new_meta_positions=zoom_new_meta_positions)
    
    # Set single title at top (only show energy range if there are metastable structures)
    if zoom_new_meta:
        title = f'{system} ($E_{{\\mathrm{{hull}}}} \\leq {e_hull_max:.2f}$ eV/atom)'
    else:
        title = f'{system}'
    fig.suptitle(title, fontsize=24, fontweight='bold', y=0.99)
    
    # Remove individual subplot titles
    ax_full.set_title('')
    ax_zoom.set_title('')
    
    # Save figure
    output_path = output_dir / f'{system.replace("-", "")}_hull.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    
    # Also save as PDF
    output_pdf = output_dir / f'{system.replace("-", "")}_hull.pdf'
    fig.savefig(output_pdf, bbox_inches='tight')
    print(f"Saved: {output_pdf}")
    
    plt.close(fig)


def _prepare_phase_data(
    entries: List[ComputedEntry],
    entry_metadata: Dict[str, Dict],
    elements: List[str],
    e_hull_max: float,
    plot_full_hull: bool
) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    """Prepare phase data for plotting to identify which phases need labels.
    
    Returns:
        (mp_stable, new_stable, new_meta) - Lists of phase data dicts
    """
    # Separate MP entries from electride entries
    mp_entries = [e for e in entries if not entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
    
    # Create phase diagrams
    pd_mp_only = PhaseDiagram(mp_entries)
    pd_all = PhaseDiagram(entries)
    
    # Determine ground states based on plot type
    if plot_full_hull:
        # Use elemental phases for full hull plot
        elem1_entry = None
        elem2_entry = None
        for entry in mp_entries:
            if len(entry.composition.elements) == 1:
                if str(entry.composition.elements[0]) == elements[0]:
                    elem1_entry = entry
                elif str(entry.composition.elements[0]) == elements[1]:
                    elem2_entry = entry
        gs1, gs2, x1, x2 = elem1_entry, elem2_entry, 0.0, 1.0
    else:
        # Smart zoom to electride region
        gs1, gs2, x1, x2 = select_ground_state_phases(entries, entry_metadata, elements, e_hull_max)
    
    # Get energies of ground states
    gs1_energy = gs1.energy / gs1.composition.num_atoms
    gs2_energy = gs2.energy / gs2.composition.num_atoms
    
    # Prepare data
    all_data = []
    for entry in entries:
        comp = entry.composition
        comp_elems = sorted([str(e) for e in comp.elements])
        if len(comp_elems) > 2:
            continue
        
        x_original = comp.get_atomic_fraction(elements[1])
        if x_original < x1 - 0.01 or x_original > x2 + 0.01:
            continue
        
        if abs(x2 - x1) < 1e-6:
            x_norm = 0.5
        else:
            x_norm = (x_original - x1) / (x2 - x1)
        
        if plot_full_hull:
            formation_energy = pd_all.get_form_energy_per_atom(entry)
        else:
            energy_per_atom = entry.energy / comp.num_atoms
            formation_energy = energy_per_atom - (x_norm * gs2_energy + (1 - x_norm) * gs1_energy)
        
        e_above_hull = pd_all.get_e_above_hull(entry)
        metadata = entry_metadata.get(entry.entry_id, {})
        is_electride = metadata.get('is_electride', False)
        on_hull = entry in pd_all.stable_entries
        on_mp_hull = entry in pd_mp_only.stable_entries
        
        all_data.append({
            'x': x_norm,
            'y': formation_energy,
            'e_above_hull': e_above_hull,
            'on_hull': on_hull,
            'on_mp_hull': on_mp_hull,
            'is_electride': is_electride,
            'formula': comp.reduced_formula,
            'entry_id': entry.entry_id,
            'metadata': metadata
        })
    
    all_data = sorted(all_data, key=lambda d: d['x'])
    
    # Identify phases
    hull_phases = [d for d in all_data if d['on_hull']]
    metastable_phases = [d for d in all_data if not d['on_hull'] and d['e_above_hull'] <= e_hull_max]
    mp_stable = [d for d in hull_phases if not d['is_electride']]
    new_stable = [d for d in hull_phases if d['is_electride']]
    new_meta = [d for d in metastable_phases if d['is_electride']]
    
    return mp_stable, new_stable, new_meta


def _plot_single_hull(
    ax,
    entries: List[ComputedEntry],
    entry_metadata: Dict[str, Dict],
    elements: List[str],
    e_hull_max: float,
    confirmed_ids: List[str],
    plot_full_hull: bool,
    y_max: float,
    mp_label_positions: List[Tuple[float, float]] = None,
    new_stable_positions: List[Tuple[float, float]] = None,
    new_meta_positions: List[Tuple[float, float]] = None
):
    """Plot a single binary hull (either full or zoomed).
    
    Args:
        ax: Matplotlib axes
        entries: All ComputedEntry objects
        entry_metadata: Metadata dict
        elements: [elem1, elem2] with metal on left
        e_hull_max: Max energy above hull to display
        confirmed_ids: List of confirmed electride IDs
        plot_full_hull: If True, plot full hull; if False, plot zoomed
        y_max: Maximum y-axis limit
        mp_label_positions: Pre-collected label positions for MP stable phases
        new_stable_positions: Pre-collected label positions for new stable structures
        new_meta_positions: Pre-collected label positions for metastable structures
    """
    
    # Separate MP entries from electride entries
    mp_entries = [e for e in entries if not entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
    
    # Create two phase diagrams: MP-only and all entries
    pd_mp_only = PhaseDiagram(mp_entries)
    pd_all = PhaseDiagram(entries)
    
    if confirmed_ids is None:
        confirmed_ids = []
    
    # Determine ground states and composition range
    if plot_full_hull:
        # Use elemental phases for full hull plot
        # Find elemental MP entries
        elem1_entry = None
        elem2_entry = None
        for entry in mp_entries:  # Only look in MP entries
            if len(entry.composition.elements) == 1:
                if str(entry.composition.elements[0]) == elements[0]:
                    elem1_entry = entry
                elif str(entry.composition.elements[0]) == elements[1]:
                    elem2_entry = entry
        gs1, gs2, x1, x2 = elem1_entry, elem2_entry, 0.0, 1.0
    else:
        # Smart zoom to electride region
        gs1, gs2, x1, x2 = select_ground_state_phases(entries, entry_metadata, elements, e_hull_max)
    
    # Get energies of ground states (per atom)
    gs1_energy = gs1.energy / gs1.composition.num_atoms
    gs2_energy = gs2.energy / gs2.composition.num_atoms
    
    # Get formulas for ground states
    if plot_full_hull:
        # For full hull, use element symbols only (no molecular formulas like N2)
        gs1_formula = str(gs1.composition.elements[0])
        gs2_formula = str(gs2.composition.elements[0])
    else:
        # For smart-zoom, use the actual formula of the ground state phases
        gs1_formula = gs1.composition.reduced_formula
        gs2_formula = gs2.composition.reduced_formula
        
        # Add parentheses around non-elemental (compound) formulas
        if len(gs1.composition.elements) > 1:
            gs1_formula = f"({gs1_formula})"
        if len(gs2.composition.elements) > 1:
            gs2_formula = f"({gs2_formula})"
    
    # Prepare data for plotting
    all_data = []
    for entry in entries:
        comp = entry.composition
        
        # Skip if not in this binary system
        comp_elems = sorted([str(e) for e in comp.elements])
        if len(comp_elems) > 2:
            continue
        
        # Calculate original composition fraction: fraction of elem2
        x_original = comp.get_atomic_fraction(elements[1])
        
        # Skip if outside ground state range
        if x_original < x1 - 0.01 or x_original > x2 + 0.01:
            continue
        
        # Normalize composition between ground states
        if abs(x2 - x1) < 1e-6:
            x_norm = 0.5
        else:
            x_norm = (x_original - x1) / (x2 - x1)
        
        # Calculate formation energy relative to ground states
        if plot_full_hull:
            # For full hull, use pymatgen's formation energy (relative to elemental phases)
            # This ensures elemental phases have exactly 0.0 formation energy
            formation_energy = pd_all.get_form_energy_per_atom(entry)
        else:
            # For smart-zoom, use manual calculation relative to selected ground states
            energy_per_atom = entry.energy / comp.num_atoms
            formation_energy = energy_per_atom - (x_norm * gs2_energy + (1 - x_norm) * gs1_energy)
        
        # Get energy above hull (using full hull with electrides)
        e_above_hull = pd_all.get_e_above_hull(entry)
        
        # Get metadata
        metadata = entry_metadata.get(entry.entry_id, {})
        is_electride = metadata.get('is_electride', False)
        
        # Get energy above MP-only hull (for electrides only, for comparison)
        # Negative value means the electride is below the original MP hull
        if is_electride:
            try:
                e_above_mp_hull = pd_mp_only.get_e_above_hull(entry)
            except ValueError:
                # Entry not in MP diagram - calculate manually
                # Use pymatgen's formation energy relative to elemental phases
                entry_form_energy = pd_all.get_form_energy_per_atom(entry)
                # Get MP decomposition
                decomp = pd_mp_only.get_decomposition(entry.composition)
                ref_energy = sum(amt * pd_mp_only.get_form_energy_per_atom(e) 
                                for e, amt in decomp.items())
                e_above_mp_hull = entry_form_energy - ref_energy
        else:
            e_above_mp_hull = 0.0  # Not applicable for MP phases
        
        # Check if on hull (using full hull)
        on_hull = entry in pd_all.stable_entries
        
        # Check if on MP-only hull
        on_mp_hull = entry in pd_mp_only.stable_entries
        
        all_data.append({
            'x': x_norm,
            'y': formation_energy,
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
    
    # Sort by x for plotting
    all_data = sorted(all_data, key=lambda d: d['x'])
    
    # Identify electride candidates
    hull_phases = [d for d in all_data if d['on_hull']]
    metastable_phases = [d for d in all_data if not d['on_hull'] and d['e_above_hull'] <= e_hull_max]
    new_stable = [d for d in hull_phases if d['is_electride']]
    new_meta = [d for d in metastable_phases if d['is_electride']]
    
    # Highlight N_excess region (only in full hull plot)
    if plot_full_hull:
        x_C_lower, x_C_upper = calculate_nexcess_boundaries_binary(
            elements[0], elements[1], 
            n_excess_max=4.0, 
            max_atoms=20
        )
        if x_C_lower is not None and x_C_upper is not None and x_C_lower < x_C_upper:
            ax.axvspan(x_C_lower, x_C_upper, alpha=0.35, color='gold', zorder=0, 
                      label='')
    
    # Plot MP-only convex hull (gray dashed line)
    mp_hull_data = [d for d in all_data if d['on_mp_hull']]
    mp_hull_x = [d['x'] for d in mp_hull_data]
    mp_hull_y = [d['y'] for d in mp_hull_data]
    if len(mp_hull_x) >= 2:
        ax.plot(mp_hull_x, mp_hull_y, 'gray', linestyle='--', linewidth=2.5, 
               alpha=0.6, zorder=1, label='MP-only hull')
    
    # Plot full convex hull line (black solid)
    hull_x = [d['x'] for d in hull_phases]
    hull_y = [d['y'] for d in hull_phases]
    ax.plot(hull_x, hull_y, 'black', linestyle='-', linewidth=3, zorder=2, label='Convex hull')
    
    # Get MP stable phases (positions should be pre-collected)
    mp_stable = [d for d in hull_phases if not d['is_electride']]
    if mp_stable:
        # Use pre-collected label positions
        if mp_label_positions is None:
            mp_label_positions = [(0.0, 30.0)] * len(mp_stable)
        
        # Plot MP stable phases (black fill, no edge)
        mp_x = [d['x'] for d in mp_stable]
        mp_y = [d['y'] for d in mp_stable]
        ax.scatter(mp_x, mp_y, c='black', s=150, marker='o', 
                  edgecolors='none', linewidth=0, zorder=3, label='MP stable')
        
        # Label MP stable phases with LaTeX formula (no MP ID)
        for i, d in enumerate(mp_stable):
            formula_latex = formula_to_latex(d['formula'])
            x_offset, y_offset = mp_label_positions[i]
            
            ax.annotate(formula_latex, xy=(d['x'], d['y']), 
                       xytext=(x_offset, y_offset), textcoords='offset points',
                       fontsize=18, ha='center', color='black', fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.4', facecolor='lightgray', 
                               alpha=0.8, edgecolor='black', linewidth=1.2),
                       arrowprops=dict(arrowstyle='-', color='black', lw=1.2, alpha=0.7))
    
    # Separate new stable structures by electride vs non-electride
    new_stable = [d for d in hull_phases if d['is_electride']]
    new_meta = [d for d in metastable_phases if d['is_electride']]
    
    # Plot stable new structures (electride and non-electride)
    if new_stable:
        # Use pre-collected label positions
        if new_stable_positions is None:
            new_stable_positions = [(0.0, 30.0)] * len(new_stable)
        
        # Plot markers (triangle for confirmed electride, square for candidate/non-electride, blue fill, red edge)
        for d in new_stable:
            is_confirmed = d['entry_id'] in confirmed_ids
            marker = '^' if is_confirmed else 's'  # Triangle for confirmed electride, square for candidate
            facecolor = 'blue'  # Always blue fill for stable structures
            marker_size = 180 if is_confirmed else 150
            ax.scatter(d['x'], d['y'], c=facecolor, s=marker_size, marker=marker, 
                      edgecolors='red', linewidth=2.5, zorder=5)
        
        # Label new stable structures (only in zoomed plot, not full hull)
        if not plot_full_hull and new_stable_positions:
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
                
                ax.annotate(label_text, xy=(d['x'], d['y']), 
                           xytext=(x_offset, y_offset), textcoords='offset points',
                           fontsize=18, ha='center', color='blue', fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.4', facecolor='lightblue', 
                                   alpha=0.9, edgecolor='blue', linewidth=2),
                           arrowprops=dict(arrowstyle='-', color='blue', lw=1.5, alpha=0.7))
    
    # Plot metastable new structures (electride and non-electride)
    if new_meta:
        # Use pre-collected label positions
        if new_meta_positions is None:
            new_meta_positions = [(0.0, 30.0)] * len(new_meta)
        
        # Plot markers (triangle for confirmed electride, square for candidate/non-electride, orange fill, red edge)
        for d in new_meta:
            is_confirmed = d['entry_id'] in confirmed_ids
            marker = '^' if is_confirmed else 's'  # Triangle for confirmed electride, square for candidate
            facecolor = 'orange'  # Always orange fill for metastable structures
            marker_size = 120 if is_confirmed else 100
            ax.scatter(d['x'], d['y'], c=facecolor, s=marker_size, marker=marker, 
                      edgecolors='red', linewidth=2, zorder=4)
        
        # Label metastable structures (only in zoomed plot, not full hull)
        if not plot_full_hull and new_meta_positions:
            for i, d in enumerate(new_meta):
                pearson = d['metadata'].get('pearson_symbol', '')
                structure_id = d['entry_id']
                formula_part = structure_id.split('_')[0] if '_' in structure_id else d['formula']
                formula_latex = formula_to_latex(formula_part)
                e_hull_val = d['e_above_hull']
                if pearson:
                    label_text = f"{pearson}-{formula_latex}"
                else:
                    label_text = f"{formula_part}"
                
                x_offset, y_offset = new_meta_positions[i]
                
                ax.annotate(label_text, xy=(d['x'], d['y']), 
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
            ax.scatter([], [], facecolors='none', s=180, marker='^', 
                      edgecolors='red', linewidth=2.5, label='Electride')
        if candidate_count > 0:
            ax.scatter([], [], facecolors='none', s=150, marker='s', 
                      edgecolors='red', linewidth=2, label='Not electride')
    
    # Configure axes
    ax.set_xlabel(f'Composition: x in {gs1_formula}$_{{1-x}}${gs2_formula}$_x$', fontsize=20, fontweight='bold')
    ax.set_ylabel('Formation Energy (eV/atom)', fontsize=20, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='lower left', fontsize=18, framealpha=0.0, edgecolor='none')
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # Set x-axis limits with some padding
    ax.set_xlim(-0.05, 1.05)
    
    # Set y-axis limits with dynamic range
    y_min = min([d['y'] for d in all_data])
    y_max_data = max([d['y'] for d in all_data if d['on_hull'] or d['e_above_hull'] <= e_hull_max])
    y_range = y_max_data - y_min
    ax.set_ylim(y_min - 0.2 * y_range, min(y_max, y_max_data) + 0.2 * y_range)


def main():
    parser = argparse.ArgumentParser(
        description='Plot binary convex hull diagrams with newly discovered electrides',
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
                       help='Chemical systems to plot (e.g., Ca-P Y-N)')
    
    # Output
    parser.add_argument('--output-dir', type=Path, default=Path('.'),
                       help='Output directory for plots (default: current directory)')
    
    # Energy above hull threshold
    parser.add_argument('--e-above-hull-max', type=float, default=0.05,
                       help='Maximum energy above hull to display (eV/atom, default: 0.05)')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load electride candidates (not just stable ones, but those close to hull)
    print("Loading electride candidates...")
    if args.csv:
        electride_candidates = load_electride_candidates_csv(args.csv, args.systems, args.e_above_hull_max)
        print(f"  Loaded {len(electride_candidates)} electride candidates from CSV")
        # Get Pearson symbols from database if it exists
        db_path = args.csv.parent / 'stable_electrides.db'
        if db_path.exists():
            structure_ids = electride_candidates['structure_id'].tolist()
            pearson_map = get_pearson_symbols_from_db(db_path, structure_ids)
            print(f"  Loaded Pearson symbols for {len(pearson_map)} structures from database")
        else:
            pearson_map = {}
    else:
        electride_candidates = load_electride_candidates_db(args.db, args.systems, args.e_above_hull_max)
        print(f"  Loaded {len(electride_candidates)} electride candidates from database")
        # Get Pearson symbols from database
        structure_ids = electride_candidates['structure_id'].tolist()
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
        
        # Create phase diagrams to identify stable and metastable candidates
        mp_entries = [e for e in entries if not entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
        pd_all = PhaseDiagram(entries)
        
        # Collect all electride candidates for confirmation prompt
        new_stable = []
        new_meta = []
        for entry in entries:
            metadata = entry_metadata.get(entry.entry_id, {})
            if metadata.get('is_electride', False):
                e_above_hull = pd_all.get_e_above_hull(entry)
                on_hull = entry in pd_all.stable_entries
                entry_dict = {
                    'entry_id': entry.entry_id,
                    'e_above_hull': e_above_hull,
                    'metadata': metadata
                }
                if on_hull:
                    new_stable.append(entry_dict)
                elif e_above_hull <= args.e_above_hull_max:
                    new_meta.append(entry_dict)
        
        # Prompt user for confirmed electrides (only once)
        confirmed_ids = prompt_for_confirmed_electrides(new_stable, new_meta)
        
        # Plot combined hull (full on top, zoomed on bottom)
        print("\n--- Plotting combined hull ---")
        plot_binary_hull_combined(system, entries, entry_metadata, elements, args.output_dir, 
                                  args.e_above_hull_max, confirmed_ids=confirmed_ids)
        
        # Print summary
        electride_data = [m for m in entry_metadata.values() if m.get('is_electride', False)]
        print(f"\n{system} System Summary:")
        print(f"  Total entries: {len(entries)}")
        print(f"  MP phases: {len(mp_entries)}")
        print(f"  Electride candidates: {len(electride_data)}")
        print(f"  New stable electrides: {len(new_stable)}")
        print(f"  New metastable electrides (E_hull ≤ {args.e_above_hull_max}): {len(new_meta)}")
        print(f"  Confirmed electrides: {len(confirmed_ids)}")
        
        if new_stable:
            print(f"\n  Stable electrides:")
            for d in new_stable:
                pearson = d['metadata'].get('pearson_symbol', 'N/A')
                confirmed = " [CONFIRMED]" if d['entry_id'] in confirmed_ids else ""
                print(f"    - {pearson} ({d['entry_id']}): E_hull = {d['e_above_hull']:.4f} eV/atom{confirmed}")
        
        if new_meta:
            print(f"\n  Metastable electrides (top 10 by stability):")
            sorted_meta = sorted(new_meta, key=lambda d: d['e_above_hull'])[:10]
            for d in sorted_meta:
                pearson = d['metadata'].get('pearson_symbol', 'N/A')
                confirmed = " [CONFIRMED]" if d['entry_id'] in confirmed_ids else ""
                print(f"    - {pearson} ({d['entry_id']}): E_hull = {d['e_above_hull']:.4f} eV/atom{confirmed}")
    
    print("\n" + "="*70)
    print("All plots completed!")
    print(f"Output directory: {args.output_dir.resolve()}")


if __name__ == '__main__':
    main()
