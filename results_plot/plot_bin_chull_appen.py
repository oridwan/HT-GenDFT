#!/usr/bin/env python3
"""
Plot binary convex hull diagrams for appendix (no interactive labels).

This script plots formation energy vs composition for binary systems with:
- Full hull only (no zoomed view)
- MP phases shown in legend with different colors
- New structures shown as rectangles in legend
- Legend positioned outside plot for clarity

Usage:
    python3 plot_bin_chull_appen.py --db Bin-Ele-HT/stable_electrides.db \
        --dft-results Bin-Ele-HT/dft_stability_results.json \
        --mp-phases Bin-Ele-HT/mp_vaspdft.json \
        --systems Ca-P Y-N \
        --output-dir hull_plots_appen \
        --e-above-hull-max 0.05
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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
    """Calculate exact composition boundaries for N_excess region in binary A-C system."""
    from math import gcd
    
    if elem_A not in VALENCE_ELECTRONS or elem_C not in VALENCE_ELECTRONS:
        return None, None
    
    val_A = VALENCE_ELECTRONS[elem_A]
    val_C = VALENCE_ELECTRONS[elem_C]
    
    valid_x_C = []
    
    for l in range(1, max_atoms):
        for n in range(1, max_atoms):
            if l + n > max_atoms:
                continue
            
            g = gcd(l, n)
            l_p, n_p = l // g, n // g
            
            if (l_p, n_p) != (l, n):
                continue
            
            n_excess = val_A * l_p - val_C * n_p
            
            if 0 < n_excess <= n_excess_max:
                x_C = n_p / (l_p + n_p)
                valid_x_C.append(x_C)
    
    if not valid_x_C:
        return None, None
    
    return min(valid_x_C), max(valid_x_C)


def load_electride_candidates_db(db_path: Path, systems: List[str], e_hull_max: float) -> pd.DataFrame:
    """Load electride candidates from ASE database filtered by chemical systems."""
    from ase.db import connect
    
    db = connect(str(db_path))
    
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
    """Create ComputedEntry objects for a given chemical system."""
    
    entries = []
    entry_metadata = {}
    temp_electride_data = []
    
    system_normalized = '-'.join(sorted(system.split('-')))
    elements = sorted(system_normalized.split('-'))
    
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
    
    # Deduplicate electrides
    deduplicated_electrides = []
    processed = set()
    
    for i, data1 in enumerate(temp_electride_data):
        if i in processed:
            continue
        
        duplicates = [data1]
        for j, data2 in enumerate(temp_electride_data):
            if j != i and j not in processed:
                same_comp = data1['composition'].reduced_formula == data2['composition'].reduced_formula
                same_sg = (data1['spacegroup'] is not None and 
                          data2['spacegroup'] is not None and 
                          data1['spacegroup'] == data2['spacegroup'])
                e_diff = abs(data1['energy_per_atom'] - data2['energy_per_atom'])
                
                if same_comp and same_sg and e_diff <= 0.05:
                    duplicates.append(data2)
                    processed.add(j)
        
        best = min(duplicates, key=lambda d: d['energy_per_atom'])
        deduplicated_electrides.append(best)
        processed.add(i)
    
    for data in deduplicated_electrides:
        entry = ComputedEntry(
            composition=data['composition'],
            energy=data['total_energy'],
            entry_id=data['structure_id']
        )
        entries.append(entry)
        
        entry_metadata[data['structure_id']] = {
            'is_electride': True,
            'spacegroup': data['spacegroup'],
            'formula': data['formula'],
            'pearson_symbol': None
        }
    
    for phase in mp_phases:
        phase_chemsys = phase['chemsys']
        comp_dict = phase['composition']
        comp = Composition(comp_dict)
        
        phase_chemsys_normalized = '-'.join(sorted(phase_chemsys.split('-')))
        
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
            
            entry_metadata[entry_id] = {
                'is_electride': False,
                'mp_id': phase.get('mp_id', entry_id.split('-')[0] if '-' in entry_id else entry_id)
            }
    
    return entries, entry_metadata, elements


def is_metal(element_symbol: str) -> bool:
    """Check if an element is a metal."""
    from pymatgen.core import Element
    try:
        elem = Element(element_symbol)
        return elem.is_metal or elem.is_alkali or elem.is_alkaline or elem.is_transition_metal or elem.is_post_transition_metal or elem.is_rare_earth_metal or elem.is_actinoid
    except:
        return False


def plot_binary_hull_appendix(
    system: str,
    entries: List[ComputedEntry],
    entry_metadata: Dict[str, Dict],
    elements: List[str],
    output_dir: Path,
    e_hull_max: float = 0.05
):
    """Plot binary convex hull for appendix (no interactive labels)."""
    
    # Ensure metal is on the left (x=0)
    if is_metal(elements[1]) and not is_metal(elements[0]):
        elements = [elements[1], elements[0]]
    
    # Separate MP entries from electride entries
    mp_entries = [e for e in entries if not entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
    
    # Create phase diagrams
    pd_mp_only = PhaseDiagram(mp_entries)
    pd_all = PhaseDiagram(entries)
    
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
    
    # Get energies of ground states
    gs1_energy = gs1.energy / gs1.composition.num_atoms
    gs2_energy = gs2.energy / gs2.composition.num_atoms
    
    # Get formulas for ground states
    gs1_formula = str(gs1.composition.elements[0])
    gs2_formula = str(gs2.composition.elements[0])
    
    # Prepare data for plotting
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
        
        formation_energy = pd_all.get_form_energy_per_atom(entry)
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
            'entry': entry,
            'metadata': metadata
        })
    
    all_data = sorted(all_data, key=lambda d: d['x'])
    
    # Identify phases
    hull_phases = [d for d in all_data if d['on_hull']]
    metastable_phases = [d for d in all_data if not d['on_hull'] and d['e_above_hull'] <= e_hull_max]
    new_stable = [d for d in hull_phases if d['is_electride']]
    new_meta = [d for d in metastable_phases if d['is_electride']]
    mp_stable = [d for d in hull_phases if not d['is_electride']]
    
    # Create figure with fixed size
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Highlight N_excess region
    x_C_lower, x_C_upper = calculate_nexcess_boundaries_binary(
        elements[0], elements[1], 
        n_excess_max=4.0, 
        max_atoms=20
    )
    if x_C_lower is not None and x_C_upper is not None and x_C_lower < x_C_upper:
        ax.axvspan(x_C_lower, x_C_upper, alpha=0.35, color='gold', zorder=0)
    
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
    
    # Define colors for MP stable phases
    mp_colors = plt.cm.tab20(np.linspace(0, 1, max(len(mp_stable), 1)))
    
    # Plot MP stable phases with different colors
    for i, d in enumerate(mp_stable):
        color = mp_colors[i % len(mp_colors)]
        formula_latex = formula_to_latex(d['formula'])
        ax.scatter(d['x'], d['y'], c=[color], s=150, marker='o', 
                  edgecolors='black', linewidth=1.5, zorder=3, label=formula_latex)
    
    # Plot stable new structures (different marker for each)
    stable_markers = ['^', 'D', 'v', 'p', 'P', 'X', 'h', 'H', '<', '>']
    for i, d in enumerate(new_stable):
        marker = stable_markers[i % len(stable_markers)]
        ax.scatter(d['x'], d['y'], c='blue', s=150, marker=marker, 
                  edgecolors='red', linewidth=2.5, zorder=5)
    
    # Add legend entry for each stable structure
    if new_stable:
        for i, d in enumerate(new_stable):
            marker = stable_markers[i % len(stable_markers)]
            pearson = d['metadata'].get('pearson_symbol', '')
            formula_latex = formula_to_latex(d['formula'])
            if pearson:
                label = f'{pearson}-{formula_latex}'
            else:
                label = formula_latex
            ax.scatter([], [], c='blue', s=150, marker=marker, 
                      edgecolors='red', linewidth=2.5, label=label)
    
    # Plot metastable new structures (orange rectangles)
    for d in new_meta:
        ax.scatter(d['x'], d['y'], c='orange', s=120, marker='s', 
                  edgecolors='red', linewidth=2, zorder=4)
    
    # Metastable structures are plotted but not added to legend
    
    # Configure axes
    ax.set_xlabel(f'Composition: x in {gs1_formula}$_{{1-x}}${gs2_formula}$_x$', fontsize=20, fontweight='bold')
    ax.set_ylabel('Formation Energy (eV/atom)', fontsize=20, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(axis='both', which='major', labelsize=18)
    
    # Set x-axis limits
    ax.set_xlim(-0.05, 1.05)
    
    # Set y-axis limits with dynamic range
    y_min = min([d['y'] for d in all_data])
    y_max_data = max([d['y'] for d in all_data if d['on_hull'] or d['e_above_hull'] <= e_hull_max])
    y_range = y_max_data - y_min
    ax.set_ylim(y_min - 0.2 * y_range, y_max_data + 0.2 * y_range)
    
    # Set title
    if new_meta:
        title = f'{system} ($E_{{\\mathrm{{hull}}}} \\leq {e_hull_max:.2f}$ eV/atom)'
    else:
        title = f'{system}'
    ax.set_title(title, fontsize=24, fontweight='bold', pad=10)
    
    # Place legend at bottom right corner
    ax.legend(loc='lower right', fontsize=18, framealpha=0.65, edgecolor='black', facecolor='white')
    
    # Use fixed margins instead of tight_layout
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.10, top=0.92)
    
    # Save figure
    output_path = output_dir / f'{system.replace("-", "")}_hull.png'
    fig.savefig(output_path, dpi=300)
    print(f"Saved: {output_path}")
    
    output_pdf = output_dir / f'{system.replace("-", "")}_hull.pdf'
    fig.savefig(output_pdf)
    print(f"Saved: {output_pdf}")
    
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Plot binary convex hull diagrams for appendix',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--db', type=Path, required=True,
                       help='Path to stable_electrides.db')
    parser.add_argument('--dft-results', type=Path, required=True,
                       help='Path to dft_stability_results.json')
    parser.add_argument('--mp-phases', type=Path, required=True,
                       help='Path to mp_vaspdft.json')
    parser.add_argument('--systems', nargs='+', required=True,
                       help='Chemical systems to plot (e.g., Ca-P Y-N)')
    parser.add_argument('--output-dir', type=Path, default=Path('.'),
                       help='Output directory for plots (default: current directory)')
    parser.add_argument('--e-above-hull-max', type=float, default=0.05,
                       help='Maximum energy above hull to display (eV/atom, default: 0.05)')
    
    args = parser.parse_args()
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading electride candidates...")
    electride_candidates = load_electride_candidates_db(args.db, args.systems, args.e_above_hull_max)
    print(f"  Loaded {len(electride_candidates)} electride candidates from database")
    
    # Load Pearson symbols from database
    structure_ids = list(electride_candidates['structure_id'])
    pearson_map = get_pearson_symbols_from_db(args.db, structure_ids)
    print(f"  Loaded Pearson symbols for {len(pearson_map)} structures")
    
    print("\nLoading DFT results...")
    dft_results = load_dft_results(args.dft_results)
    print(f"  Loaded {len(dft_results)} DFT results")
    
    print("\nLoading MP phases...")
    mp_phases = load_mp_phases(args.mp_phases)
    print(f"  Loaded {len(mp_phases)} MP phases")
    
    print("\n" + "="*70)
    for system in args.systems:
        print(f"\nProcessing {system} system...")
        print("-"*70)
        
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
        
        plot_binary_hull_appendix(system, entries, entry_metadata, elements, 
                                  args.output_dir, args.e_above_hull_max)
    
    print("\n" + "="*70)
    print("All plots completed!")
    print(f"Output directory: {args.output_dir.resolve()}")


if __name__ == '__main__':
    main()
