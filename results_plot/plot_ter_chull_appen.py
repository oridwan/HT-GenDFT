#!/usr/bin/env python3
"""
Plot ternary convex hull diagrams for appendix (no interactive labels).

This script plots formation energy vs composition for ternary systems with:
- MP phases shown in legend with different colors
- New structures shown as rectangles in legend
- Legend positioned outside plot for clarity

Usage:
    python3 plot_ter_chull_appen.py --db Ter-Ele-HT/stable_electrides.db \
        --dft-results Ter-Ele-HT/dft_stability_results.json \
        --mp-phases Ter-Ele-HT/mp_vaspdft.json \
        --systems K-B-O Cs-Al-S \
        --output-dir hull_plots_appen \
        --e-above-hull-max 0.05
"""

import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
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


def reorder_elements_for_boundary(elements: List[str]) -> Tuple[List[str], List[int]]:
    """Reorder elements to match search_ternary_electrides.py convention."""
    elem_groups = []
    for elem in elements:
        if elem in GROUP_I_II_III:
            elem_groups.append((0, elem))
        elif elem in GROUP_III_IV:
            elem_groups.append((1, elem))
        elif elem in GROUP_V_VI_VII:
            elem_groups.append((2, elem))
        else:
            elem_groups.append((3, elem))
    
    elem_groups_sorted = sorted(elem_groups, key=lambda x: x[0])
    reordered = [elem for _, elem in elem_groups_sorted]
    
    index_mapping = [elements.index(elem) for elem in reordered]
    
    return reordered, index_mapping


def calculate_nexcess_boundaries_ternary(elem_A: str, elem_B: str, elem_C: str, 
                                         n_excess_max: float = 2.0, max_atoms: int = 15) -> Tuple[List[Tuple], List[Tuple]]:
    """Calculate exact composition boundaries for N_excess region in ternary A-B-C system."""
    from math import gcd
    from functools import reduce
    
    def gcd_multiple(numbers):
        return reduce(gcd, numbers)
    
    if elem_A not in VALENCE_ELECTRONS or elem_B not in VALENCE_ELECTRONS or elem_C not in VALENCE_ELECTRONS:
        return None, None
    
    val_A = VALENCE_ELECTRONS[elem_A]
    val_B = VALENCE_ELECTRONS[elem_B]
    val_C = VALENCE_ELECTRONS[elem_C]
    
    valid_points = []
    
    for l in range(1, max_atoms):
        for m in range(1, max_atoms):
            for n in range(1, max_atoms):
                if l + m + n > max_atoms:
                    continue
                
                g = gcd_multiple([l, m, n])
                l_p, m_p, n_p = l // g, m // g, n // g
                
                if (l_p, m_p, n_p) != (l, m, n):
                    continue
                
                n_excess = val_A * l_p + val_B * m_p - val_C * n_p
                
                if 0 < n_excess <= n_excess_max:
                    x_A = l_p / (l_p + m_p + n_p)
                    x_B = m_p / (l_p + m_p + n_p)
                    valid_points.append((x_A, x_B))
    
    if not valid_points:
        return None, None
    
    from scipy.spatial import ConvexHull
    try:
        hull = ConvexHull(valid_points)
        boundary_points = [valid_points[i] for i in hull.vertices]
        center = np.mean(boundary_points, axis=0)
        angles = [np.arctan2(p[1] - center[1], p[0] - center[0]) for p in boundary_points]
        sorted_indices = np.argsort(angles)
        boundary_points = [boundary_points[i] for i in sorted_indices]
        return boundary_points, None
    except:
        return valid_points, None


def load_electride_candidates_db(db_path: Path, systems: List[str], e_hull_max: float) -> List[Dict]:
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
    
    return rows


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
    electride_candidates: List[Dict],
    dft_results: Dict,
    mp_phases: List[Dict]
) -> Tuple[List[ComputedEntry], Dict[str, Dict], List[str]]:
    """Create ComputedEntry objects for a given chemical system."""
    
    entries = []
    entry_metadata = {}
    temp_electride_data = []
    
    system_normalized = '-'.join(sorted(system.split('-')))
    elements = sorted(system_normalized.split('-'))
    
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
    
    # Deduplicate electrides
    deduplicated_electrides = []
    comp_groups = {}
    
    for data in temp_electride_data:
        comp_key = data['composition'].reduced_formula
        if comp_key not in comp_groups:
            comp_groups[comp_key] = data
        else:
            if data['energy_per_atom'] < comp_groups[comp_key]['energy_per_atom']:
                comp_groups[comp_key] = data
    
    deduplicated_electrides = list(comp_groups.values())
    
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
        
        phase_elements = sorted(list(comp_dict.keys()))
        
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
            
            entry_metadata[entry_id] = {
                'is_electride': False,
                'mp_id': phase.get('mp_id', entry_id.split('-')[0] if '-' in entry_id else entry_id)
            }
    
    return entries, entry_metadata, elements


def bary2cart(coords_bary):
    """Convert barycentric to Cartesian coordinates for ternary plot."""
    a, b, c = coords_bary
    
    sqrt3 = np.sqrt(3)
    x = 0.5 * a + c
    y = sqrt3 / 2 * a
    
    return np.array([x, y])


def plot_ternary_hull_appendix(
    system: str,
    entries: List[ComputedEntry],
    entry_metadata: Dict[str, Dict],
    elements: List[str],
    output_dir: Path,
    e_hull_max: float = 0.05
):
    """Plot ternary convex hull for appendix (no interactive labels)."""
    
    # Separate MP entries from electride entries
    mp_entries = [e for e in entries if not entry_metadata.get(e.entry_id, {}).get('is_electride', False)]
    
    # Create phase diagrams
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
        
        comp_elems = sorted([str(e) for e in comp.elements])
        if len(comp_elems) > 3:
            continue
        
        bary_coords = np.array([comp.get_atomic_fraction(elem) for elem in elements])
        cart_coords = bary2cart(bary_coords)
        
        energy_per_atom = entry.energy / comp.num_atoms
        formation_energy = energy_per_atom - sum(bary_coords[i] * ref_energies[elements[i]] 
                                                  for i in range(3))
        
        e_above_hull = pd_all.get_e_above_hull(entry)
        
        metadata = entry_metadata.get(entry.entry_id, {})
        is_electride = metadata.get('is_electride', False)
        
        on_hull = entry in pd_all.stable_entries
        on_mp_hull = entry in pd_mp_only.stable_entries
        
        all_data.append({
            'cart_x': cart_coords[0],
            'cart_y': cart_coords[1],
            'bary': bary_coords,
            'formation_energy': formation_energy,
            'e_above_hull': e_above_hull,
            'on_hull': on_hull,
            'on_mp_hull': on_mp_hull,
            'is_electride': is_electride,
            'formula': comp.reduced_formula,
            'entry_id': entry.entry_id,
            'entry': entry,
            'metadata': metadata
        })
    
    # Identify phases
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
    
    mp_stable = [d for d in hull_phases if not d['is_electride']]
    
    # Create figure with fixed size
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Prepare colormap for energy above hull (rainbow: purple=low, red=high)
    all_e_above_hull = [d['e_above_hull'] for d in all_data 
                        if d['on_hull'] or d['e_above_hull'] <= e_hull_max]
    vmin, vmax = 0.0, max(all_e_above_hull) if all_e_above_hull else 0.1
    cmap = plt.cm.rainbow
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    stable_purple = cmap(norm(0.0))
    
    # Draw ternary axes
    draw_ternary_axes(ax, elements, 'black')
    
    # Highlight N_excess region
    reordered_elems, index_mapping = reorder_elements_for_boundary(elements)
    boundary_points, _ = calculate_nexcess_boundaries_ternary(
        reordered_elems[0], reordered_elems[1], reordered_elems[2], 
        n_excess_max=2.0, 
        max_atoms=20
    )
    if boundary_points is not None and len(boundary_points) >= 3:
        polygon_points = []
        for x_reordered_0, x_reordered_1 in boundary_points:
            x_reordered_2 = 1 - x_reordered_0 - x_reordered_1
            x_reordered = [x_reordered_0, x_reordered_1, x_reordered_2]
            
            x_plot = [0, 0, 0]
            for i, idx in enumerate(index_mapping):
                x_plot[idx] = x_reordered[i]
            
            bary = np.array(x_plot)
            cart = bary2cart(bary)
            polygon_points.append(cart)
        
        polygon = Polygon(polygon_points, alpha=0.35, facecolor='gold', 
                        edgecolor='none', zorder=1)
        ax.add_patch(polygon)
    
    # Draw hull facets
    draw_hull_facets(ax, pd_all, all_data, elements, linestyle='-', color=stable_purple, 
                    linewidth=2.5, alpha=0.8, label='Convex hull', zorder=2)
    
    draw_hull_facets(ax, pd_mp_only, all_data, elements, linestyle='--', color='black', 
                    linewidth=2.5, alpha=0.8, label='MP-only hull', zorder=3)
    
    # Define colors for MP stable phases
    mp_colors = plt.cm.tab20(np.linspace(0, 1, max(len(mp_stable), 1)))
    
    # Plot MP stable phases with different colors
    for i, d in enumerate(mp_stable):
        color = mp_colors[i % len(mp_colors)]
        formula_latex = formula_to_latex(d['formula'])
        ax.scatter(d['cart_x'], d['cart_y'], color=color, s=180, marker='o', 
                  edgecolors='black', linewidth=1.5, zorder=4, label=formula_latex)
    
    # Plot stable new structures (different marker for each, purple color matching colorbar at E=0)
    stable_markers = ['^', 'D', 'v', 'p', 'P', 'X', 'h', 'H', '<', '>']
    for i, d in enumerate(new_stable):
        marker = stable_markers[i % len(stable_markers)]
        ax.scatter(d['cart_x'], d['cart_y'], color=stable_purple, s=180, marker=marker, 
                  edgecolors='red', linewidth=2.5, zorder=6)
    
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
            ax.scatter([], [], color=stable_purple, s=180, marker=marker, 
                      edgecolors='red', linewidth=2.5, label=label)
    
    # Plot metastable new structures with colormap fill (based on energy above hull)
    for d in new_meta:
        facecolor = cmap(norm(d['e_above_hull']))
        ax.scatter(d['cart_x'], d['cart_y'], c=[facecolor], s=150, marker='s', 
                  edgecolors='red', linewidth=2, zorder=5)
    
    # Metastable structures are shown with colorbar, not in legend
    
    # Configure plot
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Set title
    if new_meta:
        title = f'{system} ($E_{{\\mathrm{{hull}}}} \\leq {e_hull_max:.2f}$ eV/atom)'
    else:
        title = f'{system}'
    ax.set_title(title, fontsize=24, fontweight='bold', pad=10)
    
    # Place legend at top left corner
    ax.legend(loc='upper left', fontsize=18, framealpha=0.65, edgecolor='black', facecolor='white')
    
    # Add colorbar for metastable structures if they exist
    if new_meta:
        from matplotlib.cm import ScalarMappable
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation='vertical', pad=0.02, shrink=0.85)
        cbar.set_label('Energy above hull (eV/atom)', fontsize=18, fontweight='bold')
        cbar.ax.tick_params(labelsize=16)
    
    # Use fixed margins
    fig.subplots_adjust(left=0.03, right=1.00, bottom=0.01, top=0.95)
    
    # Save figure
    output_path = output_dir / f'{system.replace("-", "")}_ternary_hull.png'
    fig.savefig(output_path, dpi=300)
    print(f"Saved: {output_path}")
    
    output_pdf = output_dir / f'{system.replace("-", "")}_ternary_hull.pdf'
    fig.savefig(output_pdf)
    print(f"Saved: {output_pdf}")
    
    plt.close(fig)


def draw_hull_facets(ax, phase_diagram, all_data, elements, label=None, zorder=1, **kwargs):
    """Draw the convex hull facets for a ternary phase diagram."""
    entry_to_data = {d['entry']: d for d in all_data}
    
    stable_entries = phase_diagram.stable_entries
    
    stable_coords = []
    stable_entries_list = list(stable_entries)
    
    for entry in stable_entries_list:
        if entry in entry_to_data:
            d = entry_to_data[entry]
            stable_coords.append([d['cart_x'], d['cart_y']])
    
    if len(stable_coords) < 3:
        return
    
    stable_coords = np.array(stable_coords)
    
    from scipy.spatial import Delaunay
    if len(stable_coords) >= 3:
        try:
            tri = Delaunay(stable_coords)
            
            drawn_edges = set()
            
            first_line = True
            for simplex in tri.simplices:
                for i in range(3):
                    p1_idx = simplex[i]
                    p2_idx = simplex[(i+1)%3]
                    edge = tuple(sorted([p1_idx, p2_idx]))
                    
                    if edge not in drawn_edges:
                        drawn_edges.add(edge)
                        p1 = stable_coords[p1_idx]
                        p2 = stable_coords[p2_idx]
                        
                        if first_line and label:
                            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], label=label, zorder=zorder, **kwargs)
                            first_line = False
                        else:
                            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], zorder=zorder, **kwargs)
        except:
            pass


def draw_ternary_axes(ax, elements, color='black'):
    """Draw the triangular axes for a ternary plot."""
    sqrt3 = np.sqrt(3)
    vertices = np.array([
        [0.5, sqrt3/2],
        [0, 0],
        [1, 0]
    ])
    
    triangle = Polygon(vertices, fill=False, edgecolor=color, linewidth=2.5)
    ax.add_patch(triangle)
    
    offset = 0.06
    ax.text(vertices[0,0], vertices[0,1] + offset, elements[0], 
           ha='center', va='bottom', fontsize=22, fontweight='bold')
    ax.text(vertices[1,0] - offset, vertices[1,1], elements[1], 
           ha='right', va='center', fontsize=22, fontweight='bold')
    ax.text(vertices[2,0] + offset, vertices[2,1], elements[2], 
           ha='left', va='center', fontsize=22, fontweight='bold')
    
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, sqrt3/2 + 0.1)


def main():
    parser = argparse.ArgumentParser(
        description='Plot ternary convex hull diagrams for appendix',
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
                       help='Chemical systems to plot (e.g., K-B-O Cs-Al-S)')
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
    structure_ids = [row['structure_id'] for row in electride_candidates]
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
        
        if len(system.split('-')) != 3:
            print(f"  Error: {system} is not a ternary system (must have 3 elements)")
            continue
        
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
        
        plot_ternary_hull_appendix(system, entries, entry_metadata, elements, 
                                   args.output_dir, args.e_above_hull_max)
    
    print("\n" + "="*70)
    print("All plots completed!")
    print(f"Output directory: {args.output_dir.resolve()}")


if __name__ == '__main__':
    main()
