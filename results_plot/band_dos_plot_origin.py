#!/usr/bin/env python3
"""
Plot electronic band structure and density of states (DOS) for electride candidates.

This script reads VASP band structure and DOS calculations and creates combined plots
showing:
- Left: Electronic band structure with high-symmetry k-points and Fermi level
- Right: Partial DOS projected to each element

Usage:
    python3 band_dos_plot.py --db Bin-Ele-HT/stable_electrides.db --structure-id Ca3P1_s030
    python3 band_dos_plot.py --csv Bin-Ele-HT/stable_electrides.csv --structure-id Ca3P1_s030 Y9N8_s013
    python3 results_plot/band_dos_plot.py --base-dir ELECTRONIC_JOBS-Boron --all
"""

import argparse
import warnings
import re
from pathlib import Path
from typing import Dict, List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
from vasprun import vasprun


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


def parse_composition(formula: str) -> Dict[str, float]:
    """
    Parse chemical formula into element-count dictionary.
    
    Examples:
        Ca3P1 -> {'Ca': 3, 'P': 1}
        Y9N8 -> {'Y': 9, 'N': 8}
        K6BO4 -> {'K': 6, 'B': 1, 'O': 4}
    """
    pattern = r'([A-Z][a-z]?)(\d*)'
    composition = {}
    
    for match in re.finditer(pattern, formula):
        element = match.group(1)
        count_str = match.group(2)
        count = int(count_str) if count_str else 1
        
        if element:
            composition[element] = composition.get(element, 0) + count
    
    return composition


def calculate_excess_electrons(composition: str) -> int:
    """
    Calculate excess valence electrons in a composition.
    
    - Sum valence electrons from electropositive elements (Groups I, II, III)
    - Subtract valence electrons needed by electronegative elements (Groups V, VI, VII)
    
    Returns:
        Excess electrons (should be 0 < N_excess <= 4 for electride candidates)
    """
    comp = parse_composition(composition)
    
    excess = 0.0
    for elem_symbol, amount in comp.items():
        if elem_symbol in VALENCE_ELECTRONS:
            val = VALENCE_ELECTRONS[elem_symbol]
            
            # Electronegative elements subtract electrons
            if elem_symbol in ELECTRONEGATIVE_ELEMENTS:
                excess -= abs(val) * amount
            else:
                # Electropositive elements contribute electrons
                excess += val * amount
    
    return int(round(excess))


def formula_to_latex(formula: str) -> str:
    """Convert chemical formula to LaTeX format with subscripts.
    
    Ignores subscript 1 (e.g., Al1 -> Al, not Al$_{1}$).
    
    Examples:
        Ca5P3 -> Ca$_{5}$P$_{3}$
        Ca7Al1P5 -> Ca$_{7}$AlP$_{5}$
        K6BO4 -> K$_{6}$BO$_{4}$
    """
    pattern = r'([A-Z][a-z]?)(\d+)'
    
    def replace_func(m):
        element = m.group(1)
        count = m.group(2)
        if count == '1':
            return element
        else:
            return f'{element}$_{{{count}}}$'
    
    return re.sub(pattern, replace_func, formula)


def get_pearson_symbols_from_db(db_path: Path, structure_ids: List[str]) -> Dict[str, str]:
    """Get Pearson symbols for structure IDs from ASE database."""
    try:
        from ase.db import connect
    except ImportError:
        return {}
    
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


def get_pearson_symbols_from_csv(csv_path: Path, structure_ids: List[str]) -> Dict[str, str]:
    """Get Pearson symbols from CSV file."""
    import pandas as pd
    
    if not csv_path.exists():
        return {}
    
    try:
        df = pd.read_csv(csv_path)
        pearson_map = {}
        for sid in structure_ids:
            matches = df[df['structure_id'] == sid]
            if not matches.empty:
                pearson_map[sid] = matches.iloc[0].get('pearson_symbol', '')
        return pearson_map
    except Exception as e:
        return {}


def format_kpoint_label(label: str) -> str:
    """Format k-point label with proper LaTeX formatting.
    
    Examples:
        Gamma -> $\\Gamma$
        X_1 -> X$_{1}$
        GAMMA -> $\\Gamma$
    """
    # Handle backslash prefix
    if label.startswith('\\'):
        label = label[1:]
    
    # Convert common Greek letters to LaTeX
    greek_map = {
        'GAMMA': r'\Gamma',
        'Gamma': r'\Gamma',
        'DELTA': r'\Delta',
        'Delta': r'\Delta',
        'LAMBDA': r'\Lambda',
        'Lambda': r'\Lambda',
        'SIGMA': r'\Sigma',
        'Sigma': r'\Sigma',
    }
    
    # Check if it's a Greek letter
    for name, latex in greek_map.items():
        if label.upper() == name.upper():
            return f'${latex}$'
    
    # Handle subscripts (e.g., X_1, P_1)
    if '_' in label:
        parts = label.split('_')
        base = parts[0]
        subscript = '_'.join(parts[1:])
        # Check if base is Greek
        if base.upper() in [k.upper() for k in greek_map.keys()]:
            base_latex = greek_map.get(base, greek_map.get(base.upper(), greek_map.get(base.capitalize(), base)))
            return f'${base_latex}_{{{subscript}}}$'
        else:
            return f'{base}$_{{{subscript}}}$'
    
    return label


def parse_kpoints_file(kpoints_path: Path):
    """Parse KPOINTS file to extract high-symmetry point labels and positions.
    
    Handles overlapping k-points (end of one segment = start of next segment)
    and disconnected segments (labeled with " | " separator).
    
    Returns:
        labels: List of formatted k-point labels (with LaTeX)
        positions: List of k-point indices where labels should be placed
    """
    with open(kpoints_path, 'r') as f:
        lines = [l.strip() for l in f.readlines()]
    
    if len(lines) < 4:
        return [], []
    
    # Get number of divisions
    try:
        divisions = int(lines[1])
    except:
        divisions = 40  # default
    
    # Extract labels from lines with '!'
    all_labels = []
    for line in lines[4:]:
        if '!' in line:
            label = line.split('!')[-1].strip()
            all_labels.append(label)
    
    # Build k-path with labels and positions, detecting disconnected segments
    kpoint_map = {}  # position -> label
    
    current_pos = 0
    for i in range(0, len(all_labels), 2):
        if i + 1 >= len(all_labels):
            break
        
        label_start = all_labels[i]
        label_end = all_labels[i + 1]
        
        # Check if this segment is disconnected from previous one
        if current_pos in kpoint_map:
            # Position already has a label - check if it's the same or disconnected
            existing_label = kpoint_map[current_pos]
            if existing_label != label_start:
                # Disconnected segment! Merge labels with " | "
                kpoint_map[current_pos] = f"{existing_label} | {label_start}"
        else:
            # First time at this position
            kpoint_map[current_pos] = label_start
        
        # Move to end of segment
        current_pos += divisions
        
        # Add end point
        if current_pos not in kpoint_map:
            kpoint_map[current_pos] = label_end
    
    # Convert to sorted lists and format labels
    sorted_positions = sorted(kpoint_map.keys())
    labels = []
    for pos in sorted_positions:
        label_str = kpoint_map[pos]
        # Handle merged labels (with " | ")
        if ' | ' in label_str:
            parts = label_str.split(' | ')
            formatted_parts = [format_kpoint_label(part) for part in parts]
            labels.append(' | '.join(formatted_parts))
        else:
            labels.append(format_kpoint_label(label_str))
    
    positions = sorted_positions
    
    return labels, positions


def find_band_dos_dirs(structure_id: str, base_dir: Path = None) -> tuple:
    """Find BAND and DOS directories for a structure.
    
    Searches in:
    - Bin-Ele-HT/candidates/structure_id/
    - Ter-Ele-HT/candidates/structure_id/
    
    Returns:
        (band_dir, dos_dir) paths or (None, None) if not found
    """
    if base_dir is None:
        base_dir = Path.cwd()
    
    search_paths = [
        base_dir / 'Bin-Ele-HT' / 'candidates' / structure_id,
        base_dir / 'Ter-Ele-HT' / 'candidates' / structure_id,
        base_dir / structure_id,  # Direct path
    ]
    
    for search_path in search_paths:
        if search_path.exists():
            band_dir = search_path / 'BAND'
            dos_dir = search_path / 'DOS'
            
            if band_dir.exists() and dos_dir.exists():
                return band_dir, dos_dir

    # Recursive fallback for nested layouts such as:
    # ELECTRONIC_JOBS-Boron/{composition}/{structure_id}/BAND,DOS
    if base_dir.exists():
        for candidate in base_dir.glob(f"**/{structure_id}"):
            if not candidate.is_dir():
                continue
            band_dir = candidate / 'BAND'
            dos_dir = candidate / 'DOS'
            if band_dir.exists() and dos_dir.exists():
                return band_dir, dos_dir
    
    return None, None


def discover_structure_ids(base_dir: Path) -> List[str]:
    """Discover structure IDs that contain both BAND and DOS directories."""
    if not base_dir.exists():
        return []

    discovered = set()
    for band_dir in base_dir.glob("**/BAND"):
        if not band_dir.is_dir():
            continue
        struct_dir = band_dir.parent
        dos_dir = struct_dir / 'DOS'
        if dos_dir.exists() and dos_dir.is_dir():
            discovered.add(struct_dir.name)

    return sorted(discovered)


def plot_band_dos(
    structure_id: str,
    band_dir: Path,
    dos_dir: Path,
    output_dir: Path,
    pearson_symbol: str = '',
    energy_range: tuple = (-4, 2)
):
    """Plot band structure and DOS for a structure.
    
    Args:
        structure_id: Structure identifier
        band_dir: Path to BAND calculation directory
        dos_dir: Path to DOS calculation directory
        output_dir: Output directory for plots
        energy_range: Energy range relative to Fermi level (min, max) in eV
    """
    
    print(f"\nProcessing {structure_id}...")
    print(f"  BAND dir: {band_dir}")
    print(f"  DOS dir: {dos_dir}")
    
    # Read band structure
    band_xml = band_dir / 'vasprun.xml'
    if not band_xml.exists():
        print(f"  Error: vasprun.xml not found in {band_dir}")
        return
    
    print(f"  Reading band structure...")
    try:
        band_run = vasprun(str(band_xml))
    except Exception as e:
        print(f"  Error reading band structure: {e}")
        return
    
    # Read DOS
    dos_xml = dos_dir / 'vasprun.xml'
    if not dos_xml.exists():
        print(f"  Error: vasprun.xml not found in {dos_dir}")
        return
    
    print(f"  Reading DOS...")
    try:
        dos_run = vasprun(str(dos_xml))
    except Exception as e:
        print(f"  Error reading DOS: {e}")
        return
    
    # Get composition from structure_id and format with LaTeX subscripts
    composition = structure_id.split('_')[0] if '_' in structure_id else structure_id
    composition_latex = formula_to_latex(composition)
    
    # Calculate N_excess
    n_excess = calculate_excess_electrons(composition)
    
    # Create title
    if pearson_symbol:
        title = f'{pearson_symbol}-{composition_latex} ($N_{{\\mathrm{{excess}}}}={n_excess}$)'
    else:
        title = f'{composition_latex} ($N_{{\\mathrm{{excess}}}}={n_excess}$)'
    
    # Get Fermi energy
    efermi = band_run.values['calculation']['efermi']
    print(f"  Fermi energy: {efermi:.4f} eV")
    
    # Get band structure eigenvalues: shape (n_kpts, n_bands, 2)
    # Last dimension: [energy, occupation]
    eigenvalues = np.array(band_run.values['calculation']['eband_eigenvalues'])
    n_kpts, n_bands, _ = eigenvalues.shape
    bands = eigenvalues[:, :, 0]  # Extract energies
    
    print(f"  Band structure: {n_kpts} k-points, {n_bands} bands")
    
    # Shift bands relative to Fermi level
    bands_shifted = bands - efermi
    
    # Parse KPOINTS file for labels
    kpoints_file = band_dir / 'KPOINTS'
    kpt_labels, kpt_positions = parse_kpoints_file(kpoints_file)
    
    print(f"  K-path labels: {kpt_labels}")
    
    # Calculate k-point distances
    k_distance = np.arange(n_kpts)
    
    # Get DOS data
    # TDOS shape: (n_spin, n_energy, 3) where last dim is [energy, dos_up/total, dos_down/integrated]
    tdos = np.array(dos_run.values['calculation']['tdos'])
    if len(tdos.shape) == 3:
        # Take first spin channel
        tdos = tdos[0]
    dos_energies = tdos[:, 0] - efermi  # Shift relative to Fermi
    total_dos = tdos[:, 1]  # Total DOS
    
    # Get projected DOS by element
    pdos_data = dos_run.values['calculation']['pdos']
    elements = dos_run.values['elements']
    
    print(f"  Elements: {elements}")
    print(f"  DOS: {len(dos_energies)} energy points")
    
    # Process PDOS
    # pdos_data is a list of arrays, one for each atom
    # Each array has shape (n_spin, n_energy, n_orbitals+1) where first column is energy
    # We need to sum by element type
    element_pdos = {}
    if pdos_data is not None and len(pdos_data) > 0:
        composition_dict = dos_run.values['composition']
        atom_types = []
        for elem in elements:
            count = composition_dict[elem]
            atom_types.extend([elem] * count)
        
        # Sum PDOS by element
        for i, pdos_array in enumerate(pdos_data):
            if i >= len(atom_types):
                break
            elem = atom_types[i]
            pdos_array = np.array(pdos_array)
            
            # Handle shape (n_spin, n_energy, n_columns)
            if len(pdos_array.shape) == 3:
                pdos_array = pdos_array[0]  # Take first spin channel
            
            # Now shape is (n_energy, n_columns)
            # Sum over all orbital columns (skip energy column at index 0)
            pdos_total = pdos_array[:, 1:].sum(axis=1)
            
            if elem not in element_pdos:
                element_pdos[elem] = np.zeros_like(pdos_total)
            element_pdos[elem] += pdos_total
    
    # Create figure with two subplots
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[2, 1], wspace=0.08)
    ax_band = fig.add_subplot(gs[0])
    ax_dos = fig.add_subplot(gs[1], sharey=ax_band)
    
    # Plot band structure
    for i in range(n_bands):
        ax_band.plot(k_distance, bands_shifted[:, i], color='blue', 
                    linewidth=1.5, alpha=0.7, zorder=2)
    
    # Plot Fermi level
    ax_band.axhline(y=0, color='red', linestyle='--', linewidth=2.5, 
                   label='Fermi level', zorder=5)
    
    # Set high-symmetry k-point labels
    if kpt_labels and kpt_positions:
        ax_band.set_xticks(kpt_positions)
        ax_band.set_xticklabels(kpt_labels, fontsize=20)
        
        # Add vertical lines at high-symmetry points (black solid lines)
        for pos in kpt_positions:
            ax_band.axvline(x=pos, color='black', linestyle='-', 
                          linewidth=2, alpha=0.8, zorder=1)
    
    # Configure band structure axes
    ax_band.set_ylabel('Energy (eV)', fontsize=24, fontweight='bold', labelpad=10)
    ax_band.set_ylim(energy_range[0], energy_range[1])
    # Add small padding to ensure first and last k-point labels are visible
    ax_band.set_xlim(-0.001 * n_kpts, n_kpts * 1.001)
    
    # Set major tick interval to 1.0 eV
    ax_band.yaxis.set_major_locator(MultipleLocator(1.0))
    # Format y-axis tick labels to always show 1 decimal place
    ax_band.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    # Add horizontal gray dashed lines at major y-ticks (excluding y=0 where Fermi level is)
    ax_band.grid(False)  # Disable automatic grid
    for ytick in ax_band.get_yticks():
        if abs(ytick) > 0.01:  # Skip y=0 (Fermi level already shown as red dashed line)
            ax_band.axhline(y=ytick, color='gray', linestyle='--', linewidth=1.2, alpha=0.4, zorder=0)
    
    # Set background and edge colors
    ax_band.set_facecolor('white')
    for spine in ax_band.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    
    ax_band.legend(loc='upper right', fontsize=20, framealpha=0.8, edgecolor='none', facecolor='white')
    ax_band.tick_params(axis='both', which='major', labelsize=20)
    ax_band.tick_params(axis='x', which='both', bottom=False, top=False)  # Hide x-ticks
    ax_band.tick_params(axis='y', which='major', labelleft=True)  # Ensure y-tick labels are visible
    
    # Plot total DOS
    ax_dos.plot(total_dos, dos_energies, color='black', linewidth=1.5, 
               label='Total', alpha=0.8, zorder=3)
    ax_dos.fill_betweenx(dos_energies, 0, total_dos, color='gray', 
                         alpha=0.2, zorder=1)
    
    # Plot projected DOS by element
    if element_pdos:
        colors = plt.cm.Dark2(np.arange(0, 0.2 * len(element_pdos), 0.2))
        for idx, (element, pdos) in enumerate(element_pdos.items()):
            ax_dos.plot(pdos, dos_energies, linewidth=2.5, 
                       label=element, alpha=0.8, color=colors[idx], zorder=2)
    
    # Plot Fermi level on DOS
    ax_dos.axhline(y=0, color='red', linestyle='--', linewidth=2.5, zorder=5)
    
    # Calculate dynamic x-axis upper limit for DOS
    # Use median of total DOS within the energy range
    energy_mask = (dos_energies >= energy_range[0]) & (dos_energies <= energy_range[1])
    mean_dos_in_range = np.mean(total_dos[energy_mask])
    max_dos_in_range = total_dos[energy_mask].max()
    dos_xlim_max = mean_dos_in_range * 2.0
    
    # Configure DOS axes
    ax_dos.set_xlabel('DOS (states/eV)', fontsize=24, fontweight='bold')
    ax_dos.set_xlim(0, dos_xlim_max)
    
    # Set consistent number of x-axis ticks to ensure uniform plot width across all structures
    ax_dos.xaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
    
    # Sync y-ticks with band structure plot but hide labels (do this before setting ylim)
    ax_dos.set_yticks(ax_band.get_yticks())
    ax_dos.tick_params(axis='y', which='major', labelleft=False, left=True, right=False, length=5)
    
    # Set ylim AFTER syncing ticks to prevent matplotlib from auto-expanding
    ax_dos.set_ylim(energy_range[0], energy_range[1])
    ax_band.set_ylim(energy_range[0], energy_range[1])  # Re-enforce band subplot ylim too
    
    # Disable autoscaling to lock in our ylim settings
    ax_dos.autoscale(enable=False, axis='y')
    ax_band.autoscale(enable=False, axis='y')
    
    # Add horizontal gray dashed lines at major y-ticks (excluding y=0 where Fermi level is)
    for ytick in ax_dos.get_yticks():
        if abs(ytick) > 0.01:  # Skip y=0 (Fermi level already shown as red dashed line)
            ax_dos.axhline(y=ytick, color='gray', linestyle='--', linewidth=1.2, alpha=0.4, zorder=0)
    
    # Add vertical gray dashed lines at major x-ticks
    ax_dos.grid(True, alpha=0.4, linestyle='--', linewidth=1.2, axis='x', color='gray')
    ax_dos.grid(False, axis='y')  # Disable automatic y-grid (we drew manual lines above)
    
    # Set background and edge colors
    ax_dos.set_facecolor('white')
    for spine in ax_dos.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    
    ax_dos.legend(loc='upper right', fontsize=20, framealpha=0.8, edgecolor='none', facecolor='white')
    # Format x-axis tick labels to always show 1 decimal place
    ax_dos.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax_dos.tick_params(axis='x', which='major', labelsize=20)
    
    # Set main title
    fig.suptitle(title, fontsize=28, fontweight='bold', y=0.98)
    
    # Use fixed subplot margins to ensure consistent plot dimensions across all structures
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.92)
    
    # Save figure
    output_png = output_dir / f'{structure_id}_band_dos.png'
    output_pdf = output_dir / f'{structure_id}_band_dos.pdf'
    
    # Suppress tight_layout warning
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='This figure includes Axes that are not compatible with tight_layout')
        fig.savefig(output_png, dpi=300)
        print(f"  Saved: {output_png}")
        
        fig.savefig(output_pdf)
        print(f"  Saved: {output_pdf}")
    
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='Plot electronic band structure and DOS for electride candidates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Input sources
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument('--csv', type=Path,
                           help='Path to stable_electrides.csv')
    input_group.add_argument('--db', type=Path,
                           help='Path to stable_electrides.db')
    
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument('--structure-id', nargs='+',
                              help='Structure ID(s) to plot (e.g., Ca3P1_s030)')
    target_group.add_argument('--all', action='store_true',
                              help='Plot all structures found under --base-dir')
    
    parser.add_argument('--output-dir', type=Path, default=Path('band_dos_plots'),
                       help='Output directory for plots (default: band_dos_plots)')
    
    parser.add_argument('--energy-range', type=float, nargs=2, default=(-4, 2),
                       help='Energy range relative to Fermi level (min max) in eV (default: -4 2)')
    
    parser.add_argument('--base-dir', type=Path, default=None,
                       help='Base directory to search for structure data (default: current directory)')
    
    args = parser.parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("Electronic Band Structure and DOS Plotting")
    print("="*70)
    
    if args.all:
        search_base = args.base_dir if args.base_dir is not None else Path.cwd()
        structure_ids = discover_structure_ids(search_base)
        print(f"\nDiscovered {len(structure_ids)} structures under: {search_base}")
        if not structure_ids:
            print("No structures with both BAND and DOS directories were found.")
            return
    else:
        structure_ids = args.structure_id

    # Load Pearson symbols from database or CSV
    pearson_map = {}
    if args.db:
        print(f"\nLoading Pearson symbols from database: {args.db}")
        pearson_map = get_pearson_symbols_from_db(args.db, structure_ids)
        print(f"  Loaded Pearson symbols for {len([v for v in pearson_map.values() if v])} structures")
    elif args.csv:
        print(f"\nLoading Pearson symbols from CSV: {args.csv}")
        pearson_map = get_pearson_symbols_from_csv(args.csv, structure_ids)
        print(f"  Loaded Pearson symbols for {len([v for v in pearson_map.values() if v])} structures")
    else:
        print("\nNo database or CSV provided, Pearson symbols will not be shown")
    
    # Process each structure
    success_count = 0
    for structure_id in structure_ids:
        # Find BAND and DOS directories
        band_dir, dos_dir = find_band_dos_dirs(structure_id, args.base_dir)
        
        if band_dir is None or dos_dir is None:
            print(f"\nError: Could not find BAND and DOS directories for {structure_id}")
            print(f"  Searched in:")
            if args.base_dir:
                print(f"    - {args.base_dir}/Bin-Ele-HT/candidates/{structure_id}/")
                print(f"    - {args.base_dir}/Ter-Ele-HT/candidates/{structure_id}/")
                print(f"    - {args.base_dir}/{structure_id}/")
                print(f"    - {args.base_dir}/**/{structure_id}/")
            else:
                print(f"    - ./Bin-Ele-HT/candidates/{structure_id}/")
                print(f"    - ./Ter-Ele-HT/candidates/{structure_id}/")
                print(f"    - ./{structure_id}/")
                print(f"    - ./**/{structure_id}/")
            continue
        
        try:
            pearson = pearson_map.get(structure_id, '')
            plot_band_dos(structure_id, band_dir, dos_dir, args.output_dir, 
                         pearson, tuple(args.energy_range))
            success_count += 1
        except Exception as e:
            print(f"\nError processing {structure_id}: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print(f"Completed! Successfully plotted {success_count}/{len(structure_ids)} structures")
    print(f"Output directory: {args.output_dir.resolve()}")


if __name__ == '__main__':
    main()
