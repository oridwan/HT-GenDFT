#!/usr/bin/env python3
"""
Phonon Band Structure and DOS Plotting Script

Reads pre-computed phonon data files and recreates publication-quality plots 
locally without needing VASP or phonopy.

Requirements:
    - phonon_band.dat: Band structure data (REQUIRED)
    - phonon_dos.dat: DOS data (REQUIRED)
    - band_kpath.dat: K-path metadata with lattice, segments, and q-points (REQUIRED)
    - stable_electrides.db: ASE database (optional, for Pearson symbols with --db)
    - POSCAR: Structure file (optional, for Pearson symbols with --poscar)

Usage:
    python3 phonon_band_dos_plot.py --results-dir Bin-Ele-HT/candidates/ --db Bin-Ele-HT/stable_electrides.db --structure-ids Ba2P1_s001 Ca5P3_s021
    python3 phonon_band_dos_plot.py --results-dir Ter-Ele-HT/candidates/ --db Ter-Ele-HT/stable_electrides.db --structure-ids K6B1O4_s013 --output-dir ./phonon_plots
    python3 phonon_band_dos_plot.py --results-dir Bin-Ele-HT/candidates/ --structure-ids Ba2P1_s001 --poscar  # Use POSCAR for Pearson symbol
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, MaxNLocator
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import sys

# Matplotlib configuration for publication-quality plots
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5


def get_pearson_symbols_from_db(db_path: Path, structure_ids: List[str]) -> Dict[str, str]:
    """Get Pearson symbols for structure IDs from ASE database."""
    try:
        from ase.db import connect
    except ImportError:
        print("Warning: ASE not installed, cannot read Pearson symbols from database")
        return {}
    
    if not db_path.exists():
        print(f"Warning: Database file not found: {db_path}")
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


def get_pearson_symbol(structure) -> str:
    """
    Get full Pearson symbol from structure (e.g., oI14, mS26).
    
    Args:
        structure: Pymatgen Structure
    
    Returns:
        Pearson symbol with number of atoms (e.g., 'oI14', 'mS26', 'cF8')
    """
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    
    try:
        sga = SpacegroupAnalyzer(structure)
        
        # Get conventional standard structure
        conv_structure = sga.get_conventional_standard_structure()
        n_atoms = len(conv_structure)
        
        crystal_system = sga.get_crystal_system()
        
        # Get space group symbol for centering
        spg_symbol = sga.get_space_group_symbol()
        
        # Map crystal system to Pearson letter
        system_map = {
            'triclinic': 'a',
            'monoclinic': 'm',
            'orthorhombic': 'o',
            'tetragonal': 't',
            'trigonal': 'h',  # hexagonal setting
            'hexagonal': 'h',
            'cubic': 'c'
        }
        
        # Get centering from space group symbol
        # Pearson notation for centering (note: monoclinic C -> S)
        centering_map = {
            'P': 'P',  # Primitive
            'A': 'S',  # Base-centered (A-face) -> S in Pearson
            'B': 'S',  # Base-centered (B-face) -> S in Pearson  
            'C': 'S',  # Base-centered (C-face) -> S in Pearson (for monoclinic)
            'I': 'I',  # Body-centered
            'F': 'F',  # Face-centered
            'R': 'R',  # Rhombohedral
        }
        
        system_letter = system_map.get(crystal_system.lower(), 'a')
        
        # Extract centering from space group symbol
        centering = 'P'  # default
        spg_first_letter = spg_symbol[0]
        
        # Special handling for monoclinic: A, B, C all map to S
        if crystal_system.lower() == 'monoclinic':
            centering = centering_map.get(spg_first_letter, 'P')
        else:
            # For non-monoclinic, use standard mapping but keep A as A, etc.
            if spg_first_letter in ['P', 'I', 'F', 'R']:
                centering = spg_first_letter
            elif spg_first_letter in ['A', 'B', 'C']:
                # For orthorhombic, keep the specific face centering
                if crystal_system.lower() == 'orthorhombic':
                    centering = spg_first_letter
                else:
                    centering = 'S'
        
        # Full Pearson symbol: crystal system + centering + number of atoms
        pearson = f"{system_letter}{centering}{n_atoms}"
        return pearson
    
    except Exception as e:
        print(f"    Warning: Could not determine Pearson symbol: {e}")
        return None


def formula_to_latex(formula: str) -> str:
    """
    Convert chemical formula to LaTeX format with subscripts.
    
    Example: 'Cs6Al2S5' -> 'Cs$_{6}$Al$_{2}$S$_{5}$'
    """
    result = []
    i = 0
    while i < len(formula):
        char = formula[i]
        if char.isdigit():
            # Check if we're already in a subscript
            if result and result[-1].endswith('_{'):
                result.append(char)
            else:
                result.append('$_{' + char)
            # Look ahead for more digits
            i += 1
            while i < len(formula) and formula[i].isdigit():
                result.append(formula[i])
                i += 1
            result.append('}$')
        else:
            result.append(char)
            i += 1
    return ''.join(result)


def parse_band_data(band_file: Path) -> Dict:
    """
    Parse phonon_band.dat file.
    
    Returns:
        dict with keys:
            - 'distances': distances along path
            - 'qpoints': fractional coordinates of q-points
            - 'frequencies': array of shape (n_bands, n_qpoints)
            - 'n_bands': number of bands
    """
    data = np.loadtxt(band_file)
    
    # Column 0: Distance
    # Columns 1-3: qx, qy, qz
    # Columns 4+: Band frequencies
    
    distances = data[:, 0]
    qpoints = data[:, 1:4]
    frequencies = data[:, 4:]  # Shape: (n_qpoints, n_bands)
    frequencies = frequencies.T  # Transpose to (n_bands, n_qpoints)
    
    return {
        'distances': distances,
        'qpoints': qpoints,
        'frequencies': frequencies,
        'n_bands': frequencies.shape[0]
    }


def read_kpath_metadata(kpath_file: Path, qpoints: np.ndarray, distances: np.ndarray) -> Tuple[List[Tuple[int, int]], List[float], List[str]]:
    """
    Read k-path metadata from band_kpath.dat.
    
    This file contains lattice vectors, segment structure, high-symmetry tick
    labels with distances and fractional coordinates.
    
    Args:
        kpath_file: Path to band_kpath.dat file
        qpoints: Array of fractional q-point coordinates, shape (n_qpoints, 3) (unused)
        distances: Array of distances along path (from phonon_band.dat) (unused)
    
    Returns:
        segments: List of (start_idx, end_idx) tuples for each segment
        tick_positions: List of distances for high-symmetry points
        tick_labels: List of labels for high-symmetry points
    """
    segments = []
    tick_positions = []
    tick_labels = []
    
    with open(kpath_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split()
            if not parts:
                continue
            
            # Parse SEGMENT lines
            if parts[0] == 'SEGMENT' and len(parts) >= 5:
                start_idx = int(parts[2])
                end_idx = int(parts[3])
                segments.append((start_idx, end_idx + 1))
            
            # Parse TICK lines
            elif parts[0] == 'TICK' and len(parts) >= 6:
                tick_dist = float(parts[1])
                label = ' '.join(parts[5:])
                tick_positions.append(tick_dist)
                tick_labels.append(label)
    
    return segments, tick_positions, tick_labels


def parse_dos_data(dos_file: Path) -> Dict:
    """
    Parse phonon_dos.dat file.
    
    Returns:
        dict with keys:
            - 'frequencies': frequency values
            - 'total_dos': total DOS
            - 'projected_dos': dict of element -> DOS array
            - 'elements': list of element names
    """
    # Read header to get element names
    elements = []
    data_start_line = 0
    
    with open(dos_file, 'r') as f:
        for line_num, line in enumerate(f):
            if line.startswith('# Column'):
                parts = line.split(':')
                if len(parts) == 2 and 'projected DOS' in parts[1]:
                    # Extract element name
                    element = parts[1].split('projected DOS')[0].strip()
                    elements.append(element)
            elif line.startswith('#'):
                continue
            elif line.strip() and not line.strip()[0].replace('.', '').replace('-', '').isdigit():
                # This is the column header line (e.g., "Frequency  Total  Ba  P")
                # Parse element names from here if not found in comments
                if not elements:
                    cols = line.split()
                    # Skip first two columns (Frequency, Total)
                    elements = cols[2:]
                data_start_line = line_num + 1
                break
            elif line.strip():
                # First data line - no column header found
                data_start_line = line_num
                break
    
    # Read data, skipping header lines
    data = np.loadtxt(dos_file, skiprows=data_start_line)
    
    frequencies = data[:, 0]
    total_dos = data[:, 1]
    
    projected_dos = {}
    for i, element in enumerate(elements):
        projected_dos[element] = data[:, i + 2]
    
    return {
        'frequencies': frequencies,
        'total_dos': total_dos,
        'projected_dos': projected_dos,
        'elements': elements
    }


def plot_phonon_band_dos(
    band_data: Dict,
    dos_data: Dict,
    output_png: Path,
    output_pdf: Path,
    title: str = "Phonon Band Structure and DOS",
    kpath_file: Path = None
):
    """
    Create combined phonon band structure and DOS plot.
    
    Matches the style of postproc_phonon.py exactly.
    
    Args:
        band_data: Dictionary with band structure data
        dos_data: Dictionary with DOS data
        output_png: Output PNG file path
        output_pdf: Output PDF file path
        title: Plot title
        kpath_file: Path to band_kpath.dat (segment structure and tick labels)
    """
    # Enable LaTeX text rendering (matching postproc_phonon.py)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['font.family'] = 'DejaVu Sans'
    
    # Create figure with two subplots (band structure and DOS)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True,
                                     gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.08})
    
    # Extract data
    distances = band_data['distances']
    qpoints = band_data['qpoints']
    frequencies = band_data['frequencies']
    dos_freqs = dos_data['frequencies']
    total_dos = dos_data['total_dos']
    projected_dos = dos_data['projected_dos']
    elements = dos_data['elements']
    
    # Read k-path metadata (segments and tick labels)
    segments, positions, labels = read_kpath_metadata(kpath_file, qpoints, distances)
    
    # Plot band structure
    for band_idx in range(band_data['n_bands']):
        band_frequencies = frequencies[band_idx]
        ax1.plot(distances, band_frequencies, 'r-', linewidth=1.0, zorder=2)
    
    # Add high-symmetry point markers and labels
    for pos in positions:
        ax1.axvline(x=pos, color='black', linestyle='-', linewidth=0.5, zorder=1)
    
    ax1.set_xticks(positions)
    ax1.set_xticklabels(labels, fontsize=14)
    ax1.set_xlabel('Wave vector', fontsize=16, fontweight='bold')
    ax1.set_ylabel('Frequency (THz)', fontsize=16, fontweight='bold')
    ax1.axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=5)
    ax1.set_xlim(distances[0] - 0.0001, distances[-1] + 0.0001)
    
    # Set y-axis limits and ticks
    ax1.set_ylim(-2, 12)
    ax1.yaxis.set_major_locator(MultipleLocator(2.0))
    ax1.yaxis.set_minor_locator(MultipleLocator(1.0))
    ax1.tick_params(axis='y', which='major', labelsize=14)
    ax1.tick_params(axis='y', which='minor', length=3)
    
    # Plot DOS on right (ax2)
    ax2.plot(total_dos, dos_freqs, 'k-', linewidth=1.0, 
             label='Total', alpha=0.8, zorder=3)
    
    # Plot element-projected DOS with different colors
    if elements:
        colors = plt.cm.tab10(np.arange(0, 0.2 * len(elements), 0.2))
        for idx, element in enumerate(elements):
            ax2.plot(projected_dos[element], dos_freqs, 
                    linewidth=1.5, label=str(element), alpha=0.8, color=colors[idx], zorder=2)
    
    dos_mean = np.mean(total_dos)
    dos_xlim_max = dos_mean * 4.0
    
    ax2.set_xlabel('DOS', fontsize=16, fontweight='bold')
    ax2.axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=5)
    ax2.legend(loc='upper right', fontsize=14, framealpha=0.8, edgecolor='none', facecolor='white')
    ax2.set_xlim(0, dos_xlim_max)
    
    # Set consistent number of x-axis ticks (4-5 ticks)
    ax2.xaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
    ax2.tick_params(axis='x', which='major', labelsize=14)
    
    # Set y-axis limits to match band structure
    ax2.set_ylim(-2, 12)
    
    # Add title
    if title:
        fig.suptitle(title, fontsize=22, fontweight='bold', y=0.98)
    
    # Use fixed subplot margins to ensure consistent plot dimensions across all structures
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.12, top=0.92)
    
    # Save both PNG and PDF
    fig.savefig(str(output_png), dpi=300)
    fig.savefig(str(output_pdf))
    plt.close(fig)
    
    print(f"  Saved: {output_png.name} and {output_pdf.name}")


def extract_composition_from_id(structure_id: str) -> str:
    """
    Extract composition from structure ID.
    
    Example: 'Ba2P1_s001' -> 'Ba2P1'
    """
    return structure_id.split('_')[0]


def process_structure(
    structure_id: str,
    results_dir: Path,
    output_dir: Path,
    pearson_symbol: Optional[str] = None,
    use_poscar: bool = False
) -> bool:
    """
    Process a single structure and create plots.
    
    Args:
        structure_id: Structure identifier
        results_dir: Directory containing structure data
        output_dir: Output directory for plots
        pearson_symbol: Pearson symbol (if already known from database)
        use_poscar: If True, read POSCAR to determine Pearson symbol
    
    Returns:
        True if successful, False otherwise
    """
    # Locate data directory
    data_dir = results_dir / f"{structure_id}_PHON"
    
    if not data_dir.exists():
        print(f"Error: Directory not found: {data_dir}")
        return False
    
    print(f"\nProcessing: {structure_id}")
    print(f"  Data directory: {data_dir}")
    
    # Check for required files
    band_file = data_dir / "phonon_band.dat"
    dos_file = data_dir / "phonon_dos.dat"
    kpath_file = data_dir / "band_kpath.dat"
    
    if not band_file.exists():
        print(f"  Error: {band_file.name} not found")
        return False
    
    if not dos_file.exists():
        print(f"  Error: {dos_file.name} not found")
        return False
    
    if not kpath_file.exists():
        print(f"  Error: {kpath_file.name} not found")
        print(f"  This file is required for k-path structure and labels")
        print(f"  Please re-run postproc_phonon.py to generate it")
        return False
    
    # Parse data files
    print("  Reading band structure data...")
    try:
        band_data = parse_band_data(band_file)
        print(f"    Bands: {band_data['n_bands']}, Points: {len(band_data['distances'])}")
    except Exception as e:
        print(f"  Error parsing band data: {e}")
        return False
    
    print("  Reading DOS data...")
    try:
        dos_data = parse_dos_data(dos_file)
        print(f"    Elements: {', '.join(dos_data['elements']) if dos_data['elements'] else 'None'}")
    except Exception as e:
        print(f"  Error parsing DOS data: {e}")
        return False
    
    # Get Pearson symbol from POSCAR if requested and not already provided
    if use_poscar and pearson_symbol is None:
        poscar_file = data_dir / "POSCAR"
        if poscar_file.exists():
            print("  Reading POSCAR to determine Pearson symbol...")
            try:
                from pymatgen.core import Structure
                structure = Structure.from_file(poscar_file)
                pearson_symbol = get_pearson_symbol(structure)
                if pearson_symbol:
                    print(f"    Pearson symbol: {pearson_symbol}")
            except Exception as e:
                print(f"    Warning: Could not read POSCAR: {e}")
                pearson_symbol = None
        else:
            print(f"    Warning: POSCAR not found in {data_dir}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plot title
    composition = extract_composition_from_id(structure_id)
    composition_latex = formula_to_latex(composition)
    
    if pearson_symbol:
        title = f"{pearson_symbol}-{composition_latex} Phonon"
    else:
        title = f"{composition_latex} Phonon"
    
    # Create plots
    output_png = output_dir / f"PHON_{structure_id}.png"
    output_pdf = output_dir / f"PHON_{structure_id}.pdf"
    
    print("  Creating plots...")
    try:
        plot_phonon_band_dos(
            band_data=band_data,
            dos_data=dos_data,
            output_png=output_png,
            output_pdf=output_pdf,
            title=title,
            kpath_file=kpath_file
        )
    except Exception as e:
        print(f"  Error creating plots: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Plot phonon band structure and DOS from pre-computed data files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot single structure
  python3 phonon_band_dos_plot.py --results-dir Bin-Ele-HT/candidates/ --db Bin-Ele-HT/stable_electrides.db --structure-ids Ba2P1_s001
  
  # Plot multiple structures
  python3 phonon_band_dos_plot.py --results-dir Ter-Ele-HT/candidates/ --db Ter-Ele-HT/stable_electrides.db --structure-ids K6B1O4_s013 Cs2Al2S3_s013
  
  # Custom output directory
  python3 phonon_band_dos_plot.py --results-dir Bin-Ele-HT/candidates/ --db Bin-Ele-HT/stable_electrides.db --structure-ids Ba2P1_s001 --output-dir ./phonon_plots
  
  # Plot all structures in directory
  python3 phonon_band_dos_plot.py --results-dir Bin-Ele-HT/candidates/ --db Bin-Ele-HT/stable_electrides.db --all
  
  # Use POSCAR for Pearson symbol (instead of database)
  python3 phonon_band_dos_plot.py --results-dir Bin-Ele-HT/candidates/ --structure-ids Ba2P1_s001 --poscar

Directory structure expected:
  {results_dir}/{structure_id}_PHON/phonon_band.dat          (REQUIRED)
  {results_dir}/{structure_id}_PHON/phonon_dos.dat           (REQUIRED)
  {results_dir}/{structure_id}_PHON/band_kpath.dat           (REQUIRED)
  {results_dir}/{structure_id}_PHON/POSCAR                   (optional, for --poscar flag)
  
Database file:
  stable_electrides.db - ASE database with pearson_symbol field (optional)
  
Note: The .dat files are auto-generated by postproc_phonon.py.
      band_kpath.dat contains k-path segment structure and high-symmetry labels.
      If missing, re-run postproc_phonon.py to generate them.
      
Pearson symbol options:
  - Use --db to load Pearson symbols from ASE database
  - Use --poscar to determine Pearson symbols from POSCAR files
  - Use neither for plots without Pearson symbols in titles
        """
    )
    
    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help='Directory containing structure data (e.g., Bin-Ele-HT/candidates/)'
    )
    
    parser.add_argument(
        '--db',
        type=str,
        help='Path to stable_electrides.db for Pearson symbols (optional)'
    )
    
    parser.add_argument(
        '--poscar',
        action='store_true',
        help='Read POSCAR files to determine Pearson symbols (alternative to --db)'
    )
    
    parser.add_argument(
        '--structure-ids',
        type=str,
        nargs='+',
        help='Structure IDs to plot (e.g., Ba2P1_s001 Ca5P3_s021)'
    )
    
    parser.add_argument(
        '--all',
        action='store_true',
        help='Plot all structures found in results-dir'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./phonon_plots',
        help='Output directory for plots (default: ./phonon_plots)'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        sys.exit(1)
    
    output_dir = Path(args.output_dir)
    
    # Determine which structures to process
    if args.all:
        # Find all *_PHON directories
        structure_ids = []
        for phon_dir in sorted(results_dir.glob("*_PHON")):
            if phon_dir.is_dir():
                struct_id = phon_dir.name.replace('_PHON', '')
                structure_ids.append(struct_id)
        
        if not structure_ids:
            print(f"Error: No *_PHON directories found in {results_dir}")
            sys.exit(1)
        
        print(f"Found {len(structure_ids)} structures to plot")
    elif args.structure_ids:
        structure_ids = args.structure_ids
    else:
        print("Error: Must specify either --structure-ids or --all")
        parser.print_help()
        sys.exit(1)
    
    # Load Pearson symbols from database if provided
    pearson_dict = {}
    if args.db:
        db_path = Path(args.db)
        print(f"\nLoading Pearson symbols from database: {args.db}")
        pearson_dict = get_pearson_symbols_from_db(db_path, structure_ids)
        print(f"  Loaded Pearson symbols for {len([v for v in pearson_dict.values() if v])} structures")
    elif args.poscar:
        print("\nPearson symbols will be determined from POSCAR files")
    else:
        print("\nNo database or POSCAR flag provided, Pearson symbols will not be shown in titles")
    
    print("=" * 80)
    print("Phonon Band Structure and DOS Plotting")
    print("=" * 80)
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Structures to plot: {len(structure_ids)}")
    if args.db:
        print(f"Database: {args.db}")
    elif args.poscar:
        print(f"Pearson source: POSCAR files")
    print("=" * 80)
    
    # Process each structure
    successful = 0
    failed = 0
    
    for struct_id in structure_ids:
        pearson = pearson_dict.get(struct_id, '') if args.db else None
        success = process_structure(struct_id, results_dir, output_dir, pearson, use_poscar=args.poscar)
        
        if success:
            successful += 1
        else:
            failed += 1
    
    # Summary
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"Total structures: {len(structure_ids)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"\nPlots saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
