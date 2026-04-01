#!/usr/bin/env python3
"""
Phonon Post-Processing Script

Automatically post-process completed phonon calculations using phonopy and pymatgen
Processes structures marked as PHON_DONE in PHONON_JOBS/workflow.json.

Features:
- Load force sets from displacement vasprun.xml files
- Calculate force constants
- Automatic high-symmetry k-path generation (pymatgen HighSymmKpath)
- Phonon band structure calculation and plotting
- Phonon DOS calculation and plotting
- Thermal properties calculation
- Imaginary frequency detection
- All using Python APIs (no command-line tools)

Usage:
    python3 postproc_phonon.py --phonon-jobs ./PHONON_JOBS
    
    python3 postproc_phonon.py --phonon-jobs ./PHONON_JOBS --structure-ids Ba2N_s001
"""

import sys
import json
import argparse
import warnings
import re
from pathlib import Path
from typing import List, Tuple

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt

from phonopy import Phonopy
from phonopy.physical_units import get_physical_units
VaspToTHz = get_physical_units().DefaultToTHz
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.phonopy import get_phonopy_structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import PhononDos
from pymatgen.phonon.plotter import PhononBSPlotter

warnings.filterwarnings('ignore')


def formula_to_latex(formula: str) -> str:
    """Convert chemical formula to LaTeX format with subscripts.
    
    Ignores subscript 1 (e.g., Al1 -> Al, not Al$_{1}$).
    
    Examples:
        Ca5P3 -> Ca$_{5}$P$_{3}$
        Ca7Al1P5 -> Ca$_{7}$AlP$_{5}$
        K6BO4 -> K$_{6}$BO$_{4}$
        Cs6Al2S5 -> Cs$_{6}$Al$_{2}$S$_{5}$
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


# Greek letter mapping for high-symmetry point labels (LaTeX format)
GREEK_LETTERS = {
    'GAMMA': r'\Gamma',
    'DELTA': r'\Delta',
    'LAMBDA': r'\Lambda',
    'SIGMA': r'\Sigma',
    'OMEGA': r'\Omega',
    'PHI': r'\Phi',
    'PSI': r'\Psi',
    'THETA': r'\Theta',
    'ALPHA': r'\Alpha',
    'BETA': r'\Beta',
}


def convert_label_to_latex(label: str) -> str:
    """
    Convert k-point label to LaTeX format for matplotlib rendering.
    
    Args:
        label: Label from pymatgen HighSymmKpath (e.g., "GAMMA", "SIGMA_0")
    
    Returns:
        LaTeX label (e.g., r"$\Gamma$", r"$\Sigma_0$")
    
    Examples:
        "GAMMA" -> r"$\Gamma$"
        "SIGMA_0" -> r"$\Sigma_0$"
        "X" -> "X"
    """
    # Handle subscripts (e.g., SIGMA_0 -> \Sigma_0)
    if '_' in label:
        base, subscript = label.split('_', 1)
    else:
        base = label
        subscript = None
    
    # Convert base to Greek letter if applicable
    if base in GREEK_LETTERS:
        latex_label = GREEK_LETTERS[base]
        # Add subscript if present
        if subscript is not None:
            latex_label = f"{latex_label}_{{{subscript}}}"
        # Wrap in math mode
        return f"${latex_label}$"
    else:
        # Not a Greek letter, return as-is (with subscript if present)
        if subscript is not None:
            return f"${base}_{{{subscript}}}$"
        return label


def load_workflow_database(db_path: Path) -> dict:
    """Load the phonon workflow database."""
    if not db_path.exists():
        raise FileNotFoundError(f"Workflow database not found: {db_path}")
    
    with open(db_path, 'r') as f:
        return json.load(f)


def get_pearson_symbol(structure: Structure) -> str:
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
        return "Unknown"


def check_displacement_files(phon_dir: Path, n_displacements: int) -> Tuple[bool, List[Path]]:
    """
    Check if all displacement vasprun.xml files exist and are valid.
    
    Returns:
        (all_exist, vasprun_paths)
    """
    vasprun_paths = []
    
    for i in range(1, n_displacements + 1):
        disp_id = f"{i:03d}"
        disp_dir = phon_dir / disp_id
        vasprun_path = disp_dir / 'vasprun.xml'
        
        if not vasprun_path.exists():
            print(f"    Warning: Missing vasprun.xml in {disp_id}/")
            return False, []
        
        if vasprun_path.stat().st_size == 0:
            print(f"    Warning: Empty vasprun.xml in {disp_id}/")
            return False, []
        
        vasprun_paths.append(vasprun_path)
    
    return True, vasprun_paths


def extract_forces_from_vasprun(vasprun_paths: List[Path]) -> np.ndarray:
    """
    Extract forces from vasprun.xml files.
    
    Returns:
        forces: numpy array of shape (n_displacements, n_atoms, 3)
    """
    forces_list = []
    
    for vasprun_path in vasprun_paths:
        try:
            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
            forces = vr.ionic_steps[-1]['forces']  # Get forces from last ionic step
            forces_list.append(forces)
        except Exception as e:
            raise RuntimeError(f"Failed to parse {vasprun_path}: {e}")
    
    return np.array(forces_list)


def create_force_constants_with_phonopy(
    structure: Structure,
    phon_dir: Path,
    forces: np.ndarray,
    supercell_matrix: np.ndarray
) -> Phonopy:
    """
    Create force constants using phonopy Python API.
    
    Args:
        structure: Pymatgen Structure (unit cell)
        phon_dir: Path to PHON directory
        forces: Forces from displacements, shape (n_displacements, n_atoms, 3)
        supercell_matrix: Supercell matrix used for phonopy
    
    Returns:
        Phonopy object with force constants
    """
    # Convert pymatgen structure to phonopy
    phonopy_structure = get_phonopy_structure(structure)
    
    # Create Phonopy object
    phonon = Phonopy(
        phonopy_structure,
        supercell_matrix=supercell_matrix,
        primitive_matrix='auto',
        factor=VaspToTHz
    )
    
    # Generate displacements
    phonon.generate_displacements(distance=0.01)
    
    # Set forces
    phonon.forces = forces
    
    # Generate force constants
    phonon.produce_force_constants()
    
    # Save force constants
    phonon.save(filename=str(phon_dir / 'phonopy_params.yaml'), settings={'force_constants': True})
    
    return phonon


def get_band_structure_with_pymatgen(
    phonon: Phonopy,
    structure: Structure,
    npoints: int = 51
) -> PhononBandStructureSymmLine:
    """
    Calculate phonon band structure using pymatgen's automatic k-path.
    
    Args:
        phonon: Phonopy object with force constants
        structure: Pymatgen Structure
        npoints: Number of points between each high-symmetry point
    
    Returns:
        PhononBandStructureSymmLine object
    """
    # Get high-symmetry k-path using pymatgen
    kpath = HighSymmKpath(structure)
    
    # Get k-path in fractional coordinates
    kpath_dict = kpath.kpath
    path = kpath_dict['path']
    kpoints_dict = kpath_dict['kpoints']
    
    # Build paths for phonopy
    qpoints_for_phonopy = []
    for branch in path:
        branch_qpoints = []
        for i in range(len(branch) - 1):
            start_label = branch[i]
            end_label = branch[i + 1]
            
            start_coords = np.array(kpoints_dict[start_label])
            end_coords = np.array(kpoints_dict[end_label])
            
            # Generate points along this segment
            for j in range(npoints):
                frac = j / (npoints - 1)
                qpoint = start_coords + frac * (end_coords - start_coords)
                branch_qpoints.append(qpoint)
        
        qpoints_for_phonopy.append(branch_qpoints)
    
    # Run phonopy band structure calculation
    phonon.run_band_structure(qpoints_for_phonopy, is_band_connection=False)
    
    # Get results
    band_dict = phonon.get_band_structure_dict()
    
    # Flatten all results from different branches
    all_qpoints = []
    all_frequencies = []
    for qpts, freqs in zip(band_dict['qpoints'], band_dict['frequencies'], strict=False):
        all_qpoints.extend(qpts)
        all_frequencies.extend(freqs)
    
    all_qpoints = np.array(all_qpoints)
    all_frequencies = np.array(all_frequencies)  # Shape: (n_qpoints, n_bands)
    
    # Transpose frequencies to (n_bands, n_qpoints) for pymatgen
    frequencies_transposed = all_frequencies.T  # THz
    
    # Create labels_dict with converted LaTeX labels
    labels_dict = {}
    for label in set([pt for branch in path for pt in branch]):
        # Convert label to LaTeX format (GAMMA -> $\Gamma$, etc.)
        display_label = convert_label_to_latex(label)
        labels_dict[display_label] = np.array(kpoints_dict[label])
    
    # Create PhononBandStructureSymmLine
    bs = PhononBandStructureSymmLine(
        qpoints=all_qpoints,
        frequencies=frequencies_transposed,
        lattice=structure.lattice.reciprocal_lattice,
        labels_dict=labels_dict,
        structure=structure
    )
    
    return bs


def get_phonon_dos(
    phonon: Phonopy,
    structure: Structure,
    mesh: List[int] = [20, 20, 20]
) -> Tuple[PhononDos, dict]:
    """
    Calculate phonon DOS using phonopy.
    
    Args:
        phonon: Phonopy object with force constants
        structure: Pymatgen Structure
        mesh: Mesh for DOS calculation
    
    Returns:
        (total_dos, element_dos_dict) where element_dos_dict is {Element: PhononDos}
    """
    # Run mesh calculation
    phonon.run_mesh(mesh, with_eigenvectors=True, is_mesh_symmetry=False)
    
    # Get total DOS
    phonon.run_total_dos()
    total_dos_dict = phonon.get_total_dos_dict()
    
    frequencies = total_dos_dict['frequency_points']  # THz
    densities = total_dos_dict['total_dos']
    
    # Create total PhononDos object
    total_dos = PhononDos(frequencies, densities)
    
    # Get projected DOS
    phonon.run_projected_dos()
    proj_dos_dict = phonon.get_projected_dos_dict()
    
    # Get primitive cell from phonopy
    primitive = phonon.primitive
    n_atoms_primitive = len(primitive)
    
    # Group by element
    from pymatgen.core import Element
    element_dos_dict = {}
    
    for i in range(n_atoms_primitive):
        element = Element.from_Z(primitive.numbers[i])
        site_dos_array = proj_dos_dict['projected_dos'][i]  # Shape: (n_freq,)
        
        if element in element_dos_dict:
            # Add to existing element DOS
            element_dos_dict[element] += site_dos_array
        else:
            # Create new element DOS
            element_dos_dict[element] = site_dos_array.copy()
    
    # Convert to PhononDos objects
    element_dos = {}
    for element, dos_array in element_dos_dict.items():
        element_dos[element] = PhononDos(frequencies, dos_array)
    
    return total_dos, element_dos


def check_imaginary_frequencies(bs: PhononBandStructureSymmLine, tol: float = 0.1) -> Tuple[bool, float, float]:
    """
    Check for imaginary frequencies in phonon band structure.
    
    Args:
        bs: PhononBandStructureSymmLine
        tol: Tolerance in THz for considering frequency as imaginary
    
    Returns:
        (has_imaginary, min_freq, max_imaginary_freq)
    """
    has_imaginary = bs.has_imaginary_freq(tol=tol)
    min_freq, max_freq = bs.min_freq(), bs.max_freq()
    
    # Get minimum frequency value
    min_freq_value = min_freq[1] if min_freq else 0.0
    
    # Find maximum imaginary frequency (most negative)
    max_imaginary = 0.0
    if has_imaginary:
        all_freqs = []
        for branch_freqs in bs.bands:
            all_freqs.extend(branch_freqs)
        all_freqs = np.array(all_freqs).flatten()
        imaginary_freqs = all_freqs[all_freqs < -tol]
        if len(imaginary_freqs) > 0:
            max_imaginary = abs(imaginary_freqs.min())
    
    return has_imaginary, min_freq_value, max_imaginary


def plot_band_structure_and_dos(
    bs: PhononBandStructureSymmLine,
    total_dos: PhononDos,
    element_dos: dict,
    output_path: Path,
    title: str = None
):
    """Plot phonon band structure and DOS side by side with aligned y-axis."""
    # Enable LaTeX text rendering
    plt.rcParams['text.usetex'] = False
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['font.family'] = 'DejaVu Sans'
    
    # Extract band structure data
    bs_plotter = PhononBSPlotter(bs)
    bs_data = bs_plotter.bs_plot_data()
    
    # Create combined plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True,
                                     gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.08})
    
    # Plot band structure on left (ax1)
    all_distances = np.concatenate(bs_data['distances'])
    all_frequencies = np.concatenate([np.array(freqs).T for freqs in bs_data['frequency']], axis=0) # Shape (n_qpoints_total, n_bands)
    
    n_bands = all_frequencies.shape[1]
    
    for band_idx in range(n_bands):
        band_frequencies = all_frequencies[:, band_idx]
        ax1.plot(all_distances, band_frequencies, 'r-', linewidth=1.0, zorder=2)
    
    # Add high-symmetry point markers and labels
    for distance in bs_data['ticks']['distance']:
        ax1.axvline(x=distance, color='black', linestyle='-', linewidth=0.5, zorder=1)
    
    # Clean up tick labels (remove double $$)
    clean_labels = []
    
    for label in bs_data['ticks']['label']:
        # Replace double $$ with single $
        label = label.replace('$$', '$')
        clean_labels.append(label)
    
    ax1.set_xticks(bs_data['ticks']['distance'])
    ax1.set_xticklabels(clean_labels, fontsize=14)
    ax1.set_xlabel('Wave vector', fontsize=16, fontweight='bold')
    ax1.set_ylabel('Frequency (THz)', fontsize=16, fontweight='bold')
    ax1.axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=5)
    ax1.set_xlim(all_distances[0] - 0.0001, all_distances[-1] + 0.0001)
    ax1.tick_params(axis='y', which='major', labelsize=14)
    
    # Plot DOS on right (ax2)
    ax2.plot(total_dos.densities, total_dos.frequencies, 'k-', linewidth=1.0, 
             label='Total', alpha=0.8, zorder=3)
    
    # Plot element-projected DOS with different colors
    if element_dos:
        colors = plt.cm.tab10(np.arange(0, 0.2 * len(element_dos), 0.2))
        for idx, (element, el_dos) in enumerate(element_dos.items()):
            ax2.plot(el_dos.densities, el_dos.frequencies, 
                    linewidth=1.5, label=str(element), alpha=0.8, color=colors[idx], zorder=2)
    
    ax2.set_xlabel('DOS', fontsize=16, fontweight='bold')
    ax2.axhline(y=0, color='black', linestyle='--', linewidth=0.5, alpha=0.5, zorder=5)
    ax2.legend(loc='upper right', fontsize=14, framealpha=0.8, edgecolor='none', facecolor='white')
    ax2.set_xlim(left=0)
    ax2.tick_params(axis='x', which='major', labelsize=14)
    
    # Set y-axis limits to match
    ymin = min(ax1.get_ylim()[0], ax2.get_ylim()[0])
    ymax = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_ylim(ymin, ymax)
    ax2.set_ylim(ymin, ymax)
    
    # Add title
    if title:
        fig.suptitle(title, fontsize=22, fontweight='bold', y=0.98)
    
    # Use fixed subplot margins to ensure consistent plot dimensions across all structures
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.12, top=0.92)
    
    # Save both PNG and PDF
    output_png = output_path
    output_pdf = output_path.with_suffix('.pdf')
    
    fig.savefig(str(output_png), dpi=300)
    print(f"    Saved: {output_png.name}")
    
    fig.savefig(str(output_pdf))
    print(f"    Saved: {output_pdf.name}")
    
    plt.close(fig)


def save_band_structure_data(bs: PhononBandStructureSymmLine, output_path: Path):
    """Save phonon band structure data in gnuplot-friendly format."""
    # Get band structure data
    qpoints = bs.qpoints
    frequencies = bs.bands  # Shape: (n_bands, n_qpoints)
    
    # Get high-symmetry labels and segment structure
    bs_plotter = PhononBSPlotter(bs)
    bs_data = bs_plotter.bs_plot_data()
    
    # Calculate distances
    distances = np.zeros(len(qpoints))
    current_idx = 0
    for segment_distances in bs_data['distances']:
        segment_length = len(segment_distances)
        distances[current_idx:current_idx + segment_length] = segment_distances
        current_idx += segment_length
    
    # Clean and format labels
    tick_labels = []
    for label in bs_data['ticks']['label']:
        label = label.replace('$$', '$')
        
        if '|' in label:
            parts = label.split('|')
            formatted_parts = []
            for part in parts:
                part = part.strip()
                if not part.startswith('$'):
                    part = f'${part}$'
                formatted_parts.append(part)
            label = '|'.join(formatted_parts)
        else:
            if not label.startswith('$'):
                label = f'${label}$'
        
        tick_labels.append(label)
    
    tick_distances = bs_data['ticks']['distance']
    
    lattice = bs.structure.lattice
    lattice_matrix = lattice.matrix
    
    # Get fractional coordinates
    tick_qpoint_indices = []
    for tick_dist in tick_distances:
        idx = np.argmin(np.abs(distances - tick_dist))
        tick_qpoint_indices.append(idx)
    
    # Save k-path metadata
    kpath_file = output_path.parent / 'band_kpath.dat'
    with open(kpath_file, 'w') as f:
        f.write("# K-path metadata for phonon band structure\n")
        f.write("# Contains lattice, segment structure, and high-symmetry point labels\n")
        f.write("#\n")
        f.write("# LATTICE VECTORS (Angstrom):\n")
        f.write("# Format: LATTICE_A/B/C  x  y  z\n")
        f.write("#\n")
        f.write(f"LATTICE_A  {lattice_matrix[0, 0]:12.8f}  {lattice_matrix[0, 1]:12.8f}  {lattice_matrix[0, 2]:12.8f}\n")
        f.write(f"LATTICE_B  {lattice_matrix[1, 0]:12.8f}  {lattice_matrix[1, 1]:12.8f}  {lattice_matrix[1, 2]:12.8f}\n")
        f.write(f"LATTICE_C  {lattice_matrix[2, 0]:12.8f}  {lattice_matrix[2, 1]:12.8f}  {lattice_matrix[2, 2]:12.8f}\n")
        f.write("#\n")
        f.write("# SEGMENT STRUCTURE:\n")
        f.write("# Format: SEGMENT  seg_idx  start_qidx  end_qidx  n_points\n")
        f.write("#\n")
        
        # Write segment
        current_idx = 0
        for seg_idx, segment_distances in enumerate(bs_data['distances']):
            n_points = len(segment_distances)
            start_qidx = current_idx
            end_qidx = current_idx + n_points - 1
            f.write(f"SEGMENT  {seg_idx:3d}  {start_qidx:6d}  {end_qidx:6d}  {n_points:6d}\n")
            current_idx += n_points
        
        f.write("#\n")
        f.write("# HIGH-SYMMETRY POINTS:\n")
        f.write("# Format: TICK  distance  qx  qy  qz  label\n")
        f.write("# Distance is along the path, coordinates are fractional (reciprocal lattice units)\n")
        f.write("#\n")
        
        # Write tick information
        for tick_dist, idx, label in zip(tick_distances, tick_qpoint_indices, tick_labels):
            qpt = qpoints[idx]
            frac = qpt.frac_coords
            f.write(f"TICK     {tick_dist:12.8f}  {frac[0]:12.8f}  {frac[1]:12.8f}  {frac[2]:12.8f}  {label}\n")
    
    # Write band data
    with open(output_path, 'w') as f:
        # Header
        f.write("# Phonon band structure\n")
        f.write("# Column 1: Distance along path\n")
        f.write("# Column 2: q-point (fractional coordinates)\n")
        f.write(f"# Columns 3-{3+len(frequencies)}: Phonon frequencies (THz) for each band\n")
        f.write("#\n")
        f.write(f"# {'Distance':>12}  {'qx':>10} {'qy':>10} {'qz':>10}")
        for i in range(len(frequencies)):
            f.write(f"  {'Band'+str(i+1):>12}")
        f.write("\n")
        
        # Data
        for i, (dist, qpt) in enumerate(zip(distances, qpoints)):
            qcoords = qpt.frac_coords
            f.write(f"  {dist:12.6f}  {qcoords[0]:10.6f} {qcoords[1]:10.6f} {qcoords[2]:10.6f}")
            for band_freqs in frequencies:
                f.write(f"  {band_freqs[i]:12.6f}")
            f.write("\n")
    
    print(f"    Saved band structure data: {output_path.name}")
    print(f"    Saved k-path metadata: {kpath_file.name}")


def save_dos_data(total_dos: PhononDos, element_dos: dict, output_path: Path):
    """Save phonon DOS data in gnuplot-friendly format."""
    frequencies = total_dos.frequencies
    total_densities = total_dos.densities
    
    with open(output_path, 'w') as f:
        # Header
        f.write("# Phonon Density of States (DOS)\n")
        f.write("# Column 1: Frequency (THz)\n")
        f.write("# Column 2: Total DOS\n")
        if element_dos:
            for idx, element in enumerate(element_dos.keys()):
                f.write(f"# Column {idx+3}: {element} projected DOS\n")
        f.write("#\n")
        
        # Column headers
        f.write(f"{'Frequency':>12}  {'Total':>12}")
        if element_dos:
            for element in element_dos.keys():
                f.write(f"  {str(element):>12}")
        f.write("\n")
        
        # Data
        for i, freq in enumerate(frequencies):
            f.write(f"  {freq:12.6f}  {total_densities[i]:12.6f}")
            if element_dos:
                for el_dos in element_dos.values():
                    f.write(f"  {el_dos.densities[i]:12.6f}")
            f.write("\n")
    
    print(f"    Saved DOS data: {output_path.name}")


def calculate_thermal_properties(dos: PhononDos, structure: Structure, output_path: Path):
    """Calculate and save thermal properties."""
    # Temperature range
    temperatures = np.arange(0, 1001, 10)  # 0-1000 K, step 10 K
    
    # Calculate properties
    results = []
    for T in temperatures:
        if T == 0:
            # Only zero-point energy at T=0
            zpe = dos.zero_point_energy(structure)
            results.append({
                'T (K)': T,
                'F (kJ/mol)': zpe / 1000,
                'S (J/K/mol)': 0.0,
                'Cv (J/K/mol)': 0.0,
                'U (kJ/mol)': zpe / 1000
            })
        else:
            F = dos.helmholtz_free_energy(T, structure)
            S = dos.entropy(T, structure)
            Cv = dos.cv(T, structure)
            U = dos.internal_energy(T, structure)
            
            results.append({
                'T (K)': T,
                'F (kJ/mol)': F / 1000,
                'S (J/K/mol)': S,
                'Cv (J/K/mol)': Cv,
                'U (kJ/mol)': U / 1000
            })
    
    # Save to file
    with open(output_path, 'w') as f:
        f.write("# Thermal properties from phonon DOS\n")
        f.write("# T: Temperature (K)\n")
        f.write("# F: Helmholtz free energy (kJ/mol)\n")
        f.write("# S: Entropy (J/K/mol)\n")
        f.write("# Cv: Heat capacity at constant volume (J/K/mol)\n")
        f.write("# U: Internal energy (kJ/mol)\n")
        f.write("#\n")
        f.write(f"{'T (K)':>10} {'F (kJ/mol)':>15} {'S (J/K/mol)':>15} {'Cv (J/K/mol)':>15} {'U (kJ/mol)':>15}\n")
        
        for res in results:
            f.write(f"{res['T (K)']:10.1f} {res['F (kJ/mol)']:15.6f} {res['S (J/K/mol)']:15.6f} "
                   f"{res['Cv (J/K/mol)']:15.6f} {res['U (kJ/mol)']:15.6f}\n")
    
    print(f"    Saved thermal properties: {output_path.name}")


def process_structure(
    struct_id: str,
    sdata: dict,
    output_summary: dict
) -> bool:
    """
    Process a single structure's phonon calculation.
    
    Returns:
        True if successful, False otherwise
    """
    phon_dir = Path(sdata['phon_dir'])
    n_displacements = sdata['n_displacements']
    supercell_matrix = np.array(sdata.get('supercell_matrix', [[2, 0, 0], [0, 2, 0], [0, 0, 2]]))
    
    print(f"\nProcessing: {struct_id}")
    print(f"  PHON directory: {phon_dir}")
    print(f"  Displacements: {n_displacements}")
    print(f"  Supercell matrix: {np.diag(supercell_matrix).astype(int)}")
    
    # Check if all vasprun.xml files exist
    print("  Checking displacement files...")
    all_exist, vasprun_paths = check_displacement_files(phon_dir, n_displacements)
    
    if not all_exist:
        print(f"  Skipping {struct_id}: Missing or invalid displacement files")
        output_summary[struct_id] = {'status': 'FAILED', 'reason': 'Missing displacement files'}
        return False
    
    print(f"    All {n_displacements} vasprun.xml files found")
    
    # Load structure from POSCAR
    poscar_path = phon_dir / 'POSCAR'
    if not poscar_path.exists():
        print(f"  Skipping {struct_id}: POSCAR not found")
        output_summary[struct_id] = {'status': 'FAILED', 'reason': 'POSCAR not found'}
        return False
    
    structure = Structure.from_file(str(poscar_path))
    composition = structure.composition.reduced_formula
    pearson_symbol = get_pearson_symbol(structure)
    print(f"  Structure: {composition} ({len(structure)} atoms, {pearson_symbol})")
    
    # Extract forces
    print("  Extracting forces from vasprun.xml files...")
    try:
        forces = extract_forces_from_vasprun(vasprun_paths)
        print(f"    Extracted forces: shape {forces.shape}")
    except Exception as e:
        print(f"  Error extracting forces: {e}")
        output_summary[struct_id] = {'status': 'FAILED', 'reason': f'Force extraction failed: {e}'}
        return False
    
    # Create force constants with phonopy
    print("  Creating force constants...")
    try:
        phonon = create_force_constants_with_phonopy(structure, phon_dir, forces, supercell_matrix)
        print("    Force constants created successfully")
    except Exception as e:
        print(f"  Error creating force constants: {e}")
        output_summary[struct_id] = {'status': 'FAILED', 'reason': f'Force constants failed: {e}'}
        return False
    
    # Calculate band structure
    print("  Calculating phonon band structure...")
    try:
        bs = get_band_structure_with_pymatgen(phonon, structure, npoints=51)
        print("    Band structure calculated")
    except Exception as e:
        print(f"  Error calculating band structure: {e}")
        output_summary[struct_id] = {'status': 'FAILED', 'reason': f'Band structure failed: {e}'}
        return False
    
    # Check for imaginary frequencies
    has_imaginary, min_freq, max_imaginary = check_imaginary_frequencies(bs, tol=0.1)
    if has_imaginary:
        print(f"    WARNING: Imaginary frequencies detected!")
        print(f"    Minimum frequency: {min_freq:.3f} THz")
        print(f"    Maximum imaginary: {max_imaginary:.3f} THz")
    else:
        print(f"    All frequencies positive (min: {min_freq:.3f} THz)")
    
    # Save band structure data
    band_data_path = phon_dir / 'phonon_band.dat'
    save_band_structure_data(bs, band_data_path)
    
    # Calculate DOS
    print("  Calculating phonon DOS...")
    try:
        total_dos, element_dos = get_phonon_dos(phonon, structure, mesh=[20, 20, 20])
        print("    DOS calculated")
    except Exception as e:
        print(f"  Warning: DOS calculation failed: {e}")
        total_dos = None
        element_dos = None
    
    # Plot combined band structure and DOS
    if total_dos:
        # Save DOS data
        dos_data_path = phon_dir / 'phonon_dos.dat'
        save_dos_data(total_dos, element_dos, dos_data_path)
        
        # Plot combined band structure and DOS
        combined_plot_path = phon_dir / 'phonon_band_dos.png'
        composition_latex = formula_to_latex(composition)
        plot_title = f"{pearson_symbol}-{composition_latex} Phonon"
        plot_band_structure_and_dos(bs, total_dos, element_dos, combined_plot_path, 
                                     title=plot_title)
        
        # Calculate thermal properties
        thermal_path = phon_dir / 'thermal.dat'
        calculate_thermal_properties(total_dos, structure, thermal_path)
    
    # Save summary
    output_summary[struct_id] = {
        'status': 'SUCCESS',
        'has_imaginary_freq': has_imaginary,
        'min_frequency_THz': float(min_freq),
        'max_imaginary_THz': float(max_imaginary) if has_imaginary else 0.0,
        'n_displacements': n_displacements,
        'supercell_matrix': supercell_matrix.tolist()
    }
    
    print(f"  {struct_id} processed successfully")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Post-process phonon calculations using phonopy + pymatgen APIs"
    )
    parser.add_argument(
        '--phonon-jobs',
        type=str,
        default='./PHONON_JOBS',
        help="Path to PHONON_JOBS directory (default: ./PHONON_JOBS)"
    )
    parser.add_argument(
        '--structure-ids',
        type=str,
        nargs='+',
        default=None,
        help="Specific structure IDs to process (default: all PHON_DONE structures)"
    )
    parser.add_argument(
        '--output-summary',
        type=str,
        default='phonon_postproc_summary.json',
        help="Output summary file (default: phonon_postproc_summary.json)"
    )
    
    args = parser.parse_args()
    
    phonon_jobs_dir = Path(args.phonon_jobs).expanduser()
    db_path = phonon_jobs_dir / 'workflow.json'
    
    if not phonon_jobs_dir.exists():
        print(f"ERROR: PHONON_JOBS directory not found: {phonon_jobs_dir}")
        sys.exit(1)
    
    # Load workflow database
    print("="*70)
    print("Phonon Post-Processing")
    print("="*70)
    print(f"PHONON_JOBS directory: {phonon_jobs_dir}")
    
    db = load_workflow_database(db_path)
    
    # Filter structures
    if args.structure_ids:
        structures_to_process = {
            sid: sdata for sid, sdata in db['structures'].items()
            if sid in args.structure_ids and sdata['state'] == 'PHON_DONE'
        }
        print(f"User-specified structures: {len(args.structure_ids)}")
    else:
        structures_to_process = {
            sid: sdata for sid, sdata in db['structures'].items()
            if sdata['state'] == 'PHON_DONE'
        }
        print(f"All PHON_DONE structures: {len(structures_to_process)}")
    
    if not structures_to_process:
        print("\nNo PHON_DONE structures found to process!")
        sys.exit(0)
    
    print(f"Structures to process: {len(structures_to_process)}")
    print("="*70)
    
    # Process each structure
    output_summary = {}
    success_count = 0
    failed_count = 0
    
    for struct_id, sdata in structures_to_process.items():
        try:
            success = process_structure(struct_id, sdata, output_summary)
            if success:
                success_count += 1
            else:
                failed_count += 1
        except Exception as e:
            print(f"  Unexpected error processing {struct_id}: {e}")
            output_summary[struct_id] = {'status': 'FAILED', 'reason': f'Unexpected error: {e}'}
            failed_count += 1
    
    # Save summary
    summary_path = phonon_jobs_dir / args.output_summary
    with open(summary_path, 'w') as f:
        json.dump(output_summary, f, indent=2)
    
    # Final summary
    print("\n" + "="*70)
    print("Post-Processing Complete")
    print("="*70)
    print(f"Total structures: {len(structures_to_process)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {failed_count}")
    print(f"\nSummary saved to: {summary_path}")
    
    # Print structures with imaginary frequencies
    imaginary_structures = [
        sid for sid, summary in output_summary.items()
        if summary.get('status') == 'SUCCESS' and summary.get('has_imaginary_freq', False)
    ]
    
    if imaginary_structures:
        print("\nStructures with imaginary frequencies:")
        for sid in imaginary_structures:
            summary = output_summary[sid]
            print(f"  {sid}: max imaginary = {summary['max_imaginary_THz']:.3f} THz")
    
    print("="*70)


if __name__ == '__main__':
    main()