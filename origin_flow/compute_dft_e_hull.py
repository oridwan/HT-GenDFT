#!/usr/bin/env python3
"""
Compute DFT energy_above_hull for VASP-relaxed structures.

Uses:
- VASP-PBE energies from completed relaxations (vasprun.xml)
- MP DFT-PBE raw energies for competing phases (queried via legacy pymatgen MPRester)
- Cached MP entries from prescreening to minimize API calls

Important:
- Uses legacy pymatgen.ext.matproj.MPRester (not mp_api.client) for complete GGA entry coverage
- The new mp_api.client misses many stable GGA phases needed for accurate hull calculations
- Strict filtering: only accepts entries with '-GGA' or '-GGA+U' suffix in entry_id
- Uses ComputedEntry.uncorrected_energy (raw DFT, no anion/composition corrections)
- This matches VASP energies which also have no MP-style corrections
- Optionally filters to pure GGA-PBE only (excludes GGA+U) when --pure-pbe is used

Output: dft_stability_results.json with DFT-level hull analysis


python3 origin_flow/compute_dft_e_hull.py   --vasp-jobs "$PWD/VASP-out"   --db workflow.json   --prescreen-results "$PWD/prescreen/VASP_JOBS/prescreening_stability.json"   --mattersim-cache "$PWD/prescreen/VASP_JOBS/mp_mattersim.json"   --output dft_stability_results.json   --hull-threshold 0.1   --outlier-threshold 0.5 --pure-pbe
"""

import os
import sys
import json
import argparse
import time
import random
from pathlib import Path
from collections import defaultdict

from pymatgen.core import Composition
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry

# Use legacy pymatgen MPRester for complete GGA entry retrieval
# The new mp_api.client.MPRester misses many stable GGA phases
try:
    from pymatgen.ext.matproj import MPRester
except ImportError:
    print("ERROR: pymatgen package with MPRester is required for DFT hull calculations")
    print("Install with: pip install pymatgen")
    sys.exit(1)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde

# Use non-interactive backend for cluster
mpl.use('Agg')

def load_workflow_database(db_path):
    """Load workflow.json database."""
    with open(db_path, 'r') as f:
        return json.load(f)


def get_vasp_energy_from_relax(relax_dir):
    """
    Extract final energy from VASP relaxation.
    
    Tries multiple sources in order:
    1. vasprun.xml (preferred, has structure)
    2. OUTCAR (fallback if vasprun.xml missing/corrupted) + CONTCAR for structure
    
    Checks electronic convergence before returning energy (same as workflow_manager.py).
    If electronic SCF is not converged, returns (None, None).
    
    Note: Uses converged_electronic (not converged) to match workflow_manager.py behavior,
    which allows structures with NSW=30 that may not reach full ionic convergence.
    
    Returns:
        tuple: (energy_per_atom, structure) or (None, None) if failed/not converged
    """
    relax_dir = Path(relax_dir)
    vasprun_path = relax_dir / 'vasprun.xml'
    outcar_path = relax_dir / 'OUTCAR'
    contcar_path = relax_dir / 'CONTCAR'
    
    if vasprun_path.exists():
        try:
            vr = Vasprun(str(vasprun_path), parse_dos=False, parse_eigen=False)
            
            # Check electronic convergence
            if not vr.converged_electronic:
                print(f"    Warning: VASP electronic SCF not converged (from vasprun.xml)")
                return None, None
            
            final_energy = vr.final_energy  # Total energy in eV
            structure = vr.final_structure
            
            # Check for float overflow (NaN or inf values)
            if not np.isfinite(final_energy):
                print(f"    Warning: Float overflow in vasprun.xml (energy is NaN or inf)")
                return None, None
            
            n_atoms = len(structure)
            energy_per_atom = final_energy / n_atoms
            
            # Double-check energy_per_atom is also finite
            if not np.isfinite(energy_per_atom):
                print(f"    Warning: Float overflow in energy_per_atom calculation")
                return None, None
            
            return energy_per_atom, structure
            
        except Exception as e:
            print(f"    Warning: Could not parse vasprun.xml: {e}")
            print(f"    Trying OUTCAR fallback...")
    
    # Fallback: Try OUTCAR + CONTCAR
    if outcar_path.exists():
        try:
            electronic_converged = False
            with open(outcar_path, 'r') as f:
                outcar_content = f.read()
                if 'aborting loop because EDIFF is reached' in outcar_content:
                    electronic_converged = True
            
            if not electronic_converged:
                print(f"    Warning: VASP electronic SCF not converged (from OUTCAR)")
                return None, None
            
            # Extract final energy from OUTCAR
            final_energy = None
            with open(outcar_path, 'r') as f:
                for line in f:
                    if 'free energy    TOTEN' in line:
                        try:
                            parts = line.split('=')
                            if len(parts) >= 2:
                                energy_str = parts[1].split()[0]
                                # Handle overflow markers like '*******'
                                if '*' in energy_str:
                                    continue
                                final_energy = float(energy_str)
                                # Check for float overflow
                                if not np.isfinite(final_energy):
                                    final_energy = None
                        except (ValueError, IndexError):
                            continue
            
            if final_energy is None:
                print(f"    Warning: Could not extract energy from OUTCAR")
                return None, None
            
            if not contcar_path.exists():
                print(f"    Warning: CONTCAR not found, cannot get structure")
                return None, None
            
            from pymatgen.core import Structure
            structure = Structure.from_file(str(contcar_path))
            n_atoms = len(structure)
            energy_per_atom = final_energy / n_atoms
            
            # Check for float overflow in energy_per_atom
            if not np.isfinite(energy_per_atom):
                print(f"    Warning: Float overflow in energy_per_atom calculation from OUTCAR")
                return None, None
            
            print(f"    Energy extracted from OUTCAR (vasprun.xml not available)")
            return energy_per_atom, structure
            
        except Exception as e:
            print(f"    Error parsing OUTCAR: {e}")
            return None, None
    print(f"    Error: Neither vasprun.xml nor OUTCAR found/readable")
    return None, None


def save_dft_cache(new_entries_data, dft_cache_file):
    """
    Save DFT entries to single global cache file.
    
    Args:
        new_entries_data: List of dict entries to add to cache
        dft_cache_file: Path to global DFT cache file (mp_vaspdft.json)
    """
    dft_cache_file = Path(dft_cache_file)
    dft_cache_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Load existing cache
    existing_cache = {}
    if dft_cache_file.exists():
        with open(dft_cache_file, 'r') as f:
            cache_list = json.load(f)
            # Convert to dict for deduplication
            for item in cache_list:
                key = (item.get('chemsys', ''), item['entry_id'])
                existing_cache[key] = item
    
    # Add new entries (avoid duplicates)
    for new_item in new_entries_data:
        key = (new_item['chemsys'], new_item['entry_id'])
        existing_cache[key] = new_item
    
    # Save back to file
    all_cache_data = list(existing_cache.values())
    with open(dft_cache_file, 'w') as f:
        json.dump(all_cache_data, f, indent=2)
    
    print(f"    Updated DFT cache: {dft_cache_file} ({len(all_cache_data)} total entries)")


def load_dft_cache(chemsys, dft_cache_file, require_compounds=True):
    """
    Load DFT entries for a chemical system from single global cache file.
    
    Args:
        chemsys: Chemical system (e.g., 'Li-B-N')
        dft_cache_file: Path to global DFT cache file (mp_vaspdft.json)
        require_compounds: If True, check for both elementals and compounds.
                          If False, only check for elementals (for single-element systems).
    
    Returns:
        List of PDEntry objects for this chemsys, or None if cache doesn't exist or incomplete
    """
    dft_cache_file = Path(dft_cache_file)
    
    if not dft_cache_file.exists():
        return None
    
    try:
        with open(dft_cache_file, 'r') as f:
            cache_list = json.load(f)
        
        # Extract entries for this chemsys and all subsystems
        elements = chemsys.split('-')
        entries = []
        
        for item in cache_list:
            item_chemsys = item.get('chemsys', '')
            item_elements = sorted(item_chemsys.split('-'))
            
            # Include if chemsys matches OR if it's a subset
            if set(item_elements).issubset(set(elements)):
                comp = Composition(item['composition'])
                entry = PDEntry(
                    composition=comp,
                    energy=item['energy'],
                    name=item['entry_id']
                )
                entries.append(entry)
        
        if not entries:
            return None
        
        # For multi-element systems, verify we have compounds, not just elementals
        if require_compounds and len(elements) > 1:
            has_compound = any(len(e.composition.elements) > 1 for e in entries)
            if not has_compound:
                # Only have elemental phases, need to query for compounds
                return None
        
        return entries
    except Exception as e:
        print(f"    Warning: Could not load DFT cache: {e}")
        return None


def get_mp_stable_phases_dft(chemsys, mp_api_key, dft_cache_file, pure_pbe=False, use_cache=True):
    """
    Get MP phases using RAW GGA DFT energies (no MP corrections) for phase diagram construction.
    
    Strategy:
    - Use legacy pymatgen MPRester.get_entries_in_chemsys() to get ALL entries
    - Filter for GGA entries (entry_id ending with '-GGA') or GGA+U (ending with '-GGA+U')
    - Extract uncorrected_energy from ComputedEntry objects
    - Let PhaseDiagram determine stability based on GGA energies
    - Compatible with VASP raw DFT energies for hull analysis
    
    Note: The legacy pymatgen API returns more complete GGA entries compared to mp_api.client,
          which misses many stable GGA phases (e.g., mp-540703-GGA for Cs2S, mp-2654-GGA for Al2S3).
    
    Args:
        chemsys: Chemical system (e.g., 'Li-B-N')
        mp_api_key: Materials Project API key
        dft_cache_file: Path to global DFT cache file (mp_vaspdft.json)
        pure_pbe: If True, filter to GGA-PBE only (exclude PBE+U)
                  If False (default), accept both PBE and PBE+U
        use_cache: If True, check cache before querying MP. If False, force MP query (for pre-population).
    
    Returns:
        List of PDEntry objects with GGA/GGA+U energies for this chemical system
    """
    # Check DFT cache first (only if use_cache=True)
    if use_cache:
        elements_list = chemsys.split('-')
        is_single_element = (len(elements_list) == 1)
        
        # For single-element systems, don't require compounds (there are none by definition)
        cached_entries = load_dft_cache(chemsys, dft_cache_file, require_compounds=not is_single_element)
        
        if cached_entries is not None and len(cached_entries) > 0:
            # Verify we have all terminal phases (elemental phases)
            elements = set(elements_list)
            cached_elements = set()
            n_elementals = 0
            n_compounds = 0
            for entry in cached_entries:
                if len(entry.composition.elements) == 1:  # Pure elemental phase
                    cached_elements.add(str(entry.composition.elements[0]))
                    n_elementals += 1
                else:
                    n_compounds += 1
            
            if elements == cached_elements or elements.issubset(cached_elements):
                print(f"    Loaded {len(cached_entries)} stable DFT entries from cache "
                      f"({n_elementals} elementals, {n_compounds} compounds)")
                return cached_entries
            else:
                print(f"    Cache incomplete (missing elementals: {elements - cached_elements}), requerying...")
        elif not is_single_element:
            # Multi-element system but no compounds found in cache
            print(f"    Cache incomplete (no compounds found), querying...")
    
    # Query MP for phases using legacy pymatgen MPRester
    print(f"    Querying MP for phases in {chemsys} (using legacy API for complete GGA coverage)...")
    
    elements = chemsys.split('-')
    
    try:
        # Use legacy pymatgen MPRester which returns complete GGA entries
        mpr = MPRester(mp_api_key)
        
        # Get ALL entries in this chemical system
        # This automatically includes all subsystems (elementals, binaries, etc.)
        computed_entries = mpr.get_entries_in_chemsys(elements)
        
        print(f"    Retrieved {len(computed_entries)} entries from MP (filtering for GGA...)")
        
        # Process entries: filter for GGA only (entry_id ending with '-GGA' or '-GGA+U')
        entries = []
        new_entries_data = []
        seen_entries = {}  # Track by entry_id to avoid duplicates
        
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
            
            # Get uncorrected energy (raw DFT, no MP corrections)
            try:
                uncorrected_energy = comp_entry.uncorrected_energy
            except AttributeError:
                # Fallback: use corrected energy if uncorrected not available
                uncorrected_energy = comp_entry.energy
            
            composition = comp_entry.composition
            
            # Determine chemsys for this phase
            doc_elements = sorted([str(el) for el in composition.elements])
            doc_chemsys = '-'.join(doc_elements)
            
            # Extract base MP ID (e.g., 'mp-540703' from 'mp-540703-GGA')
            parts = entry_id.split('-')
            if len(parts) >= 2:
                mp_id = parts[0] + '-' + parts[1]
            else:
                mp_id = entry_id
            
            # Create PDEntry with uncorrected energy
            entry = PDEntry(
                composition=composition,
                energy=uncorrected_energy,  # Total energy (eV)
                name=entry_id
            )
            entries.append(entry)
            
            # Store for caching
            new_entries_data.append({
                'chemsys': doc_chemsys,
                'composition': {str(el): float(amt) for el, amt in composition.items()},
                'energy': uncorrected_energy,
                'entry_id': entry_id,
                'mp_id': mp_id
            })
            
            # Mark as seen
            seen_entries[entry_id] = True
        
        print(f"    Filtered to {len(entries)} GGA entries (strict '-GGA'/'-GGA+U' suffix)")
        
        # Verify we have terminal (elemental) phases
        elements_found = set()
        for entry in entries:
            if len(entry.composition.elements) == 1:  # Pure elemental phase
                elements_found.add(str(entry.composition.elements[0]))
        
        expected_elements = set(elements)
        if elements_found != expected_elements:
            missing = expected_elements - elements_found
            print(f"    WARNING: Missing terminal phases for elements: {sorted(missing)}")
        else:
            print(f"      All terminal phases present: {sorted(elements_found)}")
        
        # Save to DFT cache
        if new_entries_data:
            save_dft_cache(new_entries_data, dft_cache_file)
        
        return entries
        
    except Exception as e:
        print(f"    Error processing MP entries for {chemsys}: {e}")
        import traceback
        traceback.print_exc()
        return []


def compute_dft_hull(target_entry, mp_entries):
    """
    Compute energy_above_hull using DFT energies with PDEntry.
    
    Args:
        target_entry: PDEntry for target structure (VASP total energy in eV)
        mp_entries: List of PDEntry from MP (DFT total energies in eV)
    
    Returns:
        tuple: (decomposition_products, energy_above_hull in eV/atom)
               decomposition_products is dict of {PDEntry: fraction}
    
    Note:
        PDEntry expects total energy (eV), not energy per atom.
        PhaseDiagram internally normalizes to per-atom basis for hull calculation.
        The returned e_hull is in eV/atom.
    """
    try:
        pd = PhaseDiagram(mp_entries)
        decomp, e_hull = pd.get_decomp_and_e_above_hull(target_entry, allow_negative=True)
        return decomp, e_hull
    except Exception as e:
        # Provide detailed error info
        target_elements = set(str(el) for el in target_entry.composition.elements)
        mp_elements_found = set()
        for entry in mp_entries:
            mp_elements_found.update(str(el) for el in entry.composition.elements)
        missing_elements = target_elements - mp_elements_found
        
        print(f"    Error computing phase diagram: {e}")
        if missing_elements:
            print(f"    Missing terminal entries for elements: {sorted(missing_elements)}")
            print(f"    Target composition: {target_entry.composition.reduced_formula}")
            print(f"    Elements found in MP database: {sorted(mp_elements_found)}")
        print(f"    FAILED: Could not compute phase diagram")
        return None, None


def load_mattersim_results(prescreen_json):
    """
    Load MatterSim prescreening results for comparison.
    
    Args:
        prescreen_json: Path to prescreening_stability.json
    
    Returns:
        dict: {structure_id: {'energy_above_hull': float, 'passed_prescreening': bool}}
    """
    if not prescreen_json or not Path(prescreen_json).exists():
        return {}
    
    with open(prescreen_json, 'r') as f:
        data = json.load(f)
    
    mattersim_lookup = {}
    for result in data.get('results', []):
        if result.get('energy_above_hull') is not None:
            mattersim_lookup[result['structure_id']] = {
                'energy_above_hull': result['energy_above_hull'],
                'passed_prescreening': result.get('passed_prescreening', False)
            }
    
    return mattersim_lookup


def match_and_analyze_hulls(dft_results, mattersim_lookup, threshold=0.1, outlier_threshold=0.5):
    """
    Match DFT and MatterSim hull results and compute statistics.
    
    Filters outliers (structures with unreasonably high DFT energy_above_hull) using
    absolute threshold to avoid skewing statistics and plots with failed/unreasonable structures.
    
    Args:
        dft_results: List of DFT result dicts
        mattersim_lookup: Dict from load_mattersim_results
        threshold: E_hull threshold for stability (eV/atom)
        outlier_threshold: Absolute DFT E_hull threshold for outlier detection (eV/atom, default: 0.5)
    
    Returns:
        dict: Statistics and matched structure data (outliers filtered)
    """
    matched = []
    
    for result in dft_results:
        struct_id = result['structure_id']
        dft_e_hull = result.get('energy_above_hull')
        
        if dft_e_hull is None or struct_id not in mattersim_lookup:
            continue
        
        ms_data = mattersim_lookup[struct_id]
        
        matched.append({
            'structure_id': struct_id,
            'composition': result['composition'],
            'chemsys': result['chemsys'],
            'mattersim_e_hull': ms_data['energy_above_hull'],
            'dft_e_hull': dft_e_hull,
            'passed_prescreen': ms_data['passed_prescreening']
        })
    
    if not matched:
        return None
    
    # Filter outliers using absolute threshold on DFT energy_above_hull
    matched_filtered = []
    skipped_outliers = []
    
    for entry in matched:
        dft_e = entry['dft_e_hull']
        
        # Only use absolute threshold (structures passed prescreening should be near threshold)
        if dft_e > outlier_threshold:
            skipped_outliers.append(entry)
        else:
            matched_filtered.append(entry)
    
    if skipped_outliers:
        print(f"\n  Filtered {len(skipped_outliers)} outliers for hull comparison "
              f"(DFT E_hull > {outlier_threshold:.1f} eV/atom)")
        print(f"  Analyzing {len(matched_filtered)} structures after filtering")
    
    # Calculate statistics on filtered data
    mattersim_vals = np.array([d['mattersim_e_hull'] for d in matched_filtered])
    dft_vals = np.array([d['dft_e_hull'] for d in matched_filtered])
    passed = np.array([d['passed_prescreen'] for d in matched_filtered])
    
    # Correlation and error metrics
    correlation = np.corrcoef(mattersim_vals, dft_vals)[0, 1]
    mae = np.mean(np.abs(mattersim_vals - dft_vals))
    rmse = np.sqrt(np.mean((mattersim_vals - dft_vals)**2))
    
    # Pre-screening accuracy
    dft_stable = dft_vals < threshold
    
    true_positives = int(np.sum(passed & dft_stable))
    false_positives = int(np.sum(passed & ~dft_stable))
    true_negatives = int(np.sum(~passed & ~dft_stable))
    false_negatives = int(np.sum(~passed & dft_stable))
    
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    # Computational savings (based on original matched count before outlier filtering)
    total = len(matched)
    total_filtered = len(matched_filtered)
    passed_count = true_positives + false_positives
    filtered_count = total_filtered - passed_count
    savings = filtered_count / total_filtered if total_filtered > 0 else 0
    
    return {
        'matched_structures': matched_filtered,
        'skipped_outliers': skipped_outliers,
        'statistics': {
            'total_matched': total,
            'total_analyzed': total_filtered,
            'n_outliers_filtered': len(skipped_outliers),
            'outlier_threshold': outlier_threshold,
            'threshold': threshold,
            'correlation': correlation,
            'mae': mae,
            'rmse': rmse,
            'true_positives': true_positives,
            'false_positives': false_positives,
            'true_negatives': true_negatives,
            'false_negatives': false_negatives,
            'precision': precision,
            'recall': recall,
            'f1_score': f1_score,
            'passed_count': passed_count,
            'filtered_count': filtered_count,
            'computational_savings': savings
        }
    }


def compute_robust_limits(values, low_pct=1.0, high_pct=99.0, padding_frac=0.08, min_span=0.2):
    """
    Compute axis limits using percentiles to avoid extreme tails dominating plots.

    This only affects visualization, not the underlying statistics.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return -1.0, 1.0

    low = float(np.percentile(values, low_pct))
    high = float(np.percentile(values, high_pct))

    if high <= low:
        center = float(np.median(values))
        half = max(min_span / 2.0, 0.05)
        return center - half, center + half

    span = high - low
    padding = max(span * padding_frac, min_span * 0.05)
    return low - padding, high + padding


def plot_hull_comparison(hull_data, output_prefix='hull_comparison'):
    """
    Create scatter and residual plots comparing MatterSim vs DFT hull energies.
    
    Args:
        hull_data: Dict from match_and_analyze_hulls
        output_prefix: Prefix for output PNG files
    """
    matched = hull_data['matched_structures']
    stats = hull_data['statistics']
    
    mattersim_vals = np.array([d['mattersim_e_hull'] for d in matched])
    dft_vals = np.array([d['dft_e_hull'] for d in matched])
    
    # ===== Scatter Plot =====
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Calculate point density for color mapping
    x = dft_vals
    y = mattersim_vals
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    
    # Scale z to represent actual number of data points
    z = z * len(x)
    
    # Create custom colormap
    colors_list = ['cyan', 'dodgerblue', 'black']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list("density", colors_list, N=n_bins)
    
    # Sort points by density so densest are plotted last
    idx = z.argsort()
    x_sorted, y_sorted, z_sorted = x[idx], y[idx], z[idx]
    
    # Determine robust plot range (display-only clipping for readability)
    all_vals = np.concatenate([mattersim_vals, dft_vals])
    plot_min, plot_max = compute_robust_limits(all_vals, low_pct=1.0, high_pct=99.0)
    n_axis_clipped = int(np.sum(
        (dft_vals < plot_min) | (dft_vals > plot_max) |
        (mattersim_vals < plot_min) | (mattersim_vals > plot_max)
    ))
    
    # Create scatter plot with density-based coloring
    scatter = ax.scatter(x_sorted, y_sorted, c=z_sorted, cmap=cmap, 
                        norm=mpl.colors.LogNorm(), s=20, marker='s', 
                        edgecolors='none')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    
    # Perfect agreement line
    ax.plot([plot_min, plot_max], [plot_min, plot_max], 'r--', 
            linewidth=2, alpha=0.7, label='Perfect agreement (y=x)')
    
    # Stability threshold line
    threshold = stats['threshold']
    ax.axhline(y=threshold, color='green', linestyle=':', linewidth=2, 
               alpha=0.6, label=f'Stability threshold ({threshold} eV/atom)')
    ax.axvline(x=threshold, color='green', linestyle=':', linewidth=2, alpha=0.6)
    
    # Labels and title
    ax.set_xlabel('VASP DFT E_above_hull (eV/atom)', fontsize=16, fontweight='bold')
    ax.set_ylabel('MatterSim E_above_hull (eV/atom)', fontsize=16, fontweight='bold')
    ax.set_title('MatterSim vs VASP DFT Hull Energy Comparison', fontsize=18, fontweight='bold')
    
    # Statistics text box
    n_outliers = stats.get('n_outliers_filtered', 0)
    stats_text = (
        f"N = {stats['total_analyzed']}\n"
        f"R = {stats['correlation']:.4f}\n"
        f"MAE = {stats['mae']:.4f} eV/atom\n"
        f"RMSE = {stats['rmse']:.4f} eV/atom\n"
        f"Precision = {stats['precision']:.2%}\n"
        f"Recall = {stats['recall']:.2%}"
    )
    if n_outliers > 0:
        stats_text += f"\nOutliers excluded: {n_outliers}"
    if n_axis_clipped > 0:
        stats_text += f"\nAxis clipping: {n_axis_clipped} pts (1-99%)"
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=13,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Grid and legend
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.minorticks_on()
    ax.grid(True, which='minor', alpha=0.15, linestyle=':')
    ax.legend(loc='lower right', fontsize=12)
    
    # Equal aspect
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(plot_min, plot_max)
    ax.set_ylim(plot_min, plot_max)
    
    # Increase tick label font sizes
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    plt.tight_layout()
    scatter_file = f"{output_prefix}_scatter.png"
    plt.savefig(scatter_file, dpi=300, bbox_inches='tight')
    print(f"  Scatter plot: {scatter_file}")
    plt.close()
    
    # ===== Residual Plot =====
    fig, ax = plt.subplots(figsize=(12, 10))
    
    residuals = mattersim_vals - dft_vals
    
    # Calculate point density for color mapping
    x_res = dft_vals
    y_res = residuals
    xy_res = np.vstack([x_res, y_res])
    z_res = gaussian_kde(xy_res)(xy_res)
    z_res = z_res * len(x_res)
    
    # Create custom colormap
    colors_list_res = ['lime', 'forestgreen', 'black']
    cmap_res = LinearSegmentedColormap.from_list("density_residual", colors_list_res, N=100)
    
    # Sort points by density
    idx_res = z_res.argsort()
    x_res_sorted, y_res_sorted, z_res_sorted = x_res[idx_res], y_res[idx_res], z_res[idx_res]
    
    # Create scatter plot with density-based coloring
    scatter_res = ax.scatter(x_res_sorted, y_res_sorted, c=z_res_sorted, cmap=cmap_res, 
                            norm=mpl.colors.LogNorm(), s=20, marker='s', 
                            edgecolors='none')
    
    # Add colorbar
    cbar_res = plt.colorbar(scatter_res, ax=ax)
    
    # Reference lines
    ax.axhline(y=0, color='r', linestyle='--', linewidth=2, alpha=0.7, 
               label='Zero residual (perfect agreement)')
    
    # Mean residual
    mean_residual = np.mean(residuals)
    ax.axhline(y=mean_residual, color='orange', linestyle='-', linewidth=2, 
               alpha=0.7, label=f'Mean residual = {mean_residual:+.4f} eV/atom')
    
    # ±0.05 eV/atom bands
    ax.axhline(y=0.05, color='green', linestyle=':', linewidth=1.5, alpha=0.6, 
               label='±0.05 eV/atom')
    ax.axhline(y=-0.05, color='green', linestyle=':', linewidth=1.5, alpha=0.6)
    
    # Labels and title
    ax.set_xlabel('VASP DFT E_above_hull (eV/atom)', fontsize=16, fontweight='bold')
    ax.set_ylabel('Residual (MatterSim - VASP) (eV/atom)', fontsize=16, fontweight='bold')
    ax.set_title('Hull Energy Residuals vs VASP DFT', fontsize=18, fontweight='bold')
    
    # Keep residual panel readable around the dense region (symmetric around zero)
    res_low, res_high = compute_robust_limits(residuals, low_pct=1.0, high_pct=99.0)
    y_limit = max(abs(res_low), abs(res_high), 0.2)
    ax.set_xlim(plot_min, plot_max)
    ax.set_ylim(-y_limit, y_limit)
    
    # Increase tick label font sizes
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Grid and legend
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.minorticks_on()
    ax.grid(True, which='minor', alpha=0.15, linestyle=':')
    ax.legend(loc='best', fontsize=12)
    
    plt.tight_layout()
    residual_file = f"{output_prefix}_residuals.png"
    plt.savefig(residual_file, dpi=300, bbox_inches='tight')
    print(f"  Residual plot: {residual_file}")
    plt.close()


def print_hull_comparison_summary(hull_data):
    """Print summary of hull comparison analysis."""
    stats = hull_data['statistics']
    
    print("\n" + "="*70)
    print("MatterSim vs VASP DFT Hull Energy Comparison")
    print("="*70)
    print(f"Matched structures: {stats['total_matched']}")
    
    n_outliers = stats.get('n_outliers_filtered', 0)
    if n_outliers > 0:
        print(f"Outliers filtered: {n_outliers} (E_hull > {stats.get('outlier_threshold', 2.0):.1f} eV/atom)")
        print(f"Structures analyzed: {stats['total_analyzed']}")
    
    print(f"\nCorrelation & Error:")
    print(f"  Pearson correlation: {stats['correlation']:.4f}")
    print(f"  Mean Absolute Error: {stats['mae']:.4f} eV/atom")
    print(f"  RMSE: {stats['rmse']:.4f} eV/atom")
    
    print(f"\nPre-screening Accuracy (threshold = {stats['threshold']} eV/atom):")
    print(f"  True Positives:  {stats['true_positives']:4d} (correct: stable → passed)")
    print(f"  False Positives: {stats['false_positives']:4d} (error: unstable → passed)")
    print(f"  True Negatives:  {stats['true_negatives']:4d} (correct: unstable → filtered)")
    print(f"  False Negatives: {stats['false_negatives']:4d} (error: stable → filtered)")
    
    print(f"\nMetrics:")
    print(f"  Precision: {stats['precision']:.2%} (of passed structures, how many are stable)")
    print(f"  Recall:    {stats['recall']:.2%} (of stable structures, how many passed)")
    print(f"  F1 Score:  {stats['f1_score']:.4f}")
    
    print(f"\nComputational Savings:")
    print(f"  Total structures: {stats['total_analyzed']}")
    print(f"  Passed pre-screening: {stats['passed_count']} ({stats['passed_count']/stats['total_analyzed']:.1%})")
    print(f"  Filtered: {stats['filtered_count']} ({stats['computational_savings']:.1%} VASP computations saved)")
    
    if stats['false_negatives'] > 0:
        print(f"\n  Warning: {stats['false_negatives']} stable structures were filtered by MatterSim!")
    
    print("="*70)


def main():
    parser = argparse.ArgumentParser(
        description="Compute DFT energy_above_hull for VASP-relaxed structures"
    )
    parser.add_argument(
        '--vasp-jobs',
        type=str,
        default='./VASP_JOBS',
        help="VASP jobs directory (default: ./VASP_JOBS)"
    )
    parser.add_argument(
        '--db',
        type=str,
        default='workflow.json',
        help="Workflow database (default: workflow.json in VASP_JOBS)"
    )
    parser.add_argument(
        '--mp-api-key',
        type=str,
        default=None,
        help="Materials Project API key (default: from MP_API_KEY env)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='dft_stability_results.json',
        help="Output JSON file (default: dft_stability_results.json)"
    )
    parser.add_argument(
        '--prescreen-results',
        type=str,
        default=None,
        help="Optional: prescreening_stability.json to filter structures"
    )
    parser.add_argument(
        '--pure-pbe',
        action='store_true',
        help="Filter MP entries to pure GGA-PBE only (exclude PBE+U). "
             "Default: accept both PBE and PBE+U for accurate phase diagrams. "
             "Use this flag only if your VASP calculations use pure PBE without +U corrections."
    )
    parser.add_argument(
        '--hull-threshold',
        type=float,
        default=0.1,
        help="E_hull threshold for stability analysis (eV/atom, default: 0.1)"
    )
    parser.add_argument(
        '--outlier-threshold',
        type=float,
        default=0.5,
        help="E_hull outlier threshold for plot filtering (eV/atom, default: 0.5). "
             "Structures with E_hull above this are excluded from comparison plots to avoid skewing."
    )
    parser.add_argument(
        '--mattersim-cache',
        type=str,
        default='mp_mattersim.json',
        help="MatterSim MP cache file (default: mp_mattersim.json in VASP_JOBS). "
             "Used to extract chemical systems for consistent reference phases."
    )
    
    args = parser.parse_args()
    
    vasp_jobs = Path(args.vasp_jobs).expanduser()
    db_path = Path(args.db).expanduser()
    
    if not db_path.is_absolute():
        db_path = vasp_jobs / args.db
    
    output_path = vasp_jobs / args.output
    
    # Get MP API key
    mp_api_key = args.mp_api_key or os.environ.get('MP_API_KEY')
    if not mp_api_key:
        print("ERROR: MP_API_KEY not found in environment or arguments")
        print("Set it with: export MP_API_KEY=your_key")
        return 1
    
    if len(mp_api_key) != 32:
        print("ERROR: MP API key should be 32 characters (new API)")
        print("Get a new key from: https://next-gen.materialsproject.org/api")
        return 1
    
    # Set up cache paths (single global cache files)
    mp_mattersim_cache_path = Path(args.mattersim_cache).expanduser()
    if not mp_mattersim_cache_path.is_absolute():
        mp_mattersim_cache_path = vasp_jobs / args.mattersim_cache
    
    mp_dft_cache = vasp_jobs / "mp_vaspdft.json"
    
    print("="*70)
    print("DFT Energy Above Hull Calculation")
    print("="*70)
    print(f"VASP jobs directory: {vasp_jobs}")
    print(f"Workflow database: {db_path}")
    print(f"MP MatterSim cache: {mp_mattersim_cache_path}")
    print(f"MP DFT cache: {mp_dft_cache}")
    print(f"Output file: {output_path}")
    
    print(f"Energy source: MP GGA phases (using legacy pymatgen MPRester for complete coverage)")
    print(f"Energy corrections: NONE (using raw DFT energies)")
    print("  VASP: raw DFT energies from vasprun.xml")
    print("  MP: ComputedEntry.uncorrected_energy (no anion/composition corrections)")
    print("  Filtering: Only entries with '-GGA' or '-GGA+U' suffix (strict)")
    print("  Phase diagram stability determined by pymatgen using GGA energies")
    
    if args.pure_pbe:
        print(f"Functional filtering: Pure GGA-PBE only (PBE+U/R2SCAN/SCAN excluded)")
    else:
        print(f"Functional filtering: Mixed PBE/PBE+U (MP recommended methodology)")
    
    print("="*70 + "\n")
    
    # Load MatterSim cache to extract chemical systems for consistent reference phases
    print("="*70)
    print("Loading MatterSim Cache for Reference Phase Consistency")
    print("="*70)
    
    if not mp_mattersim_cache_path.exists():
        print(f"ERROR: MatterSim cache not found: {mp_mattersim_cache_path}")
        print("This file is required to ensure DFT and MatterSim hulls use same reference phases.")
        print("Run prescreen.py first to generate this file.")
        return 1
    
    try:
        with open(mp_mattersim_cache_path, 'r') as f:
            mattersim_cache = json.load(f)
        
        # Extract unique chemical systems from MatterSim cache
        mattersim_chemsys = set()
        for item in mattersim_cache:
            chemsys = item.get('chemsys', '')
            if chemsys:
                mattersim_chemsys.add(chemsys)
        
        print(f"Loaded MatterSim cache: {len(mattersim_cache)} reference phases")
        print(f"Unique chemical systems: {len(mattersim_chemsys)}")
        print(f"  {sorted(mattersim_chemsys)}")
        print()
        
        # Pre-populate DFT cache with GGA energies for ALL MatterSim chemical systems
        print("Pre-populating DFT cache with GGA energies for MatterSim chemical systems...")
        print("(This ensures DFT and MatterSim hulls use same reference phase set)")
        print()
        
        for chemsys in sorted(mattersim_chemsys):
            print(f"  Fetching GGA phases for {chemsys}...")
            sys.stdout.flush()
            try:
                _ = get_mp_stable_phases_dft(chemsys, mp_api_key, mp_dft_cache, 
                                            pure_pbe=args.pure_pbe, use_cache=False)
                print()
            except Exception as e:
                print(f"    ERROR: {e}\n")
                sys.stdout.flush()
        
        print("="*70)
        print("DFT cache pre-population complete")
        print("="*70 + "\n")
        
    except Exception as e:
        print(f"ERROR: Could not load MatterSim cache: {e}")
        return 1
    
    # Load workflow database
    print("Loading workflow database...")
    if not db_path.exists():
        print(f"ERROR: Workflow database not found: {db_path}")
        return 1
    
    db = load_workflow_database(db_path)
    print(f"Found {len(db['structures'])} structures in database\n")
    
    # Load pre-screening results if provided
    passed_structures = None
    if args.prescreen_results:
        prescreen_path = Path(args.prescreen_results).expanduser()
        
        # Handle relative paths: try as-is first, then try prepending vasp_jobs
        if not prescreen_path.is_absolute() and not prescreen_path.exists():
            # Try prepending vasp_jobs only if just a filename
            alt_path = vasp_jobs / args.prescreen_results
            if alt_path.exists():
                prescreen_path = alt_path
        
        if prescreen_path.exists():
            print(f"Loading pre-screening filter: {prescreen_path}")
            with open(prescreen_path, 'r') as f:
                prescreen_data = json.load(f)
            
            passed_structures = set()
            for result in prescreen_data.get('results', []):
                if result.get('passed_prescreening', False):
                    passed_structures.add(result['structure_id'])
            
            print(f"Will only process {len(passed_structures)} structures that passed pre-screening\n")
        else:
            print(f"Warning: Pre-screening file not found: {prescreen_path}")
            print("Processing all structures\n")
    
    # Scan for completed VASP relaxations
    print("Scanning for completed VASP relaxations...")
    completed_structures = []
    
    for struct_id, sdata in db['structures'].items():
        # Filter by pre-screening if provided
        if passed_structures is not None and struct_id not in passed_structures:
            continue
        
        state = sdata['state']
        
        # Check if VASP relaxation completed
        completed_states = ['RELAX_DONE', 'RELAX_TMOUT', 'SC_RUNNING', 'SC_DONE',
                          'PARCHG_RUNNING', 'PARCHG_DONE', 'PARCHG_FAILED', 'PARCHG_SKIPPED',
                          'ELF_RUNNING', 'ELF_DONE', 'ELF_FAILED']
        
        if state in completed_states:
            completed_structures.append(struct_id)
    
    print(f"Found {len(completed_structures)} completed structures\n")
    
    if not completed_structures:
        print("No completed structures found. Run VASP workflow first.")
        return 0
    
    # Group by chemical system for efficient MP queries
    print("Grouping structures by chemical system...")
    structures_by_chemsys = defaultdict(list)
    
    for struct_id in completed_structures:
        sdata = db['structures'][struct_id]
        chemsys = sdata.get('chemsys')
        
        if not chemsys:
            # Fallback: try to get from composition
            comp = Composition(sdata['composition'])
            elements = sorted([str(el) for el in comp.elements])
            chemsys = '-'.join(elements)
        
        structures_by_chemsys[chemsys].append(struct_id)
    
    print(f"Found {len(structures_by_chemsys)} unique chemical systems\n")
    
    # Process each chemical system
    results = []
    
    processed = 0
    failed = 0
    
    for chemsys, struct_ids in sorted(structures_by_chemsys.items()):
        print(f"\nProcessing {chemsys} ({len(struct_ids)} structures)...")
        
        # Get MP stable DFT entries once per chemical system (uses pre-populated cache)
        print(f"  Querying MP for stable phases...")
        try:
            mp_entries = get_mp_stable_phases_dft(chemsys, mp_api_key, mp_dft_cache, pure_pbe=args.pure_pbe, use_cache=True)
            
            if args.pure_pbe:
                print(f"  Retrieved {len(mp_entries)} GGA phases from MP (pure -GGA suffix only)")
            else:
                print(f"  Retrieved {len(mp_entries)} GGA phases from MP (-GGA and -GGA+U suffixes)")
        except Exception as e:
            print(f"  ERROR: Could not get MP entries: {e}")
            failed += len(struct_ids)
            continue
        
        # Process each structure in this chemical system
        for struct_id in struct_ids:
            sdata = db['structures'][struct_id]
            relax_dir = Path(sdata['relax_dir'])
            
            print(f"\n  {struct_id}:")
            
            # Get VASP energy
            energy_per_atom, structure = get_vasp_energy_from_relax(relax_dir)
            
            if energy_per_atom is None:
                print(f"    FAILED: Could not extract VASP energy (electronic not converged or parsing error)")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'vasp_energy_per_atom': None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'error': 'VASP electronic SCF not converged or vasprun.xml parsing failed'
                })
                continue
            
            print(f"    VASP energy: {energy_per_atom:.6f} eV/atom")
            
            # Create PDEntry for target structure
            composition = structure.composition
            total_energy = energy_per_atom * len(structure)
            target_entry = PDEntry(
                composition=composition,
                energy=total_energy,
                name=f"vasp_{struct_id}"
            )
            
            # Compute DFT hull
            decomp, e_hull = compute_dft_hull(target_entry, mp_entries)
            
            if e_hull is None:
                print(f"    FAILED: Could not compute phase diagram")
                failed += 1
                results.append({
                    'structure_id': struct_id,
                    'composition': sdata['composition'],
                    'chemsys': chemsys,
                    'vasp_energy_per_atom': float(energy_per_atom) if energy_per_atom is not None else None,
                    'energy_above_hull': None,
                    'is_stable': None,
                    'decomposition': None,
                    'error': 'Phase diagram calculation failed'
                })
                continue
            
            print(f"    DFT E_hull: {e_hull:.6f} eV/atom")
            
            # Convert decomposition products to JSON-serializable format
            decomp_list = []
            if decomp:
                for entry, fraction in decomp.items():
                    decomp_list.append({
                        'formula': entry.composition.reduced_formula,
                        'fraction': float(fraction),
                        'energy_per_atom': float(entry.energy_per_atom)
                    })
                # Print decomposition for debugging
                if e_hull > 0.001:
                    decomp_str = " + ".join([f"{d['fraction']:.3f} {d['formula']}" for d in decomp_list])
                    print(f"    Decomposes to: {decomp_str}")
            
            is_stable = bool(e_hull < 0.001)
            processed += 1
            
            results.append({
                'structure_id': struct_id,
                'composition': sdata['composition'],
                'chemsys': chemsys,
                'vasp_energy_per_atom': float(energy_per_atom),
                'energy_above_hull': float(e_hull),
                'is_stable': is_stable,
                'decomposition': decomp_list if decomp_list else None,
                'error': None
            })
    
    # Save results
    print("\n" + "="*70)
    print("Saving DFT hull results...")
    
    output_data = {
        'summary': {
            'total_structures': len(completed_structures),
            'passed_prescreening': processed,
            'failed_prescreening': failed,
            'hull_threshold': args.hull_threshold,
            'energy_reference': 'VASP-DFT-PBE',
            'mp_phase_selection': 'all_gga_entries (no is_stable filter, stability determined by PhaseDiagram)',
            'mp_functional_filter': 'pure_pbe' if args.pure_pbe else 'mixed_pbe_pbeU',
            'reference_phase_source': f'mp_mattersim.json (chemsys-based, ensures MatterSim and DFT use same reference phases)'
        },
        'results': results
    }
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"DFT results saved to: {output_path}")
    print("="*70)
    print(f"\nDFT Hull Summary:")
    print(f"  Total structures: {len(completed_structures)}")
    print(f"  Successfully processed: {processed}")
    print(f"  Failed: {failed}")
    print("="*70)
    
    # ===== Hull Comparison Analysis (MatterSim vs DFT) =====
    if args.prescreen_results:
        prescreen_path = Path(args.prescreen_results).expanduser()
        
        # Handle relative paths
        if not prescreen_path.is_absolute() and not prescreen_path.exists():
            alt_path = vasp_jobs / args.prescreen_results
            if alt_path.exists():
                prescreen_path = alt_path
        
        if prescreen_path.exists():
            print("\n" + "="*70)
            print("Performing MatterSim vs DFT Hull Comparison...")
            print("="*70)
            
            # Load MatterSim results
            mattersim_lookup = load_mattersim_results(prescreen_path)
            print(f"Loaded {len(mattersim_lookup)} MatterSim prescreening results")
            
            # Match and analyze (with outlier filtering)
            hull_data = match_and_analyze_hulls(
                results, 
                mattersim_lookup, 
                threshold=args.hull_threshold,
                outlier_threshold=args.outlier_threshold
            )
            
            if hull_data:
                # Print summary
                print_hull_comparison_summary(hull_data)
                
                # Generate plots
                print("\nGenerating comparison plots...")
                hull_plot_prefix = output_path.parent / "hull_comparison"
                plot_hull_comparison(hull_data, output_prefix=str(hull_plot_prefix))
                
                # Save hull comparison JSON
                hull_comparison_file = output_path.parent / "hull_comparison.json"
                hull_comparison_data = {
                    'summary': hull_data['statistics'],
                    'confusion_matrix': {
                        'true_positives': hull_data['statistics']['true_positives'],
                        'false_positives': hull_data['statistics']['false_positives'],
                        'true_negatives': hull_data['statistics']['true_negatives'],
                        'false_negatives': hull_data['statistics']['false_negatives']
                    },
                    'matched_structures': hull_data['matched_structures']
                }
                
                with open(hull_comparison_file, 'w') as f:
                    json.dump(hull_comparison_data, f, indent=2)
                
                print(f"\nHull comparison saved to: {hull_comparison_file}")
                print("="*70 + "\n")
            else:
                print("  Warning: No matching structures found for hull comparison")
                print("="*70 + "\n")
        else:
            print("\n" + "="*70)
            print(f"Note: Prescreening results not found at {prescreen_path}")
            print("Skipping MatterSim vs DFT hull comparison")
            print("="*70 + "\n")
    else:
        print("\n" + "="*70)
        print("Note: No prescreening results provided (--prescreen-results)")
        print("Skipping MatterSim vs DFT hull comparison")
        print("="*70 + "\n")
    
    return 0


if __name__ == '__main__':
    exit(main())
