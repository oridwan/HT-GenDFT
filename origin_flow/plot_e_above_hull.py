#!/usr/bin/env python3
"""
Plot histogram of energy_above_hull from pre-screening results.

Usage:
    python3 plot_e_above_hull.py [--input prescreening_stability.json] [--output histogram.png]
"""

import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Plot histogram of energy_above_hull from pre-screening results"
    )
    parser.add_argument(
        '--input',
        type=str,
        default='./VASP_JOBS/prescreening_stability.json',
        help="Path to prescreening_stability.json (default: ./VASP_JOBS/prescreening_stability.json)"
    )
    parser.add_argument(
        '--output',
        type=str,
        default='energy_above_hull_histogram.png',
        help="Output plot filename (default: energy_above_hull_histogram.png)"
    )
    parser.add_argument(
        '--bins',
        type=int,
        default=50,
        help="Number of histogram bins (default: 50)"
    )
    parser.add_argument(
        '--max-energy',
        type=float,
        default=None,
        help="Max energy to display (default: auto)"
    )
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"Error: File not found: {input_path}")
        return
    
    with open(input_path, 'r') as f:
        data = json.load(f)
    
    energies = []
    for result in data.get('results', []):
        e_hull = result.get('energy_above_hull')
        if e_hull is not None:
            energies.append(e_hull)
    
    if not energies:
        print("Error: No energy_above_hull values found in file")
        return
    
    energies = np.array(energies)
    
    print("="*70)
    print("Energy Above Hull Statistics")
    print("="*70)
    print(f"Total structures: {len(energies)}")
    print(f"Min: {energies.min():.4f} eV/atom")
    print(f"Max: {energies.max():.4f} eV/atom")
    print(f"Mean: {energies.mean():.4f} eV/atom")
    print(f"Median: {np.median(energies):.4f} eV/atom")
    print(f"Std: {energies.std():.4f} eV/atom")
    print("")
    
    threshold = data.get('summary', {}).get('hull_threshold', 0.1)
    passed = np.sum(energies < threshold)
    print(f"Threshold: {threshold} eV/atom")
    print(f"Structures below threshold: {passed} ({passed/len(energies)*100:.1f}%)")
    print(f"Structures above threshold: {len(energies)-passed} ({(len(energies)-passed)/len(energies)*100:.1f}%)")
    print("="*70)
    
    _, ax = plt.subplots(figsize=(10, 6))
    
    max_e = args.max_energy if args.max_energy else energies.max()
    
    ax.hist(energies, bins=args.bins, 
            range=(0, max_e),
            edgecolor='black', linewidth=0.5,
            alpha=0.7, color='steelblue')
    
    ax.axvline(threshold, color='red', linestyle='--', linewidth=2, 
               label=f'Threshold = {threshold} eV/atom')
    
    ax.set_xlabel('Energy Above Hull (eV/atom)', fontsize=12)
    ax.set_ylabel('Number of Structures', fontsize=12)
    ax.set_title('Distribution of Energy Above Hull\nMatterSim Pre-screening Results', 
                 fontsize=14, fontweight='bold')
    
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    text_str = f'Total: {len(energies)}\n'
    text_str += f'Mean: {energies.mean():.3f} eV/atom\n'
    text_str += f'Median: {np.median(energies):.3f} eV/atom\n'
    text_str += f'Passed: {passed} ({passed/len(energies)*100:.1f}%)'
    
    ax.text(0.98, 0.97, text_str,
            transform=ax.transAxes,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=10,
            family='monospace')
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {args.output}")


if __name__ == '__main__':
    main()

