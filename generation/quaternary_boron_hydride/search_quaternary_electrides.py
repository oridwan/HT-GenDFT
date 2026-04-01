#!/usr/bin/env python3
"""
Search for potential quaternary hydride compositions based on valence electron counting.
Generates a list of compositions to pass to MatterGen for structure generation.
Quaternary hydrides follow the formula A_l C_m BH_n where:
- A and C are different electropositive metals (Group I, II, or III), with A ≠ C
- B (Boron) has fixed valence of +3
- H has fixed valence of -1
"""

import json
from math import gcd
from typing import List, Dict, Tuple


VALENCE_ELECTRONS = {
    # Group I
    'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1,
    # Group II
    'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2,
    # Group III (including lanthanides, excluding actinides)
    'Al': 3, 'Ga': 3, 'In': 3, 'Tl': 3,
    'Sc': 3, 'Y': 3,
    # Lanthanides
    'La': 3, 'Ce': 3, 'Pr': 3, 'Nd': 3, 'Pm': 3, 'Sm': 3, 'Eu': 3,
    'Gd': 3, 'Tb': 3, 'Dy': 3, 'Ho': 3, 'Er': 3, 'Tm': 3, 'Yb': 3, 'Lu': 3,
    # Group IV
    'C': 4, 'Si': 4, 'Ge': 4, 'Sn': 4, 'Pb': 4,
    # Group V
    'N': 3, 'P': 3, 'As': 3, 'Sb': 3, 'Bi': 3,  # Common negative valence
    # Group VI
    'O': 2, 'S': 2, 'Se': 2, 'Te': 2, 'Po': 2,  # Common negative valence
    # Group VII
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'At': 1,  # Common negative valence
    # Transition metals
    'Fe': 2,  # Fixed valence for quaternary hydrides
}


GROUP_A_B = ['Li', 'Na', 'K', 'Rb', 'Cs',     # Group I
             'Be', 'Mg', 'Ca', 'Sr', 'Ba',     # Group II
             'Sc', 'Y', 'Gd']                   # Lanthanides

# Fixed elements for quaternary hydride system
GROUP_B = ['B']
GROUP_H = ['H']
B_VALENCE = 3
H_VALENCE = -1

def count_excess_electrons(A: str, l: int, C: str, m: int, n: int, h: int) -> float:
    """
    Count excess valence electrons in quaternary composition A_l C_m B_n H_h.
    Where B is Boron with fixed valence of +3.
    
    Positive excess = potential hydride
    Formula: n_excess = valence(A)*l + valence(C)*m + valence(B)*n + valence(H)*h
    Note: H has valence -1, so it reduces the total
    """
    valence_A = VALENCE_ELECTRONS[A]
    valence_C = VALENCE_ELECTRONS[C]
    
    excess = (valence_A * l + valence_C * m + B_VALENCE * n + H_VALENCE * h)
    
    return excess


def search_quaternary_hydrides(
    group_ab: List[str] = GROUP_A_B,
    max_atoms: int = 20,
    excess_electron_range: Tuple[float, float] = (0.1, 4.0),
    max_compositions: int = -1
) -> List[Dict]:
    """
    Search for potential quaternary hydride compositions A_l C_m B_n H_h.
    Where B is Boron with fixed valence of +3.
    
    Args:
        group_ab: List of elements for positions A and C (must be different, A ≠ C)
        max_atoms: Maximum number of atoms in unit cell
        excess_electron_range: (min, max) excess electrons for hydride
        max_compositions: Maximum number of compositions to generate (-1 for all)
        
    Returns:
        List of valid compositions as dictionaries
    """
    valid_compositions = []
    n_count = 0
    
    print("Searching for potential quaternary hydrides (A-C-B-H)...")
    print(f"Constraint: A ≠ C")
    print(f"Variable: A_l C_m B_n H_h with B valence = {B_VALENCE}, H valence = {H_VALENCE}")
    print(f"Excess electron range: {excess_electron_range[0]} to {excess_electron_range[1]}")
    print(f"Max atoms per composition: {max_atoms}")
    print("="*70)
    
    for A in group_ab:
        for C in group_ab:
            # Constraint: A must be different from C
            if A == C:
                continue
                
            for l in range(1, max_atoms):
                for m in range(1, max_atoms):
                    for n in range(1, max_atoms):
                        for h in range(1, max_atoms):
                            # Total atoms: l (A) + m (C) + n (B) + h (H)
                            if l + m + n + h > max_atoms:
                                continue
                            
                            # Check if l, m, n, h are in reduced form (gcd = 1)
                            g = gcd(gcd(gcd(l, m), n), h)
                            if g != 1:
                                continue
                            
                            excess = count_excess_electrons(A, l, C, m, n, h)
                            
                            if excess_electron_range[0] < excess <= excess_electron_range[1]:
                                composition = {
                                    A: l,
                                    C: m,
                                    'B': n,
                                    'H': h
                                }
                            
                                formula = f"{A}{l}{C}{m}B{n}H{h}"
                                
                                valid_compositions.append({
                                    'formula': formula,
                                    'composition': composition,
                                    'excess_electrons': excess,
                                    'total_atoms': l + m + n + h,
                                    'elements': [A, C, 'B', 'H']
                                })
                                
                                n_count += 1
                                
                                if n_count % 100 == 0:
                                    print(f"Found {n_count} compositions... (Latest: {formula})")
                                
                                if max_compositions > 0 and n_count >= max_compositions:
                                    print(f"\nReached maximum of {max_compositions} compositions.")
                                    return valid_compositions
    
    print(f"\nTotal valid compositions found: {len(valid_compositions)}")
    return valid_compositions


def save_compositions(compositions: List[Dict], output_file: str = "quaternary_boron_hydride_compositions.json"):
    """Save compositions to JSON file."""
    with open(output_file, 'w') as f:
        json.dump(compositions, f, indent=2)
    print(f"\nCompositions saved to: {output_file}")
    
    txt_file = output_file.replace('.json', '.txt')
    with open(txt_file, 'w') as f:
        for comp in compositions:
            f.write(f"{comp['formula']}\n")
    print(f"Formulas list saved to: {txt_file}")


def print_statistics(compositions: List[Dict]):
    """Print statistics about found compositions."""
    print("\n" + "="*70)
    print("STATISTICS")
    print("="*70)
    
    from collections import Counter
    atom_counts = Counter(c['total_atoms'] for c in compositions)
    print("\nDistribution by total atoms:")
    for atoms in sorted(atom_counts.keys()):
        print(f"  {atoms:2d} atoms: {atom_counts[atoms]:4d} compositions")
    
    excess_counts = Counter(round(c['excess_electrons'], 1) for c in compositions)
    print("\nDistribution by excess electrons:")
    for excess in sorted(excess_counts.keys()):
        print(f"  {excess:4.1f} e⁻: {excess_counts[excess]:4d} compositions")
    
    print("\nExample compositions (sorted by excess electrons):")
    sorted_comps = sorted(compositions, key=lambda x: x['excess_electrons'], reverse=True)
    for comp in sorted_comps[:10]:
        print(f"  {comp['formula']:15s} | {comp['excess_electrons']:4.1f} e⁻ | {comp['total_atoms']:2d} atoms")


def main():
    """Main execution."""
    print("="*70)
    print("QUATERNARY HYDRIDE COMPOSITION SEARCH (A-C-B-H)")
    print("="*70)
    
    compositions = search_quaternary_hydrides(
        max_atoms=20,
        excess_electron_range=(0.1, 1.1),
        max_compositions=-1
    )
    
    print_statistics(compositions)
    
    save_compositions(compositions, "quaternary_boron_hydride_compositions.json")
    
    print("\n" + "="*70)
    print(f" Search completed! Found {len(compositions)} potential quaternary hydride compositions.")
    print("="*70)
    print("\nNext steps:")
    print("1. Review quaternary_hydride_compositions.json")
    print("2. Transfer to cluster: scp quaternary_hydride_compositions.json HPC_HOST:~/SOFT/mattergen_test/")
    print("3. Run generation: sbatch generate_quaternary_csp.sh")
    print("="*70)
    print("Note: This composition uses Boron (B) with valence +3 instead of Fe with valence +2")
    print("="*70)


if __name__ == "__main__":
    main()

