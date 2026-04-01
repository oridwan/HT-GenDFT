#!/usr/bin/env python3
"""
Extract CIF files for structures that passed pre-screening.

Reads prescreening_stability.json and extracts the corresponding CIF files
from the original results directory (where generated_crystals_cif.zip files are stored).
"""

import os
import sys
import json
import zipfile
import argparse
from pathlib import Path
from collections import defaultdict


def extract_passed_structures(json_file, results_dir, output_dir):
    """
    Extract CIF files for structures that passed pre-screening.
    
    Args:
        json_file: Path to prescreening_stability.json
        results_dir: Path to directory with composition folders (e.g., Li10B1N3_structures/)
        output_dir: Path to output directory for extracted CIF files
    """
    json_file = Path(json_file).expanduser()
    results_dir = Path(results_dir).expanduser()
    output_dir = Path(output_dir).expanduser()
    
    if not json_file.exists():
        print(f"ERROR: {json_file} not found")
        return 1
    
    if not results_dir.exists():
        print(f"ERROR: Results directory {results_dir} not found")
        return 1
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("Extract Stable Structures from Pre-screening Results")
    print("="*70)
    print(f"Input JSON: {json_file}")
    print(f"Results directory: {results_dir}")
    print(f"Output directory: {output_dir}")
    print("="*70 + "\n")
    
    # Load prescreening results
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    summary = data.get('summary', {})
    results = data.get('results', [])
    
    print(f"Total structures in prescreening: {len(results)}")
    print(f"Passed pre-screening: {summary.get('passed_prescreening', 'N/A')}")
    print(f"Hull threshold: {summary.get('hull_threshold', 'N/A')} eV/atom\n")
    
    # Filter passed structures
    passed_structures = [r for r in results if r.get('passed_prescreening', False)]
    
    if not passed_structures:
        print("No structures passed pre-screening!")
        return 0
    
    print(f"Found {len(passed_structures)} structures that passed pre-screening\n")
    
    # Group by composition for efficient extraction
    by_composition = defaultdict(list)
    for struct in passed_structures:
        struct_id = struct['structure_id']
        composition = struct['composition']
        by_composition[composition].append(struct_id)
    
    print(f"Unique compositions: {len(by_composition)}\n")
    print("="*70)
    print("Extracting CIF files...")
    print("="*70 + "\n")
    
    extracted_count = 0
    missing_count = 0
    
    for comp, struct_ids in by_composition.items():
        comp_dir = results_dir / f"{comp}_structures"
        zip_path = comp_dir / "generated_crystals_cif.zip"
        
        if not zip_path.exists():
            print(f"WARNING: {comp}: ZIP file not found at {zip_path}")
            missing_count += len(struct_ids)
            continue
        
        print(f"{comp}: Extracting {len(struct_ids)} structure(s)")
        
        try:
            with zipfile.ZipFile(zip_path, 'r') as zf:
                # Get sorted list of CIF files (same order as prescreen.py)
                all_cifs = sorted([f for f in zf.namelist() if f.endswith('.cif')])
                
                for struct_id in struct_ids:
                    # Parse structure index from ID (e.g., Ce1Co5_s001 -> 1)
                    idx_str = struct_id.split('_s')[-1]
                    try:
                        idx = int(idx_str) - 1  # Convert to 0-based index
                    except ValueError:
                        print(f"  WARNING: Invalid structure ID format: {struct_id}")
                        missing_count += 1
                        continue
                    
                    if idx < 0 or idx >= len(all_cifs):
                        print(f"  WARNING: Index {idx+1} out of range for {struct_id} (ZIP has {len(all_cifs)} files)")
                        missing_count += 1
                        continue
                    
                    # Get the actual CIF filename at this index
                    original_cif_name = all_cifs[idx]
                    
                    # Extract to output directory with the structure_id as filename
                    output_file = output_dir / f"{struct_id}.cif"
                    with zf.open(original_cif_name) as src:
                        with open(output_file, 'wb') as dst:
                            dst.write(src.read())
                    
                    extracted_count += 1
                    print(f"    {struct_id}.cif (from {original_cif_name})")
        
        except Exception as e:
            print(f"  ERROR: Failed to extract from {zip_path}: {e}")
            missing_count += len(struct_ids)
    
    print("\n" + "="*70)
    print("Extraction Complete")
    print("="*70)
    print(f"Successfully extracted: {extracted_count} CIF files")
    print(f"Missing/Failed: {missing_count} CIF files")
    print(f"Output directory: {output_dir}")
    print("="*70)
    
    # Create a summary file
    summary_file = output_dir / "extraction_summary.json"
    summary_data = {
        "source_json": str(json_file),
        "results_directory": str(results_dir),
        "output_directory": str(output_dir),
        "total_passed": len(passed_structures),
        "extracted": extracted_count,
        "missing": missing_count,
        "structures": [
            {
                "structure_id": s['structure_id'],
                "composition": s['composition'],
                "energy_above_hull": s['energy_above_hull'],
                "is_stable": s['is_stable']
            }
            for s in passed_structures
        ]
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    print(f"\nSummary saved to: {summary_file}")
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Extract CIF files for structures that passed pre-screening",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python3 extract_stable_structs.py \\
    --json ./VASP_JOBS/prescreening_stability.json \\
    --results-dir ./generated_crystals \\
    --output-dir ./stable_structures
        """
    )
    
    parser.add_argument(
        '--json',
        type=str,
        default='./VASP_JOBS/prescreening_stability.json',
        help="Path to prescreening_stability.json (default: ./VASP_JOBS/prescreening_stability.json)"
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help="Path to directory containing composition folders with generated_crystals_cif.zip"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='./stable_structures',
        help="Output directory for extracted CIF files (default: ./stable_structures)"
    )
    
    args = parser.parse_args()
    
    return extract_passed_structures(args.json, args.results_dir, args.output_dir)


if __name__ == '__main__':
    sys.exit(main())

