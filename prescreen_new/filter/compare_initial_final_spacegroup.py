#!/usr/bin/env python3
"""
Compare space groups between initial CIF files and relaxed VASP CONTCAR files.

The script matches structures by structure_id:
  <structure_id>.cif  <->  <vasp_root>/<formula>/<structure_id>/Relax/CONTCAR

It computes symmetry from the actual atomic coordinates for both the CIF and the
relaxed CONTCAR using pymatgen's SpacegroupAnalyzer, then writes a CSV report.

Examples:
  python compare_initial_final_spacegroup.py \
    --cif-dir /scratch/oridwan/SuperConductorFlow/prescreen_new/filter/cif_1e_1_to_5e_2 \
    --vasp-root /scratch/oridwan/SuperConductorFlow/VASP-out-Boron \
    --output /scratch/oridwan/SuperConductorFlow/prescreen_new/filter/cif_1e_1_to_5e_2_spacegroup_compare.csv

  python compare_initial_final_spacegroup.py \
    --cif-dir /scratch/oridwan/SuperConductorFlow/prescreen_new/filter/cif_1e_1_to_5e_2-exc-Rb-Cs \
    --vasp-root /scratch/oridwan/SuperConductorFlow/VASP-out-Boron \
    --limit 20
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def build_contcar_index(vasp_root: Path, relax_subdir: str) -> Tuple[Dict[str, Path], Dict[str, List[Path]]]:
    index: Dict[str, Path] = {}
    duplicates: Dict[str, List[Path]] = {}
    pattern = relax_subdir

    for contcar in vasp_root.rglob(pattern):
        structure_id = contcar.parent.parent.name
        if structure_id in index:
            duplicates.setdefault(structure_id, [index[structure_id]]).append(contcar)
        else:
            index[structure_id] = contcar
    return index, duplicates


def resolve_contcar_path(
    vasp_root: Path,
    formula_prefix: str,
    structure_id: str,
    relax_subdir: str,
    fallback_index: Optional[Dict[str, Path]],
) -> Optional[Path]:
    direct_path = vasp_root / formula_prefix / structure_id / relax_subdir
    if direct_path.is_file():
        return direct_path
    if fallback_index is None:
        return None
    return fallback_index.get(structure_id)


def analyze_spacegroup(path: Path, symprec: float, angle_tolerance: float) -> Tuple[str, int]:
    structure = Structure.from_file(path)
    sga = SpacegroupAnalyzer(
        structure,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    )
    return sga.get_space_group_symbol(), int(sga.get_space_group_number())


def compare_folder(
    cif_dir: Path,
    vasp_root: Path,
    output: Path,
    symprec: float,
    angle_tolerance: float,
    relax_subdir: str,
    limit: Optional[int],
) -> int:
    cif_files = sorted(cif_dir.glob("*.cif"))
    if limit is not None:
        cif_files = cif_files[:limit]

    contcar_index: Optional[Dict[str, Path]] = None
    duplicates: Dict[str, List[Path]] = {}

    rows = []
    changed_count = 0
    matched_count = 0
    missing_contcar_count = 0
    error_count = 0

    for cif_path in cif_files:
        structure_id = cif_path.stem
        formula_prefix = structure_id.split("_s", 1)[0]
        contcar_path = resolve_contcar_path(
            vasp_root=vasp_root,
            formula_prefix=formula_prefix,
            structure_id=structure_id,
            relax_subdir=relax_subdir,
            fallback_index=contcar_index,
        )
        if contcar_path is None and contcar_index is None:
            contcar_index, duplicates = build_contcar_index(vasp_root, relax_subdir)
            contcar_path = resolve_contcar_path(
                vasp_root=vasp_root,
                formula_prefix=formula_prefix,
                structure_id=structure_id,
                relax_subdir=relax_subdir,
                fallback_index=contcar_index,
            )

        row = {
            "structure_id": structure_id,
            "formula": formula_prefix,
            "initial_cif": str(cif_path),
            "final_contcar": str(contcar_path) if contcar_path else "",
            "initial_sg_symbol": "",
            "initial_sg_number": "",
            "final_sg_symbol": "",
            "final_sg_number": "",
            "spacegroup_changed": "",
            "duplicate_contcar_matches": len(duplicates.get(structure_id, [])),
            "error": "",
        }

        if contcar_path is None:
            row["error"] = "Missing Relax/CONTCAR"
            missing_contcar_count += 1
            rows.append(row)
            continue

        try:
            initial_symbol, initial_number = analyze_spacegroup(
                cif_path,
                symprec=symprec,
                angle_tolerance=angle_tolerance,
            )
            final_symbol, final_number = analyze_spacegroup(
                contcar_path,
                symprec=symprec,
                angle_tolerance=angle_tolerance,
            )
        except Exception as exc:
            row["error"] = str(exc)
            error_count += 1
            rows.append(row)
            continue

        changed = (initial_symbol != final_symbol) or (initial_number != final_number)

        row["initial_sg_symbol"] = initial_symbol
        row["initial_sg_number"] = initial_number
        row["final_sg_symbol"] = final_symbol
        row["final_sg_number"] = final_number
        row["spacegroup_changed"] = changed

        rows.append(row)
        matched_count += 1
        if changed:
            changed_count += 1

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "structure_id",
        "formula",
        "initial_cif",
        "final_contcar",
        "initial_sg_symbol",
        "initial_sg_number",
        "final_sg_symbol",
        "final_sg_number",
        "spacegroup_changed",
        "duplicate_contcar_matches",
        "error",
    ]
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("=" * 70)
    print("Space Group Comparison")
    print("=" * 70)
    print(f"CIF dir: {cif_dir}")
    print(f"VASP root: {vasp_root}")
    print(f"Relax path matched: {relax_subdir}")
    print(f"symprec: {symprec}")
    print(f"angle_tolerance: {angle_tolerance}")
    print(f"CIF files scanned: {len(cif_files)}")
    print(f"Matched CONTCARs: {matched_count}")
    print(f"Missing CONTCARs: {missing_contcar_count}")
    print(f"Comparison errors: {error_count}")
    print(f"Space group changed: {changed_count}")
    print(f"Duplicate CONTCAR structure_ids in VASP root: {len(duplicates)}")
    print(f"CSV written to: {output}")

    if duplicates:
        print("\nDuplicate CONTCAR matches found for these structure_ids:")
        for structure_id in sorted(duplicates)[:25]:
            print(f"  - {structure_id}: {len(duplicates[structure_id])} matches")
        if len(duplicates) > 25:
            print(f"  ... and {len(duplicates) - 25} more")

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare initial CIF space groups with relaxed VASP CONTCAR space groups."
    )
    parser.add_argument(
        "--cif-dir",
        type=Path,
        required=True,
        help="Folder containing initial CIF files.",
    )
    parser.add_argument(
        "--vasp-root",
        type=Path,
        required=True,
        help="VASP output root, e.g. VASP-out-Boron.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("spacegroup_comparison.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--relax-subdir",
        default="Relax/CONTCAR",
        help="Relative path under each structure directory used as the final structure.",
    )
    parser.add_argument(
        "--symprec",
        type=float,
        default=0.1,
        help="Symmetry tolerance for SpacegroupAnalyzer.",
    )
    parser.add_argument(
        "--angle-tolerance",
        type=float,
        default=5.0,
        help="Angle tolerance for SpacegroupAnalyzer.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit for quick testing.",
    )
    args = parser.parse_args()

    if not args.cif_dir.is_dir():
        print(f"ERROR: CIF directory not found: {args.cif_dir}")
        return 1
    if not args.vasp_root.is_dir():
        print(f"ERROR: VASP root not found: {args.vasp_root}")
        return 1

    return compare_folder(
        cif_dir=args.cif_dir.resolve(),
        vasp_root=args.vasp_root.resolve(),
        output=args.output.resolve(),
        symprec=args.symprec,
        angle_tolerance=args.angle_tolerance,
        relax_subdir=args.relax_subdir,
        limit=args.limit,
    )


if __name__ == "__main__":
    raise SystemExit(main())
