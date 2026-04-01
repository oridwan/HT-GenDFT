#!/usr/bin/env python3
"""
Compare initial and relaxed space groups using stored metadata instead of
recomputing symmetry from structures.

Sources:
- Initial space group: filtered.db key `space_group_number`
- Final space group: electride_analysis.csv column `spacegroup`

Matching key:
- filtered.db `structure_id`
- electride_analysis.csv `formula`

Optional:
- limit comparison to the CIFs present in a specific folder
- include the expected Relax/CONTCAR path under a VASP root

Examples:
  python compare_spacegroup_db_csv.py \
    --cif-dir /scratch/oridwan/SuperConductorFlow/prescreen_new/filter/cif_1e_1_to_5e_2 \
    --initial-db /scratch/oridwan/SuperConductorFlow/prescreen_new/filter/filtered.db \
    --final-csv /scratch/oridwan/SuperConductorFlow/VASP-out-Boron/electride_analysis.csv \
    --vasp-root /scratch/oridwan/SuperConductorFlow/VASP-out-Boron \
    --output /scratch/oridwan/SuperConductorFlow/prescreen_new/filter/cif_1e_1_to_5e_2_spacegroup_compare.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def to_int(value) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except Exception:
        return None


def load_initial_db(db_path: Path) -> Dict[str, Dict[str, object]]:
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT key_value_pairs FROM systems")

    rows: Dict[str, Dict[str, object]] = {}
    for (kvp_text,) in cur:
        if not kvp_text:
            continue
        try:
            kvp = json.loads(kvp_text)
        except Exception:
            continue
        if not isinstance(kvp, dict):
            continue

        structure_id = kvp.get("structure_id")
        if structure_id is None:
            continue
        structure_id = str(structure_id).strip()
        if not structure_id:
            continue

        rows[structure_id] = {
            "structure_id": structure_id,
            "initial_space_group_number": to_int(kvp.get("space_group_number")),
            "composition": kvp.get("composition", ""),
            "e_above_hull": kvp.get("e_above_hull", ""),
            "status": kvp.get("status", ""),
        }

    conn.close()
    return rows


def load_final_csv(csv_path: Path) -> Tuple[Dict[str, Dict[str, object]], Dict[str, int]]:
    rows_by_id: Dict[str, List[Dict[str, object]]] = defaultdict(list)

    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            structure_id = row.get("formula", "")
            if structure_id is None:
                continue
            structure_id = str(structure_id).strip()
            if not structure_id:
                continue
            rows_by_id[structure_id].append(
                {
                    "structure_id": structure_id,
                    "final_space_group_number": to_int(row.get("spacegroup")),
                    "e0025": row.get("e0025", ""),
                    "e05": row.get("e05", ""),
                    "e10": row.get("e10", ""),
                    "band0": row.get("band0", ""),
                    "e_above_hull": row.get("e_above_hull", ""),
                    "is_candidate": row.get("is_candidate", ""),
                }
            )

    chosen: Dict[str, Dict[str, object]] = {}
    duplicate_counts: Dict[str, int] = {}
    for structure_id, entries in rows_by_id.items():
        chosen[structure_id] = entries[0]
        duplicate_counts[structure_id] = len(entries)
    return chosen, duplicate_counts


def load_cif_map(cif_dir: Optional[Path]) -> Dict[str, Path]:
    if cif_dir is None:
        return {}
    return {path.stem: path for path in sorted(cif_dir.glob("*.cif"))}


def build_expected_contcar_path(vasp_root: Optional[Path], structure_id: str) -> str:
    if vasp_root is None:
        return ""
    formula_prefix = structure_id.split("_s", 1)[0]
    return str(vasp_root / formula_prefix / structure_id / "Relax" / "CONTCAR")


def compare_spacegroups(
    initial_db: Path,
    final_csv: Path,
    output: Path,
    cif_dir: Optional[Path],
    vasp_root: Optional[Path],
    limit: Optional[int],
) -> int:
    initial_rows = load_initial_db(initial_db)
    final_rows, duplicate_counts = load_final_csv(final_csv)
    cif_map = load_cif_map(cif_dir)

    if cif_dir is not None:
        target_ids = sorted(cif_map)
    else:
        target_ids = sorted(initial_rows)

    if limit is not None:
        target_ids = target_ids[:limit]

    matched_target_ids = sum(1 for structure_id in target_ids if structure_id in final_rows)

    output_rows = []
    initial_missing = 0
    final_missing = 0
    changed_count = 0
    same_count = 0

    for structure_id in target_ids:
        initial = initial_rows.get(structure_id)
        final = final_rows.get(structure_id)
        cif_path = str(cif_map[structure_id]) if structure_id in cif_map else ""
        contcar_path = build_expected_contcar_path(vasp_root, structure_id)
        contcar_exists = bool(contcar_path and Path(contcar_path).is_file())

        initial_sg = initial["initial_space_group_number"] if initial else None
        final_sg = final["final_space_group_number"] if final else None

        if initial is None:
            initial_missing += 1
        if final is None:
            final_missing += 1

        changed = ""
        if initial_sg is not None and final_sg is not None:
            changed = initial_sg != final_sg
            if changed:
                changed_count += 1
            else:
                same_count += 1

        output_rows.append(
            {
                "structure_id": structure_id,
                "cif_path": cif_path,
                "relax_contcar_path": contcar_path,
                "relax_contcar_exists": contcar_exists,
                "initial_space_group_number": initial_sg if initial_sg is not None else "",
                "final_space_group_number": final_sg if final_sg is not None else "",
                "spacegroup_changed": changed,
                "initial_in_db": initial is not None,
                "final_in_csv": final is not None,
                "final_csv_duplicate_rows": duplicate_counts.get(structure_id, 0),
                "composition": initial["composition"] if initial else "",
                "initial_e_above_hull": initial["e_above_hull"] if initial else "",
                "final_e_above_hull": final["e_above_hull"] if final else "",
                "is_candidate": final["is_candidate"] if final else "",
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "structure_id",
        "cif_path",
        "relax_contcar_path",
        "relax_contcar_exists",
        "initial_space_group_number",
        "final_space_group_number",
        "spacegroup_changed",
        "initial_in_db",
        "final_in_csv",
        "final_csv_duplicate_rows",
        "composition",
        "initial_e_above_hull",
        "final_e_above_hull",
        "is_candidate",
    ]
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    print("=" * 70)
    print("Space Group Comparison From DB/CSV")
    print("=" * 70)
    print(f"Initial DB: {initial_db}")
    print(f"Final CSV: {final_csv}")
    print(f"CIF dir filter: {cif_dir if cif_dir else 'None'}")
    print(f"VASP root: {vasp_root if vasp_root else 'None'}")
    print(f"Rows written: {len(output_rows)}")
    print(f"Target IDs with final CSV matches: {matched_target_ids}")
    print(f"Missing in initial DB: {initial_missing}")
    print(f"Missing in final CSV: {final_missing}")
    print(f"Same space group: {same_count}")
    print(f"Changed space group: {changed_count}")
    print(f"CSV written to: {output}")
    if cif_dir is not None and matched_target_ids == 0:
        print("\nWARNING: The selected CIF folder has zero overlap with the final CSV.")
        print("This usually means the CIF folder and final CSV come from different subsets.")
    return 0


def main() -> int:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Compare initial and relaxed space groups using filtered.db and electride_analysis.csv."
    )
    parser.add_argument(
        "--initial-db",
        type=Path,
        default=script_dir / "filtered.db",
        help="Initial ASE DB containing structure_id and space_group_number.",
    )
    parser.add_argument(
        "--final-csv",
        type=Path,
        required=True,
        help="CSV containing final relaxed space groups, e.g. electride_analysis.csv.",
    )
    parser.add_argument(
        "--cif-dir",
        type=Path,
        default=None,
        help="Optional CIF folder to limit comparison to those structure_ids.",
    )
    parser.add_argument(
        "--vasp-root",
        type=Path,
        default=None,
        help="Optional VASP root for expected Relax/CONTCAR paths.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=script_dir / "spacegroup_compare_db_csv.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit for quick testing.",
    )
    args = parser.parse_args()

    if not args.initial_db.exists():
        print(f"ERROR: Initial DB not found: {args.initial_db}")
        return 1
    if not args.final_csv.exists():
        print(f"ERROR: Final CSV not found: {args.final_csv}")
        return 1
    if args.cif_dir is not None and not args.cif_dir.is_dir():
        print(f"ERROR: CIF dir not found: {args.cif_dir}")
        return 1
    if args.vasp_root is not None and not args.vasp_root.is_dir():
        print(f"ERROR: VASP root not found: {args.vasp_root}")
        return 1

    return compare_spacegroups(
        initial_db=args.initial_db.resolve(),
        final_csv=args.final_csv.resolve(),
        output=args.output.resolve(),
        cif_dir=args.cif_dir.resolve() if args.cif_dir else None,
        vasp_root=args.vasp_root.resolve() if args.vasp_root else None,
        limit=args.limit,
    )


if __name__ == "__main__":
    raise SystemExit(main())
