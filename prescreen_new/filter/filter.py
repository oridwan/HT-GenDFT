"""
Filter prescreening database by energy-above-hull and remove duplicates.

Notes:
- Uses e_above_hull < threshold (default: 0.1 eV/atom).
- Removes duplicates by reduced formula first, then by structure matching.
- Operates on ASE databases (e.g., prescreening_structures.db).

Example:
  python3 filter.py \
    --db ./full.db \
    --out-db ./filtered.db \
    --hull-threshold 0.1
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from ase.db import connect
from pymatgen.core import Composition
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher


def _row_to_structure(row) -> Structure:
    atoms = row.toatoms()
    return AseAtomsAdaptor.get_structure(atoms)


def _reduced_formula(row) -> str:
    # Prefer explicit composition if present, otherwise formula.
    comp = row.get("composition", None)
    if comp:
        return Composition(comp).reduced_formula
    return Composition(row.formula).reduced_formula


def filter_rows(
    db_path: Path,
    hull_threshold: float,
    require_passed: bool,
    print_every: int,
) -> Dict[str, List[Tuple[object, Structure]]]:
    by_formula: Dict[str, List[Tuple[object, Structure]]] = defaultdict(list)
    total = 0
    kept = 0
    skipped_missing = 0
    skipped_failed = 0

    with connect(db_path) as db:
        for row in db.select():
            total += 1
            if print_every > 0 and total % print_every == 0:
                print(
                    f"Scanned {total} rows | kept={kept} | missing_ehull={skipped_missing} | "
                    f"failed_prescreen={skipped_failed}"
                )
            e_hull = row.get("e_above_hull", None)
            if e_hull is None:
                skipped_missing += 1
                continue
            if e_hull >= hull_threshold:
                continue
            if require_passed:
                passed = row.get("passed_prescreening", None)
                if passed is False:
                    skipped_failed += 1
                    continue
            try:
                struct = _row_to_structure(row)
            except Exception:
                # Skip rows that cannot be converted to Structure
                continue
            rf = _reduced_formula(row)
            by_formula[rf].append((row, struct))
            kept += 1

    print(f"Total rows: {total}")
    print(f"Kept by e_above_hull < {hull_threshold}: {kept}")
    print(f"Skipped (missing e_above_hull): {skipped_missing}")
    if require_passed:
        print(f"Skipped (failed prescreening): {skipped_failed}")

    return by_formula


def dedupe_by_structure(
    by_formula: Dict[str, List[Tuple[object, Structure]]],
    matcher: StructureMatcher,
) -> List[object]:
    unique_rows: List[object] = []
    total_groups = 0
    total_structs = 0

    for idx, (rf, items) in enumerate(by_formula.items(), start=1):
        rows, structs = zip(*items)
        total_structs += len(structs)
        groups = matcher.group_structures(list(structs))
        total_groups += len(groups)
        if idx % 25 == 0:
            print(
                f"Dedupe progress: {idx}/{len(by_formula)} formulas | "
                f"current={rf} | structs={len(structs)} | groups={len(groups)}"
            )

        # Keep the first structure in each group
        for group in groups:
            first_struct = group[0]
            idx = structs.index(first_struct)
            unique_rows.append(rows[idx])

    print(f"Structures before dedupe: {total_structs}")
    print(f"Unique structure groups: {total_groups}")
    print(f"Unique structures kept: {len(unique_rows)}")
    return unique_rows


def write_filtered_db(rows: List[object], out_db: Path) -> None:
    out_db.parent.mkdir(parents=True, exist_ok=True)
    with connect(out_db) as db:
        for row in rows:
            atoms = row.toatoms()
            kv = dict(row.key_value_pairs)
            data = dict(row.data) if row.data else {}
            db.write(atoms, **kv, data=data)


def write_summary_json(rows: List[object], out_json: Path) -> None:
    out_json.parent.mkdir(parents=True, exist_ok=True)
    payload = []
    for row in rows:
        payload.append(
            {
                "id": row.id,
                "formula": row.formula,
                "composition": row.get("composition", None),
                "e_above_hull": row.get("e_above_hull", None),
                "passed_prescreening": row.get("passed_prescreening", None),
                "space_group_number": row.get("space_group_number", None),
            }
        )
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=2)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Filter prescreen database by e_above_hull and remove duplicates."
    )
    parser.add_argument(
        "--db",
        type=str,
        default="./full.db",
        help="Path to input ASE DB (default: ./full.db)",
    )
    parser.add_argument(
        "--out-db",
        type=str,
        default="./filtered.db",
        help="Path to output ASE DB (default: ./filtered.db)",
    )
    parser.add_argument(
        "--out-json",
        type=str,
        default="",
        help="Optional JSON summary output (default: disabled)",
    )
    parser.add_argument(
        "--hull-threshold",
        type=float,
        default=0.1,
        help="Energy above hull threshold (eV/atom)",
    )
    parser.add_argument(
        "--require-passed",
        action="store_true",
        help="Require passed_prescreening=True when present",
    )
    parser.add_argument(
        "--print-every",
        type=int,
        default=10000,
        help="Progress print frequency during DB scan (rows). Use 0 to disable.",
    )
    parser.add_argument(
        "--matcher-stol",
        type=float,
        default=0.2,
        help="StructureMatcher site tolerance",
    )
    parser.add_argument(
        "--matcher-ltol",
        type=float,
        default=0.2,
        help="StructureMatcher lattice tolerance",
    )
    parser.add_argument(
        "--matcher-angle-tol",
        type=float,
        default=5.0,
        help="StructureMatcher angle tolerance (degrees)",
    )

    args = parser.parse_args()
    db_path = Path(args.db).expanduser()
    out_db = Path(args.out_db).expanduser()
    out_json = Path(args.out_json).expanduser() if args.out_json else None

    if not db_path.exists():
        print(f"ERROR: Input DB not found: {db_path}")
        return 1

    by_formula = filter_rows(
        db_path=db_path,
        hull_threshold=args.hull_threshold,
        require_passed=args.require_passed,
        print_every=args.print_every,
    )

    matcher = StructureMatcher(
        stol=args.matcher_stol,
        ltol=args.matcher_ltol,
        angle_tol=args.matcher_angle_tol,
        primitive_cell=False,
        scale=True,
        attempt_supercell=False,
    )

    unique_rows = dedupe_by_structure(by_formula, matcher)
    write_filtered_db(unique_rows, out_db)
    print(f"Filtered DB written to: {out_db}")

    if out_json is not None:
        write_summary_json(unique_rows, out_json)
        print(f"Summary JSON written to: {out_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
