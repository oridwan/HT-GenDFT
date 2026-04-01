#!/usr/bin/env python3
"""
Backfill missing CIF exports from filtered.db into the expected folders.

Default behavior matches the current filter layout in this directory:
- scans filtered.db
- treats the following folders as already-exported structure sources:
  cif_less_than_1e_3
  cif_1e_3_to_5e_2
  cif_1e_1_to_5e_2
  cif_1e_1_to_5e_2-exc-Rb-Cs
  Li_Na
- targets the 5e-2 < e_above_hull < 1e-1 window
- sends structures containing Rb or Cs to cif_1e_1_to_5e_2
- sends structures without Rb/Cs to cif_1e_1_to_5e_2-exc-Rb-Cs

Examples:
  python backfill_missing_cifs.py --dry-run
  python backfill_missing_cifs.py
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path
from typing import Optional, Sequence, Set, Tuple

try:
    from ase.io import write as ase_write
except ImportError as exc:
    print(f"ERROR: ASE is required: {exc}")
    sys.exit(1)

try:
    from pyxtal.db import database_topology
except ImportError as exc:
    print(f"ERROR: PyXtal is required: {exc}")
    sys.exit(1)


DEFAULT_EXISTING_DIRS = [
    "cif_less_than_1e_3",
    "cif_1e_3_to_5e_2",
    "cif_1e_1_to_5e_2",
    "cif_1e_1_to_5e_2-exc-Rb-Cs",
    "Li_Na",
]


def _parse_structure_id(row) -> Optional[str]:
    structure_id = getattr(row, "structure_id", None)
    if structure_id:
        return structure_id

    if hasattr(row, "get"):
        try:
            structure_id = row.get("structure_id", None)
        except Exception:
            structure_id = None
        if structure_id:
            return structure_id

    kvp = getattr(row, "key_value_pairs", None)
    if isinstance(kvp, dict):
        return kvp.get("structure_id")
    if isinstance(kvp, str):
        try:
            parsed = json.loads(kvp)
        except Exception:
            return None
        if isinstance(parsed, dict):
            return parsed.get("structure_id")
    return None


def _parse_formula_elements(formula_prefix: str) -> Set[str]:
    elements: Set[str] = set()
    token = ""
    for ch in formula_prefix:
        if ch.isupper():
            if token:
                elements.add(token)
            token = ch
        elif ch.islower():
            token += ch
        elif token:
            elements.add(token)
            token = ""
    if token:
        elements.add(token)
    return elements


def _load_present_ids(base_dir: Path, existing_dirs: Sequence[str]) -> Set[str]:
    present = set()
    for name in existing_dirs:
        folder = base_dir / name
        if not folder.is_dir():
            continue
        present.update(path.stem for path in folder.glob("*.cif"))
    return present


def _in_window(value: float, min_value: float, max_value: float) -> bool:
    return min_value < value < max_value


def _export_row(db, row, row_id: int, cif_path: Path) -> Tuple[bool, str]:
    try:
        xtal = db.get_pyxtal(row_id)
        if xtal is not None:
            xtal.to_file(str(cif_path), fmt="cif")
            return True, ""
    except Exception as exc:
        pyxtal_err = f"pyxtal export failed: {exc}"
    else:
        pyxtal_err = "pyxtal export returned None"

    try:
        atoms = row.toatoms()
        ase_write(str(cif_path), atoms, format="cif")
        return True, ""
    except Exception as exc:
        return False, f"{pyxtal_err}; ASE export failed: {exc}"


def _collect_targets(
    db_path: Path,
    present_ids: Set[str],
    excluded_elements: Sequence[str],
    min_e_above_hull: float,
    max_e_above_hull: float,
) -> Tuple[int, list[Tuple[int, str, bool]]]:
    excluded_set = set(excluded_elements)
    total_rows = 0
    targets: list[Tuple[int, str, bool]] = []

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT id, key_value_pairs FROM systems")
    for row_id, kvp_text in cur:
        total_rows += 1
        if not kvp_text:
            continue
        try:
            kvp = json.loads(kvp_text)
        except Exception:
            continue
        if not isinstance(kvp, dict):
            continue

        structure_id = kvp.get("structure_id")
        if structure_id is not None:
            structure_id = str(structure_id).strip()
        if not structure_id or structure_id in present_ids:
            continue

        try:
            e_hull = float(kvp.get("e_above_hull", None))
        except Exception:
            continue
        if not _in_window(e_hull, min_e_above_hull, max_e_above_hull):
            continue

        formula_prefix = structure_id.split("_s", 1)[0]
        has_excluded = bool(_parse_formula_elements(formula_prefix) & excluded_set)
        targets.append((row_id, structure_id, has_excluded))

    conn.close()
    return total_rows, targets


def run_backfill(
    db_path: Path,
    base_dir: Path,
    existing_dirs: Sequence[str],
    regular_dir_name: str,
    excluded_dir_name: str,
    excluded_elements: Sequence[str],
    min_e_above_hull: float,
    max_e_above_hull: float,
    dry_run: bool,
    print_every: int,
) -> int:
    present_ids = _load_present_ids(base_dir, existing_dirs)
    regular_dir = base_dir / regular_dir_name
    excluded_dir = base_dir / excluded_dir_name
    regular_dir.mkdir(parents=True, exist_ok=True)
    excluded_dir.mkdir(parents=True, exist_ok=True)

    total_rows, targets = _collect_targets(
        db_path=db_path,
        present_ids=present_ids,
        excluded_elements=excluded_elements,
        min_e_above_hull=min_e_above_hull,
        max_e_above_hull=max_e_above_hull,
    )

    target_rows = len(targets)
    target_regular = sum(1 for _, _, has_excluded in targets if has_excluded)
    target_excluded = target_rows - target_regular
    exported_regular = 0
    exported_excluded = 0
    skipped_existing = 0
    failed = []
    db = None
    if not dry_run:
        db = database_topology(str(db_path))

    for row_id, structure_id, has_excluded in targets:
        outdir = regular_dir if has_excluded else excluded_dir
        cif_path = outdir / f"{structure_id}.cif"

        if cif_path.exists():
            skipped_existing += 1
            continue

        if dry_run:
            continue

        row = db.db.get(id=row_id)
        ok, err = _export_row(db, row, row_id, cif_path)
        if ok:
            if has_excluded:
                exported_regular += 1
            else:
                exported_excluded += 1
            exported_total = exported_regular + exported_excluded
            if print_every > 0 and exported_total % print_every == 0:
                print(
                    f"Exported {exported_total}/{target_rows} targets so far "
                    f"(regular={exported_regular}, exc={exported_excluded})"
                )
        else:
            failed.append((structure_id, err))

    print("=" * 70)
    print("Backfill Summary")
    print("=" * 70)
    print(f"Database: {db_path}")
    print(f"Base dir: {base_dir}")
    print(f"Existing dirs scanned: {', '.join(existing_dirs)}")
    print(f"Window: {min_e_above_hull} < e_above_hull < {max_e_above_hull}")
    print(f"Excluded elements routed to regular dir: {', '.join(excluded_elements)}")
    print(f"Regular dir: {regular_dir}")
    print(f"Excluded dir: {excluded_dir}")
    print(f"Total DB rows scanned: {total_rows}")
    print(f"Target missing rows: {target_rows}")
    print(f"Target rows to regular dir: {target_regular}")
    print(f"Target rows to excluded dir: {target_excluded}")
    print(f"Skipped because file already exists: {skipped_existing}")
    if dry_run:
        print("Dry run: no files were written")
    else:
        print(f"Exported to regular dir: {exported_regular}")
        print(f"Exported to excluded dir: {exported_excluded}")
        print(f"Failed exports: {len(failed)}")
        print(f"Final regular dir count: {sum(1 for _ in regular_dir.glob('*.cif'))}")
        print(f"Final excluded dir count: {sum(1 for _ in excluded_dir.glob('*.cif'))}")
        if failed:
            print("\nFailed exports:")
            for structure_id, err in failed[:50]:
                print(f"  - {structure_id}: {err}")
            if len(failed) > 50:
                print(f"  ... and {len(failed) - 50} more")

    return 0 if dry_run or not failed else 2


def main() -> int:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Backfill missing CIFs from filtered.db into the expected folders."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=script_dir / "filtered.db",
        help="Path to filtered ASE/PyXtal DB (default: ./filtered.db relative to this script).",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=script_dir,
        help="Directory containing the CIF folders (default: this script directory).",
    )
    parser.add_argument(
        "--existing-dirs",
        nargs="+",
        default=DEFAULT_EXISTING_DIRS,
        help="Folders scanned to decide which structure_ids are already exported.",
    )
    parser.add_argument(
        "--regular-dir",
        default="cif_1e_1_to_5e_2",
        help="Destination for structures containing excluded elements.",
    )
    parser.add_argument(
        "--excluded-dir",
        default="cif_1e_1_to_5e_2-exc-Rb-Cs",
        help="Destination for structures without the excluded elements.",
    )
    parser.add_argument(
        "--exclude-elements",
        nargs="+",
        default=["Rb", "Cs"],
        help="Elements that route a structure into the regular destination folder.",
    )
    parser.add_argument(
        "--min-e-above-hull",
        type=float,
        default=5e-2,
        help="Strict lower bound for e_above_hull (default: 5e-2).",
    )
    parser.add_argument(
        "--max-e-above-hull",
        type=float,
        default=1e-1,
        help="Strict upper bound for e_above_hull (default: 1e-1).",
    )
    parser.add_argument(
        "--print-every",
        type=int,
        default=250,
        help="Progress print frequency during export. Use 0 to disable.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report missing targets without writing CIF files.",
    )
    args = parser.parse_args()

    if not args.db.exists():
        print(f"ERROR: Database not found: {args.db}")
        return 1

    return run_backfill(
        db_path=args.db.resolve(),
        base_dir=args.base_dir.resolve(),
        existing_dirs=args.existing_dirs,
        regular_dir_name=args.regular_dir,
        excluded_dir_name=args.excluded_dir,
        excluded_elements=args.exclude_elements,
        min_e_above_hull=args.min_e_above_hull,
        max_e_above_hull=args.max_e_above_hull,
        dry_run=args.dry_run,
        print_every=args.print_every,
    )


if __name__ == "__main__":
    raise SystemExit(main())
