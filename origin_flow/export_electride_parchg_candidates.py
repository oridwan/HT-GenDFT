#!/usr/bin/env python3
"""Export PARCHG files for electride candidates and write a filtered analysis CSV.

Criteria (from ASE DB):
  - is_electride == True
  - space_group_number > 15

Then:
  1) Find each structure's SPE directory from workflow.json
  2) Extract only PARCHG-* files from PARCHG.tar.gz into an export folder
  3) Write a filtered electride_analysis.csv containing only matching structures

Example:
  python origin_flow/export_electride_parchg_candidates.py \
    --db VASP-out/electride_data.db \
    --workflow VASP-out/workflow.json \
    --csv VASP-out/electride_analysis.csv \
    --out-dir VASP-out/electride_parchg_candidates \
    --out-csv VASP-out/electride_analysis_electride_spggt15.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import tarfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    from ase.db import connect
except Exception as exc:  # pragma: no cover - runtime dependency guard
    raise SystemExit(
        "ASE is required for this script. Install/activate an environment with ase. "
        f"Import error: {exc}"
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Export electride PARCHG files and filtered CSV")
    p.add_argument("--db", required=True, help="Path to ASE database (electride_data.db)")
    p.add_argument("--workflow", default=None, help="Path to workflow.json (default: <db_dir>/workflow.json)")
    p.add_argument("--csv", dest="csv_path", default=None, help="Path to electride_analysis.csv (default: <db_dir>/electride_analysis.csv)")
    p.add_argument("--out-dir", default=None, help="Destination directory for extracted PARCHG files (default: <db_dir>/electride_parchg_candidates)")
    p.add_argument("--out-csv", default=None, help="Filtered CSV output path (default: <db_dir>/electride_analysis_electride_spggt15.csv)")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing exported PARCHG files")
    p.add_argument("--limit", type=int, default=None, help="Only process first N electride structures (debug/testing)")
    return p.parse_args()


def default_paths(args: argparse.Namespace) -> None:
    db_path = Path(args.db).resolve()
    db_dir = db_path.parent
    if args.workflow is None:
        args.workflow = str(db_dir / "workflow.json")
    if args.csv_path is None:
        args.csv_path = str(db_dir / "electride_analysis.csv")
    if args.out_dir is None:
        args.out_dir = str(db_dir / "electride_parchg_candidates")
    if args.out_csv is None:
        args.out_csv = str(db_dir / "electride_analysis_electride_spggt15.csv")


def _row_get(row, key: str, default=None):
    try:
        if hasattr(row, key):
            value = getattr(row, key)
            if value is not None:
                return value
    except Exception:
        pass
    kvp = getattr(row, "key_value_pairs", None)
    if isinstance(kvp, dict):
        return kvp.get(key, default)
    return default


def _to_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    s = str(value).strip().lower()
    return s in {"1", "true", "yes", "y", "t"}


def _to_int(value) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(float(value))
    except Exception:
        return None


def load_electride_candidates(db_path: Path, limit: Optional[int] = None) -> Dict[str, dict]:
    db = connect(str(db_path))
    selected: Dict[str, dict] = {}
    duplicate_count = 0
    total_rows = 0

    for row in db.select():
        total_rows += 1
        struct_id = _row_get(row, "structure_id")
        if not struct_id:
            continue
        is_electride = _to_bool(_row_get(row, "is_electride", False))
        sg = _to_int(_row_get(row, "space_group_number"))
        if not is_electride or sg is None or sg <= 15:
            continue

        # Deduplicate by structure_id (can happen after partial reruns).
        if struct_id in selected:
            duplicate_count += 1
            continue

        selected[struct_id] = {
            "structure_id": struct_id,
            "composition": _row_get(row, "composition"),
            "e_above_hull": _row_get(row, "e_above_hull"),
            "space_group_number": sg,
            "db_row_id": getattr(row, "id", None),
        }
        if limit is not None and len(selected) >= limit:
            break

    print(f"ASE DB rows scanned: {total_rows}")
    print(f"Electride candidates from DB (is_electride=True, space_group_number>15): {len(selected)}")
    if duplicate_count:
        print(f"Duplicate structure_id rows skipped: {duplicate_count}")
    print("")
    return selected


def load_workflow_map(workflow_path: Path) -> Dict[str, dict]:
    with workflow_path.open() as f:
        data = json.load(f)
    return data.get("structures", {})


def _safe_parchg_member_name(member: tarfile.TarInfo) -> Optional[str]:
    name = Path(member.name).name
    if not name.startswith("PARCHG-"):
        return None
    if name in {"", ".", ".."}:
        return None
    return name


def export_parchg_for_structure(struct_id: str, spe_dir: Path, out_root: Path, overwrite: bool = False) -> Tuple[bool, str, int]:
    """Return (success, reason, file_count)."""
    if not spe_dir.exists():
        return False, f"SPE directory not found: {spe_dir}", 0

    dest_dir = out_root / struct_id
    dest_dir.mkdir(parents=True, exist_ok=True)

    tar_path = spe_dir / "PARCHG.tar.gz"
    if tar_path.exists():
        extracted = 0
        try:
            with tarfile.open(tar_path, "r:gz") as tf:
                for member in tf.getmembers():
                    if not member.isfile():
                        continue
                    out_name = _safe_parchg_member_name(member)
                    if out_name is None:
                        continue
                    out_path = dest_dir / out_name
                    if out_path.exists() and not overwrite:
                        extracted += 1
                        continue
                    src_f = tf.extractfile(member)
                    if src_f is None:
                        continue
                    with src_f, out_path.open("wb") as dst_f:
                        shutil.copyfileobj(src_f, dst_f)
                    extracted += 1
        except Exception as exc:
            return False, f"Failed to extract {tar_path.name}: {exc}", 0

        if extracted == 0:
            return False, f"No PARCHG-* members found in {tar_path.name}", 0
        return True, "ok", extracted

    # Fallback: copy already-extracted files if tar.gz is absent.
    src_files = sorted([p for p in spe_dir.glob("PARCHG-*") if p.is_file()])
    if not src_files:
        return False, "No PARCHG.tar.gz and no PARCHG-* files", 0

    copied = 0
    for src in src_files:
        dst = dest_dir / src.name
        if dst.exists() and not overwrite:
            copied += 1
            continue
        shutil.copy2(src, dst)
        copied += 1
    return True, "ok", copied


def export_parchg_files(candidates: Dict[str, dict], workflow_map: Dict[str, dict], out_dir: Path, overwrite: bool = False) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    ok_count = 0
    fail_count = 0
    missing_workflow = 0
    total_files = 0

    for idx, struct_id in enumerate(candidates, start=1):
        sdata = workflow_map.get(struct_id)
        if not sdata:
            print(f"[{idx}/{len(candidates)}] {struct_id}: MISSING in workflow.json")
            missing_workflow += 1
            fail_count += 1
            continue
        spe_dir = Path(sdata.get("spe_dir") or "")
        success, reason, nfiles = export_parchg_for_structure(struct_id, spe_dir, out_dir, overwrite=overwrite)
        if success:
            ok_count += 1
            total_files += nfiles
            print(f"[{idx}/{len(candidates)}] {struct_id}: exported {nfiles} PARCHG files")
        else:
            fail_count += 1
            print(f"[{idx}/{len(candidates)}] {struct_id}: WARNING: {reason}")

    print("")
    print("PARCHG export summary")
    print("--------------------")
    print(f"Candidates requested: {len(candidates)}")
    print(f"Successful exports:   {ok_count}")
    print(f"Failed exports:       {fail_count}")
    print(f"Missing in workflow:  {missing_workflow}")
    print(f"Total PARCHG files:   {total_files}")
    print(f"Output directory:     {out_dir}")
    print("")


def _float(row: dict, key: str) -> float:
    try:
        return float(row.get(key, "0") or 0)
    except Exception:
        return 0.0


def _csv_electride_and_spg_ok(row: dict) -> bool:
    sg = _to_int(row.get("spacegroup"))
    if sg is None or sg <= 15:
        return False
    return ((_float(row, "e0025") > 0 or _float(row, "e05") > 0) and
            (_float(row, "e10") > 0 or _float(row, "band0") > 0))


def write_filtered_csv(csv_in: Path, csv_out: Path, candidate_ids: Iterable[str]) -> None:
    candidate_set = set(candidate_ids)
    if not csv_in.exists():
        print(f"WARNING: CSV not found, skipping filtered CSV write: {csv_in}")
        return

    with csv_in.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    if not fieldnames:
        raise SystemExit(f"CSV has no header: {csv_in}")

    id_col = "formula" if "formula" in fieldnames else ("structure_id" if "structure_id" in fieldnames else None)
    if id_col is None:
        raise SystemExit(f"Could not find structure-id column (expected 'formula' or 'structure_id') in {csv_in}")

    filtered = [
        row for row in rows
        if row.get(id_col) in candidate_set and _csv_electride_and_spg_ok(row)
    ]

    csv_out.parent.mkdir(parents=True, exist_ok=True)
    with csv_out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(filtered)

    print("Filtered CSV summary")
    print("--------------------")
    print(f"Input rows:          {len(rows)}")
    print(f"Filtered rows:       {len(filtered)}")
    print(f"Output CSV:          {csv_out}")
    print("")


def main() -> None:
    args = parse_args()
    default_paths(args)

    db_path = Path(args.db).resolve()
    workflow_path = Path(args.workflow).resolve()
    csv_path = Path(args.csv_path).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_csv = Path(args.out_csv).resolve()

    print("=" * 72)
    print("Export Electride PARCHG Candidates")
    print("=" * 72)
    print(f"DB:       {db_path}")
    print(f"Workflow: {workflow_path}")
    print(f"CSV:      {csv_path}")
    print(f"Out dir:  {out_dir}")
    print(f"Out CSV:  {out_csv}")
    if args.limit is not None:
        print(f"Limit:    {args.limit}")
    print("Criteria: is_electride=True and space_group_number>15")
    print("")

    candidates = load_electride_candidates(db_path, limit=args.limit)
    if not candidates:
        print("No matching candidates found in DB. Nothing to do.")
        return

    if not workflow_path.exists():
        raise SystemExit(f"workflow.json not found: {workflow_path}")
    workflow_map = load_workflow_map(workflow_path)

    export_parchg_files(candidates, workflow_map, out_dir, overwrite=args.overwrite)
    write_filtered_csv(csv_path, out_csv, candidates.keys())

    print("Done.")


if __name__ == "__main__":
    main()
