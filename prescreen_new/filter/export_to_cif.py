#!/usr/bin/env python3
"""
Convert structures from a PyXtal database to CIF files.

Examples:
  python3 export_to_cif.py --db filtered.db --output-dir cif_3e_2_to_1e_2 --overwrite  --min-e-above-hull 0.01 --max-e-above-hull 0.03

"""

import argparse
import sys
import json
from pathlib import Path

try:
    from pyxtal.db import database_topology
except ImportError as exc:
    print(f"ERROR: PyXtal is required: {exc}")
    print("Install with: conda install -c conda-forge pyxtal")
    sys.exit(1)

try:
    from ase.io import write as ase_write
except ImportError:
    ase_write = None

REQUIRED_METADATA_KEYS = {
    "chemsys",
    "composition",
    "density",
    "dof",
    "e_above_hull",
    "e_mattersim",
    "passed_prescreening",
    "pearson_symbol",
    "space_group_number",
    "status",
    "structure_id",
    "symmetrized",
    "wps",
}


def _parse_bool(value):
    """Convert common serialized truthy values to bool."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value == 1
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "t", "yes", "y"}
    return False


def _extract_metadata(row):
    """Best-effort extraction of metadata keys from a DB row."""
    metadata = {}

    kvp = getattr(row, "key_value_pairs", None)
    if isinstance(kvp, dict):
        metadata.update(kvp)
    elif isinstance(kvp, str):
        try:
            kvp_dict = json.loads(kvp)
            if isinstance(kvp_dict, dict):
                metadata.update(kvp_dict)
        except Exception:
            pass

    if hasattr(row, "get"):
        for key in REQUIRED_METADATA_KEYS:
            if key in metadata:
                continue
            try:
                val = row.get(key, None)
            except Exception:
                val = None
            if val is not None:
                metadata[key] = val

    for key in REQUIRED_METADATA_KEYS:
        if key in metadata:
            continue
        val = getattr(row, key, None)
        if val is not None:
            metadata[key] = val

    return metadata


def export_db_to_cif(
    db_path,
    output_dir,
    structure_ids=None,
    overwrite=False,
    min_e_above_hull=None,
    max_e_above_hull=None,
):
    """
    Export selected (or all) structures from a PyXtal DB to CIF files.

    Args:
        db_path: Path to .db file
        output_dir: Directory for CIF files
        structure_ids: Optional set of structure_id strings to export
        overwrite: Overwrite existing CIF files if True

    Returns:
        dict with export summary
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    db = database_topology(str(db_path))

    requested = set(structure_ids) if structure_ids else None
    found_ids = set()

    exported = 0
    skipped_existing = 0
    skipped_not_prescreened = 0
    skipped_ehull = 0
    skipped_missing_required_keys = 0
    failed = []

    for row in db.db.select():
        row_id = getattr(row, "id", None)
        structure_id = getattr(row, "structure_id", None)
        if not structure_id and hasattr(row, "get"):
            try:
                structure_id = row.get("structure_id", None)
            except Exception:
                pass
        if not structure_id:
            kvp = getattr(row, "key_value_pairs", None)
            if isinstance(kvp, dict):
                structure_id = kvp.get("structure_id")
            elif isinstance(kvp, str):
                try:
                    kvp_dict = json.loads(kvp)
                    structure_id = kvp_dict.get("structure_id")
                except Exception:
                    pass
        if not structure_id:
            structure_id = getattr(row, "unique_id", None)
        if not structure_id:
            structure_id = f"row_{row_id}"

        metadata = _extract_metadata(row)
        if metadata.get("structure_id") is None:
            metadata["structure_id"] = structure_id

        if requested is not None and structure_id not in requested:
            continue

        found_ids.add(structure_id)

        missing_keys = REQUIRED_METADATA_KEYS - set(metadata.keys())
        if missing_keys:
            skipped_missing_required_keys += 1
            continue

        if not _parse_bool(metadata.get("passed_prescreening")):
            skipped_not_prescreened += 1
            continue

        if min_e_above_hull is not None or max_e_above_hull is not None:
            e_hull = metadata.get("e_above_hull", None)
            if e_hull is None:
                skipped_ehull += 1
                continue
            try:
                e_hull = float(e_hull)
            except Exception:
                skipped_ehull += 1
                continue
            # Inclusive bounds: min <= e_hull <= max
            if min_e_above_hull is not None and e_hull < min_e_above_hull:
                skipped_ehull += 1
                continue
            if max_e_above_hull is not None and e_hull > max_e_above_hull:
                skipped_ehull += 1
                continue
        cif_path = output_dir / f"{structure_id}.cif"

        if cif_path.exists() and not overwrite:
            skipped_existing += 1
            continue

        exported_ok = False
        error_msgs = []

        # Preferred path: retrieve pyxtal object and write directly to CIF.
        if row_id is not None:
            try:
                xtal = db.get_pyxtal(row_id)
                if xtal is not None:
                    xtal.to_file(str(cif_path), fmt="cif")
                    exported_ok = True
            except Exception as exc:
                error_msgs.append(f"pyxtal export failed: {exc}")

        # Fallback path: write CIF from ASE atoms row.
        if not exported_ok:
            if ase_write is None:
                error_msgs.append("ASE not available for fallback export")
            else:
                try:
                    atoms = row.toatoms()
                    ase_write(str(cif_path), atoms, format="cif")
                    exported_ok = True
                except Exception as exc:
                    error_msgs.append(f"ASE export failed: {exc}")

        if exported_ok:
            exported += 1
        else:
            failed.append((structure_id, "; ".join(error_msgs)))

    missing = []
    if requested is not None:
        missing = sorted(requested - found_ids)

    return {
        "exported": exported,
        "skipped_existing": skipped_existing,
        "skipped_not_prescreened": skipped_not_prescreened,
        "skipped_e_above_hull": skipped_ehull,
        "skipped_missing_required_keys": skipped_missing_required_keys,
        "failed": failed,
        "missing": missing,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Export structures from a PyXtal database to CIF files."
    )
    parser.add_argument(
        "--db",
        type=Path,
        required=True,
        help="Path to PyXtal database (.db).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("./cif_from_db"),
        help="Directory to write CIF files (default: ./cif_from_db).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing CIF files.",
    )
    parser.add_argument(
        "--max-e-above-hull",
        type=float,
        default=None,
        help="Only export entries with e_above_hull <= this value (eV/atom).",
    )
    parser.add_argument(
        "--min-e-above-hull",
        type=float,
        default=None,
        help="Only export entries with e_above_hull >= this value (eV/atom).",
    )
    parser.add_argument(
        "structure_ids",
        nargs="*",
        help="Optional list of structure IDs. If omitted, exports all structures.",
    )
    args = parser.parse_args()

    if not args.db.exists():
        print(f"ERROR: Database not found: {args.db}")
        return 1

    target_ids = set(args.structure_ids) if args.structure_ids else None

    print("=" * 70)
    print("PyXtal DB -> CIF export")
    print("=" * 70)
    print(f"Database: {args.db}")
    print(f"Output dir: {args.output_dir}")
    print(f"Target IDs: {len(target_ids) if target_ids else 'all'}")
    print(f"Overwrite: {args.overwrite}")
    if args.min_e_above_hull is not None:
        print(f"Min e_above_hull: {args.min_e_above_hull} eV/atom")
    if args.max_e_above_hull is not None:
        print(f"Max e_above_hull: {args.max_e_above_hull} eV/atom")
    print("=" * 70)

    summary = export_db_to_cif(
        db_path=args.db,
        output_dir=args.output_dir,
        structure_ids=target_ids,
        overwrite=args.overwrite,
        min_e_above_hull=args.min_e_above_hull,
        max_e_above_hull=args.max_e_above_hull,
    )

    print("\nSummary")
    print("-" * 70)
    print(f"Exported: {summary['exported']}")
    print(f"Skipped existing: {summary['skipped_existing']}")
    print(f"Skipped not prescreened: {summary['skipped_not_prescreened']}")
    print(f"Skipped e_above_hull filter: {summary['skipped_e_above_hull']}")
    print(f"Skipped missing required keys: {summary['skipped_missing_required_keys']}")
    print(f"Failed: {len(summary['failed'])}")
    print(f"Missing requested IDs: {len(summary['missing'])}")

    if summary["missing"]:
        print("\nMissing IDs:")
        for sid in summary["missing"]:
            print(f"  - {sid}")

    if summary["failed"]:
        print("\nFailed exports:")
        for sid, err in summary["failed"]:
            print(f"  - {sid}: {err}")

    return 0 if not summary["failed"] else 2


if __name__ == "__main__":
    sys.exit(main())
