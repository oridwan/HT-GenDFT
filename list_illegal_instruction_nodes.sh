python3 phonon_flow_manager.py \
        --refine-jobs REFINE_VASP-out-Boron \
        --output-dir PHONON_JOBS_Boron \
        --structure-ids Cs2Rb1B2H8_s001 \
        --max-concurrent 50#!/usr/bin/env bash
set -euo pipefail

# List nodes that appear in Relax stderr logs with illegal-instruction crashes.
# Intended for building a temporary SBATCH exclude list while you migrate to a
# positive CPU-feature constraint (recommended).
#
# Usage:
#   ./list_illegal_instruction_nodes.sh [VASP-out/workflow.json]
#   ./list_illegal_instruction_nodes.sh VASP-out/workflow.json --csv-only

WORKFLOW_JSON="VASP-out/workflow.json"
CSV_ONLY=0

for arg in "$@"; do
  case "$arg" in
    --csv-only)
      CSV_ONLY=1
      ;;
    *)
      WORKFLOW_JSON="$arg"
      ;;
  esac
done

if [[ ! -f "$WORKFLOW_JSON" ]]; then
  echo "Error: workflow.json not found: $WORKFLOW_JSON" >&2
  exit 1
fi

if command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="python3"
elif command -v python >/dev/null 2>&1; then
  PYTHON_BIN="python"
else
  echo "Error: python3/python is required" >&2
  exit 1
fi

"$PYTHON_BIN" - "$WORKFLOW_JSON" "$CSV_ONLY" <<'PY'
import glob
import json
import os
import re
import sys
from collections import Counter

workflow_json = sys.argv[1]
csv_only = bool(int(sys.argv[2]))


def latest_glob(pattern):
    files = glob.glob(pattern)
    if not files:
        return ""
    try:
        return max(files, key=os.path.getmtime)
    except OSError:
        files = [f for f in files if os.path.exists(f)]
        return max(files, key=os.path.getmtime) if files else ""


def resolve_err(relax_dir, relax_job_id):
    if not relax_dir or not os.path.isdir(relax_dir):
        return ""
    if relax_job_id:
        path = os.path.join(relax_dir, f"vasp_{relax_job_id}.err")
        if os.path.isfile(path):
            return path
    return latest_glob(os.path.join(relax_dir, "vasp_*.err"))


def read_head(path, nbytes=16384):
    if not path or not os.path.isfile(path):
        return ""
    try:
        with open(path, "rb") as f:
            return f.read(nbytes).decode("utf-8", errors="ignore").lower()
    except OSError:
        return ""


def read_tail(path, nbytes=32768):
    if not path or not os.path.isfile(path):
        return ""
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as f:
            if size > nbytes:
                f.seek(-nbytes, os.SEEK_END)
            return f.read().decode("utf-8", errors="ignore").lower()
    except OSError:
        return ""


def extract_node(text):
    # Typical line: srun: error: str-c90: tasks 0-15: Exited with exit code 168
    m = re.search(r"srun:\s*error:\s*([a-z0-9._-]+):", text)
    if not m:
        return ""
    node = m.group(1)
    # Normalize FQDN to short host if present
    return node.split(".")[0]


with open(workflow_json, "r", encoding="utf-8") as f:
    data = json.load(f)

nodes = Counter()
example_struct = {}
total_illegal = 0
missing_node = 0
scanned = 0

for sid, rec in data.get("structures", {}).items():
    if rec.get("state") != "RELAX_FAILED":
        continue

    # Most illegal-instruction entries have error=None; skip known explicit non-illegal reasons.
    wf_error = (rec.get("error") or "").lower()
    if wf_error and ("time" in wf_error or "scf" in wf_error or "pmg_vasp_psp_dir" in wf_error):
        continue

    err_path = resolve_err(rec.get("relax_dir") or "", str(rec.get("relax_job_id") or ""))
    if not err_path:
        continue
    scanned += 1

    head = read_head(err_path)
    if "illegal instruction" not in head and "forrtl: severe (168)" not in head:
        continue

    total_illegal += 1
    tail = read_tail(err_path)
    node = extract_node(tail)
    if not node:
        node = "UNKNOWN"
        missing_node += 1
    nodes[node] += 1
    example_struct.setdefault(node, sid)

bad_nodes = [n for n in nodes if n != "UNKNOWN"]
bad_nodes_sorted = sorted(bad_nodes)
csv = ",".join(bad_nodes_sorted)

if csv_only:
    print(csv)
    sys.exit(0)

print(f"Workflow file: {workflow_json}")
print(f"Illegal-instruction relax failures found: {total_illegal}")
print(f"Err logs scanned: {scanned}")
print()
print("Nodes (sorted by count):")
for node, count in nodes.most_common():
    ex = example_struct.get(node, "")
    if ex:
        print(f"  {node:<16} {count:6d}    (example: {ex})")
    else:
        print(f"  {node:<16} {count:6d}")

print()
print(f"Unique bad nodes ({len(bad_nodes_sorted)}):")
if bad_nodes_sorted:
    print("  " + ", ".join(bad_nodes_sorted))
else:
    print("  <none found>")

print()
print("Ready-to-use sbatch exclude:")
if csv:
    print(f"  --exclude={csv}")
    print()
    print("For workflow manager generated VASP jobs in this repo:")
    print(f"  export VASP_JOB_EXCLUDE_NODES='{csv}'")
    print("  sbatch origin_flow/submit_workflow_manager.sh")
else:
    print("  <none>")

if missing_node:
    print()
    print(f"Note: {missing_node} illegal-instruction logs did not include a parseable srun node line.")
PY

