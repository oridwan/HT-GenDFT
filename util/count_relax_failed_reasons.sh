#!/usr/bin/env bash
set -euo pipefail

# Count RELAX_FAILED structures in workflow.json and classify likely causes.
# Uses workflow error strings first, then inspects small stderr snippets for
# generic/null-error cases (e.g., illegal instruction, timeout).
#
# Usage:
#   ./count_relax_failed_reasons.sh [VASP-out/workflow.json]
#   ./count_relax_failed_reasons.sh VASP-out/workflow.json --list

WORKFLOW_JSON="VASP-out/workflow.json"
SHOW_LISTS=0

for arg in "$@"; do
  case "$arg" in
    --list|-l)
      SHOW_LISTS=1
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
  echo "Error: python3/python is required to parse workflow.json" >&2
  exit 1
fi

"$PYTHON_BIN" - "$WORKFLOW_JSON" "$SHOW_LISTS" <<'PY'
import glob
import json
import os
import sys
from collections import Counter, defaultdict

workflow_json = sys.argv[1]
show_lists = bool(int(sys.argv[2]))


def read_head(path, nbytes=65536):
    if not path or not os.path.isfile(path):
        return ""
    try:
        with open(path, "rb") as f:
            data = f.read(nbytes)
        return data.decode("utf-8", errors="ignore").lower()
    except OSError:
        return ""


def read_tail(path, nbytes=131072):
    if not path or not os.path.isfile(path):
        return ""
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as f:
            if size > nbytes:
                f.seek(-nbytes, os.SEEK_END)
            data = f.read()
        return data.decode("utf-8", errors="ignore").lower()
    except OSError:
        return ""


def latest_glob(pattern):
    files = glob.glob(pattern)
    if not files:
        return ""
    try:
        return max(files, key=os.path.getmtime)
    except OSError:
        # Fallback if a file disappears during scan
        existing = [p for p in files if os.path.exists(p)]
        return max(existing, key=os.path.getmtime) if existing else ""


def resolve_logs(relax_dir, relax_job_id):
    err_file = ""
    out_file = ""
    if not relax_dir or not os.path.isdir(relax_dir):
        return err_file, out_file

    if relax_job_id:
        e = os.path.join(relax_dir, f"vasp_{relax_job_id}.err")
        o = os.path.join(relax_dir, f"vasp_{relax_job_id}.out")
        if os.path.isfile(e):
            err_file = e
        if os.path.isfile(o):
            out_file = o

    if not err_file:
        err_file = latest_glob(os.path.join(relax_dir, "vasp_*.err"))
    if not out_file:
        out_file = latest_glob(os.path.join(relax_dir, "vasp_*.out"))

    return err_file, out_file


TIMEOUT_PATTERNS = (
    "due to time limit",
    "time limit",
    "timelimit",
    "cancelled at",
    "died by signal 9",
)

OOM_PATTERNS = (
    "out of memory",
    "oom",
    "oom-kill",
    "killed process",
    "exceeded step memory limit",
)


def has_any(text, patterns):
    return any(p in text for p in patterns)


def classify(rec):
    wf_error = (rec.get("error") or "").strip()
    wf_error_lc = wf_error.lower()
    relax_dir = rec.get("relax_dir") or ""
    relax_job_id = str(rec.get("relax_job_id") or "")

    # Fast path: workflow.json already knows the reason.
    if "time limit" in wf_error_lc or "timeout" in wf_error_lc:
        return "timeout"
    if "out of memory" in wf_error_lc or "oom" in wf_error_lc:
        return "out_of_memory"
    if "electronic scf not converged" in wf_error_lc or (
        "scf" in wf_error_lc and "converg" in wf_error_lc
    ):
        return "scf_not_converged"
    if "pmg_vasp_psp_dir is not set" in wf_error_lc:
        return "environment_error"

    err_file, out_file = resolve_logs(relax_dir, relax_job_id)
    if not err_file and not out_file:
        return "missing_logs"

    # Illegal instruction tends to appear immediately in stderr.
    err_head = read_head(err_file, 16384)
    if "illegal instruction" in err_head or "forrtl: severe (168)" in err_head:
        return "illegal_instruction"

    # Scheduler errors usually show up near the end.
    err_tail = read_tail(err_file, 32768)
    if has_any(err_tail, TIMEOUT_PATTERNS):
        return "timeout"
    if has_any(err_tail, OOM_PATTERNS):
        return "out_of_memory"

    # Optional deeper stdout parsing is intentionally skipped for speed on large
    # workflows; explicit SCF/timeouts are typically already recorded in workflow.json.

    if "crash" in wf_error_lc or "terminated without completion marker" in wf_error_lc:
        return "crash_other"

    return "other"


with open(workflow_json, "r", encoding="utf-8") as f:
    data = json.load(f)

structures = data.get("structures", {})
counts = Counter()
examples = {}
lists = defaultdict(list)
total = 0

for sid, rec in structures.items():
    if rec.get("state") != "RELAX_FAILED":
        continue
    total += 1
    category = classify(rec)
    counts[category] += 1
    examples.setdefault(category, sid)
    if show_lists:
        lists[category].append(
            f"{sid} | job={rec.get('relax_job_id') or 'NA'} | dir={rec.get('relax_dir') or 'NA'} | workflow_error={rec.get('error') or 'NA'}"
        )

ordered = [
    "illegal_instruction",
    "timeout",
    "out_of_memory",
    "environment_error",
    "scf_not_converged",
    "crash_other",
    "missing_logs",
    "other",
]

print(f"Workflow file: {workflow_json}")
print(f"RELAX_FAILED total: {total}")
print()
print("Counts by category:")

for key in ordered:
    n = counts.get(key, 0)
    line = f"  {key:<18} {n:6d}"
    if key in examples:
        line += f"    (example: {examples[key]})"
    print(line)

for key in sorted(counts.keys()):
    if key in ordered:
        continue
    n = counts[key]
    line = f"  {key:<18} {n:6d}"
    if key in examples:
        line += f"    (example: {examples[key]})"
    print(line)

print()
print("Quick requested summary:")
print(f"  illegal_instruction: {counts.get('illegal_instruction', 0)}")
print(f"  timeout:             {counts.get('timeout', 0)}")
other_total = total - counts.get("illegal_instruction", 0) - counts.get("timeout", 0)
print(f"  other_reasons:       {other_total}")

if show_lists:
    print()
    print("Examples by category (--list enabled):")
    for key in ordered:
        items = lists.get(key, [])
        if not items:
            continue
        print()
        print(f"[{key}]")
        for line in items[:20]:
            print(line)
        if len(items) > 20:
            print(f"... ({len(items)} total entries in this category)")
PY
