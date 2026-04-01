#!/usr/bin/env bash
# Copy available structure assets (PARCHG, meta POSCARs, plots) into one export folder.
#
# Usage:
#   bash results_plot/copy_selected_assets_for_vesta.sh [OUTPUT_DIR]
#
# Environment overrides:
#   ELECTRONIC_ROOT  (default: <repo>/ELECTRONIC_JOBS-Boron)
#   PARCHG_ROOT      (default: <repo>/REFINE_VASP-out-Boron_2;
#                     falls back to ELECTRONIC_ROOT for SPE/PARCHG.tar.gz)
#   BAND_DOS_ROOT    (default: <repo>/band_dos_plots)
#   PHONON_ROOT      (default: <repo>/PHONON_Boron_check2)
#   STABILITY_JSON   (default: <PARCHG_ROOT>/mattersim_stability_results.json)
#   METADATA_CSVS    (optional, colon-separated CSV paths to prefer;
#                     used to restrict/export selected IDs when set explicitly,
#                     otherwise used only as supplemental metadata when present)
#
# Output layout:
#   OUTPUT_DIR/
#     <structure_id>_SPG-<spacegroup>_E-hull_<value>/
#       PARCHG/               # extracted PARCHG files (if tar extract succeeds)
#       plots/                # *_band_dos.png/pdf, phonon_band_dos.png/pdf when available
#       meta/                 # optional copied POSCARs
#       MANIFEST.txt          # source traceability
#
# Notes:
#   - Defaults to structure IDs from STABILITY_JSON when available.
#   - Missing files are ignored and recorded as "missing" in MANIFEST.txt.
#   - No workflow-state validation is performed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

ELECTRONIC_ROOT="${ELECTRONIC_ROOT:-$REPO_ROOT/ELECTRONIC_JOBS-Boron}"
PARCHG_ROOT="${PARCHG_ROOT:-$REPO_ROOT/REFINE_VASP-out-Boron_2}"
BAND_DOS_ROOT="${BAND_DOS_ROOT:-$REPO_ROOT/band_dos_plots}"
PHONON_ROOT="${PHONON_ROOT:-$REPO_ROOT/PHONON_Boron_check2}"
OUT_DIR="${1:-$REPO_ROOT/A-B-Boron-H_Candidates}"

mkdir -p "$OUT_DIR"

declare -A META_SPG
declare -A META_EHULL
declare -A META_SOURCE
declare -A META_SEEN
declare -A PHON_MIN_FREQ
declare -a META_IDS

echo "======================================================================"
echo "Copying selected assets for VESTA"
echo "======================================================================"
echo "Electronic root: $ELECTRONIC_ROOT"
echo "PARCHG root:     $PARCHG_ROOT"
echo "Band DOS root:   $BAND_DOS_ROOT"
echo "Phonon root:     $PHONON_ROOT"
echo "Output dir:      $OUT_DIR"
echo "======================================================================"
echo

load_metadata() {
  local csv_args=()
  local raw_csvs="${METADATA_CSVS:-}"
  local default_csv="$REPO_ROOT/VASP-out-Boron/electride_analysis_candidates-Strict.csv"
  local stability_json="${STABILITY_JSON:-$PARCHG_ROOT/mattersim_stability_results.json}"

  if [[ -n "$raw_csvs" ]]; then
    IFS=':' read -r -a csv_args <<< "$raw_csvs"
  elif [[ -f "$default_csv" ]]; then
    csv_args=("$default_csv")
  fi

  META_IDS=()
  META_SEEN=()
  while IFS=$'\t' read -r sid spg ehull source; do
    [[ -n "$sid" ]] || continue
    META_SPG["$sid"]="$spg"
    META_EHULL["$sid"]="$ehull"
    META_SOURCE["$sid"]="$source"
    if [[ -z "${META_SEEN[$sid]:-}" ]]; then
      META_IDS+=("$sid")
      META_SEEN["$sid"]=1
    fi
  done < <(METADATA_CSVS="$raw_csvs" python3 - "$stability_json" "${csv_args[@]}" <<'PY'
import csv
import json
import os
import sys
from pathlib import Path

stability_json = Path(sys.argv[1]).expanduser() if len(sys.argv) > 1 else None
csv_paths = [Path(p).expanduser() for p in sys.argv[2:] if p]
explicit_csvs = bool(os.environ.get("METADATA_CSVS", "").strip())

sid_keys = ("formula", "structure_id", "struct_id", "id")
spg_keys = ("spacegroup", "space_group_number", "space_group", "spg")
ehull_keys = ("mattersim_e_hull", "e_above_hull", "energy_above_hull", "dft_e_hull")

records = {}
order = []


def ensure_record(sid):
    if sid not in records:
        records[sid] = {"spg": "missing", "ehull": "missing", "source": "missing"}
        order.append(sid)
    return records[sid]

def parse_int(value):
    try:
        text = str(value).strip()
        return str(int(float(text))) if text else "missing"
    except Exception:
        return "missing"

def parse_float(value):
    try:
        text = str(value).strip()
        return f"{float(text):.3f}" if text else "missing"
    except Exception:
        return "missing"

def first_parsed(row, keys, parser):
    for key in keys:
        parsed = parser(row.get(key))
        if parsed != "missing":
            return parsed
    return "missing"


def load_csv(path, *, add_new, allow_ehull_override):
    if not path.exists():
        return

    try:
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                sid = next((str(row.get(k, "")).strip() for k in sid_keys if str(row.get(k, "")).strip()), "")
                if not sid:
                    continue
                if sid not in records and not add_new:
                    continue

                rec = ensure_record(sid)
                spg = first_parsed(row, spg_keys, parse_int)
                ehull = first_parsed(row, ehull_keys, parse_float)

                if spg != "missing" and rec["spg"] == "missing":
                    rec["spg"] = spg

                if ehull != "missing" and (allow_ehull_override or rec["ehull"] == "missing"):
                    rec["ehull"] = ehull
                    rec["source"] = str(path)
                elif rec["source"] == "missing":
                    rec["source"] = str(path)
    except Exception:
        return


def load_stability_json(path, *, add_new):
    if path is None or not path.exists():
        return False

    try:
        with open(path, "r") as handle:
            data = json.load(handle)
    except Exception:
        return False

    if not isinstance(data, dict):
        return False

    results = data.get("results")
    if not isinstance(results, list):
        return False

    loaded_any = False
    for entry in results:
        if not isinstance(entry, dict):
            continue

        sid = str(entry.get("structure_id", "")).strip()
        if not sid:
            continue
        if sid not in records and not add_new:
            continue

        rec = ensure_record(sid)
        spg = first_parsed(entry, spg_keys, parse_int)
        ehull = first_parsed(entry, ehull_keys, parse_float)

        if spg != "missing" and rec["spg"] == "missing":
            rec["spg"] = spg
        if ehull != "missing":
            rec["ehull"] = ehull
        rec["source"] = str(path)
        loaded_any = True

    return loaded_any


json_loaded = False
if explicit_csvs:
    for path in csv_paths:
        load_csv(path, add_new=True, allow_ehull_override=True)
    load_stability_json(stability_json, add_new=False)
else:
    json_loaded = load_stability_json(stability_json, add_new=True)
    for path in csv_paths:
        load_csv(path, add_new=not json_loaded, allow_ehull_override=not json_loaded)

for sid in order:
    rec = records[sid]
    print(f"{sid}\t{rec['spg']}\t{rec['ehull']}\t{rec['source']}")
PY
)
}

load_phonon_summary() {
  local summary_json="$PHONON_ROOT/phonon_postproc_summary.json"
  [[ -f "$summary_json" ]] || return 0

  while IFS=$'\t' read -r sid minfreq; do
    [[ -n "$sid" ]] || continue
    PHON_MIN_FREQ["$sid"]="$minfreq"
  done < <(python3 - "$summary_json" <<'PY'
import json
import sys

path = sys.argv[1]
try:
    data = json.load(open(path))
except Exception:
    sys.exit(0)

if not isinstance(data, dict):
    sys.exit(0)

for sid, meta in data.items():
    if not isinstance(meta, dict):
        continue
    val = meta.get("min_frequency_THz")
    try:
        print(f"{sid}\t{float(val):.6f}")
    except Exception:
        continue
PY
)
}

collect_structure_ids() {
  local root="$1"
  [[ -d "$root" ]] || return 0
  find "$root" -mindepth 2 -maxdepth 2 -type d 2>/dev/null | while IFS= read -r d; do
    basename "$d"
  done
}

structure_prefix_from_sid() {
  local sid="$1"
  printf '%s\n' "${sid%_s*}"
}

first_existing_file() {
  local path=""
  for path in "$@"; do
    if [[ -n "$path" && -f "$path" ]]; then
      printf '%s\n' "$path"
      return 0
    fi
  done
  return 1
}

structure_file_path() {
  local root="$1"
  local sid="$2"
  local relpath="$3"
  local prefix=""

  prefix="$(structure_prefix_from_sid "$sid")"
  first_existing_file \
    "$root/$prefix/$sid/$relpath" \
    "$root/$sid/$relpath"
}

find_plot_file() {
  local filename="$1"
  local sid="${filename%_band_dos.*}"

  first_existing_file \
    "$BAND_DOS_ROOT/$filename" \
    "$BAND_DOS_ROOT/$sid/$filename" \
    "$REPO_ROOT/results_plot/band_dos_plots/$filename" \
    "$REPO_ROOT/results_plot/band_dos_plots/$sid/$filename" \
    "$REPO_ROOT/results_plot/band_dos_plots_all/$filename" \
    "$REPO_ROOT/results_plot/band_dos_plots_all/$sid/$filename"
}

find_parchg_tar() {
  local sid="$1"
  local parchg_tar=""

  parchg_tar="$(structure_file_path "$PARCHG_ROOT" "$sid" "SPE/PARCHG.tar.gz" || true)"
  if [[ -n "$parchg_tar" ]]; then
    printf '%s\n' "$parchg_tar"
    return 0
  fi

  if [[ "$ELECTRONIC_ROOT" != "$PARCHG_ROOT" ]]; then
    parchg_tar="$(structure_file_path "$ELECTRONIC_ROOT" "$sid" "SPE/PARCHG.tar.gz" || true)"
    if [[ -n "$parchg_tar" ]]; then
      printf '%s\n' "$parchg_tar"
      return 0
    fi
  fi

  return 1
}

read_min_frequency_from_band_dat() {
  local band_dat="$1"
  [[ -f "$band_dat" ]] || return 1
  python3 - "$band_dat" <<'PY'
import sys

path = sys.argv[1]
min_freq = None

try:
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            for item in parts[4:]:
                try:
                    value = float(item)
                except Exception:
                    continue
                if min_freq is None or value < min_freq:
                    min_freq = value
except Exception:
    sys.exit(1)

if min_freq is not None:
    print(f"{min_freq:.6f}")
PY
}

copy_if_exists() {
  local src="$1"
  local dst="$2"
  if [[ -n "$src" && -f "$src" ]]; then
    cp -f "$src" "$dst"
    return 0
  fi
  return 1
}

load_metadata
load_phonon_summary

STRUCTURES=("${META_IDS[@]}")

if [[ ${#STRUCTURES[@]} -eq 0 ]]; then
  echo "No structure IDs found in metadata sources."
  exit 1
fi

for sid in "${STRUCTURES[@]}"; do
  spg="${META_SPG[$sid]:-missing}"
  ehull="${META_EHULL[$sid]:-missing}"
  meta_source="${META_SOURCE[$sid]:-missing}"
  phonon_min_freq="${PHON_MIN_FREQ[$sid]:-missing}"
  label="${sid}_SPG-${spg}_E-hull_${ehull}"

  if [[ "$phonon_min_freq" == "missing" ]]; then
    phonon_band_dat="$(structure_file_path "$PHONON_ROOT" "$sid" "PHON/phonon_band.dat" || true)"
    if [[ -n "$phonon_band_dat" ]]; then
      phonon_min_freq="$(read_min_frequency_from_band_dat "$phonon_band_dat" || true)"
      phonon_min_freq="${phonon_min_freq:-missing}"
    fi
  fi

  dst="$OUT_DIR/$label"
  parchg_dst="$dst/PARCHG"
  plots_dst="$dst/plots"
  meta_dst="$dst/meta"
  mkdir -p "$parchg_dst" "$plots_dst" "$meta_dst"

  manifest="$dst/MANIFEST.txt"
  {
    echo "label: $label"
    echo "structure_id: $sid"
    echo "space_group: $spg"
    echo "mattersim_e_above_hull: $ehull"
    echo "phonon_min_frequency_thz: $phonon_min_freq"
    echo "metadata_source: $meta_source"
    echo
  } > "$manifest"

  echo "[${sid}]"

  # 1) PARCHG
  parchg_tar="$(find_parchg_tar "$sid" || true)"
  if [[ -n "$parchg_tar" && -f "$parchg_tar" ]]; then
    echo "  - extracted PARCHG tar: $parchg_tar"
    echo "parchg_tar: $parchg_tar" >> "$manifest"
    tar -xzf "$parchg_tar" -C "$parchg_dst" >/dev/null 2>&1 || true
  else
    echo "  - missing PARCHG tar"
    echo "parchg_tar: missing" >> "$manifest"
  fi

  # Optional POSCARs for reference in VESTA
  spe_poscar="$(structure_file_path "$ELECTRONIC_ROOT" "$sid" "SPE/POSCAR" || true)"
  band_poscar="$(structure_file_path "$ELECTRONIC_ROOT" "$sid" "BAND/POSCAR" || true)"
  copy_if_exists "$spe_poscar" "$meta_dst/POSCAR_SPE" || true
  copy_if_exists "$band_poscar" "$meta_dst/POSCAR_BAND" || true
  echo "poscar_spe: ${spe_poscar:-missing}" >> "$manifest"
  echo "poscar_band: ${band_poscar:-missing}" >> "$manifest"

  # 2) Phonon plot
  phonon_png="$(structure_file_path "$PHONON_ROOT" "$sid" "PHON/phonon_band_dos.png" || true)"
  phonon_pdf="$(structure_file_path "$PHONON_ROOT" "$sid" "PHON/phonon_band_dos.pdf" || true)"
  if copy_if_exists "$phonon_png" "$plots_dst/"; then
    echo "  - copied phonon PNG: $phonon_png"
    echo "phonon_png: $phonon_png" >> "$manifest"
  else
    echo "  - missing phonon PNG"
    echo "phonon_png: missing" >> "$manifest"
  fi
  copy_if_exists "$phonon_pdf" "$plots_dst/" && echo "phonon_pdf: $phonon_pdf" >> "$manifest" || true

  # 3) Band + DOS plot
  band_png="$(find_plot_file "${sid}_band_dos.png" || true)"
  band_pdf="$(find_plot_file "${sid}_band_dos.pdf" || true)"

  if copy_if_exists "$band_png" "$plots_dst/"; then
    echo "  - copied band+DOS PNG: $band_png"
    echo "band_dos_png: $band_png" >> "$manifest"
  else
    echo "  - missing band+DOS PNG"
    echo "band_dos_png: missing" >> "$manifest"
  fi
  copy_if_exists "$band_pdf" "$plots_dst/" && echo "band_dos_pdf: $band_pdf" >> "$manifest" || true

  echo
done

echo "Done. Export folder: $OUT_DIR"
