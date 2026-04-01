#!/usr/bin/env python3
"""
Build a combined electride report PDF with:
1) PARCHG isosurface-style visuals (from PARCHG.tar.gz)
2) Electronic band + DOS plot
3) Phonon band + DOS plot
4) 3D structure view

Typical usage:
    python results_plot/build_multimodal_report.py \
        --electronic-root ELECTRONIC_JOBS-Boron \
        --phonon-root PHONON_JOBS_Boron \
        --structure-ids Cs1Rb1B1H4_s009 \
        --parchg-percentile 97 \
        --iso-backend auto \
        --output-pdf Final_Results/multimodal_report.pdf \
        --assets-dir Final_Results/report_assets

    python results_plot/build_multimodal_report.py \
        --electronic-root ELECTRONIC_JOBS-Boron \
        --phonon-root PHONON_JOBS_Boron \
        --all \
        --parchg-percentile 97 \
        --iso-backend auto \
        --output-pdf results_plot/multimodal_report_all.pdf \
        --assets-dir results_plot/report_assets_all

    /users/oridwan/miniconda3/envs/mattersim/bin/python results_plot/build_multimodal_report.py \
        --electronic-root ELECTRONIC_JOBS-Boron \
        --phonon-root PHONON_JOBS_Boron \
        --all \
        --output-pdf results_plot/multimodal_report_all.pdf
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import importlib.util
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Chgcar
from scipy.ndimage import gaussian_filter

ISO_PLOTLY_WIDTH = 2400
ISO_PLOTLY_HEIGHT = 900
ISO_PLOTLY_SCALE = 3
ISO_MPL_DPI = 320
REPORT_PDF_DPI = 100


def discover_structure_ids(electronic_root: Path) -> List[str]:
    """Discover structure IDs under electronic_root that have SPE/BAND/DOS folders."""
    if not electronic_root.exists():
        return []

    found = set()
    for spe_dir in electronic_root.glob("**/SPE"):
        if not spe_dir.is_dir():
            continue
        struct_dir = spe_dir.parent
        if (struct_dir / "BAND").is_dir() and (struct_dir / "DOS").is_dir():
            found.add(struct_dir.name)
    return sorted(found)


def find_electronic_structure_dir(electronic_root: Path, structure_id: str) -> Optional[Path]:
    """Find ELECTRONIC_JOBS structure directory containing SPE/BAND/DOS for structure_id."""
    for candidate in electronic_root.glob(f"**/{structure_id}"):
        if not candidate.is_dir():
            continue
        if (candidate / "SPE").is_dir() and (candidate / "BAND").is_dir() and (candidate / "DOS").is_dir():
            return candidate
    return None


def find_phonon_dir(phonon_root: Path, structure_id: str) -> Optional[Path]:
    """Find PHON directory for structure_id under phonon_root."""
    for candidate in phonon_root.glob(f"**/{structure_id}/PHON"):
        if candidate.is_dir():
            return candidate
    return None


def load_phonon_postproc_summary(phonon_root: Path) -> Dict[str, dict]:
    """Load phonon post-processing summary JSON if present."""
    summary_path = phonon_root / "phonon_postproc_summary.json"
    if not summary_path.exists():
        return {}
    try:
        with open(summary_path, "r") as f:
            data = json.load(f)
        if isinstance(data, dict):
            return data
    except Exception as exc:
        print(f"  Warning: Failed to read {summary_path}: {exc}")
    return {}


def _parse_int(value: object) -> Optional[int]:
    if value is None:
        return None
    try:
        text = str(value).strip()
        if not text:
            return None
        return int(float(text))
    except Exception:
        return None


def _parse_float(value: object) -> Optional[float]:
    if value is None:
        return None
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except Exception:
        return None


def load_structure_metadata(metadata_csv: Optional[Path]) -> Dict[str, Dict[str, object]]:
    """
    Load per-structure metadata from CSV.
    Expected common columns:
      - structure id: formula (or structure_id/struct_id/id)
      - space group: spacegroup (or space_group_number/space_group/spg)
      - hull energy: e_above_hull (or energy_above_hull/mattersim_e_hull/dft_e_hull)
    """
    if metadata_csv is None or not metadata_csv.exists():
        return {}

    sid_keys = ("formula", "structure_id", "struct_id", "id")
    spg_keys = ("spacegroup", "space_group_number", "space_group", "spg")
    ehull_keys = ("e_above_hull", "energy_above_hull", "mattersim_e_hull", "dft_e_hull")

    out: Dict[str, Dict[str, object]] = {}
    try:
        with open(metadata_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sid = None
                for k in sid_keys:
                    v = row.get(k)
                    if v and str(v).strip():
                        sid = str(v).strip()
                        break
                if not sid:
                    continue

                spg = None
                for k in spg_keys:
                    spg = _parse_int(row.get(k))
                    if spg is not None:
                        break

                ehull = None
                for k in ehull_keys:
                    ehull = _parse_float(row.get(k))
                    if ehull is not None:
                        break

                out[sid] = {"spacegroup": spg, "e_above_hull": ehull}
    except Exception as exc:
        print(f"  Warning: Failed to read metadata CSV {metadata_csv}: {exc}")
        return {}

    return out


def resolve_metadata_csv(preferred_csv: Optional[Path]) -> Optional[Path]:
    """Resolve metadata CSV path; auto-detect when explicit path is not provided."""
    if preferred_csv is not None:
        path = preferred_csv.expanduser().resolve()
        if path.exists():
            return path
        print(f"  Warning: --metadata-csv not found: {path}")
        return None

    candidates = [
        Path("./VASP-out-Boron/electride_analysis.csv"),
        Path("./VASP-out-Boron/electride_analysis_electride_spg_gt15.csv"),
        Path("./REFINE_VASP-out-Boron/stable_electrides.csv"),
    ]
    for c in candidates:
        rc = c.expanduser().resolve()
        if rc.exists():
            return rc
    return None


def read_min_frequency_from_band_dat(phonon_root: Path, structure_id: str) -> Optional[float]:
    """
    Fallback: read minimum phonon frequency (THz) from PHON/phonon_band.dat.
    Expected numeric rows are:
      distance qx qy qz band1 band2 ...
    """
    phon_dir = find_phonon_dir(phonon_root, structure_id)
    if phon_dir is None:
        return None

    band_dat = phon_dir / "phonon_band.dat"
    if not band_dat.exists():
        return None

    min_freq: Optional[float] = None
    try:
        with open(band_dat, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) <= 4:
                    continue
                try:
                    freqs = [float(x) for x in parts[4:]]
                except ValueError:
                    continue
                if not freqs:
                    continue
                row_min = min(freqs)
                if min_freq is None or row_min < min_freq:
                    min_freq = row_min
    except Exception:
        return None

    return float(min_freq) if min_freq is not None else None


def ensure_band_dos_plot(
    structure_id: str,
    electronic_root: Path,
    out_dir: Path,
    force: bool = False,
) -> Optional[Path]:
    """Ensure band+DOS PNG exists for structure_id by calling band_dos_plot.py if needed."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{structure_id}_band_dos.png"
    if out_png.exists() and not force:
        return out_png

    script = Path(__file__).with_name("band_dos_plot.py")
    if not script.exists():
        print(f"  Warning: Missing script {script}")
        return out_png if out_png.exists() else None

    cmd = [
        sys.executable,
        str(script),
        "--base-dir",
        str(electronic_root),
        "--structure-id",
        structure_id,
        "--output-dir",
        str(out_dir),
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"  Warning: band_dos_plot.py failed for {structure_id}: {exc}")
        return out_png if out_png.exists() else None

    return out_png if out_png.exists() else None


def ensure_phonon_plot(
    structure_id: str,
    phonon_root: Path,
    run_postproc_if_missing: bool = False,
) -> Optional[Path]:
    """
    Return phonon_band_dos.png path if available.
    Optionally run refined_flow/postproc_phonon.py for this structure if missing.
    """
    phon_dir = find_phonon_dir(phonon_root, structure_id)
    if phon_dir is None:
        return None

    out_png = phon_dir / "phonon_band_dos.png"
    if out_png.exists():
        return out_png

    if not run_postproc_if_missing:
        return None

    postproc_script = Path(__file__).resolve().parent.parent / "refined_flow" / "postproc_phonon.py"
    if not postproc_script.exists():
        print(f"  Warning: Missing postproc script {postproc_script}")
        return None

    cmd = [
        sys.executable,
        str(postproc_script),
        "--phonon-jobs",
        str(phonon_root),
        "--structure-ids",
        structure_id,
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"  Warning: postproc_phonon.py failed for {structure_id}: {exc}")
        return out_png if out_png.exists() else None

    return out_png if out_png.exists() else None


def _load_chgcar_from_tar_member(tar_path: Path, member_name: str) -> Chgcar:
    """Extract one member from tar.gz to temp file and parse as CHGCAR-like volumetric file."""
    with tarfile.open(tar_path, "r:gz") as tf:
        member = tf.getmember(member_name)
        extracted = tf.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"Could not extract {member_name} from {tar_path}")
        data = extracted.read()

    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp_path = Path(tmp.name)
        tmp.write(data)

    try:
        chg = Chgcar.from_file(str(tmp_path))
    finally:
        tmp_path.unlink(missing_ok=True)
    return chg


def _render_parchg_with_plotly(
    ds: np.ndarray,
    iso: float,
    title_text: str,
    out_png: Path,
    atom_frac_coords: np.ndarray,
    atom_species: List[str],
    show_annotations: bool = True,
) -> bool:
    """Isosurface rendering via plotly + kaleido static export."""
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        cell_len = max(1.0, float(np.ceil(np.max(atom_frac_coords) + 1e-9)))
        gx = np.linspace(0.0, cell_len, ds.shape[0])
        gy = np.linspace(0.0, cell_len, ds.shape[1])
        gz = np.linspace(0.0, cell_len, ds.shape[2])
        x, y, z = np.meshgrid(gx, gy, gz, indexing="ij")
        values = ds.ravel()
        iso_max = float(values.max())
        if not np.isfinite(iso_max) or iso_max <= iso:
            iso_max = float(iso + 1e-6)

        views = [
            dict(eye=dict(x=2.2, y=0.0, z=0.0), up=dict(x=0, y=0, z=1), projection=dict(type="orthographic")),
            dict(eye=dict(x=0.0, y=2.2, z=0.0), up=dict(x=0, y=0, z=1), projection=dict(type="orthographic")),
            dict(eye=dict(x=0.0, y=0.0, z=2.2), up=dict(x=0, y=1, z=0), projection=dict(type="orthographic")),
        ]
        fig = make_subplots(
            rows=1,
            cols=3,
            specs=[[{"type": "scene"}, {"type": "scene"}, {"type": "scene"}]],
            horizontal_spacing=0.02,
        )

        color_map = _element_colors(atom_species)
        frame_edges = [
            ((0.0, 0.0, 0.0), (cell_len, 0.0, 0.0)),
            ((0.0, 0.0, 0.0), (0.0, cell_len, 0.0)),
            ((0.0, 0.0, 0.0), (0.0, 0.0, cell_len)),
            ((cell_len, cell_len, cell_len), (0.0, cell_len, cell_len)),
            ((cell_len, cell_len, cell_len), (cell_len, 0.0, cell_len)),
            ((cell_len, cell_len, cell_len), (cell_len, cell_len, 0.0)),
            ((cell_len, 0.0, 0.0), (cell_len, cell_len, 0.0)),
            ((cell_len, 0.0, 0.0), (cell_len, 0.0, cell_len)),
            ((0.0, cell_len, 0.0), (cell_len, cell_len, 0.0)),
            ((0.0, cell_len, 0.0), (0.0, cell_len, cell_len)),
            ((0.0, 0.0, cell_len), (cell_len, 0.0, cell_len)),
            ((0.0, 0.0, cell_len), (0.0, cell_len, cell_len)),
        ]

        for idx, cam in enumerate(views, start=1):
            fig.add_trace(
                go.Isosurface(
                    x=x.ravel(),
                    y=y.ravel(),
                    z=z.ravel(),
                    value=values,
                    isomin=float(iso),
                    isomax=iso_max,
                    surface_count=1,
                    opacity=0.72,
                    colorscale="Turbo",
                    caps=dict(x_show=False, y_show=False, z_show=False),
                    showscale=False,
                    showlegend=False,
                    name="PARCHG",
                )
                ,
                row=1,
                col=idx,
            )
            for el in sorted(set(atom_species)):
                sel = np.array([s == el for s in atom_species], dtype=bool)
                pts = atom_frac_coords[sel]
                if len(pts) == 0:
                    continue
                fig.add_trace(
                    go.Scatter3d(
                        x=pts[:, 0],
                        y=pts[:, 1],
                        z=pts[:, 2],
                        mode="markers",
                        name=el,
                        legendgroup=el,
                        showlegend=False,
                        marker=dict(
                            size=5,
                            color=mcolors.to_hex(color_map[el]),
                            line=dict(color="black", width=0.6),
                            opacity=0.95,
                        ),
                    ),
                    row=1,
                    col=idx,
                )
            for i, (a, b) in enumerate(frame_edges):
                fig.add_trace(
                    go.Scatter3d(
                        x=[a[0], b[0]],
                        y=[a[1], b[1]],
                        z=[a[2], b[2]],
                        mode="lines",
                        line=dict(color="black", width=2),
                        showlegend=False,
                        hoverinfo="skip",
                        name=f"edge-{i}",
                    ),
                    row=1,
                    col=idx,
                )
            scene_key = "scene" if idx == 1 else f"scene{idx}"
            fig.update_layout(
                **{
                    scene_key: dict(
                        xaxis=dict(visible=False, range=[0.0, cell_len]),
                        yaxis=dict(visible=False, range=[0.0, cell_len]),
                        zaxis=dict(visible=False, range=[0.0, cell_len]),
                        camera=cam,
                        aspectmode="cube",
                        bgcolor="white",
                    )
                }
            )

        fig.update_layout(
            title=(dict(text=title_text, x=0.5, xanchor="center") if show_annotations else None),
            paper_bgcolor="white",
            margin=dict(l=0, r=0, t=(56 if show_annotations else 0), b=0),
            showlegend=False,
        )
        fig.write_image(
            str(out_png),
            width=ISO_PLOTLY_WIDTH,
            height=ISO_PLOTLY_HEIGHT,
            scale=ISO_PLOTLY_SCALE,
        )
        return True
    except Exception:
        return False


def _render_parchg_with_matplotlib(
    ds: np.ndarray,
    mask: np.ndarray,
    title_text: str,
    out_png: Path,
    atom_frac_coords: np.ndarray,
    atom_species: List[str],
    show_annotations: bool = True,
) -> bool:
    """Fallback point-cloud renderer for environments without 3D mesh backends."""
    points = np.argwhere(mask)
    if len(points) == 0:
        return False
    max_points = 18000
    if len(points) > max_points:
        rng = np.random.default_rng(0)
        keep = rng.choice(len(points), size=max_points, replace=False)
        points = points[keep]

    cell_len = max(1.0, float(np.ceil(np.max(atom_frac_coords) + 1e-9)))
    fx = max(ds.shape[0] - 1, 1)
    fy = max(ds.shape[1] - 1, 1)
    fz = max(ds.shape[2] - 1, 1)
    coords = np.column_stack(
        [
            points[:, 0] / fx * cell_len,
            points[:, 1] / fy * cell_len,
            points[:, 2] / fz * cell_len,
        ]
    )

    vals = ds[points[:, 0], points[:, 1], points[:, 2]]
    val_norm = (vals - vals.min()) / (vals.max() - vals.min() + 1e-12)

    fig = plt.figure(figsize=(12.6, 4.4))
    views = [
        (0, 0),
        (0, 90),
        (90, -90),
    ]
    color_map = _element_colors(atom_species)
    edges = [
        ((0.0, 0.0, 0.0), (cell_len, 0.0, 0.0)),
        ((0.0, 0.0, 0.0), (0.0, cell_len, 0.0)),
        ((0.0, 0.0, 0.0), (0.0, 0.0, cell_len)),
        ((cell_len, cell_len, cell_len), (0.0, cell_len, cell_len)),
        ((cell_len, cell_len, cell_len), (cell_len, 0.0, cell_len)),
        ((cell_len, cell_len, cell_len), (cell_len, cell_len, 0.0)),
        ((cell_len, 0.0, 0.0), (cell_len, cell_len, 0.0)),
        ((cell_len, 0.0, 0.0), (cell_len, 0.0, cell_len)),
        ((0.0, cell_len, 0.0), (cell_len, cell_len, 0.0)),
        ((0.0, cell_len, 0.0), (0.0, cell_len, cell_len)),
        ((0.0, 0.0, cell_len), (cell_len, 0.0, cell_len)),
        ((0.0, 0.0, cell_len), (0.0, cell_len, cell_len)),
    ]

    for i, (elev, azim) in enumerate(views, start=1):
        ax = fig.add_subplot(1, 3, i, projection="3d")
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=val_norm, cmap="turbo", s=2.4, alpha=0.52, linewidths=0)
        for el in sorted(set(atom_species)):
            sel = np.array([s == el for s in atom_species], dtype=bool)
            pts = atom_frac_coords[sel]
            if len(pts) == 0:
                continue
            ax.scatter(
                pts[:, 0],
                pts[:, 1],
                pts[:, 2],
                s=36,
                c=[color_map[el]],
                edgecolors="black",
                linewidths=0.32,
                alpha=0.98,
                label=el if i == 1 else None,
            )
        for a, b in edges:
            ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], color="black", alpha=0.26, linewidth=0.8)
        ax.set_xlim(0.0, cell_len)
        ax.set_ylim(0.0, cell_len)
        ax.set_zlim(0.0, cell_len)
        ax.set_proj_type("ortho")
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
    if show_annotations:
        fig.suptitle(title_text, fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout(rect=[0, 0, 1, 1])
    fig.savefig(out_png, dpi=ISO_MPL_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return True


def render_parchg_voxel_isosurface(
    tar_path: Path,
    label: str,
    out_png: Path,
    percentile: float = 99.5,
    max_grid: int = 64,
    iso_backend: str = "auto",
    abs_isovalue_ea3: Optional[float] = None,
    show_annotations: bool = True,
) -> bool:
    """
    Render an isosurface-style voxel view from a selected PARCHG-* file in PARCHG.tar.gz.
    Backend order is configurable (`auto`: plotly -> matplotlib).
    """
    candidate_names = [label, f"PARCHG-{label}"] if not label.startswith("PARCHG-") else [label]

    with tarfile.open(tar_path, "r:gz") as tf:
        members = {Path(m.name).name: m for m in tf.getmembers() if m.isfile()}

    member_name = None
    for cand in candidate_names:
        if cand in members:
            member_name = cand
            break
    if member_name is None:
        print(f"  Warning: {label} not found in {tar_path.name}")
        return False

    chg = _load_chgcar_from_tar_member(tar_path, member_name)
    data = chg.data.get("total")
    if data is None:
        # Fallback to first available volumetric key
        data = next(iter(chg.data.values()))

    # Convert to charge density in e/Ang^3.
    volume = float(chg.structure.lattice.volume)
    if volume <= 0:
        return False
    data = data.astype(float, copy=False) / volume

    # Downsample unit-cell grid, then tile to 2x2x2 for better connectivity view.
    nx, ny, nz = data.shape
    target_uc = max(8, int(max_grid // 2))
    step = max(1, int(np.ceil(max(nx, ny, nz) / target_uc)))
    ds_uc = data[::step, ::step, ::step].astype(float, copy=False)
    parchg_supercell = (2, 2, 2)
    ds = np.tile(ds_uc, parchg_supercell)

    if ds.size == 0:
        return False

    # Mild smoothing gives much cleaner surfaces than raw noisy voxels.
    ds = gaussian_filter(ds, sigma=0.8)

    use_absolute_iso = abs_isovalue_ea3 is not None
    if use_absolute_iso:
        iso = float(abs_isovalue_ea3)
        iso_label = f"{iso:.4g} e/Å³"
    else:
        iso = float(np.percentile(ds, percentile))
        iso_label = f"p{percentile:.1f}"
    mask = ds >= iso

    # If percentile threshold too strict, relax a bit.
    if not use_absolute_iso and np.count_nonzero(mask) < 50:
        iso = float(np.percentile(ds, max(90.0, percentile - 5.0)))
        mask = ds >= iso
        iso_label = f"p~{percentile:.1f}"

    if np.count_nonzero(mask) == 0:
        if use_absolute_iso:
            print(
                f"  Warning: No voxels above absolute isovalue {abs_isovalue_ea3:.4g} e/Å³ "
                f"(range: {ds.min():.3g}..{ds.max():.3g} e/Å³)"
            )
        return False

    out_png.parent.mkdir(parents=True, exist_ok=True)
    atom_frac_uc = np.mod(chg.structure.frac_coords, 1.0)
    atom_species_uc = [str(site.specie.symbol) for site in chg.structure.sites]
    shifts = np.array([[i, j, k] for i in range(2) for j in range(2) for k in range(2)], dtype=float)
    atom_frac_coords = np.vstack([atom_frac_uc + s for s in shifts])
    atom_species = [sp for _ in range(len(shifts)) for sp in atom_species_uc]
    title_text = f"{member_name} (iso {iso_label}, 2x2x2)"

    backend = iso_backend.lower()
    if backend not in {"auto", "plotly", "matplotlib"}:
        print(f"  Warning: Unknown iso backend '{iso_backend}', using auto")
        backend = "auto"

    if backend == "auto":
        order = ["plotly", "matplotlib"]
    else:
        order = [backend]

    for name in order:
        if name == "plotly":
            if _render_parchg_with_plotly(
                ds,
                iso,
                title_text,
                out_png,
                atom_frac_coords,
                atom_species,
                show_annotations=show_annotations,
            ):
                return True
        else:
            if _render_parchg_with_matplotlib(
                ds,
                mask,
                title_text,
                out_png,
                atom_frac_coords,
                atom_species,
                show_annotations=show_annotations,
            ):
                return True

    return False


def render_missing_panel(ax, title: str, msg: str, title_color: str = "black") -> None:
    ax.set_title(title, fontsize=10, pad=2, color=title_color)
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=9)
    ax.axis("off")


def _trim_white_border(img: np.ndarray, white_threshold: float = 0.985, pad: int = 4) -> np.ndarray:
    """Trim near-white border from image arrays."""
    arr = img
    if arr.ndim == 2:
        mask = arr < white_threshold
    else:
        rgb = arr[..., :3]
        mask = np.any(rgb < white_threshold, axis=-1)
    if not np.any(mask):
        return img

    ys, xs = np.where(mask)
    y0, y1 = int(ys.min()), int(ys.max()) + 1
    x0, x1 = int(xs.min()), int(xs.max()) + 1
    y0 = max(0, y0 - pad)
    x0 = max(0, x0 - pad)
    y1 = min(arr.shape[0], y1 + pad)
    x1 = min(arr.shape[1], x1 + pad)
    return arr[y0:y1, x0:x1]


def render_image_panel(
    ax,
    title: str,
    image_path: Path,
    trim_white: bool = False,
    title_color: str = "black",
) -> None:
    ax.set_title(title, fontsize=10, pad=2, color=title_color)
    try:
        img = plt.imread(image_path)
        if trim_white:
            img = _trim_white_border(img)
        ax.imshow(img)
        ax.axis("off")
    except Exception as exc:
        render_missing_panel(ax, title, f"Failed to load image\n{image_path}\n{exc}", title_color=title_color)


def _element_colors(species: List[str]) -> Dict[str, Tuple[float, float, float, float]]:
    """Assign deterministic colors per element symbol."""
    unique = sorted(set(species))
    cmap = plt.colormaps["tab20"]
    if len(unique) == 1:
        return {unique[0]: cmap(0.0)}
    return {el: cmap(i / (len(unique) - 1)) for i, el in enumerate(unique)}


def render_structure_3d_panel(ax, structure: Structure, title: str) -> None:
    """Render a lightweight 3D structure view from a pymatgen Structure."""
    coords = structure.cart_coords
    species = [str(site.specie.symbol) for site in structure.sites]
    color_map = _element_colors(species)

    # Size points by approximate atomic radius where available.
    sizes = []
    for site in structure.sites:
        r = getattr(site.specie, "atomic_radius", None) or getattr(site.specie, "atomic_radius_calculated", None)
        if r is None:
            r = 1.0
        sizes.append(float(r) * 70.0)

    colors = [color_map[s] for s in species]
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=sizes, c=colors, alpha=0.95, edgecolors="k", linewidths=0.2)

    # Equal-ish axis scaling
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    center = (mins + maxs) / 2.0
    half = max(np.max(maxs - mins) / 2.0, 1.0)
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)

    ax.set_title(title, fontsize=11)
    ax.view_init(elev=22, azim=38)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    # Compact legend (unique elements only)
    handles = []
    labels = []
    for el in sorted(set(species)):
        handles.append(plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[el], markeredgecolor="k", markersize=6))
        labels.append(el)
    if handles:
        ax.legend(handles, labels, loc="upper right", fontsize=7, framealpha=0.7)


def write_structure_py3dmol_html(
    structure: Structure,
    out_html: Path,
    js_source: Optional[str] = None,
) -> bool:
    """
    Export an interactive py3Dmol HTML viewer for a pymatgen Structure.
    Returns True on success.
    """
    if importlib.util.find_spec("py3Dmol") is None:
        return False

    try:
        import py3Dmol
    except Exception:
        return False

    try:
        cif_str = structure.to(fmt="cif")
        if js_source:
            view = py3Dmol.view(width=920, height=720, js=js_source)
        else:
            view = py3Dmol.view(width=920, height=720)
        view.addModel(cif_str, "cif")
        view.setStyle({"sphere": {"scale": 0.30}, "stick": {"radius": 0.13}})
        view.addUnitCell()
        view.zoomTo()
        out_html.parent.mkdir(parents=True, exist_ok=True)
        view.write_html(str(out_html))
        return True
    except Exception:
        return False


def build_report(
    structure_ids: List[str],
    electronic_root: Path,
    phonon_root: Path,
    output_pdf: Path,
    assets_dir: Path,
    parchg_labels: List[str],
    parchg_percentile: float,
    parchg_abs_isovalue: Optional[float],
    max_grid: int,
    iso_backend: str,
    py3dmol_js_path: Optional[Path],
    structure_supercell: Tuple[int, int, int],
    run_postproc_phonon: bool,
    force_regen_band_dos: bool,
    structure_metadata: Optional[Dict[str, Dict[str, object]]] = None,
) -> Dict[str, Dict[str, str]]:
    """Generate combined report PDF and return per-structure status summary."""
    status: Dict[str, Dict[str, str]] = {}
    phonon_summary = load_phonon_postproc_summary(phonon_root)
    structure_metadata = structure_metadata or {}

    band_dos_assets = assets_dir / "band_dos"
    parchg_assets = assets_dir / "parchg"
    structure_assets = assets_dir / "structure_3dmol"
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    structure_assets.mkdir(parents=True, exist_ok=True)

    py3dmol_js_source: Optional[str] = None
    if py3dmol_js_path is not None:
        local_js = structure_assets / "3Dmol-min.js"
        try:
            shutil.copy2(py3dmol_js_path, local_js)
            py3dmol_js_source = "./3Dmol-min.js"
        except Exception as exc:
            print(f"  Warning: failed to copy py3Dmol JS from {py3dmol_js_path}: {exc}")

    with PdfPages(output_pdf) as pdf:
        for sid in structure_ids:
            print(f"\nProcessing {sid} ...")
            status[sid] = {}

            meta = structure_metadata.get(sid, {})
            sid_spg = _parse_int(meta.get("spacegroup")) if isinstance(meta, dict) else None
            sid_ehull = _parse_float(meta.get("e_above_hull")) if isinstance(meta, dict) else None
            status[sid]["spacegroup"] = str(sid_spg) if sid_spg is not None else "missing"
            status[sid]["e_above_hull"] = f"{sid_ehull:.6f}" if sid_ehull is not None else "missing"

            e_dir = find_electronic_structure_dir(electronic_root, sid)
            if e_dir is None:
                print("  Missing electronic directory")
                status[sid]["electronic"] = "missing"
                continue
            status[sid]["electronic"] = str(e_dir)

            spe_dir = e_dir / "SPE"
            tar_path = spe_dir / "PARCHG.tar.gz"

            # Structure 3D view source (prefer SPE POSCAR)
            structure = None
            structure_poscar = None
            for poscar_candidate in [e_dir / "SPE" / "POSCAR", e_dir / "BAND" / "POSCAR", e_dir / "DOS" / "POSCAR"]:
                if poscar_candidate.exists():
                    try:
                        structure = Structure.from_file(str(poscar_candidate))
                        structure_poscar = poscar_candidate
                        break
                    except Exception:
                        continue
            status[sid]["structure_poscar"] = str(structure_poscar) if structure_poscar else "missing"
            status[sid]["structure_supercell"] = "x".join(str(i) for i in structure_supercell)

            structure_for_plot = None
            if structure is not None:
                structure_for_plot = structure.copy()
                try:
                    structure_for_plot.make_supercell(list(structure_supercell))
                except Exception as exc:
                    status[sid]["structure_supercell_error"] = str(exc)
                    structure_for_plot = structure

            # Interactive 3D structure viewer (py3Dmol HTML)
            struct_html = None
            if structure_for_plot is not None:
                candidate_html = structure_assets / f"{sid}_structure_3dmol.html"
                if write_structure_py3dmol_html(
                    structure_for_plot,
                    candidate_html,
                    js_source=py3dmol_js_source,
                ):
                    struct_html = candidate_html
            status[sid]["structure_py3dmol_html"] = str(struct_html) if struct_html else "missing"
            status[sid]["structure_py3dmol_js"] = str(py3dmol_js_path) if py3dmol_js_path else "CDN"

            # Band+DOS plot
            band_png = ensure_band_dos_plot(
                sid,
                electronic_root=electronic_root,
                out_dir=band_dos_assets,
                force=force_regen_band_dos,
            )
            status[sid]["band_dos"] = str(band_png) if band_png else "missing"

            # Phonon plot
            phonon_png = ensure_phonon_plot(
                sid,
                phonon_root=phonon_root,
                run_postproc_if_missing=run_postproc_phonon,
            )
            status[sid]["phonon"] = str(phonon_png) if phonon_png else "missing"
            phon_sum = phonon_summary.get(sid, {}) if isinstance(phonon_summary.get(sid), dict) else {}
            min_freq_thz = phon_sum.get("min_frequency_THz")
            min_freq_source = "summary"
            if min_freq_thz is None:
                min_freq_thz = read_min_frequency_from_band_dat(phonon_root, sid)
                min_freq_source = "phonon_band.dat" if min_freq_thz is not None else "missing"

            min_freq_val: Optional[float] = None
            try:
                if min_freq_thz is not None:
                    min_freq_val = float(min_freq_thz)
            except (TypeError, ValueError):
                min_freq_val = None

            if min_freq_val is not None:
                status[sid]["min_frequency_THz"] = f"{min_freq_val:.6f}"
                status[sid]["min_frequency_source"] = min_freq_source
                phonon_title = f"Phonon Band + DOS (min freq: {min_freq_val:.3f} THz)"
                phonon_title_color = "red" if min_freq_val < -0.5 else "green"
            else:
                status[sid]["min_frequency_THz"] = "missing"
                status[sid]["min_frequency_source"] = "missing"
                phonon_title = "Phonon Band + DOS"
                phonon_title_color = "black"

            # PARCHG visuals (up to first two labels for report panel)
            parchg_pngs: List[Tuple[str, Path]] = []
            if tar_path.exists():
                for label in parchg_labels[:2]:
                    out_png = parchg_assets / f"{sid}_{label}.png"
                    ok = render_parchg_voxel_isosurface(
                        tar_path=tar_path,
                        label=label,
                        out_png=out_png,
                        percentile=parchg_percentile,
                        max_grid=max_grid,
                        iso_backend=iso_backend,
                        abs_isovalue_ea3=parchg_abs_isovalue,
                        show_annotations=False,
                    )
                    if ok:
                        parchg_pngs.append((label, out_png))
            status[sid]["parchg_tar"] = str(tar_path) if tar_path.exists() else "missing"
            status[sid]["parchg_images"] = ", ".join(str(p) for _, p in parchg_pngs) if parchg_pngs else "missing"

            # Build one page per structure
            fig = plt.figure(figsize=(13.2, 16.2), constrained_layout=False)
            gs = fig.add_gridspec(
                4,
                3,
                width_ratios=[1.0, 1.0, 1.0],
                height_ratios=[1.15, 1.15, 1.0, 1.45],
                wspace=0.04,
                hspace=0.07,
            )
            ax_p0 = fig.add_subplot(gs[0, :])
            ax_p1 = fig.add_subplot(gs[1, :])
            ax_struct = fig.add_subplot(gs[2, 0], projection="3d")
            ax_band = fig.add_subplot(gs[2, 1:3])
            ax_phon = fig.add_subplot(gs[3, :])
            title_chunks = [sid]
            if sid_spg is not None:
                title_chunks.append(f"SPG-{sid_spg}")
            if sid_ehull is not None:
                title_chunks.append(f"E-hull_{sid_ehull:.3f}")
            fig.suptitle("_".join(title_chunks), fontsize=14, fontweight="bold", y=0.992)
            fig.subplots_adjust(left=0.02, right=0.995, top=0.968, bottom=0.018, wspace=0.04, hspace=0.08)

            # Structure panel
            if structure_for_plot is not None:
                render_structure_3d_panel(
                    ax_struct,
                    structure_for_plot,
                    f"3D Structure ({structure_supercell[0]}x{structure_supercell[1]}x{structure_supercell[2]})",
                )
                if struct_html is not None:
                    ax_struct.text2D(
                        0.02,
                        0.02,
                        f"Interactive py3Dmol: {struct_html.name}",
                        transform=ax_struct.transAxes,
                        fontsize=7,
                        ha="left",
                    )
            else:
                render_missing_panel(ax_struct, "3D Structure", "POSCAR not found")

            # PARCHG panels
            if len(parchg_pngs) > 0:
                render_image_panel(ax_p0, f"PARCHG Isosurface ({parchg_pngs[0][0]})", parchg_pngs[0][1], trim_white=True)
            else:
                render_missing_panel(ax_p0, "PARCHG Isosurface", "No PARCHG visual generated")

            if len(parchg_pngs) > 1:
                render_image_panel(ax_p1, f"PARCHG Isosurface ({parchg_pngs[1][0]})", parchg_pngs[1][1], trim_white=True)
            else:
                render_missing_panel(ax_p1, "PARCHG Isosurface", "Second PARCHG panel not available")

            # Band + DOS panel
            if band_png is not None and band_png.exists():
                render_image_panel(ax_band, "Electronic Band + DOS", band_png)
            else:
                render_missing_panel(ax_band, "Electronic Band + DOS", "Plot missing")

            # Phonon panel
            if phonon_png is not None and phonon_png.exists():
                render_image_panel(ax_phon, phonon_title, phonon_png, title_color=phonon_title_color)
            else:
                render_missing_panel(
                    ax_phon,
                    phonon_title,
                    "Plot missing (run postproc_phonon.py after PHON_DONE)",
                    title_color=phonon_title_color,
                )

            pdf.savefig(fig, dpi=REPORT_PDF_DPI)
            plt.close(fig)

    return status


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build combined PDF report with PARCHG visuals, band+DOS, and phonon plots.",
    )
    parser.add_argument(
        "--electronic-root",
        type=Path,
        default=Path("./ELECTRONIC_JOBS-Boron"),
        help="Path to electronic jobs root (default: ./ELECTRONIC_JOBS-Boron)",
    )
    parser.add_argument(
        "--phonon-root",
        type=Path,
        default=Path("./PHONON_JOBS_Boron"),
        help="Path to phonon jobs root (default: ./PHONON_JOBS_Boron)",
    )
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--structure-ids", nargs="+", help="Structure IDs to process")
    target.add_argument("--all", action="store_true", help="Process all structures found under --electronic-root")
    parser.add_argument(
        "--output-pdf",
        type=Path,
        default=Path("./multimodal_report.pdf"),
        help="Output report PDF path",
    )
    parser.add_argument(
        "--assets-dir",
        type=Path,
        default=Path("./report_assets"),
        help="Directory for generated intermediate images",
    )
    parser.add_argument(
        "--parchg-labels",
        nargs="+",
        default=["band0", "band1"],
        help="PARCHG labels to visualize (default: band0 band1)",
    )
    parser.add_argument(
        "--parchg-percentile",
        type=float,
        default=99.5,
        help="Percentile threshold for isosurface mask (default: 99.5)",
    )
    parser.add_argument(
        "--parchg-abs-isovalue",
        type=float,
        default=None,
        help="Absolute isovalue in e/Å^3 (overrides --parchg-percentile when set)",
    )
    parser.add_argument(
        "--max-grid",
        type=int,
        default=64,
        help="Maximum volumetric grid dimension after downsampling (default: 64)",
    )
    parser.add_argument(
        "--iso-backend",
        choices=["auto", "plotly", "matplotlib"],
        default="auto",
        help="Backend for PARCHG isosurface rendering (default: auto)",
    )
    parser.add_argument(
        "--py3dmol-js-path",
        type=Path,
        default=None,
        help="Local path to 3Dmol-min.js for offline py3Dmol HTML rendering (optional)",
    )
    parser.add_argument(
        "--structure-supercell",
        nargs=3,
        type=int,
        default=[2, 2, 2],
        metavar=("NA", "NB", "NC"),
        help="Supercell for 3D structure view (default: 2 2 2)",
    )
    parser.add_argument(
        "--run-postproc-phonon",
        action="store_true",
        help="Run refined_flow/postproc_phonon.py for missing phonon plots",
    )
    parser.add_argument(
        "--force-regen-band-dos",
        action="store_true",
        help="Regenerate band+DOS plots even if cached images exist",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="Optional path to write per-structure generation summary JSON",
    )
    parser.add_argument(
        "--metadata-csv",
        type=Path,
        default=None,
        help=(
            "Optional CSV containing structure metadata columns like "
            "formula/spacegroup/e_above_hull. "
            "If omitted, script auto-detects common files."
        ),
    )

    args = parser.parse_args()

    electronic_root = args.electronic_root.expanduser().resolve()
    phonon_root = args.phonon_root.expanduser().resolve()
    output_pdf = args.output_pdf.expanduser().resolve()
    assets_dir = args.assets_dir.expanduser().resolve()

    if not electronic_root.exists():
        raise FileNotFoundError(f"Electronic root not found: {electronic_root}")

    if args.all:
        structure_ids = discover_structure_ids(electronic_root)
        print(f"Discovered {len(structure_ids)} structures under {electronic_root}")
    else:
        structure_ids = args.structure_ids or []

    if not structure_ids:
        print("No structures to process.")
        return

    metadata_csv = resolve_metadata_csv(args.metadata_csv)
    structure_metadata = load_structure_metadata(metadata_csv)
    if metadata_csv is not None:
        print(f"Using metadata CSV: {metadata_csv} ({len(structure_metadata)} rows parsed)")
    else:
        print("No metadata CSV found; title will use structure ID only.")

    if any(i <= 0 for i in args.structure_supercell):
        raise ValueError(f"--structure-supercell values must be positive: {args.structure_supercell}")
    if args.parchg_abs_isovalue is not None and args.parchg_abs_isovalue <= 0:
        raise ValueError(f"--parchg-abs-isovalue must be positive, got {args.parchg_abs_isovalue}")
    py3dmol_js_path = None
    if args.py3dmol_js_path is not None:
        py3dmol_js_path = args.py3dmol_js_path.expanduser().resolve()
        if not py3dmol_js_path.exists():
            raise FileNotFoundError(f"--py3dmol-js-path not found: {py3dmol_js_path}")

    status = build_report(
        structure_ids=structure_ids,
        electronic_root=electronic_root,
        phonon_root=phonon_root,
        output_pdf=output_pdf,
        assets_dir=assets_dir,
        parchg_labels=args.parchg_labels,
        parchg_percentile=args.parchg_percentile,
        parchg_abs_isovalue=args.parchg_abs_isovalue,
        max_grid=args.max_grid,
        iso_backend=args.iso_backend,
        py3dmol_js_path=py3dmol_js_path,
        structure_supercell=tuple(args.structure_supercell),
        run_postproc_phonon=args.run_postproc_phonon,
        force_regen_band_dos=args.force_regen_band_dos,
        structure_metadata=structure_metadata,
    )

    if args.summary_json:
        summary_path = args.summary_json.expanduser().resolve()
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as f:
            json.dump(status, f, indent=2)
        print(f"Saved summary: {summary_path}")

    print(f"\nSaved report: {output_pdf}")


if __name__ == "__main__":
    main()
