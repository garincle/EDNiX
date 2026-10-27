"""
Supplementary Table 1 of the paper: acquisition parameters of every dataset (Clément,
2026-10-27: "une figure supplémentaire avec les caractéristiques des données ? champ IRM ?
TR ?").

Read from the raw BIDS tree only (never from EDNiX outputs): for every image,
  the JSON sidecar (BIDS inheritance: <bids>/<suffix>.json, <bids>/task-<x>_<suffix>.json,
  then the sidecar next to the image, the latter wins) and the NIfTI header (dimensions,
  voxel size, pixdim[4] with its time unit). The image data are not loaded.

  a  resting-state fMRI, datasets of Figures 4-5 (paper_data.FUNC_DATASETS); the BOLD
     files are those the run selection reads (endfmri of config_ednix.BIDS_CONFIGS),
     the arousal state is the one given there.
  b  anatomy, datasets of Figure 2 (paper_data.ANAT_DATASETS), T1w and T2w images.

Per dataset: n subjects / sessions / images, field strength, scanner, TR, TE, flip
angle, voxel size, matrix, volumes and duration (fMRI): median and [min-max] across
images when they differ.

Values are reported, never corrected. A TR is flagged (column TR_check of the per-image
table, and 'TR flagged' in the dataset row) when the sidecar value is missing, outside
0.1-10 s, or differs by more than 5 % from pixdim[4] of the NIfTI header; the TR printed
is then the header value when that one is plausible. Example: BIDS_Grandjean (rat)
sidecars give RepetitionTime = 0.001, i.e. the TR in ms written in a field in s.
"""
import os
import re
import glob
import json

import numpy as np
import pandas as pd

from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis import paper_style as ps

opj = os.path.join
MM = 1.0 / 25.4

TR_RANGE = (0.1, 10.0)          # s; a TR outside is a unit error in the sidecar
TR_TOLERANCE = 0.05             # sidecar vs NIfTI header
ANAT_SUFFIXES = ("T1w", "T2w")


# ── arousal state and BOLD pattern of each fMRI dataset (run selection config) ─────────
def _func_config():
    """{BIDS folder: (endfmri pattern, state text)} from config_ednix.BIDS_CONFIGS."""
    try:
        from Statistics.Group_fMRI.run_selection.config_ednix import BIDS_CONFIGS
    except Exception:
        return {}
    out = {}
    for c in BIDS_CONFIGS:
        st = c.get("state")
        if isinstance(st, dict):
            m = st.get("map", {})
            vals = m.values() if isinstance(m, dict) else [v for _, v in m]
            st = " / ".join(dict.fromkeys(vals)) + f" (from {st.get('from', '?')})"
        out[pdat.bids_name(c["bids_dir"])] = (c.get("endfmri", "*_bold.nii*"), st or "unknown")
    return out


# ── one image ──────────────────────────────────────────────────────────────────────────
def _entity(path, key):
    m = re.search(rf"(?:^|_){key}-([A-Za-z0-9]+)", os.path.basename(path))
    return m.group(1) if m else None


def _load_json(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except Exception:
        return {}


def _sidecar(img, bids_dir, suffix):
    """Merged sidecar, BIDS inheritance from the dataset root down to the image."""
    info = {}
    task = _entity(img, "task")
    for p in [opj(bids_dir, f"{suffix}.json")] + \
             ([opj(bids_dir, f"task-{task}_{suffix}.json")] if task else []):
        info.update(_load_json(p))
    stem = re.sub(r"\.nii(\.gz)?$", "", img)
    info.update(_load_json(stem + ".json"))
    return info


def _header(img):
    """(shape, voxel size mm, TR from pixdim[4] in s or nan)."""
    import nibabel as nib
    try:
        h = nib.load(img).header
    except Exception:
        return None, None, np.nan
    shape = tuple(int(s) for s in h.get_data_shape())
    zooms = tuple(float(z) for z in h.get_zooms())
    tr = np.nan
    if len(zooms) > 3 and zooms[3] > 0:
        unit = h.get_xyzt_units()[1]
        tr = zooms[3] * {"msec": 1e-3, "usec": 1e-6}.get(unit, 1.0)
    return shape, zooms[:3], tr


def _ok(tr):
    return tr is not None and np.isfinite(tr) and TR_RANGE[0] <= tr <= TR_RANGE[1]


def _image_row(img, bids_dir, species, suffix):
    info = _sidecar(img, bids_dir, suffix)
    shape, vox, tr_hdr = _header(img)
    row = dict(species=species, dataset=pdat.bids_name(bids_dir), image=img, suffix=suffix,
               subject=_entity(img, "sub"), session=_entity(img, "ses") or "1",
               run=_entity(img, "run"), task=_entity(img, "task"),
               field_T=info.get("MagneticFieldStrength"),
               manufacturer=info.get("Manufacturer"),
               model=info.get("ManufacturersModelName"),
               TE_ms=1e3 * float(info["EchoTime"]) if "EchoTime" in info else np.nan,
               flip_deg=float(info["FlipAngle"]) if "FlipAngle" in info else np.nan,
               multiband=info.get("MultibandAccelerationFactor"),
               vox_x=vox[0] if vox else np.nan, vox_y=vox[1] if vox else np.nan,
               vox_z=vox[2] if vox else np.nan,
               matrix="x".join(str(s) for s in shape[:3]) if shape else None,
               n_volumes=shape[3] if shape and len(shape) > 3 else 1)
    if suffix == "bold":
        tr_js = float(info["RepetitionTime"]) if "RepetitionTime" in info else np.nan
        if not np.isfinite(tr_js):
            check = "missing in sidecar"
        elif not _ok(tr_js):
            check = f"sidecar {tr_js:g} s outside {TR_RANGE[0]}-{TR_RANGE[1]} s"
        elif _ok(tr_hdr) and abs(tr_js - tr_hdr) / tr_hdr > TR_TOLERANCE:
            check = f"sidecar {tr_js:g} s vs header {tr_hdr:g} s"
        else:
            check = ""
        tr = tr_js if not check else (tr_hdr if _ok(tr_hdr) else np.nan)
        row.update(TR_sidecar_s=tr_js, TR_header_s=tr_hdr, TR_s=tr, TR_check=check,
                   duration_min=tr * row["n_volumes"] / 60 if np.isfinite(tr) else np.nan)
    return row


# ── datasets ───────────────────────────────────────────────────────────────────────────
def _images(bids_dir, kind, pattern):
    out = []
    for lev in (("sub-*", "ses-*"), ("sub-*",)):
        out += glob.glob(opj(bids_dir, *lev, kind, pattern))
    return sorted(p for p in set(out) if re.search(r"\.nii(\.gz)?$", p))


def scan(func_datasets=None, anat_datasets=None):
    """Per-image table (fMRI and anatomy)."""
    func_datasets = pdat.FUNC_DATASETS if func_datasets is None else func_datasets
    anat_datasets = pdat.ANAT_DATASETS if anat_datasets is None else anat_datasets
    cfg = _func_config()
    rows = []
    for sp, dirs in func_datasets.items():
        for bd in dirs:
            pat = cfg.get(pdat.bids_name(bd), ("*_bold.nii*", None))[0]
            imgs = [p for p in _images(bd, "func", pat) if "_bold." in os.path.basename(p)]
            print(f"  [acq] {sp:11s} {pdat.bids_name(bd):18s} {len(imgs)} BOLD")
            rows += [_image_row(p, bd, sp, "bold") for p in imgs]
    for sp, dirs in anat_datasets.items():
        for bd in dirs:
            n0 = len(rows)
            for sfx in ANAT_SUFFIXES:
                rows += [_image_row(p, bd, sp, sfx) for p in _images(bd, "anat", f"*_{sfx}.nii*")]
            print(f"  [acq] {sp:11s} {pdat.bids_name(bd):18s} {len(rows) - n0} T1w/T2w")
    return pd.DataFrame(rows)


def load(refresh=False):
    """scan() cached as <CACHE_DIR>/acquisition.pkl."""
    cache = opj(pdat.CACHE_DIR, "acquisition.pkl")
    if not refresh and os.path.isfile(cache):
        print(f"  [cache] {cache}")
        return pd.read_pickle(cache)
    df = scan()
    os.makedirs(pdat.CACHE_DIR, exist_ok=True)
    df.to_pickle(cache)
    return df


# ── summary ────────────────────────────────────────────────────────────────────────────
def _num(v, fmt):
    v = pd.to_numeric(pd.Series(v), errors="coerce").dropna()
    if v.empty:
        return "n.r."
    lo, med, hi = v.min(), v.median(), v.max()
    if np.isclose(lo, hi, rtol=1e-3):
        return format(med, fmt)
    return f"{format(med, fmt)} [{format(lo, fmt)}-{format(hi, fmt)}]"


def _cat(v):
    v = [format(x, ".3g") if isinstance(x, (int, float, np.number)) else str(x)
         for x in pd.Series(v).dropna() if str(x).strip()]
    return ", ".join(dict.fromkeys(v)) if v else "n.r."


def _voxel(g):
    """Most frequent voxel size; '*' when the dataset has several."""
    v = g[["vox_x", "vox_y", "vox_z"]].round(3).dropna()
    if v.empty:
        return "n.r."
    counts = v.value_counts()
    return "x".join(format(x, ".3g") for x in counts.index[0]) + ("*" if len(counts) > 1 else "")


def summary(df):
    """One row per (modality, dataset), in figure order. '*' after the voxel size:
    several voxel sizes in the dataset (the most frequent is given)."""
    if df is None or df.empty:
        return pd.DataFrame()
    cfg = _func_config()
    out = []
    df = df.assign(modality=np.where(df["suffix"].eq("bold"), "fMRI", "anat"))
    for (bold, ds), g in df.groupby(["modality", "dataset"], sort=False):
        sp = g["species"].iloc[0]
        r = dict(modality=bold, species=sp, dataset=ds,
                 n_subjects=g["subject"].nunique(),
                 n_sessions=g[["subject", "session"]].drop_duplicates().shape[0],
                 n_images=len(g),
                 contrast=_cat(g["suffix"]) if bold == "anat" else "BOLD",
                 field_T=_cat(g["field_T"]),
                 scanner=_cat((g["manufacturer"].fillna("") + " " + g["model"].fillna("")).str.strip()
                              .replace("", np.nan)),
                 TE_ms=_num(g["TE_ms"], ".3g"), flip_deg=_num(g["flip_deg"], ".3g"),
                 voxel_mm=_voxel(g), matrix=_cat(g["matrix"].value_counts().index[:1]))
        if bold == "fMRI":
            r.update(state=cfg.get(ds, (None, "unknown"))[1],
                     TR_s=_num(g["TR_s"], ".3g"),
                     n_volumes=_num(g["n_volumes"], ".0f"),
                     duration_min=_num(g["duration_min"], ".3g"),
                     multiband=_cat(g["multiband"]),
                     TR_flagged=int((g["TR_check"] != "").sum()))
        out.append(r)
    s = pd.DataFrame(out)
    order = ps.dataset_sort(list(s["dataset"].unique()), dict(zip(s["dataset"], s["species"])))
    s["_k"] = s["dataset"].map({d: i for i, d in enumerate(order)})
    s["_m"] = s["modality"].map({"fMRI": 0, "anat": 1})
    return s.sort_values(["_m", "_k"]).drop(columns=["_k", "_m"]).reset_index(drop=True)


def lines(s, df):
    out = ["n.r. = not reported in the BIDS sidecars; values = median [min-max] across images"]
    for r in s.itertuples():
        txt = (f"{r.modality} {r.species} {r.dataset}: {r.n_subjects} subjects, "
               f"{r.n_sessions} sessions, {r.n_images} images; field {r.field_T} T; "
               f"{r.scanner}; TE {r.TE_ms} ms; flip {r.flip_deg} deg; voxel {r.voxel_mm} mm")
        if r.modality == "fMRI":
            txt += (f"; TR {r.TR_s} s; {r.n_volumes} volumes; {r.duration_min} min; "
                    f"state {r.state}")
            if r.TR_flagged:
                ex = df[(df["dataset"] == r.dataset) & (df["TR_check"] != "")]["TR_check"].iloc[0]
                txt += f"; TR flagged in {int(r.TR_flagged)} runs (e.g. {ex})"
        out.append(txt)
    return out


# ── figure ─────────────────────────────────────────────────────────────────────────────
FUNC_COLS = [("species", "Species"), ("dataset", "Dataset"), ("state", "State"), ("n", "Subjects /\nsessions / runs"),
             ("field_T", "B0 (T)"), ("scanner", "Scanner"), ("TR_s", "TR (s)"),
             ("TE_ms", "TE (ms)"), ("flip_deg", "Flip (°)"), ("voxel_mm", "Voxel (mm)"),
             ("n_volumes", "Volumes"), ("duration_min", "Duration\n(min)")]
ANAT_COLS = [("species", "Species"), ("dataset", "Dataset"), ("contrast", "Contrast"),
             ("n", "Subjects /\nsessions / images"), ("field_T", "B0 (T)"),
             ("scanner", "Scanner"), ("TE_ms", "TE (ms)"), ("flip_deg", "Flip (°)"),
             ("voxel_mm", "Voxel (mm)"), ("matrix", "Matrix")]
FONT_PT, LINE_MM, GAP_MM = 5.0, 2.45, 2.2


def _tcol(sp):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import supp_panels as fsp
    return fsp._tcol(sp)


def _cells(s, cols):
    """Text of every cell (wrapped), the species of every row."""
    import textwrap
    s = s.copy()
    s["n"] = s.n_subjects.astype(str) + " / " + s.n_sessions.astype(str) + " / " + \
        s.n_images.astype(str)
    species = list(s.species)
    s["dataset"] = [ps.dataset_short(d, sp) for d, sp in zip(s.dataset, s.species)]
    s["species"] = [sp if i == 0 or sp != species[i - 1] else "" for i, sp in enumerate(species)]
    if any(c == "TR_s" for c, _ in cols):
        s["TR_s"] = [t + (" †" if f else "") for t, f in zip(s.TR_s, s.TR_flagged)]
    if "state" in s:
        s["state"] = s["state"].str.replace(" (from ", "\n(from ", regex=False)
    body = [["\n".join(textwrap.fill(p, 24) for p in str(r[c]).split("\n")) for c, _ in cols]
            for _, r in s.iterrows()]
    return body, species


def _measure(fig, texts, bold=False):
    """Width (mm) of each text at FONT_PT."""
    r = fig.canvas.get_renderer()
    out = []
    for t in texts:
        a = fig.text(0, 0, t, fontsize=FONT_PT, fontweight="bold" if bold else "normal")
        out.append(a.get_window_extent(r).width / fig.dpi * 25.4)
        a.remove()
    return out


def _layout(fig, body, cols):
    heads = [h for _, h in cols]
    w = [max(_measure(fig, [h], bold=True) + _measure(fig, [row[j] for row in body]))
         for j, h in enumerate(heads)]
    x = np.r_[0, np.cumsum(np.array(w) + GAP_MM)][:-1]
    head_lines = max(h.count("\n") + 1 for h in heads)
    row_lines = [max(c.count("\n") + 1 for c in row) for row in body]
    return x, float(x[-1] + w[-1]), head_lines, row_lines


def _draw(fig, W, H, x0, y0, s, cols, title, letter):
    """Draw one table, its top-left corner at (x0, y0) mm from the top-left of the figure;
    returns its height (mm)."""
    body, species = _cells(s, cols)
    x, width, hl, rl = _layout(fig, body, cols)
    X = lambda mm: (x0 + mm) / W
    Y = lambda mm: 1 - (y0 + mm) / H
    fig.text((x0 - 1.5) / W, Y(0), letter, ha="right", va="top", fontsize=ps.PANEL_LABEL_SIZE,
             fontweight="bold")
    fig.text(X(0), Y(0.3), title, ha="left", va="top", fontsize=6, fontweight="bold")
    y = 4.2
    fig.add_artist(_hline(X(0), X(width), Y(y), 0.5))
    y += 0.6
    for j, (_, h) in enumerate(cols):
        fig.text(X(x[j]), Y(y), h, ha="left", va="top", fontsize=FONT_PT, fontweight="bold",
                 linespacing=1.15)
    y += hl * LINE_MM + 0.4
    fig.add_artist(_hline(X(0), X(width), Y(y), 0.3))
    y += 0.6
    for i, row in enumerate(body):
        if i and species[i] != species[i - 1]:
            fig.add_artist(_hline(X(0), X(width), Y(y - 0.3), 0.15, color="0.75"))
        for j, cell in enumerate(row):
            fig.text(X(x[j]), Y(y), cell, ha="left", va="top", fontsize=FONT_PT,
                     linespacing=1.15, color=_tcol(species[i]) if cols[j][0] in ("species", "dataset")
                     else "black",
                     fontweight="bold" if cols[j][0] in ("species", "dataset") else "normal")
        y += rl[i] * LINE_MM + 0.6
    fig.add_artist(_hline(X(0), X(width), Y(y - 0.2), 0.5))
    return y


def _hline(x0, x1, y, lw, color="black"):
    from matplotlib.lines import Line2D
    return Line2D([x0, x1], [y, y], lw=lw, color=color, solid_capstyle="butt")


def _height(fig, s, cols):
    body, _ = _cells(s, cols)
    _, width, hl, rl = _layout(fig, body, cols)
    return 4.2 + 0.6 + hl * LINE_MM + 0.4 + 0.6 + sum(r * LINE_MM + 0.6 for r in rl), width


def plot_table(s, out_base, width_mm=183):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if s is None or s.empty:
        return None
    parts = [(s[s.modality == "fMRI"], FUNC_COLS,
              "Resting-state fMRI (Figs 4-5, Supplementary Figs 5, 7-11)"),
             (s[s.modality == "anat"], ANAT_COLS,
              "Anatomical MRI (Figs 2-3, Supplementary Figs 4, 6)")]
    parts = [p for p in parts if len(p[0])]
    note = ("Values read from the BIDS JSON sidecars and NIfTI headers: median [min-max] "
            "across images when they differ; n.r., not reported in the sidecars; *, several "
            "voxel sizes in the dataset (most frequent shown); †, sidecar TR implausible or "
            "inconsistent with the NIfTI header (pixdim[4]), header value shown. State: "
            "arousal state during the resting-state scans.")
    with plt.rc_context(ps.PAPER_RC):
        probe = plt.figure(figsize=(width_mm * MM, 1))
        dims = [_height(probe, p[0], p[1]) for p in parts]
        plt.close(probe)
        x0 = 4.0
        H = sum(h for h, _ in dims) + 5.0 * (len(parts) - 1) + 12.0
        fig = plt.figure(figsize=(width_mm * MM, H * MM))
        y = 1.0
        for (sub, cols, title), letter in zip(parts, "ab"):
            y += _draw(fig, width_mm, H, x0, y, sub, cols, title, letter) + 5.0
        too_wide = [w for _, w in dims if w > width_mm - x0 - 1]
        if too_wide:
            print(f"  [WARN] table wider than the page: {max(too_wide):.0f} mm")
        fig.text(x0 / width_mm, 1 - (y - 3.0) / H, note, ha="left", va="top", fontsize=4.5,
                 wrap=True)
        paths = ps.save(fig, out_base)
        plt.close(fig)
    return paths
