#!/usr/bin/env python3
"""
EDNiX — Paper Figure Orchestrator
=================================
ONE entry point that produces every figure of the EDNiX paper, organized
strictly by the figure numbering of the manuscript. It calls the project's
real functions (no re-implementation):

  Fig 1  — (not script-generated: pipeline / architecture overview)
  Fig 2  — Morphometry & allometric scaling
             2A  morphometry by cortical division          → make_all_figures
             2B  PGLS body weight + neuron number          → pgls_connect
  Fig 3  — Developmental & aging trajectories
             3A  macaque surface trajectory                → longitudinal
             3B  macaque thickness (non-monotonic)         → longitudinal
             3C  association vs primary growth (lvl3)       → longitudinal (lvl3)
             3D  human aging (BIDS_Park)                    → longitudinal
  Fig 4  — Functional network quality
             4A  homotopic-specificity threshold method    → threshold_explorer
             4B  category proportions by species/anaesth   → threshold_explorer
  Fig 5  — Seed-based connectivity + FC fingerprint
             5   retrosplenial seed (niftotoWBsurface PNG) → mounted here
             5c  cross-species fingerprint similarity      → pgls_connect
  Fig 6  — FDG-PET glucose metabolism (mounted when ready)

Each figure is a method `fig2()`, `fig3()`, … so you can run one at a time:
    python ednix_paper_figures.py --only 3
or all:
    python ednix_paper_figures.py

CONFIG — everything dataset/path-specific is in the CONFIG block below.
Nothing is left blank for you to fill except genuinely missing data
(PET renders, seed PNGs) which are clearly flagged.
"""

import os
import sys
import argparse
import warnings
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

sys.path.insert(0, "/home/cgarin/PycharmProjects/EDNiX")
sys.path.insert(0, "/home/cgarin/PycharmProjects/EDNiX/Exemples/Study_EDNiX/analysis")

opj = os.path.join

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG — single source of truth
# ═══════════════════════════════════════════════════════════════════════════════

ATLAS_LIB      = "/home/cgarin/PycharmProjects/EDNiX/Atlases_library"
ATLAS_NAME     = "EDNIxCSC"
NEWICK_FILE    = "/scratch2/EDNiX3/Liste des espèces.nwk"
HERCULANO_XLSX = "/scratch2/EDNiX3/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx"

# results roots (everything is in EDNiX3, NOT EDNiX)
RESULTS_DIR    = "/scratch2/EDNiX3/results/multispecies_analysis"
CSV_DIR        = opj(RESULTS_DIR, "csv")
COMBO_DIR      = opj(RESULTS_DIR, "figures", "cross_species", "combo")   # make_all_figures output
FC_RESULTS_DIR = "/scratch2/EDNiX3/results/cross_species/qc_recap"
ARTICLE_DIR    = "/scratch2/EDNiX3/results/article_figures"

# species → BIDS  (single + multi)
SPECIES_BIDS = {
    "Rat":        "/scratch2/EDNiX3/Rat/BIDS_Grandjean",
    "Mouse":      "/scratch2/EDNiX3/Mouse/BIDS_Grandjean2",
    "Dog":        "/scratch2/EDNiX3/Dog/BIDS_Boch_K9",
    "Marmoset":   "/scratch2/EDNiX3/Marmoset/BIDS_Tian",
    "Mouselemur": "/scratch2/EDNiX3/Mouselemur/BIDS_Garin",
}
SPECIES_MULTI_BIDS = {
    "Macaque": ["/scratch2/EDNiX3/Macaque/BIDS_BenHamed",
                "/scratch2/EDNiX3/Macaque/BIDS_Zhu_Garin"],
    "Human":   ["/scratch2/EDNiX3/Human/BIDS_Merida",
                "/scratch2/EDNiX3/Human/BIDS_Park",
                "/scratch2/EDNiX3/Human/BIDS_Castrillon"],
}

# Longitudinal datasets (Fig 3)
MACAQUE_LONGI_BIDS  = ["/scratch2/EDNiX3/Macaque/BIDS_Zhu_Garin"]
MACAQUE_AGE_XLSX    = "/scratch2/EDNiX3/Macaque/BIDS_Zhu_Garin/Garin_macaque.xlsx"
HUMAN_AGING_BIDS    = ["/scratch2/EDNiX3/Human/BIDS_Park"]
HUMAN_AGE_PREFIX    = "AgeMRI_W"

# the 15 cortical ROIs used for PGLS / fingerprint / lvl3 development (Fig 2B,3C,5c)
REGIONS_LVL3 = [
    'Auditory cortex (Superior temporal )', 'Insula and others in lateral sulcus',
    'Middle Temporal, Inferior temporal , Temporal pole  (MIPT)',
    'Motor and premotor', 'Olfactory cortex', 'Orbital PFC (oPFC)',
    'Orbital frontal cortex (oFC)', 'Periarchicortex',
    'Posterior medial cortex (PMC)', 'Posterior parietal cortex',
    'Prefrontal cortex (PFC)', 'Somatosensory cortex',
    'Ventral areas of the temporal lobe (vent Temp)',
    'Visual pre and extra striate cortex', 'Visual striate cortex']

# Primary vs association split for Fig 3D (paper: assoc. grows longer than primary)
PRIMARY_REGIONS = ['Somatosensory cortex', 'Visual striate cortex',
                   'Auditory cortex (Superior temporal )', 'Motor and premotor']
ASSOCIATION_REGIONS = ['Prefrontal cortex (PFC)', 'Orbital PFC (oPFC)',
                       'Posterior parietal cortex', 'Posterior medial cortex (PMC)',
                       'Middle Temporal, Inferior temporal , Temporal pole  (MIPT)']

# FC thresholds — None = auto (suggest_thresholds 25th/33rd pct). Fill after first run.
THRESH_INTRA = None
THRESH_DELTA = None

# Fig 4 layout: set False to drop panel D (network metrics) — "D is not interesting"
FIG4_INCLUDE_D = True

# Fig 5/6 external renders (PNG produced by _3dLMEr → niftotoWBsurface, and PET)
SEED_PNG = {
    # "Macaque": "/scratch2/EDNiX/.../Macaque_retrosplenial.png", ...
}
PET_PNG = {
    # "Human":   "/scratch2/EDNiX/.../Human_FDG.png",
    # "Macaque": "/scratch2/EDNiX/.../Macaque_FDG.png",
}


# ═══════════════════════════════════════════════════════════════════════════════
# SHARED HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def _dirs():
    # mounted figures live at the article root; only these subdirs are needed
    for sub in ("panels", "stats", "supplementary"):
        os.makedirs(opj(ARTICLE_DIR, sub), exist_ok=True)

# Master accumulator for every reported statistic (R, p, CI, slope, …)
_PAPER_VALUES = []

def _log_value(figure, measure, **kw):
    """Append one reported statistic to the master paper_values.csv."""
    row = dict(figure=figure, measure=measure)
    row.update(kw)
    _PAPER_VALUES.append(row)
    # also echo to console so you see it live
    kvs = "  ".join(f"{k}={v}" for k, v in kw.items())
    print(f"  [VALUE] {figure} | {measure} | {kvs}")

def _flush_paper_values():
    if not _PAPER_VALUES:
        return
    out = opj(ARTICLE_DIR, "stats", "paper_values.csv")
    pd.DataFrame(_PAPER_VALUES).to_csv(out, index=False)
    print(f"\n  [PAPER VALUES] {len(_PAPER_VALUES)} statistics → {out}")


# Module-level species ordering / colours (used by combo + supplementary)
_PHYLO_ORDER = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
                "Macaque", "Human", "Dog", "Cat", "Pig"]
_SPECIES_COLORS = {
    "Mouse": "#1B4F72", "Rat": "#2980B9", "Mouselemur": "#1D6B3F",
    "Marmoset": "#27AE60", "Macaque": "#82E0AA", "Human": "#196F3D",
    "Dog": "#E67E22", "Cat": "#F39C12", "Bat": "#8E44AD", "Pig": "#A04000"}


def _combo_plot(data_dict, regions, atlas_level, hemisphere, out_path):
    """
    3 rows (surface / volume / thickness) × N regions grid.
    One dot per subject coloured by species + median bar, log-y, no title.
    Used for Fig 2A (lvl1 Iso/Allo/Periallo) and Suppl 4-6 (other lvls).
    """
    import numpy as np
    dfs = {}
    for key in ("surface", "volume", "thickness"):
        d = data_dict.get(key)
        if d is None or d.empty:
            continue
        if "atlas_level" in d.columns:
            d = d[d["atlas_level"] == atlas_level]
        if "hemisphere" in d.columns and hemisphere != "all":
            h = d[d["hemisphere"] == hemisphere]
            d = h if not h.empty else d[d["hemisphere"] == "bilateral"]
        if not d.empty:
            dfs[key] = d
    if not dfs:
        print(f"  [combo] no data at lvl{atlas_level} ({hemisphere}) → skipped")
        return None

    present = set()
    for d in dfs.values():
        present |= set(d["region"].unique())
    if regions:
        regs = [r for r in regions if r in present]
    else:
        regs = sorted(present)
    if not regs:
        print(f"  [combo lvl{atlas_level}] none of the requested regions present")
        return None

    species = set()
    for d in dfs.values():
        species |= set(d["species"].unique())
    sp_list = [s for s in _PHYLO_ORDER if s in species] + \
              [s for s in species if s not in _PHYLO_ORDER]

    n_reg = len(regs)
    fig, axes = plt.subplots(3, n_reg, figsize=(max(2.4 * n_reg, 5), 8),
                             squeeze=False)
    rng = np.random.default_rng(0)
    mods = [("surface", "surface_area_mm2", "Surface (mm²)"),
            ("volume",  "volume_mm3",        "Volume (mm³)"),
            ("thickness", "thickness_mm",    "Thickness (mm)")]
    for ri, (mod_key, col, mod_label) in enumerate(mods):
        d = dfs.get(mod_key)
        for ci, region in enumerate(regs):
            ax = axes[ri][ci]
            if d is None or d.empty:
                ax.axis("off"); continue
            sub = d[d["region"] == region]
            for si, sp in enumerate(sp_list):
                vals = sub.loc[sub["species"] == sp, col].dropna().values
                vals = vals[np.isfinite(vals) & (vals > 0)]
                if len(vals) == 0:
                    continue
                color = _SPECIES_COLORS.get(sp, "#888888")
                jit = rng.uniform(-0.22, 0.22, len(vals))
                ax.scatter(np.full(len(vals), si) + jit, vals, s=10, alpha=0.5,
                           color=color, edgecolors="none", zorder=3)
                med = float(np.median(vals))
                ax.plot([si - 0.3, si + 0.3], [med, med],
                        color=color, lw=2.2, zorder=4)
            ax.set_yscale("log")
            ax.set_xticks(range(len(sp_list)))
            if ri == 2:
                ax.set_xticklabels(sp_list, rotation=55, ha="right", fontsize=7)
            else:
                ax.set_xticklabels([])
            if ci == 0:
                ax.set_ylabel(mod_label, fontsize=9)
            if ri == 0:
                short = region.split("(")[0].strip()[:22]
                ax.set_title(short, fontsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=7)

    plt.tight_layout(pad=0.4, h_pad=0.2, w_pad=0.2)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [combo] {out_path}")
    return out_path


def _letterbox(img, target_ar):
    """Pad an RGBA/RGB image array with white to reach `target_ar` (= W/H),
    keeping the original centred. This makes every panel occupy an identical
    box so the mount renders them all the SAME visible size, centred — no
    distortion (unlike aspect='auto') and no size mismatch (unlike raw imshow).
    """
    h, w = img.shape[:2]
    ar = w / h
    if img.ndim == 2:                      # grayscale → RGB
        img = np.dstack([img] * 3)
    if img.shape[2] == 3:                  # add opaque alpha
        img = np.dstack([img, np.ones((h, w), img.dtype)])
    fill = 1.0 if img.dtype.kind == "f" else 255
    if ar < target_ar:                     # too tall → pad left/right
        new_w = int(round(h * target_ar))
        pad = (new_w - w) // 2
        out = np.full((h, new_w, img.shape[2]), fill, dtype=img.dtype)
        out[:, pad:pad + w] = img
    else:                                  # too wide → pad top/bottom
        new_h = int(round(w / target_ar))
        pad = (new_h - h) // 2
        out = np.full((new_h, w, img.shape[2]), fill, dtype=img.dtype)
        out[pad:pad + h, :] = img
    return out


def _mount(panel_paths, layout, out_name, labels=None, title=None, figsize=None):
    """Assemble PNG panels at a UNIFORM size, each centred in its cell, with a
    centred trailing row when the panel count doesn't fill the grid."""
    paths = [p for p in panel_paths if p and os.path.exists(p)]
    if not paths:
        print(f"  [mount] no panels for {out_name}")
        return None
    nrows, ncols = layout
    imgs = [mpimg.imread(p) for p in paths[:nrows * ncols]]
    # common aspect = widest panel, so nothing is cropped; all letterboxed to it
    target_ar = max(im.shape[1] / im.shape[0] for im in imgs)
    imgs = [_letterbox(im, target_ar) for im in imgs]

    if figsize is None:
        cell_w = 5.6
        figsize = (ncols * cell_w, nrows * cell_w / target_ar)
    fig = plt.figure(figsize=figsize)

    n = len(imgs)
    full_rows = n // ncols
    margin = 0.012                          # uniform gutter (figure fraction)
    cw, ch = 1.0 / ncols, 1.0 / nrows
    for i, im in enumerate(imgs):
        r, c = divmod(i, ncols)
        row_n = ncols if r < full_rows else (n - full_rows * ncols)
        x_off = (1.0 - row_n * cw) / 2.0    # centre a partial trailing row
        x0 = x_off + c * cw + margin
        y0 = 1.0 - (r + 1) * ch + margin
        ax = fig.add_axes([x0, y0, cw - 2 * margin, ch - 2 * margin])
        ax.imshow(im); ax.axis("off")
        if labels and i < len(labels) and labels[i]:
            ax.text(0.0, 1.0, labels[i], transform=ax.transAxes,
                    fontsize=18, fontweight="bold", va="top", ha="left")
    out = opj(ARTICLE_DIR, out_name)        # ← article root (matches output tree)
    fig.savefig(out, dpi=300, facecolor="white", pad_inches=0.02)
    plt.close(fig)
    print(f"  [MOUNTED] {out}")
    return out

def _build_species_config():
    cfg, label_paths, frags = {}, {}, {}
    from ednix_bids_tools import get_atlas_label_path, find_species_path
    for sp, bd in SPECIES_BIDS.items():
        cfg[sp] = {"bids_dirs": [bd], "list_to_keep": [], "list_to_remove": []}
    for sp, bds in SPECIES_MULTI_BIDS.items():
        cfg[sp] = {"bids_dirs": bds, "list_to_keep": [], "list_to_remove": []}
    for sp in cfg:
        try:
            frag = find_species_path(ATLAS_LIB, sp)
            frags[sp] = frag
            label_paths[sp] = get_atlas_label_path(ATLAS_LIB, frag, ATLAS_NAME,
                                                   prefer_statslut=False)
        except Exception as e:
            print(f"  [WARN] atlas {sp}: {e}")
    return cfg, label_paths, frags


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 2 — MORPHOMETRY & ALLOMETRIC SCALING
# ═══════════════════════════════════════════════════════════════════════════════

def fig2(data_plot=None):
    """
    2A morphometry by division (lvl1)  → from make_all_figures output (reuse)
    2B PGLS body weight + neuron number → pgls_connect machinery
    """
    print("\n" + "=" * 64 + "\n  FIG 2 — Morphometry & allometric scaling\n" + "=" * 64)
    from Statistics.Evo.PGLS import load_newick_vcv, run_pgls_all_region_pairs
    from pgls_connect import (load_herculano, run_morphometry_scaling,
                              _build_species_config as pgls_species_config)
    from ednix_bids_tools import collect_multi_species

    species_config, atlas_label_paths, atlas_frags = _build_species_config()

    # ── 2A: combo plot at level 1 (Iso/Allo/Periallo, bilateral, log y) ──────
    # If make_all_figures already produced combo_bilateral_log.png, REUSE it.
    existing_2a = opj(COMBO_DIR, "combo_bilateral_log.png")
    if os.path.exists(existing_2a):
        p2a = existing_2a
        print(f"  [Fig2A] reusing existing: {existing_2a}")
    else:
        data_l1 = collect_multi_species(
            species_config, regions_of_interest=["Isocortex", "Allocortex", "Periallocortex"],
            extract=("surface", "volume", "thickness"), atlas_name=ATLAS_NAME,
            atlas_label_paths=atlas_label_paths, atlas_library_root=ATLAS_LIB,
            species_atlas_fragments=atlas_frags)
        # ── Render 2A with the ORIGINAL plot_combo (grey violins + anesthesia
        # scatter, PAPER_RC) — "strictly what was done before". Only fall back
        # to the internal _combo_plot if the original cannot be imported/run.
        def _as_df(v):
            if isinstance(v, list):
                v = [x for x in v if x is not None and not getattr(x, "empty", True)]
                return pd.concat(v, ignore_index=True) if v else pd.DataFrame()
            return v if v is not None else pd.DataFrame()
        p2a = opj(ARTICLE_DIR, "panels", "Fig2A_combo_bilateral_log.png")
        try:
            from EDNiX_figures import plot_combo
            plot_combo(_as_df(data_l1.get("surface")),
                       _as_df(data_l1.get("volume")),
                       _as_df(data_l1.get("thickness")),
                       p2a, regions=["Isocortex", "Allocortex", "Periallocortex"],
                       hemisphere="bilateral", log_scale=True, atlas_level=1)
            print(f"  [Fig2A] rendered with original plot_combo (violin style)")
        except Exception as e:
            print(f"  [Fig2A] plot_combo unavailable ({e}); using _combo_plot fallback")
            p2a = _combo_plot(
                data_l1, regions=["Isocortex", "Allocortex", "Periallocortex"],
                atlas_level=1, hemisphere="bilateral", out_path=p2a)

    # ── 2B: PGLS scaling at level 3 (15 cortical ROIs) ───────────────────────
    data = collect_multi_species(
        species_config, regions_of_interest=REGIONS_LVL3,
        extract=("surface", "volume", "thickness"), atlas_name=ATLAS_NAME,
        atlas_label_paths=atlas_label_paths, atlas_library_root=ATLAS_LIB,
        species_atlas_fragments=atlas_frags)

    vcv_df, _ = load_newick_vcv(NEWICK_FILE)
    hh = load_herculano(HERCULANO_XLSX)

    panels = []
    all_pairs = []
    scaling_dir = opj(ARTICLE_DIR, "panels", "fig2_scaling")
    for key, col, label in [("surface", "surface_area_mm2", "Surface area (mm²)"),
                            ("volume", "volume_mm3", "Volume (mm³)"),
                            ("thickness", "thickness_mm", "Thickness (mm)")]:
        df = data.get(key)
        if df is None or df.empty:
            continue
        # 2B: allometric scaling (body/brain/neurons) — pgls_connect's function
        run_morphometry_scaling(df, col, label, key, vcv_df, hh,
                                out_dir=opj(scaling_dir, key), atlas_level=3)
        comp = opj(scaling_dir, key, f"scaling_{key}_composite.png")
        if os.path.exists(comp):
            panels.append(comp)
        # stats: all region-pairs PGLS
        dfr = run_pgls_all_region_pairs(
            df_morph=df, metric_col=col, vcv_df=vcv_df, regions=REGIONS_LVL3,
            atlas_level=3, hemisphere="bilateral", min_species=4,
            log_x=True, log_y=True)
        if dfr is not None and not dfr.empty:
            dfr["metric"] = key
            all_pairs.append(dfr)
            # log the headline scaling stats (slope, p, R2) for the paper
            for _, r in dfr.iterrows():
                if r.get("significant_fdr05", False):
                    _log_value("Fig2B", f"{key}_{r.get('region_x','?')}_vs_{r.get('region_y','?')}",
                               slope=round(float(r.get("slope", float("nan"))), 3),
                               p=round(float(r.get("p_value", r.get("pval", float("nan")))), 4),
                               R2=round(float(r.get("r_squared", r.get("r2", float("nan")))), 3))

    if all_pairs:
        out = opj(ARTICLE_DIR, "stats", "fig2_pgls_region_pairs.csv")
        pd.concat(all_pairs, ignore_index=True).to_csv(out, index=False)
        print(f"  [stats] {out}")

    # Final Fig 2 montage: 2A (combo lvl1) on top, then 2B scaling panels below
    final = [p2a] + panels if p2a else panels
    _mount(final, layout=(len(final), 1), out_name="Fig2_morphometry.png",
           figsize=(10, len(final) * 5),
           labels=["A"] + ["B"] * max(len(panels), 0))


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 3 — DEVELOPMENT & AGING  (lvl1 for A/B/D, lvl3 for C)
# ═══════════════════════════════════════════════════════════════════════════════

def fig3():
    """
    Main Fig 3 = 2×3 grid:
        top row    Macaque (volume, surface, thickness)  — development
        bottom row Human   (volume, surface, thickness)  — aging
    Fig 3D (separate) = association vs primary cortex growth, lvl3, macaque.
    Plus supplementary: every region trajectory + aging stats CSV.
    """
    print("\n" + "=" * 64 + "\n  FIG 3 — Development & aging\n" + "=" * 64)
    from longitudinal import (load_morphometry, format_macaque_age,
                              format_human_age, _build_atlas_paths)
    from Statistics.Group_fMRI.EDNiX_longitudinal import run_longitudinal

    all_stats = []
    MODS = [("volume", "volume_mm3", "Volume (mm³)"),
            ("surface", "surface_area_mm2", "Surface area (mm²)"),
            ("thickness", "thickness_mm", "Cortical thickness (mm)")]
    REGS_L1 = ["Isocortex", "Allocortex", "Periallocortex"]

    # ── Macaque (development) ────────────────────────────────────────────────
    mac_paths, _ = _build_atlas_paths(["Macaque"])
    mac_morph = load_morphometry({"Macaque": MACAQUE_LONGI_BIDS}, mac_paths)
    mac_age = format_macaque_age(MACAQUE_AGE_XLSX)
    mac_panels = {}
    for mod, col, label in MODS:
        df = mac_morph.get(mod, pd.DataFrame())
        if df is None or df.empty:
            continue
        out = opj(ARTICLE_DIR, "panels", "Macaque", mod)
        stats = run_longitudinal(
            label=f"Macaque_{mod}", df_morph=df, df_age=mac_age,
            metric_col=col, metric_label=label,
            age_col="age_years", age_label="Age (years)",
            out_dir=out, atlas_level=1, hemisphere="left",
            regions=REGS_L1, min_sessions=3)
        if stats is not None and not stats.empty:
            all_stats.append(stats)
        traj = opj(out, "Isocortex_trajectory.png")   # isocortex = main panel
        if os.path.exists(traj):
            mac_panels[mod] = traj

    # ── Human (aging) ────────────────────────────────────────────────────────
    hum_paths, _ = _build_atlas_paths(["Human"])
    hum_morph = load_morphometry({"Human": HUMAN_AGING_BIDS}, hum_paths)
    hum_age = format_human_age(HUMAN_AGING_BIDS, age_prefix=HUMAN_AGE_PREFIX)
    hum_panels = {}
    for mod, col, label in MODS:
        df = hum_morph.get(mod, pd.DataFrame())
        if df is None or df.empty:
            continue
        out = opj(ARTICLE_DIR, "panels", "Human", mod)
        stats = run_longitudinal(
            label=f"Human_{mod}", df_morph=df, df_age=hum_age,
            metric_col=col, metric_label=label,
            age_col="age", age_label="Age (years)",
            out_dir=out, atlas_level=1, hemisphere="left",
            regions=REGS_L1, min_sessions=3)
        if stats is not None and not stats.empty:
            all_stats.append(stats)
        traj = opj(out, "Isocortex_trajectory.png")
        if os.path.exists(traj):
            hum_panels[mod] = traj

    # ── 2×3 mount: row1 macaque, row2 human (vol, surf, thick) ───────────────
    order = ["volume", "surface", "thickness"]
    grid_panels = ([mac_panels.get(m) for m in order]
                   + [hum_panels.get(m) for m in order])
    _mount(grid_panels, layout=(2, 3),
           out_name="Fig3_development_aging.png",
           labels=["A", "B", "C", "D", "E", "F"],
           figsize=(18, 10),
           title="Figure 3 | Development (macaque, top) and aging (human, bottom)")

    # ── Fig 3D: association vs primary growth, lvl3 ──────────────────────────
    mac_morph_l3 = load_morphometry({"Macaque": MACAQUE_LONGI_BIDS}, mac_paths)
    _fig3d_assoc_vs_primary(mac_morph_l3, mac_age,
                            opj(ARTICLE_DIR, "panels", "Fig3D_assoc_vs_primary.png"),
                            all_stats)

    # ── Human aging quantification (volume reduction, surf/thick preserved) ──
    _fig3_human_aging_stats(hum_morph, hum_age, all_stats)

    # ── stats summary ────────────────────────────────────────────────────────
    if all_stats:
        out = opj(ARTICLE_DIR, "stats", "fig3_longitudinal_summary.csv")
        pd.concat([s for s in all_stats if not s.empty],
                  ignore_index=True).to_csv(out, index=False)
        print(f"  [stats] {out}")


def _fig3_human_aging_stats(hum_morph, hum_age, all_stats):
    """
    Quantify human aging: regression of each modality on age (Iso/Allo/Periallo).
    Paper claim: volume reduces, surface & thickness relatively preserved.
    Writes a dedicated CSV with slope, p, % change per decade.
    """
    import numpy as np
    from longitudinal import _filter_df
    rows = []
    for mod, col in [("volume", "volume_mm3"), ("surface", "surface_area_mm2"),
                     ("thickness", "thickness_mm")]:
        df = hum_morph.get(mod, pd.DataFrame())
        if df is None or df.empty:
            continue
        df = _filter_df(df, atlas_level=1, hemisphere="left")
        df["subject"] = df["subject"].astype(str)
        age = hum_age.copy(); age["subject"] = age["subject"].astype(str)
        m = df.merge(age, on="subject", how="inner")
        for region in ["Isocortex", "Allocortex", "Periallocortex"]:
            rdf = m[m["region"] == region][["age", col]].dropna()
            if len(rdf) < 5:
                continue
            sl, ic, r, p, se = _stats_linregress(rdf["age"].values, rdf[col].values)
            mean_val = rdf[col].mean()
            pct_decade = 100 * sl * 10 / mean_val if mean_val else np.nan
            rows.append(dict(analysis="human_aging_regression", modality=mod,
                             region=region, n=len(rdf), slope=sl, p=p, r=r,
                             pct_change_per_decade=pct_decade))
            _log_value("Fig3_human_aging", f"{mod}_{region}",
                       n=len(rdf), R=round(r, 3), p=round(p, 4),
                       slope=round(sl, 4),
                       pct_per_decade=round(pct_decade, 2))
    if rows:
        out = opj(ARTICLE_DIR, "stats", "fig3_human_aging_regression.csv")
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f"  [stats] {out}")
        all_stats.append(pd.DataFrame(rows))


def _stats_linregress(x, y):
    from scipy import stats as _st
    r = _st.linregress(x, y)
    return r.slope, r.intercept, r.rvalue, r.pvalue, r.stderr


def _fig3d_assoc_vs_primary(mac_morph, mac_age, out_path, all_stats):
    """
    Fig 3D: mean % growth of association vs primary cortex (lvl3) over development.
    Uses surface area. Robust region matching (substring, case-insensitive) so
    it works even if the atlas region names differ slightly from the config.
    """
    import numpy as np
    from longitudinal import _filter_df
    df = mac_morph.get("surface", pd.DataFrame())
    if df is None or df.empty:
        print("  [Fig3D] no macaque surface data at all")
        return None

    df = _filter_df(df, atlas_level=3, hemisphere="left")
    if df.empty:
        print("  [Fig3D] no lvl3 data after filter — check atlas_level=3 exists")
        return None

    # Diagnostic: show what region names actually exist
    available = sorted(df["region"].unique())
    print(f"  [Fig3D] {len(available)} lvl3 regions available, e.g.: {available[:5]}")

    def _match(region_list):
        """Match config regions to actual names by case-insensitive substring."""
        matched = []
        for want in region_list:
            wl = want.lower().strip()
            for avail in available:
                al = avail.lower().strip()
                # match if either contains the other's key words
                if wl in al or al in wl or _key_overlap(wl, al):
                    matched.append(avail)
                    break
        return list(dict.fromkeys(matched))   # dedupe, keep order

    def _key_overlap(a, b):
        # match on distinctive keyword (e.g. 'somatosensory', 'prefrontal')
        keys = ['somatosensory', 'visual striate', 'auditory', 'motor',
                'prefrontal', 'orbital', 'parietal', 'posterior medial',
                'temporal']
        for k in keys:
            if k in a and k in b:
                return True
        return False

    prim_regions  = _match(PRIMARY_REGIONS)
    assoc_regions = _match(ASSOCIATION_REGIONS)
    print(f"  [Fig3D] matched PRIMARY: {prim_regions}")
    print(f"  [Fig3D] matched ASSOCIATION: {assoc_regions}")

    if not prim_regions or not assoc_regions:
        print("  [Fig3D] region matching failed — check PRIMARY/ASSOCIATION_REGIONS "
              "against the available names printed above")
        return None

    df["subject"] = df["subject"].astype(str)
    mac_age = mac_age.copy(); mac_age["subject"] = mac_age["subject"].astype(str)
    if "session" in df.columns and "session" in mac_age.columns:
        df["session"] = df["session"].astype(str)
        mac_age["session"] = mac_age["session"].astype(str)
        merged = df.merge(mac_age, on=["subject", "session"], how="inner")
    else:
        merged = df.merge(mac_age, on="subject", how="inner")
    if merged.empty:
        print("  [Fig3D] merge with age table empty"); return None

    def _grp_traj(regions):
        sub = merged[merged["region"].isin(regions)]
        g = sub.groupby(["subject", "age_years"])["surface_area_mm2"].sum().reset_index()
        rows = []
        for s, gs in g.groupby("subject"):
            gs = gs.sort_values("age_years")
            base = gs["surface_area_mm2"].iloc[0]
            if base <= 0 or len(gs) < 2:
                continue
            gs = gs.copy()
            gs["pct"] = 100 * (gs["surface_area_mm2"] - base) / base
            rows.append(gs)
        return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    assoc = _grp_traj(assoc_regions)
    prim  = _grp_traj(prim_regions)
    if assoc.empty or prim.empty:
        print(f"  [Fig3D] insufficient trajectory data "
              f"(assoc={len(assoc)}, prim={len(prim)} rows)")
        return None

    fig, ax = plt.subplots(figsize=(7, 5))
    for traj, color, lbl in [(assoc, "#196F3D", "Association"),
                             (prim, "#E67E22", "Primary")]:
        for s, gs in traj.groupby("subject"):
            ax.plot(gs["age_years"], gs["pct"], "-", color=color, alpha=0.25, lw=0.8)
        grid = np.linspace(traj["age_years"].min(), traj["age_years"].max(), 50)
        mat = []
        for s, gs in traj.groupby("subject"):
            gs = gs.sort_values("age_years")
            if len(gs) >= 2:
                mat.append(np.interp(grid, gs["age_years"], gs["pct"],
                                     left=np.nan, right=np.nan))
        if mat:
            mat = np.array(mat)
            mean = np.nanmean(mat, 0)
            sem = np.nanstd(mat, 0) / np.sqrt(np.maximum(np.sum(~np.isnan(mat), 0), 1))
            ax.plot(grid, mean, color=color, lw=2.6, label=lbl, zorder=5)
            ax.fill_between(grid, mean - sem, mean + sem, color=color, alpha=0.18)

    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("Age (years)", fontsize=12)
    ax.set_ylabel("% surface-area change from baseline", fontsize=12)
    ax.legend(fontsize=11, frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [panel] {out_path}")

    # stats: final % + t-test assoc vs primary
    assoc_final = assoc.groupby("subject")["pct"].last()
    prim_final  = prim.groupby("subject")["pct"].last()
    from scipy import stats as _st
    if len(assoc_final) >= 2 and len(prim_final) >= 2:
        t, p = _st.ttest_ind(assoc_final, prim_final, equal_var=False)
    else:
        t, p = float("nan"), float("nan")
    _log_value("Fig3D", "assoc_vs_primary_final_growth",
               assoc_pct=round(float(assoc_final.mean()), 2),
               primary_pct=round(float(prim_final.mean()), 2),
               t=round(float(t), 3), p=round(float(p), 4),
               n_assoc=len(assoc_final), n_primary=len(prim_final))
    return out_path


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 4 — FUNCTIONAL NETWORK QUALITY
# ═══════════════════════════════════════════════════════════════════════════════

def fig4():
    """
    Fig 4 — Functional network quality (5 panels):
      4A  homotopic-specificity threshold method (threshold_explorer fig 01)
      4B  category proportions by species        (fig 03 species)
      4C  mouse anaesthesia regimes              (fig 03 by-bids, mouse subset)
      4D  network metrics: eigenvalue ratio / silhouette by category
      4E  surrogate / perturbation sensitivity   (ednix_fc_perturbation)
    Writes 00_per_run_summary.csv (consumed elsewhere) + logs values.
    """
    print("\n" + "=" * 64 + "\n  FIG 4 — Functional network quality\n" + "=" * 64)
    import numpy as np
    from scipy import stats as _st
    from ednix_threshold_explorer import (
        run_pipeline_multilevel, classify_runs, suggest_thresholds)

    qc_csv = opj(CSV_DIR, "qc.csv")
    if not os.path.exists(qc_csv):
        print(f"  [Fig4] {qc_csv} not found — run multispecies_analysis first")
        return
    df_qc = pd.read_csv(qc_csv)

    # Mixed atlas variants: humans use 'EDNIxCSC' (no LR), macaques use 'EDNIxCSCLR'.
    # Try both, keep whichever loads runs successfully.
    res = None
    last_err = None
    for use_lr_try in (True, False):
        try:
            res = run_pipeline_multilevel(
                df_qc, fit_kind="correlation",
                bids_root_template="/scratch2/EDNiX3/{species}/{bids_dir}",
                atlas_levels=(2, 3, 4),
                subsets=("all", "primates", "primates_rodents"),
                thresh_intra=THRESH_INTRA, thresh_delta=THRESH_DELTA, stringency=0.0,
                n_top=10, atlas_name=ATLAS_NAME, use_lr=use_lr_try,
                fp_test="species_mean", fdr_alpha=0.05,
                fig_dir=FC_RESULTS_DIR, verbose=True)
            print(f"  [Fig4] use_lr={use_lr_try} succeeded")
            break
        except ValueError as e:
            last_err = e
            print(f"  [Fig4] use_lr={use_lr_try} failed: {e}")
    if res is None:
        print(f"  [Fig4] both use_lr attempts failed — last error: {last_err}")
        print("  [Fig4] Hint: check that matrices exist under "
              "/scratch2/EDNiX3/{species}/{bids_dir}/sub-X/ses-Y/func/acpc-func/"
              "Stats/Correl_matrix/EDNIxCSC[LR]/correlation/*_run_*_correlation_matrix.csv")
        return

    base = opj(FC_RESULTS_DIR, "lvl2", "all")
    summary_csv = opj(base, "00_per_run_summary.csv")

    # ── log overall classification counts (4B values) ────────────────────────
    df_scored = None
    if os.path.exists(summary_csv):
        df_scored = pd.read_csv(summary_csv)
        if "fp_category" in df_scored.columns:
            n_tot = len(df_scored)
            for cat in ["Specific", "Unspecific", "No", "Spurious"]:
                nc = int((df_scored["fp_category"] == cat).sum())
                _log_value("Fig4B", f"overall_{cat}",
                           n=nc, pct=round(100 * nc / max(n_tot, 1), 1))
            _log_value("Fig4B", "total_runs", n=n_tot,
                       n_species=int(df_scored["species"].nunique())
                       if "species" in df_scored.columns else None)
            # per-species rates
            if "species" in df_scored.columns:
                for sp in df_scored["species"].unique():
                    sub = df_scored[df_scored["species"] == sp]
                    n = len(sub)
                    pct_spec = 100 * (sub["fp_category"] == "Specific").sum() / max(n, 1)
                    _log_value("Fig4B", f"species_{sp}_Specific",
                               n=n, pct=round(pct_spec, 1))

    # ── 4C: mouse anaesthesia regimes (by bids_dir within Mouse) ─────────────
    p4c = _fig4c_mouse_anaesthesia(df_scored,
                                   opj(ARTICLE_DIR, "panels", "Fig4C_mouse_anaesthesia.png"))

    # ── 4D: network metrics (eigenvalue ratio / silhouette) by category ──────
    p4d = None
    if FIG4_INCLUDE_D:
        p4d = _fig4d_network_metrics(df_qc, df_scored,
                                     opj(ARTICLE_DIR, "panels", "Fig4D_network_metrics.png"))

    # ── 4E: surrogate / perturbation sensitivity ─────────────────────────────
    p4e = _fig4e_surrogate(df_scored,
                           opj(ARTICLE_DIR, "panels", "Fig4E_surrogate.png"))

    # ── mount 4A-4E ───────────────────────────────────────────────────────────
    # Panel B: prefer per-BIDS proportions; fall back to per-species.
    p4b_bids = opj(base, "03_category_proportions_by_bids.png")
    p4b_sp   = opj(base, "03_category_proportions_by_species.png")
    p4b = p4b_bids if os.path.exists(p4b_bids) else p4b_sp
    panels = [
        opj(base, "01_threshold_selection_distributions.png"),   # 4A
        p4b,                                                      # 4B (per-BIDS)
        p4c,                                                      # 4C
        p4d,                                                      # 4D
        p4e,                                                      # 4E
    ]
    labels = ["A", "B", "C", "D", "E"]
    if not FIG4_INCLUDE_D:                  # drop "D is not interesting"
        panels.pop(3); labels = ["A", "B", "C", "D"]
    n = len([p for p in panels if p and os.path.exists(p)])
    layout = (2, 3) if n >= 5 else ((2, 2) if n >= 3 else (1, n))
    _mount(panels, layout=layout, out_name="Fig4_network_quality.png",
           labels=labels)


def _fig4c_mouse_anaesthesia(df_scored, out_path):
    """4C: classification proportions for the mouse anaesthesia regimes."""
    if df_scored is None or "species" not in df_scored.columns:
        print("  [Fig4C] no scored data")
        return None
    mouse = df_scored[df_scored["species"] == "Mouse"].copy()
    if mouse.empty or "fp_category" not in mouse.columns:
        print("  [Fig4C] no mouse runs")
        return None

    # group by bids_dir (each anaesthesia regime is a separate BIDS in your data)
    groups = sorted(mouse["bids_dir"].dropna().unique()) if "bids_dir" in mouse else []
    if not groups:
        return None

    CATS = ["Specific", "Unspecific", "No", "Spurious"]
    CAT_COLORS = {"Specific": "#009E73", "Unspecific": "#E69F00",
                  "No": "#D55E00", "Spurious": "#CC79A7"}
    import numpy as np
    fig, ax = plt.subplots(figsize=(max(5, len(groups) * 1.5), 4.5))
    x = np.arange(len(groups)); bottoms = np.zeros(len(groups))
    for cat in CATS:
        vals = []
        for g in groups:
            gd = mouse[mouse["bids_dir"] == g]["fp_category"]
            n = max(len(gd), 1)
            pct = 100 * (gd == cat).sum() / n
            vals.append(pct)
            _log_value("Fig4C", f"mouse_{g}_{cat}", pct=round(pct, 1), n=len(gd))
        vals = np.array(vals)
        ax.bar(x, vals, 0.7, bottom=bottoms, color=CAT_COLORS[cat],
               label=cat, edgecolor="white", lw=0.5)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 5:
                ax.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                        fontsize=8, color="white", fontweight="bold")
        bottoms += vals
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace("BIDS_", "") for g in groups],
                       rotation=25, ha="right", fontsize=9)
    ax.set_ylabel("% runs", fontsize=11); ax.set_ylim(0, 105)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4,
              fontsize=8, frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [panel] {out_path}")
    return out_path


def _fig4d_network_metrics(df_qc, df_scored, out_path):
    """
    4D: network-structure metrics (eigenvalue ratio, silhouette) by category.
    Tests Specific vs Spurious to confirm the classification captures real
    segregation differences.
    """
    import numpy as np
    from scipy import stats as _st
    if df_scored is None or "fp_category" not in df_scored.columns:
        print("  [Fig4D] no scored data")
        return None

    # merge network QC metrics from df_qc onto scored runs
    metrics = [("net_eigenvalue_ratio", "Eigenvalue ratio (λ₁/λ₂)"),
               ("net_silhouette", "Silhouette score")]
    have = [(m, lbl) for m, lbl in metrics
            if m in df_qc.columns or m in df_scored.columns]
    if not have:
        print("  [Fig4D] no network metrics in qc.csv or scored runs")
        return None

    key = [c for c in ("species", "subject", "session") if c in df_qc.columns
           and c in df_scored.columns]
    # only bring in metrics NOT already present in df_scored — avoids the
    # _x/_y suffix collision that drops the bare column name
    need = [m for m, _ in have if m not in df_scored.columns and m in df_qc.columns]
    if need and key:
        merged = df_scored.merge(df_qc[key + need], on=key, how="left")
    else:
        merged = df_scored.copy()
    have = [(m, lbl) for m, lbl in have if m in merged.columns]   # final guard
    if not have:
        print("  [Fig4D] network metrics unresolved after merge — skipping")
        return None

    CAT_COLORS = {"Specific": "#009E73", "Unspecific": "#E69F00",
                  "No": "#D55E00", "Spurious": "#CC79A7"}
    CATS = ["Specific", "Unspecific", "No", "Spurious"]
    fig, axes = plt.subplots(1, len(have), figsize=(len(have) * 5, 4.5))
    if len(have) == 1:
        axes = [axes]
    rng = np.random.default_rng(0)
    for ax, (m, lbl) in zip(axes, have):
        for ci, cat in enumerate(CATS):
            vals = merged.loc[merged["fp_category"] == cat, m].dropna().values
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                continue
            jit = rng.uniform(-0.18, 0.18, len(vals))
            ax.scatter(np.full(len(vals), ci) + jit, vals, s=14, alpha=0.5,
                       color=CAT_COLORS[cat], edgecolors="none")
            med = np.median(vals)
            ax.plot([ci - 0.25, ci + 0.25], [med, med],
                    color=CAT_COLORS[cat], lw=2.5, zorder=5)
        ax.set_xticks(range(len(CATS))); ax.set_xticklabels(CATS, fontsize=8,
                                                            rotation=20, ha="right")
        ax.set_ylabel(lbl, fontsize=10)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        # Specific vs Spurious test
        spec = merged.loc[merged["fp_category"] == "Specific", m].dropna()
        spur = merged.loc[merged["fp_category"] == "Spurious", m].dropna()
        if len(spec) >= 2 and len(spur) >= 2:
            t, p = _st.ttest_ind(spec, spur, equal_var=False)
            _log_value("Fig4D", f"{m}_Specific_vs_Spurious",
                       specific_median=round(float(spec.median()), 3),
                       spurious_median=round(float(spur.median()), 3),
                       t=round(float(t), 3), p=round(float(p), 4),
                       n_spec=len(spec), n_spur=len(spur))
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [panel] {out_path}")
    return out_path


def _fig4e_surrogate(df_scored, out_path):
    """
    4E: surrogate sensitivity. threshold_explorer already produces
    08_human_surrogate_sensitivity.png — reuse it. Fall back to the
    perturbation module only if that file is absent.
    """
    pre = opj(FC_RESULTS_DIR, "lvl2", "all", "08_human_surrogate_sensitivity.png")
    if os.path.exists(pre):
        print(f"  [Fig4E] using threshold_explorer surrogate: {pre}")
        return pre
    if df_scored is None:
        print("  [Fig4E] no scored data and no surrogate PNG")
        return None
    try:
        from ednix_fc_perturbation import run_perturbation_analysis
    except ImportError:
        print("  [Fig4E] no surrogate PNG and ednix_fc_perturbation not importable")
        return None
    cats = df_scored["fp_category"] if "fp_category" in df_scored.columns else None
    if cats is None:
        return None
    t_in = THRESH_INTRA if THRESH_INTRA is not None else 0.24
    t_de = THRESH_DELTA if THRESH_DELTA is not None else 0.10
    try:
        run_perturbation_analysis(
            df_scored, cats, thresh_intra=t_in, thresh_delta=t_de,
            species="Human", bids_root_template="/scratch2/EDNiX3/{species}/{bids_dir}",
            output_dir=os.path.dirname(out_path), max_runs=10)
        produced = opj(os.path.dirname(out_path), "perturbation_Human.png")
        if os.path.exists(produced):
            return produced
    except Exception as e:
        print(f"  [Fig4E] surrogate failed: {e}")
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 5 — SEED-BASED + FINGERPRINT
# ═══════════════════════════════════════════════════════════════════════════════

def fig5():
    """
    5  retrosplenial seed maps (mounts existing niftotoWBsurface PNGs)
    5c cross-species fingerprint similarity (pgls_connect machinery)
    """
    print("\n" + "=" * 64 + "\n  FIG 5 — Seed-based connectivity + fingerprint\n" + "=" * 64)

    # 5c — fingerprint similarity (always runnable from corr matrices)
    from pgls_connect import (_build_species_config as pgls_cfg,
                              _collect_species_mean_matrices,
                              compute_fingerprint_similarity,
                              plot_fingerprint_heatmap, N_PERMUTATIONS)
    species_config, _, _ = pgls_cfg()
    species_data = _collect_species_mean_matrices(species_config)
    fp_dir = opj(ARTICLE_DIR, "panels", "fig5c_fingerprint")
    if len(species_data) >= 2:
        df_sim, df_p = compute_fingerprint_similarity(species_data,
                                                      n_perm=N_PERMUTATIONS)
        os.makedirs(fp_dir, exist_ok=True)
        plot_fingerprint_heatmap(df_sim, df_p, fp_dir)
    else:
        print("  [Fig5c] <2 species with matrices — fingerprint skipped")

    # 5 — seed-based surface renders (external PNGs)
    if SEED_PNG:
        from Statistics.Evo.PGLS import _phylo_sort
        order = _phylo_sort(list(SEED_PNG.keys()))
        paths = [SEED_PNG[s] for s in order if s in SEED_PNG]
        n = len(paths); ncols = min(4, n); nrows = (n + ncols - 1) // ncols
        _mount(paths, layout=(nrows, ncols), out_name="Fig5_seed_based.png",
               figsize=(ncols * 4, nrows * 4),
               title="Figure 5 | Seed-based connectivity from retrosplenial cortex")
    else:
        print("  [Fig5] no SEED_PNG supplied — add niftotoWBsurface renders to CONFIG")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 6 — FDG-PET
# ═══════════════════════════════════════════════════════════════════════════════

def fig6():
    print("\n" + "=" * 64 + "\n  FIG 6 — FDG-PET glucose metabolism\n" + "=" * 64)
    if not PET_PNG:
        print("  [Fig6] no PET_PNG supplied — add to CONFIG when PET maps ready")
        return
    paths = list(PET_PNG.values())
    n = len(paths); ncols = min(2, n); nrows = (n + ncols - 1) // ncols
    _mount(paths, layout=(nrows, ncols), out_name="Fig6_pet.png",
           labels=["A", "B", "C", "D"][:n],
           title="Figure 6 | Cross-species glucose metabolism (FDG-PET)")


# ═══════════════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY — correlations + per-region trajectory mounts + QC densities
# ═══════════════════════════════════════════════════════════════════════════════

def supplementary():
    """
    Supplementary figures, numbered in the paper's logical order:

      Suppl 1-3   QC distributions (anatomical / functional / network)
      Suppl 4-6   morphometry by species & dataset (surface / volume / thickness)
      Suppl 7-9   regional outlier rates (surface / volume / thickness)
      Suppl 10    allocortex / periallocortex scaling
      Suppl 11-13 FC classification matrices by category (all / primates / +rodents)
      Suppl 14    mouse anaesthesia regimes

    All reported correlations land in stats/paper_values.csv via _log_value.
    """
    print("\n" + "=" * 64 + "\n  SUPPLEMENTARY (paper order)\n" + "=" * 64)
    import numpy as np
    from scipy import stats as _st
    from ednix_bids_tools import collect_multi_species
    sup = opj(ARTICLE_DIR, "supplementary")
    os.makedirs(sup, exist_ok=True)

    species_config, atlas_label_paths, atlas_frags = _build_species_config()
    data = collect_multi_species(
        species_config, regions_of_interest=[],
        extract=("surface", "volume", "thickness", "qc"), atlas_name=ATLAS_NAME,
        atlas_label_paths=atlas_label_paths, atlas_library_root=ATLAS_LIB,
        species_atlas_fragments=atlas_frags)
    qc = data.get("qc")

    # ── Suppl 1-3: QC distributions (anat / func / network) ──────────────────
    if qc is not None and not qc.empty:
        qcols = set(qc.columns)
        def _present(metrics):
            return [(m, lbl, unit) for m, lbl, unit in metrics if m in qcols]
        try:
            from ednix_fc_perturbation import plot_qc_densities_per_metric
            # Suppl 1 — anatomical QC
            anat_m = _present([
                ("anat_avg_snr_gray", "Gray-matter SNR", "SNR"),
                ("anat_cnr", "Contrast-to-noise", "CNR"),
                ("anat_cortical_contrast", "Cortical contrast", "ratio"),
                ("anat_template_correlation", "Template registration", "r"),
                ("anat_efc", "Entropy focus", "EFC"),
                ("anat_fwhm_avg", "Smoothness", "FWHM")])
            if anat_m:
                plot_qc_densities_per_metric(qc, metrics=anat_m,
                                             output_dir=opj(sup, "Suppl1_anat_QC"))
            # Suppl 2 — functional QC (real column names from qc.csv)
            func_m = _present([
                ("func_TSNR", "Temporal SNR", "tSNR"),
                ("func_avg_snr_gray", "Gray-matter SNR", "SNR"),
                ("func_cnr", "Contrast-to-noise", "CNR"),
                ("func_mean_fd", "Mean framewise displacement", "mm"),
                ("func_mean_dvars", "Mean DVARS", "a.u."),
                ("func_gcor", "Global correlation", "GCOR"),
                ("func_ghost_ratio", "Ghost-to-signal", "ratio")])
            if func_m:
                plot_qc_densities_per_metric(qc, metrics=func_m,
                                             output_dir=opj(sup, "Suppl2_func_QC"))
            # Suppl 3 — network QC
            net_m = _present([
                ("net_eigenvalue_ratio", "Eigenvalue ratio", "λ₁/λ₂"),
                ("net_top_eigenvalue", "Top eigenvalue", "λ₁"),
                ("net_silhouette", "Silhouette", "score"),
                ("net_davies_bouldin", "Davies–Bouldin", "index")])
            if net_m:
                plot_qc_densities_per_metric(qc, metrics=net_m,
                                             output_dir=opj(sup, "Suppl3_network_QC"))
        except ImportError:
            print("  [Suppl1-3] ednix_fc_perturbation not importable")

    # ── Suppl 4-6: combo plots for all regions at every atlas level ──────────
    # Suppl 4 → lvl 2 (all regions at that level)
    # Suppl 5 → lvl 3
    # Suppl 6 → lvl 4
    # Lvl 1 (Iso/Allo/Periallo) is in main Fig 2A, so not duplicated.
    level_to_suppl = {2: 4, 3: 5, 4: 6}
    for lvl, suppl_n in level_to_suppl.items():
        try:
            _combo_plot(
                data, regions=None,   # None = use all regions present at that level
                atlas_level=lvl, hemisphere="bilateral",
                out_path=opj(sup, f"Suppl{suppl_n}_combo_lvl{lvl}_bilateral_log.png"))
        except Exception as e:
            print(f"  [Suppl{suppl_n}] combo lvl{lvl} failed: {e}")

    # ── Suppl 7-9: regional outlier rates ────────────────────────────────────
    for n, (key, col) in zip(
            [7, 8, 9],
            [("surface", "surface_area_mm2"), ("volume", "volume_mm3"),
             ("thickness", "thickness_mm")]):
        df = data.get(key)
        if df is None or df.empty:
            continue
        _suppl_outlier_rates(df, col, opj(sup, f"Suppl{n}_{key}_outliers.csv"))

    # ── Suppl 11-13: FC classification matrices by subset ────────────────────
    # threshold_explorer writes per-subset folders; mount their category-matrix
    # figure under explicit Suppl11/12/13 names.
    subset_map = [(11, "all"), (12, "primates"), (13, "primates_rodents")]
    for n, subset in subset_map:
        cand = [
            opj(FC_RESULTS_DIR, "lvl2", subset, "04_fc_matrices_by_category_species.png"),
            opj(FC_RESULTS_DIR, "lvl2", subset, "02_classification_model_per_species.png"),
            opj(FC_RESULTS_DIR, "lvl2", subset, "06_fingerprint_heatmap_by_bids.png"),
        ]
        found = [p for p in cand if os.path.exists(p)]
        if found:
            _mount(found, layout=(len(found), 1),
                   out_name=f"../supplementary/Suppl{n}_FC_{subset}.png",
                   figsize=(7, len(found) * 5))
        else:
            print(f"  [Suppl{n}] no FC panels for subset '{subset}' "
                  f"(run Fig 4 first to generate threshold_explorer outputs)")

    # ── Suppl 14: mouse anaesthesia (copy the Fig4C panel) ───────────────────
    src14 = opj(ARTICLE_DIR, "panels", "Fig4C_mouse_anaesthesia.png")
    if os.path.exists(src14):
        import shutil
        dst14 = opj(sup, "Suppl14_mouse_anaesthesia.png")
        shutil.copy(src14, dst14)
        print(f"  [supp] {dst14}")
    else:
        print("  [Suppl14] Fig4C not found — run Fig 4 first")

    # ── correlations into paper_values (inter-modality) ──────────────────────
    surf, vol, thk = data.get("surface"), data.get("volume"), data.get("thickness")
    if all(d is not None and not d.empty for d in (surf, vol, thk)):
        kcols = ["species", "subject", "region"]
        m = (surf[kcols + ["surface_area_mm2"]]
             .merge(vol[kcols + ["volume_mm3"]], on=kcols)
             .merge(thk[kcols + ["thickness_mm"]], on=kcols))
        for sp in m["species"].unique():
            sm = m[m["species"] == sp]
            for a, b in [("surface_area_mm2", "volume_mm3"),
                         ("surface_area_mm2", "thickness_mm"),
                         ("volume_mm3", "thickness_mm")]:
                d = sm[[a, b]].dropna()
                if len(d) < 5:
                    continue
                r, p = _st.pearsonr(np.log10(d[a].clip(lower=1e-9)),
                                    np.log10(d[b].clip(lower=1e-9)))
                z = np.arctanh(r); se = 1 / np.sqrt(len(d) - 3)
                _log_value("Suppl_corr", f"{sp}_{a}_vs_{b}",
                           n=len(d), R=round(r, 3), p=round(p, 4),
                           ci_low=round(float(np.tanh(z - 1.96 * se)), 3),
                           ci_high=round(float(np.tanh(z + 1.96 * se)), 3))

    print(f"  Supplementary outputs → {sup}/")


def _suppl_morphometry_by_species(df, col, label, out_path):
    """Strip plot of one morphometry metric across species & datasets."""
    import numpy as np
    sp_list = [s for s in _PHYLO_ORDER if s in df["species"].unique()]
    fig, ax = plt.subplots(figsize=(max(6, len(sp_list) * 1.3), 4.5))
    rng = np.random.default_rng(0)
    for si, sp in enumerate(sp_list):
        vals = df.loc[df["species"] == sp, col].dropna().values
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        jit = rng.uniform(-0.2, 0.2, len(vals))
        ax.scatter(np.full(len(vals), si) + jit, vals, s=10, alpha=0.4,
                   color=_SPECIES_COLORS.get(sp, "#888"), edgecolors="none")
        ax.plot([si - 0.25, si + 0.25], [np.median(vals)] * 2,
                color=_SPECIES_COLORS.get(sp, "#888"), lw=2.4, zorder=5)
    ax.set_yscale("log")
    ax.set_xticks(range(len(sp_list)))
    ax.set_xticklabels(sp_list, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel(label, fontsize=10)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [supp] {out_path}")


def _suppl_outlier_rates(df, col, out_path):
    """Per-region outlier rate (>2.5 MAD from species median) → CSV."""
    import numpy as np
    rows = []
    for (sp, region), g in df.groupby(["species", "region"]):
        vals = g[col].dropna().values
        vals = vals[np.isfinite(vals)]
        if len(vals) < 5:
            continue
        med = np.median(vals)
        mad = np.median(np.abs(vals - med)) + 1e-9
        n_out = int(np.sum(np.abs(vals - med) / (1.4826 * mad) > 2.5))
        rows.append(dict(species=sp, region=region, n=len(vals),
                         n_outliers=n_out, outlier_rate=round(n_out / len(vals), 3)))
    if rows:
        pd.DataFrame(rows).sort_values("outlier_rate", ascending=False)\
            .to_csv(out_path, index=False)
        print(f"  [supp] {out_path}")




FIGURES = {2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6}

def main():
    ap = argparse.ArgumentParser(description="EDNiX paper figure orchestrator")
    ap.add_argument("--only", type=int, nargs="+", default=None,
                    help="figure numbers to run (e.g. --only 3 4). Default: all")
    ap.add_argument("--supp", action="store_true",
                    help="also run supplementary figures + correlations")
    ap.add_argument("--supp-only", action="store_true",
                    help="run ONLY supplementary (skip main figures)")
    args = ap.parse_args()
    _dirs()

    if args.supp_only:
        supplementary()
        _flush_paper_values()
        print(f"\n  Done. Outputs → {ARTICLE_DIR}/")
        return

    to_run = args.only or sorted(FIGURES.keys())
    for n in to_run:
        if n in FIGURES:
            try:
                FIGURES[n]()
            except Exception as e:
                import traceback
                print(f"\n  [Fig {n} FAILED] {e}")
                traceback.print_exc()
        else:
            print(f"  [skip] no Fig {n}")

    if args.supp:
        try:
            supplementary()
        except Exception as e:
            import traceback
            print(f"\n  [Supplementary FAILED] {e}")
            traceback.print_exc()

    _flush_paper_values()
    print(f"\n  Done. Outputs → {ARTICLE_DIR}/{{panels,mounted,stats,supplementary}}/")


if __name__ == "__main__":
    main()