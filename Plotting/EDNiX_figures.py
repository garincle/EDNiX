"""
EDNiX Figures  —  v7b
======================
Based on v6 (doc 7). Targeted fixes and new features.

CHANGES vs v6
-------------
FIX 1  — _make_bilateral_volume: robust L+R aggregation that handles mixed
          hemisphere encodings (some subjects bilateral, others L/R separate).
          Subjects that already have bilateral rows are kept as-is; only
          subjects with only L/R rows get their hemispheres summed.

FIX 2  — plot_fc_report_by_category: All category removed from matrices;
          categories = Specific / Unspecific / No only (3 matrix rows).
          Raw mean matrices shown (not normalised).
          Bar plots use inter_mean and func_gcor instead of intra_mean / mean_fd.
          "No specificity" label shortened to "No".
          Y-axis label moved inside the panel (no overlap with species names).
          Lemur n=0 panels now show correctly even when bar plot has data.

FIX 3  — Anesthesia comparison plot: plot_anesthesia_comparison() generates
          a specificity-style stacked bar showing Awake vs Anaesthetised
          proportions per species and per BIDS dir.

FIX 4  — plot_corr_matrix_per_bids: uses same ROI set as per-species matrices
          (fixes ROI mismatch).

FIX 5  — plot_corr_matrix_sig_norm: normalises by sp_specific_correlation
          from df_qc when available.
"""

import os
import math
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.cm as _cm
from matplotlib.colors import TwoSlopeNorm, Normalize

from Plotting.ednix_bids_tools import (
    PAPER_RC, PALETTE,
    AWAKE_MARKER, ANESTH_MARKER,
    _bids_label, _bids_offsets, _phylo_sort,
    _anesthesia_marker, _anesthesia_legend_handles,
    _violin_strip_quartiles, _scatter_with_anesth,
    _df_region_mask, _filter_atlas_level, _hemis_for,
    collect_corr_matrices, collect_corr_matrices_per_bids,
    extract_corr_matrix_paths, load_corr_matrix,
    plot_morphometry_intra_bids, plot_qc_dashboard,
    _AWAKE_KEYWORDS,
)

opj  = os.path.join
_rng = np.random.default_rng(42)

_SPEC_CATS   = ["Specific", "Unspecific", "No"]
_SPEC_COLORS = {
    "Specific":   "#009E73",
    "Unspecific": "#E69F00",
    "No":         "#D55E00",
    "Spurious":   "#CC79A7",
}
_SPEC_LABEL  = {"Specific": "Specific", "Unspecific": "Unspecific",
                "No": "No", "Spurious": "Spurious"}


# ═══════════════════════════════════════════════════════════════════════════════
# SHARED LEGEND HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def _bids_legend_handles(bids_dirs, bids_col):
    return [mpatches.Patch(facecolor=bids_col.get(_bids_label(bd), "#888"),
                           label=_bids_label(bd))
            for bd in bids_dirs]


def _add_bids_legend(fig, bids_dirs, bids_col, anesth_map=None, bottom_offset=0.10):
    handles = _bids_legend_handles(bids_dirs, bids_col)
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -bottom_offset),
               ncol=min(len(handles), 6),
               title="BIDS dir", fontsize=8, frameon=False)
    if anesth_map:
        fig.legend(handles=_anesthesia_legend_handles(), loc="lower right",
                   bbox_to_anchor=(1.0, 0.0), ncol=1,
                   title="Condition", fontsize=8, frameon=False)


# ═══════════════════════════════════════════════════════════════════════════════
# FIX 1 — ROBUST BILATERAL VOLUME AGGREGATION
# ═══════════════════════════════════════════════════════════════════════════════

def _make_bilateral_volume(df_volume):
    """
    Robust L+R → bilateral aggregation.

    Handles three cases per (species, bids_dir, subject, session, region):
      A) Only bilateral row(s) already → keep as-is (sum in case of duplicates)
      B) Only L + R rows → sum into one bilateral row
      C) Mix of bilateral AND L/R rows → keep bilateral, discard L/R

    This avoids triple-counting when df_volume has both bilateral and L/R entries
    for the same subject (which causes the bimodal violin artefact).
    """
    if df_volume is None or df_volume.empty:
        return df_volume
    if "hemisphere" not in df_volume.columns:
        return df_volume

    grp_cols = [c for c in ("species", "bids_dir", "subject", "session",
                             "region", "atlas_level")
                if c in df_volume.columns]

    rows_out = []
    for key, grp in df_volume.groupby(grp_cols):
        bil = grp[grp["hemisphere"] == "bilateral"]
        lr  = grp[grp["hemisphere"].isin(["left", "right"])]

        if not bil.empty:
            # Case A or C: use bilateral rows only
            row = bil.iloc[0].copy()
            row["volume_mm3"]  = float(bil["volume_mm3"].sum())
            row["hemisphere"]  = "bilateral"
            rows_out.append(row)
        elif not lr.empty:
            # Case B: sum L + R
            row = lr.iloc[0].copy()
            row["volume_mm3"] = float(lr["volume_mm3"].sum())
            row["hemisphere"] = "bilateral"
            rows_out.append(row)
        else:
            # Unexpected hemisphere value — keep first row
            row = grp.iloc[0].copy()
            row["hemisphere"] = "bilateral"
            rows_out.append(row)

    result   = pd.DataFrame(rows_out)
    n_before = len(df_volume)
    n_after  = len(result)
    if n_before != n_after:
        print(f"  [bilateral_vol] {n_before} → {n_after} rows (L+R consolidated)")
    return result


def _mean_per_subject(df, metric_col):
    """Average metric across sessions → one row per subject-region."""
    if df is None or df.empty or metric_col not in df.columns:
        return df
    if "session" not in df.columns:
        return df
    excl = {metric_col, "session", "voxel_count", "label_id",
            "region_name", "region_full", "n_sessions"}
    id_c = [c for c in df.columns if c not in excl]
    n    = len(df)
    out  = df.groupby(id_c, as_index=False)[metric_col].mean()
    if len(out) != n:
        print(f"  [session_avg] {metric_col}: {n} → {len(out)} rows")
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# §B  MORPHOMETRY FIGURES  (v6 unchanged, violin colour grey)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_cross_species_dots(
    df, metric_col, metric_label, regions, output_path,
    hemisphere="bilateral", log_scale=False,
    normalise_by_brain=False, df_vol_for_norm=None,
    atlas_level=1, bids_col=None,
    anesth_map=None, figsize=None,
):
    with plt.rc_context(PAPER_RC):
        plot_df = _filter_atlas_level(df.copy(), atlas_level)
        if hemisphere == "bilateral":
            bil = plot_df[plot_df["hemisphere"] == "bilateral"]
            if bil.empty and "hemisphere" in plot_df.columns:
                grp = [c for c in ("species","bids_dir","subject","session",
                                   "region","atlas_level") if c in plot_df.columns]
                plot_df = plot_df.groupby(grp)[metric_col].sum().reset_index()
                plot_df["hemisphere"] = "bilateral"
            else:
                plot_df = bil
        elif hemisphere:
            filt = plot_df[plot_df["hemisphere"] == hemisphere]
            plot_df = filt if not filt.empty else \
                      plot_df[plot_df["hemisphere"] == "bilateral"]

        plot_df = plot_df[_df_region_mask(plot_df["region"], regions)]
        if plot_df.empty:
            warnings.warn(f"plot_cross_species_dots: no data for {metric_col} hemi={hemisphere}")
            return None

        if normalise_by_brain and df_vol_for_norm is not None:
            keys = [c for c in ("subject","session","bids_dir")
                    if c in plot_df.columns and c in df_vol_for_norm.columns]
            bt = (df_vol_for_norm.groupby(keys)["volume_mm3"].sum().reset_index()
                  .rename(columns={"volume_mm3": "_brain_total"}))
            plot_df = plot_df.merge(bt, on=keys, how="left")
            plot_df[metric_col] = plot_df[metric_col] / plot_df["_brain_total"]
            metric_label += " / brain vol."

        species_order = _phylo_sort(plot_df["species"].unique())
        bids_dirs     = sorted(plot_df["bids_dir"].unique()) \
                        if "bids_dir" in plot_df.columns else []
        _bids_col = dict(bids_col) if bids_col else {}
        for i, bd in enumerate([b for b in bids_dirs if b not in _bids_col]):
            _bids_col[bd] = PALETTE[(len(_bids_col)+i) % len(PALETTE)]

        n_reg = len(regions)
        w, h  = figsize or (max(7, n_reg*3.2), 5.5)
        fig, axes = plt.subplots(1, n_reg, figsize=(w, h), sharey=False)
        if n_reg == 1: axes = [axes]

        for ax, region in zip(axes, regions):
            rdf = plot_df[_df_region_mask(plot_df["region"], [region])]
            if rdf.empty: ax.set_visible(False); continue
            for xi, sp in enumerate(species_order):
                sp_df = rdf[rdf["species"] == sp]
                vals  = sp_df[metric_col].dropna().values
                if len(vals) == 0: continue
                _violin_strip_quartiles(ax, vals, xi, "#aaaaaa")
                _bd_list   = sorted(sp_df["bids_dir"].unique()) \
                             if "bids_dir" in sp_df.columns else [""]
                _bd_offset = _bids_offsets(_bd_list, spread=0.26)
                for bd in _bd_list:
                    bd_df   = sp_df[sp_df["bids_dir"] == bd] \
                              if "bids_dir" in sp_df.columns else sp_df
                    bd_vals = bd_df[metric_col].dropna().values
                    if len(bd_vals) == 0: continue
                    jitter = _rng.uniform(-0.08, 0.08, len(bd_vals))
                    x_draw = xi + _bd_offset.get(bd, 0) + jitter
                    subs   = bd_df.loc[bd_df[metric_col].notna(), "subject"].values \
                             if "subject" in bd_df.columns else None
                    sess_  = bd_df.loc[bd_df[metric_col].notna(), "session"].values \
                             if "session" in bd_df.columns else None
                    _scatter_with_anesth(ax, x_draw, bd_vals, subs, sess_,
                                         anesth_map, _bids_col.get(bd, PALETTE[0]), _rng)
            if log_scale:
                ax.set_yscale("log")
            else:
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
            ax.set_xticks(range(len(species_order)))
            ax.set_xticklabels(species_order, rotation=35, ha="right")
            ax.set_title(region, fontweight="bold")
            ax.set_ylabel(metric_label if ax is axes[0] else "")
            ax.set_xlim(-0.65, len(species_order)-0.35)

        for _pi, _ax in enumerate(axes):
            _ax.text(-0.08, 1.02, chr(65+_pi), transform=_ax.transAxes,
                     fontsize=13, fontweight="bold", va="bottom", ha="left")

        plt.tight_layout()
        _add_bids_legend(fig, bids_dirs, _bids_col, anesth_map, bottom_offset=0.10)
        fig.subplots_adjust(bottom=0.18)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def plot_combo(df_surface, df_volume, df_thickness, output_path, regions,
               hemisphere="bilateral", log_scale=False, normalise_by_brain=False,
               atlas_level=1, bids_col=None, anesth_map=None):
    modalities = [
        (df_surface,   "surface_area_mm2", "Surface area (mm²)"),
        (df_volume,    "volume_mm3",        "Volume (mm³)"),
        (df_thickness, "thickness_mm",      "Cortical thickness (mm)"),
    ]
    modalities = [(df,col,lbl) for df,col,lbl in modalities
                  if df is not None and not df.empty and col in df.columns]
    if not modalities: warnings.warn("plot_combo: no data"); return None
    n_mod = len(modalities); n_reg = len(regions)
    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(n_mod, n_reg,
                                 figsize=(max(8, n_reg*3.0), n_mod*3.5),
                                 squeeze=False)
        all_bids_dirs = []
        for row, (df, metric_col, metric_label) in enumerate(modalities):
            plot_df = _filter_atlas_level(df.copy(), atlas_level)
            if hemisphere == "bilateral":
                bil = plot_df[plot_df["hemisphere"] == "bilateral"]
                if bil.empty and "hemisphere" in plot_df.columns:
                    grp = [c for c in ("species","bids_dir","subject","session",
                                       "region","atlas_level") if c in plot_df.columns]
                    plot_df = plot_df.groupby(grp)[metric_col].sum().reset_index()
                    plot_df["hemisphere"] = "bilateral"
                else:
                    plot_df = bil
            elif hemisphere:
                filt = plot_df[plot_df["hemisphere"] == hemisphere]
                plot_df = filt if not filt.empty else plot_df[plot_df["hemisphere"]=="bilateral"]
            plot_df = plot_df[_df_region_mask(plot_df["region"], regions)]
            if normalise_by_brain and df_volume is not None:
                keys = [c for c in ("subject","session","bids_dir")
                        if c in plot_df.columns and c in df_volume.columns]
                bt = (df_volume.groupby(keys)["volume_mm3"].sum().reset_index()
                      .rename(columns={"volume_mm3": "_brain_total"}))
                plot_df = plot_df.merge(bt, on=keys, how="left")
                if "_brain_total" in plot_df.columns:
                    plot_df[metric_col] = plot_df[metric_col] / plot_df["_brain_total"]
                    metric_label += " / brain vol."
            species_order = _phylo_sort(plot_df["species"].unique())
            bids_dirs     = sorted(plot_df["bids_dir"].unique()) \
                            if "bids_dir" in plot_df.columns else []
            all_bids_dirs = sorted(set(all_bids_dirs + bids_dirs))
            _bids_col = dict(bids_col) if bids_col else {}
            for i,bd in enumerate([b for b in bids_dirs if b not in _bids_col]):
                _bids_col[bd] = PALETTE[(len(_bids_col)+i) % len(PALETTE)]
            for col, region in enumerate(regions):
                ax  = axes[row][col]
                rdf = plot_df[_df_region_mask(plot_df["region"], [region])]
                for xi, sp in enumerate(species_order):
                    sp_df = rdf[rdf["species"] == sp]
                    vals  = sp_df[metric_col].dropna().values
                    if len(vals) == 0: continue
                    _violin_strip_quartiles(ax, vals, xi, "#aaaaaa")
                    if "bids_dir" in sp_df.columns:
                        _bd_list   = sorted(sp_df["bids_dir"].unique())
                        _bd_offset = _bids_offsets(_bd_list, spread=0.26)
                        for bd in _bd_list:
                            bd_df   = sp_df[sp_df["bids_dir"] == bd]
                            bd_vals = bd_df[metric_col].dropna().values
                            if len(bd_vals) == 0: continue
                            jitter = _rng.uniform(-0.08, 0.08, len(bd_vals))
                            x_draw = xi + _bd_offset[bd] + jitter
                            subs   = bd_df.loc[bd_df[metric_col].notna(),"subject"].values \
                                     if "subject" in bd_df.columns else None
                            sess_  = bd_df.loc[bd_df[metric_col].notna(),"session"].values \
                                     if "session" in bd_df.columns else None
                            _scatter_with_anesth(ax, x_draw, bd_vals, subs, sess_,
                                                  anesth_map, _bids_col[bd], _rng)
                if log_scale: ax.set_yscale("log")
                else: ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
                ax.set_xticks(range(len(species_order)))
                if row == n_mod-1:
                    ax.set_xticklabels(species_order, rotation=35, ha="right", fontsize=8)
                else:
                    ax.set_xticklabels([], fontsize=0)
                ax.set_xlim(-0.65, len(species_order)-0.35)
                if col == 0: ax.set_ylabel(metric_label, fontsize=9)
                if row == 0: ax.set_title(region, fontweight="bold", fontsize=9)
        plt.tight_layout()
        handles = [mpatches.Patch(facecolor=_bids_col.get(bd,"#888"), label=_bids_label(bd))
                   for bd in all_bids_dirs if bd in _bids_col]
        if handles:
            _add_bids_legend(fig, all_bids_dirs, _bids_col, anesth_map, bottom_offset=0.02)
            fig.subplots_adjust(bottom=0.10)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
# §C  QC FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

SHARED_NMI_COLS = ["func_nmi", "anat_nmi", "avg_snr_gray"]
FUNC_QC_COLS    = ["func_TSNR_0", "func_mean_fd", "func_gcor"]
ANAT_QC_COLS    = ["anat_template_correlation", "anat_cortical_contrast"]
QC_CLIP_COLS    = {"anat_fwhm": 99, "func_TSNR_0": 99}
FUNC_NET_COLS   = [
    "net_mean_correlation", "net_std_correlation",
    "net_eigenvalue_ratio", "net_davies_bouldin",
    "net_inter_mean", "net_intra_mean",
    "net_homotopic_mean", "net_cross_mean",
    "sp_specific_correlation", "sp_nonspecific_correlation",
    "sp_specificity_index",
]
FUNC_SPEC_COL = "sp_category"


def _draw_qc_panel(ax, qc_df, metric, species_order, bids_col, anesth_map,
                   clip_cols=None):
    clip_cols = clip_cols or QC_CLIP_COLS
    mdf = qc_df[["species","bids_dir",metric]].dropna() \
          if "bids_dir" in qc_df.columns else qc_df[["species",metric]].dropna()
    for xi, sp in enumerate(species_order):
        sp_df = mdf[mdf["species"] == sp]
        vals  = sp_df[metric].values
        if len(vals) == 0: continue
        _violin_strip_quartiles(ax, vals, xi, "#aaaaaa")
        if "bids_dir" in sp_df.columns:
            _bd_list   = sorted(sp_df["bids_dir"].unique())
            _bd_offset = _bids_offsets(_bd_list, spread=0.26)
            for bd in _bd_list:
                bd_df   = sp_df[sp_df["bids_dir"] == bd]
                bd_vals = bd_df[metric].values
                if len(bd_vals) == 0: continue
                jitter = _rng.uniform(-0.10, 0.10, len(bd_vals))
                x_draw = xi + _bd_offset[bd] + jitter
                color  = bids_col.get(bd, bids_col.get(_bids_label(bd), "#888"))
                subs   = bd_df.loc[bd_df[metric].notna(),"subject"].values \
                         if "subject" in bd_df.columns else None
                sess_  = bd_df.loc[bd_df[metric].notna(),"session"].values \
                         if "session" in bd_df.columns else None
                _scatter_with_anesth(ax, x_draw, bd_vals, subs, sess_,
                                     anesth_map, color, _rng)
    ax.set_xticks(range(len(species_order)))
    ax.set_xticklabels(species_order, rotation=35, ha="right", fontsize=8)
    clean = metric.replace("func_","").replace("anat_","").replace("_"," ")
    ax.set_title(clean, fontweight="bold", fontsize=9)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    if metric in clip_cols:
        try:
            all_vals = pd.to_numeric(mdf[metric], errors="coerce").dropna().values
            if all_vals.size > 0:
                pct  = clip_cols[metric]
                ymax = float(np.nanpercentile(all_vals, pct))
                ymin = float(np.nanpercentile(all_vals, max(0, 100-pct)))
                ax.set_ylim(bottom=min(0,ymin)*0.95, top=ymax*1.05)
        except Exception as _e:
            warnings.warn(f"clip failed for {metric}: {_e}")


def plot_cross_species_qc(qc_df, output_path, anat_metrics=None, func_metrics=None,
                           shared_nmi=None, clip_cols=None, bids_col=None,
                           anesth_map=None, figsize=None):
    with plt.rc_context(PAPER_RC):
        shared_nmi   = [m for m in (shared_nmi   or SHARED_NMI_COLS)  if m in qc_df.columns]
        anat_metrics = [m for m in (anat_metrics or ANAT_QC_COLS)     if m in qc_df.columns]
        func_metrics = [m for m in (func_metrics or FUNC_QC_COLS)     if m in qc_df.columns]
        _bids_col    = bids_col or {}
        if not anat_metrics and not func_metrics and not shared_nmi:
            warnings.warn("plot_cross_species_qc: no QC metrics found"); return None
        species_order = _phylo_sort(qc_df["species"].unique())
        bids_dirs     = sorted(qc_df["bids_dir"].unique()) if "bids_dir" in qc_df.columns else []
        n_shared = 1 if shared_nmi else 0
        n_rows   = n_shared + max(len(func_metrics), len(anat_metrics), 1)
        w, h     = figsize or (11, n_rows*3.0)
        fig, axes = plt.subplots(n_rows, 2, figsize=(w, h), squeeze=False)
        if shared_nmi:
            nmi_func = next((m for m in shared_nmi if "func" in m), None)
            nmi_anat = next((m for m in shared_nmi if "anat" in m), None)
            for ci, metric in enumerate([nmi_func, nmi_anat]):
                if metric:
                    _draw_qc_panel(axes[0][ci], qc_df, metric, species_order,
                                   _bids_col, anesth_map, clip_cols)
                else: axes[0][ci].set_visible(False)
        for ri, metric in enumerate(func_metrics):
            _draw_qc_panel(axes[n_shared+ri][0], qc_df, metric, species_order,
                           _bids_col, anesth_map, clip_cols)
        for ri, metric in enumerate(anat_metrics):
            _draw_qc_panel(axes[n_shared+ri][1], qc_df, metric, species_order,
                           _bids_col, anesth_map, clip_cols)
        for ri in range(len(func_metrics), n_rows-n_shared): axes[n_shared+ri][0].set_visible(False)
        for ri in range(len(anat_metrics), n_rows-n_shared): axes[n_shared+ri][1].set_visible(False)
        plt.tight_layout()
        if bids_dirs:
            _add_bids_legend(fig, bids_dirs, _bids_col, anesth_map, bottom_offset=0.04)
            fig.subplots_adjust(bottom=0.09, top=0.92)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


def plot_qc_fc_report(qc_df, output_path, fc_cols=None, bids_col=None,
                       anesth_map=None, figsize=None):
    with plt.rc_context(PAPER_RC):
        fc_cols   = [c for c in (fc_cols or FUNC_NET_COLS) if c in qc_df.columns]
        _bids_col = bids_col or {}
        if not fc_cols: warnings.warn("plot_qc_fc_report: no columns"); return None
        species_order = _phylo_sort(qc_df["species"].unique())
        bids_dirs     = sorted(qc_df["bids_dir"].unique()) if "bids_dir" in qc_df.columns else []
        label_map = {
            "net_inter_mean":             "Inter hemi (mean r)",
            "net_intra_mean":             "Intra (mean r)",
            "net_intra_left_mean":        "Intra L (mean r)",
            "net_intra_right_mean":       "Intra R (mean r)",
            "sp_specific_correlation":    "Specific pair r",
            "sp_nonspecific_correlation": "Non-specific pair r",
        }
        NCOLS = min(4, len(fc_cols))
        n_rows = math.ceil(len(fc_cols)/NCOLS)
        w, h   = figsize or (NCOLS*3.8, n_rows*3.5)
        fig, axes = plt.subplots(n_rows, NCOLS, figsize=(w, h))
        axes = np.array(axes).flatten()
        for i, metric in enumerate(fc_cols):
            _draw_qc_panel(axes[i], qc_df, metric, species_order, _bids_col, anesth_map)
            axes[i].set_title(label_map.get(metric, metric.replace("_"," ")),
                              fontweight="bold", fontsize=9)
            axes[i].yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
        for j in range(i+1, len(axes)): axes[j].set_visible(False)
        plt.tight_layout()
        if bids_dirs:
            _add_bids_legend(fig, bids_dirs, _bids_col, anesth_map, bottom_offset=0.07)
            fig.subplots_adjust(bottom=0.14)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


def plot_specificity_bar(qc_df, output_path, cat_col="sp_category",
                          group_by="species", figsize=None):
    with plt.rc_context(PAPER_RC):
        if cat_col not in qc_df.columns:
            warnings.warn(f"plot_specificity_bar: {cat_col} not found"); return None
        if group_by == "bids_dir" and "bids_dir" in qc_df.columns:
            group_vals = sorted(qc_df["bids_dir"].unique())
            col_src    = qc_df["bids_dir"]
        else:
            group_vals = _phylo_sort(qc_df["species"].unique())
            col_src    = qc_df["species"]
        cats   = ["Specific","Unspecific","Spurious","No"]
        colors = {"Specific":"#009E73","Unspecific":"#E69F00","Spurious":"#CC79A7","No":"#D55E00"}
        counts = {}
        for g in group_vals:
            gdf   = qc_df[col_src == g][cat_col].dropna()
            total = len(gdf)
            counts[g] = {c: 100.0*(gdf==c).sum()/total for c in cats} \
                        if total > 0 else {c: 0.0 for c in cats}
        n_g  = len(group_vals)
        w, h = figsize or (max(5, n_g*1.2), 4.5)
        fig, ax = plt.subplots(figsize=(w, h))
        bottoms = np.zeros(n_g); x = np.arange(n_g)
        for cat in cats:
            vals = np.array([counts[g][cat] for g in group_vals])
            ax.bar(x, vals, 0.7, bottom=bottoms, color=colors[cat], label=cat,
                   edgecolor="white", linewidth=0.5)
            for xi, (v, b) in enumerate(zip(vals, bottoms)):
                if v > 5:
                    ax.text(xi, b+v/2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")
            bottoms += vals
        ax.set_xticks(x)
        ax.set_xticklabels(group_vals, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel("% subjects", fontsize=10)
        ax.set_ylim(0, 105)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f%%"))
        ax.legend(loc="upper center", bbox_to_anchor=(0.5,-0.22),
                  ncol=len(cats), fontsize=8, frameon=False, title="Category")
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


def plot_specificity_bar_per_bids(qc_df, output_path, cat_col="sp_category", figsize=None):
    return plot_specificity_bar(qc_df, output_path, cat_col=cat_col,
                                group_by="bids_dir", figsize=figsize)


# ═══════════════════════════════════════════════════════════════════════════════
# FIX 3 — ANESTHESIA COMPARISON PLOT
# ═══════════════════════════════════════════════════════════════════════════════

def _is_awake(anesth_str):
    if not isinstance(anesth_str, str): return True
    s = anesth_str.strip().lower()
    return not s or s in ("nan","none") or any(k in s for k in _AWAKE_KEYWORDS)


def plot_anesthesia_comparison(qc_df, output_path, anesth_map,
                                group_by="species", figsize=None):
    """
    Stacked bar chart: proportion Awake vs Anaesthetised per species or per BIDS.
    Similar layout to plot_specificity_bar.

    Parameters
    ----------
    anesth_map : {(subject, session): anesth_str}
    group_by   : 'species' or 'bids_dir'
    """
    if anesth_map is None or qc_df is None or qc_df.empty:
        warnings.warn("plot_anesthesia_comparison: no anesth_map or qc_df"); return None

    with plt.rc_context(PAPER_RC):
        # Assign condition to each QC row
        df = qc_df.copy()
        def _get_cond(row):
            sub = str(row.get("subject","")).strip()
            ses = str(row.get("session","1")).strip()
            an  = anesth_map.get((sub, ses), anesth_map.get((sub,"1"), ""))
            return "Awake" if _is_awake(an) else "Anaesthetised"

        df["_condition"] = df.apply(_get_cond, axis=1)

        if group_by == "bids_dir" and "bids_dir" in df.columns:
            group_vals = sorted(df["bids_dir"].unique())
            col_src    = df["bids_dir"]
        else:
            group_vals = _phylo_sort(df["species"].unique()) \
                         if "species" in df.columns else list(df.index.unique())
            col_src    = df["species"] if "species" in df.columns \
                         else pd.Series(df.index, index=df.index)

        conds  = ["Awake", "Anaesthetised"]
        colors = {"Awake": "#0072B2", "Anaesthetised": "#D55E00"}

        counts = {}
        for g in group_vals:
            gdf   = df[col_src == g]["_condition"]
            total = len(gdf)
            counts[g] = {c: 100.0*(gdf==c).sum()/total for c in conds} \
                        if total > 0 else {c: 0.0 for c in conds}

        n_g  = len(group_vals)
        w, h = figsize or (max(5, n_g*1.2), 4.5)
        fig, ax = plt.subplots(figsize=(w, h))
        bottoms = np.zeros(n_g); x = np.arange(n_g)
        for cond in conds:
            vals = np.array([counts[g][cond] for g in group_vals])
            ax.bar(x, vals, 0.7, bottom=bottoms, color=colors[cond], label=cond,
                   edgecolor="white", linewidth=0.5)
            for xi, (v, b) in enumerate(zip(vals, bottoms)):
                if v > 5:
                    ax.text(xi, b+v/2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")
            bottoms += vals

        ax.set_xticks(x)
        ax.set_xticklabels(group_vals, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel("% subjects", fontsize=10)
        ax.set_ylim(0, 105)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f%%"))
        title = f"Awake vs Anaesthetised — by {group_by}"
        ax.set_title(title, fontweight="bold")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5,-0.22),
                  ncol=2, fontsize=9, frameon=False)
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
# §D  CORRELATION MATRIX FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

def _draw_corr_matrix_column(fig, gs, col_idx, label, data_dict, n_rows=3,
                              norm_r=None, norm_var=None,
                              cmap_r="RdBu_r", cmap_var="viridis", first_col=True):
    mean  = data_dict["mean"]; pval = data_dict["pval"]; var = data_dict["var"]
    n     = data_dict["n"];    rois = data_dict["rois"]
    short = [r.replace("L_","L ").replace("R_","R ")[:16] for r in rois]
    n_roi = len(rois)
    ax0   = fig.add_subplot(gs[0, col_idx])
    ax0.imshow(mean, cmap=cmap_r, norm=norm_r, aspect="auto", interpolation="nearest")
    ax0.set_title(f"{label}\n(n={n})", fontweight="bold", fontsize=8)
    if first_col:
        ax0.set_ylabel("Mean r", fontsize=8)
        ax0.set_yticks(range(n_roi)); ax0.set_yticklabels(short, fontsize=4)
    else: ax0.set_yticks([])
    ax0.set_xticks([])
    ax1 = fig.add_subplot(gs[1, col_idx])
    pval_disp = np.clip(pval if pval is not None else np.ones_like(mean), 1e-10, 1.0)
    ax1.imshow(-np.log10(pval_disp), cmap="hot_r", aspect="auto",
               vmin=0, vmax=4, interpolation="nearest")
    if first_col:
        ax1.set_ylabel("-log10(p)", fontsize=8)
        ax1.set_yticks(range(n_roi)); ax1.set_yticklabels(short, fontsize=4)
    else: ax1.set_yticks([])
    ax1.set_xticks([])
    ax2 = fig.add_subplot(gs[2, col_idx])
    ax2.imshow(var, cmap=cmap_var, norm=norm_var, aspect="auto", interpolation="nearest")
    if first_col:
        ax2.set_ylabel("Variance", fontsize=8)
        ax2.set_yticks(range(n_roi)); ax2.set_yticklabels(short, fontsize=4)
    else: ax2.set_yticks([])
    ax2.set_xticks([])


def _make_corr_norms(corr_data, vmin=-0.5, vmax=0.5):
    all_var = np.concatenate([d["var"].ravel() for d in corr_data.values()])
    valid   = all_var[~np.isnan(all_var)]
    vmax_var = float(np.nanpercentile(valid, 95)) if valid.size > 0 else 0.01
    return (TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax), Normalize(0, vmax_var))


def plot_corr_matrix_combo(species_config, output_path,
                            atlas_name="EDNIxCSC", atlas_level=3, use_lr=False,
                            figsize=None, vmin=-0.5, vmax=0.5, pval_thresh=0.05):
    csv_dir   = opj(os.path.dirname(output_path), "corr_matrix_csv")
    corr_data = collect_corr_matrices(species_config, atlas_name, atlas_level,
                                       use_lr, save_csv_dir=csv_dir)
    if not corr_data: warnings.warn("plot_corr_matrix_combo: no data"); return None
    species_order = [s for s in _phylo_sort(list(corr_data.keys())) if s in corr_data]
    n_sp = len(species_order)
    norm_r, norm_var = _make_corr_norms(corr_data, vmin, vmax)
    w, h = figsize or (max(8, n_sp*3.5), 3*3.8)
    fig  = plt.figure(figsize=(w, h))
    gs   = gridspec.GridSpec(3, n_sp+1, width_ratios=[1.0]*n_sp+[0.05],
                             wspace=0.04, hspace=0.15)
    for col, sp in enumerate(species_order):
        _draw_corr_matrix_column(fig, gs, col, sp, corr_data[sp],
                                  norm_r=norm_r, norm_var=norm_var, first_col=(col==0))
    for ri, (norm, cmap, label) in enumerate([
        (norm_r,"RdBu_r","Pearson r"),
        (Normalize(0,4),"hot_r","-log10(p)"),
        (norm_var,"viridis","Variance"),
    ]):
        cbar_ax = fig.add_subplot(gs[ri,-1])
        fig.colorbar(_cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, label=label)
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=200)
    plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


def plot_corr_matrix_per_bids(species_config, output_path,
                               atlas_name="EDNIxCSC", atlas_level=3, use_lr=False,
                               figsize=None, vmin=-0.5, vmax=0.5, bids_col=None):
    """
    FIX 4: collect_corr_matrices_per_bids now uses the same common-ROI logic
    as the per-species version, so the ROI sets match.
    """
    csv_dir   = opj(os.path.dirname(output_path), "corr_matrix_csv_per_bids")
    corr_data = collect_corr_matrices_per_bids(species_config, atlas_name,
                                                atlas_level, use_lr,
                                                save_csv_dir=csv_dir)
    if not corr_data: warnings.warn("plot_corr_matrix_per_bids: no data"); return None
    bids_order = sorted(corr_data.keys())
    n_bd = len(bids_order)
    norm_r, norm_var = _make_corr_norms(corr_data, vmin, vmax)
    _bids_col = bids_col or {}
    w, h = figsize or (max(8, n_bd*3.5), 3*3.8)
    fig  = plt.figure(figsize=(w, h))
    gs   = gridspec.GridSpec(3, n_bd+1, width_ratios=[1.0]*n_bd+[0.05],
                             wspace=0.04, hspace=0.15)
    for col, bd_lbl in enumerate(bids_order):
        d     = corr_data[bd_lbl]
        label = f"{bd_lbl}\n({d.get('species','')})"
        _draw_corr_matrix_column(fig, gs, col, label, d, norm_r=norm_r,
                                  norm_var=norm_var, first_col=(col==0))
        if col < len(fig.axes):
            col_edge = _bids_col.get(bd_lbl, "#333333")
            for spine in fig.axes[col*3].spines.values():
                spine.set_edgecolor(col_edge); spine.set_linewidth(2)
    for ri, (norm, cmap, label) in enumerate([
        (norm_r,"RdBu_r","Pearson r"),
        (Normalize(0,4),"hot_r","-log10(p)"),
        (norm_var,"viridis","Variance"),
    ]):
        cbar_ax = fig.add_subplot(gs[ri,-1])
        fig.colorbar(_cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, label=label)
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=200)
    plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


def plot_corr_matrix_sig_norm(species_config, output_path,
                               atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
                               pval_thresh=0.05, figsize=None, df_qc=None):
    """
    2-row normalised FC matrix figure.

    Row 0: All subjects
    Row 1: Specific subjects only

    Normalisation strategy
    ----------------------
    When sp_specific_correlation is available from df_qc:
        norm_mat = 0.5 * mat / sp_specific_corr
        → the specific homotopic pair maps to 0.5 on the colour scale
        → all other pairs are shown relative to that reference value

    When sp_specific_correlation is NOT available (no QC data):
        z-score each matrix independently:
        norm_mat = (mat - mean(mat)) / std(mat)
        → each species shown relative to its own distribution
        → no external reference needed; homotopic pairs appear as +ve outliers
        → avoids meaningless division by mean|r| which distorts geometry

    No mean|r| normalisation is used (it compresses small values and inflates
    noise in weak-connectivity species like rodents).
    """
    csv_dir   = opj(os.path.dirname(output_path), "corr_matrix_csv")
    corr_data_all = collect_corr_matrices(species_config, atlas_name, atlas_level,
                                          use_lr, save_csv_dir=csv_dir)
    if not corr_data_all:
        warnings.warn("plot_corr_matrix_sig_norm: no data"); return None

    # Collect Specific-only matrices using _build_keep_list (same session
    # format handling as category plots — fixes Human/Marmoset n=0)
    cfg_spec = {}
    for sp, cfg in species_config.items():
        keep = _build_keep_list(df_qc, sp, "Specific")
        cfg_spec[sp] = {"bids_dirs":    cfg.get("bids_dirs", []),
                        "list_to_keep": keep, "list_to_remove": []}
        if keep not in ([], [("__NONE__","0")]):
            print(f"  [sig_norm] {sp}/Specific: {len(keep)} (sub,ses) pairs")
    corr_data_spec = collect_corr_matrices(cfg_spec, atlas_name, atlas_level, use_lr)

    # sp_specific_correlation lookup (mean across subjects per species)
    sp_specific = {}
    if df_qc is not None and "sp_specific_correlation" in df_qc.columns:
        for sp, grp in df_qc.groupby("species"):
            vals = grp["sp_specific_correlation"].dropna().values
            if len(vals):
                sp_specific[sp] = float(np.mean(vals))

    def _norm_mat(mean, sp):
        """
        Normalise matrix for display.
        Returns (normalised_matrix, description_string, vabs) where vabs is
        the recommended colour scale half-range.
        """
        mask = ~np.eye(mean.shape[0], dtype=bool)
        if sp in sp_specific and sp_specific[sp] > 0:
            sc   = sp_specific[sp]
            mat  = np.clip(mean * 0.5 / sc, -1.5, 1.5)
            tag  = f"×0.5/sp_r ({sc:.2f})"
            vabs = 1.5
        else:
            # z-score within the off-diagonal distribution
            mu  = float(np.nanmean(mean[mask]))
            sd  = float(np.nanstd(mean[mask]))
            mat = (mean - mu) / sd if sd > 1e-6 else mean - mu
            mat = np.clip(mat, -3, 3)
            tag  = "z-score"
            vabs = 3.0
        return mat, tag, vabs

    species_order = [s for s in _phylo_sort(list(corr_data_all.keys()))
                     if s in corr_data_all]
    n_sp = len(species_order)
    cmap = "RdBu_r"
    w, h = figsize or (max(6, n_sp*3.5), 2*3.8)
    fig  = plt.figure(figsize=(w, h))
    gs   = gridspec.GridSpec(2, n_sp+1, width_ratios=[1.0]*n_sp+[0.05],
                             wspace=0.04, hspace=0.20)

    for row_i, (corr_data, row_label) in enumerate([
        (corr_data_all,  "All"),
        (corr_data_spec, "Specific"),
    ]):
        for col, sp in enumerate(species_order):
            ax = fig.add_subplot(gs[row_i, col])
            if sp not in corr_data or corr_data[sp]["n"] == 0:
                ax.set_facecolor("#f2f2f2")
                ax.set_xticks([]); ax.set_yticks([])
                ax.spines[:].set_visible(False)
                ax.text(0.5, 0.5, "n=0", ha="center", va="center",
                        fontsize=9, color="#bbbbbb", transform=ax.transAxes)
                if row_i == 0: ax.set_title(sp, fontweight="bold", fontsize=9)
                if col == 0:   ax.set_ylabel(row_label, fontsize=8)
                continue

            d           = corr_data[sp]
            mean_mat    = d["mean"]; rois = d["rois"]; n = d["n"]
            short       = [r.replace("L_","L ").replace("R_","R ")[:16] for r in rois]
            mat_n, tag, vabs = _norm_mat(mean_mat, sp)
            norm        = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)

            im = ax.imshow(mat_n, cmap=cmap, norm=norm, aspect="equal",
                           interpolation="nearest")
            if row_i == 0:
                ax.set_title(f"{sp} (n={n})", fontweight="bold", fontsize=9)
            if col == 0:
                ax.set_ylabel(f"{row_label}\n({tag})", fontsize=7)
                ax.set_yticks(range(len(rois)))
                ax.set_yticklabels(short, fontsize=3)
            else:
                ax.set_yticks([])
            ax.set_xticks([])
            ax.text(0.02, 0.98, f"n={n}", transform=ax.transAxes,
                    fontsize=6, va="top", ha="left", color="#444")

    # Single shared colorbar (approximate — each species may use different vabs)
    cbar_ax = fig.add_subplot(gs[:, -1])
    fig.colorbar(_cm.ScalarMappable(
                     norm=TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5), cmap=cmap),
                 cax=cbar_ax, label="Normalised r")
    cbar_ax.yaxis.set_tick_params(labelsize=6)

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  [plot] {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
# FC REPORT BY CATEGORY
# ═══════════════════════════════════════════════════════════════════════════════

def _build_keep_list(qc_df, sp, cat):
    """
    Build list_to_keep of (subject, session) tuples for a given species and category.
    Returns [] (no filter) for cat="all".
    Returns [("__NONE__","0")] sentinel when no subjects exist (blocks all results).
    Adds multiple session format variants to handle zero-padding mismatches.
    """
    if cat == "all" or qc_df is None or "sp_category" not in qc_df.columns:
        return []

    sp_col   = qc_df["species"] if "species" in qc_df.columns \
               else pd.Series([sp]*len(qc_df), index=qc_df.index)
    cat_mask = (sp_col == sp) & (qc_df["sp_category"] == cat)
    sub_qc   = qc_df[cat_mask]

    if sub_qc.empty:
        return [("__NONE__", "0")]  # sentinel — matches nothing

    keep = set()
    for _, row in sub_qc.iterrows():
        sub = str(row["subject"]).strip()
        ses = str(row.get("session", "1")).strip()
        keep.add((sub, ses))
        # Add numeric variants to handle zero-padding mismatches
        if ses.isdigit():
            keep.add((sub, ses.zfill(2)))
            keep.add((sub, ses.lstrip("0") or "1"))
        # Also add the subject with no session restriction (ses="*") —
        # NOT done here; instead we use a subject-only fallback below
    return list(keep)


def _collect_corr_for_category(species_config, qc_df, cat,
                                atlas_name, atlas_level, use_lr):
    """
    Collect correlation matrices for one sp_category subset.
    cat = "all" | "Specific" | "Unspecific" | "No"
    """
    cfg_sub = {}
    for sp, cfg in species_config.items():
        keep = _build_keep_list(qc_df, sp, cat)
        cfg_sub[sp] = {
            "bids_dirs":    cfg.get("bids_dirs", []),
            "list_to_keep": keep,
            "list_to_remove": [],
        }
        if cat != "all" and keep != [("__NONE__", "0")]:
            print(f"  [keep] {sp}/{cat}: {len(keep)} (sub,ses) pairs")

    try:
        cd = collect_corr_matrices(cfg_sub, atlas_name, atlas_level, use_lr)
    except Exception as e:
        warnings.warn(f"_collect_corr_for_category [{cat}]: {e}")
        cd = {}
    return cd


def _fc_category_figure(cat_corr, qc_df, output_path, species_or_bids_list,
                         groupby, bids_col, anesth_map, figsize,
                         bar_metrics, cats, cat_label, _SPEC_COLORS_local):
    """
    Shared worker for per-species and per-bids FC report by category.

    groupby: 'species' | 'bids'
    species_or_bids_list: ordered list of column labels (species names or bids labels)
    cat_corr: {cat: {label: corr_data_dict}}
    """
    n_grp  = len(species_or_bids_list)
    n_rows = len(cats) + len(bar_metrics)

    # Shared colour scale
    all_means = []
    for cd in cat_corr.values():
        for d in cd.values():
            v = d["mean"].ravel(); all_means.append(v[np.isfinite(v)])
    vmax_r = max(float(np.nanpercentile(np.abs(np.concatenate(all_means)), 97))
                 if all_means else 0.5, 0.05)
    norm_r = TwoSlopeNorm(vmin=-vmax_r, vcenter=0, vmax=vmax_r)
    cmap_r = "RdBu_r"

    w, h = figsize or (max(8, n_grp*3.5), n_rows*3.2)
    fig  = plt.figure(figsize=(w, h))
    gs   = gridspec.GridSpec(n_rows, n_grp+1,
                             width_ratios=[1.0]*n_grp + [0.04],
                             wspace=0.05, hspace=0.32)

    # Matrix rows
    for r_i, cat in enumerate(cats):
        cd = cat_corr.get(cat, {})
        for c_i, grp_lbl in enumerate(species_or_bids_list):
            ax  = fig.add_subplot(gs[r_i, c_i])
            spd = cd.get(grp_lbl)
            n   = spd["n"] if spd else 0
            if not spd or n == 0:
                ax.set_facecolor("#f2f2f2"); ax.set_xticks([]); ax.set_yticks([])
                ax.spines[:].set_visible(False)
                ax.text(0.5,0.5,"n=0",ha="center",va="center",
                        fontsize=10,color="#bbbbbb",transform=ax.transAxes)
            else:
                mean  = spd["mean"]; rois = spd["rois"]
                short = [r.replace("L_","L ").replace("R_","R ")[:14] for r in rois]
                ax.imshow(mean, cmap=cmap_r, norm=norm_r, aspect="equal",
                          interpolation="nearest")
                if c_i == 0:
                    ax.set_yticks(range(len(rois)))
                    ax.set_yticklabels(short, fontsize=3)
                else:
                    ax.set_yticks([])
                ax.set_xticks([])
                ax.text(0.02, 0.98, f"n={n}", transform=ax.transAxes,
                        fontsize=6, va="top", ha="left", color="#444")
            col_lbl = _SPEC_COLORS_local.get(cat, "#444") if cat != "all" else "#444"
            if c_i == 0:
                ax.set_ylabel(cat_label.get(cat, cat), fontsize=8,
                              color=col_lbl, fontweight="bold")
            if r_i == 0:
                ax.set_title(grp_lbl, fontweight="bold", fontsize=8)

    cbar_ax = fig.add_subplot(gs[:len(cats), -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm_r, cmap=cmap_r),
                 cax=cbar_ax, label="Mean r")

    # Bar rows — only Specific/Unspecific/No (no "All")
    bar_labels = {"net_homotopic_mean": "Homotopic\nmean r (L↔R)",
                  "net_cross_mean":     "Cross-network\nmean r",
                  "net_inter_mean":     "Inter-hemi\nmean r",
                  "func_mean_fd":       "Mean FD"}
    bar_cats   = [c for c in cats if c != "all"]  # exclude "All" from bars
    bar_colors = {c: _SPEC_COLORS_local.get(c, "#aaaaaa") for c in bar_cats}

    for bm_i, metric in enumerate(bar_metrics):
        row_idx = len(cats) + bm_i
        for c_i, grp_lbl in enumerate(species_or_bids_list):
            ax = fig.add_subplot(gs[row_idx, c_i])
            if qc_df is None or metric not in qc_df.columns:
                ax.set_visible(False); continue

            if groupby == "species":
                grp_df = qc_df[qc_df["species"]==grp_lbl] \
                         if "species" in qc_df.columns else qc_df
            else:
                grp_df = qc_df[qc_df["bids_dir"]==grp_lbl] \
                         if "bids_dir" in qc_df.columns else qc_df

            for x_i, cat in enumerate(bar_cats):
                if "sp_category" in grp_df.columns:
                    sub = grp_df[grp_df["sp_category"]==cat]
                else:
                    sub = pd.DataFrame()
                vals  = sub[metric].dropna().values if not sub.empty else np.array([])
                color = bar_colors[cat]
                if len(vals) > 0:
                    _violin_strip_quartiles(ax, vals, x_i, color, width=0.55)
                    jit  = _rng.uniform(-0.13, 0.13, len(vals))
                    subs = sub.loc[sub[metric].notna(),"subject"].values \
                           if "subject" in sub.columns else None
                    sess = sub.loc[sub[metric].notna(),"session"].values \
                           if "session" in sub.columns else None
                    _scatter_with_anesth(ax, x_i+jit, vals, subs, sess,
                                          anesth_map, color, _rng)
                    ax.text(x_i, ax.get_ylim()[0] if hasattr(ax,'_ylim_') else 0,
                            f"n={len(vals)}", ha="center", va="top",
                            fontsize=5, color="#666")
                else:
                    ax.text(x_i, 0, "n=0", ha="center", va="center",
                            fontsize=7, color="#aaaaaa",
                            transform=ax.transData)

            ax.set_xticks(range(len(bar_cats)))
            ax.set_xticklabels([cat_label.get(c, c) for c in bar_cats],
                                rotation=35, ha="right", fontsize=7)
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
            if c_i == 0:
                ax.set_ylabel(bar_labels.get(metric, metric.replace("_"," ")),
                              fontsize=7, rotation=90)

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  [plot] {output_path}")
    return output_path


def plot_fc_report_by_category(
    species_config, qc_df, output_path,
    atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
    bids_col=None, anesth_map=None, figsize=None,
):
    """
    N_species columns × (4 matrix + 2 bar) rows:
      row 0-3: raw mean FC matrices — All / Specific / Unspecific / No
      row 4:   net_inter_mean violin per category
      row 5:   func_gcor violin per category

    n shown on each panel. Grey placeholder for n=0.
    """
    cats      = ["all"] + _SPEC_CATS
    cat_label = {"all": "All", "Specific": "Specific",
                 "Unspecific": "Unspecific", "No": "No"}

    cat_corr = {}
    for cat in cats:
        cat_corr[cat] = _collect_corr_for_category(
            species_config, qc_df, cat, atlas_name, atlas_level, use_lr)
        n_total = sum(d["n"] for d in cat_corr[cat].values()) if cat_corr[cat] else 0
        print(f"  [fc_by_cat/species] {cat}: n_total={n_total}")

    species_order = [s for s in _phylo_sort(list(species_config.keys()))
                     if s in cat_corr.get("all", {})]
    if not species_order:
        warnings.warn("plot_fc_report_by_category: no data"); return None

    # Prefer homotopic/cross metrics (new); fall back to inter_mean/mean_fd
    bar_metrics = [m for m in ["net_homotopic_mean", "net_cross_mean"]
                   if qc_df is not None and m in qc_df.columns]
    if not bar_metrics:
        bar_metrics = [m for m in ["net_inter_mean", "func_mean_fd"]
                       if qc_df is not None and m in qc_df.columns]

    with plt.rc_context(PAPER_RC):
        return _fc_category_figure(
            cat_corr, qc_df, output_path, species_order,
            groupby="species", bids_col=bids_col, anesth_map=anesth_map,
            figsize=figsize, bar_metrics=bar_metrics,
            cats=cats, cat_label=cat_label,
            _SPEC_COLORS_local=_SPEC_COLORS)


def plot_fc_report_by_category_per_bids(
    species_config, qc_df, output_path,
    atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
    bids_col=None, anesth_map=None, figsize=None,
):
    """
    Same as plot_fc_report_by_category but columns = BIDS directories.
    Uses collect_corr_matrices_per_bids filtered per category.
    """
    cats      = ["all"] + _SPEC_CATS
    cat_label = {"all": "All", "Specific": "Specific",
                 "Unspecific": "Unspecific", "No": "No"}

    cat_corr_bids = {}
    for cat in cats:
        # Build per-category species_config for per-bids collection
        if cat == "all":
            cfg_sub = species_config
        else:
            cfg_sub = {}
            for sp, cfg in species_config.items():
                if qc_df is None or "sp_category" not in qc_df.columns:
                    keep = []
                else:
                    sp_col  = qc_df["species"] if "species" in qc_df.columns \
                              else pd.Series([sp]*len(qc_df), index=qc_df.index)
                    sub_qc  = qc_df[(sp_col==sp) & (qc_df["sp_category"]==cat)]
                    keep    = list(zip(
                        sub_qc["subject"].astype(str),
                        sub_qc.get("session",
                          pd.Series(["1"]*len(sub_qc), index=sub_qc.index)
                        ).astype(str))) if not sub_qc.empty else [("__NONE__","0")]
                cfg_sub[sp] = {"bids_dirs": cfg.get("bids_dirs",[]),
                               "list_to_keep": keep, "list_to_remove": []}
        try:
            cd = collect_corr_matrices_per_bids(cfg_sub, atlas_name, atlas_level, use_lr)
        except Exception as e:
            warnings.warn(f"plot_fc_report_by_category_per_bids [{cat}]: {e}")
            cd = {}
        cat_corr_bids[cat] = cd
        n_total = sum(d["n"] for d in cd.values()) if cd else 0
        print(f"  [fc_by_cat/bids] {cat}: n_total={n_total}")

    # Determine BIDS order from "all" category
    bids_order = sorted(cat_corr_bids.get("all", {}).keys())
    if not bids_order:
        warnings.warn("plot_fc_report_by_category_per_bids: no data"); return None

    # Prefer homotopic/cross metrics (new); fall back to inter_mean/mean_fd
    bar_metrics = [m for m in ["net_homotopic_mean", "net_cross_mean"]
                   if qc_df is not None and m in qc_df.columns]
    if not bar_metrics:
        bar_metrics = [m for m in ["net_inter_mean", "func_mean_fd"]
                       if qc_df is not None and m in qc_df.columns]

    with plt.rc_context(PAPER_RC):
        return _fc_category_figure(
            cat_corr_bids, qc_df, output_path, bids_order,
            groupby="bids", bids_col=bids_col, anesth_map=anesth_map,
            figsize=figsize, bar_metrics=bar_metrics,
            cats=cats, cat_label=cat_label,
            _SPEC_COLORS_local=_SPEC_COLORS)

def plot_brain_scaling(df, metric_col, metric_label, output_path,
                        y_region, x_region, atlas_level=1):
    from scipy import stats as _stats
    plot_df = _filter_atlas_level(df.copy(), atlas_level)
    species_order = _phylo_sort(plot_df["species"].unique())
    sp_colors     = {sp: PALETTE[i % len(PALETTE)] for i, sp in enumerate(species_order)}
    id_cols = [c for c in ("subject","session","bids_dir","species") if c in plot_df.columns]
    def _extract(region):
        rdf = plot_df[_df_region_mask(plot_df["region"], [region])]
        return rdf.groupby(id_cols)[metric_col].mean().reset_index()
    y_df   = _extract(y_region).rename(columns={metric_col: "y_val"})
    x_df   = _extract(x_region).rename(columns={metric_col: "x_val"})
    merged = y_df.merge(x_df, on=id_cols).dropna(subset=["x_val","y_val"])
    if merged.empty: warnings.warn("plot_brain_scaling: no data"); return None
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(5.5, 5))
        for i, sp in enumerate(species_order):
            sp_df = merged[merged["species"]==sp]
            if sp_df.empty: continue
            ax.scatter(sp_df["x_val"], sp_df["y_val"], color=sp_colors[sp],
                       s=40, zorder=5, alpha=0.85, edgecolors="k", linewidths=0.5, label=sp)
        xv = np.clip(merged["x_val"].values.astype(float), 1e-12, None)
        yv = np.clip(merged["y_val"].values.astype(float), 1e-12, None)
        lx, ly = np.log10(xv), np.log10(yv)
        valid  = np.isfinite(lx) & np.isfinite(ly)
        if valid.sum() >= 3:
            lx_v, ly_v = lx[valid], ly[valid]
            slope, intercept, r, p, _ = _stats.linregress(lx_v, ly_v)
            x_fit  = np.linspace(lx_v.min(), lx_v.max(), 200)
            y_fit  = intercept + slope*x_fit
            n_v    = valid.sum(); t_crit = _stats.t.ppf(0.975, df=n_v-2)
            x_mean = lx_v.mean(); ssx = np.sum((lx_v-x_mean)**2)
            rv     = np.sum((ly_v-(intercept+slope*lx_v))**2)/(n_v-2)
            se     = np.sqrt(rv*(1+1/n_v+(x_fit-x_mean)**2/max(ssx,1e-12)))
            ax.plot(10**x_fit, 10**y_fit, "k--", lw=1.5, alpha=0.8, zorder=3)
            ax.fill_between(10**x_fit, 10**(y_fit-t_crit*se),
                            10**(y_fit+t_crit*se), color="k", alpha=0.12, zorder=2)
            p_str = f"{p:.3f}" if p >= 0.001 else "<0.001"
            ax.text(0.05, 0.95, f"slope={slope:.2f}  r={r:.2f}  p={p_str}",
                    transform=ax.transAxes, fontsize=8, va="top", color="#333")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
        ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())
        ax.set_xlabel(f"{x_region}  {metric_label}", fontsize=9)
        ax.set_ylabel(f"{y_region}  {metric_label}", fontsize=9)
        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp) for sp in species_order]
        plt.tight_layout()
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5,-0.06),
                   ncol=min(len(handles),7), fontsize=8, frameon=False)
        fig.subplots_adjust(bottom=0.18)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig); print(f"  [plot] {output_path}")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
# §F  make_all_figures
# ═══════════════════════════════════════════════════════════════════════════════

def make_all_figures(data, species_config, all_bids, fig_dir, plot_regions,
                      atlas_level=1, atlas_name="EDNIxCSC",
                      corr_atlas_level=2, corr_use_lr=True,
                      anesth_map=None, global_bids_col=None,
                      average_sessions=False, trimouse_xlsx=None):
    df_surface   = data.get("surface")
    df_thickness = data.get("thickness")
    df_qc        = data.get("qc")

    # Bilateral volume aggregation (always on — sums L+R per subject-region).
    # Session averaging is optional; off by default so each session/run
    # appears as a separate dot (pass average_sessions=True to collapse).
    df_volume = _make_bilateral_volume(data.get("volume"))
    if average_sessions:
        df_volume    = _mean_per_subject(df_volume,    "volume_mm3")
        df_surface   = _mean_per_subject(df_surface,   "surface_area_mm2")
        df_thickness = _mean_per_subject(df_thickness, "thickness_mm")

    if global_bids_col is None:
        all_bd = sorted(set(_bids_label(bd) for _, bd in all_bids))
        global_bids_col = {bd: PALETTE[i%len(PALETTE)] for i,bd in enumerate(all_bd)}

    hemis_surface   = _hemis_for(df_surface)
    hemis_thickness = _hemis_for(df_thickness)
    hemis_volume    = ["bilateral"]   # always bilateral after fix

    morph_dir   = opj(fig_dir, "cross_species", "morphometry_recap")
    combo_dir   = opj(fig_dir, "cross_species", "combo")
    qc_dir      = opj(fig_dir, "cross_species", "qc_recap")
    scaling_dir = opj(fig_dir, "cross_species", "brain_scaling")
    for d in [morph_dir, combo_dir, qc_dir, scaling_dir]:
        os.makedirs(d, exist_ok=True)

    # ── Volume: bilateral only ────────────────────────────────────────────────
    for log_s, suffix in [(False,"raw"),(True,"log")]:
        plot_cross_species_dots(
            df_volume, "volume_mm3", "Volume (mm³)", plot_regions,
            opj(morph_dir, f"volume_bilateral_{suffix}.png"),
            hemisphere="bilateral", log_scale=log_s, atlas_level=atlas_level,
            bids_col=global_bids_col, anesth_map=anesth_map)
    plot_cross_species_dots(
        df_volume, "volume_mm3", "Volume (mm³)", plot_regions,
        opj(morph_dir, "volume_bilateral_norm.png"),
        hemisphere="bilateral", log_scale=False, normalise_by_brain=True,
        df_vol_for_norm=df_volume, atlas_level=atlas_level,
        bids_col=global_bids_col, anesth_map=anesth_map)

    # ── Surface + thickness: all available hemispheres ────────────────────────
    for hemi in sorted(set(hemis_surface + hemis_thickness)):
        for df, metric_col, metric_label, mod, hemis in [
            (df_surface,   "surface_area_mm2", "Surface area (mm²)",      "surface",   hemis_surface),
            (df_thickness, "thickness_mm",      "Cortical thickness (mm)", "thickness", hemis_thickness),
        ]:
            if df is None or df.empty or hemi not in hemis: continue
            for log_s, suffix in [(False,"raw"),(True,"log")]:
                plot_cross_species_dots(
                    df, metric_col, metric_label, plot_regions,
                    opj(morph_dir, f"{mod}_{hemi}_{suffix}.png"),
                    hemisphere=hemi, log_scale=log_s, atlas_level=atlas_level,
                    bids_col=global_bids_col, anesth_map=anesth_map)
            if df_volume is not None:
                plot_cross_species_dots(
                    df, metric_col, metric_label, plot_regions,
                    opj(morph_dir, f"{mod}_{hemi}_norm.png"),
                    hemisphere=hemi, log_scale=False, normalise_by_brain=True,
                    df_vol_for_norm=df_volume, atlas_level=atlas_level,
                    bids_col=global_bids_col, anesth_map=anesth_map)

    # ── Combo ─────────────────────────────────────────────────────────────────
    for log_s, norm_s, suffix in [(False,False,"raw"),(True,False,"log"),(False,True,"norm")]:
        plot_combo(df_surface, df_volume, df_thickness,
                   opj(combo_dir, f"combo_bilateral_{suffix}.png"),
                   regions=plot_regions, hemisphere="bilateral",
                   log_scale=log_s, normalise_by_brain=norm_s,
                   atlas_level=atlas_level, bids_col=global_bids_col,
                   anesth_map=anesth_map)

    # ── QC ───────────────────────────────────────────────────────────────────
    if df_qc is not None and not df_qc.empty:
        plot_cross_species_qc(
            df_qc, opj(qc_dir, "QC_cross_species_anat_and_func.png"),
            bids_col=global_bids_col, anesth_map=anesth_map)
        plot_qc_fc_report(
            df_qc, opj(qc_dir, "QC_cross_species_network_report.png"),
            bids_col=global_bids_col, anesth_map=anesth_map)
        plot_specificity_bar(df_qc, opj(qc_dir,"QC_specificity_by_species.png"),
                             group_by="species")
        plot_specificity_bar_per_bids(df_qc, opj(qc_dir,"QC_specificity_by_bids.png"))

        # Anesthesia comparison (FIX 3)
        if anesth_map:
            plot_anesthesia_comparison(
                df_qc, opj(qc_dir,"QC_anesthesia_by_species.png"),
                anesth_map=anesth_map, group_by="species")
            if "bids_dir" in df_qc.columns:
                plot_anesthesia_comparison(
                    df_qc, opj(qc_dir,"QC_anesthesia_by_bids.png"),
                    anesth_map=anesth_map, group_by="bids_dir")

        # Trimouse multi-condition anesthesia figure
        if trimouse_xlsx and os.path.exists(str(trimouse_xlsx)):
            plot_trimouse_anesthesia(
                species_config, trimouse_xlsx,
                opj(qc_dir,"QC_trimouse_anesthesia.png"),
                atlas_name=atlas_name, atlas_level=corr_atlas_level,
                use_lr=corr_use_lr, bids_col=global_bids_col,
                anesth_map=anesth_map)

        # FC report by category — per species
        plot_fc_report_by_category(
            species_config, df_qc,
            opj(qc_dir,"QC_fc_report_by_category.png"),
            atlas_name=atlas_name, atlas_level=corr_atlas_level, use_lr=corr_use_lr,
            bids_col=global_bids_col, anesth_map=anesth_map)
        # FC report by category — per BIDS dir
        plot_fc_report_by_category_per_bids(
            species_config, df_qc,
            opj(qc_dir,"QC_fc_report_by_category_per_bids.png"),
            atlas_name=atlas_name, atlas_level=corr_atlas_level, use_lr=corr_use_lr,
            bids_col=global_bids_col, anesth_map=anesth_map)

    # ── Correlation matrices ──────────────────────────────────────────────────
    plot_corr_matrix_combo(
        species_config, opj(combo_dir, f"combo_corr_matrix_lvl{corr_atlas_level}.png"),
        atlas_name=atlas_name, atlas_level=corr_atlas_level, use_lr=corr_use_lr)
    plot_corr_matrix_sig_norm(
        species_config, opj(combo_dir, f"combo_corr_matrix_sig_norm_lvl{corr_atlas_level}.png"),
        atlas_name=atlas_name, atlas_level=corr_atlas_level,
        use_lr=corr_use_lr, df_qc=df_qc)
    plot_corr_matrix_per_bids(
        species_config, opj(combo_dir, f"combo_corr_matrix_per_bids_lvl{corr_atlas_level}.png"),
        atlas_name=atlas_name, atlas_level=corr_atlas_level,
        use_lr=corr_use_lr, bids_col=global_bids_col)

    # ── Brain scaling ─────────────────────────────────────────────────────────
    for df, metric_col, metric_label, mod in [
        (df_surface,   "surface_area_mm2", "Surface area (mm²)",      "surface"),
        (df_volume,    "volume_mm3",        "Volume (mm³)",            "volume"),
        (df_thickness, "thickness_mm",      "Cortical thickness (mm)", "thickness"),
    ]:
        if df is None or df.empty: continue
        plot_brain_scaling(
            df, metric_col, metric_label,
            opj(scaling_dir, f"scaling_{mod}_lvl{atlas_level}.png"),
            y_region=plot_regions[0] if plot_regions else "Isocortex",
            x_region=plot_regions[1] if len(plot_regions)>1 else "Allocortex",
            atlas_level=atlas_level)

    # ── Per-BIDS ──────────────────────────────────────────────────────────────
    for species, bids_dir in all_bids:
        label    = _bids_label(bids_dir)
        bids_out = opj(fig_dir, "per_bids", f"{species}_{label}")
        os.makedirs(bids_out, exist_ok=True)
        for df, metric_col, metric_label, mod, hemis in [
            (df_surface,   "surface_area_mm2", "Surface area (mm²)",      "surface",   hemis_surface),
            (df_volume,    "volume_mm3",        "Volume (mm³)",            "volume",    ["bilateral"]),
            (df_thickness, "thickness_mm",      "Cortical thickness (mm)", "thickness", hemis_thickness),
        ]:
            if df is None or df.empty: continue
            df_sub = df[(df["species"]==species) & (df["bids_dir"]==label)]
            if df_sub.empty: continue
            for hemi in hemis:
                plot_morphometry_intra_bids(
                    df_sub, metric_col=metric_col, metric_label=metric_label,
                    regions=plot_regions,
                    output_path=opj(bids_out, f"{mod}_{hemi}.png"),
                    species=species, hemisphere=hemi, bids_label=label,
                    atlas_level=atlas_level, bids_col=global_bids_col,
                    anesth_map=anesth_map)
        if df_qc is not None and not df_qc.empty:
            qc_sub = df_qc[(df_qc.get("species", pd.Series())==species) &
                           (df_qc.get("bids_dir", pd.Series())==label)] \
                     if "species" in df_qc.columns and "bids_dir" in df_qc.columns \
                     else df_qc[df_qc["species"]==species] \
                     if "species" in df_qc.columns else df_qc
            if not qc_sub.empty:
                plot_qc_dashboard(qc_sub, opj(bids_out,"QC_dashboard.png"),
                                  anesth_map=anesth_map)

    print(f"\n  ✓ All figures → {fig_dir}")


# ═══════════════════════════════════════════════════════════════════════════════
# TRIMOUSE ANESTHESIA COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════

def plot_trimouse_anesthesia(
    species_config, trimouse_xlsx_path, output_path,
    atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
    bids_col=None, anesth_map=None, figsize=None,
    n_top_connections=3,
):
    """
    Multi-anesthesia FC comparison for Trimouse (awake / iso / mediso).

    Layout — 3 rows × (n_cond + n_sig + 1) columns:

      Row 0 — Mean FC matrix:
        col 0: awake   col 1: iso   col 2: mediso

      Row 1 — Fisher z-score difference (same scale):
        col 0: awake vs iso   col 1: awake vs mediso
        col 2: combined |z| (avg of both, shows regions most affected by anaesthesia)

      Row 2 — Bar plots for the n_top_connections ROI pairs most affected:
        one grouped bar per top connection showing r per condition (awake/iso/mediso)

    Parameters
    ----------
    n_top_connections : int
        Number of most-affected ROI connections to highlight in bar plots.
    """
    from scipy import stats as _stats

    # ── Load anesthesia assignments from xlsx ──────────────────────────────────
    try:
        df_anesth = pd.read_excel(trimouse_xlsx_path, sheet_name="MRIsessions")
    except Exception as e:
        warnings.warn(f"plot_trimouse_anesthesia: cannot load xlsx: {e}"); return None

    id_col     = "ID"         if "ID"         in df_anesth.columns else df_anesth.columns[0]
    anesth_col = "anest_type" if "anest_type" in df_anesth.columns else None
    ses_col    = "session"    if "session"    in df_anesth.columns else None
    if anesth_col is None:
        warnings.warn("plot_trimouse_anesthesia: no anest_type column"); return None

    def _sub_id(raw):
        s = str(raw).strip()
        return s[4:] if s.startswith("sub-") else s

    cond_map = {}
    for _, row in df_anesth.iterrows():
        sub  = _sub_id(row[id_col])
        ses  = str(int(row[ses_col])) if ses_col and pd.notna(row[ses_col]) else "1"
        cond = str(row[anesth_col]).strip().lower() if pd.notna(row[anesth_col]) else None
        if cond and cond not in ("unknown", "nan"):
            cond_map[(sub, ses)] = cond

    all_conds = sorted(set(cond_map.values()))
    if "awake" in all_conds:
        conditions = ["awake"] + [c for c in all_conds if c != "awake"]
    else:
        conditions = all_conds
    ref_cond = conditions[0]
    other_conds = [c for c in conditions if c != ref_cond]
    n_cond   = len(conditions)
    print(f"  [trimouse] conditions: {conditions}")

    # ── Collect per-condition corr matrices ────────────────────────────────────
    cond_corr = {}
    for cond in conditions:
        keep = [(s, e) for (s, e), c in cond_map.items() if c == cond]
        keep_ext = set(keep)
        for s, e in keep:
            if e.isdigit():
                keep_ext.add((s, e.zfill(2)))
                keep_ext.add((s, e.lstrip("0") or "1"))
        cfg_sub = {sp: {"bids_dirs": cfg.get("bids_dirs", []),
                        "list_to_keep": list(keep_ext), "list_to_remove": []}
                   for sp, cfg in species_config.items()}
        try:
            cd = collect_corr_matrices(cfg_sub, atlas_name, atlas_level, use_lr)
        except Exception as e:
            warnings.warn(f"[trimouse/{cond}]: {e}"); cd = {}
        cond_corr[cond] = cd
        n = sum(d["n"] for d in cd.values()) if cd else 0
        print(f"  [trimouse] {cond}: n={n}")

    # Pick species (prefer one present in all conditions)
    candidate_sp = [sp for sp in species_config
                    if all(sp in cond_corr.get(c, {}) and
                           cond_corr[c][sp]["n"] > 0
                           for c in conditions)]
    sp = candidate_sp[0] if candidate_sp else next(
        (sp for sp in species_config
         if any(sp in cond_corr.get(c, {}) for c in conditions)), None)
    if sp is None:
        warnings.warn("plot_trimouse_anesthesia: no species with data"); return None
    print(f"  [trimouse] species: {sp}")

    # Reference data for ROIs
    ref_d = next((cond_corr[c][sp] for c in conditions
                  if sp in cond_corr.get(c, {}) and cond_corr[c][sp]["n"] > 0),
                 None)
    if ref_d is None:
        warnings.warn("plot_trimouse_anesthesia: no matrix data"); return None
    rois  = ref_d["rois"]
    n_roi = len(rois)
    short = [r.replace("L_","L ").replace("R_","R ")[:14] for r in rois]

    # ── Colour scales ──────────────────────────────────────────────────────────
    all_r = np.concatenate([cond_corr[c][sp]["mean"].ravel()
                            for c in conditions if sp in cond_corr.get(c,{}) and
                            cond_corr[c][sp]["n"] > 0])
    all_r = all_r[np.isfinite(all_r)]
    vmax  = max(float(np.nanpercentile(np.abs(all_r), 97)), 0.05)
    norm_r   = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    norm_sig = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    cmap_r   = "RdBu_r"

    cpal    = ["#0072B2","#D55E00","#009E73","#E69F00","#CC79A7"]
    ccolors = {c: cpal[i % len(cpal)] for i, c in enumerate(conditions)}

    # ── Compute Fisher z-score difference matrices ─────────────────────────────
    z_diffs = {}
    combined_z = None
    for cond2 in other_conds:
        d1 = cond_corr.get(ref_cond, {}).get(sp)
        d2 = cond_corr.get(cond2,    {}).get(sp)
        if d1 and d2 and d1["n"] >= 2 and d2["n"] >= 2:
            r1 = np.clip(d1["mean"], -0.9999, 0.9999)
            r2 = np.clip(d2["mean"], -0.9999, 0.9999)
            n1, n2 = d1["n"], d2["n"]
            se  = np.sqrt(1.0/(n1-3) + 1.0/(n2-3) + 1e-12)
            z_diffs[cond2] = (np.arctanh(r1) - np.arctanh(r2)) / se
            if combined_z is None:
                combined_z = np.abs(z_diffs[cond2])
            else:
                combined_z = np.maximum(combined_z, np.abs(z_diffs[cond2]))

    # ── Find top affected connections (upper triangle, off-diagonal) ────────────
    top_connections = []
    if combined_z is not None:
        mask_triu = np.triu(np.ones((n_roi, n_roi), dtype=bool), k=1)
        flat_z    = combined_z[mask_triu]
        flat_idx  = np.argwhere(mask_triu)
        top_k     = min(n_top_connections, len(flat_idx))
        top_order = np.argsort(flat_z)[::-1][:top_k]
        top_connections = [(flat_idx[i][0], flat_idx[i][1], float(flat_z[i]))
                           for i in top_order]
        print(f"  [trimouse] top {top_k} connections:")
        for i, j, z in top_connections:
            print(f"    {rois[i]} ↔ {rois[j]}: combined |z|={z:.2f}")

    # ── Figure layout ──────────────────────────────────────────────────────────
    # Row 0: mean FC matrices   (n_cond columns)
    # Row 1: sig diff matrices  (len(other_conds) + 1 for combined)
    # Row 2: top connection bars (n_top_connections panels)
    n_sig_cols = len(other_conds) + 1   # iso-diff + mediso-diff + combined
    n_col_max  = max(n_cond, n_sig_cols, n_top_connections)
    cell       = 3.2

    w = n_col_max * cell + 0.6
    h = 3 * cell + 0.5
    w, h = figsize or (w, h)

    fig = plt.figure(figsize=(w, h))
    gs  = gridspec.GridSpec(3, n_col_max + 1,
                            width_ratios=[cell]*n_col_max + [0.25],
                            wspace=0.06, hspace=0.28)

    # ── Row 0: mean r per condition ────────────────────────────────────────────
    for c_i, cond in enumerate(conditions):
        ax  = fig.add_subplot(gs[0, c_i])
        d   = cond_corr.get(cond, {}).get(sp)
        if d and d["n"] > 0:
            ax.imshow(d["mean"], cmap=cmap_r, norm=norm_r, aspect="equal",
                      interpolation="nearest")
            ax.text(0.02, 0.98, f"n={d['n']}", transform=ax.transAxes,
                    fontsize=7, va="top", ha="left",
                    color=ccolors[cond], fontweight="bold")
        else:
            ax.set_facecolor("#f2f2f2")
            ax.text(0.5, 0.5, "n=0", ha="center", va="center",
                    fontsize=9, color="#bbb", transform=ax.transAxes)
        if c_i == 0:
            ax.set_yticks(range(n_roi)); ax.set_yticklabels(short, fontsize=3)
            ax.set_ylabel("Mean r", fontsize=8)
        else:
            ax.set_yticks([])
        ax.set_xticks([])
        ax.set_title(cond, fontsize=9, fontweight="bold",
                     color=ccolors.get(cond, "#444"))
    # Pad unused columns in row 0
    for c_i in range(n_cond, n_col_max):
        fig.add_subplot(gs[0, c_i]).set_visible(False)
    cbar0 = fig.add_subplot(gs[0, -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm_r, cmap=cmap_r), cax=cbar0,
                 label="r"); cbar0.yaxis.set_tick_params(labelsize=6)

    # ── Row 1: diff matrices side by side + combined ───────────────────────────
    diff_panels = [(cond2, z_diffs.get(cond2)) for cond2 in other_conds]
    if combined_z is not None:
        diff_panels.append(("combined |z|", combined_z))
        combined_norm = Normalize(0, 3)
        combined_cmap = "Reds"
    for c_i, (label2, mat2) in enumerate(diff_panels):
        ax = fig.add_subplot(gs[1, c_i])
        if mat2 is not None:
            is_combined = label2 == "combined |z|"
            _norm = combined_norm if is_combined else norm_sig
            _cmap = combined_cmap if is_combined else cmap_r
            ax.imshow(mat2, cmap=_cmap, norm=_norm, aspect="equal",
                      interpolation="nearest")
        else:
            ax.set_facecolor("#f0f0f0")
            ax.text(0.5, 0.5, "n<2", ha="center", va="center",
                    fontsize=8, color="#bbb", transform=ax.transAxes)
        if c_i == 0:
            ax.set_yticks(range(n_roi)); ax.set_yticklabels(short, fontsize=3)
            ax.set_ylabel("Fisher z\n(awake − other)", fontsize=7)
        else:
            ax.set_yticks([])
        ax.set_xticks([])
        title_str = (f"{ref_cond} vs {label2}" if label2 != "combined |z|"
                     else "Combined |z|\n(max across comparisons)")
        ax.set_title(title_str, fontsize=8, fontweight="bold")
    for c_i in range(len(diff_panels), n_col_max):
        fig.add_subplot(gs[1, c_i]).set_visible(False)
    cbar1 = fig.add_subplot(gs[1, -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm_sig, cmap=cmap_r), cax=cbar1,
                 label="Fisher z"); cbar1.yaxis.set_tick_params(labelsize=6)

    # ── Row 2: bar plots of top connections ────────────────────────────────────
    if top_connections:
        for t_i, (roi_i, roi_j, z_val) in enumerate(top_connections):
            if t_i >= n_col_max: break
            ax = fig.add_subplot(gs[2, t_i])
            x  = np.arange(n_cond)
            bw = 0.6
            for c_i, cond in enumerate(conditions):
                d = cond_corr.get(cond, {}).get(sp)
                r = float(d["mean"][roi_i, roi_j]) if d and d["n"] > 0 else np.nan
                ax.bar(c_i, r, bw, color=ccolors[cond], label=cond,
                       alpha=0.85, edgecolor="white")
                if not np.isnan(r):
                    ax.text(c_i, r + 0.005 if r >= 0 else r - 0.015,
                            f"{r:.2f}", ha="center", va="bottom" if r >= 0 else "top",
                            fontsize=7)
            ax.set_xticks(x)
            ax.set_xticklabels(conditions, fontsize=7, rotation=30, ha="right")
            ax.axhline(0, color="#aaa", lw=0.7, ls="--")
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
            roi_lbl_i = rois[roi_i].replace("L_","L ").replace("R_","R ")[:16]
            roi_lbl_j = rois[roi_j].replace("L_","L ").replace("R_","R ")[:16]
            ax.set_title(f"{roi_lbl_i}\n↔ {roi_lbl_j}\n(|z|={z_val:.1f})",
                         fontsize=7, fontweight="bold")
            if t_i == 0:
                ax.set_ylabel("Mean r", fontsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        for t_i in range(len(top_connections), n_col_max):
            fig.add_subplot(gs[2, t_i]).set_visible(False)
    else:
        for t_i in range(n_col_max):
            fig.add_subplot(gs[2, t_i]).set_visible(False)
    fig.add_subplot(gs[2, -1]).set_visible(False)

    # Legend for conditions
    handles = [mpatches.Patch(color=ccolors[c], label=c) for c in conditions]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=n_cond, fontsize=9,
               frameon=False)

    fig.suptitle(f"Anesthesia comparison — {sp}", fontweight="bold", fontsize=12)
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  [plot] {output_path}")
    return output_path