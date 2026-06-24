#!/usr/bin/env python3
"""
ednix_perturbation_and_qc.py
============================
Two analysis modules for the EDNiX paper:

  MODULE 1: Per-BIDS time-series perturbation
  --------------------------------------------
  Runs the FC perturbation analysis (noise, smooth, shift, phase) for EACH
  BIDS dataset independently. Produces:
    - Per-BIDS perturbation matrix + proportion figures
    - A combined bar plot comparing original vs perturbed specificity
      rates across all BIDS datasets side by side

  MODULE 2: QC density plots
  ---------------------------
  Reads per-subject QC metrics (tSNR, SNR anat, SNR func) from df_qc and
  plots pooled KDE distributions following the same style as the threshold
  explorer Fig 01 middle panel (faded per-BIDS curves, bold global density,
  percentile markers).

USAGE
-----
    from ednix_perturbation_and_qc import (
        run_perturbation_all_bids,
        plot_qc_densities,
    )

    # Module 1
    run_perturbation_all_bids(
        df_scored, cats,
        thresh_intra=0.24, thresh_delta=0.1,
        bids_root_template="/scratch2/EDNiX/{species}/{bids_dir}",
        output_dir="figures/perturbation",
    )

    # Module 2
    plot_qc_densities(
        df_qc,
        output_dir="figures/qc_densities",
    )
"""

import os
import warnings
import re
import glob
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.cm as _mpl_cm
from matplotlib.colors import TwoSlopeNorm
from scipy import stats as _stats

opj = os.path.join
ope = os.path.exists

# ─────────────────────────────────────────────────────────────────────────────
# Phylogenetic color palette (must match threshold_explorer)
# ─────────────────────────────────────────────────────────────────────────────
_SPECIES_COLORS = {
    "Mouse":      "#1B4F72",
    "Rat":        "#2980B9",
    "Mouselemur": "#1D6B3F",
    "Marmoset":   "#27AE60",
    "Macaque":    "#82E0AA",
    "Human":      "#196F3D",
    "Dog":        "#E67E22",
    "Cat":        "#F39C12",
    "Bat":        "#8E44AD",
    "Pig":        "#A04000",
}
_PHYLO = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
          "Macaque", "Human", "Dog", "Cat", "Pig"]
_CAT_COLORS = {"Specific": "#009E73", "Unspecific": "#E69F00",
               "No": "#D55E00", "Spurious": "#CC79A7"}
CATS = ["Specific", "Unspecific", "No", "Spurious"]

def _sp_order(df):
    vals = df["species"].unique().tolist() if "species" in df.columns else []
    return [s for s in _PHYLO if s in vals] + [s for s in vals if s not in _PHYLO]

def _sp_color(sp):
    return _SPECIES_COLORS.get(sp, "#888888")

def _save(fig, path):
    if path:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [saved] {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 1 — PER-BIDS PERTURBATION
# ═══════════════════════════════════════════════════════════════════════════════

def _perturb_ts(data, kind, rng, **params):
    """Apply perturbation to 4D fMRI array. See ednix_fc_perturbation.py."""
    if kind is None:
        return data.copy()
    out = data.copy().astype(np.float64)
    mask = np.any(out != 0, axis=-1)
    n_t = out.shape[-1]

    if kind == "noise":
        sf = params.get("sigma_frac", 0.3)
        for idx in zip(*np.where(mask)):
            s = np.std(out[idx])
            if s > 0:
                out[idx] += rng.normal(0, sf * s, n_t)
    elif kind == "smooth":
        w = params.get("window", 5)
        kernel = np.ones(w) / w
        for idx in zip(*np.where(mask)):
            out[idx] = np.convolve(out[idx], kernel, mode='same')
    elif kind == "shift":
        ss = params.get("shift_std", 0.5)
        gs = rng.normal(0, ss, n_t)
        out[mask] += gs[np.newaxis, :]
    elif kind == "phase":
        for idx in zip(*np.where(mask)):
            ft = np.fft.rfft(out[idx])
            ph = rng.uniform(0, 2 * np.pi, len(ft))
            ph[0] = 0
            if n_t % 2 == 0:
                ph[-1] = 0
            out[idx] = np.fft.irfft(np.abs(ft) * np.exp(1j * ph), n=n_t)
    return out


def _reconstruct_matrix(ts_4d, atlas_path, atlas_level, tr, fit_kind="correlation"):
    """Reconstruct FC matrix from 4D data via Nilearn (mirrors EDNiX step 10)."""
    from nilearn.input_data import NiftiLabelsMasker
    from nilearn.connectome import ConnectivityMeasure

    atlas_img = nib.load(atlas_path)
    ad = atlas_img.get_fdata()
    level_data = ad[:, :, :, atlas_level] if ad.ndim == 4 else ad
    atlas_3d = nib.Nifti1Image(level_data, atlas_img.affine, atlas_img.header)
    func_img = nib.Nifti1Image(ts_4d, atlas_img.affine)

    masker = NiftiLabelsMasker(
        labels_img=atlas_3d, detrend=False, smoothing_fwhm=None,
        low_pass=None, high_pass=None, t_r=tr,
        standardize='zscore_sample', memory=None, verbose=0)
    ts = masker.fit_transform(func_img)
    roi_names = [str(l) for l in masker.labels_]

    conn = ConnectivityMeasure(kind=fit_kind, standardize='zscore_sample')
    mat = conn.fit_transform([ts])[0]
    np.fill_diagonal(mat, 0)
    return roi_names, mat


def _find_func_paths(df_scored, species, bids_dir_name, bids_root_template):
    """Find residual fMRI + atlas paths for a specific species/BIDS combo."""
    entries = []
    mask = pd.Series(True, index=df_scored.index)
    if "species" in df_scored.columns:
        mask &= df_scored["species"] == species
    if "bids_dir" in df_scored.columns:
        mask &= df_scored["bids_dir"] == bids_dir_name

    for _, row in df_scored[mask].iterrows():
        sp  = str(row.get("species", ""))
        bd  = str(row.get("bids_dir", ""))
        sub = str(row.get("subject", ""))
        ses = str(row.get("session", ""))
        root = bids_root_template.format(species=sp, bids_dir=bd)
        ses_dir = opj(root, f"sub-{sub}", f"ses-{ses}") if ses else opj(root, f"sub-{sub}")

        res_files = sorted(glob.glob(
            opj(ses_dir, "func", "acpc-func", "postprocessed_rs",
                "*_desc-fMRI_residual.nii.gz")))
        atlas_files = sorted(glob.glob(
            opj(ses_dir, "func", "acpc-func", "labels",
                f"{sub}_seg-EDNIxCSC*_dseg.nii.gz")))

        if res_files and atlas_files:
            entries.append(dict(
                species=sp, bids_dir=bd, subject=sub, session=ses,
                residual_path=res_files[0], atlas_path=atlas_files[0],
                category=str(row.get("fp_category", "")),
            ))
    return entries


def run_perturbation_all_bids(
    df_scored, cats, thresh_intra, thresh_delta,
    bids_root_template="/scratch2/EDNiX/{species}/{bids_dir}",
    atlas_level=2,
    fit_kind="correlation",
    max_runs_per_bids=10,
    seed=42,
    output_dir=None,
    verbose=True,
):
    """
    Run perturbation analysis for EACH BIDS dataset and produce:
      1. Per-BIDS perturbation figure (matrix + proportions)
      2. Combined comparison bar plot (original vs perturbed, all BIDS)

    Returns dict of per-BIDS results.
    """
    try:
        from ednix_threshold_explorer import (
            compute_primary_score_from_matrix, _classify_one)
    except ImportError:
        from Exemples.Study_EDNiX.analysis.ednix_threshold_explorer import (
            compute_primary_score_from_matrix, _classify_one)

    df = df_scored.copy()
    df["fp_category"] = cats.values if hasattr(cats, "values") else cats

    perturbations = [
        ("Original",              None,     {}),
        ("Noise (σ=0.3)",         "noise",  dict(sigma_frac=0.3)),
        ("Smooth (w=5)",          "smooth", dict(window=5)),
        ("Global shift (σ=0.5)",  "shift",  dict(shift_std=0.5)),
        ("Phase randomize",       "phase",  {}),
    ]

    # Identify all (species, bids_dir) combinations
    if "species" not in df.columns or "bids_dir" not in df.columns:
        print("  [perturbation] need species + bids_dir columns")
        return {}

    groups = df.groupby(["species", "bids_dir"]).size().reset_index()
    groups.columns = ["species", "bids_dir", "n_runs"]
    groups = groups.sort_values("n_runs", ascending=False)

    all_bids_results = {}
    rng = np.random.default_rng(seed)

    for _, grp_row in groups.iterrows():
        sp = grp_row["species"]
        bd = grp_row["bids_dir"]
        label = f"{bd} ({sp})"
        print(f"\n{'='*60}")
        print(f"  Perturbation: {label}")
        print(f"{'='*60}")

        entries = _find_func_paths(df, sp, bd, bids_root_template)
        if not entries:
            print(f"  [skip] no residual fMRI found")
            continue
        if len(entries) > max_runs_per_bids:
            entries = entries[:max_runs_per_bids]
            print(f"  capped to {max_runs_per_bids} runs")

        results = {lab: [] for lab, _, _ in perturbations}

        for entry in entries:
            if verbose:
                print(f"    {entry['subject']}/{entry['session']}...", end=" ", flush=True)
            try:
                img = nib.load(entry["residual_path"])
                data_4d = img.get_fdata()
                # squeeze singleton dims (e.g. shape 64,64,18,1,1000)
                while data_4d.ndim > 4:
                    data_4d = np.squeeze(data_4d, axis=3)
                run_tr = float(img.header.get_zooms()[-1])
                if run_tr < 0.01 or run_tr > 30:
                    run_tr = 2.0
            except Exception as e:
                if verbose:
                    print(f"SKIP ({e})")
                continue

            for lab, kind, params in perturbations:
                perturbed = _perturb_ts(data_4d, kind, rng, **params)
                try:
                    rois, mat = _reconstruct_matrix(
                        perturbed, entry["atlas_path"], atlas_level, run_tr, fit_kind)
                    ps = compute_primary_score_from_matrix(mat, rois)
                    cat = _classify_one(ps["mean_intra"], ps["mean_inter_hetero"],
                                        thresh_intra, thresh_delta)
                    results[lab].append(dict(
                        intra=ps["mean_intra"], delta=ps["delta"],
                        category=cat, matrix=mat, rois=rois))
                except Exception:
                    pass
            if verbose:
                print("done")

        all_bids_results[label] = results

        # Per-BIDS figure
        if output_dir:
            bids_dir_safe = re.sub(r'[^a-zA-Z0-9_]', '_', f"{sp}_{bd}")
            _plot_single_bids_perturbation(
                results, perturbations, label, thresh_intra, thresh_delta,
                compute_primary_score_from_matrix, _classify_one,
                opj(output_dir, bids_dir_safe, "perturbation.png"))

    # Combined comparison figure
    if output_dir and all_bids_results:
        _plot_combined_perturbation_comparison(
            all_bids_results, perturbations,
            opj(output_dir, "combined_perturbation_comparison.png"))

    return all_bids_results


def _plot_single_bids_perturbation(results, perturbations, title,
                                    thresh_intra, thresh_delta,
                                    score_fn, classify_fn, output_path):
    """Per-BIDS perturbation figure: matrices (top) + proportions (bottom)."""
    n_pert = len(perturbations)
    fig = plt.figure(figsize=(n_pert * 3.2, 7.5))
    gs = gridspec.GridSpec(2, n_pert + 1,
                           height_ratios=[1.6, 1.0],
                           width_ratios=[1.0] * n_pert + [0.04],
                           wspace=0.10, hspace=0.40)

    all_v = []
    for lab, _, _ in perturbations:
        for r in results.get(lab, []):
            all_v.extend(r["matrix"].ravel())
    all_v = np.array(all_v)
    all_v = all_v[np.isfinite(all_v)]
    vmax = max(float(np.nanpercentile(np.abs(all_v), 97)), 0.05) if len(all_v) else 0.5
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    for ci, (lab, _, _) in enumerate(perturbations):
        ax = fig.add_subplot(gs[0, ci])
        runs = results.get(lab, [])
        if runs:
            med = np.nanmedian(np.stack([r["matrix"] for r in runs], 0), 0)
            ax.imshow(med, cmap="RdBu_r", norm=norm, aspect="equal",
                      interpolation="nearest")
            ax.text(0.02, 0.98, f"n={len(runs)}", transform=ax.transAxes,
                    fontsize=6, va="top")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(lab, fontsize=8, fontweight="bold")

    cbar_ax = fig.add_subplot(gs[0, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm, cmap="RdBu_r"),
                 cax=cbar_ax, label="r")

    ax_bar = fig.add_subplot(gs[1, :n_pert])
    x = np.arange(n_pert)
    bottoms = np.zeros(n_pert)
    for cat in CATS:
        vals = []
        for lab, _, _ in perturbations:
            runs = results.get(lab, [])
            n = max(len(runs), 1)
            vals.append(100.0 * sum(1 for r in runs if r["category"] == cat) / n)
        vals = np.array(vals)
        ax_bar.bar(x, vals, 0.65, bottom=bottoms, color=_CAT_COLORS[cat],
                   label=cat, edgecolor="white", lw=0.5)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 5:
                ax_bar.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold")
        bottoms += vals
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([l for l, _, _ in perturbations], fontsize=7, rotation=20, ha="right")
    ax_bar.set_ylabel("% runs"); ax_bar.set_ylim(0, 105)
    ax_bar.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
                  ncol=4, fontsize=7, frameon=False)
    ax_bar.spines["top"].set_visible(False); ax_bar.spines["right"].set_visible(False)

    fig.suptitle(title, fontsize=11, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    _save(fig, output_path)


def _plot_combined_perturbation_comparison(all_bids_results, perturbations,
                                            output_path):
    """
    Combined bar plot: for each BIDS dataset, show % Specific for
    Original vs each perturbation type, side by side.
    """
    bids_labels = list(all_bids_results.keys())
    n_bids = len(bids_labels)
    pert_labels = [lab for lab, _, _ in perturbations]
    n_pert = len(pert_labels)

    # Build data: specificity rate per (BIDS, perturbation)
    spec_rates = np.zeros((n_bids, n_pert))
    for bi, bl in enumerate(bids_labels):
        res = all_bids_results[bl]
        for pi, pl in enumerate(pert_labels):
            runs = res.get(pl, [])
            n = max(len(runs), 1)
            spec_rates[bi, pi] = 100.0 * sum(1 for r in runs if r["category"] == "Specific") / n

    fig, ax = plt.subplots(figsize=(max(8, n_bids * 1.8), 5))
    bar_w = 0.8 / n_pert
    x = np.arange(n_bids)
    pert_colors = ["#333333", "#E69F00", "#56B4E9", "#CC79A7", "#009E73"]

    for pi, pl in enumerate(pert_labels):
        offset = (pi - n_pert / 2 + 0.5) * bar_w
        bars = ax.bar(x + offset, spec_rates[:, pi], bar_w,
                      label=pl, color=pert_colors[pi % len(pert_colors)],
                      edgecolor="white", lw=0.5)
        for xi, v in enumerate(spec_rates[:, pi]):
            if v > 3:
                ax.text(x[xi] + offset, v + 1, f"{v:.0f}",
                        ha="center", va="bottom", fontsize=5.5)

    ax.set_xticks(x)
    ax.set_xticklabels(bids_labels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("% Specific runs", fontsize=10)
    ax.set_ylim(0, 110)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25),
              ncol=n_pert, fontsize=7, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Specificity rate: original vs perturbations (per BIDS)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    _save(fig, output_path)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 2 — QC DENSITY PLOTS (tSNR, SNR)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_qc_densities(df_qc, metrics=None, output_dir=None):
    """
    Plot KDE density distributions for QC metrics, following the style of
    threshold_explorer Fig 01 middle panel:
      - faded per-BIDS KDE curves (colored by species)
      - bold global pooled density on top
      - P25 / P50 / P75 percentile markers

    Parameters
    ----------
    df_qc     : DataFrame with QC columns (func_TSNR_0, func_avg_snr_gray,
                anat_template_correlation, etc.) + species + bids_dir
    metrics   : list of (column_name, display_title, x_label) tuples,
                or None for defaults
    output_dir: directory for output PNGs
    """
    if metrics is None:
        metrics = [
            ("func_TSNR_0",             "Functional tSNR",            "tSNR"),
            ("func_avg_snr_gray",       "Functional SNR (gray)",      "SNR"),
            ("anat_template_correlation","Anat registration quality",  "Template correlation"),
            ("anat_cortical_contrast",  "Anat cortical contrast",     "CNR"),
        ]
        # keep only metrics actually present in df_qc
        metrics = [(c, t, xl) for c, t, xl in metrics if c in df_qc.columns]

    if not metrics:
        print("  [QC density] no QC metrics found in df_qc columns")
        print(f"  Available: {[c for c in df_qc.columns if 'snr' in c.lower() or 'tsnr' in c.lower() or 'qc' in c.lower()]}")
        return

    n_met = len(metrics)
    fig, axes = plt.subplots(1, n_met, figsize=(n_met * 5.5, 5))
    if n_met == 1:
        axes = [axes]

    sp_list = _sp_order(df_qc)
    bids_list = sorted(df_qc["bids_dir"].dropna().unique()) if "bids_dir" in df_qc.columns else []

    # map BIDS → dominant species
    bids_species = {}
    if "species" in df_qc.columns:
        for bd in bids_list:
            sp_counts = df_qc.loc[df_qc["bids_dir"] == bd, "species"].value_counts()
            bids_species[bd] = sp_counts.index[0] if len(sp_counts) else None

    for ax, (col, title, xlabel) in zip(axes, metrics):
        # Per-BIDS faded KDEs
        for bd in bids_list:
            v = df_qc.loc[df_qc["bids_dir"] == bd, col].dropna().values
            v = v[np.isfinite(v)]
            if len(v) < 5:
                continue
            xr = np.linspace(v.min(), v.max(), 300)
            try:
                k = _stats.gaussian_kde(v, bw_method="silverman")
                sp = bids_species.get(bd)
                c = _sp_color(sp) if sp else "#888888"
                ax.plot(xr, k(xr), color=c, lw=1.0, alpha=0.30, zorder=1)
                ax.fill_between(xr, k(xr), color=c, alpha=0.05, zorder=1)
            except Exception:
                pass

        # Global pooled density
        g = df_qc[col].dropna().values
        g = g[np.isfinite(g)]
        if len(g) < 5:
            ax.set_title(title + " (insufficient data)", fontsize=10)
            continue
        xg = np.linspace(g.min(), g.max(), 500)
        kg = _stats.gaussian_kde(g, bw_method="silverman")
        ax.fill_between(xg, kg(xg), color="#0072B2", alpha=0.12, zorder=2)
        ax.plot(xg, kg(xg), color="#0072B2", lw=2.2, zorder=3, label="global (pooled)")

        # Percentiles
        trans = ax.get_xaxis_transform()
        for p, ls in [(25, ":"), (50, "--"), (75, ":")]:
            xp = np.percentile(g, p)
            ax.axvline(xp, color="#555555", lw=1.0, ls=ls, zorder=4)
            ax.text(xp, 0.92, f"P{p}", rotation=90, fontsize=6,
                    color="#555555", ha="right", va="top", transform=trans)

        # Median label
        med = np.median(g)
        ax.text(med, 0.98, f" median={med:.2f}", color="#0072B2", fontsize=7,
                ha="left", va="top", fontweight="bold", transform=trans)

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Density", fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Legends
        h0, l0 = ax.get_legend_handles_labels()
        leg1 = ax.legend(h0, l0, fontsize=7, frameon=False, loc="upper right")
        ax.add_artist(leg1)
        if sp_list:
            sp_handles = [mpatches.Patch(color=_sp_color(sp), alpha=0.5, label=sp)
                          for sp in sp_list]
            ax.legend(handles=sp_handles, fontsize=6, frameon=False,
                      loc="upper left", ncol=1, title="species", title_fontsize=6)

    fig.suptitle("Cross-species QC metric distributions", fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = opj(output_dir, "qc_density_distributions.png") if output_dir else None
    _save(fig, out)
    if not out:
        plt.show()
    return fig


def plot_qc_densities_per_metric(df_qc, metrics=None, output_dir=None):
    """
    Same as plot_qc_densities but one figure per metric (larger, more readable).
    Each figure has two panels:
      A: KDE density (per-BIDS + global)
      B: Strip plot per species (individual subjects visible)
    """
    if metrics is None:
        metrics = [
            ("func_TSNR_0",             "Functional tSNR",            "tSNR"),
            ("func_avg_snr_gray",       "Functional SNR (gray)",      "SNR"),
            ("anat_template_correlation","Anat registration quality",  "Template corr."),
            ("anat_cortical_contrast",  "Anat cortical contrast",     "CNR"),
        ]
        metrics = [(c, t, xl) for c, t, xl in metrics if c in df_qc.columns]

    sp_list = _sp_order(df_qc)
    bids_list = sorted(df_qc["bids_dir"].dropna().unique()) if "bids_dir" in df_qc.columns else []
    bids_species = {}
    if "species" in df_qc.columns:
        for bd in bids_list:
            sc = df_qc.loc[df_qc["bids_dir"] == bd, "species"].value_counts()
            bids_species[bd] = sc.index[0] if len(sc) else None

    rng = np.random.default_rng(42)

    for col, title, xlabel in metrics:
        if col not in df_qc.columns:
            continue

        fig, (ax_kde, ax_strip) = plt.subplots(1, 2, figsize=(14, 5),
                                                gridspec_kw={"width_ratios": [1.2, 1]})

        # ── Panel A: KDE density ─────────────────────────────────────────────
        for bd in bids_list:
            v = df_qc.loc[df_qc["bids_dir"] == bd, col].dropna().values
            v = v[np.isfinite(v)]
            if len(v) < 5:
                continue
            xr = np.linspace(v.min(), v.max(), 300)
            try:
                k = _stats.gaussian_kde(v, bw_method="silverman")
                sp = bids_species.get(bd)
                c = _sp_color(sp) if sp else "#888"
                ax_kde.plot(xr, k(xr), color=c, lw=1.0, alpha=0.30)
                ax_kde.fill_between(xr, k(xr), color=c, alpha=0.05)
            except Exception:
                pass

        g = df_qc[col].dropna().values
        g = g[np.isfinite(g)]
        if len(g) >= 5:
            xg = np.linspace(g.min(), g.max(), 500)
            kg = _stats.gaussian_kde(g, bw_method="silverman")
            ax_kde.fill_between(xg, kg(xg), color="#0072B2", alpha=0.12)
            ax_kde.plot(xg, kg(xg), color="#0072B2", lw=2.2, label="global")
            trans = ax_kde.get_xaxis_transform()
            for p, ls in [(25, ":"), (50, "--"), (75, ":")]:
                xp = np.percentile(g, p)
                ax_kde.axvline(xp, color="#555", lw=1.0, ls=ls)
                ax_kde.text(xp, 0.92, f"P{p}", rotation=90, fontsize=6,
                            color="#555", ha="right", va="top", transform=trans)

        ax_kde.set_xlabel(xlabel); ax_kde.set_ylabel("Density")
        ax_kde.set_title(f"A: {title} — pooled distribution", fontsize=10, fontweight="bold")
        ax_kde.set_ylim(bottom=0)
        ax_kde.spines["top"].set_visible(False); ax_kde.spines["right"].set_visible(False)
        if sp_list:
            ax_kde.legend(
                handles=[mpatches.Patch(color=_sp_color(s), alpha=0.5, label=s) for s in sp_list],
                fontsize=6, frameon=False, loc="upper left", ncol=1)

        # ── Panel B: Strip plot per species ──────────────────────────────────
        for si, sp in enumerate(sp_list):
            v = df_qc.loc[df_qc["species"] == sp, col].dropna().values
            v = v[np.isfinite(v)]
            if len(v) == 0:
                continue
            jit = rng.uniform(-0.2, 0.2, len(v))
            ax_strip.scatter(v, si + jit, s=12, alpha=0.5, color=_sp_color(sp),
                             edgecolors="none")
            # median bar
            med = np.median(v)
            ax_strip.plot([med, med], [si - 0.3, si + 0.3],
                          color=_sp_color(sp), lw=2.5, zorder=5)

        ax_strip.set_yticks(range(len(sp_list)))
        ax_strip.set_yticklabels(sp_list, fontsize=9)
        ax_strip.set_xlabel(xlabel)
        ax_strip.set_title(f"B: {title} — per species", fontsize=10, fontweight="bold")
        ax_strip.spines["top"].set_visible(False); ax_strip.spines["right"].set_visible(False)

        plt.tight_layout()
        out = opj(output_dir, f"qc_{col}.png") if output_dir else None
        _save(fig, out)
        if not out:
            plt.show()
