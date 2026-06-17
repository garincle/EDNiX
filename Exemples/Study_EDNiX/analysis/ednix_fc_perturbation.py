#!/usr/bin/env python3
"""
ednix_fc_perturbation.py — Time-series level sensitivity analysis
=================================================================

Unlike Fig 08 in ednix_threshold_explorer.py (which perturbs the finished
correlation MATRIX), this script applies perturbations directly to the
preprocessed fMRI RESIDUAL TIME-SERIES, then reconstructs the correlation
matrices from scratch using the same Nilearn pipeline as EDNiX step 10.

This is a more rigorous test because it preserves the statistical structure
of the time-series (autocorrelation, spectral content) before perturbation,
and reconstructs the matrix through the same masking/erosion/extraction
steps as the real pipeline.

PERTURBATION TYPES
------------------
  original : no perturbation (reference)
  noise    : add Gaussian noise to each voxel's time-series (σ = fraction of
             voxel-wise temporal std)
  smooth   : temporally smooth the time-series (moving average window)
  shift    : add a spatially uniform signal shift (global signal injection)
  phase    : randomize the phase of each voxel's Fourier spectrum while
             preserving amplitude (destroys temporal correlations)

USAGE
-----
  from ednix_fc_perturbation import run_perturbation_analysis
  run_perturbation_analysis(
      df_scored,                           # from threshold_explorer
      cats,                                # classification labels
      thresh_intra=0.24, thresh_delta=0.1,
      species="Human",
      bids_root_template="/scratch2/EDNiX/{species}/{bids_dir}",
      output_dir="figures/perturbation",
  )
"""

import os
import warnings
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as _mpl_cm
from matplotlib.colors import TwoSlopeNorm
from scipy import stats as _stats

# Import scoring/classification from the threshold explorer
try:
    from ednix_threshold_explorer import (
        compute_primary_score_from_matrix, _classify_one,
        _load_matrix, _save_or_show, _CAT_COLORS, CATS,
        _SPECIES_COLORS, _sp_order
    )
except ImportError:
    from Exemples.Study_EDNiX.analysis.ednix_threshold_explorer import (
        compute_primary_score_from_matrix, _classify_one,
        _load_matrix, _save_or_show, _CAT_COLORS, CATS,
    )

opj = os.path.join
ope = os.path.exists


# =============================================================================
# TIME-SERIES PERTURBATIONS
# =============================================================================

def perturb_timeseries(ts_data, kind, rng, **params):
    """
    Apply a perturbation to a 4D fMRI array (x, y, z, t).

    Parameters
    ----------
    ts_data : np.ndarray, shape (x, y, z, t)
        The preprocessed fMRI residual volumes.
    kind : str
        'noise', 'smooth', 'shift', 'phase', or None (original).
    rng : np.random.Generator
    params : perturbation-specific parameters

    Returns
    -------
    perturbed : np.ndarray, same shape as ts_data
    """
    if kind is None:
        return ts_data.copy()

    data = ts_data.copy().astype(np.float64)
    mask = np.any(data != 0, axis=-1)  # brain mask from non-zero voxels
    n_t = data.shape[-1]

    if kind == "noise":
        sigma_frac = params.get("sigma_frac", 0.3)
        # Add Gaussian noise proportional to each voxel's temporal std
        for idx in zip(*np.where(mask)):
            vox_std = np.std(data[idx])
            if vox_std > 0:
                data[idx] += rng.normal(0, sigma_frac * vox_std, n_t)

    elif kind == "smooth":
        window = params.get("window", 5)
        # Temporal moving average (reduces high-frequency signal, inflates correlations)
        kernel = np.ones(window) / window
        for idx in zip(*np.where(mask)):
            data[idx] = np.convolve(data[idx], kernel, mode='same')

    elif kind == "shift":
        shift_std = params.get("shift_std", 0.5)
        # Add a spatially uniform signal (simulates global signal contamination)
        global_signal = rng.normal(0, shift_std, n_t)
        data[mask] += global_signal[np.newaxis, :]

    elif kind == "phase":
        # Phase randomization: preserves power spectrum, destroys correlations
        for idx in zip(*np.where(mask)):
            ts = data[idx]
            ft = np.fft.rfft(ts)
            phases = rng.uniform(0, 2 * np.pi, len(ft))
            phases[0] = 0  # keep DC component
            if n_t % 2 == 0:
                phases[-1] = 0  # keep Nyquist
            ft_rand = np.abs(ft) * np.exp(1j * phases)
            data[idx] = np.fft.irfft(ft_rand, n=n_t)

    return data


# =============================================================================
# MATRIX RECONSTRUCTION (mirrors EDNiX step 10)
# =============================================================================

def reconstruct_matrix(ts_4d, atlas_path, atlas_level, tr,
                       fit_kinds=("correlation",)):
    """
    Reconstruct FC matrices from perturbed 4D data using the same
    Nilearn pipeline as EDNiX step 10.

    Parameters
    ----------
    ts_4d : np.ndarray (x, y, z, t)
    atlas_path : str, path to the atlas NIfTI (4D, level in dim4)
    atlas_level : int, which level to extract from the 4D atlas
    tr : float, repetition time in seconds
    fit_kinds : tuple of str, e.g. ("correlation", "partial correlation")

    Returns
    -------
    dict : {fit_kind: (rois_list, matrix_2d)}
    """
    from nilearn.input_data import NiftiLabelsMasker
    from nilearn.connectome import ConnectivityMeasure

    # Load atlas and extract the requested level
    atlas_img = nib.load(atlas_path)
    atlas_data = atlas_img.get_fdata()

    if atlas_data.ndim == 4:
        atlas_level_data = atlas_data[:, :, :, atlas_level]
    else:
        atlas_level_data = atlas_data

    # Create a 3D atlas NIfTI in memory
    atlas_3d = nib.Nifti1Image(atlas_level_data, atlas_img.affine, atlas_img.header)

    # Create a temporary 4D NIfTI from the perturbed data
    # Use the atlas header for spatial reference
    func_img = nib.Nifti1Image(ts_4d, atlas_img.affine)

    # Extract ROI labels
    labels = sorted(set(int(v) for v in np.unique(atlas_level_data) if v != 0))

    masker = NiftiLabelsMasker(
        labels_img=atlas_3d,
        detrend=False,
        smoothing_fwhm=None,
        low_pass=None,
        high_pass=None,
        t_r=tr,
        standardize='zscore_sample',
        memory=None, verbose=0
    )

    try:
        time_series = masker.fit_transform(func_img)
    except Exception as e:
        warnings.warn(f"Masker failed: {e}")
        return {}

    # Get the actual labels that were extracted
    extracted_labels = masker.labels_
    roi_names = [str(l) for l in extracted_labels]

    results = {}
    for fit_kind in fit_kinds:
        conn = ConnectivityMeasure(kind=fit_kind, standardize='zscore_sample')
        matrix = conn.fit_transform([time_series])[0]
        np.fill_diagonal(matrix, 0)
        results[fit_kind] = (roi_names, matrix)

    return results


# =============================================================================
# FIND PATHS FOR A SPECIES
# =============================================================================

def _find_func_paths(df_scored, species, bids_root_template):
    """
    For each run in df_scored matching species, find:
      - the residual fMRI path
      - the atlas path
      - the existing correlation matrix path (for ROI name reference)
    """
    import re, glob

    entries = []
    sp_df = df_scored[df_scored["species"] == species] if "species" in df_scored.columns else df_scored

    for _, row in sp_df.iterrows():
        sp = str(row.get("species", ""))
        bd = str(row.get("bids_dir", ""))
        sub = str(row.get("subject", ""))
        ses = str(row.get("session", ""))
        corr_path = str(row.get("corr_matrix_path", ""))

        root = bids_root_template.format(species=sp, bids_dir=bd)
        sub_dir = opj(root, f"sub-{sub}")
        ses_dir = opj(sub_dir, f"ses-{ses}") if ses else sub_dir

        # Find the residual fMRI
        res_pattern = opj(ses_dir, "func", "acpc-func", "postprocessed_rs",
                          "*_desc-fMRI_residual.nii.gz")
        res_files = sorted(glob.glob(res_pattern))

        # Find the atlas in labels dir
        atlas_pattern = opj(ses_dir, "func", "acpc-func", "labels",
                            f"{sub}_seg-EDNIxCSC*_dseg.nii.gz")
        atlas_files = sorted(glob.glob(atlas_pattern))

        if res_files and atlas_files:
            entries.append({
                "species": sp, "bids_dir": bd, "subject": sub, "session": ses,
                "residual_path": res_files[0],
                "atlas_path": atlas_files[0],
                "corr_matrix_path": corr_path,
                "category": str(row.get("fp_category", "")),
            })

    return entries


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_perturbation_analysis(
    df_scored, cats, thresh_intra, thresh_delta,
    species="Human",
    bids_root_template="/scratch2/EDNiX/{species}/{bids_dir}",
    atlas_level=2,
    fit_kind="correlation",
    tr=None,          # auto-detected from NIfTI header if None
    max_runs=10,      # cap for speed
    seed=42,
    output_dir=None,
    verbose=True,
):
    """
    Full perturbation analysis for one species:
      1. Load residual fMRI 4D + atlas for each run
      2. For each perturbation type, perturb the time-series
      3. Reconstruct the FC matrix via Nilearn
      4. Score and classify the reconstructed matrix
      5. Plot comparison (original vs perturbed)
    """
    df = df_scored.copy()
    df["fp_category"] = cats.values if hasattr(cats, "values") else cats

    entries = _find_func_paths(df, species, bids_root_template)
    if not entries:
        print(f"  [perturbation] no runs found for {species}")
        return None

    if len(entries) > max_runs:
        print(f"  [perturbation] capping to {max_runs}/{len(entries)} runs")
        entries = entries[:max_runs]

    print(f"  [perturbation] {species}: {len(entries)} runs")

    perturbations = [
        ("Original",              None,     {}),
        ("Noise (σ=0.3×std)",     "noise",  dict(sigma_frac=0.3)),
        ("Temporal smooth (w=5)", "smooth", dict(window=5)),
        ("Global shift (σ=0.5)",  "shift",  dict(shift_std=0.5)),
        ("Phase randomize",       "phase",  {}),
    ]

    rng = np.random.default_rng(seed)

    # Results: for each perturbation, collect per-run scores
    results = {label: [] for label, _, _ in perturbations}

    for entry in entries:
        if verbose:
            print(f"    {entry['subject']}/{entry['session']}...", end=" ", flush=True)

        # Load 4D residual
        try:
            img = nib.load(entry["residual_path"])
            data_4d = img.get_fdata()
            run_tr = tr or float(img.header.get_zooms()[-1])
        except Exception as e:
            if verbose:
                print(f"SKIP ({e})")
            continue

        for label, kind, params in perturbations:
            perturbed = perturb_timeseries(data_4d, kind, rng, **params)

            # Reconstruct matrix
            try:
                mat_dict = reconstruct_matrix(
                    perturbed, entry["atlas_path"], atlas_level, run_tr,
                    fit_kinds=(fit_kind,))
                if fit_kind not in mat_dict:
                    continue
                rois, mat = mat_dict[fit_kind]
            except Exception as e:
                if verbose:
                    print(f"[{label}: matrix failed]", end=" ")
                continue

            # Score
            ps = compute_primary_score_from_matrix(mat, rois)
            cat = _classify_one(ps["mean_intra"], ps["mean_inter_hetero"],
                                thresh_intra, thresh_delta)
            results[label].append({
                "subject": entry["subject"],
                "session": entry["session"],
                "intra": ps["mean_intra"],
                "inter_hetero": ps["mean_inter_hetero"],
                "delta": ps["delta"],
                "category": cat,
                "matrix": mat,
                "rois": rois,
            })

        if verbose:
            print("done")

    # ── PLOT ─────────────────────────────────────────────────────────────────
    n_pert = len(perturbations)
    fig = plt.figure(figsize=(n_pert * 3.4, 8))
    gs = gridspec.GridSpec(2, n_pert + 1,
                           height_ratios=[1.6, 1.0],
                           width_ratios=[1.0] * n_pert + [0.04],
                           wspace=0.10, hspace=0.40)

    # Compute unified vmax
    all_mats = []
    for label, _, _ in perturbations:
        for r in results[label]:
            all_mats.append(r["matrix"])
    if all_mats:
        allv = np.concatenate([m.ravel() for m in all_mats])
        allv = allv[np.isfinite(allv)]
        vmax = max(float(np.nanpercentile(np.abs(allv), 97)), 0.05)
    else:
        vmax = 0.5
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    for c_i, (label, _, _) in enumerate(perturbations):
        runs = results[label]
        ax = fig.add_subplot(gs[0, c_i])

        if runs:
            med_mat = np.nanmedian(np.stack([r["matrix"] for r in runs], 0), 0)
            med_ps = compute_primary_score_from_matrix(med_mat, runs[0]["rois"])
            med_cat = _classify_one(med_ps["mean_intra"], med_ps["mean_inter_hetero"],
                                    thresh_intra, thresh_delta)
            ax.imshow(med_mat, cmap="RdBu_r", norm=norm,
                      aspect="equal", interpolation="nearest")
            ax.text(0.02, 0.98, f"n={len(runs)}", transform=ax.transAxes,
                    fontsize=6, va="top", ha="left")
            ax.text(0.5, -0.06,
                    f"intra={med_ps['mean_intra']:.3f}  Δ={med_ps['delta']:.3f}\n"
                    f"→ {med_cat}",
                    transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                    color=_CAT_COLORS.get(med_cat, "#444"), fontweight="bold")
        else:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    fontsize=9, transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(label, fontsize=9, fontweight="bold")

    cbar_ax = fig.add_subplot(gs[0, -1])
    fig.colorbar(_mpl_cm.ScalarMappable(norm=norm, cmap="RdBu_r"),
                 cax=cbar_ax, label="r")

    # Bottom: stacked proportions
    ax_bar = fig.add_subplot(gs[1, :n_pert])
    x = np.arange(n_pert)
    bottoms = np.zeros(n_pert)
    for cat in CATS:
        vals = []
        for label, _, _ in perturbations:
            runs = results[label]
            n_total = len(runs)
            n_cat = sum(1 for r in runs if r["category"] == cat)
            vals.append(100.0 * n_cat / max(n_total, 1))
        vals = np.array(vals)
        ax_bar.bar(x, vals, 0.65, bottom=bottoms, color=_CAT_COLORS[cat],
                   label=cat, edgecolor="white", lw=0.5)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 5:
                ax_bar.text(xi, b + v / 2, f"{v:.0f}%", ha="center", va="center",
                            fontsize=8, color="white", fontweight="bold")
        bottoms += vals
    pert_labels = [f"{label}\n(n={len(results[label])})"
                   for label, _, _ in perturbations]
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(pert_labels, fontsize=8)
    ax_bar.set_ylabel(f"% of {species} runs", fontsize=9)
    ax_bar.set_ylim(0, 105)
    ax_bar.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18),
                  ncol=len(CATS), fontsize=8, frameon=False)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    fig.suptitle(f"{species}: time-series perturbation sensitivity\n"
                 f"(thresholds intra>{thresh_intra:.3f}, Δ>{thresh_delta:.3f})",
                 fontsize=12, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.92])

    out_path = opj(output_dir, f"perturbation_{species}.png") if output_dir else None
    _save_or_show(fig, out_path)

    # Summary table
    summary = []
    for label, _, _ in perturbations:
        runs = results[label]
        n = len(runs)
        props = {c: sum(1 for r in runs if r["category"] == c) for c in CATS}
        med_intra = float(np.nanmedian([r["intra"] for r in runs])) if runs else np.nan
        med_delta = float(np.nanmedian([r["delta"] for r in runs])) if runs else np.nan
        summary.append({"perturbation": label, "n_runs": n,
                        "median_intra": med_intra, "median_delta": med_delta,
                        **{f"n_{c}": props[c] for c in CATS}})
    df_summary = pd.DataFrame(summary)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        df_summary.to_csv(opj(output_dir, f"perturbation_{species}_summary.csv"),
                          index=False)
    if verbose:
        print(f"\n  Perturbation summary ({species}):")
        print(df_summary.to_string(index=False))

    return dict(results=results, summary=df_summary, fig=fig)
