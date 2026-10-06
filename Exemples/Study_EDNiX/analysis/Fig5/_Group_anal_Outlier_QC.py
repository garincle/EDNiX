"""
_Group_anal_Outlier_QC.py
──────────────────────────
Per-subject outlier QC against a group-level effect (3dttest++ or 3dLMEr).

For a given group result (the "overall effect" map produced by
`_3dttest_EDNiX` or `_3dLMEr_EDNiX`), this module:

  1. Gathers the individual subject maps that went *into* that group model.
  2. Computes, for every subject, a distance-to-group-effect score
     (spatial correlation distance by default).
  3. Flags outliers with a robust (MAD-based) statistic — no assumption
     of normality, resistant to the very outliers you're trying to find.
  4. Plots a per-subject QC figure (distance vs. subject, outliers in red).
  5. Optionally renders a subject-vs-group *difference* map, always with
     a fixed colour scale (vmin=-1, vmax=+1) — never a data-driven
     vmax mirrored into vmin. Both the subject map and the group map are
     first rescaled onto the same [-1, 1] range so the difference is
     meaningful and comparable across seeds/species.

Design notes
------------
* The group/overall effect is the *reference* the distances are measured
  against — it is deliberately never plotted as one of the points in the
  per-subject QC figure. It only appears as a background/threshold
  concept, never as a 2-D scatter entry.
* This module does not re-run AFNI. It reads the outputs already produced
  by `_3dttest_EDNiX` (Grp_SBA_3dTTEST/<seed>/) or `_3dLMEr_EDNiX`
  (Grp_SBA_3dLME_network/<base_name>/), so it can be run as a lightweight
  QC pass after either group analysis.
"""

import os
import glob
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from nilearn import plotting

opj = os.path.join
opb = os.path.basename
ope = os.path.exists


# =============================================================================
# 1) Core distance computation
# =============================================================================

def _squeeze_trailing_singleton(data):
    """Drop trailing singleton axes, e.g. (48, 82, 48, 1) -> (48, 82, 48).
    AFNI-written volumes routinely carry a stray 4th dim of size 1 (a
    sub-brick axis with a single element) — that's not a real grid
    mismatch, just a shape bookkeeping artifact, and should never trigger
    the "resample onto the group grid" error below."""
    while data.ndim > 3 and data.shape[-1] == 1:
        data = data[..., 0]
    return data


def _masked_vector(img_path, mask_data):
    """Load a NIfTI and return its data restricted to voxels where mask_data != 0."""
    data = nib.load(img_path).get_fdata()
    data = _squeeze_trailing_singleton(data)
    mask_data = _squeeze_trailing_singleton(mask_data)
    if data.shape != mask_data.shape:
        raise ValueError(
            f"Grid mismatch for {img_path}: {data.shape} vs mask {mask_data.shape}. "
            "Resample onto the group grid before computing distances."
        )
    return data[mask_data != 0]


def _distance(subj_vec, group_vec, metric="correlation"):
    """Return (distance, pearson_r) between one subject's vector and the group vector."""
    finite = np.isfinite(subj_vec) & np.isfinite(group_vec)
    subj_vec, group_vec = subj_vec[finite], group_vec[finite]

    if subj_vec.size < 2 or np.std(subj_vec) == 0 or np.std(group_vec) == 0:
        return np.nan, np.nan

    r = np.corrcoef(subj_vec, group_vec)[0, 1]

    if metric == "correlation":
        dist = 1.0 - r                      # 0 = identical spatial pattern, 2 = anti-correlated
    elif metric == "cosine":
        cos = np.dot(subj_vec, group_vec) / (np.linalg.norm(subj_vec) * np.linalg.norm(group_vec))
        dist = 1.0 - cos
    elif metric == "euclidean":
        # z-score each map within the mask first so amplitude differences across
        # subjects/species don't dominate a purely spatial-pattern comparison.
        s = (subj_vec - subj_vec.mean()) / subj_vec.std()
        g = (group_vec - group_vec.mean()) / group_vec.std()
        dist = np.sqrt(np.mean((s - g) ** 2))
    else:
        raise ValueError(f"Unknown metric: {metric}")

    return dist, r


def compute_group_distance(group_map_path, subject_maps, mask_path, metric="correlation"):
    """
    Compute each subject's distance to the group/overall effect map.

    Parameters
    ----------
    group_map_path : str
        Path to the overall-effect map (e.g. *_ttest-stat_fisher_cmap.nii.gz
        or a 3dLMEr *_coef.nii.gz contrast map). This map is the reference
        only — it is never itself treated as an observation.
    subject_maps : dict {subject_id: path} or list[(subject_id, path)]
        The individual maps that fed into the group model.
    mask_path : str
        Mask restricting the comparison (the same group_overlap mask used
        for the group analysis is strongly recommended).
    metric : 'correlation' | 'cosine' | 'euclidean'

    Returns
    -------
    pd.DataFrame with columns: subject, distance, pearson_r, n_voxels
    """
    if isinstance(subject_maps, dict):
        subject_maps = list(subject_maps.items())

    mask_data  = nib.load(mask_path).get_fdata()
    group_vec  = _masked_vector(group_map_path, mask_data)

    rows = []
    for subj_id, subj_path in subject_maps:
        if not ope(subj_path):
            print(f"  [WARN] missing map for subject={subj_id}: {subj_path}")
            continue
        try:
            subj_vec = _masked_vector(subj_path, mask_data)
        except ValueError as e:
            print(f"  [WARN] {e}")
            continue
        dist, r = _distance(subj_vec, group_vec, metric=metric)
        rows.append({"subject": str(subj_id), "distance": dist, "pearson_r": r,
                      "n_voxels": int(mask_data.sum())})

    df = pd.DataFrame(rows).sort_values("distance", ascending=False).reset_index(drop=True)
    return df


# =============================================================================
# 2) Robust outlier flagging (median absolute deviation)
# =============================================================================

def flag_outliers(df, column="distance", threshold=3.5, direction="two-sided"):
    """
    Flag outliers using the modified z-score (Iglewicz & Hoaglin, 1993):
        Mi = 0.6745 * (xi - median) / MAD
    Robust to the very outliers it is meant to detect, unlike mean/SD.
    threshold=3.5 is the conventional cutoff.

    direction:
      - 'two-sided' (default, original behavior): flag |Mi| > threshold —
        unusually high OR unusually low values.
      - 'upper': flag only Mi > threshold — unusually high values. Use this
        when only "more distant / less correlated" is actually bad and an
        unusually LOW value (very well correlated) should never be
        penalized, e.g. distance-to-group-effect QC.
      - 'lower': flag only Mi < -threshold — the mirror image of 'upper'.
    """
    df = df.copy()
    x = df[column].values.astype(float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))

    if mad == 0:
        # Degenerate case (all values identical or MAD collapses) — fall back to SD.
        sd = np.nanstd(x)
        z = np.zeros_like(x) if sd == 0 else (x - np.nanmean(x)) / sd
    else:
        z = 0.6745 * (x - med) / mad

    df["modified_z"] = z
    if direction == "upper":
        df["is_outlier"] = z > threshold
    elif direction == "lower":
        df["is_outlier"] = z < -threshold
    elif direction == "two-sided":
        df["is_outlier"] = np.abs(z) > threshold
    else:
        raise ValueError(f"Unknown direction: {direction!r} (expected 'two-sided', 'upper', or 'lower')")
    return df


# =============================================================================
# 3) Per-subject QC plot (distance vs subject — the group effect is NOT a point)
# =============================================================================

def plot_distance_qc(df, output_path, title, threshold=3.5, column="distance",
                      higher_is_worse=True, dpi=300):
    """
    Bar plot of each subject's distance to the group effect, sorted worst→best.
    Outliers (flagged by `flag_outliers`) are highlighted in red with a
    threshold line. The group/overall effect itself is never plotted as a
    bar — it is the implicit reference (distance = 0) the bars are measured
    from, so it never appears as a 2-D data point.
    """
    df_sorted = df.sort_values(column, ascending=not higher_is_worse).reset_index(drop=True)
    colors = np.where(df_sorted["is_outlier"], "#d62728", "#4472c4")

    fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(df_sorted)), 5))
    ax.bar(df_sorted["subject"], df_sorted[column], color=colors, edgecolor="black", linewidth=0.5)

    # Convert the modified-z threshold back into distance units for a reference line.
    med = np.nanmedian(df[column])
    mad = np.nanmedian(np.abs(df[column] - med))
    if mad > 0:
        cutoff = med + (threshold / 0.6745) * mad
        ax.axhline(cutoff, color="#d62728", linestyle="--", linewidth=1,
                   label=f"outlier cutoff (MAD, z>{threshold})")
        ax.legend(loc="upper right", fontsize=8)

    ax.set_ylabel(column.replace("_", " "))
    ax.set_xlabel("subject")
    ax.set_title(title)
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    print(f"  [saved] {output_path}")
    return output_path


# =============================================================================
# 4) Fixed-scale [-1, 1] normalization + subject-vs-group difference map
# =============================================================================

def _normalize_to_unit_range(data, target_percentile=95):
    """Same robust-scale normalization used elsewhere in the pipeline:
    the target_percentile-th percentile of |data| -> ±1, then clipped."""
    nz = np.isfinite(data) & (data != 0)
    if not np.any(nz):
        return None
    scale = np.percentile(np.abs(data[nz]), target_percentile)
    if scale <= 0:
        return None
    out = np.zeros_like(data)
    out[nz] = np.clip(data[nz] / scale, -1, 1)
    return out


def plot_fixed_scale_map(map_path, output_folder, template, title, cut_coords, dpi=300):
    """
    Plot a NIfTI map with a HARD-CODED vmin=-1 / vmax=+1 colour scale.
    Unlike `_safe_plot` elsewhere in the pipeline (which derives vmax from
    the data and lets nilearn mirror it into vmin=-vmax), this never lets
    the data set the scale — vmin is always -1 and vmax is always +1, so
    maps stay directly comparable across subjects/seeds/species.
    """
    if not ope(map_path):
        print(f"  [WARN] {map_path} not found, skipping plot")
        return False
    try:
        data = nib.load(map_path).get_fdata()
        if np.all(data == 0):
            print(f"  [WARN] {map_path} contains only zeros, skipping plot")
            return False
        display = plotting.plot_stat_map(
            map_path, bg_img=template, display_mode="y", cut_coords=cut_coords,
            title=title, cmap="cold_hot", colorbar=True,
            vmin=-1.0, vmax=1.0,
        )
        out_png = opj(output_folder, f"{opb(map_path).replace('.nii.gz', '')}.png")
        display.savefig(out_png, dpi=dpi)
        display.close()
        print(f"  [saved] {out_png}")
        return True
    except Exception as e:
        print(f"  [WARN] Could not save plot for {map_path}: {e}")
        return False


def compute_and_plot_difference_map(subject_map_path, group_map_path, mask_path,
                                     output_folder, subject_id, template=None,
                                     cut_coords=None, target_percentile=95, dpi=300):
    """
    Normalize both the subject map and the group map onto [-1, 1] (robust,
    percentile-based scaling — see `_normalize_to_unit_range`), then save
    and plot their voxelwise difference (subject - group) with a fixed
    vmin=-1 / vmax=+1 scale. Useful to *see where* an outlier subject
    departs from the overall effect, not just that it does.
    """
    mask_data  = nib.load(mask_path).get_fdata()
    group_img  = nib.load(group_map_path)
    subj_img   = nib.load(subject_map_path)

    group_norm = _normalize_to_unit_range(group_img.get_fdata(), target_percentile)
    subj_norm  = _normalize_to_unit_range(subj_img.get_fdata(),  target_percentile)
    if group_norm is None or subj_norm is None:
        print(f"  [WARN] normalization failed for subject={subject_id}")
        return None

    diff = np.zeros_like(group_norm)
    diff[mask_data != 0] = subj_norm[mask_data != 0] - group_norm[mask_data != 0]
    diff = np.clip(diff, -1, 1)

    diff_path = opj(output_folder, f"{subject_id}_vs_group_diff.nii.gz")
    nib.save(nib.Nifti1Image(diff, group_img.affine, group_img.header), diff_path)

    if template is not None and cut_coords is not None:
        plot_fixed_scale_map(
            diff_path, output_folder, template,
            f"{subject_id} minus group effect", cut_coords, dpi=dpi,
        )
    return diff_path


# =============================================================================
# 5) Orchestration — one seed / one contrast at a time
# =============================================================================

def run_outlier_qc(
    group_map_path,
    subject_maps,
    mask_path,
    output_folder,
    label,                      # e.g. seed_name or contrast name, used in filenames/titles
    metric="correlation",
    outlier_threshold=3.5,
    plot_diff_for_outliers=True,
    template=None,
    cut_coords=None,
    target_percentile=95,
):
    """
    Full QC pass for one group result: compute distances, flag outliers,
    save a CSV, save the QC bar plot, and (optionally) save fixed-scale
    difference maps for every flagged outlier subject.

    Returns the results DataFrame.
    """
    os.makedirs(output_folder, exist_ok=True)

    df = compute_group_distance(group_map_path, subject_maps, mask_path, metric=metric)
    if df.empty:
        print(f"  [SKIP] {label}: no valid subject/group comparisons")
        return df

    df = flag_outliers(df, column="distance", threshold=outlier_threshold)

    csv_path = opj(output_folder, f"{label}_outlier_QC.csv")
    df.to_csv(csv_path, index=False)
    print(f"  [saved] {csv_path}")

    plot_distance_qc(
        df, opj(output_folder, f"{label}_outlier_QC.png"),
        title=f"{label}: subject distance to overall effect ({metric})",
        threshold=outlier_threshold,
    )

    if plot_diff_for_outliers:
        subj_lookup = dict(subject_maps.items() if isinstance(subject_maps, dict) else subject_maps)
        for _, row in df[df["is_outlier"]].iterrows():
            subj_path = subj_lookup.get(row["subject"])
            if subj_path is None or not ope(subj_path):
                continue
            compute_and_plot_difference_map(
                subj_path, group_map_path, mask_path, output_folder, row["subject"],
                template=template, cut_coords=cut_coords, target_percentile=target_percentile,
            )

    return df


# =============================================================================
# 6) Convenience gatherers matching the existing pipeline's output layout
# =============================================================================

def gather_3dttest_subject_maps(seed_output_folder):
    """
    Collect {subject_id: path} for the per-subject average Fisher maps saved
    by `_3dttest_EDNiX` inside Grp_SBA_3dTTEST/<seed_name>/.
    """
    subject_maps = {}
    for f in glob.glob(opj(seed_output_folder, "*_avg_fisher_map.nii.gz")):
        subj_id = opb(f).replace("_avg_fisher_map.nii.gz", "")
        subject_maps[subj_id] = f
    return subject_maps


def run_outlier_qc_3dttest(seed_output_folder, seed_name, tag="", mask_path=None,
                            metric="correlation", outlier_threshold=3.5,
                            template=None, cut_coords=None):
    """
    QC wrapper for a `_3dttest_EDNiX` seed result.
    group map = {seed_name}{tag}_ttest-stat_fisher_cmap.nii.gz
    subjects  = *_avg_fisher_map.nii.gz in the same folder
    """
    group_map = opj(seed_output_folder, f"{seed_name}{tag}_ttest-stat_fisher_cmap.nii.gz")
    subject_maps = gather_3dttest_subject_maps(seed_output_folder)

    if mask_path is None:
        # Fall back to "wherever the group map is non-zero" if no explicit
        # group_overlap mask is passed in.
        mask_path = opj(seed_output_folder, f"{seed_name}{tag}_auto_mask.nii.gz")
        gimg = nib.load(group_map)
        auto_mask = (gimg.get_fdata() != 0).astype(np.int16)
        nib.save(nib.Nifti1Image(auto_mask, gimg.affine, gimg.header), mask_path)

    return run_outlier_qc(
        group_map, subject_maps, mask_path, seed_output_folder, label=f"{seed_name}{tag}",
        metric=metric, outlier_threshold=outlier_threshold,
        template=template, cut_coords=cut_coords,
    )


def gather_3dLMEr_subject_maps(design_matrix_path, sing_afni=None):
    """
    Collect one map per subject from a 3dLMEr design_matrix.txt, averaging
    across that subject's sessions/runs/hemispheres (left+mirrored-right)
    in memory. Returns {subject_id: path_to_averaged_nifti}, writing the
    averaged files next to the design matrix.
    """
    dm_df = pd.read_csv(design_matrix_path, sep="\t")
    out_dir = opj(os.path.dirname(design_matrix_path), "subject_avg_for_QC")
    os.makedirs(out_dir, exist_ok=True)

    subject_maps = {}
    for subj_id, rows in dm_df.groupby("Subj"):
        paths = [p for p in rows["InputFile"] if ope(p)]
        if not paths:
            continue
        ref_img = nib.load(paths[0])
        stack = []
        for p in paths:
            img = nib.load(p)
            if img.shape != ref_img.shape:
                print(f"  [WARN] grid mismatch, skipping {p} for subject {subj_id}")
                continue
            stack.append(img.get_fdata())
        if not stack:
            continue
        avg_data = np.mean(np.stack(stack, axis=0), axis=0)
        avg_path = opj(out_dir, f"{subj_id}_avg_for_QC.nii.gz")
        nib.save(nib.Nifti1Image(avg_data, ref_img.affine, ref_img.header), avg_path)
        subject_maps[str(subj_id)] = avg_path

    return subject_maps


def run_outlier_qc_3dLMEr(contrast_output_folder, design_matrix_path, contrast_name, tag="",
                           coef_or_zstat="coef", mask_path=None, metric="correlation",
                           outlier_threshold=3.5, template=None, cut_coords=None):
    """
    QC wrapper for a `_3dLMEr_EDNiX` contrast result.
    group map = {contrast_name}{tag}_coef.nii.gz  (or *_zstat.nii.gz)
    subjects  = per-subject average of the InputFile column in design_matrix.txt
    """
    group_map = opj(contrast_output_folder, f"{contrast_name}{tag}_{coef_or_zstat}.nii.gz")
    subject_maps = gather_3dLMEr_subject_maps(design_matrix_path)

    if mask_path is None:
        mask_path = opj(contrast_output_folder, f"{contrast_name}{tag}_auto_mask.nii.gz")
        gimg = nib.load(group_map)
        auto_mask = (gimg.get_fdata() != 0).astype(np.int16)
        nib.save(nib.Nifti1Image(auto_mask, gimg.affine, gimg.header), mask_path)

    return run_outlier_qc(
        group_map, subject_maps, mask_path, contrast_output_folder,
        label=f"{contrast_name}{tag}_{coef_or_zstat}",
        metric=metric, outlier_threshold=outlier_threshold,
        template=template, cut_coords=cut_coords,
    )


# =============================================================================
# 7) Cross-seed / cross-contrast summary — spot subjects that are
#    systematically off across many seeds, not just one
# =============================================================================

def summarize_across_labels(results_by_label, output_folder, dpi=300):
    """
    results_by_label : dict {label: df} — the per-label output of run_outlier_qc.
    Builds a subject x label distance matrix, saves it as CSV, and plots a
    heatmap so subjects that are outliers across many seeds/contrasts stand
    out from subjects that are only borderline on one.
    """
    os.makedirs(output_folder, exist_ok=True)
    wide = pd.concat(
        [df.set_index("subject")["distance"].rename(label) for label, df in results_by_label.items()
         if not df.empty],
        axis=1,
    )
    wide["n_outlier_flags"] = pd.concat(
        [df.set_index("subject")["is_outlier"].rename(label) for label, df in results_by_label.items()
         if not df.empty],
        axis=1,
    ).sum(axis=1)
    wide = wide.sort_values("n_outlier_flags", ascending=False)

    csv_path = opj(output_folder, "outlier_QC_summary.csv")
    wide.to_csv(csv_path)
    print(f"  [saved] {csv_path}")

    heat_cols = [c for c in wide.columns if c != "n_outlier_flags"]
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(heat_cols)), max(4, 0.3 * len(wide))))
    im = ax.imshow(wide[heat_cols].values, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(heat_cols)))
    ax.set_xticklabels(heat_cols, rotation=90, fontsize=7)
    ax.set_yticks(range(len(wide)))
    ax.set_yticklabels(
        [f"{s}  ({int(n)} flags)" for s, n in zip(wide.index, wide["n_outlier_flags"])],
        fontsize=7,
    )
    fig.colorbar(im, ax=ax, label="distance to group effect")
    ax.set_title("Subject distance to overall effect, across seeds/contrasts")
    plt.tight_layout()
    png_path = opj(output_folder, "outlier_QC_summary.png")
    fig.savefig(png_path, dpi=dpi)
    plt.close(fig)
    print(f"  [saved] {png_path}")

    return wide
