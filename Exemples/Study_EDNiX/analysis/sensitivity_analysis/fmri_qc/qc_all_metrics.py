"""
Former Supp_Figure_4-6_QC.py (every QC metric, every dataset, exploratory plots).
The paper QC figures (Supp. 5-6) are drawn by library/qc_figures.py (run_paper.py).
OUT_DIR is unchanged: run_selection/config_ednix.QC_EXCEL_PATH reads
<OUT_DIR>/qc_raw_per_individual.xlsx written here.

EDNiX BIDS Tools v6  —  QC Figures only
=========================================
Generates QC figures using the same plotting functions as EDNiX_figures.py.
"""

import os
import sys
import warnings
import math
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from scipy import stats as _stats

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
    PALETTE, get_atlas_label_path, find_species_path,
    collect_multi_species, _phylo_sort,
)
from Exemples.Study_EDNiX.analysis.EDNiX_figures import (
    plot_cross_species_qc,
    plot_qc_fc_report,
)

opj = os.path.join

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

ATLAS_LIB = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
ATLAS_NAME = 'EDNIxCSC'
from Exemples.Study_EDNiX.analysis import paper_data as pdat
OUT_DIR = pdat.QC_DIR                    # created below before the first write

# ═══════════════════════════════════════════════════════════════════════════════
# SPECIES ORDERING & COLOURS (same as main script)
# ═══════════════════════════════════════════════════════════════════════════════

PHYLO_ORDER = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
               "Macaque", "Human", "Dog", "Cat", "Pig"]

# species colours of the paper
from Exemples.Study_EDNiX.analysis.paper_style import SPECIES_COLORS  # noqa: E402

# ═══════════════════════════════════════════════════════════════════════════════
# QC COLUMN DEFINITIONS — SCIENTIFICALLY LOGICAL ORDER
# ═══════════════════════════════════════════════════════════════════════════════

# Anatomical QC
ANAT_METRICS = [
    "anat_template_correlation",
    "anat_global_snr",
    "anat_cortical_contrast",
    "anat_avg_snr_gray",
    "anat_avg_snr_white",
    "anat_fwhm",
    "anat_noise_estimate",
    "anat_regional_snr_Cerebellum_Left",
    "anat_regional_snr_Cerebellum_Right",
    "anat_regional_snr_Cerebellum White_Left",
    "anat_regional_snr_Cerebellum White_Right",
    "anat_regional_snr_Isocortex_Left",
    "anat_regional_snr_Isocortex_Right",
    "anat_regional_snr_Allocortex_Left",
    "anat_regional_snr_Allocortex_Right",
    "anat_regional_snr_Periallocortex_Left",
    "anat_regional_snr_Periallocortex_Right",
    "anat_regional_snr_Subcortical areas_Left",
    "anat_regional_snr_Subcortical areas_Right",
    "anat_regional_snr_Diencephalon_Left",
    "anat_regional_snr_Diencephalon_Right",
    "anat_regional_snr_Cortical White matter_Left",
    "anat_regional_snr_Cortical White matter_Right",
    "anat_regional_snr_4th ventricles_Left",
    "anat_regional_snr_4th ventricles_Right",
    "anat_regional_snr_Brain stem_Left",
    "anat_regional_snr_Brain stem_Right",
    "anat_nmi",
]

# Functional QC
FUNC_METRICS = [
    "func_avg_snr_gray",
    "func_avg_snr_white",
    "func_cnr",
    "func_cortical_contrast",
    "func_TSNR",
    "func_TSNRcvarinv",
    "func_mean_fd",
    "func_mean_dvars",
    "func_gcor",
    "func_pct_anat_uncovered",
    "func_pct_fmri_uncovered",
    "func_mi",
    "func_cc",
    "func_jaccard",
    "func_nmi",
    "func_avg_velocity_0",
    "func_avg_velocity_1",
    "func_avg_velocity_2",
    "func_avg_velocity_3",
    "func_avg_velocity_4",
    "func_avg_velocity_5",
    "func_avg_enorm",
    "func_avg_outcount",
    "func_fwhm",
    "func_stdev",
]

# Network metrics
NET_METRICS = [
    "net_mean_correlation",
    "net_std_correlation",
    "net_top_eigenvalue",
    "net_eigenvalue_ratio",
    "net_p_intra_vs_inter",
    "net_p_left_vs_right",
    "net_silhouette",
    "net_davies_bouldin",

    "homotopic_mean_primary",
    "cross_mean_primary",

    "homotopic_mean_all",
    "cross_mean_all",

    "sp_specificity_index_primary",
    "sp_specificity_index_EDNiX",

    "sp_specificity_index_all",
    "contam_global_median_r",
    "contam_global_mean_abs_r",
    "contam_frac_positive",
    "contam_frac_high",
    "contam_iqr_r"
]

# ═══════════════════════════════════════════════════════════════════════════════
# SPECIES / BIDS DIRECTORIES
# ═══════════════════════════════════════════════════════════════════════════════

species_bids_dict = {
    'Rat':        '/scratch2/EDNiX/Rat/BIDS_Grandjean',
    'Mouse':      '/scratch2/EDNiX/Mouse/BIDS_Grandjean2',
    'Dog':        '/scratch2/EDNiX/Dog/BIDS_Boch_K9',
    'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin',
}
species_multi_bids = {
    'Macaque': [
        '/scratch2/EDNiX/Macaque/BIDS_Tremblay',
        '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
        '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',
    ],
    'Human': [
        '/scratch2/EDNiX/Human/BIDS_Castrillon',
        '/scratch2/EDNiX/Human/BIDS_Babayan',
    ],
    'Marmoset': [
        '/scratch2/EDNiX/Marmoset/BIDS_Tian',
        '/scratch2/EDNiX/Marmoset/BIDS_Hori',
    ],
}

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Build species_config and resolve atlas paths
# ═══════════════════════════════════════════════════════════════════════════════

species_config = {}
atlas_label_paths = {}
species_atlas_fragments = {}

for species, bids_dir in species_bids_dict.items():
    species_config[species] = {
        'bids_dirs': [bids_dir], 'list_to_keep': [], 'list_to_remove': [],
    }
for species, bids_dirs in species_multi_bids.items():
    species_config[species] = {
        'bids_dirs': bids_dirs, 'list_to_keep': [], 'list_to_remove': [],
    }

for species in species_config:
    try:
        fragment = find_species_path(ATLAS_LIB, species)
        species_atlas_fragments[species] = fragment
        atlas_label_paths[species] = get_atlas_label_path(
            ATLAS_LIB, fragment, ATLAS_NAME, prefer_statslut=False)
        print(f"  [{species}] {atlas_label_paths[species]}")
    except (FileNotFoundError, ValueError) as e:
        print(f"  [WARN] {species}: {e}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Collect QC data
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 64)
print("  Collecting QC data...")
print("=" * 64)

data = collect_multi_species(
    species_config,
    regions_of_interest=None,
    extract=('qc',),
    atlas_name=ATLAS_NAME,
    atlas_label_paths=atlas_label_paths,
    fit_kind='correlation',
    atlas_library_root=ATLAS_LIB,
    species_atlas_fragments=species_atlas_fragments,
)

df_qc = data.get('qc')

print(f"  QC data collected: {len(df_qc)} rows, {len(df_qc.columns)} columns")

# ── Save RAW per-individual QC (before any session averaging) ────────────────
raw_qc_path = opj(OUT_DIR, "qc_raw_per_individual.xlsx")
os.makedirs(OUT_DIR, exist_ok=True)
df_qc.to_excel(raw_qc_path, index=False)
print(f"  [saved] {raw_qc_path}")

# ── Plot granularity ──────────────────────────────────────────────────────────
# PREVIOUSLY averaged away 'session' here (_qc_id deliberately excluded it),
# collapsing every session of a subject into one point before plotting --
# correct for the morphometry figures this pattern was written for ("one dot
# per animal"), wrong for QC: a subject with one bad session among several
# good ones would show a clean subject-level average, hiding exactly the
# variability QC figures exist to surface. 'run' is real granularity too
# where extract_qc_paths finds more than one QC file per session (empty/
# constant here -- confirmed this pipeline's QC json is one-per-session, not
# one-per-run -- but included for when it isn't, elsewhere).
#
# QC_PLOT_UNIT: "session" (default) -- one point per (subject, session), no
# averaging at all. "subject" restores the old behavior for comparison.
QC_PLOT_UNIT = "session"

if QC_PLOT_UNIT == "subject" and 'session' in df_qc.columns:
    _qc_id = [c for c in ('species', 'bids_dir', 'subject') if c in df_qc.columns]
    _qc_skip = set(_qc_id) | {'session', 'run'}
    _qc_num = [c for c in df_qc.columns
               if c not in _qc_skip and pd.api.types.is_numeric_dtype(df_qc[c])]
    _qc_cat = [c for c in df_qc.columns
               if c not in _qc_skip and c not in _qc_num]

    df_qc_plot = df_qc.groupby(_qc_id, as_index=False)[_qc_num].mean()
    if _qc_cat:
        _cat_first = (df_qc.sort_values('session')
                      .groupby(_qc_id, as_index=False)[_qc_cat]
                      .first())
        df_qc_plot = df_qc_plot.merge(_cat_first, on=_qc_id, how='left')
    print(f"  [session_avg] qc: {len(df_qc)} -> {len(df_qc_plot)} rows "
          f"(QC_PLOT_UNIT='subject' -- one point per animal, averaged "
          f"across session/run)")
else:
    df_qc_plot = df_qc
    n_multi = 0
    if {'subject', 'session'} <= set(df_qc.columns):
        n_multi = int((df_qc.groupby(['species', 'bids_dir', 'subject'] if
                                     'bids_dir' in df_qc.columns else ['species', 'subject'])
                      ['session'].transform('nunique') > 1).sum())
    print(f"  [no_avg] qc: {len(df_qc)} row(s) plotted as-is, one point per "
          f"(subject, session{', run' if 'run' in df_qc.columns else ''}) "
          f"({n_multi} row(s) belong to a subject with >1 session)")


# ── Filter to only columns that exist ──────────────────────────────────────────
def _filter_existing_metrics(df, metric_list):
    """Return only metrics that exist in the DataFrame."""
    return [m for m in metric_list if m in df.columns]


def _split_metrics(metric_list):
    """Split a list of metrics into two roughly equal halves."""
    n = len(metric_list)
    half = (n + 1) // 2
    return metric_list[:half], metric_list[half:]


ANAT_METRICS_EXISTING = _filter_existing_metrics(df_qc_plot, ANAT_METRICS)
FUNC_METRICS_EXISTING = _filter_existing_metrics(df_qc_plot, FUNC_METRICS)
NET_METRICS_EXISTING = _filter_existing_metrics(df_qc_plot, NET_METRICS)

print(f"\n  Existing anatomical metrics: {len(ANAT_METRICS_EXISTING)}/{len(ANAT_METRICS)}")
print(f"  Existing functional metrics: {len(FUNC_METRICS_EXISTING)}/{len(FUNC_METRICS)}")
print(f"  Existing network metrics: {len(NET_METRICS_EXISTING)}/{len(NET_METRICS)}")

# ── Split metrics evenly between left and right columns ──────────────────────
ANAT_LEFT, ANAT_RIGHT = _split_metrics(ANAT_METRICS_EXISTING)
FUNC_LEFT, FUNC_RIGHT = _split_metrics(FUNC_METRICS_EXISTING)
NET_LEFT, NET_RIGHT = _split_metrics(NET_METRICS_EXISTING)

print(f"\n  Anatomical split: {len(ANAT_LEFT)} left, {len(ANAT_RIGHT)} right")
print(f"  Functional split: {len(FUNC_LEFT)} left, {len(FUNC_RIGHT)} right")
print(f"  Network split: {len(NET_LEFT)} left, {len(NET_RIGHT)} right")

# ── Build global BIDS colour map ──────────────────────────────────────────────
_bd_sources = []
if df_qc_plot is not None and not df_qc_plot.empty and 'bids_dir' in df_qc_plot.columns:
    _bd_sources.extend(df_qc_plot['bids_dir'].unique().tolist())

_all_bids_dirs = sorted(set(
    os.path.basename(str(bd).rstrip('/')) for bd in _bd_sources))
GLOBAL_BIDS_COL = {
    bd: PALETTE[i % len(PALETTE)] for i, bd in enumerate(_all_bids_dirs)
}
print(f"\n  Global BIDS colour map: {GLOBAL_BIDS_COL}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Generate QC figures
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 64)
print("  Generating QC figures...")
print("=" * 64)

os.makedirs(OUT_DIR, exist_ok=True)

out_paths = []

# ── Fig QC-1: Anatomical QC ─────────────────────────────────────────────────
print("\n  Generating Fig QC-1: Anatomical QC...")
if ANAT_METRICS_EXISTING:
    p = plot_cross_species_qc(
        df_qc_plot,
        output_path=opj(OUT_DIR, "QC_Fig1_anatomical.png"),
        func_metrics=ANAT_LEFT,
        anat_metrics=ANAT_RIGHT,
        bids_col=GLOBAL_BIDS_COL,
    )
    if p:
        out_paths.append(p)
else:
    print("  No anatomical metrics found - skipping.")

# ── Fig QC-2: Functional QC ─────────────────────────────────────────────────
print("\n  Generating Fig QC-2: Functional QC...")
if FUNC_METRICS_EXISTING:
    p = plot_cross_species_qc(
        df_qc_plot,
        output_path=opj(OUT_DIR, "QC_Fig2_functional.png"),
        func_metrics=FUNC_LEFT,
        anat_metrics=FUNC_RIGHT,
        bids_col=GLOBAL_BIDS_COL,
    )
    if p:
        out_paths.append(p)
else:
    print("  No functional metrics found - skipping.")

# ── Fig QC-3: Network metrics ──────────────────────────────────────────────
print("\n  Generating Fig QC-3: Network metrics...")
if NET_METRICS_EXISTING:
    p = plot_cross_species_qc(
        df_qc_plot,
        output_path=opj(OUT_DIR, "QC_Fig3_network.png"),
        func_metrics=NET_LEFT,
        anat_metrics=NET_RIGHT,
        bids_col=GLOBAL_BIDS_COL,
    )
    if p:
        out_paths.append(p)
else:
    print("  No network metrics found - skipping.")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Summary Report
# ═══════════════════════════════════════════════════════════════════════════════

def generate_qc_summary(df_qc, output_path):
    """Generate a summary table of QC metrics per species."""
    if df_qc is None or df_qc.empty:
        return

    # Select all numeric columns
    summary_cols = []
    for col in df_qc.columns:
        if pd.api.types.is_numeric_dtype(df_qc[col]):
            summary_cols.append(col)

    if not summary_cols:
        return

    # Group by species and compute mean ± std
    summary = df_qc.groupby('species')[summary_cols].agg(['mean', 'std', 'count']).round(3)

    # Save to CSV
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    summary.to_csv(output_path)
    print(f"  QC summary → {output_path}")
    return summary


summary_path = opj(OUT_DIR, 'qc_summary.csv')
generate_qc_summary(df_qc_plot, summary_path)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Correlation Analysis: sp_specificity_index vs Functional Metrics
# ═══════════════════════════════════════════════════════════════════════════════
# USING METHOD 1: Partial Correlation controlling for BIDS effects
# ORDERED BY R² (coefficient of determination)
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 64)
print("  Correlation Analysis: sp_specificity_index vs Functional Metrics")
print("  Method: Partial Correlation (controlling for BIDS effects)")
print("  Ordering: By R² (highest first)")
print("=" * 64)


def find_best_partial_correlations(df, target_col="sp_specificity_index",
                                   metric_list=None, n_top=3):
    """
    Find the best partial correlations between sp_specificity_index and metrics.
    Controls for BIDS effects using residualization.
    Orders by R² (coefficient of determination).

    Returns:
    --------
    DataFrame with metrics sorted by R² (highest first)
    """
    if target_col not in df.columns:
        print(f"  Target column '{target_col}' not found in data!")
        return None

    # Filter metrics that exist in the DataFrame
    metrics = [m for m in metric_list if m in df.columns]
    print(f"  Testing {len(metrics)} functional metrics...")

    results = []
    for col in metrics:
        if col == target_col:
            continue

        # Remove NaN values
        valid = df[[target_col, col, 'bids_dir']].dropna()
        if len(valid) < 5:
            print(f"    {col}: insufficient data (n={len(valid)})")
            continue

        # Check if there's variance
        if valid[col].nunique() <= 1:
            print(f"    {col}: constant value, skipping")
            continue

        # ─── Partial correlation controlling for BIDS ──────────────────────
        # Create dummy variables for BIDS
        bids_dummies = pd.get_dummies(valid['bids_dir'], prefix='bids', drop_first=True)

        # Prepare design matrix
        X = bids_dummies.values
        X = np.column_stack([np.ones(len(X)), X])

        x = valid[col].values
        y = valid[target_col].values

        try:
            # Residual of x (metric) controlling for BIDS
            beta_x = np.linalg.lstsq(X, x, rcond=None)[0]
            resid_x = x - X @ beta_x

            # Residual of y (target) controlling for BIDS
            beta_y = np.linalg.lstsq(X, y, rcond=None)[0]
            resid_y = y - X @ beta_y

            # Correlation of residuals (partial correlation). pearsonr() on residuals
            # uses n-2 df; a partial correlation controlling for k covariates has
            # n-2-k df (k = number of BIDS dummies), so the p-value is recomputed.
            r_partial = float(np.corrcoef(resid_x, resid_y)[0, 1])
            k_cov = X.shape[1] - 1
            df_partial = len(valid) - 2 - k_cov
            if df_partial < 1 or abs(r_partial) >= 1:
                p_partial = np.nan
            else:
                t_partial = r_partial * np.sqrt(df_partial / (1.0 - r_partial ** 2))
                p_partial = float(2 * _stats.t.sf(abs(t_partial), df_partial))
            r2_partial = r_partial ** 2

            # Also compute regular Spearman for comparison
            r_spearman, p_spearman = _stats.spearmanr(x, y)
            r2_spearman = r_spearman ** 2

            results.append({
                'metric': col,
                'n': len(valid),
                'partial_r': r_partial,
                'partial_r2': r2_partial,
                'partial_p': p_partial,
                'spearman_r': r_spearman,
                'spearman_r2': r2_spearman,
                'spearman_p': p_spearman,
                'df': df_partial,
                'significant': bool(p_partial < 0.05) if not np.isnan(p_partial) else False
            })

            print(f"    {col}: partial r={r_partial:.3f}, R²={r2_partial:.3f}, p={p_partial:.4f}")

        except Exception as e:
            print(f"    {col}: partial correlation failed: {e}")
            continue

    if not results:
        print("  No correlations found!")
        return None

    results_df = pd.DataFrame(results)

    # Sort by R² (highest first)
    results_df = results_df.sort_values('partial_r2', ascending=False)

    # Keep only top N
    results_df = results_df.head(n_top)

    return results_df


def plot_partial_correlation_dashboard(df, target_col="sp_specificity_index_EDNiX",
                                       metric_list=None, output_path=None,
                                       n_top=3, bids_col=None):
    """
    Create a dashboard with two rows:
    Row 1: Top 3 partial correlations across species (ordered by R²)
    Row 2: Top 3 partial correlations across BIDS (ordered by R²)
    """
    # Find partial correlations across all data
    corr_results = find_best_partial_correlations(df, target_col, metric_list, n_top)

    if corr_results is None or corr_results.empty:
        print("  No correlations found!")
        return None

    # Get top metrics
    top_metrics = corr_results['metric'].tolist()

    print(f"\n  Top {len(top_metrics)} correlations (ordered by R²):")
    print(corr_results[['metric', 'partial_r', 'partial_r2', 'partial_p', 'n']].to_string(index=False))

    # ── Calculate per-BIDS partial correlations for top metrics ──────────────
    bids_partial_data = {}
    metric_avg_r2 = {}

    if 'bids_dir' in df.columns:
        bids_dirs = sorted(df['bids_dir'].dropna().unique())

        for metric in top_metrics:
            bids_partial_data[metric] = {}
            r2_values = []
            r_values = []
            p_values = []

            for bd in bids_dirs:
                bd_df = df[df['bids_dir'] == bd].dropna(subset=[target_col, metric])

                # Need enough data within BIDS to compute partial correlation
                if len(bd_df) >= 5:
                    try:
                        # For within-BIDS, we can't control for BIDS (only one BIDS)
                        # So we use regular correlation within BIDS
                        x_sp = bd_df[metric].values
                        y_sp = bd_df[target_col].values
                        # Both must have variance, otherwise spearmanr returns nan (and warns)
                        if np.ptp(x_sp) == 0 or np.ptp(y_sp) == 0:
                            continue
                        r_sp, p_sp = _stats.spearmanr(x_sp, y_sp)
                        if np.isnan(r_sp):
                            continue
                        r2 = r_sp ** 2

                        bids_partial_data[metric][bd] = {
                            'r': r_sp,
                            'r2': r2,
                            'p': p_sp,
                            'n': len(bd_df),
                            'significant': p_sp < 0.05
                        }
                        r2_values.append(r2)
                        r_values.append(r_sp)
                        p_values.append(p_sp)
                    except:
                        pass

            # Average R² across BIDS (this is the effect size metric)
            if r2_values:
                # Fixed-effect meta-analysis of the within-dataset Spearman rho
                # (Fisher z, weights (n-3)/1.06). Replaces the mean of r^2 (biased
                # upward in small datasets: E[r^2] = 1/(n-1) under H0) and the
                # geometric mean of p-values (not a valid combined test).
                ns = np.array([bids_partial_data[metric][bd]['n'] for bd in bids_partial_data[metric]])
                rs = np.clip(np.array(r_values), -0.999999, 0.999999)
                w = (ns - 3) / 1.06
                keep = w > 0
                if keep.any():
                    z_bar = np.sum(w[keep] * np.arctanh(rs[keep])) / np.sum(w[keep])
                    se = 1.0 / np.sqrt(np.sum(w[keep]))
                    meta_r = float(np.tanh(z_bar))
                    meta_p = float(2 * _stats.norm.sf(abs(z_bar) / se))
                else:
                    meta_r, meta_p = np.nan, np.nan
                metric_avg_r2[metric] = {
                    'avg_r2': meta_r ** 2 if not np.isnan(meta_r) else 0.0,
                    'avg_r': meta_r,
                    'avg_p': meta_p,
                    'n_bids': len(r_values)
                }

    # Sort metrics by average R² across BIDS
    sorted_bids_metrics = sorted(metric_avg_r2.items(),
                                key=lambda x: x[1]['avg_r2'],
                                reverse=True)
    top_bids_metrics = [m[0] for m in sorted_bids_metrics[:n_top]]

    print(f"\n  Top {len(top_bids_metrics)} across BIDS (by average R²):")
    for metric, stats in sorted_bids_metrics[:n_top]:
        print(f"    {metric}: avg R²={stats['avg_r2']:.3f}, avg ρ={stats['avg_r']:.3f}")

    # ── Create figure with 2 rows ──────────────────────────────────────────────
    n_cols = min(n_top, 3)
    n_rows = 2

    fig = plt.figure(figsize=(15, 10))

    # ── Row 1: Top correlations across all data (partial, ordered by R²) ──
    for i, metric in enumerate(top_metrics):
        ax = plt.subplot2grid((n_rows, n_cols), (0, i), rowspan=1)

        row = corr_results[corr_results['metric'] == metric].iloc[0]
        r_partial = row['partial_r']
        r2_partial = row['partial_r2']
        p_partial = row['partial_p']
        n = row['n']

        # Plot data by species
        valid = df[[target_col, metric, 'species']].dropna()
        x = valid[metric].values
        y = valid[target_col].values

        species_order = _phylo_sort(valid['species'].unique())
        for sp in species_order:
            sp_df = valid[valid['species'] == sp]
            color = SPECIES_COLORS.get(sp, "#888")
            ax.scatter(sp_df[metric], sp_df[target_col],
                       color=color, s=40, alpha=0.8, label=sp, edgecolors='k', linewidths=0.5)

        # Overall regression line (for visualization)
        x = bd_df[metric].values
        y = bd_df[target_col].values
        # Need at least 3 distinct x values and some y variance for a regression line
        if len(x) > 2 and np.ptp(x) > 0 and np.ptp(y) > 0:
            slope, intercept, r_val, p_val, se = _stats.linregress(x, y)
            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = intercept + slope * x_line
            ax.plot(x_line, y_line, '-', color=color, lw=1.5, alpha=0.6)

        # Title with partial correlation stats (R and R²)
        signif = '***' if p_partial < 0.001 else '**' if p_partial < 0.01 else '*' if p_partial < 0.05 else ''
        ax.set_title(f"{metric}\nPartial ρ={r_partial:.3f}  R²={r2_partial:.3f}  p={p_partial:.4f} {signif}\nn={n}",
                     fontsize=9, fontweight='bold')
        ax.set_xlabel(metric, fontsize=8)
        if i == 0:
            ax.set_ylabel(target_col, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add species legend (only once)
        if i == 0:
            sp_handles = [mpatches.Patch(facecolor=SPECIES_COLORS.get(sp, "#888"),
                                         label=sp, edgecolor='k', linewidth=0.3)
                          for sp in species_order]
            ax.legend(handles=sp_handles, fontsize=7, frameon=False,
                      bbox_to_anchor=(1.01, 1), loc='upper left', title='Species')

    # ── Row 2: Top correlations across BIDS (ordered by R²) ──────────────────
    if bids_col is None:
        bids_colors = {bd: PALETTE[i % len(PALETTE)] for i, bd in enumerate(sorted(df['bids_dir'].dropna().unique()))}
    else:
        bids_colors = bids_col

    all_bids = sorted(df['bids_dir'].dropna().unique())

    for i, metric in enumerate(top_bids_metrics[:n_cols]):
        ax = plt.subplot2grid((n_rows, n_cols), (1, i), rowspan=1)

        # Get average stats for this metric
        stats = metric_avg_r2.get(metric, {})
        avg_r2 = stats.get('avg_r2', 0)
        avg_r = stats.get('avg_r', 0)

        # Calculate overall partial correlation for this metric
        valid_all = df[[target_col, metric]].dropna()
        r_all, p_all = _stats.spearmanr(valid_all[metric].values, valid_all[target_col].values)
        r2_all = r_all ** 2

        # Plot data by BIDS
        valid = df[[target_col, metric, 'bids_dir']].dropna()

        for bd in all_bids:
            bd_df = valid[valid['bids_dir'] == bd]
            if len(bd_df) >= 3:
                color = bids_colors.get(bd, "#888")

                # Scatter points
                ax.scatter(bd_df[metric], bd_df[target_col],
                           color=color, s=35, alpha=0.7, edgecolors='k', linewidths=0.4)

                # Per-BIDS regression line
                x = bd_df[metric].values
                y = bd_df[target_col].values
                if len(x) > 2 and np.ptp(x) > 0 and np.ptp(y) > 0:
                    slope, intercept, r_val, p_val, se = _stats.linregress(x, y)
                    x_line = np.linspace(x.min(), x.max(), 100)
                    y_line = intercept + slope * x_line
                    ax.plot(x_line, y_line, '-', color=color, lw=1.5, alpha=0.6)

        # Title with R² (the selection criterion)
        signif = '***' if p_all < 0.001 else '**' if p_all < 0.01 else '*' if p_all < 0.05 else ''
        ax.set_title(f"{metric}\nAll: ρ={r_all:.3f}  R²={r2_all:.3f}  p={p_all:.4f}\nAvg R²={avg_r2:.3f}",
                     fontsize=9, fontweight='bold')

        # Add per-BIDS correlation info
        if metric in bids_partial_data and bids_partial_data[metric]:
            corr_text = []
            for bd in all_bids:
                if bd in bids_partial_data[metric]:
                    c = bids_partial_data[metric][bd]
                    sig = '***' if c['p'] < 0.001 else '**' if c['p'] < 0.01 else '*' if c['p'] < 0.05 else ''
                    corr_text.append(f"{bd}: ρ={c['r']:.3f}{sig}")

            if corr_text:
                ax.text(0.02, 0.98, '\n'.join(corr_text), transform=ax.transAxes,
                        va='top', fontsize=6, bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

        ax.set_xlabel(metric, fontsize=8)
        if i == 0:
            ax.set_ylabel(target_col, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add BIDS legend (only once)
        if i == 0:
            bd_handles = [mpatches.Patch(facecolor=bids_colors.get(bd, "#888"),
                                         label=bd, edgecolor='k', linewidth=0.3)
                          for bd in all_bids]
            ax.legend(handles=bd_handles, fontsize=6, frameon=False,
                      bbox_to_anchor=(1.01, 1), loc='upper left', title='BIDS dir')

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close(fig)
        print(f"\n  [plot] {output_path}")
    else:
        plt.show()

    return corr_results, bids_partial_data, metric_avg_r2


# ── Run correlation analysis using Method 1 ───────────────────────────────────
if df_qc_plot is not None and not df_qc_plot.empty:
    # Use functional metrics only
    corr_results, bids_partial_data, metric_avg_r2 = plot_partial_correlation_dashboard(
        df_qc_plot,
        target_col="sp_specificity_index_EDNiX",
        metric_list=FUNC_METRICS_EXISTING,
        output_path=opj(OUT_DIR, "sp_specificity_correlation_dashboard.png"),
        n_top=3,
        bids_col=GLOBAL_BIDS_COL,
    )

    if corr_results is not None:
        # Save to CSV
        corr_path = opj(OUT_DIR, 'sp_specificity_correlations.csv')
        corr_results.to_csv(corr_path, index=False)
        print(f"\n  Correlation results saved to: {corr_path}")

        # Save per-BIDS correlations
        if bids_partial_data:
            bids_corr_list = []
            for metric, bd_data in bids_partial_data.items():
                for bd, stats in bd_data.items():
                    bids_corr_list.append({
                        'metric': metric,
                        'bids_dir': bd,
                        'r_spearman': stats['r'],
                        'r2_spearman': stats['r2'],
                        'p_spearman': stats['p'],
                        'n': stats['n'],
                        'significant': stats['significant']
                    })

            if bids_corr_list:
                bids_corr_df = pd.DataFrame(bids_corr_list)
                bids_corr_path = opj(OUT_DIR, 'sp_specificity_correlations_per_bids.csv')
                bids_corr_df.to_csv(bids_corr_path, index=False)
                print(f"  Per-BIDS correlations saved to: {bids_corr_path}")

        # Save metric ranking by R²
        if metric_avg_r2:
            ranking_df = pd.DataFrame([
                {'metric': m, 'avg_r2': s['avg_r2'], 'avg_r': s['avg_r'], 'avg_p': s['avg_p'], 'n_bids': s['n_bids']}
                for m, s in metric_avg_r2.items()
            ]).sort_values('avg_r2', ascending=False)
            ranking_path = opj(OUT_DIR, 'sp_specificity_metric_ranking_by_r2.csv')
            ranking_df.to_csv(ranking_path, index=False)
            print(f"  Metric ranking by R² saved to: {ranking_path}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Report
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 64)
print("  QC Figures complete!")
print("=" * 64)
print(f"  Output directory: {OUT_DIR}")
print(f"  QC summary: {summary_path}")
print("\n  Generated files:")
for p in out_paths:
    print(f"    {os.path.basename(p)}")
print(f"    sp_specificity_correlations.csv")
print(f"    sp_specificity_correlation_dashboard.png")
print(f"    sp_specificity_metric_ranking_by_r2.csv")
print("\n  Done.")