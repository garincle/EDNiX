import os
import warnings
import json
import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats
import seaborn as sns
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from fMRI.extract_filename import extract_filename

# =============================================================================
# CONSTANTS
# =============================================================================

_SENSORY_PATTERNS = ['Somatosensory', 'Motor', 'Visual', 'Auditory']

_CATEGORY_COLORS = {
    'Specific':   '#2ca02c',
    'Unspecific': '#ff7f0e',
    'No':         '#1f77b4',
    'Spurious':   '#d62728',
}

_CATEGORIES = ['Specific', 'Unspecific', 'No', 'Spurious']

# =============================================================================
# MATRIX VALIDATION & LOADING
# =============================================================================

def validate_matrix(matrix, roi_names):
    """
    Ensure matrix is square and consistent with roi_names.
    Attempts trimming if shape is off by one axis.
    """
    if matrix.shape[0] != matrix.shape[1]:
        if isinstance(roi_names, (list, np.ndarray)):
            if len(roi_names) == matrix.shape[0]:
                matrix = matrix[:, :len(roi_names)]
            elif len(roi_names) == matrix.shape[1]:
                matrix = matrix[:len(roi_names), :]
            else:
                raise ValueError(
                    f"ROI names count ({len(roi_names)}) doesn't match "
                    f"matrix dimensions ({matrix.shape[0]}x{matrix.shape[1]})"
                )

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Failed to make matrix square. Final shape: {matrix.shape}")

    if isinstance(roi_names, (list, np.ndarray)) and len(roi_names) != matrix.shape[0]:
        raise ValueError(
            f"ROI names count ({len(roi_names)}) doesn't match "
            f"matrix dimension ({matrix.shape[0]}) after adjustment"
        )

    return matrix, roi_names


def filter_bilateral_rois(matrix, roi_names):
    """
    Keep only ROIs that have both L_ and R_ counterparts.
    """
    left_rois  = {name[2:]: i for i, name in enumerate(roi_names) if name.startswith("L_")}
    right_rois = {name[2:]: i for i, name in enumerate(roi_names) if name.startswith("R_")}
    bilateral  = set(left_rois) & set(right_rois)

    keep = [i for i, name in enumerate(roi_names)
            if name.startswith(("L_", "R_")) and name[2:] in bilateral]

    return matrix[np.ix_(keep, keep)], [roi_names[i] for i in keep]


def load_and_validate_matrix(matrix_path):
    """
    Load a correlation matrix CSV, filter to bilateral ROIs,
    and return (matrix, roi_names).
    """
    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"Matrix file not found: {matrix_path}")

    df = pd.read_csv(matrix_path, index_col=0)
    df.columns = df.columns.astype(str).str.strip()
    df.index   = df.index.astype(str).str.strip()

    # Filter to bilateral pairs
    all_rois = set(df.columns)
    seen, bilateral_rois = set(), []
    for roi in df.columns:
        if roi not in seen:
            counterpart = ('R_' + roi[2:] if roi.startswith('L_')
                           else 'L_' + roi[2:] if roi.startswith('R_')
                           else None)
            if counterpart and counterpart in all_rois:
                for r in (roi, counterpart):
                    if r not in seen:
                        bilateral_rois.append(r)
                        seen.add(r)

    if bilateral_rois:
        df = df.loc[bilateral_rois, bilateral_rois]

    if list(df.columns) != list(df.index):
        raise ValueError("Header and index mismatch after bilateral filtering")

    return df.values.astype(float), list(df.columns)

# =============================================================================
# SPECIFICITY SCORING
# =============================================================================

def _find_homotopic_pairs(roi_names):
    """
    Return list of (pattern, l_idx, r_idx) for each primary sensory network
    found in roi_names.  Falls back to positional pairing if no L_/R_ prefix.
    """
    pairs = []
    for pattern in _SENSORY_PATTERNS:
        l = [i for i, n in enumerate(roi_names) if n.startswith('L_') and pattern in n]
        r = [i for i, n in enumerate(roi_names) if n.startswith('R_') and pattern in n]
        if l and r:
            pairs.append((pattern, l[0], r[0]))

    if not pairs:                          # fallback: no hemisphere prefix
        for pattern in _SENSORY_PATTERNS:
            m = [i for i, n in enumerate(roi_names) if pattern in n]
            if len(m) >= 2:
                pairs.append((pattern, m[0], m[1]))

    return pairs


def analyze_specificity(correlation_matrix, roi_names,
                        intra_thresh=0.240,
                        delta_thresh=0.050):
    """
    Score and classify one FC matrix.

    Scoring (all aggregations use median for robustness)
    -------
    intra        : median homotopic same-network r   (L_Motor ↔ R_Motor, …)
    inter        : median ALL cross-network pairs    (reference only)
    inter_hetero : median contralateral cross-network pairs only
                   L_netA ↔ R_netB  and  R_netA ↔ L_netB
                   excludes ipsilateral (L_netA ↔ L_netB, R_netA ↔ R_netB)
    delta        : intra − inter_hetero

    Classification
    --------------
    Spurious   : inter_hetero > intra          (delta < 0)
    No         : intra <= intra_thresh
    Unspecific : intra > intra_thresh  AND  delta <= delta_thresh
    Specific   : intra > intra_thresh  AND  delta >  delta_thresh
    """
    try:
        corr_matrix, roi_names = validate_matrix(correlation_matrix, roi_names)
    except Exception as e:
        return {'error': str(e)}

    base = {
        'target_specificity': None,
        'matrix_shape':       corr_matrix.shape,
        'n_rois':             len(roi_names),
        'intra_threshold':    intra_thresh,
        'delta_threshold':    delta_thresh,
    }

    homotopic_pairs = _find_homotopic_pairs(roi_names)

    if not homotopic_pairs:
        base['target_specificity'] = {
            'Specific_ROI_pair':        'Not found',
            'Intra_Correlation':        np.nan,
            'Inter_Correlation':        np.nan,
            'Inter_Correlation_Hetero': np.nan,
            'Delta':                    np.nan,
            'Intra_Threshold':          intra_thresh,
            'Delta_Threshold':          delta_thresh,
            'Category':                 'No',
            'n_homotopic_pairs':        0,
            'error':                    'No primary sensory ROI pairs found',
        }
        return base

    # ── intra ─────────────────────────────────────────────────────────────────
    intra_corrs  = [corr_matrix[l, r] for _, l, r in homotopic_pairs]
    median_intra = float(np.nanmedian(intra_corrs))

    # ── inter ─────────────────────────────────────────────────────────────────
    inter_all    = []
    inter_hetero = []

    n = len(homotopic_pairs)
    for i in range(n):
        for j in range(i + 1, n):
            _, ia_L, ia_R = homotopic_pairs[i]
            _, ib_L, ib_R = homotopic_pairs[j]

            for ia, ib in [(ia_L, ib_L), (ia_L, ib_R), (ia_R, ib_L), (ia_R, ib_R)]:
                r = float(corr_matrix[ia, ib])
                if np.isfinite(r):
                    inter_all.append(r)

            # contralateral cross-network only
            for ia, ib in [(ia_L, ib_R), (ia_R, ib_L)]:
                r = float(corr_matrix[ia, ib])
                if np.isfinite(r):
                    inter_hetero.append(r)

    median_inter        = float(np.nanmedian(inter_all))    if inter_all    else np.nan
    median_inter_hetero = float(np.nanmedian(inter_hetero)) if inter_hetero else np.nan
    delta               = (float(median_intra - median_inter_hetero)
                           if np.isfinite(median_inter_hetero) else np.nan)

    # ── classify ──────────────────────────────────────────────────────────────
    if not np.isfinite(median_intra) or not np.isfinite(median_inter_hetero):
        category = 'No'
    elif median_inter_hetero > median_intra:
        category = 'Spurious'
    elif median_intra <= intra_thresh:
        category = 'No'
    elif delta <= delta_thresh:
        category = 'Unspecific'
    elif delta > delta_thresh and median_intra > intra_thresh:
        category = 'Specific'
    else:
        warnings.warn(
            f"Unhandled classification state: intra={median_intra:.3f}, "
            f"inter_hetero={median_inter_hetero:.3f}, delta={delta:.3f} "
            f"— defaulting to 'Unclassified'"
        )
        category = 'Unclassified'

    pair_string = " + ".join(
        f"{roi_names[l][2:]}↔{roi_names[r][2:]}" for _, l, r in homotopic_pairs[:3]
    )
    if len(homotopic_pairs) > 3:
        pair_string += f" + {len(homotopic_pairs) - 3} more"

    base['target_specificity'] = {
        'Specific_ROI_pair':        pair_string,
        'Intra_Correlation':        round(median_intra,        3),
        'Inter_Correlation':        round(median_inter,        3),
        'Inter_Correlation_Hetero': round(median_inter_hetero, 3),
        'Delta':                    round(delta,               3),
        'Intra_Threshold':          intra_thresh,
        'Delta_Threshold':          delta_thresh,
        'Category':                 category,
        'n_homotopic_pairs':        len(homotopic_pairs),
    }
    return base

# =============================================================================
# NETWORK METRICS
# =============================================================================

def calculate_network_metrics(correlation_matrix):
    """
    Eigenvalue, clustering, and global correlation statistics.
    Mean/Median/Std are all reported explicitly; no aggregation choice needed.
    """
    metrics  = {}
    n_nodes  = correlation_matrix.shape[0]

    # Eigenvalue analysis
    try:
        eigs = np.sort(np.abs(np.linalg.eig(correlation_matrix)[0]))[::-1]
        metrics['Top_Eigenvalue']   = round(float(eigs[0]), 3)
        metrics['Eigenvalue_Ratio'] = round(float(eigs[0] / eigs[1]), 3) if len(eigs) > 1 else np.nan
    except Exception:
        metrics.update({'Top_Eigenvalue': np.nan, 'Eigenvalue_Ratio': np.nan})

    # Clustering metrics
    if n_nodes > 1:
        try:
            labels = SpectralClustering(
                n_clusters=min(7, n_nodes - 1),
                affinity='precomputed',
                assign_labels='kmeans',
                random_state=42,
            ).fit_predict((correlation_matrix + 1) / 2)
            metrics.update({
                'Silhouette_Score':   round(silhouette_score(correlation_matrix, labels),   3),
                'Davies_Bouldin':     round(davies_bouldin_score(correlation_matrix, labels), 3),
                'Calinski_Harabasz':  round(calinski_harabasz_score(correlation_matrix, labels), 3),
            })
        except Exception:
            metrics.update({
                'Silhouette_Score':  np.nan,
                'Davies_Bouldin':    np.nan,
                'Calinski_Harabasz': np.nan,
            })

    # Global correlation statistics (lower triangle)
    mask = np.tri(n_nodes, k=-1).astype(bool)
    if mask.any():
        vals = correlation_matrix[mask]
        metrics.update({
            'Mean_Correlation':   round(float(np.nanmean(vals)),   3),
            'Median_Correlation': round(float(np.nanmedian(vals)), 3),
            'Std_Correlation':    round(float(np.nanstd(vals)),    3),
        })
    else:
        metrics.update({
            'Mean_Correlation':   np.nan,
            'Median_Correlation': np.nan,
            'Std_Correlation':    np.nan,
        })

    return metrics

# =============================================================================
# HEMISPHERE COMPARISON
# =============================================================================

def analyze_hemisphere_comparisons(correlation_matrix, roi_names):
    """
    Compare intra-left, intra-right, and inter-hemisphere connectivity.
    Summary scalars use median; raw arrays passed to scipy for t-tests.
    """
    left_idx  = [i for i, n in enumerate(roi_names) if n.startswith('L_')]
    right_idx = [i for i, n in enumerate(roi_names) if n.startswith('R_')]

    if not left_idx or not right_idx:
        return None

    intra_left_vals  = correlation_matrix[np.ix_(left_idx,  left_idx)
                                          ][np.triu_indices(len(left_idx),  k=1)]
    intra_right_vals = correlation_matrix[np.ix_(right_idx, right_idx)
                                          ][np.triu_indices(len(right_idx), k=1)]

    inter_lr_vals = [
        correlation_matrix[roi_names.index(n), roi_names.index('R_' + n[2:])]
        for n in roi_names
        if n.startswith('L_') and 'R_' + n[2:] in roi_names
    ]

    t_lr,    p_lr    = stats.ttest_rel(intra_left_vals, intra_right_vals)
    t_intra, p_intra = stats.ttest_ind(
        np.concatenate([intra_left_vals, intra_right_vals]),
        inter_lr_vals,
    )

    return {
        'intra_left':          intra_left_vals,
        'intra_right':         intra_right_vals,
        'inter':               inter_lr_vals,
        'intra_left_median':   float(np.nanmedian(intra_left_vals)),
        'intra_right_median':  float(np.nanmedian(intra_right_vals)),
        'inter_median':        float(np.nanmedian(inter_lr_vals)),
        'p_left_vs_right':     p_lr,
        'p_intra_vs_inter':    p_intra,
        't_left_vs_right':     t_lr,
        't_intra_vs_inter':    t_intra,
    }

# =============================================================================
# QC PLOTS
# =============================================================================

def generate_qc_plots(corr_matrix, roi_names, output_dir, prefix,
                      intra_thresh=0.2, delta_thresh=0.05):
    """
    Four-panel QC figure:
      A — correlation matrix heatmap
      B — correlation value distribution
      C — hemisphere connectivity violins
      D — summary statistics table
    """
    fig = plt.figure(figsize=(18, 10))
    gs  = GridSpec(2, 3, figure=fig,
                   width_ratios=[1, 1, 1.5],
                   height_ratios=[1, 1.2],
                   wspace=0.3, hspace=0.4)

    # ── A: heatmap ────────────────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    img = ax1.imshow(corr_matrix, cmap='RdBu_r', vmin=-0.8, vmax=0.8)
    plt.colorbar(img, ax=ax1, fraction=0.04, pad=0.01)
    ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # ── B: distribution ───────────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    tril_mask   = np.tri(corr_matrix.shape[0], k=-1).astype(bool)
    corr_values = corr_matrix[tril_mask]
    ax2.hist(corr_values, bins=50, color='#4e79a7', edgecolor='white', alpha=0.8)
    ax2.axvline(float(np.nanmedian(corr_values)), color='#d62728',
                lw=1.5, ls='--', label=f'median={np.nanmedian(corr_values):.3f}')
    ax2.set_xlabel('Correlation value')
    ax2.set_ylabel('Count')
    ax2.legend(fontsize=8, frameon=False)
    ax2.text(-0.1, 1.1, 'B', transform=ax2.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # ── C: hemisphere violins ─────────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    hemi = analyze_hemisphere_comparisons(corr_matrix, roi_names) or {}

    left_data  = hemi.get('intra_left',  [np.nan])
    right_data = hemi.get('intra_right', [np.nan])
    inter_data = hemi.get('inter',       [np.nan])

    plot_df = pd.DataFrame({
        'Connectivity': np.concatenate([left_data, right_data, inter_data]),
        'Type': (  ['Intra-Left']  * len(left_data)
                 + ['Intra-Right'] * len(right_data)
                 + ['Inter-Hemi']  * len(inter_data)),
    })

    palette = ['#4e79a7', '#f28e2b', '#59a14f']
    sns.violinplot(x='Type', y='Connectivity', data=plot_df,
                   palette=palette, ax=ax3, inner=None, saturation=0.8)
    sns.stripplot(x='Type', y='Connectivity', data=plot_df,
                  palette=palette, ax=ax3, alpha=0.7, size=5, jitter=0.1)

    if len(left_data) == len(right_data):
        for lv, rv in zip(left_data, right_data):
            jit = np.random.uniform(-0.1, 0.1)
            ax3.plot([jit, 1 + jit], [lv, rv],
                     color='#7f7f7f', alpha=0.3, lw=1, ls='--')

    y_max = max(np.nanmax(left_data), np.nanmax(right_data), np.nanmax(inter_data))
    y_min = min(np.nanmin(left_data), np.nanmin(right_data), np.nanmin(inter_data))
    if 'p_left_vs_right' in hemi:
        ax3.text(0.5, y_max * 1.05, f"p={hemi['p_left_vs_right']:.3f}",
                 ha='center', fontsize=10)
    if 'p_intra_vs_inter' in hemi:
        ax3.text(1.5, y_max * 1.10, f"p={hemi['p_intra_vs_inter']:.3f}",
                 ha='center', fontsize=10)
    ax3.set_ylim(y_min * 1.1, y_max * 1.15)
    ax3.text(-0.1, 1.1, 'C', transform=ax3.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # ── D: table ──────────────────────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, :])
    ax4.axis('off')

    metrics     = calculate_network_metrics(corr_matrix) or {}
    spec        = analyze_specificity(corr_matrix, roi_names,
                                      intra_thresh, delta_thresh) or {}
    td          = spec.get('target_specificity', {})

    pos_frac    = float(np.mean(corr_values > 0))
    neg_frac    = float(np.mean(corr_values < 0))
    iqr         = float(np.nanpercentile(corr_values, 75)
                        - np.nanpercentile(corr_values, 25))

    def _fmt(v, decimals=3):
        return f"{v:.{decimals}f}" if np.isfinite(float(v)) else 'n/a'

    table_content = [
        # row 0 — section headers
        ["Network metrics", "", "Hemisphere results", "", "Matrix stats", ""],
        # row 1
        ["Top eigenvalue",      _fmt(metrics.get('Top_Eigenvalue',   np.nan), 2),
         "Left vs right p",     _fmt(hemi.get('p_left_vs_right',     np.nan)),
         "Min",                 _fmt(np.nanmin(corr_values))],
        # row 2
        ["Eigenvalue ratio",    _fmt(metrics.get('Eigenvalue_Ratio', np.nan)),
         "Intra vs inter p",    _fmt(hemi.get('p_intra_vs_inter',    np.nan)),
         "Max",                 _fmt(np.nanmax(corr_values))],
        # row 3
        ["Silhouette score",    _fmt(metrics.get('Silhouette_Score', np.nan)),
         "Median intra-left",   _fmt(np.nanmedian(left_data)),
         "Mean",                _fmt(np.nanmean(corr_values))],
        # row 4
        ["Davies-Bouldin",      _fmt(metrics.get('Davies_Bouldin',   np.nan)),
         "Median intra-right",  _fmt(np.nanmedian(right_data)),
         "Median",              _fmt(np.nanmedian(corr_values))],
        # row 5
        ["Calinski-Harabasz",   _fmt(metrics.get('Calinski_Harabasz', np.nan), 0),
         "Median inter-hemi",   _fmt(np.nanmedian(inter_data)),
         "Std",                 _fmt(np.nanstd(corr_values))],
        # row 6
        ["Positive corr.",      f"{pos_frac:.1%}",
         "Negative corr.",      f"{neg_frac:.1%}",
         "IQR",                 _fmt(iqr)],
        # row 7 — section header
        ["Specificity results", "", "", "", "", ""],
        # row 8
        ["Specific pair",       "homotopic same-network",
         "Intra (median)",      _fmt(td.get('Intra_Correlation',        np.nan)),
         "Intra threshold",     _fmt(td.get('Intra_Threshold', intra_thresh))],
        # row 9
        ["Non-specific pair",   "heterotopic contralateral",
         "Inter hetero (med.)", _fmt(td.get('Inter_Correlation_Hetero', np.nan)),
         "Delta",               _fmt(td.get('Delta',                    np.nan))],
        # row 10
        ["", "", "", "", "Classification", td.get('Category', 'N/A')],
    ]

    table = ax4.table(cellText=table_content, loc='center', cellLoc='center',
                      colWidths=[0.20, 0.15, 0.20, 0.15, 0.15, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    for (row, col), cell in table.get_celld().items():
        if row in [0, 7]:
            cell.set_facecolor('#404040')
            cell.set_text_props(color='white', weight='bold')
        elif row == 1 and col == 3:
            p = hemi.get('p_left_vs_right', 1.0)
            cell.set_facecolor('#d62728' if p < 0.05 else '#2ca02c')
            cell.set_text_props(color='white')
        elif row == 2 and col == 3:
            p = hemi.get('p_intra_vs_inter', 1.0)
            cell.set_facecolor('#2ca02c' if p < 0.05 else '#d62728')
            cell.set_text_props(color='white')
        elif row == 10 and col == 5:
            cat = td.get('Category', 'N/A')
            if cat in _CATEGORY_COLORS:
                cell.set_facecolor(_CATEGORY_COLORS[cat])
                cell.set_text_props(
                    color='white' if cat != 'No' else 'black',
                    weight='bold',
                )
        elif row > 0 and col in [0, 2, 4]:
            cell.set_facecolor('#f0f0f0')

    ax4.text(-0.05, 1.05, 'D', transform=ax4.transAxes,
             fontsize=16, fontweight='bold', va='top')

    plt.savefig(opj(output_dir, f'{prefix}_qc_report.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
    return metrics, hemi

# =============================================================================
# MAIN QC PIPELINE
# =============================================================================

def fMRI_QC_matrix(path_func, dir_prepro_orig, intra_thresh, delta_thresh,
                   RS, nb_run, diary_file):
    """
    Run per-subject, per-run FC QC classification.

    Parameters
    ----------
    intra_thresh : float
        Minimum median homotopic correlation to pass (want > this).
    delta_thresh : float
        Minimum intra − inter_hetero gap to be Specific (want > this).
    """
    run_cmd.msg('## Working on fMRI QC matrix analysis ##', diary_file, 'HEADER')

    out_results = opj(path_func, 'QC')
    os.makedirs(out_results, exist_ok=True)

    all_classifications = []

    for i in range(int(nb_run)):
        root_RS     = extract_filename(RS[i])
        matrix_file = opj(dir_prepro_orig, 'Stats', 'Correl_matrix',
                          f'EDNIxCSCLR_2_run_{i}_correlation_matrix.csv')

        if not ope(matrix_file):
            run_cmd.msg(f'WARNING: Missing matrix file {matrix_file}',
                        diary_file, 'WARNING')
            continue

        try:
            corr, roi_names = load_and_validate_matrix(matrix_file)
            corr, roi_names = filter_bilateral_rois(corr, roi_names)

            spec  = analyze_specificity(corr, roi_names, intra_thresh, delta_thresh)
            td    = spec.get('target_specificity', {})
            cat   = td.get('Category', 'No')

            record = {
                'run':                      root_RS,
                'category':                 cat,
                'intra_correlation':        td.get('Intra_Correlation',        np.nan),
                'inter_correlation':        td.get('Inter_Correlation',        np.nan),
                'inter_correlation_hetero': td.get('Inter_Correlation_Hetero', np.nan),
                'delta':                    td.get('Delta',                    np.nan),
                'intra_threshold':          intra_thresh,
                'delta_threshold':          delta_thresh,
                'specific_pair':            td.get('Specific_ROI_pair', 'Not found'),
            }
            all_classifications.append(record)

            _emoji  = {'Specific': '🟢', 'Unspecific': '🟠',
                       'No': '🔵', 'Spurious': '🔴'}.get(cat, '⚪')
            _colour = {'Specific': 'OKGREEN', 'Unspecific': 'WARNING',
                       'No': 'ENDC',          'Spurious': 'FAIL'}.get(cat, 'ENDC')

            run_cmd.msg(
                f"  → {_emoji} Run {root_RS}: {cat} "
                f"(intra={record['intra_correlation']:.3f}, "
                f"inter_hetero={record['inter_correlation_hetero']:.3f}, "
                f"delta={record['delta']:.3f})",
                diary_file, _colour,
            )

            metrics, hemi = generate_qc_plots(
                corr, roi_names, out_results, root_RS,
                intra_thresh, delta_thresh,
            )

            # ── save JSON ─────────────────────────────────────────────────────
            def _serialise(obj):
                if isinstance(obj, dict):    return {k: _serialise(v) for k, v in obj.items()}
                if isinstance(obj, list):    return [_serialise(v) for v in obj]
                if isinstance(obj, np.ndarray): return obj.tolist()
                return obj

            with open(opj(out_results, f'{root_RS}_full_results.json'), 'w') as fh:
                json.dump(_serialise({
                    'network_metrics':     metrics,
                    'hemisphere_results':  hemi,
                    'specificity_results': spec,
                    'timestamp':           str(datetime.datetime.now()),
                }), fh, indent=2)

            pd.DataFrame(corr, index=roi_names, columns=roi_names).to_csv(
                opj(out_results, f'{root_RS}_correlation_matrix.csv'))

            run_cmd.msg(f'\n{root_RS} analysis complete\n', diary_file, 'OKGREEN')

        except Exception as e:
            run_cmd.msg(f'Error processing {matrix_file}: {e}', diary_file, 'FAIL')
            continue

    # ── summary ───────────────────────────────────────────────────────────────
    if not all_classifications:
        run_cmd.msg('No runs classified.', diary_file, 'WARNING')
        return

    pd.DataFrame(all_classifications).to_csv(
        opj(out_results, 'all_runs_classification.csv'), index=False)

    _bar  = '=' * 60
    _dash = '-' * 40
    run_cmd.msg(f'\n{_bar}', diary_file, 'OKGREEN')
    run_cmd.msg('CLASSIFICATION SUMMARY', diary_file, 'OKGREEN')
    run_cmd.msg(f'  Thresholds: intra > {intra_thresh},  delta > {delta_thresh}',
                diary_file, 'OKGREEN')
    run_cmd.msg(_dash, diary_file, 'OKGREEN')

    n_total = len(all_classifications)
    _emoji_map  = {'Specific': '🟢', 'Unspecific': '🟠', 'No': '🔵', 'Spurious': '🔴'}
    _colour_map = {'Specific': 'OKGREEN', 'Unspecific': 'WARNING',
                   'No': 'ENDC',          'Spurious': 'FAIL'}

    for cat in _CATEGORIES:
        n   = sum(1 for c in all_classifications if c['category'] == cat)
        pct = 100.0 * n / n_total
        run_cmd.msg(f"  {_emoji_map[cat]} {cat}: {n} runs ({pct:.1f}%)",
                    diary_file, _colour_map[cat])

    run_cmd.msg(_bar, diary_file, 'OKGREEN')
    run_cmd.msg('QC analysis completed successfully', diary_file, 'OKGREEN')