#!/usr/bin/env python
# ═══════════════════════════════════════════════════════════════════════════════
# fMRI PARAMETERS & SP_SPECIFICITY DASHBOARD
# =============================================================================
# This script:
#   1. Extracts fMRI parameters (TR, n_timepoints, TE) from BIDS datasets
#   2. Generates a dashboard with:
#      - Row 1: Violin plots of fMRI parameters by BIDS directory (colored by BIDS)
#      - Row 2: Correlations of fMRI parameters vs sp_specificity_index
#        (partial correlation controlling for BIDS effects, colored by BIDS)
# ═══════════════════════════════════════════════════════════════════════════════

import os
import sys
import glob
import json
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats as _stats

# Add paths
sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/Archives/')
sys.path.insert(0, '//')

from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
    _parse_bids_entity, _parse_run, _spco, PALETTE, _phylo_sort,
    get_atlas_label_path, find_species_path, collect_multi_species
)

opj = os.path.join

# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────

SPECIES_COLORS = {
    "Mouse": "#1A5276", "Rat": "#2E86C1", "Bat": "#00A896",
    "Mouselemur": "#58D68D", "Marmoset": "#F4D03F", "Macaque": "#F5B041",
    "Human": "#D35400", "Dog": "#8E44AD", "Cat": "#C0392B", "Pig": "#6E2C00",
}

# BIDS colors - consistent mapping
BIDS_COLORS = {
    'BIDS_Grandjean': '#0072B2',
    'BIDS_Grandjean2': '#E69F00',
    'BIDS_Boch_K9': '#009E73',
    'BIDS_Tian': '#CC79A7',
    'BIDS_Garin': '#56B4E9',
    'BIDS_Tremblay': '#D55E00',
    'BIDS_BenHamed': '#F0E442',
    'BIDS_Zhu_Garin': '#000000',
    'BIDS_Merida': '#1A5276',
    'BIDS_Park': '#2E86C1',
    'BIDS_Castrillon': '#58D68D',
}

PHYLO_ORDER = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
               "Macaque", "Human", "Dog", "Cat", "Pig"]

# Macaque/BenHamed default parameters
MACAQUE_BENHAMED_PARAMS = {
    'TR': 2.0,
    'TE': 0.030,
    'TotalReadoutTime': 0.045,
    'PhaseEncodingDirection': 'j-',
    'SliceEncodingDirection': 'k',
}


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: fMRI PARAMETER EXTRACTION
# ═══════════════════════════════════════════════════════════════════════════════

def extract_fmri_parameters(
        bids_dir,
        endfmri="*_bold.nii.gz",
        endjson="*_bold.json",
        list_to_keep=None,
        list_to_remove=None,
        ntimepoint_threshold=100,
        species=None,
        bids_label=None
):
    """Extract fMRI parameters from BIDS dataset with auto-fixes."""
    list_to_keep = list_to_keep or []
    list_to_remove = list_to_remove or []
    records = []

    print(f"  [extract] Scanning: {bids_dir}")

    for sub_dir in sorted(glob.glob(opj(bids_dir, "sub-*"))):
        sub = _parse_bids_entity(sub_dir, "sub")
        if sub is None:
            continue

        ses_dirs = sorted(glob.glob(opj(sub_dir, "ses-*"))) or [sub_dir]

        for ses_dir in ses_dirs:
            ses = _parse_bids_entity(ses_dir, "ses") or "1"

            if list_to_keep and (sub, ses) not in list_to_keep:
                continue
            if list_to_remove and (sub, ses) in list_to_remove:
                continue

            func_dir = opj(ses_dir, "func")
            if not os.path.exists(func_dir):
                continue

            bold_files = sorted(glob.glob(opj(func_dir, endfmri)))
            json_files = sorted(glob.glob(opj(func_dir, endjson)))

            if not bold_files:
                continue

            for bold_path in bold_files:
                base_name = bold_path.replace('.nii.gz', '').replace('.nii', '')
                json_path = base_name + '.json'

                if not os.path.exists(json_path):
                    json_path = None
                    for jf in json_files:
                        if os.path.basename(jf).replace('.json', '') in os.path.basename(bold_path):
                            json_path = jf
                            break

                run = _parse_run(bold_path)

                try:
                    img = nib.load(bold_path)
                    n_timepoints = img.shape[3] if len(img.shape) >= 4 else 1
                except Exception as e:
                    print(f"    Warning: Could not load {bold_path}: {e}")
                    continue

                if n_timepoints < ntimepoint_threshold:
                    print(f"    [skip] {sub}: {n_timepoints} TRs < threshold")
                    continue

                info = {}
                if json_path and os.path.exists(json_path):
                    try:
                        with open(json_path, 'r') as f:
                            info = json.load(f)
                    except Exception as e:
                        print(f"    Warning: Could not parse JSON: {e}")

                # Extract TR with auto-fix for Rat data
                tr = None
                if "RepetitionTime" in info:
                    tr = float(info["RepetitionTime"])
                    # Fix Rat TR=0.001s → 1.0s
                    if tr is not None and tr < 0.01:
                        is_rat = species == 'Rat' or 'Rat' in str(bids_dir)
                        if is_rat and n_timepoints >= 100:
                            print(f"    [fix] {sub} ses-{ses}: TR {tr}s → 1.0s")
                            tr = 1.0
                elif "SliceTiming" in info:
                    slice_timing = sorted(info["SliceTiming"])
                    if len(slice_timing) > 1:
                        intervals = [slice_timing[i + 1] - slice_timing[i]
                                     for i in range(len(slice_timing) - 1)]
                        tr = ((sum(intervals) / (len(slice_timing) - 1)) * len(slice_timing))

                # Extract other parameters
                te = float(info["EchoTime"]) if "EchoTime" in info else None
                total_readout = float(info["TotalReadoutTime"]) if "TotalReadoutTime" in info else None
                phase_encoding = info.get("PhaseEncodingDirection")
                slice_encoding = info.get("SliceEncodingDirection")
                slice_timing = info.get("SliceTiming")
                task = info.get("TaskName")

                # Fill Macaque/BenHamed missing params
                is_benhamed = bids_label == 'BIDS_BenHamed' or 'BIDS_BenHamed' in str(bids_dir)
                if is_benhamed and species == 'Macaque':
                    if tr is None or np.isnan(tr) or tr < 0.01:
                        tr = MACAQUE_BENHAMED_PARAMS['TR']
                        print(f"    [fill] {sub} ses-{ses}: TR set to {tr}s")
                    if te is None or np.isnan(te):
                        te = MACAQUE_BENHAMED_PARAMS['TE']
                    if total_readout is None or np.isnan(total_readout):
                        total_readout = MACAQUE_BENHAMED_PARAMS['TotalReadoutTime']
                    if phase_encoding is None or phase_encoding == '':
                        phase_encoding = MACAQUE_BENHAMED_PARAMS['PhaseEncodingDirection']
                    if slice_encoding is None or slice_encoding == '':
                        slice_encoding = MACAQUE_BENHAMED_PARAMS['SliceEncodingDirection']

                duration = n_timepoints * tr if tr else None

                record = {
                    'subject': sub,
                    'session': ses,
                    'run': run,
                    'bold_path': bold_path,
                    'TR': tr,
                    'n_timepoints': n_timepoints,
                    'duration_seconds': duration,
                    'duration_minutes': duration / 60 if duration else None,
                    'TE': te,
                    'EchoTime': te,
                    'TotalReadoutTime': total_readout,
                    'PhaseEncodingDirection': phase_encoding,
                    'SliceEncodingDirection': slice_encoding,
                    'SliceTiming': str(slice_timing) if slice_timing else None,
                    'n_slices': len(slice_timing) if slice_timing else None,
                    'AcquisitionTime': info.get("AcquisitionTime"),
                    'TaskName': task,
                    'run_label': run,
                }

                if species:
                    record['species'] = species
                if bids_label:
                    record['bids_dir'] = bids_label

                records.append(record)

    df = pd.DataFrame(records)

    if not df.empty:
        print(f"\n  [extract] {len(df)} runs from {bids_dir}")
        if 'TR' in df.columns and not df['TR'].isna().all():
            valid_tr = df[df['TR'] > 0.01]['TR'] if 'TR' in df.columns else df['TR']
            if not valid_tr.empty:
                print(f"    TR: {valid_tr.min():.3f}s - {valid_tr.max():.3f}s (mean: {valid_tr.mean():.3f}s)")
        print(f"    n_timepoints: {df['n_timepoints'].min()} - {df['n_timepoints'].max()}")

    return df


def collect_fmri_parameters_from_config(species_config, ntimepoint_threshold=100):
    """Collect fMRI parameters across multiple species."""
    all_dfs = []

    for species, cfg in species_config.items():
        bids_dirs = cfg.get('bids_dirs', [])
        if isinstance(bids_dirs, str):
            bids_dirs = [bids_dirs]

        lk = cfg.get('list_to_keep', [])
        lr = cfg.get('list_to_remove', [])

        print(f"\n{'=' * 60}\n  Species: {species}\n{'=' * 60}")

        for bids_dir in bids_dirs:
            bids_label = os.path.basename(bids_dir.rstrip('/'))
            print(f"\n  Processing: {bids_label}")

            df = extract_fmri_parameters(
                bids_dir,
                list_to_keep=lk,
                list_to_remove=lr,
                ntimepoint_threshold=ntimepoint_threshold,
                species=species,
                bids_label=bids_label
            )
            if not df.empty:
                all_dfs.append(df)

    if all_dfs:
        df = pd.concat(all_dfs, ignore_index=True)
        print(f"\n{'=' * 60}")
        print(f"  Total: {len(df)} runs across {df['species'].nunique()} species")
        if 'TR' in df.columns:
            valid_tr = df[df['TR'] > 0.01]['TR'] if 'TR' in df.columns else df['TR']
            if not valid_tr.empty:
                print(f"  Global TR: {valid_tr.mean():.3f}s ± {valid_tr.std():.3f}s")
        print(f"{'=' * 60}")
        return df

    print("\n  No runs found.")
    return pd.DataFrame()


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: VIOLIN PLOT HELPER
# ═══════════════════════════════════════════════════════════════════════════════

def _violin_strip_quartiles(ax, data, x_pos, color, width=0.4):
    """Violin + quartile plot matching plot_morphometry_intra_bids style."""
    data = np.asarray(data, dtype=float)
    data = data[~np.isnan(data)]
    if len(data) == 0:
        return
    if np.std(data) == 0 or len(data) < 2:
        ax.hlines(data[0], x_pos - width / 2, x_pos + width / 2, colors=color, linewidths=2, alpha=0.6)
        return
    from scipy.stats import gaussian_kde
    try:
        kde = gaussian_kde(data, bw_method="scott")
    except Exception:
        return
    y_g = np.linspace(data.min(), data.max(), 200)
    dens = kde(y_g)
    dens = dens / dens.max() * (width / 2)
    ax.fill_betweenx(y_g, x_pos - dens, x_pos + dens, color=color, alpha=0.25, linewidth=0)
    ax.plot(x_pos - dens, y_g, color=color, lw=0.8, alpha=0.5)
    ax.plot(x_pos + dens, y_g, color=color, lw=0.8, alpha=0.5)
    for q, lw in zip(np.percentile(data, [25, 50, 75]), [1, 2, 1]):
        dq = float(kde(q)) / kde(y_g).max() * (width / 2)
        ax.hlines(q, x_pos - dq, x_pos + dq, colors=color, linewidths=lw, zorder=3)


# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: MAIN DASHBOARD - fMRI PARAMETERS vs SP_SPECIFICITY
# ═══════════════════════════════════════════════════════════════════════════════
def plot_fmri_params_with_specificity_dashboard(
        df_qc,
        df_fmri,
        output_dir,
        fmri_params=None,
        target_col="sp_specificity_index",
        bids_col=None,
        figsize=(18, 12)
):
    """
    Create a dashboard with:
    - Row 1: Violin plots of fMRI parameters by BIDS directory (colored by BIDS)
    - Row 2: Correlations of each fMRI parameter vs sp_specificity_index
      (partial correlation controlling for BIDS effects, colored by BIDS)
    """
    if df_fmri.empty:
        print("  No fMRI parameters to plot.")
        return None

    if fmri_params is None:
        fmri_params = ['TR', 'n_timepoints', 'TE']

    # Filter to only params that exist
    existing_params = [p for p in fmri_params if p in df_fmri.columns and not df_fmri[p].isna().all()]
    if not existing_params:
        print("  No valid fMRI parameters found.")
        return None

    # ─── Merge QC and fMRI data ──────────────────────────────────────────────
    # Determine merge columns that exist in both DataFrames
    possible_merge_cols = ['subject', 'session', 'species', 'bids_dir']
    merge_cols = []
    for col in possible_merge_cols:
        if col in df_qc.columns and col in df_fmri.columns:
            merge_cols.append(col)

    if not merge_cols:
        print("  No common columns to merge on!")
        return None

    print(f"  Merging on: {merge_cols}")

    df_merged = pd.merge(df_qc, df_fmri, on=merge_cols, how='inner')

    if df_merged.empty:
        print("  No matching data after merge!")
        return None

    print(f"  Merged data: {len(df_merged)} rows")

    # Filter suspicious TR values
    if 'TR' in df_merged.columns:
        before = len(df_merged)
        df_merged = df_merged[df_merged['TR'] > 0.01]
        if len(df_merged) < before:
            print(f"  [filter] Removed {before - len(df_merged)} runs with TR < 0.01s")

    # Check if target exists
    if target_col not in df_merged.columns:
        print(f"  Target '{target_col}' not found in merged data!")
        return None

    # ─── ORDER BIDS DIRECTORIES BY SPECIES ──────────────────────────────────
    # Get species for each BIDS directory
    bids_species = {}
    for bd in df_merged['bids_dir'].unique():
        bids_species[bd] = df_merged[df_merged['bids_dir'] == bd]['species'].iloc[0]

    # Sort BIDS directories: first by species order, then by BIDS name
    def sort_key(bd):
        sp = bids_species.get(bd, '')
        sp_idx = PHYLO_ORDER.index(sp) if sp in PHYLO_ORDER else 999
        return (sp_idx, bd)

    bids_dirs = sorted(df_merged['bids_dir'].unique(), key=sort_key)
    print(f"\n  BIDS directories (ordered by species):")
    for bd in bids_dirs:
        print(f"    {bd} ({bids_species.get(bd, '')})")

    # ─── BIDS colors ──────────────────────────────────────────────────────────
    if bids_col is None:
        # Use consistent colors from BIDS_COLORS dict
        bids_colors = {}
        available_colors = list(PALETTE)
        for i, bd in enumerate(bids_dirs):
            if bd in BIDS_COLORS:
                bids_colors[bd] = BIDS_COLORS[bd]
            else:
                bids_colors[bd] = available_colors[i % len(available_colors)]
    else:
        bids_colors = bids_col

    # ─── CREATE FIGURE: 2 ROWS, 3 COLUMNS ─────────────────────────────────────
    n_cols = len(existing_params)
    n_rows = 2

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # ─── ROW 1: VIOLIN PLOTS OF fMRI PARAMETERS (colored by BIDS) ────────────
    for idx, param in enumerate(existing_params):
        ax = axes[idx]

        # Collect data by BIDS directory
        bids_positions = {}
        for i, bd in enumerate(bids_dirs):
            vals = df_merged[df_merged['bids_dir'] == bd][param].dropna().values
            if len(vals) > 0:
                bids_positions[i] = vals

        # Plot violin for each BIDS directory
        for i, bd in enumerate(bids_dirs):
            if i not in bids_positions:
                continue
            vals = bids_positions[i]
            color = bids_colors.get(bd, "#888")

            _violin_strip_quartiles(ax, vals, i, color, width=0.4)

            jitter = np.random.normal(0, 0.04, len(vals))
            ax.scatter([i] * len(vals) + jitter, vals,
                       color=color, s=18, alpha=0.5, edgecolors='k', linewidths=0.3)

        # X-axis labels
        x_labels = []
        for bd in bids_dirs:
            sp = bids_species.get(bd, '')
            short_bd = bd.replace('BIDS_', '').replace('_', ' ')
            x_labels.append(f"{short_bd}\n({sp})")

        ax.set_xticks(range(len(bids_dirs)))
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)

        # Y-axis formatting
        if param == 'TR':
            ylabel = 'TR (s)'
            ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
        elif param == 'n_timepoints':
            ylabel = 'Number of TRs'
            ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.0f'))
        elif param == 'TE':
            ylabel = 'TE (s)'
            ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.4f'))
        else:
            ylabel = param
            ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))

        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(param, fontweight='bold', fontsize=11, pad=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.6, len(bids_dirs) - 0.4)

        # ─── BIDS legend (only first subplot) ──────────────────────────────
        if idx == 0:
            bd_handles = []
            for bd in bids_dirs:
                color = bids_colors.get(bd, "#888")
                short_bd = bd.replace('BIDS_', '').replace('_', ' ')
                bd_handles.append(mpatches.Patch(facecolor=color, label=short_bd,
                                                 edgecolor='k', linewidth=0.3))
            ax.legend(handles=bd_handles, fontsize=6, frameon=False,
                      bbox_to_anchor=(1.05, 1), loc='upper left', title='BIDS dir')

    # ─── ROW 2: PARTIAL CORRELATIONS (controlling for BIDS, colored by BIDS) ──
    corr_results = {}

    for param in existing_params:
        if param not in df_merged.columns:
            continue

        # Remove NaN
        valid = df_merged[[param, target_col, 'bids_dir']].dropna()

        # Check if we have enough data
        if len(valid) < 5:
            print(f"  {param}: insufficient data (n={len(valid)})")
            corr_results[param] = {'r': np.nan, 'p': np.nan, 'n': len(valid), 'r2': np.nan}
            continue

        # Check if there's variance
        if valid[param].nunique() <= 1:
            print(f"  {param}: constant value, skipping correlation")
            corr_results[param] = {'r': np.nan, 'p': np.nan, 'n': len(valid), 'r2': np.nan}
            continue

        # ─── Partial correlation controlling for BIDS ──────────────────────
        # Create dummy variables for BIDS
        bids_dummies = pd.get_dummies(valid['bids_dir'], prefix='bids', drop_first=True)

        # Get residuals
        X = bids_dummies.values
        X = np.column_stack([np.ones(len(X)), X])

        x = valid[param].values
        y = valid[target_col].values

        try:
            # Residual of x
            beta_x = np.linalg.lstsq(X, x, rcond=None)[0]
            resid_x = x - X @ beta_x

            # Residual of y
            beta_y = np.linalg.lstsq(X, y, rcond=None)[0]
            resid_y = y - X @ beta_y

            # Correlation of residuals
            r_partial, p_partial = _stats.pearsonr(resid_x, resid_y)

            # Calculate R² from the partial correlation
            r2_partial = r_partial ** 2

            # Also calculate the effect size (partial η² equivalent)
            # This gives the proportion of variance explained by the parameter
            # after accounting for BIDS effects
            try:
                # Alternative: Calculate from residuals
                ss_total = np.sum((resid_y - np.mean(resid_y)) ** 2)
                ss_residual = np.sum((resid_y - (resid_x * r_partial)) ** 2)
                r2_from_ss = 1 - (ss_residual / ss_total) if ss_total > 0 else 0
                # Use the simpler R² from correlation
                r2_partial = r_partial ** 2
            except:
                r2_partial = r_partial ** 2

        except Exception as e:
            print(f"  {param}: partial correlation failed: {e}")
            r_partial, p_partial = np.nan, np.nan
            r2_partial = np.nan

        # ─── Regular correlation (for comparison) ──────────────────────────
        r_regular, p_regular = _stats.spearmanr(x, y)
        r2_regular = r_regular ** 2 if not np.isnan(r_regular) else np.nan

        corr_results[param] = {
            'r_partial': r_partial,
            'p_partial': p_partial,
            'r2_partial': r2_partial,
            'r_regular': r_regular,
            'p_regular': p_regular,
            'r2_regular': r2_regular,
            'n': len(valid)
        }

        print(f"  {param}: partial r={r_partial:.3f} (R²={r2_partial:.3f}, p={p_partial:.4f}), regular r={r_regular:.3f} (R²={r2_regular:.3f})")

    # Plot correlations in second row
    for idx, param in enumerate(existing_params):
        ax = axes[idx + n_cols]  # Second row

        if param not in corr_results:
            ax.set_visible(False)
            continue

        # Get data
        valid = df_merged[[param, target_col, 'species', 'bids_dir']].dropna()

        # ─── Plot by BIDS (colored by BIDS) ──────────────────────────────────
        for bd in bids_dirs:
            bd_df = valid[valid['bids_dir'] == bd]
            if not bd_df.empty:
                color = bids_colors.get(bd, "#888")
                ax.scatter(bd_df[param], bd_df[target_col],
                           color=color, s=40, alpha=0.8, label=bd if idx == 0 else "",
                           edgecolors='k', linewidths=0.5)

        # ─── Per-BIDS regression lines ──────────────────────────────────────
        per_bids_r2 = {}
        for bd in bids_dirs:
            bd_df = valid[valid['bids_dir'] == bd]
            if len(bd_df) >= 3 and bd_df[param].nunique() > 1 and bd_df[target_col].nunique() > 1:
                color = bids_colors.get(bd, "#888")
                x = bd_df[param].values
                y = bd_df[target_col].values
                try:
                    slope, intercept, r_val, p_val, se = _stats.linregress(x, y)
                    x_line = np.linspace(x.min(), x.max(), 100)
                    y_line = intercept + slope * x_line
                    ax.plot(x_line, y_line, '-', color=color, lw=1.5, alpha=0.6)
                    # Store R² for per-BIDS regression
                    per_bids_r2[bd] = r_val ** 2
                except:
                    pass

        # ─── Overall regression line ────────────────────────────────────────
        x_all = valid[param].values
        y_all = valid[target_col].values
        overall_r2 = np.nan
        if len(x_all) > 2 and valid[param].nunique() > 1 and valid[target_col].nunique() > 1:
            try:
                slope, intercept, r_val, p_val, se = _stats.linregress(x_all, y_all)
                x_line = np.linspace(x_all.min(), x_all.max(), 100)
                y_line = intercept + slope * x_line
                ax.plot(x_line, y_line, 'k--', lw=1.5, alpha=0.6)
                overall_r2 = r_val ** 2
            except:
                pass

        # ─── Stats: Show PARTIAL correlation with R² ────────────────────────
        r = corr_results[param]['r_partial']
        p = corr_results[param]['p_partial']
        r2 = corr_results[param]['r2_partial']
        n = corr_results[param]['n']

        if not np.isnan(r) and not np.isnan(p):
            signif = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            # Show both ρ (correlation) and R²
            ax.set_title(f"{target_col} vs {param}\nρ={r:.3f}  R²={r2:.3f}  p={p:.4f}{signif}\nn={n}",
                         fontsize=9, fontweight='bold')
        else:
            ax.set_title(f"{target_col} vs {param}\nNo valid correlation",
                         fontsize=10, fontweight='bold')

        # ─── Per-BIDS correlation text with R² ──────────────────────────────
        corr_text = []
        for bd in bids_dirs:
            bd_df = valid[valid['bids_dir'] == bd]
            if len(bd_df) >= 3 and bd_df[param].nunique() > 1 and bd_df[target_col].nunique() > 1:
                try:
                    r_sp, p_sp = _stats.spearmanr(bd_df[param].values, bd_df[target_col].values)
                    r2_sp = r_sp ** 2 if not np.isnan(r_sp) else np.nan
                    if not np.isnan(r_sp):
                        sig = '***' if p_sp < 0.001 else '**' if p_sp < 0.01 else '*' if p_sp < 0.05 else ''
                        short_bd = bd.replace('BIDS_', '').replace('_', ' ')
                        corr_text.append(f"{short_bd}: ρ={r_sp:.3f} R²={r2_sp:.3f}{sig}")
                except:
                    pass

        if corr_text:
            # Add overall R² to the text box
            if not np.isnan(overall_r2):
                corr_text.insert(0, f"Overall R²={overall_r2:.3f}")
                corr_text.insert(1, "---")

            # Truncate if too many
            if len(corr_text) > 8:
                corr_text = corr_text[:8] + [f"... ({len(corr_text)-8} more)"]

            ax.text(0.02, 0.98, '\n'.join(corr_text), transform=ax.transAxes,
                    va='top', fontsize=6, bbox=dict(boxstyle='round',
                                                    facecolor='white', alpha=0.9))

        # ─── Formatting ──────────────────────────────────────────────────────
        if param == 'TR':
            xlabel = 'TR (s)'
        elif param == 'n_timepoints':
            xlabel = 'Number of TRs'
        elif param == 'TE':
            xlabel = 'TE (s)'
        else:
            xlabel = param

        ax.set_xlabel(xlabel, fontsize=8)
        if idx == 0:
            ax.set_ylabel(target_col, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # ─── BIDS legend (only first subplot in second row) ────────────────
        if idx == 0:
            bd_handles = []
            for bd in bids_dirs:
                color = bids_colors.get(bd, "#888")
                short_bd = bd.replace('BIDS_', '').replace('_', ' ')
                bd_handles.append(mpatches.Patch(facecolor=color, label=short_bd,
                                                 edgecolor='k', linewidth=0.3))
            ax.legend(handles=bd_handles, fontsize=6, frameon=False,
                      bbox_to_anchor=(1.05, 1), loc='upper left', title='BIDS dir')

    plt.tight_layout()

    output_path = opj(output_dir, "FMRI_params_with_specificity_dashboard.png")
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"\n  [plot] {output_path}")

    # ─── Save correlation results with R² ─────────────────────────────────────
    corr_list = []
    for param, info in corr_results.items():
        corr_list.append({
            'fmri_param': param,
            'partial_r': info.get('r_partial', np.nan),
            'partial_r2': info.get('r2_partial', np.nan),  # Added R²
            'partial_p': info.get('p_partial', np.nan),
            'regular_r': info.get('r_regular', np.nan),
            'regular_r2': info.get('r2_regular', np.nan),  # Added R²
            'regular_p': info.get('p_regular', np.nan),
            'n': info.get('n', 0)
        })

    if corr_list:
        corr_df = pd.DataFrame(corr_list)
        corr_path = opj(output_dir, "fmri_specificity_correlations.csv")
        corr_df.to_csv(corr_path, index=False)
        print(f"  [export] {corr_path}")

    return output_path
# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":

    # ─── CONFIGURATION ──────────────────────────────────────────────────────────
    ATLAS_LIB = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
    ATLAS_NAME = 'EDNIxCSC'
    OUT_DIR = '/scratch2/EDNiX/results/EDNiX_Figures/Suppl_Figure_4-5-6_QC'

    # ─── SPECIES CONFIGURATION ─────────────────────────────────────────────────
    species_bids_dict = {
        'Rat': '/scratch2/EDNiX/Rat/BIDS_Grandjean',
        'Mouse': '/scratch2/EDNiX/Mouse/BIDS_Grandjean2',
        'Dog': '/scratch2/EDNiX/Dog/BIDS_Boch_K9',
        'Marmoset': '/scratch2/EDNiX/Marmoset/BIDS_Tian',
        'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin',
    }
    species_multi_bids = {
        'Macaque': [
            '/scratch2/EDNiX/Macaque/BIDS_Tremblay',
            '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
            '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',
        ],
        'Human': [
            '/scratch2/EDNiX/Human/BIDS_Merida',
            '/scratch2/EDNiX/Human/BIDS_Park',
            '/scratch2/EDNiX/Human/BIDS_Castrillon',
        ],
    }

    species_config = {}
    for species, bids_dir in species_bids_dict.items():
        species_config[species] = {
            'bids_dirs': [bids_dir],
            'list_to_keep': [],
            'list_to_remove': [],
        }
    for species, bids_dirs in species_multi_bids.items():
        species_config[species] = {
            'bids_dirs': bids_dirs,
            'list_to_keep': [],
            'list_to_remove': [],
        }

    # ─── RESOLVE ATLAS PATHS ───────────────────────────────────────────────────
    atlas_label_paths = {}
    species_atlas_fragments = {}

    for species in species_config:
        try:
            fragment = find_species_path(ATLAS_LIB, species)
            species_atlas_fragments[species] = fragment
            atlas_label_paths[species] = get_atlas_label_path(
                ATLAS_LIB, fragment, ATLAS_NAME, prefer_statslut=False)
            print(f"  [{species}] {atlas_label_paths[species]}")
        except (FileNotFoundError, ValueError) as e:
            print(f"  [WARN] {species}: {e}")

    # ─── STEP 1: COLLECT QC DATA ──────────────────────────────────────────────
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

    if df_qc is None or df_qc.empty:
        print("  No QC data collected - exiting.")
        sys.exit(1)

    print(f"  QC data collected: {len(df_qc)} rows, {len(df_qc.columns)} columns")

    # ─── PLOT GRANULARITY ──────────────────────────────────────────────────────
    # Same fix as Supp_Figure_4-6_QC.py: this used to average every session
    # of a subject into one point (_qc_id excluded 'session') before plotting
    # -- masks exactly the session/run-level variability a QC dashboard exists
    # to show. QC_PLOT_UNIT="session" (default): no averaging, one point per
    # (subject, session). "subject" restores the old behavior.
    QC_PLOT_UNIT = "session"

    if QC_PLOT_UNIT == "subject" and 'session' in df_qc.columns:
        _qc_id = [c for c in ('species', 'bids_dir', 'subject') if c in df_qc.columns]
        _qc_skip = set(_qc_id) | {'session', 'run'}
        _qc_num = [c for c in df_qc.columns
                   if c not in _qc_skip and pd.api.types.is_numeric_dtype(df_qc[c])]
        _qc_cat = [c for c in df_qc.columns
                   if c not in _qc_skip and c not in _qc_num]

        if _qc_id and _qc_num:
            df_qc_plot = df_qc.groupby(_qc_id, as_index=False)[_qc_num].mean()
            if _qc_cat:
                _cat_first = (df_qc.sort_values('session')
                              .groupby(_qc_id, as_index=False)[_qc_cat]
                              .first())
                df_qc_plot = df_qc_plot.merge(_cat_first, on=_qc_id, how='left')
            print(f"  [session_avg] qc: {len(df_qc)} -> {len(df_qc_plot)} rows "
                  f"(QC_PLOT_UNIT='subject')")
        else:
            df_qc_plot = df_qc
    else:
        df_qc_plot = df_qc
        print(f"  [no_avg] qc: {len(df_qc)} row(s) plotted as-is, one point per "
              f"(subject, session{', run' if 'run' in df_qc.columns else ''})")


    # ─── STEP 2: EXTRACT fMRI PARAMETERS ──────────────────────────────────────
    print("\n" + "=" * 64)
    print("  Extracting fMRI parameters...")
    print("=" * 64)

    df_fmri = collect_fmri_parameters_from_config(
        species_config,
        ntimepoint_threshold=100
    )

    if df_fmri.empty:
        print("  No fMRI parameters extracted.")
        sys.exit(1)

    # ─── STEP 3: GENERATE DASHBOARD ───────────────────────────────────────────
    print("\n" + "=" * 64)
    print("  Generating fMRI Parameters vs sp_specificity Dashboard...")
    print("=" * 64)

    # Generate the dashboard
    plot_fmri_params_with_specificity_dashboard(
        df_qc_plot,
        df_fmri,
        OUT_DIR,
        fmri_params=['TR', 'n_timepoints', 'TE'],
        target_col="sp_specificity_index",
        figsize=(18, 12)
    )

    # ─── STEP 4: EXPORT fMRI PARAMETERS ───────────────────────────────────────
    excel_path = opj(OUT_DIR, "fmri_parameters.xlsx")
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        df_fmri.to_excel(writer, sheet_name='fMRI_params', index=False)
    print(f"  [export] fMRI parameters → {excel_path}")

    # ─── FINAL REPORT ──────────────────────────────────────────────────────────
    print("\n" + "=" * 64)
    print("  ANALYSIS COMPLETE!")
    print("=" * 64)
    print(f"  Output directory: {OUT_DIR}")
    print("\n  Generated files:")
    print("    FMRI_params_with_specificity_dashboard.png")
    print("    fmri_parameters.xlsx")
    print("    fmri_specificity_correlations.csv")
    print("\n  Done.")