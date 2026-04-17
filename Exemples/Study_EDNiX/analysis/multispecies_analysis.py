"""
EDNiX BIDS Tools v3  —  Main analysis script
=============================================
Processes all species: surface, thickness, volume, QC
Produces DataFrames, exports and paper-quality figures.

Changes v3
----------
- Hemisphere detection from data (not from atlas name)
- atlas_level=1 filter applied in all plot calls (configurable via ATLAS_LEVEL)
- QC exported to two separate Excel files (anatomical / functional)
- QC figure separates anat and func on distinct row blocks
- plot_combo updated with hemisphere + atlas_level awareness

Figure organisation
-------------------
figures/
├── cross_species/
│   ├── morphometry_recap/
│   │   <mod>_<hemi>_<raw|log|norm>.png
│   ├── combo/
│   │   combo_<hemi>_<raw|log|norm>.png
│   └── qc_recap/
│       QC_cross_species_anat_and_func.png
└── per_bids/
    └── <species>_<bids_label>/
        <mod>_<hemi>.png  +  QC_dashboard.png
"""

import os
import sys
import math
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')
from Plotting.ednix_bids_tools import (
    get_atlas_label_path, find_species_path,
    extract_and_process_surfaces_wb, extract_and_process_thickness_wb,
    extract_corr_matrix_paths, collect_corr_matrices, load_corr_matrix,
    collect_multi_species,
    export_to_excel, export_to_csv, export_summary_stats,
    plot_morphometry_intra_bids,
    plot_qc_dashboard,
    _violin_strip_quartiles,
    PAPER_RC, PALETTE,
)

opj = os.path.join
rng = np.random.default_rng(42)

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

ATLAS_LIB   = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
ATLAS_NAME  = 'EDNIxCSC'   # bilateral atlas (no L_/R_ prefix)
ATLAS_LEVEL = 1            # which hierarchy level to show in figures (1-4)

OUT_DIR  = '/scratch2/EDNiX/results/multispecies_analysis'
FIG_DIR  = opj(OUT_DIR, 'figures')

SING_WB = (
    'vglrun singularity run '
    '--bind /srv/projects/,/srv/projects/easymribrain,/scratch2/,/scratch/ '
    '/home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/'
    'connectome_workbench_1.5.0-freesurfer-update.sif'
)

# Regions shown in every figure panel — must match exact names in xlsx/label files
PLOT_REGIONS = [
    'Isocortex', 'Allocortex', 'Periallocortex',
]

REGIONS_OF_INTEREST = PLOT_REGIONS  # filter applied during data collection

# Phylogenetic order — from most distant to closest to Human
PHYLO_ORDER = [
    'Bat', 'Rat', 'Mouse', 'Mouse_lemur',
    'Marmoset', 'Macaque', 'Human',
    'Dog', 'Pig', 'Cat',
]

def _phylo_sort(species_list):
    """Sort species by phylogenetic proximity (defined in PHYLO_ORDER)."""
    known   = [s for s in PHYLO_ORDER if s in species_list]
    unknown = sorted([s for s in species_list if s not in PHYLO_ORDER])
    return known + unknown

# ── QC column groups ─────────────────────────────────────────────────────────
# NMI is shown in its own shared row (func + anat side by side)
SHARED_NMI_COLS = ['func_nmi', 'anat_nmi', 'avg_snr_gray']
# Func QC (below NMI row): SNR gray, TSNR, Mean FD, gcor
FUNC_QC_COLS = [
    'func_TSNR_0',
    'func_mean_fd',
    'func_gcor',
]
# Anat QC (below NMI row): template correlation, SNR gray, noise estimate
ANAT_QC_COLS = [
    'anat_template_correlation',
    'anat_cortical_contrast',
]
# Clip percentile for noisy QC metrics
QC_CLIP_COLS = {
    'anat_fwhm':   99,
    'func_TSNR_0': 99,
}
# Network recap: scalars from full_results.json — in display order
FUNC_NET_COLS = [
    'net_mean_correlation', 'net_std_correlation',
    'net_eigenvalue_ratio', 'net_davies_bouldin',
    'net_inter_mean', 'net_intra_mean',
    'sp_specific_correlation', 'sp_nonspecific_correlation',
]
# Specificity categories (for bar plot)
FUNC_SPEC_COL = 'sp_category'

# Atlas config for correlation matrix plots
CORR_ATLAS_LEVEL = 2   # lvl2 for matrix plots
CORR_USE_LR      = True   # use LR version for lvl2

# Brain scaling config
SCALE_Y_REGION   = 'Allocortex'     # y-axis region for morphometry scaling
SCALE_X_REGION   = 'Isocortex'      # internal scaling denominator (plot3)
# FC-specific allometric scaling
FC_Y_ROI_A = 'L_Medial_prefrontal_cortex_(mPFC)'  # first ROI of the pair
FC_Y_ROI_B = 'L_Posterior_parietal_cortex'         # second ROI of the pair
FC_X_ROI_A   = 'L_Somatosensory_cortex'              # x-axis internal scaling
FC_X_ROI_B   = 'R_Somatosensory_cortex'              # x-axis internal scaling
# Body and brain weights per species (grams) — update as needed
SPECIES_BODY_WEIGHT_G = {
    'Mouse':       25,
    'Rat':         300,
    'Bat':         20,
    'Mouse_lemur': 60,
    'Marmoset':    350,
    'Macaque':     7000,
    'Human':       70000,
    'Dog':         15000,
    'Cat':         4000,
    'Pig':         50000,
}
SPECIES_BRAIN_WEIGHT_G = {
    'Mouse':       0.4,
    'Rat':         2.0,
    'Bat':         0.8,
    'Mouse_lemur': 1.8,
    'Marmoset':    7.7,
    'Macaque':     70,
    'Human':       1232,
    'Dog':         72,
    'Cat':         30,
    'Pig':         180,
}

# ─────────────────────────────────────────────────────────────────────────────
# SPECIES / BIDS DIRECTORIES
# ─────────────────────────────────────────────────────────────────────────────

species_bids_dict = {
    'Rat':         '/scratch2/EDNiX/Rat/BIDS_Gd',
    'Mouse':       '/scratch2/EDNiX/Mouse/BIDS_Gd',
    'Dog':         '/scratch2/EDNiX/Dog/BIDS_k9',
    'Marmoset':    '/scratch2/EDNiX/Marmoset/BIDS_NIH_MBM',
    'Mouse_lemur': '/scratch2/EDNiX/Mouse_lemur/BIDS_Garin',
}
species_multi_bids = {
    'Macaque': [
        '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
        '/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin',
    ],
    'Human': [
        '/scratch2/EDNiX/Human/ds004513-download',
    ],
}

# Flat list of (species, bids_dir) for per-BIDS figures
ALL_BIDS = (
    [(sp, bd) for sp, bd in species_bids_dict.items()] +
    [(sp, bd) for sp, bds in species_multi_bids.items() for bd in bds]
)


def _bids_label(bids_dir):
    return os.path.basename(bids_dir.rstrip('/'))


# ─────────────────────────────────────────────────────────────────────────────
# STEP 0  —  Pre-compute surfaces & thickness with wb_command (run once)
# ─────────────────────────────────────────────────────────────────────────────

for species, bids_dir in ALL_BIDS:
    print(f"\n  wb_command: {species} <- {_bids_label(bids_dir)}")
    extract_and_process_thickness_wb(
        SING_WB, bids_dir, regions_to_process=PLOT_REGIONS, overwrite=True)
    extract_and_process_surfaces_wb(
        SING_WB, bids_dir, regions_to_process=PLOT_REGIONS,
        should_average_hemispheres=False, overwrite=False)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1  —  Build species_config and resolve atlas paths
# ─────────────────────────────────────────────────────────────────────────────

species_config          = {}
atlas_label_paths       = {}
species_atlas_fragments = {}

for species, bids_dir in species_bids_dict.items():
    species_config[species] = {
        'bids_dirs':      [bids_dir],
        'list_to_keep':   [],
        'list_to_remove': [],
    }
for species, bids_dirs in species_multi_bids.items():
    species_config[species] = {
        'bids_dirs':      bids_dirs,
        'list_to_keep':   [],
        'list_to_remove': [],
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

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2  —  Collect all data
# ─────────────────────────────────────────────────────────────────────────────

data = collect_multi_species(
    species_config,
    regions_of_interest     = REGIONS_OF_INTEREST,
    extract                 = ('surface', 'volume', 'thickness', 'qc'),
    atlas_name              = ATLAS_NAME,
    atlas_label_paths       = atlas_label_paths,
    atlas_library_root      = ATLAS_LIB,
    species_atlas_fragments = species_atlas_fragments,
)

df_surface   = data.get('surface')
df_volume    = data.get('volume')
df_thickness = data.get('thickness')
df_qc        = data.get('qc')

for name, df in [('surface', df_surface), ('volume', df_volume),
                  ('thickness', df_thickness), ('qc', df_qc)]:
    n = len(df) if df is not None and not df.empty else 0
    print(f"  {name:10s}: {n} rows")

# ─────────────────────────────────────────────────────────────────────────────
# Global BIDS colour mapping — consistent across ALL modalities and figures
# ─────────────────────────────────────────────────────────────────────────────

_all_bids_dirs = sorted(set(
    list(df_surface['bids_dir'].unique()   if df_surface  is not None and not df_surface.empty  else []) +
    list(df_volume['bids_dir'].unique()    if df_volume   is not None and not df_volume.empty   else []) +
    list(df_thickness['bids_dir'].unique() if df_thickness is not None and not df_thickness.empty else []) +
    list(df_qc['bids_dir'].unique()        if df_qc       is not None and not df_qc.empty and
                                              'bids_dir' in df_qc.columns else [])
))
# Normalize all keys to basename so QC lookups always match
_all_bids_dirs = sorted(set(os.path.basename(str(bd).rstrip('/'))
                             for bd in _all_bids_dirs))
GLOBAL_BIDS_COL = {bd: PALETTE[i % len(PALETTE)] for i, bd in enumerate(_all_bids_dirs)}
print(f"\n  Global BIDS colour map: {GLOBAL_BIDS_COL}")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3  —  Detect available hemispheres per modality
# ─────────────────────────────────────────────────────────────────────────────

def _hemis_for(df):
    """
    Return list of hemispheres to iterate for this modality.
    Always includes 'bilateral' if present.
    If both left/right AND bilateral exist, returns all three.
    If only bilateral exists (e.g. EDNIxCSC volume atlas), returns ['bilateral'].
    """
    if df is None or df.empty or "hemisphere" not in df.columns:
        return ["bilateral"]
    hemis = sorted(df["hemisphere"].unique())
    return hemis  # may be ['bilateral'], ['left','right'], or all three


HEMIS_SURFACE   = _hemis_for(df_surface)
HEMIS_THICKNESS = _hemis_for(df_thickness)
HEMIS_VOLUME    = _hemis_for(df_volume)

print(f"\n  Hemispheres — surface:{HEMIS_SURFACE}  "
      f"thickness:{HEMIS_THICKNESS}  volume:{HEMIS_VOLUME}")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4  —  Export
# ─────────────────────────────────────────────────────────────────────────────

os.makedirs(OUT_DIR, exist_ok=True)

# Full results + CSVs
export_to_excel(data, opj(OUT_DIR, 'ednix_results.xlsx'))
export_to_csv(data,   opj(OUT_DIR, 'csv'))

# Summary stats — all modalities, hemisphere-aware
export_summary_stats(data, opj(OUT_DIR, 'ednix_summary_stats.xlsx'),
                     groupby=('species', 'bids_dir'))

# QC — single Excel file with anat / func / network tabs
if df_qc is not None and not df_qc.empty:
    id_cols = [c for c in ('species', 'subject', 'session', 'bids_dir')
               if c in df_qc.columns]

    anat_cols_present = [c for c in ANAT_QC_COLS  if c in df_qc.columns]
    func_cols_present = [c for c in FUNC_QC_COLS  if c in df_qc.columns]
    net_cols_present  = [c for c in FUNC_NET_COLS if c in df_qc.columns]

    qc_xlsx = opj(OUT_DIR, 'qc_report.xlsx')
    with pd.ExcelWriter(qc_xlsx, engine='openpyxl') as writer:
        if func_cols_present:
            (df_qc[id_cols + func_cols_present]
             .dropna(subset=func_cols_present, how='all')
             .to_excel(writer, sheet_name='func_qc', index=False))
        if anat_cols_present:
            (df_qc[id_cols + anat_cols_present]
             .dropna(subset=anat_cols_present, how='all')
             .to_excel(writer, sheet_name='anat_qc', index=False))
        if net_cols_present:
            (df_qc[id_cols + net_cols_present]
             .dropna(subset=net_cols_present, how='all')
             .to_excel(writer, sheet_name='network_qc', index=False))
    print(f"  QC report → qc_report.xlsx  "
          f"(func:{len(func_cols_present)}  anat:{len(anat_cols_present)}  "
          f"net:{len(net_cols_present)} cols)")

# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def _bids_colors(all_bids_pairs):
    unique = list(dict.fromkeys(bd for _, bd in all_bids_pairs))
    return {bd: PALETTE[i % len(PALETTE)] for i, bd in enumerate(unique)}


def _total_brain_volume(df_vol):
    keys   = [c for c in ('subject', 'session', 'bids_dir') if c in df_vol.columns]
    totals = df_vol.groupby(keys)['volume_mm3'].transform('sum')
    return totals


def _filter_atlas_level(df, atlas_level):
    """Filter by atlas_level if column present, otherwise return as-is."""
    if df is None or df.empty:
        return df
    if "atlas_level" in df.columns:
        return df[df["atlas_level"] == atlas_level]
    return df


def plot_cross_species_dots(
    df, metric_col, metric_label, regions, output_path,
    hemisphere='left', log_scale=False,
    normalise_by_brain=False, df_vol_for_norm=None,
    atlas_level=1, bids_col=None, figsize=None,
):
    """
    Cross-species comparison: one panel per region.
    x-axis = species, dots coloured by BIDS directory.
    bids_col : global colour dict {bids_dir: colour} — if None, built locally.
    """
    import matplotlib.ticker as ticker
    with plt.rc_context(PAPER_RC):
        plot_df = _filter_atlas_level(df.copy(), atlas_level)
        if hemisphere and hemisphere != 'bilateral':
            filtered = plot_df[plot_df['hemisphere'] == hemisphere]
            # If no left/right data found, fall back to bilateral
            plot_df = filtered if not filtered.empty else \
                      plot_df[plot_df['hemisphere'] == 'bilateral']
        elif hemisphere == 'bilateral':
            plot_df = plot_df[plot_df['hemisphere'] == 'bilateral']
        # if hemisphere is None: keep all
        plot_df = plot_df[plot_df['region'].isin(regions)]
        if plot_df.empty:
            warnings.warn(f"plot_cross_species_dots: no data for {metric_col} "
                          f"hemi={hemisphere}"); return None

        if normalise_by_brain and df_vol_for_norm is not None:
            keys = [c for c in ('subject', 'session', 'bids_dir')
                    if c in plot_df.columns and c in df_vol_for_norm.columns]
            brain_totals = (
                df_vol_for_norm.groupby(keys)['volume_mm3']
                .sum().reset_index().rename(columns={'volume_mm3': '_brain_total'})
            )
            plot_df = plot_df.merge(brain_totals, on=keys, how='left')
            plot_df[metric_col] = plot_df[metric_col] / plot_df['_brain_total']
            metric_label = metric_label + ' / brain vol.'

        species_order = _phylo_sort(plot_df['species'].unique())
        bids_dirs     = sorted(plot_df['bids_dir'].unique())
        # Use global mapping if provided, fill missing keys locally
        _bids_col = dict(bids_col) if bids_col else {}
        for i, bd in enumerate([b for b in bids_dirs if b not in _bids_col]):
            _bids_col[bd] = PALETTE[(len(_bids_col) + i) % len(PALETTE)]

        n_reg = len(regions)
        w, h  = figsize or (max(7, n_reg * 3.2), 5.5)

        fig, axes = plt.subplots(1, n_reg, figsize=(w, h), sharey=False)
        if n_reg == 1:
            axes = [axes]

        for ax, region in zip(axes, regions):
            rdf = plot_df[plot_df['region'] == region]
            if rdf.empty:
                ax.set_visible(False); continue

            for xi, sp in enumerate(species_order):
                sp_df = rdf[rdf['species'] == sp]
                vals  = sp_df[metric_col].dropna().values
                if len(vals) == 0:
                    continue
                _violin_strip_quartiles(ax, vals, xi, '#888888')
                for bd in sp_df['bids_dir'].unique():
                    bd_vals = sp_df[sp_df['bids_dir'] == bd][metric_col].dropna().values
                    if len(bd_vals) == 0:
                        continue
                    jitter = rng.uniform(-0.18, 0.18, len(bd_vals))
                    ax.scatter(xi + jitter, bd_vals,
                               color=_bids_col[bd], s=20, zorder=5,
                               alpha=0.80, linewidths=0.5, edgecolors='k')

            if log_scale:
                ax.set_yscale('log')
            else:
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

            # KW and n= removed per user request

            ax.set_xticks(range(len(species_order)))
            ax.set_xticklabels(species_order, rotation=35, ha='right')
            ax.set_title(region, fontweight='bold')
            ax.set_ylabel(metric_label if ax is axes[0] else '')
            ax.set_xlim(-0.65, len(species_order) - 0.35)

        handles = [mpatches.Patch(facecolor=_bids_col[bd], label=_bids_label(bd))
                   for bd in bids_dirs]
        hemi_str = f' [{hemisphere}]' if hemisphere else ''
        norm_str = ' (norm.)'   if normalise_by_brain else ''
        log_str  = ' [log]'     if log_scale          else ''
        # No figure title — panels labeled A/B/C
        # Panel labels A/B/C/D
        for _pi, _ax in enumerate(axes):
            _ax.text(-0.08, 1.02, chr(65+_pi), transform=_ax.transAxes,
                     fontsize=13, fontweight='bold', va='bottom', ha='left')
        plt.tight_layout()
        n_leg = len(bids_dirs)
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, -0.08),
                   ncol=min(n_leg, 6),
                   title='BIDS dir', fontsize=8,
                   frameon=False)
        fig.subplots_adjust(bottom=0.16)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path


def plot_cross_species_qc(qc_df, output_path, anat_metrics=None, func_metrics=None,
                           shared_nmi=None, clip_cols=None, figsize=None):
    """
    Cross-species QC figure.
    Layout:
      Row 0        : NMI side-by-side (func left | anat right)
      Rows 1..N    : func-only panels (left) | anat-only panels (right)
    clip_cols : dict {col_name: percentile} — caps y-axis at that percentile.
    """
    import matplotlib.ticker as ticker
    with plt.rc_context(PAPER_RC):
        shared_nmi  = [m for m in (shared_nmi  or SHARED_NMI_COLS) if m in qc_df.columns]
        anat_metrics= [m for m in (anat_metrics or ANAT_QC_COLS)   if m in qc_df.columns]
        func_metrics= [m for m in (func_metrics or FUNC_QC_COLS)   if m in qc_df.columns]
        clip_cols   = clip_cols or QC_CLIP_COLS

        if not anat_metrics and not func_metrics and not shared_nmi:
            warnings.warn('plot_cross_species_qc: no QC metrics found'); return None

        species_order = _phylo_sort(qc_df['species'].unique())
        bids_dirs     = sorted(qc_df['bids_dir'].unique())                         if 'bids_dir' in qc_df.columns else []
        # Build bids_col keyed by basename — normalize all entries
        all_bd_in_data = sorted(qc_df['bids_dir'].unique())                          if 'bids_dir' in qc_df.columns else []
        bids_col = {os.path.basename(str(k).rstrip('/')): v
                    for k, v in GLOBAL_BIDS_COL.items()}
        for i, bd in enumerate(all_bd_in_data):
            key = os.path.basename(str(bd).rstrip('/'))
            if key not in bids_col:
                bids_col[key] = PALETTE[(len(bids_col) + i) % len(PALETTE)]

        # Layout: 1 NMI row + max(func, anat) rows; 2 columns
        n_shared = 1 if shared_nmi else 0
        n_rows   = n_shared + max(len(func_metrics), len(anat_metrics), 1)
        n_cols   = 2

        w, h = figsize or (n_cols * 5.5, n_rows * 3.0)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(w, h), squeeze=False)

        def _draw_panel(ax, metric):
            """Draw one violin+dot panel, clip y-axis if requested."""
            mdf = qc_df[['species', 'bids_dir', metric]].dropna()                   if 'bids_dir' in qc_df.columns                   else qc_df[['species', metric]].dropna()

            for xi, sp in enumerate(species_order):
                sp_df = mdf[mdf['species'] == sp]
                vals  = sp_df[metric].values
                if len(vals) == 0:
                    continue
                _violin_strip_quartiles(ax, vals, xi, '#888888')
                if 'bids_dir' in sp_df.columns:
                    for bd in sp_df['bids_dir'].unique():
                        bd_v = sp_df[sp_df['bids_dir'] == bd][metric].values
                        if len(bd_v) == 0: continue
                        jitter = rng.uniform(-0.18, 0.18, len(bd_v))
                        # Normalize to basename for color lookup
                        bd_key = os.path.basename(str(bd).rstrip('/'))
                        # Try basename first, then full key, then species-based fallback
                        if bd_key in bids_col:
                            dot_color = bids_col[bd_key]
                        elif bd in bids_col:
                            dot_color = bids_col[bd]
                        else:
                            # Fallback: use species index in PALETTE
                            dot_color = PALETTE[
                                list(species_order).index(sp)
                                % len(PALETTE)] if sp in species_order                                 else PALETTE[xi % len(PALETTE)]
                        ax.scatter(xi + jitter, bd_v,
                                   color=dot_color, s=18, zorder=5,
                                   alpha=0.80, linewidths=0.5, edgecolors='k')

            ax.set_xticks(range(len(species_order)))
            ax.set_xticklabels(species_order, rotation=35, ha='right', fontsize=8)
            clean = metric.replace('func_','').replace('anat_','').replace('_',' ')
            ax.set_title(clean, fontweight='bold', fontsize=9)
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

            # Clip y-axis at given percentile if metric is noisy
            if metric in clip_cols:
                try:
                    all_vals = pd.to_numeric(mdf[metric], errors='coerce').dropna().values
                    if all_vals.size > 0:
                        pct  = clip_cols[metric]
                        ymax = float(np.nanpercentile(all_vals, pct))
                        ymin = float(np.nanpercentile(all_vals, max(0, 100 - pct)))
                        ax.set_ylim(bottom=min(0, ymin) * 0.95, top=ymax * 1.05)
                except Exception as _e:
                    warnings.warn(f'clip failed for {metric}: {_e}')

        # Row 0: shared NMI (func left, anat right)
        if shared_nmi:
            nmi_func = next((m for m in shared_nmi if 'func' in m), None)
            nmi_anat = next((m for m in shared_nmi if 'anat' in m), None)
            if nmi_func:
                _draw_panel(axes[0][0], nmi_func)
                axes[0][0].set_title('NMI (func)', fontweight='bold',
                                     fontsize=10, color='black')
            else:
                axes[0][0].set_visible(False)
            if nmi_anat:
                _draw_panel(axes[0][1], nmi_anat)
                axes[0][1].set_title('NMI (anat)', fontweight='bold',
                                     fontsize=10, color='black')
            else:
                axes[0][1].set_visible(False)

        # Column headers via fig.text (placed above first row of panels)
        # These are drawn after tight_layout so we use supxlabel-style positioning
        _col_headers = [('Functional QC', '#0072B2'), ('Anatomical QC', '#E69F00')]

        # Remaining rows: func (col 0) and anat (col 1)
        for row_i, metric in enumerate(func_metrics):
            row = n_shared + row_i
            _draw_panel(axes[row][0], metric)

        for row_i, metric in enumerate(anat_metrics):
            row = n_shared + row_i
            _draw_panel(axes[row][1], metric)

        # Hide empty cells
        for row_i in range(len(func_metrics), n_rows - n_shared):
            axes[n_shared + row_i][0].set_visible(False)
        for row_i in range(len(anat_metrics), n_rows - n_shared):
            axes[n_shared + row_i][1].set_visible(False)

        # No suptitle — column headers serve as titles
        plt.tight_layout()
        # Add column headers after tight_layout (positions in figure fraction)
        for col_i, (label, _unused_color) in enumerate(_col_headers):
            ax_ref = axes[0][col_i]
            bbox   = ax_ref.get_position()
            x_mid  = (bbox.x0 + bbox.x1) / 2
            fig.text(x_mid, 0.995, label, ha='center', va='top',
                     fontsize=12, fontweight='bold', color='black',
                     transform=fig.transFigure)
        if bids_dirs:
            handles = [mpatches.Patch(
                            facecolor=bids_col.get(
                                os.path.basename(str(bd).rstrip('/')),
                                PALETTE[i % len(PALETTE)]),
                            label=_bids_label(bd))
                       for i, bd in enumerate(bids_dirs)]
            fig.legend(handles=handles, loc='lower center',
                       bbox_to_anchor=(0.5, -0.02),
                       ncol=min(len(handles), 6),
                       title='BIDS dir', fontsize=8, frameon=False)
            fig.subplots_adjust(bottom=0.08, top=0.92)

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path


def plot_qc_fc_report(qc_df, output_path, fc_cols=None, figsize=None):
    """
    Network / FC quality report: spreading of classification metrics per species.
    Shows intra/inter hemisphere correlations, pct pos, specificity scores.
    """
    import matplotlib.ticker as ticker
    with plt.rc_context(PAPER_RC):
        fc_cols = [c for c in (fc_cols or FUNC_NET_COLS) if c in qc_df.columns]
        if not fc_cols:
            warnings.warn('plot_qc_fc_report: no network columns found'); return None

        species_order = _phylo_sort(qc_df['species'].unique())
        bids_dirs     = sorted(qc_df['bids_dir'].unique())                         if 'bids_dir' in qc_df.columns else []
        _bids_col = dict(GLOBAL_BIDS_COL)
        for i, bd in enumerate([b for b in bids_dirs if b not in _bids_col]):
            _bids_col[bd] = PALETTE[(len(_bids_col) + i) % len(PALETTE)]

        label_map = {
            'net_intra_left_mean':         'Intra L (mean r)',
            'net_intra_right_mean':        'Intra R (mean r)',
            'net_inter_mean':              'Inter hemi (mean r)',
            'net_intra_left_pct_pos':      'Intra L (% pos)',
            'net_intra_right_pct_pos':     'Intra R (% pos)',
            'net_inter_pct_pos':           'Inter (% pos)',
            'net_p_intra_vs_inter':        'p intra vs inter',
            'sp_specific_correlation':     'Specific pair r',
            'sp_nonspecific_correlation':  'Non-specific pair r',
            'sp_specific':                 'Specificity (yes=1)',
        }

        NCOLS  = min(4, len(fc_cols))
        n_rows = math.ceil(len(fc_cols) / NCOLS)
        w, h   = figsize or (NCOLS * 3.8, n_rows * 3.5)
        fig, axes = plt.subplots(n_rows, NCOLS, figsize=(w, h))
        axes = np.array(axes).flatten()

        for i, metric in enumerate(fc_cols):
            ax  = axes[i]
            mdf = qc_df[['species', 'bids_dir', metric]].dropna()                   if 'bids_dir' in qc_df.columns                   else qc_df[['species', metric]].dropna()

            for xi, sp in enumerate(species_order):
                sp_df = mdf[mdf['species'] == sp]
                vals  = sp_df[metric].values
                if len(vals) == 0:
                    continue
                _violin_strip_quartiles(ax, vals, xi, '#888888')
                if 'bids_dir' in sp_df.columns:
                    for bd in sp_df['bids_dir'].unique():
                        bd_vals = sp_df[sp_df['bids_dir'] == bd][metric].values
                        if len(bd_vals) == 0:
                            continue
                        jitter = rng.uniform(-0.18, 0.18, len(bd_vals))
                        ax.scatter(xi + jitter, bd_vals,
                                   color=_bids_col.get(bd, '#888888'), s=18, zorder=5,
                                   alpha=0.80, linewidths=0.5, edgecolors='k')

            ax.set_xticks(range(len(species_order)))
            ax.set_xticklabels(species_order, rotation=35, ha='right', fontsize=8)
            ax.set_title(label_map.get(metric, metric.replace('_', ' ')),
                         fontweight='bold', fontsize=9)
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        # No suptitle for network QC
        plt.tight_layout()
        if bids_dirs:
            handles = [mpatches.Patch(facecolor=_bids_col.get(bd, '#888'), label=_bids_label(bd))
                       for bd in bids_dirs]
            fig.legend(handles=handles, loc='lower center',
                       bbox_to_anchor=(0.5, -0.05),
                       ncol=min(len(handles), 6),
                       title='BIDS dir', fontsize=8, frameon=False)
            fig.subplots_adjust(bottom=0.12)

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path

def plot_combo(df_surface, df_volume, df_thickness, output_path, regions,
               hemisphere='left', log_scale=False, normalise_by_brain=False,
               atlas_level=1, bids_col=None):
    """
    Combined figure: surface / volume / thickness — 3 rows × N regions.
    For each modality the hemisphere filter falls back to 'bilateral' if the
    requested hemisphere is absent (e.g. bilateral atlas for volumes).
    bids_col : global colour dict {bids_dir: colour}.
    """
    import matplotlib.ticker as ticker

    modalities = [
        (df_surface,   'surface_area_mm2', 'Surface area (mm²)'),
        (df_volume,    'volume_mm3',        'Volume (mm³)'),
        (df_thickness, 'thickness_mm',      'Cortical thickness (mm)'),
    ]
    modalities = [(df, col, lbl) for df, col, lbl in modalities
                  if df is not None and not df.empty and col in df.columns]
    if not modalities:
        warnings.warn('plot_combo: no data'); return None

    n_mod = len(modalities)
    n_reg = len(regions)

    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(n_mod, n_reg,
                                 figsize=(max(8, n_reg * 3.0), n_mod * 3.5),
                                 squeeze=False)
        all_bids_dirs = []

        for row, (df, metric_col, metric_label) in enumerate(modalities):
            plot_df = _filter_atlas_level(df.copy(), atlas_level)

            # Hemisphere filter with bilateral fallback
            if hemisphere and hemisphere != 'bilateral':
                filtered = plot_df[plot_df['hemisphere'] == hemisphere]
                plot_df  = filtered if not filtered.empty else \
                           plot_df[plot_df['hemisphere'] == 'bilateral']
            elif hemisphere == 'bilateral':
                bil = plot_df[plot_df['hemisphere'] == 'bilateral']
                # If bilateral not found, sum L+R per subject/session/region
                if bil.empty and 'hemisphere' in plot_df.columns:
                    grp_cols = [c for c in ('species','bids_dir','subject','session',
                                            'region','atlas_level') if c in plot_df.columns]
                    plot_df = (plot_df.groupby(grp_cols)[metric_col]
                               .sum().reset_index())
                else:
                    plot_df = bil

            plot_df = plot_df[plot_df['region'].isin(regions)]

            if normalise_by_brain and df_volume is not None:
                keys = [c for c in ('subject', 'session', 'bids_dir')
                        if c in plot_df.columns and c in df_volume.columns]
                brain_totals = (
                    df_volume.groupby(keys)['volume_mm3']
                    .sum().reset_index().rename(columns={'volume_mm3': '_brain_total'})
                )
                plot_df = plot_df.merge(brain_totals, on=keys, how='left')
                if '_brain_total' in plot_df.columns:
                    plot_df[metric_col] = plot_df[metric_col] / plot_df['_brain_total']
                    metric_label += ' / brain vol.'

            species_order = _phylo_sort(plot_df['species'].unique())
            bids_dirs     = sorted(plot_df['bids_dir'].unique()) \
                            if 'bids_dir' in plot_df.columns else []
            all_bids_dirs = sorted(set(all_bids_dirs + bids_dirs))
            _bids_col = dict(bids_col) if bids_col else {}
            for i, bd in enumerate([b for b in bids_dirs if b not in _bids_col]):
                _bids_col[bd] = PALETTE[(len(_bids_col) + i) % len(PALETTE)]

            for col, region in enumerate(regions):
                ax  = axes[row][col]
                rdf = plot_df[plot_df['region'] == region]

                for xi, sp in enumerate(species_order):
                    sp_df = rdf[rdf['species'] == sp]
                    vals  = sp_df[metric_col].dropna().values
                    if len(vals) == 0:
                        continue
                    _violin_strip_quartiles(ax, vals, xi, '#888888')
                    if 'bids_dir' in sp_df.columns:
                        for bd in sp_df['bids_dir'].unique():
                            bd_v = sp_df[sp_df['bids_dir'] == bd][metric_col].dropna().values
                            if len(bd_v) == 0:
                                continue
                            jitter = rng.uniform(-0.18, 0.18, len(bd_v))
                            ax.scatter(xi + jitter, bd_v,
                                       color=_bids_col[bd], s=16, zorder=5,
                                       alpha=0.8, linewidths=0.4, edgecolors='k')

                if log_scale:
                    ax.set_yscale('log')
                else:
                    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

                ax.set_xticks(range(len(species_order)))
                if row == n_mod - 1:
                    ax.set_xticklabels(species_order, rotation=35, ha='right', fontsize=8)
                else:
                    ax.set_xticklabels([], fontsize=0)
                ax.set_xlim(-0.65, len(species_order) - 0.35)

                if col == 0:
                    ax.set_ylabel(metric_label, fontsize=9)
                if row == 0:
                    ax.set_title(region, fontweight='bold', fontsize=9)

        _bids_col_final = dict(bids_col) if bids_col else {}
        handles = [mpatches.Patch(facecolor=_bids_col_final.get(bd, '#888'),
                                  label=_bids_label(bd))
                   for bd in all_bids_dirs if bd in _bids_col_final]
        hemi_str = f' [{hemisphere}]' if hemisphere else ''
        norm_str = ' (normalised)' if normalise_by_brain else ''
        log_str  = ' [log scale]'  if log_scale          else ''
        # No figure title
        plt.tight_layout()
        if handles:
            fig.legend(handles=handles, loc='lower center',
                       bbox_to_anchor=(0.5, -0.04),
                       ncol=min(len(handles), 6),
                       title='BIDS dir', fontsize=8, frameon=False)
            fig.subplots_adjust(bottom=0.12)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path


def plot_specificity_bar(qc_df, output_path, cat_col='sp_category', figsize=None):
    """
    Bar chart showing the proportion of each specificity category per species.
    Categories: Specific / Unspecific / Spurious / No
    """
    with plt.rc_context(PAPER_RC):
        if cat_col not in qc_df.columns:
            warnings.warn(f'plot_specificity_bar: column {cat_col} not found'); return None

        species_order = _phylo_sort(qc_df['species'].unique())
        cats    = ['Specific', 'Unspecific', 'Spurious', 'No']
        colors  = {'Specific':   '#009E73',
                   'Unspecific': '#E69F00',
                   'Spurious':   '#CC79A7',
                   'No':         '#D55E00'}

        # Count per species
        counts = {}
        for sp in species_order:
            sp_df = qc_df[qc_df['species'] == sp][cat_col].dropna()
            total = len(sp_df)
            if total == 0:
                counts[sp] = {c: 0.0 for c in cats}
            else:
                counts[sp] = {c: 100.0 * (sp_df == c).sum() / total for c in cats}

        n_sp = len(species_order)
        w, h = figsize or (max(5, n_sp * 1.2), 4.5)
        fig, ax = plt.subplots(figsize=(w, h))

        bar_w   = 0.7
        bottoms = np.zeros(n_sp)
        x       = np.arange(n_sp)
        for cat in cats:
            vals = np.array([counts[sp][cat] for sp in species_order])
            ax.bar(x, vals, bar_w, bottom=bottoms,
                   color=colors[cat], label=cat, edgecolor='white', linewidth=0.5)
            # Add % label if > 5%
            for xi, (v, b) in enumerate(zip(vals, bottoms)):
                if v > 5:
                    ax.text(xi, b + v / 2, f'{v:.0f}%',
                            ha='center', va='center', fontsize=7, color='white',
                            fontweight='bold')
            bottoms += vals

        ax.set_xticks(x)
        ax.set_xticklabels(species_order, rotation=35, ha='right', fontsize=9)
        ax.set_ylabel('% subjects', fontsize=10)
        ax.set_ylim(0, 105)
        ax.yaxis.set_major_formatter(
            __import__('matplotlib.ticker', fromlist=['FormatStrFormatter'])
            .FormatStrFormatter('%.0f%%'))
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.22),
                  ncol=len(cats), fontsize=8, frameon=False, title='Category')
        # No figure title
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path


def plot_corr_matrix_combo(
    species_config,
    output_path,
    atlas_name="EDNIxCSC",
    atlas_level=3,
    use_lr=False,
    figsize=None,
    vmin=-0.5, vmax=0.5,
    pval_thresh=0.05,
):
    """
    Combined correlation matrix figure.
    3 rows (mean r, p<thresh mask, variance) × N species columns.
    Only ROIs common to all species are shown.
    Species ordered phylogenetically.
    """
    import matplotlib.ticker as ticker
    from matplotlib.colors import TwoSlopeNorm

    csv_dir   = os.path.join(os.path.dirname(output_path), 'corr_matrix_csv')
    corr_data = collect_corr_matrices(
        species_config, atlas_name, atlas_level, use_lr,
        save_csv_dir=csv_dir)
    if not corr_data:
        warnings.warn('plot_corr_matrix_combo: no matrix data found'); return None

    species_order = [s for s in _phylo_sort(list(corr_data.keys()))
                     if s in corr_data]
    n_sp = len(species_order)

    import matplotlib.cm as _cm
    from matplotlib.colors import Normalize

    norm_r   = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    cmap_r   = 'RdBu_r'
    cmap_var = 'viridis'

    # Compute global variance range for consistent colorbar
    all_var_vals = np.concatenate([corr_data[sp]['var'].ravel()
                                   for sp in species_order])
    valid_var = all_var_vals[~np.isnan(all_var_vals)]
    vmax_var  = float(np.nanpercentile(valid_var, 95)) if valid_var.size > 0 else 0.01
    norm_var  = Normalize(0, vmax_var)

    n_rows = 3
    # Reserve 0.12 of width on the right for colorbars
    w, h   = figsize or (max(8, n_sp * 3.5), n_rows * 3.8)
    fig    = plt.figure(figsize=(w, h))
    # Grid: n_rows rows × (n_sp main + 1 cbar) columns
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(n_rows, n_sp + 1,
                           width_ratios=[1.0] * n_sp + [0.05],
                           wspace=0.04, hspace=0.15)

    for col, sp in enumerate(species_order):
        d     = corr_data[sp]
        mean  = d['mean']
        pval  = d['pval']
        var   = d['var']
        n     = d['n']
        rois  = d['rois']
        short = [r.replace('L_','L ').replace('R_','R ')[:16] for r in rois]
        n_roi = len(rois)

        # Row 0: mean correlation
        ax0 = fig.add_subplot(gs[0, col])
        ax0.imshow(mean, cmap=cmap_r, norm=norm_r, aspect='auto',
                   interpolation='nearest')
        ax0.set_title(f'{sp} (n={n})', fontweight='bold', fontsize=8)
        if col == 0:
            ax0.set_ylabel('Mean r', fontsize=8)
            ax0.set_yticks(range(n_roi))
            ax0.set_yticklabels(short, fontsize=4)
        else:
            ax0.set_yticks([])
        ax0.set_xticks([])

        # Row 1: p-value heatmap (log scale, capped at 1)
        ax1 = fig.add_subplot(gs[1, col])
        pval_disp = np.clip(pval, 1e-10, 1.0)   # avoid log(0)
        ax1.imshow(-np.log10(pval_disp), cmap='hot_r', aspect='auto',
                   vmin=0, vmax=4, interpolation='nearest')
        if col == 0:
            ax1.set_ylabel('-log10(p)', fontsize=8)
            ax1.set_yticks(range(n_roi))
            ax1.set_yticklabels(short, fontsize=4)
        else:
            ax1.set_yticks([])
        ax1.set_xticks([])

        # Row 2: variance
        ax2 = fig.add_subplot(gs[2, col])
        ax2.imshow(var, cmap=cmap_var, norm=norm_var, aspect='auto',
                   interpolation='nearest')
        if col == 0:
            ax2.set_ylabel('Variance', fontsize=8)
            ax2.set_yticks(range(n_roi))
            ax2.set_yticklabels(short, fontsize=4)
        else:
            ax2.set_yticks([])
        ax2.set_xticks([])

    # Colorbars in the last column
    # Row 0: mean r
    cbar_ax_r = fig.add_subplot(gs[0, -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm_r, cmap=cmap_r),
                 cax=cbar_ax_r, label='Pearson r')
    cbar_ax_r.yaxis.set_label_position('right')

    # Row 1: -log10(p)
    norm_p = Normalize(vmin=0, vmax=4)
    cbar_ax_p = fig.add_subplot(gs[1, -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm_p, cmap='hot_r'),
                 cax=cbar_ax_p, label='-log10(p)')
    cbar_ax_p.yaxis.set_label_position('right')

    # Row 2: variance
    cbar_ax_v = fig.add_subplot(gs[2, -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm_var, cmap=cmap_var),
                 cax=cbar_ax_v, label='Variance')
    cbar_ax_v.yaxis.set_label_position('right')

    lr_str = 'LR' if use_lr else 'bilateral'
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches='tight', dpi=200)
    plt.close(fig)
    print(f'  [plot] {output_path}')
    return output_path


def plot_corr_matrix_sig_norm(
        species_config,
        output_path,
        atlas_name="EDNIxCSC",
        atlas_level=2,
        use_lr=True,
        pval_thresh=0.05,
        figsize=None,
):
    """
    Normalized significant FC matrix — one panel per species.
    Shows only cells with p < pval_thresh; values normalized per species
    (z-score of significant cells) so cross-species comparison is meaningful.
    """
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.cm as _cm

    csv_dir = os.path.join(os.path.dirname(output_path), 'corr_matrix_csv')
    corr_data = collect_corr_matrices(
        species_config, atlas_name, atlas_level, use_lr,
        save_csv_dir=csv_dir)
    if not corr_data:
        warnings.warn('plot_corr_matrix_sig_norm: no data');
        return None

    species_order = [s for s in _phylo_sort(list(corr_data.keys()))
                     if s in corr_data]
    n_sp = len(species_order)

    w, h = figsize or (max(6, n_sp * 3.5), 5.5)  # taller for x-axis labels
    fig = plt.figure(figsize=(w, h))
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(1, n_sp + 1,
                           width_ratios=[1.0] * n_sp + [0.05],
                           wspace=0.04)

    # [-1, 1] normalization: each species normalized independently by its own max abs value
    # so the scale is comparable across species without z-score distortion
    norm = TwoSlopeNorm(vmin=-3.0, vcenter=0, vmax=3.0)  # r / mean(|r|) scale
    cmap = 'RdBu_r'

    for col, sp in enumerate(species_order):
        d = corr_data[sp]
        mean = d['mean']
        pval = d['pval']
        rois = d['rois']
        n = d['n']
        # Short labels — strip hemisphere prefix, keep max 16 chars
        short = [r.replace('L_', 'L ').replace('R_', 'R ')[:16] for r in rois]
        n_roi = len(rois)

        # Normalize: each cell / mean(abs(r)) of the species matrix
        # Values > 1 = stronger than average FC, < 1 = weaker, negative = anti-correlated
        mat_mean_r = float(np.nanmean(np.abs(mean)))
        sig_norm = mean / mat_mean_r if mat_mean_r > 1e-6 else mean.copy()

        ax = fig.add_subplot(gs[0, col])
        ax.imshow(sig_norm, cmap=cmap, norm=norm, aspect='equal',
                  interpolation='nearest')
        ax.set_title(f'{sp} (n={n})', fontweight='bold', fontsize=8)
        # y-axis: ROI labels on first panel only
        if col == 0:
            ax.set_yticks(range(n_roi))
            ax.set_yticklabels(short, fontsize=4)
            ax.set_ylabel('r / mean(|r|)', fontsize=8)
        else:
            ax.set_yticks([])
        # No x-axis ROI labels on sig_norm (too crowded, y-axis already shows them)
        ax.set_xticks([])

    # Colorbar
    import matplotlib.cm as _cm
    from matplotlib.colors import Normalize
    cbar_ax = fig.add_subplot(gs[0, -1])
    fig.colorbar(_cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=cbar_ax, label='r / mean(|r|)')
    cbar_ax.yaxis.set_label_position('right')

    lr_str = 'LR' if use_lr else 'bilateral'
    # No figure title
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches='tight', dpi=200)
    plt.close(fig)
    print(f'  [plot] {output_path}')
    return output_path


def _plot_fc_scaling(df, y_col, y_label, output_path, x_cols=None, x_labels=None):
    """
    Standalone FC allometric scaling plot.
    Works directly on dataframe columns (not region-based).
    x_cols / x_labels : list of column names and axis labels.
                        None entries are skipped.
    Adds 95% prediction interval and A/B/C panel labels.
    """
    import matplotlib.ticker as ticker
    from scipy import stats as _stats

    x_cols   = [c for c in (x_cols   or []) if c is not None]
    x_labels = [l for c, l in zip(x_cols or [], x_labels or []) if c is not None]
    # Filter to columns that exist and have data
    pairs = [(c, l) for c, l in zip(x_cols, x_labels)
             if c in df.columns and df[c].notna().any()]
    if not pairs:
        warnings.warn(f'_plot_fc_scaling: no x columns with data for {y_col}'); return None

    species_order = _phylo_sort(df['species'].unique())
    sp_colors     = {sp: PALETTE[i % len(PALETTE)] for i, sp in enumerate(species_order)}
    PANEL_LABELS  = list('ABCDEFGH')

    n_panels = len(pairs)
    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, n_panels, figsize=(n_panels * 4.5, 4.5))
        if n_panels == 1:
            axes = [axes]

        for panel_i, (ax, (x_col, x_label)) in enumerate(zip(axes, pairs)):
            sub = df[[y_col, x_col, 'species']].dropna()
            if sub.empty:
                ax.set_visible(False); continue

            # Scatter per species
            for sp in species_order:
                sp_sub = sub[sub['species'] == sp]
                if sp_sub.empty: continue
                ax.scatter(sp_sub[x_col], sp_sub[y_col],
                           color=sp_colors[sp], s=40, zorder=5,
                           alpha=0.85, edgecolors='k', linewidths=0.5,
                           label=sp)

            # OLS + 95% prediction interval on log-log
            xv = np.clip(sub[x_col].values.astype(float), 1e-12, None)
            yv = sub[y_col].values.astype(float)
            # Allow negative FC values → use linear scale for FC, log for weight
            if x_col in ('body_g', 'brain_g'):
                lx = np.log10(xv)
                ly = yv          # FC r values can be negative — keep linear
                use_log_x = True
            else:
                lx = xv          # already FC values
                ly = yv
                use_log_x = False

            valid = np.isfinite(lx) & np.isfinite(ly)
            if valid.sum() >= 3:
                lx_v, ly_v = lx[valid], ly[valid]
                slope, intercept, r, p, _ = _stats.linregress(lx_v, ly_v)
                x_fit     = np.linspace(lx_v.min(), lx_v.max(), 200)
                y_fit     = intercept + slope * x_fit
                n_v       = valid.sum()
                t_crit    = _stats.t.ppf(0.975, df=n_v - 2)
                x_mean    = lx_v.mean()
                ssx       = np.sum((lx_v - x_mean)**2)
                resid_var = np.sum((ly_v - (intercept + slope*lx_v))**2) / (n_v - 2)
                se_pred   = np.sqrt(resid_var * (1 + 1/n_v +
                                    (x_fit - x_mean)**2 / (ssx if ssx > 1e-12 else 1e-12)))
                x_plot = 10**x_fit if use_log_x else x_fit
                ax.plot(x_plot, y_fit, color='k', lw=1.5, linestyle='--',
                        alpha=0.8, zorder=3)
                ax.fill_between(x_plot,
                                y_fit - t_crit * se_pred,
                                y_fit + t_crit * se_pred,
                                color='k', alpha=0.15, zorder=2)
                p_str = f'{p:.3f}' if p >= 0.001 else '<0.001'
                ax.text(0.05, 0.95, f'slope={slope:.2f}  r={r:.2f}  p={p_str}',
                        transform=ax.transAxes, fontsize=8, va='top', color='#333333')

            if use_log_x:
                ax.set_xscale('log')
                ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())

            ax.set_xlabel(x_label, fontsize=9)
            ax.set_ylabel(y_label if panel_i == 0 else '', fontsize=9)
            # Panel label A/B/C
            ax.text(-0.12, 1.02, PANEL_LABELS[panel_i], transform=ax.transAxes,
                    fontsize=13, fontweight='bold', va='bottom', ha='left')

        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                   for sp in species_order]
        plt.tight_layout()
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, -0.06),
                   ncol=min(len(handles), 7),
                   fontsize=8, frameon=False, title='Species')
        fig.subplots_adjust(bottom=0.18)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path


def plot_brain_scaling(
    df_morph,
    metric_col,
    metric_label,
    output_path,
    y_region=None,
    x_region=None,
    atlas_level=1,
    hemisphere='bilateral',
    figsize=None,
):
    """
    Allometric brain scaling — 3-panel figure.
    For each subject/session compute the mean of y_region metric, then plot:
      Panel 1 : y vs body weight  (log-log)
      Panel 2 : y vs brain weight (log-log)
      Panel 3 : y vs x_region metric (log-log, internal scaling)
    Each dot = one subject, coloured by species (phylogenetic order).
    Regression line (OLS on log-log) per panel.
    """
    import matplotlib.ticker as ticker
    from scipy import stats as _stats

    y_region = y_region or SCALE_Y_REGION
    x_region = x_region or SCALE_X_REGION

    with plt.rc_context(PAPER_RC):
        plot_df = _filter_atlas_level(df_morph.copy(), atlas_level)

        # hemisphere filter
        if hemisphere and hemisphere != 'bilateral':
            filt = plot_df[plot_df['hemisphere'] == hemisphere]
            plot_df = filt if not filt.empty else plot_df[plot_df['hemisphere'] == 'bilateral']
        elif hemisphere == 'bilateral':
            bil = plot_df[plot_df['hemisphere'] == 'bilateral']
            if bil.empty:
                # sum L+R
                grp = [c for c in ('species','bids_dir','subject','session',
                                   'atlas_level') if c in plot_df.columns]
                plot_df = (plot_df[plot_df['region'].isin([y_region, x_region])]
                           .groupby(grp + ['region'])[metric_col]
                           .sum().reset_index())
            else:
                plot_df = bil

        # Extract per-subject means for y_region and x_region
        id_cols = [c for c in ('species','bids_dir','subject','session')
                   if c in plot_df.columns]

        def _region_mean(reg):
            sub = plot_df[plot_df['region'] == reg]
            if sub.empty:
                return pd.DataFrame()
            return sub.groupby(id_cols)[metric_col].mean().reset_index()                      .rename(columns={metric_col: reg})

        df_y = _region_mean(y_region)
        df_x = _region_mean(x_region)

        if df_y.empty:
            warnings.warn(f'plot_brain_scaling: no data for y_region={y_region}')
            return None

        df_merged = df_y.copy()
        if not df_x.empty:
            df_merged = df_merged.merge(df_x, on=id_cols, how='inner')

        # Add weight columns
        df_merged['body_g']  = df_merged['species'].map(SPECIES_BODY_WEIGHT_G)
        df_merged['brain_g'] = df_merged['species'].map(SPECIES_BRAIN_WEIGHT_G)

        species_order = _phylo_sort(df_merged['species'].unique())
        sp_colors     = {sp: PALETTE[i % len(PALETTE)]
                         for i, sp in enumerate(species_order)}

        panels = [
            ('body_g',  'Body weight (g)',  'body'),
            ('brain_g', 'Brain weight (g)', 'brain'),
        ]
        if x_region and x_region in df_merged.columns:
            panels.append((x_region, f'{x_region[:20]}', 'internal'))
        # Only keep panels with data
        panels = [(xc, xl, xt) for xc, xl, xt in panels
                  if xc in df_merged.columns and df_merged[xc].notna().any()]

        n_panels = len(panels)
        w, h = figsize or (n_panels * 4.5, 4.5)
        fig, axes = plt.subplots(1, n_panels, figsize=(w, h))
        if n_panels == 1:
            axes = [axes]

        for ax, (x_col, x_label, _) in zip(axes, panels):
            sub = df_merged[[y_region, x_col, 'species']].dropna()
            if sub.empty:
                ax.set_visible(False); continue

            for sp in species_order:
                sp_sub = sub[sub['species'] == sp]
                if sp_sub.empty: continue
                ax.scatter(sp_sub[x_col], sp_sub[y_region],
                           color=sp_colors[sp], s=40, zorder=5,
                           alpha=0.85, edgecolors='k', linewidths=0.5,
                           label=sp)

            # OLS regression on log-log + 95% confidence interval
            log_x = np.log10(np.clip(sub[x_col].values.astype(float), 1e-9, None))
            log_y = np.log10(np.clip(sub[y_region].values.astype(float), 1e-9, None))
            valid = np.isfinite(log_x) & np.isfinite(log_y)
            if valid.sum() >= 3:
                lx, ly = log_x[valid], log_y[valid]
                slope, intercept, r, p, se = _stats.linregress(lx, ly)
                x_fit   = np.linspace(lx.min(), lx.max(), 200)
                y_fit   = intercept + slope * x_fit
                # 95% prediction interval via t-distribution
                n_v    = valid.sum()
                t_crit = _stats.t.ppf(0.975, df=n_v - 2)
                x_mean = lx.mean()
                ssx    = np.sum((lx - x_mean)**2)
                resid_var = np.sum((ly - (intercept + slope*lx))**2) / (n_v - 2)
                if ssx > 1e-12:
                    # Prediction interval (includes individual variability)
                    se_pred = np.sqrt(resid_var * (1 + 1/n_v +
                                      (x_fit - x_mean)**2 / ssx))
                else:
                    # All x-values identical (same weight per species) —
                    # use residual std as width
                    se_pred = np.full_like(x_fit, np.sqrt(resid_var))
                ax.plot(10**x_fit, 10**y_fit,
                        color='k', lw=1.5, linestyle='--', alpha=0.8, zorder=3)
                ax.fill_between(10**x_fit,
                                10**(y_fit - t_crit * se_pred),
                                10**(y_fit + t_crit * se_pred),
                                color='k', alpha=0.15, zorder=2)
                p_str = f'{p:.3f}' if p >= 0.001 else '<0.001'
                ax.text(0.05, 0.95,
                        f'slope={slope:.2f}  r={r:.2f}  p={p_str}',
                        transform=ax.transAxes, fontsize=8,
                        va='top', color='#333333')

            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_xlabel(x_label, fontsize=9)
            ax.set_ylabel(f'{y_region} ({metric_label})' if ax is axes[0] else '',
                          fontsize=9)
            ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())

        # Panel labels A/B/C
        _pl = list('ABCDEFGH')
        for _pi, _ax in enumerate(axes):
            _ax.text(-0.12, 1.02, _pl[_pi], transform=_ax.transAxes,
                     fontsize=13, fontweight='bold', va='bottom', ha='left')
        # Shared legend
        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                   for sp in species_order]
        plt.tight_layout()
        fig.legend(handles=handles, loc='lower center',
                   bbox_to_anchor=(0.5, -0.06),
                   ncol=min(len(handles), 7),
                   fontsize=8, frameon=False, title='Species')
        fig.subplots_adjust(bottom=0.18)

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5A  —  CROSS-SPECIES FIGURES
# ─────────────────────────────────────────────────────────────────────────────

# Collect all hemispheres across modalities for the main loop
_all_hemis = sorted(set(HEMIS_SURFACE + HEMIS_THICKNESS + HEMIS_VOLUME))

for hemi in _all_hemis:
    morph_dir = opj(FIG_DIR, 'cross_species', 'morphometry_recap')
    os.makedirs(morph_dir, exist_ok=True)

    for df, metric_col, metric_label, mod, hemis in [
        (df_surface,   'surface_area_mm2', 'Surface area (mm²)',      'surface',   HEMIS_SURFACE),
        (df_volume,    'volume_mm3',        'Volume (mm³)',            'volume',    HEMIS_VOLUME),
        (df_thickness, 'thickness_mm',      'Cortical thickness (mm)', 'thickness', HEMIS_THICKNESS),
    ]:
        if df is None or df.empty:
            continue
        if hemi not in hemis:
            continue   # skip hemispheres not in this modality

        for log_s, suffix in [(False, 'raw'), (True, 'log')]:
            plot_cross_species_dots(
                df, metric_col, metric_label, PLOT_REGIONS,
                opj(morph_dir, f'{mod}_{hemi}_{suffix}.png'),
                hemisphere=hemi, log_scale=log_s, atlas_level=ATLAS_LEVEL,
                bids_col=GLOBAL_BIDS_COL,
            )

        if df_volume is not None and not df_volume.empty and metric_col != 'volume_mm3':
            plot_cross_species_dots(
                df, metric_col, metric_label, PLOT_REGIONS,
                opj(morph_dir, f'{mod}_{hemi}_norm.png'),
                hemisphere=hemi, log_scale=False,
                normalise_by_brain=True, df_vol_for_norm=df_volume,
                atlas_level=ATLAS_LEVEL, bids_col=GLOBAL_BIDS_COL,
            )

    # Combined modalities figure
    combo_dir = opj(FIG_DIR, 'cross_species', 'combo')
    os.makedirs(combo_dir, exist_ok=True)
    for log_s, norm_s, suffix in [
        (False, False, 'raw'), (True, False, 'log'), (False, True, 'norm')
    ]:
        plot_combo(df_surface, df_volume, df_thickness,
                   opj(combo_dir, f'combo_{hemi}_{suffix}.png'),
                   regions=PLOT_REGIONS, hemisphere=hemi,
                   log_scale=log_s, normalise_by_brain=norm_s,
                   atlas_level=ATLAS_LEVEL, bids_col=GLOBAL_BIDS_COL)

# QC cross-species — func QC + anat QC aligned side by side
qc_dir = opj(FIG_DIR, 'cross_species', 'qc_recap')
os.makedirs(qc_dir, exist_ok=True)
if df_qc is not None and not df_qc.empty:
    plot_cross_species_qc(
        df_qc,
        opj(qc_dir, 'QC_cross_species_anat_and_func.png'),
        anat_metrics=ANAT_QC_COLS,
        func_metrics=FUNC_QC_COLS,
    )
    plot_qc_fc_report(
        df_qc,
        opj(qc_dir, 'QC_cross_species_network_report.png'),
        fc_cols=FUNC_NET_COLS,
    )
    # Specificity category bar chart (separate figure)
    plot_specificity_bar(
        df_qc,
        opj(qc_dir, 'QC_cross_species_specificity_bar.png'),
        cat_col=FUNC_SPEC_COL,
    )

# Correlation matrix combo (mean / p<0.05 / variance per species)
plot_corr_matrix_combo(
    species_config,
    opj(combo_dir, f'combo_corr_matrix_lvl{CORR_ATLAS_LEVEL}.png'),
    atlas_name=ATLAS_NAME,
    atlas_level=CORR_ATLAS_LEVEL,
    use_lr=CORR_USE_LR,
)
# Normalized significant FC — companion figure
plot_corr_matrix_sig_norm(
    species_config,
    opj(combo_dir, f'combo_corr_matrix_sig_norm_lvl{CORR_ATLAS_LEVEL}.png'),
    atlas_name=ATLAS_NAME,
    atlas_level=CORR_ATLAS_LEVEL,
    use_lr=CORR_USE_LR,
    pval_thresh=1.0,   # show ALL cells normalized [-1,1], no significance masking
)

# ── Brain scaling / allometric plots ─────────────────────────────────────────
scaling_dir = opj(FIG_DIR, 'cross_species', 'brain_scaling')
os.makedirs(scaling_dir, exist_ok=True)

for df, metric_col, metric_label, mod in [
    (df_surface,   'surface_area_mm2', 'Surface area (mm²)',      'surface'),
    (df_volume,    'volume_mm3',        'Volume (mm³)',            'volume'),
    (df_thickness, 'thickness_mm',      'Cortical thickness (mm)', 'thickness'),
]:
    if df is None or df.empty:
        continue
    plot_brain_scaling(
        df, metric_col, metric_label,
        opj(scaling_dir, f'scaling_{mod}_lvl{ATLAS_LEVEL}.png'),
        y_region=SCALE_Y_REGION,
        x_region=SCALE_X_REGION,
        atlas_level=ATLAS_LEVEL,
    )

# ── FC allometric scaling — mPFC↔PPC pair vs body weight / brain weight / somato FC
# Runs independently — does not require _corr_lvl0 to succeed.
# ROI matching uses partial name matching (case-insensitive).

def _find_roi_exact_then_partial(exact_name, partial_pattern, rois):
    """
    ROI lookup: try exact match first (FC_Y_ROI_A/B config names),
    then fall back to case-insensitive partial match.
    Returns matched ROI name or None.
    """
    # 1. Exact match
    if exact_name in rois:
        return exact_name
    # 2. Case-insensitive exact
    exact_lo = exact_name.lower()
    for r in rois:
        if r.lower() == exact_lo:
            return r
    # 3. Partial fallback
    pat = partial_pattern.lower()
    return next((r for r in rois if pat in r.lower()), None)

print("\n  [FC scaling] extracting mPFC↔PPC ROI pair per subject...")
print(f"  Target ROIs: y_a={FC_Y_ROI_A!r}  y_b={FC_Y_ROI_B!r}")
print(f"               x_a={FC_X_ROI_A!r}  x_b={FC_X_ROI_B!r}")
_fc_pair_rows = []
for _sp, _cfg in species_config.items():
    for _bd in _cfg.get('bids_dirs', []):
        _recs = extract_corr_matrix_paths(
            _bd, ATLAS_NAME, 2, False,
            _cfg.get('list_to_keep', []), _cfg.get('list_to_remove', []))
        _done = set()
        for _rec in _recs:
            _k = (_rec['subject'], _rec['session'])
            if _k in _done: continue
            _done.add(_k)
            try:
                _rois, _mat = load_corr_matrix(_rec['path'])
                # Use exact config names first, partial pattern as fallback
                _roi_a = _find_roi_exact_then_partial(FC_Y_ROI_A, 'prefrontal',       _rois)
                _roi_b = _find_roi_exact_then_partial(FC_Y_ROI_B, 'posterior_parietal', _rois)
                _roi_sa = _find_roi_exact_then_partial(FC_X_ROI_A, 'somatosensory',    _rois)
                _roi_sb = _find_roi_exact_then_partial(FC_X_ROI_B, 'somatosensory',    _rois)
                _row = dict(species=_sp,
                            bids_dir=os.path.basename(_bd),
                            subject=_rec['subject'],
                            session=_rec['session'])
                if _roi_a and _roi_b:
                    _ia = _rois.index(_roi_a)
                    _ib = _rois.index(_roi_b)
                    _row['FC_mpfc_ppc'] = float(_mat[_ia, _ib])
                    print(f"    {_sp} sub-{_rec['subject']}: "
                          f"{_roi_a} ↔ {_roi_b} = {_row['FC_mpfc_ppc']:.3f}")
                else:
                    print(f"    [WARN] {_sp} sub-{_rec['subject']}: "
                          f"mPFC={_roi_a}  PPC={_roi_b}  "
                          f"(available: {_rois[:5]}...)")
                # Somatosensory: average L and R if both found, else whichever is found
                _somato_vals = []
                for _roi_s in [_roi_sa, _roi_sb]:
                    if _roi_s and _roi_s in _rois:
                        _is = _rois.index(_roi_s)
                        _somato_vals.append(float(np.nanmean(_mat[_is, :])))
                if _somato_vals:
                    _row['FC_somato_mean'] = float(np.mean(_somato_vals))
                _fc_pair_rows.append(_row)
            except Exception as _e:
                print(f"    [ERROR] {_rec['path']}: {_e}")

print(f"  [FC scaling] {len(_fc_pair_rows)} rows, "
      f"valid FC_mpfc_ppc: "
      f"{sum('FC_mpfc_ppc' in r for r in _fc_pair_rows)}")

if _fc_pair_rows:
    _df_fc2 = pd.DataFrame(_fc_pair_rows)
    _df_fc2.to_csv(opj(scaling_dir, 'fc_roi_pair_per_subject.csv'), index=False)

    if 'FC_mpfc_ppc' in _df_fc2.columns and _df_fc2['FC_mpfc_ppc'].notna().any():
        _df_fc2['body_g']  = _df_fc2['species'].map(SPECIES_BODY_WEIGHT_G)
        _df_fc2['brain_g'] = _df_fc2['species'].map(SPECIES_BRAIN_WEIGHT_G)

        # 3 panels: mPFC↔PPC vs (A) body weight, (B) brain weight, (C) somato FC
        _x_cols   = ['body_g', 'brain_g']
        _x_labels = ['Body weight (g)', 'Brain weight (g)']
        if 'FC_somato_mean' in _df_fc2.columns and _df_fc2['FC_somato_mean'].notna().any():
            _x_cols.append('FC_somato_mean')
            _x_labels.append('Somatosensory mean FC (r)')

        _plot_fc_scaling(
            _df_fc2,
            y_col='FC_mpfc_ppc',
            y_label='mPFC ↔ PPC  (r)',
            output_path=opj(scaling_dir, 'scaling_fc_lvl1.png'),
            x_cols=_x_cols,
            x_labels=_x_labels,
        )
    else:
        print("  [FC scaling] WARNING: no valid FC_mpfc_ppc values — plot skipped")
else:
    print("  [FC scaling] WARNING: no subjects found — check ROI names in lvl2 matrices")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5B  —  PER-BIDS FIGURES
# ─────────────────────────────────────────────────────────────────────────────

for species, bids_dir in ALL_BIDS:
    label    = _bids_label(bids_dir)
    bids_out = opj(FIG_DIR, 'per_bids', f'{species}_{label}')
    os.makedirs(bids_out, exist_ok=True)

    for df, metric_col, metric_label, mod, hemis in [
        (df_surface,   'surface_area_mm2', 'Surface area (mm²)',      'surface',   HEMIS_SURFACE),
        (df_volume,    'volume_mm3',        'Volume (mm³)',            'volume',    HEMIS_VOLUME),
        (df_thickness, 'thickness_mm',      'Cortical thickness (mm)', 'thickness', HEMIS_THICKNESS),
    ]:
        if df is None or df.empty:
            continue
        # Match bids_dir by label (basename)
        df_sub = df[(df['species'] == species) & (df['bids_dir'] == label)]
        if df_sub.empty:
            continue

        for hemi in hemis:
            plot_morphometry_intra_bids(
                df_sub, metric_col=metric_col, metric_label=metric_label,
                regions=PLOT_REGIONS,
                output_path=opj(bids_out, f'{mod}_{hemi}.png'),
                species=species, hemisphere=hemi, bids_label=label,
                atlas_level=ATLAS_LEVEL, bids_col=GLOBAL_BIDS_COL,
            )

    # QC per BIDS
    if df_qc is not None and not df_qc.empty:
        qc_sub = df_qc[(df_qc['species'] == species) & (df_qc['bids_dir'] == label)] \
                 if 'bids_dir' in df_qc.columns \
                 else df_qc[df_qc['species'] == species]
        if not qc_sub.empty:
            plot_qc_dashboard(qc_sub, opj(bids_out, 'QC_dashboard.png'))

print(f"\n  All figures → {FIG_DIR}")
print("""
  Folder structure:
  figures/
  ├── cross_species/
  │   ├── morphometry_recap/   (raw / log / norm  x  modality  x  hemisphere)
  │   ├── combo/               (3-row combined  x  raw|log|norm  x  hemisphere)
  │   └── qc_recap/            (QC_cross_species_anat_and_func.png)
  └── per_bids/
      └── <species>_<bids>/    (surface|volume|thickness  x  hemi + QC_dashboard)
""")
print("  Done.")