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
    PAPER_RC, PALETTE, extract_regions_from_legend
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

regionliste = extract_regions_from_legend('/srv/projects/easymribrain/data/Atlas/Classiff/Legende_EDNiX.xlsx')

# Regions shown in every figure panel — must match exact names in xlsx/label files
PLOT_REGIONS = ['Auditory cortex (Superior temporal )', 'Insula and others in lateral sulcus', 'Middle Temporal, Inferior temporal , Temporal pole  (MIPT)', 'Motor and premotor', 'Olfactory cortex', 'Orbital PFC (oPFC)', 'Orbital frontal cortex (oFC)', 'Periarchicortex', 'Posterior medial cortex (PMC)', 'Posterior parietal cortex', 'Prefrontal cortex (PFC)', 'Somatosensory cortex', 'Ventral areas of the temporal lobe (vent Temp)', 'Visual pre and extra striate cortex', 'Visual striate cortex']

REGIONS_OF_INTEREST = PLOT_REGIONS  # filter applied during data collection

# Phylogenetic order — from most distant to closest to Human
PHYLO_ORDER = [
    'Bat', 'Rat', 'Mouse', 'Mouselemur',
    'Marmoset', 'Macaque', 'Human',
    'Dog', 'Pig', 'Cat',
]

# ── anesthesia marker map ─────────────────────────────────────────────────
_ANESTH_MARKER_MAP = {
    'awake': 'o',
    'isoflurane': 's',
    'iso': 's',
    'ketamine': '^',
    'propofol': 'D',
    'medetomidine': 'P',
    'dexmed': 'P',
}


def _anesthesia_marker(anesth_str):
    if not isinstance(anesth_str, str):
        return 'o'
    s = anesth_str.strip().lower()
    for key, m in _ANESTH_MARKER_MAP.items():
        if key in s:
            return m
    return 'o'


def _bids_offsets(bids_dirs, spread=0.22):
    n = len(bids_dirs)
    if n <= 1:
        return {b: 0.0 for b in bids_dirs}
    offs = np.linspace(-spread / 2, spread / 2, n)
    return {b: float(o) for b, o in zip(bids_dirs, offs)}


# Build per-(subject, session) anesthesia map from the study database
anesth_map = {}
_db = pd.read_excel(
    '/home/common/benhalab/CASCAD/EDNiX/databse_ednix_study.xlsx')
_db.columns = _db.columns.str.strip()
for _, _row in _db.iterrows():
    _sub = str(_row.get('subject', '')).strip()
    _ses = str(_row.get('session', '1')).strip()
    _an = str(_row.get('Anesthesia', '')).strip()
    if _sub:
        anesth_map[(_sub, _ses)] = _an
        anesth_map[(_sub, '1')] = anesth_map.get((_sub, '1'), _an)
print(f"  anesth_map: {len(anesth_map)} entries")

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

# ─────────────────────────────────────────────────────────────────────────────
# SPECIES / BIDS DIRECTORIES
# ─────────────────────────────────────────────────────────────────────────────

species_bids_dict = {
    'Rat':         '/scratch2/EDNiX/Rat/BIDS_Grandjean',
    'Mouse':       '/scratch2/EDNiX/Mouse/BIDS_Grandjean',
    'Dog':         '/scratch2/EDNiX/Dog/BIDS_k9',
    'Marmoset':    '/scratch2/EDNiX/Marmoset/BIDS_NIH_MBM',
    'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin',
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
                _bd_list = sorted(sp_df['bids_dir'].unique())
                _bd_offset = _bids_offsets(_bd_list)
                for bd in _bd_list:
                    bd_df = sp_df[sp_df['bids_dir'] == bd]
                    bd_vals = bd_df[metric_col].dropna().values
                    if len(bd_vals) == 0:
                        continue
                    jitter = rng.uniform(-0.10, 0.10, len(bd_vals))
                    x_draw = xi + _bd_offset[bd] + jitter
                    bd_subs = (bd_df.loc[bd_df[metric_col].notna(), 'subject'].values
                               if 'subject' in bd_df.columns else None)
                    bd_sess = (bd_df.loc[bd_df[metric_col].notna(), 'session'].values
                               if 'session' in bd_df.columns else None)
                    # Group by anesthesia marker
                    if anesth_map and bd_subs is not None:
                        from collections import defaultdict
                        _grp = defaultdict(lambda: ([], []))
                        for xv, yv, sub, ses in zip(x_draw, bd_vals, bd_subs, bd_sess):
                            m = _anesthesia_marker(
                                anesth_map.get((str(sub), str(ses)),
                                               anesth_map.get((str(sub), '1'), '')))
                            _grp[m][0].append(xv);
                            _grp[m][1].append(yv)
                        for _m, (_xs, _ys) in _grp.items():
                            ax.scatter(_xs, _ys, color=_bids_col[bd], marker=_m,
                                       s=20, zorder=5, alpha=0.80,
                                       linewidths=0.5, edgecolors='k')
                    else:
                        ax.scatter(x_draw, bd_vals, color=_bids_col[bd], s=20,
                                   zorder=5, alpha=0.80,
                                   linewidths=0.5, edgecolors='k')

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
    Network / FC quality report — violin+dot per species, coloured by BIDS dir.
    Anesthesia regime shown via marker shape (circle=awake, square=iso, etc.).
    """
    import matplotlib.ticker as ticker
    from collections import defaultdict
    with plt.rc_context(PAPER_RC):
        fc_cols = [c for c in (fc_cols or FUNC_NET_COLS) if c in qc_df.columns]
        if not fc_cols:
            warnings.warn('plot_qc_fc_report: no network columns found'); return None

        species_order = _phylo_sort(qc_df['species'].unique())
        bids_dirs = sorted(qc_df['bids_dir'].unique()) \
                    if 'bids_dir' in qc_df.columns else []
        _bids_col = dict(GLOBAL_BIDS_COL)
        for _i, bd in enumerate([b for b in bids_dirs if b not in _bids_col]):
            _bids_col[bd] = PALETTE[(len(_bids_col) + _i) % len(PALETTE)]

        label_map = {
            'net_intra_left_mean':        'Intra L (mean r)',
            'net_intra_right_mean':       'Intra R (mean r)',
            'net_inter_mean':             'Inter hemi (mean r)',
            'net_intra_left_pct_pos':     'Intra L (% pos)',
            'net_intra_right_pct_pos':    'Intra R (% pos)',
            'net_inter_pct_pos':          'Inter (% pos)',
            'net_p_intra_vs_inter':       'p intra vs inter',
            'sp_specific_correlation':    'Specific pair r',
            'sp_nonspecific_correlation': 'Non-specific pair r',
            'sp_specific':               'Specificity (yes=1)',
        }

        NCOLS  = min(4, len(fc_cols))
        n_rows = math.ceil(len(fc_cols) / NCOLS)
        w, h   = figsize or (NCOLS * 3.8, n_rows * 3.5)
        fig, axes = plt.subplots(n_rows, NCOLS, figsize=(w, h))
        axes = np.array(axes).flatten()

        for i, metric in enumerate(fc_cols):
            ax = axes[i]
            cols_needed = ['species', metric]
            if 'bids_dir'  in qc_df.columns: cols_needed.append('bids_dir')
            if 'subject'   in qc_df.columns: cols_needed.append('subject')
            if 'session'   in qc_df.columns: cols_needed.append('session')
            mdf = qc_df[cols_needed].dropna(subset=[metric])

            for xi, sp in enumerate(species_order):
                sp_df = mdf[mdf['species'] == sp]
                vals  = sp_df[metric].values
                if len(vals) == 0:
                    continue
                # Violin coloured by first BIDS of this species
                _bds = sorted(sp_df['bids_dir'].unique()) if 'bids_dir' in sp_df.columns else []
                _vc  = _bids_col.get(_bds[0], '#888888') if _bds else '#888888'
                _violin_strip_quartiles(ax, vals, xi, _vc)

                if 'bids_dir' in sp_df.columns:
                    _bd_list   = sorted(sp_df['bids_dir'].unique())
                    _bd_offset = _bids_offsets(_bd_list)
                    for bd in _bd_list:
                        bd_df   = sp_df[sp_df['bids_dir'] == bd]
                        bd_vals = bd_df[metric].dropna().values
                        if len(bd_vals) == 0:
                            continue
                        jitter = rng.uniform(-0.10, 0.10, len(bd_vals))
                        x_draw = xi + _bd_offset[bd] + jitter
                        dot_color = _bids_col.get(bd, '#888888')
                        # Anesthesia markers for network / functional report
                        if anesth_map and 'subject' in bd_df.columns:
                            bd_subs = bd_df['subject'].values
                            bd_sess = bd_df['session'].values \
                                      if 'session' in bd_df.columns else ['1'] * len(bd_subs)
                            _grp = defaultdict(lambda: ([], []))
                            for xv, yv, sub, ses in zip(x_draw, bd_vals, bd_subs, bd_sess):
                                _m = _anesthesia_marker(
                                    anesth_map.get((str(sub), str(ses)),
                                                   anesth_map.get((str(sub), '1'), '')))
                                _grp[_m][0].append(xv); _grp[_m][1].append(yv)
                            for _m, (_xs, _ys) in _grp.items():
                                ax.scatter(_xs, _ys, color=dot_color, marker=_m,
                                           s=20, zorder=5, alpha=0.80,
                                           linewidths=0.5, edgecolors='k')
                        else:
                            ax.scatter(x_draw, bd_vals, color=dot_color, s=20,
                                       zorder=5, alpha=0.80,
                                       linewidths=0.5, edgecolors='k')

            ax.set_xticks(range(len(species_order)))
            ax.set_xticklabels(species_order, rotation=35, ha='right', fontsize=8)
            ax.set_title(label_map.get(metric, metric.replace('_', ' ')),
                         fontweight='bold', fontsize=9)
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        plt.tight_layout()
        # BIDS colour legend
        if bids_dirs:
            bids_handles = [mpatches.Patch(facecolor=_bids_col.get(bd, '#888'),
                                           label=_bids_label(bd))
                            for bd in bids_dirs]
            fig.legend(handles=bids_handles, loc='lower center',
                       bbox_to_anchor=(0.5, -0.06),
                       ncol=min(len(bids_handles), 6),
                       title='BIDS dir', fontsize=8, frameon=False)
        # Anesthesia marker legend
        if anesth_map:
            _anesth_seen = {}
            for v in anesth_map.values():
                if isinstance(v, str) and v.strip():
                    _anesth_seen.setdefault(v.strip(), _anesthesia_marker(v))
            _ahandles = [
                plt.Line2D([0],[0], marker=_m, color='#333', linestyle='None',
                           markerfacecolor='#333', markersize=7, label=_lbl)
                for _lbl, _m in sorted(_anesth_seen.items())
            ]
            fig.legend(handles=_ahandles, loc='lower right',
                       bbox_to_anchor=(0.99, 0.0), ncol=1,
                       title='Anesthesia', fontsize=7, frameon=True,
                       edgecolor='#ccc')
        fig.subplots_adjust(bottom=0.15)
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
                        _bd_list = sorted(sp_df['bids_dir'].unique())
                        _bd_offset = _bids_offsets(_bd_list)
                        for bd in _bd_list:
                            bd_df = sp_df[sp_df['bids_dir'] == bd]
                            bd_vals = bd_df[metric_col].dropna().values
                            if len(bd_vals) == 0:
                                continue
                            jitter = rng.uniform(-0.10, 0.10, len(bd_vals))
                            x_draw = xi + _bd_offset[bd] + jitter
                            bd_subs = (bd_df.loc[bd_df[metric_col].notna(), 'subject'].values
                                       if 'subject' in bd_df.columns else None)
                            bd_sess = (bd_df.loc[bd_df[metric_col].notna(), 'session'].values
                                       if 'session' in bd_df.columns else None)
                            # Group by anesthesia marker
                            if anesth_map and bd_subs is not None:
                                from collections import defaultdict
                                _grp = defaultdict(lambda: ([], []))
                                for xv, yv, sub, ses in zip(x_draw, bd_vals, bd_subs, bd_sess):
                                    m = _anesthesia_marker(
                                        anesth_map.get((str(sub), str(ses)),
                                                       anesth_map.get((str(sub), '1'), '')))
                                    _grp[m][0].append(xv);
                                    _grp[m][1].append(yv)
                                for _m, (_xs, _ys) in _grp.items():
                                    ax.scatter(_xs, _ys, color=_bids_col[bd], marker=_m,
                                               s=20, zorder=5, alpha=0.80,
                                               linewidths=0.5, edgecolors='k')
                            else:
                                ax.scatter(x_draw, bd_vals, color=_bids_col[bd], s=20,
                                           zorder=5, alpha=0.80,
                                           linewidths=0.5, edgecolors='k')

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



def _anesth_legend_handles_from_map(anesth_map):
    """Build marker legend handles for the anesthesia map used globally."""
    seen = {}
    for v in (anesth_map or {}).values():
        if isinstance(v, str) and v.strip():
            seen.setdefault(v.strip(), _anesthesia_marker(v))
    return [
        plt.Line2D([0],[0], marker=_m, color='#333', linestyle='None',
                   markerfacecolor='#333', markersize=7, label=_lbl)
        for _lbl, _m in sorted(seen.items())
    ]


def _bids_anesth_label(bids_lbl, qc_df):
    """
    Return dominant anesthesia protocol for a BIDS directory,
    inferred from anesth_map + subject/session columns in qc_df.
    Returns a short string like 'awake' or 'iso' or 'mixed'.
    """
    if 'bids_dir' not in qc_df.columns or 'subject' not in qc_df.columns:
        return ''
    bd_df = qc_df[qc_df['bids_dir'] == bids_lbl]
    protocols = set()
    for _, row in bd_df.iterrows():
        sub = str(row.get('subject', ''))
        ses = str(row.get('session', '1'))
        an  = anesth_map.get((sub, ses), anesth_map.get((sub, '1'), ''))
        if an.strip():
            protocols.add(an.strip().lower().split('/')[0].strip())
    if not protocols:
        return ''
    if len(protocols) == 1:
        return list(protocols)[0]
    return 'mixed'


def plot_specificity_bar_per_bids(qc_df, output_path, cat_col='sp_category', figsize=None):
    """
    Stacked bar chart of specificity categories — one bar per BIDS directory.
    BIDS dirs grouped/coloured by species.
    Anesthesia regime shown in legend.
    """
    with plt.rc_context(PAPER_RC):
        if cat_col not in qc_df.columns or 'bids_dir' not in qc_df.columns:
            warnings.warn('plot_specificity_bar_per_bids: missing columns'); return None

        cats   = ['Specific', 'Unspecific', 'Spurious', 'No']
        colors = {'Specific':'#009E73','Unspecific':'#E69F00',
                  'Spurious':'#CC79A7','No':'#D55E00'}

        # Build ordered list: BIDS dirs sorted by species (phylo), then by name
        sp_of_bids = {}
        if 'species' in qc_df.columns:
            for bd in qc_df['bids_dir'].unique():
                sps = qc_df[qc_df['bids_dir'] == bd]['species'].dropna().unique()
                sp_of_bids[bd] = sps[0] if len(sps) > 0 else 'Unknown'
        else:
            for bd in qc_df['bids_dir'].unique():
                sp_of_bids[bd] = 'Unknown'

        phylo_rank = {sp: i for i, sp in enumerate(_phylo_sort(list(set(sp_of_bids.values()))))}
        bids_order = sorted(sp_of_bids.keys(),
                            key=lambda b: (phylo_rank.get(sp_of_bids[b], 99), b))

        counts = {}
        for bd in bids_order:
            sub = qc_df[qc_df['bids_dir'] == bd][cat_col].dropna()
            total = len(sub)
            counts[bd] = {c: 100.0 * (sub == c).sum() / total
                          if total > 0 else 0.0
                          for c in cats}

        n_bd = len(bids_order)
        w, h = figsize or (max(5, n_bd * 1.3), 5.5)
        fig, ax = plt.subplots(figsize=(w, h))

        bar_w   = 0.70
        x       = np.arange(n_bd)
        bottoms = np.zeros(n_bd)
        for cat in cats:
            vals = np.array([counts[bd][cat] for bd in bids_order])
            ax.bar(x, vals, bar_w, bottom=bottoms,
                   color=colors[cat], label=cat,
                   edgecolor='white', linewidth=0.5)
            for xi, (v, b) in enumerate(zip(vals, bottoms)):
                if v > 5:
                    ax.text(xi, b + v/2, f'{v:.0f}%',
                            ha='center', va='center',
                            fontsize=7, color='white', fontweight='bold')
            bottoms += vals

        # x-tick labels: "BidsLabel\n(Species)" with colour per species
        sp_colors_map = {sp: PALETTE[i % len(PALETTE)]
                         for i, sp in enumerate(_phylo_sort(list(set(sp_of_bids.values()))))}
        ax.set_xticks(x)
        tick_labels = [f'{_bids_label(bd)}\n({sp_of_bids[bd]})' for bd in bids_order]
        ax.set_xticklabels(tick_labels, rotation=35, ha='right', fontsize=8)
        for tick, bd in zip(ax.get_xticklabels(), bids_order):
            tick.set_color(sp_colors_map.get(sp_of_bids[bd], '#333'))

        ax.set_ylabel('% subjects', fontsize=10)
        ax.set_ylim(0, 108)
        import matplotlib.ticker as _tk
        ax.yaxis.set_major_formatter(_tk.FormatStrFormatter('%.0f%%'))

        # Category legend
        cat_handles = [mpatches.Patch(facecolor=colors[c], label=c) for c in cats]
        ax.legend(handles=cat_handles, loc='upper center',
                  bbox_to_anchor=(0.5, -0.28), ncol=len(cats),
                  fontsize=8, frameon=False, title='Category')

        # Anesthesia marker legend
        if anesth_map:
            ah = _anesth_legend_handles_from_map(anesth_map)
            if ah:
                fig.legend(handles=ah, loc='upper right',
                           bbox_to_anchor=(0.99, 0.99), ncol=1,
                           title='Anesthesia', fontsize=7,
                           frameon=True, edgecolor='#ccc')

        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')
    return output_path


def _collect_corr_per_bids(species_config, atlas_name, atlas_level, use_lr):
    """
    Like collect_corr_matrices but returns one entry per (species, bids_dir)
    instead of pooling all BIDS dirs of the same species.
    Returns dict: { (species, bids_lbl): {mean, var, pval, rois, n, species, bids_lbl} }
    """
    from scipy import stats as _stats
    result = {}
    for species, cfg in species_config.items():
        bids_dirs = cfg.get('bids_dirs', [])
        if isinstance(bids_dirs, str):
            bids_dirs = [bids_dirs]
        lk = cfg.get('list_to_keep',   [])
        lr = cfg.get('list_to_remove', [])

        for bids_dir in bids_dirs:
            bids_lbl = os.path.basename(bids_dir.rstrip('/'))
            col_key  = (species, bids_lbl)
            records  = extract_corr_matrix_paths(
                bids_dir, atlas_name, atlas_level, use_lr, lk, lr)
            if not records:
                continue

            sub_ses_runs = {}
            rois_ref = None
            for rec in records:
                try:
                    rois, mat = load_corr_matrix(rec['path'])
                    valid = mat[~np.isnan(mat)]
                    if valid.size == 0 or np.all(valid == 0):
                        continue
                    if rois_ref is None:
                        rois_ref = rois
                    key = (rec['subject'], rec['session'])
                    sub_ses_runs.setdefault(key, []).append((rois, mat))
                except Exception as _e:
                    print(f'    [ERROR] {rec["path"]}: {_e}')

            if not sub_ses_runs:
                continue

            subject_means = []
            for (sub, ses), run_list in sorted(sub_ses_runs.items()):
                rois_sets = [set(r) for r, _ in run_list]
                common    = rois_sets[0].intersection(*rois_sets[1:])
                ref_rois  = [r for r in run_list[0][0] if r in common]
                if not ref_rois: continue
                mats = []
                for rois_r, mat_r in run_list:
                    try:
                        idx = [list(rois_r).index(r) for r in ref_rois]
                        mats.append(mat_r[np.ix_(idx, idx)])
                    except ValueError:
                        pass
                if mats:
                    subject_means.append(
                        (ref_rois,
                         np.nanmean(np.stack(mats, axis=0), axis=0)))

            if not subject_means:
                continue

            all_roi_sets = [set(r) for r, _ in subject_means]
            common_sp    = all_roi_sets[0].intersection(*all_roi_sets[1:])
            rois_ref     = [r for r in subject_means[0][0] if r in common_sp]

            aligned = []
            for rois_s, mat_s in subject_means:
                idx = [list(rois_s).index(r) for r in rois_ref if r in rois_s]
                if len(idx) == len(rois_ref):
                    aligned.append(mat_s[np.ix_(idx, idx)])

            if not aligned:
                continue

            stack = np.stack(aligned, axis=0)
            n     = stack.shape[0]
            mean  = np.nanmean(stack, axis=0)
            var   = np.nanvar(stack,  axis=0)
            _, pv = _stats.ttest_1samp(stack, 0, axis=0, nan_policy='omit')                     if n >= 2 else (None, np.full(mean.shape, np.nan))

            result[col_key] = dict(mean=mean, var=var, pval=pv,
                                   rois=rois_ref, n=n,
                                   species=species, bids_lbl=bids_lbl)
            print(f'  [per_bids_corr] {species}/{bids_lbl}: n={n} rois={len(rois_ref)}')

    # Restrict all entries to common ROIs
    if result:
        all_roi_sets = [set(v['rois']) for v in result.values()]
        common_rois_set = all_roi_sets[0].intersection(*all_roi_sets[1:])
        ref_order = next(iter(result.values()))['rois']
        common_rois = [r for r in ref_order if r in common_rois_set]
        for col_key, d in result.items():
            lr_rois = d['rois']
            idx = [lr_rois.index(r) for r in common_rois]
            d['mean'] = d['mean'][np.ix_(idx, idx)]
            d['var']  = d['var'][np.ix_(idx, idx)]
            d['pval'] = d['pval'][np.ix_(idx, idx)]
            d['rois'] = common_rois

    return result


def plot_corr_matrix_combo_per_bids(
    species_config, output_path,
    atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
    figsize=None, vmin=-0.5, vmax=0.5,
):
    """
    Correlation matrix (mean r | -log10 p | variance) — one column per BIDS dir.
    Title of each column: "Species / BIDSdir (n=X, anesthesia)".
    Anesthesia marker legend added from global anesth_map.
    """
    import matplotlib.gridspec as gridspec
    import matplotlib.cm as _cm
    from matplotlib.colors import TwoSlopeNorm, Normalize

    bids_data = _collect_corr_per_bids(species_config, atlas_name, atlas_level, use_lr)
    if not bids_data:
        warnings.warn('plot_corr_matrix_combo_per_bids: no data'); return None

    # Order columns: phylogenetic species order, then bids_lbl alphabetically
    phylo_rank = {sp: i for i, sp in enumerate(_phylo_sort(
        [d['species'] for d in bids_data.values()]))}
    cols = sorted(bids_data.keys(),
                  key=lambda k: (phylo_rank.get(bids_data[k]['species'], 99),
                                 bids_data[k]['bids_lbl']))

    common_rois = bids_data[cols[0]]['rois']
    short_rois  = [r.replace('L_','L ').replace('R_','R ')[:16] for r in common_rois]
    n_roi = len(common_rois)
    n_cols = len(cols)
    n_rows = 3

    norm_r   = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    cmap_r   = 'RdBu_r'
    cmap_var = 'viridis'
    all_var  = np.concatenate([bids_data[k]['var'].ravel() for k in cols])
    vmax_var = float(np.nanpercentile(all_var[~np.isnan(all_var)], 95))                if all_var[~np.isnan(all_var)].size > 0 else 0.01
    norm_var = Normalize(0, vmax_var)

    w, h = figsize or (max(8, n_cols * 3.5), n_rows * 3.8)
    fig  = plt.figure(figsize=(w, h), facecolor='white')
    gs   = gridspec.GridSpec(n_rows, n_cols + 1,
                             width_ratios=[1.0]*n_cols + [0.05],
                             wspace=0.04, hspace=0.15)

    for c_i, col_key in enumerate(cols):
        d   = bids_data[col_key]
        # Build title with anesthesia info
        an_str = _bids_anesth_label(d['bids_lbl'],
                                    df_qc if df_qc is not None else pd.DataFrame())
        title  = f"{d['species']} / {d['bids_lbl']} (n={d['n']}" + (f", {an_str})" if an_str else ")")

        ax0 = fig.add_subplot(gs[0, c_i])
        ax0.imshow(d['mean'], cmap=cmap_r, norm=norm_r,
                   aspect='auto', interpolation='nearest')
        ax0.set_title(title, fontweight='bold', fontsize=7)
        if c_i == 0:
            ax0.set_ylabel('Mean r', fontsize=8)
            ax0.set_yticks(range(n_roi)); ax0.set_yticklabels(short_rois, fontsize=4)
        else:
            ax0.set_yticks([])
        ax0.set_xticks([])

        ax1 = fig.add_subplot(gs[1, c_i])
        ax1.imshow(-np.log10(np.clip(d['pval'], 1e-10, 1.0)),
                   cmap='hot_r', aspect='auto', vmin=0, vmax=4,
                   interpolation='nearest')
        if c_i == 0:
            ax1.set_ylabel('-log10(p)', fontsize=8)
            ax1.set_yticks(range(n_roi)); ax1.set_yticklabels(short_rois, fontsize=4)
        else:
            ax1.set_yticks([])
        ax1.set_xticks([])

        ax2 = fig.add_subplot(gs[2, c_i])
        ax2.imshow(d['var'], cmap=cmap_var, norm=norm_var,
                   aspect='auto', interpolation='nearest')
        if c_i == 0:
            ax2.set_ylabel('Variance', fontsize=8)
            ax2.set_yticks(range(n_roi)); ax2.set_yticklabels(short_rois, fontsize=4)
        else:
            ax2.set_yticks([])
        ax2.set_xticks([])

    fig.colorbar(_cm.ScalarMappable(norm=norm_r, cmap=cmap_r),
                 cax=fig.add_subplot(gs[0,-1]), label='Pearson r').ax.yaxis.set_label_position('right')
    fig.colorbar(_cm.ScalarMappable(norm=Normalize(0,4), cmap='hot_r'),
                 cax=fig.add_subplot(gs[1,-1]), label='-log10(p)').ax.yaxis.set_label_position('right')
    fig.colorbar(_cm.ScalarMappable(norm=norm_var, cmap=cmap_var),
                 cax=fig.add_subplot(gs[2,-1]), label='Variance').ax.yaxis.set_label_position('right')

    plt.tight_layout()
    # Anesthesia legend
    ah = _anesth_legend_handles_from_map(anesth_map)
    if ah:
        fig.legend(handles=ah, loc='lower right', bbox_to_anchor=(0.99, 0.0),
                   ncol=1, title='Anesthesia', fontsize=7, frameon=True, edgecolor='#ccc')

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, bbox_inches='tight', dpi=200)
    plt.close(fig)
    print(f'  [plot] {output_path}')
    return output_path


def plot_corr_matrix_sig_norm_per_bids(
    species_config, output_path,
    atlas_name="EDNIxCSC", atlas_level=2, use_lr=True,
    figsize=None,
):
    """
    Normalised FC matrix (r / mean|r|) — one panel per BIDS directory.
    Anesthesia regime shown in legend.
    """
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.cm as _cm
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import Normalize

    bids_data = _collect_corr_per_bids(species_config, atlas_name, atlas_level, use_lr)
    if not bids_data:
        warnings.warn('plot_corr_matrix_sig_norm_per_bids: no data'); return None

    phylo_rank = {sp: i for i, sp in enumerate(_phylo_sort(
        [d['species'] for d in bids_data.values()]))}
    cols = sorted(bids_data.keys(),
                  key=lambda k: (phylo_rank.get(bids_data[k]['species'], 99),
                                 bids_data[k]['bids_lbl']))

    common_rois = bids_data[cols[0]]['rois']
    short_rois  = [r.replace('L_','L ').replace('R_','R ')[:16] for r in common_rois]
    n_roi  = len(common_rois)
    n_cols = len(cols)

    norm = TwoSlopeNorm(vmin=-3.0, vcenter=0, vmax=3.0)
    cmap = 'RdBu_r'

    w, h = figsize or (max(6, n_cols * 3.5), 5.5)
    fig  = plt.figure(figsize=(w, h), facecolor='white')
    gs   = gridspec.GridSpec(1, n_cols + 1,
                             width_ratios=[1.0]*n_cols + [0.05],
                             wspace=0.04)

    for c_i, col_key in enumerate(cols):
        d = bids_data[col_key]
        mat_mean_r = float(np.nanmean(np.abs(d['mean'])))
        sig_norm   = d['mean'] / mat_mean_r if mat_mean_r > 1e-6 else d['mean'].copy()

        an_str = _bids_anesth_label(d['bids_lbl'],
                                    df_qc if df_qc is not None else pd.DataFrame())
        title  = f"{d['species']} / {d['bids_lbl']} (n={d['n']}" + (f", {an_str})" if an_str else ")")

        ax = fig.add_subplot(gs[0, c_i])
        ax.imshow(sig_norm, cmap=cmap, norm=norm,
                  aspect='equal', interpolation='nearest')
        ax.set_title(title, fontweight='bold', fontsize=7)
        if c_i == 0:
            ax.set_yticks(range(n_roi))
            ax.set_yticklabels(short_rois, fontsize=4)
            ax.set_ylabel('r / mean(|r|)', fontsize=8)
        else:
            ax.set_yticks([])
        ax.set_xticks([])

    fig.colorbar(_cm.ScalarMappable(norm=norm, cmap=cmap),
                 cax=fig.add_subplot(gs[0,-1]),
                 label='r / mean(|r|)').ax.yaxis.set_label_position('right')

    plt.tight_layout()
    ah = _anesth_legend_handles_from_map(anesth_map)
    if ah:
        fig.legend(handles=ah, loc='lower right', bbox_to_anchor=(0.99, 0.0),
                   ncol=1, title='Anesthesia', fontsize=7, frameon=True, edgecolor='#ccc')

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
            # Allow negative FC values → use linear scale for FC, log for weight*
            print(xv)
            print(yv)
            print(x_col)
            print(y_col)
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
    # Specificity category bar chart — by species (original)
    plot_specificity_bar(
        df_qc,
        opj(qc_dir, 'QC_cross_species_specificity_bar.png'),
        cat_col=FUNC_SPEC_COL,
    )
    # Specificity category bar chart — by BIDS dir (extra figure)
    plot_specificity_bar_per_bids(
        df_qc,
        opj(qc_dir, 'QC_cross_species_specificity_bar_per_bids.png'),
        cat_col=FUNC_SPEC_COL,
    )

# Correlation matrix combo — by species (original)
plot_corr_matrix_combo(
    species_config,
    opj(combo_dir, f'combo_corr_matrix_lvl{CORR_ATLAS_LEVEL}.png'),
    atlas_name=ATLAS_NAME,
    atlas_level=CORR_ATLAS_LEVEL,
    use_lr=CORR_USE_LR,
)
# Correlation matrix combo — by BIDS dir (extra figure)
plot_corr_matrix_combo_per_bids(
    species_config,
    opj(combo_dir, f'combo_corr_matrix_per_bids_lvl{CORR_ATLAS_LEVEL}.png'),
    atlas_name=ATLAS_NAME,
    atlas_level=CORR_ATLAS_LEVEL,
    use_lr=CORR_USE_LR,
)

# Normalized significant FC — by species (original)
plot_corr_matrix_sig_norm(
    species_config,
    opj(combo_dir, f'combo_corr_matrix_sig_norm_lvl{CORR_ATLAS_LEVEL}.png'),
    atlas_name=ATLAS_NAME,
    atlas_level=CORR_ATLAS_LEVEL,
    use_lr=CORR_USE_LR,
    pval_thresh=1.0,
)
# Normalized significant FC — by BIDS dir (extra figure)
plot_corr_matrix_sig_norm_per_bids(
    species_config,
    opj(combo_dir, f'combo_corr_matrix_sig_norm_per_bids_lvl{CORR_ATLAS_LEVEL}.png'),
    atlas_name=ATLAS_NAME,
    atlas_level=CORR_ATLAS_LEVEL,
    use_lr=CORR_USE_LR,
)

def _ols_on_ax(ax, xv, yv, log_x):
    """
    Fit OLS on (xv, yv) and draw a dashed regression line + 95% prediction
    band on `ax`.

    Parameters
    ----------
    ax     : matplotlib Axes
    xv, yv : 1-D float arrays (may contain NaN / inf)
    log_x  : bool ? if True, x is log10-transformed before regression;
             the axis scale must be set to 'log' by the caller.

    Returns
    -------
    (slope, r, p) on success, or None when degenerate (< 3 finite points,
    or all x identical).  A horizontal median line is drawn in the degenerate
    case so the panel is never blank.
    """
    from scipy import stats as _st

    lx    = np.log10(np.clip(xv, 1e-9, None)) if log_x else np.array(xv, dtype=float)
    valid = np.isfinite(lx) & np.isfinite(yv)
    if valid.sum() < 3:
        return None

    lxv, lyv = lx[valid], np.array(yv, dtype=float)[valid]
    ssx = float(np.sum((lxv - lxv.mean()) ** 2))

    if ssx < 1e-12:
        # Degenerate: all x-values are identical (e.g. one species, one weight)
        ax.axhline(float(np.median(lyv)), color='k',
                   lw=1.2, linestyle='--', alpha=0.55)
        return None

    slope, intercept, r, p, _ = _st.linregress(lxv, lyv)
    n_v    = int(valid.sum())
    t_crit = _st.t.ppf(0.975, df=n_v - 2)
    x_fit  = np.linspace(lxv.min(), lxv.max(), 200)
    y_fit  = intercept + slope * x_fit
    resid  = np.sum((lyv - (intercept + slope * lxv)) ** 2) / (n_v - 2)
    se     = np.sqrt(resid * (1 + 1 / n_v + (x_fit - lxv.mean()) ** 2 / ssx))
    x_plot = 10 ** x_fit if log_x else x_fit

    ax.plot(x_plot, y_fit,
            color='k', lw=1.5, linestyle='--', alpha=0.8, zorder=3)
    ax.fill_between(x_plot,
                    y_fit - t_crit * se,
                    y_fit + t_crit * se,
                    color='k', alpha=0.15, zorder=2)
    p_str = f'{p:.3f}' if p >= 0.001 else '<0.001'
    ax.text(0.05, 0.95,
            f'slope={slope:.2f}  r={r:.2f}  p={p_str}',
            transform=ax.transAxes, fontsize=8,
            va='top', color='#333333')
    return slope, r, p



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PART 5  ?  Figure 2: summary heatmap (all connections × all species)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def plot_fc_connections_summary(df_fc, connections, output_path, figsize=None):
    """
    Heatmap where:
      rows    = FC connections (in the order listed in `connections`)
      columns = species (phylogenetic order)
      cell    = median normalised FC (r / mean|r|) across all subjects

    Cells are annotated with the numeric median.
    A diverging RdBu_r colormap is centred at 0; the colour range is set
    to the 95th percentile of |median| across all cells.
    """
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.cm as _cm

    species_order = _phylo_sort(df_fc['species'].unique())

    # Restrict to connections that actually produced data
    keys   = [k for k, *_ in connections
              if k in df_fc.columns and df_fc[k].notna().any()]
    labels = {k: lbl for k, *_, lbl in connections}

    if not keys:
        warnings.warn('plot_fc_connections_summary: no data')
        return None

    # Build median matrix  (n_connections × n_species)
    mat = np.full((len(keys), len(species_order)), np.nan)
    for ki, k in enumerate(keys):
        for si, sp in enumerate(species_order):
            vals = df_fc.loc[df_fc['species'] == sp, k].dropna().values
            if len(vals):
                mat[ki, si] = float(np.median(vals))

    finite = mat[np.isfinite(mat)]
    vext   = float(np.nanpercentile(np.abs(finite), 95)) if finite.size > 0 else 1.0
    norm   = TwoSlopeNorm(vmin=-vext, vcenter=0, vmax=vext)
    cmap   = 'RdBu_r'

    n_r, n_c = mat.shape
    w, h = figsize or (max(5, n_c * 1.0 + 3.5), max(3, n_r * 0.6 + 1.5))

    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(w, h))
        ax.imshow(mat, cmap=cmap, norm=norm,
                  aspect='auto', interpolation='nearest')

        # Annotate each cell with its median value
        for ki in range(n_r):
            for si in range(n_c):
                v = mat[ki, si]
                if np.isfinite(v):
                    text_color = 'white' if abs(v) > 0.45 * vext else '#222'
                    ax.text(si, ki, f'{v:.2f}',
                            ha='center', va='center',
                            fontsize=7, color=text_color)

        ax.set_xticks(range(n_c))
        ax.set_xticklabels(species_order, rotation=30, ha='right', fontsize=9)
        ax.set_yticks(range(n_r))
        ax.set_yticklabels([labels.get(k, k) for k in keys], fontsize=8)
        ax.set_title('FC connections of interest  -  median  r / mean|r|',
                     fontweight='bold', fontsize=10)

        cb = fig.colorbar(_cm.ScalarMappable(norm=norm, cmap=cmap),
                          ax=ax, fraction=0.03, pad=0.03)
        cb.set_label('r / mean|r|', fontsize=8)
        cb.ax.yaxis.set_label_position('right')

        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches='tight', dpi=200)
        plt.close(fig)
        print(f'  [plot] {output_path}')

    return output_path

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