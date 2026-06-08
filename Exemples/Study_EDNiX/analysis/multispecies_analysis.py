"""
EDNiX BIDS Tools v6  —  Main analysis script
=============================================
Processes all species: surface, thickness, volume, QC
Produces DataFrames, exports and paper-quality figures.

FIXES vs previous version
--------------------------
★ STEP 0: run_wb_extraction uses n_atlas_levels=N_ATLAS_LEVELS (was 0 → nothing written)
★ STEP 0: removed the redundant extract_and_process_*_from_rois calls that
          followed run_wb_extraction — those are only needed when wb_command
          already ran but xlsx files are missing/corrupt.
★ STEP 4: per-subject session averaging applied to all modalities BEFORE
          make_all_figures() — one dot per animal on the plots.
          Raw multi-session data is still exported to CSV/xlsx unchanged.
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Plotting.ednix_bids_tools import (
    PAPER_RC, PALETTE,
    get_atlas_label_path, find_species_path,
    extract_corr_matrix_paths, load_corr_matrix,
    collect_multi_species,
    run_wb_extraction,
    extract_and_process_surfaces_from_rois,
    extract_and_process_thickness_from_rois,
    extract_regions_from_legend,
    export_to_excel, export_to_csv, export_summary_stats,
    _bids_label,
)
from Plotting.EDNiX_figures import make_all_figures

opj = os.path.join

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

ATLAS_LIB      = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
ATLAS_NAME     = 'EDNIxCSC'
ATLAS_LEVEL    = 1   # hierarchy level shown in figures (1–4)
N_ATLAS_LEVELS = 1   # how many levels Stage 0 extracts — must be >= 1

OUT_DIR = '/scratch2/EDNiX/results/multispecies_analysis'
FIG_DIR = opj(OUT_DIR, 'figures')

SING_WB = (
    'vglrun singularity run '
    '--bind /srv/projects/,/srv/projects/easymribrain,/scratch2/,/scratch/ '
    '/home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/'
    'connectome_workbench_1.5.0-freesurfer-update.sif'
)

PLOT_REGIONS        = ['Isocortex', 'Allocortex', 'Periallocortex']
REGIONS_OF_INTEREST = PLOT_REGIONS

CORR_ATLAS_LEVEL = 2
CORR_USE_LR      = True

# ═══════════════════════════════════════════════════════════════════════════════
# SPECIES / BIDS DIRECTORIES
# ═══════════════════════════════════════════════════════════════════════════════

species_bids_dict = {
    'Rat':        '/scratch2/EDNiX/Rat/BIDS_Grandjean',
    'Mouse':      '/scratch2/EDNiX/Mouse/BIDS_Grandjean2',
    'Dog':        '/scratch2/EDNiX/Dog/BIDS_Boch_K9',
    'Marmoset':   '/scratch2/EDNiX/Marmoset/BIDS_Tian',
    'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin',
}
species_multi_bids = {
    'Macaque': [

        '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
        '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin',
    ],
    'Human': [
        '/scratch2/EDNiX/Human/BIDS_Merida',
        '/scratch2/EDNiX/Human/BIDS_Park',
        '/scratch2/EDNiX/Human/BIDS_Castrillon',
    ],
}

ALL_BIDS = (
    [(sp, bd) for sp, bd in species_bids_dict.items()] +
    [(sp, bd) for sp, bds in species_multi_bids.items() for bd in bds]
)

# ═══════════════════════════════════════════════════════════════════════════════
# ANESTHESIA MAP
# ═══════════════════════════════════════════════════════════════════════════════

anesth_map = {}
_db = pd.read_excel('/home/common/benhalab/CASCAD/EDNiX/databse_ednix_study.xlsx')
_db.columns = _db.columns.str.strip()
for _, _row in _db.iterrows():
    _sub = str(_row.get('subject', '')).strip()
    _ses = str(_row.get('session', '1')).strip()
    _an  = str(_row.get('Anesthesia', '')).strip()
    if _sub:
        anesth_map[(_sub, _ses)] = _an
        anesth_map.setdefault((_sub, '1'), _an)
print(f"  anesth_map: {len(anesth_map)} entries")

# ═══════════════════════════════════════════════════════════════════════════════
# QC COLUMN GROUPS
# ═══════════════════════════════════════════════════════════════════════════════

FUNC_QC_COLS = ['func_TSNR_0', 'func_mean_fd', 'func_gcor']
ANAT_QC_COLS = ['anat_template_correlation', 'anat_cortical_contrast']
FUNC_NET_COLS = [
    'net_mean_correlation', 'net_std_correlation',
    'net_eigenvalue_ratio', 'net_davies_bouldin',
    'net_inter_mean', 'net_intra_mean']

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 0  —  Stage 0: wb_command surface area + thickness extraction
# ═══════════════════════════════════════════════════════════════════════════════
# Run once per dataset.  After a successful run set both overwrite flags to
# False to skip already-computed subjects.
#
# KEY POINTS:
#   • n_atlas_levels must be >= 1  (was accidentally 0 → empty xlsx)
#   • Do NOT call extract_and_process_*_from_rois after this — that helper
#     is only for rebuilding xlsx when ROI files exist but xlsx are missing.
'''
run_wb_extraction(
    SING_WB,
    ALL_BIDS,
    atlas_name='EDNIx',
    regions_to_process  = REGIONS_OF_INTEREST,
    overwrite_thickness = True,   # ← set False after first successful run
    overwrite_surface   = True,   # ← set False after first successful run
    n_atlas_levels      = N_ATLAS_LEVELS,
)
'''

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1  —  Build species_config and resolve atlas paths
# ═══════════════════════════════════════════════════════════════════════════════

species_config          = {}
atlas_label_paths       = {}
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
# STEP 2  —  Collect all data  (reads xlsx / NIfTI / JSON — no wb_command)
# ═══════════════════════════════════════════════════════════════════════════════

data = collect_multi_species(
    species_config,
    regions_of_interest     = REGIONS_OF_INTEREST,
    extract                 = ('surface', 'volume', 'thickness', 'qc'),
    atlas_name              = ATLAS_NAME,
    atlas_label_paths       = atlas_label_paths,
    fit_kind='correlation',
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
    print(f"  {name:10s}: {n} rows (raw, multi-session)")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3  —  Export raw multi-session data (unchanged)
# ═══════════════════════════════════════════════════════════════════════════════

os.makedirs(OUT_DIR, exist_ok=True)

export_to_excel(data, opj(OUT_DIR, 'ednix_results.xlsx'))
export_to_csv(data,   opj(OUT_DIR, 'csv'))
export_summary_stats(data, opj(OUT_DIR, 'ednix_summary_stats.xlsx'),
                     groupby=('species', 'bids_dir'))

if df_qc is not None and not df_qc.empty:
    id_cols = [c for c in ('species', 'subject', 'session', 'bids_dir')
               if c in df_qc.columns]
    anat_p = [c for c in ANAT_QC_COLS  if c in df_qc.columns]
    func_p = [c for c in FUNC_QC_COLS  if c in df_qc.columns]
    net_p  = [c for c in FUNC_NET_COLS if c in df_qc.columns]
    with pd.ExcelWriter(opj(OUT_DIR, 'qc_report.xlsx'), engine='openpyxl') as w:
        if func_p:
            (df_qc[id_cols + func_p].dropna(subset=func_p, how='all')
             .to_excel(w, sheet_name='func_qc', index=False))
        if anat_p:
            (df_qc[id_cols + anat_p].dropna(subset=anat_p, how='all')
             .to_excel(w, sheet_name='anat_qc', index=False))
        if net_p:
            (df_qc[id_cols + net_p].dropna(subset=net_p, how='all')
             .to_excel(w, sheet_name='network_qc', index=False))
    print(f"  QC report → qc_report.xlsx "
          f"(func:{len(func_p)} anat:{len(anat_p)} net:{len(net_p)} cols)")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4  —  Average sessions per subject for plotting
#
# Multiple sessions per subject → multiple dots per species tick → visual
# bimodality.  We collapse to one value per (species, bids_dir, subject,
# region, hemisphere, atlas_level) by taking the mean across sessions.
#
# The raw data exported in STEP 3 is untouched.
# ═══════════════════════════════════════════════════════════════════════════════

def _mean_per_subject(df, metric_col):
    """Average metric_col across sessions → one row per subject-region."""
    if df is None or df.empty or metric_col not in df.columns:
        return df
    if 'session' not in df.columns:
        return df

    # Identity columns: everything except session, metric, and raw-count cols
    exclude = {metric_col, 'session', 'voxel_count', 'label_id',
               'region_name', 'region_full'}
    id_cols = [c for c in df.columns if c not in exclude]

    agg      = df.groupby(id_cols, as_index=False)[metric_col].mean()
    n_before = len(df)
    n_after  = len(agg)
    if n_before != n_after:
        print(f"  [session_avg] {metric_col}: "
              f"{n_before} → {n_after} rows "
              f"({n_before - n_after} collapsed from multi-session subjects)")
    return agg


df_surface_plot   = _mean_per_subject(df_surface,   'surface_area_mm2')
df_volume_plot    = _mean_per_subject(df_volume,     'volume_mm3')
df_thickness_plot = _mean_per_subject(df_thickness,  'thickness_mm')

# QC: average numeric metrics per subject, keep first-occurrence categoricals
if df_qc is not None and not df_qc.empty and 'session' in df_qc.columns:
    _qc_id   = [c for c in ('species', 'bids_dir', 'subject') if c in df_qc.columns]
    _qc_skip = set(_qc_id) | {'session'}
    _qc_num  = [c for c in df_qc.columns
                if c not in _qc_skip
                and pd.api.types.is_numeric_dtype(df_qc[c])]
    _qc_cat  = [c for c in df_qc.columns
                if c not in _qc_skip and c not in _qc_num]
    df_qc_plot = df_qc.groupby(_qc_id, as_index=False)[_qc_num].mean()
    if _qc_cat:
        _cat_first = (df_qc.sort_values('session')
                           .groupby(_qc_id, as_index=False)[_qc_cat]
                           .first())
        df_qc_plot = df_qc_plot.merge(_cat_first, on=_qc_id, how='left')
    print(f"  [session_avg] qc: {len(df_qc)} → {len(df_qc_plot)} rows")
else:
    df_qc_plot = df_qc

data_plot = {
    'surface':   df_surface_plot,
    'volume':    df_volume_plot,
    'thickness': df_thickness_plot,
    'qc':        df_qc_plot,
}

for name, df in [('surface', df_surface_plot), ('volume', df_volume_plot),
                 ('thickness', df_thickness_plot), ('qc', df_qc_plot)]:
    n = len(df) if df is not None and not df.empty else 0
    print(f"  {name:10s}: {n} rows (plot, one per subject)")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5  —  Build global BIDS colour map
# ═══════════════════════════════════════════════════════════════════════════════

_bd_sources = []
for df in (df_volume_plot, df_surface_plot, df_thickness_plot, df_qc_plot):
    if df is not None and not df.empty and 'bids_dir' in df.columns:
        _bd_sources.extend(df['bids_dir'].unique().tolist())

_all_bids_dirs = sorted(set(
    os.path.basename(str(bd).rstrip('/')) for bd in _bd_sources))
GLOBAL_BIDS_COL = {
    bd: PALETTE[i % len(PALETTE)] for i, bd in enumerate(_all_bids_dirs)
}
print(f"\n  Global BIDS colour map: {GLOBAL_BIDS_COL}")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6  —  Generate all figures (per-subject averaged data)
# ═══════════════════════════════════════════════════════════════════════════════

make_all_figures(
    data             = data_plot,
    species_config   = species_config,
    all_bids         = ALL_BIDS,
    fig_dir          = FIG_DIR,
    plot_regions     = PLOT_REGIONS,
    atlas_level      = ATLAS_LEVEL,
    atlas_name       = ATLAS_NAME,
    corr_atlas_level = CORR_ATLAS_LEVEL,
    corr_use_lr      = CORR_USE_LR,
    anesth_map       = anesth_map,
    global_bids_col  = GLOBAL_BIDS_COL,
    trimouse_xlsx='/home/common/benhalab/CASCAD/EDNiX/Mouse/BIDS_Grandjean/Trimouse-Grandjean.xlsx'
)

print(f"\n  Done.  All outputs → {OUT_DIR}")
print("""
  Folder structure:
  figures/
  ├── cross_species/
  │   ├── morphometry_recap/   (raw / log / norm  ×  modality  ×  hemisphere)
  │   ├── combo/               (3-row combined + corr matrices per species & BIDS)
  │   ├── qc_recap/            (anat+func QC, network QC, specificity by species & BIDS)
  │   └── brain_scaling/       (allometric scaling per modality)
  └── per_bids/
      └── <species>_<bids>/    (surface|volume|thickness × hemi + QC_dashboard)
""")

from Plotting.ednix_threshold_explorer import run_pipeline_multilevel, run_pipeline
'''
result = run_pipeline(
    data["qc"],
    bids_root_template = "/scratch2/EDNiX/{species}/{bids_dir}",
    n_top    = 10,
    stringency=0,
    fig_dir  = "/scratch2/EDNiX/results/cross_species/qc_recap")
'''
result = run_pipeline_multilevel(data["qc"], fit_kind="correlation",
                            bids_root_template="/scratch2/EDNiX/{species}/{bids_dir}",
                            atlas_levels=(2, 3, 4),
                            subsets=("all", "primates", "primates_rodents"),
                            thresh_intra=None, thresh_delta=None, stringency=0.0,
                            n_top=10, atlas_name="EDNIxCSC", use_lr=True,
                            fp_test="species_mean", fdr_alpha=0.05,
                            fig_dir  = "/scratch2/EDNiX/results/cross_species/qc_recap",
                            verbose=True)