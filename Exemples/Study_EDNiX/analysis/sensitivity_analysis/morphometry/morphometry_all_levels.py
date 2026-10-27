"""
Morphometry -- every level, every hemisphere, every scaling (sensitivity analysis).

Former Figure_2_anat.py, without the paper figure and without the PGLS part:
  * the paper Figure 2 (level 2, bilateral, log10) is Figure_2_morphometry.py;
  * the outlier rates of the paper (Supp. Figs 6-8, within-dataset allometric residuals)
    are Supplementary/Supp_Fig_06_08_outliers.py;
  * PGLS (frontal scaling vs body/brain/neurons, frontal vs occipital) is
    phylo_preliminary/PGLS_frontal_scaling.py.

What stays here, all written under OUT_DIR (unchanged, other scripts read it):
  figures/Figure2|supplementary/...  dots + combo plots for levels 2-4 (raw/log/norm,
                                     left/right/bilateral)
  stats/correlations_lvl2.csv        surface/volume/thickness species-mean correlations
                                     (log10, level 2, BH-FDR)
  supplementary/outliers_raw_pooled_<mod>_lvl<k>.png  MAD outlier rates with the
                                     dataset pooled within species (former Suppl7-9;
                                     mixes scanners / ages, hence a sensitivity check)
  csv / xlsx exports of the collected tables.
"""

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
    PALETTE, get_atlas_label_path, find_species_path,
    collect_multi_species,
    export_to_excel, export_to_csv, export_summary_stats,)

opj = os.path.join

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIG — single source of truth
# ═══════════════════════════════════════════════════════════════════════════════

ATLAS_LIB      = "/home/cgarin/PycharmProjects/EDNiX/Atlases_library"
ATLAS_NAME     = "EDNIxCSC"
NEWICK_FILE    = "/scratch2/EDNiX/Liste des espèces.nwk"
HERCULANO_XLSX = "/scratch2/EDNiX/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx"

# results roots (everything is in EDNiX, NOT EDNiX)
RESULTS_DIR    = "/scratch2/EDNiX/results/multispecies_analysis"
CSV_DIR        = opj(RESULTS_DIR, "csv")
COMBO_DIR      = opj(RESULTS_DIR, "figures", "cross_species", "combo")   # make_all_figures output
FC_RESULTS_DIR = "/scratch2/EDNiX/results/cross_species/qc_recap"
ARTICLE_DIR    = "/scratch2/EDNiX/results/article_figures"

# Longitudinal datasets (Fig 3)
MACAQUE_LONGI_BIDS  = ["/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin"]
MACAQUE_AGE_XLSX    = "/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin/Garin_macaque.xlsx"
HUMAN_AGING_BIDS    = ["/scratch2/EDNiX/Human/BIDS_Park"]
HUMAN_AGE_PREFIX    = "AgeMRI_W"

ATLAS_LEVEL    = 1   # hierarchy level shown in figures (1–4)
N_ATLAS_LEVELS = 4   # how many levels Stage 0 extracts — must be >= 1

OUT_DIR = '/scratch2/EDNiX/results/EDNiX_Figures/Figure_2_anat'

SING_WB = (
    'vglrun singularity run '
    '--bind /srv/projects/,/srv/projects/easymribrain,/scratch2/,/scratch/ '
    '/home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/'
    'connectome_workbench_1.5.0-freesurfer-update.sif'
)

PLOT_REGIONS        = []
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
        '/scratch2/EDNiX/Human/BIDS_Babayan',
    ],
    'Marmoset': [
        '/scratch2/EDNiX/Marmoset/BIDS_Tian',
        '/scratch2/EDNiX/Marmoset/BIDS_Hori',
    ],
}

ALL_BIDS = (
    [(sp, bd) for sp, bd in species_bids_dict.items()] +
    [(sp, bd) for sp, bds in species_multi_bids.items() for bd in bds]
)

# Module-level species ordering / colours
PHYLO_ORDER = ["Mouse", "Rat", "Mouselemur", "Bat", "Marmoset",
               "Macaque", "Human", "Dog", "Cat", "Pig"]

# Phylogenetic gradient - EACH major clade has a DISTINCT color
SPECIES_COLORS = {
    # Rodents (cool blues)
    "Mouse": "#1A5276",  # Dark blue
    "Rat": "#2E86C1",  # Medium blue

    # Bat (sister to primates - teal/cyan)
    "Bat": "#00A896",  # Teal (completely different from blue)

    # Primates (warm - green → yellow → orange → warm)
    "Mouselemur": "#58D68D",  # Green (basal primates)
    "Marmoset": "#F4D03F",  # Yellow (New World)
    "Macaque": "#F5B041",  # Orange (Old World)
    "Human": "#D35400",  # Burnt orange (Great Apes)

    # Carnivores (cool purple → magenta)
    "Dog": "#8E44AD",  # Purple (Canidae)
    "Cat": "#C0392B",  # Red (Felidae) - NOT orange!

    # Ungulates (brown/earth)
    "Pig": "#6E2C00",  # Dark brown
}
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
    regions_of_interest     = None,
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
# STEP 6  —  Generate all figures
# ═══════════════════════════════════════════════════════════════════════════════

from Exemples.Study_EDNiX.analysis import EDNiX_figures
atlas_name = "EDNIxCSC",
corr_atlas_level = 2
data = data_plot
species_config = species_config
all_bids = ALL_BIDS
global_bids_col = None
atlas_level = ATLAS_LEVEL
atlas_name = ATLAS_NAME
average_sessions = False


df_surface = data.get("surface")
df_thickness = data.get("thickness")
df_qc = data.get("qc")

# Bilateral volume aggregation (always on — sums L+R per subject-region).
# Session averaging is optional; off by default so each session/run
# appears as a separate dot (pass average_sessions=True to collapse).
df_volume = EDNiX_figures._make_bilateral_volume(data.get("volume"))
if average_sessions:
    df_volume = _mean_per_subject(df_volume, "volume_mm3")
    df_surface = _mean_per_subject(df_surface, "surface_area_mm2")
    df_thickness = _mean_per_subject(df_thickness, "thickness_mm")

if global_bids_col is None:
    all_bd = sorted(set(EDNiX_figures._bids_label(bd) for _, bd in all_bids))
    global_bids_col = {bd: PALETTE[i % len(PALETTE)] for i, bd in enumerate(all_bd)}

hemis_surface = EDNiX_figures._hemis_for(df_surface)
hemis_thickness = EDNiX_figures._hemis_for(df_thickness)
hemis_volume = ["bilateral"]  # always bilateral after fix

fig_dir = opj(OUT_DIR, 'figures')
morph_dir = opj(fig_dir, "cross_species", "morphometry_recap")
combo_dir = opj(fig_dir, "cross_species", "combo")
qc_dir = opj(fig_dir, "cross_species", "qc_recap")
scaling_dir = opj(fig_dir, "cross_species", "brain_scaling")
for d in [morph_dir, combo_dir, qc_dir, scaling_dir]:
    os.makedirs(d, exist_ok=True)
# ── Ground truth from dlabel vertex data (Rat sub-300101) ────────────────────
# Only regions that actually have vertices on the pial surface are included.
# Subcortical spillover (Brainstem 1vtx, Claustrum 2vtx, Thalamus 8vtx etc.)
# kept/excluded based on biological relevance, not just vertex count.

level1_cortical = [
    "Isocortex",
    #"Allocortex",
    #"Periallocortex",
    # EXCLUDED: Cortical_White_matter (2 vtx, artifact)
    # EXCLUDED: Subcortical_areas (201 vtx, not cortical surface)
]

level2_cortical = [
    "Frontal_lobe",
    "Occipital_lobe",
    "Temporal_lobe",
    "Parietal_lobe",
    "Insular_lobe",
    #"Paleocortex",
    #"Archicortex",
    #"Periarchicortex",
    # EXCLUDED: Brainstem (388 vtx — deep structure spillover)
    # EXCLUDED: Cerebral_nuclei (30 vtx — subcortical spillover)
    # EXCLUDED: Cortical_subplate (106 vtx — amygdala/claustrum, deep)
    # NOTE: Parietal_lobe absent from Rat map 2 vertex data
]

level3_cortical = [
    "Motor_premotor",
    "PFC",
    "Somatosensory_Cx",
    "Post_parietal_Cx",
    "Striate_Cx",
    "Extra_striate_Cx",
    "Auditory_Cx",
    "Insula",
    "Olfactory_Cx",
    #"Hippocampal_form",
    #"Subiculum",
    #"Periarchicortex",   # still present at map 3 (465 vtx)
    # EXCLUDED: oFC, PMC, MIPT, Vent_temporal — absent from Rat map 3
    # EXCLUDED: Brainstem (1 vtx), Basal_forebrain (26), Pallidum (5),
    #           Claustrum (2), Amygdala (79 — deep), Thalamus (8)
    # UNKNOWN_51, UNKNOWN_65, UNKNOWN_70 — need alias table update
]

level4_cortical = [
    "BA_6",
    "BA_4_M1",
    "BA_1_3_S1",
    "S2",
    "BA_18_V2",
    "Belt_A2",
    "Insular_Cx",
    "BA_27",
    "BA_35_36",
    #"Presubiculum",
    #"BA_26_29_30",
    # EXCLUDED: CA (54 vtx — deep hippocampus)
    # EXCLUDED: DG (292 vtx — deep hippocampus)
    # EXCLUDED: Subiculum (248 vtx — deep)
    # EXCLUDED: Amygdala (109 vtx — deep)
    # EXCLUDED: Claustrum (2 vtx — deep)
    # EXCLUDED: Subthalamus (11 vtx — subcortical)
    # EXCLUDED: Globus_pallidus (6 vtx — subcortical)
]

########################## Figure 1 ##########################
atlas_level_l  = [2, 3, 4]
plot_regions_l = [level2_cortical, level3_cortical, level4_cortical]

for d in [morph_dir, combo_dir]:
    os.makedirs(d, exist_ok=True)

fig_dir_l = ['figures', 'supplementary', 'supplementary']
morph_dir_l = ["other", "Suppl1", "other"]
combo_dir_l = ["Figure2", "Suppl1",  "other"]

for atlas_level, plot_regions, fig_dir_, morph_dir_, combo_dir_ in zip(atlas_level_l, plot_regions_l, fig_dir_l, morph_dir_l, combo_dir_l):

    fig_dir = opj(OUT_DIR, fig_dir_)
    morph_dir = opj(fig_dir, morph_dir_)
    combo_dir = opj(fig_dir, combo_dir_)

    # ── Volume: bilateral only ────────────────────────────────────────────────
    for log_s, suffix in [(False, "raw"), (True, "log")]:
        EDNiX_figures.plot_cross_species_dots(
            df_volume, "volume_mm3", "Volume (mm³)", plot_regions,
            opj(morph_dir, f"lvl_{str(atlas_level)}_volume_bilateral_{suffix}.png"),
            hemisphere="bilateral", log_scale=log_s, atlas_level=atlas_level,
            bids_col=global_bids_col, anesth_map=anesth_map)
    EDNiX_figures.plot_cross_species_dots(
        df_volume, "volume_mm3", "Volume (mm³)", plot_regions,
        opj(morph_dir, f"lvl_{str(atlas_level)}_volume_bilateral_norm.png"),
        hemisphere="bilateral", log_scale=False, normalise_by_brain=True,
        df_vol_for_norm=df_volume, atlas_level=atlas_level,
        bids_col=global_bids_col, anesth_map=anesth_map)

    # ── Surface + thickness: all available hemispheres ────────────────────────
    for hemi in sorted(set(hemis_surface + hemis_thickness)):
        for df, metric_col, metric_label, mod, hemis in [
            (df_surface, "surface_area_mm2", "Surface area (mm²)", "surface", hemis_surface),
            (df_thickness, "thickness_mm", "Cortical thickness (mm)", "thickness", hemis_thickness),
        ]:
            if df is None or df.empty or hemi not in hemis: continue
            for log_s, suffix in [(False, "raw"), (True, "log")]:
                EDNiX_figures.plot_cross_species_dots(
                    df, metric_col, metric_label, plot_regions,
                    opj(morph_dir, f"lvl_{str(atlas_level)}_{mod}_{hemi}_{suffix}.png"),
                    hemisphere=hemi, log_scale=log_s, atlas_level=atlas_level,
                    bids_col=global_bids_col, anesth_map=anesth_map)
            if df_volume is not None:
                EDNiX_figures.plot_cross_species_dots(
                    df, metric_col, metric_label, plot_regions,
                    opj(morph_dir, f"lvl_{str(atlas_level)}_{mod}_{hemi}_norm.png"),
                    hemisphere=hemi, log_scale=False, normalise_by_brain=True,
                    df_vol_for_norm=df_volume, atlas_level=atlas_level,
                    bids_col=global_bids_col, anesth_map=anesth_map)

    # ── Combo ─────────────────────────────────────────────────────────────────
    for log_s, norm_s, suffix in [(False, False, "raw"), (True, False, "log"), (False, True, "norm")]:
        EDNiX_figures.plot_combo(df_surface, df_volume, df_thickness,
                   opj(combo_dir, f"lvl_{str(atlas_level)}_combo_bilateral_{suffix}.png"),
                   regions=plot_regions, hemisphere="bilateral",
                   log_scale=log_s, normalise_by_brain=norm_s,
                   atlas_level=atlas_level, bids_col=global_bids_col,
                   anesth_map=anesth_map)
# ── Cross-modality correlations (Surface vs Volume, Thickness vs Surface/Volume)
# One mean per species → Pearson R + 95% CI (Fisher z) logged to paper_values.csv
from scipy import stats as _st
import numpy as np

def _pearson_ci(x, y):
    # Species means span several orders of magnitude: on raw values Pearson r is
    # driven by the largest brains (leverage). Correlate log10 values instead
    # (allometric convention). Species are not independent: see the PGLS section.
    x = np.log10(np.asarray(x, float)); y = np.log10(np.asarray(y, float))
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 4: return None
    r, p = _st.pearsonr(x, y)
    z, se = np.arctanh(r), 1.0 / np.sqrt(n - 3)
    zc = _st.norm.ppf(0.975)
    return dict(r=round(float(r),4), p=round(float(p),8),
                ci_lo=round(float(np.tanh(z - zc*se)),4),
                ci_hi=round(float(np.tanh(z + zc*se)),4), n=n)

corr_rows = []
for region in level2_cortical:
    def _reg_mean(df, col):
        """Species mean of the bilateral value of `region` at CORR_ATLAS_LEVEL."""
        if df is None or df.empty: return {}
        d = df[df["region"] == region]
        if "atlas_level" in d.columns:
            d = d[d["atlas_level"] == CORR_ATLAS_LEVEL]
        if "hemisphere" in d.columns:
            d = d[d["hemisphere"] == "bilateral"]
        return d.groupby("species")[col].mean().to_dict()

    ms = _reg_mean(df_surface,   "surface_area_mm2")
    mv = _reg_mean(df_volume,    "volume_mm3")
    mt = _reg_mean(df_thickness, "thickness_mm")

    for (mx, my, xcol, ycol, pair) in [
        (ms, mv, "surface_area_mm2", "volume_mm3",    "surf_vs_vol"),
        (ms, mt, "surface_area_mm2", "thickness_mm",  "surf_vs_thick"),
        (mv, mt, "volume_mm3",       "thickness_mm",  "vol_vs_thick"),
    ]:
        common = sorted(set(mx) & set(my))
        if len(common) < 4: continue
        xs = np.array([mx[s] for s in common])
        ys = np.array([my[s] for s in common])
        res = _pearson_ci(xs, ys)
        if res:
            print(f"  [corr] lvl{CORR_ATLAS_LEVEL} {region} {pair}: "
                  f"R={res['r']}  p={res['p']}  "
                  f"95%CI[{res['ci_lo']}, {res['ci_hi']}]  n={res['n']}")
            corr_rows.append(dict(atlas_level=CORR_ATLAS_LEVEL, region=region,
                                 pair=pair, **res))

if corr_rows:
    import pandas as pd
    corr_df = pd.DataFrame(corr_rows)
    # Benjamini-Hochberg over all region x pair tests
    from statsmodels.stats.multitest import multipletests as _mt
    corr_df["p_fdr"] = _mt(corr_df["p"].values, method="fdr_bh")[1]
    corr_df["scale"] = "log10"
    corr_out = opj(OUT_DIR, "stats", f"correlations_lvl{CORR_ATLAS_LEVEL}.csv")
    os.makedirs(os.path.dirname(corr_out), exist_ok=True)
    corr_df.to_csv(corr_out, index=False)
    print(f"  [corr] → {corr_out}")


# ── Regional outlier rates, datasets pooled within species (sensitivity check of
#    Supplementary Figs 6-8, which use within-dataset allometric residuals)
for atlas_level, plot_regions, fig_dir_, morph_dir_, combo_dir_ in zip(atlas_level_l, plot_regions_l, fig_dir_l,
                                                                           morph_dir_l, combo_dir_l):
    # MAD-based: outlier = |value - median| / (1.4826 * MAD) > 2.5
    # One rate per region per species, aggregated across species for the paper sentence.
    # ── Regional outlier rates (Supplementary Fig. 7-9)
    # MAD-based: outlier = |value - median| / (1.4826 * MAD) > 2.5
    # One rate per region per species, aggregated across species for the paper sentence.
    outlier_rows = []
    for df, col, mod in [
        (df_surface, "surface_area_mm2", "surface"),
        (df_volume, "volume_mm3", "volume"),
        (df_thickness, "thickness_mm", "thickness"),
    ]:
        if df is None or df.empty or col not in df.columns:
            continue
        d = df[df["atlas_level"] == atlas_level].copy() if "atlas_level" in df.columns else df.copy()
        if "hemisphere" in d.columns:
            bil = d[d["hemisphere"] == "bilateral"]
            d = bil if not bil.empty else d

        for (sp, region), grp in d.groupby(["species", "region"]):
            vals = grp[col].dropna().values
            vals = vals[np.isfinite(vals)]
            if len(vals) < 5:
                continue
            med = np.median(vals)
            mad = np.median(np.abs(vals - med))
            if mad == 0:
                continue
            n_out = int(np.sum(np.abs(vals - med) / (1.4826 * mad) > 2.5))
            outlier_rows.append(dict(
                atlas_level=atlas_level,
                modality=mod,
                species=sp,
                region=region,
                n=len(vals),
                n_outliers=n_out,
                outlier_rate=round(n_out / len(vals), 4),
            ))

    if outlier_rows:
        out_df = pd.DataFrame(outlier_rows)

        # ── Per-modality supplementary figure ────────────────────────────────
        # ── Per-modality supplementary figure ────────────────────────────────
        for mod in ("surface", "volume", "thickness"):
            sub = out_df[out_df["modality"] == mod]
            if sub.empty:
                continue

            # Regions sorted by mean outlier rate across species
            reg_order = (sub.groupby("region")["outlier_rate"]
                         .mean()
                         .sort_values(ascending=False)
                         .index.tolist())
            species_list = sorted(sub["species"].unique(),
                                  key=lambda s: PHYLO_ORDER.index(s)
                                  if s in PHYLO_ORDER else 99)
            n_reg = len(reg_order)
            n_sp = len(species_list)

            fig, ax = plt.subplots(figsize=(max(6, n_reg * 0.65), max(4, n_sp * 0.55) + 1.5))
            x = np.arange(n_reg)
            bar_w = 0.8 / max(n_sp, 1)

            for si, sp in enumerate(species_list):
                sp_sub = sub[sub["species"] == sp].set_index("region")
                rates = [sp_sub.loc[r, "outlier_rate"] * 100
                         if r in sp_sub.index else 0.0
                         for r in reg_order]
                offset = (si - n_sp / 2 + 0.5) * bar_w
                ax.bar(x + offset, rates, bar_w,
                       color=SPECIES_COLORS.get(sp, "#888"),
                       edgecolor="none", alpha=0.85, label=sp)

            ax.set_xticks(x)
            ax.set_xticklabels(reg_order, rotation=45, ha="right", fontsize=7)
            ax.set_ylabel("Outlier rate (%)", fontsize=10)
            ax.set_title(f"Regional outlier rates — {mod}  (atlas level {atlas_level})",
                         fontsize=10, fontweight="bold")
            ax.axhline(5, color="#e67e22", lw=0.8, ls="--", alpha=0.6)
            ax.axhline(10, color="#c0392b", lw=0.8, ls="--", alpha=0.6)
            ax.text(n_reg - 0.5, 5.3, ">5%", fontsize=7, color="#e67e22", ha="right")
            ax.text(n_reg - 0.5, 10.3, ">10%", fontsize=7, color="#c0392b", ha="right")
            ax.legend(fontsize=7, frameon=False,
                      bbox_to_anchor=(1.01, 1), loc="upper left", title="Species")
            for sp_ in ("top", "right"):
                ax.spines[sp_].set_visible(False)
            plt.tight_layout()

            fig_out = opj(OUT_DIR, "supplementary",
                          f"outliers_raw_pooled_{mod}_lvl{atlas_level}.png")
            os.makedirs(os.path.dirname(fig_out), exist_ok=True)
            fig.savefig(fig_out, dpi=200, bbox_inches="tight")
            plt.close(fig)
            print(f"  [outliers] → {fig_out}")
