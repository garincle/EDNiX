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

def _species_mean(df, col):
    if df is None or df.empty: return {}
    d = df[df["atlas_level"] == atlas_level] if "atlas_level" in df.columns else df
    h = d[d["hemisphere"] == "bilateral"] if "hemisphere" in d.columns else d
    return h.groupby("species")[col].mean().to_dict()

corr_rows = []
for region in level2_cortical:
    sv = {sp: v for sp, v in _species_mean(df_surface, "surface_area_mm2").items()
          if region in df_surface[df_surface["region"] == region]["region"].values}
    # simpler: filter per region
    def _reg_mean(df, col):
        if df is None or df.empty: return {}
        d = df[(df["atlas_level"] == 2) & (df["region"] == region)] \
            if "atlas_level" in df.columns else df[df["region"] == region]
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
            print(f"  [corr] lvl{atlas_level} {region} {pair}: "
                  f"R={res['r']}  p={res['p']}  "
                  f"95%CI[{res['ci_lo']}, {res['ci_hi']}]  n={res['n']}")
            corr_rows.append(dict(atlas_level=atlas_level, region=region,
                                 pair=pair, **res))

if corr_rows:
    import pandas as pd
    corr_df = pd.DataFrame(corr_rows)
    # Benjamini-Hochberg over all region x pair tests
    from statsmodels.stats.multitest import multipletests as _mt
    corr_df["p_fdr"] = _mt(corr_df["p"].values, method="fdr_bh")[1]
    corr_df["scale"] = "log10"
    corr_out = opj(OUT_DIR, "stats", f"correlations_lvl{atlas_level}.csv")
    os.makedirs(os.path.dirname(corr_out), exist_ok=True)
    corr_df.to_csv(corr_out, index=False)
    print(f"  [corr] → {corr_out}")


# ── Regional outlier rates (Supplementary Fig. 7-9)
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

            suppl_n = {"surface": 7, "volume": 8, "thickness": 9}[mod]
            fig_out = opj(OUT_DIR, "supplementary",
                          f"Suppl{suppl_n}_outliers_{mod}_lvl{atlas_level}.png")
            os.makedirs(os.path.dirname(fig_out), exist_ok=True)
            fig.savefig(fig_out, dpi=200, bbox_inches="tight")
            plt.close(fig)
            print(f"  [outliers] → {fig_out}")

"""
PGLS — Evolutionary comparative analyses
=========================================
Two analyses:

  A) Frontal lobe scaling vs body weight, brain weight, neuron count (Fig. 2B)
     → 3 separate figures: surface, volume, thickness (each with 3 panels)

  B) Frontal vs Occipital lobe (Fig. 2C)
     → 1 figure with 3 panels: surface, volume, thickness

For each modality (surface, volume, thickness):
  - species-mean values at level 2 (Frontal / Occipital)
  - PGLS with Pagel's λ estimated by ML
  - Plots + stats exported
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as tck

import sys

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/Exemples/Study_EDNiX/analysis')

from Statistics.Evo.PGLS import (
    load_newick_vcv,
    pgls_regression,
    pgls_on_ax,
    PAPER_RC,
    _phylo_sort,
    _species_colors,
    _pgls_on_ax_with_prediction
)

# ── Output directories ────────────────────────────────────────────────────────
pgls_dir = opj(OUT_DIR, 'figures', 'cross_species', 'brain_scaling')
pgls_stats_dir = opj(OUT_DIR, 'stats')
os.makedirs(pgls_dir, exist_ok=True)
os.makedirs(pgls_stats_dir, exist_ok=True)

print("\n" + "=" * 64)
print("  PGLS — Evolutionary comparative analyses")
print("=" * 64)

# ── Load phylogeny ──────────────────────────────────────────────────────────
vcv_df, _ = load_newick_vcv(NEWICK_FILE)


# ── Helpers ─────────────────────────────────────────────────────────────────
def _species_mean(df, metric_col, region, atlas_level, hemisphere="bilateral"):
    """Return {species: mean} for a region."""
    if df is None or df.empty:
        return {}
    d = df.copy()
    if "atlas_level" in d.columns:
        d = d[d["atlas_level"] == atlas_level]
    if "hemisphere" in d.columns and hemisphere != "all":
        h = d[d["hemisphere"] == hemisphere]
        d = h if not h.empty else d[d["hemisphere"] == "bilateral"]
    d = d[d["region"] == region]
    return {} if d.empty else d.groupby("species")[metric_col].mean().to_dict()


def _load_herculano_from_excel(xlsx_path):
    """
    Load body weight, brain weight, neuron count from Herculano-Houzel Excel file.
    Uses T5 – Whole Brain sheet primarily.
    """
    # EDNiX species to Excel species mapping
    SPECIES_MAP = {
        "Mouse": ["Mus musculus"],
        "Rat": ["Rattus norvegicus"],
        "Mouselemur": ["Microcebus murinus"],
        "Marmoset": ["Callithrix jacchus"],
        "Macaque": ["Macaca mulatta", "Macaca fascicularis", "Macaca radiata"],
        "Human": ["Homo sapiens"],
        "Dog": ["Canis familiaris"],
        "Cat": ["Felis catus"],
        "Bat": ["Pteropus vampyrus"],
        "Pig": ["Sus scrofa domesticus"],
    }

    body_g, brain_g, neuron_M = {}, {}, {}

    if not os.path.exists(xlsx_path):
        print(f"  [Herculano] File not found: {xlsx_path}")
        return body_g, brain_g, neuron_M

    # ── Read T5 – Whole Brain ──────────────────────────────────────────────────
    try:
        df_raw = pd.read_excel(xlsx_path, sheet_name="T5 – Whole Brain", header=None)

        header_row = None
        for idx, row in df_raw.iterrows():
            row_str = ' '.join([str(x).strip() for x in row.values if pd.notna(x)])
            if 'Species' in row_str and 'Body mass' in row_str and 'Brain mass' in row_str:
                header_row = idx
                break

        if header_row is not None:
            df_t5 = pd.read_excel(xlsx_path, sheet_name="T5 – Whole Brain", header=header_row)
            df_t5.columns = [str(c).replace('\n', ' ').strip() for c in df_t5.columns]

            # Print columns for debugging
            print(f"  [Herculano] T5 columns: {df_t5.columns.tolist()}")

            sp_col = None
            body_col = None
            brain_col = None
            neuron_col = None

            for col in df_t5.columns:
                col_clean = col.lower().strip()
                if 'species' in col_clean:
                    sp_col = col
                elif 'body mass' in col_clean:
                    body_col = col
                elif 'brain mass' in col_clean:
                    brain_col = col
                elif 'neurons' in col_clean and 'mg' not in col_clean and '%' not in col_clean:
                    neuron_col = col

            print(f"  [Herculano] neuron_col = '{neuron_col}'")

            if sp_col and body_col and brain_col and neuron_col:
                for idx, row in df_t5.iterrows():
                    sp_name = str(row[sp_col]).strip().lower() if pd.notna(row[sp_col]) else ""
                    if not sp_name:
                        continue

                    for ednix_name, excel_names in SPECIES_MAP.items():
                        if any(excel_name.lower() in sp_name for excel_name in excel_names):
                            # Body mass
                            if body_col and pd.notna(row[body_col]) and ednix_name not in body_g:
                                try:
                                    val_str = str(row[body_col]).replace(',', '').strip()
                                    if '±' in val_str:
                                        val_str = val_str.split('±')[0]
                                    val = float(val_str)
                                    body_g[ednix_name] = val
                                except:
                                    pass

                            # Brain mass
                            if brain_col and pd.notna(row[brain_col]) and ednix_name not in brain_g:
                                try:
                                    val_str = str(row[brain_col]).replace(',', '').strip()
                                    if '±' in val_str:
                                        val_str = val_str.split('±')[0]
                                    val = float(val_str)
                                    brain_g[ednix_name] = val
                                except:
                                    pass

                            # Raw neuron count - THIS IS THE KEY
                            if neuron_col and pd.notna(row[neuron_col]) and ednix_name not in neuron_M:
                                try:
                                    val_str = str(row[neuron_col]).replace(',', '').strip()
                                    print(f"  [DEBUG] {ednix_name}: neuron_str = '{val_str}'")
                                    if '±' in val_str:
                                        val_str = val_str.split('±')[0]
                                    # Remove ×10⁶ if present
                                    if '×10⁶' in val_str:
                                        val_str = val_str.replace('×10⁶', '')
                                    val = float(val_str)
                                    # If value is > 1e6, it's raw count, convert to millions
                                    if val > 1e6:
                                        val = val / 1e6
                                    neuron_M[ednix_name] = val
                                except Exception as e:
                                    print(f"  [DEBUG] Failed to parse {ednix_name}: {e}")
                            break

    except Exception as e:
        print(f"  [Herculano] Error reading T5: {e}")
        import traceback
        traceback.print_exc()

    # ── Read Carnivora sheet for Dog/Cat ──────────────────────────────────────
    try:
        df_carn = pd.read_excel(xlsx_path, sheet_name="Carnivora (2017)", header=1)
        df_carn.columns = [str(c).replace('\n', ' ').strip() for c in df_carn.columns]

        sp_col = None
        body_col = None
        brain_col = None
        neuron_col = None

        for col in df_carn.columns:
            col_clean = col.lower().strip()
            if 'species' in col_clean:
                sp_col = col
            elif 'body mass' in col_clean:
                body_col = col
            elif 'brain mass' in col_clean:
                brain_col = col
            elif 'n whole brain' in col_clean:
                neuron_col = col

        if sp_col and body_col and brain_col and neuron_col:
            for idx, row in df_carn.iterrows():
                sp_name = str(row[sp_col]).strip().lower() if pd.notna(row[sp_col]) else ""
                if 'dog' in sp_name and 'Dog' not in body_g:
                    body_g['Dog'] = float(str(row[body_col]).replace(',', '').strip()) * 1000
                    brain_g['Dog'] = float(str(row[brain_col]).replace(',', '').strip())
                    val_str = str(row[neuron_col]).replace('×10⁶', '').replace(',', '').strip()
                    neuron_M['Dog'] = float(val_str)
                elif 'cat' in sp_name and 'Cat' not in body_g:
                    body_g['Cat'] = float(str(row[body_col]).replace(',', '').strip()) * 1000
                    brain_g['Cat'] = float(str(row[brain_col]).replace(',', '').strip())
                    val_str = str(row[neuron_col]).replace('×10⁶', '').replace(',', '').strip()
                    neuron_M['Cat'] = float(val_str)
    except Exception as e:
        print(f"  [Herculano] Carnivora sheet error: {e}")

    # ── Read Chiroptera sheet for Bat ─────────────────────────────────────────
    try:
        df_bat = pd.read_excel(xlsx_path, sheet_name="Chiroptera - Bats (2020)", header=1)
        df_bat.columns = [str(c).replace('\n', ' ').strip() for c in df_bat.columns]

        sp_col = None
        body_col = None
        brain_col = None
        neuron_col = None

        for col in df_bat.columns:
            col_clean = col.lower().strip()
            if 'species' in col_clean:
                sp_col = col
            elif 'body mass' in col_clean:
                body_col = col
            elif 'brain mass' in col_clean:
                brain_col = col
            elif 'n whole brain' in col_clean:
                neuron_col = col

        if sp_col and body_col and brain_col and neuron_col:
            for idx, row in df_bat.iterrows():
                sp_name = str(row[sp_col]).strip().lower() if pd.notna(row[sp_col]) else ""
                if 'pteropus vampyrus' in sp_name and 'Bat' not in body_g:
                    body_g['Bat'] = float(str(row[body_col]).replace(',', '').strip())
                    brain_g['Bat'] = float(str(row[brain_col]).replace(',', '').strip())
                    val_str = str(row[neuron_col]).replace('×10⁶', '').replace(',', '').strip()
                    neuron_M['Bat'] = float(val_str)
    except Exception as e:
        print(f"  [Herculano] Chiroptera sheet error: {e}")

    # ── Report ─────────────────────────────────────────────────────────────────
    all_species = ["Mouse", "Rat", "Mouselemur", "Marmoset", "Macaque", "Human", "Dog", "Cat", "Bat", "Pig"]

    print(f"\n  [Herculano] Loaded from Excel (raw neuron counts in millions):")
    print(f"    body_g:    {len(body_g)} species")
    print(f"    brain_g:   {len(brain_g)} species")
    print(f"    neuron_M:  {len(neuron_M)} species")

    missing_body = [s for s in all_species if s not in body_g]
    missing_brain = [s for s in all_species if s not in brain_g]
    missing_neuron = [s for s in all_species if s not in neuron_M]

    if missing_body:
        print(f"\n  ⚠️  MISSING body_g: {missing_body}")
    if missing_brain:
        print(f"  ⚠️  MISSING brain_g: {missing_brain}")
    if missing_neuron:
        print(f"  ⚠️  MISSING neuron_M: {missing_neuron}")

    print("\n  [Herculano] Values (whole brain neurons in millions):")
    print(f"    {'Species':12s} {'Body (g)':>10s} {'Brain (g)':>10s} {'Neurons (M)':>12s}")
    print(f"    {'-' * 12} {'-' * 10} {'-' * 10} {'-' * 12}")
    for s in all_species:
        b = body_g.get(s, 'MISSING')
        br = brain_g.get(s, 'MISSING')
        n = neuron_M.get(s, 'MISSING')
        b_str = f"{b:.0f}" if isinstance(b, (int, float)) else str(b)
        br_str = f"{br:.2f}" if isinstance(br, (int, float)) else str(br)
        n_str = f"{n:.0f}" if isinstance(n, (int, float)) else str(n)
        print(f"    {s:12s} {b_str:>10s} {br_str:>10s} {n_str:>12s}")

    return body_g, brain_g, neuron_M

# ── Load Herculano-Houzel data from Excel ──────────────────────────────────
HERCULANO_XLSX = "/scratch2/EDNiX/Herculano-Houzel_mammalian_brain_dataset_complete.xlsx"
body_g, brain_g, neuron_M = _load_herculano_from_excel(HERCULANO_XLSX)

# ── Modalities ──────────────────────────────────────────────────────────────
MODALITIES = [
    ("surface", df_surface_plot, "surface_area_mm2", "Surface area (mm²)"),
    ("volume", df_volume_plot, "volume_mm3", "Volume (mm³)"),
    ("thickness", df_thickness_plot, "thickness_mm", "Thickness (mm)"),
]

# ── Results storage ─────────────────────────────────────────────────────────
all_results = []

# ═══════════════════════════════════════════════════════════════════════════
# A) Frontal lobe scaling — 3 separate figures (surface, volume, thickness)
#    Each figure: 1x3 panels (body, brain, neuron)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 64)
print("  A) Frontal lobe scaling — body / brain / neuron (Fig. 2B)")
print("─" * 64)

predictors = [
    ("body_weight", "Body weight (g)", body_g),
    ("brain_weight", "Brain weight (g)", brain_g),
    ("neuron_count", "Neuron count (×10⁶)", neuron_M),
]

for mod, df_mod, metric_col, metric_label in MODALITIES:
    if df_mod is None or df_mod.empty:
        print(f"  [{mod}] no data — skipped")
        continue

    # Get Frontal lobe means
    means = _species_mean(df_mod, metric_col, "Frontal_lobe",
                          atlas_level=2, hemisphere="bilateral")
    if len(means) < 3:
        print(f"  [{mod}] <3 species with Frontal data — skipped")
        continue

    sp_list = sorted(means.keys())
    y_vals = np.array([means[s] for s in sp_list])
    print(f"\n  [{mod}] n={len(sp_list)} species: {sp_list}")

    # Create figure: 1 row, 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Use the SPECIES_COLORS defined in the main script
    sp_col = SPECIES_COLORS

    for col, (pred_name, x_label, xdata) in enumerate(predictors):
        ax = axes[col]
        x_vals = np.array([xdata.get(s, np.nan) for s in sp_list])

        # Filter finite values
        valid = [(s, xv, yv) for s, xv, yv in zip(sp_list, x_vals, y_vals)
                 if np.isfinite(xv) and np.isfinite(yv)]
        if len(valid) < 3:
            ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10)
            ax.set_xscale("log")
            ax.set_yscale("log")
            continue

        sp_sub, x_sub, y_sub = zip(*valid)
        x_arr = np.array(x_sub)
        y_arr = np.array(y_sub)

        # Scatter
        for s, xv, yv in zip(sp_sub, x_arr, y_arr):
            ax.scatter(xv, yv, color=sp_col.get(s, "#888"), s=55, zorder=5,
                       edgecolors="k", linewidths=0.5)

        # PGLS line with PREDICTION intervals (wider than confidence intervals)
        res = pgls_regression(x_arr, y_arr, vcv_df, list(sp_sub),
                              log_x=True, log_y=True)
        if res:
            # Use the modified pgls_on_ax with prediction intervals
            _pgls_on_ax_with_prediction(ax, res, log_x=True, log_y=True)

            all_results.append(dict(
                modality=mod, predictor=pred_name, **res,
                region="Frontal_lobe", atlas_level=2
            ))

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(x_label, fontsize=10)
        ax.xaxis.set_major_formatter(tck.LogFormatterSciNotation())
        ax.yaxis.set_major_formatter(tck.LogFormatterSciNotation())

        # Only first column gets y-label
        if col == 0:
            ax.set_ylabel(f"Frontal lobe {metric_label}", fontsize=10)

    # Legend
    handles = [mpatches.Patch(facecolor=sp_col.get(s, "#888"), label=s)
               for s in _phylo_sort(sp_list)]
    fig.legend(handles=handles, fontsize=7,
               bbox_to_anchor=(0.5, -0.02), loc="upper center", ncol=5)

    plt.tight_layout()
    out_png = opj(pgls_dir, f"Fig2B_frontal_scaling_{mod}.png")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] {out_png}")

# ═══════════════════════════════════════════════════════════════════════════
# B) Frontal vs Occipital — 1 figure with 3 panels (surface, volume, thickness)
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 64)
print("  B) Frontal vs Occipital (Fig. 2C)")
print("─" * 64)

# Collect data for all modalities
occipital_data = {}
for mod, df_mod, metric_col, metric_label in MODALITIES:
    if df_mod is None or df_mod.empty:
        continue

    means_frontal = _species_mean(df_mod, metric_col, "Frontal_lobe",
                                  atlas_level=2, hemisphere="bilateral")
    means_occipital = _species_mean(df_mod, metric_col, "Occipital_lobe",
                                    atlas_level=2, hemisphere="bilateral")

    common_sp = sorted(set(means_frontal) & set(means_occipital) & set(vcv_df.index))
    if len(common_sp) < 3:
        print(f"  [{mod}] <3 species with both Frontal & Occipital — skipped")
        continue

    occipital_data[mod] = {
        'x': means_occipital,
        'y': means_frontal,
        'species': common_sp,
        'label': metric_label
    }

if occipital_data:
    # Create 1x3 grid
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sp_col = SPECIES_COLORS

    for idx, (mod, data) in enumerate(occipital_data.items()):
        ax = axes[idx]
        sp_list = data['species']
        x_occ = np.array([data['x'][s] for s in sp_list])
        y_frontal = np.array([data['y'][s] for s in sp_list])

        # Scatter
        for s, xv, yv in zip(sp_list, x_occ, y_frontal):
            ax.scatter(xv, yv, color=sp_col.get(s, "#888"), s=55, zorder=5,
                       edgecolors="k", linewidths=0.5)

        # PGLS with prediction intervals
        res = pgls_regression(x_occ, y_frontal, vcv_df, sp_list,
                              log_x=True, log_y=True)
        if res:
            _pgls_on_ax_with_prediction(ax, res, log_x=True, log_y=True)
            all_results.append(dict(
                modality=mod, predictor="occipital_lobe", **res,
                region_x="Occipital_lobe", region_y="Frontal_lobe", atlas_level=2
            ))

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(f"Occipital lobe {data['label']}", fontsize=10)
        ax.xaxis.set_major_formatter(tck.LogFormatterSciNotation())
        ax.yaxis.set_major_formatter(tck.LogFormatterSciNotation())

        # Only first column gets y-label
        if idx == 0:
            ax.set_ylabel(f"Frontal lobe {data['label']}", fontsize=10)

    # Legend
    first_sp = list(occipital_data.values())[0]['species']
    handles = [mpatches.Patch(facecolor=sp_col.get(s, "#888"), label=s)
               for s in _phylo_sort(first_sp)]
    fig.legend(handles=handles, fontsize=7,
               bbox_to_anchor=(0.5, -0.02), loc="upper center", ncol=5)

    plt.tight_layout()
    out_png = opj(pgls_dir, "Fig2C_frontal_vs_occipital_all_modalities.png")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [plot] {out_png}")

# ═══════════════════════════════════════════════════════════════════════════
# Export stats
# ═══════════════════════════════════════════════════════════════════════════

if all_results:
    df_results = pd.DataFrame(all_results)
    out_csv = opj(pgls_stats_dir, "pgls_frontal_scaling.csv")
    df_results.to_csv(out_csv, index=False)

    # Print summary table
    print("\n" + "─" * 64)
    print("  SUMMARY — PGLS results")
    print("─" * 64)
    summary_cols = ["modality", "predictor", "slope", "r_pgls", "p", "lambda_pagel", "n_species"]
    summary = df_results[summary_cols]
    print(summary.to_string(index=False))
    print(f"\n  [stats] {out_csv}")

print("\n  PGLS done.")
print(f"  Figures → {pgls_dir}")
print(f"  Stats   → {pgls_stats_dir}")
