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

import pandas as pd
import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
    PALETTE, get_atlas_label_path, find_species_path, run_wb_extraction,
)

opj = os.path.join

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

ATLAS_LIB      = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
ATLAS_NAME     = 'EDNIxCSC'
ATLAS_LEVEL    = 1   # hierarchy level shown in figures (1–4)
N_ATLAS_LEVELS = 4   # how many levels Stage 0 extracts — must be >= 1

OUT_DIR = '/scratch2/EDNiX/results/multispecies_analysis'
FIG_DIR = opj(OUT_DIR, 'figures')

SING_WB = (
    'vglrun singularity run '
    '--bind /srv/projects/,/srv/projects/easymribrain,/scratch2/,/scratch/ '
    '/home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/'
    'connectome_workbench_1.5.0-freesurfer-update.sif'
)

level1_cortical = [
    "Isocortex",
    # "Allocortex",
    # "Periallocortex",
    # EXCLUDED: Cortical_White_matter (2 vtx, artifact)
    # EXCLUDED: Subcortical_areas (201 vtx, not cortical surface)
]

level2_cortical = [
    "Frontal_lobe",
    "Occipital_lobe",
    "Temporal_lobe",
    "Parietal_lobe",
    "Insular_lobe",
    # "Paleocortex",
    # "Archicortex",
    # "Periarchicortex",
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
    # "Hippocampal_form",
    # "Subiculum",
    # "Periarchicortex",   # still present at map 3 (465 vtx)
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
    "Striate_Cx",  # still present at map 4
    "Auditory_Cx",  # still present at map 4
    "BA_18_V2",
    "Belt_A2",
    "Insular_Cx",
    "BA_27",
    "Olfactory_bulb",  # 20 vtx — present at map 4 for Rat
    "BA_35_36",
    # "Presubiculum",
    # "BA_26_29_30",
    # EXCLUDED: CA (54 vtx — deep hippocampus)
    # EXCLUDED: DG (292 vtx — deep hippocampus)
    # EXCLUDED: Subiculum (248 vtx — deep)
    # EXCLUDED: Amygdala (109 vtx — deep)
    # EXCLUDED: Claustrum (2 vtx — deep)
    # EXCLUDED: Subthalamus (11 vtx — subcortical)
    # EXCLUDED: Globus_pallidus (6 vtx — subcortical)
]

########################## Figure 1 ##########################
atlas_level_l = [2, 3, 4]


PLOT_REGIONS        = level2_cortical
REGIONS_OF_INTEREST = [level2_cortical, level3_cortical, level4_cortical]

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
# STEP 0  —  Stage 0: wb_command surface area + thickness extraction
# ═══════════════════════════════════════════════════════════════════════════════
# Run once per dataset.  After a successful run set both overwrite flags to
# False to skip already-computed subjects.
#

from Exemples.Study_EDNiX.atlases.Atlas_nomenclature import build_flat_atlas_config

run_wb_extraction(
    SING_WB,
    ALL_BIDS,
    atlas_name            = 'EDNIx',
    overwrite_thickness   = True,
    overwrite_surface     = True,
    n_atlas_levels        = N_ATLAS_LEVELS,
    atlas_config_override = build_flat_atlas_config(),
)