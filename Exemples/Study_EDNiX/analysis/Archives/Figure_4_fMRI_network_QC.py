"""
ARCHIVED (former Figure 4, network QC): not in the paper any more. Kept because
Fig5/fig5_pipeline.CLASSIFICATION_PATH still reads its output
(Figure_4_rsfMRI_QC/lvl2/all/classification_results.xlsx) -- do not delete that file.

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

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis.ednix_bids_tools import (
    get_atlas_label_path, find_species_path,
    collect_multi_species,
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

ALL_BIDS = (
    [(sp, bd) for sp, bd in species_bids_dict.items()] +
    [(sp, bd) for sp, bds in species_multi_bids.items() for bd in bds]
)

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

from Exemples.Study_EDNiX.analysis.Archives.ednix_threshold_explorer import run_pipeline_multilevel
data = collect_multi_species(
    species_config,
    regions_of_interest     = None,
    extract                 = ('qc'),
    atlas_name              = ATLAS_NAME,
    atlas_label_paths       = atlas_label_paths,
    fit_kind='correlation',
    atlas_library_root      = ATLAS_LIB,
    species_atlas_fragments = species_atlas_fragments,
)

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
                            thresh_intra=None, thresh_inter=None, stringency=0.0,
                            n_top=10, atlas_name="EDNIxCSC", use_lr=True,
                            fp_test="species_mean", fdr_alpha=0.05,
                            fig_dir  = "/scratch2/EDNiX/results/EDNiX_Figures/Figure_4_rsfMRI_QC",
                            use_inter_all=False,
                            normalize_method="translate",
                            verbose=True)