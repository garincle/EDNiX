import os
import sys
import glob

# ── EDNiX environment ─────────────────────────────────────────────────────────
MAIN_PATH        = '/home/cgarin/PycharmProjects/EDNiX/'
ATLAS_LIB        = os.path.join(MAIN_PATH, 'Atlases_library')
sys.path.insert(0, MAIN_PATH)
sys.path.insert(0, os.path.join(MAIN_PATH, 'code'))

# ── Species / atlas config ────────────────────────────────────────────────────
SPECIES_NAME     = 'Macaque'           # used in surface filenames
SPECIES_FRAGMENT = 'atlas/mammals/primates/Macaque'

# ── BIDS dataset ──────────────────────────────────────────────────────────────
BIDS_DIR = os.path.join('/scratch2/EDNiX/Macaque/BIDS_BenHamed')
RESULTS_DIR      = os.path.join(BIDS_DIR, 'Results')

# ── Atlas paths ───────────────────────────────────────────────────────────────
ATLAS_EDNIX      = os.path.join(ATLAS_LIB, SPECIES_FRAGMENT, 'EDNiX')
BASE_SS          = os.path.join(ATLAS_EDNIX, 'volumes', f'{SPECIES_NAME}_space-acpc_desc-SS_T1w.nii.gz')
MASK_FUNC        = os.path.join(ATLAS_EDNIX, 'volumes', 'masks', f'{SPECIES_NAME}_desc-Cerebral_Gray_mask.nii.gz')
SCENE_TEMPLATE   = os.path.join(ATLAS_EDNIX, 'surfaces', 'Native_resol', f'{SPECIES_NAME}_thickness.scene')
TEMPLATE_LOW     = '/scratch2/EDNiX/Macaque/BIDS_BenHamed/sub-Scooby/ses-1/func/templates/EDNiX/preprocessing/BASE_SS_fMRI.nii.gz'

# ── Singularity / wb_command ─────────────────────────────────────────────────
SING_WB = (
    'vglrun singularity exec '
    f'--bind {ATLAS_LIB}/,/scratch2/,/scratch/ '
    f'{MAIN_PATH}/Singularity_library/workbench_2.1.0.sif'
)

# ── Surface rendering config ─────────────────────────────────────────────────
SURFACE_CFG = dict(
    delete_surf_file = False,
    species_name     = SPECIES_NAME,
    atlas_lib_path   = ATLAS_LIB,
    species_fragment = SPECIES_FRAGMENT,
    sing_wb          = SING_WB,
    surface_type     = 'pial',
    palette_name     = 'RY-BC-BL',
    pos_max          = 90,
    neg_max          = 90,
    colorbar_mode    = 'percent',
    scene_template   = SCENE_TEMPLATE,
    scene_label      = 'New Scene 1',
)

# ── Statistical parameters ────────────────────────────────────────────────────
STAT_CFG = dict(
    alpha        = 0.05,
    lower_cutoff = 0.90,
    upper_cutoff = 0.95,
    cut_coords   = 5,
)

REGIONS_OF_INTEREST = [
    'L_BA_26_29_30__retrosplenial',
    'R_BA_26_29_30__retrosplenial',
]

# ── Imports ───────────────────────────────────────────────────────────────────
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
from Tools import Load_subject_with_BIDS, load_bids
from Plotting.ednix_bids_tools import extract_bold_paths, parse_label_file, get_atlas_label_path
import Statistics.Group_fMRI._Group_anal_3dTtest       as _ttest
import Statistics.Group_fMRI._Group_anal_3dLMEr_SBA    as _lmer_sba
import Statistics.Group_fMRI._Group_anal_3dLMEr_Mirror as _lmer_mirror
from Tools import niftotoWBsurface

# ── BIDS loading ──────────────────────────────────────────────────────────────
allinfo_study_c = load_bids.Load_BIDS_to_pandas(
    BIDS_DIR, modalities=['anat'], suffixes=['T1w'], extensions=['.nii.gz']
)
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

bold_paths = extract_bold_paths(
    bids_dir       = BIDS_DIR,
    list_to_keep   = [],
    list_to_remove = [],
    task           = 'rest',
)
print(f"  {len(bold_paths['bold_path'])} BOLD files found.")

# ── Mean functional images ────────────────────────────────────────────────────
POSTPRO_SUBDIR = os.path.join('templates', 'EDNiX', 'postprocessed_rs')
seen_sessions  = set()
mean_imgs      = []

for sub, ses, bold_p in zip(bold_paths['subject'], bold_paths['session'], bold_paths['bold_path']):
    key = (sub, ses)
    if key in seen_sessions:
        continue
    seen_sessions.add(key)
    ses_dir    = os.path.join(BIDS_DIR, f'sub-{sub}', f'ses-{ses}')
    postpro    = os.path.join(ses_dir, 'func', POSTPRO_SUBDIR)
    candidates = sorted(glob.glob(os.path.join(postpro, '*fMRI_Mean_Image_SS.nii.gz')))
    if candidates:
        mean_imgs.append(candidates[0])
    else:
        print(f"  [WARN] No mean image found in {postpro}")

print(f"  {len(mean_imgs)} mean images collected.")

# ── Atlas label file ──────────────────────────────────────────────────────────
label_path = get_atlas_label_path(ATLAS_LIB, SPECIES_FRAGMENT, 'EDNIxCSCLR')
label_df   = parse_label_file(label_path)
print(f"  {len(label_df)} regions in label file.")

# ── Shared kwargs for all analyses ────────────────────────────────────────────
SHARED = dict(
    bids_dir            = BIDS_DIR,
    templatehigh        = BASE_SS,
    templatelow         = TEMPLATE_LOW,
    oversample_map      = False,
    mask_func           = MASK_FUNC,
    label_df            = label_df,
    MAIN_PATH           = MAIN_PATH,
    bold_paths          = bold_paths,
    mean_imgs           = mean_imgs,
    regions_of_interest = REGIONS_OF_INTEREST,
    method_mask_func    = 'onlyprovidedmask',
    **STAT_CFG,
)


# ── ANALYSIS 1: 3dTtest++ ─────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  ANALYSIS 1  —  3dTtest++ seed-based")
print("=" * 70)
'''
_ttest._3dttest_EDNiX(**SHARED)

_lmer_mirror._3dLMEr_EDNiX(
    **SHARED,
    model          = "(1|Subj)*Hemisphere",
    glt_spec       = [
        (1, 'Overall_Effect',    'Hemisphere : 0.5*L +0.5*R'),
        (2, 'Hemisphere_Effect', 'Hemisphere : 1*L -1*R'),
    ],
    contrast_names = ['Overall_Effect', 'Hemisphere_Effect'],
    midline_x      = 0.0,
    visualize      = 'threshold',
    percent        = 5,
)
'''
niftotoWBsurface.render_stat_maps(
    results_dir         = RESULTS_DIR,
    results_subdir      = 'Grp_SBA_3dTTEST',
    stat_suffixes       = ['ttest-stat_fisher_zmap', 'ttest-stat_fisher_cmap'],
    regions_of_interest = REGIONS_OF_INTEREST,
    surface_cfg         = SURFACE_CFG,
)

# ── ANALYSIS 3: 3dLMEr Mirror ─────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  ANALYSIS 3  —  3dLMEr seed-based with hemisphere mirroring")
print("=" * 70)


niftotoWBsurface.render_stat_maps(
    results_dir         = RESULTS_DIR,
    results_subdir      = 'Grp_SBA_3dLME_network',
    stat_suffixes       = ['Overall_Effect_zstat', 'Overall_Effect_coef'],
    regions_of_interest = ['BA_26_29_30__retrosplenial'],
    surface_cfg         = SURFACE_CFG,
)