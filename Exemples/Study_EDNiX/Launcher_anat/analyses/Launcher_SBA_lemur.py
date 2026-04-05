import os
import sys
import glob
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
from Tools import Load_subject_with_BIDS, load_bids

# ─────────────────────────────────────────────────────────────────────────────
# EDNiX environment

MAIN_PATH = os.path.join('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, os.path.join(MAIN_PATH))
# ─────────────────────────────────────────────────────────────────────────────
# BIDS dataset
# ─────────────────────────────────────────────────────────────────────────────
species  = 'Mouse_lemur'
bids_dir = os.path.join('/scratch2/EDNiX', species, 'BIDS_Garin')

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['anat'], suffixes= ['T2w'], extensions=['.nii.gz'])
### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# ─────────────────────────────────────────────────────────────────────────────
# Import EDNiX tools
# ─────────────────────────────────────────────────────────────────────────────
sys.path.insert(0, os.path.join(MAIN_PATH, 'code'))
from Plotting.ednix_bids_tools import extract_bold_paths, parse_label_file, get_atlas_label_path

import Statistics.Group_fMRI._Group_anal_3dTtest      as _ttest
import Statistics.Group_fMRI._Group_anal_3dLMEr_SBA   as _lmer_sba
import Statistics.Group_fMRI._Group_anal_3dLMEr_Mirror as _lmer_mirror

# ─────────────────────────────────────────────────────────────────────────────
# Subjects / sessions to exclude
# ─────────────────────────────────────────────────────────────────────────────
LIST_TO_REMOVE = []
LIST_TO_KEEP = []   # empty = keep all (after removals)

# ─────────────────────────────────────────────────────────────────────────────
# BOLD paths  — extract_bold_paths() handles all path resolution
# Returns dict: subject, session, run, bold_path
# ─────────────────────────────────────────────────────────────────────────────
bold_paths = extract_bold_paths(
    bids_dir       = bids_dir,
    list_to_keep   = LIST_TO_KEEP,
    list_to_remove = LIST_TO_REMOVE,
    task           = 'rest',
)
print(f"  {len(bold_paths['bold_path'])} BOLD files found.")

# ─────────────────────────────────────────────────────────────────────────────
# Mean functional images — one per session
# Glob for *mean*.nii.gz in each session's postprocessed_rs/ folder
# ─────────────────────────────────────────────────────────────────────────────
POSTPRO_SUBDIR = os.path.join('templates', 'EDNiX', 'postprocessed_rs')

seen_sessions = set()
mean_imgs     = []

for sub, ses, bold_p in zip(bold_paths['subject'], bold_paths['session'], bold_paths['bold_path']):
    key = (sub, ses)
    if key in seen_sessions:
        continue
    seen_sessions.add(key)

    ses_dir    = os.path.join(bids_dir, f'sub-{sub}', f'ses-{ses}')
    postpro    = os.path.join(ses_dir, 'func', POSTPRO_SUBDIR)
    candidates = sorted(glob.glob(os.path.join(postpro, '*fMRI_Mean_Image_SS.nii.gz')))
    if candidates:
        mean_imgs.append(candidates[0])
    else:
        print(f"  [WARN] No mean image found in {postpro}")

print(f"  {len(mean_imgs)} mean images collected.")

# ─────────────────────────────────────────────────────────────────────────────
# Atlas label file  →  label_df
# parse_label_file() returns a DataFrame:
#   region_name, base_region, label_id, hemisphere, R, G, B, A
# ─────────────────────────────────────────────────────────────────────────────
ATLAS_LIB  = os.path.join(MAIN_PATH, 'Atlases_library')
SPECIES_FRAGMENT = 'atlas/mammals/primates/Mouselemur'
label_path = get_atlas_label_path(ATLAS_LIB, SPECIES_FRAGMENT, 'EDNIxCSCLR')
label_df   = parse_label_file(label_path)
print(f"  {len(label_df)} regions in label file.")

# Regions to process (base_region substrings); None = all regions
regions_of_interest = [
    'L_BA_26_29_30__retrosplenial',
    'R_BA_26_29_30__retrosplenial']

# ─────────────────────────────────────────────────────────────────────────────
# Study template & masks
# ─────────────────────────────────────────────────────────────────────────────
BASE_SS        = os.path.join(ATLAS_LIB, SPECIES_FRAGMENT, 'EDNiX/volumes/', 'Mouselemur_space-acpc_desc-SS_T1w.nii.gz')
mask_func      = os.path.join(ATLAS_LIB, SPECIES_FRAGMENT, 'EDNiX/volumes/masks/Mouselemur_desc-Cerebral_Gray_mask.nii.gz')
templatehigh   = BASE_SS
templatelow    = '/scratch2/EDNiX/Mouse_lemur/BIDS_Garin/sub-147BCBB/ses-01/func/templates/EDNiX/preprocessing/BASE_SS_fMRI.nii.gz'
oversample_map = False

# ─────────────────────────────────────────────────────────────────────────────
# Shared statistical parameters
# ─────────────────────────────────────────────────────────────────────────────
alpha        = 0.001
lower_cutoff = 0.90
upper_cutoff = 0.95
cut_coords   = [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6]

# ─────────────────────────────────────────────────────────────────────────────
# ── ANALYSIS 1 :  3dTtest++ (one-sample t-test, SBA) ─────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  ANALYSIS 1  —  3dTtest++ seed-based")
print("=" * 70)

_ttest._3dttest_EDNiX(
    bids_dir            = bids_dir,
    templatehigh        = templatehigh,
    templatelow         = templatelow,
    oversample_map      = oversample_map,
    mask_func           = mask_func,
    cut_coords          = cut_coords,
    label_df            = label_df,
    lower_cutoff        = lower_cutoff,
    upper_cutoff        = upper_cutoff,
    MAIN_PATH           = MAIN_PATH,
    alpha               = alpha,
    bold_paths          = bold_paths,
    mean_imgs           = mean_imgs,
    regions_of_interest = regions_of_interest)

# ─────────────────────────────────────────────────────────────────────────────
# ── ANALYSIS 2 :  3dLMEr SBA (mixed-effects, no mirroring) ───────────────────
# model   : random subject intercept
# gltCode : group effect averaged across sessions
# stat_img: one label per -gltCode, in the same order
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  ANALYSIS 2  —  3dLMEr seed-based (no mirroring)")
print("=" * 70)

_lmer_sba._3dLMEr_EDNiX(
    bids_dir            = bids_dir,
    templatehigh        = templatehigh,
    templatelow         = templatelow,
    oversample_map      = oversample_map,
    mask_func           = mask_func,
    cut_coords          = cut_coords,
    label_df            = label_df,
    lower_cutoff        = lower_cutoff,
    upper_cutoff        = upper_cutoff,
    MAIN_PATH           = MAIN_PATH,
    alpha               = alpha,
    bold_paths          = bold_paths,
    mean_imgs           = mean_imgs,
    model               = '"Sess*(1|Subj)"',
    gltCode             = "-gltCode groupeffect 'Sess : 0.5*Sess_01 +0.5*Sess_02'",
    stat_img            = ['groupeffect'],
    regions_of_interest = regions_of_interest,
)


# ─────────────────────────────────────────────────────────────────────────────
# ── ANALYSIS 3 :  3dLMEr Mirror (mixed-effects + hemisphere mirroring) ────────
# glt_spec      : list of (index, label, contrast) tuples
# contrast_names: matching labels used for sub-brick extraction and file naming
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  ANALYSIS 3  —  3dLMEr seed-based with hemisphere mirroring")
print("=" * 70)

_lmer_mirror._3dLMEr_EDNiX(
    bids_dir            = bids_dir,
    templatehigh        = templatehigh,
    templatelow         = templatelow,
    oversample_map      = oversample_map,
    mask_func           = mask_func,
    cut_coords          = cut_coords,
    label_df            = label_df,
    lower_cutoff        = lower_cutoff,
    upper_cutoff        = upper_cutoff,
    MAIN_PATH           = MAIN_PATH,
    alpha               = alpha,
    bold_paths          = bold_paths,
    mean_imgs           = mean_imgs,
    model               = "(1|Subj)*Hemisphere",
    glt_spec            = [(1, 'Overall_Effect', 'Hemisphere : 0.5*L +0.5*R'),
    (2, 'Hemisphere_Effect', 'Hemisphere : 1*L -1*R')],
    contrast_names      = ['Overall_Effect', 'Hemisphere_Effect'],
    midline_x           = 0.0,
    visualize           = 'threshold',   # 'threshold' | 'percentile' | None
    percent             = 5,
    regions_of_interest = regions_of_interest,)