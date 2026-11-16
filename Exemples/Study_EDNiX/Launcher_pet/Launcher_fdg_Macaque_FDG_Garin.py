import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))  # the EDNiX folder
opn = os.path.normpath
opj = os.path.join

# MAIN_PATH: this EDNiX folder; DATA_ROOT: folder of the BIDS datasets (environment variable EDNIX_DATA)
from Tools.ednix_paths import MAIN_PATH, DATA_ROOT
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS, load_bids
from modalities.pet import _0_Pipeline_launcher

######################################################################################################################
# [18F]FDG dynamic PET of the macaques rebuilt from the split_FDG frames (Vanderbilt, Philips Vereos PET/CT): 7 animals,
# 9 sessions (Flipper and Olvia x2), 31 frames on a 256x256x164 grid at 1 mm, 2-161 min post-injection (group average).
# The real frame timing is UNKNOWN: fill_macaque_json.py wrote a NOMINAL timing (31 equal frames from 120 to 9660 s,
# 'FrameTimingNominal' in the sidecar). Equal durations -> every duration-weighted mean is an unweighted mean, so the
# results depend only on which frames are selected. Valid here: motion correction + bratio of the late window
# SUV_window=[84, 99] min (3 late frames, ~17-19 of 31, plateau of the brain signal). NOT valid: kinetic models (srtm),
# SUV, time-activity curves.
# Alternative without dynamic processing: endpet='*_trc-FDG_rec-late3_pet.nii.gz' (static mean of the last 3 frames).
# No MRI in this dataset: scenario 3v2 (PET study template = common space of the sessions); pet_template_bridge=True
# (default) registers the EDNiX macaque template to it and brings the atlases back to every session.
# Figure 6 uses the baseline (earliest, pre-implantation) session of each animal: multimodal_figure.SESSION_RULE.
# Manual brain masks (native PET space, one per session): pet/raw/masks/ (read by EDNiX since PR #87: brain mask of bratio).
######################################################################################################################
species = 'Macaque'
bids_dir = Load_subject_with_BIDS.linux_path(opj(DATA_ROOT, 'Macaque/BIDS_Macaque_FDG_Garin/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['pet'], suffixes=['pet'], extensions=['.nii.gz'])

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
list_to_keep = []
list_to_remove = []

### PET images
endpet = '*_trc-FDG_pet.nii.gz'            # dynamic image (not the rec-late3 static image)
tracer = 'FDG'
isotope = 'F-18'

### orientation: unknown (the 1 mm grid may already have been reoriented): check on one session, then fix
humanPosition = ['']
animalPosition = ['']
orientation = ''
Need_check_orient = False

### no anatomy, no CT
type_norm = 'T1w'
T_PET = 'T1w'
creat_study_template = False
coregistration_longitudinal = False
IhaveanANAT = False
anat_pet_same_space = False
IhaveaCT = False
creat_study_template_pet = True   # common space of 3v2 + crop drawn once (step 3) and brought to every session

### motion correction (dynamic, 31 frames of ~5 min from 2 min post-injection)
nb_void = 3                        # first 3 frames (low brain uptake) left out of the motion estimation
crop_for_MoCo = 'mask'
redo_crop = False

### parametric images
SUV = False                       # per-animal dose / weight unknown (group values only)
bratio = True
SUV_window = [84, 99]            # late window: 3 late frames (~17-19 of 31) on the nominal timing, plateau of the brain signal
Do_cut = False
duration = 10

Skip_step = [1,'pvc', 9]
_0_Pipeline_launcher.preprocess_pet(
    Skip_step, MAIN_PATH, bids_dir, species, allinfo_study_c, endpet, tracer,
    animalPosition, humanPosition, orientation,
    T_PET, type_norm, creat_study_template, coregistration_longitudinal,
    IhaveanANAT=IhaveanANAT, anat_pet_same_space=anat_pet_same_space,
    IhaveaCT=IhaveaCT, creat_study_template_pet=creat_study_template_pet, studyT_transform='SyNRA',
    isotope=isotope, nb_void=nb_void, crop_for_MoCo=crop_for_MoCo, redo_crop=redo_crop,
    SUV=SUV, bratio=bratio, SUV_window=SUV_window, Do_cut=Do_cut, duration=duration,
    export_resolution='pet', selected_atlases='all', coreg_metric='mattes',
    atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
    reference='EDNiX', list_to_keep=list_to_keep, list_to_remove=list_to_remove,
    Need_check_orient=Need_check_orient, Need_check_MC=False, Need_check_coreg=False)
