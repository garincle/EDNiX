import os
import sys
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS, load_bids
from modalities.pet import _0_Pipeline_launcher

######################################################################################################################
# Castrillon et al. 2023 (Sci. Adv. 9:eadi7632): simultaneous [18F]FDG PET / BOLD fMRI on an integrated 3T PET/MR
# (Siemens Biograph mMR), list-mode DYNAMIC acquisition (33 frames over 45 min: 10x12 s, 8x30 s, 8x60 s, 2x180 s,
# 5x300 s), bolus injection, eyes open (second TUM.exp session: closed). Sessions here: 'open' / 'closed'.
# T1w MPRAGE processed by the anatomical launcher of BIDS_Castrillon (must have run first).
# Scenario 1v2a: rigid PET -> T1w registration, then the anatomy -> EDNiX template chain of the anatomical pipeline.
# The PET is DYNAMIC (unlike the static 10-min frame of CERMEP): the ratio image is the mean of the LATE frames
# (20-45 min, the window of the Patlak fit of the paper), not of the whole acquisition (dominated by the blood pool).
######################################################################################################################
species = 'Human'
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Human/BIDS_Castrillon/'))
# the PET files are listed (as in the CERMEP and Tremblay launchers), not the T1w
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['pet'], suffixes=['pet'], extensions=['.nii.gz'])

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
list_to_keep = []
list_to_remove = [('s020', 'open'), ('s020', 'closed'), ('s023', 'closed'), ('s037', 'open')]

### PET images: CHECK the file names first (ls BIDS_Castrillon/sub-s003/ses-open/pet); the pattern must match ONE 4D PET
### file per session (not a parametric map such as a CMRglc or Ki image)
endpet = '*pet.nii.gz'
tracer = 'FDG'
isotope = 'F-18'

### orientation (humans: nothing to change)
humanPosition = ['humanlike']
animalPosition = ['humanlike']
orientation = ''

### anatomy (same values as in the anatomical launcher of BIDS_Castrillon)
type_norm = 'T1w'
T_PET = 'T1w'
creat_study_template = False
coregistration_longitudinal = False
IhaveanANAT = True
anat_pet_same_space = False       # hybrid PET/MR, but the PET is registered anyway (1v2a works whether or not the images
                                  # are already aligned; True would trust the scanner coordinates without checking)
IhaveaCT = False                  # no CT: attenuation correction from the MR (pseudo-CT)
ct_pet_same_space = False

### motion correction (dynamic run)
nb_void = 10                      # first frames not registered (low counts, blood pool): the 10 x 12 s frames = 2 min.
                                  # CHECK against FrameTimesStart of the sidecar if the frames were rebinned
crop_for_MoCo = 'auto'            # head box drawn in freeview on the session mean the first time; 'none' to run without
                                  # display

### parametric images
SUV = True                        # needs InjectedRadioactivity (sidecar) and body weight (participants.tsv); set False if
                                  # they are missing: the ratio to the brain mean (bratio) does not need them
bratio = True
SUV_window = [20, 45]             # minutes after injection: late frames, as the Patlak window of the paper
Do_cut = False
duration = 10

### Reset before relaunch (modalities/pet/_000_Data_Reset.py): deletes the EDNiX outputs of the selected sessions and
### keeps the raw BIDS files and the crops drawn by hand (desc-crop). RESET_KEEP_REGISTRATION keeps orig/matrices.
### RESET_BEFORE_RUN = True with RESET_DRY_RUN = True previews; RESET_DRY_RUN = False deletes after ONE 'YES'.
RESET_BEFORE_RUN = False
RESET_DRY_RUN = True
RESET_KEEP_REGISTRATION = False
if RESET_BEFORE_RUN:
    from modalities.pet import _000_Data_Reset
    _000_Data_Reset.reset_all(allinfo_study_c, bids_dir, list_to_keep=list_to_keep, list_to_remove=list_to_remove,
                              dry_run=RESET_DRY_RUN, keep_registration=RESET_KEEP_REGISTRATION)

Skip_step = ['pvc', 9]            # no partial-volume correction in this first pass (the paper applied iterative Yang PVC
                                  # with GM / WM / CSF masks: a later sensitivity analysis), no PET during fMRI
_0_Pipeline_launcher.preprocess_pet(
    Skip_step, MAIN_PATH, bids_dir, species, allinfo_study_c, endpet, tracer,
    animalPosition, humanPosition, orientation,
    T_PET, type_norm, creat_study_template, coregistration_longitudinal,
    IhaveanANAT=IhaveanANAT, anat_pet_same_space=anat_pet_same_space,
    IhaveaCT=IhaveaCT, ct_pet_same_space=ct_pet_same_space,
    isotope=isotope, nb_void=nb_void, crop_for_MoCo=crop_for_MoCo,
    SUV=SUV, bratio=bratio, SUV_window=SUV_window, Do_cut=Do_cut, duration=duration,
    export_resolution='pet', selected_atlases='all',
    atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
    reference='EDNiX', list_to_keep=list_to_keep, list_to_remove=list_to_remove,
    Need_check_orient=False, Need_check_MC=False, Need_check_coreg=False)
