import os
import sys
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS, load_bids
from modalities.pet import _0_Pipeline_launcher

######################################################################################################################
# CERMEP-IDB-MRXFDG (Merida et al. 2021, EJNMMI Res 11:91): [18F]FDG PET/CT (Biograph mCT64, one 10-min frame 50 min
# after injection), T1w MPRAGE (1.5 T) processed by Launcher_anat_Human_CERMEP.py.
# Scenario 1v2a (default): rigid PET -> T1w registration. The authors distribute the NIfTI images coregistered to
# the PET (SPM12, Merida et al. 2021): with anat_pet_same_space = True (scenario 1v1) the PET follows the T1w without
# registration and relies on that coregistration (check the QC image); 1v2a re-estimates it and works in both cases. Through the CT (scenario 1v2b): run ../Launcher_ct/Launcher_ct_CERMEP.py first, then set IhaveaCT = True.
######################################################################################################################
species = 'Human'
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Human/BIDS_Merida/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['pet'], suffixes=['pet'], extensions=['.nii.gz'])

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
list_to_keep = []
list_to_remove = []

### PET images
endpet = '*pet.nii.gz'
tracer = 'FDG'
isotope = 'F-18'

### orientation (humans: nothing to change)
humanPosition = ['humanlike']
animalPosition = ['humanlike']
orientation = ''

### anatomy (same values as in Launcher_anat_Human_CERMEP.py)
type_norm = 'T1w'
T_PET = 'T1w'
creat_study_template = False
coregistration_longitudinal = True
IhaveanANAT = True
anat_pet_same_space = False       # True with the NIfTI set coregistered to the PET (scenario 1v1)
IhaveaCT = False                  # True after the CT pipeline (scenario 1v2b)
ct_pet_same_space = True          # PET/CT scanner

### parametric images (static acquisition: no motion correction, no kinetic model)
SUV = True                        # images in Bq/mL: needs InjectedRadioactivity and body_weight; images in SUV: kept
bratio = True

Skip_step = ['pvc', 8, 9]
_0_Pipeline_launcher.preprocess_pet(
    Skip_step, MAIN_PATH, bids_dir, species, allinfo_study_c, endpet, tracer,
    animalPosition, humanPosition, orientation,
    T_PET, type_norm, creat_study_template, coregistration_longitudinal,
    IhaveanANAT=IhaveanANAT, anat_pet_same_space=anat_pet_same_space,
    IhaveaCT=IhaveaCT, ct_pet_same_space=ct_pet_same_space,
    isotope=isotope, SUV=SUV, bratio=bratio,
    export_resolution='pet', selected_atlases='all',
    atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
    reference='EDNiX', list_to_keep=list_to_keep, list_to_remove=list_to_remove,
    Need_check_orient=False, Need_check_coreg=False)
