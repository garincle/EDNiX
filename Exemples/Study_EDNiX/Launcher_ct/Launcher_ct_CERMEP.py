import os
import sys
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS, load_bids
from modalities.ct import _0_Pipeline_launcher

######################################################################################################################
# Low-dose CT of the CERMEP-IDB-MRXFDG PET/CT: rigid CT -> T1w registration (scenario 1v2b of Illustration/pet4.svg).
# Run after Launcher_anat_Human_CERMEP.py and before ../Launcher_pet/Launcher_fdg_CERMEP.py with IhaveaCT = True.
######################################################################################################################
species = 'Human'
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Human/BIDS_Merida/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['ct'], suffixes=['ct'], extensions=['.nii.gz'])

Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
list_to_keep = []
list_to_remove = []

endct = '*_ct.nii.gz'             # raw CT in sub-*/ses-*/ct
humanPosition = ['humanlike']     # same values as the PET of the PET/CT
animalPosition = ['humanlike']
orientation = ''

type_norm = 'T1w'
T_CT = 'T1w'
creat_study_template = False
coregistration_longitudinal = True
IhaveanANAT = True

### Reset before relaunch (modalities/ct/_000_Data_Reset.py): deletes the EDNiX outputs of the selected sessions and
### keeps the raw BIDS files. RESET_KEEP_REGISTRATION keeps orig/matrices.
### RESET_BEFORE_RUN = True with RESET_DRY_RUN = True previews; RESET_DRY_RUN = False deletes after ONE 'YES'.
RESET_BEFORE_RUN = False
RESET_DRY_RUN = True
RESET_KEEP_REGISTRATION = False
if RESET_BEFORE_RUN:
    from modalities.ct import _000_Data_Reset
    _000_Data_Reset.reset_all(allinfo_study_c, bids_dir, list_to_keep=list_to_keep, list_to_remove=list_to_remove,
                              dry_run=RESET_DRY_RUN, keep_registration=RESET_KEEP_REGISTRATION)

Skip_step = [1,2]
_0_Pipeline_launcher.preprocess_ct(
    Skip_step, MAIN_PATH, bids_dir, species, allinfo_study_c, endct,
    animalPosition, humanPosition, orientation,
    T_CT, type_norm, creat_study_template, coregistration_longitudinal,
    IhaveanANAT=IhaveanANAT, reference='EDNiX', coreg_metric='mattes', coreg_init='auto',
    atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
    list_to_keep=list_to_keep, list_to_remove=list_to_remove,
    Need_check_orient=False, Need_check_coreg=False)
