import os
import sys
import Tools.Read_atlas
import pet._0_Pipeline_launcher
from Tools import Load_subject_with_BIDS, load_bids
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

species    = 'Macaque'
# Override os.path.join to always return Linux-style paths
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Macaque/imagina/BIDS/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['pet'], suffixes= ['pet'], extensions=['.nii.gz'])

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = []
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

endpet = '*_trc-FDG_*.nii.gz' # string
endjson = '*_trc-FDG_*.json' # string

humanPosition     = ['']
orientation       = 'RPS' # "LPI" or ''
animalPosition    = [''] # valid only for species smaller than humans
## prior anat processing
coregistration_longitudinal = False #True or False
type_norm = 'T1w' # T1 or T2
### co-registration func to anat to template to with T1 ? T2? use the correct  suffix as in the BIDS
T_PET = 'T1w' # string
### if you don't have any anat image you will need to put several image in the folderforTemplate_Anat (refer to the doc)
Method_mask_pet = '3dSkullStrip_monkeynodil' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"
#### ANTs function of the co-registration HammingWindowedSinc is advised
IhaveanANAT = True # True or False
anat_pet_same_space = False # True or False
type_of_transform = 'BOLDAffine'
aff_metric_ants_Transl = 'mattes' # string
aff_metric_ants = 'mattes'
do_pet_to_func = True # True or False

##### if you don't have an anat then template will be the same as anat...
#creat_study_template was created with the anat type_norm img, and you want to use it as standart space
creat_study_template = False # True or Fals
blur = 0# float
#Dilate the functional brain mask by n layers
dilate_mask = 5 # int

############ Right in a list format the steps that you want to skip
doWARPonfunc = 'header'

Skip_step = ['itk_1', 'itk_2', 'Clean']
pet._0_Pipeline_launcher.preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endfmri, endjson,
                    animalPosition, humanPosition, orientation,
                    T_PET, type_norm, creat_study_template,
                    anat_pet_same_space, coregistration_longitudinal,
                    Method_mask_pet, do_pet_to_func, folderforTemplate_Anat='', IhaveanANAT=True,
                    overwrite_option=True,
                    doWARPonfunc=doWARPonfunc, registration_fast=False, type_of_transform=type_of_transform, n_for_ANTS='lanczosWindowedSinc', aff_metric_ants=aff_metric_ants, aff_metric_ants_Transl=aff_metric_ants_Transl, dilate_mask=dilate_mask,
                    list_to_keep=list_to_keep, list_to_remove=list_to_remove, atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
                    reference='EDNiX', use_erode_WM_func_masks = True, use_erode_V_func_masks=True,
                    blur=0,
                    selected_atlases_matrix='all', wanted_level_matrix='all',
                    selected_atlases_SBA='default', panda_files_SBA='default')