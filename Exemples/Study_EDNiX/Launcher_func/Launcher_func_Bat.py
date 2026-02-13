import os
import sys
import Tools.Read_atlas
import fMRI._0_Pipeline_launcher
import pandas as pd
from Tools import Load_subject_with_BIDS, load_bids
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

species    = 'Bat'
# Override os.path.join to always return Linux-style paths
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Bat/BIDS_bat/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['func'], suffixes= ['bold'], extensions=['.nii.gz'])

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = [('10', '1')]
#
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

#### fMRI pre-treatment
T1_eq = 5 # int
REF_int = 0 # int
ntimepoint_treshold = 100
endfmri = '*_task-restingstate_*.nii.gz' # string
endjson = '*_task-restingstate_*.json' # string
endmap = '*fMRI_epi.nii.gz' # string
humanPosition     = ['']
orientation       = 'LAI' # "LPI" or ''
animalPosition    = [''] # valid only for species smaller than humans
## prior anat processing
coregistration_longitudinal = False #True or False
type_norm = 'T2w' # T1 or T2
### co-registration func to anat to template to with T1 ? T2? use the correct  suffix as in the BIDS
TfMRI = 'T2w' # string
### if you don't have any anat image you will need to put several image in the folderforTemplate_Anat (refer to the doc)
Method_mask_func = 'NoSkullStrip' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"
#### ANTs function of the co-registration HammingWindowedSinc is advised
IhaveanANAT = True # True or False
anat_func_same_space = False # True or False
type_of_transform = 'SyNCC'
aff_metric_ants_Transl = 'CC' # string
aff_metric_ants = 'CC'
do_anat_to_func = True # True or False
Slice_timing_info = 'Auto'
##### if you don't have an anat then template will be the same as anat...
#creat_study_template was created with the anat type_norm img, and you want to use it as standart space
creat_study_template = True # True or Fals
blur = 0# float
#Dilate the functional brain mask by n layers
dilate_mask = 0 # int
smoothSBA = 0.4
# for the seed base analysis, you need to provide the names and the labels of the regions you want to use as "seeds"
selected_atlases = [['pcc', 0]]  # Using NEW VERSION format (single atlas)
############ Right in a list format the steps that you want to skip
doWARPonfunc = False
resting_or_task = 'resting'  # 'resting' or 'task'

Skip_step = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,'itk_1', 'Clean']
fMRI._0_Pipeline_launcher.preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endfmri, endjson, endmap, resting_or_task,
                    animalPosition, humanPosition, orientation,
                    Slice_timing_info,
                    TfMRI, type_norm, creat_study_template,
                    anat_func_same_space, coregistration_longitudinal,
                    Method_mask_func, extra_erode=0, do_anat_to_func=do_anat_to_func, folderforTemplate_Anat='', IhaveanANAT=True,
                    ntimepoint_treshold=100, REF_int=0, T1_eq=5, correction_direction='Auto', overwrite_option=True,
                    DwellT='Auto', SED='Auto', TR=2, TRT='Auto',
                    nb_ICA_run=20, ICA_cleaning='Skip',
                    costAllin='lpa',
                    doWARPonfunc=doWARPonfunc, registration_fast=False, type_of_transform=type_of_transform, n_for_ANTS='lanczosWindowedSinc', aff_metric_ants=aff_metric_ants, aff_metric_ants_Transl=aff_metric_ants_Transl, dilate_mask=dilate_mask,
                    list_to_keep=list_to_keep, list_to_remove=list_to_remove, atlas_followers=[['Stuart', 'pcc'], ['ctab', 'txt'], [1, 1], [1, 1]],
                    reference='EDNiX', post_treatment_method='Grandjean',
                    band='0.01 0.1', blur=0, do_not_correct_signal = False, extract_exterior_CSF = False, extract_WM=True, extract_Vc = False, extract_GS = False,
                    use_erode_WM_func_masks = True, use_erode_V_func_masks=True, normalize='Skip',
                    selected_atlases_matrix='all', wanted_level_matrix='all',
                    selected_atlases_SBA=['pcc', 0], panda_files_SBA=[pd.DataFrame({'region':['pcc'],'label':[1]})],
                    SBAspace=['atlas'], erod_seed=True, smoothSBA=smoothSBA,
                    specific_roi_tresh=0.2, delta_thresh=0.1,
                    oversample_map=False, use_cortical_mask_func=False, n_cut=10, threshold_val=10)