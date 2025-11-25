import pandas as pd
import os
from bids import BIDSLayout
from bids.reports import BIDSReport
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import fMRI._0_Pipeline_launcher
from Tools import Load_subject_with_BIDS
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
bids_dir = Load_subject_with_BIDS.linux_path(opj('/home/cgarin/Documents/Snake/'))

########### Subject loader with BIDS  ##############
layout = BIDSLayout(bids_dir, validate=False)
report = BIDSReport(layout)
df = layout.to_df()
df.head()

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = df[(df['suffix'] == 'bold') & (df['extension'] == '.nii.gz')]
### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# choose if you want to select or remove ID from you analysis:
list_to_keep   = []
list_to_remove = []
species    = 'Snake'
reference  = 'EDNiX'
addatlas   = ''

all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

#### functional images paramters definition
Slice_timing_info = 'Auto'
ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii.gz' # string
endjson = '*_task-rest_*.json' # string
endmap = '*_map.nii.gz' # string
# "AHF" stands for "Animal Head First",  "AFF" stands for "Animal Feet First", "humanlike" means no change to be done.
humanPosition     = ['']
orientation       = 'LIP' # "LPI" or ''
animalPosition    = [''] # valid only for species smaller than humans
## prior anat processing
coregistration_longitudinal = False #True or False
type_norm = 'T2w' # T1 or T2
TfMRI = 'T2w' # string

## masking
Method_mask_func = 'Vol_sammba_400' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"
creat_study_template = True # True or False

doWARPonfunc = True
resting_or_task = 'resting'  # 'resting' or 'task'
anat_func_same_space = True # True or False
type_of_transform = 'BOLDAffine'
aff_metric_ants_Transl = 'mattes' # string
aff_metric_ants = 'mattes'
do_anat_to_func = False # True or False
dilate_mask=0

Skip_step = [1,2,3,4,6,7,8,9,10,12,13,14,'Clean']
fMRI._0_Pipeline_launcher.preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endfmri, endjson, endmap, resting_or_task,
                    animalPosition, humanPosition, orientation,
                    Slice_timing_info,
                    TfMRI, type_norm, creat_study_template,
                    anat_func_same_space, coregistration_longitudinal,
                    Method_mask_func, do_anat_to_func, folderforTemplate_Anat='', IhaveanANAT=True,
                    ntimepoint_treshold=100, REF_int=0, T1_eq=5, correction_direction='Auto', overwrite_option=True,
                    DwellT='Auto', SED='Auto', TR='Auto', TRT='Auto',
                    nb_ICA_run=20, ICA_cleaning='Skip',
                    costAllin='lpa',
                    doWARPonfunc=doWARPonfunc, registration_fast=False, type_of_transform=type_of_transform, n_for_ANTS='lanczosWindowedSinc', aff_metric_ants=aff_metric_ants, aff_metric_ants_Transl=aff_metric_ants_Transl, dilate_mask=dilate_mask,
                    list_to_keep=list_to_keep, list_to_remove=list_to_remove, atlas_followers=[['Seedsrandom'], ['ctab'], [1], [0]],
                    reference='EDNiX', post_treatment_method='AFNI',
                    band='0.01 0.1', blur=0, do_not_correct_signal = False, extract_exterior_CSF = False, extract_WM=False, extract_Vc = False, extract_GS = True,
                    use_erode_WM_func_masks = False, use_erode_V_func_masks=False, normalize='Skip',
                    selected_atlases_matrix='all', wanted_level_matrix='all',
                    selected_atlases_SBA=['Seedsrandom', 0], panda_files_SBA=[pd.DataFrame({'region':['1', '2', '3', '4','5', '6',],'label':[1,2,3,4,5,6]})],
                    SBAspace=['func'], erod_seed=True, smoothSBA=0.7,
                    specific_roi_tresh=0.2, delta_thresh=0.1,
                    oversample_map=False, use_cortical_mask_func=False, n_cut=10, threshold_val=10)
