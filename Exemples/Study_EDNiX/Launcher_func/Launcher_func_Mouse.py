import pandas as pd
import os
from bids import BIDSLayout
opn = os.path.normpath
opj = os.path.join
MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import fonctions._0_Pipeline_launcher

species = 'Mouse'
## linux ##FS_dir
# Override os.path.join to always return Linux-style paths
bids_dir = '/srv/projects/easymribrain/data/MRI/Mouse/BIDS_Gd'

########### Subject loader with BIDS##############
layout= BIDSLayout(bids_dir,  validate=False)
df = layout.to_df()
df.head()

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = df[(df['suffix'] == 'T2w') & (df['extension'] == '.nii')]

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = [('jgrAesAWc11R', '1'),
    ('jgrAesAWc11R1L', '0'),
    ('jgrAesAWc12R', '1'),
    ('jgrAesAWc1NT', '1')]

#### functional images paramters definition
Slice_timing_info = '-tpattern alt+z'
TR = '1.2'  # 'value du calculate in s', 'Auto', 'None'

ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii' # string
endjson = '*_task-rest_*.json' # string
endmap = '*_map.nii.gz' # string
orientation = 'LSP' # string
humanPosition     = ['']
animalPosition    = [''] # valid only for species smaller than humans

## prior anatomical processing
coregistration_longitudinal = False #True or False
type_norm = 'acq-RARE_T2w' # T1 or T2
### co-registration func to anat to template to with T1 ? T2? use the correct  suffix as in the BIDS
TfMRI = 'acq-RARE_T2w' # string
### if you don't have any anatomical image you will need to put several image in the folderforTemplate_Anat (refer to the doc)
folderforTemplate_Anat = ''

## masking
doMaskingfMRI = True # True or False
Method_mask_func = 'Vol_sammba_400' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"

#### ANTs function of the co-registration HammingWindowedSinc is advised
IhaveanANAT = True # True or False
anat_func_same_space = True # True or False
use_master_for_Allineate = False
n_for_ANTS = 'hammingWindowedSinc' # string
registration_fast = False
type_of_transform = 'BOLDAffine'
aff_metric_ants_Transl = 'mattes' # string
aff_metric_ants = 'mattes'
do_anat_to_func = True # True or False

##### if you don't have an anat then template will be the same as anat...
#creat_study_template was created with the anat type_norm img, and you want to use it as standart space
creat_study_template = True # True or False
#folder where you stored the stdy template
study_template_atlas_forlder = '/scratch/cgarin/Mouse/BIDS_Gd/sty_template/'  # string
stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_mask.nii.gz') # sting
stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz') # sting
GM_mask_studyT = opj(study_template_atlas_forlder, 'atlases', 'Gmask.nii.gz') # sting
resting_or_task = 'resting'  # 'resting' or 'task'
doWARPonfunc = False
#Smooth
blur = 0 # float (changed from 0.2 to 0 as requested)
#Dilate the functional brain mask by n layers
dilate_mask = 4 # int
SBAspace = ['func', 'atlas'] #list containing at least on of the string 'func', 'anat', 'atlas'
smoothSBA = 0.2

############ Right in a list format the steps that you want to skip
Skip_step = [4,100,200] # Changed to match requested values
fonctions._0_Pipeline_launcher.preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endfmri, endjson, endmap, resting_or_task,
                    animalPosition, humanPosition, orientation,
                    Slice_timing_info,
                    TfMRI, type_norm, creat_study_template,
                    anat_func_same_space, coregistration_longitudinal,
                    Method_mask_func, do_anat_to_func, folderforTemplate_Anat='', IhaveanANAT=True,
                    ntimepoint_treshold=100, REF_int=0, T1_eq=5, correction_direction='Auto', overwrite_option=True,
                    DwellT='Auto', SED='Auto', TR=TR, TRT='Auto',
                    nb_ICA_run=20, ICA_cleaning='Skip',
                    costAllin='lpa',
                    doWARPonfunc=doWARPonfunc, registration_fast=False, type_of_transform=type_of_transform, n_for_ANTS='lanczosWindowedSinc', aff_metric_ants=aff_metric_ants, aff_metric_ants_Transl=aff_metric_ants_Transl, dilate_mask=dilate_mask,
                    list_to_keep=[], list_to_remove=[], atlas_followers=[[], [], [], []],
                    reference='EDNiX', post_treatment_method='Grandjean',
                    band='0.01 0.1', blur=0, do_not_correct_signal = False, extract_exterior_CSF = False, extract_WM=True, extract_Vc = False, extract_GS = False,
                    use_erode_WM_func_masks = True, use_erode_V_func_masks=True, normalize='Skip',
                    selected_atlases_matrix='all', wanted_level_matrix='all',
                    selected_atlases_SBA='default', panda_files_SBA='default',
                    SBAspace=['atlas'], erod_seed=True, smoothSBA=smoothSBA,
                    specific_roi_tresh=0.2, delta_thresh=0.1,
                    oversample_map=False, use_cortical_mask_func=False, n_cut=10, threshold_val=10)