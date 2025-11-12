from Statistics.Group_fMRI import _Group_anal__func_DicLearn
import pandas as pd
import os
from bids import BIDSLayout
from bids.reports import BIDSReport
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
from Tools import Load_BIDS_data_for_analysis
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
list_to_remove = [('JACQUELINE', '1')]
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
## prior anatomical processing
coregistration_longitudinal = False #True or False
type_norm = 'T2w' # T1 or T2
TfMRI = 'T2w' # string

## masking
Method_mask_func = 'Vol_sammba_400' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"
creat_study_template = True # True or False

doWARPonfunc = True
resting_or_task = 'resting'  # 'resting' or 'task'
anat_func_same_space = False # True or False
type_of_transform = 'BOLDAffine'
aff_metric_ants_Transl = 'mattes' # string
aff_metric_ants = 'meansquares'
do_anat_to_func = True # True or False
dilate_mask=0
selected_atlases = []  # Using NEW VERSION format (single atlas)

oversample_map = False
oversample_dictionary = False
min_size = 10
cut_coords = 10
alpha = 0.0001
alpha_dic = 9
component_list = [4,7,12,16,20]
lower_cutoff = 0.1
upper_cutoff = 0.95
templatelow = opj('/home/cgarin/Documents/Snake/sub-SILAS/ses-1/func/templates/EDNiX/preprocessing/BASE_SS_fMRI.nii.gz')
#templatehigh = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 's648*tudy_template.nii.gz') # sting
templatehigh = opj('/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/reptiles/Snake/EDNiX/volumes/Snake_space-acpc_desc-template_T2w.nii.gz')
TR = '1'  # 'value du calculate in s', 'Auto', 'None'
smoothing = None
mask_func     = opj('/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/reptiles/Snake/EDNiX/volumes/masks/Snake_desc-Gray_mask.nii.gz') #
redo = False
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max, mean_imgs, images_dir =  Load_BIDS_data_for_analysis.reverse_load_data_bids(bids_dir, 'all', file_pattern="-acpc-template_desc-fMRI_residual.nii.gz")
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')
ratio_n_voxel = 0.5

_Group_anal__func_DicLearn.dicstat(oversample_map, mask_func, cut_coords, alpha_dic, component_list, oversample_dictionary,
              bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, templatelow, templatehigh, TR, smoothing, ratio_n_voxel,redo)
