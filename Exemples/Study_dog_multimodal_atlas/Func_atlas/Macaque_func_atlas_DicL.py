import os
import sys
import Tools.Read_atlas
from Statistics.Group_fMRI import _Group_anal__func_DicLearn
from Tools import Load_subject_with_BIDS, load_bids, Load_BIDS_data_for_analysis
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

species    = 'Macaque'
# Override os.path.join to always return Linux-style paths
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Macaque/BIDS_BenHamed/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['anat'], suffixes= ['T1w'], extensions=['.nii.gz'])

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = []

all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii.gz' # string
endjson = '*_task-rest_*.json' # string
endmap = '*_map.nii.gz' # string


oversample_map = True
oversample_dictionary = False
min_size = 10
cut_coords = 10
alpha = 0.0001
alpha_dic = 10
component_list = [7,12,16,20]
lower_cutoff = 0.2
upper_cutoff = 0.85
templatelow = opj('/scratch2/EDNiX/Macaque/BIDS_BenHamed/sub-Elak/ses-1/func/templates/EDNiX/preprocessing/BASE_SS_fMRI.nii.gz')
#templatehigh = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 's648*tudy_template.nii.gz') # sting
templatehigh = opj('/scratch2/EDNiX/Macaque/BIDS_BenHamed/sub-Elak/ses-1/func/templates/EDNiX/preprocessing/BASE_SS.nii.gz')
TR = '2'  # 'value du calculate in s', 'Auto', 'None'
smoothing = 4
mask_func     = opj('/scratch2/EDNiX/Macaque/BIDS_BenHamed/sub-Elak/ses-1/func/templates/EDNiX/masks/Gmask.nii.gz') #
redo = True
#'Spurious', 'Specific', 'No', 'Unspecific'
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max, mean_imgs, images_dir =  Load_BIDS_data_for_analysis.reverse_load_data_bids(bids_dir, 'all', file_pattern="_space-template_desc-fMRI_residual.nii.gz")
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')
ratio_n_voxel = 0.5
print(len(all_ID))
_Group_anal__func_DicLearn.dicstat(oversample_map, mask_func, cut_coords, alpha_dic, component_list, oversample_dictionary,
              bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, templatelow, templatehigh, TR, smoothing, ratio_n_voxel,redo)
