import os
import sys
import Tools.Read_atlas
from Statistics.Group_fMRI import _Group_anal__func_DicLearn, Select_parameters_DicL, Select_cpmnt_number
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
component_list = [7,12,17]
lower_cutoff = 0.2
upper_cutoff = 0.85
templatelow = opj('/scratch2/EDNiX/Macaque/BIDS_BenHamed/sub-Elak/ses-1/func/templates/EDNiX/preprocessing/BASE_SS_fMRI.nii.gz')
#templatehigh = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 's648*tudy_template.nii.gz') # sting
templatehigh = opj('/scratch2/EDNiX/Macaque/BIDS_BenHamed/sub-Elak/ses-1/func/templates/EDNiX/preprocessing/BASE_SS.nii.gz')
TR = '2'  # 'value du calculate in s', 'Auto', 'None'
smoothing = 3
mask_func     = opj('/scratch2/EDNiX/Macaque/BIDS_BenHamed/sub-Elak/ses-1/func/templates/EDNiX/masks/Gmask.nii.gz') #
redo = True
#'Spurious', 'Specific', 'No', 'Unspecific'
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max, mean_imgs, images_dir =  Load_BIDS_data_for_analysis.reverse_load_data_bids(bids_dir, 'Specific', file_pattern="_space-template_desc-fMRI_residual.nii.gz")
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')
ratio_n_voxel = 0.5
print(len(all_ID))

specific = 'Specific'
#'Spurious', 'Specific', 'No', 'Unspecific'
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max, mean_imgs, images_dir =  Load_BIDS_data_for_analysis.reverse_load_data_bids(bids_dir, specific, file_pattern="_space-template_desc-fMRI_residual.nii.gz")
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')
print(len(all_ID))
#### DL analysis and functional atlas building
alpha_dic = 10
smoothing = 3
method_mask_func = 'onlyprovidedmask'

'''
# Run optimization with stability
results_df, best_params = Select_parameters_DicL.optimize_parameters_with_stability(
    bids_dir=bids_dir,
    images_dir=images_dir,
    mask_func=mask_func,
    component_range=[7, 16, 20],
    smoothing_range=[3, 4],
    alpha_range=[10, 11],
    TR=1.0,
    n_repeats=10,  # 10 split-half repeats per config
    test_size=0.5,
    random_state=42,
    label_specific=specific)

# Usage example

# Define parameters
specificity_groups = ['Specific', 'Unspecific', 'Spurious', 'No', 'all']
component_list = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
alpha_dic = 10
smoothing = 3
output_dir = os.path.join(bids_dir, 'DicL', 'specificity_comparison')

# Run comparison
results = Select_cpmnt_number.compare_specificity_groups(
    bids_dir=bids_dir,
    specificity_groups=specificity_groups,
    component_list=component_list,
    alpha_dic=alpha_dic,
    smoothing=smoothing,
    TR=TR,
    mask_func=mask_func,
    n_repeats=10,
    output_dir=output_dir)

'''
_Group_anal__func_DicLearn.dicstat(oversample_map, mask_func, cut_coords, alpha_dic, component_list, oversample_dictionary, bids_dir, images_dir, mean_imgs, min_size, lower_cutoff,
            upper_cutoff, MAIN_PATH, templatelow, templatehigh, TR, smoothing, redo, method_mask_func, specific)

