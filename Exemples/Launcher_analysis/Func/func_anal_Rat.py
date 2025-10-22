import pandas as pd
import os
from bids import BIDSLayout
opn = os.path.normpath
opj = os.path.join
MAIN_PATH = opj('/srv/projects/easymribrain/code/EDNiX/')
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import analyses.Group_fMRI._Group_anal__func_DicLearn
import analyses.Group_fMRI._Group_anal_3dLMEr_Mirror
import Tools.Load_BIDS_data_for_analysis
species = 'RatWHS'

## linux ##
# Override os.path.join to always return Linux-style paths
bids_dir = '/srv/projects/easymribrain/data/MRI/Rat/BIDS_Gd_shared/'
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')
atlas_dir = opj(MAIN_PATH, "Atlas_library", "Atlases_V2", species)
Lut_dir = opj(MAIN_PATH, "Atlas_library", "LUT_files")

# Define your path variables
path_vars = {'FS_dir': FS_dir,
    'atlas_dir': atlas_dir,
    'Lut_dir': Lut_dir}
# Load and process config.
config_file_path = opj(MAIN_PATH, "Atlas_library", "Atlases_V2", "atlas_config_V2.json")
config = Tools.Read_atlas.load_config(Tools.Load_subject_with_BIDS.linux_path(config_file_path), path_vars)

BASE_SS = config["paths"]["BASE_SS"]
BASE_mask = config["paths"]["BASE_mask"]
GM_mask = config["paths"]["GM_mask"]
Aseg_ref = config["paths"]["Aseg_ref"]
Aseg_refLR = config["paths"]["Aseg_refLR"]

########### Subject loader with BIDS##############
layout= BIDSLayout(bids_dir,  validate=False)
df = layout.to_df()
df.head()

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = df[(df['task'] == 'rest') & (df['extension'] == '.nii.gz')]
list_of_ones = [1] * len(allinfo_study_c)
allinfo_study_c['session'] = list_of_ones

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = [
        ('300301', '1'),
        ('300302', '1'),
        ('300303', '1'),
        ('300304', '1'),
        ('300305', '1'),
        ('300306', '1'),
        ('300307', '1'),
        ('300308', '1'),
        ('300309', '1'),
        ('300600', '1'),
        ('300601', '1'),
        ('300602', '1'),
        ('300603', '1'),
        ('300604', '1'),
        ('300605', '1'),
        ('300606', '1'),
        ('300607', '1'),
        ('300608', '1'),
        ('300609', '1'),
        ('300800', '1'),
        ('300801', '1'),
        ('300802', '1'),
        ('300803', '1'),
        ('300804', '1'),
        ('300805', '1'),
        ('300806', '1'),
        ('300807', '1'),
        ('300808', '1'),
        ('300809', '1')]


all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

atlas_dfs = Tools.Read_atlas.extract_atlas_definitions(config)
(lvl1, lvl1LR, lvl2, lvl2LR,
    lvl3, lvl3LR, lvl4, lvl4LR) = (
    atlas_dfs['lvl1'], atlas_dfs['lvl1LR'],
    atlas_dfs['lvl2'], atlas_dfs['lvl2LR'],
    atlas_dfs['lvl3'], atlas_dfs['lvl3LR'],
    atlas_dfs['lvl4'], atlas_dfs['lvl4LR'])
# Get all atlas file paths
(lvl1_file, lvl1LR_file, lvl2_file, lvl2LR_file,
    lvl3_file, lvl3LR_file, lvl4_file, lvl4LR_file) = (
    config["atlas_definitions"]["lvl1"]["atlas_file"],
    config["atlas_definitions"]["lvl1LR"]["atlas_file"],
    config["atlas_definitions"]["lvl2"]["atlas_file"],
    config["atlas_definitions"]["lvl2LR"]["atlas_file"],
    config["atlas_definitions"]["lvl3"]["atlas_file"],
    config["atlas_definitions"]["lvl3LR"]["atlas_file"],
    config["atlas_definitions"]["lvl4"]["atlas_file"],
    config["atlas_definitions"]["lvl4LR"]["atlas_file"])

# Create combined lists
list_atlases = [lvl1_file, lvl2_file, lvl3_file, lvl4_file,
    lvl1LR_file, lvl2LR_file, lvl3LR_file, lvl4LR_file]

overwrite_option = True #True or False overwrite previous analysis if in BIDS

#### functional images paramters definition
Slice_timing_info = 'Auto'
correction_direction = 'y' # 'x', 'x-', 'y', 'y-', 'Auto', 'None'
### Dwell Time (necessery only if you want to appply TOPUP)
DwellT = 'None' # 'value du calculate', 'Auto', 'None'
### TotalReadoutTime (necessery only if you want to appply TOPUP)
TRT = 'None'
### Slice encoding direction (SED) (necessery only if you want to restrict the transfo for anat to func)
SED = 'Auto' #  "i", "i-", "j", "j-", "k", "k-", 'Auto', 'None'
### YOU NEED TO PROVIDE A TR if not in .json, otherwise it will fail
TR = 'Auto'  # 'value du calculate in s', 'Auto', 'None'

#### fMRI pre-treatment

endfmri = '*_task-rest_bold.nii.gz' # string
endjson = '*_task-rest_bold.json' # string
endmap = '*_map.nii.gz' # string
orientation = 'RIP' # string
deoblique='WARP_without_3drefit' #header or WARP

## prior anatomical processing
coregistration_longitudinal = False #True or False
type_norm = 'T2w' # T1 or T2
### co-registration func to anat to template to with T1 ? T2? use the correct  suffix as in the BIDS
TfMRI = 'T2w' # string
### if you don't have any anatomical image you will need to put several image in the folderforTemplate_Anat (refer to the doc)
folderforTemplate_Anat = ''

## masking
doMaskingfMRI = True # True or False
Method_mask_func = 'nilearn' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"
costAllin = 'ls' # string

#### ANTs function of the co-registration HammingWindowedSinc is advised
IhaveanANAT = True # True or False
anat_func_same_space = True # True or False
use_master_for_Allineate = False
n_for_ANTS = 'hammingWindowedSinc' # string
registration_fast = False
type_of_transform = 'BOLDAffine'
aff_metric_ants_Transl = 'mattes' # string
aff_metric_ants = 'mattes'
do_anat_to_func = False # True or False

##### if you don't have an anat then template will be the same as anat...
#creat_study_template was created with the anat type_norm img, and you want to use it as standart space
creat_study_template = False # True or False
#folder where you stored the stdy template
study_template_atlas_forlder = bids_dir + '/sty_template'
# sting
stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_mask.nii.gz') # sting
stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz') # sting
GM_mask_studyT = opj(study_template_atlas_forlder, 'atlases', 'Gmask.nii.gz') # sting

#######for melodic cleaning (step 6)
ICA_cleaning = 'Skip'
nb_ICA_run = 20 # int
nb_ICA     = 30 # int

#######for 3dDeconvolve cleaning (step 7)
#If you don't want to correct the signal at all!  do_not_correct_signal  = True
do_not_correct_signal  = False # True or False
### you can use the CSF (first layer outside the brain) as regressor
extract_exterior_CSF = False # True or False
### you can use the White Matter as regressor
extract_WM = True # True or False
#use the eroded  White Matter functional mask (produced during the anat processing)
use_erode_WM_func_masks  = True # True or False
### you can use the Ventricules as regressor (not advised for small species as often not enough voxels)
extract_Vc = False # True or False
#use the eroded ventricular functional mask (produced during the anat processing)
use_erode_V_func_masks = False # True or False
#Global signal regression ?
extract_GS = False # True or False

### Band path filtering
band = '0.01 0.1' # string
normalize = 'Skip'
#Smooth
blur = 0.5 # float
#Dilate the functional brain mask by n layers
dilate_mask = 2 # int
#retrain the analysis to the gray matter instate of the brain
use_cortical_mask_func = False # True or False

#### coordinate of the template plot in list form, each number will be a slice (plotting.plot_stat_map = cut_coords)
cut_coordsX = [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6] #list of int
cut_coordsY = [-7, -6, -5, -3, -2, 0, 1, 3, 4, 5] #list of int
cut_coordsZ = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8] #list of int

SBAspace = ['func', 'atlas'] #list containing at least on of the string 'func', 'anat', 'atlas'
erod_seed  = True
smoothSBA = 0

#######for matrix analysis (step 11)
#### name of the atlases  you want to use for the matrix analysis
selected_atlases_matrix = list_atlases.copy()
segmentation_name_list = [lvl1, lvl2, lvl3, lvl4, lvl1LR, lvl2LR, lvl3LR, lvl4LR]

#######for seed analysis (step 11)
# threshold_val is the percentage of the correlation image that will be removed
threshold_val = 10 # int
##use high quality anat image as background for figures
oversample_map = True # True or False
# for the seed base analysis, you need to provide the names and the labels of the regions you want to use as "seeds"
selected_atlases = ['atlaslvl4_LR.nii.gz']  # Your selected atlases for SBA
panda_files = [pd.DataFrame({'region':['retrosplenial'],'label':[162]})]  # DataFrames for levels 3 and 4

#For QC value to define specific and non-spe correlation
specific_roi_tresh = 0.2
delta_thresh = 0.1

#### specific to analysis ###

oversample_map = True
oversample_dictionary = False
min_size = 10
cut_coords = 10
alpha = 0.0001
alpha_dic = 11
component_list = [6]
lower_cutoff = 0.1
upper_cutoff = 0.95
network_quality =  'all'
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max, mean_imgs, images_dir =  Tools.Load_BIDS_data_for_analysis.reverse_load_data_bids(bids_dir, network_quality, file_pattern="residual_in_template.nii.gz")

treshold_or_stat = 'stat'
templatehigh = '/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/RatWHS/templateT2_SS.nii.gz'
mask_func     = opj('/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/RatWHS/', 'Gmask.nii.gz') # sting
templatelow = '/srv/projects/easymribrain/data/MRI/Rat/BASE_SS_fMRI.nii.gz'
artefactpt = ''

TR = None  # 'value du calculate in s', 'Auto', 'None'
smoothing = 0.5
# Define model and GLT specifications
model = "(1|Subj)*Hemisphere"
glt_spec = [(1, 'Overall_Effect', 'Hemisphere : 0.5*L +0.5*R'),
    (2, 'Hemisphere_Effect', 'Hemisphere : 1*L -1*R')]
contrast_names = ['Overall_Effect', 'Hemisphere_Effect']
midline_x=0
visualize = 'percentile'
percent = 10
ratio_n_voxel = 0.5


analyses.Group_fMRI._Group_anal__func_DicLearn.dicstat(oversample_map, mask_func, cut_coords, alpha_dic, component_list, oversample_dictionary,
                                                       bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, templatelow, templatehigh, TR, smoothing, ratio_n_voxel)

analyses.Group_fMRI._Group_anal_3dTtest._3dttest_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords, panda_files, selected_atlases,
                                                       lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha, all_ID, all_Session, all_data_path, endfmri, mean_imgs)

analyses.Group_fMRI._Group_anal_3dLMEr_Mirror._3dLMEr_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords,
                                                            panda_files, selected_atlases, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha,
                                                            all_ID, all_Session, all_data_path, endfmri, mean_imgs,
                                                            model, glt_spec, contrast_names, midline_x, visualize, percent)