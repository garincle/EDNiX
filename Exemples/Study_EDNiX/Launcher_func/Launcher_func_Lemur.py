import pandas as pd
import os
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport
opn = os.path.normpath
opj = os.path.join
MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import fonctions._0_Pipeline_launcher
from anatomical import set_launcher
from Tools import Load_subject_with_BIDS

bids_dir = Load_subject_with_BIDS.linux_path(opj('/srv/projects/easymribrain/scratch/Mouse_lemur/BIDS_Garin/'))

# which format ?
BIDStype = 1

########### Subject loader with BIDS  ##############
layout = BIDSLayout(bids_dir, validate=False)
report = BIDSReport(layout)
df = layout.to_df()
df.head()

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = df[(df['suffix'] == 'T2w') & (df['extension'] == '.nii.gz')]

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# choose if you want to select or remove ID from you analysis:
list_to_keep   = []
list_to_remove = []

species    = 'Mouselemur'
reference  = 'EDNiX'
addatlas   = ''
# choose if you want to select or remove ID from you analysis
list_to_keep = []
#####################################################################################XX should be 01 check that it work
list_to_remove = [('300BA', '2'), ('314CA', '1')]
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

overwrite_option = True #True or False overwrite previous analysis if in BIDS

#### functional images paramters definition
Slice_timing_info = '-tpattern alt+z'
correction_direction = 'Auto' # 'x', 'x-', 'y', 'y-', 'Auto', 'None'
### Dwell Time (necessery only if you want to appply TOPUP)
DwellT = 'Auto' # 'value du calculate', 'Auto', 'None'
### TotalReadoutTime (necessery only if you want to appply TOPUP)
TRT = 'Auto'
### Slice encoding direction (SED) (necessery only if you want to restrict the transfo for anat to func)
SED = 'Auto' #  "i", "i-", "j", "j-", "k", "k-", 'Auto', 'None'
### YOU NEED TO PROVIDE A TR if not in .json, otherwise it will fail
TR = 'Auto'  # 'value du calculate in s', 'Auto', 'None'

#### fMRI pre-treatment
T1_eq = 5 # int
REF_int = 0 # int
ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii.gz' # string
endjson = '*_task-rest_*.json' # string
endmap = '*_map.nii.gz' # string
humanPosition     = ['']
animalPosition    = [''] # valid only for species smaller than humans
orientation = 'LSP' # string
deoblique='header' #header or WARP

## prior anatomical processing
coregistration_longitudinal = False #True or False
type_norm = 'T2w' # T1 or T2
### co-registration func to anat to template to with T1 ? T2? use the correct  suffix as in the BIDS
TfMRI = 'T2w' # string
### if you don't have any anatomical image you will need to put several image in the folderforTemplate_Anat (refer to the doc)
folderforTemplate_Anat = ''

## masking
doMaskingfMRI = True # True or False
Method_mask_func = '3dSkullStrip_marmoset' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"
costAllin = '' # string

#### ANTs function of the co-registration HammingWindowedSinc is advised
IhaveanANAT = True # True or False
anat_func_same_space = False # True or False
use_master_for_Allineate = False
n_for_ANTS = 'hammingWindowedSinc' # string
registration_fast = False
type_of_transform = 'BOLDAffine'
aff_metric_ants_Transl = 'mattes' # string
aff_metric_ants = 'meansquares'
do_anat_to_func = True # True or False

##### if you don't have an anat then template will be the same as anat...
#creat_study_template was created with the anat type_norm img, and you want to use it as standart space
creat_study_template = True # True or False
#folder where you stored the stdy template
study_template_atlas_forlder = '/scratch/cgarin/Mouse_lemur/BIDS_Garin/sty_template/' # sting
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
blur = 0 # float
#Dilate the functional brain mask by n layers
dilate_mask = 0 # int
#retrain the analysis to the gray matter instate of the brain
use_cortical_mask_func = False # True or False

#### coordinate of the template plot in list form, each number will be a slice (plotting.plot_stat_map = cut_coords)
cut_coordsX = [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6] #list of int
cut_coordsY = [-7, -6, -5, -3, -2, 0, 1, 3, 4, 5] #list of int
cut_coordsZ = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8] #list of int

SBAspace = ['func'] #list containing at least on of the string 'func', 'anat', 'atlas'
erod_seed  = True
smoothSBA = 0.5

atlas_followers =  [[], [], [], []]
path_ATLAS = "/home/cgarin/PycharmProjects/Atlases_library/"
(FS_refs, template_dir, reference,balsa_folder, BALSAname, balsa_brainT1,BASE_atlas_folder, BASE_template, BASE_SS,
 BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask,CSF, GM, WM, Aseg_ref,list_atlas, path_label_code,all_ID,
 all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
 fs_tools,reftemplate_path,MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH,bids_dir,allinfo_study_c,species,list_to_keep,
                                                               list_to_remove,reference,type_norm,'', '', atlas_followers)

#######for matrix analysis (step 10)
#### name of the atlases  you want to use for the matrix analysis
selected_atlases_matrix = list_atlas.copy()
segmentation_name_list = []

#######for seed analysis (step 11)
# threshold_val is the percentage of the correlation image that will be removed
threshold_val = 10 # int
##use high quality anat image as background for figures
oversample_map = False # True or False
# for the seed base analysis, you need to provide the names and the labels of the regions you want to use as "seeds"
selected_atlases = [['EDNIxCSC', 3]]  # Using NEW VERSION format (single atlas)
panda_files = [pd.DataFrame({'region':['retrosplenial'],'label':[162]})]  # Using NEW VERSION format (single DataFrame)
doWARPonfunc = True
reference = 'EDNiX'
resting_or_task = 'resting'  # 'resting' or 'task'
#For QC value to define specific and non-spe correlation
specific_roi_tresh = 0.2
delta_thresh = 0.1
post_treatment_method = 'AFNI'
Skip_step = [1,2,3,4,5,6,7,8,12,13,14,15,16,100,200]
fonctions._0_Pipeline_launcher.preprocess_data(species, all_ID, all_Session, all_data_path, all_Session_max,
                    BASE_SS, BASE_mask, T1_eq, Slice_timing_info, anat_func_same_space,
                    correction_direction, REF_int, SBAspace, erod_seed, smoothSBA, deoblique, orientation,
                    TfMRI, GM_mask_studyT, GM, creat_study_template, type_norm, coregistration_longitudinal,
                    dilate_mask, overwrite_option, nb_ICA_run, blur, ICA_cleaning, extract_exterior_CSF, extract_WM,
                    n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl, list_atlas, selected_atlases, panda_files,
                    endfmri, endjson, endmap, oversample_map, use_cortical_mask_func, cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val,
                    Skip_step, bids_dir, costAllin, use_erode_WM_func_masks, do_not_correct_signal, use_erode_V_func_masks,
                    folderforTemplate_Anat, IhaveanANAT, do_anat_to_func, Method_mask_func, segmentation_name_list,
                    band, animalPosition, humanPosition, doWARPonfunc, resting_or_task,
                    extract_Vc, selected_atlases_matrix, specific_roi_tresh, delta_thresh, extract_GS, MAIN_PATH,
                    DwellT, SED, TR, TRT, type_of_transform, ntimepoint_treshold, registration_fast, normalize,
                    reftemplate_path, reference, BASE_atlas_folder, post_treatment_method)