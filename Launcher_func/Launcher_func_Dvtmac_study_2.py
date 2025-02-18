import os
import subprocess
import numpy as np
import pandas as pd
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

##############################################################  TO DO !! ##############################################################

MAIN_PATH = opj('/','srv','projects','easymribrain')
s_bind = ' --bind ' + opj('/', 'scratch', 'cgarin/') + ',' + MAIN_PATH
s_path = opj(MAIN_PATH, 'code', 'singularity')
# Freesurfer set up
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')
sys.path.append(opj(MAIN_PATH,'code','EasyMRI_brain-master'))
import fonctions._0_Pipeline_launcher


species = 'Macaque'
bids_dir = opj('/scratch/cgarin/Macaque/BIDS_Cdt_Garin')
##########################################
########### Subject loader################
##########################################

layout= BIDSLayout(bids_dir,  validate=True)
print(layout)
subject = layout.get_subjects()
print(subject)
tasks = layout.get_tasks()
print(tasks)

###report
report = BIDSReport(layout)
#counter = report.generate()
#main_report = counter.most_common()[0][0]
#print(main_report)


# Ask get() to return the ids of subjects that have T1w files #return_type='filename
T1 = layout.get(return_type='filename', target='subject', suffix='T1', extension='nii.gz')
print(T1)
###question


# Ask get() to return the ids of subjects that have T2w files
#T2 = layout.get(return_type='filename', target='subject', suffix='T2w', extension='nii.gz')
#print(T2)

# Ask get() to return the ids of subjects that have T1w files
Bold = layout.get(return_type='filename', target='subject', suffix='bold', extension='nii')

# Ask get() to return the ids of subjects that have T1w files
#topup_dir = layout.get(return_type='filename', target='subject', suffix='epi', extension='nii.gz')

# Convert the layout to a pandas dataframe
df = layout.to_df()
df.head()

##############################################################  TO DO !! ##############################################################
#### Create a pandas sheet for the dataset (I like it, it help to know what you are about to process
allinfo_study_c = df[(df['suffix'] == 'bold') & (df['extension'] == '.nii')]
allinfo_study_c.rename(columns={'session': 'Session'}, inplace=True)
allinfo_study_c.rename(columns={'subject': 'ID'}, inplace=True)
allinfo_study_c.rename(columns={'path': 'DICOMdir'}, inplace=True)

##############Select all the monkey of the study
filter1 = allinfo_study_c["Session"].isin([95])
allinfo_study_c_formax = allinfo_study_c[~filter1]

#########creat lists of indiv and usefull variables
all_data_path = []
all_ID = []
all_Session = []
all_data_path_max = []
all_ID_max = []
max_session = []
max_sessionlist = []

for ID in pd.unique(allinfo_study_c_formax.ID):
    list_session = allinfo_study_c_formax.loc[allinfo_study_c_formax['ID'] == ID].Session.dropna()
    listereverse = list(map(int, list_session))  # Convert strings to integers
    listereverse.sort(reverse=True)  # Sort in descending order
    if 95 in listereverse:
        listereverse.remove(95)
    max_session.append(np.array(listereverse).max())

    for Session in listereverse:
        print('session numuber ' + str(Session))
        # Organization of the folders
        data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
        all_data_path.append(data_path)
        all_Session.append(Session)
        all_ID.append(ID)
        max_sessionlist.append(np.array(listereverse).max())

for ID, Session in zip(pd.unique(allinfo_study_c_formax.ID), max_session):
    # Organization of the folders
    data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
    all_data_path_max.append(data_path)
    all_ID_max.append(ID)

######## select animals that have not been analyzed yet
removelist = []
###quantum no T2Flair ses 3
######### select the indiv you want to analyse!!!
for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, max_sessionlist)):
    if ID in ["Trinity"] and Session in [6,3]:
        removelist.append(num)
    elif ID in ["Sonic"] and Session in [4]:
        removelist.append(num)

all_ID =  [item for i, item in enumerate(all_ID) if i in removelist]
all_Session =  [item for i, item in enumerate(all_Session) if i in removelist]
all_data_path =  [item for i, item in enumerate(all_data_path) if i in removelist]
max_sessionlist =  [item for i, item in enumerate(max_sessionlist) if i in removelist]

##############################################################  TO DO !!!! ##############################################################
                                ##############################################################
                                ########### 2. variable specific of you dataset ##############
                                ##############################################################

#### Choose the type of analysis you want to to do
coregistration_longitudinal = True #True or False
overwrite_option = True #True or False overwrite previous analysis if in BIDS

#### if you don't have an anat image!!  not great but sometine you don't have the choice
IhaveanANAT = True # True or False

#####if IhaveanANAT = False
### you will need to put in the folderforTemplate_Anat

folderforTemplate_Anat = ''
anat_subject = opj(folderforTemplate_Anat,'template.nii.gz') # string
brainmask     = opj(folderforTemplate_Anat,'brainmask.nii.gz') # string
V_mask        = opj(folderforTemplate_Anat,'V_mask.nii.gz') # string
W_mask = opj(folderforTemplate_Anat,'W_mask.nii.gz') # string
G_mask = opj(folderforTemplate_Anat,'G_mask.nii.gz') # string

#### find the good fmri image: as it is not always standart, look in you BIDS and help use to know how you fmri dataset end by ?:
### specify the suffix to be used by glob.glob to select all fmri image (or map) in their respective folders
endfmri = '*_task-rest_*.nii' # string
endjson = '*_task-rest_*.json' # string

####find on image in the opposite direction of th BOLD aquistion either one per run or one per session or none !!!
### if the pipeline doesn't find the image it will continue anyway so be carefull!
endmap = '*_map.nii' # string

##### is your anat and func in the same same space ? iff they are you can put anat_func_same_space = True  and it will use the mask of the anat to help
# with the co-registration. It also add other problem, so even if they are in the same space you can put anat_func_same_space = False
anat_func_same_space = False # True or False

### co-registration func to anat to template to with T1 ? T2? use the correct  suffix as in the BIDS
TfMRI = 'T1' # string

#### Specify if you have a T1 and T2 image in the same space
IgotbothT1T2 = False # True or False

#### fMRI pre-treatment
### number of TR to remove at the begining
T1_eq = 5 # int
#### Choose which image you want to use as ref (0 the first one, 1 the second run, ect...)
REF_int = 0 # int
Slice_timing_info = 'Auto'

#### ==> you need to check in the json or on the image map what is the encoding direction
# For exemple if  phase encoding direction is "j-" you should put physical coordinates y- (ijk = xyz) (info: if image orientation is LPI it means that your image has been acquired from P to A)
###################for fudge and anat to func!!!
# from FSL
####--unwarpdir=dir
#specifies the direction of the unwarping/warping - i.e. phase-encode direction - with dir being one of x,y,z,x-,y-,z-).
#Note: the conventions for x, y and z are voxel axes (in the input image), but with the x-axis being swapped if the qform or sform determinant is positive
# (to match the internal FSL voxel coordinate axes), though due to other possible changes in the pipeline for reconstructing/analysing both the EPI
# and fieldmap images, it is usually necessary to determine the sign of the direction by trial-and-error once
# for each scanning/reconstruction/analysis protocol (each dataset using the same protocol will need the same unwarpdir).
### CHATGPT
#To determine the correct unwarp direction:
#Check the Phase Encoding Direction: Determine the direction in which the phase encoding was applied during acquisition.
#This is often specified in the metadata of the image or the acquisition protocol.
#Refer to the Image Orientation: Understand the image orientation and how the voxel axes (i, j, k) map to the physical axes (x, y, z).
#For example, if you know:
#Your image was acquired with phase encoding direction from posterior to anterior (P to A).
#The image orientation is LPA (Left-Posterior-Anterior).
#Then:
#The phase encoding direction corresponds to the y-axis.
#The distortions occurred in the y direction, so you would use --unwarpdir=y.

##### if their is an image in the  opposite direction, then we will apply fudge but for this you need to know
####### so here you need to provide the path for the se_map and b02b0.cnf
#b02b0.cnf will be provided, but I never undestood if it should be really modified

# However, you need to change se map file!!!!!!
# se_map need 4 numbers in two lines:
# First is the line for the map img and second for the func
# In our example, it should be the -1 if the map image is j- and 1 if j (then for the second line it is the same but for the func img).
# The last number is used to correct for different in acquistion between func and fmap.
#I don't know how it works, but refer to doc if necessary (it should be TotalReadoutTime).
#https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide#A--datain
# If same image, it should be the same number (whatever it is)
# Ones you have created the topup file you should save it somewhere

### correction_direction (necessery only if you want to appply TOPUP)
correction_direction = 'Auto' # 'x', 'x-', 'y', 'y-', 'Auto', 'None'

### Dwell Time (necessery only if you want to appply TOPUP)
DwellT = 'Auto' # 'value du calculate', 'Auto', 'None'
### TotalReadoutTime (necessery only if you want to appply TOPUP)
TRT = 'Auto'

#### where are stored the file for topup??
study_fMRI_Refth = opj(MAIN_PATH,'code','4topup.txt') #string (path)

########################################################################################################################

### Slice encoding direction (SED) (necessery only if you want to restrict the transfo for anat to func)
SED = 'Auto' #  "i", "i-", "j", "j-", "k", "k-", 'Auto', 'None'

### YOU NEED TO PROVIDE A TR if not in .json, otherwise it will fail
TR = 'Auto'  # 'value du calculate in s', 'Auto', 'None'
ntimepoint_treshold = 100
##### masking steps SUPER IMPORTANT!!
# you can choose to not do it (not advised)
doMaskingfMRI = True # True or False

# if anat_func_same_space == True
# then it will automatically take the mask of the anat image (be sure that they are in the same space!)

# if no
### Method_mask_func=="3dAllineate" or Method_mask_func = 'nilearn'
#### 3dAllineate is based ont the linerar alignment of the anat to the func to send the anat mask to the func

#### nilearn is a theshold based method (you can play with the threshold level)
Method_mask_func = '3dSkullStrip_monkey' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"

### if Method_mask_func=="3dAllineate" choose a method a alignment
costAllin = '' # string

### if Method_mask_func=="nilearn" choose a cutoff
lower_cutoff = 0.05 # int
upper_cutoff = 0.5 # int

###############################################################################################################
############################################### coregistration steps ##########################################
###############################################################################################################
########### define orientation############
orientation = 'RAI' # string
###############################################################################################################
######################## Probably the most "obscure part of the script" ########################
###############################################################################################################
#### If your anat and func are not in the same space, then, it is easy, just put deoblique='header'
# another option is to put all of several animals in deoblique_exeption1, and no deoblique will be applied

#### if they are, and you want to make your life easy by using the anatomical mask, you need to have the anat and func in the space after deoblique
#### their is two way to do it deoblique='header': if the oblique is not too crazy and the center of the FOV similar it should work, but often, it is not the case.
## deoblique='WARP' BOLD images will be warped to fix the deoblique difference, you can use it ONLY if you have done the same thing of the anatomical images.
# It will have for unfortunate consequence to warp the func multiple times to go in the atlas space and one time in the original space.
# However, this function can help to solve common space problem.....

deoblique='WARP_without_3drefit' #header or WARP

#### ANTs function of the co-registration HammingWindowedSinc is advised
n_for_ANTS = 'hammingWindowedSinc' # string
type_of_transform = 'SyNBoldAff'
aff_metric_ants = 'MI'
registration_fast = False
####Choose to normalize using T1 or T2 or T2w as in you anat file!!!!!
### define the acronyme/suffix of the anat as in the BIDS
type_norm = 'T1' # T1 or T2
### define the acronyme/suffix of the other anat as in the BIDS
otheranat = 'T2FLAIR' # sting

###### sometime, the functional quality is so poor that co-registering the anat to the functional image will creat mistakes
###### if it is the case and if you !!!! FUNC IS IN THE ANAT SPACE !!!!! you may try do_anat_to_func = False,
# it will assume that no coregistration between the anatomical image and the func is necessary
do_anat_to_func = True # True or False

#######################################################################
######################### study template??? ###########################
#######################################################################
##### if you don't have an anat then template will be the same as anat...
#creat_study_template was created with the anat type_norm img, and you want to use it as standart space
creat_study_template = True # True or False

########## if creat_study_template = True ##########

######no need to answer this question if you are not doing a study template
#folder where you stored the stdy template
study_template_atlas_forlder = '/scratch/cgarin/Macaque/BIDS_Cdt_Garin/Study_template_test/'  # sting
stdy_template_mask = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_mask.nii.gz') # sting
stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz') # sting
GM_mask_studyT = opj('/scratch/cgarin/Macaque/BIDS_Cdt_Garin/Study_template_test/Atlases_ref_in_stdy_template/Gmask.nii.gz') # sting

########## if creat_study_template = False ##########
diratlas_orig = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','Atlases_V2', species)

# if creat_study_template== False you need to provide this
BASE_SS     = opj(diratlas_orig, 'template.nii.gz') # sting
BASE_mask   = opj(diratlas_orig, 'brain_mask.nii.gz') # sting
GM_mask     =opj(diratlas_orig, 'Gmask.nii.gz') # sting

    ##########################################################
    ##### define atlases that are in template space ##########
    ##########################################################

####put all atlases and template to process in the same folder named: ...
list_atlases = [opj(diratlas_orig, 'atlaslvl1.nii.gz'),
opj(diratlas_orig, 'atlaslvl2.nii.gz'),
opj(diratlas_orig, 'atlaslvl3.nii.gz'),
opj(diratlas_orig, 'atlaslvl4.nii.gz'),
opj(diratlas_orig, 'atlaslvl1_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl2_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl3_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl4_LR.nii.gz')]

#######for melodic cleaning (step 4)
melodic_prior_post_TTT = False # True or False
nb_ICA_run = 20 # int
nb_ICA     = 30 # int

#######for 3dDeconvolve cleaning (step 5)
#If you dont want to correct the signal at all!  do_not_correct_signal  = True
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
#Smooth
blur = 3 # float
#Dilate the functional brain mask by n layers
dilate_mask = 0 # int
#retrain the analysis to the gray matter
use_cortical_mask_func = False # True or False

#######for seed analysis (step 11)
#### name of the atlases  you want to use for the seed base analysis
selected_atlases = ['atlaslvl3_LR.nii.gz', 'atlaslvl4_LR.nii.gz'] #liste

# for the seed base analysis, you need to provide the names and the labels of the regions you want to use as "seeds"
panda_files = [pd.DataFrame({'region':[
'Somatosensory cortex',
'Posterior parietal cortex',
'Visual pre and extra striate cortex',
'Visual striate cortex',
'Auditory cortex (Superior temporal)',
'Insula and others in lateral sulcus',
'Septum',
'Hippocampal formation',
'Periarchicortex',
'Striatum',
'Basal forebrain',
'Amygdala',
'Hypothalamus',
'Thalamus'],'label':[58,59,61,62,64,67,68,71,74,75,76,79,80,81]}), pd.DataFrame({'region':[
'retrosplenial',
'BA 23',
'BA 24',
'BA 32',
'BA 9',
'OB'],'label':[162,128,114,112,107,153]})] # liste of pandas dataframe


#### coordinate of the template plot in list form, each number will be a slice (plotting.plot_stat_map = cut_coords)
cut_coordsX = [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6] #list of int
cut_coordsY = [-7, -6, -5, -3, -2, 0, 1, 3, 4, 5] #list of int
cut_coordsZ = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8] #list of int

SBAspace = ['func', 'atlas'] #list containing at least on of the string 'func', 'anat', 'atlas'
erod_seed  = True

#Threshold the correlation image np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], threshold_val)
# threshold_val is the percentage of the correlation image that will be removed
threshold_val = 10 # int

##use high quality anat image as background for figures
oversample_map = False # True or False

#######for matrix analysis (step 10)
#### name of the atlases  you want to use for the matrix analysis
# Read the Excel file into a DataFrame
file_path = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','Classiff','Legende_VDualvf2_formatrix.xlsx')
legendPNAS = pd.read_excel(file_path, 'Legend_2023')

selected_atlases_matrix = [opj(diratlas_orig, 'atlaslvl1.nii.gz'),
opj(diratlas_orig, 'atlaslvl2.nii.gz'),
opj(diratlas_orig, 'atlaslvl3.nii.gz'),
opj(diratlas_orig, 'atlaslvl4.nii.gz'),
opj(diratlas_orig, 'atlaslvl1_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl2_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl3_LR.nii.gz'),
opj(diratlas_orig, 'atlaslvl4_LR.nii.gz')]

# Select the desired columns and rename them
pandas1 = legendPNAS[['NEWlvl1_label', 'NEWLVL1']].rename(columns={'NEWlvl1_label': 'label', 'NEWLVL1': 'region'})
pandas1 = pd.DataFrame(data={'label': pandas1['label'].unique(), 'region': pandas1['region'].unique()}).dropna()
pandas1['label'] = pandas1['label'].astype(int)
# Select the desired columns and rename them
pandas2 = legendPNAS[['NEWlvl2_label', 'NEWLVL2']].rename(columns={'NEWlvl2_label': 'label', 'NEWLVL2': 'region'})
pandas2 =  pd.DataFrame(data={'label': pandas2['label'].unique(), 'region': pandas2['region'].unique()}).dropna()
pandas2['label'] = pandas2['label'].astype(int)
# Select the desired columns and rename them
pandas3 = legendPNAS[['NEWlvl3_label', 'NEWLVL3']].rename(columns={'NEWlvl3_label': 'label', 'NEWLVL3': 'region'})
pandas3 =  pd.DataFrame(data={'label': pandas3['label'].unique(), 'region': pandas3['region'].unique()}).dropna()
pandas3['label'] = pandas3['label'].astype(int)
# Select the desired columns and rename them
pandas4 = legendPNAS[['NEWlvl4_label', 'NEWLVL4']].rename(columns={'NEWlvl4_label': 'label', 'NEWLVL4': 'region'})
pandas4 =  pd.DataFrame(data={'label': pandas4['label'].unique(), 'region': pandas4['region'].unique()}).dropna()
pandas4['label'] = pandas4['label'].astype(int)

#### name of the regions and labels  you want to use for the matrix analysis
segmentation_name_list = [pandas1, pandas2, pandas3, pandas4] # liste of pandas dataframe

specific_roi_tresh = 0.4
unspecific_ROI_thresh = 0.2
Seed_name = 'Periarchicortex'

############ Right in a list format the steps that you want to skip
Skip_step = [4,100,200]

    ############################################################
    ######################## START de pipeline #################
    ############################################################

fonctions._0_Pipeline_launcher.preprocess_data(all_ID, all_Session, all_data_path, max_sessionlist, stdy_template, stdy_template_mask, BASE_SS, BASE_mask, T1_eq, Slice_timing_info, anat_func_same_space,
    correction_direction, REF_int, study_fMRI_Refth, SBAspace, erod_seed, deoblique, orientation,
    TfMRI, GM_mask_studyT, GM_mask, creat_study_template, type_norm, coregistration_longitudinal, dilate_mask, overwrite_option, nb_ICA_run, blur, melodic_prior_post_TTT,
    extract_exterior_CSF, extract_WM, n_for_ANTS, aff_metric_ants, list_atlases, selected_atlases, panda_files, endfmri, endjson, endmap, oversample_map, use_cortical_mask_func,
    cut_coordsX, cut_coordsY, cut_coordsZ, threshold_val, Skip_step, bids_dir, costAllin, use_erode_WM_func_masks, do_not_correct_signal, use_erode_V_func_masks,
    folderforTemplate_Anat, IhaveanANAT, doMaskingfMRI, do_anat_to_func, Method_mask_func, segmentation_name_list, band, extract_Vc, lower_cutoff, upper_cutoff, selected_atlases_matrix,
    specific_roi_tresh, unspecific_ROI_thresh, Seed_name, extract_GS, MAIN_PATH, DwellT, SED, TR, TRT, type_of_transform, ntimepoint_treshold, registration_fast, s_bind, s_path)