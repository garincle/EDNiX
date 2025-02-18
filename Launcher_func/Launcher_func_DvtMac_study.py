import os
import subprocess
import numpy as np
import pandas as pd
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport
import glob
import nibabel as nib
from fonctions.extract_filename import extract_filename
import analyses
import analyses._Groupe_anal__func_DicLearn
import analyses._Groupe_anal_func_network_torch
import analyses._Group_anal_3dMEMA
import analyses._Groupe_anal_3dTtest
import analyses._Group_anal_3dLME_SBA

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
    if ID in ["Unity"] and Session in [5]: ### failed! proably header
        removelist.append(num)
    elif ID in ["Quantum"] and Session in [3]: ### data problem
        removelist.append(num)

all_ID =  [item for i, item in enumerate(all_ID) if i not in removelist]
all_Session =  [item for i, item in enumerate(all_Session) if not i in removelist]
all_data_path =  [item for i, item in enumerate(all_data_path) if not i in removelist]
max_sessionlist =  [item for i, item in enumerate(max_sessionlist) if not i in removelist]

############################################################## NOTHING TO DO HERE ##############################################################

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### now, we need to creat the list of variable to feed to function "_0_Pipeline_launcher" #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#########creat lists of indiv and usefull variables
# just a list of all datapath
all_data_path = []
# just a list of all subject ID
all_ID = []
# just a list of all Session
all_Session = []

### for longitudinal datasets, you may want to know which session is the last, this is the way to handle it
### for longitudinal datasets, you may want to know which session is the last, this is the way to handle it
all_data_path_max = []
all_ID_max = []
max_session = []
max_sessionlist = []

### this is just to list all subject in a ID + 'ses-' + str(Session) way (usefull for the variable "deoblique_exeption")
animal_ID = []

#let's add all the the string to those lists
for ID in pd.unique(allinfo_study_c_formax.ID):
    list_session = allinfo_study_c_formax.loc[allinfo_study_c_formax['ID'] == ID].Session.dropna()
    listereverse = list(list_session)
    listereverse.reverse()
    if len(list(list_session))>1:
        max_session.append(np.array(listereverse).max())
    else:
        max_session.append(np.array(listereverse))

    for Session in listereverse:
        print('session numuber ' + str(Session))

        # Organization of the folders
        data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
        all_data_path.append(data_path)
        all_Session.append(Session)
        all_ID.append(ID)
        if len(list(list_session)) > 1:
            max_sessionlist.append(np.array(listereverse).max())
        else:
            max_sessionlist.append(np.array(listereverse))
        animal_ID.append(ID + 'ses-' + str(Session))

for ID, Session in zip(pd.unique(allinfo_study_c_formax.ID), max_session):
    # Organization of the folders
    data_path = opj(bids_dir,'sub-' + ID,'ses-' + str(Session))
    all_data_path_max.append(data_path)
    all_ID_max.append(ID)

##############################################################  TO DO !!!! ##############################################################

######## sometime you aleady analysed some subjects, and you don't want to do it again,
# this can not be done in the previous step, as if you are working with the longitudinal dataset it can mess up your analysis
######## select animals that have been analyzed already

removelist = []
######### select the indiv you want to remove !!!
for num, (ID, Session, data_path, max_ses) in enumerate(zip(all_ID, all_Session, all_data_path, max_sessionlist)):
    if ID in []:
        removelist.append(num)

############################################################## NOTHING TO DO HERE ##############################################################

#### apply this to the already created lists
all_ID =  [item for i, item in enumerate(all_ID) if i not in removelist]
all_Session =  [item for i, item in enumerate(all_Session) if i not in removelist]
all_data_path =  [item for i, item in enumerate(all_data_path) if i not in removelist]
max_sessionlist =  [item for i, item in enumerate(max_sessionlist) if i not in removelist]

ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii.gz' # string

images_dir = []
mean_imgs = []
seed_to_voxel_correlations_all_fish = []
for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):

    dir_fMRI_Refth_RS = opj(data_path, 'func')
    dir_fMRI_Refth_RS_prepro = opj(dir_fMRI_Refth_RS, '01_prepro')
    dir_fMRI_Refth_RS_prepro1 = opj(dir_fMRI_Refth_RS_prepro, '01_funcspace')
    dir_fMRI_Refth_RS_prepro2 = opj(dir_fMRI_Refth_RS_prepro, '02_anatspace')
    dir_fMRI_Refth_RS_prepro3 = opj(dir_fMRI_Refth_RS_prepro, '03_atlas_space')

    # Check the func runs
    list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))
    RS = [os.path.basename(i) for i in list_RS]

    if len(list_RS) == 0:
        nl = 'ERROR : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endfmri) and here it is ' + str(
            opj(dir_fMRI_Refth_RS, endfmri)) + ' I would check how you define "endfmri"'
        raise ValueError(nl)
    list_RS_list = list_RS.copy()
    list_pop_index = []

    for imageF in list_RS_list:
        # Load the fMRI NIfTI image
        fmri_image = nib.load(imageF)
        # Get the shape of the image (x, y, z, t)
        image_shape = fmri_image.shape
        # Check the number of time points (4th dimension)
        if len(image_shape) == 4:
            ntimepoint = image_shape[3]  # The 4th dimension represents time

            if int(ntimepoint) < ntimepoint_treshold:
                index_of_imageF = list_RS.index(imageF)
                list_RS.pop(index_of_imageF)
                list_pop_index.append(index_of_imageF)

    nb_run = len(list_RS)
    # Setup for distortion correction using Fieldmaps

    if ope(opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz')):
        mean_imgs.append(opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz'))


    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])
        ### here we want to select all the image preTTT in a common space
        if ope(opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')):
            images_dir.append(opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz'))

oversample_map = True

########## if creat_study_template = False ##########
diratlas_orig = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','Atlases_V2', species)
# if creat_study_template== False you need to provide this
BASE_SS     = opj(diratlas_orig, 'template.nii.gz') # sting
BASE_mask   = opj(diratlas_orig, 'brain_mask.nii.gz') # sting
mask_func     = opj(diratlas_orig, 'Gmask.nii.gz') # sting

folder_atlases = diratlas_orig

min_size = 10
cut_coordsX = [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6] #list of int
cut_coordsY = [-7, -6, -5, -3, -2, 0, 1, 3, 4, 5] #list of int
cut_coordsZ = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8] #list of int
print(images_dir)
alpha = 0.01
alpha_dic = 9
component_list = [7, 17]
lower_cutoff = 0.1
upper_cutoff = 0.95



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

treshold_or_stat = 'stat'
templatelow = "/srv/projects/easymribrain/data/MRI/Dog/BIDS_k9/sub-28/ses-1/func/01_prepro/03_atlas_space/BASE_SS_fMRI.nii.gz"  # Low-resolution atlas
templatehigh = mask_func   # High-resolution anatomical image

#analyses._Groupe_anal_3dTtest._3dttest_EDNiX(bids_dir, BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, cut_coordsZ, panda_files, selected_atlases,
#              lower_cutoff, upper_cutoff, s_bind, afni_sif, alpha ,all_ID, all_Session, all_data_path, max_sessionlist, endfmri, mean_imgs, ntimepoint_treshold)

#analyses._Group_anal_3dMEMA._3dMEMA_EDNiX(bids_dir, BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, cut_coordsZ, panda_files, selected_atlases,
#              lower_cutoff, upper_cutoff, s_bind, afni_sif, alpha ,all_ID, all_Session, all_data_path, max_sessionlist, endfmri, mean_imgs, ntimepoint_treshold)

#analyses._Groupe_anal__func_DicLearn.dicstat(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, alpha_dic, component_list,
#              cut_coordsZ, bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, afni_sif, s_bind, templatelow, templatehigh)

analyses._Group_anal_3dLME_SBA._3dLME_EDNiX(bids_dir, BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, cut_coordsZ, panda_files, selected_atlases,
              lower_cutoff, upper_cutoff, sresh, Seed_name, extract_GS, MAIN_PATH, DwellT, SED, TR, TRT, type_of_transform, ntimepoint_treshold, registration_fast, s_bind, s_path)