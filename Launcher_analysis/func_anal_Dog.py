import os
import subprocess
import numpy as np
import pandas as pd
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport
import glob
import nibabel as nib

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
sys.path.append(os.path.join(MAIN_PATH,'code','EasyMRI_brain-master'))
import fonctions
from fonctions.extract_filename import extract_filename
import analyses
import analyses._Groupe_anal__func_DicLearn
import analyses._Groupe_anal_func_network_torch

##############################################################  TO DO !! ##############################################################

s_bind = ' --bind ' + opj('/', 'scratch', 'cgarin/') + ',' + MAIN_PATH
s_path = opj(MAIN_PATH, 'code', 'singularity')
afni_sif    = ' ' + opj(s_path , 'afni_make_build_AFNI_23.1.10.sif') + ' '
# Freesurfer set up
FS_dir    = opj(MAIN_PATH,'FS_Dir_tmp')

##### where is EasyMRI_brain?
sys.path.append(opj(MAIN_PATH,'code','EasyMRI_brain-master'))
import fonctions._0_Pipeline_launcher

species = 'CatinDog'
bids_dir = opj('/scratch/cgarin/','Dog','BIDS_k9')
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
T1 = layout.get(return_type='filename', target='subject', suffix='T1w', extension='nii.gz')
print(T1)
###question


# Ask get() to return the ids of subjects that have T2w files
#T2 = layout.get(return_type='filename', target='subject', suffix='T2w', extension='nii.gz')
#print(T2)

# Ask get() to return the ids of subjects that have T1w files
Bold = layout.get(return_type='filename', target='subject', suffix='bold', extension='nii.gz')

# Ask get() to return the ids of subjects that have T1w files
#topup_dir = layout.get(return_type='filename', target='subject', suffix='epi', extension='nii.gz')

# Convert the layout to a pandas dataframe
df = layout.to_df()
df.head()


allinfo_study_c = df[(df['suffix'] == 'bold') & (df['extension'] == '.nii.gz')]

allinfo_study_c.rename(columns={'session': 'Session'}, inplace=True)
allinfo_study_c.rename(columns={'subject': 'ID'}, inplace=True)
allinfo_study_c.rename(columns={'path': 'DICOMdir'}, inplace=True)

folders = glob.glob(opj(bids_dir,'sub*'))

##############################################################  TO DO !! ##############################################################


allinfo_study_c_formax = allinfo_study_c.copy()

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
mask_func     =opj(diratlas_orig, 'Gmask.nii.gz') # sting

folder_atlases = diratlas_orig

min_size = 10
cut_coordsX = [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6] #list of int
cut_coordsY = [-7, -6, -5, -3, -2, 0, 1, 3, 4, 5] #list of int
cut_coordsZ = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8] #list of int
print(images_dir)
alpha = 10
component_list = [10, 20]
lower_cutoff = 0.1
upper_cutoff = 0.95



analyses._Groupe_anal_func_network_torch.dicstat_torch(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, alpha, component_list,
              cut_coordsZ, bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, afni_sif, s_bind)






#analyses._Groupe_anal__func_DicLearn.dicstat(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, alpha, component_list,
#              cut_coordsZ, bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, afni_sif, s_bind)




'''
#fMRI 3dTest

oversample_map = True
########## if creat_study_template = False ##########
diratlas_orig = opj(MAIN_PATH,'data','Atlas','13_Atlas_project','Atlases_V2', species)
# if creat_study_template== False you need to provide this
BASE_SS     = opj(diratlas_orig, 'template.nii.gz') # sting
BASE_mask   = opj(diratlas_orig, 'brain_mask.nii.gz') # sting
mask_func     =opj(diratlas_orig, 'Gmask.nii.gz') # sting

folder_atlases = diratlas_orig

min_size = 10
cut_coordsX = [-6, -5, -4, -2, -1, 1, 3, 4, 5, 6] #list of int
cut_coordsY = [-7, -6, -5, -3, -2, 0, 1, 3, 4, 5] #list of int
cut_coordsZ = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8] #list of int
print(images_dir)
alpha = 10.5
component_list = [10, 20]
lower_cutoff = 0.1
upper_cutoff = 0.95

output_dir_process = opn('/media/cgarin/Clement_1/4_Chimpanzee/1_fMRI/2_post_process/1_Seed_Base/')
output_dir = '/media/cgarin/Clement_1/4_Chimpanzee/1_fMRI/1_Preproc_Celine_data'
template_filename = opj('/home/cgarin/Documents/8_Multispecies/13_Atlas_project/0_Atlas_modify/Atlas/Chimpanzee/templateTLRC.nii.gz')
atlas_filename = '/home/cgarin/Documents/8_Multispecies/13_Atlas_project/New_atlas_PET/Chimpanzee/PMC_func_reso.nii.gz'
cortical_mask_func = '/home/cgarin/Documents/8_Multispecies/13_Atlas_project/New_atlas_PET/Chimpanzee/cortical_mask_func_reso.nii.gz'
size_or_stat = 'stat'
treshold_or_stat = 70
dolvl1 = False

selected_atlases = ['atlaslvl3_LR.nii.gz', 'atlaslvl4_LR.nii.gz'] #liste

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


analyses.LME_EDNiX(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, panda_files, selected_atlases,
              cut_coordsZ,  cortical_mask_func, bids_dir, mean_imgs, lower_cutoff, upper_cutoff, s_bind, afni_sif, treshold_or_stat,
              seed_to_voxel_correlations_all_fish)




##### do 3dLME #####
panda_disign_matrix = allinfo_study[['ID', regressor, 'Sexe']]
panda_disign_matrix.rename(columns={'ID': 'Subj'}, inplace=True)
panda_disign_matrix.rename(columns={'Sexe': 'Subj'}, inplace=True)
panda_disign_matrix['InputFile'] = list_of_img_to_analyze

filename_disign_matrix = 'disign_matrix.txt'
disign_matrix_txt = file_results + filename_disign_matrix
### remove the disign matrix in case it exists
if os.path.exists(disign_matrix_txt):
    os.remove(disign_matrix_txt)

# creat the path
if not os.path.exists(file_results): os.mkdir(file_results)

if not os.path.exists(disign_matrix_txt):
    panda_disign_matrix.to_csv(disign_matrix_txt, index=None, sep='\t', mode='a')

if os.path.exists(file_results + '3dLME_glt.nii.gz'):
    os.remove(file_results + '3dLME_glt.nii.gz')
    os.remove(file_results + 'resid.nii.gz')
'''


