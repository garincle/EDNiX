import os
import subprocess
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
sys.path.append(opj('/home/cgarin/PycharmProjects/EasyMRIbrain_sing/'))

import analyses._Group_anal__func_DicLearn
import analyses._Group_anal_3dTtest
import analyses._Group_anal_3dLMEr_SBA

s_bind = ' --bind ' + opj('/', 'scratch', 'cgarin/') + ',' + MAIN_PATH
s_path = opj(MAIN_PATH, 'code', 'singularity')
afni_sif    = ' ' + opj(s_path , 'afni_make_build_AFNI_23.1.10.sif') + ' '

##### where is EasyMRI_brain?
sys.path.append(opj(MAIN_PATH,'code','EasyMRI_brain-master'))
import fonctions._0_Pipeline_launcher

species = 'CatinDog'
bids_dir = opj('/','srv','projects','easymribrain','data', 'MRI', 'Dog', 'BIDS_k9')
layout = BIDSLayout(bids_dir, validate=True)
print(layout)
subject = layout.get_subjects()
print(subject)
tasks = layout.get_tasks()
print(tasks)

###report
report = BIDSReport(layout)
# counter = report.generate()
# main_report = counter.most_common()[0][0]
# print(main_report)

# Ask get() to return the ids of subjects that have T1w files #return_type='filename
T1 = layout.get(return_type='filename', target='subject', suffix='T1w', extension='nii.gz')
print(T1)
###question

# Ask get() to return the ids of subjects that have T2w files
# T2 = layout.get(return_type='filename', target='subject', suffix='T2w', extension='nii.gz')
# print(T2)

# Ask get() to return the ids of subjects that have T1w files
Bold = layout.get(return_type='filename', target='subject', suffix='bold', extension='nii.gz')

# Convert the layout to a pandas dataframe
df = layout.to_df()
df.head()

ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii.gz' # string

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
alpha = 0.001
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
'Amygdala','Hypothalamus',
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

analyses._Group_anal_3dTtest._3dttest_EDNiX(bids_dir, BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, cut_coordsZ, panda_files, selected_atlases,
              lower_cutoff, upper_cutoff, s_bind, afni_sif, alpha, all_ID, all_Session, all_data_path, max_sessionlist, endfmri, mean_imgs, ntimepoint_treshold)

analyses._Group_anal__func_DicLearn.dicstat(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, alpha_dic, component_list,
             cut_coordsZ, bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, afni_sif, s_bind, templatelow, templatehigh)
