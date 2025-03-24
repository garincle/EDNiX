import os
import subprocess
import pandas as pd
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport
MAIN_PATH = r'/mnt/c/Users/cgarin/Documents/EDNiX'
sys.path.append('/mnt/c/Users/cgarin/PycharmProjects/EDNiX')
import analyses._Group_anal__func_DicLearn
import analyses._Group_anal_3dTtest
import analyses._Group_anal_3dLMEr_SBA
import Launcher_analysis.Load_BIDS_data_for_analysis
import Tools.Load_subject_with_BIDS

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

MAIN_PATH = r'/mnt/c/Users/cgarin/Documents/EDNiX'
sys.path.append('/mnt/c/Users/cgarin/PycharmProjects/EDNiX')

species = 'CatinDog'
# Override os.path.join to always return Linux-style paths
bids_dir = Tools.Load_subject_with_BIDS.linux_path(opj(r"C:\Users\cgarin\Desktop\BIDS_k9"))
FS_dir    = Tools.Load_subject_with_BIDS.linux_path(opj(MAIN_PATH,'FS_Dir_tmp'))
atlas_dir = Tools.Load_subject_with_BIDS.linux_path(opj(r"C:\Users\cgarin\Documents\EDNiX\Atlas_library\Atlases_V2", species))
Lut_dir = Tools.Load_subject_with_BIDS.linux_path(opj(r"C:\Users\cgarin\Documents\EDNiX\Atlas_library\LUT_files"))

########### Subject loader with BIDS##############
layout= BIDSLayout(bids_dir,  validate=True)
###report
report = BIDSReport(layout)
# Ask get() to return the ids of subjects that have T1w files #return_type='filename
T1 = layout.get(return_type='filename', target='subject', suffix='T1w', extension='nii.gz')
print(T1)
# Ask get() to return the ids of subjects that have T1w files
Bold = layout.get(return_type='filename', target='subject', suffix='epi', extension='nii.gz')
# Convert the layout to a pandas dataframe
df = layout.to_df()
df.head()

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = df[(df['suffix'] == 'bold') & (df['extension'] == '.nii.gz')]
list_of_ones = [1] * len(allinfo_study_c)
allinfo_study_c['session'] = list_of_ones

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = []

ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii.gz' # string
images_dir, all_ID, all_Session, all_data_path, max_sessionlist, mean_imgs, templatelow =  Launcher_analysis.Load_BIDS_data_for_analysis.load_data(bids_dir, df, ntimepoint_treshold, list_to_keep, list_to_remove, endfmri)

oversample_map = True
oversample_dictionary = False
min_size = 10
cut_coords = 10
alpha = 0.0001
alpha_dic = 10
component_list = [7, 17]
lower_cutoff = 0.1
upper_cutoff = 0.95

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
'Amygdala','Hypothalamus',
'Thalamus'],'label':[58,59,61,62,64,67,68,71,74,75,76,79,80,81]}), pd.DataFrame({'region':[
'retrosplenial',
'BA 23',
'BA 24',
'BA 32',
'BA 9',
'OB'],'label':[162,128,114,112,107,153]})] # liste of pandas dataframe

treshold_or_stat = 'stat'
study_template_atlas_forlder = bids_dir + '/sty_template'
folder_atlases = opj(study_template_atlas_forlder, 'atlases/') # sting
mask_func     = opj(folder_atlases, 'Gmask.nii.gz') # sting
type_norm = 'T1w' # T1 or T2
#templatehigh = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz') # sting
templatehigh = opj(r'/mnt/c/Users/cgarin/Documents/EDNiX/Atlas_library/Atlases_V2/CatinDog/template_SS.nii.gz')
TR = '1'  # 'value du calculate in s', 'Auto', 'None'
smoothing = 4

analyses._Group_anal__func_DicLearn.dicstat(oversample_map, mask_func, cut_coords, alpha_dic, component_list, oversample_dictionary,
              bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, templatelow, templatehigh, TR, smoothing)

analyses._Group_anal_3dTtest._3dttest_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords, panda_files, selected_atlases,
              lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha ,all_ID, all_Session, all_data_path, endfmri, mean_imgs, ntimepoint_treshold, smoothing)
