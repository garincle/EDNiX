import os
import pandas as pd
from bids import BIDSLayout
from bids.reports import BIDSReport
opn = os.path.normpath
opj = os.path.join
MAIN_PATH = opj('/srv/projects/easymribrain/code/EDNiX/')
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import analyses._Group_anal__func_DicLearn
import analyses._Group_anal_3dTtest
import analyses._Group_anal_3dLMEr_Mirror
import Launcher_analysis.Load_BIDS_data_for_analysis

species = 'CatinDog'
## linux ##
# Override os.path.join to always return Linux-style paths
bids_dir = '/srv/projects/easymribrain/data/MRI/Dog/BIDS_k9'
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
layout= BIDSLayout(bids_dir,  validate=True)
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
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = (
    Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove))


BASE_SS = config["paths"]["BASE_SS"]
BASE_mask = config["paths"]["BASE_mask"]
GM_mask = config["paths"]["GM_mask"]
Aseg_ref = config["paths"]["Aseg_ref"]
Aseg_refLR = config["paths"]["Aseg_refLR"]

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
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = (
    Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove))

# Extract all atlas definitions as DataFrames and get their file paths
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

# for the seed base analysis, you need to provide the names and the labels of the regions you want to use as "seeds"
selected_atlases = [lvl3LR_file, lvl4LR_file]  # Your selected atlases for SBA
panda_files = [lvl3LR, lvl4LR]  # DataFrames for levels 3 and 4

import pandas as pd
selected_atlases = ['atlaslvl4_LR.nii.gz', 'registered_bodiesfaces_mask.nii', 'registered_bodiesfaces_LR.nii']  # Your selected atlases for SBA
panda_files = [pd.DataFrame({'region': ['retrosplenial'], 'label': [162]}), pd.DataFrame({'region': ['bodies-faces'], 'label': [1]}),
    pd.DataFrame({'region': ['R_bodies-faces', 'L_bodies-faces'], 'label': [2, 1]})]

ntimepoint_treshold = 100
endfmri = '*_task-rest_*.nii.gz' # string
images_dir, all_ID, all_Session, all_data_path, max_sessionlist, mean_imgs, templatelow =  Launcher_analysis.Load_BIDS_data_for_analysis.load_data(bids_dir, df, ntimepoint_treshold, list_to_keep, list_to_remove, endfmri)

oversample_map = True
oversample_dictionary = False
min_size = 10
cut_coords = 10
alpha = 0.0001
alpha_dic = 11
component_list = [7, 17]
lower_cutoff = 0.1
upper_cutoff = 0.95

treshold_or_stat = 'stat'
study_template_atlas_forlder = bids_dir + '/sty_template'
folder_atlases = opj(study_template_atlas_forlder, 'atlases/') # sting
mask_func     = opj(folder_atlases, 'Gmask.nii.gz') # sting
type_norm = 'T1w' # T1 or T2
#templatehigh = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz') # sting
templatehigh = opj(study_template_atlas_forlder, 'studytemplate2_T1w/study_template.nii.gz')
TR = '1'  # 'value du calculate in s', 'Auto', 'None'
smoothing = 3.5

# Define model and GLT specifications
model = "(1|Subj)*Hemisphere"
glt_spec = [(1, 'Overall_Effect', 'Hemisphere : 0.5*L +0.5*R'),
    (2, 'Hemisphere_Effect', 'Hemisphere : 1*L -1*R')]
contrast_names = ['Overall_Effect', 'Hemisphere_Effect']
midline_x=1
visualize = 'percentile'
percent = 10

templatehigh = opj(r'/mnt/c/Users/cgarin/Documents/EDNiX/Atlas_library/Atlases_V2/CatinDog/template_SS.nii.gz')
#analyses._Group_anal__func_DicLearn.dicstat(oversample_map, mask_func, cut_coords, alpha_dic, component_list, oversample_dictionary,
#              bids_dir, images_dir, mean_imgs, min_size, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, templatelow, templatehigh, TR, smoothing)

templatehigh = opj(study_template_atlas_forlder, 'studytemplate2_T1w/study_template.nii.gz')
analyses._Group_anal_3dTtest._3dttest_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords, panda_files, selected_atlases,
              lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha ,all_ID, all_Session, all_data_path, endfmri, mean_imgs, ntimepoint_treshold)

analyses._Group_anal_3dLMEr_Mirror._3dLMEr_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords,
                  panda_files, selected_atlases, lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha,
                  all_ID, all_Session, all_data_path, endfmri, mean_imgs,
                  ntimepoint_treshold, model, glt_spec, contrast_names, midline_x, visualize, percent)