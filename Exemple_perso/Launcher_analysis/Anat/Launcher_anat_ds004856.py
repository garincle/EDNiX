import os
from bids import BIDSLayout
from bids.reports import BIDSReport
import pandas as pd
opn = os.path.normpath
opj = os.path.join
MAIN_PATH = opj('/srv/projects/easymribrain/code/EDNiX/')
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import Exemples.Study_Zhu_Garin.old.volstat

species = 'Human'
# Override os.path.join to always return Linux-style paths
bids_dir = Tools.Load_subject_with_BIDS.linux_path(opj("/scratch2/EDNiX/Human/ds004856/"))  # Will be replaced from OLD
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
FreeSlabel_ctab = config["lookup_tables"]["FreeSlabel_ctab_list"]

########### Subject loader with BIDS##############
layout= BIDSLayout(bids_dir,  validate=True)
report = BIDSReport(layout)
df = layout.to_df()
df.head()

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = df[(df['suffix'] == 'T1w') & (df['extension'] == '.nii.gz')]

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = []
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

# Load data (assuming allinfo_study_c is already loaded)
participants_path = bids_dir + "/participants.tsv"
participants_df = pd.read_csv(participants_path, sep='\t')

# First convert subject IDs to strings in both datasets
participants_df['numeric_id'] = participants_df['participant_id'].str.extract(r'sub-(\d+)').astype(str)
allinfo_study_c['subject_str'] = allinfo_study_c['subject'].astype(str)

# Rebuild age_map with string keys
age_map = {}
for _, row in participants_df.iterrows():
    subj_id = row['numeric_id']
    age_data = {
        1: row.get('AgeMRI_W1'),
        2: row.get('AgeMRI_W2'),
        3: row.get('AgeMRI_W3')
    }
    age_map[subj_id] = age_data

# Convert session numbers to integers explicitly
allinfo_study_c['age'] = allinfo_study_c.apply(
    lambda row: age_map.get(str(row['subject']), {}).get(int(row['session'])),
    axis=1)
# Verification
print(f"Added age to {allinfo_study_c['age'].notna().sum()} rows")
print("Sample ages:")
print(allinfo_study_c[['subject', 'session', 'age']].head(10))

regressor_list = ['age']
segmentation_name_list = [lvl1, lvl1LR, lvl2, lvl2LR,
    lvl3, lvl3LR, lvl4, lvl4LR]
segmentation_ID_list = ['lvl1', 'lvl1LR', 'lvl2', 'lvl2LR',
    'lvl3', 'lvl3LR', 'lvl4', 'lvl4LR']

atlas_names_Seg_list = [os.path.basename(path) for path in list_atlases]
type_norm = 'T1w' # T1 or T2

Exemples.Stat_study_Zhu_Garin.old.volstat.extractVol(MAIN_PATH, FS_dir, allinfo_study_c, regressor_list, all_ID, all_Session, all_data_path, type_norm,
                                                     segmentation_name_list, segmentation_ID_list,
                                                     atlas_names_Seg_list, list_atlases, bids_dir)
