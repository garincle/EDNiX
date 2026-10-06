import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import Exemples.Stat_study_Zhu_Garin.old.volstat

import os
from bids import BIDSLayout
from bids.reports import BIDSReport
opn = os.path.normpath
opj = os.path.join
MAIN_PATH = opj('/srv/projects/easymribrain/code/EDNiX/')
import Tools.Read_atlas

species = 'Macaque'
# Override os.path.join to always return Linux-style paths
bids_dir = Tools.Load_subject_with_BIDS.linux_path(opj('/srv/projects/easymribrain/data/MRI/Macaque/BIDS_Cdt_Garin'))
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
allinfo_study_c = df[(df['suffix'] == 'bold') & (df['extension'] == '.nii.gz')]

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = []
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
#### for 14 (surfaces) ####
list_atlases_2 = [lvl1_file, lvl2_file, lvl3_file, lvl4_file,
    lvl1LR_file, lvl2LR_file, lvl3LR_file, lvl4LR_file]

### file to standardize space
FS_buckner40_TIF = opj(FS_dir,'MacaqueYerkes19')
FS_buckner40_GCS = opj(FS_dir,'MacaqueYerkes19')
FreeSlabel_ctab_list = [FreeSlabel_ctab[0], FreeSlabel_ctab[0], FreeSlabel_ctab[0],FreeSlabel_ctab[0],
                        FreeSlabel_ctab[1], FreeSlabel_ctab[1], FreeSlabel_ctab[1],FreeSlabel_ctab[1]]
Lut_file = opj(FreeSlabel_ctab[0])
coregistration_longitudinal = True
check_visualy_each_img = False
do_manual_crop = False
overwrite_option = True

check_visualy_final_mask = False
deoblique='WARP_without_3drefit'

n_for_ANTS='hammingWindowedSinc'
type_of_transform = 'SyN'
aff_metric_ants_Transl = 'mattes'
aff_metric_ants = 'MI'

####Choose to normalize using T1 or T2
type_norm = 'T1' # T1 or T2
otheranat = 'T2FLAIR'
orientation = 'RAP'
BIDStype = 1

###masking
masking_img = 'T1'
brain_skullstrip_1 ='Custum_Macaque2' # bet2_ANTS or MachinL
#precise
brain_skullstrip_2 ='Custum_QWARP' # bet2_ANTS or MachinL
do_fMRImasks = True
do_fMRImasks = True
fMRImasks = 'aseg' #must be aseg or custom, if custom  please add a ventricle and whitte matter mask in the template space named such as Vmask, Wmask
Align_img_to_template = '3dAllineate' #3dAllineate or No or @Align_Centers
cost3dAllineate = 'hel'

#creat_study_template with type_norm img
creat_study_template = True
#do you want to use all the data or only the last one of each subject (for longitudinal co-registration)
which_on = 'max' # all or max
type_of_transform_stdyT = 'SyN'
aff_metric_ants_Transl_template = 'mattes'
Atemplate_to_Stemplate = 'T1'
template_skullstrip = 'Manual'
do_surfacewith = 'T1'


regressor_list = ['age']
segmentation_name_list = [lvl1, lvl1LR, lvl2, lvl2LR,
    lvl3, lvl3LR, lvl4, lvl4LR]
segmentation_ID_list = ['lvl1', 'lvl1LR', 'lvl2', 'lvl2LR',
    'lvl3', 'lvl3LR', 'lvl4', 'lvl4LR']

atlas_names_Seg_list = [os.path.basename(path) for path in list_atlases]
Exemples.Stat_study_Zhu_Garin.old.volstat.extractVol(MAIN_PATH, FS_dir, allinfo_study_c, regressor_list, all_ID, all_Session, all_data_path, type_norm, segmentation_name_list,
                                                     segmentation_ID_list, atlas_names_Seg_list, list_atlases, bids_dir)
