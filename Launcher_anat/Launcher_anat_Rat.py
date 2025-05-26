import os
from bids import BIDSLayout
from bids.reports import BIDSReport
opn = os.path.normpath
opj = os.path.join
MAIN_PATH = opj('/srv/projects/easymribrain/code/EDNiX/')
import Tools.Load_subject_with_BIDS
import Tools.Read_atlas
import anatomical._0_Pipeline_launcher

species = 'RatWHS'
# Override os.path.join to always return Linux-style paths
bids_dir = Tools.Load_subject_with_BIDS.linux_path(opj(r"C:\Users\cgarin\Desktop\BIDS_k9"))
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
allinfo_study_c = df[(df['suffix'] == 'T2w') & (df['extension'] == '.nii.gz')]

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


coregistration_longitudinal = False
check_visualy_each_img = False
do_manual_crop = False
overwrite_option = True

check_visualy_final_mask = False
deoblique='WARP_without_3drefit'

n_for_ANTS='hammingWindowedSinc'
type_of_transform = 'SyNCC'
aff_metric_ants = 'MI'

####Choose to normalize using T1 or T2
type_norm = 'T2w' # T1 or T2
otheranat = ''
orientation = 'RIP'
BIDStype = 2

###masking
masking_img = 'T2w'
brain_skullstrip_1 ='3dSkullStrip_Rat' # bet2_ANTS or MachinL
#precise
brain_skullstrip_2 ='Custum_QWARPT2' # bet2_ANTS or MachinL
do_fMRImasks = True
fMRImasks = 'aseg' #must be aseg or custom, if custom  please add a ventricle and whitte matter mask in the template space named such as Vmask, Wmask
Align_img_to_template = '3dAllineate' #3dAllineate or No or @Align_Centers
cost3dAllineate = 'lpa'

#creat_study_template with type_norm img
creat_study_template = False
#do you want to use all the data or only the last one of each subject (for longitud inal co-registration)
which_on = 'all' # all or max
type_of_transform_stdyT = 'SyN'
Atemplate_to_Stemplate = ''
template_skullstrip = ''


do_surfacewith = 'T2w'
### file for standardize space
FS_buckner40_TIF = opj(FS_dir,'MacaqueYerkes19')
FS_buckner40_GCS = opj(FS_dir,'MacaqueYerkes19')

list_atlases = [opj(atlas_dir, 'atlaslvl1.nii.gz'),
opj(atlas_dir, 'atlaslvl2.nii.gz'),
opj(atlas_dir, 'atlaslvl3.nii.gz'),
opj(atlas_dir, 'atlaslvl4.nii.gz'),
opj(atlas_dir, 'atlaslvl1_LR.nii.gz'),
opj(atlas_dir, 'atlaslvl2_LR.nii.gz'),
opj(atlas_dir, 'atlaslvl3_LR.nii.gz'),
opj(atlas_dir, 'atlaslvl4_LR.nii.gz'),
opj(atlas_dir,'Gmask.nii.gz'),
opj(atlas_dir, 'Wmask.nii.gz')]

BASE_SS     = opj(atlas_dir,'templateT2.nii.gz') # sting
BASE_mask   = opj(atlas_dir, 'brain_mask.nii.gz') # sting
Aseg_ref    = opj(atlas_dir, 'atlas_forSEG_final.nii.gz')
Aseg_refLR  = opj(atlas_dir, 'atlas_forSEG_final_LR.nii.gz')

#### for 14 ####
list_atlases_2 = [opj(atlas_dir, 'atlaslvl1.nii.gz'),
opj(atlas_dir, 'atlaslvl2.nii.gz'),
opj(atlas_dir, 'atlaslvl3.nii.gz'),
opj(atlas_dir, 'atlaslvl4.nii.gz')]

FreeSlabel_ctab_list = [opj(Lut_dir,'Multispecies_LUT_Dual.txt'),
opj(Lut_dir,'Multispecies_LUT_Dual.txt'),
opj(Lut_dir,'Multispecies_LUT_Dual.txt'),
opj(Lut_dir,'Multispecies_LUT_Dual.txt')]

Lut_file = opj(Lut_dir,'Multispecies_LUT_Dual.txt')

### Block1: step 1,2 (orienting, cleaning images)
### Block1: step 3 (study template)
### Block2: step 4 (study template mask QC itksnap)
### Block2: step 5 (template brain)
### Block2: step 6 (registration stdy template)
### Block2: step 7 (registration anat to template)
### Block3: step 8,9 (altases, masks, fmri masks)
### Block3: step 10,11,12,13,14,15 (surfaces)
### Block3: step 16 (QC)
### Block3: step 100 (Clean)
### Block3: step 200 (QC itksnap)
Skip_step = [100,200]

anatomical._0_Pipeline_launcher.preprocess_anat(BIDStype, deoblique, BASE_mask, coregistration_longitudinal, creat_study_template,
    orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, n_for_ANTS, aff_metric_ants, Skip_step,
    check_visualy_each_img, do_fMRImasks, BASE_SS, which_on, all_ID_max, all_data_path_max, all_ID,
    all_Session, all_data_path, template_skullstrip, list_atlases, Aseg_ref, Aseg_refLR, FS_dir,
    do_surfacewith, Atemplate_to_Stemplate, FS_buckner40_TIF,FS_buckner40_GCS, Lut_file, otheranat,
    type_norm, all_Session_max, bids_dir, check_visualy_final_mask, FreeSlabel_ctab_list,
    list_atlases_2, cost3dAllineate, Align_img_to_template, species, type_of_transform,
    type_of_transform_stdyT, fMRImasks, overwrite_option, MAIN_PATH)