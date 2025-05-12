#import
import os
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = r'/srv/projects/easymribrain'
sys.path.append('/home/cgarin/PycharmProjects/EasyMRIbrain_sing/')

import fonctions
import Tools.Load_subject_with_BIDS

species = 'RatWHS'
# Override os.path.join to always return Linux-style paths
bids_dir = Tools.Load_subject_with_BIDS.linux_path(opj('/srv/projects/easymribrain/data/MRI/Rat/BIDS_Gd/'))
FS_dir    = Tools.Load_subject_with_BIDS.linux_path(opj(MAIN_PATH,'FS_Dir_tmp'))
atlas_dir = Tools.Load_subject_with_BIDS.linux_path(opj(MAIN_PATH,'data','Atlas','13_Atlas_project','Atlases_V2', species))
Lut_dir = Tools.Load_subject_with_BIDS.linux_path(opj(MAIN_PATH,'data','Atlas','13_Atlas_project','LUT_files'))

########### Subject loader with BIDS##############
layout= BIDSLayout(bids_dir,  validate=True)
###report
report = BIDSReport(layout)
# Ask get() to return the ids of subjects that have T1w files #return_type='filename
T1 = layout.get(return_type='filename', target='subject', suffix='T1w', extension='nii.gz')
print(T1)
# Ask get() to return the ids of subjects that have T1w files
Bold = layout.get(return_type='filename', target='subject', suffix='bold', extension='nii.gz')
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
list_to_remove = [
        ('300301', 0),
        ('300302', 0),
        ('300303', 0),
        ('300304', 0),
        ('300305', 0),
        ('300306', 0),
        ('300307', 0),
        ('300308', 0),
        ('300309', 0),
        ('300600', 0),
        ('300601', 0),
        ('300602', 0),
        ('300603', 0),
        ('300604', 0),
        ('300605', 0),
        ('300606', 0),
        ('300607', 0),
        ('300608', 0),
        ('300609', 0),
        ('300800', 0),
        ('300801', 0),
        ('300802', 0),
        ('300803', 0),
        ('300804', 0),
        ('300805', 0),
        ('300806', 0),
        ('300807', 0),
        ('300808', 0),
        ('300809', 0)]
all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

coregistration_longitudinal = False
check_visualy_each_img = False
do_manual_crop = False
overwrite_option = True

check_visualy_final_mask = False
deoblique='WARP_without_3drefit'

n_for_ANTS='hammingWindowedSinc'
type_of_transform = 'SyNOnly'
aff_metric_ants = 'CC'

####Choose to normalize using T1 or T2
type_norm = 'T2w' # T1 or T2
otheranat = ''
orientation = 'RIP'
BIDStype = 2

###masking
masking_img = 'T2w'
brain_skullstrip_1 ='XXX' # bet2_ANTS or MachinL
#precise
brain_skullstrip_2 ='XXX' # bet2_ANTS or MachinL
do_fMRImasks = True
fMRImasks = 'XXX' #must be aseg or custom, if custom  please add a ventricle and whitte matter mask in the template space named such as Vmask, Wmask
Align_img_to_template = 'XXX' #3dAllineate or No or @Align_Centers
cost3dAllineate = 'XXX'

#creat_study_template with type_norm img
creat_study_template = False
#do you want to use all the data or only the last one of each subject (for longitud inal co-registration)
which_on = 'all' # all or max
type_of_transform_stdyT = 'XXX'
Atemplate_to_Stemplate = 'XXX'
template_skullstrip = 'XXX'


do_surfacewith = 'XXX'
### file for standardize space
FS_buckner40_TIF = opj(FS_dir,'XXX')
FS_buckner40_GCS = opj(FS_dir,'XXX')

list_atlases = [opj(atlas_dir, 'atlaslvl1.nii.gz'),
opj(atlas_dir, 'atlaslvl2.nii.gz'),
opj(atlas_dir, 'atlaslvl3.nii.gz'),
opj(atlas_dir, 'atlaslvl4.nii.gz'),
opj(atlas_dir, 'atlaslvl1_LR.nii.gz'),
opj(atlas_dir, 'atlaslvl2_LR.nii.gz'),
opj(atlas_dir, 'atlaslvl3_LR.nii.gz'),
opj(atlas_dir, 'atlaslvl4_LR.nii.gz'),
opj(atlas_dir,'XXX'),
opj(atlas_dir, 'XXX')]

BASE_SS     = opj(atlas_dir, 'templateT2.nii.gz') # sting
BASE_mask   = opj(atlas_dir, 'brain_mask.nii.gz') # sting
Aseg_ref    = opj(atlas_dir, 'XXX')
Aseg_refLR  = opj(atlas_dir, 'XXX')

#### for 14 ####
list_atlases_2 = [opj(atlas_dir, 'atlaslvl1.nii.gz'),
opj(atlas_dir, 'atlaslvl2.nii.gz'),
opj(atlas_dir, 'atlaslvl3.nii.gz'),
opj(atlas_dir, 'atlaslvl4.nii.gz')]

FreeSlabel_ctab_list = [opj(Lut_dir,'XXX'),
opj(Lut_dir,'XXX'),
opj(Lut_dir,'XXX'),
opj(Lut_dir,'XXX')]

Lut_file = opj(Lut_dir,'XXX')

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
Skip_step = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,100]

fonctions._0_Pipeline_launcher.preprocess_anat(BIDStype, deoblique, BASE_mask, coregistration_longitudinal, creat_study_template,
                                               orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, n_for_ANTS, aff_metric_ants, Skip_step,
                                               check_visualy_each_img, do_fMRImasks, BASE_SS, which_on, all_ID_max, all_data_path_max, all_ID,
                                               all_Session, all_data_path, template_skullstrip, list_atlases, Aseg_ref, Aseg_refLR, FS_dir,
                                               do_surfacewith, Atemplate_to_Stemplate, FS_buckner40_TIF, FS_buckner40_GCS, Lut_file, otheranat,
                                               type_norm, all_Session_max, bids_dir, check_visualy_final_mask, FreeSlabel_ctab_list,
                                               list_atlases_2, cost3dAllineate, Align_img_to_template, species, type_of_transform,
                                               type_of_transform_stdyT, fMRImasks, overwrite_option, MAIN_PATH)