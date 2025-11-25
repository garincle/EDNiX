import os
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport

opj = os.path.join
opb = os.path.basename

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS
from anat import _0_Pipeline_launcher
from anat.load_transfo_parameters import build_transfos
from Tools import Load_EDNiX_requirement
from Plotting import Plot_BIDS_surface_for_QC

########################################################################################################################
#                                       Set the pipeline parameters (see manual.....)                                  #
########################################################################################################################

# Where are the data
# Override os.path.join to always return Linux-style paths
bids_dir = Load_subject_with_BIDS.linux_path(opj('/srv/projects/easymribrain/scratch/EDNiX/Bat/BIDS_bat/'))

# which format ?
BIDStype = 1

########### Subject loader with BIDS  ##############
layout = BIDSLayout(bids_dir, validate=False)
report = BIDSReport(layout)
df = layout.to_df()
df.head()

#### Create a pandas sheet for the dataset (I like it, it helps to know what you are about to process)
allinfo_study_c = df[(df['suffix'] == 'T2w') & (df['extension'] == '.nii.gz')]

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# choose if you want to select or remove ID from you analysis:
list_to_keep   = []
list_to_remove = []

species    = 'Bat'

# is it a longitudinal study ?
coregistration_longitudinal = False
#do you want to use all the data or only the last one of each subject
which_on  = 'all'                       # "all" or "max"

# create  a study template
creat_study_template  = True

#### Choose to normalize using T1 or T2

type_norm = 'acq-coronal_T2w'                   # T1w or T2w
otheranat = ''                         # '' if none otherwise T1w or T2w
masking_img = 'acq-coronal_T2w'  # could be T1w or T2w (if left empty it will be set to "type_norm")

# to get the correct header orientation (valid only with non-human studies) --------------------------------------------
# if you don't know anything about it : leave it empty
# "AHF" stands for "Animal Head First",  "AFF" stands for "Animal Feet First", "humanlike" means no change to be done.
humanPosition     = ['']
orientation       = 'LAI' # "LPI" or ''
animalPosition    = [''] # valid only for species smaller than humans

### masking and skull stripping ----------------------------------------------------------------------------------------
# step 1 : coarse method (use for cropping and acpc setting)
brain_skullstrip_1  = '_bet20.5'            # bet2_ANTS or MachinL see skullstrip method script for mmore information
# step 2 : precise method
brain_skullstrip_2  = 'NoSkullStrip'            # bet2_ANTS or MachinL
# step 3 : valid only for study or session template :
template_skullstrip = 'Manual'

# valid for resting-state Statistics -------------------------------------------------------------------------------------
do_fMRImasks  = True
fMRImasks     = 'custom' # must be aseg or custom, if custom  please add a ventricle and white matter mask in the template space, named such as Vmask and Wmask

### Coregistration parameters    ---------------------------------------------------------------------------------------
Align_img_to_template = 'Ants'
MNIBcorrect_indiv               = ''                      # 'N4' by default. could be set as 'N3'

### Coregistration parameters    ---------------------------------------------------------------------------------------
Align_img_to_template = 'Ants'
list_transfo = build_transfos(
    align={'type_of_transform': 'Translation', 'affmetric': 'mattes', 'affmetricT': 'mattes'},
    coreg={'type_of_transform': 'SyNBold', 'affmetric': 'MI', 'affmetricT': 'MI'})

########################################################################################################################
#                                       Select the preprocessing steps                                                 #
#                                                                                                                      #
#   Block1: step 1,2 (orienting, cleaning images)                                                                      #
#   Block1: step 3 (study template)                                                                                    #
#   Block2: step 4 (study template mask QC itksnap)                                                                    #
#   Block2: step 5 (template brain)                                                                                    #
#   Block2: step 6 (registration stdy template)                                                                        #
#   Block2: step 7 (registration anat to template)                                                                     #
#   Block3: step 8,9 (atlases, masks, fmri masks)                                                                      #
#   Block3: step 10,11,12,13,14,15 (surfaces)                                                                          #
#   Block3: step 16 (QC)                                                                                               #
#   Block3: step 100 (Clean)                                                                                           #
#   Block3: step 200 (QC itksnap)                                                                                      #
########################################################################################################################

Skip_step = [1,2,3,4,5,6,10,11,12,13,14,15,16,'itk_2', 'flat_map', 'Clean']

########################################################################################################################
#                                       Run the preprocessing steps                                                    #
########################################################################################################################

_0_Pipeline_launcher.preprocess_anat(Skip_step,
                     MAIN_PATH, bids_dir, BIDStype, species,
                     allinfo_study_c, list_to_keep, list_to_remove,
                     type_norm, otheranat, masking_img,
                     orientation, animalPosition, humanPosition,
                     coregistration_longitudinal, creat_study_template, which_on,
                     brain_skullstrip_1, brain_skullstrip_2, template_skullstrip,
                     list_transfo, Align_img_to_template, MNIBcorrect_indiv,
                     fMRImasks, reference='EDNiX', do_fMRImasks=True, atlas_followers=[['Stuart', 'pcc'], ['ctab', 'txt'], [1, 1], [1, 1]], addatlas='',
                     transfo_message='do_as_I_said', force_myelin_same_space=False,
                     check_visualy_final_mask=False, check_visualy_each_img=False, overwrite_option=True, preftool='ITK')
'''
### Surface QC summary creation --------------------------------------------------------------------------------
# Function 1: Load EDNiX requirements
sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
    MAIN_PATH, '', bids_dir, 'yes')
# Function 2: Create QC summary for all subjects
Plot_BIDS_surface_for_QC.create_surface_qc_summary(
    sing_wb=sing_wb,
    bids_root=bids_dir,
    output_dir=bids_dir + "/QC/Surface",
    template_scene=bids_dir + "/sub-jgrAesMEDISOc22r/ses-2/anat/native/surfaces/Native_resol/Exemple1.scene",
    scene_ID_name="jgrAesMEDISOc22r",
    scene_name="Exemple1")
'''