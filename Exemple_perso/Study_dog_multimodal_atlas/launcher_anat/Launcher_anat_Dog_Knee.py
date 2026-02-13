import os
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport

opj = os.path.join
opb = os.path.basename

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS, load_bids
from anat import _0_Pipeline_launcher
from anat.load_transfo_parameters import build_transfos
from Tools import Load_EDNiX_requirement
from Plotting import Plot_BIDS_surface_for_QC
########################################################################################################################
#                                       Set the pipeline parameters (see manual.....)                                  #
########################################################################################################################

# Where are the data

# Override os.path.join to always return Linux-style paths
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Dog/BIDS_knee/'))
# which format ?
BIDStype = 2

allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['anat'], suffixes= ['T1w'], extensions=['.nii.gz'])

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# choose if you want to select or remove ID from you analysis:
list_to_keep   = []
list_to_remove = []

species    = 'Dog'
# is it a longitudinal study ?
coregistration_longitudinal = False
#do you want to use all the data or only the last one of each subject
which_on  = 'all'                       # "all" or "max"

# create  a study template
creat_study_template  = False

#### Choose to normalize using T1 or T2

type_norm = 'T1w'                      # T1w or T2w
otheranat = ''                         # '' if none otherwise T1w or T2w
masking_img = 'T1w' # could be T1w or T2w (if left empty it will be set to "type_norm")

# to get the correct header orientation (valid only with non-human studies) --------------------------------------------

# if you don't know anything about it : leave it empty
# "AHF" stands for "Animal Head First",  "AFF" stands for "Animal Feet First", "humanlike" means no change to be done.
humanPosition     = ['humanlike']
orientation       = '' # "LPI" or ''
animalPosition    = ['humanlike'] # valid only for species smaller than humans

### masking and skull stripping ----------------------------------------------------------------------------------------


# step 1 : coarse method (use for cropping and acpc setting)
brain_skullstrip_1  = 'muSkullStrip_cross_species'            # bet2_ANTS or MachinL see skullstrip method script for mmore information
# step 2 : precise method
brain_skullstrip_2  = 'NoSkullStrip'            # bet2_ANTS or MachinL
# step 3 : valid only for study or session template :
template_skullstrip = 'NoSkullStrip'

# valid for resting-state Statistics -------------------------------------------------------------------------------------
fMRImasks     = 'aseg' # must be aseg or custom
                       # if custom  please add a ventricle and white matter mask in the template space
                       # named such as Vmask and Wmask

### Coregistration parameters    ---------------------------------------------------------------------------------------
#
Align_img_to_template = 'Ants'

list_transfo = build_transfos(
    align={'type_of_transform': 'Rigid', 'affmetric': 'MI', 'affmetricT': 'MI'},
    coreg={'type_of_transform': 'SyNCC', 'affmetric': 'MI', 'affmetricT': 'MI'})

MNIBcorrect_indiv               = ''                      # 'N4' by default. could be set as 'N3'

########################################################################################################################
#                                                                                                                      #
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
#                                                                                                                      #
########################################################################################################################


########################################################################################################################
#                                       Run the preprocessing steps                                                    #
########################################################################################################################
Skip_step = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,'itk_2','itk_1','flat_map', 'Clean']
_0_Pipeline_launcher.preprocess_anat(Skip_step,
                     MAIN_PATH, bids_dir, BIDStype, species,
                     allinfo_study_c, list_to_keep, list_to_remove,
                     type_norm, otheranat, masking_img,
                     orientation, animalPosition, humanPosition,
                     coregistration_longitudinal, creat_study_template, which_on,
                     brain_skullstrip_1, brain_skullstrip_2, template_skullstrip,
                     list_transfo, Align_img_to_template, MNIBcorrect_indiv,
                     fMRImasks, reference='EDNiX', do_fMRImasks=True, atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]], addatlas='',
                     transfo_message='do_as_I_said', force_myelin_same_space=False,
                     check_visualy_final_mask=True, check_visualy_each_img=False, overwrite_option=True, preftool='ITK')

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
    template_scene=bids_dir + "/sub-01/ses-1/anat/native/surfaces/Native_resol/Exemple1.scene",
    scene_ID_name="01",
    scene_name="Exemple1")
'''
