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
bids_dir = Load_subject_with_BIDS.linux_path(opj('/srv/projects/easymribrain/scratch/EDNiX/Rat/BIDS_Gd/'))
# which format ?
BIDStype = 2

allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['anat'], suffixes= ['T2w'], extensions=['.nii.gz'])

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# choose if you want to select or remove ID from you analysis:
list_to_keep   = []
list_to_remove = [('300101', '1'), ('300102', '1'), ('300103', '1'), ('300104', '1'), ('300105', '1'), ('300106', '1'), ('300107', '1'), ('300108', '1'), ('300109', '1'), ('300200', '1'), ('300201', '1'), ('300202', '1'), ('300203', '1'), ('300204', '1'), ('300205', '1'), ('300206', '1'), ('300207', '1'), ('300208', '1'), ('300209', '1'), ('300400', '1'), ('300401', '1'), ('300402', '1'), ('300403', '1'), ('300404', '1'), ('300405', '1'), ('300406', '1'), ('300407', '1'), ('300408', '1'), ('300409', '1'), ('300500', '1'), ('300501', '1'), ('300502', '1'), ('300503', '1'), ('300504', '1'), ('300505', '1'), ('300506', '1'), ('300507', '1'), ('300508', '1'), ('300509', '1'), ('300700', '1'), ('300701', '1'), ('300702', '1'), ('300703', '1'), ('300704', '1'), ('300705', '1'), ('300706', '1'), ('300707', '1'), ('300708', '1'), ('300709', '1'), ('300900', '1'), ('300901', '1'), ('300902', '1'), ('300903', '1'), ('300904', '1'), ('300905', '1'), ('300906', '1'), ('300907', '1'), ('300908', '1'), ('300909', '1'), ('301000', '1'), ('301001', '1'), ('301002', '1'), ('301003', '1'), ('301004', '1'), ('301005', '1'), ('301006', '1'), ('301007', '1'), ('301008', '1'), ('301009', '1'), ('301100', '1'), ('301101', '1'), ('301102', '1'), ('301103', '1'), ('301104', '1'), ('301105', '1'), ('301106', '1'), ('301107', '1'), ('301108', '1'), ('301109', '1'), ('301200', '1'), ('301201', '1'), ('301202', '1'), ('301203', '1'), ('301204', '1'), ('301205', '1'), ('301206', '1'), ('301207', '1'), ('301208', '1'), ('301209', '1'), ('301300', '1'), ('301301', '1'), ('301302', '1'), ('301303', '1'), ('301304', '1'), ('301305', '1'), ('301306', '1'), ('301307', '1'), ('301308', '1'), ('301309', '1'), ('301400', '1'), ('301401', '1'), ('301402', '1'), ('301403', '1'), ('301404', '1'), ('301405', '1'), ('301406', '1'), ('301407', '1'), ('301408', '1'), ('301409', '1'), ('301500', '1'), ('301501', '1'), ('301502', '1'), ('301503', '1'), ('301504', '1'), ('301505', '1'), ('301506', '1'), ('301507', '1'), ('301508', '1'), ('301509', '1'), ('301600', '1'), ('301601', '1'), ('301602', '1'), ('301603', '1'), ('301604', '1'), ('301605', '1')]

species    = 'Rat'
# is it a longitudinal study ?
coregistration_longitudinal = False
#do you want to use all the data or only the last one of each subject
which_on  = 'all'                       # "all" or "max"

# create  a study template
creat_study_template  = False

#### Choose to normalize using T1 or T2

type_norm = 'T2w'                      # T1w or T2w
otheranat = ''                         # '' if none otherwise T1w or T2w
# to get the correct header orientation (valid only with non-human studies) --------------------------------------------

# if you don't know anything about it : leave it empty
# "AHF" stands for "Animal Head First",  "AFF" stands for "Animal Feet First", "humanlike" means no change to be done.
humanPosition     = ['humanlike']
orientation       = '' # "LPI" or ''
animalPosition    = ['humanlike'] # valid only for species smaller than humans

### masking and skull stripping ----------------------------------------------------------------------------------------

masking_img = 'T2w' # could be T1w or T2w (if left empty it will be set to "type_norm")

# step 1 : coarse method (use for cropping and acpc setting)
brain_skullstrip_1  = '3dSkullStrip_Rat'            # bet2_ANTS or MachinL see skullstrip method script for mmore information
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
    coreg={'type_of_transform': 'BOLDAffine', 'affmetric': '', 'affmetricT': ''})

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

Skip_step = ['itk_2', 'flat_map', 'Clean']

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
                     fMRImasks, reference='EDNiX', do_fMRImasks=True, atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]], addatlas='',
                     transfo_message='do_as_I_said', force_myelin_same_space=False,
                     check_visualy_final_mask=False, check_visualy_each_img=False, overwrite_option=True, preftool='ITK')

### Surface QC summary creation --------------------------------------------------------------------------------
# Function 1: Load EDNiX requirements
sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
    MAIN_PATH, '', bids_dir, 'yes')
# Function 2: Create QC summary for all subjects
Plot_BIDS_surface_for_QC.create_surface_qc_summary(
    sing_wb=sing_wb,
    bids_root=bids_dir,
    output_dir=bids_dir + "/QC/Surface",
    template_scene=bids_dir + "/sub-301105/ses-1/anat/native/surfaces/Native_resol/Exemple1.scene",
    scene_ID_name="301105",
    scene_name="Exemple1")