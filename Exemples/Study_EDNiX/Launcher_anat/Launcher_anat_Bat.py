import os
import sys
from bids import BIDSLayout
from bids.reports import BIDSReport

opj = os.path.join
opb = os.path.basename

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS
from anatomical import _0_Pipeline_launcher
from anatomical import set_launcher

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

allinfo_study_c = df[(df['suffix'] == 'T2w') & (df['extension'] == '.nii.gz')]

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)

# choose if you want to select or remove ID from you analysis:
list_to_keep   = []
list_to_remove = []

species    = 'Bat'
reference  = 'EDNiX'
addatlas   = ''

# is it a longitudinal study ?
coregistration_longitudinal = False
#do you want to use all the data or only the last one of each subject
which_on  = 'all'                       # "all" or "max"

# create  a study template
creat_study_template  = True

#### Choose to normalize using T1 or T2

type_norm = 'acq-coronal_T2w' # T1 or T2
otheranat = ''                   # '' if none otherwise T1w or T2w
force_myelin_same_space = False
# to get the correct header orientation (valid only with non-human studies) --------------------------------------------

# if you don't know anything about it : leave it empty
# "AHF" stands for "Animal Head First",  "AFF" stands for "Animal Feet First", "humanlike" means no change to be done.
humanPosition     = ['']
orientation       = 'LAI' # "LPI" or ''
animalPosition    = [''] # valid only for species smaller than humans

### masking and skull stripping ----------------------------------------------------------------------------------------

masking_img = 'acq-coronal_T2w' # could be T1w or T2w (if left empty it will be set to "type_norm")

# step 1 : coarse method (use for cropping and acpc setting)
brain_skullstrip_1  = 'bet2'            # bet2_ANTS or MachinL see skullstrip method script for mmore information
# step 2 : precise method
brain_skullstrip_2  = 'NoSkullStrip'            # bet2_ANTS or MachinL
# step 3 : valid only for study or session template :
template_skullstrip = 'Manual'


# valid for resting-state Statistics -------------------------------------------------------------------------------------
do_fMRImasks  = True
fMRImasks     = 'custom' # must be aseg or custom
                       # if custom  please add a ventricle and white matter mask in the template space
                       # named such as Vmask and Wmask


### Coregistration parameters    ---------------------------------------------------------------------------------------
#
transfo_message       = 'do_as_I_said'                           # 'do_it_for_me' or 'do_as_I_said'

Align_img_to_template = 'Ants'
# coregistration between native anatomy and template for acpc setting
transfo_align_dict = {"name": 'align',
                              "interpol": '',               # should be 'nearestNeighbour','linear','hammingWindowedSinc'
                              "type_of_transform": 'Translation',      # should be "SyN","SyNCC","Rigid" or "Affine"
                              "affmetric": 'mattes',              # should be "mattes" (means MI for mutual information),"GC" (global correlation) or "meansquare"
                              "affmetricT": 'mattes'              # should be "mattes" (means MI for mutual information),"GC" (global correlation) or "meansquare"
                              }

transfo_between_dict = {"name": 'between',
                              "interpol": '',
                              "type_of_transform": '',
                              "affmetric": '',
                              "affmetricT": ''
                              }
# coregistration between native anatomy and templates
transfo_coreg_dict = {"name": 'coreg',
                              "interpol": 'hammingWindowedSinc',
                              "type_of_transform": 'SyNBold',
                              "affmetric": 'MI',
                              "affmetricT": 'MI'
                              }
# skull stripping step1
transfo_SS1_dict = {"name": 'SS1',
                              "interpol": '',
                              "type_of_transform": '',
                              "affmetric": '',
                              "affmetricT": ''
                              }
# skull stripping step2
transfo_SS2_dict = {"name": 'SS2',
                              "interpol": '',
                              "type_of_transform": '',
                              "affmetric": '',
                              "affmetricT": ''
                              }
# skull stripping study or session template
transfo_SS3_dict = {"name": 'SS3',
                              "interpol": '',
                              "type_of_transform": '',
                              "affmetric": '',
                              "affmetricT": ''
                              }
# coregistration between study template and reference
transfo_stdyT_dict = {"name": 'stdyT',
                              "interpol": '',
                              "type_of_transform": 'SyNBold',
                              "affmetric": 'mattes',
                              "affmetricT": 'mattes'
                              }

list_transfo = [transfo_align_dict,transfo_between_dict,transfo_coreg_dict,
                transfo_SS1_dict,transfo_SS2_dict,transfo_SS3_dict,transfo_stdyT_dict]


MNIBcorrect_indiv               = ''                      # 'N4' by default. could be set as 'N3'

# set the behaviour of the preprocessing pipeline  ---------------------------------------------------------------------
overwrite_option            = True
check_visualy_final_mask    = False
check_visualy_each_img      = False
do_manual_crop              = False
preftool                    = 'ITK'                 # 'freeview' or 'ITK'



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

Skip_step = [1,2,3,4,5,6,7,8,9,100,'Clean']

########################################################################################################################
#                                       Run the preprocessing steps                                                    #
########################################################################################################################

# get some crucial parameters

atlas_followers =  [[], [], [], []]

(FS_refs, template_dir, reference,balsa_folder, BALSAname, balsa_brainT1,BASE_atlas_folder, BASE_template, BASE_SS,
 BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask,CSF, GM, WM, Aseg_ref,list_atlas, path_label_code,all_ID,
 all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
 fs_tools,reftemplate_path,MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH,bids_dir,allinfo_study_c,species,list_to_keep,
                                                               list_to_remove,reference,type_norm,MNIBcorrect_indiv, masking_img, atlas_followers)

_0_Pipeline_launcher.preprocess_anat(BIDStype, BASE_mask, coregistration_longitudinal, creat_study_template,
    orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, Skip_step,
    check_visualy_each_img, do_fMRImasks, BASE_SS, which_on, all_ID_max, all_data_path_max, all_ID,
    all_Session, all_data_path, template_skullstrip, list_atlas, otheranat, force_myelin_same_space, reference,BALSAname,
    type_norm, all_Session_max, bids_dir, check_visualy_final_mask,
    transfo_message,Align_img_to_template, list_transfo,species, fMRImasks, overwrite_option, MAIN_PATH, fs_tools,reftemplate_path,preftool,FS_refs,path_label_code,
                    BASE_atlas_folder,MNIBcorrect_indiv, animalPosition, humanPosition,balsa_folder,
                    balsa_brainT1,addatlas)

from Tools import Load_EDNiX_requirement
from Plotting import Plot_BIDS_surface_for_QC

sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
    MAIN_PATH, reftemplate_path, bids_dir, 'yes')
# Function 2: Create QC summary for all subjects
Plot_BIDS_surface_for_QC.create_surface_qc_summary(
    sing_wb=sing_wb,
    bids_root=bids_dir,
    output_dir=bids_dir + "/QC/Surface",
    template_scene=bids_dir + "/sub-967HACA/ses-01/anat/native/surfaces/Native_resol/Exemple1.scene",
    scene_ID_name="967HACA",
    scene_name="Exemple1")