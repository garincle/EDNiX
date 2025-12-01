import os
import sys
import Tools.Read_atlas
import fMRI._0_Pipeline_launcher
from Tools import Load_subject_with_BIDS, load_bids
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

species    = 'Rat'
# Override os.path.join to always return Linux-style paths
bids_dir = Load_subject_with_BIDS.linux_path(opj('/srv/projects/easymribrain/scratch/EDNiX/Rat/BIDS_Gd/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['anat'], suffixes= ['T2w'], extensions=['.nii.gz'])

### select the subject, session to process
Tools.Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
# choose if you want to select or remove ID from you analysis
list_to_keep = []
list_to_remove = [('300101', '1'), ('300102', '1'), ('300103', '1'), ('300104', '1'), ('300105', '1'), ('300106', '1'), ('300107', '1'), ('300108', '1'), ('300109', '1'),
                  ('300200', '1'), ('300201', '1'), ('300202', '1'), ('300203', '1'), ('300204', '1'), ('300205', '1'), ('300206', '1'), ('300207', '1'),
                  ('300208', '1'), ('300209', '1'), ('300400', '1'), ('300401', '1'), ('300402', '1'), ('300403', '1'), ('300404', '1'), ('300405', '1'),
                  ('300406', '1'), ('300407', '1'), ('300408', '1'), ('300409', '1'), ('300500', '1'), ('300501', '1'), ('300502', '1'), ('300503', '1'),
                  ('300504', '1'), ('300505', '1'), ('300506', '1'), ('300507', '1'), ('300508', '1'), ('300509', '1'), ('300700', '1'), ('300701', '1'),
                  ('300702', '1'), ('300703', '1'), ('300704', '1'), ('300705', '1'), ('300706', '1'), ('300707', '1'), ('300708', '1'), ('300709', '1'),
                  ('300900', '1'), ('300901', '1'), ('300902', '1'), ('300903', '1'), ('300904', '1'), ('300905', '1'), ('300906', '1'), ('300907', '1'),
                  ('300908', '1'), ('300909', '1'), ('301000', '1'), ('301001', '1'), ('301002', '1'), ('301003', '1'), ('301004', '1'), ('301005', '1'),
                  ('301006', '1'), ('301007', '1'), ('301008', '1'), ('301009', '1'), ('301100', '1'), ('301101', '1'), ('301102', '1'), ('301103', '1'),
                  ('301104', '1'), ('301105', '1'), ('301106', '1'), ('301107', '1'), ('301108', '1'), ('301109', '1'), ('301200', '1'), ('301201', '1'),
                  ('301202', '1'), ('301203', '1'), ('301204', '1'), ('301205', '1'), ('301206', '1'), ('301207', '1'), ('301208', '1'), ('301209', '1'),
                  ('301300', '1'), ('301301', '1'), ('301302', '1'), ('301303', '1'), ('301304', '1'), ('301305', '1'), ('301306', '1'), ('301307', '1'),
                  ('301308', '1'), ('301309', '1'), ('301400', '1'), ('301401', '1'), ('301402', '1'), ('301403', '1'), ('301404', '1'), ('301405', '1'),
                  ('301406', '1'), ('301407', '1'), ('301408', '1'), ('301409', '1'), ('301500', '1'), ('301501', '1'), ('301502', '1'), ('301503', '1'),
                  ('301504', '1'), ('301505', '1'), ('301506', '1'), ('301507', '1'), ('301508', '1'), ('301509', '1'), ('301600', '1'), ('301601', '1'),
                  ('301602', '1'), ('301603', '1'), ('301604', '1'), ('301605', '1'), ('301606', '1'), ('301607', '1'), ('301608', '1'), ('301609', '1')]

all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Tools.Load_subject_with_BIDS.load_data_bids(allinfo_study_c, bids_dir, list_to_keep, list_to_remove)
#### fMRI pre-treatment
T1_eq = 5 # int
REF_int = 0 # int
ntimepoint_treshold = 100
endfmri = '*_task-rest_bold.nii.gz' # string
endjson = '*_task-rest_bold.json' # string
endmap = '*_map.nii.gz' # string
humanPosition     = ['humanlike']
orientation       = '' # "LPI" or ''
animalPosition    = ['humanlike'] # valid only for species smaller than humans
## prior anat processing
coregistration_longitudinal = False #True or False
type_norm = 'T2w' # T1 or T2
### co-registration func to anat to template to with T1 ? T2? use the correct  suffix as in the BIDS
TfMRI = 'T2w' # string
### if you don't have any anat image you will need to put several image in the folderforTemplate_Anat (refer to the doc)
Method_mask_func = 'nilearn' # string 3dAllineate or nilearn or creat a manual mask in the funcsapce folder name "manual_mask.nii.gz"
#### ANTs function of the co-registration HammingWindowedSinc is advised
IhaveanANAT = True # True or False
anat_func_same_space = True # True or False
type_of_transform = 'BOLDAffine'
aff_metric_ants_Transl = 'mattes' # string
aff_metric_ants = 'mattes'
do_anat_to_func = False # True or False
Slice_timing_info = 'Auto'
##### if you don't have an anat then template will be the same as anat...
#creat_study_template was created with the anat type_norm img, and you want to use it as standart space
creat_study_template = False # True or Fals
blur = 0.5 # float
#Dilate the functional brain mask by n layers
dilate_mask = 2 # int
smoothSBA = 0.5
# for the seed base analysis, you need to provide the names and the labels of the regions you want to use as "seeds"
selected_atlases = [['EDNIxCSC', 3]]  # Using NEW VERSION format (single atlas)
############ Right in a list format the steps that you want to skip
doWARPonfunc = False
resting_or_task = 'resting'  # 'resting' or 'task'

Skip_step = ['itk_1', 'itk_2', 'Clean']
fMRI._0_Pipeline_launcher.preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endfmri, endjson, endmap, resting_or_task,
                    animalPosition, humanPosition, orientation,
                    Slice_timing_info,
                    TfMRI, type_norm, creat_study_template,
                    anat_func_same_space, coregistration_longitudinal,
                    Method_mask_func, do_anat_to_func, folderforTemplate_Anat='', IhaveanANAT=True,
                    ntimepoint_treshold=100, REF_int=0, T1_eq=5, correction_direction='Auto', overwrite_option=True,
                    DwellT='Auto', SED='Auto', TR='Auto', TRT='Auto',
                    nb_ICA_run=20, ICA_cleaning='Skip',
                    costAllin='lpa',
                    doWARPonfunc=doWARPonfunc, registration_fast=False, type_of_transform=type_of_transform, n_for_ANTS='lanczosWindowedSinc', aff_metric_ants=aff_metric_ants, aff_metric_ants_Transl=aff_metric_ants_Transl, dilate_mask=dilate_mask,
                    list_to_keep=list_to_keep, list_to_remove=list_to_remove, atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
                    reference='EDNiX', post_treatment_method='Grandjean',
                    band='0.01 0.1', blur=0, do_not_correct_signal = False, extract_exterior_CSF = False, extract_WM=True, extract_Vc = False, extract_GS = False,
                    use_erode_WM_func_masks = True, use_erode_V_func_masks=True, normalize='Skip',
                    selected_atlases_matrix='all', wanted_level_matrix='all',
                    selected_atlases_SBA='default', panda_files_SBA='default',
                    SBAspace=['atlas'], erod_seed=True, smoothSBA=smoothSBA,
                    specific_roi_tresh=0.2, delta_thresh=0.1,
                    oversample_map=False, use_cortical_mask_func=False, n_cut=10, threshold_val=10)