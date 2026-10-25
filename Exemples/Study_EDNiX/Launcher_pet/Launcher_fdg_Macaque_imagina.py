import os
import sys
opn = os.path.normpath
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS, load_bids
from modalities.pet import _0_Pipeline_launcher

######################################################################################################################
# [18F]FDG PET of the macaques of BIDS_Tremblay, with the T1w of the anatomical pipeline
# (Launcher_anat_Macaque_imagina.py must have run first). Scenario 1v2a of Illustration/pet4.svg: the PET is
# registered rigidly to the anatomy, then follows the anatomy -> EDNiX template chain of the anatomical pipeline.
######################################################################################################################
species = 'Macaque'
bids_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Macaque/BIDS_Tremblay/'))
allinfo_study_c = load_bids.Load_BIDS_to_pandas(bids_dir, modalities=['pet'], suffixes=['pet'], extensions=['.nii.gz'])

### select the subject, session to process
Load_subject_with_BIDS.print_included_tuples(allinfo_study_c)
list_to_keep = []
list_to_remove = []

### PET images
endpet = '*_trc-FDG_*pet.nii.gz'  # raw PET images in sub-*/ses-*/pet
tracer = 'FDG'
isotope = 'unknown'               # 'unknown': TracerRadionuclide of the sidecar, else F-18 for FDG

### orientation of the PET header (as for the anatomy)
humanPosition = ['']              # '': PatientPosition of the sidecar
animalPosition = ['']             # '' = humanlike
orientation = 'RIP'               # code forced on the header ('' to use the positions): CHECK it on one subject
Need_check_orient = False         # True: freeview after the reorientation

### anatomy: DECLARATIONS of what the anatomical pipeline did (same values as in the anatomical launcher); they only
### tell where the anatomy -> template transforms are, they build nothing for the PET
type_norm = 'T1w'
T_PET = 'T1w'                     # anatomical image the PET is registered to
creat_study_template = False      # the ANATOMICAL pipeline used an anatomical study template
coregistration_longitudinal = False  # the ANATOMICAL pipeline registered the sessions to a reference session
# PET sessions are not registered to each other: each one goes to the anatomy of its session (1v2a here) and they
# meet in the AC-PC / template spaces through the anatomical chain. A PET study template (creat_study_template_pet,
# not set here: False) only serves the crop and the scenario without anatomy (3v2); it is not registered to the anatomy.
IhaveanANAT = True
anat_pet_same_space = False       # True only for a hybrid PET-MR acquisition (scenario 1v1)

### motion correction
nb_void = 7                       # first frames not registered (low counts)
crop_for_MoCo = 'auto'            # head box: without PET study template, drawn in freeview on the session mean at the
                                  # first of steps 4/5/6 that runs, saved in pet/orig/preprocessing/forMC/
                                  # <session>_desc-crop_pet.nii.gz and reused; QC in pet/QC/*_desc-cropbox_pet.png;
                                  # 'none' to run without display
redo_crop = False                 # True: draw every crop again (session crops and study-template crop)

### parametric images
SUV = True                        # needs the injected radioactivity (sidecar) and the body weight (participants.tsv)
bratio = True                     # ratio to the mean over the brain mask
SUV_window = None                 # [start, end] min of the frames averaged for a dynamic acquisition; None: all
Do_cut = False
duration = 10

### reference-region kinetic models (dynamic receptor tracers only)
srtm = False
reference_signal = 'unknown'
DV_time = 0

Skip_step = [1,2,3,4,5,'pvc', 9]
_0_Pipeline_launcher.preprocess_pet(
    Skip_step, MAIN_PATH, bids_dir, species, allinfo_study_c, endpet, tracer,
    animalPosition, humanPosition, orientation,
    T_PET, type_norm, creat_study_template, coregistration_longitudinal,
    IhaveanANAT=IhaveanANAT, anat_pet_same_space=anat_pet_same_space,
    isotope=isotope, nb_void=nb_void, crop_for_MoCo=crop_for_MoCo, redo_crop=redo_crop,
    SUV=SUV, bratio=bratio, SUV_window=SUV_window, Do_cut=Do_cut, duration=duration,
    srtm=srtm, reference_signal=reference_signal, DV_time=DV_time,
    export_resolution='pet', selected_atlases='all',
    coreg_metric='mattes', type_of_transform_template='SyNRA',
    # PET -> T1w: rigid (same subject; 'Similarity' adds a scale that can shrink the PET brain into the MRI brain),
    # multi-start initialisation, PET restricted to the brain of the anatomy in a second pass (extracerebral FDG
    # uptake of the macaque: temporal muscles, Harderian glands, eye muscles)
    coreg_transform='Rigid', coreg_init='auto', pet_coreg_mask='anat', pet_mask_dilate_mm=3.0,
    coreg_target='anat',
    atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
    reference='EDNiX', list_to_keep=list_to_keep, list_to_remove=list_to_remove,
    Need_check_orient=Need_check_orient, Need_check_MC=False, Need_check_coreg=False)
