#import
import os

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

######### to do
### add a verbose feature
### add an option to not redo if it already exists
### improve overwrite option
### QC
### add script to say what have been done per animal

##########################################
########### Subject loader################
##########################################

#https://bids-standard.github.io/pybids/reports/index.html

from Tools import Load_EDNiX_requirement, check_nii, getpath
from anatomical import transfoparams
from anatomical import _loop1
from anatomical import _loop2
from anatomical import _loop3
from anatomical import set_launcher
from anatomical.studytemplate import studytemplate

def preprocess_anat(Skip_step,
                     MAIN_PATH, bids_dir, BIDStype, species,
                     allinfo_study_c, list_to_keep, list_to_remove,
                     type_norm, otheranat, masking_img,
                     orientation, animalPosition, humanPosition,
                     coregistration_longitudinal, creat_study_template, which_on,
                     brain_skullstrip_1, brain_skullstrip_2, template_skullstrip,
                     list_transfo, Align_img_to_template, MNIBcorrect_indiv,
                     fMRImasks, reference='EDNiX', do_fMRImasks=True, atlas_followers=[[], [], [], []], addatlas='',
                     transfo_message='do_as_I_said', force_myelin_same_space=False,
                     check_visualy_final_mask=False, check_visualy_each_img=False, overwrite_option=True, preftool='ITK'):

    """
    Preprocess anatomical MRI data for a BIDS-compliant study, including orientation,
    skull-stripping, normalization, template creation, coregistration, and mask generation.

    The function is designed for both human and non-human datasets (e.g., rodents), supporting
    longitudinal studies, multi-modal anatomical acquisitions, and fMRI-related masks. Preprocessing
    is organized in modular steps that can be selectively skipped.

    Parameters
    ----------
    Skip_step : list
        Steps of the pipeline to skip. Can include integers or descriptive strings (e.g., 'itk_1', 'flat_map').
    MAIN_PATH : str
        Base path for EDNiX installation and resources.
    bids_dir : str
        Path to BIDS-formatted dataset.
    BIDStype : int
        Format identifier for BIDS dataset structure.
    species : str
        Subject species ('Rat', 'Mouse', 'Human', etc.).
    allinfo_study_c : pandas.DataFrame
        DataFrame containing BIDS subject/session information for the study.
    list_to_keep : list
        List of (subject, session) tuples to include.
    list_to_remove : list
        List of (subject, session) tuples to exclude.
    type_norm : str
        Anatomical contrast used for normalization ('T1w' or 'T2w').
    otheranat : str
        Optional second anatomical modality ('T1w' or 'T2w') if available.
    masking_img : str
        Image contrast used for masking/skull-stripping (usually same as type_norm).
    orientation : str
        Target orientation ('LPI' or '').
    animalPosition : list
        Orientation of the animal during scanning (e.g., ['AHF', 'AFF', 'humanlike']).
    humanPosition : list
        Orientation of human subjects ('humanlike' or other).
    coregistration_longitudinal : bool
        Whether to perform longitudinal coregistration across sessions.
    creat_study_template : bool
        Flag to create a study-specific template.
    which_on : str use all the data or only the last one of each subject
        Indicates which data to use for the study-specific template: 'all' for all sessions, 'max' for last session only.
    brain_skullstrip_1 : str
        First-step skull-stripping method for coarse brain extraction.
    brain_skullstrip_2 : str
        Second-step skull-stripping method for fine extraction.
    template_skullstrip : str
        Skull-stripping method for study template or session template.
    list_transfo : dict
        Transformation parameters for rigid and nonlinear registrations.
    Align_img_to_template : str
        Coregistration tool/method (should be "3dAllineate" / "No" / "@Align_Centers" / " Ants").
    MNIBcorrect_indiv : str
        Bias field correction method ('N4', 'N3', or '').
    fMRImasks : str
        Type of fMRI masks to generate ('aseg' or 'custom').
    reference : str, optional
        Reference BIDS dataset or session ('EDNiX' by default).
    do_fMRImasks : bool, optional
        Whether to generate fMRI masks (default: True).
    atlas_followers : list of lists, optional [[], [], [], []]
        Specifies atlas-following regions for multi-level processing.
    addatlas : str, optional
        Additional atlas to include.
    transfo_message : str, optional
        Control behavior for transformations ('do_it_for_me' or 'do_as_I_said').
    force_myelin_same_space : bool, optional
        Force myelin maps to stay in the same space.
    check_visualy_final_mask : bool, optional
        Open QC viewer for the final mask (default: False).
    check_visualy_each_img : bool, optional
        Open QC viewer for each anatomical image (default: False).
    overwrite_option : bool, optional
        Overwrite existing files (default: True).
    preftool : str, optional
        Preferred visualization tool for QC ('ITK' or 'freeview').

    Workflow
    --------
    1. Load BIDS dataset and create DataFrame overview.
    2. Select or exclude subjects and sessions.
    3. Apply orientation correction and coarse cleaning of images.
    4. Optionally create a study template and template mask.
    5. Skull-strip anatomical images using stepwise methods.
    6. Normalize images to template space using ANTs or other tools.
    7. Generate fMRI-related masks and apply atlas transformations.
    8. Optional QC checks on masks and intermediate results.
    9. Optional cleaning and final QC visualization.

    Notes
    -----
    - Steps are modular; users can skip any by specifying Skip_step.
    - Supports both animal and human datasets, including longitudinal designs.
    - Designed for integration with BIDS-formatted datasets for reproducibility.
    """

    print(masking_img)
    (FS_refs, template_dir, reference,balsa_folder, BALSAname, balsa_brainT1,BASE_atlas_folder, BASE_template, BASE_SS,
     BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask,CSF, GM, WM, Aseg_ref,list_atlases, path_label_code,all_ID,
     all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
     fs_tools,reftemplate_path,MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH,bids_dir,allinfo_study_c,species,list_to_keep,
                                                                   list_to_remove,reference,type_norm,MNIBcorrect_indiv, masking_img, atlas_followers)
    ### singularity set up
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _,sing_synstrip,Unetpath =  Load_EDNiX_requirement.load_requirement(MAIN_PATH,reftemplate_path,bids_dir,'yes')


    # Set the environment variable for the current process
    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"

    ######### define automatically other useful parameters  (touch as you own risque) #########

    overwrite    = ''
    if overwrite_option == True:
        overwrite = ' -overwrite'

    type_norm_orig = type_norm
    otheranat_orig = otheranat
    # Usage
    type_norm = check_nii.normalize_anat_type(type_norm)
    if otheranat not in ['', None]:
        otheranat = check_nii.normalize_anat_type(otheranat)
    masking_img = check_nii.normalize_anat_type(masking_img)
    listTimage   = [type_norm]
    IgotbothT1T2 = bool(otheranat.strip())

    listTimage_orig   = [type_norm_orig]
    if IgotbothT1T2 == True:
        listTimage_orig.append(otheranat_orig)

    if IgotbothT1T2 == True:
        listTimage.append(otheranat)

    BASE_SS_coregistr = BASE_SS
    BASE_SS_mask      = BASE_mask
    anat_ref_path = [all_ID[0], all_Session[0], all_data_path[0]]
    # for study template

    study_template_atlas_folder = opj(bids_dir, 'sty_template')
    studyacpc_dir, studyprepro_dir, studytransfo_dir, studyvolume_dir, studylabels_dir, studymasks_dir,_,_,_,_,_,_,_,_ = getpath.stytemplate(
        study_template_atlas_folder, reference, BALSAname)

    Align_img_to_template,list_transfo = transfoparams.set(transfo_message, IgotbothT1T2, creat_study_template, brain_skullstrip_1,
                                     brain_skullstrip_2, template_skullstrip, Align_img_to_template, list_transfo)


    ####################################################################################
    ########################## Start the pipeline !!!!!!!!!!!!!!!!!!!!!!   #############
    ####################################################################################

    # first loop : (contains the first two steps : reorientation and skull stripping)

    for ID, Session, data_path, in zip(all_ID, all_Session, all_data_path):
        ### warning change data_path
        _loop1.run(ID, Session, data_path, data_path,BIDStype, listTimage,otheranat, listTimage_orig, type_norm, orientation,IgotbothT1T2, force_myelin_same_space,
                   Align_img_to_template, list_transfo,masking_img,brain_skullstrip_1, anat_ref_path,
                   BASE_SS_coregistr, BASE_SS_mask, BASE_SS,check_visualy_each_img,check_visualy_final_mask, overwrite,
                   bids_dir,Skip_step,MNIBcorrect_indiv, animalPosition, humanPosition,sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb,sing_synstrip, Unetpath,preftool)

    ####will creat problems???
    if creat_study_template == True:
        # break (contains steps 3 and 4: make study template and skull stripping of that template)

        studytemplate.create(study_template_atlas_folder,Skip_step,which_on, all_ID_max, all_Session_max,
                             all_data_path_max, all_ID,all_Session, all_data_path,type_norm,BASE_SS_coregistr, BASE_SS_mask,
                             list_transfo, template_skullstrip,preftool,check_visualy_final_mask,
                             sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip,reference, Unetpath, BALSAname)


    # second loop (adapt to the longitudinal co-registration and study template

    ref_suffix  = '_space-acpc_desc-SS-step2_'
    mask_suffix = '_desc-step2_'

    stdy_template      = opj(studyprepro_dir, 'studyTemplate_' + type_norm + '.nii.gz')
    stdy_template_mask = opj(studyprepro_dir, 'studyTemplate_mask.nii.gz')

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):

        _loop2.run(ID, Session, data_path, max_ses,type_norm,ref_suffix,coregistration_longitudinal,bids_dir, mask_suffix,
        creat_study_template,stdy_template,stdy_template_mask,BASE_SS,BASE_mask,Skip_step,listTimage,list_transfo,
        brain_skullstrip_2,check_visualy_final_mask,MNIBcorrect_indiv, sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip,preftool,sing_wb,reference,masking_img,Unetpath)


    if creat_study_template==True:
        # second break : adapt the atlases
         BASE_SS, BASE_mask, Aseg_ref = studytemplate.use(study_template_atlas_folder, Skip_step,
                                                                       list_transfo, list_atlases,
                                                                       BASE_SS, BASE_mask,BASE_atlas_folder, species,
                                                                       stdy_template, fMRImasks, reference, BALSAname,
                                                                       path_label_code, sing_afni, sing_wb, which_on,
                                                                       all_data_path_max,all_data_path, listTimage,
                                                                       creat_study_template,
                                                                       coregistration_longitudinal,
                                                                       'native_space-studyTemplate_')


    ###### LOOP AGAIN ######
    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        refnb = 0
        for i, j in enumerate(list_transfo):
            if list_transfo[i]["name"] == 'coreg':
                refnb = i
        _loop3.run(all_ID, all_Session, all_data_path,ID, Session, data_path, max_ses,creat_study_template,
                   study_template_atlas_folder,BASE_SS,BASE_mask,BASE_atlas_folder,coregistration_longitudinal,
                   bids_dir,type_norm,Skip_step,check_visualy_each_img,do_fMRImasks,
                   list_transfo[refnb]["interpol"], fs_tools,sing_afni, sing_itk,sing_fs,sing_wb,
                   fMRImasks,list_atlases,listTimage,species,path_label_code,FS_refs,reference,
        balsa_folder,BALSAname,balsa_brainT1,[0, 0],addatlas, list_transfo)

