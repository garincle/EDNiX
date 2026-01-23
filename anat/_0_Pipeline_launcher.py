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
from anat import transfoparams
from anat import _loop1
from anat import _loop2
from anat import _loop3
from anat import set_launcher
from anat.studytemplate import studytemplate

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
    Preprocess anatomical MRI data for BIDS-compliant studies, supporting human and non-human datasets.

    Performs modular preprocessing including orientation correction, skull-stripping, normalization,
    study template creation, longitudinal coregistration, and fMRI mask generation.

    :param Skip_step: Steps of the pipeline to skip (integers or descriptive strings)
    :type Skip_step: list
    :param MAIN_PATH: Base path for EDNiX installation and resources
    :type MAIN_PATH: str
    :param bids_dir: Path to BIDS-formatted dataset
    :type bids_dir: str
    :param BIDStype: Format identifier for BIDS dataset structure
    :type BIDStype: int
    :param species: Subject species ('Human', 'Rat', 'Mouse', etc. as in Atlases_library/atlas)
    :type species: str
    :param allinfo_study_c: DataFrame with subject/session metadata
    :type allinfo_study_c: pandas.DataFrame
    :param list_to_keep: List of (subject, session) tuples to include
    :type list_to_keep: list
    :param list_to_remove: List of (subject, session) tuples to exclude
    :type list_to_remove: list
    :param type_norm: Anatomical contrast for normalization ('T1w' or 'T2w')
    :type type_norm: str
    :param otheranat: Optional secondary anatomical contrast
    :type otheranat: str
    :param masking_img: Image contrast for masking/skull-stripping
    :type masking_img: str
    :param orientation: Target image orientation ('LPI', etc.)
    :type orientation: str
    :param animalPosition: Animal positioning during scanning
    :type animalPosition: list
    :param humanPosition: Human positioning
    :type humanPosition: list
    :param coregistration_longitudinal: Perform longitudinal coregistration
    :type coregistration_longitudinal: bool
    :param creat_study_template: Flag to create study-specific template
    :type creat_study_template: bool
    :param which_on: Template session selection: 'all' or 'max'
    :type which_on: str
    :param brain_skullstrip_1: Coarse skull-stripping method
    :type brain_skullstrip_1: str
    :param brain_skullstrip_2: Fine skull-stripping method
    :type brain_skullstrip_2: str
    :param template_skullstrip: Skull-stripping method for study/session template
    :type template_skullstrip: str
    :param list_transfo: Parameters for rigid and nonlinear transformations
    :type list_transfo: dict
    :param Align_img_to_template: Coregistration tool/method
    :type Align_img_to_template: str
    :param MNIBcorrect_indiv: Bias field correction method ('N4', 'N3', or '')
    :type MNIBcorrect_indiv: str
    :param fMRImasks: Type of fMRI masks to generate ('aseg' or 'custom')
    :type fMRImasks: str
    :param reference: Reference dataset or session
    :type reference: str
    :param do_fMRImasks: Generate fMRI masks
    :type do_fMRImasks: bool
    :param atlas_followers: Multi-level atlas-following regions
    :type atlas_followers: list of lists
    :param addatlas: Additional atlas to include
    :type addatlas: str
    :param transfo_message: Transformation control ('do_it_for_me' or 'do_as_I_said')
    :type transfo_message: str
    :param force_myelin_same_space: Maintain myelin maps in same space
    :type force_myelin_same_space: bool
    :param check_visualy_final_mask: Open QC viewer for final mask
    :type check_visualy_final_mask: bool
    :param check_visualy_each_img: Open QC viewer for each anatomical image
    :type check_visualy_each_img: bool
    :param overwrite_option: Overwrite existing files
    :type overwrite_option: bool
    :param preftool: Preferred QC visualization tool ('ITK' or 'freeview')
    :type preftool: str

    :return: None
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
         stdy_template_SS = opj(studyprepro_dir, 'studyTemplate_SS_' + type_norm + '.nii.gz')
         _, _, _ = studytemplate.use(study_template_atlas_folder, Skip_step,
                                                                       list_transfo, list_atlases,
                                                                       BASE_SS, BASE_mask,BASE_atlas_folder, species,
                                                                       stdy_template_SS, fMRImasks, reference, BALSAname,
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

