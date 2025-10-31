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

from anatomical.studytemplate import studytemplate

def preprocess_anat(BIDStype, BASE_mask, coregistration_longitudinal, creat_study_template,
    orientation, masking_img, brain_skullstrip_1, brain_skullstrip_2, Skip_step,
    check_visualy_each_img, do_fMRImasks, BASE_SS, which_on, all_ID_max, all_data_path_max, all_ID,
    all_Session, all_data_path, template_skullstrip, list_atlases,otheranat,force_myelin_same_space, reference,BALSAname,
    type_norm, all_Session_max, bids_dir, check_visualy_final_mask,
    transfo_message,Align_img_to_template, list_transfo,species, fMRImasks, overwrite_option, MAIN_PATH, fs_tools,reftemplate_path,preftool,FS_refs,path_label_code,
                    BASE_atlas_folder,MNIBcorrect_indiv, animalPosition, humanPosition,balsa_folder,
                    balsa_brainT1,addatlas):

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
                             all_data_path_max, all_ID,all_Session, all_data_path,type_norm,BASE_SS_coregistr,BASE_SS_mask,
                             list_transfo,template_skullstrip,preftool,check_visualy_final_mask,
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
                                                                       stdy_template, type_norm, reference, BALSAname,
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

