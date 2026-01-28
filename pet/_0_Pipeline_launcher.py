# import
import os
import subprocess
import glob
import json
import nibabel as nib
import numpy as np
import pandas as pd
from anat import set_launcher
from atlases import atlas4func
import inspect

opj = os.path.join
opb = os.path.basename
opr = os.path.relpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
opi = os.path.isfile

spgo = subprocess.getoutput

from Tools import run_cmd, getpath, diaryfile, Load_EDNiX_requirement, check_nii
from anat.loop3 import backtonative
from fMRI import chooseanat
##### from pet import

def preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endpet, endjson,
                    animalPosition, humanPosition, orientation,
                    T_PET, type_norm, creat_study_template,
                    anat_pet_same_space, coregistration_longitudinal,
                    Method_mask_pet, do_pet_to_func, folderforTemplate_Anat='', IhaveanANAT=True,
                    overwrite_option=True,
                    doWARPonfunc='WARP', registration_fast=False, type_of_transform='BOLDAffine', n_for_ANTS='lanczosWindowedSinc', aff_metric_ants=aff_metric_ants, aff_metric_ants_Transl=aff_metric_ants_Transl, dilate_mask=dilate_mask,
                    list_to_keep=[], list_to_remove=[], atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
                    reference='EDNiX', use_erode_WM_func_masks = True, use_erode_V_func_masks=True,
                    blur=0,
                    selected_atlases_matrix='all', wanted_level_matrix='all',
                    selected_atlases_SBA='default', panda_files_SBA='default'):


    (FS_refs, template_dir, reference, balsa_folder, BALSAname, balsa_brainT1, BASE_atlas_folder, BASE_template, BASE_SS,
    BASE_mask, BASE_Gmask, BASE_WBGmask, BASE_Wmask, BASE_Vmask, CSF, GM, WM, Aseg_ref, list_atlas, path_label_code, all_ID,
    all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
    fs_tools, reftemplate_path, MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH, bids_dir, allinfo_study_c,
                                                                                   species, list_to_keep,
                                                                                   list_to_remove, reference, type_norm,
                                                                                   '', '', atlas_followers)

    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')

    #######for matrix analysis (step 10)
    #### name of the atlases  you want to use for the matrix analysis
    Lut_dir, selected_atlases_matrix, wanted_level, segmentation_name_list, selected_atlases, panda_files = (
        atlas4func.setup(MAIN_PATH, species, reference, selected_atlases_matrix, wanted_level_matrix,
                         selected_atlases_SBA, panda_files_SBA, template_dir))

    if species in ['Human', 'Chimpanzee']:
        config_f = opj(MAIN_PATH, 'Tool_library', 'config', 'b02b0Human.cnf')
    else:
        config_f = opj(MAIN_PATH, 'Tool_library', 'config', 'b02b0Macaque.cnf')

    # Set the environment variable for the current process
    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
    overwrite = ''
    if overwrite_option == True:
        overwrite = ' -overwrite'
    # Usage
    type_norm = check_nii.normalize_anat_type(type_norm)
    T_PET = check_nii.normalize_anat_type(T_PET)

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session)
        run_cmd.printcolor(nl, 'HEADER')

        # Resting data organization folders ###########################################################################
        (dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_raw_matrices, path_pet,
            dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
            dir_prepro_orig_process, dir_prepro_orig_postprocessed,
            dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
            dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
            dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
            dir_prepro_template_process, dir_prepro_template_postprocessed, dir_prepro_acpc_matrices, dir_prepro_orig_matrices) \
            = getpath.pet(data_path, reference)

        for path in [dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_raw_matrices, path_pet,
            dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
            dir_prepro_orig_process, dir_prepro_orig_postprocessed,
            dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
            dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
            dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
            dir_prepro_template_process, dir_prepro_template_postprocessed, dir_prepro_acpc_matrices, dir_prepro_orig_matrices]:

            if not os.path.exists(path):
                os.makedirs(path)

        diary_file = diaryfile.create(opj(path_pet, str(ID) + ' session ' + str(Session)), nl)
        launcher_parameters = diaryfile.create(opj(path_pet, str(ID) + 'launcher_parameters session ' + str(Session)), nl)

        frame = inspect.currentframe()
        args_dict = inspect.getargvalues(frame).locals

        with open(launcher_parameters, "w") as f:
            f.write("Function: preprocess_data\n")
            f.write("Effective parameters at runtime:\n\n")
            for k, v in args_dict.items():
                if k == "kwargs":
                    for kk, vv in v.items():
                        f.write(f"- {kk} = {vv}\n")
                else:
                    f.write(f"- {k} = {v}\n")
        # ----------------------------------------

        diary_WARNING = opj(path_pet, 'MAJOR_WARNING.txt')

        # link with the individual anat template ################################################################
        if IhaveanANAT == False:
            if folderforTemplate_Anat == '':
                nl = 'ERROR: You need to provide either an anat image or a folder with template images to perform the registration'
                raise ValueError(run_cmd.error(nl, diary_file))
            else:
                (anat_subject, brainmask, G_mask, V_mask, W_mask, dir_transfo, FS_dir,
                 dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.create(folderforTemplate_Anat, diary_file)
        else:
            (anat_subject, brainmask, G_mask, WBG_mask, V_mask, W_mask, dir_transfo, FS_dir,
             dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.retrieve(ID, data_path,
                                                                                   Session, anat_pet_same_space,
                                                                                   use_erode_V_func_masks,
                                                                                   use_erode_WM_func_masks,
                                                                                   T_PET, diary_file)
        study_template_atlas_folder = opj(bids_dir, 'sty_template')
        targetsuffix = 'space-acpc_desc-SS'
        info = backtonative.get(ID,data_path,bids_dir,Session,max_ses,targetsuffix,type_norm,BASE_SS, BASE_mask,BASE_atlas_folder,study_template_atlas_folder,
            creat_study_template,coregistration_longitudinal,reference,'Final',species,diary_file)

        # Check the func runs
        list_RS = sorted(glob.glob(opj(path_pet, endpet)))
        RS = [opb(i) for i in list_RS]

        if len(list_RS) == 0:
            nl = ('ERROR : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endpet) and here it is ' +
                        str(opj(path_pet, endpet)) + ' I would check how you define "endpet"')

            raise ValueError(run_cmd.error(nl, diary_file))


            ##maybe load some PET parameters with .json file??

            DIR = os.getcwd()
            nl = 'Working path : ' + DIR
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            # Run the preprocessing ####################################################################################

            if 1 in Skip_step:
                nl = 'skip step ' + str(1)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:


            if 2 in Skip_step:
                nl = 'skip step ' + str(2)
                run_cmd.msg(nl, diary_file, 'OKGREEN')

            else:
