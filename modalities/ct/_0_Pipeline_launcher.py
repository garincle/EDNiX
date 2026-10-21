# import
import os
import subprocess
import glob
import json
import ants
import json
import pandas as pd
from modalities.anat import set_launcher
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

from modalities.ct import reorient, coreg2anat


from modalities.fMRI import chooseanat




def preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endct,Need_check_orient,Need_check_coreg,
                    animalPosition, humanPosition,
                    T_CT, type_norm, creat_study_template,
                    anat_ct_same_space, coregistration_longitudinal,doBALSAspace,
        folderforTemplate_Anat, IhaveanANAT,
                    overwrite_option,
                    list_to_keep, list_to_remove,atlas_followers,
                    reference)


    (FS_refs, template_dir, reference, _, _, _, _, _, _,
    _, _, _k, _, _, _, _, _, _, _, _, all_ID,
    all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
    fs_tools, reftemplate_path, MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH, bids_dir, allinfo_study_c,
                                                                                   species, list_to_keep,
                                                                                   list_to_remove, reference, type_norm,
                                                                                   '', '', atlas_followers, doBALSAspace)

    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, _, _ = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')


    # Set the environment variable for the current process
    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
    overwrite = ''
    if overwrite_option == True:
        overwrite = ' -overwrite'
    # Usage
    type_norm = check_nii.normalize_anat_type(type_norm)


    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session)
        run_cmd.printcolor(nl, 'HEADER')

        # PET organization folders ###########################################################################
        (dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_raw_matrices, path_ct,
            dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
            dir_prepro_orig_process, dir_prepro_orig_postprocessed,
            dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
            dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
            dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
            dir_prepro_template_process, dir_prepro_template_postprocessed, dir_prepro_acpc_matrices, dir_prepro_orig_matrices) \
            = getpath.ct(data_path, reference)


        diary_file = diaryfile.create(opj(path_ct, str(ID) + ' session ' + str(Session)), nl)
        launcher_parameters = diaryfile.create(opj(path_ct, str(ID) + 'launcher_parameters session ' + str(Session)), nl)

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

        diary_WARNING = opj(path_ct, 'MAJOR_WARNING.txt')


        list_img = sorted(glob.glob(opj(dir_prepro_raw, ID + '*_'+ endct+ '.nii.gz')))
        print(opj(dir_prepro_raw, ID +  '*_' + endct + '.nii.gz'))
        print(list_img)

        for paths in [path_ct,dir_prepro_raw_process,dir_prepro_orig,dir_prepro_orig_process,
                     dir_prepro_acpc,dir_prepro_acpc_process]:
            os.makedirs(paths, exist_ok=True)

        if not IhaveanANAT:
            if not folderforTemplate_Anat :
                nl = 'ERROR: You need to provide either an anat image or a folder with template images to perform the registration'
                raise ValueError(run_cmd.error(nl, diary_file))
            else:
                (anat_subject, brainmask, G_mask, WBG_mask, V_mask, W_mask, dir_transfo, FS_dir,
                 dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.create(folderforTemplate_Anat,
                                                                                     False,
                                                                                     True, species,
                                                                                     type_norm, diary_file)
                anat_raw = anat_subject
        else:
            (anat_raw,anat_subject,_,_, _, _,_,dir_transfo) = chooseanat.retrieve4pet(ID, data_path,
                                                                                   Session, anat_ct_same_space,
                                                                                   True,True,
                                                                                   type_norm, diary_file)

            if not anat_subject or not ope(anat_subject):
                nl = (f'INFO: No anat found for {ID} session {Session}, '
                      f'searching across all sessions for this subject...')
                run_cmd.msg(nl, diary_file, 'WARNING')

                anat_found = False
                for ses in sorted(set(all_Session)):
                    candidate_path = opj(bids_dir, f'sub-{ID}', f'ses-{ses}')
                    if not ope(candidate_path):
                        continue
                    try:
                        (anat_raw,anat_subject,brainmask,G_mask, WBG_mask, W_mask,V_mask,dir_transfo) \
                            = chooseanat.retrieve4pet(ID, candidate_path,ses,anat_ct_same_space,True,True,
                                                      type_norm, diary_file)
                        if anat_subject  and ope(anat_subject):
                            nl = f'INFO: Found anat for {ID} in session {ses}, using it as fallback.'
                            run_cmd.msg(nl, diary_file, 'OKGREEN')
                            anat_session_used = ses  # ← on mémorise la session anat réelle
                            anat_data_path_used = candidate_path
                            anat_found = True
                            break
                    except Exception:
                        continue

                if not anat_found:
                    nl = f'ERROR: No anat found for {ID} in any session. Cannot proceed.'
                    raise FileNotFoundError(run_cmd.error(nl, diary_file))


        if 1 in Skip_step:
            run_cmd.msg('INFO: skip step ' + str(1), diary_file, 'OKGREEN')
        else:
            for img in list_img:
                jsonfile = img.replace('.nii.gz', '.json')
                if opi(jsonfile):
                    f = open(jsonfile)
                    info_ct = json.load(f)
                    if humanPosition == 'unknown':
                        if "PatientPosition" in info_ct:
                            humanPosition = info_ct["PatientPosition"]

                reorient.doit(dir_prepro_orig_process, img, humanPosition, animalPosition, anat_raw,
                              Need_check_orient,diary_file, sing_afni, sing_fs)


        if 2 in Skip_step:
            nl = 'skip step ' + str(1)
            run_cmd.msg(nl, diary_file, 'OKGREEN')
        else:
            if IhaveanANAT:
                for img in list_img:
                    coreg2anat.doit(img, anat_raw, T_CT, anat_subject,
                     dir_transfo, dir_prepro_orig_process, dir_prepro_acpc,
                     dir_prepro_acpc_process,
                     dir_prepro_orig_matrices, anat_ct_same_space, False,
                     Need_check_coreg, diary_file, sing_afni, sing_fs)






