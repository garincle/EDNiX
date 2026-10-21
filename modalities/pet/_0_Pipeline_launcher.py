# import
import os
import subprocess
import glob
import json
import nibabel as nib
import numpy as np
import pandas as pd

from Exemples.Study_EDNiX.Launcher_pet.Launcher_fdg_CERMEP import baseline
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

from modalities.pet import loop1, studyT, loop2, baseline_def

from atlases import getlabels
from modalities.fMRI import chooseanat

known_tracor = {"FDG":"F-18","RACLO":"C-11","BUPRE":"C-11","2FNQ1P":"F-18","H2O":"O-15","PE2I":"C-11"}


def preprocess_data(
                    Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endpet, endjson,dcmfolder,endfiledcm,
                    animalPosition, humanPosition, orientation,tracor,
                    T_PET, type_norm, creat_study_template,nbvoid,DV_time, Do_cut, duration, suv, bratio, strm,
                    anat_pet_same_space, coregistration_longitudinal,IhaveaCT,ct_pet_same_space,
                    Method_mask_pet, do_pet_to_func, studyT_transform, doBALSAspace,selected_atlases, reference_signal,isotope,folderforTemplate_Anat, IhaveanANAT,
                    overwrite_option,Need_check_MC,Need_check_coreg,do_manual,
                    doWARPonfunc, registration_fast, type_of_transform,n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl,
                    list_to_keep, list_to_remove,atlas_followers,
                    reference, use_erode_WM_func_masks ,
                    selected_atlases_matrix, wanted_level_matrix):


    (FS_refs, template_dir, reference, balsa_folder, BALSAname, balsa_brainT1, BASE_atlas_folder, BASE_template, BASE_SS,
    BASE_mask, BASE_Gmask, BASE_WBGmask, BASE_Wmask, BASE_Vmask, CSF, GM, WM, Aseg_ref, list_atlas, path_label_code, all_ID,
    all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
    fs_tools, reftemplate_path, MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH, bids_dir, allinfo_study_c,
                                                                                   species, list_to_keep,
                                                                                   list_to_remove, reference, type_norm,
                                                                                   '', '', atlas_followers, doBALSAspace)

    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, sing_turku, _, Unetpath = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')

    template_atlas_file = opj(BASE_atlas_folder, reference + '_seg-' + selected_atlases["name"] + '_dseg.nii.gz')

    if strm ==1:
        if reference_signal['type'] == 'masks':
            reference_signal = getlabels.gettype(template_atlas_file, selected_atlases['level'], '')

        if reference_signal['type'] == 'atlas' and reference_signal['label']=='unknown':

            if not 'EDNIx' in selected_atlases["name"]:
                if not opi(opj(path_label_code,selected_atlases["name"] + '_gray.csv')):

                    file = opj(path_label_code, selected_atlases["name"] + '_StatsLUT.txt')
                    df = pd.read_csv(file,sep=' ',header=['ROIs', 'label','R','G','B','A'])
                    df.dropna(subset=['R','G','B','A'], how='any', inplace=True)
                    os.makedirs(opj(opd(template_atlas_file), 'for_pet'), exist_ok=True)
                    df.to_csv(
                        opj(opd(template_atlas_file), 'for_pet', opb(template_atlas_file).replace('.nii.gz', '_gray.csv')),
                        index=False,
                        sep='\t')
                    atlas_label = opj(opd(template_atlas_file), 'for_pet',
                                      opb(template_atlas_file).replace('.nii.gz', '_gray.csv'))
                reference_signal = getlabels.gettype(template_atlas_file,
                                                     selected_atlases['level'],
                                                     opj(path_label_code,selected_atlases["name"] + '_gray.csv'))

                atlas_label = opj(path_label_code,selected_atlases["name"] + '_gray.csv')
            else:
                json_file = opj(MAIN_PATH, "Atlases_library", 'atlas_config_V3.json')
                f = open(json_file)
                info = json.load(f)

                l1, l2, newline,newfile = baseline_def.check_atlasjson(baseline, selected_atlases['name'], selected_atlases['level'],
                                                               template_atlas_file,info)

                if l1 != None:
                    toget = info['atlas_definitions'][selected_atlases['name'] + '[' + str(selected_atlases['level']) + ']'][
                        'regions']
                    df = pd.DataFrame(toget)
                    df.rename(columns={'label': 'ROIs', 'region': 'labels'}, inplace=True)

                    if not newline.empty:
                        df = pd.concat([df, newline], ignore_index=True)
                    os.makedirs(opj(opd(template_atlas_file), 'for_pet'), exist_ok=True)
                    df.to_csv(opj(opd(template_atlas_file), 'for_pet', opb(template_atlas_file).replace('.nii.gz', '_gray.csv')), index=False,
                              sep='\t')
                    if newfile:
                        atlas_file = newfile
                    atlas_label = opj(opd(template_atlas_file), 'for_pet', opb(template_atlas_file).replace('.nii.gz', '_gray.csv'))
                    reference_signal['label'] = l1
        else:
            reference_signal=''

    study_template_atlas_folder = opj(bids_dir,'derivatives','sty_template')


    # Set the environment variable for the current process
    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
    overwrite = ''
    if overwrite_option == True:
        overwrite = ' -overwrite'
    # Usage
    type_norm = check_nii.normalize_anat_type(type_norm)
    T_PET = check_nii.normalize_anat_type(T_PET)


    if isotope == 'unknown':
        if tracor in known_tracor:
            isotope = known_tracor[tracor]


    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session)
        run_cmd.printcolor(nl, 'HEADER')

        # PET organization folders ###########################################################################
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

            os.makedirs(path,exist_ok=True)


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

        # required for quantitative image (SUV)
        # if unknown (or useless) leave it
        Bq_inj = 'unknown'  # injection of FDG in MBq ex: 145
        weight = 'unknown'  # subject's weight in kg  ex: 7.98
        bloodc = 'unknown'  # in g/L # contraction of substance in blood

        info = ['weight_kg', 'tracor_qte_MBq', 'blood_concentration']

        if opi(opj(bids_dir, 'participants.tsv')):
            name = ID.split('-')[1]
            sesname = Session.split('-')[1]
            infotsv = pd.read_csv(opj(bids_dir, 'participants.tsv'), sep='\t')

            mask = (infotsv['id'] == name) & (infotsv['session_name'] == sesname)
            values = {}
            for col in info:
                if col in infotsv.columns:
                    values[col] = float(infotsv.loc[mask, col].values[0])
                else:
                    values[col] = 'unknown'
            weight = values.get(info[0])
            Bq_inj = values.get(info[1])
            bloodc = values.get(info[2])

        if not Bq_inj == 'unknown' and not weight == 'unknown':
            if not bloodc == 'unknown':
                ratio = (Bq_inj / weight * 1000) / bloodc
            else:
                ratio = Bq_inj / weight * 1000
        else:
            ratio = 1

        if isotope == 'unknown':
            if tracor in known_tracor:
                isotope = known_tracor[tracor]

        dictionary = {"Subject": ID,
                      "Session": Session,
                      "Template": reference,
                      "atlas": selected_atlases['name'],
                      "scannerPosition": humanPosition,
                      "realPosition": animalPosition,
                      "Tracor": tracor,
                      "isotope": isotope,
                      "nb of ignored frames for motion correction": str(nbvoid),
                      "Dose injected (in MBq)": str(Bq_inj),
                      "Subject Weight (in Kg)": str(weight),
                      "Blood_sugar (in g/L)": str(bloodc),
                      "Ratio ": str(ratio),
                      }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(path_pet, 'Prepro_Settings.json'), "w") as outfile:
            outfile.write(json_object)


        loop1(ID, tracor, data_path, reference, humanPosition, animalPosition, anat_raw, Need_check_orient, ref_CT,
              Skip_step, diary_file, sing_afni, sing_fs)



    if 2 in Skip_step:
        nl = 'skip step ' + str(1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')

    else:
        diary_file = diaryfile.create(opj(study_template_atlas_folder, 'Study_template_pet'),
                                  'Create the study template for PET')
        studyTfile = studyT.make_template(study_template_atlas_folder, studyT_transform, tracor, all_ID, all_data_path, reference,
                  diary_file)

    if 3 in Skip_step:
        nl = 'skip step ' + str(1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')

    else:
        studyT.manual_crop(study_template_atlas_folder,studyTfile,'desc-crop',diary_file,sing_fs)


    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session)
        run_cmd.printcolor(nl, 'HEADER')

        # PET organization folders ###########################################################################
        (dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_raw_matrices, path_pet,
         dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
         dir_prepro_orig_process, dir_prepro_orig_postprocessed,
         dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
         dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
         dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
         dir_prepro_template_process, dir_prepro_template_postprocessed, dir_prepro_acpc_matrices,
         dir_prepro_orig_matrices) \
            = getpath.pet(data_path, reference)

        diary_file = diaryfile.create(opj(path_pet, str(ID) + ' session ' + str(Session) + 'block2'), nl)

        # link with the individual anat template ################################################################
        anat_session_used = Session
        anat_data_path_used = data_path

        if not IhaveanANAT:
            if not folderforTemplate_Anat :
                nl = 'ERROR: You need to provide either an anat image or a folder with template images to perform the registration'
                raise ValueError(run_cmd.error(nl, diary_file))
            else:
                (anat_subject, brainmask, G_mask, WBG_mask, V_mask, W_mask, dir_transfo, FS_dir,
                 dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.create(folderforTemplate_Anat,
                                                                                     False,
                                                                                     True, species,
                                                                                     T_PET, diary_file)

        else:
            (anat_raw,anat_subject,brainmask,G_mask, WBG_mask, W_mask,V_mask,dir_transfo) = chooseanat.retrieve4pet(ID, data_path,
                                                                                   Session, anat_pet_same_space,
                                                                                   True,True,
                                                                                   T_PET, diary_file)

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
                            = chooseanat.retrieve4pet(ID, candidate_path,ses,anat_pet_same_space,True,True,
                                                      T_PET, diary_file)
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


        loop2(ID,tracor,data_path, reference, anat_pet_same_space,IhaveanANAT,anat_raw,anat_subject,brainmask,G_mask, W_mask,V_mask,dir_transfo,
              templateT1w, templatemask,T_PET,do_manual, atlas_file,atlas_label,Need_check_MC,Need_check_coreg,nbvoid,studyTfile,'Linear',Skip_step,reference_signal, isotope, DV_time,Do_cut, duration,suv,bratio,ratio,strm,diary_file,sing_afni, sing_fs,sing_wb,sing_turku)




