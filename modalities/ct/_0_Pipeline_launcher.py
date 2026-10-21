# import
import os
import inspect

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd, getpath, diaryfile, Load_EDNiX_requirement, check_nii
from modalities.anat import set_launcher
from modalities.anat.loop3 import backtonative
from modalities.ct import reorient, coreg2anat
from modalities.pet import petutils as pu
from modalities.pet import studyT
from modalities.pet._0_Pipeline_launcher import _anat_files


def preprocess_ct(Skip_step, MAIN_PATH, bids_dir, species, allinfo_study_c, endct,
                  animalPosition, humanPosition, orientation,
                  T_CT, type_norm, creat_study_template, coregistration_longitudinal,
                  IhaveanANAT=True, anat_ct_same_space=False,
                  template_CT='', creat_study_template_ct=False, studyT_transform='SyN',
                  type_of_transform_template='SyN', coreg_metric='mattes',
                  CT_header_wrong=False, deoblique='WARP', do_manual=False,
                  Need_check_orient=False, Need_check_coreg=False,
                  reference='EDNiX', atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
                  use_erode_WM_func_masks=True, use_erode_V_func_masks=True,
                  list_to_keep=[], list_to_remove=[], overwrite_option=True, doBALSAspace=False, **kwargs):
    """
    Preprocess the CT images of a BIDS study for the PET pipeline (Illustration/pet4.svg): reorientation
    and registration to the anatomy (scenario 1v2b), to a CT template (2v1) or to a CT study template (2v2).
    The registration record written in ct/orig/matrices is read by the PET pipeline (IhaveaCT=True). Run it
    before the PET pipeline.

    Steps (Skip_step): 1 reorientation, 2 CT study template, 3 registration and export (CT in the AC-PC /
    template space, anatomy in the CT space, quality control).

    Parameters
    ----------
    Skip_step : list
        Steps NOT to run.
    MAIN_PATH : str
        Folder of the EDNiX installation (the one that contains Atlases_library).
    bids_dir : str
        Path to the BIDS dataset.
    species : str
        Species of the dataset, as in the anatomical pipeline.
    allinfo_study_c : pandas.DataFrame
        BIDS table of the study (Tools.load_bids.Load_BIDS_to_pandas with modalities=['ct']).
    endct : str
        Glob pattern of the raw CT images in sub-*/ses-*/ct, e.g. '*_ct.nii.gz'.
    animalPosition, humanPosition : list
        Position of the animal and position declared to the scanner, as in the anatomical pipeline; use the
        same values as the PET of the same PET/CT acquisition.
    orientation : str
        Orientation code forced on the CT header; '' to derive it from the positions.
    T_CT : str
        Anatomical image the CT is registered to ('T1w' or 'T2w').
    type_norm : str
        Anatomical contrast of the normalisation of the anatomical pipeline (same value as there).
    creat_study_template : bool
        Whether the anatomical pipeline used a study template (same value as there).
    coregistration_longitudinal : bool
        Whether the anatomical pipeline used a longitudinal registration (same value as there).
    IhaveanANAT : bool, optional
        The subjects have an anatomical MRI processed by the anatomical pipeline (scenario 1v2b).
    anat_ct_same_space : bool, optional
        CT and anatomy acquired in the same scanner space: the CT follows the AC-PC matrix of the anatomy.
    template_CT : str, optional
        CT template in the space of the reference template (scenario 2v1, no anatomy).
    creat_study_template_ct : bool, optional
        Build a CT study template (step 2) and register the CTs to it (scenario 2v2, no anatomy, no CT template).
    studyT_transform : str, optional
        ANTs transform of the study template registrations ('SyN').
    type_of_transform_template : str, optional
        ANTs transform of the CT -> CT template registration (2v1), e.g. 'SyN'.
    coreg_metric : str, optional
        Metric of the registrations ('mattes': mutual information).
    CT_header_wrong : bool, optional
        The reconstruction wrote a wrong orientation: the orientation code of the raw anatomy of the same scanner
        space is copied before the reorientation (with anat_ct_same_space).
    deoblique : str, optional
        'WARP' resamples an oblique CT on a cardinal grid (3dWarp, as the anatomy); 'header' drops the obliquity.
    do_manual : bool, optional
        Manual correction of the CT -> anatomy registration in freeview.
    Need_check_orient, Need_check_coreg : bool, optional
        Open freeview to check the orientation and the registration.
    reference : str, optional
        Reference template ('EDNiX').
    atlas_followers : list
        Atlases of the anatomical pipeline, as there.
    use_erode_WM_func_masks, use_erode_V_func_masks : bool, optional
        Same values as in the PET pipeline (masks of the anatomy).
    list_to_keep, list_to_remove : list
        (subject, session) pairs to process / exclude.
    overwrite_option : bool, optional
        Overwrite the AFNI outputs.
    doBALSAspace : bool, optional
        Same value as in the anatomical pipeline.
    **kwargs : dict
        Unknown keyword arguments are reported and ignored.
    """
    if kwargs:
        run_cmd.printcolor('WARNING: unknown argument(s) ignored by preprocess_ct: ' + ', '.join(sorted(kwargs)) +
                           ' — check the spelling in your launcher.', 'WARNING')
    Skip_step = list(Skip_step)
    (FS_refs, template_dir, reference, balsa_folder, BALSAname, balsa_brainT1, BASE_atlas_folder, BASE_template,
     BASE_SS, BASE_mask, BASE_Gmask, BASE_WBGmask, BASE_Wmask, BASE_Vmask, CSF, GM, WM, Aseg_ref, list_atlas,
     path_label_code, all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
     fs_tools, reftemplate_path, MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH, bids_dir, allinfo_study_c,
                                                                                    species, list_to_keep,
                                                                                    list_to_remove, reference, type_norm,
                                                                                    '', '', atlas_followers, doBALSAspace)
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, _, _ = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')
    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
    type_norm = check_nii.normalize_anat_type(type_norm)
    T_CT = check_nii.normalize_anat_type(T_CT)
    if IhaveanANAT:
        sc = 'anat_same_space' if anat_ct_same_space else '1v2b'
    elif template_CT:
        sc = '2v1'
    elif creat_study_template_ct:
        sc = '2v2'
    else:
        raise ValueError('no anatomy, no CT template and no CT study template: nothing to register the CT to')

    frame = inspect.currentframe()
    args_dict = inspect.getargvalues(frame).locals.copy()
    study_dir = opj(bids_dir, 'sty_template', 'ct')
    study_template = studyT.template_name(study_dir, '', 'ct')

    sessions = []
    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session) + ' (CT)'
        run_cmd.printcolor(nl, 'HEADER')
        paths = getpath.ct(data_path, reference)
        os.makedirs(paths['path'], exist_ok=True)
        diary_file = diaryfile.create(opj(paths['path'], str(ID) + ' session ' + str(Session)), nl)
        launcher_parameters = diaryfile.create(opj(paths['path'], str(ID) + 'launcher_parameters session ' + str(Session)), nl)
        with open(launcher_parameters, "w") as f:
            f.write("Function: preprocess_ct\nEffective parameters at runtime:\n\n")
            for k, v in args_dict.items():
                if k not in ('frame', 'allinfo_study_c'):
                    f.write(f"- {k} = {v}\n")
        raw = pu.list_raw_images(paths['path'], endct, 'ct')
        if not raw:
            run_cmd.msg('WARNING: no CT found with ' + opj(paths['path'], endct) + ': session skipped', diary_file,
                        'WARNING')
            continue
        if len(raw) > 1:
            run_cmd.msg('WARNING: several CTs in this session; the first is used: ' + opb(raw[0]), diary_file, 'WARNING')
        raw = raw[0]
        base = pu.base_name(raw, 'ct')
        anat = anat_path = anat_ses = None
        if IhaveanANAT:
            anat, anat_path, anat_ses = _anat_files(ID, data_path, Session, all_Session, bids_dir, anat_ct_same_space,
                                                    use_erode_V_func_masks, use_erode_WM_func_masks, T_CT, diary_file)
        ct_img = pu.out_name(paths['orig_process'], base, desc='reorient', suffix='ct')
        if 1 in Skip_step:
            run_cmd.msg('INFO: skip step 1', diary_file, 'OKGREEN')
        else:
            run_cmd.msg('##  Working on step 1 (CT pipeline).  ##', diary_file, 'HEADER')
            info = pu.raw_sidecar(raw)
            hp = (humanPosition[0] if isinstance(humanPosition, (list, tuple)) and humanPosition else humanPosition) \
                or info.get('PatientPosition', '')
            ap = animalPosition[0] if isinstance(animalPosition, (list, tuple)) and animalPosition else animalPosition
            check = (anat['reorient'] if anat_ct_same_space else anat['acpc']) if anat else ''
            reorient.doit(paths['orig_process'], raw, hp, ap, orientation, check, Need_check_orient, CT_header_wrong,
                          diary_file, sing_afni, sing_fs, deoblique=deoblique,
                          header_ref=anat['raw'] if (anat and anat_ct_same_space) else '')
        sessions.append(dict(ID=ID, Session=Session, max_ses=max_ses, paths=paths, base=base, raw=raw, ct=ct_img,
                             anat=anat, anat_path=anat_path, anat_ses=anat_ses, diary=diary_file))

    if sc == '2v2':
        os.makedirs(study_dir, exist_ok=True)
        sdiary = diaryfile.create(opj(study_dir, 'Study_template_ct'), 'Study template of the CT images')
        if 2 in Skip_step:
            run_cmd.msg('INFO: skip step 2', sdiary, 'OKGREEN')
        else:
            studyT.make_template(study_dir, [(s['ct'], s['paths']['orig_matrices'], s['base']) for s in sessions], '',
                                 studyT_transform, sdiary, suffix='ct')

    study_template_atlas_folder = opj(bids_dir, 'sty_template')
    for s in sessions:
        if 3 in Skip_step:
            run_cmd.msg('INFO: skip step 3', s['diary'], 'OKGREEN')
            continue
        run_cmd.msg('##  Working on step 3 (CT pipeline).  ##', s['diary'], 'HEADER')
        packages = None
        if s['anat']:
            packages = backtonative.get(s['ID'], s['anat_path'], bids_dir, s['anat_ses'], s['max_ses'],
                                        'space-acpc_desc-SS', type_norm, BASE_SS, BASE_mask, BASE_atlas_folder,
                                        study_template_atlas_folder, creat_study_template, coregistration_longitudinal,
                                        reference, 'Final', species, s['diary'])
        coreg2anat.register(sc, s['ct'], s['raw'], s['base'], s['paths'], s['anat'], packages, template_CT, reference,
                            study_template, T_CT, coreg_metric, type_of_transform_template, do_manual,
                            Need_check_coreg, bids_dir, s['diary'], sing_fs)
