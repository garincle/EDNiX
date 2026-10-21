# import
import os
import glob
import json
import inspect

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd, getpath, diaryfile, Load_EDNiX_requirement, check_nii
from modalities.anat import set_launcher
from modalities.anat.loop3 import backtonative
from modalities.fMRI import chooseanat
from modalities.pet import petutils as pu
from modalities.pet import loop1, loop2, studyT, coreg2anat, transfo
from atlases import getlabels

# kept for the launchers of the PET engineer (tracer -> isotope)
known_tracor = pu.KNOWN_TRACER


def _anat_files(ID, data_path, Session, all_Session, bids_dir, anat_pet_same_space, use_erode_V, use_erode_WM, T_PET,
                diary_file):
    """Anatomical files of the subject in AC-PC space (chooseanat.retrieve, as the functional pipeline), with the
    same fallback on the other sessions of the subject. Returns (anat dict, anat data_path, anat session)."""
    def get(path, ses):
        (anat_subject, brainmask, G_mask, WBG_mask, V_mask, W_mask, dir_transfo, FS_dir, dir_prepro, volumes_dir,
         labels_dir, masks_dir) = chooseanat.retrieve(ID, path, ses, anat_pet_same_space, use_erode_V, use_erode_WM,
                                                      T_PET, diary_file)
        return dict(ID=ID, acpc=anat_subject,
                    brain=opj(volumes_dir, ID + '_space-acpc_desc-SS_' + T_PET + '.nii.gz'),
                    mask=brainmask, gray=G_mask, white=W_mask, csf=V_mask, labels_dir=labels_dir,
                    dir_transfo=dir_transfo,
                    reorient=opj(dir_prepro, ID + '_space-raw_desc-reorient_' + T_PET + '.nii.gz'),
                    raw=opj(dir_prepro, ID + '_space-raw_desc-template_' + T_PET + '.nii.gz'),
                    surfaces=sorted(glob.glob(opj(opd(volumes_dir), 'surfaces', 'Native_resol', ID + '.[lr].*.surf.gii'))))
    try:
        anat = get(data_path, Session)
        if opi(anat['acpc']):
            return anat, data_path, Session
    except Exception:
        pass
    if anat_pet_same_space:
        raise FileNotFoundError(run_cmd.error(f'ERROR: anat_pet_same_space=True but no anatomy for {ID} in session '
                                              f'{Session}', diary_file))
    run_cmd.msg(f'INFO: no anat found for {ID} session {Session}, searching the other sessions of this subject...',
                diary_file, 'WARNING')
    for ses in sorted(set(all_Session)):
        candidate = opj(bids_dir, f'sub-{ID}', f'ses-{ses}')
        if not ope(candidate):
            continue
        try:
            anat = get(candidate, ses)
        except Exception:
            continue
        if opi(anat['acpc']):
            run_cmd.msg(f'INFO: anat of {ID} found in session {ses}, used for this PET session.', diary_file, 'OKGREEN')
            return anat, candidate, ses
    raise FileNotFoundError(run_cmd.error(f'ERROR: no anat found for {ID} in any session. Cannot proceed.', diary_file))


def _ct_record(ID, data_path, bids_dir, all_Session, reference, same_session_only):
    """Registration record of the CT pipeline for this subject (this session first). Returns (record, path)."""
    sessions = [data_path] + ([] if same_session_only else
                              [opj(bids_dir, f'sub-{ID}', f'ses-{s}') for s in sorted(set(all_Session))])
    for path in sessions:
        recs = sorted(glob.glob(opj(getpath.ct(path, reference)['orig_matrices'], '*_desc-registration_ct.json')))
        if recs:
            return transfo.read_record(recs[0]), recs[0]
    return None, None


def _ct_raw(ID, data_path, reference):
    """Raw CT of the session (for the orientation check / header copy), or ''."""
    recs = sorted(glob.glob(opj(getpath.ct(data_path, reference)['orig_matrices'], '*_desc-registration_ct.json')))
    if recs:
        return transfo.read_record(recs[0]).get('raw', '')
    return ''


def _selected_atlases(selected_atlases, list_atlas):
    """{atlas: [levels] or 'all'} from 'all', [name, level] or [[name, level|'all'], ...]."""
    if selected_atlases in ('all', None, ''):
        return {a: 'all' for a in list_atlas[0]}
    if isinstance(selected_atlases, (list, tuple)) and selected_atlases and not isinstance(selected_atlases[0], (list, tuple)):
        selected_atlases = [selected_atlases]
    out = {}
    for a in selected_atlases:
        name, level = a[0], (a[1] if len(a) > 1 else 'all')
        if level == 'all':
            out[name] = 'all'
        else:
            if out.get(name) == 'all':
                continue
            out.setdefault(name, []).append(int(level))
    return out


def preprocess_pet(Skip_step, MAIN_PATH, bids_dir, species, allinfo_study_c, endpet, tracer,
                   animalPosition, humanPosition, orientation,
                   T_PET, type_norm, creat_study_template, coregistration_longitudinal,
                   IhaveanANAT=True, anat_pet_same_space=False,
                   IhaveaCT=False, ct_pet_same_space=True,
                   template_PET='', creat_study_template_pet=False, studyT_transform='SyN',
                   isotope='unknown', PET_header_wrong=False, deoblique='WARP', dcmfolder='', endfiledcm='dcm',
                   nb_void=0, crop_for_MoCo='auto', redo_crop=False, MoCo_transform='Rigid',
                   coreg_metric='mattes', type_of_transform_template='SyN', do_manual=False,
                   SUV=True, bratio=True, SUV_window=None, Do_cut=False, duration=10,
                   srtm=False, reference_signal='unknown', DV_time=0, kinetic_end=None,
                   pvc_fwhm=None, pvc_method='RBV',
                   do_pet_to_func=False, endfmri='*_bold.nii.gz', T1_eq=0,
                   export_resolution='pet', export_dynamic=False,
                   selected_atlases='all', atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]],
                   reference='EDNiX', use_erode_WM_func_masks=True, use_erode_V_func_masks=True,
                   list_to_keep=[], list_to_remove=[], overwrite_option=True,
                   Need_check_orient=False, Need_check_MC=False, Need_check_coreg=False,
                   doBALSAspace=False, **kwargs):
    """
    Preprocess PET images of a BIDS study: reorientation, motion correction within and between runs,
    registration to the anatomy / CT / a template (scenarios of Illustration/pet4.svg), parametric images
    (SUV, ratio to the brain mean, time blocks), regional analysis, reference-region kinetic models
    (SRTM, Logan) and export to the AC-PC space of the anatomy and to the template spaces. Every
    computation is done in the PET space (orig); the outputs follow Illustration/pet2b.svg.

    Steps (Skip_step): 1 reorientation, 2 PET study template, 3 manual crop of the study template,
    4 within-run motion correction, 5 between-runs motion correction, 6 registration and atlases / masks
    in the PET space, 7 parametric images and regional values, 8 kinetic models, 'pvc' partial-volume
    correction, 9 PET during the fMRI runs, 10 export to the anatomical and template spaces.

    Parameters
    ----------
    Skip_step : list
        Steps NOT to run (numbers above, or 'pvc').
    MAIN_PATH : str
        Folder of the EDNiX installation (the one that contains Atlases_library).
    bids_dir : str
        Path to the BIDS dataset.
    species : str
        Species of the dataset ('Macaque', 'Human'...), as in the anatomical pipeline.
    allinfo_study_c : pandas.DataFrame
        BIDS table of the study (Tools.load_bids.Load_BIDS_to_pandas with modalities=['pet']).
    endpet : str
        Glob pattern of the raw PET images in sub-*/ses-*/pet, e.g. '*_trc-FDG_*pet.nii.gz'.
    tracer : str
        Tracer name, as in the BIDS trc- entity (FDG, RACLO, PE2I...); names the study template and the
        session files.
    animalPosition, humanPosition : list
        Position of the animal and position declared to the scanner, as in the anatomical pipeline
        (['']: humanlike / PatientPosition of the sidecar). One value, or one per run.
    orientation : str
        Orientation code forced on the PET header (as `orientation` of the anatomical pipeline); '' to
        derive it from the positions.
    T_PET : str
        Anatomical image the PET is registered to ('T1w' or 'T2w').
    type_norm : str
        Anatomical contrast of the normalisation of the anatomical pipeline (same value as there).
    creat_study_template : bool
        Whether the anatomical pipeline used a study template (same value as there; anatomy -> template bridge).
    coregistration_longitudinal : bool
        Whether the anatomical pipeline used a longitudinal registration (same value as there).
    IhaveanANAT : bool, optional
        The subjects have an anatomical MRI processed by the anatomical pipeline (scenario 1 of pet4.svg).
    anat_pet_same_space : bool, optional
        PET and anatomy acquired in the same session and scanner space (hybrid PET-MR, scenario 1v1): the PET
        is reoriented like the anatomy and follows its AC-PC matrix, without registration.
    IhaveaCT : bool, optional
        A CT processed by the CT pipeline (modalities/ct) is used: PET -> CT -> anatomy (1v2b) or
        CT template (2v1, 2v2).
    ct_pet_same_space : bool, optional
        PET and CT share the scanner space (PET/CT scanner): no PET -> CT registration.
    template_PET : str, optional
        PET template in the space of the reference template (scenario 3v1, no anatomy, no CT).
    creat_study_template_pet : bool, optional
        Build a PET study template (step 2): used to crop the images for the motion correction and, without
        anatomy, CT or PET template, as the common space (3v2, no atlas).
    studyT_transform : str, optional
        ANTs transform of the PET study template registrations ('SyN').
    isotope : str, optional
        Radionuclide ('F-18', 'C-11'...); 'unknown': TracerRadionuclide of the sidecar, else from the tracer name.
    PET_header_wrong : bool, optional
        The reconstruction wrote a wrong orientation: the orientation code of the raw anatomy (or CT) of the same
        scanner space is copied before the reorientation.
    deoblique : str, optional
        'WARP' resamples an oblique PET on a cardinal grid (3dWarp, as the anatomy); 'header' drops the obliquity.
    dcmfolder, endfiledcm : str, optional
        DICOM folder of the PET and file extension, to read the isotope and scan start for the SIF (optional).
    nb_void : int, optional
        Number of first frames not registered by the motion correction (they get the transform of the first
        registered frame).
    crop_for_MoCo : str, optional
        Head box used by the motion correction: 'auto' (crop of the study template if it exists, else drawn on the
        session), 'manual' (drawn on every session, freeview), 'none' (whole field of view).
    redo_crop : bool, optional
        Draw the crops again (they are kept from one run to the next otherwise).
    MoCo_transform : str, optional
        ANTs transform of the motion correction ('Rigid').
    coreg_metric : str, optional
        Metric of the rigid registrations (motion, PET -> anatomy / CT): 'mattes' (mutual information).
    type_of_transform_template : str, optional
        ANTs transform of the PET -> PET template registration (scenario 3v1), e.g. 'SyN'.
    do_manual : bool, optional
        Manual correction of the PET -> anatomy registration in freeview (scenario 1v2a).
    SUV : bool, optional
        Compute the SUV (body weight and injected radioactivity from the sidecar / sessions.tsv /
        participants.tsv) and, with the blood glucose, the glucose-corrected SUV.
    bratio : bool, optional
        Compute the ratio of the PET to its mean over the brain mask.
    SUV_window : list, optional
        [start, end] in minutes of the frames averaged for the SUV and ratio images of a dynamic run (frames whose
        middle is in the window); None = all frames.
    Do_cut : bool, optional
        Also write the mean of every time block of `duration` minutes (divided by its brain mean).
    duration : float, optional
        Duration (min) of the time blocks of Do_cut.
    srtm : bool, optional
        Reference-region kinetic models of the dynamic runs (step 8): SRTM (regional bfmsrtm, voxel-wise imgbfbp,
        imgsrtm) and Logan DVR (imgdv) with the Turku PET Centre tools.
    reference_signal : dict or str, optional
        Reference region (atlases/getlabels.py): {'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0,
        'labels': [7, 1007] or a name pattern or 'default', 'name': 'Cer'}, {'type': 'mask', 'name': 'white'},
        or 'unknown' for the interactive choice (saved for the next subjects).
    DV_time : float, optional
        Start (min) of the linear part of the Logan plot (imgdv); 0 = no DVR image.
    kinetic_end : float, optional
        End (min) of the kinetic fits; None = all frames.
    pvc_fwhm : float or list, optional
        Spatial resolution of the scanner (PSF FWHM, mm): partial-volume correction (PETPVC) of the static images;
        None = no correction.
    pvc_method : str, optional
        PETPVC method ('RBV', 'GTM', 'MG'...).
    do_pet_to_func : bool, optional
        Simultaneous PET-MR: mean PET during every fMRI run of the session (step 9).
    endfmri : str, optional
        Glob pattern of the raw fMRI runs in sub-*/ses-*/func (step 9).
    T1_eq : int, optional
        Number of first fMRI volumes removed by the functional pipeline (step 9).
    export_resolution : str, optional
        Grid of the exported images: 'pet' (target space at the PET voxel size) or 'target' (anatomy / template grid).
    export_dynamic : bool, optional
        Also export the motion-corrected dynamic images (large files).
    selected_atlases : str or list, optional
        Atlases of the regional analysis: 'all' (every level of every atlas of atlas_followers) or
        [[name, level or 'all'], ...].
    atlas_followers : list
        Atlases of the anatomical pipeline ([[names], [types], [levels], [...]]), as there.
    reference : str, optional
        Reference template ('EDNiX').
    use_erode_WM_func_masks, use_erode_V_func_masks : bool, optional
        Use the eroded white-matter / ventricle masks (reference region 'white', 'csf').
    list_to_keep, list_to_remove : list
        (subject, session) pairs to process / exclude.
    overwrite_option : bool, optional
        Overwrite the AFNI outputs.
    Need_check_orient, Need_check_MC, Need_check_coreg : bool, optional
        Open freeview to check the orientation, the motion correction and the registrations.
    doBALSAspace : bool, optional
        Same value as in the anatomical pipeline.
    **kwargs : dict
        Unknown keyword arguments are reported and ignored.
    """
    if kwargs:
        run_cmd.printcolor('WARNING: unknown argument(s) ignored by preprocess_pet: ' + ', '.join(sorted(kwargs)) +
                           ' — check the spelling in your launcher.', 'WARNING')
    Skip_step = list(Skip_step)

    (FS_refs, template_dir, reference, balsa_folder, BALSAname, balsa_brainT1, BASE_atlas_folder, BASE_template,
     BASE_SS, BASE_mask, BASE_Gmask, BASE_WBGmask, BASE_Wmask, BASE_Vmask, CSF, GM, WM, Aseg_ref, list_atlas,
     path_label_code, all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
     fs_tools, reftemplate_path, MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH, bids_dir, allinfo_study_c,
                                                                                    species, list_to_keep,
                                                                                    list_to_remove, reference, type_norm,
                                                                                    '', '', atlas_followers, doBALSAspace)

    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, sing_turku, _, _ = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')

    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
    type_norm = check_nii.normalize_anat_type(type_norm)
    T_PET = check_nii.normalize_anat_type(T_PET)
    if crop_for_MoCo not in ('auto', 'manual', 'none'):
        raise ValueError(f"crop_for_MoCo={crop_for_MoCo!r}: expected 'auto', 'manual' or 'none'")
    if export_resolution not in ('pet', 'target'):
        raise ValueError(f"export_resolution={export_resolution!r}: expected 'pet' or 'target'")

    # atlases of the regional analysis and their label tables
    levels = _selected_atlases(selected_atlases, list_atlas)
    atlas_names = list(levels)
    label_tables = {}
    for a in atlas_names:
        try:
            label_tables[a] = getlabels.label_table(path_label_code, a)
        except FileNotFoundError as e:
            run_cmd.printcolor('WARNING: ' + str(e) + ' (regions named label_<id>)', 'WARNING')
            label_tables[a] = None
    ref_atlas = {'name': (reference_signal.get('atlas') if isinstance(reference_signal, dict) else None) or
                 (atlas_names[0] if atlas_names else ''),
                 'level': (reference_signal.get('level', 0) if isinstance(reference_signal, dict) else 0)}
    if ref_atlas['name'] and ref_atlas['name'] not in label_tables:
        # atlas of the reference region not in selected_atlases: brought into the PET space too
        label_tables[ref_atlas['name']] = getlabels.label_table(path_label_code, ref_atlas['name'])
        levels[ref_atlas['name']] = [ref_atlas['level']]
        atlas_names = list(levels)

    study_dir = opj(bids_dir, 'sty_template', 'pet')
    study_template = studyT.template_name(study_dir, tracer)
    reference_json = opj(study_dir, f'reference_region_trc-{tracer}.json')

    opts = dict(Skip_step=Skip_step, T_PET=T_PET, ct_pet_same_space=ct_pet_same_space, nb_void=int(nb_void),
                crop_for_MoCo=crop_for_MoCo, redo_crop=redo_crop, MoCo_transform=MoCo_transform,
                coreg_metric=coreg_metric, type_of_transform_template=type_of_transform_template,
                studyT_transform=studyT_transform, do_manual=do_manual, Need_check_MC=Need_check_MC,
                Need_check_coreg=Need_check_coreg, SUV=SUV, bratio=bratio, SUV_window=SUV_window, Do_cut=Do_cut,
                duration=duration, srtm=srtm, reference_signal=reference_signal, DV_time=DV_time,
                kinetic_end=kinetic_end, pvc_fwhm=pvc_fwhm, pvc_method=pvc_method, do_pet_to_func=do_pet_to_func,
                endfmri=endfmri, T1_eq=int(T1_eq), export_resolution=export_resolution, export_dynamic=export_dynamic,
                isotope=isotope)

    frame = inspect.currentframe()
    _loc = inspect.getargvalues(frame).locals
    args_dict = {k: _loc[k] for k in inspect.signature(preprocess_pet).parameters if k in _loc}

    # ======================================================================================== block 1 (step 1)
    sessions = []
    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session) + ' (PET block 1)'
        run_cmd.printcolor(nl, 'HEADER')
        paths = getpath.pet(data_path, reference)
        os.makedirs(paths['path'], exist_ok=True)
        diary_file = diaryfile.create(opj(paths['path'], str(ID) + ' session ' + str(Session)), nl)
        launcher_parameters = diaryfile.create(opj(paths['path'], str(ID) + 'launcher_parameters session ' + str(Session)), nl)
        with open(launcher_parameters, "w") as f:
            f.write("Function: preprocess_pet\nEffective parameters at runtime:\n\n")
            for k, v in args_dict.items():
                if k not in ('frame', 'allinfo_study_c'):
                    f.write(f"- {k} = {v}\n")

        raw_list = pu.list_raw_images(paths['path'], endpet, 'pet')
        if not raw_list:
            run_cmd.msg('WARNING: no PET image found with ' + opj(paths['path'], endpet) + ': session skipped',
                        diary_file, 'WARNING')
            continue
        info0 = pu.raw_sidecar(raw_list[0])
        iso = isotope
        if iso in ('unknown', '', None):
            iso = pu.normalise_isotope(info0.get('TracerRadionuclide')) or pu.KNOWN_TRACER.get(str(tracer).upper())
            if iso is None and srtm:
                raise ValueError(run_cmd.error('ERROR: isotope unknown (no TracerRadionuclide in the sidecar, tracer '
                                               f'{tracer} not in the table): give isotope in the launcher', diary_file))
        opts['isotope'] = iso
        session = pu.session_base(pu.base_name(raw_list[0]))

        # images for the orientation check / header copy
        check_ref, header_ref = '', ''
        if IhaveanANAT:
            try:
                anat, _, _ = _anat_files(ID, data_path, Session, all_Session, bids_dir, anat_pet_same_space,
                                         use_erode_V_func_masks, use_erode_WM_func_masks, T_PET, diary_file)
                check_ref = anat['reorient'] if anat_pet_same_space else anat['acpc']
                header_ref = anat['raw'] if anat_pet_same_space else ''
            except FileNotFoundError:
                if not (1 in Skip_step):
                    raise
        if IhaveaCT and ct_pet_same_space and not header_ref:
            header_ref = _ct_raw(ID, data_path, reference)

        values = pu.subject_values(bids_dir, ID, Session, info0)
        settings = {"Subject": ID, "Session": Session, "Template": reference, "Atlases": levels,
                    "Tracer": tracer, "Isotope": iso, "humanPosition": humanPosition, "animalPosition": animalPosition,
                    "orientation": orientation, "Number of frames not registered (nb_void)": int(nb_void),
                    "Injected radioactivity (MBq)": values['dose_MBq'], "Body weight (kg)": values['weight_kg'],
                    "Blood glucose (mmol/L)": values['glucose_mmol_L'], "Sources": values['sources'],
                    "Runs": [opb(r) for r in raw_list]}
        with open(opj(paths['path'], 'Prepro_Settings.json'), 'w') as f:
            json.dump(settings, f, indent=2, default=str)

        if 1 not in Skip_step:
            run_cmd.msg(pu.step_header(1), diary_file, 'HEADER')
        reoriented = loop1.run(raw_list, paths, humanPosition, animalPosition, orientation, check_ref, header_ref,
                               Need_check_orient, PET_header_wrong, deoblique, iso, srtm, (dcmfolder, endfiledcm),
                               Skip_step, diary_file, sing_afni, sing_fs, sing_turku)
        sessions.append(dict(ID=ID, Session=Session, data_path=data_path, max_ses=max_ses, paths=paths,
                             session=session, raw=raw_list, reoriented=reoriented, diary=diary_file, isotope=iso))

    if not sessions:
        raise ValueError('no PET session to process: check endpet, list_to_keep and list_to_remove')

    # ======================================================================================== steps 2 and 3
    os.makedirs(study_dir, exist_ok=True)
    study_diary = diaryfile.create(opj(study_dir, 'Study_template_pet'), 'Study template of the PET images')
    if creat_study_template_pet:
        if 2 in Skip_step:
            run_cmd.msg('INFO: skip step 2', study_diary, 'OKGREEN')
        else:
            run_cmd.msg(pu.step_header(2), study_diary, 'HEADER')
            items = []
            for s in sessions:
                first = s['reoriented'][0]
                mean = first.replace('_desc-reorient_', '_desc-reorient-mean_') if pu.is_dynamic(first) else first
                items.append((mean, s['paths']['orig_matrices'], s['session']))
            studyT.make_template(study_dir, items, tracer, studyT_transform, study_diary)
        if 3 in Skip_step or crop_for_MoCo == 'none':
            run_cmd.msg('INFO: skip step 3', study_diary, 'OKGREEN')
        else:
            run_cmd.msg(pu.step_header(3), study_diary, 'HEADER')
            studyT.manual_crop(study_template, study_diary, sing_fs, redo=redo_crop)
    study_crop = studyT.crop_name(study_template) if opi(studyT.crop_name(study_template)) else ''

    # ======================================================================================== block 2 (steps 4-10)
    study_template_atlas_folder = opj(bids_dir, 'sty_template')
    for s in sessions:
        ID, Session, data_path, paths, diary_file = s['ID'], s['Session'], s['data_path'], s['paths'], s['diary']
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session) + ' (PET block 2)'
        run_cmd.printcolor(nl, 'HEADER')
        opts['isotope'] = s['isotope']

        anat, packages, ct_record = None, None, None
        if IhaveanANAT:
            anat, anat_path, anat_ses = _anat_files(ID, data_path, Session, all_Session, bids_dir, anat_pet_same_space,
                                                    use_erode_V_func_masks, use_erode_WM_func_masks, T_PET, diary_file)
            packages = backtonative.get(ID, anat_path, bids_dir, anat_ses, s['max_ses'], 'space-acpc_desc-SS',
                                        type_norm, BASE_SS, BASE_mask, BASE_atlas_folder, study_template_atlas_folder,
                                        creat_study_template, coregistration_longitudinal, reference, 'Final', species,
                                        diary_file)
        if IhaveaCT:
            ct_record, _ = _ct_record(ID, data_path, bids_dir, all_Session, reference, ct_pet_same_space)
        sc = coreg2anat.scenario(IhaveanANAT, anat_pet_same_space, IhaveaCT,
                                 None if ct_record is None else sorted(ct_record['targets']),
                                 template_PET, creat_study_template_pet and opi(study_template))
        run_cmd.msg('INFO: registration scenario ' + sc + ' (Illustration/pet4.svg)', diary_file, 'OKGREEN')

        template_files = dict(name=reference, label_base=opj(BASE_atlas_folder, species + '_seg-'),
                              atlas_dir=BASE_atlas_folder, species=species, brain=BASE_mask, gray=BASE_Gmask,
                              white=BASE_Wmask if use_erode_WM_func_masks else BASE_Wmask.replace('_desc-erod-', '_desc-'),
                              csf=BASE_Vmask if use_erode_V_func_masks else BASE_Vmask.replace('_desc-erod-', '_desc-'))
        target_files = {}
        for pkg in packages or []:
            target_files[pkg[0]] = dict(label_base=pkg[3], list_atlas=atlas_names, masks=[pkg[2]])
        if sc in ('2v1', '3v1'):
            target_files[reference] = dict(label_base=template_files['label_base'], list_atlas=atlas_names,
                                           masks=[BASE_mask])
        reg = dict(scenario=sc, anat=anat, packages=packages, ct_record=ct_record, template_PET=template_PET,
                   reference=reference, study_template=study_template, study_crop=study_crop,
                   template_files=template_files, target_files=target_files, list_atlas=atlas_names,
                   label_tables=label_tables, levels=levels, ref_atlas=ref_atlas, reference_json=reference_json)
        loop2.run(ID, Session, data_path, bids_dir, paths, s['session'], s['reoriented'], s['raw'], opts, reg,
                  diary_file, sing_afni, sing_fs, sing_wb, sing_turku)


# name of the PET engineer's entry point
preprocess_data = preprocess_pet
