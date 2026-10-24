"""
Block 2 of the PET pipeline for one subject / session (steps 4 to 10, see _0_Pipeline_launcher).
"""
import os
import glob
import json
import shutil

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet import studyT, MoCo, coreg2anat, preproPET, norm2template, roi_signal, transfo
from modalities.pet import cut_func, pvc as pvcmod
from modalities.pet import loop1
from atlases import getlabels


def _skip(step, Skip_step, diary_file):
    if step in Skip_step:
        run_cmd.msg('INFO: skip step ' + str(step), diary_file, 'OKGREEN')
        return True
    run_cmd.msg(pu.step_header(step), diary_file, 'HEADER')
    return False


def _mean_of(img, desc_mean):
    """Mean image written next to a step output (dynamic), or the image itself (static)."""
    m = img.replace('_desc-' + desc_mean.rsplit('-mean', 1)[0] + '_', '_desc-' + desc_mean + '_')
    return m if pu.is_dynamic(img) else img


def crop_for_session(crop_for_MoCo, study_crop, session, paths, first_mean, diary_file, sing_fs, redo_crop):
    """Crop box (file) in the space of the first run of the session, or None."""
    if crop_for_MoCo == 'none':
        return None
    prefix = opj(paths['orig_matrices'], session + '_to_studyTemplate_')
    out = opj(paths['orig_formc'], session + '_desc-cropbox_mask.nii.gz')
    os.makedirs(paths['orig_formc'], exist_ok=True)
    if crop_for_MoCo == 'auto' and study_crop and opi(study_crop) and opi(prefix + transfo.ANTS_AFF):
        return studyT.crop_mask(study_crop, first_mean, out, transfo.registration(prefix))
    crop = opj(paths['orig_formc'], session + '_desc-crop_pet.nii.gz')
    studyT.manual_crop_session(first_mean, crop, diary_file, sing_fs, redo=redo_crop)
    return studyT.crop_mask(crop, first_mean, out)


def run(ID, Session, data_path, bids_dir, paths, session, reoriented, raw_list, opts, reg, diary_file,
        sing_afni, sing_fs, sing_wb, sing_turku):
    """
    opts: options of the launcher (dict); reg: registration inputs (scenario, anat, packages, ct_record,
    template_PET, reference, study_template, template_files, list_atlas, label_tables, levels).
    """
    Skip_step = opts['Skip_step']
    timings = [pu.frame_timing(pu.raw_sidecar(r), pu.n_frames(r)) for r in raw_list]
    infos = [pu.raw_sidecar(r) for r in raw_list]
    sifs = [loop1.sif_file(paths, r) for r in raw_list]
    sifs = [s if opi(s) else '' for s in sifs]
    for img in reoriented:
        if not opi(img):
            raise FileNotFoundError(f'{img} not found: run step 1 first')
    first_mean = _mean_of(reoriented[0], 'reorient-mean')

    # ---------------------------------------------------------------------------------------- 4 within-run MoCo
    within = []
    crop = None
    skip4 = _skip(4, Skip_step, diary_file) if any(pu.is_dynamic(i) for i in reoriented) else True
    for img, timing, sif in zip(reoriented, timings, sifs):
        base = pu.base_name(img)
        out = pu.out_name(paths['orig_process'], base, desc='MoCo-within')
        if not pu.is_dynamic(img):
            within.append(img)
            continue
        if skip4:
            if opi(out):
                within.append(out)
            else:
                # step 4 skipped and never run: the frames are used without within-run correction
                run_cmd.msg('WARNING: ' + opb(out) + ' not found (step 4 skipped): the reoriented frames are used '
                            'without within-run motion correction', diary_file, 'WARNING')
                within.append(img)
            continue
        within.append(out)
        if crop is None:
            crop = crop_for_session(opts['crop_for_MoCo'], reg.get('study_crop'), session, paths, first_mean,
                                    diary_file, sing_fs, opts['redo_crop']) or ''
        # the box of the session (first run) is resampled on the grid of each run by MoCo (world coordinates)
        MoCo.dynamic(img, paths['orig_process'], paths['orig_formc'], timing, opts['nb_void'], crop or None,
                     opts['Need_check_MC'], diary_file, sing_fs, opts['MoCo_transform'], opts['coreg_metric'],
                     sif=img.replace('.nii.gz', '.sif'))
        # motion plot in the QC folder of the subject (read by the interface)
        os.makedirs(paths['qc'], exist_ok=True)
        png = opj(paths['orig_formc'], base + '_desc-MoCo-within_motion.png')
        if opi(png):
            shutil.copyfile(png, opj(paths['qc'], opb(png)))

    # ---------------------------------------------------------------------------------------- 5 between-runs MoCo
    between = [pu.out_name(paths['orig_process'], pu.base_name(i), desc='MoCo-between') for i in within]
    session_mean = pu.out_name(paths['orig_process'], session, desc='MoCo-between-mean')
    if not _skip(5, Skip_step, diary_file):
        means = [_mean_of(i, 'MoCo-within-mean') if 'MoCo-within' in i else _mean_of(i, 'reorient-mean')
                 for i in within]
        if crop is None and len(within) > 1:
            crop = crop_for_session(opts['crop_for_MoCo'], reg.get('study_crop'), session, paths, first_mean,
                                    diary_file, sing_fs, opts['redo_crop']) or ''
        between, session_mean = MoCo.static(within, means, timings,
                                            [i.replace('.nii.gz', '.sif') if opi(i.replace('.nii.gz', '.sif')) else ''
                                             for i in within],
                                            session, paths['orig_process'], paths['orig_formc'], crop or None,
                                            diary_file, opts['MoCo_transform'], opts['coreg_metric'])

    # ---------------------------------------------------------------------------------------- 6 registration
    rec_file = coreg2anat.record_path(paths, session)
    if not _skip(6, Skip_step, diary_file):
        if not opi(session_mean):
            raise FileNotFoundError(f'{session_mean} not found: run step 5 first')
        study_template = reg.get('study_template')

        # The crop box of steps 4-5 is used in the MoCo metric only: the MoCo outputs are full FOV. ANTs
        # initialises a rigid registration from the centre of mass of the *whole* moving image, so a
        # full-FOV PET (bladder, liver, injection site, scanner bed) starts far from the head even when
        # the metric is restricted to a brain mask. Crop the session mean for the registration itself;
        # the transforms produced are world-to-world and remain valid for the full session mean, which
        # stays the reference of the record (used by to_pet_space, the QC and norm2template.export).
        # crop is set in step 4 (dynamic runs) or step 5 (multiple static runs); it is None / '' when
        # crop_for_MoCo='none' or when a single static run gives no crop. In those cases the full mean is
        # used, exactly as before.
        session_mean_for_reg = session_mean
        if not (crop and opi(crop)):
            # steps 4-5 skipped (or crop_for_MoCo='none'): the crop was not computed earlier in this call.
            # Recompute it here so the registration always sees a head-only moving image. crop_for_session
            # reuses the manual crop / study-template crop already on disk (no freeview if they exist) and
            # is idempotent for the cropbox mask.
            crop = crop_for_session(opts['crop_for_MoCo'], reg.get('study_crop'), session, paths, first_mean,
                                    diary_file, sing_fs, opts['redo_crop']) or ''
            if crop and opi(crop):
                run_cmd.msg('INFO: crop computed for step 6 alone (steps 4-5 not run): ' + opb(crop),
                            diary_file, 'OKGREEN')
        if crop and opi(crop):
            session_mean_for_reg = MoCo.crop_mean(
                session_mean, crop, session_mean.replace('.nii.gz', '_desc-forReg.nii.gz'))
            run_cmd.msg('INFO: registration estimated on the cropped session mean: ' +
                        opb(session_mean_for_reg) + ' (crop box: ' + opb(crop) + ')', diary_file, 'OKGREEN')
        else:
            run_cmd.msg('WARNING: no crop box available: the full session mean is used for the '
                        'registration (ANTs centre-of-mass initialisation may be biased by the scanner '
                        'bed / shoulders)', diary_file, 'WARNING')

        if reg['scenario'] == '3v2':
            # registration of the motion-corrected session mean to the study template. Done here (not in
            # coreg2anat.register, which just reads the prefix) so it can also use the cropped mean.
            transfo.register(study_template, session_mean_for_reg, coreg2anat.study_prefix(paths, session),
                             opts['studyT_transform'])

        anchor = None
        if reg['scenario'] == '1v1' or (reg['ct_record'] is not None and opts['ct_pet_same_space']):
            # same scanner space as the anatomy / CT: the motion-corrected session is brought back to the head
            # position of the frame acquired with it
            anchor, how = MoCo.scanner_anchor(reoriented, within, timings, paths['orig_formc'], opts['nb_void'],
                                              opts.get('same_space_window'))
            run_cmd.msg('INFO: same scanner space as the anatomy / CT: ' + how, diary_file, 'OKGREEN')

        coreg2anat.register(reg['scenario'], session_mean, session, paths, reg['anat'], reg['packages'],
                            reg['ct_record'], reg['template_PET'], reg['reference'], study_template,
                            opts['T_PET'], opts['ct_pet_same_space'], opts['coreg_metric'],
                            opts['type_of_transform_template'], opts['do_manual'], opts['Need_check_coreg'],
                            bids_dir, diary_file, sing_fs,
                            template_mask=(reg.get('template_files') or {}).get('brain', ''), anchor=anchor,
                            moving_for_estimate=session_mean_for_reg, coreg_transform=opts.get('coreg_transform', 'Rigid'),
                            coreg_init=opts.get('coreg_init', 'auto'),
                            pet_coreg_mask=opts.get('pet_coreg_mask', 'none'),
                            pet_mask_percentile=opts.get('pet_mask_percentile', 30.0),
                            pet_mask_dilate_mm=opts.get('pet_mask_dilate_mm', 3.0),
                            coreg_target=opts.get('coreg_target', 'anat'),
                            synthetic_fwhm=opts.get('synthetic_fwhm'))

        record = transfo.read_record(rec_file)
        record['pet_space'] = coreg2anat.to_pet_space(record, session, paths, reg['anat'], reg['template_files'],
                                                      reg['list_atlas'], diary_file, opts['T_PET'])
        transfo.write_record(rec_file, record)
    record = transfo.read_record(rec_file) if opi(rec_file) else None

    if record is not None and opi(session_mean) and os.path.getmtime(session_mean) > os.path.getmtime(rec_file):
        run_cmd.msg('WARNING: the session mean is newer than the registration record (steps 4-5 run again): run step 6 '
                    'again, the registration and the atlases in the PET space belong to the previous correction',
                    diary_file, 'WARNING')
    pet_space = (record or {}).get('pet_space', dict(atlases={}, masks={}, anat=''))
    atlases, masks = pet_space['atlases'], pet_space['masks']

    # ---------------------------------------------------------------------------------------- 7 parametric images
    missing = [b for b in between if not opi(b)]
    if not _skip(7, Skip_step, diary_file):
        if missing:
            raise FileNotFoundError(f'{missing[0]} not found: run steps 4 and 5 first')
        for img, timing, info in zip(between, timings, infos):
            # injected radioactivity of this run (sidecar), body weight / glucose of the session
            values = pu.subject_values(bids_dir, ID, Session, info)
            run_cmd.msg('INFO: ' + opb(img) + ': subject data for the SUV: ' + json.dumps(values, default=str),
                        diary_file, 'OKGREEN')
            preproPET.suv(img, paths, timing, info, values, opts['isotope'], masks, atlases, reg['label_tables'],
                          reg['levels'], opts['SUV'], opts['bratio'], opts['SUV_window'], opts['Do_cut'],
                          opts['duration'], diary_file)

    # ---------------------------------------------------------------------------------------- 8 kinetic models
    if opts['srtm'] and not _skip(8, Skip_step, diary_file):
        for img, timing, info in zip(between, timings, infos):
            if not pu.is_dynamic(img):
                continue
            atlas0 = reg['ref_atlas']
            present = []
            if atlas0['name'] in atlases:
                import numpy as np
                present = sorted(int(x) for x in np.unique(roi_signal.labels_volume(atlases[atlas0['name']],
                                                                                    atlas0.get('level', 0))) if x > 0)
            ref = getlabels.resolve(opts['reference_signal'], atlas0, present,
                                    reg['label_tables'].get(atlas0['name']), reg['reference_json'])
            if ref['type'] == 'atlas' and ref['atlas'] not in atlases:
                raise ValueError(f"reference region in atlas {ref['atlas']}, which is not in the PET space "
                                 f"(available: {sorted(atlases)})")
            sif = img.replace('.nii.gz', '.sif')
            preproPET.srtm(img, paths, timing, sif, info, opts['isotope'], masks, atlases, reg['label_tables'],
                           reg['levels'], ref, opts['DV_time'], diary_file, sing_turku, opts['kinetic_end'],
                           opts.get('DV_k2ref'))

    # ---------------------------------------------------------------------------------------- pvc
    if opts['pvc_fwhm'] is not None and not _skip('pvc', Skip_step, diary_file):
        for info in infos:
            if str(info.get('ReconFilterType', 'none')).lower() not in ('none', '', '["none"]') and \
                    info.get('ReconFilterSize'):
                run_cmd.msg(f"INFO: post-reconstruction filter {info.get('ReconFilterType')} "
                            f"{info.get('ReconFilterSize')} mm in the sidecar: pvc_fwhm must include it",
                            diary_file, 'WARNING')
                break
        for img in between:
            base = pu.base_name(img)
            for desc in ('SUV-mean', 'SUV', 'activity-mean', 'activity'):
                src = pu.out_name(paths['orig_post'], base, desc=desc)
                if opi(src) and not pu.is_dynamic(src):
                    out = pvcmod.correct(src, pet_space.get('probseg', {}), opts['pvc_fwhm'],
                                         pu.out_name(paths['orig_post'], base, desc=desc + '-PVC'), diary_file,
                                         opts['pvc_method'])
                    if atlases:
                        roi_signal.read(out, atlases, reg['levels'], reg['label_tables'], paths['orig_roi'], None,
                                        desc + '-PVC', diary_file)
                    break

    # ---------------------------------------------------------------------------------------- 9 PET during fMRI
    if opts['do_pet_to_func'] and not _skip(9, Skip_step, diary_file):
        for img, timing, info in zip(between, timings, infos):
            if pu.is_dynamic(img):
                out = cut_func.cut(data_path, opts['endfmri'], img, timing, info, opts['T1_eq'], masks.get('brain'),
                                   paths['orig_post'], diary_file)
                if atlases:
                    for o in out:
                        roi_signal.read(o, atlases, reg['levels'], reg['label_tables'], paths['orig_roi'], None,
                                        pu.strip_ext(o).split('_desc-')[1].rsplit('_', 1)[0], diary_file)

    # ---------------------------------------------------------------------------------------- 10 export
    if not _skip(10, Skip_step, diary_file):
        if record is None:
            raise FileNotFoundError(f'{rec_file} not found: run step 6 first')
        images = []
        for img in between:
            base = pu.base_name(img)
            for f in sorted(glob.glob(opj(paths['orig_post'], base + '_desc-*_pet.nii.gz'))):
                if pu.is_dynamic(f) and not opts['export_dynamic']:
                    continue
                # PET-space skull-stripped images (redone with the mask of each target), working files and
                # integer QC maps (fit at a bound, number of Logan points: not interpolated)
                if any(k in opb(f) for k in ('-SS_', 'desc-decaycorr', 'seg4D', '-err_', '-npoints_', 'dynamic-SS')):
                    continue
                images.append(f)
            if opts['export_dynamic']:
                images.append(img)
        maxvals = {'SUV': 2.5, 'SUV-mean': 2.5, 'bratio': 1.5, 'bratio-mean': 1.5}
        norm2template.export(images, record, session, paths, bids_dir, opts['export_resolution'], diary_file,
                             sing_wb, maxvals, reg['target_files'])
