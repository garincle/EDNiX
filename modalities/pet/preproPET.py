"""
Steps 7 and 8 of the PET pipeline for one run, in the PET space (orig/postprocessing):

suv  (step 7)  activity image of the run (static image, or duration-weighted mean of the frames of
               SUV_window for a dynamic run) -> SUV, SUVglu, ratio to the brain mean (bratio), their
               skull-stripped versions, time blocks (Do_cut), regional values and TACs.
srtm (step 8)  dynamic run -> reference-region models (roi_signal.dft): regional and voxel-wise BPnd,
               Logan DVR.
Images are not resampled here (PET space); step 9 exports them to the anatomical / template spaces.
"""
import os
import shutil


opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet import paramImg, cut_image, roi_signal


def suv(file, paths, timing, pet_info, values, isotope, masks, atlases, label_tables, levels, SUV, bratio,
        SUV_window, Do_cut, duration, diary_name):
    """Returns (3D parametric images written, activity image used, ROI tables)."""
    post, roi_dir, pict_dir = paths['orig_post'], paths['orig_roi'], paths['orig_pict']
    os.makedirs(post, exist_ok=True)
    base = pu.base_name(file)
    dynamic = pu.is_dynamic(file)
    brain = masks.get('brain')
    produced, tables = [], {}

    # 1. activity image (decay-corrected to the injection time when the data are not decay-corrected)
    src = file
    if dynamic and pet_info.get('ImageDecayCorrected') is False:
        if isotope in pu.HALF_LIFE_S and pet_info.get('InjectionStart') is not None:
            img, d = pu.load(file)
            src = pu.out_name(post, base, desc='decaycorr')
            pu.save_like(pu.decay_correct_frames(d, timing, isotope, float(pet_info['InjectionStart'])), img, src)
            pu.write_json(src, file, f'Decay correction to the injection time ({isotope}).')
        else:
            run_cmd.msg('WARNING: data not decay-corrected and isotope / injection time unknown: the SUV is biased',
                        diary_name, 'WARNING')
    if dynamic:
        act = pu.out_name(post, base, desc='activity-mean')
        pu.write_mean(src, act, timing, SUV_window,
                      f'Duration-weighted mean of the frames of the window {SUV_window} min (None = all).')
        tag = '-mean'
    else:
        act = src
        tag = ''

    # 2. SUV, SUVglu
    if SUV:
        f, report = pu.suv_scale(pet_info, values, isotope)
        if f is None:
            run_cmd.msg('WARNING: ' + report, diary_name, 'WARNING')
        else:
            for k, v in report.items():
                if k.startswith('warning') or k == 'assumption':
                    run_cmd.msg('WARNING: SUV: ' + str(v), diary_name, 'WARNING')
            out = paramImg.suv(act, f, pu.out_name(post, base, desc='SUV' + tag), report)
            produced.append(out)
            if dynamic:
                paramImg.suv(src, f, pu.out_name(post, base, desc='SUV'), report)     # dynamic SUV (not exported)
            if values.get('glucose_mmol_L'):
                produced.append(paramImg.suvglu(out, values['glucose_mmol_L'], pu.out_name(post, base, desc='SUVglu' + tag)))

    # 3. ratio to the brain mean, time blocks
    if bratio:
        if brain:
            out, m = paramImg.bratio(act, brain, pu.out_name(post, base, desc='bratio' + tag))
            run_cmd.msg(f'INFO: mean over the brain mask: {m:.5g}', diary_name, 'OKGREEN')
            produced.append(out)
        else:
            run_cmd.msg('WARNING: no brain mask in the PET space: bratio not computed', diary_name, 'WARNING')
    if Do_cut and dynamic:
        produced += cut_image.doit(src, timing, duration, brain, post, base, 'bratio' if brain else 'activity')

    # 4. skull-stripped versions
    if brain:
        for img in list(produced):
            d = pu.strip_ext(img).split('_desc-')[1].rsplit('_', 1)[0]
            produced.append(paramImg.skullstrip(img, brain, pu.out_name(post, base, desc=d + '-SS')))

    # 5. regional values and TACs
    if atlases:
        for img in produced:
            if img.endswith('-SS_pet.nii.gz'):
                continue
            d = pu.strip_ext(img).split('_desc-')[1].rsplit('_', 1)[0]
            tables.update({(k, d): v for k, v in roi_signal.read(img, atlases, levels, label_tables, roi_dir, None, d,
                                                                 diary_name).items()})
        if dynamic:
            marks = [t[0] / 60.0 for t in cut_image.setTime(timing, duration)] if Do_cut else None
            tables.update({(k, 'Kinetic'): v for k, v in roi_signal.read(src, atlases, levels, label_tables, roi_dir,
                                                                         pict_dir, 'Kinetic', diary_name, timing,
                                                                         marks).items()})
    return produced, src, tables


def srtm(file, paths, timing, sif, pet_info, isotope, masks, atlases, label_tables, levels, ref, DV_time,
         diary_name, sing_turku, end_min=None):
    """Reference-region models of a dynamic run (skull-stripped data for the voxel-wise models)."""
    if not pu.is_dynamic(file):
        run_cmd.msg('INFO: ' + opb(file) + ' is static: no kinetic model', diary_name, 'OKGREEN')
        return []
    if not sif or not opi(sif):
        raise FileNotFoundError(f'SIF of {file} not found ({sif}): step 1 must create it (srtm=True)')
    post = paths['orig_post']
    base = pu.base_name(file)
    img = file
    if masks.get('brain'):
        from modalities.pet import paramImg
        img = paramImg.skullstrip(file, masks['brain'], pu.out_name(post, base, desc='dynamic-SS'))
        shutil.copyfile(sif, img[:-7] + '.sif')
    unit = pet_info.get('Units') or 'kBq/mL'
    decay_corrected = pet_info.get('ImageDecayCorrected', True) is not False
    return roi_signal.dft(img, atlases, ref, masks, label_tables, timing, img[:-7] + '.sif' if img != file else sif,
                          unit, decay_corrected, isotope, DV_time, paths['orig_roi'], paths['orig_pict'], post,
                          diary_name, sing_turku, levels, end_min)
