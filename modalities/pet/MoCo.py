"""
Motion correction of the PET images (steps 4 and 5 of the PET pipeline).

dynamic (step 4, within a dynamic run)
    Reference: duration-weighted mean of the frames from `nb_void` on (the first frames, low in counts
    and with another distribution of the tracer, are not registered; they receive the transform of the
    first registered frame, the nearest in time); two passes, the second one registering to the mean of
    the frames corrected by the first (a mean of misaligned frames is a blurred reference). Every frame is registered rigidly to the reference on
    images cropped to the head box (crop drawn on the study template, or on the session), and the
    transform is applied to the full frame (linear interpolation, which preserves the local mean
    concentration). Motion parameters and mean displacement are written for the quality control.
static (step 5, between the runs of one session)
    Every run (its mean) is registered rigidly to the mean of the runs and resampled on the grid of the
    first run, so that all the runs of the session share one space. With one run, the images are copied.

Outputs (orig/preprocessing): <base>_desc-MoCo-within_pet.nii.gz, <base>_desc-MoCo-between_pet.nii.gz and
their means (-mean), the session mean <session>_desc-MoCo-between-mean_pet.nii.gz; matrices and motion
tables in preprocessing/forMC.
"""
import os
import shutil

import numpy as np
import pandas as pd
import ants
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

opj = os.path.join
opd = os.path.dirname
opb = os.path.basename
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet import sifFile

petinterp = 'linear'


def _register(fixed, moving, prefix, transform, metric):
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform=transform, aff_metric=metric,
                            outprefix=prefix)
    return reg['fwdtransforms'][-1]


def motion_report(mats, timing, radius_mm, tsv, png, title):
    rows = []
    for i, m in enumerate(mats):
        t, r, d = pu.rigid_parameters(m, radius_mm)
        rows.append(dict(frame=i, start_s=timing['start'][i] if timing else np.nan,
                         tx_mm=t[0], ty_mm=t[1], tz_mm=t[2], rx_deg=r[0], ry_deg=r[1], rz_deg=r[2],
                         mean_displacement_mm=d))
    df = pd.DataFrame(rows)
    df.to_csv(tsv, sep='\t', index=False, float_format='%.4f')
    x = df['start_s'] / 60.0 if timing else df['frame']
    fig, ax = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
    for c in ('tx_mm', 'ty_mm', 'tz_mm'):
        ax[0].plot(x, df[c], 'o-', ms=3, label=c[:2])
    for c in ('rx_deg', 'ry_deg', 'rz_deg'):
        ax[1].plot(x, df[c], 'o-', ms=3, label=c[:2])
    ax[2].plot(x, df['mean_displacement_mm'], 'ko-', ms=3)
    ax[0].set_ylabel('translation (mm)')
    ax[1].set_ylabel('rotation (deg)')
    ax[2].set_ylabel(f'mean displ. (mm)\n(sphere r={radius_mm:.0f} mm)')
    ax[2].set_xlabel('frame start (min)' if timing else 'frame')
    ax[0].legend(fontsize=7, ncol=3)
    ax[1].legend(fontsize=7, ncol=3)
    ax[0].set_title(title, fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=120)
    plt.close(fig)
    return df


def dynamic(img, out_dir, formc_dir, timing, nb_void, crop_mask, Need_check, diary_name, sing_fs,
            transform='Rigid', metric='mattes', sif='', passes=2):
    """Within-run motion correction of one dynamic image (see module docstring). Returns the corrected image."""
    nl = 'Motion correction within the run ' + opb(img)
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(formc_dir, exist_ok=True)

    base = pu.base_name(img)
    mc4D = pu.out_name(out_dir, base, desc='MoCo-within')
    mc3D = pu.out_name(out_dir, base, desc='MoCo-within-mean')
    ref3D = pu.out_name(formc_dir, base, desc='MoCo-ref')

    full4d = ants.image_read(img)
    if full4d.dimension != 4:
        raise ValueError(f'{img} is not a dynamic image')
    frames_img = ants.ndimage_to_list(full4d)
    nfr = len(frames_img)
    if nb_void >= nfr:
        raise ValueError(f'nb_void = {nb_void} but the run has only {nfr} frames')
    frames = np.arange(nb_void, nfr)
    w = timing['duration'][frames] if timing is not None else np.ones(len(frames))
    w = w / w.sum()
    ref_full = frames_img[frames[0]] * float(w[0])
    for k, i in enumerate(frames[1:], start=1):
        ref_full = ref_full + frames_img[i] * float(w[k])
    ants.image_write(ref_full, ref3D, ri=False)
    pu.write_json(ref3D, img, f'Reference of the motion correction: duration-weighted mean of frames {nb_void}-{nfr - 1}.')

    if crop_mask:
        cmask = ants.resample_image_to_target(ants.image_read(crop_mask), ref_full, interp_type='nearestNeighbor')
        ref_crop = ants.crop_image(ref_full, label_image=cmask, label=1)
        ants.image_write(ref_crop, pu.out_name(formc_dir, base, desc='crop-mean'), ri=False)
    else:
        cmask = None
        ref_crop = ref_full

    mats, corrected = [None] * nfr, [None] * nfr
    for p in range(passes):
        # pass 1: reference = mean of the uncorrected frames; pass 2: mean of the frames corrected by pass 1
        # (sharper reference; two-pass registration to the mean, as SPM realign / AFNI 3dvolreg -base mean)
        for i in frames:
            fr = frames_img[i]
            fr_crop = ants.crop_image(fr, label_image=cmask, label=1) if cmask is not None else fr
            mats[i] = _register(ref_crop, fr_crop, opj(formc_dir, base + '_mc' + f'{i:03d}' + '_'), transform, metric)
            corrected[i] = ants.apply_transforms(fixed=ref_full, moving=fr, transformlist=[mats[i]],
                                                 interpolator=petinterp)
        if p < passes - 1:
            ref_full = corrected[frames[0]] * float(w[0])
            for k, i in enumerate(frames[1:], start=1):
                ref_full = ref_full + corrected[i] * float(w[k])
            ants.image_write(ref_full, ref3D, ri=False)
            ref_crop = ants.crop_image(ref_full, label_image=cmask, label=1) if cmask is not None else ref_full
    for i in range(nb_void):
        # nearest registered frame in time (the engineer left these frames uncorrected)
        dst = opj(formc_dir, base + '_mc' + f'{i:03d}' + '_0GenericAffine.mat')
        shutil.copyfile(mats[nb_void], dst)
        mats[i] = dst
        corrected[i] = ants.apply_transforms(fixed=ref_full, moving=frames_img[i], transformlist=[dst],
                                             interpolator=petinterp)
    ants.image_write(ants.list_to_ndimage(full4d, corrected), mc4D, ri=False)
    pu.write_json(mc4D, [img, ref3D, crop_mask or 'no crop'],
                  f'Motion correction within the run ({transform}, {metric}, ANTsPy); frames < {nb_void} '
                  f'received the transform of frame {nb_void}; linear interpolation.',
                  Matrices=[opb(m) for m in mats])
    sifFile.copy(sif, mc4D)
    pu.write_mean(mc4D, mc3D, timing, None, 'Duration-weighted mean of the motion-corrected frames.')

    radius = pu.head_radius_mm(ref_full.numpy(), ref_full.spacing)
    motion_report(mats, timing, radius, opj(formc_dir, base + '_desc-MoCo-within_motion.tsv'),
                  opj(formc_dir, base + '_desc-MoCo-within_motion.png'), base)

    if Need_check:
        run_cmd.do(sing_fs + 'freeview -v ' + mc4D + ' -subtitle "Check-residual-head-motion!"', diary_name)
    return mc4D


def _grid4d(g3, like4d):
    """Empty 4D image with the spatial grid of g3 and the frames of like4d (template of list_to_ndimage)."""
    n = like4d.shape[3]
    d = np.eye(4)
    d[:3, :3] = g3.direction
    return ants.from_numpy(np.zeros(tuple(g3.shape) + (n,), np.float32), origin=tuple(g3.origin) + (0.0,),
                           spacing=tuple(g3.spacing) + (like4d.spacing[3],), direction=d)


def static(list_img, list_mean, list_timing, list_sif, session, out_dir, formc_dir, crop_mask, diary_name,
           transform='Rigid', metric='mattes'):
    """
    Between-run motion correction of the runs of one session (images after step 4 for dynamic runs, after
    step 1 for static ones). list_mean: mean image of every run. Returns (corrected images, session mean).
    """
    nl = 'Between-runs motion correction (' + str(len(list_img)) + ' run(s))'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(formc_dir, exist_ok=True)
    session_mean = pu.out_name(out_dir, session, desc='MoCo-between-mean')
    outs = []

    if len(list_img) == 1:
        base = pu.base_name(list_img[0])
        out = pu.out_name(out_dir, base, desc='MoCo-between')
        shutil.copyfile(list_img[0], out)
        pu.write_json(out, list_img[0], 'Single run: no between-runs correction (copy).')
        sifFile.copy(list_sif[0], out)
        shutil.copyfile(list_mean[0], session_mean)
        if pu.is_dynamic(out):
            shutil.copyfile(list_mean[0], pu.out_name(out_dir, base, desc='MoCo-between-mean'))
        pu.write_json(session_mean, list_mean[0], 'Mean image of the session (single run).')
        return [out], session_mean

    first = ants.image_read(list_mean[0])
    means = [ants.resample_image_to_target(ants.image_read(m), first, interp_type='linear') for m in list_mean]
    ref = means[0] * 0
    for m in means:
        ref = ref + m / len(means)
    ref_file = pu.out_name(formc_dir, session, desc='between-ref')
    ants.image_write(ref, ref_file, ri=False)
    cmask = ants.resample_image_to_target(ants.image_read(crop_mask), first, interp_type='nearestNeighbor') \
        if crop_mask else None
    ref_crop = ants.crop_image(ref, label_image=cmask, label=1) if cmask is not None else ref

    corrected_means = []
    for img, mean, timing, sif in zip(list_img, list_mean, list_timing, list_sif):
        base = pu.base_name(img)
        out = pu.out_name(out_dir, base, desc='MoCo-between')
        mv = ants.image_read(mean)
        mv_crop = ants.crop_image(mv, label_image=ants.resample_image_to_target(cmask, mv, interp_type='nearestNeighbor'),
                                  label=1) if cmask is not None else mv
        mat = _register(ref_crop, mv_crop, opj(formc_dir, base + '_between_'), transform, metric)
        side = pu.read_json(pu.sidecar(img))
        src = side.get('Sources', [None])[0] if isinstance(side.get('Sources'), list) else None
        if side.get('Matrices') and src and opi(src):
            # dynamic run corrected at step 4: the frames of the reoriented run are resampled once, with the
            # between-runs and the per-frame transforms composed (no second interpolation)
            raw4d = ants.image_read(src)
            frames_raw = ants.ndimage_to_list(raw4d)
            if len(frames_raw) != len(side['Matrices']):
                raise ValueError(f'{img}: {len(side["Matrices"])} frame matrices for {len(frames_raw)} frames')
            moved_frames = [ants.apply_transforms(fixed=first, moving=fr, interpolator=petinterp,
                                                  transformlist=[mat, opj(formc_dir, m)])
                            for fr, m in zip(frames_raw, side['Matrices'])]
            moved = ants.list_to_ndimage(_grid4d(first, raw4d), moved_frames)
            how = 'per-frame and between-runs transforms composed, one resampling of the reoriented frames'
        else:
            full = ants.image_read(img)
            moved = ants.apply_transforms(fixed=first, moving=full, transformlist=[mat], interpolator=petinterp,
                                          imagetype=3 if full.dimension == 4 else 0)
            how = 'one resampling'
        ants.image_write(moved, out, ri=False)
        pu.write_json(out, [img, ref_file], f'Between-runs motion correction ({transform}, {metric}, ANTsPy) on the '
                      f'grid of the first run; linear interpolation; {how}.', Matrix=opb(mat))
        sifFile.copy(sif, out)
        mean_out = pu.out_name(out_dir, base, desc='MoCo-between-mean')
        if pu.is_dynamic(out):
            pu.write_mean(out, mean_out, timing, None, 'Duration-weighted mean after the between-runs correction.')
        else:
            shutil.copyfile(out, mean_out)
        corrected_means.append(mean_out)
        outs.append(out)

    acc = None
    for m in corrected_means:
        a = ants.image_read(m)
        acc = a / len(corrected_means) if acc is None else acc + a / len(corrected_means)
    ants.image_write(acc, session_mean, ri=False)
    pu.write_json(session_mean, corrected_means, 'Mean image of the session (mean of the run means after correction).')
    return outs, session_mean
