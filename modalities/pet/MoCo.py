"""
Motion correction of the PET images (steps 4 and 5 of the PET pipeline).

dynamic (step 4, within a dynamic run)
    Reference: duration-weighted mean of the frames from `nb_void` on (the first frames, low in counts
    and with another distribution of the tracer, are not registered; they receive the transform of the
    first registered frame, the nearest in time); two passes, the second one registering to the mean of
    the frames corrected by the first (a mean of misaligned frames is a blurred reference); metric computed in a
    head mask of the reference (MOCO_REG, head_mask). Every frame is registered rigidly to the reference on
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

# Registration settings of the motion correction. The metric is computed in a head mask of the reference: without
# it the background (most of a PET field of view) dominates the mutual-information histogram and the frames of a
# still phantom received 0.4-2 mm of spurious motion (pairwise error 1.5-2 mm, max 3-4 mm, against 0.2-0.3 mm with
# the mask). Three levels ending at full resolution with all voxels sampled (ANTsPy's default Rigid ends with 10
# iterations at full resolution, sampling 20 % of the voxels).
MOCO_REG = dict(aff_iterations=(1000, 500, 250), aff_shrink_factors=(4, 2, 1), aff_smoothing_sigmas=(2, 1, 0),
                aff_random_sampling_rate=1.0)

def crop_mean(mean_file, crop_mask, out_file):
    """Cropped version of the session mean (same world coords, smaller FOV), for the registration to the
    anatomy / template. The crop box of step 4-5 is only used in the MoCo metric; the MoCo outputs are
    full FOV, so ANTs' center-of-mass initialisation (bladder, injection site, ...) is biased when the
    full image is passed to the registration. The transforms produced on the cropped image are
    world-to-world and remain valid for the full session mean downstream."""
    ref = ants.image_read(mean_file)
    msk = ants.resample_image_to_target(ants.image_read(crop_mask), ref, interp_type='nearestNeighbor')
    ants.image_write(ants.crop_image(ref, label_image=msk, label=1), out_file, ri=False)
    pu.write_json(out_file, [mean_file, crop_mask],
                  'Session mean cropped to the MoCo crop box (same world coords, smaller FOV), for the registration '
                  'to the anatomy / template.')
    return out_file

def head_mask(img, frac=0.1, dilate_vox=2):
    """Head of a PET image: voxels above frac x its 99th percentile, dilated (ANTsImage)."""
    a = img.numpy()
    pos = a[np.isfinite(a) & (a > 0)]
    thr = frac * float(np.percentile(pos, 99)) if pos.size else 0.0
    return ants.iMath(ants.threshold_image(img, thr, 1e30, 1, 0), 'MD', dilate_vox)


def _identity_matrix(path):
    """Writes an identity rigid transform (3D) in ITK .mat format and returns its path."""
    ants.write_transform(ants.new_ants_transform(dimension=3, transform_type='AffineTransform'), path)
    return path


def _has_signal(img):
    """True if the image has a positive, finite total mass (ANTs' centre-of-mass initialisation needs it)."""
    a = img.numpy()
    a = a[np.isfinite(a)]
    return bool(a.size) and float(np.sum(a[a > 0])) > 0.0


def _register(fixed, moving, prefix, transform, metric, fallback=None):
    """
    Rigid registration moving -> fixed, returns the path of the matrix. A frame without signal in the cropped
    box (all zeros) cannot be registered: ANTs' initialisation (image moments) aborts with 'Total Mass of the
    image was zero'. Such a frame receives `fallback` (matrix of the previous registered frame), or the identity.
    """
    for suffix in ('0GenericAffine.mat',):
        if opi(prefix + suffix):
            os.remove(prefix + suffix)
    if not _has_signal(moving):
        if fallback and opi(fallback):
            shutil.copyfile(fallback, prefix + '0GenericAffine.mat')
            return prefix + '0GenericAffine.mat'
        return _identity_matrix(prefix + '0GenericAffine.mat')
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform=transform, aff_metric=metric,
                            outprefix=prefix, mask=head_mask(fixed), **MOCO_REG)
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

    if not _has_signal(ref_full):
        raise ValueError(f'{img}: the reference (mean of frames {nb_void}-{nfr - 1}) has no signal; check the image and nb_void')

    if crop_mask:
        cmask = ants.resample_image_to_target(ants.image_read(crop_mask), ref_full, interp_type='nearestNeighbor')
        ref_crop = ants.crop_image(ref_full, label_image=cmask, label=1)
        ants.image_write(ref_crop, pu.out_name(formc_dir, base, desc='crop-mean'), ri=False)
        if not _has_signal(ref_crop):
            raise ValueError(f'{crop_mask}: the crop mask does not overlap the head of {img} (empty box); check its space')

    else:
        cmask = None
        ref_crop = ref_full

    mats, corrected = [None] * nfr, [None] * nfr
    for p in range(passes):
        # pass 1: reference = mean of the uncorrected frames; pass 2: mean of the frames corrected by pass 1
        # (sharper reference; two-pass registration to the mean, as SPM realign / AFNI 3dvolreg -base mean)
        prev = None
        for i in frames:
            fr = frames_img[i]
            fr_crop = ants.crop_image(fr, label_image=cmask, label=1) if cmask is not None else fr
            if not _has_signal(fr_crop):
                run_cmd.msg(f'WARNING: frame {i} has no signal in the registration box (all zeros): it receives the '
                            f'transform of the previous frame (identity if none)', diary_name, 'WARNING')
            mats[i] = _register(ref_crop, fr_crop, opj(formc_dir, base + '_mc' + f'{i:03d}' + '_'), transform, metric,
                                fallback=prev)
            prev = mats[i]
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
           transform='Rigid', metric='mattes', passes=2):
    """
    Between-run motion correction of the runs of one session (images after step 4 for dynamic runs, after
    step 1 for static ones). list_mean: mean image of every run. Reference: mean of the run means on the grid of
    the first run; second pass on the mean of the means corrected by the first (as within a run). Returns
    (corrected images, session mean).
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
    cmask = ants.resample_image_to_target(ants.image_read(crop_mask), first, interp_type='nearestNeighbor') \
        if crop_mask else None

    def crop(im):
        if cmask is None:
            return im
        return ants.crop_image(im, label_image=ants.resample_image_to_target(cmask, im, interp_type='nearestNeighbor'),
                               label=1)

    mats = []
    for p in range(passes):
        ref_crop = crop(ref)
        mats = [_register(ref_crop, crop(ants.image_read(m)), opj(formc_dir, pu.base_name(i) + '_between_'),
                          transform, metric) for i, m in zip(list_img, list_mean)]
        if p < passes - 1:
            moved = [ants.apply_transforms(fixed=first, moving=ants.image_read(m), transformlist=[t],
                                           interpolator=petinterp) for m, t in zip(list_mean, mats)]
            ref = moved[0] * 0
            for m in moved:
                ref = ref + m / len(moved)
    ants.image_write(ref, ref_file, ri=False)

    corrected_means = []
    for img, mean, timing, sif, mat in zip(list_img, list_mean, list_timing, list_sif, mats):
        base = pu.base_name(img)
        out = pu.out_name(out_dir, base, desc='MoCo-between')
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



def mean_rigid(mats, out):
    """
    Mean of rigid ANTs transforms (x' = A (x - c) + c + t): chordal L2 mean of the rotations (SVD projection of the
    mean matrix on SO(3)) and mean of the translations, written as one ANTs affine (centre 0) to `out`.
    """
    Rs, ts = [], []
    for m in mats:
        tr = ants.read_transform(m)
        p = np.asarray(tr.parameters, dtype=float)
        c = np.asarray(tr.fixed_parameters, dtype=float)[:3]
        A, t = p[:9].reshape(3, 3), p[9:12]
        Rs.append(A)
        ts.append(t + c - A @ c)
    U, _, Vt = np.linalg.svd(np.mean(Rs, axis=0))
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt
    tr = ants.create_ants_transform(transform_type='AffineTransform', precision='double', dimension=3, matrix=R,
                                    translation=np.mean(ts, axis=0), center=(0.0, 0.0, 0.0))
    ants.write_transform(tr, out)
    return out


SAME_SPACE_WINDOW_MIN = 10.0


def scanner_anchor(runs, corrected, timings, formc_dir, nb_void=0, window_min=None):
    """
    Chains between the session space (after the motion correction: average head position of the session) and the
    scanner space of the PET frames acquired with the anatomy (hybrid PET-MR, 1v1) or the CT (PET/CT): without it, a
    session whose head moved would be shifted from the anatomy / CT by the mean motion, with no registration to
    catch it.
    runs: reoriented runs (step 1); corrected: outputs of step 4 (MoCo-within, or the run itself when not corrected);
    window_min: [t0, t1] (min, PET TimeZero) of the anatomy / CT acquisition; None = the first 10 min of the first
    run from its first registered frame (CT acquired before the PET). The head position of the window is the mean
    (mean_rigid) of the motion transforms of its registered frames (frame nb_void when the window has none), the
    run is the one that overlaps the window most.
    Returns ({'to': ..., 'from': ...}, description): 'to' resamples a session image into the scanner space.
    """
    r = 0
    if window_min is not None:
        t0, t1 = [60.0 * float(x) for x in window_min]
        overlap = [0.0 if tm is None else float(np.clip(np.minimum(tm['end'], t1) - np.maximum(tm['start'], t0),
                                                         0, None).sum()) for tm in timings]
        if max(overlap) > 0:
            r = int(np.argmax(overlap))
    chain_from = []
    if len(runs) > 1:
        between = opj(formc_dir, pu.base_name(runs[r]) + '_between_0GenericAffine.mat')
        if not opi(between):
            raise FileNotFoundError(f'{between} not found: run step 5 first')
        chain_from.append([between, False])
    mats = pu.read_json(pu.sidecar(corrected[r])).get('Matrices') if 'MoCo-within' in opb(corrected[r]) else None
    desc = 'scanner position of run ' + opb(runs[r])
    if mats:
        tm = timings[r]
        if window_min is None:
            t0 = float(tm['start'][min(nb_void, len(mats) - 1)])
            t1 = t0 + 60.0 * SAME_SPACE_WINDOW_MIN
        frames = [k for k in range(nb_void, len(mats)) if t0 <= tm['mid'][k] <= t1]
        if not frames:
            centre = (t0 + t1) / 2.0
            frames = [min(range(nb_void, len(mats)), key=lambda k: abs(tm['mid'][k] - centre))]
        if len(frames) == 1:
            pos = opj(formc_dir, mats[frames[0]])
        else:
            pos = mean_rigid([opj(formc_dir, mats[k]) for k in frames],
                             opj(formc_dir, pu.base_name(runs[r]) + '_desc-scanner-position_0GenericAffine.mat'))
        chain_from.append([pos, False])
        desc += f', frames {frames} ({t0 / 60:.1f}-{t1 / 60:.1f} min)'
    return {'to': [[f, not inv] for f, inv in reversed(chain_from)], 'from': chain_from}, desc
