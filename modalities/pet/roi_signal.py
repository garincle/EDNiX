"""
Regional analysis in the PET space (steps 7 and 8): regional values of the parametric images,
time-activity curves (TACs) of the dynamic images, and the reference-region kinetic models (Turku tools).

Files (orig/postprocessing/ROI_analysis), for every atlas and level:
    <base>_seg-<atlas>-lvl<k>_vol_ROIs.tsv            label, region, voxels, volume (mm3)
    <base>_seg-<atlas>-lvl<k>_desc-<desc>_ROIs.tsv    label, region, voxels, mean, SD, median (3D image)
    <base>_seg-<atlas>-lvl<k>_desc-Kinetic_ROIs.tsv   frame times + one TAC per region (4D image)
    .dft versions of the TACs and of the reference TAC (Turku format), BPnd results (bfmsrtm)
    pictures/<base>_seg-..._time_activity_curve.png/.svg, ..._mean-time_activity_curve.png/.svg
The PET engineer's version called 3dROIstats and rebuilt the order of the labels; values are now attached
to their label (Tools/roiStats.py). TAC tables carry the frame times, which the DFT files need.
"""
import os

import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd, roiStats
from modalities.pet import petutils as pu
from modalities.pet import dftFile
from modalities.pet.turkutools import turkubfbp, turkustrm, turkudv

dftname = 'Kinetic'


def atlas_levels(atlas_file):
    shape = nib.load(atlas_file).shape
    return shape[3] if len(shape) > 3 else 1


def labels_volume(atlas_file, level):
    return roiStats._labels_volume(atlas_file, level)


def prefix(roi_dir, base, atlas, level):
    return opj(roi_dir, f'{base}_seg-{atlas}-lvl{level}')


def read(img, atlases, levels, label_tables, roi_dir, pict_dir, desc, diary_name, timing=None, cut_marks=None,
         n_plot=16):
    """
    Regional values of `img` (PET space) for every atlas (dict name -> atlas file in PET space) and level
    (levels: dict name -> list of levels, or 'all'). Returns {(atlas, level): table file}.
    """
    run_cmd.msg('extract the regional values of ' + opb(img), diary_name, 'OKGREEN')
    os.makedirs(roi_dir, exist_ok=True)
    nii, data = pu.load(img)
    vox_mm3 = float(np.prod(nii.header.get_zooms()[:3]))
    base = pu.base_name(img)
    out = {}
    for atlas, afile in atlases.items():
        if nib.load(afile).shape[:3] != data.shape[:3]:
            raise ValueError(f'{afile} is not on the grid of {img}')
        lv = levels.get(atlas, 'all') if isinstance(levels, dict) else levels
        lv = list(range(atlas_levels(afile))) if lv == 'all' else list(lv)
        for level in lv:
            lab = labels_volume(afile, level)
            pre = prefix(roi_dir, base, atlas, level)
            vol_file = pre + '_vol_ROIs.tsv'
            if not opi(vol_file):
                roiStats.regional_table(lab, label_tables.get(atlas), vox_mm3).to_csv(vol_file, sep='\t', index=False)
            if data.ndim == 3:
                tab = roiStats.stats(data, lab, label_tables.get(atlas), vox_mm3)
                f = pre + f'_desc-{desc}_ROIs.tsv'
                tab.to_csv(f, sep='\t', index=False, float_format='%.6g')
            else:
                tac = roiStats.tacs(data, lab, label_tables.get(atlas), vox_mm3)
                f = pre + f'_desc-{dftname}_ROIs.tsv'
                tab = tac.copy()
                tab.columns = [str(c) for c in tab.columns]
                if timing is not None:
                    tab.insert(0, 'frame_end_s', timing['end'])
                    tab.insert(0, 'frame_start_s', timing['start'])
                tab.to_csv(f, sep='\t', index=False, float_format='%.6g')
                with open(f.replace('.tsv', '_names.tsv'), 'w') as fh:
                    fh.write('label\tregion\tvolume_mm3\n')
                    for l, n, v in zip(tac.columns, tac.attrs['names'], tac.attrs['volumes_mm3']):
                        fh.write(f'{l}\t{n}\t{v:.6g}\n')
                if pict_dir:
                    plot_tacs(tac, timing, opj(pict_dir, opb(pre)), cut_marks, n_plot)
            out[(atlas, level)] = f
    return out


def plot_tacs(tac, timing, pre, cut_marks=None, n_plot=16):
    os.makedirs(opd(pre), exist_ok=True)
    x = timing['mid'] / 60.0 if timing is not None else np.arange(tac.shape[0])
    vols = np.asarray(tac.attrs['volumes_mm3'])
    order = np.argsort(vols)[::-1][:n_plot]
    fig, ax = plt.subplots(figsize=(8, 5))
    for i in order:
        ax.plot(x, tac.iloc[:, i], 'o-', lw=1, ms=3, label=tac.attrs['names'][i])
    ax.set_xlabel('time (min)')
    ax.set_ylabel('activity concentration')
    ax.set_title(f'{len(order)} largest regions of {tac.shape[1]}', fontsize=8)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6)
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(pre + '_time_activity_curve.' + ext)
    plt.close(fig)
    M, S = tac.mean(axis=1), tac.std(axis=1)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(x, M, yerr=S, fmt='o-', lw=2, color='red', ms=3)
    for t in (cut_marks or []):
        ax.axvline(t, color='k', lw=0.8)
    ax.set_xlabel('time (min)')
    ax.set_ylabel('mean +- SD across regions')
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(pre + '_mean-time_activity_curve.' + ext)
    plt.close(fig)


def dft(img, atlases, ref, masks, label_tables, timing, sif, unit, decay_corrected, isotope, DV_time,
        roi_dir, pict_dir, out_dir, diary_name, sing_turku, levels='all', end_min=None, study='EDNiX'):
    """
    Reference-region models of one dynamic image (step 8):
      regional BPnd (bfmsrtm) for every atlas level, parametric BPnd (imgbfbp, imgsrtm) and DVR (imgdv,
      if DV_time > 0) images. ref: resolved reference region (atlases.getlabels.resolve).
    The voxel-wise tools need decay-corrected data: they are skipped otherwise (bfmsrtm then gets -i=).
    Returns the parametric images written.
    """
    nl = 'Reference-region kinetic models of ' + opb(img)
    run_cmd.msg(nl, diary_name, 'HEADER')
    os.makedirs(roi_dir, exist_ok=True)
    nii, data = pu.load(img)
    if data.ndim != 4:
        raise ValueError(f'{img} is not dynamic')
    vox_mm3 = float(np.prod(nii.header.get_zooms()[:3]))
    base = pu.base_name(img)

    # reference TAC
    ref_dft = opj(roi_dir, f"{base}_{ref['name']}_{dftname}.dft")
    if ref['type'] == 'atlas':
        lab = labels_volume(atlases[ref['atlas']], int(ref['level']))
        _, nvox = dftFile.ref(ref_dft, data, lab, ref['labels'], vox_mm3, timing, unit, ref['name'], study)
    else:
        name = ref['name']
        if name == 'white+csf':
            msk = pu.mask_bool(masks['white'], img) | pu.mask_bool(masks['csf'], img)   # union (the engineer's a*b was empty)
        else:
            msk = pu.mask_bool(masks[name], img)
        _, nvox = dftFile.mask(ref_dft, data, msk, vox_mm3, timing, unit, name, study)
    run_cmd.msg(f"INFO: reference region {ref['name']}: {nvox} voxels ({nvox * vox_mm3:.1f} mm3)", diary_name, 'OKGREEN')

    # regional BPnd
    for atlas, afile in atlases.items():
        lv = levels.get(atlas, 'all') if isinstance(levels, dict) else levels
        lv = list(range(atlas_levels(afile))) if lv == 'all' else list(lv)
        for level in lv:
            labv = labels_volume(afile, level)
            tac = roiStats.tacs(data, labv, label_tables.get(atlas), vox_mm3)
            pre = prefix(roi_dir, base, atlas, level)
            ttac = dftFile.roi(pre + f'_desc-{dftname}_ROIs.dft', tac, timing, unit, study)
            turkustrm.bfm(ttac, ref_dft, sif, pre + f'_desc-{dftname}_ROIs_bp.res',
                          pre + f'_desc-{dftname}_ROIs_fit.dft', None if decay_corrected else isotope,
                          opj(pict_dir, opb(pre)) if pict_dir else None, ref['name'], diary_name, sing_turku,
                          end_min)

    out = []
    if not decay_corrected:
        run_cmd.msg('WARNING: the image is not decay-corrected: the parametric images (imgbfbp, imgsrtm, imgdv) '
                    'need decay-corrected data and are not computed', diary_name, 'WARNING')
        return out
    bp = pu.out_name(out_dir, base, desc='BP')
    out.append(turkubfbp.run(img, ref_dft, sif, bp, pu.out_name(out_dir, base, desc='BP-err'), diary_name,
                             sing_turku, end_min))
    out.append(turkustrm.img(img, ref_dft, sif, pu.out_name(out_dir, base, desc='SRTM'), diary_name, sing_turku,
                             end_min))
    if DV_time and float(DV_time) > 0:
        out.append(turkudv.run(img, ref_dft, sif, DV_time, pu.out_name(out_dir, base, desc='DVR'),
                               pu.out_name(out_dir, base, desc='DVR-npoints'), diary_name, sing_turku, end_min))
    return out
