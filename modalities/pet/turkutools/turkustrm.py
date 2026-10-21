"""
SRTM with tpcclib: regional fit with basis functions (bfmsrtm) and parametric image (imgsrtm).
(File name kept from the PET engineer; the model is SRTM, Lammertsma & Hume 1996.)
"""
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet import dftFile
from modalities.pet.turkutools._nifti import TurkuWorkspace


def bfm(ttac_dft, ref_dft, sif, res_file, fit_dft, isotope_if_not_corrected, pict_dir, ref_name,
        diary_name, sing_turku, end_min=None, n_plot=16):
    """
    dftweigh ttac sif, then bfmsrtm [-fit=file] [-i=isotope] [-end=min] TTAC RTAC results
    http://www.turkupetcentre.net/petanalysis/tpcclib/doc/bfmsrtm.html
    -i= is given only when the data are NOT decay-corrected (bfmsrtm then handles the decay; giving it for
    decay-corrected data, as the engineer's version did for every image, is wrong).
    Figures: Before_fit (TACs of the n_plot largest regions + reference) and one fit per region.
    """
    nl = '# Regional BPnd: SRTM with basis functions (bfmsrtm)'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    for f in (ttac_dft, ref_dft, sif):
        if not opi(f):
            raise FileNotFoundError(f"Input file not found: {f}")
    run_cmd.run(sing_turku + 'dftweigh ' + ttac_dft + ' ' + sif, diary_name)
    opts = '-fit=' + fit_dft
    if isotope_if_not_corrected:
        opts += ' -i=' + isotope_if_not_corrected
    if end_min:
        opts += f' -end={end_min:g}'
    run_cmd.run(sing_turku + 'bfmsrtm ' + opts + ' ' + ttac_dft + ' ' + ref_dft + ' ' + res_file, diary_name)

    if pict_dir:
        os.makedirs(pict_dir, exist_ok=True)
        t_start, t_end, values, names, vols = dftFile.read(ttac_dft)
        _, _, ref_val, _, _ = dftFile.read(ref_dft)
        x = (t_start + t_end) / 2.0 / 60.0
        order = np.argsort(vols)[::-1][:n_plot]
        fig, ax = plt.subplots(figsize=(8, 5))
        for i in order:
            ax.plot(x, values[:, i], 'o-', lw=1, ms=3, label=names[i])
        ax.plot(x, ref_val[:, 0], 'ko-', lw=2.5, ms=3, label='reference: ' + ref_name)
        ax.set_xlabel('time (min)')
        ax.set_ylabel('activity')
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6)
        fig.tight_layout()
        for ext in ('png', 'svg'):
            fig.savefig(opj(pict_dir, 'Before_fit.' + ext))
        plt.close(fig)
        if opi(fit_dft):
            fs, fe, fit_val, fit_names, _ = dftFile.read(fit_dft)
            xf = (fs + fe) / 2.0 / 60.0
            for i, name in enumerate(names):
                j = fit_names.index(name) if name in fit_names else i
                if j >= fit_val.shape[1]:
                    continue
                fig, ax = plt.subplots(figsize=(4, 3))
                ax.plot(x, values[:, i], 'o', ms=3, label='data')
                ax.plot(xf, fit_val[:, j], '-', lw=1.5, label='SRTM fit')
                ax.set_title(name, fontsize=8)
                ax.set_xlabel('time (min)')
                ax.legend(fontsize=6)
                fig.tight_layout()
                fig.savefig(opj(pict_dir, name + '_fit.png'))
                plt.close(fig)
    return res_file


def img(img, ref_dft, sif, out, diary_name, sing_turku, end_min=None):
    """
    imgsrtm [-end=min] imgfile rtacfile bpfile: parametric BPnd image, SRTM (linearised / NNLS solution).
    http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgsrtm.html
    """
    nl = '# Parametric image of BPnd: SRTM (imgsrtm)'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    for f in (img, ref_dft):
        if not opi(f):
            raise FileNotFoundError(f"Input file not found: {f}")
    with TurkuWorkspace(opd(out)) as ws:
        nii = ws.nii(img, sif)
        bp = ws.path('srtm.nii')
        opts = (f'-end={end_min:g} ' if end_min else '')
        run_cmd.run(sing_turku + 'imgsrtm ' + opts + nii + ' ' + ref_dft + ' ' + bp, diary_name)
        ws.back(bp, img, out)
    pu.write_json(out, [img, ref_dft], 'BPnd, SRTM (tpcclib imgsrtm; Lammertsma & Hume 1996).')
    return out
