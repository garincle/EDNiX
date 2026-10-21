"""imgdv (tpcclib): parametric DVR image, Logan plot with a reference region."""
import os

opi = os.path.isfile
opd = os.path.dirname

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet.turkutools._nifti import TurkuWorkspace


def run(img, ref_dft, sif, DV_time, dvr_out, n_out, diary_name, sing_turku, end_min=None):
    """
    imgdv [-n=file] [-end=min] input image starttime resultimage
    http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgdv.html
    With a reference-region TAC as input the result is the DVR (Logan et al. 1996); DV_time is the start
    of the linear part of the plot (t*, minutes). n_out: number of points of the fitted line.
    """
    nl = ('# Parametric image of DVR: multiple-time graphical analysis (Logan plot, reference region; '
          'Logan et al. 1996), imgdv, t* = ' + str(DV_time) + ' min')
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    for f in (img, ref_dft):
        if not opi(f):
            raise FileNotFoundError(f"Input file not found: {f}")
    with TurkuWorkspace(opd(dvr_out)) as ws:
        nii = ws.nii(img, sif)
        dvr, npts = ws.path('dvr.nii'), ws.path('npts.nii')
        opts = '-n=' + npts + (f' -end={end_min:g}' if end_min else '')
        run_cmd.run(sing_turku + 'imgdv ' + opts + ' ' + ref_dft + ' ' + nii + ' ' + str(DV_time) + ' ' + dvr,
                    diary_name)
        ws.back(dvr, img, dvr_out)
        ws.back(npts, img, n_out)
    pu.write_json(dvr_out, [img, ref_dft], 'DVR, Logan plot with a reference region (tpcclib imgdv; Logan et al. 1996).',
                  StartTime_min=DV_time)
    pu.write_json(n_out, [img, ref_dft], 'Number of points of the Logan line (imgdv -n).')
    return dvr_out
