"""imgbfbp (tpcclib): parametric BPnd image, SRTM fitted with basis functions."""
import os

opi = os.path.isfile
opd = os.path.dirname

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet.turkutools._nifti import TurkuWorkspace


def run(img, ref_dft, sif, bp_out, err_out, diary_name, sing_turku, end_min=None):
    """
    imgbfbp [-err=file] [-end=min] imgfile rtacfile bpfile [sif]
    http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgbfbp.html
    The image and the reference TAC must be decay-corrected (checked by the caller).
    err_out: 1 / 2 where theta3 reached its lower / upper limit (fit not reliable).
    """
    nl = ('# Parametric image of BPnd: simplified reference tissue model (SRTM, Lammertsma & Hume 1996) '
          'with basis functions (Gunn et al. 1997), imgbfbp')
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    for f in (img, ref_dft, sif):
        if not opi(f):
            raise FileNotFoundError(f"Input file not found: {f}")
    with TurkuWorkspace(opd(bp_out)) as ws:
        nii = ws.nii(img, sif)
        bp, err = ws.path('bp.nii'), ws.path('bperr.nii')
        opts = '-err=' + err + (f' -end={end_min:g}' if end_min else '')
        run_cmd.run(sing_turku + 'imgbfbp ' + opts + ' ' + nii + ' ' + ref_dft + ' ' + bp + ' ' + nii[:-4] + '.sif',
                    diary_name)
        ws.back(bp, img, bp_out)
        ws.back(err, img, err_out)
    desc = ('BPnd, SRTM with basis functions (tpcclib imgbfbp; Lammertsma & Hume 1996, Gunn et al. 1997).')
    pu.write_json(bp_out, [img, ref_dft, sif], desc)
    pu.write_json(err_out, [img, ref_dft, sif], 'theta3 at its lower (1) or upper (2) limit (imgbfbp -err).')
    return bp_out
