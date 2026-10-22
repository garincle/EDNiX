"""imgweigh (tpcclib): counts of each frame of a dynamic NIfTI image, written in its SIF (fit weights)."""
import os
import shutil

opi = os.path.isfile
opd = os.path.dirname

from Tools import run_cmd
from modalities.pet.turkutools._nifti import TurkuWorkspace


def weight(file, sif, diary_name, sing_turku):
    """
    imgweigh petfile sif: reads the frame times of `sif` (written by sifFile.write) and rewrites it with
    the counts estimated from the image. http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgweigh.html
    """
    for f in (file, sif):
        if not opi(f):
            raise FileNotFoundError(f"Input file not found: {f}")
    with TurkuWorkspace(opd(sif)) as ws:
        nii = ws.nii(file, sif)
        tmp_sif = nii[:-4] + '.sif'
        run_cmd.run(sing_turku + 'imgweigh ' + nii + ' ' + tmp_sif, diary_name)
        shutil.copyfile(tmp_sif, sif)
    return sif
