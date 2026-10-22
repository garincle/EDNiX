"""
Partial-volume correction of a static PET image (optional, PETPVC: Thomas et al. 2016, Phys Med Biol
61:7975).

The PET engineer passed the voxel size as the point-spread function: the PSF is the spatial resolution
of the scanner and reconstruction (FWHM, mm), a property of the acquisition that must be given
(pvc_fwhm in the launcher); the correction is skipped without it.

Tissue classes given to PETPVC: a 4D image of fractions summing to one in every voxel,
    GM   whole-brain grey matter (cortex, deep nuclei, cerebellum: Whole_Brain_Gray mask)
    WM   white matter (non-eroded mask), outside the grey matter
    CSF  rest of the brain mask (sulcal CSF and ventricles)
    rest outside the brain mask
computed from the fractions of the PET voxels covered by each mask (coreg2anat.tissue_fraction). An eroded
white-matter mask or a ventricle-only CSF class (the masks of the functional pipeline) would leave tissue in
the 'rest' class and bias the correction.
"""
import os

import numpy as np

opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu


def tissue_classes(probseg, img):
    """[GM, WM, CSF, rest] fractions (list of 3D arrays) on the grid of `img`."""
    for k in ('GM', 'WM', 'brain'):
        if not probseg.get(k) or not opi(probseg[k]):
            raise FileNotFoundError(f'tissue fraction {k} in the PET space needed for the partial-volume correction '
                                    '(run step 6 again)')
    ref, _ = pu.load(img)
    fr = {}
    for k in ('GM', 'WM', 'brain'):
        im, a = pu.load(probseg[k])
        if im.shape[:3] != ref.shape[:3] or not np.allclose(im.affine, ref.affine, atol=1e-3):
            raise ValueError(f'{probseg[k]} is not on the grid of {img}')
        fr[k] = np.clip(a if a.ndim == 3 else a[..., 0], 0.0, 1.0)
    gm = fr['GM']
    wm = np.minimum(fr['WM'], 1.0 - gm)
    csf = np.clip(fr['brain'] - gm - wm, 0.0, 1.0)
    rest = np.clip(1.0 - gm - wm - csf, 0.0, 1.0)
    return [gm, wm, csf, rest]


def correct(img, probseg, fwhm_mm, out, diary_name, algo='RBV'):
    """
    algo: MG, GTM, LABBE, MTC, RBV, IY... (PETPVC names). probseg: dict GM / WM / brain of tissue fractions in
    the PET space (registration record, pet_space['probseg']). Needs nipype and the petpvc program.
    """
    try:
        from nipype.interfaces.petpvc import PETPVC
    except ImportError as e:
        raise ImportError('partial-volume correction needs nipype and PETPVC (pip install nipype; petpvc on the '
                          'PATH)') from e
    if fwhm_mm is None:
        raise ValueError('pvc_fwhm (scanner resolution, FWHM in mm) is needed for the partial-volume correction')
    fw = np.atleast_1d(np.asarray(fwhm_mm, dtype=float))
    fw = np.repeat(fw, 3) if fw.size == 1 else fw
    if fw.size != 3 or np.any(fw <= 0):
        raise ValueError(f'pvc_fwhm={fwhm_mm!r}: one value or three positive values (mm) expected')
    nii, d = pu.load(img)
    seg_file = out.replace('.nii.gz', '_seg4D.nii.gz')
    pu.save_like(np.stack(tissue_classes(probseg, img), axis=3), nii, seg_file)
    run_cmd.msg(f'Partial-volume correction ({algo}, PSF FWHM {fw.tolist()} mm)', diary_name, 'OKGREEN')
    p = PETPVC()
    p.inputs.in_file = img
    p.inputs.mask_file = seg_file
    p.inputs.out_file = out
    p.inputs.pvc = algo
    p.inputs.fwhm_x, p.inputs.fwhm_y, p.inputs.fwhm_z = [float(x) for x in fw]
    p.run()
    pu.write_json(out, [img, seg_file], f'Partial-volume correction (PETPVC {algo}, PSF FWHM {fw.tolist()} mm; '
                                        'classes GM, WM, CSF, rest).')
    return out
