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

Two engines (pvc_engine):
    'petpvc'  the PETPVC program through nipype (every PETPVC method);
    'python'  RBV and MG written here with numpy / scipy, same inputs and same definitions as PETPVC, no external
              program (validated on a phantom of known uptake: tests/unit/test_pet.py);
    'auto'    PETPVC when the program `petpvc` is on the PATH, else Python (RBV, MG only).
The engine used is written in the JSON sidecar of the result.
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


PYTHON_METHODS = ('RBV', 'MG')
ENGINES = ('auto', 'python', 'petpvc')


def _fwhm3(fwhm_mm):
    fw = np.atleast_1d(np.asarray(fwhm_mm, dtype=float))
    fw = np.repeat(fw, 3) if fw.size == 1 else fw
    if fw.size != 3 or np.any(fw <= 0):
        raise ValueError(f'pvc_fwhm={fwhm_mm!r}: one value or three positive values (mm) expected')
    return fw


def blur(a, fwhm_mm, spacing):
    """Gaussian point-spread function of FWHM fwhm_mm (3 values, mm, along the image axes); zero outside the grid."""
    from scipy import ndimage
    sigma = [f / (2.0 * np.sqrt(2.0 * np.log(2.0))) / float(v) for f, v in zip(fwhm_mm, spacing)]
    return ndimage.gaussian_filter(np.asarray(a, dtype=np.float64), sigma, mode='constant')


def gtm(pet, classes, fwhm_mm, spacing):
    """
    Geometric transfer matrix (Rousset et al. 1998, J Nucl Med 39:904) with tissue fractions as weights: the observed
    mean of class j is sum_i G[j, i] t_i, G[j, i] = sum_v f_j(v) PSF*f_i(v) / sum_v f_j(v). Returns the true class
    means t (least squares) and the observed ones.
    """
    w = [np.asarray(c, dtype=np.float64) for c in classes]
    blurred = [blur(c, fwhm_mm, spacing) for c in w]
    sw = np.array([max(c.sum(), 1e-12) for c in w])
    g = np.array([[float((wj * bi).sum()) / sj for bi in blurred] for wj, sj in zip(w, sw)])
    obs = np.array([float((wj * pet).sum()) / sj for wj, sj in zip(w, sw)])
    t, *_ = np.linalg.lstsq(g, obs, rcond=None)
    return t, obs, blurred


def rbv(pet, classes, fwhm_mm, spacing):
    """
    Region-based voxel-wise correction (Thomas et al. 2011, Eur J Nucl Med Mol Imaging 38:1104): the GTM means
    painted on the tissue fractions give a synthetic image S = sum_i t_i f_i; corrected = PET * S / (PSF*S).
    Voxels where PSF*S is ~0 (outside the object) are left as they are.
    """
    t, _, _ = gtm(pet, classes, fwhm_mm, spacing)
    syn = sum(ti * np.asarray(c, dtype=np.float64) for ti, c in zip(t, classes))
    bs = blur(syn, fwhm_mm, spacing)
    ok = bs > 1e-3 * max(float(np.abs(bs).max()), 1e-12)
    out = np.where(ok, pet * np.divide(syn, bs, out=np.ones_like(syn), where=ok), pet)
    return out, t


def mg(pet, classes, fwhm_mm, spacing, min_gm=0.1):
    """
    Muller-Gartner et al. 1992 (J Cereb Blood Flow Metab 12:571): the spill-in of the other classes (their GTM means
    times their smoothed fractions) is removed and the grey matter rescaled by its smoothed fraction,
    GM = (PET - sum_{i != GM} t_i PSF*f_i) / PSF*f_GM, where PSF*f_GM >= min_gm; 0 elsewhere. Class 0 is the grey
    matter. The other-class means come from the GTM (as the 'MG' of PETPVC uses region means).
    """
    t, _, blurred = gtm(pet, classes, fwhm_mm, spacing)
    other = sum(ti * bi for ti, bi in zip(t[1:], blurred[1:]))
    bgm = blurred[0]
    ok = bgm >= min_gm
    return np.where(ok, (pet - other) / np.where(ok, bgm, 1.0), 0.0), t


def engine_for(algo, engine='auto'):
    """'python' or 'petpvc' for this method and request; raises when the request cannot be met."""
    import shutil
    if engine not in ENGINES:
        raise ValueError(f'pvc_engine={engine!r}: expected one of {ENGINES}')
    have = shutil.which('petpvc') is not None
    if engine == 'petpvc' and not have:
        raise FileNotFoundError("pvc_engine='petpvc' but the program petpvc is not on the PATH (conda install -c "
                                "aramislab petpvc, or build https://github.com/UCL/PETPVC); pvc_engine='python' "
                                f"runs {PYTHON_METHODS} without it")
    if engine == 'petpvc' or (engine == 'auto' and have):
        return 'petpvc'
    if algo not in PYTHON_METHODS:
        raise ValueError(f'pvc_method={algo!r} needs the PETPVC program (not on the PATH); without it: '
                         f'{PYTHON_METHODS}')
    return 'python'


def correct(img, probseg, fwhm_mm, out, diary_name, algo='RBV', engine='auto'):
    """
    algo: MG, GTM, LABBE, MTC, RBV, IY... (PETPVC names; RBV and MG without PETPVC). probseg: dict GM / WM / brain of
    tissue fractions in the PET space (registration record, pet_space['probseg']). engine: see the module docstring.
    """
    if fwhm_mm is None:
        raise ValueError('pvc_fwhm (scanner resolution, FWHM in mm) is needed for the partial-volume correction')
    fw = _fwhm3(fwhm_mm)
    eng = engine_for(algo, engine)
    nii, d = pu.load(img)
    classes = tissue_classes(probseg, img)
    seg_file = out.replace('.nii.gz', '_seg4D.nii.gz')
    pu.save_like(np.stack(classes, axis=3), nii, seg_file)
    run_cmd.msg(f'Partial-volume correction ({algo}, PSF FWHM {fw.tolist()} mm, engine {eng})', diary_name, 'OKGREEN')
    extra = {}
    if eng == 'petpvc':
        from nipype.interfaces.petpvc import PETPVC
        p = PETPVC()
        p.inputs.in_file = img
        p.inputs.mask_file = seg_file
        p.inputs.out_file = out
        p.inputs.pvc = algo
        p.inputs.fwhm_x, p.inputs.fwhm_y, p.inputs.fwhm_z = [float(x) for x in fw]
        p.run()
    else:
        fn = rbv if algo == 'RBV' else mg
        res, t = fn(d.astype(np.float64), classes, fw, nii.header.get_zooms()[:3])
        pu.save_like(res, nii, out)
        extra['ClassMeansGTM'] = dict(zip(('GM', 'WM', 'CSF', 'rest'), [round(float(v), 6) for v in t]))
    pu.write_json(out, [img, seg_file], f'Partial-volume correction ({algo}, PSF FWHM {fw.tolist()} mm; classes GM, '
                                        f'WM, CSF, rest; engine {eng}).', PVCEngine=eng, **extra)
    return out
