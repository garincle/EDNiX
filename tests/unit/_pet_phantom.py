"""
Synthetic data for the PET tests: an analytic head (scalp, skull, grey and white matter, ventricles,
cerebellum) sampled on an anatomical grid (AC-PC space) and on a PET grid placed with a known rigid pose,
a fake anatomical derivative of EDNiX (images, masks, 4D atlas, AC-PC and template matrices) and the
AFNI / Workbench / Turku commands emulated in numpy.
"""
import os
import re
import json
import shlex
import shutil

import numpy as np
import nibabel as nib
from scipy.ndimage import gaussian_filter
from scipy.spatial.transform import Rotation

opj = os.path.join

F = np.diag([-1.0, -1.0, 1.0])     # RAS <-> LPS (ANTs / ITK physical space)


def grid(n, vox, center=(0.0, 0.0, 0.0)):
    A = np.diag([vox, vox, vox, 1.0])
    A[:3, 3] = np.asarray(center) - vox * (n - 1) / 2.0
    return A


def world_points(shape, A):
    ijk = np.stack(np.meshgrid(*[np.arange(s) for s in shape], indexing='ij'), -1).reshape(-1, 3)
    return ijk @ A[:3, :3].T + A[:3, 3]


def _ell(p, c, r):
    return (((p - c) / r) ** 2).sum(-1) <= 1.0


def tissues(p):
    """Tissue labels of world points p (AC-PC space, RAS, mm): 0 air, 1 scalp, 2 skull, 3 CSF, 4 GM, 5 WM,
    6 cerebellum GM."""
    lab = np.zeros(len(p), np.int16)
    lab[_ell(p, (0, 0, 0), (30, 36, 32))] = 1
    lab[_ell(p, (0, 0, 0), (27, 33, 29))] = 2
    brain = _ell(p, (0, 0, 0), (24, 30, 26))
    lab[brain] = 4
    lab[_ell(p, (0, 0, 0), (20, 26, 22))] = 5
    lab[_ell(p, (0, 0, 2), (4, 8, 4))] = 3
    lab[brain & _ell(p, (0, -22, -14), (12, 7, 7))] = 6
    return lab


T1 = {0: 0, 1: 300, 2: 60, 3: 80, 4: 500, 5: 800, 6: 500}
PET = {0: 0, 1: 800, 2: 300, 3: 400, 4: 4000, 5: 1500, 6: 3600}   # Bq/mL, FDG-like


def pose(rot_deg=(8.0, -5.0, 12.0), t_mm=(6.0, -9.0, 4.0)):
    """x_acpc = R x_pet + t (RAS)."""
    return Rotation.from_euler('xyz', rot_deg, degrees=True).as_matrix(), np.asarray(t_mm, float)


def pet_frame(shape, A, R, t, scale=1.0, motion=None, smooth_mm=3.0):
    p = world_points(shape, A)
    if motion is not None:
        Rm, tm = motion
        p = p @ Rm.T + tm
    q = p @ R.T + t
    vals = np.vectorize(PET.get)(tissues(q)).astype(np.float32).reshape(shape) * scale
    return gaussian_filter(vals, smooth_mm / A[0, 0])


def save(data, A, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    nib.save(nib.Nifti1Image(np.asarray(data, np.float32), A), path)
    return path


def write_affine_mat(path, M_ras, t_ras):
    """ANTs matrix file of the RAS map x -> M x + t (fixed points -> moving points)."""
    import ants
    M = F @ M_ras @ F
    t = F @ t_ras
    tr = ants.create_ants_transform(transform_type='AffineTransform', dimension=3,
                                    matrix=M, translation=t, center=(0.0, 0.0, 0.0))
    ants.write_transform(tr, path)
    return path


def anat_derivative(bids, ID='01', ses='1', template_scale=1.08, template_shift=(2.0, -3.0, 1.5), tdir=None,
                    species='Test'):
    """EDNiX anatomical outputs of one subject + a 'template' space related to AC-PC by an affine map
    (template files in tdir, default <bids>/template, named <species>_...)."""
    dp = opj(bids, 'sub-' + ID, 'ses-' + ses)
    vol = opj(dp, 'anat', 'native', 'volumes')
    msk, lab, mat, prep = opj(vol, 'masks'), opj(vol, 'labels'), opj(dp, 'anat', 'matrices'), opj(dp, 'anat', 'preprocessing')
    A = grid(72, 1.0)
    shape = (72, 72, 72)
    t = tissues(world_points(shape, A)).reshape(shape)
    t1 = gaussian_filter(np.vectorize(T1.get)(t).astype(np.float32), 0.7)
    files = dict(acpc=save(t1, A, opj(vol, ID + '_space-acpc_desc-template_T1w.nii.gz')),
                 brain=save(t1 * (t >= 3), A, opj(vol, ID + '_space-acpc_desc-SS_T1w.nii.gz')),
                 mask=save((t >= 3), A, opj(msk, ID + '_space-acpc_mask.nii.gz')),
                 gray=save(np.isin(t, (4, 6)), A, opj(msk, ID + '_desc-Cerebral_Gray_mask.nii.gz')),
                 white=save(t == 5, A, opj(msk, ID + '_desc-erod-White_mask.nii.gz')),
                 csf=save(t == 3, A, opj(msk, ID + '_desc-erod-Vent_mask.nii.gz')))
    # 4D atlas, level 0: cerebrum (1) / cerebellum (7); level 1: left / right hemispheres of the cortex, WM, CSF
    p = world_points(shape, A)
    lvl0 = np.where(t == 6, 7, np.where(t >= 3, 1, 0)).astype(np.int16)
    left = (p[:, 0] < 0).reshape(shape)
    lvl1 = np.zeros(shape, np.int16)
    lvl1[(t == 4) & left] = 10
    lvl1[(t == 4) & ~left] = 1010
    lvl1[t == 5] = 20
    lvl1[t == 3] = 30
    lvl1[t == 6] = 7
    atlas = save(np.stack([lvl0, lvl1], -1), A, opj(lab, ID + '_seg-TEST_dseg.nii.gz'))
    files.update(labels_dir=lab, dir_transfo=mat, ID=ID, atlas=atlas,
                 reorient=save(t1, A, opj(prep, ID + '_space-raw_desc-reorient_T1w.nii.gz')),
                 raw=save(t1, A, opj(prep, ID + '_space-raw_desc-template_T1w.nii.gz')), surfaces=[])
    os.makedirs(mat, exist_ok=True)
    write_affine_mat(opj(mat, 'acpc_0GenericAffine.mat'), np.eye(3), np.zeros(3))
    # template: x_native = S x_template + s  (fixed = template, moving = native)
    S = np.eye(3) / template_scale
    s = np.asarray(template_shift)
    write_affine_mat(opj(mat, 'native_to_EDNiX_Final_0GenericAffine.mat'), S, s)
    tdir = tdir or opj(bids, 'template')
    At = grid(80, 1.0)
    qt = world_points((80,) * 3, At) @ S.T + s
    tt = tissues(qt).reshape((80,) * 3)
    tpl = save(gaussian_filter(np.vectorize(T1.get)(tt).astype(np.float32), 0.7) * (tt >= 3), At,
               opj(tdir, species + '_space-acpc_desc-SS_T1w.nii.gz'))
    save(gaussian_filter(np.vectorize(T1.get)(tt).astype(np.float32), 0.7), At,
         opj(tdir, species + '_space-acpc_desc-template_T1w.nii.gz'))
    tmask = save(tt >= 3, At, opj(tdir, 'masks', species + '_mask.nii.gz'))
    save(np.isin(tt, (4, 6)), At, opj(tdir, 'masks', species + '_desc-Cerebral_Gray_mask.nii.gz'))
    save(tt == 5, At, opj(tdir, 'masks', species + '_desc-erod-White_mask.nii.gz'))
    save(tt == 3, At, opj(tdir, 'masks', species + '_desc-erod-Vent_mask.nii.gz'))
    # template atlas (for the export of the labels)
    pt = world_points((80,) * 3, At)
    l0 = np.where(tt == 6, 7, np.where(tt >= 3, 1, 0)).astype(np.int16)
    l1 = np.where(tt == 6, 7, np.where(tt == 5, 20, np.where(tt == 3, 30, np.where(tt == 4, np.where(
        pt[:, 0].reshape((80,) * 3) < 0, 10, 1010), 0)))).astype(np.int16)
    save(np.stack([l0, l1], -1), At, opj(tdir, 'labels', species + '_seg-TEST_dseg.nii.gz'))
    files.update(template=tpl, template_mask=tmask, template_S=S, template_s=s,
                 template_label_base=opj(tdir, 'labels', species + '_seg-'))
    package = ['EDNiX', tpl, tmask, files['template_label_base'], '', '', '',
               [opj(mat, 'native_to_EDNiX_Final_0GenericAffine.mat')], [True],
               [opj(mat, 'native_to_EDNiX_Final_0GenericAffine.mat')], [False]]
    return files, package, dp


def dynamic_pet(path, R, t, n=40, vox=2.0, frames=((0, 30), (30, 30), (60, 60), (120, 120), (240, 300),
                                                     (540, 600), (1140, 600), (1740, 600)),
                motions=None, units='Bq/mL', dose=100.0, decay_corrected=True, extra=None):
    """Dynamic PET (activity rising as 1 - exp(-t/10 min)) with optional per-frame rigid motions (PET space)."""
    A = grid(n, vox)
    data = []
    for k, (s, d) in enumerate(frames):
        mid = s + d / 2.0
        scale = 1.0 - np.exp(-mid / 600.0)
        data.append(pet_frame((n,) * 3, A, R, t, scale, None if motions is None else motions[k]))
    save(np.stack(data, -1), A, path)
    info = dict(FrameTimesStart=[f[0] for f in frames], FrameDuration=[f[1] for f in frames], Units=units,
                InjectedRadioactivity=dose, InjectedRadioactivityUnits='MBq', ImageDecayCorrected=decay_corrected,
                ImageDecayCorrectionTime=0, InjectionStart=0, TimeZero='10:00:00', ScanStart=0,
                TracerName='FDG', TracerRadionuclide='F18')
    info.update(extra or {})
    with open(path.replace('.nii.gz', '.json'), 'w') as f:
        json.dump(info, f)
    return path, A, info


# ---------------------------------------------------------------------------------------------- command emulator
CALLS = []


def _axcodes(A):
    return nib.orientations.aff2axcodes(A)


OPP = {'R': 'L', 'L': 'R', 'A': 'P', 'P': 'A', 'S': 'I', 'I': 'S'}
DIR = {'R': (1, 0, 0), 'L': (-1, 0, 0), 'A': (0, 1, 0), 'P': (0, -1, 0), 'S': (0, 0, 1), 'I': (0, 0, -1)}


def _afni_code(A):
    return ''.join(OPP[c] for c in _axcodes(A))


def emulate(cmd, diary=None):
    tok = shlex.split(cmd)
    prog = tok[0]
    CALLS.append(prog)
    if prog == '3dinfo':
        img = nib.load(tok[-1])
        if tok[1] == '-is_oblique':
            M = img.affine[:3, :3] / np.linalg.norm(img.affine[:3, :3], axis=0)
            return (b'1\n' if np.abs(np.abs(M).max(0) - 1).max() > 1e-3 else b'0\n'), b''
        if tok[1] == '-orient':
            return (_afni_code(img.affine) + '\n').encode(), b''
    if prog == '3drefit':
        f = tok[-1]
        img = nib.load(f)
        if '-orient' in tok:
            new = tok[tok.index('-orient') + 1]
            old = _afni_code(img.affine)
            P = sum(np.outer(DIR[OPP[n]], DIR[OPP[o]]) for n, o in zip(new, old))
            A = img.affine.copy()
            A[:3, :] = P @ A[:3, :]
            nib.save(nib.Nifti1Image(np.asarray(img.dataobj), A, img.header), f)
        return b'', b''
    if prog == '3dWarp':
        shutil.copyfile(tok[-1], tok[tok.index('-prefix') + 1])
        return b'', b''
    if prog in ('wb_command', 'dftweigh', 'freeview'):
        return b'', b''
    if prog == 'imgweigh':
        return b'', b''
    if prog == 'bfmsrtm':
        args = [a for a in tok[1:] if not a.startswith('-')]
        fit = next(a.split('=', 1)[1] for a in tok if a.startswith('-fit='))
        shutil.copyfile(args[0], fit)
        with open(args[2], 'w') as fh:
            fh.write('emulated bfmsrtm\n')
        return b'', b''
    if prog in ('imgbfbp', 'imgsrtm', 'imgdv'):
        args = [a for a in tok[1:] if not a.startswith('-')]
        if prog == 'imgdv':
            img_in, out = args[1], args[3]
        else:
            img_in, out = args[0], args[2]
        im = nib.load(img_in)
        d = np.asarray(im.dataobj)[..., -1]
        nib.save(nib.Nifti1Image((d > 0).astype(np.float32), im.affine), out)
        for a in tok:
            m = re.match(r'-(err|n)=(.+)', a)
            if m:
                nib.save(nib.Nifti1Image(np.zeros(d.shape, np.float32), im.affine), m.group(2))
        return b'', b''
    raise RuntimeError('command not emulated: ' + cmd)


def patch_run_cmd(monkeypatch):
    from Tools import run_cmd
    monkeypatch.setattr(run_cmd, 'run', emulate)
    monkeypatch.setattr(run_cmd, 'get', emulate)
    monkeypatch.setattr(run_cmd, 'do', lambda cmd, diary=None: emulate(cmd, diary)[0].decode())


LABELS = [(1, 'Cerebrum'), (7, 'Cerebellum'), (10, 'L_Ctx'), (1010, 'R_Ctx'), (20, 'WM'), (30, 'CSF')]


def atlas_library(root, species='Testmonkey', reference='EDNiX'):
    """Minimal Atlases_library: <root>/Atlases_library/atlas/Primate/<species>/<reference>/volumes + label_code.
    Returns (MAIN_PATH, template folder)."""
    lib = opj(root, 'Atlases_library', 'atlas', 'Primate', species)
    tdir = opj(lib, reference, 'volumes')
    os.makedirs(opj(lib, 'label_code'), exist_ok=True)
    with open(opj(lib, 'label_code', 'TEST_label.txt'), 'w') as f:
        for lid, name in LABELS:
            f.write(f'{name}\n{lid} 100 100 100 255\n')
    main = opj(root, 'EDNiX')
    os.makedirs(main, exist_ok=True)
    return main, tdir


def no_containers(monkeypatch):
    from Tools import Load_EDNiX_requirement
    monkeypatch.setattr(Load_EDNiX_requirement, 'load_requirement', lambda *a, **k: ('', '', '', '', '', '', '', ''))
    from Tools import QC_plot
    monkeypatch.setattr(QC_plot, 'mosaic', lambda a, b, png: open(png, 'w').close())


def bids_table(entries):
    import pandas as pd
    return pd.DataFrame([dict(subject=s, session=e, modality='pet') for s, e in entries])
