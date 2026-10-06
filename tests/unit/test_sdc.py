import numpy as np
import nibabel as nib
import ants
from scipy import ndimage

from modalities.fMRI import sdc_field as S
import modalities.fMRI._2_coregistration_to_norm as C


def _phantom(tmp_path):
    x, y, z = np.mgrid[0:40, 0:48, 0:32].astype(float)
    true = np.exp(-((x - 20) ** 2 / 70 + (y - 24) ** 2 / 120 + (z - 16) ** 2 / 50)) * 100 \
        + 30 * np.exp(-((x - 26) ** 2 / 10 + (y - 16) ** 2 / 12 + (z - 18) ** 2 / 10))
    shift = 3.0 * np.exp(-((x - 17) ** 2 / 50 + (y - 15) ** 2 / 35 + (z - 16) ** 2 / 50))
    dist = ndimage.map_coordinates(true, np.array([x, y - shift, z]), order=3)   # distortion along y
    A = np.diag([2.0, 2.0, 2.0, 1.0])
    A[:3, 3] = [-40, -48, -32]
    paths = {}
    for n, a in (('true', true), ('dist', dist), ('shift', shift)):
        paths[n] = str(tmp_path / f'{n}.nii.gz')
        nib.save(nib.Nifti1Image(a.astype('f4'), A), paths[n])
    return paths


def _corr(a, b):
    return np.corrcoef(a.numpy().ravel(), b.numpy().ravel())[0, 1]


def test_sign_is_found_and_the_field_undistorts(tmp_path):
    p = _phantom(tmp_path)
    info = S.determine_sign(p['dist'], p['true'], p['shift'], axis=1)
    assert info['status'] == 'ok'
    S.write_field(S.shift_to_field(p['shift'], 1, info['sign']), str(tmp_path / 'f.nii.gz'))
    D, T = ants.image_read(p['dist']), ants.image_read(p['true'])
    out = ants.apply_transforms(fixed=D, moving=D, transformlist=[str(tmp_path / 'f.nii.gz')])
    assert _corr(out, T) > _corr(D, T) and _corr(out, T) > 0.999


def test_negligible_distortion_is_not_applied(tmp_path):
    p = _phantom(tmp_path)
    nib.save(nib.Nifti1Image(np.zeros((40, 48, 32), 'f4'), nib.load(p['shift']).affine), p['shift'])
    assert S.determine_sign(p['dist'], p['true'], p['shift'], axis=1)['status'] == 'negligible'


def test_one_shot_chain_order_with_explicit_whichtoinvert(tmp_path):
    # [run_inRef (affine), SDC field]: with whichtoinvert=None ANTsPy would invert the affine
    p = _phantom(tmp_path)
    info = S.determine_sign(p['dist'], p['true'], p['shift'], axis=1)
    tx = ants.create_ants_transform(transform_type='AffineTransform', dimension=3)
    tx.set_parameters([1, 0, 0, 0, 1, 0, 0, 0, 1, 4.0, 0, 0])
    tx.set_fixed_parameters([0, 0, 0])
    ants.write_transform(tx, str(tmp_path / 'rigid.mat'))
    D, T, SH = ants.image_read(p['dist']), ants.image_read(p['true']), ants.image_read(p['shift'])
    run = ants.apply_transforms(fixed=D, moving=D, transformlist=[str(tmp_path / 'rigid.mat')], whichtoinvert=[True])
    shift_run = ants.apply_transforms(fixed=SH, moving=SH, transformlist=[str(tmp_path / 'rigid.mat')], whichtoinvert=[True])
    S.write_field(S.shift_to_field(shift_run, 1, info['sign']), str(tmp_path / 'f_run.nii.gz'))
    chain = [str(tmp_path / 'rigid.mat'), str(tmp_path / 'f_run.nii.gz')]
    out = ants.apply_transforms(fixed=T, moving=run, transformlist=chain, whichtoinvert=[False, False])
    rigid_only = ants.apply_transforms(fixed=T, moving=run, transformlist=chain[:1])
    assert _corr(out, T) > 0.999 and _corr(out, T) > _corr(rigid_only, T)


def test_build_sdc_field_writes_field_in_the_reoriented_space(tmp_path, monkeypatch):
    p = _phantom(tmp_path)
    monkeypatch.setattr(C._2b_fix_orient, 'fix_orient',
                        lambda out, inp, *a, **k: nib.save(nib.load(inp), out))   # header-only identity
    out = C._build_sdc_field('sub-01_run-1', 0, p['dist'], p['true'], p['shift'], 'y', str(tmp_path),
                             str(tmp_path), [], '', '', '', 'header', '', str(tmp_path / 'diary.txt'))
    assert out and ants.image_read(out).components == 3
    monkeypatch.setenv('EDNIX_SDC', '0')
    assert C._build_sdc_field('sub-01_run-1', 0, p['dist'], p['true'], p['shift'], 'y', str(tmp_path),
                              str(tmp_path), [], '', '', '', 'header', '', str(tmp_path / 'diary.txt')) is None


def test_which_field_each_run_uses(tmp_path):
    f = lambda r: tmp_path / f'{r}_space-func_desc-sdc_1Warp.nii.gz'
    f('run1').write_bytes(b'x')                          # reference run: field accepted
    d = str(tmp_path)
    assert S.run_field(d, 'run1', 'run1', 'new') == str(f('run1'))
    assert S.run_field(d, 'run2', 'run1', 'new') is None          # own field rejected: no borrowing
    assert S.run_field(d, 'run2', 'run1', 'old') == str(f('run1'))    # shared field map
    assert S.run_field(d, 'run2', 'run1', '2_mapdir') == str(f('run1'))
    assert S.run_field(d, 'run1', 'run1', 'very_old') is None


def test_small_distortion_sign_is_found_from_opposite_effects(tmp_path):
    p = _phantom(tmp_path)
    sh = nib.load(p['shift'])
    small = np.asarray(sh.dataobj) * 0.15                     # max 0.45 voxel
    x, y, z = np.mgrid[0:40, 0:48, 0:32].astype(float)
    true = np.asarray(nib.load(p['true']).dataobj)
    dist = ndimage.map_coordinates(true, np.array([x, y - small, z]), order=3)
    nib.save(nib.Nifti1Image(small.astype('f4'), sh.affine), p['shift'])
    nib.save(nib.Nifti1Image(dist.astype('f4'), sh.affine), p['dist'])
    info = S.determine_sign(p['dist'], p['true'], p['shift'], axis=1)
    assert info['status'] == 'ok' and info['sign'] == 1


def test_registration_restriction_follows_the_phase_encoding_axis():
    I = np.eye(3)
    assert S.pe_restriction('j-', I) == (0.1, 1, 0.1)
    assert S.pe_restriction('y', I) == (0.1, 1, 0.1)
    assert S.pe_restriction('Auto', I) == (1, 1, 1)
    # animal data reoriented by a header relabel: voxel axis j now points along physical z
    D = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], float)
    assert S.pe_restriction('j', D) == (0.1, 0.1, 1)
