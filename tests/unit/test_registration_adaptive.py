"""SyNadaptive: candidate choice, standard output names, inheritance (synthetic volumes, no download)."""
import os
import json
import numpy as np
import pytest

ants = pytest.importorskip('ants')

from modalities.anat import registration_adaptive, transfoparams, norm2template


def _brain(center, radii, shape=(32, 32, 32)):
    """Three nested ellipsoids (CSF-like rim, GM-like shell, WM-like core) + a mask."""
    z, y, x = np.mgrid[:shape[0], :shape[1], :shape[2]]
    d = ((z - center[0]) / radii[0]) ** 2 + ((y - center[1]) / radii[1]) ** 2 + ((x - center[2]) / radii[2]) ** 2
    img = np.where(d < 1, 30., 0.) + np.where(d < 0.6, 40., 0.) + np.where(d < 0.25, 50., 0.)
    img = ants.smooth_image(ants.from_numpy(img.astype('float32')), 0.7)
    return img, ants.from_numpy((d < 1).astype('float32'))


@pytest.fixture(scope='module')
def pair():
    fixed, fmask = _brain((16, 16, 16), (11, 10, 9))
    moving, mmask = _brain((17, 15, 17), (9, 11, 8))
    return fixed, fmask, moving, mmask


def _translation(fixed, moving, prefix):
    return ants.registration(fixed, moving, 'Translation', outprefix=prefix)['fwdtransforms']


def test_A_kept_when_checks_pass(tmp_path, pair):
    fixed, fmask, moving, mmask = pair
    pre = str(tmp_path / 'xfm_')
    res = registration_adaptive.register(fixed, moving, _translation(fixed, moving, str(tmp_path / 'T_')), pre,
                                         fmask, mmask, thresholds=dict(tissue_dice=0.0, mask_dice=0.0, max_folding=1.0))
    rep = json.load(open(pre + 'adaptive_QC.json'))
    assert rep['chosen'] == 'A_SyNCC' and 'B_refine' not in rep
    assert res['fwdtransforms'] == [pre + '1Warp.nii.gz', pre + '0GenericAffine.mat']
    assert res['invtransforms'] == [pre + '0GenericAffine.mat', pre + '1InverseWarp.nii.gz']
    assert all(os.path.isfile(f) for f in res['fwdtransforms'] + res['invtransforms'])
    assert rep['A_SyNCC']['mask_dice'] > 0.9


def test_B_tried_when_checks_fail(tmp_path, pair):
    fixed, fmask, moving, mmask = pair
    pre = str(tmp_path / 'xfm_')
    res = registration_adaptive.register(fixed, moving, _translation(fixed, moving, str(tmp_path / 'T_')), pre,
                                         fmask, mmask, thresholds=dict(tissue_dice=1.01))
    rep = json.load(open(pre + 'adaptive_QC.json'))
    assert 'B_refine' in rep and rep['chosen'] in ('A_SyNCC', 'B_refine')
    assert all(os.path.isfile(f) for f in res['fwdtransforms'] + res['invtransforms'])
    assert not any(n.startswith('xfm_refine_') for n in os.listdir(tmp_path))  # B's own files are cleaned


def test_norm_with_SyNadaptive_writes_standard_outputs(tmp_path, pair):
    fixed, fmask, moving, mmask = pair
    for name, img in (('tpl', fixed), ('tplmask', fmask), ('sub', moving), ('submask', mmask)):
        ants.image_write(img, str(tmp_path / f'{name}.nii.gz'))
    out = tmp_path / 'xfm'
    mTx = norm2template.norm('sub', str(tmp_path / 'sub.nii.gz'), str(tmp_path / 'submask.nii.gz'), str(tmp_path),
                             'tpl', str(tmp_path / 'tpl.nii.gz'), str(tmp_path / 'tplmask.nii.gz'),
                             str(out), 'SyNadaptive', 'sub_T', 'sub_S', 'mattes', 'mattes', 'linear',
                             str(tmp_path / 'diary.txt'), '', '_desc-test', 0)
    for f in ('sub_S_0GenericAffine.mat', 'sub_S_1Warp.nii.gz', 'sub_S_1InverseWarp.nii.gz', 'sub_S_adaptive_QC.json'):
        assert (out / f).is_file()
    assert mTx['fwdtransforms'][1].endswith('sub_S_0GenericAffine.mat')  # FS_2_WB reads the affine at index 1
    assert (tmp_path / 'sub_space-tpl_desc-test.nii.gz').is_file()


def test_inherited_blocks_use_SyNCC():
    _, blocks = transfoparams.set('do_as_I_said', False, False, '', '', '', '',
                                  [{'name': 'coreg', 'type_of_transform': 'SyNadaptive'}])
    types = {b['name']: b['type_of_transform'] for b in blocks}
    assert types['coreg'] == 'SyNadaptive'
    assert types['SS1'] == types['SS2'] == types['SS3'] == types['stdyT'] == 'SyNCC'
    assert registration_adaptive.ants_type('SyNadaptive') == 'SyNCC'
    assert registration_adaptive.ants_type('Rigid') == 'Rigid'
