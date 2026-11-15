import os
import numpy as np
import ants
import pytest
from modalities.anat import norm2template as N


def _images(tmp_path):
    a = np.zeros((40, 48, 36), dtype='float32')
    zz, yy, xx = np.meshgrid(*[np.arange(s) for s in a.shape], indexing='ij')
    a[((zz - 20) / 12) ** 2 + ((yy - 24) / 16) ** 2 + ((xx - 18) / 10) ** 2 < 1] = 100
    a[((zz - 20) / 4) ** 2 + ((yy - 30) / 5) ** 2 + ((xx - 18) / 3) ** 2 < 1] = 40   # asymmetric feature
    fixed = ants.from_numpy(a, spacing=(0.2, 0.2, 0.2))
    moving = ants.from_numpy(np.roll(a, (3, -4, 2), (0, 1, 2)), spacing=(0.2, 0.2, 0.2), origin=(1.0, -0.5, 0.3))
    paths = {}
    for k, im in (('fixed', fixed), ('moving', moving)):
        paths[k] = str(tmp_path / (k + '.nii.gz'))
        ants.image_write(im, paths[k])
    return paths


def _quiet(monkeypatch):
    monkeypatch.setattr(N.run_cmd, 'msg', lambda *a, **k: None)


def test_acpc_keeps_default_when_brain_in_grid(tmp_path, monkeypatch):
    _quiet(monkeypatch)
    p = _images(tmp_path)
    calls = []
    real = N._acpc_candidate
    monkeypatch.setattr(N, '_acpc_candidate', lambda kind, seed, *a: calls.append((kind, seed)) or real(kind, seed, *a))
    out = N.acpc('sub', p['moving'], '', str(tmp_path), p['fixed'], str(tmp_path / 'mat'), 'Rigid', 'MI', 'MI',
                 'linear', None, '_desc-64_T2w')
    assert calls == [('default', 1)]
    assert os.path.isfile(out[0]) and os.path.isfile(tmp_path / 'sub_space-acpc_desc-64_T2w.nii.gz')
    kept, dice = N._brain_kept_and_dice(ants.image_read(p['fixed']), ants.image_read(p['moving']), out)
    assert kept > 0.99 and dice > 0.9
    assert not [f for f in os.listdir(tmp_path / 'mat') if f.startswith('acpc_try')]
    assert 'ANTS_RANDOM_SEED' not in os.environ


def test_acpc_retries_and_keeps_best_when_default_pushes_brain_out(tmp_path, monkeypatch):
    _quiet(monkeypatch)
    p = _images(tmp_path)
    scores = iter([(0.33, 0.2), (0.85, 0.6), (0.97, 0.8), (0.95, 0.7), (0.99, 0.75), (0.96, 0.7)])
    calls = []

    def cand(kind, seed, fixed, moving, mmask, prefix, *a):
        calls.append((kind, seed))
        tx = prefix + 'R_0GenericAffine.mat'
        ants.write_transform(ants.create_ants_transform(transform_type='Euler3DTransform', dimension=3), tx)
        return [tx]
    monkeypatch.setattr(N, '_acpc_candidate', cand)
    monkeypatch.setattr(N, '_brain_kept_and_dice', lambda f, m, tx: next(scores))
    N.acpc('sub', p['moving'], '', str(tmp_path), p['fixed'], str(tmp_path / 'mat'), 'Rigid', 'MI', 'MI',
           'linear', None, '_desc-64_T2w')
    assert calls == [('default', 1), ('rotsearch', 1), ('default', 2), ('rotsearch', 2), ('default', 3), ('rotsearch', 3)]
    desc = open(tmp_path / 'sub_space-acpc_desc-64_T2w.json').read()
    assert 'default seed 2' in desc          # highest Dice


def test_acpc_stops_after_rotsearch_when_it_fixes_it(tmp_path, monkeypatch):
    _quiet(monkeypatch)
    p = _images(tmp_path)
    scores = iter([(0.33, 0.2), (0.98, 0.8)])
    calls = []

    def cand(kind, seed, fixed, moving, mmask, prefix, *a):
        calls.append((kind, seed))
        tx = prefix + 'R_0GenericAffine.mat'
        ants.write_transform(ants.create_ants_transform(transform_type='Euler3DTransform', dimension=3), tx)
        return [tx]
    monkeypatch.setattr(N, '_acpc_candidate', cand)
    monkeypatch.setattr(N, '_brain_kept_and_dice', lambda f, m, tx: next(scores))
    N.acpc('sub', p['moving'], '', str(tmp_path), p['fixed'], str(tmp_path / 'mat'), 'Rigid', 'MI', 'MI',
           'linear', None, '_desc-64_T2w')
    assert calls == [('default', 1), ('rotsearch', 1)]
    assert 'rotsearch seed 1' in open(tmp_path / 'sub_space-acpc_desc-64_T2w.json').read()


def test_brain_kept_detects_brain_out_of_grid(tmp_path):
    p = _images(tmp_path)
    f, m = ants.image_read(p['fixed']), ants.image_read(p['moving'])
    t = ants.create_ants_transform(transform_type='Euler3DTransform', dimension=3, translation=(5.0, 0, 0))
    ants.write_transform(t, str(tmp_path / 'shift.mat'))
    kept, _ = N._brain_kept_and_dice(f, m, [str(tmp_path / 'shift.mat')])
    assert kept < N.ACPC_MIN_KEPT
