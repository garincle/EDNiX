"""Manual brain mask of the PET (pet/raw/masks): reorientation, crop box, brain mask of scenario 3v2."""
import glob
import os

import numpy as np
import nibabel as nib
import pytest

import _pet_phantom as P
from modalities.pet import petutils as pu
from modalities.pet import manualmask, studyT, transfo, coreg2anat

opj = os.path.join


def _brain_mask(A, shape, R, t):
    """Brain (CSF, GM, WM, cerebellum) of the phantom on the PET grid A: x_acpc = R x_pet + t."""
    p = P.world_points(shape, A)
    lab = P.tissues(p @ R.T + t).reshape(shape)
    return np.isin(lab, (3, 4, 5, 6))


def _save_mask(m, A, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    img = nib.Nifti1Image(m.astype(np.uint8), A)
    img.set_sform(A, code=1)
    img.set_qform(A, code=1)
    nib.save(img, path)
    return path


def test_header_only_reorientation_keeps_the_voxels(tmp_path):
    A = P.grid(20, 2.0)
    raw = P.save(np.random.rand(20, 20, 20).astype(np.float32), A, str(tmp_path / 'raw_pet.nii.gz'))
    m = np.zeros((20, 20, 20), bool)
    m[5:12, 6:15, 8:13] = True
    mask = _save_mask(m, A, str(tmp_path / 'raw' / 'masks' / 's_desc-manual_mask.nii.gz'))
    # relabelled axes (3drefit -orient): same voxel array, another header
    A2 = A.copy()
    A2[:3, :3] = A[:3, :3] @ np.diag([1, -1, -1])
    reo = P.save(np.asarray(nib.load(raw).dataobj), A2, str(tmp_path / 'reo_pet.nii.gz'))
    pu.write_json(reo, raw, 'axes relabelled RAI -> LPI (3drefit -orient)', WorldTransform=(A2 @ np.linalg.inv(A)).tolist())
    out = manualmask.to_orig(mask, raw, reo, str(tmp_path / 'orig.nii.gz'), str(tmp_path / 'diary.txt'))
    o = nib.load(out)
    assert np.array_equal(np.asarray(o.dataobj) > 0, m) and np.allclose(o.affine, A2)


def test_resampled_reorientation_follows_the_world_transform(tmp_path):
    A = P.grid(24, 2.0)
    raw = P.save(np.random.rand(24, 24, 24).astype(np.float32), A, str(tmp_path / 'raw_pet.nii.gz'))
    m = np.zeros((24, 24, 24), bool)
    m[4:10, 8:20, 10:16] = True
    mask = _save_mask(m, A, str(tmp_path / 'm.nii.gz'))
    # step 1 resampled the image (3dWarp) on another grid, then moved it by W (90 degrees about z + shift)
    W = np.eye(4)
    W[:3, :3] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    W[:3, 3] = [4.0, -6.0, 2.0]
    A_reo = P.grid(30, 2.0, center=(4.0, -6.0, 2.0))
    reo = P.save(np.zeros((30, 30, 30), np.float32), A_reo, str(tmp_path / 'reo_pet.nii.gz'))
    pu.write_json(reo, raw, 'obliquity corrected by resampling (3dWarp -deoblique -linear)', WorldTransform=W.tolist())
    o = nib.load(manualmask.to_orig(mask, raw, reo, str(tmp_path / 'orig.nii.gz'), str(tmp_path / 'diary.txt')))
    got = np.asarray(o.dataobj) > 0
    # centre of the mask in world coordinates moved by W
    c_raw = (A @ np.r_[np.argwhere(m).mean(0), 1.0])[:3]
    c_out = (A_reo @ np.r_[np.argwhere(got).mean(0), 1.0])[:3]
    assert np.linalg.norm((W @ np.r_[c_raw, 1.0])[:3] - c_out) < 1.0
    assert abs(int(got.sum()) - int(m.sum())) <= 0.05 * m.sum()


def test_box_and_dilation(tmp_path):
    A = P.grid(30, 2.0)
    m = np.zeros((30, 30, 30), bool)
    m[10:15, 12:18, 14:16] = True
    mask = _save_mask(m, A, str(tmp_path / 'm.nii.gz'))
    img = P.save(np.ones((30, 30, 30), np.float32), A, str(tmp_path / 'img.nii.gz'))
    b = np.asarray(nib.load(manualmask.box(mask, img, str(tmp_path / 'box.nii.gz'), margin_mm=4.0)).dataobj) > 0
    idx = np.argwhere(b)
    assert idx.min(0).tolist() == [8, 10, 12] and idx.max(0).tolist() == [16, 19, 17]
    d = np.asarray(nib.load(manualmask.dilated(mask, 2.0, str(tmp_path / 'd.nii.gz'))).dataobj) > 0
    assert d.sum() > m.sum() and d[m].all()


def test_masked_image_keeps_the_world_coordinates(tmp_path):
    A = P.grid(30, 2.0)
    d = np.random.rand(30, 30, 30).astype(np.float32) + 1.0
    img = P.save(d, A, str(tmp_path / 'img.nii.gz'))
    m = np.zeros((30, 30, 30), bool)
    m[10:15, 12:18, 14:16] = True
    mask = _save_mask(m, A, str(tmp_path / 'm.nii.gz'))
    o = nib.load(manualmask.masked_image(mask, img, str(tmp_path / 'brain.nii.gz'), dilate_mm=0.0, margin_mm=4.0))
    od = np.asarray(o.dataobj)
    assert od.shape == (9, 10, 6)                                   # box of the mask + 2 voxels on each side
    # voxel (i, j, k) of the crop is voxel (i, j, k) + (8, 10, 12) of the image, at the same world position
    assert np.allclose(o.affine @ np.r_[0, 0, 0, 1.0], A @ np.r_[8, 10, 12, 1.0])
    assert np.allclose(od[2:7, 2:8, 2:4], d[10:15, 12:18, 14:16]) and od.sum() == pytest.approx(d[m].sum(), rel=1e-5)
    assert pu.read_json(pu.sidecar(str(tmp_path / 'brain.nii.gz')))['BrainMasked'] is True


@pytest.mark.slow
def test_manual_masks_give_the_brain_and_the_crop_in_3v2(tmp_path, monkeypatch):
    """Scenario 3v2 (no anatomy): with the masks in pet/raw/masks the study template is built on the brain of the
    session means, crop_for_MoCo='mask' replaces the freeview crop, the bratio is computed over the manual brain (it
    replaces the brain of the reference template brought by the bridge), and the brain template is bridged to the
    reference template."""
    from test_pet import _study
    from modalities.pet import _0_Pipeline_launcher as L
    from Tools import getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch, subjects=('01', '02'), anat=False)
    for sub in ('01', '02'):
        shape = nib.load(s[sub]['raw']).shape[:3]
        _save_mask(_brain_mask(s[sub]['A'], shape, s[sub]['R'], s[sub]['t']), s[sub]['A'],
                   opj(s[sub]['dp'], 'pet', 'raw', 'masks', f'sub-{sub}_ses-1_trc-FDG_desc-manual_mask.nii.gz'))
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1'), ('02', '1')]), '*_trc-FDG_*pet.nii.gz',
                     'FDG', [''], [''], '', 'T1w', 'T1w', False, False, IhaveanANAT=False,
                     creat_study_template_pet=True, studyT_transform='TRSAA', crop_for_MoCo='mask', bratio=True,
                     selected_atlases=[['TEST', 'all']], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    st = studyT.template_name(opj(bids, 'sty_template', 'pet'), 'FDG')
    assert os.path.isfile(st) and not os.path.isfile(studyT.crop_name(st))      # step 3 (freeview crop) skipped
    assert studyT.is_brain_masked(st)                                          # template built on the brain
    for sub in ('01', '02'):
        paths = getpath.pet(s[sub]['dp'], 'EDNiX')
        session = f'sub-{sub}_ses-1_trc-FDG'
        rec = transfo.read_record(coreg2anat.record_path(paths, session))
        assert rec['scenario'] == '3v2' and sorted(rec['targets']) == ['EDNiX', 'studyTemplate']
        assert rec['pet_space']['atlases']                                      # atlases through the bridge
        assert os.path.isfile(opj(paths['orig_formc'], session + '_desc-reorient-mean-brain_pet.nii.gz'))
        assert os.path.isfile(opj(paths['orig_formc'], session + '_desc-MoCo-mean-brain_pet.nii.gz'))
        orig = manualmask.orig_file(paths, session)
        assert os.path.isfile(orig)
        assert os.path.isfile(opj(paths['orig_formc'], session + '_desc-cropbox_mask.nii.gz'))
        br = sorted(glob.glob(opj(paths['orig_post'], session + '*_desc-bratio*_pet.nii.gz')))
        assert br, 'bratio not computed with the manual mask'
        b = nib.load(br[0])
        bd = np.asarray(b.dataobj, dtype=np.float32)
        bd = bd[..., -1] if bd.ndim == 4 else bd
        m = pu.mask_bool(manualmask.on_grid(orig, br[0], str(tmp_path / f'chk{sub}.nii.gz')))
        assert float(bd[m].mean()) == pytest.approx(1.0, abs=0.05)
