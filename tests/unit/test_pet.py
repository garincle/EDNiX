"""PET / CT pipelines: helpers, transforms, motion correction, registration and an end-to-end run on a phantom."""
import os
import json
import math

import numpy as np
import nibabel as nib
import pandas as pd
import pytest

import ants
import _pet_phantom as P
from modalities.pet import petutils as pu
from modalities.pet import transfo, dftFile, cut_image, run_cut, MoCo, coreg2anat, loop1, loop2, studyT
from modalities.pet.turkutools._nifti import TurkuWorkspace
from atlases import getlabels

opj = os.path.join
opb = os.path.basename

# reproducible registrations (ANTs random sampling, one thread): the pose tolerances below are not flaky
if hasattr(ants, 'config'):                    # antspyx >= 0.6
    ants.config.set_ants_deterministic(True, 123)
else:                                          # antspyx 0.5: same settings by hand
    os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = '1'
    np.random.seed(123)


# ------------------------------------------------------------------------------------------------ helpers
def test_isotopes_and_units():
    for code in ('F18', '18F', 'F-18', '^18^Fluorine', 'Fluorine-18'):
        assert pu.normalise_isotope(code) == 'F-18'
    assert pu.normalise_isotope('C11') == 'C-11'
    assert pu.isotope_from_half_life(6586.0) == 'F-18' and pu.isotope_from_half_life(1223.0) == 'C-11'
    assert pu.isotope_from_half_life(3000.0) is None
    assert pu.concentration_to_bq_per_ml('kBq/mL') == 1000.0 and pu.concentration_to_bq_per_ml('Bq/cc') == 1.0
    assert pu.concentration_to_bq_per_ml('nCi/mL') == pytest.approx(37.0)
    assert pu.activity_to_mbq(5, 'mCi') == pytest.approx(185.0)
    assert pu.activity_to_mbq(2e8, 'Bq') == pytest.approx(200.0)


def test_suv_factor_and_decay_reference():
    info = dict(Units='kBq/mL', ImageDecayCorrected=True, ImageDecayCorrectionTime=0, InjectionStart=0)
    f, rep = pu.suv_scale(info, dict(weight_kg=10.0, dose_MBq=100.0), 'F-18')
    # SUV = C(kBq/mL) * 1000 * 10 000 g / 1e8 Bq
    assert f == pytest.approx(1000 * 10000 / 1e8)
    # image decay-corrected to 10 min after injection: dose decayed to that time
    info2 = dict(info, ImageDecayCorrectionTime=600, InjectionStart=0)
    f2, _ = pu.suv_scale(info2, dict(weight_kg=10.0, dose_MBq=100.0), 'F-18')
    assert f2 == pytest.approx(f * math.exp(math.log(2) * 600 / pu.HALF_LIFE_S['F-18']))
    assert pu.suv_scale(info, dict(weight_kg=None, dose_MBq=100.0), 'F-18')[0] is None
    assert pu.suv_scale(dict(Units='SUV'), dict(weight_kg=None, dose_MBq=None), 'F-18')[0] == 1.0


def test_decay_correction_of_frames_restores_a_constant_activity():
    lam = pu.decay_constant('C-11')
    t = dict(start=np.array([0., 60., 300.]), duration=np.array([60., 240., 600.]))
    # frame average of A0 exp(-lam t) over [s, s + d]
    mean_decayed = np.array([(np.exp(-lam * s) - np.exp(-lam * (s + d))) / (lam * d) for s, d in zip(t['start'], t['duration'])])
    d4 = np.ones((2, 2, 2, 3), np.float32) * mean_decayed.astype(np.float32)
    out = pu.decay_correct_frames(d4, t, 'C-11', 0.0)
    assert np.allclose(out, 1.0, rtol=1e-5)


def test_frame_timing_window_and_weighted_mean():
    info = dict(FrameTimesStart=[0, 60, 120, 600], FrameDuration=[60, 60, 480, 600])
    t = pu.frame_timing(info, 4)
    assert np.allclose(t['mid'], [30, 90, 360, 900])
    assert list(pu.window_frames(t, [5, 20])) == [2, 3]
    d = np.stack([np.full((2, 2, 2), v, np.float32) for v in (1, 2, 3, 4)], -1)
    m = pu.weighted_mean(d, t, np.array([2, 3]))
    assert np.allclose(m, (3 * 480 + 4 * 600) / 1080)
    with pytest.raises(ValueError):
        pu.frame_timing(info, 5)


def test_names():
    raw = '/x/sub-01_ses-1_trc-FDG_run-2_pet.nii.gz'
    assert pu.base_name(raw) == 'sub-01_ses-1_trc-FDG_run-2'
    assert pu.base_name('/x/sub-01_ses-1_trc-FDG_run-2_desc-MoCo-within_pet.nii.gz') == 'sub-01_ses-1_trc-FDG_run-2'
    assert pu.session_base('sub-01_ses-1_trc-FDG_run-2') == 'sub-01_ses-1_trc-FDG'
    assert pu.out_name('/o', 'b', desc='SUV', space='EDNiX').endswith('/o/b_space-EDNiX_desc-SUV_pet.nii.gz')


def test_subject_values_from_bids_tables(tmp_path):
    (tmp_path / 'participants.tsv').write_text('participant_id\tbody_weight\tblood_concentration\nsub-01\t7.5\t0.9008\n')
    v = pu.subject_values(str(tmp_path), '01', '1', dict(InjectedRadioactivity=3.0, InjectedRadioactivityUnits='mCi'))
    assert v['weight_kg'] == 7.5 and v['dose_MBq'] == pytest.approx(111.0)
    assert v['glucose_mmol_L'] == pytest.approx(5.0, rel=1e-3)


def test_dft_roundtrip(tmp_path):
    f = dftFile.write(str(tmp_path / 'a.dft'), [0, 60], [60, 180], [[1.0, 2.0], [3.0, np.nan]], ['cer', 'put'],
                      [100.0, 50.0], 'Bq/mL', 'study')
    lines = open(f).read().splitlines()
    assert lines[0] == 'DFT\tcer\tput' and lines[3].startswith('Times (sec)\t100') and lines[5].endswith('.')
    s, e, v, n, vol = dftFile.read(f)
    assert list(n) == ['cer', 'put'] and np.allclose(vol, [100, 50]) and np.allclose(e, [60, 180])
    assert v[0, 1] == 2.0 and np.isnan(v[1, 1])


def test_reference_tac_is_the_mean_over_the_union_of_the_labels(tmp_path):
    lab = np.zeros((4, 4, 1), int)
    lab[:3, :, 0] = 7        # 12 voxels
    lab[3, 0, 0] = 1007      # 1 voxel
    data = np.zeros((4, 4, 1, 2), np.float32)
    data[lab == 7] = [1.0, 1.0]
    data[lab == 1007] = [14.0, 14.0]
    t = dict(start=np.array([0., 60.]), end=np.array([60., 120.]))
    f, n = dftFile.ref(str(tmp_path / 'r.dft'), data, lab, [7, 1007], 1.0, t, 'Bq/mL', 'Cer')
    _, _, v, _, vol = dftFile.read(f)
    assert n == 13 and np.allclose(v[:, 0], 2.0) and vol[0] == 13     # (12 x 1 + 14) / 13, not (1 + 14) / 2


def test_time_blocks_and_fmri_overlap():
    t = pu.frame_timing(dict(FrameTimesStart=[0, 300, 600, 900, 1500], FrameDuration=[300, 300, 300, 600, 600]), 5)
    blocks = cut_image.setTime(t, 10)
    assert [list(b[2]) for b in blocks] == [[0, 1], [2], [3], [4]]
    func = dict(AcquisitionTime='10:05:00.000000', RepetitionTime=2.0)
    ov, (a, b) = run_cut.getframes(dict(TimeZero='10:00:00'), t, func, 160, 10)
    assert (a, b) == (320.0, 620.0)
    assert np.allclose(ov, [0, 280, 20, 0, 0])


def test_reference_region_resolution(tmp_path):
    df = pd.DataFrame(dict(label_id=[7, 1007, 10], region_name=['L_Cerebellum', 'R_Cerebellum', 'L_Cortex']))
    atlas = dict(name='EDNIxCSCLR', level=0)
    r = getlabels.resolve(dict(type='atlas', labels='Cerebell', name='Cer'), atlas, [7, 1007, 10], df, None)
    assert r['labels'] == [7, 1007] and r['atlas'] == 'EDNIxCSCLR' and r['level'] == 0
    assert getlabels.resolve(dict(type='atlas', labels='default'), atlas, [7], df, None)['labels'] == [7, 1007]
    saved = tmp_path / 'ref.json'
    saved.write_text(json.dumps(dict(type='mask', name='white')))
    assert getlabels.resolve('unknown', atlas, [7], df, str(saved)) == dict(type='mask', name='white')
    with pytest.raises(ValueError):
        getlabels.resolve(dict(type='mask', name='bone'), atlas, [], df, None)
    assert getlabels.cerebellum(dict(name='BAL2024')) == [29, 30]


def test_scenarios():
    assert coreg2anat.scenario(True, True, False, None, '', False) == '1v1'
    assert coreg2anat.scenario(True, False, False, None, '', False) == '1v2a'
    assert coreg2anat.scenario(True, False, True, ['acpc'], '', False) == '1v2b'
    assert coreg2anat.scenario(False, False, True, ['EDNiX'], '', False) == '2v1'
    assert coreg2anat.scenario(False, False, True, ['studyTemplate'], '', False) == '2v2'
    assert coreg2anat.scenario(False, False, False, None, 'tpl.nii.gz', False) == '3v1'
    assert coreg2anat.scenario(False, False, False, None, '', True) == '3v2'
    with pytest.raises(ValueError):
        coreg2anat.scenario(False, False, False, None, '', False)


def test_turku_outputs_get_the_input_geometry_back(tmp_path):
    A = P.grid(6, 2.0)
    ref = P.save(np.random.rand(6, 6, 6, 3), A, str(tmp_path / 'in.nii.gz'))
    flipped = np.diag([-2.0, 2.0, 2.0, 1.0])
    flipped[:3, 3] = A[:3, 3] + [10.0, 0, 0]
    data = np.arange(216, dtype=np.float32).reshape(6, 6, 6)
    nib.save(nib.Nifti1Image(data[::-1], flipped), str(tmp_path / 'out.nii'))
    out = TurkuWorkspace.back(str(tmp_path / 'out.nii'), ref, str(tmp_path / 'res.nii.gz'))
    assert np.allclose(nib.load(out).affine, A) and np.allclose(nib.load(out).get_fdata(), data)


# ------------------------------------------------------------------------------------------------ transforms
def test_chain_composition_matches_the_direct_map(tmp_path):
    # M -> A: x_A -> x_M = R1 x_A + t1 ; A -> T: x_T -> x_A = R2 x_T + t2
    R1, t1 = P.pose((5, 0, 10), (3, -2, 1))
    R2, t2 = P.pose((0, -7, 4), (-4, 1, 2))
    m1 = P.write_affine_mat(str(tmp_path / 'm1.mat'), R1, t1)
    m2 = P.write_affine_mat(str(tmp_path / 'm2.mat'), R2, t2)
    chain = transfo.compose(transfo.affine(m1), transfo.affine(m2))
    m12 = P.write_affine_mat(str(tmp_path / 'm12.mat'), R1 @ R2, R1 @ t2 + t1)
    A = P.grid(40, 2.0)
    img = P.save(P.pet_frame((40,) * 3, A, np.eye(3), np.zeros(3)), A, str(tmp_path / 'img.nii.gz'))
    a = transfo.apply(img, img, chain, 'to').numpy()
    b = transfo.apply(img, img, transfo.affine(m12), 'to').numpy()
    assert np.corrcoef(a.ravel(), b.ravel())[0, 1] > 0.999
    back = transfo.apply(img, transfo.apply(img, img, chain, 'to'), chain, 'from').numpy()
    orig = nib.load(img).get_fdata()
    inner = orig > 0
    assert np.corrcoef(back[inner], orig[inner])[0, 1] > 0.97


# ------------------------------------------------------------------------------------------------ motion / registration
@pytest.mark.slow
def test_within_run_motion_correction_recovers_the_motion(tmp_path):
    R, t = P.pose((0, 0, 0), (0, 0, 0))
    # frames 0-1 (not registered) are in the position of frame 2, the assumption of the correction
    moves = [((4, 0, 6), (3, -2, 1))] * 3 + [((-5, 3, 0), (-4, 2, 2)), ((0, -6, 4), (2, 4, -3)),
                                            ((3, 3, -5), (-3, -3, 2)), ((6, 0, 0), (0, 3, 3)), ((-4, -4, 4), (4, 0, -2))]
    motions = [(P.pose(r, m)[0], np.asarray(m, float)) for r, m in moves]
    img, A, info = P.dynamic_pet(str(tmp_path / 'sub-01_ses-1_trc-FDG_pet.nii.gz'), R, t, motions=motions)
    timing = pu.frame_timing(info, 8)
    out = MoCo.dynamic(img, str(tmp_path / 'pre'), str(tmp_path / 'formc'), timing, 2, None, False,
                       str(tmp_path / 'diary.txt'), '')
    mc = nib.load(out).get_fdata()
    raw = nib.load(img).get_fdata()

    def r(a, b):
        return np.corrcoef((a / a.max()).ravel(), (b / b.max()).ravel())[0, 1]
    pairs = [(2, 7), (3, 6), (4, 5), (2, 3)]
    r_raw = min(r(raw[..., i], raw[..., j]) for i, j in pairs)
    r_mc = min(r(mc[..., i], mc[..., j]) for i, j in pairs)
    assert r_mc > 0.995 and r_raw < 0.97, (r_raw, r_mc)
    # frames 0-1 (not registered) received the transform of frame 2
    assert r(mc[..., 1], mc[..., 2]) > 0.99
    tab = pd.read_csv(opj(str(tmp_path / 'formc'), 'sub-01_ses-1_trc-FDG_desc-MoCo-within_motion.tsv'), sep='\t')
    assert len(tab) == 8 and tab['mean_displacement_mm'].max() > 2.0


@pytest.mark.slow
def test_rigid_pet_to_anatomy_registration_recovers_the_pose(tmp_path):
    import ants
    files, _, _ = P.anat_derivative(str(tmp_path / 'bids'))
    R, t = P.pose()
    A = P.grid(40, 2.0)
    pet = P.save(P.pet_frame((40,) * 3, A, R, t), A, str(tmp_path / 'pet.nii.gz'))
    chain = coreg2anat.rigid_to(files['acpc'], pet, str(tmp_path / 'r_'), 'mattes', mask_file=files['mask'])
    tr = ants.read_transform(chain['to'][0][0])
    err = []
    for q in ([0, 0, 0], [15, 10, -5], [-12, -20, 10], [5, -25, -12]):
        q = np.asarray(q, float)
        p_true = R.T @ (q - t)                            # acpc point -> PET point (RAS)
        p_est = P.F @ np.asarray(tr.apply_to_point(tuple(P.F @ q)))
        err.append(np.linalg.norm(p_true - p_est))
    assert max(err) < 1.0, err                           # < half a PET voxel


def _mstart_error(chain, R, t):
    """Largest error (mm) of a PET -> AC-PC chain at points of the brain of the phantom."""
    trs = [ants.read_transform(f) for f, _ in chain['to']]
    err = []
    for q in ([0, 0, 0], [15, 10, -5], [-12, -20, 10], [5, -25, -12]):
        q = np.asarray(q, float)
        p = tuple(P.F @ q)
        for tr in trs:
            p = tr.apply_to_point(p)
        err.append(np.linalg.norm(R.T @ (q - t) - P.F @ np.asarray(p)))
    return max(err)


def _render(shape, A, R, t, lut, extra=(), smooth_mm=3.0, ss=2):
    """Phantom on the grid A (pose x_acpc = R x + t), tissue values `lut`, plus ellipsoids (centre, radii, value) in
    AC-PC coordinates painted over the head (extracerebral uptake, torso, bed)."""
    from scipy.ndimage import gaussian_filter
    Af = A.copy()
    Af[:3, :3] = A[:3, :3] / ss
    Af[:3, 3] = A[:3, 3] - (ss - 1) / (2.0 * ss) * np.diag(A[:3, :3])
    fine = tuple(n * ss for n in shape)
    q = P.world_points(fine, Af) @ R.T + t
    v = np.asarray(lut, np.float32)[P.tissues(q)]
    for c, r, val in extra:
        v[P._ell(q, c, r)] = val
    v = v.reshape(fine).reshape(shape[0], ss, shape[1], ss, shape[2], ss).mean(axis=(1, 3, 5))
    return gaussian_filter(v, smooth_mm / A[0, 0]).astype(np.float32) if smooth_mm else v.astype(np.float32)


@pytest.mark.slow
def test_multistart_registration_of_a_pet_with_the_torso_in_the_field_of_view(tmp_path):
    # 160 mm field of view: head + torso + heart (FDG); the centre-of-mass start of ANTs is pulled to the torso
    files, _, _ = P.anat_derivative(str(tmp_path / 'bids'))
    torso = [((0, -10, -80), (40, 35, 40), 900.0), ((0, 0, -85), (14, 14, 14), 9000.0)]
    A = P.grid(80, 2.0, center=(0, 0, -30))
    for rot, tr in (((8, -5, 12), (6, -9, 4)), ((5, 5, -25), (-40, 30, 70))):
        R, t = P.pose(rot, tr)
        pet = P.save(_render((80,) * 3, A, R, t, [P.PET[k] for k in range(7)], torso), A,
                     str(tmp_path / 'pet.nii.gz'))
        prefix = str(tmp_path / 'r_')
        chain = coreg2anat.rigid_to(files['acpc'], pet, prefix, 'mattes', mask_file=files['mask'])
        assert _mstart_error(chain, R, t) < 1.0, (rot, tr)
        with open(prefix + 'candidates.json') as f:
            starts = json.load(f)
        assert {s['init'] for s in starts['starts']} >= {'com', 'fov'}
        # the kept start has the best mutual information among the starts that cover the brain
        ok = [s for s in starts['starts'] if s['coverage'] >= 0.9]
        assert starts['kept'] == min(ok, key=lambda s: s['metric'])['init']
    # the AC-PC matrix of the anatomy is a start when given
    cands = coreg2anat.init_candidates(ants.image_read(files['acpc']), None, ants.image_read(pet),
                                       str(tmp_path / 'c_'), 'auto', opj(files['dir_transfo'], 'acpc_0GenericAffine.mat'))
    assert [c[0] for c in cands] == ['com', 'fov', 'acpc', 'identity', 'search']   # 'brain' needs both brains
    with pytest.raises(ValueError):
        coreg2anat.init_candidates(ants.image_read(files['acpc']), None, ants.image_read(pet), str(tmp_path / 'c_'),
                                   'bad')


@pytest.mark.slow
def test_ct_to_anatomy_centre_of_mass_and_bed(tmp_path):
    # CT in HU (air -1000) with the bed under the head, off-centre field of view: the centre of mass of the raw CT
    # is far from the head; the CT shifted to >= 0 and restricted to the head is registered within half a voxel
    files, _, _ = P.anat_derivative(str(tmp_path / 'bids'))
    hu = [-1000.0, 40.0, 1200.0, 5.0, 38.0, 30.0, 38.0]
    bed = [((0, -45, -40), (60, 4, 80), 400.0)]
    R, t = P.pose((10, -6, 8), (12, 20, -15))
    A = P.grid(70, 2.0, center=(20, -25, 10))
    ct = _render((70,) * 3, A, R, t, hu, bed, smooth_mm=1.0)
    ctf = P.save(ct, A, str(tmp_path / 'ct.nii.gz'))
    img = ants.image_read(ctf)
    head = coreg2anat.ct_head_mask(img)
    q = P.world_points((70,) * 3, A) @ R.T + t
    in_bed = P._ell(q, *bed[0][:2]).reshape((70,) * 3)
    assert head.numpy()[in_bed].mean() < 0.05                        # the bed is not in the head mask
    com_raw = np.asarray(ants.get_center_of_mass(img))
    com_pos = np.asarray(ants.get_center_of_mass(coreg2anat._nonnegative(img) * head))
    head_c = P.F @ (R.T @ (np.zeros(3) - t))                          # centre of the head in the CT (LPS)
    assert np.linalg.norm(com_raw - head_c) > 2 * np.linalg.norm(com_pos - head_c)
    chain = coreg2anat.rigid_to(files['acpc'], ctf, str(tmp_path / 'ct_'), 'mattes', mask_file=files['mask'],
                                dilate_mm=10.0, moving_mask=head)
    assert _mstart_error(chain, R, t) < 1.0


@pytest.mark.slow
def test_ct_large_rotation_search_and_brain_window_pass(tmp_path):
    # low-dose CT (noise 25 HU) of a head flexed by 80 degrees, with the neck, the cervical spine, the head holder
    # and the bed: the centre-of-mass start fails, the search over the rotations finds it; the second pass in the
    # brain window keeps it within half a voxel
    from modalities.ct import coreg2anat as ctc
    files, _, _ = P.anat_derivative(str(tmp_path / 'bids'))
    hu = [-1000.0, 40.0, 1200.0, 5.0, 38.0, 30.0, 38.0]
    extra = [((0, -5, -55), (18, 18, 30), 40.0), ((0, -12, -50), (5, 5, 28), 700.0),
             ((0, -38, -5), (34, 4, 36), 150.0), ((0, -48, -40), (60, 4, 80), 400.0)]
    A = P.grid(80, 2.0, center=(15, -20, -10))
    R, t = P.pose((80, 0, 0), (12, 20, -15))
    ct = _render((80,) * 3, A, R, t, hu, extra, smooth_mm=1.0)
    ct = ct + np.random.default_rng(0).normal(0, 25, ct.shape).astype(np.float32)
    ctf = P.save(ct, A, str(tmp_path / 'ct.nii.gz'))
    head = coreg2anat.ct_head_mask(ants.image_read(ctf))
    prefix = str(tmp_path / 'ct_')
    chain = coreg2anat.rigid_to(files['acpc'], ctf, prefix, 'mattes', mask_file=files['mask'], dilate_mm=10.0,
                                moving_mask=head)
    with open(prefix + 'candidates.json') as f:
        starts = {s['init']: s for s in json.load(f)['starts']}
    assert starts['search']['metric'] < starts['com']['metric'] and _mstart_error(chain, R, t) < 1.0
    tr = ants.read_transform(str(tmp_path / 'ct_init-search.mat'))
    lin = np.asarray(tr.parameters[:9]).reshape(3, 3)
    assert np.allclose(lin @ lin.T, np.eye(3), atol=1e-4)                # the start is a rotation (no scale)
    win = ctc.brain_window(ctf, str(tmp_path / 'win.nii.gz'))
    w = nib.load(win).get_fdata()
    lab = P.tissues(P.world_points((80,) * 3, A) @ R.T + t).reshape((80,) * 3)
    assert w[lab == 4].mean() > w[lab == 5].mean() > w[lab == 3].mean()  # grey > white > CSF
    assert w[lab == 2].mean() < 0.1 and w[lab == 0].mean() < 0.05       # bone and air dark
    second = coreg2anat.refine_in_brain(files['acpc'], files['mask'], win, chain, str(tmp_path / 'b_'), 'mattes',
                                        3.0, 'Rigid', str(tmp_path / 'bm.nii.gz'))
    assert _mstart_error(second, R, t) < 1.0


@pytest.mark.slow
def test_ct_with_neck_and_shoulders_intracranial_start(tmp_path):
    # CT field of view down to the shoulders (CERMEP-like): the centres of mass and of the field of view are pulled
    # to the body and those starts converge with half of the brain outside; the intracranial cavity (skull ring
    # filled slice by slice) gives the 'brain' start
    files, _, _ = P.anat_derivative(str(tmp_path / 'bids'))
    hu = [-1000.0, 40.0, 1200.0, 5.0, 38.0, 30.0, 38.0]
    body = [((0, -5, -60), (18, 18, 30), 40.0), ((0, -12, -55), (5, 5, 30), 700.0),
            ((0, -5, -115), (70, 30, 35), 40.0), ((0, -38, -5), (34, 4, 36), 150.0),
            ((0, -48, -60), (80, 4, 100), 400.0)]
    A = P.grid(96, 2.0, center=(10, -15, -45))
    R, t = P.pose((-25, 10, 0), (12, 20, -15))
    ct = _render((96,) * 3, A, R, t, hu, body, smooth_mm=1.0)
    ct = ct + np.random.default_rng(0).normal(0, 25, ct.shape).astype(np.float32)
    ctf = P.save(ct, A, str(tmp_path / 'ct.nii.gz'))
    img = ants.image_read(ctf)
    cav = coreg2anat.ct_intracranial_mask(img)
    lab = P.tissues(P.world_points((96,) * 3, A) @ R.T + t).reshape((96,) * 3)
    c = cav.numpy() > 0
    assert c[lab <= 2].mean() < 0.01 and c[lab >= 3].sum() > 1000            # inside the skull only
    prefix = str(tmp_path / 'ct_')
    chain = coreg2anat.rigid_to(files['acpc'], ctf, prefix, 'mattes', mask_file=files['mask'], dilate_mm=10.0,
                                moving_mask=coreg2anat.ct_head_mask(img), moving_brain=cav)
    with open(prefix + 'candidates.json') as f:
        starts = {s['init']: s for s in json.load(f)['starts']}
    assert starts['com']['coverage'] < 0.9 and starts['brain']['coverage'] >= 0.9
    assert _mstart_error(chain, R, t) < 1.0
    # a start that leaves the CT outside the mask of the anatomy (AC-PC matrix of other scanner coordinates) fails
    # with an explicit reason instead of the ITK 'zero-valued spacing' error
    far = coreg2anat._translation_file(str(tmp_path / 'far.mat'), (400.0, 0.0, 0.0))
    with pytest.raises(RuntimeError, match='no overlap'):
        coreg2anat.rigid_to(files['acpc'], ctf, str(tmp_path / 'f_'), 'mattes', mask_file=files['mask'],
                            dilate_mm=10.0, moving_mask=coreg2anat.ct_head_mask(img), init='acpc', acpc=far)
    assert coreg2anat.ct_intracranial_mask(img * 0 - 1000.0) is None          # no skull: no cavity


@pytest.mark.slow
def test_pet_masks_second_pass_and_synthetic_target(tmp_path):
    # head PET with hot extracerebral structures next to the cortex (temporal muscle, glands, neck)
    files, _, _ = P.anat_derivative(str(tmp_path / 'bids'))
    hot = [((29, 0, -6), (5, 16, 14), 4500.0), ((-12, 33, -10), (6, 5, 5), 6000.0),
           ((12, 33, -10), (6, 5, 5), 6000.0), ((0, -20, -33), (14, 8, 6), 3500.0)]
    R, t = P.pose((-6, 10, -8), (-5, 7, -3))
    A = P.grid(40, 2.5)
    pet = P.save(_render((40,) * 3, A, R, t, [P.PET[k] for k in range(7)], hot, smooth_mm=2.1), A,
                 str(tmp_path / 'pet.nii.gz'))
    # threshold mask: one component, covers the brain, excludes the air
    m = coreg2anat.threshold_mask(pet, 50.0, 0.0).numpy() > 0
    truth = P.tissues(P.world_points((40,) * 3, A) @ R.T + t).reshape((40,) * 3)
    assert m[truth >= 4].mean() > 0.9 and m[truth == 0].mean() < 0.1
    out = coreg2anat.pet_mask('threshold', pet, str(tmp_path / 'm.nii.gz'), 50.0, 2.0)
    assert os.path.isfile(out) and coreg2anat.pet_mask('none', pet, '') is None
    with pytest.raises(ValueError):
        coreg2anat.pet_mask('bad', pet, '')
    first = coreg2anat.rigid_to(files['acpc'], pet, str(tmp_path / 'p1_'), 'mattes', mask_file=files['mask'])
    second = coreg2anat.refine_in_brain(files['acpc'], files['mask'], pet, first, str(tmp_path / 'p2_'), 'mattes',
                                        3.0, 'Rigid', str(tmp_path / 'brain_pet.nii.gz'))
    assert _mstart_error(second, R, t) < 1.0
    bm = nib.load(str(tmp_path / 'brain_pet.nii.gz')).get_fdata() > 0
    assert bm[truth >= 3].mean() > 0.95 and bm[truth == 0].mean() < 0.05
    # FDG-like target from the tissue masks: grey > white > CSF, registration within half a PET voxel
    syn = coreg2anat.synthetic_pet(files, (2.5, 2.5, 2.5), str(tmp_path / 'syn.nii.gz'))
    s = nib.load(syn).get_fdata()
    tt = P.tissues(P.world_points(s.shape, nib.load(syn).affine)).reshape(s.shape)
    assert s[tt == 4].mean() > s[tt == 5].mean() > s[tt == 0].mean()
    chain = coreg2anat.rigid_to(syn, pet, str(tmp_path / 's_'), 'mattes', mask_file=files['mask'])
    assert _mstart_error(chain, R, t) < 1.25


# ------------------------------------------------------------------------------------------------ end to end
@pytest.mark.slow
def test_pipeline_1v2a_end_to_end_on_a_phantom(tmp_path, monkeypatch):
    P.patch_run_cmd(monkeypatch)
    bids = str(tmp_path / 'bids')
    files, package, dp = P.anat_derivative(bids)
    (tmp_path / 'bids' / 'participants.tsv').write_text('participant_id\tbody_weight\nsub-01\t10\n')
    R, t = P.pose()
    motions = [None] * 4 + [(P.pose((0, 0, 2), (1, 0, 0))[0], np.array([1.0, 0, 0]))] * 4
    raw, A, info = P.dynamic_pet(opj(dp, 'pet', 'sub-01_ses-1_trc-FDG_pet.nii.gz'), R, t, motions=motions)
    from Tools import getpath
    paths = getpath.pet(dp, 'EDNiX')
    diary = str(tmp_path / 'diary.txt')
    Skip = []
    reor = loop1.run([raw], paths, [''], [''], '', files['acpc'], '', False, False, 'WARP', 'F-18', True, ('', ''),
                     Skip, diary, '', '', '')
    assert os.path.isfile(reor[0]) and os.path.isfile(loop1.sif_file(paths, raw))
    df = pd.DataFrame(dict(label_id=[1, 7, 10, 1010, 20, 30],
                           region_name=['Cerebrum', 'Cerebellum', 'L_Ctx', 'R_Ctx', 'WM', 'CSF']))
    opts = dict(Skip_step=Skip, T_PET='T1w', ct_pet_same_space=True, nb_void=2, crop_for_MoCo='none', redo_crop=False,
                MoCo_transform='Rigid', coreg_metric='mattes', coreg_transform='Rigid', type_of_transform_template='SyN', studyT_transform='SyN',
                do_manual=False, Need_check_MC=False, Need_check_coreg=False, SUV=True, bratio=True,
                SUV_window=[10, 40], Do_cut=True, duration=10, srtm=True,
                reference_signal=dict(type='atlas', atlas='TEST', level=0, labels=[7], name='Cer'), DV_time=10,
                kinetic_end=None, pvc_fwhm=None, pvc_method='RBV', do_pet_to_func=False, endfmri='*_bold.nii.gz',
                T1_eq=0, export_resolution='pet', export_dynamic=False, isotope='F-18')
    reg = dict(scenario='1v2a', anat=files, packages=[package], ct_record=None, template_PET='', reference='EDNiX',
               study_template='', study_crop='', template_files=None,
               target_files={'EDNiX': dict(label_base=files['template_label_base'], list_atlas=['TEST'],
                                           masks=[files['template_mask']])},
               list_atlas=['TEST'], label_tables={'TEST': df}, levels={'TEST': 'all'},
               ref_atlas=dict(name='TEST', level=0), reference_json=str(tmp_path / 'ref.json'))
    loop2.run('01', '1', dp, bids, paths, 'sub-01_ses-1_trc-FDG', reor, [raw], opts, reg, diary, '', '', '', '')

    rec = transfo.read_record(coreg2anat.record_path(paths, 'sub-01_ses-1_trc-FDG'))
    assert sorted(rec['targets']) == ['EDNiX', 'acpc'] and rec['scenario'] == '1v2a'
    # atlas in the PET space vs the true labels of the PET grid
    atlas_pet = nib.load(rec['pet_space']['atlases']['TEST']).get_fdata()[..., 0]
    q = P.world_points((40,) * 3, A) @ R.T + t
    truth = P.tissues(q).reshape((40,) * 3)
    dice = 2 * ((atlas_pet == 7) & (truth == 6)).sum() / ((atlas_pet == 7).sum() + (truth == 6).sum())
    assert dice > 0.8, dice
    # SUV of the white matter at the end of the scan: 1500 Bq/mL x (1 - exp(-t/10min)) x 10 kg / 100 MBq
    suv = nib.load(pu.out_name(paths['orig_post'], 'sub-01_ses-1_trc-FDG', desc='SUV-mean')).get_fdata()
    act = nib.load(pu.out_name(paths['orig_post'], 'sub-01_ses-1_trc-FDG', desc='activity-mean')).get_fdata()
    assert np.allclose(suv, act * 1e-4, rtol=1e-5, atol=1e-7)
    # brain ratio has mean 1 over the brain mask
    br = nib.load(pu.out_name(paths['orig_post'], 'sub-01_ses-1_trc-FDG', desc='bratio-mean')).get_fdata()
    brain = nib.load(rec['pet_space']['masks']['brain']).get_fdata() > 0.5
    assert br[brain].mean() == pytest.approx(1.0, rel=1e-4)
    # regional tables, TACs, DFT files and parametric images of the kinetic models
    roi = paths['orig_roi']
    tab = pd.read_csv(opj(roi, 'sub-01_ses-1_trc-FDG_seg-TEST-lvl1_desc-SUV-mean_ROIs.tsv'), sep='\t')
    gm = tab.loc[tab['region'] == 'L_Ctx', 'mean'].iloc[0]
    wm = tab.loc[tab['region'] == 'WM', 'mean'].iloc[0]
    assert gm > wm > 0
    assert os.path.isfile(opj(roi, 'sub-01_ses-1_trc-FDG_Cer_Kinetic.dft'))
    assert os.path.isfile(opj(roi, 'sub-01_ses-1_trc-FDG_seg-TEST-lvl0_desc-Kinetic_ROIs_bp.res'))
    for d in ('BP', 'SRTM', 'DVR'):
        assert os.path.isfile(pu.out_name(paths['orig_post'], 'sub-01_ses-1_trc-FDG', desc=d))
    # export: SUV in the template space, aligned with the template
    exp = pu.out_name(opj(paths['templates'], 'EDNiX'), 'sub-01_ses-1_trc-FDG', desc='bratio-mean', space='EDNiX')
    assert os.path.isfile(exp) and os.path.isfile(exp.replace('bratio-mean', 'bratio-mean-SS'))
    e = nib.load(exp)
    assert os.path.isfile(opj(paths['templates'], 'EDNiX', 'EDNiX_desc-grid-pet-resolution.nii.gz'))
    pts = P.world_points(e.shape, e.affine) @ files['template_S'].T + files['template_s']
    from scipy.ndimage import gaussian_filter
    truth_t = gaussian_filter(np.vectorize(P.PET.get)(P.tissues(pts)).reshape(e.shape).astype(float),
                              3.0 / abs(e.affine[0, 0]))          # the phantom PET is smoothed by 3 mm
    x = e.get_fdata()
    r0 = np.corrcoef(x.ravel(), truth_t.ravel())[0, 1]
    assert r0 > 0.98
    # aligned: a shift of one voxel in any direction decreases the agreement
    assert all(np.corrcoef(np.roll(x, sh, ax).ravel(), truth_t.ravel())[0, 1] < r0 for ax in range(3) for sh in (-1, 1))
    assert os.path.isfile(opj(paths['templates'], 'EDNiX', 'labels', 'Test_seg-TEST_dseg.nii.gz'))
    assert os.path.isfile(opj(paths['qc'], 'sub-01_ses-1_trc-FDG_space-acpc_pet.png'))


# ------------------------------------------------------------------------------------------------ launchers
CT_HU = {0: -1000, 1: 40, 2: 1200, 3: 10, 4: 38, 5: 30, 6: 38}


def _study(tmp_path, monkeypatch, subjects=('01',), anat=True, species='Testmonkey'):
    P.patch_run_cmd(monkeypatch)
    P.no_containers(monkeypatch)
    main, tdir = P.atlas_library(str(tmp_path), species)
    bids = str(tmp_path / 'bids')
    out = {}
    for k, sub in enumerate(subjects):
        if anat:
            files, package, dp = P.anat_derivative(bids, ID=sub, tdir=tdir, species=species)
        else:
            dp = opj(bids, 'sub-' + sub, 'ses-1')
            files = None
            P.anat_derivative(str(tmp_path / 'unused'), ID=sub, tdir=tdir, species=species)
        R, t = P.pose((8 - 4 * k, -5, 12), (6, -9 + 3 * k, 4))
        raw, A, info = P.dynamic_pet(opj(dp, 'pet', f'sub-{sub}_ses-1_trc-FDG_pet.nii.gz'), R, t)
        out[sub] = dict(files=files, dp=dp, R=R, t=t, A=A, raw=raw)
    with open(opj(bids, 'participants.tsv'), 'w') as f:
        f.write('participant_id\tbody_weight\n' + ''.join(f'sub-{s}\t10\n' for s in subjects))
    return main, bids, tdir, out


def _chain_error(chain, truth):
    """Max error (mm) of a linear 'to' chain (points of the fixed space through the matrices in list order, inverted
    when flagged) against truth(q), the true point of the moving space (RAS) of the fixed point q."""
    trs = []
    for f, inv in chain:
        tr = ants.read_transform(f)
        trs.append(tr.invert() if inv else tr)
    err = []
    for q in ([0, 0, 0], [15, 10, -5], [-12, -20, 10], [5, -25, -12]):
        q = np.asarray(q, float)
        p = tuple(P.F @ q)
        for tr in trs:
            p = tr.apply_to_point(p)
        err.append(np.linalg.norm(truth(q) - P.F @ np.asarray(p)))
    return max(err)


def _chain_pose_error(chain, R, t):
    """_chain_error of a PET -> AC-PC chain against the true pose x_PET = R^T (x_acpc - t)."""
    return _chain_error(chain, lambda q: R.T @ (q - t))


def _template_truth(R, t):
    """True PET point of a template point (template x -> native S x + s of _template_image, then the PET pose)."""
    S, sh = np.eye(3) / 1.08, np.array([2.0, -3.0, 1.5])
    return lambda q: R.T @ (S @ q + sh - t)


def _pose_error(mat_chain, R, t):
    import ants
    tr = ants.read_transform(mat_chain[0][0])
    err = []
    for q in ([0, 0, 0], [15, 10, -5], [-12, -20, 10], [5, -25, -12]):
        q = np.asarray(q, float)
        p_est = P.F @ np.asarray(tr.apply_to_point(tuple(P.F @ q)))
        err.append(np.linalg.norm(R.T @ (q - t) - p_est))
    return max(err)


@pytest.mark.slow
def test_preprocess_pet_launcher_scenario_1v2a(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     [''], [''], '', 'T1w', 'T1w', False, False, nb_void=2, SUV_window=[10, 40], srtm=True,
                     reference_signal=dict(type='atlas', atlas='TEST', level=0, labels='Cerebell', name='Cer'),
                     DV_time=10, selected_atlases=[['TEST', 'all']], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    from Tools import getpath
    paths = getpath.pet(s['01']['dp'], 'EDNiX')
    rec = transfo.read_record(coreg2anat.record_path(paths, 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '1v2a' and sorted(rec['targets']) == ['EDNiX', 'acpc']
    assert _pose_error(rec['targets']['acpc']['to'], s['01']['R'], s['01']['t']) < 1.0
    st = json.load(open(opj(paths['path'], 'Prepro_Settings.json')))
    assert st['Isotope'] == 'F-18' and st['Body weight (kg)'] == 10.0
    tab = pd.read_csv(opj(paths['orig_roi'], 'sub-01_ses-1_trc-FDG_seg-TEST-lvl1_desc-bratio-mean_ROIs.tsv'), sep='\t')
    assert set(tab['region']) == {name for _, name in P.LABELS}
    assert os.path.isfile(pu.out_name(opj(paths['templates'], 'EDNiX'), 'sub-01_ses-1_trc-FDG', desc='SUV-mean',
                                      space='EDNiX'))
    assert os.path.isfile(opj(paths['templates'], 'EDNiX', 'labels', 'Testmonkey_seg-TEST_dseg.nii.gz'))
    assert any('session' in f and f.endswith('.txt') for f in os.listdir(paths['path']))


@pytest.mark.slow
def test_ct_pipeline_then_pet_through_the_ct_1v2b(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    from modalities.ct import _0_Pipeline_launcher as LC
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    R, t = s['01']['R'], s['01']['t']
    Act = P.grid(80, 1.0)
    q = P.world_points((80,) * 3, Act) @ R.T + t                 # PET/CT scanner: CT in the space of the PET
    ct = np.vectorize(CT_HU.get)(P.tissues(q)).reshape((80,) * 3).astype(np.float32)
    P.save(ct, Act, opj(s['01']['dp'], 'ct', 'sub-01_ses-1_ct.nii.gz'))
    LC.preprocess_ct([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_ct.nii.gz', [''], [''], '',
                     'T1w', 'T1w', False, False, atlas_followers=[['TEST'], ['txt'], [2], [1]])
    from Tools import getpath
    ctp = getpath.ct(s['01']['dp'], 'EDNiX')
    crec = transfo.read_record(opj(ctp['orig_matrices'], 'sub-01_ses-1_desc-registration_ct.json'))
    assert crec['scenario'] == '1v2b' and _pose_error(crec['targets']['acpc']['to'], R, t) < 1.0
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     [''], [''], '', 'T1w', 'T1w', False, False, IhaveaCT=True, ct_pet_same_space=True,
                     selected_atlases=[['TEST', 'all']], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    rec = transfo.read_record(coreg2anat.record_path(getpath.pet(s['01']['dp'], 'EDNiX'), 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '1v2b'
    to = rec['targets']['acpc']['to']
    # PET follows the CT (chain of the CT), from the scanner position of its first registered frame (MoCo matrix)
    assert to[:len(crec['targets']['acpc']['to'])] == crec['targets']['acpc']['to']
    assert all('forMC' in f for f, _ in to[len(crec['targets']['acpc']['to']):])
    assert _chain_pose_error(to, R, t) < 1.0


@pytest.mark.slow
def test_pet_study_template_without_anatomy_3v2(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    main, bids, tdir, s = _study(tmp_path, monkeypatch, subjects=('01', '02'), anat=False)
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1'), ('02', '1')]), '*_trc-FDG_*pet.nii.gz',
                     'FDG', [''], [''], '', 'T1w', 'T1w', False, False, IhaveanANAT=False,
                     creat_study_template_pet=True, studyT_transform='TRSAA', crop_for_MoCo='auto',
                     selected_atlases=[['TEST', 'all']], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    st = studyT.template_name(opj(bids, 'sty_template', 'pet'), 'FDG')
    assert os.path.isfile(st) and os.path.isfile(studyT.crop_name(st))
    from Tools import getpath
    for sub in ('01', '02'):
        paths = getpath.pet(s[sub]['dp'], 'EDNiX')
        rec = transfo.read_record(coreg2anat.record_path(paths, f'sub-{sub}_ses-1_trc-FDG'))
        assert rec['scenario'] == '3v2' and list(rec['targets']) == ['studyTemplate'] and not rec['pet_space']['atlases']
        assert os.path.isfile(pu.out_name(paths['orig_post'], f'sub-{sub}_ses-1_trc-FDG', desc='SUV-mean'))
        assert os.path.isfile(pu.out_name(opj(paths['templates'], 'studyTemplate'), f'sub-{sub}_ses-1_trc-FDG',
                                          desc='SUV-mean', space='studyTemplate'))


@pytest.mark.slow
def test_hybrid_pet_mr_1v1_follows_the_reorientation_of_the_anatomy(tmp_path, monkeypatch):
    """Sphinx animal (AHF declared HFS) in a hybrid PET-MR; T1w and PET share the scanner space, not their grids, and
    their fields of view lie past the coordinate origin, where 3drefit -orient re-signs the origins (the emulator
    follows 3drefit.c). The PET receives the world transform of the reorientation of the T1w (an independent
    relabelling would shift it by 7 mm), and the motion-corrected session is brought back to the head position of
    the first 10 min (the head moved later in the run)."""
    from modalities.pet import _0_Pipeline_launcher as L
    from Tools import get_orientation, getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    dp, A0, files = s['01']['dp'], s['01']['A'], s['01']['files']
    A0 = nib.load(files['reorient']).affine                         # reoriented T1w = AC-PC grid (identity AC-PC)
    A_anat_raw = P.raw_header(A0, 'HFS', 'AHF')                     # raw T1w in the scanner space
    t1 = nib.load(files['reorient']).get_fdata().astype(np.float32)
    P.save(t1, A_anat_raw, files['raw'])
    W = A0 @ np.linalg.inv(A_anat_raw)                              # scanner -> AC-PC world
    Ap = P.grid(40, 2.0)                                            # PET grid once in the AC-PC world
    A_pet_raw = np.linalg.inv(W) @ Ap
    indep = P.refit_orient(A_pet_raw, get_orientation.getreal('HFS', 'AHF', P._afni_code(A_pet_raw)))
    assert np.abs(indep - Ap)[:3, 3].max() > 5.0                    # the hazard of an independent relabelling
    still = (np.eye(3), np.zeros(3))
    moved = (P.pose((4, -3, 5), (4, -3, 3))[0], np.array([4.0, -3.0, 3.0]))
    raw = opj(dp, 'pet', 'sub-01_ses-1_trc-FDG_pet.nii.gz')
    P.dynamic_pet(raw, np.eye(3), np.zeros(3), motions=[still] * 5 + [moved] * 3)
    img = nib.load(raw)
    nib.save(nib.Nifti1Image(np.asarray(img.dataobj), A_pet_raw), raw)
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     ['AHF'], ['HFS'], '', 'T1w', 'T1w', False, False, anat_pet_same_space=True, SUV_window=[0, 10],
                     crop_for_MoCo='none', selected_atlases=[['TEST', 0]], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    paths = getpath.pet(dp, 'EDNiX')
    reo = pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG', desc='reorient')
    assert np.allclose(nib.load(reo).affine, Ap, atol=1e-4)        # PET in the world of the reoriented T1w
    rec = transfo.read_record(coreg2anat.record_path(paths, 'sub-01_ses-1_trc-FDG'))
    to = rec['targets']['acpc']['to']
    assert rec['scenario'] == '1v1' and to[0][0].endswith('acpc_0GenericAffine.mat') and 'scanner-position' in to[1][0]
    # frames of 0-10 min (still) exported to AC-PC: on the truth with the anchor, shifted without it
    exp = pu.out_name(paths['acpc'], 'sub-01_ses-1_trc-FDG', desc='activity-mean', space='acpc')
    grid = nib.load(exp)
    truth = P.pet_frame(grid.shape[:3], grid.affine, np.eye(3), np.zeros(3))
    act = pu.out_name(paths['orig_post'], 'sub-01_ses-1_trc-FDG', desc='activity-mean')
    no_anchor = transfo.apply(exp, act, dict(to=to[:1]), 'to').numpy()
    r = lambda a: np.corrcoef(a.ravel(), truth.ravel())[0, 1]
    assert r(grid.get_fdata()) > 0.99 and r(grid.get_fdata()) > r(no_anchor) + 0.01, (r(grid.get_fdata()), r(no_anchor))


@pytest.mark.slow
def test_pet_ct_same_scanner_space_with_off_centre_fields_of_view(tmp_path, monkeypatch):
    """PET/CT of a sphinx animal: CT (1 mm, 80 mm) and PET (2 mm, 80 mm) grids past the coordinate origin. The CT is
    reoriented from the positions, the PET receives the CT's world transform, and the PET -> AC-PC chain (CT
    registration + head position of the first 10 min) recovers the pose."""
    from modalities.pet import _0_Pipeline_launcher as L
    from modalities.ct import _0_Pipeline_launcher as LC
    from Tools import getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    dp, R, t = s['01']['dp'], s['01']['R'], s['01']['t']
    Act, Ap = P.grid(80, 1.0), P.grid(40, 2.0)                      # grids once reoriented (x_acpc = R x + t)
    A_ct_raw = P.raw_header(Act, 'HFS', 'AHF')
    W = Act @ np.linalg.inv(A_ct_raw)
    q = P.world_points((80,) * 3, Act) @ R.T + t
    P.save(np.vectorize(CT_HU.get)(P.tissues(q)).reshape((80,) * 3).astype(np.float32), A_ct_raw,
           opj(dp, 'ct', 'sub-01_ses-1_ct.nii.gz'))
    raw = opj(dp, 'pet', 'sub-01_ses-1_trc-FDG_pet.nii.gz')
    P.dynamic_pet(raw, R, t)
    img = nib.load(raw)
    nib.save(nib.Nifti1Image(np.asarray(img.dataobj), np.linalg.inv(W) @ Ap), raw)
    LC.preprocess_ct([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_ct.nii.gz', ['AHF'], ['HFS'], '',
                     'T1w', 'T1w', False, False, atlas_followers=[['TEST'], ['txt'], [2], [1]])
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     ['AHF'], ['HFS'], '', 'T1w', 'T1w', False, False, IhaveaCT=True, ct_pet_same_space=True,
                     crop_for_MoCo='none', selected_atlases=[['TEST', 0]], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    paths = getpath.pet(dp, 'EDNiX')
    reo = pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG', desc='reorient')
    assert np.allclose(nib.load(reo).affine, Ap, atol=1e-4)
    rec = transfo.read_record(coreg2anat.record_path(paths, 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '1v2b' and _chain_pose_error(rec['targets']['acpc']['to'], R, t) < 1.0


@pytest.mark.slow
def test_ct_study_template_without_anatomy_2v2(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    from modalities.ct import _0_Pipeline_launcher as LC
    from Tools import getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch, subjects=('01', '02'), anat=False)
    for sub in ('01', '02'):
        R, t = s[sub]['R'], s[sub]['t']
        Act = P.grid(80, 1.0)
        q = P.world_points((80,) * 3, Act) @ R.T + t
        P.save(np.vectorize(CT_HU.get)(P.tissues(q)).reshape((80,) * 3).astype(np.float32), Act,
               opj(s[sub]['dp'], 'ct', f'sub-{sub}_ses-1_ct.nii.gz'))
    table = P.bids_table([('01', '1'), ('02', '1')])
    LC.preprocess_ct([], main, bids, 'Testmonkey', table, '*_ct.nii.gz', [''], [''], '', 'T1w', 'T1w', False, False,
                     IhaveanANAT=False, creat_study_template_ct=True, studyT_transform='TRSAA',
                     atlas_followers=[['TEST'], ['txt'], [2], [1]])
    L.preprocess_pet([], main, bids, 'Testmonkey', table, '*_trc-FDG_*pet.nii.gz', 'FDG', [''], [''], '', 'T1w', 'T1w',
                     False, False, IhaveanANAT=False, IhaveaCT=True, crop_for_MoCo='none',
                     selected_atlases=[['TEST', 'all']], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    from modalities.pet import bridge
    assert os.path.isfile(bridge.record_path(opj(bids, 'sty_template', 'pet'), 'EDNiX'))
    for sub in ('01', '02'):
        paths = getpath.pet(s[sub]['dp'], 'EDNiX')
        rec = transfo.read_record(coreg2anat.record_path(paths, f'sub-{sub}_ses-1_trc-FDG'))
        # step 'bridge': the CT study template is linked to the reference template through the PET study template
        assert rec['scenario'] == '2v2' and sorted(rec['targets']) == ['EDNiX', 'studyTemplate']
        err = _chain_error(rec['targets']['EDNiX']['to'], _template_truth(s[sub]['R'], s[sub]['t']))
        assert err < 2.0, err
        atlas = nib.load(rec['pet_space']['atlases']['TEST']).get_fdata()[..., 0]
        A = nib.load(rec['reference']).affine
        truth = P.tissues(P.world_points(atlas.shape, A) @ s[sub]['R'].T + s[sub]['t']).reshape(atlas.shape)
        dice = 2 * ((atlas == 7) & (truth == 6)).sum() / ((atlas == 7).sum() + (truth == 6).sum())
        assert dice > 0.7, dice
        assert os.path.isfile(pu.out_name(opj(paths['templates'], 'studyTemplate'), f'sub-{sub}_ses-1_trc-FDG',
                                          desc='activity-mean', space='studyTemplate'))
    # both subjects in the study template: their exported means agree
    a, b = [nib.load(pu.out_name(opj(getpath.pet(s[x]['dp'], 'EDNiX')['templates'], 'studyTemplate'),
                                 f'sub-{x}_ses-1_trc-FDG', desc='activity-mean', space='studyTemplate')).get_fdata()
            for x in ('01', '02')]
    assert np.corrcoef(a.ravel(), b.ravel())[0, 1] > 0.97


def _template_image(tdir, values, species, name):
    """Image of the tissues of the template space (template x -> native S x + s) with the given values."""
    from scipy.ndimage import gaussian_filter
    At = P.grid(80, 1.0)
    S = np.eye(3) / 1.08
    q = P.world_points((80,) * 3, At) @ S.T + np.array([2.0, -3.0, 1.5])
    v = np.vectorize(values.get)(P.tissues(q)).reshape((80,) * 3).astype(np.float32)
    return P.save(gaussian_filter(v, 1.5), At, opj(tdir, species + name))


@pytest.mark.slow
def test_pet_template_without_anatomy_3v1(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    from Tools import getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch, anat=False)
    tpl = _template_image(tdir, P.PET, 'Testmonkey', '_desc-FDG_pet.nii.gz')
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     [''], [''], '', 'T1w', 'T1w', False, False, IhaveanANAT=False, template_PET=tpl,
                     type_of_transform_template='TRSAA', selected_atlases=[['TEST', 'all']],
                     atlas_followers=[['TEST'], ['txt'], [2], [1]])
    paths = getpath.pet(s['01']['dp'], 'EDNiX')
    rec = transfo.read_record(coreg2anat.record_path(paths, 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '3v1' and list(rec['targets']) == ['EDNiX']
    atlas_pet = nib.load(rec['pet_space']['atlases']['TEST']).get_fdata()[..., 0]
    q = P.world_points((40,) * 3, s['01']['A']) @ s['01']['R'].T + s['01']['t']
    truth = P.tissues(q).reshape((40,) * 3)
    brain = 2 * ((atlas_pet > 0) & (truth >= 3)).sum() / ((atlas_pet > 0).sum() + (truth >= 3).sum())
    assert brain > 0.9, brain
    assert _chain_error(rec['targets']['EDNiX']['to'], _template_truth(s['01']['R'], s['01']['t'])) < 1.5
    assert set(rec['pet_space']['masks']) == {'brain', 'gray', 'white', 'csf'}


@pytest.mark.slow
def test_ct_template_without_anatomy_2v1(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    from modalities.ct import _0_Pipeline_launcher as LC
    from Tools import getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch, anat=False)
    R, t = s['01']['R'], s['01']['t']
    Act = P.grid(80, 1.0)
    q = P.world_points((80,) * 3, Act) @ R.T + t
    P.save(np.vectorize(CT_HU.get)(P.tissues(q)).reshape((80,) * 3).astype(np.float32), Act,
           opj(s['01']['dp'], 'ct', 'sub-01_ses-1_ct.nii.gz'))
    tct = _template_image(tdir, CT_HU, 'Testmonkey', '_desc-template_ct.nii.gz')
    LC.preprocess_ct([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_ct.nii.gz', [''], [''], '',
                     'T1w', 'T1w', False, False, IhaveanANAT=False, template_CT=tct,
                     type_of_transform_template='TRSAA', atlas_followers=[['TEST'], ['txt'], [2], [1]])
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     [''], [''], '', 'T1w', 'T1w', False, False, IhaveanANAT=False, IhaveaCT=True,
                     selected_atlases=[['TEST', 'all']], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    rec = transfo.read_record(coreg2anat.record_path(getpath.pet(s['01']['dp'], 'EDNiX'), 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '2v1' and list(rec['targets']) == ['EDNiX'] and rec['pet_space']['atlases']
    assert _chain_error(rec['targets']['EDNiX']['to'], _template_truth(R, t)) < 1.5
    assert os.path.isfile(opj(getpath.ct(s['01']['dp'], 'EDNiX')['qc'], 'sub-01_ses-1_space-EDNiX_ct.png'))


@pytest.mark.slow
def test_between_runs_correction_aligns_the_runs_with_one_resampling(tmp_path, monkeypatch):
    P.patch_run_cmd(monkeypatch)
    from Tools import getpath
    dp = str(tmp_path / 'sub-01' / 'ses-1')
    R0, t0 = np.eye(3), np.zeros(3)
    raw1, A, _ = P.dynamic_pet(opj(dp, 'pet', 'sub-01_ses-1_trc-FDG_run-1_pet.nii.gz'), R0, t0)
    m = (P.pose((5, -4, 6), (4, -3, 2))[0], np.array([4.0, -3.0, 2.0]))
    raw2, _, _ = P.dynamic_pet(opj(dp, 'pet', 'sub-01_ses-1_trc-FDG_run-2_pet.nii.gz'), R0, t0, motions=[m] * 8)
    paths = getpath.pet(dp, 'EDNiX')
    diary = str(tmp_path / 'd.txt')
    reor = loop1.run([raw1, raw2], paths, [''], [''], '', '', '', False, False, 'WARP', 'F-18', False, ('', ''), [],
                     diary, '', '', '')
    opts = dict(Skip_step=[6, 7, 8, 9, 10], nb_void=2, crop_for_MoCo='none', redo_crop=False, MoCo_transform='Rigid',
                coreg_metric='mattes', Need_check_MC=False, srtm=False, pvc_fwhm=None, do_pet_to_func=False)
    loop2.run('01', '1', dp, str(tmp_path), paths, 'sub-01_ses-1_trc-FDG', reor, [raw1, raw2], opts, dict(), diary,
              '', '', '', '')
    b1 = nib.load(pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG_run-1', desc='MoCo-between-mean')).get_fdata()
    b2 = nib.load(pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG_run-2', desc='MoCo-between-mean')).get_fdata()
    r1 = nib.load(pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG_run-1', desc='MoCo-within-mean')).get_fdata()
    r2 = nib.load(pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG_run-2', desc='MoCo-within-mean')).get_fdata()
    assert np.corrcoef(b1.ravel(), b2.ravel())[0, 1] > 0.995 > np.corrcoef(r1.ravel(), r2.ravel())[0, 1]
    side = json.load(open(pu.sidecar(pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG_run-2', desc='MoCo-between'))))
    assert 'composed' in side['Description']
    assert os.path.isfile(pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG', desc='MoCo-between-mean'))


# ------------------------------------------------------------------------------------------------ deep-review additions
def test_frame_times_are_post_injection_and_fmri_matching_stays_in_bids_time():
    # TimeZero = ScanStart, injection 30 min before the scan
    info = dict(FrameTimesStart=[0, 600, 1200], FrameDuration=[600, 600, 600], InjectionStart=-1800)
    tm = pu.frame_timing(info, 3)
    assert np.allclose(tm['start'], [1800, 2400, 3000]) and tm['offset'] == -1800
    assert list(pu.window_frames(tm, [40, 60])) == [1, 2]          # 40-60 min after the injection
    func = dict(AcquisitionTime='10:20:00', RepetitionTime=2.0)       # 1200-1500 s after TimeZero: third PET frame
    ov, _ = run_cut.getframes(dict(TimeZero='10:00:00'), tm, func, 150, 0)
    assert np.allclose(ov, [0, 0, 300])
    # midnight: TimeZero 23:50, fMRI at 00:05 the next day
    ov, (a, b) = run_cut.getframes(dict(TimeZero='23:50:00'), pu.frame_timing(
        dict(FrameTimesStart=[0, 600, 1200], FrameDuration=[600, 600, 600]), 3), dict(AcquisitionTime='00:05:00',
                                                                                     RepetitionTime=2.0), 150, 0)
    assert a == 900 and ov[1] == 300


def test_units_isotopes_and_timing_are_parsed_robustly():
    assert pu.normalise_isotope('[18F]') == 'F-18' and pu.normalise_isotope('C11') == 'C-11'
    assert pu.activity_to_mbq(10, 'mci') == pytest.approx(370.0) and pu.activity_to_mbq(5, 'mbq') == 5.0
    assert pu.concentration_to_bq_per_ml('kBq/cm3') == 1e3 and pu.concentration_to_bq_per_ml('BQML') == 1.0
    assert pu.concentration_to_bq_per_ml('counts') is None and pu.as_bool('false') is False
    tm = pu.frame_timing(dict(FrameTimesStart=[0], FrameDuration=[600], FrameReferenceTime=300), 1)
    assert tm['mid'][0] == 300
    tm = pu.frame_timing(dict(FrameTimesStart=[0, 60], FrameDuration=[60, 60], FrameReferenceTime=[30000, 90000]), 2)
    assert np.allclose(tm['mid'], [30, 90])                           # ms reference times ignored
    with pytest.raises(ValueError):
        pu.frame_timing(dict(FrameTimesStart=[0, 30], FrameDuration=[60, 60]), 2)   # overlapping frames
    # unknown isotope and image decay-corrected after the injection: the dose cannot be decayed
    f, why = pu.suv_scale(dict(Units='Bq/mL', ImageDecayCorrectionTime=3600, InjectionStart=0),
                          dict(weight_kg=70, dose_MBq=200), None)
    assert f is None and 'isotope' in why
    f, rep = pu.suv_scale(dict(Units='Bq/mL', ImageDecayCorrectionTime=3600, InjectionStart=0),
                          dict(weight_kg=70, dose_MBq=200), 'F-18')
    assert f == pytest.approx(70e3 / (200e6 * math.exp(-math.log(2) / 6586.2 * 3600)))


def test_mean_rigid_and_scanner_anchor(tmp_path):
    from scipy.spatial.transform import Rotation
    d = str(tmp_path)
    R1 = Rotation.from_euler('xyz', (3, -2, 5), degrees=True).as_matrix()
    m1 = P.write_affine_mat(opj(d, 'm1.mat'), R1, np.array([1.0, 2, 3]))
    same = MoCo.mean_rigid([m1, m1], opj(d, 'mean.mat'))
    a, b = ants.read_transform(m1), ants.read_transform(same)
    assert max(np.linalg.norm(np.subtract(a.apply_to_point(p), b.apply_to_point(p)))
               for p in [(0, 0, 0), (20, -10, 5)]) < 1e-4
    # two runs: session -> run (between) -> frames of the window (within), and its exact inverse
    runs = [opj(d, f'sub-01_trc-FDG_run-{k}_desc-reorient_pet.nii.gz') for k in (1, 2)]
    corrected = [opj(d, f'sub-01_trc-FDG_run-{k}_desc-MoCo-within_pet.nii.gz') for k in (1, 2)]
    for k in (1, 2):
        P.write_affine_mat(opj(d, f'sub-01_trc-FDG_run-{k}_between_0GenericAffine.mat'), np.eye(3), np.array([k, 0, 0.]))
    mats = []
    for i in range(3):
        mats.append(f'f{i}.mat')
        P.write_affine_mat(opj(d, mats[-1]), np.eye(3), np.array([0, i, 0.]))
    for c in corrected:
        pu.write_json(c, '', '', Matrices=mats)
    tms = [pu.frame_timing(dict(FrameTimesStart=[0, 300, 900], FrameDuration=[300, 600, 600]), 3)] * 2
    ch, how = MoCo.scanner_anchor(runs, corrected, tms, d, nb_void=0, window_min=[0, 6])
    assert [opb(f) for f, _ in ch['from']] == ['sub-01_trc-FDG_run-1_between_0GenericAffine.mat', 'f0.mat']
    assert [inv for _, inv in ch['to']] == [True, True] and ch['to'][0][0] == ch['from'][1][0]
    ch, _ = MoCo.scanner_anchor(runs, corrected, tms, d, nb_void=0)          # first 10 min: frames 0 and 1
    assert 'scanner-position' in ch['from'][1][0]


def test_world_transform_of_the_reference_is_applied_to_the_header(tmp_path, monkeypatch):
    from modalities.pet import reorient
    P.patch_run_cmd(monkeypatch)
    A0 = P.grid(30, 1.0)
    raw = P.save(np.random.rand(30, 30, 30), P.raw_header(A0, 'HFS', 'AHF'), str(tmp_path / 'ref_raw.nii.gz'))
    reo = str(tmp_path / 'ref_reo.nii.gz')
    P.save(nib.load(raw).get_fdata(), A0, reo)
    W = reorient.reference_world_map(raw, reo, str(tmp_path), str(tmp_path / 'd.txt'), '')
    other = P.save(np.random.rand(10, 12, 14), np.linalg.inv(W) @ P.grid(10, 3.0), str(tmp_path / 'o.nii.gz'))
    out = reorient.apply_world_map(other, W, str(tmp_path / 'o2.nii.gz'))
    assert np.allclose(nib.load(out).affine, P.grid(10, 3.0)) and np.allclose(nib.load(out).get_fdata(),
                                                                              nib.load(other).get_fdata(), atol=1e-6)
    pu.write_json(reo, raw, 'x', WorldTransform=np.eye(4).tolist())
    assert np.allclose(reorient.reference_world_map(raw, reo, str(tmp_path), str(tmp_path / 'd.txt'), ''), np.eye(4))


def test_tissue_fractions_and_pvc_classes(tmp_path):
    from modalities.pet import pvc
    A1, A2 = P.grid(40, 1.0), P.grid(10, 4.0)
    x = P.world_points((40,) * 3, A1)[:, 0].reshape((40,) * 3)
    half = P.save((x < 0).astype(np.float32), A1, str(tmp_path / 'half.nii.gz'))
    ref = P.save(np.ones((10,) * 3), A2, str(tmp_path / 'ref.nii.gz'))
    ident = transfo.affine(P.write_affine_mat(str(tmp_path / 'id.mat'), np.eye(3), np.zeros(3)))
    fr = nib.load(coreg2anat.tissue_fraction(half, ref, ident, (4.0, 4.0, 4.0), str(tmp_path / 'f.nii.gz'))).get_fdata()
    xp = P.world_points((10,) * 3, A2)[:, 0].reshape((10,) * 3)
    assert np.allclose(fr[xp < -6][..., ], 1.0, atol=0.02) and np.allclose(fr[xp > 6], 0.0, atol=0.02)
    assert 0.25 < fr[np.abs(xp) == 2].mean() < 0.75                  # voxels on the edge: partial fraction
    gm = P.save(np.full((10,) * 3, 0.5), A2, str(tmp_path / 'gm.nii.gz'))
    wm = P.save(np.full((10,) * 3, 0.7), A2, str(tmp_path / 'wm.nii.gz'))
    br = P.save(np.full((10,) * 3, 0.9), A2, str(tmp_path / 'br.nii.gz'))
    cls = pvc.tissue_classes(dict(GM=gm, WM=wm, brain=br), ref)
    assert np.allclose(sum(cls), 1.0) and np.allclose(cls[1], 0.5) and np.allclose(cls[2], 0.0)


@pytest.mark.slow
def test_static_single_frame_not_decay_corrected_gives_the_analytic_suv(tmp_path, monkeypatch):
    """One 10-min frame 50 min after the injection (TimeZero = scan start), not decay-corrected, stored as 4D with one
    frame: the SUV of the white matter must be C x W / D with C decay-corrected to the injection."""
    from modalities.pet import _0_Pipeline_launcher as L
    from Tools import getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    dp, R, t = s['01']['dp'], s['01']['R'], s['01']['t']
    raw = opj(dp, 'pet', 'sub-01_ses-1_trc-FDG_pet.nii.gz')
    A = P.grid(40, 2.0)
    lam = math.log(2) / 6586.2
    t0, dt = 3000.0, 600.0
    meas = math.exp(-lam * t0) * (1 - math.exp(-lam * dt)) / (lam * dt)    # mean decay over the frame
    P.save(P.pet_frame((40,) * 3, A, R, t, smooth_mm=0.1)[..., None] * meas, A, raw)
    with open(raw.replace('.nii.gz', '.json'), 'w') as f:
        json.dump(dict(FrameTimesStart=[0], FrameDuration=[dt], Units='Bq/mL', InjectedRadioactivity=100.0,
                       InjectedRadioactivityUnits='MBq', ImageDecayCorrected='false', InjectionStart=-t0,
                       TimeZero='10:00:00', TracerRadionuclide='[18F]'), f)
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     [''], [''], '', 'T1w', 'T1w', False, False, selected_atlases=[['TEST', 1]],
                     atlas_followers=[['TEST'], ['txt'], [2], [1]])
    paths = getpath.pet(dp, 'EDNiX')
    assert nib.load(pu.out_name(paths['orig_process'], 'sub-01_ses-1_trc-FDG', desc='reorient')).ndim == 3
    tab = pd.read_csv(opj(paths['orig_roi'], 'sub-01_ses-1_trc-FDG_seg-TEST-lvl1_desc-SUV_ROIs.tsv'), sep='\t')
    wm = float(tab.loc[tab['region'] == 'WM', 'median'].iloc[0])
    assert wm == pytest.approx(P.PET[5] * 10e3 / 100e6, rel=0.02)          # 1500 Bq/mL x 10 kg / 100 MBq
    assert os.path.isfile(pu.out_name(paths['acpc'], 'sub-01_ses-1_trc-FDG', desc='activity', space='acpc'))


def test_anat_of_another_session_and_reference_session(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    bids = str(tmp_path / 'bids')
    files, package, dp1 = P.anat_derivative(bids, ses='1')
    dp2 = opj(bids, 'sub-01', 'ses-2')
    os.makedirs(opj(dp2, 'pet'))
    anat, apath, ases = L._anat_files('01', dp2, '2', ['2'], bids, False, True, True, 'T1w', str(tmp_path / 'd.txt'))
    assert ases == '1' and os.path.normpath(apath) == os.path.normpath(dp1) and anat['acpc'] == files['acpc']
    assert L._anat_max_session(bids, '01', 'EDNiX', '2') == '1'


@pytest.mark.slow
def test_rerun_of_the_late_steps_and_skipped_registration(tmp_path, monkeypatch):
    """Steps 7-10 run again alone; the scenario of the existing record is kept when step 6 is skipped (the CT record,
    needed only by step 6, is then not required)."""
    from modalities.pet import _0_Pipeline_launcher as L
    from Tools import getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    args = ([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG', [''], [''], '',
            'T1w', 'T1w', False, False)
    kw = dict(crop_for_MoCo='none', selected_atlases=[['TEST', 0]], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    L.preprocess_pet(*args, **kw)
    paths = getpath.pet(s['01']['dp'], 'EDNiX')
    exp = pu.out_name(paths['acpc'], 'sub-01_ses-1_trc-FDG', desc='SUV-mean', space='acpc')
    first = os.path.getmtime(exp)
    L.preprocess_pet([1, 2, 3, 4, 5, 6], *args[1:], IhaveaCT=True, **kw)    # no CT record: fine, step 6 skipped
    assert os.path.getmtime(exp) > first
    with pytest.raises(FileNotFoundError):
        L.preprocess_pet([1, 2, 3, 4, 5], *args[1:], IhaveaCT=True, **kw)


def test_launcher_rejects_old_entry_point_and_bad_options(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    with pytest.raises(TypeError):
        L.preprocess_data([], 'x')
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    args = ([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG', [''], [''], '',
            'T1w', 'T1w', False, False)
    with pytest.raises(ValueError):
        L.preprocess_pet(*args, isotope='F-19', atlas_followers=[['TEST'], ['txt'], [2], [1]])
    with pytest.raises(ValueError):
        L.preprocess_pet(*args, pvc_fwhm=5, pvc_method='GTM', atlas_followers=[['TEST'], ['txt'], [2], [1]])
    with pytest.raises(ValueError):
        L.preprocess_pet(*args, MoCo_transform='SyN', atlas_followers=[['TEST'], ['txt'], [2], [1]])


@pytest.mark.slow
def test_within_run_motion_estimates_are_accurate(tmp_path):
    """Relative motion between frames recovered within 0.6 mm (mean) on a noisy phantom (3-6 mm, 4-6 degrees)."""
    R, t = np.eye(3), np.zeros(3)
    moves = [((4, 0, 6), (3, -2, 1))] * 3 + [((-5, 3, 0), (-4, 2, 2)), ((0, -6, 4), (2, 4, -3)),
                                            ((3, 3, -5), (-3, -3, 2)), ((6, 0, 0), (0, 3, 3)), ((-4, -4, 4), (4, 0, -2))]
    motions = [(P.pose(r, m)[0], np.asarray(m, float)) for r, m in moves]
    img, A, info = P.dynamic_pet(str(tmp_path / 'sub-01_ses-1_trc-FDG_pet.nii.gz'), R, t, motions=motions)
    im = nib.load(img)
    x = np.asarray(im.dataobj, np.float32)
    rng = np.random.default_rng(0)
    dur = np.asarray(info['FrameDuration'], float)
    x = x + rng.normal(size=x.shape).astype(np.float32) * np.sqrt(np.clip(x, 0, None) * 60 / dur)
    nib.save(nib.Nifti1Image(x, im.affine), img)
    MoCo.dynamic(img, str(tmp_path / 'pre'), str(tmp_path / 'formc'), pu.frame_timing(info, 8), 2, None, False,
                 str(tmp_path / 'diary.txt'), '')
    mats = json.load(open(str(tmp_path / 'pre' / 'sub-01_ses-1_trc-FDG_desc-MoCo-within_pet.json')))['Matrices']
    T = [ants.read_transform(str(tmp_path / 'formc' / m)) for m in mats]
    pts = [np.array(p, float) for p in ([0, 0, 0], [15, 10, -5], [-12, -20, 10], [5, -25, -12])]
    errs = []
    for i in range(2, 8):
        for j in range(i + 1, 8):
            Ti = T[i].invert()
            for p in pts:
                true = motions[j][0].T @ (motions[i][0] @ p + motions[i][1] - motions[j][1])
                est = P.F @ np.asarray(T[j].apply_to_point(Ti.apply_to_point(tuple(P.F @ p))))
                errs.append(np.linalg.norm(true - est))
    assert np.mean(errs) < 0.6 and np.max(errs) < 1.5, (np.mean(errs), np.max(errs))


def test_session_crop_redrawn_when_the_session_image_changed(tmp_path, monkeypatch):
    # two sessions of one subject with different head positions: each keeps its own crop; a crop drawn before step 1
    # was run again with another orientation (session 1) no longer matches its image and is drawn again
    from Tools import getpath
    drawn = []
    diary = str(tmp_path / 'diary.txt')

    def draw(image, crop, diary, sing_fs, redo=False):               # freeview replaced by a box around the head
        if os.path.isfile(crop) and not redo:
            return crop
        img = ants.image_read(image)
        head = MoCo.head_mask(img, 0.3, 0)
        ants.image_write(ants.crop_image(img, label_image=head, label=1), crop, ri=False)
        drawn.append(os.path.basename(crop))
        return crop
    monkeypatch.setattr(studyT, 'manual_crop_session', draw)
    lut = [P.PET[k] for k in range(7)]
    A = P.grid(48, 2.0)
    boxes = {}
    for ses, pose in (('1', ((5, 0, 0), (8, -6, 4))), ('2', ((-8, 4, 0), (-10, 5, -3)))):
        dp = str(tmp_path / 'bids' / 'sub-01' / f'ses-{ses}')
        paths = getpath.pet(dp, 'EDNiX')
        os.makedirs(paths['orig_formc'], exist_ok=True)
        mean = P.save(_render((48,) * 3, A, *P.pose(*pose), lut), A, os.path.join(paths['orig_formc'], 'mean.nii.gz'))
        base = f'sub-01_ses-{ses}_trc-FDG'
        boxes[ses] = (loop2.crop_for_session('auto', '', base, paths, mean, diary, '', False), mean, paths, base)
        assert os.path.isfile(os.path.join(paths['qc'], base + '_desc-cropbox_pet.png'))
    assert drawn == ['sub-01_ses-1_trc-FDG_desc-crop_pet.nii.gz', 'sub-01_ses-2_trc-FDG_desc-crop_pet.nii.gz']
    assert os.path.isfile(str(tmp_path / 'bids' / 'QC' / 'pet_crop' / 'sub-01_ses-1_trc-FDG_desc-cropbox_pet.png'))
    # same images: the crops are reused
    for ses in ('1', '2'):
        _, mean, paths, base = boxes[ses]
        loop2.crop_for_session('auto', '', base, paths, mean, diary, '', False)
    assert len(drawn) == 2
    # session 1 reoriented again (axes flipped: other world coordinates for the same voxels): drawn again
    _, mean, paths, base = boxes['1']
    img = nib.load(mean)
    nib.save(nib.Nifti1Image(np.asarray(img.dataobj)[::-1, :, ::-1], img.affine), mean)
    box = loop2.crop_for_session('auto', '', base, paths, mean, diary, '', False)
    assert drawn[-1] == 'sub-01_ses-1_trc-FDG_desc-crop_pet.nii.gz' and len(drawn) == 3
    head = MoCo.head_mask(ants.image_read(mean), 0.3, 0).numpy() > 0
    assert (ants.image_read(box).numpy() > 0)[head].mean() > 0.99


def _folded_brain(out, sp=0.5, fwhm=3.0, gap_mm=2.0, muscles='lateral'):
    """Macaque-like brain (ellipsoid 60 x 50 x 44 mm) with a 2 mm cortex folded by radial sulci, its T1, its masks
    (names of the anatomical pipeline) and an FDG PET (grey 1, white 0.35, CSF 0.05, temporal muscles 0.6 at
    gap_mm outside the brain) smoothed to fwhm and sampled at 1 mm."""
    from scipy import ndimage
    from modalities.pet import boundary
    n = (170, 150, 140)
    g = (np.indices(n).T - np.array(n) / 2).T * sp
    r = np.sqrt((g[0] / 30) ** 2 + (g[1] / 25) ** 2 + (g[2] / 22) ** 2)
    brain, head = r < 1, r < 1.35
    theta, phi = np.arctan2(g[1], g[0]), np.arctan2(g[2], np.hypot(g[0], g[1]))
    sulc = np.zeros(n, bool)
    for k in range(0, 10, 2):
        sulc |= (np.abs(np.sin(5 * theta + k)) < 0.04) & (r > 0.7)
    sulc = (sulc | ((np.abs(np.sin(4 * phi)) < 0.04) & (r > 0.7))) & brain
    gm = brain & ~sulc & (ndimage.distance_transform_edt(~(sulc | ~brain), sampling=sp) < 2.0)
    wm = brain & ~sulc & ~gm
    muscle = head & (ndimage.distance_transform_edt(~brain, sampling=sp) > gap_mm)
    muscle &= (np.abs(g[0]) > 22) & (np.abs(g[2]) < 12) if muscles == 'lateral' else g[2] < 15   # 'all': not dorsal
    def img(a):
        return ants.from_numpy(a.astype('float32'), spacing=(sp,) * 3)
    f = dict(mask=opj(out, 'M_space-acpc_mask.nii.gz'), gm=opj(out, 'M_desc-Whole_Brain_Gray_mask.nii.gz'),
             wm=opj(out, 'M_desc-White_mask.nii.gz'), t1=opj(out, 'M_space-acpc_T1w.nii.gz'), pet=opj(out, 'pet.nii.gz'))
    for k, a in (('mask', brain), ('gm', gm), ('wm', wm), ('t1', 0.4 * gm + 0.8 * wm + 0.6 * (head & ~brain))):
        ants.image_write(img(a), f[k])
    s = np.where(brain, 0.05, 0.0)
    s[wm], s[gm], s[muscle] = 0.35, 1.0, 0.6
    pet = ants.resample_image(boundary._gauss(img(s), fwhm), (1.0, 1.0, 1.0), use_voxels=False, interp_type=0)
    ants.image_write(pet, f['pet'])
    centre = list(np.array(n) / 2 * sp)
    return f, dict(brain=brain, gm=gm, sp=sp), centre


def test_boundary_check_separates_partial_volume_from_misregistration(tmp_path):
    """A dim outer cortex with a correct registration is partial volume (verdict: boundary on the MRI boundary); a
    PET 4 % too large or shifted 1.5 mm is reported as a scale / registration error. Through check_record."""
    from scipy import ndimage
    from modalities.pet import boundary
    f, truth, centre = _folded_brain(str(tmp_path))
    mat = opj(str(tmp_path), 'sub-M', 'ses-1', 'pet', 'orig', 'matrices')
    os.makedirs(mat)

    def record(name, matrix, t):
        p = opj(mat, name + '.mat')
        ants.write_transform(ants.create_ants_transform(transform_type='AffineTransform', dimension=3, matrix=matrix,
                                                        translation=t, center=centre), p)
        rec = opj(mat, name + '_desc-registration_pet.json')
        transfo.write_record(rec, dict(reference=f['pet'], targets=dict(acpc=dict(
            image=f['t1'], mask=f['mask'], to=[[p, False]], **{'from': [[p, True]]}))))
        return boundary.check_record(rec, 3.0, str(tmp_path))

    ok = record('S1', np.eye(3), [0, 0, 0])
    assert ok['model'] == 'grey / white matter'
    assert ok['fwhm_estimated_mm'] == pytest.approx(3.0, abs=0.5)        # resolution of the phantom PET
    assert abs(ok['dilation_mm']) < 0.3 and all(abs(v) < 0.45 for v in ok['shifts_mm'].values())
    assert ok['verdict'].startswith('PET boundary on the MRI boundary')
    assert os.path.isfile(opj(str(tmp_path), 'QC', 'pet_boundary', 'S1_desc-registration_boundary.png'))
    assert os.path.isfile(opj(str(tmp_path), 'sub-M', 'ses-1', 'pet', 'QC', 'S1_desc-registration_boundary.png'))
    # partial volume alone: with the correct registration the gyral crowns are far below the sulcal banks
    pet_a = ants.resample_image_to_target(ants.image_read(f['pet']), ants.image_read(f['mask'])).numpy()
    depth = ndimage.distance_transform_edt(truth['brain'], sampling=truth['sp'])
    crown, bank = truth['gm'] & (depth < 1.0), truth['gm'] & (depth > 3) & (depth < 6)
    assert pet_a[crown].mean() / pet_a[bank].mean() < 0.75

    assert set(ok['clean_sectors']) >= {'+y', '-y', '+z', '-z'}             # no muscle there
    big = record('S2', np.eye(3) * 0.96, [0, 0, 0])      # samples the PET 4 % closer to the centre: PET 4 % larger
    assert big['dilation_mm'] > 0.5 and 'larger' in big['verdict']
    moved = record('S3', np.eye(3), [0, 1.5, 0])         # PET content 1.5 mm towards -y
    assert moved['translation_mm']['y'] < -0.7 and 'along y' in moved['verdict']


def test_boundary_check_leaves_out_sectors_with_extracerebral_uptake(tmp_path):
    """Muscles beside every sector but the dorsal one and a 6.5-mm PET, exact registration: the uptake pulls the
    PET edge about 1 mm inward (a false 'PET smaller than the MRI'); those sectors are left out and the verdict is
    'not measurable', with the dorsal shift near 0."""
    from modalities.pet import boundary
    f, _, centre = _folded_brain(str(tmp_path), fwhm=6.5, gap_mm=2.5, muscles='all')
    p = str(tmp_path / 'id.mat')
    ants.write_transform(ants.create_ants_transform(transform_type='AffineTransform', dimension=3, matrix=np.eye(3),
                                                    translation=[0, 0, 0], center=centre), p)
    res = boundary.check(f['t1'], f['mask'], f['pet'], {'to': [[p, False]], 'from': [[p, True]]},
                         str(tmp_path / 'S_'), 6.5, anat=dict(gray_wb=f['gm'], white_full=f['wm']))
    assert min(res['shifts_mm'][s] for s in ('+x', '-x', '+y', '-y', '-z')) < -0.7      # the bias, left out
    assert res['clean_sectors'] == ['+z'] and abs(res['shifts_mm']['+z']) < 0.3
    assert res['dilation_mm'] is None and res['verdict'].startswith('not measurable')


def test_pvc_resolution_auto_reads_the_estimate_and_keeps_the_filter(tmp_path):
    rec = str(tmp_path / 'S_desc-registration_pet.json')
    diary = str(tmp_path / 'diary.txt')
    open(diary, 'w').close()
    assert loop2.pvc_resolution('auto', rec, [{}], diary) is None                 # no estimate: skipped
    with open(rec.replace('_pet.json', '_boundary.json'), 'w') as f:
        json.dump(dict(fwhm_estimated_mm=3.5), f)
    assert loop2.pvc_resolution('auto', rec, [{}], diary) == 3.5
    filt = [{'ReconFilterType': 'Gaussian', 'ReconFilterSize': 5.0}]
    assert loop2.pvc_resolution('auto', rec, filt, diary) == 5.0                  # never below the filter
    assert loop2.pvc_resolution(4.0, rec, filt, diary) == 4.0                     # a value given is kept
