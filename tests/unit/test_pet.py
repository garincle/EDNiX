"""PET / CT pipelines: helpers, transforms, motion correction, registration and an end-to-end run on a phantom."""
import os
import json
import math

import numpy as np
import nibabel as nib
import pandas as pd
import pytest

import _pet_phantom as P
from modalities.pet import petutils as pu
from modalities.pet import transfo, dftFile, cut_image, run_cut, MoCo, coreg2anat, loop1, loop2, studyT
from modalities.pet.turkutools._nifti import TurkuWorkspace
from atlases import getlabels

opj = os.path.join


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
                MoCo_transform='Rigid', coreg_metric='mattes', type_of_transform_template='SyN', studyT_transform='SyN',
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
    assert rec['targets']['acpc']['to'] == crec['targets']['acpc']['to']        # PET follows the CT


@pytest.mark.slow
def test_pet_study_template_without_anatomy_3v2(tmp_path, monkeypatch):
    from modalities.pet import _0_Pipeline_launcher as L
    main, bids, tdir, s = _study(tmp_path, monkeypatch, subjects=('01', '02'), anat=False)
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1'), ('02', '1')]), '*_trc-FDG_*pet.nii.gz',
                     'FDG', [''], [''], '', 'T1w', 'T1w', False, False, IhaveanANAT=False,
                     creat_study_template_pet=True, studyT_transform='Affine', crop_for_MoCo='auto',
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
    """Raw PET in the scanner space of a sphinx animal (AHF declared HFS): after the relabelling of step 1 it lies in
    the space of the reoriented anatomy, and the AC-PC matrix of the anatomy brings it to AC-PC."""
    from modalities.pet import _0_Pipeline_launcher as L
    from Tools import get_orientation, getpath
    main, bids, tdir, s = _study(tmp_path, monkeypatch)
    dp, A0 = s['01']['dp'], s['01']['A']
    raw, _, info = P.dynamic_pet(opj(dp, 'pet', 'sub-01_ses-1_trc-FDG_pet.nii.gz'), np.eye(3), np.zeros(3))
    probe = str(tmp_path / 'probe.nii.gz')
    P.save(np.zeros((4, 4, 4)), A0, probe)
    new = get_orientation.getreal('HFS', 'AHF', P._afni_code(A0))
    P.emulate('3drefit -orient ' + new + ' ' + probe)
    Pm = nib.load(probe).affine @ np.linalg.inv(A0)               # relabelling of step 1 (world map)
    img = nib.load(raw)
    nib.save(nib.Nifti1Image(np.asarray(img.dataobj), np.linalg.inv(Pm) @ A0), raw)   # scanner space
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     ['AHF'], ['HFS'], '', 'T1w', 'T1w', False, False, anat_pet_same_space=True,
                     selected_atlases=[['TEST', 0]], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    paths = getpath.pet(dp, 'EDNiX')
    rec = transfo.read_record(coreg2anat.record_path(paths, 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '1v1' and rec['targets']['acpc']['to'][0][0].endswith('acpc_0GenericAffine.mat')
    assert np.allclose(nib.load(rec['reference']).affine, A0, atol=1e-4)
    atlas_pet = nib.load(rec['pet_space']['atlases']['TEST']).get_fdata()[..., 0]
    truth = P.tissues(P.world_points((40,) * 3, A0)).reshape((40,) * 3)
    dice = 2 * ((atlas_pet == 7) & (truth == 6)).sum() / ((atlas_pet == 7).sum() + (truth == 6).sum())
    assert dice > 0.85, dice


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
                     type_of_transform_template='Affine', selected_atlases=[['TEST', 'all']],
                     atlas_followers=[['TEST'], ['txt'], [2], [1]])
    paths = getpath.pet(s['01']['dp'], 'EDNiX')
    rec = transfo.read_record(coreg2anat.record_path(paths, 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '3v1' and list(rec['targets']) == ['EDNiX']
    atlas_pet = nib.load(rec['pet_space']['atlases']['TEST']).get_fdata()[..., 0]
    q = P.world_points((40,) * 3, s['01']['A']) @ s['01']['R'].T + s['01']['t']
    truth = P.tissues(q).reshape((40,) * 3)
    brain = 2 * ((atlas_pet > 0) & (truth >= 3)).sum() / ((atlas_pet > 0).sum() + (truth >= 3).sum())
    assert brain > 0.85, brain
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
                     type_of_transform_template='Affine', atlas_followers=[['TEST'], ['txt'], [2], [1]])
    L.preprocess_pet([], main, bids, 'Testmonkey', P.bids_table([('01', '1')]), '*_trc-FDG_*pet.nii.gz', 'FDG',
                     [''], [''], '', 'T1w', 'T1w', False, False, IhaveanANAT=False, IhaveaCT=True,
                     selected_atlases=[['TEST', 'all']], atlas_followers=[['TEST'], ['txt'], [2], [1]])
    rec = transfo.read_record(coreg2anat.record_path(getpath.pet(s['01']['dp'], 'EDNiX'), 'sub-01_ses-1_trc-FDG'))
    assert rec['scenario'] == '2v1' and list(rec['targets']) == ['EDNiX'] and rec['pet_space']['atlases']
