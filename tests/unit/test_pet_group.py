"""Group analysis of PET images (Statistics/Group_PET) on synthetic BIDS studies."""
import json
import os

import numpy as np
import pandas as pd
import nibabel as nib
import pytest
from scipy import stats

from Statistics.Group_PET import _glm
from Statistics.Group_PET import pet_group as G
from Statistics.Group_PET import pet_roi_group as R

SHAPE = (24, 24, 24)
VOX = 2.0
AFF = np.diag([VOX, VOX, VOX, 1.0])
AFF[:3, 3] = -VOX * (np.array(SHAPE) - 1) / 2
TEMPLATE = 'EDNiX'


def _phantom():
    """Brain sphere: cortex-like shell (label 1, 'blob' part label 3), white-matter core (2), 'cerebellum' (7),
    CSF ventricle (low uptake) in the middle."""
    ijk = np.indices(SHAPE).transpose(1, 2, 3, 0) - (np.array(SHAPE) - 1) / 2
    r = np.linalg.norm(ijk, axis=-1)
    brain = r <= 10
    lab = np.zeros(SHAPE, np.int16)
    lab[brain & (r > 6.5)] = 1
    lab[brain & (r <= 6.5)] = 2
    lab[brain & (ijk[..., 2] < -6)] = 7
    blob = (np.linalg.norm(ijk - np.array([0, 8, 2]), axis=-1) <= 3) & (lab == 1)
    lab[blob] = 3
    csf = r <= 2.5
    lab[csf] = 0
    return brain, lab, blob, csf


def _image(rng, lab, csf, scale, effect=0.0):
    uptake = np.zeros(SHAPE, np.float32)
    uptake[lab == 1] = 8.0
    uptake[lab == 3] = 8.0 * (1 - effect)
    uptake[lab == 2] = 4.5
    uptake[lab == 7] = 7.0
    uptake[csf] = 0.6
    img = scale * (uptake + rng.normal(0, 0.35, SHAPE).astype(np.float32))
    img[np.linalg.norm(np.indices(SHAPE).transpose(1, 2, 3, 0) - (np.array(SHAPE) - 1) / 2, axis=-1) > 11] = \
        0.2 * scale                                                        # scalp / background, outside the brain
    return img


def make_study(root, n_per_group=8, effect=0.15, seed=0, with_roi=True, dynamic=()):
    """dynamic: subject numbers exported as a dynamic run (desc-SUV-mean) whose file names carry no ses- entity
    (sub-X_trc-FDG_..., as BIDS_Tremblay), the session being only the folder name."""
    rng = np.random.default_rng(seed)
    brain, lab, blob, csf = _phantom()
    rows = []
    k = 0
    for g in ('patient', 'control'):
        for _ in range(n_per_group):
            k += 1
            sub, ses = f'{k:02d}', '1'
            scale = rng.uniform(0.6, 1.6)                                  # SUV differs between subjects
            img = _image(rng, lab, csf, scale, effect if g == 'patient' else 0.0)
            d = root / f'sub-{sub}' / f'ses-{ses}' / 'pet'
            t = d / 'templates' / TEMPLATE
            (t / 'labels').mkdir(parents=True)
            (t / 'masks').mkdir(parents=True)
            dyn = k in dynamic
            base = f'sub-{sub}_trc-FDG' if dyn else f'sub-{sub}_ses-{ses}_trc-FDG'
            dsc = 'SUV-mean' if dyn else 'SUV'
            nib.save(nib.Nifti1Image(img, AFF), str(t / f'{base}_space-{TEMPLATE}_desc-{dsc}_pet.nii.gz'))
            nib.save(nib.Nifti1Image(img * brain, AFF), str(t / f'{base}_space-{TEMPLATE}_desc-{dsc}-SS_pet.nii.gz'))
            nib.save(nib.Nifti1Image(lab, AFF), str(t / 'labels' / 'Phantom_seg-PH_dseg.nii.gz'))
            nib.save(nib.Nifti1Image(brain.astype(np.uint8), AFF), str(t / 'masks' / 'Phantom_mask.nii.gz'))
            if with_roi:
                roi = d / 'orig' / 'postprocessing' / 'ROI_analysis'
                roi.mkdir(parents=True)
                tab = []
                for l, name in ((1, 'cortex'), (2, 'white'), (3, 'blob'), (7, 'cerebellum'), (9, 'tiny')):
                    v = img[lab == l]
                    if l == 9:
                        n = 4 if k == 3 else 30                           # too small in one subject only
                        v = rng.normal(8 * scale, 0.3, n)
                    tab.append(dict(label=l, region=name, n_voxels=v.size, volume_mm3=v.size * 8.0, n_valid=v.size,
                                    mean=float(v.mean()), sd=float(v.std(ddof=1)), median=float(np.median(v))))
                pd.DataFrame(tab).to_csv(roi / f'{base}_seg-PH-lvl0_desc-{dsc}_ROIs.tsv', sep='\t', index=False)
            rows.append(dict(participant_id=f'sub-{sub}', group=g, age=float(rng.uniform(20, 60))))
    pd.DataFrame(rows).to_csv(root / 'participants.tsv', sep='\t', index=False)
    return dict(brain=brain, lab=lab, blob=blob, csf=csf)


LABELS = {'PH': pd.DataFrame(dict(label_id=[1, 2, 3, 7], region_name=['cortex', 'white', 'blob', 'cerebellum']))}


# ------------------------------------------------------------------------------------------------- GLM
def test_glm_matches_scipy_two_sample_and_regression():
    rng = np.random.default_rng(1)
    tab = pd.DataFrame(dict(group=['a'] * 7 + ['b'] * 9, age=rng.uniform(20, 70, 16)))
    Y = rng.normal(size=(16, 50))
    Y[:7] += 0.8
    X, c, _, _ = _glm.design_matrix(tab, 'two-sample', 'group', ['a', 'b'])
    cope, t, _ = _glm.GLM(X, c).fit(Y)
    ref = stats.ttest_ind(Y[:7], Y[7:], axis=0)
    assert np.allclose(t, ref.statistic) and np.allclose(cope, Y[:7].mean(0) - Y[7:].mean(0))
    X, c, _, _ = _glm.design_matrix(tab, 'regression', covariates=['age'])
    _, t, _ = _glm.GLM(X, c).fit(Y)
    r = np.array([stats.pearsonr(tab.age, Y[:, v])[0] for v in range(50)])
    assert np.allclose(t, r * np.sqrt(14 / (1 - r ** 2)))
    with pytest.raises(ValueError):
        _glm.design_matrix(tab.assign(age2=tab.age * 2), 'regression', covariates=['age', 'age2'])


def test_t_to_z_is_exact_and_finite_in_the_tails():
    t = np.array([-3.0, 0.0, 2.5, 40.0, 1e3])
    z = _glm.t_to_z(t, 10)
    assert np.allclose(z[:3], np.sign(t[:3]) * stats.norm.isf(stats.t.sf(np.abs(t[:3]), 10)))
    assert np.all(np.isfinite(z)) and z[-1] > z[-2] > 5


def test_bh_fdr_matches_statsmodels():
    sm = pytest.importorskip('statsmodels.stats.multitest')
    p = np.random.default_rng(0).uniform(size=40) ** 3
    assert np.allclose(_glm.bh_fdr(p), sm.multipletests(p, method='fdr_bh')[1])


def test_critical_values_of_the_null_distribution():
    null = np.arange(1, 101)                    # 100 null maxima
    k = _glm.critical_value(null, 0.05, integer=True)
    assert k == 96 and _glm.fwe_p([k], null)[0] <= 0.05 and _glm.fwe_p([k - 1], null)[0] > 0.05


@pytest.mark.slow
def test_permutation_controls_the_family_wise_error():
    """Null data (two groups, smooth noise): clusters and voxels FWE-significant in about alpha of the datasets."""
    from scipy import ndimage
    rng = np.random.default_rng(5)
    mask = np.ones((14, 14, 14), bool)
    tab = pd.DataFrame(dict(group=['a'] * 6 + ['b'] * 6))
    X, c, _, _ = _glm.design_matrix(tab, 'two-sample', 'group', ['a', 'b'])
    glm = _glm.GLM(X, c)
    hits_c = hits_v = 0
    n_sim = 60
    for s in range(n_sim):
        Y = np.stack([ndimage.gaussian_filter(rng.normal(size=mask.shape), 1.2).ravel() for _ in range(12)])
        perm = _glm.permutation_test(glm, Y, mask, voxel_p=0.01, n_perm=200, seed=s)
        thr = perm['threshold']
        vol = perm['t'].reshape(mask.shape)
        obs = _glm._cluster_sizes_signed(vol, thr)
        hits_c += _glm.fwe_p([obs], perm['max_cluster'])[0] <= 0.05
        hits_v += _glm.fwe_p([np.abs(perm['t']).max()], perm['max_t'])[0] <= 0.05
    assert hits_c / n_sim < 0.15 and hits_v / n_sim < 0.15


# ------------------------------------------------------------------------------------------------- pieces
def test_smoothing_inside_the_mask_keeps_the_outside_signal_out():
    mask = np.zeros((20, 20, 20), bool)
    mask[5:15, 5:15, 5:15] = True
    d = np.where(mask, 1.0, 50.0)[..., None].astype(np.float32)     # very high signal around the mask (scalp)
    out = G.smooth_in_mask(d, mask, 6.0, (2, 2, 2))
    assert np.allclose(out[mask], 1.0, atol=1e-5) and np.all(out[~mask] == 0)


def test_reference_normalisation_removes_the_global_scale(tmp_path):
    make_study(tmp_path, n_per_group=3, effect=0.0, with_roi=False)
    imgs = G.collect_images(str(tmp_path), TEMPLATE)
    assert imgs.desc.iloc[0] == 'SUV' and len(imgs) == 6 and not imgs.image.str.contains('-SS_').any()
    data, ref = G.load_on_grid(list(imgs.image))
    files = G.template_files(imgs)
    brain = G.mask_on_grid(files['brain'], ref)
    cer = G.reference_mask({'type': 'atlas', 'atlas': 'PH', 'level': 0, 'labels': [7], 'name': 'Cer'}, ref, files)
    norm, _ = G.normalise(data, 'reference', brain, cer)
    ctx = (G.labels_on_grid(files['atlases']['PH'], 0, ref) == 1)
    suvr = [norm[..., k][ctx].mean() for k in range(6)]
    assert np.std(suvr) < 0.01 and abs(np.mean(suvr) - 8 / 7) < 0.02
    assert np.std([data[..., k][ctx].mean() for k in range(6)]) > 0.5      # SUV itself varied


def test_analysis_mask_excludes_csf_and_background(tmp_path):
    ph = make_study(tmp_path, n_per_group=3, effect=0.0, with_roi=False)
    imgs = G.collect_images(str(tmp_path), TEMPLATE)
    data, ref = G.load_on_grid(list(imgs.image))
    brain = G.mask_on_grid(G.template_files(imgs)['brain'], ref)
    norm, _ = G.normalise(data, 'global', brain)
    m, _ = G.analysis_mask(norm, brain, 0.8)
    assert not (m & ph['csf']).any() and not (m & ~ph['brain']).any()
    assert (m & (ph['lab'] == 1)).sum() > 0.95 * (ph['lab'] == 1).sum()


def test_images_on_different_grids_are_resampled(tmp_path):
    make_study(tmp_path, n_per_group=2, effect=0.0, with_roi=False)
    imgs = G.collect_images(str(tmp_path), TEMPLATE)
    f = imgs.image.iloc[1]
    im = nib.load(f)
    orig = np.asarray(im.dataobj).copy()
    aff = im.affine.copy()
    aff[:3, :3] *= 0.5                                         # finer voxels (another scanner)
    from nibabel.processing import resample_from_to
    rs = resample_from_to(im, (tuple(2 * s for s in im.shape), aff), order=1)
    nib.save(rs, f)
    data, _ = G.load_on_grid(list(imgs.image))
    assert data.shape == SHAPE + (4,)
    assert np.corrcoef(data[..., 1].ravel(), orig.ravel())[0, 1] > 0.95


# ------------------------------------------------------------------------------------------------- end to end
def test_two_group_analysis_finds_the_hypometabolism(tmp_path):
    ph = make_study(tmp_path, n_per_group=8, effect=0.15)
    res = G.run_group_analysis(str(tmp_path), TEMPLATE, 'two-sample', group_col='group',
                               groups=['patient', 'control'], covariates=['age'],
                               normalization={'type': 'atlas', 'atlas': 'PH', 'level': 0, 'labels': [7],
                                              'name': 'Cer'},
                               label_tables=LABELS, fwhm=4.0, voxel_p=0.001, n_perm=500, seed=1)
    out = res['output_dir']
    assert out.endswith('Grp_PET_two-sample_desc-SUV_norm-Cer')
    cl = res['clusters']
    sig = cl[cl.significant]
    assert len(sig) >= 1 and (sig.sign == '-').all()                  # hypometabolism only
    z = np.asarray(nib.load(res['maps']['zstat_thresholded']).dataobj)
    assert (z < 0).sum() > 0 and ((z != 0) & ~ph['blob']).sum() < 0.5 * (z != 0).sum() + 30
    assert sig.iloc[0].peak_region == 'blob' and sig.iloc[0].p_FWE_cluster <= 0.05
    pc = np.asarray(nib.load(res['maps']['percent_difference']).dataobj)
    assert -20 < np.median(pc[ph['blob'] & (pc != 0)]) < -9                      # 15 % effect, smoothed
    s = json.load(open(os.path.join(out, 'patient-vs-control_summary.json')))
    assert s['n'] == 16 and s['df'] == 13 and s['normalising_value_test'] is not None
    ins = pd.read_csv(os.path.join(out, 'inputs.tsv'), sep='\t')
    assert len(ins) == 16 and ins.normalising_value.notna().all()
    for f in ('mask_group.nii.gz', 'inputs_4D.nii.gz', 'group-patient_mean.nii.gz',
              'patient-vs-control_logp-FWE-voxel.nii.gz', 'patient-vs-control_clusters.tsv'):
        assert os.path.isfile(os.path.join(out, f)), f


def test_null_study_gives_no_significant_cluster(tmp_path):
    make_study(tmp_path, n_per_group=6, effect=0.0, seed=3, with_roi=False)
    res = G.run_group_analysis(str(tmp_path), TEMPLATE, 'two-sample', group_col='group',
                               groups=['patient', 'control'], fwhm=4.0, n_perm=300, seed=2)
    assert res['summary']['n_significant'] == 0
    assert not np.any(np.asarray(nib.load(res['maps']['zstat_thresholded']).dataobj))


def test_paired_sessions_and_selection(tmp_path):
    make_study(tmp_path, n_per_group=4, effect=0.0, with_roi=False)
    rng = np.random.default_rng(0)
    _, lab, _, csf = _phantom()
    for sub in [f'{k:02d}' for k in range(1, 9)]:                     # second session with a 20 % blob decrease
        t = tmp_path / f'sub-{sub}' / 'ses-2' / 'pet' / 'templates' / TEMPLATE
        (t / 'labels').mkdir(parents=True)
        (t / 'masks').mkdir(parents=True)
        img = _image(rng, lab, csf, rng.uniform(0.6, 1.6), 0.2)
        nib.save(nib.Nifti1Image(img, AFF), str(t / f'sub-{sub}_ses-2_trc-FDG_space-{TEMPLATE}_desc-SUV_pet.nii.gz'))
    res = G.run_group_analysis(str(tmp_path), TEMPLATE, 'paired', paired_levels=['2', '1'], fwhm=4.0,
                               n_perm=1000, list_to_remove=[('08', '1')])
    s = res['summary']
    assert s['n'] == 7 and s['n_perm'] == 128 and s['exhaustive']          # 2^7 sign flips
    assert res['contrast'] == 'paired-difference'
    with pytest.raises(ValueError):
        G.run_group_analysis(str(tmp_path), TEMPLATE, 'two-sample', group_col='group',
                             groups=['patient', 'control'], aggregate='session')


# ------------------------------------------------------------------------------------------------- regions
def test_roi_analysis_filters_small_regions_and_normalises_to_the_reference(tmp_path):
    make_study(tmp_path, n_per_group=8, effect=0.15)
    res = R.run_roi_analysis(str(tmp_path), 'PH', 0, 'SUV', 'two-sample', group_col='group',
                             groups=['patient', 'control'], normalization=[7], min_voxels=10)
    rep = res['report'].set_index('label')
    assert not rep.loc[9, 'kept'] and rep.loc[[1, 2, 3, 7], 'kept'].all()   # 4 voxels in one subject
    st = res['stats'].set_index('region')
    assert st.loc['blob', 'p_FDR'] < 0.001 and st.loc['blob', 'cope'] < 0
    assert -18 < st.loc['blob', 'percent_difference'] < -12
    assert st.loc['cerebellum', 'cope'] == pytest.approx(0, abs=1e-9)        # the reference is 1 by definition
    assert st.loc['cortex', 'p_FDR'] > 0.05
    vals = res['values']
    assert np.allclose(vals[7], 1.0)
    out = res['output_dir']
    for f in ('values_raw.tsv', 'values_normalised.tsv', 'regions_voxel_counts.tsv',
              'patient-vs-control_regions.tsv'):
        assert os.path.isfile(os.path.join(out, f)), f


def test_roi_reference_is_the_voxel_weighted_mean_of_its_labels():
    ses = pd.DataFrame(dict(subject=['1', '1', '1'], session=['1'] * 3, label=[6, 1006, 5],
                            region=['L', 'R', 'x'], n_valid=[100, 300, 50], n_voxels=[100, 300, 50],
                            volume_mm3=[1, 1, 1], mean=[2.0, 4.0, 7.0]))
    _, ref = R.normalise_regions(ses, [6, 1006])
    assert ref.normalising_value.iloc[0] == pytest.approx((100 * 2 + 300 * 4) / 400)


def test_surface_rendering_uses_the_ribbon_projection_and_skips_empty_maps(tmp_path, monkeypatch):
    from Tools import niftotoWBsurface
    from Statistics.Group_PET import pet_surface
    calls = []
    monkeypatch.setattr(niftotoWBsurface, 'nifti_to_surface_wb', lambda **k: calls.append(k))
    z = np.zeros(SHAPE, np.float32)
    z[5:8, 5:8, 5:8] = -4
    nib.save(nib.Nifti1Image(z, AFF), str(tmp_path / 'c_zstat_thresholded.nii.gz'))
    nib.save(nib.Nifti1Image(np.zeros(SHAPE, np.float32), AFF), str(tmp_path / 'c_zstat.nii.gz'))
    res = {'output_dir': str(tmp_path), 'maps': {'zstat_thresholded': str(tmp_path / 'c_zstat_thresholded.nii.gz'),
                                                 'zstat': str(tmp_path / 'c_zstat.nii.gz')}}
    done = pet_surface.render_pet_surfaces(res, 'Dog', '/atlas', 'atlas/mammals/carnivorans/Dog', 'wb ')
    assert list(done) == ['zstat_thresholded'] and len(calls) == 1
    k = calls[0]
    assert k['projection'] == 'ribbon' and k['pos_max'] == pytest.approx(4) and k['neg_max'] == pytest.approx(-4)
    assert k['output_dir'] == str(tmp_path / 'surfaces') and k['display_zero'] is False


# ------------------------------------------------------------------------------------------------- naming
def test_dynamic_naming_and_session_from_the_folder(tmp_path, capsys):
    """desc-SUV finds the desc-SUV-mean images of dynamic runs (and the reverse); the session comes from the
    ses- folder when the file name has no ses- entity."""
    make_study(tmp_path, n_per_group=3, dynamic={1, 2, 3, 4, 5, 6})
    imgs = G.collect_images(str(tmp_path), TEMPLATE, 'SUV')
    assert len(imgs) == 6 and set(imgs.desc) == {'SUV-mean'} and set(imgs.session) == {'1'}
    assert 'resolved to desc-SUV-mean' in capsys.readouterr().out
    long = R.collect_roi_tables(str(tmp_path), 'PH', 0, 'SUV')
    assert set(long.session) == {'1'} and long.subject.nunique() == 6
    # mixed study: static and dynamic runs side by side, each session uses the one it has
    mixed = tmp_path / 'mixed'
    make_study(mixed, n_per_group=3, dynamic={2, 5})
    for d in ('SUV', 'SUV-mean'):
        imgs = G.collect_images(str(mixed), TEMPLATE, d)
        assert len(imgs) == 6 and imgs.groupby('desc').size().to_dict() == {'SUV': 4, 'SUV-mean': 2}
    with pytest.raises(FileNotFoundError, match='desc-bratio-mean'):
        G.collect_images(str(mixed), TEMPLATE, 'bratio')
    res = G.run_group_analysis(str(mixed), TEMPLATE, 'two-sample', desc='SUV', group_col='group',
                               groups=['patient', 'control'], normalization='global', n_perm=50,
                               output_dir=str(mixed / 'out'))
    assert json.load(open(os.path.join(res['output_dir'], 'patient-vs-control_summary.json')))['desc'] == 'SUV'


# ------------------------------------------------------------------------------------------------- subject tables
def test_no_covariate_needs_no_participants_table(tmp_path):
    """One-sample without covariate and two-sample with group_assignment: participants.tsv is not read."""
    make_study(tmp_path, n_per_group=4, effect=0.3)
    (tmp_path / 'participants.tsv').write_text('broken\tfile\nwithout\tids\n')
    res = G.run_group_analysis(str(tmp_path), TEMPLATE, 'one-sample', desc='SUV', covariates=None,
                               normalization='global', n_perm=50, output_dir=str(tmp_path / 'one'))
    assert res['summary']['n'] == 8
    groups = {'patient': ['sub-01', '02', '03', '04'], 'control': ['05', '06', '07', '08']}
    res = G.run_group_analysis(str(tmp_path), TEMPLATE, 'two-sample', desc='SUV', group_assignment=groups,
                               normalization=dict(type='atlas', atlas='PH', level=0, labels=[7], name='Cer'),
                               label_tables=LABELS, n_perm=200, output_dir=str(tmp_path / 'two'))
    assert res['contrast'] == 'patient-vs-control' and res['summary']['n'] == 8
    roi = R.run_roi_analysis(str(tmp_path), 'PH', 0, 'SUV', 'two-sample', group_assignment=groups,
                             normalization=[7], output_dir=str(tmp_path / 'roi'))
    st = roi['stats'].set_index('region')
    assert st.loc['blob', 'cope'] < 0 and st.loc['blob', 'p_FDR'] < 0.05
    with pytest.raises(ValueError, match='no image for subject'):
        G.assign_groups(G.collect_images(str(tmp_path), TEMPLATE, 'SUV'), {'a': ['99'], 'b': ['01']})
    with pytest.raises(KeyError, match='no ID column.*columns found'):         # a covariate needs the table
        G.run_group_analysis(str(tmp_path), TEMPLATE, 'regression', desc='SUV', covariates=['age'],
                             output_dir=str(tmp_path / 'reg'))


def test_participants_table_with_another_id_column(tmp_path):
    make_study(tmp_path, n_per_group=4)
    p = pd.read_csv(tmp_path / 'participants.tsv', sep='\t')
    p = p.rename(columns={'participant_id': 'Subject'})
    p['Subject'] = p['Subject'].str.replace('sub-', '')
    p.loc[0, 'age'] = np.nan
    p.to_csv(tmp_path / 'participants.tsv', sep='\t', index=False, na_rep='n/a')
    imgs = G.collect_images(str(tmp_path), TEMPLATE, 'SUV')
    t = G.subject_table(imgs, str(tmp_path), ['age', 'group'])
    assert t.age.dtype.kind == 'f' and t.age.isna().sum() == 1 and set(t.group) == {'patient', 'control'}
    with pytest.raises(KeyError, match="available: .*'age'"):
        G.subject_table(imgs, str(tmp_path), ['weight'])


def test_resampling_does_not_depend_on_the_nilearn_version(monkeypatch):
    """Resampling through nibabel.processing: nilearn's resample_to_img (keywords changed in 0.11) is not used.
    Trilinear resampling of a linear ramp is exact; nearest keeps the labels."""
    import nilearn.image

    def _fail(*a, **k):
        raise AssertionError('nilearn.image.resample_to_img called')
    monkeypatch.setattr(nilearn.image, 'resample_to_img', _fail)
    ijk = np.indices(SHAPE).astype(np.float32)
    ramp = nib.Nifti1Image((2 * ijk[0] + ijk[1] - ijk[2]).astype(np.float32), AFF)
    aff2 = AFF.copy()
    aff2[:3, 3] += VOX / 2                                                 # half-voxel shift
    ref = nib.Nifti1Image(np.zeros(SHAPE, np.float32), aff2)
    out = np.asarray(G.resample_to_ref(ramp, ref, order=1).dataobj)
    inner = (slice(1, -2),) * 3
    assert np.allclose(out[inner], (2 * ijk[0] + ijk[1] - ijk[2] + 0.5 * (2 + 1 - 1))[inner], atol=1e-4)
    lab = nib.Nifti1Image(np.rint(ijk[0] / 6).astype(np.int16), AFF)
    labs = np.unique(np.asarray(G.resample_to_ref(lab, ref, order=0).dataobj))
    assert set(labs) <= set(np.unique(np.asarray(lab.dataobj)))


# ------------------------------------------------------------------------------------------------- group mean
@pytest.mark.parametrize('n', [1, 2])
def test_group_mean_for_small_groups_with_surfaces(tmp_path, monkeypatch, n):
    """n = 1 or 2 (no GLM possible): mean / SD / CV maps of the ratio image, regional means, surfaces of the
    mean; design names with '_' are accepted by the GLM functions."""
    from Tools import niftotoWBsurface
    from Statistics.Group_PET import pet_surface
    ph = make_study(tmp_path, n_per_group=1, effect=0.0, dynamic={1, 2})     # 2 subjects, desc-SUV-mean
    keep = [('01', '1')] if n == 1 else None
    res = G.run_group_mean(str(tmp_path), TEMPLATE, 'SUV', list_to_keep=keep, output_dir=str(tmp_path / 'mean'))
    assert res['n'] == n and ('sd_all' in res['maps']) == (n > 1) and ('cv_all' in res['maps']) == (n > 1)
    m = np.asarray(nib.load(res['maps']['mean_all']).dataobj)
    mask = np.asarray(nib.load(os.path.join(res['output_dir'], 'mask_group.nii.gz')).dataobj) > 0
    assert not (mask & ph['csf']).any() and np.all(m[~mask] == 0)
    assert m[mask & (ph['lab'] == 1)].mean() > m[mask & (ph['lab'] == 7)].mean()   # cortex 8 > cerebellum 7
    calls = []
    monkeypatch.setattr(niftotoWBsurface, 'nifti_to_surface_wb', lambda **k: calls.append(k))
    done = pet_surface.render_pet_surfaces(res, 'Macaque', '/atlas', 'atlas/mammals/primates/Macaque', 'wb ',
                                           maps=['mean_all', 'cv_all'])
    assert 'mean_all' in done and calls[0]['palette_name'] == 'videen_style' and not calls[0]['display_negative']
    roi = R.run_roi_mean(str(tmp_path), 'PH', 0, 'SUV', list_to_keep=keep, output_dir=str(tmp_path / 'roi'))
    st = roi['stats'].set_index('region')
    assert set(st.n) == {n} and 'tiny' in st.index and st.loc['cortex', 'mean'] > st.loc['white', 'mean']
    assert os.path.isfile(os.path.join(roi['output_dir'], 'regions_mean.tsv'))
    assert G._design_name('one_sample') == 'one-sample'
    with pytest.raises(ValueError, match='expected one of'):
        G._design_name('mean')
