"""Figure 6 data path (Exemples/Study_EDNiX/analysis/library/multimodal_profiles.py) on synthetic BIDS trees."""
import os

import numpy as np
import pandas as pd
import nibabel as nib
import pytest

pytest.importorskip('seaborn')
from Exemples.Study_EDNiX.analysis.library import multimodal_profiles as MP
from Plotting import Plot_interspecies_regional_profiles as PIRP

opj = os.path.join
REGIONS = ['Frontal_lobe', 'Temporal_lobe', 'Parietal_lobe', 'Occipital_lobe', 'Insular_lobe']
EXTRA = ['Lateral_ventricules', 'Brainstem']          # not cortical; ventricles excluded from the brain volume
NAMES = REGIONS + EXTRA


def _labels(atlas_lib, sp):
    d = opj(atlas_lib, 'atlas', 'mammals', 'primates', sp, 'label_code')
    os.makedirs(d, exist_ok=True)
    with open(opj(d, 'EDNIxCSCLR_StatsLUT.txt'), 'w') as f:
        for i, n in enumerate(NAMES):
            f.write(f'{10 + i} L_{n} 10 10 10 0\n{1010 + i} R_{n} 10 10 10 0\n')


def _dseg(bids, sub, sizes):
    """4D dseg: level 1 (index 0) one label; level 2 (index 1) the regions with `sizes` voxels per hemisphere."""
    d = opj(bids, f'sub-{sub}', 'ses-1', 'anat', 'native', 'volumes', 'labels')
    os.makedirs(d, exist_ok=True)
    n = 2 * sum(sizes)
    lab = np.zeros((n, 1, 1, 2), np.int32)
    lab[..., 0] = 1
    k = 0
    for i, s in enumerate(sizes):
        for off in (0, 1000):
            lab[k:k + s, 0, 0, 1] = 10 + i + off
            k += s
    nib.save(nib.Nifti1Image(lab, np.eye(4)), opj(d, f'sub-{sub}_seg-EDNIxCSCLR_dseg.nii.gz'))


def _matrix(bids, sub, strength, rng):
    d = opj(bids, f'sub-{sub}', 'ses-1', 'func', 'acpc-func', 'Stats', 'Correl_matrix', 'EDNIxCSCLR', 'correlation')
    os.makedirs(d, exist_ok=True)
    names = [f'{h}_{n}' for h in 'LR' for n in NAMES]
    s = np.array([strength[n] for h in 'LR' for n in NAMES])
    r = np.clip(np.sqrt(np.outer(s, s)) + rng.normal(0, 0.01, (len(s), len(s))), -0.95, 0.95)
    r = (r + r.T) / 2
    np.fill_diagonal(r, 0)
    for run in (0, 1):
        pd.DataFrame(r, index=names, columns=names).to_csv(opj(d, f'EDNIxCSCLR_1_run_{run}_correlation_matrix.csv'))


def _pet(bids, sub, values):
    d = opj(bids, f'sub-{sub}', 'ses-1', PIRP.ROI_SUBDIR)
    os.makedirs(d, exist_ok=True)
    rows = [dict(label=10 + i + off, region=f'{h}_{n}', n_voxels=100, n_valid=100, volume_mm3=800.0,
                 mean=values[n], sd=0.1, median=values[n])
            for i, n in enumerate(NAMES) for h, off in (('L', 0), ('R', 1000))]
    pd.DataFrame(rows).to_csv(opj(d, f'sub-{sub}_ses-1_trc-FDG_seg-EDNIxCSCLR-lvl1_desc-bratio_ROIs.tsv'),
                              sep='\t', index=False)


@pytest.fixture
def cfg(tmp_path):
    rng = np.random.default_rng(0)
    lib = str(tmp_path / 'lib')
    src = {'morphometry': {}, 'FC': {}, 'metabolism': {}}
    for sp, nsub, occ in (('Human', 6, 1.10), ('Macaque', 4, 1.40)):
        _labels(lib, sp)
        bids = str(tmp_path / sp)
        for s in range(nsub):
            _dseg(bids, f'{sp[0]}{s}', [40, 30, 20, 10, 5, 50, 60])
            _matrix(bids, f'{sp[0]}{s}', dict(zip(NAMES, [0.3, 0.2, 0.25, 0.4, 0.1, 0.05, 0.05])), rng)
            v = dict(zip(NAMES, [1.1, 0.95, 1.05, occ, 1.0, 0.5, 0.8]))
            _pet(bids, f'{sp[0]}{s}', {k: x + rng.normal(0, 0.02) for k, x in v.items()})
        for m in src:
            src[m][sp] = [bids]
    return dict(SOURCES=src, ATLAS='EDNIxCSCLR', LEVEL=2, ATLAS_LIB=lib, DESC=('bratio', 'bratio-mean'))


def test_modalities_are_read_at_level_2_with_canonical_names(cfg):
    t = MP.collect_modalities(cfg)
    assert set(t) == {'morphometry', 'FC', 'metabolism'}
    for df in t.values():
        assert set(REGIONS) <= set(df['region'])
    m = t['morphometry']
    h0 = m[(m.subject == 'H0')].set_index('region')['value']
    # ventricles excluded from the brain: 2 x (40+30+20+10+5+60) = 330 voxels
    assert 'Lateral_ventricules' not in h0.index
    assert np.isclose(h0['Frontal_lobe'], 80 / 330 * 100)
    assert m.groupby('species')['subject'].nunique().to_dict() == {'Human': 6, 'Macaque': 4}


def test_node_strength_is_mean_fisher_z_without_diagonal(tmp_path):
    r = np.array([[0, 0.5, 0.2], [0.5, 0, 0.1], [0.2, 0.1, 0]])
    p = tmp_path / 'm.csv'
    pd.DataFrame(r, index=list('abc'), columns=list('abc')).to_csv(p)
    s = MP.node_strength(str(p))
    assert np.isclose(s['a'], np.mean(np.arctanh([0.5, 0.2])))


def test_profiles_and_figure(cfg, tmp_path):
    from Exemples.Study_EDNiX.analysis.library import multimodal_figure as F
    t = MP.collect_modalities(cfg)
    F.REGIONS, F.LEVEL = None, 2
    cortical = ['Frontal_lobe', 'Insular_lobe', 'Temporal_lobe', 'Parietal_lobe', 'Occipital_lobe']
    F.CORTICAL_ONLY = True
    assert F.select_regions(t) == cortical                                # Atlas_nomenclature order, no Brainstem
    F.CORTICAL_ONLY = False
    allr = F.select_regions(t)                                            # every brain region (PR #61), no ventricles
    assert allr[:5] == cortical and 'Brainstem' in allr and 'Lateral_ventricules' not in allr
    prof_b, stats_b, cov_b = F.metabolism_stats(t['metabolism'], REGIONS)
    assert (cov_b.n_missing == 0).all()
    assert list(prof_b.columns) == REGIONS
    occ = stats_b.set_index('region').loc['Occipital_lobe']
    assert occ.q_fdr < 0.05 and occ.difference > 0.2
    prof_c, regions_c = MP.multimodal_profiles(t, REGIONS, F.SPECIES)
    assert regions_c == REGIONS
    z = prof_c['Human']['FC']
    assert np.allclose(z.mean(axis=1), 0) and np.allclose(z.std(axis=1, ddof=1), 1)
    cmc = MP.cross_modal_correlations(prof_c)
    assert len(cmc) == 6
    pet = _pet_stats(F, cfg, tmp_path, REGIONS)
    summ = MP.summary_stats(prof_c, regions_c, F.SPECIES, n=200)
    fig = F.build_figure({}, (0.5, 1.5), pet, prof_c, regions_c, summ, None, {}, (-0.3, 0.3))
    fig.savefig(tmp_path / 'fig.png')
    assert "Figure 6" in F.legend_text(pet, prof_c, regions_c, (0.5, 1.5), cmc, summ, (-0.3, 0.3))


def _pet_human2(tmp_path, rng, nsub=5):
    """Second human FDG dataset, two sessions per subject (eyes open / closed), occipital higher when open."""
    bids = str(tmp_path / 'Human2')
    for s in range(nsub):
        for ses, occ in (('open', 1.16), ('closed', 1.08)):
            d = opj(bids, f'sub-C{s}', f'ses-{ses}', PIRP.ROI_SUBDIR)
            os.makedirs(d, exist_ok=True)
            v = dict(zip(NAMES, [1.1, 0.95, 1.05, occ, 1.0, 0.5, 0.8]))
            rows = [dict(label=10 + i + off, region=f'{h}_{n}', n_voxels=100, n_valid=100, volume_mm3=800.0,
                         mean=v[n] + rng.normal(0, 0.02), sd=0.1, median=v[n])
                    for i, n in enumerate(NAMES) for h, off in (('L', 0), ('R', 1000))]
            pd.DataFrame(rows).to_csv(opj(d, f'sub-C{s}_ses-{ses}_trc-FDG_seg-EDNIxCSCLR-lvl1_desc-bratio-mean_ROIs.tsv'),
                                      sep='\t', index=False)
    return bids


def _pet_stats(F, cfg, tmp_path, regions, monkeypatch=None):
    h2 = _pet_human2(tmp_path, np.random.default_rng(3))
    ds = [('Human_Merida', 'Human', cfg['SOURCES']['metabolism']['Human'][0]), ('Human_Castrillon', 'Human', h2),
          ('Macaque_Tremblay', 'Macaque', cfg['SOURCES']['metabolism']['Macaque'][0])]
    if monkeypatch is not None:
        monkeypatch.setattr(F, 'PET_DATASETS', ds)
        monkeypatch.setattr(F, 'N_PERM', 200)
        return None
    F.PET_DATASETS, F.N_PERM = ds, 200
    return F.pet_dataset_stats(regions)


def test_region_missing_in_one_species_is_dropped_not_the_species(cfg):
    """Regression: a region absent from every human table used to drop all humans (listwise deletion)."""
    from Exemples.Study_EDNiX.analysis.library import multimodal_figure as F
    met = MP.collect_modalities(cfg)['metabolism']
    met = met[~((met.species == 'Human') & (met.region == 'Insular_lobe'))]          # macaque-only region
    met.loc[(met.subject == 'M0') & (met.region == 'Parietal_lobe'), 'value'] = np.nan  # one NaN subject
    prof_b, stats_b, cov_b = F.metabolism_stats(met, REGIONS)
    assert set(prof_b.index.get_level_values(0)) == {'Human', 'Macaque'}
    assert list(prof_b.columns) == ['Frontal_lobe', 'Temporal_lobe', 'Occipital_lobe']
    miss = cov_b.set_index(['region', 'species'])['n_missing']
    assert miss[('Insular_lobe', 'Human')] == 6 and miss[('Parietal_lobe', 'Macaque')] == 1
    with pytest.raises(ValueError, match='species found'):
        PIRP.compare_species(prof_b.loc[['Macaque']], 'Human', 'permutation', 10)


def test_fingerprint_names_are_kept_and_classified_at_level_3():
    names = ['L_Orbital_frontal_cortex_(oFC)', 'R_Prefrontal_cortex_(PFC)', 'L_Thalamus', 'R_Striate_Cx']
    keys = [MP.region_key(n) for n in names]
    assert keys[0] == 'Orbital_frontal_cortex_(oFC)'                     # base name, as the Fig 5 profiles
    assert MP.cortical_regions(keys, 3) == ['Orbital_frontal_cortex_(oFC)', 'Prefrontal_cortex_(PFC)',
                                            'Striate_Cx']


def test_regions_follow_the_lobe_order_of_the_atlas_figure():
    from Exemples.Study_EDNiX.atlases import Atlas_nomenclature as AN
    names = ['R_Visual_striate_cortex', 'L_Somatosensory_cortex', 'Prefrontal_cortex_(PFC)', 'Thalamus',
             'Auditory_cortex_(Superior_temporal)', 'Insula_and_others_in_lateral_sulcus', 'unknown_x',
             'Orbital_frontal_cortex_(oFC)']
    assert [AN.lobe_of(n) for n in AN.atlas_order(names)] == [
        'Frontal_lobe', 'Frontal_lobe', 'Insular_lobe', 'Temporal_lobe', 'Parietal_lobe', 'Occipital_lobe',
        'Diencephalon', None]
    assert AN.atlas_order(names)[:2] == ['Orbital_frontal_cortex_(oFC)', 'Prefrontal_cortex_(PFC)']
    # Figure 6 and the cross-species fingerprints use the same order
    keys = [MP.region_key(n) for n in names]
    assert MP.cortical_regions(keys, 3) == [k for k in PIRP.atlas_sorted(keys) if k not in ('Thalamus', 'unknown_x')]


def test_radar_starts_at_top_counter_clockwise_with_lobe_ring():
    import matplotlib.pyplot as plt
    regions = ['Orbital_frontal_cortex_(oFC)', 'Insula_and_others_in_lateral_sulcus',
               'Auditory_cortex_(Superior_temporal)', 'Visual_striate_cortex']
    ax = plt.figure().add_subplot(projection='polar')
    ang = PIRP.setup_radar(ax, regions, (-1, 1))
    assert np.isclose(ax.get_theta_offset(), np.pi / 2) and ax.get_theta_direction() == 1
    assert len(ax.patches) == 4                                  # one rim segment per region (lobe colour)
    PIRP.radar_labels(ax, ang, regions)
    assert [t.get_ha() for t in ax.get_xticklabels()] == ['center', 'right', 'center', 'left']
    plt.close('all')


def test_bars_can_start_from_zero():
    import matplotlib.pyplot as plt
    idx = pd.MultiIndex.from_tuples([(sp, f'{sp}{i}') for sp in ('Human', 'Macaque') for i in range(3)])
    prof = pd.DataFrame({'Frontal_lobe': [0.9, 0.95, 0.92, 1.1, 1.2, 1.15]}, index=idx)
    st = PIRP.compare_species(prof, 'Human', 'permutation', 50)
    ax = plt.figure().add_subplot()
    PIRP.draw_bars(ax, prof, st, 'Human', 'bratio', bar_from=0)
    bars = [p for p in ax.patches if hasattr(p, 'get_height')]
    assert all(np.isclose(p.get_y(), 0) for p in bars) and ax.get_ylim()[0] == 0
    plt.close('all')


def test_figure6_paper_stats_and_figure(cfg, tmp_path, monkeypatch):
    """Figure 6 through paper_stats / paper_figures: same numbers as the functions of the
    figure module, stats files and figure written (panel A renders stubbed)."""
    from Exemples.Study_EDNiX.analysis.library import multimodal_figure as F
    from Exemples.Study_EDNiX.analysis import paper_stats as pst, paper_figures as pfg
    for k, v in cfg.items():
        monkeypatch.setattr(F, k, v)
    monkeypatch.setattr(F, 'REGIONS', None)
    monkeypatch.setattr(F, 'FC_SOURCE', 'matrix')           # synthetic trees have no fMRI selection
    monkeypatch.setattr(F, 'surface_pngs', lambda diffs=None: ({}, (0.5, 1.5), {}, (-0.3, 0.3)))
    _pet_stats(F, cfg, tmp_path, None, monkeypatch)
    res = pst.run(('figure6',), pst.Inputs(), paper_dir=str(tmp_path))
    st = res['figure6']['draw']['stats']
    pet = F.pet_dataset_stats()
    pd.testing.assert_frame_equal(st['pet']['stats'], pet['stats'])
    base = tmp_path / 'Figure_6' / 'Figure_6'
    assert (tmp_path / 'Figure_6' / 'Figure_6_stats.xlsx').is_file()
    sheets = pd.ExcelFile(tmp_path / 'Figure_6' / 'Figure_6_stats.xlsx').sheet_names
    assert {'B_vs_macaque', 'C_agreement', 'C_pair_similarity'} <= set(sheets)
    pfg.run(res, None, paper_dir=str(tmp_path))
    assert (tmp_path / 'Figure_6' / 'Figure_6.png').is_file()
    assert 'Figure 6 |' in open(str(base) + '_legend.txt').read()


def test_brain_regions_keep_subcortical_and_drop_ventricles():
    keys = ['Prefrontal_cortex_(PFC)', 'Thalamus', 'Lateral_ventricle', 'unknown_x']
    out = MP.brain_regions(keys, 3)
    assert 'Prefrontal_cortex_(PFC)' in out and 'Thalamus' in out
    assert 'Lateral_ventricle' not in out and 'unknown_x' not in out


def test_fc_strength_kinds(tmp_path):
    m = pd.DataFrame([[1, .5, -.5], [.5, 1, .2], [-.5, .2, 1]], index=list('abc'), columns=list('abc'))
    f = tmp_path / 'm.csv'
    m.to_csv(f)
    z = np.arctanh
    pos = MP.node_strength(str(f), 'positive')
    assert np.isclose(pos['a'], (z(.5) + 0) / 2) and np.isclose(pos['c'], (0 + z(.2)) / 2)
    assert np.isclose(MP.node_strength(str(f), 'absolute')['a'], (z(.5) + z(.5)) / 2)
    assert np.isclose(MP.node_strength(str(f), 'signed')['a'], 0.0)          # +z and -z cancel
    with pytest.raises(ValueError):
        MP.node_strength(str(f), 'mean')


def test_sba_tables_from_the_selected_runs():
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    common = ['Frontal_lobe|q80', 'Frontal_lobe|q20', 'Temporal_lobe|q80', 'Temporal_lobe|q20']

    def prof(a, b, c, d):
        return pd.DataFrame({'region': common, 'median_r': [a, b, c, d]})
    cdict = dict(cfg=dict(VALUE_COL='median_r', ZSCORE=None), common=common,
                 species_of_label={'H1': 'Human', 'M1': 'Macaque', 'R1': 'Rat'},
                 profiles_by_label={'H1': {'r1': prof(.5, -.1, .4, -.2), 'r2': prof(.6, -.2, .3, -.1)},
                                    'M1': {'r1': prof(.4, -.1, .5, -.3)}, 'R1': {'r1': prof(.1, .1, .1, .1)}},
                 selected_units={'H1': {'r1': prof(.5, -.1, .4, -.2)}, 'M1': {'r1': prof(.4, -.1, .5, -.3)},
                                 'R1': {}})
    orig = fa._ordered_common_values, fa._normalize_values
    try:
        fa._ordered_common_values = lambda p, common, col: p[col].values.astype(float)
        fa._normalize_values = lambda v, norm: v
        t = MP.sba_tables(cdict, ['Human', 'Macaque'], 'selected', 'z')
        assert set(t) == {'Q80', 'Q20'}
        assert set(t['Q80'].species) == {'Human', 'Macaque'} and len(t['Q80']) == 4   # 2 selected runs x 2 regions
        assert t['Q20'].query("species == 'Human' and region == 'Frontal_lobe'").value.iloc[0] == -.1
        assert t['Q80'].subject.nunique() == 2
        assert len(MP.sba_tables(cdict, ['Human'], 'all', 'z')['Q80']) == 4     # both Human runs, 2 regions
    finally:
        fa._ordered_common_values, fa._normalize_values = orig


def test_lower_case_fingerprint_names_are_matched_to_the_atlas_spelling():
    """Regression (Clement's run, 2026-11-03): SBA names are lower case, PET/volume names are not: no region
    matched, every region was reported missing in every run."""
    def tab(names):
        return pd.DataFrame([dict(species='Human', bids_dir='d', subject=f's{i}', region=n, value=float(i))
                             for i in range(3) for n in names])
    t = {'metabolism': tab(['Prefrontal_cortex_(PFC)', 'Thalamus']), 'morphometry': tab(['Prefrontal_cortex_(PFC)']),
         'Q80': tab(['prefrontal_cortex_(pfc)', 'thalamus', 'only_sba'])}
    h = MP.harmonize_region_names(t)
    assert set(h['Q80'].region) == {'Prefrontal_cortex_(PFC)', 'Thalamus', 'only_sba'}
    kept, cov = MP.complete_regions(h['Q80'], ['Prefrontal_cortex_(PFC)', 'Thalamus'], ['Human'], 'Q80')
    assert kept == ['Prefrontal_cortex_(PFC)', 'Thalamus'] and (cov.n_missing == 0).all()


def test_sba_tails_in_r_units():
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    common = ['frontal|q80', 'frontal|q20', 'temporal|q80', 'temporal|q20']
    raw = [.4, -.1, .3, .05]                                   # temporal: even the lowest fifth is positive
    prof = pd.DataFrame({'base_name': common, 'q20q80_r': [1, -1, 1, -1], 'q20q80_r_raw': raw})
    c = dict(cfg=dict(VALUE_COL='q20q80_r', ZSCORE=None), common=common, species_of_label={'H': 'Human'},
             profiles_by_label={'H': {'r': prof}}, selected_units={'H': {'r': prof}})
    t = MP.sba_tables(c, ['Human'], 'selected')
    pos = t['pos'].set_index('region').value
    neg = t['neg'].set_index('region').value
    assert pos['frontal'] == .4 and pos['temporal'] == .3
    assert np.isclose(neg['frontal'], .1) and neg['temporal'] == 0.0           # q20 > 0 is not a negative index
    with pytest.raises(ValueError, match='_raw'):
        c['profiles_by_label']['H']['r'] = prof.drop(columns='q20q80_r_raw')
        c['selected_units']['H'] = c['profiles_by_label']['H']
        MP.sba_tables(c, ['Human'], 'selected')


def test_summary_stats_and_superior_inferior_position(tmp_path):
    rng = np.random.default_rng(1)
    regs = [f'R{i}' for i in range(12)]
    base = np.linspace(-1.5, 1.5, 12)

    def prof(noise, flip=1):
        return pd.DataFrame(flip * base + rng.normal(0, noise, (6, 12)), columns=regs)
    prof_c = {'Human': {'morphometry': prof(.2), 'metabolism': prof(.2)},
              'Macaque': {'morphometry': prof(.2), 'metabolism': prof(.2, flip=-1)}}
    si = pd.DataFrame([dict(species=sp, region=r, si=v) for sp in prof_c for r, v in zip(regs, np.linspace(-1, 1, 12))])
    out = MP.summary_stats(prof_c, regs, ['Human', 'Macaque'], si, n=300)
    c = out[out.analysis == 'conservation'].set_index('modality')
    assert c.loc['morphometry', 'r'] > .9 and c.loc['morphometry', 'p_holm'] < .05
    assert c.loc['metabolism', 'r'] < -.9                                  # mirrored profile
    s_ = out[out.analysis == 'superior-inferior'].set_index(['species', 'modality'])
    assert s_.loc[('Human', 'metabolism'), 'r'] > .9 and s_.loc[('Macaque', 'metabolism'), 'r'] < -.9
    assert (out.p_holm >= out.p).all()
    d = out[out.analysis == 'difference vs FDG difference'].set_index('modality')
    assert 'morphometry' in d.index and np.isfinite(d.loc['morphometry', 'r'])


def test_centroid_si_orientation(tmp_path):
    lab = np.zeros((4, 4, 8), np.int16)
    lab[:, :, :4] = 1                      # inferior half
    lab[:, :, 4:] = 2                      # superior half
    img = nib.Nifti1Image(np.stack([lab, lab], axis=3), np.eye(4))     # RAS+: +z = superior
    f = str(tmp_path / 'd.nii.gz')
    nib.save(img, f)
    ldf = pd.DataFrame({'label_id': [1, 2], 'region_name': ['L_Inf', 'L_Sup']})
    v = MP.centroid_si(f, ldf, 2)
    assert v['Sup'] > 0 > v['Inf'] and abs(v['Sup'] + v['Inf']) < 1e-9 and v['Sup'] <= 1
