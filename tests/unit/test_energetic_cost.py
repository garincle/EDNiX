import numpy as np
import pandas as pd
import nibabel as nib

from Exemples.Study_EDNiX.analysis.library import energetic_cost as EC


def _smooth(shape, seed, sigma=3):
    from scipy import ndimage
    x = ndimage.gaussian_filter(np.random.default_rng(seed).standard_normal(shape), sigma)
    return (x - x.mean()) / x.std()


def test_degree_centrality_matches_brute_force():
    rng = np.random.default_rng(0)
    base = rng.standard_normal((2, 200))
    ts = np.vstack([base[0] + .8 * rng.standard_normal((30, 200)), base[1] + .8 * rng.standard_normal((20, 200)),
                    rng.standard_normal((10, 200))])
    d = EC.degree_centrality(ts, chunk=7)
    r = np.corrcoef(ts)
    np.fill_diagonal(r, 0)
    thr = EC.r_threshold(200)
    ref = np.where(r > thr, r, 0).sum(1)
    assert np.allclose(d, ref, atol=1e-4)
    assert d[:30].mean() > d[50:].mean()                      # voxels of a large network have higher degree


def test_degree_map_and_group(tmp_path):
    rng = np.random.default_rng(1)
    d = rng.standard_normal((5, 5, 4, 60)).astype(np.float32)
    d[..., :] += rng.standard_normal(60)[None, None, None, :] * 2
    d[0, 0, 0] = 0                                             # constant voxel -> outside
    p = tmp_path / 'b.nii.gz'
    nib.save(nib.Nifti1Image(d, np.eye(4)), str(p))
    m, aff = EC.degree_map(str(p))
    assert np.isnan(m[0, 0, 0]) and abs(np.nanmean(m)) < 1e-5
    g = EC.group_degree_map([m, m, m])
    assert np.allclose(g[np.isfinite(g)], m[np.isfinite(m)])


def test_resample_shift_and_smooth():
    a = np.zeros((10, 10, 10), np.float32)
    a[3:7, 3:7, 3:7] = 1
    aff = np.eye(4)
    same = EC.resample_to_grid(a, aff, a.shape, aff)
    assert np.allclose(same[3:7, 3:7, 3:7], 1) and np.isnan(same[0, 0, 0]) is np.False_
    shifted = np.eye(4)
    shifted[0, 3] = 1
    out = EC.resample_to_grid(a, aff, a.shape, shifted)       # ref voxel i is world i+1 -> source i+1
    assert out[2, 4, 4] == 1 and out[6, 4, 4] == 0
    a2 = np.where(a > 0, 1.0, np.nan).astype(np.float32)
    s = EC.smooth_masked(a2, 2.0, aff)
    assert np.allclose(s[np.isfinite(s)], 1.0, atol=1e-5) and np.isnan(s[0, 0, 0])


def test_coupling_analysis_detects_alignment_and_not_chance():
    shape = (24, 24, 24)
    mask = np.ones(shape, bool)
    deg = _smooth(shape, 1)
    pet = .7 * deg + .7 * _smooth(shape, 2)
    st, res = EC.coupling_analysis(deg, pet, mask, n_surr=200, seed=0)
    assert st['r'] > .5 and st['p_surrogate'] < .05
    assert abs(np.nanmean(res)) < 1e-4
    other = _smooth(shape, 3)
    st2, _ = EC.coupling_analysis(deg, other, mask, n_surr=200, seed=0)
    assert st2['p_surrogate'] > .05                            # smooth unrelated maps: spatial null keeps the P honest


def test_expansion_table_and_analysis():
    regs = [f'R{i}' for i in range(12)]
    exp = np.linspace(-1, 1, 12)
    rows_v, rows_f = [], []
    for sp, sign in (('Human', 1), ('Macaque', 0)):
        for s in range(3):
            for r, e in zip(regs, exp):
                rows_v.append(dict(species=sp, subject=f's{s}', region=r, value=2 ** (sign * e)))
                rows_f.append(dict(species=sp, subject=f's{s}', region=r, value=2 ** (sign * .5 * e)))
    tab = EC.expansion_table(pd.DataFrame(rows_v), pd.DataFrame(rows_f), regions=regs)
    assert np.allclose(tab['log2_expansion'], exp) and np.allclose(tab['log2_fdg_ratio'], .5 * exp)
    out = EC.expansion_analysis(tab, cortical=regs[:8], n=300)
    assert out.loc[0, 'r'] > .99 and set(out['regions']) == {'all regions', 'cortical regions'}
