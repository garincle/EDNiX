import numpy as np
import pandas as pd
import pytest

from Statistics.Group_fMRI.EDNiX_longitudinal import fit_gamm, _ml_random_intercept
from Statistics.Group_fMRI._Group_anal_3dLMEr_SBA import _cluster_mask
from Statistics.Group_fMRI.Select_cpmnt_number import match_components
import Statistics.Evo.PGLS as P


def _long_data(seed, effect):
    rng = np.random.default_rng(seed)
    rows = []
    for s in range(12):
        b = rng.normal(0, 1)
        for a in np.sort(rng.uniform(0, 10, 4)):
            rows.append(dict(subject=s, age=a, y=b + effect * np.sin(a / 3) + rng.normal(0, 1)))
    return pd.DataFrame(rows)


def test_random_intercept_ml_matches_statsmodels():
    sm = pytest.importorskip('statsmodels.api')
    rng = np.random.default_rng(3)
    rows = [dict(subject=s, age=a, y=b + 0.3 * a + rng.normal(0, 1))
            for s in range(30) for b in [rng.normal(0, 2)] for a in rng.uniform(0, 10, 5)]
    d = pd.DataFrame(rows)
    X = np.column_stack([np.ones(len(d)), d.age])
    m = sm.MixedLM(d.y, X, groups=d.subject).fit(reml=False)
    llf, beta, cov, _ = _ml_random_intercept(d.y.values, X, d.subject.values)
    assert abs(llf - m.llf) < 1e-3 and np.allclose(beta, m.fe_params.values, atol=1e-4)


@pytest.mark.slow
def test_age_effect_type_one_error_is_controlled():
    p = np.array([fit_gamm(_long_data(s, 0.0), 'y', 'age')['p_age'] for s in range(300)])
    assert (p < 0.05).mean() < 0.10          # was ~0.19 with the minimum coefficient p-value


def test_cluster_extent_is_applied():
    z = np.zeros((20, 20, 20))
    z[2:5, 2:5, 2:5] = 5
    z[10, 10, 10] = -5
    assert _cluster_mask(z, 0.001, 5).sum() == 27


def test_component_matching_ignores_sign():
    from Statistics.Group_fMRI.Select_parameters_DicL import match_components as match_dicl
    rng = np.random.default_rng(0)
    A = rng.normal(size=(4, 500))
    B = A[[2, 0, 3, 1]] * np.array([1, -1, 1, -1])[:, None] + 0.1 * rng.normal(size=(4, 500))
    assert match_components(A, B)['mean_correlation'] > 0.99
    assert match_dicl(A, B)['mean'] > 0.99


def test_pgls_covariance_from_newick(tmp_path):
    nwk = tmp_path / 't.nwk'
    nwk.write_text('((A:1,B:1)n1:2,(C:0.5,D:0.5):2.5)root;')
    v, _ = P.load_newick_vcv(str(nwk), {k: k for k in 'ABCD'})
    assert v.loc['A', 'A'] == 3 and v.loc['A', 'B'] == 2 and v.loc['A', 'C'] == 0 and v.loc['C', 'D'] == 2.5
