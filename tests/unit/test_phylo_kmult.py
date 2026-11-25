"""phylo_kmult: K_mult (Adams 2014) calibration and the Fig. 4e pools on synthetic profiles."""
import numpy as np
import pytest

from Paper.analysis.fmri_fingerprints import phylo, phylo_kmult as pk

SPS = phylo.order_species(['Mouse', 'Rat', 'Mouselemur', 'Marmoset', 'Macaque', 'Dog', 'Human'])


def _C():
    return phylo.phylo_matrices(SPS)[1]


def test_k_is_one_under_brownian_motion_and_lower_without_structure():
    C = _C()
    T = pk.KTest(C)
    rng = np.random.default_rng(0)
    L = np.linalg.cholesky(C)
    kb = [T.k(L @ rng.standard_normal((7, 30))) for _ in range(300)]
    kr = [T.k(rng.standard_normal((7, 30))) for _ in range(300)]
    assert np.median(kb) == pytest.approx(1.0, abs=0.08)
    assert np.median(kr) < 0.75


def test_k_matches_relabelled_tree_and_stack():
    C = _C()
    T = pk.KTest(C)
    Y = np.random.default_rng(1).standard_normal((7, 30))
    p = np.array([3, 0, 6, 1, 5, 2, 4])
    T2 = pk.KTest(C[np.ix_(p, p)])
    assert T2.k(Y[p]) == pytest.approx(T.k(Y))          # same species, same tree, reordered
    stack = T.null(Y)
    assert stack.shape == (5040,)
    assert stack[0] == pytest.approx(T.k(Y))            # identity permutation first
    k, pval = T.test(Y)
    assert 1 / 5040 <= pval <= 1


def test_delta_zero_for_identical_pools():
    T = pk.KTest(_C())
    Y = np.random.default_rng(2).standard_normal((7, 30))
    d = T.delta(Y, Y)
    assert d["delta"] == pytest.approx(0.0) and d["p"] == pytest.approx(1.0)


def test_loso_residuals_remove_the_shared_profile():
    rng = np.random.default_rng(3)
    common = rng.standard_normal(30)
    sp_of_ds = np.array(SPS)
    groups = {d: pk._z(5 * common + 0.3 * rng.standard_normal(30)) for d in range(7)}
    res = pk.loso_residuals(groups, sp_of_ds, SPS)
    for d in groups:
        assert abs(np.corrcoef(groups[d], common)[0, 1]) > 0.95
        assert abs(np.corrcoef(res[d], common)[0, 1]) < 0.3


def test_compute_detects_brownian_species_signal():
    rng = np.random.default_rng(4)
    C = _C()
    sig = np.linalg.cholesky(C / C.max()) @ rng.standard_normal((7, 30))
    common = rng.standard_normal(30)
    vecs, ds, sel = [], [], []
    for d, s in enumerate(SPS):
        for r in range(20):
            vecs.append(common + sig[d] + 0.5 * rng.standard_normal(30))
            ds.append(d)
            sel.append(r < 10)
    out = pk.compute(np.array(vecs), np.array(ds), np.array(sel), np.array(SPS), n_draws=10)
    for k in pk.KINDS:
        assert out[k]["all"][1] < 0.05 and out[k]["selected"][1] < 0.05
        assert len(out[k]["random"]) == 10
    lo = out["loso_species"]
    assert len(lo) == 7 * 2 * 2 and (lo.n_perm == 720).all()
    assert any("leave-one-species-out" in ln for ln in pk.result_lines(out))
    assert set(pk.table(out).pool) == {"all runs", "selected", "random"}
