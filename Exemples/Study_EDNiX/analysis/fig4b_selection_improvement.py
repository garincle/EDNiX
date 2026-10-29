"""
fig4b_selection_improvement.py -- is the gain of the run selection in Fig. 4b significant?

Fig. 4b shows two separate tests (all runs: no species effect; selected runs: species effect).
Two tests, one non-significant and one significant, do not show that the two situations differ
(Gelman & Stern 2006), and the pair-level Wilcoxon ignores that pairs share runs and datasets.
This script tests the improvement itself, with the dataset as the unit of dependence.

Statistic (same as Fig. 4b): T = median r(same species, other dataset) - median r(different
species), r = Pearson correlation between the regional profiles of two runs.

  1. SELECTION vs RANDOM CHOICE (run-level). The selected runs are compared with 2,000 random
     draws of as many runs per dataset, taken from all the runs of that dataset. P = share of
     draws with T >= T(selected). Answers: does the criterion do better than keeping the same
     number of runs at random?
  2. IMPROVEMENT vs SPECIES-LABEL PERMUTATION (dataset-level). Delta = T(selected) - T(all runs).
     Species labels are permuted across datasets (the same permutation for both pools), as for
     the P in Fig. 4b; P = share of permutations with Delta >= observed.
  3. CLUSTER BOOTSTRAP over datasets (95% interval of Delta, percentile; datasets resampled
     with replacement, all their runs kept).

    python fig4b_selection_improvement.py

Reads the saved run selection like run_paper.py (needs its 07_state); writes
<PAPER_DIR>/Figure_4/Figure_4_selection_improvement.txt
"""
import os
import sys

import numpy as np

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

N_PERM = 2000
SEED = 0


def _corr_matrix(vecs):
    """Pearson r between rows, pairwise complete (>= 3 shared finite values)."""
    n = len(vecs)
    r = np.full((n, n), np.nan)
    fin = np.isfinite(vecs)
    for i in range(n):
        for j in range(i + 1, n):
            both = fin[i] & fin[j]
            if both.sum() >= 3:
                r[i, j] = r[j, i] = np.corrcoef(vecs[i][both], vecs[j][both])[0, 1]
    return r


def _t_stat(r, ds, sp_of_ds, idx):
    """median r(same species, other dataset) - median r(different species), runs `idx`."""
    idx = np.asarray(idx)
    d = ds[idx]
    sub = r[np.ix_(idx, idx)]
    iu = np.triu_indices(len(idx), k=1)
    rr, da, db = sub[iu], d[iu[0]], d[iu[1]]
    ok = np.isfinite(rr)
    rr, da, db = rr[ok], da[ok], db[ok]
    same_sp = sp_of_ds[da] == sp_of_ds[db]
    within = rr[same_sp & (da != db)]
    across = rr[~same_sp]
    if len(within) < 1 or len(across) < 1:
        return np.nan
    return float(np.median(within) - np.median(across))


def analyse(r, ds, selected, sp_of_ds, n_perm=N_PERM, seed=SEED):
    """r: run x run correlation matrix; ds: dataset index per run; selected: boolean per run;
    sp_of_ds: species per dataset index. Returns a dict of results."""
    rng = np.random.default_rng(seed)
    all_idx, sel_idx = np.arange(len(ds)), np.flatnonzero(selected)
    t_all, t_sel = _t_stat(r, ds, sp_of_ds, all_idx), _t_stat(r, ds, sp_of_ds, sel_idx)
    delta = t_sel - t_all
    # 1. random choice of the same number of runs per dataset
    per_ds = {d: np.flatnonzero(ds == d) for d in np.unique(ds)}
    n_sel = {d: int(selected[i].sum()) for d, i in per_ds.items()}
    rand = np.array([_t_stat(r, ds, sp_of_ds, np.concatenate(
        [rng.choice(per_ds[d], n_sel[d], replace=False) for d in per_ds if n_sel[d]]))
        for _ in range(n_perm)])
    rand = rand[np.isfinite(rand)]
    p_random = (1 + int(np.sum(rand >= t_sel - 1e-12))) / (len(rand) + 1)
    # 2. species labels permuted across datasets, same permutation for both pools
    null = []
    for _ in range(n_perm):
        spp = rng.permutation(sp_of_ds)
        a, s = _t_stat(r, ds, spp, all_idx), _t_stat(r, ds, spp, sel_idx)
        if np.isfinite(a) and np.isfinite(s):
            null.append(s - a)
    null = np.array(null)
    p_perm = (1 + int(np.sum(null >= delta - 1e-12))) / (len(null) + 1)
    # 3. cluster bootstrap over datasets
    boots = []
    dsets = np.unique(ds)
    for _ in range(n_perm):
        pick = rng.choice(dsets, len(dsets), replace=True)
        idx_all = np.concatenate([per_ds[d] for d in pick])
        # a resampled dataset drawn twice is a second copy: give copies new dataset ids
        new_ds = np.concatenate([np.full(len(per_ds[d]), k) for k, d in enumerate(pick)])
        sp_b = np.array([sp_of_ds[d] for d in pick])
        sel_b = selected[idx_all]
        rb = r[np.ix_(idx_all, idx_all)]
        a = _t_stat(rb, new_ds, sp_b, np.arange(len(idx_all)))
        s = _t_stat(rb, new_ds, sp_b, np.flatnonzero(sel_b))
        if np.isfinite(a) and np.isfinite(s):
            boots.append(s - a)
    boots = np.array(boots)
    lo, hi = (np.percentile(boots, [2.5, 97.5]) if len(boots) else (np.nan, np.nan))
    return dict(t_all=t_all, t_sel=t_sel, delta=delta, p_random=p_random,
                random_median=float(np.median(rand)), random_q95=float(np.percentile(rand, 95)),
                p_perm=p_perm, boot_lo=float(lo), boot_hi=float(hi), n_boot=len(boots),
                n_runs=len(ds), n_selected=int(selected.sum()), n_datasets=len(per_ds),
                n_perm=n_perm)


def report(res):
    return [
        "Fig. 4b: is the gain from selection significant?  T = median r(same species, other "
        "dataset) - median r(different species)",
        f"  T, all runs ({res['n_runs']} runs, {res['n_datasets']} datasets) = {res['t_all']:+.3f}",
        f"  T, selected runs ({res['n_selected']}) = {res['t_sel']:+.3f}; "
        f"Delta = T(selected) - T(all) = {res['delta']:+.3f}",
        f"  1. selection vs random choice of as many runs per dataset ({res['n_perm']} draws): "
        f"random median T = {res['random_median']:+.3f}, 95th percentile "
        f"{res['random_q95']:+.3f}; one-sided P = {res['p_random']:.4g}",
        f"  2. Delta vs species labels permuted across datasets ({res['n_perm']} permutations): "
        f"one-sided P = {res['p_perm']:.4g}",
        f"  3. cluster bootstrap over datasets ({res['n_boot']} resamples): Delta 95% interval "
        f"[{res['boot_lo']:+.3f}, {res['boot_hi']:+.3f}]",
    ]


def main():
    from Exemples.Study_EDNiX.analysis import paper_stats as pst
    from Exemples.Study_EDNiX.analysis import paper_data as pd_
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    c = pst.Inputs().fmri()
    cfg = c["cfg"]
    value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
    common, sol = c["common"], c["species_of_label"]
    labels = sorted(c["profiles_by_label"])
    vecs, ds, keys, selected = [], [], [], []
    for i, lab in enumerate(labels):
        kept = set(c["selected_units"].get(lab, {}))
        for key, prof in c["profiles_by_label"][lab].items():
            v = fa._ordered_common_values(prof, common, value_col)
            if not np.isfinite(v).any():
                continue
            vecs.append(fa._normalize_values(v, norm))
            ds.append(i)
            keys.append((lab, key))
            selected.append(key in kept)
    vecs, ds, selected = np.vstack(vecs), np.array(ds), np.array(selected)
    if selected.sum() == 0:
        raise SystemExit("no selected run matches a run key of profiles_by_label: check the keys")
    sp_of_ds = np.array([sol.get(l, l) for l in labels])
    res = analyse(_corr_matrix(vecs), ds, selected, sp_of_ds)
    lines = report(res)
    out = os.path.join(pd_.PAPER_DIR, "Figure_4")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "Figure_4_selection_improvement.txt")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"[saved] {path}")


if __name__ == "__main__":
    main()
