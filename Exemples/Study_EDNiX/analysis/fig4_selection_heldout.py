"""
fig4_selection_heldout.py -- non-circular check of the run selection: select on half of the
regions, evaluate on the other half.

Any criterion that keeps runs resembling a target makes the kept runs resemble each other on
the values it scored. Here the regional values are split at random into two halves (50 splits):
  selection half   leave-one-species-out consensus (median over the other species of the
                   species medians of the dataset means, all runs) and, per dataset, the 10 runs
                   with the highest Pearson r with it (simplified version of the criterion:
                   one round, no invariance weighting)
  held-out half    the kept runs are evaluated on values that played no part in their choice:
                   E1  median r between runs of different datasets (same and different species)
                   E2  median over datasets of r(group profile, LOSO consensus), held-out values
                   G   median r same dataset - median r across datasets
Compared with 10 runs drawn at random per dataset (200 draws per split). P = share of random
draws whose split-averaged statistic is at least as good as the selection.

    python fig4_selection_heldout.py

Writes <PAPER_DIR>/Figure_4/Figure_4_selection_heldout.txt
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

N_SPLITS = 50
N_DRAWS = 200
N_KEEP = 10
SEED = 0


def _z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def _corr_rows(m):
    return pd.DataFrame(m.T).corr(min_periods=3).to_numpy()


def _corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() >= 3 else np.nan


def _loso(vecs, ds, sp_of_ds, d):
    """Median over the other species of the species medians of the dataset means (all runs)."""
    means = {e: np.nanmean(vecs[ds == e], axis=0) for e in np.unique(ds)}
    sps = [s for s in set(sp_of_ds) if s != sp_of_ds[d]]
    return np.nanmedian([np.nanmedian([means[e] for e in means if sp_of_ds[e] == s], axis=0)
                         for s in sps], axis=0)


def _evaluate(r, vecs, ds, sp_of_ds, idx):
    d = ds[idx]
    iu = np.triu_indices(len(idx), k=1)
    rr = r[np.ix_(idx, idx)][iu]
    da, db = d[iu[0]], d[iu[1]]
    ok = np.isfinite(rr)
    rr, da, db = rr[ok], da[ok], db[ok]
    cross = da != db
    e1 = float(np.median(rr[cross]))
    g = float(np.median(rr[~cross]) - e1)
    e2 = float(np.median([_corr(_z(np.nanmean(vecs[idx[d == e]], axis=0)),
                                _z(_loso(vecs, ds, sp_of_ds, e))) for e in np.unique(d)]))
    return e1, e2, g


def analyse(vecs, ds, sp_of_ds, n_splits=N_SPLITS, n_draws=N_DRAWS, n_keep=N_KEEP, seed=SEED):
    rng = np.random.default_rng(seed)
    nreg = vecs.shape[1]
    per_ds = {d: np.flatnonzero(ds == d) for d in np.unique(ds)}
    sel_stats, rand_stats = [], []
    for _ in range(n_splits):
        perm = rng.permutation(nreg)
        a, b = perm[: nreg // 2], perm[nreg // 2:]
        va, vb = vecs[:, a], vecs[:, b]
        keep = []
        for d, idx in per_ds.items():
            target = _loso(va, ds, sp_of_ds, d)
            sc = np.array([_corr(va[i], target) for i in idx])
            sc[~np.isfinite(sc)] = -np.inf
            keep.append(idx[np.argsort(-sc)[: min(n_keep, len(idx))]])
        keep = np.concatenate(keep)
        rb = _corr_rows(vb)
        sel_stats.append(_evaluate(rb, vb, ds, sp_of_ds, keep))
        rand_stats.append([_evaluate(rb, vb, ds, sp_of_ds, np.concatenate(
            [rng.choice(idx, min(n_keep, len(idx)), replace=False) for idx in per_ds.values()]))
            for _ in range(n_draws)])
    sel = np.mean(sel_stats, axis=0)                       # (E1, E2, G) averaged over splits
    rnd = np.mean(np.array(rand_stats), axis=0)            # draws x 3, averaged over splits
    p = [(1 + np.sum(rnd[:, 0] >= sel[0])) / (n_draws + 1),
         (1 + np.sum(rnd[:, 1] >= sel[1])) / (n_draws + 1),
         (1 + np.sum(rnd[:, 2] <= sel[2])) / (n_draws + 1)]
    names = ("E1 median r between runs of different datasets",
             "E2 median r(dataset group profile, LOSO consensus)",
             "G  same dataset - cross-dataset (lower = less dataset-specific)")
    lines = [f"Fig. 4 held-out check: select on half of the {nreg} regional values, evaluate on "
             f"the other half ({n_splits} random splits; {n_keep} runs per dataset; "
             f"{n_draws} random draws per split)"]
    for k in range(3):
        lines.append(f"  {names[k]}: selected {sel[k]:+.3f}, random median "
                     f"{np.median(rnd[:, k]):+.3f}; one-sided P = {p[k]:.4g}")
    lines.append("  criterion simplified (one round, no invariance weighting): tests the principle, "
                 "not the exact selection of the paper")
    return lines


def main():
    from Exemples.Study_EDNiX.analysis import paper_data as pd_
    from Exemples.Study_EDNiX.analysis.fig4b_selection_improvement import load_runs
    vecs, ds, selected, sp_of_ds, labels = load_runs()
    lines = analyse(vecs, ds, sp_of_ds)
    out = os.path.join(pd_.PAPER_DIR, "Figure_4")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "Figure_4_selection_heldout.txt")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"[saved] {path}")


if __name__ == "__main__":
    main()
