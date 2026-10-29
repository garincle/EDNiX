"""
fig4_selection_benefit.py -- does keeping the 10 selected runs make a dataset more comparable
with the other species than not selecting?

Unit = dataset (paired design: each dataset is measured under every condition).
Group profile of a dataset = mean of the (normalised) profiles of the runs of a given pool,
z-scored across regions. Three pools per dataset:
    all       every run
    selected  the runs kept by the final round
    random    as many runs as selected, drawn at random from the dataset (2,000 draws; one
              draw = one random subset for every dataset)

Two outcomes:
  A  r with the cross-species consensus (median across the OTHER species of the species
     profiles; leave-one-species-out, as in Fig. 4a). Partly circular: it is the quantity the
     selection maximises, so A shows that the criterion works, not that it generalises.
  B  r with the other dataset of the same species (datasets with a sibling only). Independent of
     the criterion for the sibling's runs: asks whether selected profiles of two independent
     datasets of one species agree better than unselected ones (Fig. 4c, Supp. 10).

Tests (one-sided, selected > comparison):
  paired Wilcoxon signed-rank and paired t-test on Fisher z of r, selected vs all, across datasets;
  selected vs random: share of the 2,000 random draws whose mean r over datasets is >= selected
  (and, per dataset, the share of draws with r >= selected).

    python fig4_selection_benefit.py

Writes <PAPER_DIR>/Figure_4/Figure_4_selection_benefit.txt|xlsx
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

N_DRAWS = 2000
SEED = 0


def _z(v):
    v = np.asarray(v, float)
    return (v - np.nanmean(v)) / np.nanstd(v)


def _group(vecs, idx):
    return _z(np.nanmean(vecs[idx], axis=0))


def _corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() >= 3 else np.nan


def outcomes(groups, sp_of_ds):
    """groups: {dataset: group profile}. Returns (A, B) dicts dataset -> r.
    A: r with the median over other species of the species-median profiles.
    B: r with the mean profile of the other datasets of the same species (absent if none)."""
    species = {s: np.nanmedian([groups[d] for d in groups if sp_of_ds[d] == s], axis=0)
               for s in set(sp_of_ds[list(groups)])}
    A, B = {}, {}
    for d, g in groups.items():
        others = [species[s] for s in species if s != sp_of_ds[d]]
        A[d] = _corr(g, _z(np.nanmedian(others, axis=0))) if others else np.nan
        sib = [groups[e] for e in groups if e != d and sp_of_ds[e] == sp_of_ds[d]]
        if sib:
            B[d] = _corr(g, _z(np.nanmean(sib, axis=0)))
    return A, B


def _fz(x):
    return np.arctanh(np.clip(np.asarray(x, float), -0.999, 0.999))


def _paired(sel, base, name):
    ds = [d for d in sel if d in base and np.isfinite(sel[d]) and np.isfinite(base[d])]
    a, b = _fz([sel[d] for d in ds]), _fz([base[d] for d in ds])
    diff = a - b
    w = stats.wilcoxon(diff, alternative="greater") if np.any(diff != 0) else None
    t = stats.ttest_1samp(diff, 0.0, alternative="greater")
    return (f"{name}: n = {len(ds)} datasets; median r {np.median([sel[d] for d in ds]):.3f} "
            f"vs {np.median([base[d] for d in ds]):.3f}; datasets improved "
            f"{int((diff > 0).sum())}/{len(ds)}; Wilcoxon signed-rank one-sided P = "
            f"{(w.pvalue if w is not None else np.nan):.4g}; paired t (Fisher z) "
            f"t({len(ds) - 1}) = {t.statistic:.2f}, one-sided P = {t.pvalue:.4g}")


def analyse(vecs, ds, selected, sp_of_ds, n_draws=N_DRAWS, seed=SEED):
    rng = np.random.default_rng(seed)
    sets = np.unique(ds)
    idx_all = {d: np.flatnonzero(ds == d) for d in sets}
    idx_sel = {d: np.flatnonzero((ds == d) & selected) for d in sets}
    g_all = {d: _group(vecs, idx_all[d]) for d in sets}
    g_sel = {d: _group(vecs, idx_sel[d]) for d in sets if len(idx_sel[d]) >= 2}
    A_all, B_all = outcomes(g_all, sp_of_ds)
    A_sel, B_sel = outcomes(g_sel, sp_of_ds)
    drawsA, drawsB = [], []
    for _ in range(n_draws):
        g = {d: _group(vecs, rng.choice(idx_all[d], len(idx_sel[d]), replace=False))
             for d in g_sel}
        a, b = outcomes(g, sp_of_ds)
        drawsA.append(a)
        drawsB.append(b)
    tab = []
    for d in sorted(g_sel):
        ra = np.array([x.get(d, np.nan) for x in drawsA])
        rb = np.array([x.get(d, np.nan) for x in drawsB])
        tab.append(dict(dataset=d, species=sp_of_ds[d], n_runs=len(idx_all[d]),
                        n_selected=len(idx_sel[d]),
                        A_all=A_all.get(d), A_selected=A_sel.get(d),
                        A_random_median=float(np.nanmedian(ra)),
                        A_p_vs_random=(1 + np.sum(ra >= A_sel[d] - 1e-12)) / (n_draws + 1),
                        B_all=B_all.get(d), B_selected=B_sel.get(d),
                        B_random_median=float(np.nanmedian(rb)) if np.isfinite(rb).any() else np.nan,
                        B_p_vs_random=((1 + np.sum(rb >= B_sel[d] - 1e-12)) / (n_draws + 1)
                                       if d in B_sel else np.nan)))
    tab = pd.DataFrame(tab)
    lines = ["Fig. 4: does selection make a dataset more comparable with other species?",
             f"  pools: all runs, selected runs, {n_draws} random subsets of the same size "
             "(group profile = mean of the run profiles, z-scored across regions)",
             "  note: 'selected vs all' favours 'all' (a mean of many runs is less noisy than a "
             "mean of 10); the fair comparison is 'selected vs random' (same number of runs)"]
    lines.append("  A  r with the leave-one-species-out consensus (partly circular)")
    lines.append("    " + _paired(A_sel, A_all, "selected vs all"))
    mA = np.array([np.nanmean(list(a.values())) for a in drawsA])
    obsA = float(np.nanmean(list(A_sel.values())))
    lines.append(f"    selected vs random: mean r {obsA:.3f} vs random median "
                 f"{np.median(mA):.3f} (95th percentile {np.percentile(mA, 95):.3f}); "
                 f"one-sided P = {(1 + np.sum(mA >= obsA - 1e-12)) / (n_draws + 1):.4g}")
    lines.append("  B  r with the other dataset of the same species (independent of the criterion)")
    lines.append("    " + _paired(B_sel, B_all, "selected vs all"))
    mB = np.array([np.nanmean(list(b.values())) for b in drawsB])
    obsB = float(np.nanmean(list(B_sel.values())))
    lines.append(f"    selected vs random: mean r {obsB:.3f} vs random median "
                 f"{np.median(mB):.3f} (95th percentile {np.percentile(mB, 95):.3f}); "
                 f"one-sided P = {(1 + np.sum(mB >= obsB - 1e-12)) / (n_draws + 1):.4g}")
    return lines, tab


def main():
    from Exemples.Study_EDNiX.analysis import paper_data as pd_
    from Exemples.Study_EDNiX.analysis.fig4b_selection_improvement import load_runs
    vecs, ds, selected, sp_of_ds, labels = load_runs()
    lines, tab = analyse(vecs, ds, selected, sp_of_ds)
    tab.insert(0, "label", [labels[i] for i in tab["dataset"]])
    out = os.path.join(pd_.PAPER_DIR, "Figure_4")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "Figure_4_selection_benefit.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    tab.to_excel(os.path.join(out, "Figure_4_selection_benefit.xlsx"), index=False)
    print("\n".join(lines))
    print(tab.round(3).to_string(index=False))
    print(f"[saved] {out}/Figure_4_selection_benefit.txt|xlsx")


if __name__ == "__main__":
    main()
