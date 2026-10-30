"""
fig4_species_structure.py -- does the run selection remove noise WITHOUT flattening the
differences between species?

The selection raises every correlation (runs, datasets, profiles: they all move towards the
cross-species common profile), so a raw correlation cannot tell "noise removed" from "all
profiles pulled onto the same curve". The species signal is what is left once that common
profile is taken out. Per dataset:

    group profile  g_d  = mean of the (normalised) profiles of its runs, z-scored across regions
    common profile C    = median across species of the species medians of g_d (Fig. 4a, black)
    residual       R_d  = g_d - beta_d * C   (beta_d = least-squares coefficient: the part of g_d
                                               along the common profile is removed)

Three pools: all runs, selected runs, 200 random subsets of as many runs per dataset.

Questions (the same statistic is computed for the three pools, selected vs random is tested):
  S1  Do the residuals of two datasets of the SAME species agree more than those of different
      species?  gap = mean r(R_d, R_e) same species - mean r different species;
      P: species labels permuted across datasets.
  S2  Is the residual structure phylogenetic?  species-level residual distance (1 - r between
      the mean residual of each species) vs MRCA age of the species pair, Spearman Mantel
      (exact permutation, phylo.py).  Evidence that what remains is biology, ordered by the
      tree, not noise.
  Both are also given on the raw group profiles (without removing C) for comparison.

    python fig4_species_structure.py

Writes <PAPER_DIR>/Figure_4/Figure_4_species_structure.txt
"""
import os
import sys

import numpy as np

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

N_DRAWS = 200
N_PERM = 2000
SEED = 0


def _z(v):
    v = np.asarray(v, float)
    return (v - np.nanmean(v)) / np.nanstd(v)


def _corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() >= 3 else np.nan


def group_profiles(vecs, ds, idx_by_ds):
    return {d: _z(np.nanmean(vecs[idx], axis=0)) for d, idx in idx_by_ds.items()}


def common_profile(groups, sp_of_ds):
    sps = sorted(set(sp_of_ds[list(groups)]))
    spm = [np.nanmedian([groups[d] for d in groups if sp_of_ds[d] == s], axis=0) for s in sps]
    return _z(np.nanmedian(spm, axis=0))


def residuals(groups, C):
    out = {}
    for d, g in groups.items():
        ok = np.isfinite(g) & np.isfinite(C)
        beta = np.cov(g[ok], C[ok])[0, 1] / np.var(C[ok], ddof=1)
        out[d] = g - beta * C
    return out


def gap(profiles, sp_of_ds, labels=None):
    """mean r same species - mean r different species over dataset pairs."""
    ks = sorted(profiles)
    sp = sp_of_ds if labels is None else labels
    same, diff = [], []
    for i, d in enumerate(ks):
        for e in ks[i + 1:]:
            (same if sp[d] == sp[e] else diff).append(_corr(profiles[d], profiles[e]))
    return float(np.nanmean(same) - np.nanmean(diff)) if same and diff else np.nan


def gap_perm_p(profiles, sp_of_ds, rng, n_perm=N_PERM):
    obs = gap(profiles, sp_of_ds)
    sets = np.array(sorted(profiles))
    null = []
    for _ in range(n_perm):
        sp = sp_of_ds.copy()
        sp[sets] = sp_of_ds[rng.permutation(sets)]
        g = gap(profiles, sp_of_ds, labels=sp)
        if np.isfinite(g):
            null.append(g)
    null = np.array(null)
    return obs, (1 + int(np.sum(null >= obs - 1e-12))) / (len(null) + 1)


def species_mantel(profiles, sp_of_ds, with_p=True):
    """Spearman Mantel between the species-level distance (1 - r of the species mean
    profiles) and the MRCA age. Returns (rho, p_greater) or (rho, None)."""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import phylo
    sps = phylo.order_species(sorted(set(sp_of_ds[list(profiles)])))
    prof = {s: np.nanmean([profiles[d] for d in profiles if sp_of_ds[d] == s], axis=0)
            for s in sps}
    n = len(sps)
    Df = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            Df[i, j] = Df[j, i] = 1.0 - _corr(prof[sps[i]], prof[sps[j]])
    Dp, _ = phylo.phylo_matrices(sps)
    if with_p:
        m = phylo.mantel_test(Df, Dp, method="spearman")
        return m["r"], m["p_greater"]
    iu = np.triu_indices(n, 1)
    return phylo._corr(Df[iu], Dp[iu], "spearman"), None


def analyse(vecs, ds, selected, sp_of_ds, n_draws=N_DRAWS, n_perm=N_PERM, seed=SEED):
    rng = np.random.default_rng(seed)
    sets = np.unique(ds)
    idx_all = {d: np.flatnonzero(ds == d) for d in sets}
    idx_sel = {d: np.flatnonzero((ds == d) & selected) for d in sets}
    n_sel = {d: len(i) for d, i in idx_sel.items()}

    def stats_for(idx_by_ds, with_p):
        g = group_profiles(vecs, ds, idx_by_ds)
        res = residuals(g, common_profile(g, sp_of_ds))
        out = {}
        for name, prof in (("raw", g), ("residual", res)):
            gp = gap_perm_p(prof, sp_of_ds, rng, n_perm) if with_p else (gap(prof, sp_of_ds), None)
            rho, p_m = species_mantel(prof, sp_of_ds, with_p=with_p)
            out[name] = dict(gap=gp[0], p_gap=gp[1], rho=rho, p_rho=p_m)
        return out

    a, s = stats_for(idx_all, True), stats_for(idx_sel, True)
    draws = [stats_for({d: rng.choice(idx_all[d], n_sel[d], replace=False) for d in sets},
                       False) for _ in range(n_draws)]
    lines = ["Fig. 4: species structure after the run selection (noise removed, species "
             "differences kept?)",
             "  group profile per dataset = mean of its runs; residual = group profile minus "
             "its projection on the common profile (Fig. 4a)",
             f"  random = as many runs per dataset as selected, {n_draws} draws; "
             "P(selected vs random) = share of draws at least as large as the selection"]
    for name, title in (("raw", "raw group profiles"),
                        ("residual", "residuals (common profile removed)")):
        gd = np.array([x[name]["gap"] for x in draws])
        rd = np.array([x[name]["rho"] for x in draws])
        lines.append(f"  {title}")
        lines.append("    S1 same-species minus different-species r between datasets: "
                     f"all runs {a[name]['gap']:+.3f} (P = {a[name]['p_gap']:.4g}); random "
                     f"median {np.nanmedian(gd):+.3f} (95th pct {np.nanpercentile(gd, 95):+.3f}); "
                     f"selected {s[name]['gap']:+.3f} (permutation P = {s[name]['p_gap']:.4g}); "
                     f"P(selected vs random) = "
                     f"{(1 + np.sum(gd >= s[name]['gap'] - 1e-12)) / (n_draws + 1):.4g}")
        lines.append("    S2 Mantel (Spearman) species distance vs MRCA age: "
                     f"all runs rho = {a[name]['rho']:+.3f} (P = {a[name]['p_rho']:.4g}); random "
                     f"median {np.nanmedian(rd):+.3f} (95th pct {np.nanpercentile(rd, 95):+.3f}); "
                     f"selected rho = {s[name]['rho']:+.3f} (P = {s[name]['p_rho']:.4g}); "
                     f"P(selected vs random) = "
                     f"{(1 + np.sum(rd >= s[name]['rho'] - 1e-12)) / (n_draws + 1):.4g}")
    lines.append("  reading: S1/S2 on the residuals >= random and > 0 means the species signal "
                 "survives (and is cleaner); a collapse towards 0 would mean the selection "
                 "flattens species differences")
    return lines


def main():
    from Exemples.Study_EDNiX.analysis import paper_data as pd_
    from Exemples.Study_EDNiX.analysis.fig4b_selection_improvement import load_runs
    vecs, ds, selected, sp_of_ds, labels = load_runs()
    lines = analyse(vecs, ds, selected, sp_of_ds)
    out = os.path.join(pd_.PAPER_DIR, "Figure_4")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "Figure_4_species_structure.txt")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"[saved] {path}")


if __name__ == "__main__":
    main()
