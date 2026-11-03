"""
paper_extra_stats.py -- two analyses that the manuscript needs and that the Figure stats files did not carry.

  1. Fig. 2: isocortical surface area vs grey-matter volume across species (the "R = [XX], p = [XX], 95 % CI"
     of the Results). Input: the sheet `values_by_species` of Figure_2_stats.xlsx (median per species and lobe).
       - 35 points (7 species x 5 lobes), log10-log10, Pearson R; 95 % CI and P by bootstrap / permutation of
         SPECIES (lobes of one species are not independent, so species are resampled as blocks);
       - 7 points (sum of the five lobe medians per species), Pearson R with Fisher CI, exact Spearman rho.
       Species are not independent either (shared ancestry): this is a descriptive correlation, not a PGLS.
  2. Fig. 6c: Spearman r between the mean profiles of two modalities over the 15 regions, with a permutation P
     (regions permuted, 10,000) and a 95 % bootstrap CI (regions resampled, 10,000); Holm correction over the 6
     correlations. Input: sheets C_<species>_<modality> of Figure_6_stats.xlsx.

    python paper_extra_stats.py [Figure_2_stats.xlsx] [Figure_6_stats.xlsx]
Writes <PAPER_DIR>/Supplementary_stats/extra_stats.txt when PAPER_DIR is available.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

N_RES = 10000
SEED = 0


def _pearson_ci_fisher(r, n, z=1.959964):
    zr = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    return float(np.tanh(zr - z * se)), float(np.tanh(zr + z * se))


def fig2_scaling(xlsx):
    v = pd.read_excel(xlsx, sheet_name="values_by_species")
    s = v[v.modality == "surface"].set_index(["species", "region"])["median"]
    g = v[v.modality == "volume"].set_index(["species", "region"])["median"]
    d = pd.concat([np.log10(s).rename("s"), np.log10(g).rename("v")], axis=1).dropna().reset_index()
    sp = list(d.species.unique())
    rng = np.random.default_rng(SEED)
    r = stats.pearsonr(d.s, d.v)[0]
    boot = []
    for _ in range(N_RES):
        pick = rng.choice(sp, len(sp))
        b = pd.concat([d[d.species == k] for k in pick])
        if b.s.std() > 0 and b.v.std() > 0:
            boot.append(np.corrcoef(b.s, b.v)[0, 1])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    # permutation: species-level blocks of (s) re-assigned to species-level blocks of (v)
    blocks_s = [d[d.species == k].s.values for k in sp]
    blocks_v = [d[d.species == k].v.values for k in sp]
    null = []
    for _ in range(N_RES):
        pm = rng.permutation(len(sp))
        null.append(np.corrcoef(np.concatenate(blocks_s), np.concatenate([blocks_v[i] for i in pm]))[0, 1])
    p_perm = float((np.sum(np.array(null) >= r) + 1) / (N_RES + 1))
    # species totals
    tot = d.assign(s=10 ** d.s, v=10 ** d.v).groupby("species")[["s", "v"]].sum()
    ls, lv = np.log10(tot.s), np.log10(tot.v)
    r7 = stats.pearsonr(ls, lv)[0]
    lo7, hi7 = _pearson_ci_fisher(r7, len(tot))
    rho7 = stats.spearmanr(ls, lv)[0]
    # exact one-sided permutation P of the 7-species Spearman
    nullr = [stats.spearmanr(ls.values, lv.values[list(pm)])[0] for pm in itertools.permutations(range(len(tot)))]
    p7 = float(np.mean(np.array(nullr) >= rho7 - 1e-12))
    slope = stats.linregress(d.s, d.v).slope
    return [
        "Fig. 2: isocortical surface area vs grey-matter volume across species (log10 of the median per species and lobe)",
        f"  {len(d)} points ({len(sp)} species x {d.region.nunique()} lobes): Pearson R = {r:.3f}; species-block bootstrap "
        f"95% CI [{lo:.3f}, {hi:.3f}] ({N_RES} resamples); species-block permutation one-sided P = {p_perm:.2g}; "
        f"log-log slope (volume on surface) = {slope:.2f}",
        f"  {len(tot)} species totals (sum of the five lobes): Pearson R = {r7:.3f}, Fisher 95% CI [{lo7:.3f}, {hi7:.3f}]; "
        f"Spearman rho = {rho7:.3f}, exact one-sided P = {p7:.2g} ({len(nullr)} permutations)",
        "  descriptive: species share ancestry (no phylogenetic correction); lobes of a species are not independent",
    ]


def _holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for rank, i in enumerate(o):
        run = max(run, (len(p) - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj


def fig6c_correlations(xlsx):
    rng = np.random.default_rng(SEED)
    x = pd.read_excel(xlsx, sheet_name=None)
    rows = []
    import itertools as it
    for sp in ("Human", "Macaque"):
        prof = {m: x[f"C_{sp}_{m}"[:31]].drop(columns="subject").mean()
                for m in ("morphometry", "FC", "Q80", "Q20", "metabolism") if f"C_{sp}_{m}"[:31] in x}
        for a, b in it.combinations(list(prof), 2):
            u, w = prof[a].values, prof[b].values
            r = stats.spearmanr(u, w)[0]
            null = np.array([stats.spearmanr(u, rng.permutation(w))[0] for _ in range(N_RES)])
            p2 = float((np.sum(np.abs(null) >= abs(r) - 1e-12) + 1) / (N_RES + 1))
            bs = []
            for _ in range(N_RES):
                i = rng.integers(0, len(u), len(u))
                if len(set(u[i])) > 1 and len(set(w[i])) > 1:
                    bs.append(stats.spearmanr(u[i], w[i])[0])
            lo, hi = np.percentile(bs, [2.5, 97.5])
            rows.append((sp, a, b, len(u), r, lo, hi, p2))
    adj = _holm([r[-1] for r in rows])
    out = ["Fig. 6c: Spearman r between mean regional profiles (descriptive; regions permuted for P, resampled for the CI; "
           f"{N_RES} each; two-sided P; Holm over all the correlations)"]
    for (sp, a, b, n, r, lo, hi, p2), pa in zip(rows, adj):
        out.append(f"  {sp} {a}-{b}: r = {r:+.2f} [{lo:+.2f}, {hi:+.2f}], n = {n} regions, P = {p2:.3f}, Holm P = {pa:.3f}")
    return out


def main(argv):
    try:
        from Exemples.Study_EDNiX.analysis import paper_data as pd_
        base = pd_.PAPER_DIR
    except Exception:                                   # noqa: BLE001
        base = None
    f2 = argv[1] if len(argv) > 1 else os.path.join(base, "Figure_2", "Figure_2_stats.xlsx")
    f6 = argv[2] if len(argv) > 2 else os.path.join(base, "Figure_6", "Figure_6_stats.xlsx")
    lines = fig2_scaling(f2) + [""] + fig6c_correlations(f6)
    print("\n".join(lines))
    if base:
        out = os.path.join(base, "Supplementary_stats")
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "extra_stats.txt"), "w") as f:
            f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv)
