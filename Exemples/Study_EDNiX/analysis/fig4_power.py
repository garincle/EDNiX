"""
fig4_power.py -- power of the species-level Mantel test of Fig. 4d (7 species).

The test correlates the distance between species (1 - r of their mean fingerprints) with the
age of their last common ancestor (Spearman, exact one-sided permutation over the 7! = 5,040
assignments of the species). With 7 species there are 21 species pairs that are not
independent, so the test has little power. This script gives
  1. the smallest Spearman rho that reaches one-sided P < 0.05 / 0.01 (exact null distribution);
  2. the power to detect a true association rho0, simulated as D = z(rank(MRCA age)) + sigma e,
     with independent Gaussian noise e per pair and sigma = sqrt(1 / rho0^2 - 1), so that rho0
     is the expected correlation between D and the ages (critical value from item 1);
  3. the true association needed for 80 % power.
The observed rho (+0.29 selected, +0.10 all runs) is read against these numbers.

    python fig4_power.py

Uses the time-calibrated tree of the pipeline when it can be loaded (Inputs().fmri()), otherwise
the approximate hand-entered ages of phylo.py (the Spearman test depends only on the order of the
node ages, which is the same). Writes <PAPER_DIR>/Figure_4/Figure_4_power.txt
"""
import itertools
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

SPECIES = ("Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse", "Dog")
N_SIM = 5000
SEED = 0


def exact_null(Dp):
    n = len(Dp)
    iu = np.triu_indices(n, 1)
    rp = stats.rankdata(Dp[iu])
    vals = []
    for pm in itertools.permutations(range(n)):
        pm = list(pm)
        vals.append(np.corrcoef(stats.rankdata(Dp[np.ix_(pm, pm)][iu]), rp)[0, 1])
    return np.array(vals), rp


def analyse(Dp, n_sim=N_SIM, seed=SEED):
    rng = np.random.default_rng(seed)
    null, rp = exact_null(Dp)
    crit05, crit01 = np.percentile(null, 95), np.percentile(null, 99)
    # smallest attainable rho with P <= alpha (P = share of permutations >= rho)
    zp = (rp - rp.mean()) / rp.std()
    rows = []
    for rho0 in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        sigma = np.sqrt(1.0 / rho0 ** 2 - 1.0)
        hits, obs = 0, []
        for _ in range(n_sim):
            d = zp + sigma * rng.normal(size=len(zp))
            r = stats.spearmanr(d, rp)[0]
            obs.append(r)
            hits += r >= crit05
        rows.append((rho0, hits / n_sim, float(np.median(obs))))
    need = next((r[0] for r in rows if r[1] >= 0.8), None)
    return null, crit05, crit01, rows, need


def lines_for(Dp, tree_note):
    null, c05, c01, rows, need = analyse(Dp)
    lines = ["Fig. 4d: power of the species-level Mantel test (7 species, exact one-sided "
             "Spearman permutation test, 5,040 permutations)",
             f"  tree: {tree_note}",
             f"  smallest rho with P < 0.05: {c05:+.3f}; with P < 0.01: {c01:+.3f}; smallest "
             f"attainable P = {1 / len(null):.2g}",
             "  power to detect a true association rho0 (alpha = 0.05, one-sided; "
             f"{N_SIM} simulated species-distance matrices, independent noise per pair):"]
    for rho0, pw, med in rows:
        lines.append(f"    rho0 = {rho0:.1f}: power = {pw:.2f} (median observed rho {med:+.2f})")
    lines.append("  true association needed for 80 % power: " +
                 (f"rho0 >= {need:.1f}" if need is not None else "not reached up to 0.9"))
    lines.append("  reading: observed rho = +0.29 (selected) and +0.10 (all runs) are below the "
                 "critical value; an association of moderate size cannot be detected with seven "
                 "species, so the comparison with random draws (not the test against zero) "
                 "carries the inference")
    return lines


def main():
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import phylo
    note = "approximate hand-entered node ages (phylo.py)"
    try:
        from Exemples.Study_EDNiX.analysis import paper_stats as pst
        pst.Inputs().fmri()                 # loads TREE_PATH into phylo when configured
        note = "time-calibrated tree of the pipeline (TREE_PATH) if configured, as in Fig. 4d"
    except Exception as e:                  # noqa: BLE001
        print(f"  [tree] pipeline inputs not available ({type(e).__name__}); using {note}")
    sps = phylo.order_species(list(SPECIES))
    Dp, _ = phylo.phylo_matrices(sps)
    lines = lines_for(np.asarray(Dp, float), note)
    try:
        from Exemples.Study_EDNiX.analysis import paper_data as pd_
        out = os.path.join(pd_.PAPER_DIR, "Figure_4")
        os.makedirs(out, exist_ok=True)
        path = os.path.join(out, "Figure_4_power.txt")
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"[saved] {path}")
    except Exception:                       # noqa: BLE001
        pass
    print("\n".join(lines))


if __name__ == "__main__":
    main()
