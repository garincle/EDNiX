"""
phylo_power_species.py -- how many species are needed for the Spearman Mantel test
(functional distance between species vs MRCA age) to reach P < 0.05?

Model (Brownian motion on a time-calibrated tree, standard for comparative data;
Felsenstein 1985; Harmon & Glor 2010 for Mantel power on phylogenetic data):
    x_s = sqrt(a) * c + sqrt(1 - a) * ( sqrt(w) * p_s + sqrt(1 - w) * e_s )
  c    common profile shared by all species (the consensus network)
  p_s  phylogenetic component: Brownian motion on the tree, cov(p_s, p_t) = shared
       branch length / root age
  e_s  species-specific, non-phylogenetic part (biology + dataset/measurement noise)
  a    share of the common profile (A_COMMON, set so the mean between-species r is ~0.5)
  w    phylogenetic share of the rest, CALIBRATED so that the median Spearman Mantel rho
       on the paper's 7-species tree equals the observed rho (TARGETS)
Functional distance = 1 - Pearson r between species profiles of N_VALUES regional values
(as in Fig. 4). Mantel: one-sided, species-label permutations (exact for 7 species).

Larger samples: random pure-birth (Yule) trees scaled to the same root age (94 Mya),
i.e. N mammal species sampled without any preference for close relatives.

Writes, next to this script or in <PAPER_DIR>/Validation when it exists:
    phylo_power_species.csv   power per N and target rho
    phylo_power_species.png   power curves
    phylo_power_species.txt   sentences for the paper

    python phylo_power_species.py
"""
import itertools
import math
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

opj = os.path.join

SEED = 20261102
N_VALUES = 30                 # regional values per profile (Fig. 4: 30 kept values)
A_COMMON = 0.5                # share of the common profile
ROOT_AGE = 94.0               # Mya, Boreoeutheria (paper tree)
TARGETS = {                   # observed Spearman Mantel rho (Fig. 4)
    "selected runs, all values (rho = 0.29)": 0.29,
    "selected runs, held-out values (rho = 0.15)": 0.15,
    "optimistic (rho = 0.40)": 0.40,
}
N_SPECIES = [7, 9, 12, 15, 20, 25, 30, 40]
N_SIM = 1000                  # simulated datasets per (N, target)
N_PERM = 999                  # Mantel permutations (N > 7)
N_CALIB = 2000                # simulations per grid point for the calibration
ALPHA = 0.05

# Paper tree (fmri_fingerprints.phylo.SPECIES_CLADES; Spearman Mantel uses only the order)
PAPER_SPECIES = ["Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse", "Dog"]
PAPER_CLADES = [
    (("Rat", "Mouse"), 13.0),
    (("Human", "Macaque"), 29.0),
    (("Human", "Macaque", "Marmoset"), 43.0),
    (("Human", "Macaque", "Marmoset", "Mouselemur"), 74.0),
    (("Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse"), 87.0),
    (("Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse", "Dog"), 94.0),
]


def out_dir():
    paper = "/scratch2/EDNiX/results/EDNiX_Figures/Paper/Validation"
    if os.path.isdir(os.path.dirname(paper)):
        os.makedirs(paper, exist_ok=True)
        return paper
    return os.path.dirname(os.path.abspath(__file__))


def paper_mrca():
    n = len(PAPER_SPECIES)
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        a, b = PAPER_SPECIES[i], PAPER_SPECIES[j]
        age = min(t for members, t in PAPER_CLADES if a in members and b in members)
        D[i, j] = D[j, i] = age
    return D


def yule_mrca(n, rng):
    """MRCA-age matrix of a random pure-birth (Yule) tree with n tips, root scaled to ROOT_AGE."""
    parent, time = [None], [0.0]          # node 0 = root, split at t = 0
    leaves = []
    for _ in range(2):
        parent.append(0); time.append(None); leaves.append(len(parent) - 1)
    t = 0.0
    while len(leaves) < n:
        t += rng.exponential(1.0 / len(leaves))
        k = int(rng.integers(len(leaves)))
        node = leaves.pop(k)
        time[node] = t                      # this lineage splits now
        for _ in range(2):
            parent.append(node); time.append(None); leaves.append(len(parent) - 1)
    t_end = t + rng.exponential(1.0 / len(leaves))

    def ancestors(i):
        out = []
        while i is not None:
            out.append(i); i = parent[i]
        return out

    anc = [ancestors(l) for l in leaves]
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        sj = set(anc[j])
        m = next(x for x in anc[i][1:] if x in sj)
        D[i, j] = D[j, i] = t_end - time[m]
    return D * (ROOT_AGE / t_end)


def bm_cholesky(D):
    """Cholesky factor of the Brownian covariance (shared branch length / root age)."""
    root = D.max()
    C = (root - D) / root
    np.fill_diagonal(C, 1.0)
    return np.linalg.cholesky(C + 1e-10 * np.eye(len(C)))


def simulate_profiles(L, w, rng, n_sets):
    """n_sets x n_species x N_VALUES profiles."""
    n = L.shape[0]
    c = rng.standard_normal((n_sets, 1, N_VALUES))
    p = np.einsum("ij,sjk->sik", L, rng.standard_normal((n_sets, n, N_VALUES)))
    e = rng.standard_normal((n_sets, n, N_VALUES))
    return math.sqrt(A_COMMON) * c + math.sqrt(1 - A_COMMON) * (
        math.sqrt(w) * p + math.sqrt(1 - w) * e)


def func_dist(X):
    """1 - Pearson r between rows, for each set: n_sets x n x n."""
    Z = X - X.mean(axis=2, keepdims=True)
    Z /= np.linalg.norm(Z, axis=2, keepdims=True)
    return 1.0 - np.einsum("sik,sjk->sij", Z, Z)


def rank_rows(M):
    """Average ranks (ties averaged) along the last axis."""
    return rankdata(M, axis=-1)


def _std(v):
    v = v - v.mean(axis=-1, keepdims=True)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def rank_matrix(D):
    """Matrix whose upper-triangle entries are the ranks of D's upper triangle. Permuting
    species labels only moves the values between pairs, so the ranks of a permuted matrix
    are the permuted ranks: ranking once is enough."""
    n = D.shape[0]
    iu = np.triu_indices(n, 1)
    R = np.zeros_like(D, dtype=float)
    R[iu] = rankdata(D[iu])
    return R + R.T


def mantel_power(Dtrees, Dfuncs, rng, exact):
    """One-sided Spearman Mantel for each (tree, functional matrix) pair. Returns rho, P."""
    n = Dtrees.shape[1]
    iu = np.triu_indices(n, 1)
    rho, pval = [], []
    if exact:
        perms = np.array(list(itertools.permutations(range(n))))
    for Dt, Df in zip(Dtrees, Dfuncs):
        x = _std(rankdata(Df[iu]))
        if not exact:
            perms = np.vstack([np.arange(n)[None],
                               np.array([rng.permutation(n) for _ in range(N_PERM)])])
        R = rank_matrix(Dt)
        yr = _std(R[perms[:, :, None], perms[:, None, :]][:, iu[0], iu[1]])
        null = yr @ x
        obs = null[0]
        if exact:
            p = float(np.mean(null >= obs - 1e-12))
        else:
            p = (1 + int(np.sum(null[1:] >= obs - 1e-12))) / (N_PERM + 1)
        rho.append(obs); pval.append(p)
    return np.array(rho), np.array(pval)


def spearman_only(Dt, Dfuncs):
    n = Dt.shape[0]
    iu = np.triu_indices(n, 1)
    y = _std(rank_rows(Dt[iu][None])[0])
    x = _std(rank_rows(Dfuncs[:, iu[0], iu[1]]))
    return x @ y


def calibrate(target, D7, rng):
    L = bm_cholesky(D7)
    grid = np.linspace(0.0, 0.99, 34)
    med = []
    for w in grid:
        X = simulate_profiles(L, w, rng, N_CALIB)
        med.append(np.median(spearman_only(D7, func_dist(X))))
    med = np.array(med)
    if target > med.max():
        return grid[-1], med.max()
    k = np.searchsorted(med, target)          # med increases with w
    k = min(max(k, 1), len(grid) - 1)
    w = np.interp(target, med[k - 1:k + 1], grid[k - 1:k + 1])
    return float(w), float(target)


def main():
    rng = np.random.default_rng(SEED)
    D7 = paper_mrca()
    rows = []
    for label, target in TARGETS.items():
        w, achieved = calibrate(target, D7, rng)
        print(f"[{label}] phylogenetic share w = {w:.3f} (median rho on paper tree {achieved:+.3f})")
        for n in N_SPECIES:
            if n == 7:
                Dtrees = np.repeat(D7[None], N_SIM, axis=0)
            else:
                Dtrees = np.stack([yule_mrca(n, rng) for _ in range(N_SIM)])
            X = np.stack([simulate_profiles(bm_cholesky(Dt), w, rng, 1)[0] for Dt in Dtrees])
            rho, p = mantel_power(Dtrees, func_dist(X), rng, exact=(n == 7))
            rows.append(dict(scenario=label, target_rho_7sp=target, phylo_share_w=w,
                             n_species=n, power=float(np.mean(p < ALPHA)),
                             median_rho=float(np.median(rho)), median_p=float(np.median(p))))
            print(f"   N = {n:2d}: power {rows[-1]['power']:.2f}, median rho {rows[-1]['median_rho']:+.2f}")
    df = pd.DataFrame(rows)

    od = out_dir()
    df.to_csv(opj(od, "phylo_power_species.csv"), index=False)

    fig, ax = plt.subplots(figsize=(3.5, 2.6), dpi=300)
    for label, g in df.groupby("scenario", sort=False):
        ax.plot(g["n_species"], g["power"], marker="o", ms=3, lw=1, label=label)
    ax.axhline(0.8, color="0.5", lw=0.6, ls="--")
    ax.set_xlabel("Number of species", fontsize=7)
    ax.set_ylabel("Power (one-sided P < 0.05)", fontsize=7)
    ax.set_ylim(0, 1)
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=5, frameon=False, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(opj(od, "phylo_power_species.png"))

    lines = [f"Spearman Mantel power, Brownian-motion model calibrated on the 7-species tree "
             f"({N_SIM} simulations per point, {N_VALUES} regional values, alpha = {ALPHA}).", ""]
    for label, g in df.groupby("scenario", sort=False):
        ok = g[g["power"] >= 0.8]
        n80 = int(ok["n_species"].min()) if len(ok) else None
        p7 = float(g.loc[g["n_species"] == 7, "power"].iloc[0])
        lines.append(f"{label}: power with 7 species = {p7:.2f}; 80% power reached with "
                     + (f"{n80} species." if n80 else f"> {max(N_SPECIES)} species."))
    with open(opj(od, "phylo_power_species.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"-> {od}")


if __name__ == "__main__":
    main()
