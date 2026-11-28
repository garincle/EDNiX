"""Statistics of the Fig. 4 map analyses (no I/O).

Profiles are dataset vectors (one value per P80/P20 channel of the common
regions, or per region / edge for the amplitude analyses); `groups` is always
{dataset index: vector} and `sp_of_ds[d]` the species of dataset d.

Phylogenetic tests on the species means (one row per species, `sps` order):
  K_mult   multivariate Blomberg K (Adams 2014, Syst Biol 63:685), ratio to its
           Brownian expectation; P = share of the 7! = 5040 species
           permutations with K >= observed (exact, one-sided).
  Mantel   Spearman rho between a species distance matrix and the MRCA age
           matrix; exact P over the same permutations.
Residuals are taken against the GENERAL common profile (median of the species
medians of the all-runs maps), leave-one-species-out (LOSO) for the species of
the dataset, fixed for every group of runs (selected, random draw, noise bin)
so that a group is never compared with a reference built from itself.
"""
import itertools

import numpy as np
from scipy import stats as _st


def z(v):
    v = np.asarray(v, float)
    return (v - np.nanmean(v)) / np.nanstd(v)


def corr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def resid(g, C):
    """Residual of g after regression on the common profile C (slope only: both are
    z-scored, or C is the reference pattern)."""
    g, C = np.asarray(g, float), np.asarray(C, float)
    ok = np.isfinite(g) & np.isfinite(C)
    return g - np.cov(g[ok], C[ok])[0, 1] / np.var(C[ok], ddof=1) * C


def species_means(groups, sp_of_ds, sps):
    return np.vstack([np.nanmean([groups[d] for d in groups if sp_of_ds[d] == s], axis=0)
                      for s in sps])


def common(groups, sp_of_ds, species, zscore=True):
    """Median over `species` of the species medians of the dataset vectors."""
    c = np.nanmedian([np.nanmedian([groups[d] for d in groups if sp_of_ds[d] == s], axis=0)
                      for s in species], axis=0)
    return z(c) if zscore else c


def commons(groups_all, sp_of_ds, sps, zscore=True):
    """General common profile Cg and its LOSO versions Cs[s] (species s left out)."""
    Cg = common(groups_all, sp_of_ds, sps, zscore)
    Cs = {s: common(groups_all, sp_of_ds, [t for t in sps if t != s], zscore) for s in sps}
    return Cg, Cs


def species_profiles(groups, sp_of_ds, sps, Cg, Cs):
    """Species x value matrices: raw, residual vs Cg, residual vs the LOSO common of the
    dataset's species (resLOSO), and the alignment r(species mean, LOSO common)."""
    res = {d: resid(g, Cg) for d, g in groups.items()}
    loso = {d: resid(g, Cs[sp_of_ds[d]]) for d, g in groups.items()}
    align = [corr(species_means({d: g for d, g in groups.items() if sp_of_ds[d] == s},
                                sp_of_ds, [s])[0], Cs[s]) for s in sps]
    return dict(raw=species_means(groups, sp_of_ds, sps), res=species_means(res, sp_of_ds, sps),
                resLOSO=species_means(loso, sp_of_ds, sps), align=np.array(align))


def dist_corr(Y):
    """1 - Pearson r between the rows of Y."""
    n = len(Y)
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        D[i, j] = D[j, i] = 1.0 - corr(Y[i], Y[j])
    return D


def dist_euc(Y):
    Y = np.where(np.isfinite(Y), Y, 0.0)
    return np.sqrt(((Y[:, None, :] - Y[None, :, :]) ** 2).sum(-1))


def p_high(null, obs):
    null = np.asarray(null, float)
    return float((1 + np.sum(null >= obs - 1e-12)) / (null.size + 1))


def p_low(null, obs):
    null = np.asarray(null, float)
    return float((1 + np.sum(null <= obs + 1e-12)) / (null.size + 1))


class Tests:
    """Exact permutation tests over the species (n! permutations of the rows).
    Dp: MRCA age matrix, Cphy: phylogenetic covariance (phylo.phylo_matrices)."""

    def __init__(self, Dp, Cphy):
        self.n = len(Dp)
        self.iu = np.triu_indices(self.n, 1)
        self.Dp = np.asarray(Dp, float)
        self.rp = _st.rankdata(self.Dp[self.iu])
        self.perms = [list(p) for p in itertools.permutations(range(self.n))]
        self.Ci = np.linalg.inv(Cphy)
        one = np.ones(self.n)
        self.k_exp = (np.trace(Cphy) - self.n / (one @ self.Ci @ one)) / (self.n - 1)

    def _rho(self, D, p):
        return float(np.corrcoef(_st.rankdata(D[np.ix_(p, p)][self.iu]), self.rp)[0, 1])

    def mantel(self, D, with_p=True, keep_null=False):
        obs = self._rho(D, list(range(self.n)))
        if not with_p:
            return obs, None
        null = np.array([self._rho(D, p) for p in self.perms])
        out = (obs, float(np.mean(null >= obs - 1e-12)))
        return out + (null,) if keep_null else out

    def _k(self, Y):
        one = np.ones(self.n)
        a = (one @ self.Ci @ Y) / (one @ self.Ci @ one)
        E = Y - np.outer(one, a)
        return (np.trace(E.T @ E) / np.trace(E.T @ self.Ci @ E)) / self.k_exp

    def kmult(self, Y, with_p=True, keep_null=False):
        Y = np.where(np.isfinite(Y), np.asarray(Y, float), 0.0)
        if Y.ndim == 1:
            Y = Y[:, None]
        obs = self._k(Y)
        if not with_p:
            return obs, None
        null = np.array([self._k(Y[p]) for p in self.perms])
        out = (obs, float(np.mean(null >= obs - 1e-12)))
        return out + (null,) if keep_null else out

    def pattern(self, M, with_p=True, keep_null=False):
        """Measures of the P80/P20 (z-scored) profiles: Mantel on 1 - r of the raw,
        residual and LOSO-residual species profiles, Mantel on the alignment
        difference, K_mult of the raw and LOSO-residual profiles."""
        a = M["align"]
        kw = dict(with_p=with_p, keep_null=keep_null)
        return dict(K_raw=self.kmult(M["raw"], **kw), K_resLOSO=self.kmult(M["resLOSO"], **kw),
                    M_raw=self.mantel(dist_corr(M["raw"]), **kw),
                    M_res=self.mantel(dist_corr(M["res"]), **kw),
                    M_resLOSO=self.mantel(dist_corr(M["resLOSO"]), **kw),
                    M_align=self.mantel(np.abs(a[:, None] - a[None, :]), **kw))

    def amplitude(self, G, sp_of_ds, sps, Cs, with_p=True):
        """Measures of unnormalised profiles (regional mean Fisher z, or edges):
        K_raw (amplitude and pattern), K_resSub (g - LOSO common: keeps amplitude),
        K_resReg (regression residual: removes the scaled common pattern), K_amp
        (Blomberg K of the RMS alone), Mantel Euclidean (raw, g - LOSO) and 1 - r."""
        raw = species_means(G, sp_of_ds, sps)
        sub = species_means({d: g - Cs[sp_of_ds[d]] for d, g in G.items()}, sp_of_ds, sps)
        reg = species_means({d: resid(g, Cs[sp_of_ds[d]]) for d, g in G.items()}, sp_of_ds, sps)
        amp = species_means({d: np.array([rms(g)]) for d, g in G.items()}, sp_of_ds, sps)
        return dict(K_raw=self.kmult(raw, with_p), K_resSub=self.kmult(sub, with_p),
                    K_resReg=self.kmult(reg, with_p), K_amp=self.kmult(amp, with_p),
                    M_euc=self.mantel(dist_euc(raw), with_p),
                    M_eucRes=self.mantel(dist_euc(sub), with_p),
                    M_corr=self.mantel(dist_corr(raw), with_p))


def rms(g):
    g = np.asarray(g, float)
    return float(np.sqrt(np.nanmean(g ** 2)))


def amp_ratio(G, sp_of_ds):
    """Mean |amplitude difference| between two datasets of the SAME species divided by
    the mean between datasets of DIFFERENT species (amplitude = RMS of the vector).
    >= 1: two datasets of one species differ as much as two species (amplitude is a
    property of the protocol); << 1: amplitude is species-like."""
    a = {d: rms(g) for d, g in G.items()}
    ds_ = sorted(a)
    pairs = list(itertools.combinations(ds_, 2))
    w = [abs(a[i] - a[j]) for i, j in pairs if sp_of_ds[i] == sp_of_ds[j]]
    b = [abs(a[i] - a[j]) for i, j in pairs if sp_of_ds[i] != sp_of_ds[j]]
    if not w or not b:
        return np.nan
    return float(np.mean(w) / np.mean(b))


def identifiability(prof, sp_of_ds, sps, Cs):
    """prof {d: z profile}. Mean r same species / other dataset and between species,
    their gap, sibling hits (datasets whose most similar dataset is the other dataset
    of their species), on the raw profiles and on the residuals vs the fixed LOSO
    common; Dsp = species distance matrix (1 - r of the species means)."""
    ds_ = sorted(prof)
    res = {d: resid(prof[d], Cs[sp_of_ds[d]]) for d in ds_}
    out = {}
    n = 0
    for tag, P in (("raw", prof), ("res", res)):
        same, betw = [], []
        for a, b in itertools.combinations(ds_, 2):
            (same if sp_of_ds[a] == sp_of_ds[b] else betw).append(corr(P[a], P[b]))
        out[f"r_same_{tag}"] = float(np.mean(same)) if same else np.nan
        out[f"r_between_{tag}"] = float(np.mean(betw))
        out[f"gap_{tag}"] = out[f"r_same_{tag}"] - out[f"r_between_{tag}"]
        hits, n = 0, 0
        for a in ds_:
            sib = [b for b in ds_ if b != a and sp_of_ds[b] == sp_of_ds[a]]
            if not sib:
                continue
            best = max((b for b in ds_ if b != a), key=lambda b: corr(P[a], P[b]))
            hits += best in sib
            n += 1
        out[f"sibling_{tag}"] = hits
    out["n_sibling"] = n
    out["Dsp"] = dist_corr(species_means(prof, sp_of_ds, sps))
    return out


IDENT_SCALARS = ("r_same_raw", "r_between_raw", "gap_raw", "sibling_raw",
                 "r_same_res", "r_between_res", "gap_res", "sibling_res")


def trend(values):
    """Spearman rho and least-squares slope of values on their index (1..B)."""
    b = np.arange(1, len(values) + 1)
    return float(_st.spearmanr(b, values)[0]), float(np.polyfit(b, values, 1)[0])
