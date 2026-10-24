"""
General linear model of the PET group analyses (voxel-wise and regional), in numpy.

One image (or one regional value) per subject: y = X b + e, a t contrast c, and non-parametric inference
by permutation (Freedman & Lane 1983, recommended by Winkler et al. 2014 and Anderson & Robinson 2001):
the model is split into the effect of interest and the nuisance part Z = X C0 (C0: basis of the null space
of c, Ridgway 2009); the residuals of the nuisance-only model are permuted (two-sample, regression) or
sign-flipped (one-sample, paired differences), added back to its fit, and the full model is refitted.
The maximum statistic over the brain gives family-wise error rates (Nichols & Holmes 2002): maximum |t|
for voxels, maximum cluster size for clusters (cluster-forming threshold two-sided, face connectivity, as
the AFNI NN1 2-sided tables of the fMRI group analyses). Eklund et al. (2016) showed that this
permutation test controls the family-wise error rate where parametric cluster inference did not.
"""
import numpy as np
import pandas as pd
from scipy import ndimage, stats

DESIGNS = ('one-sample', 'two-sample', 'regression', 'paired')


# ---------------------------------------------------------------------------------------------- design
def design_matrix(table, design, group_col=None, groups=None, covariates=None, test_covariate=None):
    """
    Design matrix of one row per observation of `table` (DataFrame).

    design          'one-sample'  intercept (+ covariates); contrast: intercept
                    'two-sample'  one column per group of `groups` (cell means, = groups[0] - groups[1])
                                  (+ covariates); contrast: groups[0] - groups[1]
                    'regression'  intercept + covariates; contrast: `test_covariate` (first covariate)
                    'paired'      rows are within-subject differences: as 'one-sample'
    covariates      columns of `table` (numeric, or text = categorical, coded with one column fewer than
                    its levels); numeric covariates are centred (the intercept / group means are then the
                    values at the sample mean of the covariates).
    Returns X (n x p), contrast (p,), column names, and a mask of the rows used (rows with NaN in a used
    column are dropped).
    """
    if design not in DESIGNS:
        raise ValueError(f'design={design!r}: expected one of {DESIGNS}')
    covariates = list(covariates or [])
    used = list(covariates) + ([group_col] if design == 'two-sample' else [])
    missing = [c for c in used if c not in table.columns]
    if missing:
        raise KeyError(f'column(s) {missing} not in the table of the subjects (columns: {list(table.columns)})')
    keep = np.ones(len(table), bool)
    for c in used:
        keep &= table[c].notna().to_numpy()
    if design == 'two-sample':
        if groups is None or len(groups) != 2:
            raise ValueError("design='two-sample' needs groups=[group A, group B] (contrast A - B)")
        keep &= table[group_col].astype(str).isin([str(g) for g in groups]).to_numpy()
    t = table[keep].reset_index(drop=True)
    cols, names = [], []
    if design == 'two-sample':
        for g in groups:
            cols.append((t[group_col].astype(str) == str(g)).to_numpy(float))
            names.append(f'{group_col}={g}')
        for g, col in zip(groups, cols[:2]):
            if col.sum() < 2:
                raise ValueError(f'group {g!r} has {int(col.sum())} observation(s): at least 2 are needed')
    else:
        cols.append(np.ones(len(t)))
        names.append('intercept')
    for c in covariates:
        v = t[c]
        if pd.api.types.is_numeric_dtype(v):
            cols.append(v.to_numpy(float) - v.to_numpy(float).mean())
            names.append(c)
        else:
            levels = sorted(v.astype(str).unique())
            for lev in levels[1:]:
                cols.append((v.astype(str) == lev).to_numpy(float))
                names.append(f'{c}={lev}')
    X = np.column_stack(cols)
    contrast = np.zeros(X.shape[1])
    if design == 'two-sample':
        contrast[:2] = [1, -1]
    elif design == 'regression':
        if not covariates:
            raise ValueError("design='regression' needs at least one covariate")
        tc = test_covariate or covariates[0]
        if tc not in names:
            raise ValueError(f'test_covariate={tc!r} must be a numeric covariate (columns: {names})')
        contrast[names.index(tc)] = 1
    else:
        contrast[0] = 1
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError(f'rank-deficient design (columns {names}): remove collinear covariates')
    if X.shape[0] - X.shape[1] < 2:
        raise ValueError(f'{X.shape[0]} observations for {X.shape[1]} parameters: not enough degrees of freedom')
    return X, contrast, names, keep


def contrast_name(design, contrast, names, groups=None):
    if design == 'two-sample':
        return f'{groups[0]}-vs-{groups[1]}'
    if design == 'regression':
        return f'effect-of-{names[int(np.flatnonzero(contrast)[0])]}'
    return 'mean' if design == 'one-sample' else 'paired-difference'


# ---------------------------------------------------------------------------------------------- fit
class GLM:
    """OLS fit of Y (n x V) for one t contrast, with the matrices reused by the permutations."""

    def __init__(self, X, contrast):
        self.X = np.asarray(X, float)
        self.c = np.asarray(contrast, float)
        n, p = self.X.shape
        self.df = n - p
        self.pinv = np.linalg.pinv(self.X)
        xtx_inv = self.pinv @ self.pinv.T
        self.var_c = float(self.c @ xtx_inv @ self.c)
        # nuisance part of the model: Z = X C0, C0 = basis of the null space of c (Ridgway 2009)
        vt = np.linalg.svd(self.c[None, :])[2]
        C0 = vt[1:].T
        Z = self.X @ C0 if C0.size else np.zeros((n, 0))
        self.Hz = Z @ np.linalg.pinv(Z) if Z.shape[1] else np.zeros((n, n))
        self.Rz = np.eye(n) - self.Hz

    def fit(self, Y):
        """cope (V,), t (V,), residuals (n x V)."""
        beta = self.pinv @ Y
        resid = Y - self.X @ beta
        sigma2 = (resid ** 2).sum(0) / self.df
        cope = self.c @ beta
        with np.errstate(divide='ignore', invalid='ignore'):
            t = cope / np.sqrt(sigma2 * self.var_c)
        t[~np.isfinite(t)] = 0.0
        return cope, t, resid

    def t_only(self, Y):
        return self.fit(Y)[1]


def t_to_z(t, df):
    """Signed z of the same one-sided p as t (log scale: exact far in the tails)."""
    t = np.asarray(t, float)
    logp = stats.t.logsf(np.abs(t), df)
    z = -stats.norm.ppf(np.exp(np.clip(logp, -745, 0)))
    big = logp < -700
    if np.any(big):
        # asymptotic inverse of the normal tail for p below the double precision range
        z[big] = np.sqrt(-2 * logp[big] - np.log(-2 * logp[big]) - np.log(2 * np.pi))
    return np.sign(t) * z


def t_threshold(voxel_p, df):
    """|t| of a two-sided voxel-wise p."""
    return float(stats.t.isf(voxel_p / 2, df))


# ---------------------------------------------------------------------------------------------- clusters
def _cluster_sizes_signed(t3d, thr):
    """Max cluster size with positive and negative voxels clustered separately (sign-consistent clusters)."""
    best = 0
    for sup in (t3d >= thr, t3d <= -thr):
        lab, n = ndimage.label(sup)
        if n:
            best = max(best, int(np.bincount(lab.ravel())[1:].max()))
    return best


# ---------------------------------------------------------------------------------------------- permutations
def permutation_test(glm, Y, mask3d, voxel_p=0.001, n_perm=5000, sign_flip=False, seed=0, progress=None):
    """
    Freedman-Lane permutation (or sign flipping) of the t contrast of `glm` on Y (n x V, V = voxels of
    mask3d in C order). Returns dict: t (observed), max_t (n_perm,), max_cluster (n_perm,), threshold
    (|t| forming the clusters), n_perm, exhaustive (bool).
    The first permutation is the identity (Phipson & Smyth 2010: p = (b + 1) / (m + 1) counted with it).
    Clusters: |t| >= t(voxel_p / 2), positive and negative voxels in separate clusters.
    """
    n = Y.shape[0]
    rng = np.random.default_rng(seed)
    fitted_z = glm.Hz @ Y
    resid_z = glm.Rz @ Y
    thr = t_threshold(voxel_p, glm.df)
    vol = np.zeros(mask3d.shape, np.float32)
    idx = np.flatnonzero(mask3d.ravel())

    shuffles, exhaustive = _shuffles(n, n_perm, sign_flip, rng)
    max_t = np.empty(len(shuffles))
    max_c = np.empty(len(shuffles))
    t_obs = None
    for k, s in enumerate(shuffles):
        Ystar = (resid_z * s[:, None] if sign_flip else resid_z[s]) + fitted_z
        t = glm.t_only(Ystar)
        if k == 0:
            t_obs = t
        max_t[k] = np.abs(t).max() if t.size else 0
        vol.ravel()[idx] = t
        max_c[k] = _cluster_sizes_signed(vol, thr)
        if progress and k % 500 == 0:
            progress(k, len(shuffles))
    return dict(t=t_obs, max_t=max_t, max_cluster=max_c, threshold=thr, n_perm=len(shuffles),
                exhaustive=exhaustive)


def _shuffles(n, n_perm, sign_flip, rng):
    """Identity first, then random permutations / sign flips; every one of them when there are fewer than
    n_perm (exhaustive test)."""
    from math import factorial
    total = 2 ** n if sign_flip else factorial(n)
    if total <= n_perm:
        import itertools
        if sign_flip:
            out = [np.array(s, float) for s in itertools.product([1.0, -1.0], repeat=n)]
        else:
            out = [np.array(p) for p in itertools.permutations(range(n))]
        return out, True
    out = [np.ones(n) if sign_flip else np.arange(n)]
    for _ in range(n_perm - 1):
        out.append(rng.choice([-1.0, 1.0], n) if sign_flip else rng.permutation(n))
    return out, False


def fwe_p(observed, null_max):
    """FWE-corrected p of every observed value against the null maximum distribution (identity included)."""
    null_max = np.sort(np.asarray(null_max))
    observed = np.asarray(observed, float)
    ge = len(null_max) - np.searchsorted(null_max, observed, side='left')
    return ge / len(null_max)


def critical_value(null_max, alpha, integer=False):
    """Smallest value v with FWE p(v) = #(null maxima >= v) / m <= alpha: every value above the
    (m - floor(alpha m))-th smallest null maximum. integer=True (cluster sizes): that maximum + 1."""
    s = np.sort(np.asarray(null_max, float))
    m = len(s)
    allowed = int(np.floor(alpha * m))
    if allowed == 0:
        return np.inf
    edge = s[m - allowed - 1]
    return edge + 1 if integer else float(np.nextafter(edge, np.inf))


# ---------------------------------------------------------------------------------------------- tables
def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (NaN kept)."""
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    q = p[ok]
    order = np.argsort(q)
    m = len(q)
    adj = q[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.minimum(adj, 1)
    out[ok] = res
    return out


def parametric_p(t, df):
    """Two-sided p of t."""
    return 2 * stats.t.sf(np.abs(t), df)
