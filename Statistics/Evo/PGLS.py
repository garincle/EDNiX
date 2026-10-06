"""
ednix_pgls_functions.py  —  TRUE PGLS with Pagel λ estimation
==============================================================
Phylogenetic Generalised Least Squares following Freckleton et al. (2002)
and Grafen (1989).  No R / rpy2 required.

Key difference from the previous version
-----------------------------------------
The old code fixed V = V_phylo (pure Brownian motion, λ=1 implicit).
This version estimates **Pagel's λ** by maximum likelihood, then fits GLS
at the optimal λ.  That means:

  V(λ) = λ · V_phylo  +  (1−λ) · diag(V_phylo)

  λ = 0  →  OLS   (no phylogenetic signal)
  λ = 1  →  pure Brownian-motion GLS
  0<λ<1  →  intermediate phylogenetic signal

The likelihood is profiled over σ² (concentrated out analytically), so
the optimisation is only 1-D over λ ∈ [0, 1].

Public API (backward-compatible with old ednix_pgls_functions.py)
------------------------------------------------------------------
load_newick_vcv(newick_path, species_map=None)
    → (vcv_df, tip_names)

pgls_intercept(species_vals, vcv_df)
    → dict  (intercept, se, t, p, lambda_pagel, n, df_residual)

pgls_regression(x_vals, y_vals, vcv_df, species,
                log_x=False, log_y=False, n_pred=200)
    → dict  (slope, intercept, r, p, se_slope, lambda_pagel,
             n_species, df_residual, x_fit, y_fit, ci_lo, ci_hi,
             species_used, log_x, log_y, ll, ll_ols)

run_pgls_all_region_pairs(df_morph, metric_col, vcv_df, ...)
    → DataFrame

plot_pgls_region_pair(result_row, df_morph, ...)
    → output_path

plot_pgls_summary_heatmap(df_results, regions, output_path, ...)
    → output_path

plot_pgls_dotplot(df_results, output_path, ...)
    → output_path

plot_brain_scaling_pgls(df_morph, metric_col, ..., vcv_df, ...)
    → output_path

ols_on_ax(ax, xv, yv, log_x=False)        [kept for back-compat]
pgls_on_ax(ax, xv, yv, vcv_df, species,
           log_x=False, log_y=False)       [PGLS-aware drop-in]
"""

import os
import itertools
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats as _stats
from scipy.optimize import minimize_scalar


# ─────────────────────────────────────────────────────────────────────────────
# Paper style
# ─────────────────────────────────────────────────────────────────────────────

PAPER_RC = {
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "font.family":       "sans-serif",
    "font.size":         11,
    "axes.titlesize":    12,
    "axes.labelsize":    11,
    "legend.frameon":    False,
    "figure.dpi":        150,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
}

PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#56B4E9", "#D55E00", "#F0E442", "#000000",
]

PHYLO_ORDER = [
    "Bat", "Rat", "Mouse", "Mouselemur",
    "Marmoset", "Macaque", "Human",
    "Dog", "Pig", "Cat",
]

DEFAULT_SPECIES_MAP = {
    "Mouse":      "Mus_musculus",
    "Rat":        "Rattus_norvegicus",
    "Bat":        "Phyllostomus_discolor",
    "Mouselemur": "Microcebus_murinus",
    "Marmoset":   "Callithrix_jacchus",
    "Macaque":    "Macaca_mulatta",
    "Human":      "Homo_sapiens",
    "Dog":        "Canis_lupus",
    "Pig":        "Sus_scrofa",
    "Rabbit":     "Oryctolagus_cuniculus",
    "Baboon":     "Papio_anubis",
    "Chimp":      "Pan_troglodytes",
}


# ─────────────────────────────────────────────────────────────────────────────
# PART 1 — Newick parser & VCV builder  (unchanged from original)
# ─────────────────────────────────────────────────────────────────────────────

class _Node:
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = float(length) if length else 0.0
        self.children = []

    def is_leaf(self):
        return len(self.children) == 0


def _parse_newick(s: str) -> _Node:
    s = s.strip().rstrip(";")

    def _parse(tok, pos):
        node = _Node()
        if pos < len(tok) and tok[pos] == "(":
            pos += 1
            while True:
                child, pos = _parse(tok, pos)
                node.children.append(child)
                if pos >= len(tok):
                    break
                if tok[pos] == ",":
                    pos += 1; continue
                if tok[pos] == ")":
                    pos += 1; break
        label_start = pos
        while pos < len(tok) and tok[pos] not in (",", ")", ":", ";"):
            pos += 1
        node.name = tok[label_start:pos].strip().strip("'\"")
        if pos < len(tok) and tok[pos] == ":":
            pos += 1
            bl_start = pos
            while pos < len(tok) and tok[pos] not in (",", ")", ";"):
                pos += 1
            try:
                node.length = float(tok[bl_start:pos])
            except ValueError:
                node.length = 0.0
        return node, pos

    root, _ = _parse(s, 0)
    return root


def _collect_leaves(node: _Node, depth: float = 0.0):
    depth += node.length
    if node.is_leaf():
        return [(node.name, depth)]
    out = []
    for child in node.children:
        out.extend(_collect_leaves(child, depth))
    return out


def _pairwise_shared(node: _Node, depth: float = 0.0):
    depth += node.length
    if node.is_leaf():
        return {}, {node.name: depth}
    child_results = [_pairwise_shared(c, depth) for c in node.children]
    all_pairs, all_leaves = {}, {}
    for pairs, leaves in child_results:
        all_pairs.update(pairs)
        all_leaves.update(leaves)
    for (_, leaves_a), (_, leaves_b) in itertools.combinations(
        [(p, l) for p, l in child_results], 2
    ):
        for la in leaves_a:
            for lb in leaves_b:
                all_pairs[frozenset([la, lb])] = depth
    return all_pairs, all_leaves


def load_newick_vcv(newick_path: str, species_map: dict = None) -> tuple:
    """
    Parse a Newick file → phylogenetic VCV matrix (Brownian-motion covariance).

    Returns
    -------
    vcv_df    : pd.DataFrame (species × species), Cov(i,j) = shared branch length
    tip_names : list of EDNiX species labels found in the tree
    """
    species_map = species_map or DEFAULT_SPECIES_MAP
    with open(newick_path, encoding="utf-8") as fh:
        nwk_str = fh.read().strip()

    root = _parse_newick(nwk_str)
    pairs, leaf_depths = _pairwise_shared(root)
    tip_names = list(leaf_depths.keys())

    tip_to_short = {}
    for short, pattern in species_map.items():
        for tip in tip_names:
            if tip.lower() == pattern.lower() or pattern.lower() in tip.lower():
                tip_to_short[tip] = short
                break

    ednix_names = [tip_to_short.get(t, t) for t in tip_names]
    n = len(ednix_names)
    vcv = np.zeros((n, n))
    for i, ti in enumerate(tip_names):
        vcv[i, i] = leaf_depths[ti]
        for j, tj in enumerate(tip_names):
            if i == j:
                continue
            vcv[i, j] = pairs.get(frozenset([ti, tj]), 0.0)

    vcv_df = pd.DataFrame(vcv, index=ednix_names, columns=ednix_names)
    print(f"  [VCV] {n} tips loaded from {Path(newick_path).name}")
    print(f"  [VCV] Species: {ednix_names}")
    return vcv_df, ednix_names


# ─────────────────────────────────────────────────────────────────────────────
# PART 2 — TRUE PGLS ENGINE with Pagel λ estimation
# ─────────────────────────────────────────────────────────────────────────────

def _make_vpagel(V: np.ndarray, lam: float) -> np.ndarray:
    """
    Build Pagel's λ-scaled covariance matrix.

      V(λ) = λ · V  +  (1−λ) · diag(V)

    diag(V) preserves variances (diagonal of V), zeroes off-diagonal covariances.
    """
    return lam * V + (1.0 - lam) * np.diag(np.diag(V))


def _gls_at_lambda(X: np.ndarray, y: np.ndarray,
                   V: np.ndarray, lam: float) -> dict | None:
    """
    Fit GLS at a fixed λ.  σ² is concentrated out analytically (profile ML).

    Returns a dict with keys:
        beta, se, t, p, sigma2, ll, df, resid, Vi, cov_beta
    or None on numerical failure.

    Log-likelihood uses the profile form:
        ℓ(λ) = −n/2 · log(σ̂²) − 1/2 · log|V(λ)|  +  const

    This is the standard used by ape::pgls and nlme::gls.
    """
    n, k = X.shape
    Vp = _make_vpagel(V, lam)

    # Cholesky for stability; fall back to direct inv on failure
    try:
        L = np.linalg.cholesky(Vp)
        Vi = np.linalg.solve(L.T, np.linalg.solve(L, np.eye(n)))
        logdet = 2.0 * np.sum(np.log(np.diag(L)))
    except np.linalg.LinAlgError:
        try:
            Vi = np.linalg.inv(Vp)
            sign, logdet = np.linalg.slogdet(Vp)
            if sign <= 0:
                return None
        except np.linalg.LinAlgError:
            return None

    XtVi = X.T @ Vi
    XtViX = XtVi @ X

    try:
        XtViX_inv = np.linalg.inv(XtViX)
    except np.linalg.LinAlgError:
        return None

    beta = XtViX_inv @ (XtVi @ y)
    resid = y - X @ beta
    df = n - k

    if df <= 0:
        return None

    # Concentrated (profile) σ²
    sigma2 = max(float(resid @ Vi @ resid) / df, 1e-15)

    # Profile log-likelihood  (constant −n/2·log(2π) omitted for optimisation)
    ll = -0.5 * (n * np.log(sigma2) + logdet + n)

    cov_beta = sigma2 * XtViX_inv
    se = np.sqrt(np.clip(np.diag(cov_beta), 0, None))
    t = np.where(se > 0, beta / se, np.nan)
    p = 2.0 * _stats.t.sf(np.abs(t), df=df)

    return dict(beta=beta, se=se, t=t, p=p,
                sigma2=sigma2, ll=ll, df=df,
                resid=resid, Vi=Vi, cov_beta=cov_beta,
                logdet=logdet)


def _estimate_lambda(X: np.ndarray, y: np.ndarray, V: np.ndarray,
                     n_grid: int = 20) -> tuple[float, float]:
    """
    Estimate Pagel's λ by maximising the profile log-likelihood over [0, 1].

    Strategy: coarse grid search to find a bracket, then bounded Brent
    minimisation.  This avoids getting stuck in flat regions near λ=0 or λ=1.

    Returns
    -------
    lambda_hat : float  — ML estimate of λ
    ll_max     : float  — maximised log-likelihood at λ_hat
    """
    def neg_ll(lam):
        r = _gls_at_lambda(X, y, V, lam)
        return -r["ll"] if r is not None else 1e10

    # Grid search for initialisation
    grid = np.linspace(0.0, 1.0, n_grid)
    ll_grid = np.array([neg_ll(l) for l in grid])
    best_idx = int(np.argmin(ll_grid))

    lo = max(0.0, grid[best_idx - 1] if best_idx > 0 else 0.0)
    hi = min(1.0, grid[best_idx + 1] if best_idx < n_grid - 1 else 1.0)

    res = minimize_scalar(neg_ll, bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-6, "maxiter": 500})

    # keep whichever is better: the bounded optimum or the best grid point
    if res.fun <= ll_grid[best_idx]:
        lam_hat = float(np.clip(res.x, 0.0, 1.0))
        ll_max = -float(res.fun)
    else:
        lam_hat = float(grid[best_idx])
        ll_max = -float(ll_grid[best_idx])

    return lam_hat, ll_max


def _subset_vcv(vcv_df: pd.DataFrame, species: list) -> tuple:
    """Return sub-matrix of vcv_df for the given species, warn on missing."""
    present = [s for s in species if s in vcv_df.index]
    missing = [s for s in species if s not in vcv_df.index]
    if missing:
        warnings.warn(f"  [PGLS] species not in VCV (dropped): {missing}")
    return vcv_df.loc[present, present], present


def _phylo_sort(species_list):
    known   = [s for s in PHYLO_ORDER if s in species_list]
    unknown = sorted(s for s in species_list if s not in PHYLO_ORDER)
    return known + unknown


# ─────────────────────────────────────────────────────────────────────────────
# PART 3 — PGLS intercept-only model
# ─────────────────────────────────────────────────────────────────────────────

def pgls_intercept(species_vals: dict, vcv_df: pd.DataFrame) -> dict | None:
    """
    PGLS intercept-only model:  y_species ~ μ

    Tests whether the phylogenetically corrected mean differs from zero.
    λ is estimated by ML.

    Parameters
    ----------
    species_vals : {species: value}  (NaN values dropped)
    vcv_df       : phylogenetic VCV DataFrame

    Returns
    -------
    dict: intercept, se, t, p, lambda_pagel, ll, n_species, df_residual
    or None when insufficient data.
    """
    clean = {sp: v for sp, v in species_vals.items()
             if pd.notna(v) and np.isfinite(float(v))}
    if len(clean) < 3:
        return None

    sp_list = sorted(clean.keys())
    vcv_sub, present = _subset_vcv(vcv_df, sp_list)
    if len(present) < 3:
        return None

    y = np.array([clean[s] for s in present], dtype=float)
    X = np.ones((len(present), 1))
    V = vcv_sub.values.astype(float)

    lam_hat, ll_max = _estimate_lambda(X, y, V)
    fit = _gls_at_lambda(X, y, V, lam_hat)
    if fit is None:
        return None

    return dict(
        intercept    = float(fit["beta"][0]),
        se           = float(fit["se"][0]),
        t            = float(fit["t"][0]),
        p            = float(fit["p"][0]),
        lambda_pagel = lam_hat,
        ll           = ll_max,
        sigma2       = fit["sigma2"],
        n_species    = len(present),
        df_residual  = fit["df"],
        species      = present,
    )


# ─────────────────────────────────────────────────────────────────────────────
# PART 4 — PGLS regression  (y ~ x, with optional log transforms)
# ─────────────────────────────────────────────────────────────────────────────

def pgls_regression(
    x_vals,
    y_vals,
    vcv_df: pd.DataFrame,
    species: list,
    log_x: bool = False,
    log_y: bool = False,
    n_pred: int = 200,
) -> dict | None:
    """
    PGLS simple regression with ML estimation of Pagel's λ.

      y ~ β₀ + β₁·x   (under Brownian-motion covariance scaled by λ)

    Parameters
    ----------
    x_vals, y_vals : arrays aligned with `species`
    vcv_df         : phylogenetic VCV DataFrame
    species        : list of species names aligned with x/y arrays
    log_x, log_y   : apply log₁₀ before fitting (standard allometric scaling)
    n_pred         : number of points in the fitted curve for plotting

    Returns
    -------
    dict:
        slope, intercept, r_pgls, p, se_slope, se_intercept,
        lambda_pagel, ll, ll_ols,      ← ll_ols = log-lik at λ=0 (OLS)
        n_species, df_residual,
        x_fit, y_fit, ci_lo, ci_hi,   ← in original (back-transformed) scale
        species_used, log_x, log_y
    or None on failure / degeneracy.

    Notes on r_pgls
    ---------------
    We report the square root of the phylogenetic R² = 1 − RSS_PGLS/TSS_PGLS
    where both RSS and TSS are computed in the GLS metric (Vi-weighted), which
    is the correct phylogenetic analogue of the Pearson r.
    """
    x_arr = np.asarray(x_vals, dtype=float)
    y_arr = np.asarray(y_vals, dtype=float)
    sp_arr = np.asarray(species)

    # ── mask invalid values ──────────────────────────────────────────────────
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if log_x:
        mask &= (x_arr > 0)
    if log_y:
        mask &= (y_arr > 0)

    sp_ok = sp_arr[mask].tolist()
    x_ok  = x_arr[mask]
    y_ok  = y_arr[mask]

    if len(sp_ok) < 3:
        return None

    vcv_sub, present = _subset_vcv(vcv_df, sp_ok)
    if len(present) < 3:
        return None

    idx = [sp_ok.index(s) for s in present]
    xv  = x_ok[idx]
    yv  = y_ok[idx]

    lx = np.log10(xv) if log_x else xv
    ly = np.log10(yv) if log_y else yv

    if np.std(lx) < 1e-12:
        warnings.warn("  [PGLS] degenerate x (zero variance) — skipped")
        return None

    V = vcv_sub.values.astype(float)
    X = np.column_stack([np.ones(len(lx)), lx])

    # ── estimate λ by ML ─────────────────────────────────────────────────────
    lam_hat, ll_max = _estimate_lambda(X, ly, V)

    # ── fit GLS at optimal λ ─────────────────────────────────────────────────
    fit = _gls_at_lambda(X, ly, V, lam_hat)
    if fit is None:
        return None

    # ── OLS log-likelihood (λ=0) for comparison ─────────────────────────────
    fit_ols = _gls_at_lambda(X, ly, V, 0.0)
    ll_ols  = fit_ols["ll"] if fit_ols is not None else np.nan
    fit_bm  = _gls_at_lambda(X, ly, V, 1.0)
    ll_bm   = fit_bm["ll"] if fit_bm is not None else np.nan
    # likelihood-ratio tests of lambda against its bounds (chi2, 1 df; caper convention).
    # With ~10 species lambda is weakly identified: report these, not lambda alone.
    p_lambda0 = float(_stats.chi2.sf(max(2 * (ll_max - ll_ols), 0.0), 1)) if np.isfinite(ll_ols) else np.nan
    p_lambda1 = float(_stats.chi2.sf(max(2 * (ll_max - ll_bm), 0.0), 1)) if np.isfinite(ll_bm) else np.nan

    intercept = float(fit["beta"][0])
    slope     = float(fit["beta"][1])
    se_int    = float(fit["se"][0])
    se_sl     = float(fit["se"][1])
    p_slope   = float(fit["p"][1])
    df        = fit["df"]
    Vi        = fit["Vi"]

    # ── phylogenetic R²  (GLS metric) ────────────────────────────────────────
    ly_fit = X @ fit["beta"]
    ly_mean_gls = float(np.ones(len(ly)) @ Vi @ ly / (np.ones(len(ly)) @ Vi @ np.ones(len(ly))))
    tss_gls = float((ly - ly_mean_gls) @ Vi @ (ly - ly_mean_gls))
    rss_gls = float(fit["resid"] @ Vi @ fit["resid"])
    r2_pgls = max(0.0, 1.0 - rss_gls / tss_gls) if tss_gls > 1e-15 else 0.0
    r_pgls  = float(np.sign(slope) * np.sqrt(r2_pgls))

    # ── prediction curve + 95 % CI (in log space, then back-transformed) ─────
    x_fit_l = np.linspace(lx.min(), lx.max(), n_pred)
    y_fit_l = intercept + slope * x_fit_l

    t_crit   = _stats.t.ppf(0.975, df=df)
    n_v      = len(lx)
    x_mean_l = float(lx.mean())
    ssx      = max(float(np.sum((lx - x_mean_l) ** 2)), 1e-12)

    # 95 % confidence band of the GLS regression line: se = sqrt(x' Cov(beta) x).
    # (The previous band was an OLS *prediction* interval, 1 + 1/n + (x-x̄)²/Sxx,
    # which ignores the phylogenetic covariance.)
    Xg = np.column_stack([np.ones_like(x_fit_l), x_fit_l])
    var_mean = np.clip(np.einsum('ij,jk,ik->i', Xg, fit["cov_beta"], Xg), 0, None)
    se_pred = np.sqrt(var_mean)
    # 95 % prediction interval for a new species of average phylogenetic
    # position: Var = x' Cov(beta) x + sigma2 * mean(diag(V_lambda)).
    tip_var = float(np.mean(np.diag(_make_vpagel(V, lam_hat))))
    se_pi = np.sqrt(var_mean + fit["sigma2"] * tip_var)

    if log_x:
        x_plot = 10.0 ** x_fit_l
    else:
        x_plot = x_fit_l

    if log_y:
        y_plot = 10.0 ** y_fit_l
        ci_lo  = 10.0 ** (y_fit_l - t_crit * se_pred)
        ci_hi  = 10.0 ** (y_fit_l + t_crit * se_pred)
        pi_lo  = 10.0 ** (y_fit_l - t_crit * se_pi)
        pi_hi  = 10.0 ** (y_fit_l + t_crit * se_pi)
    else:
        y_plot = y_fit_l
        ci_lo  = y_fit_l - t_crit * se_pred
        ci_hi  = y_fit_l + t_crit * se_pred
        pi_lo  = y_fit_l - t_crit * se_pi
        pi_hi  = y_fit_l + t_crit * se_pi

    return dict(
        slope        = slope,
        intercept    = intercept,
        r_pgls       = r_pgls,          # phylogenetic R (GLS metric)
        r            = r_pgls,          # alias for back-compat
        p            = p_slope,
        se_slope     = se_sl,
        se_intercept = se_int,
        lambda_pagel = lam_hat,
        ll           = ll_max,
        ll_ols       = ll_ols,
        ll_bm        = ll_bm,
        p_lambda0    = p_lambda0,
        p_lambda1    = p_lambda1,
        n_species    = n_v,
        df_residual  = df,
        x_fit        = x_plot,
        y_fit        = y_plot,
        ci_lo        = ci_lo,
        ci_hi        = ci_hi,
        pi_lo        = pi_lo,
        pi_hi        = pi_hi,
        species_used = present,
        log_x        = log_x,
        log_y        = log_y,
    )


def _pgls_on_ax_with_prediction(ax, res, log_x=False, log_y=False):
    """
    Draw PGLS regression line with 95% PREDICTION interval (not confidence).
    Prediction intervals are wider and more appropriate for individual points.
    """
    import matplotlib.ticker as tck

    # Plot the line
    ax.plot(
        res["x_fit"], res["y_fit"],
        color="k", lw=1.5, ls="--", alpha=0.8, zorder=3
    )

    # Bands computed in pgls_regression: GLS confidence band of the line and
    # prediction interval for a new species. (The previous version widened the
    # confidence band by an arbitrary factor sqrt(1+1/n+1/df) * 1.5.)
    if 'pi_lo' in res:
        ax.fill_between(res["x_fit"], res["pi_lo"], res["pi_hi"],
                        color="k", alpha=0.08, zorder=1, label="_nolegend_")
    ax.fill_between(res["x_fit"], res["ci_lo"], res["ci_hi"],
                    color="k", alpha=0.15, zorder=2, label="_nolegend_")

    # Add annotation with PGLS stats
    p_str = f"p={res['p']:.3f}" if res['p'] >= 0.001 else "p<0.001"
    ax.text(
        0.05, 0.85,
        f"PGLS  slope={res['slope']:.2f}  "
        f"r={res['r_pgls']:.2f}  "
        f"{p_str}  "
        f"λ={res['lambda_pagel']:.2f}",
        transform=ax.transAxes,
        fontsize=7,
        va="top",
        ha="left",
        color="#333",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.85)
    )

    if log_x:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(tck.LogFormatterSciNotation())

    if log_y:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(tck.LogFormatterSciNotation())


# ─────────────────────────────────────────────────────────────────────────────
# PART 5 — Axes-level helpers
# ─────────────────────────────────────────────────────────────────────────────

def ols_on_ax(ax, xv, yv, log_x=False):
    """Plain OLS line + 95 % band.  Kept for back-compatibility."""
    import matplotlib.ticker as tck
    xv = np.asarray(xv, dtype=float)
    yv = np.asarray(yv, dtype=float)
    lx = np.log10(np.clip(xv, 1e-9, None)) if log_x else xv
    ok = np.isfinite(lx) & np.isfinite(yv)
    if ok.sum() < 3:
        return None
    lxv, lyv = lx[ok], yv[ok]
    ssx = float(np.sum((lxv - lxv.mean()) ** 2))
    if ssx < 1e-12:
        ax.axhline(float(np.median(lyv)), color="k", lw=1.2, ls="--", alpha=0.55)
        return None
    slope, intercept, r, p, _ = _stats.linregress(lxv, lyv)
    n_v = ok.sum()
    t_c = _stats.t.ppf(0.975, df=n_v - 2)
    xf  = np.linspace(lxv.min(), lxv.max(), 200)
    yf  = intercept + slope * xf
    res = np.sum((lyv - (intercept + slope * lxv)) ** 2) / (n_v - 2)
    se  = np.sqrt(res * (1 + 1 / n_v + (xf - lxv.mean()) ** 2 / ssx))
    xp  = 10.0 ** xf if log_x else xf
    ax.plot(xp, yf, color="k", lw=1.5, ls="--", alpha=0.8, zorder=3)
    ax.fill_between(xp, yf - t_c * se, yf + t_c * se,
                    color="k", alpha=0.15, zorder=2)
    p_str = f"{p:.3f}" if p >= 0.001 else "<0.001"
    ax.text(0.05, 0.95, f"OLS slope={slope:.2f}  r={r:.2f}  p={p_str}",
            transform=ax.transAxes, fontsize=8, va="top", color="#555")
    return slope, r, p


def pgls_on_ax(ax, xv, yv, vcv_df, species,
               log_x=False, log_y=False):
    """
    PGLS regression line + 95 % prediction band drawn directly on `ax`.
    Annotates with slope, phylogenetic r, p-value, and Pagel's λ.

    Returns (slope, r_pgls, p, lambda_pagel) on success, None on failure.
    """
    import matplotlib.ticker as tck

    res = pgls_regression(xv, yv, vcv_df, species, log_x=log_x, log_y=log_y)
    if res is None:
        return None

    ax.plot(res["x_fit"], res["y_fit"],
            color="k", lw=1.5, ls="--", alpha=0.8, zorder=3)
    ax.fill_between(res["x_fit"], res["ci_lo"], res["ci_hi"],
                    color="k", alpha=0.15, zorder=2)

    p_str = f"{res['p']:.3f}" if res["p"] >= 0.001 else "<0.001"
    lam   = res["lambda_pagel"]
    ax.text(0.05, 0.95,
            f"PGLS  slope={res['slope']:.2f}  r={res['r_pgls']:.2f}"
            f"  p={p_str}  λ={lam:.2f}",
            transform=ax.transAxes, fontsize=8, va="top", color="#333")

    if log_x:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(tck.LogFormatterSciNotation())
    if log_y:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(tck.LogFormatterSciNotation())

    return res["slope"], res["r_pgls"], res["p"], lam


# ─────────────────────────────────────────────────────────────────────────────
# PART 6 — All-pairs PGLS across cortical regions
# ─────────────────────────────────────────────────────────────────────────────

def run_pgls_all_region_pairs(
    df_morph: pd.DataFrame,
    metric_col: str,
    vcv_df: pd.DataFrame,
    regions: list = None,
    atlas_level: int = 1,
    hemisphere: str = "bilateral",
    min_species: int = 4,
    log_x: bool = True,
    log_y: bool = True,
) -> pd.DataFrame:
    """
    For every ordered pair (region_x, region_y), fit:

        log₁₀(mean_y_per_species)  ~  log₁₀(mean_x_per_species)

    using PGLS with ML-estimated Pagel λ.

    Returns
    -------
    DataFrame with columns:
        region_x, region_y, slope, intercept, r_pgls, p, p_fdr,
        se_slope, lambda_pagel, ll, ll_ols,
        n_species, df_residual, significant_fdr05, mean_x, mean_y
    Sorted by p (ascending).
    """
    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == atlas_level]
    if hemisphere and hemisphere != "all" and "hemisphere" in df.columns:
        h = df[df["hemisphere"] == hemisphere]
        df = h if not h.empty else df[df["hemisphere"] == "bilateral"]

    if regions is None:
        regions = sorted(df["region"].dropna().unique().tolist())

    # Per-species means per region
    region_means: dict = {}
    for reg in regions:
        sub = df[df["region"] == reg]
        if sub.empty:
            continue
        means = sub.groupby("species")[metric_col].mean()
        if isinstance(means.index, pd.MultiIndex):
            region_means[reg] = {k[0]: v for k, v in means.items()}
        else:
            region_means[reg] = means.to_dict()

    avail = [r for r in regions
             if r in region_means and len(region_means[r]) >= min_species]
    print(f"  [PGLS pairs] regions ≥{min_species} species: {len(avail)}/{len(regions)}")

    results = []
    pairs = list(itertools.combinations(avail, 2))
    print(f"  [PGLS pairs] testing {len(pairs)} pairs ...")

    for reg_x, reg_y in pairs:
        mx, my = region_means[reg_x], region_means[reg_y]
        common = [s for s in vcv_df.index
                  if s in mx and s in my
                  and pd.notna(mx[s]) and pd.notna(my[s])]
        if len(common) < min_species:
            continue

        xv = np.array([mx[s] for s in common])
        yv = np.array([my[s] for s in common])

        res = pgls_regression(xv, yv, vcv_df, common,
                              log_x=log_x, log_y=log_y)
        if res is None:
            continue

        results.append(dict(
            region_x     = reg_x,
            region_y     = reg_y,
            slope        = res["slope"],
            intercept    = res["intercept"],
            r_pgls       = res["r_pgls"],
            r            = res["r_pgls"],   # back-compat alias
            p            = res["p"],
            se_slope     = res["se_slope"],
            lambda_pagel = res["lambda_pagel"],
            ll           = res["ll"],
            ll_ols       = res["ll_ols"],
            n_species    = res["n_species"],
            df_residual  = res["df_residual"],
            mean_x       = float(np.mean(xv)),
            mean_y       = float(np.mean(yv)),
            species_used = ", ".join(res["species_used"]),
        ))

    if not results:
        warnings.warn("  [PGLS pairs] no results returned.")
        return pd.DataFrame()

    df_res = pd.DataFrame(results).sort_values("p").reset_index(drop=True)

    # Benjamini-Hochberg FDR over the finite p-values (a single NaN previously
    # propagated through the cumulative minimum and turned every p_fdr into NaN)
    p_vals = df_res["p"].to_numpy(dtype=float)
    ok = np.isfinite(p_vals)
    n_tests = int(ok.sum())
    p_fdr = np.full(len(p_vals), np.nan)
    if n_tests:
        pv = p_vals[ok]                              # already sorted ascending
        ranks = np.arange(1, n_tests + 1)
        q = np.minimum.accumulate((pv * n_tests / ranks)[::-1])[::-1]
        p_fdr[ok] = np.minimum(q, 1.0)
    df_res["p_fdr"] = p_fdr
    df_res["significant_fdr05"] = df_res["p_fdr"] < 0.05

    n_sig = df_res["significant_fdr05"].sum()
    print(f"  [PGLS pairs] {n_tests} pairs | {n_sig} FDR<0.05")
    return df_res


# ─────────────────────────────────────────────────────────────────────────────
# PART 7 — Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────

def _species_colors(species_list):
    u = _phylo_sort(list(set(species_list)))
    return {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(u)}


def plot_pgls_region_pair(
    result_row: pd.Series,
    df_morph: pd.DataFrame,
    metric_col: str,
    metric_label: str,
    vcv_df: pd.DataFrame,
    output_path: str,
    atlas_level: int = 1,
    hemisphere: str = "bilateral",
    log_x: bool = True,
    log_y: bool = True,
):
    """
    Scatter + PGLS line for one (region_x, region_y) pair.
    Each dot = one species (per-species mean).  λ is shown in the annotation.
    """
    import matplotlib.ticker as tck

    reg_x = result_row["region_x"]
    reg_y = result_row["region_y"]

    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == atlas_level]
    if hemisphere != "all" and "hemisphere" in df.columns:
        h = df[df["hemisphere"] == hemisphere]
        df = h if not h.empty else df[df["hemisphere"] == "bilateral"]

    mx = df[df["region"] == reg_x].groupby("species")[metric_col].mean()
    my = df[df["region"] == reg_y].groupby("species")[metric_col].mean()
    common = sorted(set(mx.index) & set(my.index) & set(vcv_df.index))
    if len(common) < 3:
        return None

    xv = mx.loc[common].values
    yv = my.loc[common].values
    sp_colors = _species_colors(common)

    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(5.5, 5))
        for sp, x, y in zip(common, xv, yv):
            ax.scatter(x, y, color=sp_colors[sp], s=55, zorder=5,
                       edgecolors="k", linewidths=0.5, label=sp)

        pgls_on_ax(ax, xv, yv, vcv_df, common, log_x=log_x, log_y=log_y)

        if log_x:
            ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")

        ax.set_xlabel(f"{reg_x}  ({metric_label})", fontsize=9)
        ax.set_ylabel(f"{reg_y}  ({metric_label})", fontsize=9)

        p_fdr = float(result_row.get("p_fdr", np.nan))
        sig = ("***" if p_fdr < 0.001 else "**" if p_fdr < 0.01
               else "*" if p_fdr < 0.05 else "n.s.")
        lam = float(result_row.get("lambda_pagel", np.nan))
        ax.set_title(f"{reg_x}  ×  {reg_y}\n"
                     f"FDR={p_fdr:.3f} {sig}   λ={lam:.2f}",
                     fontweight="bold", fontsize=9)

        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                   for sp in _phylo_sort(common)]
        ax.legend(handles=handles, fontsize=7,
                  bbox_to_anchor=(1.01, 1), loc="upper left")
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def plot_pgls_summary_heatmap(
    df_results: pd.DataFrame,
    regions: list,
    output_path: str,
    value_col: str = "slope",
    title: str = "PGLS allometric slopes",
):
    """Symmetric heatmap of `value_col`.  Cells with FDR<0.05 marked *."""
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.cm as _cm

    if df_results.empty:
        return None

    n = len(regions)
    mat = np.full((n, n), np.nan)
    mat_sig = np.zeros((n, n), bool)
    idx = {r: i for i, r in enumerate(regions)}

    for _, row in df_results.iterrows():
        ri = idx.get(row["region_x"])
        ci = idx.get(row["region_y"])
        if ri is None or ci is None:
            continue
        v = float(row[value_col]) if pd.notna(row[value_col]) else np.nan
        mat[ri, ci] = mat[ci, ri] = v
        sig = bool(row.get("significant_fdr05", False))
        mat_sig[ri, ci] = mat_sig[ci, ri] = sig

    finite = mat[np.isfinite(mat)]
    if finite.size == 0:
        return None

    vext = float(np.nanpercentile(np.abs(finite), 95))
    norm = TwoSlopeNorm(vmin=-vext, vcenter=0, vmax=vext) if vext > 0 else None

    short = [r[:20] for r in regions]
    w = max(5, n * 0.9 + 2.5)
    h = max(4, n * 0.65 + 1.5)

    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(w, h))
        im = ax.imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto",
                       interpolation="nearest")
        for ri in range(n):
            for ci in range(n):
                v = mat[ri, ci]
                if not np.isfinite(v):
                    continue
                txt = f"{v:.2f}" + (" *" if mat_sig[ri, ci] else "")
                tc  = "white" if (vext > 0 and abs(v) > 0.4 * vext) else "#222"
                ax.text(ci, ri, txt, ha="center", va="center",
                        fontsize=6, color=tc)
        ax.set_xticks(range(n))
        ax.set_xticklabels(short, rotation=40, ha="right", fontsize=8)
        ax.set_yticks(range(n)); ax.set_yticklabels(short, fontsize=8)
        ax.set_title(f"{title}  (FDR<0.05 marked *)",
                     fontweight="bold", fontsize=10)
        cb = fig.colorbar(_cm.ScalarMappable(norm=norm, cmap="RdBu_r"),
                          ax=ax, fraction=0.03, pad=0.04)
        cb.set_label(value_col, fontsize=8)
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


def plot_pgls_dotplot(
    df_results: pd.DataFrame,
    output_path: str,
    value_col: str = "slope",
    top_n: int = 40,
    title: str = "PGLS allometric slopes — top pairs",
):
    """Dot-plot of top PGLS pairs; dot size ∝ −log₁₀(p_fdr); colour = significance."""
    if df_results.empty:
        return None

    df_plot = df_results.sort_values("p").head(top_n).iloc[::-1].reset_index(drop=True)
    labels = [f"{r['region_x'][:16]} ↔ {r['region_y'][:16]}"
              for _, r in df_plot.iterrows()]
    colors = ["#D55E00" if s else "#999" for s in df_plot["significant_fdr05"]]
    sizes  = np.clip(-np.log10(df_plot["p_fdr"].clip(lower=1e-10)) * 30, 20, 300)

    h = max(4, len(df_plot) * 0.32 + 1.5)
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(7.5, h))
        y_pos = np.arange(len(df_plot))
        ax.errorbar(df_plot[value_col], y_pos, xerr=df_plot["se_slope"],
                    fmt="none", color="#666", lw=0.8, zorder=2)
        ax.scatter(df_plot[value_col], y_pos, c=colors, s=sizes,
                   zorder=3, edgecolors="k", linewidths=0.4, alpha=0.9)
        ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.5)
        ax.axvline(1, color="#aaa", lw=0.6, ls=":", alpha=0.4)
        ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel(f"PGLS {value_col} ± SE", fontsize=9)
        ax.set_title(title, fontweight="bold", fontsize=10)
        ax.legend(handles=[
            mpatches.Patch(facecolor="#D55E00", label="FDR<0.05"),
            mpatches.Patch(facecolor="#999",    label="n.s."),
        ], fontsize=8, loc="lower right")

        for i, row in df_plot.iterrows():
            p_str = f"p={row['p_fdr']:.3f}" if row["p_fdr"] >= 0.001 else "p<0.001"
            lam   = row.get("lambda_pagel", np.nan)
            note  = f"{p_str}  λ={lam:.2f}" if np.isfinite(float(lam)) else p_str
            ax.text(ax.get_xlim()[1] * 0.98, i, note,
                    va="center", ha="right", fontsize=6, color="#444")

        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# PART 8 — Brain scaling (allometric body/brain weight regressions)
# ─────────────────────────────────────────────────────────────────────────────

def plot_brain_scaling_pgls(
    df_morph: pd.DataFrame,
    metric_col: str,
    metric_label: str,
    output_path: str,
    vcv_df: pd.DataFrame,
    y_region: str,
    x_region: str = None,
    atlas_level: int = 1,
    hemisphere: str = "bilateral",
    species_body_weight_g: dict = None,
    species_brain_weight_g: dict = None,
    figsize=None,
):
    """
    3-panel allometric figure (body weight / brain weight / internal region)
    with PGLS lines and Pagel λ annotations.
    """
    import matplotlib.ticker as tck

    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == atlas_level]
    if hemisphere != "all" and "hemisphere" in df.columns:
        h = df[df["hemisphere"] == hemisphere]
        df = h if not h.empty else df[df["hemisphere"] == "bilateral"]

    df_y = df[df["region"] == y_region].groupby("species")[metric_col].mean()
    if df_y.empty:
        warnings.warn(f"plot_brain_scaling_pgls: no data for y_region={y_region}")
        return None

    sp_list   = list(df_y.index)
    y_vals    = df_y.values
    sp_colors = _species_colors(sp_list)

    panels = []
    if species_body_weight_g:
        bw = np.array([species_body_weight_g.get(s, np.nan) for s in sp_list])
        panels.append((bw, "Body weight (g)", True, False))
    if species_brain_weight_g:
        brw = np.array([species_brain_weight_g.get(s, np.nan) for s in sp_list])
        panels.append((brw, "Brain weight (g)", True, False))
    if x_region:
        df_x = df[df["region"] == x_region].groupby("species")[metric_col].mean()
        xrv  = np.array([float(df_x.loc[s]) if s in df_x.index else np.nan
                         for s in sp_list])
        panels.append((xrv, f"{x_region[:24]} ({metric_label})", True, True))

    panels = [(xv, xl, lx, ly) for xv, xl, lx, ly in panels
              if np.isfinite(xv).sum() >= 3]
    if not panels:
        warnings.warn(f"plot_brain_scaling_pgls: no panels with data for {y_region}")
        return None

    n_panels = len(panels)
    w, h = figsize or (n_panels * 4.5, 4.5)

    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, n_panels, figsize=(w, h))
        if n_panels == 1:
            axes = [axes]

        for pi, (ax, (xv, xlabel, log_x, log_y)) in enumerate(zip(axes, panels)):
            for sp, x, y in zip(sp_list, xv, y_vals):
                if not (np.isfinite(x) and np.isfinite(y)):
                    continue
                ax.scatter(x, y, color=sp_colors[sp], s=40, zorder=5,
                           alpha=0.85, edgecolors="k", linewidths=0.5, label=sp)
            pgls_on_ax(ax, xv, y_vals, vcv_df, sp_list,
                       log_x=log_x, log_y=log_y)
            ax.set_xlabel(xlabel, fontsize=9)
            if pi == 0:
                ax.set_ylabel(f"{y_region} ({metric_label})", fontsize=9)
            ax.text(-0.12, 1.02, chr(65 + pi), transform=ax.transAxes,
                    fontsize=13, fontweight="bold", va="bottom", ha="left")

        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                   for sp in _phylo_sort(sp_list)]
        plt.tight_layout()
        fig.legend(handles=handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.06),
                   ncol=min(len(handles), 7), fontsize=8,
                   frameon=False, title="Species")
        fig.subplots_adjust(bottom=0.18)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# PART 9 — FC connections PGLS  (unchanged logic, λ now estimated)
# ─────────────────────────────────────────────────────────────────────────────

def run_pgls_fc_connections(
    df_fc: pd.DataFrame,
    fc_columns: list,
    vcv_df: pd.DataFrame,
    x_columns: list = None,
    x_labels: list = None,
    x_log: list = None,
    min_species: int = 4,
) -> pd.DataFrame:
    """PGLS for each FC column vs each x column. λ estimated by ML."""
    x_columns = x_columns or []
    x_labels  = x_labels  or x_columns
    x_log     = x_log     or [False] * len(x_columns)

    results = []
    for fc_col in fc_columns:
        if fc_col not in df_fc.columns:
            continue
        sp_means_y = df_fc.groupby("species")[fc_col].mean().dropna()
        for x_col, x_lbl, lx in zip(x_columns, x_labels, x_log):
            if x_col not in df_fc.columns:
                continue
            sp_means_x = df_fc.groupby("species")[x_col].mean().dropna()
            common = sorted(set(sp_means_x.index) & set(sp_means_y.index)
                            & set(vcv_df.index))
            if len(common) < min_species:
                continue
            xv = sp_means_x.loc[common].values
            yv = sp_means_y.loc[common].values
            res = pgls_regression(xv, yv, vcv_df, common, log_x=lx, log_y=False)
            if res is None:
                continue
            results.append(dict(
                fc_col       = fc_col,
                x_col        = x_col,
                x_label      = x_lbl,
                slope        = res["slope"],
                r_pgls       = res["r_pgls"],
                r            = res["r_pgls"],
                p            = res["p"],
                se_slope     = res["se_slope"],
                lambda_pagel = res["lambda_pagel"],
                n_species    = res["n_species"],
                species_used = ", ".join(res["species_used"]),
            ))

    if not results:
        return pd.DataFrame()

    df_r = pd.DataFrame(results).sort_values("p").reset_index(drop=True)
    n = len(df_r)
    ranks = np.arange(1, n + 1)
    p_fdr = np.minimum.accumulate(
        df_r["p"].values[::-1] * n / ranks[::-1])[::-1]
    df_r["p_fdr"]             = np.minimum(p_fdr, 1.0)
    df_r["significant_fdr05"] = df_r["p_fdr"] < 0.05
    return df_r


def plot_fc_connections_pgls(
    df_fc: pd.DataFrame,
    connections: list,
    vcv_df: pd.DataFrame,
    output_path_template: str,
    species_body_weight_g: dict = None,
    species_brain_weight_g: dict = None,
    figsize=None,
):
    """One PGLS figure per FC connection (body / brain / internal reference)."""
    import matplotlib.ticker as tck

    df_fc = df_fc.copy()
    if species_body_weight_g:
        df_fc["body_g"]  = df_fc["species"].map(species_body_weight_g)
    if species_brain_weight_g:
        df_fc["brain_g"] = df_fc["species"].map(species_brain_weight_g)

    sp_colors = _species_colors(df_fc["species"].unique())
    out_paths = []

    for conn in connections:
        col_key = conn[0]
        y_label = conn[-1]
        if col_key not in df_fc.columns or df_fc[col_key].notna().sum() == 0:
            continue

        x_specs = []
        if "body_g"  in df_fc.columns:
            x_specs.append(("body_g",  "Body weight (g)",  True))
        if "brain_g" in df_fc.columns:
            x_specs.append(("brain_g", "Brain weight (g)", True))
        if (col_key != "FC_somato_lr"
                and "FC_somato_lr" in df_fc.columns
                and df_fc["FC_somato_lr"].notna().any()):
            x_specs.append(("FC_somato_lr", "Somato L-R (r/mean|r|)", False))

        x_specs = [(xc, xl, lx) for xc, xl, lx in x_specs
                   if xc in df_fc.columns and df_fc[xc].notna().any()]
        if not x_specs:
            continue

        n_panels = len(x_specs)
        w, h = figsize or (n_panels * 4.5, 4.5)

        with plt.rc_context(PAPER_RC):
            fig, axes = plt.subplots(1, n_panels, figsize=(w, h))
            if n_panels == 1:
                axes = [axes]

            for pi, (ax, (x_col, x_label, log_x)) in enumerate(zip(axes, x_specs)):
                sub = df_fc[[col_key, x_col, "species"]].dropna()
                if sub.empty:
                    ax.set_visible(False); continue
                for sp in _phylo_sort(sub["species"].unique()):
                    ss = sub[sub["species"] == sp]
                    if ss.empty: continue
                    ax.scatter(ss[x_col], ss[col_key],
                               color=sp_colors.get(sp, "#888"),
                               s=40, zorder=5, alpha=0.85,
                               edgecolors="k", linewidths=0.5)
                pgls_on_ax(ax, sub[x_col].values.astype(float),
                           sub[col_key].values.astype(float),
                           vcv_df, sub["species"].tolist(),
                           log_x=log_x, log_y=False)
                ax.axhline(0, color="#aaa", lw=0.8, ls=":", zorder=1)
                ax.set_xlabel(x_label, fontsize=9)
                if pi == 0:
                    ax.set_ylabel(y_label, fontsize=9)
                ax.text(-0.12, 1.02, chr(65 + pi), transform=ax.transAxes,
                        fontsize=13, fontweight="bold", va="bottom")

            handles = [mpatches.Patch(facecolor=sp_colors.get(sp, "#888"), label=sp)
                       for sp in _phylo_sort(df_fc["species"].unique())]
            plt.tight_layout()
            fig.legend(handles=handles, loc="lower center",
                       bbox_to_anchor=(0.5, -0.06), ncol=min(len(handles), 7),
                       fontsize=8, frameon=False, title="Species")
            fig.subplots_adjust(bottom=0.18)
            safe = col_key.replace("FC_", "")
            out  = output_path_template.replace(".png", f"_{safe}.png")
            os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
            fig.savefig(out, bbox_inches="tight", dpi=200)
            plt.close(fig)
            print(f"  [plot] {out}")
            out_paths.append(out)

    return out_paths
