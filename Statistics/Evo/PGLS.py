"""
ednix_pgls_functions.py
=======================
Phylogenetic Generalised Least Squares (PGLS) functions for EDNiX.

Replaces all plain OLS regression (anatomical allometry, FC scaling) with
PGLS using a Brownian-motion covariance structure derived from a Newick tree.
No R / rpy2 required: implemented via NumPy and SciPy GLS.

Public API
----------
load_newick_vcv(newick_path, species_map=None)
    ? (vcv_df, tree_species_list)

pgls_intercept(species_vals, vcv_df)
    ? dict (intercept, se, t, p, lambda_pagel, n, df_residual)

pgls_regression(x_vals, y_vals, vcv_df, log_x=False, log_y=False)
    ? dict (slope, intercept, r, p, se_slope, lambda_pagel, n, x_fit, y_fit, ci_lo, ci_hi)

run_pgls_all_region_pairs(df_morph, metric_col, vcv_df,
                           regions=None, atlas_level=1,
                           hemisphere='bilateral', min_species=4)
    ? df_results (one row per region pair, sorted by p_fdr)

plot_pgls_region_pair(result_row, df_morph, metric_col, metric_label,
                       vcv_df, output_path, atlas_level=1)
    ? output_path

plot_pgls_summary_heatmap(df_results, regions, output_path, value_col='slope')
    ? output_path

plot_pgls_dotplot(df_results, output_path)
    ? output_path

ols_on_ax(ax, xv, yv, log_x=False)           [unchanged helper, kept for back-compat]
pgls_on_ax(ax, xv, yv, vcv_df, log_x=False)  [PGLS-aware drop-in replacement]
"""

import os
import re
import warnings
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats as _stats

# ?????????????????????????????????????????????????????????????????????????????
# PAPER STYLE  (mirrors ednix_bids_tools)
# ?????????????????????????????????????????????????????????????????????????????

PAPER_RC = {
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     1.2,
    "xtick.major.width":  1.2,
    "ytick.major.width":  1.2,
    "font.family":        "sans-serif",
    "font.size":          11,
    "axes.titlesize":     12,
    "axes.labelsize":     11,
    "legend.frameon":     False,
    "figure.dpi":         150,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
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

# Mapping: EDNiX species name ? Newick tip label substring
# (case-insensitive partial match used if exact not found)
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


# ?????????????????????????????????????????????????????????????????????????????
# PART 1 ? NEWICK PARSER & VCV BUILDER
# ?????????????????????????????????????????????????????????????????????????????

class _Node:
    """Minimal binary tree node for Newick parsing."""
    def __init__(self, name="", length=0.0):
        self.name     = name
        self.length   = float(length) if length else 0.0
        self.children = []

    def is_leaf(self):
        return len(self.children) == 0


def _parse_newick(s: str) -> _Node:
    """Recursive-descent Newick parser. Handles branch lengths and labels."""
    s = s.strip().rstrip(";")

    def _parse(tok: str, pos: int):
        node = _Node()
        if pos < len(tok) and tok[pos] == "(":
            pos += 1                          # consume '('
            while True:
                child, pos = _parse(tok, pos)
                node.children.append(child)
                if pos >= len(tok):
                    break
                if tok[pos] == ",":
                    pos += 1
                    continue
                if tok[pos] == ")":
                    pos += 1
                    break
        # read optional label
        label_start = pos
        while pos < len(tok) and tok[pos] not in (",", ")", ":", ";"):
            pos += 1
        node.name = tok[label_start:pos].strip().strip("'\"")
        # read optional branch length
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
    """Return list of (leaf_name, root_to_tip_distance)."""
    depth += node.length
    if node.is_leaf():
        return [(node.name, depth)]
    out = []
    for child in node.children:
        out.extend(_collect_leaves(child, depth))
    return out


def _pairwise_shared(node: _Node, depth: float = 0.0):
    """
    Return dict {frozenset({leaf_a, leaf_b}): shared_branch_length}.
    Shared branch length = distance from root to the MRCA.
    Under Brownian motion: Cov(i,j) = t_mrca(i,j).
    """
    depth += node.length
    if node.is_leaf():
        return {}, {node.name: depth}

    child_results = [_pairwise_shared(c, depth) for c in node.children]
    # Merge all pairwise dicts and leaf-depth dicts
    all_pairs  = {}
    all_leaves = {}
    for pairs, leaves in child_results:
        all_pairs.update(pairs)
        all_leaves.update(leaves)

    # Cross-child pairs: their MRCA is the current node at `depth`
    for (pairs_a, leaves_a), (pairs_b, leaves_b) in itertools.combinations(
        [(p, l) for p, l in child_results], 2
    ):
        for la in leaves_a:
            for lb in leaves_b:
                all_pairs[frozenset([la, lb])] = depth

    return all_pairs, all_leaves


def load_newick_vcv(newick_path: str, species_map: dict = None) -> tuple:
    """
    Parse a Newick file and build the phylogenetic VCV matrix (Brownian motion).

    Parameters
    ----------
    newick_path  : path to .nwk file
    species_map  : dict mapping EDNiX short names ? Newick tip-label substrings.
                   Defaults to DEFAULT_SPECIES_MAP.

    Returns
    -------
    vcv_df       : pd.DataFrame  (species × species), Cov = shared branch length
    tip_names    : list of Newick tip labels found in the tree
    """
    species_map = species_map or DEFAULT_SPECIES_MAP

    with open(newick_path, encoding="utf-8") as fh:
        nwk_str = fh.read().strip()

    root = _parse_newick(nwk_str)
    pairs, leaf_depths = _pairwise_shared(root)
    tip_names = list(leaf_depths.keys())

    # Build reverse map: newick_tip ? EDNiX_short_name
    tip_to_short = {}
    for short, pattern in species_map.items():
        for tip in tip_names:
            if tip.lower() == pattern.lower() or pattern.lower() in tip.lower():
                tip_to_short[tip] = short
                break

    # Build EDNiX-named list of species present in tree
    ednix_names = [tip_to_short.get(t, t) for t in tip_names]
    n = len(ednix_names)

    vcv = np.zeros((n, n))
    for i, ti in enumerate(tip_names):
        vcv[i, i] = leaf_depths[ti]           # variance = root-to-tip
        for j, tj in enumerate(tip_names):
            if i == j:
                continue
            shared = pairs.get(frozenset([ti, tj]), 0.0)
            vcv[i, j] = shared
            vcv[j, i] = shared

    vcv_df = pd.DataFrame(vcv, index=ednix_names, columns=ednix_names)
    print(f"  [VCV] Loaded {n} tips from {Path(newick_path).name}")
    print(f"  [VCV] Species in VCV: {ednix_names}")
    return vcv_df, ednix_names


# ?????????????????????????????????????????????????????????????????????????????
# PART 2 ? CORE GLS ENGINE
# ?????????????????????????????????????????????????????????????????????????????

def _gls_core(X: np.ndarray, y: np.ndarray, V: np.ndarray):
    """
    Pure-numpy GLS.  Solves  ? = (X'V?¹X)?¹ X'V?¹y.

    Parameters
    ----------
    X  : (n, k) design matrix (already includes intercept column)
    y  : (n,)   response
    V  : (n, n) covariance matrix (assumed SPD)

    Returns
    -------
    beta       : (k,) coefficient vector
    se         : (k,) standard errors
    residuals  : (n,) y - X @ beta
    sigma2     : float, residual variance
    t_stats    : (k,) t-statistics
    p_vals     : (k,) two-tailed p-values
    df_resid   : int
    """
    n, k = X.shape
    try:
        V_inv = np.linalg.inv(V)
    except np.linalg.LinAlgError:
        V_inv = np.linalg.pinv(V)

    XtVi   = X.T @ V_inv
    XtViX  = XtVi @ X
    try:
        XtViX_inv = np.linalg.inv(XtViX)
    except np.linalg.LinAlgError:
        XtViX_inv = np.linalg.pinv(XtViX)

    beta      = XtViX_inv @ (XtVi @ y)
    residuals = y - X @ beta
    df_resid  = n - k

    if df_resid <= 0:
        raise ValueError(f"Too few observations (n={n}) for k={k} parameters.")

    sigma2 = float(residuals @ V_inv @ residuals) / df_resid
    cov_b  = sigma2 * XtViX_inv
    se     = np.sqrt(np.diag(cov_b))

    t_stats = beta / np.where(se > 0, se, np.nan)
    p_vals  = 2 * _stats.t.sf(np.abs(t_stats), df=df_resid)

    return beta, se, residuals, sigma2, t_stats, p_vals, df_resid


def _subset_vcv(vcv_df: pd.DataFrame, species: list):
    """Return the sub-matrix of vcv_df for the given species list."""
    present = [s for s in species if s in vcv_df.index]
    missing = [s for s in species if s not in vcv_df.index]
    if missing:
        warnings.warn(f"  [PGLS] Species not in VCV: {missing} ? excluded.")
    return vcv_df.loc[present, present], present


def _phylo_sort(species_list):
    known   = [s for s in PHYLO_ORDER if s in species_list]
    unknown = sorted([s for s in species_list if s not in PHYLO_ORDER])
    return known + unknown


# ?????????????????????????????????????????????????????????????????????????????
# PART 3 ? PGLS INTERCEPT (test whether mean ? 0)
# ?????????????????????????????????????????????????????????????????????????????

def pgls_intercept(species_vals: dict, vcv_df: pd.DataFrame) -> dict:
    """
    PGLS intercept-only model:  y_species ~ 1

    Tests whether the phylogenetically corrected mean is significantly
    different from zero.

    Parameters
    ----------
    species_vals : {species: value}   (NaN values are dropped)
    vcv_df       : phylogenetic VCV DataFrame

    Returns
    -------
    dict with: intercept, se, t, p, n_species, df_residual
               or None when insufficient data.
    """
    clean = {sp: v for sp, v in species_vals.items()
             if pd.notna(v) and np.isfinite(v)}
    if len(clean) < 3:
        return None

    sp_list = sorted(clean.keys())
    vcv_sub, present = _subset_vcv(vcv_df, sp_list)
    y = np.array([clean[s] for s in present], dtype=float)
    X = np.ones((len(present), 1))

    try:
        beta, se, _, _, t_stats, p_vals, df_r = _gls_core(X, y, vcv_sub.values)
    except Exception as exc:
        warnings.warn(f"  [PGLS intercept] GLS failed: {exc}")
        return None

    return dict(
        intercept   = float(beta[0]),
        se          = float(se[0]),
        t           = float(t_stats[0]),
        p           = float(p_vals[0]),
        n_species   = len(present),
        df_residual = df_r,
        species     = present,
    )


# ?????????????????????????????????????????????????????????????????????????????
# PART 4 ? PGLS REGRESSION  (y ~ x, with optional log transforms)
# ?????????????????????????????????????????????????????????????????????????????

def pgls_regression(
    x_vals: np.ndarray,
    y_vals: np.ndarray,
    vcv_df: pd.DataFrame,
    species: list,
    log_x: bool = False,
    log_y: bool = False,
    n_pred: int = 200,
) -> dict:
    """
    PGLS simple regression:  y ~ x   (one predictor + intercept).

    Parameters
    ----------
    x_vals, y_vals : arrays aligned with `species`
    vcv_df         : phylogenetic VCV
    species        : list of species names aligned with x_vals / y_vals
    log_x, log_y   : apply log10 transform before fitting
    n_pred         : number of points in the fitted curve

    Returns
    -------
    dict: slope, intercept, r, p, se_slope, se_intercept,
          n_species, df_residual,
          x_fit, y_fit, ci_lo, ci_hi  (for plotting, in original scale)
          species_used
    or None on failure.
    """
    x_arr = np.asarray(x_vals, dtype=float)
    y_arr = np.asarray(y_vals, dtype=float)
    sp_arr = np.asarray(species)

    # mask NaN and non-positive (for log)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if log_x:
        mask &= (x_arr > 0)
    if log_y:
        mask &= (y_arr > 0)

    sp_ok  = sp_arr[mask].tolist()
    x_ok   = x_arr[mask]
    y_ok   = y_arr[mask]

    if len(sp_ok) < 3:
        return None

    vcv_sub, present = _subset_vcv(vcv_df, sp_ok)
    idx = [sp_ok.index(s) for s in present]
    xv  = x_ok[idx]
    yv  = y_ok[idx]

    lx = np.log10(xv) if log_x else xv
    ly = np.log10(yv) if log_y else yv

    # Check degeneracy
    if np.std(lx) < 1e-12:
        return None

    X = np.column_stack([np.ones(len(lx)), lx])

    try:
        beta, se, resid, sigma2, t_stats, p_vals, df_r = _gls_core(
            X, ly, vcv_sub.values)
    except Exception as exc:
        warnings.warn(f"  [PGLS regression] GLS failed: {exc}")
        return None

    intercept, slope = float(beta[0]), float(beta[1])
    se_int,  se_sl   = float(se[0]),   float(se[1])
    t_slope,  p_slope = float(t_stats[1]), float(p_vals[1])

    # Pearson r on GLS residuals vs fitted (indicative, not phylo-corrected)
    ly_fit   = X @ beta
    r_val, _ = _stats.pearsonr(ly, ly_fit)

    # Prediction curve + 95% prediction interval
    x_fit_l  = np.linspace(lx.min(), lx.max(), n_pred)
    y_fit_l  = intercept + slope * x_fit_l
    n_v      = len(lx)
    t_crit   = _stats.t.ppf(0.975, df=df_r)
    x_mean_l = lx.mean()
    ssx      = float(np.sum((lx - x_mean_l) ** 2))
    se_pred  = np.sqrt(sigma2 * (1 + 1 / n_v + (x_fit_l - x_mean_l) ** 2
                                  / max(ssx, 1e-12)))

    # Back-transform
    x_plot = 10 ** x_fit_l if log_x else x_fit_l
    if log_y:
        y_plot  = 10 ** y_fit_l
        ci_lo   = 10 ** (y_fit_l - t_crit * se_pred)
        ci_hi   = 10 ** (y_fit_l + t_crit * se_pred)
    else:
        y_plot  = y_fit_l
        ci_lo   = y_fit_l - t_crit * se_pred
        ci_hi   = y_fit_l + t_crit * se_pred

    return dict(
        slope        = slope,
        intercept    = intercept,
        r            = r_val,
        p            = p_slope,
        se_slope     = se_sl,
        se_intercept = se_int,
        n_species    = n_v,
        df_residual  = df_r,
        x_fit        = x_plot,
        y_fit        = y_plot,
        ci_lo        = ci_lo,
        ci_hi        = ci_hi,
        species_used = present,
        log_x        = log_x,
        log_y        = log_y,
    )


# ?????????????????????????????????????????????????????????????????????????????
# PART 5 ? AX-LEVEL HELPERS (drop-in replacements)
# ?????????????????????????????????????????????????????????????????????????????

def ols_on_ax(ax, xv, yv, log_x=False):
    """
    [DEPRECATED ? use pgls_on_ax]
    Plain OLS regression line + 95 % prediction band on ax.
    Kept for backward-compatibility.
    """
    import matplotlib.ticker as ticker
    xv = np.asarray(xv, dtype=float)
    yv = np.asarray(yv, dtype=float)
    lx = np.log10(np.clip(xv, 1e-9, None)) if log_x else xv
    valid = np.isfinite(lx) & np.isfinite(yv)
    if valid.sum() < 3:
        return None
    lxv, lyv = lx[valid], yv[valid]
    ssx = float(np.sum((lxv - lxv.mean()) ** 2))
    if ssx < 1e-12:
        ax.axhline(float(np.median(lyv)), color="k", lw=1.2,
                   linestyle="--", alpha=0.55)
        return None
    slope, intercept, r, p, _ = _stats.linregress(lxv, lyv)
    n_v    = valid.sum()
    t_crit = _stats.t.ppf(0.975, df=n_v - 2)
    x_fit  = np.linspace(lxv.min(), lxv.max(), 200)
    y_fit  = intercept + slope * x_fit
    resid  = np.sum((lyv - (intercept + slope * lxv)) ** 2) / (n_v - 2)
    se     = np.sqrt(resid * (1 + 1 / n_v + (x_fit - lxv.mean()) ** 2
                               / ssx))
    x_plot = 10 ** x_fit if log_x else x_fit
    ax.plot(x_plot, y_fit, color="k", lw=1.5, linestyle="--",
            alpha=0.8, zorder=3)
    ax.fill_between(x_plot, y_fit - t_crit * se,
                    y_fit + t_crit * se,
                    color="k", alpha=0.15, zorder=2)
    p_str = f"{p:.3f}" if p >= 0.001 else "<0.001"
    ax.text(0.05, 0.95,
            f"slope={slope:.2f}  r={r:.2f}  p={p_str}",
            transform=ax.transAxes, fontsize=8, va="top", color="#333")
    return slope, r, p


def pgls_on_ax(ax, xv, yv, vcv_df, species, log_x=False, log_y=False):
    """
    PGLS regression line + 95 % prediction band drawn directly on `ax`.

    Parameters
    ----------
    ax           : matplotlib Axes
    xv, yv       : 1-D arrays aligned with `species`
    vcv_df       : phylogenetic VCV DataFrame
    species      : list of species names
    log_x, log_y : log10-transform before fitting

    Returns
    -------
    (slope, r, p) on success, or None on failure / degeneracy.
    """
    import matplotlib.ticker as ticker

    res = pgls_regression(xv, yv, vcv_df, species,
                          log_x=log_x, log_y=log_y)
    if res is None:
        return None

    ax.plot(res["x_fit"], res["y_fit"],
            color="k", lw=1.5, linestyle="--", alpha=0.8, zorder=3)
    ax.fill_between(res["x_fit"], res["ci_lo"], res["ci_hi"],
                    color="k", alpha=0.15, zorder=2)
    p_str = f"{res['p']:.3f}" if res["p"] >= 0.001 else "<0.001"
    ax.text(0.05, 0.95,
            f"PGLS slope={res['slope']:.2f}  r={res['r']:.2f}  p={p_str}",
            transform=ax.transAxes, fontsize=8, va="top", color="#333")

    if log_x:
        import matplotlib.ticker as tck
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(tck.LogFormatterSciNotation())
    if log_y:
        ax.set_yscale("log")

    return res["slope"], res["r"], res["p"]


# ?????????????????????????????????????????????????????????????????????????????
# PART 6 ? ALL-PAIRS PGLS ACROSS REGIONS
# ?????????????????????????????????????????????????????????????????????????????

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
    For every ordered pair (region_x, region_y) of cortical regions, fit:

        log10(mean_y_per_species)  ~  log10(mean_x_per_species)

    using PGLS with phylogenetic covariance from `vcv_df`.

    Parameters
    ----------
    df_morph    : long-format morphometry DataFrame with columns
                  [species, region, hemisphere, atlas_level, <metric_col>]
    metric_col  : 'surface_area_mm2', 'volume_mm3', or 'thickness_mm'
    vcv_df      : phylogenetic VCV (from load_newick_vcv)
    regions     : list of region names to include; None = use all
    atlas_level : int, filter if atlas_level column present
    hemisphere  : 'bilateral', 'left', or 'right'
    min_species : minimum number of species with data to attempt PGLS
    log_x/log_y : whether to log-transform before fitting (default: both True,
                  standard allometric scaling)

    Returns
    -------
    DataFrame with columns:
        region_x, region_y, slope, intercept, r, p, p_fdr, se_slope,
        n_species, df_residual, significant_fdr05, mean_x, mean_y
    Sorted by p (ascending).
    """
    # ?? filter atlas level & hemisphere ?????????????????????????????????????
    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == atlas_level]

    if hemisphere and hemisphere != "all":
        if hemisphere in df.get("hemisphere", pd.Series()).unique():
            df = df[df["hemisphere"] == hemisphere]
        else:
            # Fall back to bilateral if requested hemi absent
            bil = df[df["hemisphere"] == "bilateral"] if "hemisphere" in df.columns else df
            df  = bil if not bil.empty else df

    if regions is None:
        regions = sorted(df["region"].dropna().unique().tolist())

    # ?? compute per-species mean for each region ??????????????????????????
    id_cols = [c for c in ("species",) if c in df.columns]
    region_means: dict = {}
    for reg in regions:
        sub = df[df["region"] == reg]
        if sub.empty:
            continue
        means = sub.groupby(id_cols)[metric_col].mean()
        region_means[reg] = means.to_dict() if isinstance(means, pd.Series) else {}
        # Flatten if multi-level index
        if isinstance(means.index, pd.MultiIndex):
            region_means[reg] = {k[0]: v for k, v in means.items()}
        else:
            region_means[reg] = means.to_dict()

    available_regions = [r for r in regions if r in region_means
                         and len(region_means[r]) >= min_species]

    print(f"  [PGLS pairs] regions with ?{min_species} species: "
          f"{len(available_regions)} / {len(regions)}")

    # ?? run all ordered pairs ?????????????????????????????????????????????
    results = []
    pairs   = list(itertools.combinations(available_regions, 2))
    print(f"  [PGLS pairs] testing {len(pairs)} pairs ...")

    for reg_x, reg_y in pairs:
        means_x = region_means[reg_x]
        means_y = region_means[reg_y]

        # Common species
        common_sp = [s for s in vcv_df.index
                     if s in means_x and s in means_y
                     and pd.notna(means_x[s]) and pd.notna(means_y[s])]
        if len(common_sp) < min_species:
            continue

        xv = np.array([means_x[s] for s in common_sp])
        yv = np.array([means_y[s] for s in common_sp])

        res = pgls_regression(xv, yv, vcv_df, common_sp,
                              log_x=log_x, log_y=log_y)
        if res is None:
            continue

        results.append(dict(
            region_x     = reg_x,
            region_y     = reg_y,
            slope        = res["slope"],
            intercept    = res["intercept"],
            r            = res["r"],
            p            = res["p"],
            se_slope     = res["se_slope"],
            n_species    = res["n_species"],
            df_residual  = res["df_residual"],
            mean_x       = float(np.mean(xv)),
            mean_y       = float(np.mean(yv)),
            species_used = ", ".join(res["species_used"]),
        ))

    if not results:
        warnings.warn("  [PGLS pairs] No results returned.")
        return pd.DataFrame()

    df_res = pd.DataFrame(results).sort_values("p").reset_index(drop=True)

    # FDR (Benjamini-Hochberg)
    n   = len(df_res)
    rnk = _stats.rankdata(df_res["p"].values, method="ordinal")
    df_res["p_fdr"] = np.minimum(1.0, df_res["p"].values * n / rnk)
    df_res["significant_fdr05"] = df_res["p_fdr"] < 0.05

    n_sig = df_res["significant_fdr05"].sum()
    print(f"  [PGLS pairs] {n} pairs tested | {n_sig} FDR<0.05")
    return df_res


# ?????????????????????????????????????????????????????????????????????????????
# PART 7 ? PLOTTING
# ?????????????????????????????????????????????????????????????????????????????

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
    Scatter plot for one (region_x, region_y) PGLS result.

    Each dot = one species (per-species mean), coloured by PHYLO_ORDER.
    The PGLS regression line + 95 % prediction band are overlaid.

    Parameters
    ----------
    result_row   : one row from the DataFrame returned by run_pgls_all_region_pairs
    df_morph     : same morphometry DataFrame used in run_pgls_all_region_pairs
    metric_col   : e.g. 'surface_area_mm2'
    metric_label : y-axis label, e.g. 'Surface area (mm²)'
    vcv_df       : phylogenetic VCV
    output_path  : PNG save path

    Returns
    -------
    output_path on success, None on failure.
    """
    import matplotlib.ticker as ticker

    reg_x = result_row["region_x"]
    reg_y = result_row["region_y"]

    # ?? filter data ??????????????????????????????????????????????????????
    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == atlas_level]
    if hemisphere != "all" and "hemisphere" in df.columns:
        h_df = df[df["hemisphere"] == hemisphere]
        df   = h_df if not h_df.empty else df[df["hemisphere"] == "bilateral"]

    id_cols = ["species"]
    means_x = df[df["region"] == reg_x].groupby(id_cols)[metric_col].mean()
    means_y = df[df["region"] == reg_y].groupby(id_cols)[metric_col].mean()

    common_sp = sorted(set(means_x.index) & set(means_y.index)
                       & set(vcv_df.index))
    if len(common_sp) < 3:
        warnings.warn(f"  [plot pair] Not enough species for {reg_x} vs {reg_y}")
        return None

    xv = means_x.loc[common_sp].values
    yv = means_y.loc[common_sp].values
    sp_colors = _species_colors(common_sp)

    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(5.5, 5))

        # Scatter
        for sp, x, y in zip(common_sp, xv, yv):
            ax.scatter(x, y, color=sp_colors[sp], s=55, zorder=5,
                       edgecolors="k", linewidths=0.5, label=sp)

        # PGLS line
        pgls_on_ax(ax, xv, yv, vcv_df, common_sp,
                   log_x=log_x, log_y=log_y)

        if log_x:
            ax.set_xscale("log")
            ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
        if log_y:
            ax.set_yscale("log")

        ax.set_xlabel(f"{reg_x}  ({metric_label})", fontsize=9)
        ax.set_ylabel(f"{reg_y}  ({metric_label})", fontsize=9)

        p_fdr = float(result_row["p_fdr"]) if "p_fdr" in result_row.index else np.nan
        sig_star = "***" if p_fdr < 0.001 else (
                   "**"  if p_fdr < 0.01  else (
                   "*"   if p_fdr < 0.05  else "n.s."))
        ax.set_title(f"{reg_x}  ×  {reg_y}\nFDR {p_fdr:.3f} {sig_star}",
                     fontweight="bold", fontsize=9)

        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                   for sp in _phylo_sort(common_sp)]
        ax.legend(handles=handles, fontsize=7,
                  bbox_to_anchor=(1.01, 1), loc="upper left",
                  borderaxespad=0.)
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
    """
    Symmetric heatmap of `value_col` for all region pairs.
    Cells with FDR<0.05 are marked with *.
    """
    from matplotlib.colors import TwoSlopeNorm
    import matplotlib.cm as _cm

    if df_results.empty:
        warnings.warn("plot_pgls_summary_heatmap: empty results")
        return None

    n = len(regions)
    mat      = np.full((n, n), np.nan)
    mat_sig  = np.full((n, n), False)
    idx_map  = {r: i for i, r in enumerate(regions)}

    for _, row in df_results.iterrows():
        ri = idx_map.get(row["region_x"])
        ci = idx_map.get(row["region_y"])
        if ri is None or ci is None:
            continue
        val = float(row[value_col]) if pd.notna(row[value_col]) else np.nan
        mat[ri, ci]   = val
        mat[ci, ri]   = val          # symmetric
        sig = bool(row.get("significant_fdr05", False))
        mat_sig[ri, ci] = sig
        mat_sig[ci, ri] = sig

    finite = mat[np.isfinite(mat)]
    if finite.size == 0:
        warnings.warn("plot_pgls_summary_heatmap: no finite values")
        return None

    vext = float(np.nanpercentile(np.abs(finite), 95))
    norm = TwoSlopeNorm(vmin=-vext, vcenter=0, vmax=vext) if vext > 0 else None
    cmap = "RdBu_r"

    short = [r[:20] for r in regions]
    w = max(5, n * 0.9 + 2.5)
    h = max(4, n * 0.65 + 1.5)

    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(w, h))
        im = ax.imshow(mat, cmap=cmap, norm=norm,
                       aspect="auto", interpolation="nearest")
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
        ax.set_yticks(range(n))
        ax.set_yticklabels(short, fontsize=8)
        ax.set_title(title + "  (FDR<0.05 marked *)",
                     fontweight="bold", fontsize=10)
        cb = fig.colorbar(
            _cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
            fraction=0.03, pad=0.04)
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
    title: str = "PGLS allometric slopes ? top pairs",
):
    """
    Dot-plot of PGLS results, sorted by p-value.
    Dot size ? -log10(p_fdr), colour = red if significant.

    Parameters
    ----------
    df_results : from run_pgls_all_region_pairs
    output_path: PNG path
    top_n      : show only top_n pairs (by p-value)
    """
    if df_results.empty:
        warnings.warn("plot_pgls_dotplot: empty results")
        return None

    df_plot = df_results.sort_values("p").head(top_n).reset_index(drop=True)
    df_plot = df_plot.iloc[::-1].reset_index(drop=True)   # ascending p on y

    labels = [f"{r['region_x'][:16]} ? {r['region_y'][:16]}"
              for _, r in df_plot.iterrows()]
    colors = ["#D55E00" if s else "#999999"
              for s in df_plot["significant_fdr05"]]
    sizes  = np.clip(-np.log10(df_plot["p_fdr"].clip(lower=1e-10)) * 30,
                     20, 300)

    h = max(4, len(df_plot) * 0.32 + 1.5)
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(7.5, h))
        y_pos = np.arange(len(df_plot))

        ax.errorbar(df_plot[value_col], y_pos,
                    xerr=df_plot["se_slope"],
                    fmt="none", color="#666", lw=0.8, zorder=2)
        ax.scatter(df_plot[value_col], y_pos,
                   c=colors, s=sizes, zorder=3,
                   edgecolors="k", linewidths=0.4, alpha=0.9)

        ax.axvline(0, color="k", lw=0.8, linestyle="--", alpha=0.5)
        ax.axvline(1, color="#aaa", lw=0.6, linestyle=":", alpha=0.4)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel(f"PGLS {value_col}  ± SE", fontsize=9)
        ax.set_title(title, fontweight="bold", fontsize=10)

        sig_p  = mpatches.Patch(facecolor="#D55E00", label="FDR < 0.05")
        nosig_p = mpatches.Patch(facecolor="#999999", label="n.s.")
        ax.legend(handles=[sig_p, nosig_p], fontsize=8, loc="lower right")

        # p_fdr annotations
        for i, row in df_plot.iterrows():
            p_str = (f"p={row['p_fdr']:.3f}" if row["p_fdr"] >= 0.001
                     else "p<0.001")
            ax.text(ax.get_xlim()[1] * 0.98, i, p_str,
                    va="center", ha="right", fontsize=6, color="#444")

        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


# ?????????????????????????????????????????????????????????????????????????????
# PART 8 ? PGLS for brain-scaling (allometric body/brain weight regressions)
# ?????????????????????????????????????????????????????????????????????????????

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
    Allometric brain-scaling figure with PGLS (3 panels):
      A) metric vs body weight  (log-log)
      B) metric vs brain weight (log-log)
      C) metric_y_region vs metric_x_region (log-log, internal scaling)

    Replaces the plain-OLS version in the main script.
    """
    import matplotlib.ticker as ticker

    # filter
    df = df_morph.copy()
    if "atlas_level" in df.columns:
        df = df[df["atlas_level"] == atlas_level]
    if hemisphere != "all" and "hemisphere" in df.columns:
        h_df = df[df["hemisphere"] == hemisphere]
        df   = h_df if not h_df.empty else df[df["hemisphere"] == "bilateral"]

    id_cols = ["species"]
    df_y_reg = df[df["region"] == y_region].groupby(id_cols)[metric_col].mean()
    if df_y_reg.empty:
        warnings.warn(f"plot_brain_scaling_pgls: no data for y_region={y_region}")
        return None

    sp_list   = list(df_y_reg.index)
    y_vals    = df_y_reg.values
    sp_colors = _species_colors(sp_list)

    panels = []
    if species_body_weight_g:
        bw = np.array([species_body_weight_g.get(s, np.nan) for s in sp_list])
        panels.append((bw, "Body weight (g)", True, False))
    if species_brain_weight_g:
        brw = np.array([species_brain_weight_g.get(s, np.nan) for s in sp_list])
        panels.append((brw, "Brain weight (g)", True, False))
    if x_region:
        df_x_reg = df[df["region"] == x_region].groupby(id_cols)[metric_col].mean()
        x_reg_vals = np.array([
            float(df_x_reg.loc[s]) if s in df_x_reg.index else np.nan
            for s in sp_list])
        panels.append((x_reg_vals, f"{x_region[:24]} ({metric_label})", True, True))

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

        for pi, (ax, (xv, xlabel, log_x, log_y)) in enumerate(
                zip(axes, panels)):
            for sp, x, y in zip(sp_list, xv, y_vals):
                if not (np.isfinite(x) and np.isfinite(y)):
                    continue
                ax.scatter(x, y, color=sp_colors[sp], s=40, zorder=5,
                           alpha=0.85, edgecolors="k", linewidths=0.5, label=sp)

            pgls_on_ax(ax, xv, y_vals, vcv_df, sp_list,
                       log_x=log_x, log_y=log_y)

            ax.set_xlabel(xlabel, fontsize=9)
            ax.set_ylabel(f"{y_region} ({metric_label})" if pi == 0 else "",
                          fontsize=9)
            ax.text(-0.12, 1.02, chr(65 + pi), transform=ax.transAxes,
                    fontsize=13, fontweight="bold", va="bottom", ha="left")

        handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                   for sp in _phylo_sort(sp_list)]
        plt.tight_layout()
        fig.legend(handles=handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.06),
                   ncol=min(len(handles), 7),
                   fontsize=8, frameon=False, title="Species")
        fig.subplots_adjust(bottom=0.18)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        print(f"  [plot] {output_path}")
    return output_path


# ?????????????????????????????????????????????????????????????????????????????
# PART 9 ? FC CONNECTIONS PGLS
# ?????????????????????????????????????????????????????????????????????????????

def run_pgls_fc_connections(
    df_fc: pd.DataFrame,
    fc_columns: list,
    vcv_df: pd.DataFrame,
    x_columns: list = None,
    x_labels: list = None,
    x_log: list = None,
    min_species: int = 4,
) -> pd.DataFrame:
    """
    Run PGLS for each FC column vs each x column (body/brain weight, internal FC).

    Parameters
    ----------
    df_fc      : DataFrame with 'species' and FC metric columns
    fc_columns : list of y column names (e.g. ['FC_mpfc_pmc', ...])
    vcv_df     : phylogenetic VCV
    x_columns  : list of x column names
    x_labels   : display labels for x columns
    x_log      : list of bool ? whether to log-transform each x
    min_species: minimum species required

    Returns
    -------
    DataFrame: fc_col, x_col, slope, r, p, p_fdr, se_slope, n_species, significant_fdr05
    """
    x_columns = x_columns or []
    x_labels  = x_labels  or x_columns
    x_log     = x_log     or [False] * len(x_columns)

    results = []
    for fc_col in fc_columns:
        if fc_col not in df_fc.columns:
            continue
        # compute per-species mean
        sp_means_y = df_fc.groupby("species")[fc_col].mean().dropna()

        for x_col, x_lbl, lx in zip(x_columns, x_labels, x_log):
            if x_col not in df_fc.columns:
                continue
            sp_means_x = df_fc.groupby("species")[x_col].mean().dropna()
            common_sp  = sorted(set(sp_means_x.index) & set(sp_means_y.index)
                                 & set(vcv_df.index))
            if len(common_sp) < min_species:
                continue
            xv = sp_means_x.loc[common_sp].values
            yv = sp_means_y.loc[common_sp].values

            res = pgls_regression(xv, yv, vcv_df, common_sp,
                                  log_x=lx, log_y=False)
            if res is None:
                continue
            results.append(dict(
                fc_col      = fc_col,
                x_col       = x_col,
                x_label     = x_lbl,
                slope       = res["slope"],
                r           = res["r"],
                p           = res["p"],
                se_slope    = res["se_slope"],
                n_species   = res["n_species"],
                species_used = ", ".join(res["species_used"]),
            ))

    if not results:
        return pd.DataFrame()

    df_r = pd.DataFrame(results).sort_values("p").reset_index(drop=True)
    n    = len(df_r)
    rnk  = _stats.rankdata(df_r["p"].values, method="ordinal")
    df_r["p_fdr"]            = np.minimum(1.0, df_r["p"].values * n / rnk)
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
    """
    One figure per FC connection.  Panels: body weight, brain weight,
    and optionally somato L-R FC.  All regressions use PGLS.

    Parameters
    ----------
    connections          : list of (col_key, ..., y_label) tuples
                           (same format as FC_CONNECTIONS_OF_INTEREST in main)
    output_path_template : e.g. '/path/to/scaling_fc_norm.png'
                           actual files will be '_mpfc_pmc.png' etc.
    """
    import matplotlib.ticker as ticker

    df_fc = df_fc.copy()
    if species_body_weight_g:
        df_fc["body_g"]  = df_fc["species"].map(species_body_weight_g)
    if species_brain_weight_g:
        df_fc["brain_g"] = df_fc["species"].map(species_brain_weight_g)

    species_order = _phylo_sort(df_fc["species"].unique())
    sp_colors     = _species_colors(species_order)
    out_paths     = []

    for conn in connections:
        col_key = conn[0]
        y_label = conn[-1]

        if col_key not in df_fc.columns or df_fc[col_key].notna().sum() == 0:
            print(f"  [FC PGLS] {col_key}: no data ? skipped")
            continue

        x_specs = []
        if "body_g" in df_fc.columns:
            x_specs.append(("body_g", "Body weight (g)", True))
        if "brain_g" in df_fc.columns:
            x_specs.append(("brain_g", "Brain weight (g)", True))
        # internal reference (somato) ? skip if plotting that pair itself
        if (col_key != "FC_somato_lr"
                and "FC_somato_lr" in df_fc.columns
                and df_fc["FC_somato_lr"].notna().any()):
            x_specs.append(("FC_somato_lr", "Somato L-R (r/mean|r|)", False))

        x_specs = [(xc, xl, lx) for xc, xl, lx in x_specs
                   if xc in df_fc.columns and df_fc[xc].notna().any()]
        if not x_specs:
            continue

        n_panels = len(x_specs)
        w, h     = figsize or (n_panels * 4.5, 4.5)

        with plt.rc_context(PAPER_RC):
            fig, axes = plt.subplots(1, n_panels, figsize=(w, h))
            if n_panels == 1:
                axes = [axes]

            for pi, (ax, (x_col, x_label, log_x)) in enumerate(
                    zip(axes, x_specs)):
                sub = df_fc[[col_key, x_col, "species"]].dropna()
                if sub.empty:
                    ax.set_visible(False)
                    continue

                sp_list = sub["species"].tolist()
                for sp in _phylo_sort(sub["species"].unique()):
                    sp_sub = sub[sub["species"] == sp]
                    if sp_sub.empty:
                        continue
                    ax.scatter(sp_sub[x_col], sp_sub[col_key],
                               color=sp_colors[sp], s=40, zorder=5,
                               alpha=0.85, edgecolors="k", linewidths=0.5,
                               label=sp)

                pgls_on_ax(ax, sub[x_col].values.astype(float),
                           sub[col_key].values.astype(float),
                           vcv_df, sp_list, log_x=log_x, log_y=False)

                ax.axhline(0, color="#aaa", lw=0.8, linestyle=":", zorder=1)
                ax.set_xlabel(x_label, fontsize=9)
                ax.set_ylabel(y_label if pi == 0 else "", fontsize=9)
                ax.text(-0.12, 1.02, chr(65 + pi), transform=ax.transAxes,
                        fontsize=13, fontweight="bold", va="bottom", ha="left")

            handles = [mpatches.Patch(facecolor=sp_colors[sp], label=sp)
                       for sp in _phylo_sort(df_fc["species"].unique())]
            plt.tight_layout()
            fig.legend(handles=handles, loc="lower center",
                       bbox_to_anchor=(0.5, -0.06),
                       ncol=min(len(handles), 7),
                       fontsize=8, frameon=False, title="Species")
            fig.subplots_adjust(bottom=0.18)

            safe_key = col_key.replace("FC_", "")
            out = output_path_template.replace(".png", f"_{safe_key}.png")
            os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
            fig.savefig(out, bbox_inches="tight", dpi=200)
            plt.close(fig)
            print(f"  [plot] {out}")
            out_paths.append(out)

    return out_paths