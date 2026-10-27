"""
phylo.py
─────────────
Phylogenetic comparative tests for Figure 5, replacing the 1-D
"distance to consensus vs divergence from human" regression as the
primary phylogenetic test.

Why: species are not independent data points (Felsenstein 1985), and
"divergence from human" only uses one row of the tree (Rat and Mouse
share the same abscissa). The tests here use the FULL time-calibrated
tree:

1. MANTEL (primary): pairwise functional distance between every species
   pair vs pairwise phylogenetic distance (MRCA age). With <= 8 species
   the permutation null is enumerated EXACTLY (7! = 5040), no Monte
   Carlo error. Known caveat (Harmon & Glor 2010; Guillot & Rousset
   2013): Mantel tests on phylogenetically structured data can have
   inflated type-I error -- report alongside PGLS/K, never alone.

2. BLOMBERG'S K per region (Blomberg, Garland & Ives 2003): is the
   between-species variation of each region's connectivity value
   phylogenetically structured? K = 1 is Brownian-motion expectation,
   K << 1 = less signal than BM. Exact permutation p + BH-FDR across
   regions.

3. PGLS (Grafen 1989; Pagel's lambda, Pagel 1999): generalized least
   squares with the tree's Brownian covariance, lambda by ML. Used for
   (a) phylogenetic signal of the F6 distance trait (intercept-only),
   (b) confound checks: distance ~ covariate (tSNR, motion, n runs, any
   user-supplied species covariate such as brain volume). With n = 7
   species, lambda is poorly identified (Freckleton et al. 2002;
   Munkemuller et al. 2012) -- every model is ALSO reported at the fixed
   Brownian value lambda = 1 so results never hinge on that estimate.

Tree: the pipeline loads the time-calibrated Newick tree given in
cfg["TREE_PATH"] (use_newick_tree) and every MRCA age comes from it. The
hand-entered SPECIES_CLADES below are only a fallback when TREE_PATH=None
(approximate TimeTree medians, not for publication). The Spearman Mantel
test depends only on the order of node ages; Pearson Mantel, Blomberg's K
and PGLS depend on their values.
"""

import itertools
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as _stats


# ═══════════════════════════════════════════════════════════════════════════
# Tree
# ═══════════════════════════════════════════════════════════════════════════

# Topology: (((((Human, Macaque) Catarrhini, Marmoset) Anthropoidea,
# Mouselemur) Primates, (Rat, Mouse) Muridae) Euarchontoglires, Dog)
# Boreoeutheria. Node ages: same provenance and caveat as
# profile_stats.DIVERGENCE_TIME_MYA (hand-entered approximate TimeTree
# medians, not queried). The Rat-Mouse age (13 Mya) is the only node not
# shared with that dict. The Spearman Mantel test depends only on the ORDER
# of node ages; Pearson Mantel, Blomberg's K and PGLS depend on their values.
SPECIES_CLADES = [
    (("Rat", "Mouse"), 13.0),
    (("Human", "Macaque"), 29.0),
    (("Human", "Macaque", "Marmoset"), 43.0),
    (("Human", "Macaque", "Marmoset", "Mouselemur"), 74.0),
    (("Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse"), 87.0),
    (("Human", "Macaque", "Marmoset", "Mouselemur", "Rat", "Mouse", "Dog"), 94.0),
]

# Leaf order for plotting (tree order, so blocks show up on heatmaps)
TREE_ORDER = ["Human", "Macaque", "Marmoset", "Mouselemur", "Mouse", "Rat", "Dog"]


# ═══════════════════════════════════════════════════════════════════════════
# Time-calibrated tree from a Newick file (replaces SPECIES_CLADES once loaded)
# ═══════════════════════════════════════════════════════════════════════════

_TREE_PAIR_AGES = None      # {frozenset({sp_a, sp_b}): MRCA age}, set by use_newick_tree
_TREE_SOURCE = None


def parse_newick(text):
    """
    Minimal Newick parser (no external dependency). Handles nested clades,
    tip and internal labels (quoted or not), branch lengths, [comments].
    Returns (nodes, root) with nodes[i] = dict(name, length, parent, children).
    """
    text = re.sub(r"\[[^\]]*\]", "", text.strip())
    if text.endswith(";"):
        text = text[:-1]
    nodes, stack, i = [], [], 0

    def new(parent):
        nodes.append(dict(name="", length=0.0, parent=parent, children=[]))
        if parent is not None:
            nodes[parent]["children"].append(len(nodes) - 1)
        return len(nodes) - 1

    def read_label(j):
        if j < len(text) and text[j] in "'\"":
            q = text[j]; k = text.index(q, j + 1)
            return text[j + 1:k], k + 1
        k = j
        while k < len(text) and text[k] not in ",():;":
            k += 1
        return text[j:k].strip(), k

    root = new(None)
    cur = root
    while i < len(text):
        c = text[i]
        if c == "(":
            cur = new(cur); stack.append(cur); i += 1
        elif c == ",":
            parent = nodes[cur]["parent"]
            cur = new(parent); i += 1
        elif c == ")":
            cur = nodes[cur]["parent"]; i += 1
            lab, i = read_label(i)
            if lab:
                nodes[cur]["name"] = lab
        elif c == ":":
            k = i + 1
            while k < len(text) and text[k] not in ",();":
                k += 1
            nodes[cur]["length"] = float(text[i + 1:k]); i = k
        elif c.isspace():
            i += 1
        else:
            lab, i = read_label(i)
            nodes[cur]["name"] = lab
    # the outer "(...)" created a child of the synthetic root: collapse it
    if len(nodes[root]["children"]) == 1 and not nodes[root]["name"]:
        root = nodes[root]["children"][0]
        nodes[root]["parent"] = None
    return nodes, root


def newick_pair_ages(path, species_to_tip, verbose=True):
    """
    MRCA age for every pair of the species in species_to_tip
    ({pipeline species: tip name in the tree}), read from a time-calibrated
    Newick file. Tips not listed are ignored (the tree can contain more
    species). Age of MRCA(a, b) = mean over the two tips of (root-to-tip
    depth - root-to-MRCA depth), so a slightly non-ultrametric tree (rounding
    in the source) gives one consistent value per pair; the maximum
    deviation from ultrametricity among the used tips is reported.
    """
    with open(path, encoding="utf-8") as f:
        nodes, root = parse_newick(f.read())
    tip_idx = {n["name"]: i for i, n in enumerate(nodes) if not n["children"]}
    missing = {sp: t for sp, t in species_to_tip.items() if t not in tip_idx}
    if missing:
        raise KeyError(f"tip(s) not found in {path}: {missing}. Tips available: "
                       f"{sorted(tip_idx)}")
    depth = {root: 0.0}
    order = [root]
    for i in order:
        for c in nodes[i]["children"]:
            depth[c] = depth[i] + nodes[c]["length"]
            order.append(c)

    def ancestors(i):
        out = []
        while i is not None:
            out.append(i); i = nodes[i]["parent"]
        return out

    ages = {}
    sp = sorted(species_to_tip)
    for a, b in itertools.combinations(sp, 2):
        ta, tb = tip_idx[species_to_tip[a]], tip_idx[species_to_tip[b]]
        anc_b = set(ancestors(tb))
        m = next(x for x in ancestors(ta) if x in anc_b)
        ages[frozenset((a, b))] = 0.5 * ((depth[ta] - depth[m]) + (depth[tb] - depth[m]))
    bad = {tuple(sorted(k)): v for k, v in ages.items() if not v > 0}
    if bad:
        raise ValueError(f"{path}: non-positive MRCA age for {len(bad)} species pair(s) "
                         f"(e.g. {list(bad.items())[:3]}) -- the tree has no (or zero) branch "
                         f"lengths, or two species map to the same tip; a time-calibrated "
                         f"tree is required")
    tip_depths = [depth[tip_idx[t]] for t in species_to_tip.values()]
    if verbose:
        print(f"  [tree] {path}: {len(tip_idx)} tips, {len(sp)} used; max deviation from "
              f"ultrametricity among used tips = {max(tip_depths) - min(tip_depths):.4f} Mya")
    return ages


def use_newick_tree(path, species_to_tip, verbose=True):
    """Make mrca_age / phylo_matrices use the Newick tree instead of the
    hand-entered SPECIES_CLADES. Returns the pair-age dict."""
    global _TREE_PAIR_AGES, _TREE_SOURCE
    _TREE_PAIR_AGES = newick_pair_ages(path, species_to_tip, verbose=verbose)
    _TREE_SOURCE = path
    return _TREE_PAIR_AGES


def mrca_age(a, b, clades=SPECIES_CLADES):
    """Age (Mya) of the most recent common ancestor of species a and b:
    from the Newick tree when use_newick_tree() was called, otherwise from
    the hand-entered SPECIES_CLADES."""
    if a == b:
        return 0.0
    if _TREE_PAIR_AGES is not None:
        key = frozenset((a, b))
        if key not in _TREE_PAIR_AGES:
            raise KeyError(f"no age for {a!r}-{b!r} in the loaded tree {_TREE_SOURCE}")
        return float(_TREE_PAIR_AGES[key])
    ages = [age for members, age in clades if a in members and b in members]
    if not ages:
        raise KeyError(f"no clade in the tree contains both {a!r} and {b!r} -- "
                       f"add them to SPECIES_CLADES")
    return float(min(ages))


def check_tree_against_divergence_times(divergence_times, verbose=True):
    """Every 'Human vs X' MRCA age must equal DIVERGENCE_TIME_MYA[X],
    otherwise F6 and the tests here would silently use different trees."""
    bad = []
    for sp, t in divergence_times.items():
        if sp == "Human":
            continue
        try:
            age = mrca_age("Human", sp)
        except KeyError:
            bad.append((sp, t, None))
            continue
        if abs(age - t) > 1e-9:
            bad.append((sp, t, age))
    if verbose:
        if bad:
            print(f"  [WARN] tree / DIVERGENCE_TIME_MYA mismatch: "
                  f"{[(s, f'dict={t}', f'tree={a}') for s, t, a in bad]}")
        else:
            print(f"  [tree] SPECIES_CLADES consistent with DIVERGENCE_TIME_MYA "
                  f"for all {len(divergence_times) - 1} Human-X pairs")
    return not bad


def order_species(species):
    known = [s for s in TREE_ORDER if s in species]
    return known + sorted(s for s in species if s not in TREE_ORDER)


def phylo_matrices(species):
    """
    Returns (D_mrca, C) for the given species order:
      D_mrca[i, j] = MRCA age (0 on diagonal) -- used for Mantel.
      C[i, j]      = shared branch length from the root = root - MRCA(i, j),
                     C[i, i] = root -- the Brownian-motion covariance
                     (ultrametric tree, so C is positive definite).
    root = MRCA age of the whole species SET given (not always 94 Mya:
    at an atlas level where Dog drops out, the root moves to 87).
    """
    n = len(species)
    D = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        D[i, j] = D[j, i] = mrca_age(species[i], species[j])
    root = float(D.max())
    C = root - D
    np.fill_diagonal(C, root)
    return D, C


# ═══════════════════════════════════════════════════════════════════════════
# Functional distances
# ═══════════════════════════════════════════════════════════════════════════

def _ordered_values(profile, common, value_col):
    """Same contract as profile_stats._ordered_common_values (copied so
    this module has no package-path dependency and can be unit-tested
    standalone)."""
    if profile is None or "base_name" not in profile.columns:
        return np.full(len(common), np.nan)
    lut = (profile.assign(_b=profile["base_name"].str.lower().str.strip())
                  .groupby("_b")[value_col].mean())
    return np.array([lut.get(b, np.nan) for b in common], dtype=float)


def _zscore(v):
    v = np.asarray(v, dtype=float)
    f = np.isfinite(v)
    if f.sum() < 2 or np.nanstd(v[f]) == 0:
        return np.full_like(v, np.nan)
    out = np.full_like(v, np.nan)
    out[f] = (v[f] - v[f].mean()) / v[f].std()
    return out


def _vec(prof, common, value_col, zscore):
    """zscore follows cfg["ZSCORE"]: True/"zscore", "scale" (RMS, no
    centering, sign preserved) or False/None/"none" -- same contract as
    profile_stats._normalize_values."""
    v = _ordered_values(prof, common, value_col)
    if zscore is True or zscore == "zscore":
        return _zscore(v)
    if zscore == "scale":
        f = np.isfinite(v)
        rms = float(np.sqrt(np.mean(v[f] ** 2))) if f.sum() >= 2 else 0.0
        if rms == 0:
            return np.full_like(v, np.nan)
        out = np.full_like(v, np.nan)
        out[f] = v[f] / rms
        return out
    return v


def functional_distance_matrix(profiles_by_label, common, value_col="median_r",
                               zscore=True, min_regions=3):
    """
    profiles_by_label: {label: (species, profile_df)} (group-map profiles).

    Species-pair distance = MEAN over every cross-species BIDS pair of
    (1 - Pearson r) -- NOT the distance between species-collapsed
    profiles, so a 2-BIDS species' internal disagreement (Macaque
    BenHamed vs ZhuGarin) is carried into the distance rather than
    hidden by averaging the two profiles first.

    Returns (D_func, species_order, bids_pairs_df). bids_pairs_df has
    one row per BIDS pair, including WITHIN-species pairs (same_species
    True) -- those are reported, never used in the Mantel test.
    """
    vecs = {lab: (sp, _vec(p, common, value_col, zscore))
            for lab, (sp, p) in profiles_by_label.items()}
    rows = []
    for (la, (sa, va)), (lb, (sb, vb)) in itertools.combinations(sorted(vecs.items()), 2):
        ok = np.isfinite(va) & np.isfinite(vb)
        if ok.sum() < min_regions:
            continue
        r = float(np.corrcoef(va[ok], vb[ok])[0, 1])
        rows.append(dict(label_a=la, label_b=lb, species_a=sa, species_b=sb,
                         same_species=sa == sb, pearson_r=r, distance=1 - r,
                         n_regions=int(ok.sum())))
    pairs = pd.DataFrame(rows)
    species = order_species(sorted({sp for sp, _ in vecs.values()}))
    n = len(species)
    D = np.full((n, n), np.nan)
    np.fill_diagonal(D, 0.0)
    if not pairs.empty:
        cross = pairs[~pairs["same_species"]]
        for i, j in itertools.combinations(range(n), 2):
            m = cross[((cross.species_a == species[i]) & (cross.species_b == species[j])) |
                      ((cross.species_a == species[j]) & (cross.species_b == species[i]))]
            if not m.empty:
                D[i, j] = D[j, i] = float(m["distance"].mean())
    return D, species, pairs


# ═══════════════════════════════════════════════════════════════════════════
# Mantel
# ═══════════════════════════════════════════════════════════════════════════

def _corr(x, y, method):
    if method == "spearman":
        x, y = _stats.rankdata(x), _stats.rankdata(y)
    return float(np.corrcoef(x, y)[0, 1])


def mantel_test(D_x, D_y, method="pearson", max_exact=8, n_perm=10000, seed=0):
    """
    Correlation between the upper triangles of D_x and D_y; null built
    by permuting species labels of D_y (rows and columns together).
    n <= max_exact: all n! permutations enumerated (exact p, includes
    the identity). Otherwise n_perm random permutations with the
    (1 + count) / (n_perm + 1) correction.

    Ties in D_y (e.g. Human-Rat and Human-Mouse both 87 Mya) make some
    permutations reproduce the observed matrix -- that is correct exact
    behaviour, and why the smallest achievable p with this tree is
    > 1/5040.

    Returns dict(r, p_greater, p_two_sided, n_species, n_pairs,
    n_perm, exact, null).
    """
    n = D_x.shape[0]
    iu = np.triu_indices(n, 1)
    x, y = D_x[iu], D_y[iu]
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return None
    obs = _corr(x[ok], y[ok], method)

    exact = n <= max_exact
    if exact:
        perms = itertools.permutations(range(n))
    else:
        rng = np.random.default_rng(seed)
        perms = (rng.permutation(n) for _ in range(n_perm))

    null = []
    for p in perms:
        p = np.asarray(p)
        yp = D_y[np.ix_(p, p)][iu]
        okp = np.isfinite(x) & np.isfinite(yp)
        null.append(_corr(x[okp], yp[okp], method))
    null = np.asarray(null)
    tol = 1e-12
    if exact:
        p_g = float(np.mean(null >= obs - tol))
        p_2 = float(np.mean(np.abs(null) >= abs(obs) - tol))
    else:
        p_g = (1 + int(np.sum(null >= obs - tol))) / (len(null) + 1)
        p_2 = (1 + int(np.sum(np.abs(null) >= abs(obs) - tol))) / (len(null) + 1)
    return dict(r=obs, p_greater=p_g, p_two_sided=p_2, method=method,
                n_species=n, n_pairs=int(ok.sum()), n_perm=len(null),
                exact=exact, null=null)


def run_mantel(profiles_by_label, common, value_col="median_r", zscore=True,
               verbose=True, label=""):
    """Convenience wrapper: functional matrix + tree matrix + Mantel
    (Pearson and Spearman). Returns dict with everything needed to plot
    and export."""
    D_func, species, pairs = functional_distance_matrix(
        profiles_by_label, common, value_col=value_col, zscore=zscore)
    if len(species) < 4:
        if verbose:
            print(f"  [SKIP] Mantel{label}: only {len(species)} species (need >= 4)")
        return None
    D_phylo, C = phylo_matrices(species)
    res_p = mantel_test(D_func, D_phylo, method="pearson")
    res_s = mantel_test(D_func, D_phylo, method="spearman")
    if res_p is None:
        if verbose:
            print(f"  [SKIP] Mantel{label}: too few finite species pairs")
        return None
    if verbose:
        print(f"  [Mantel{label}] {len(species)} species, {res_p['n_pairs']} pairs, "
              f"{'exact' if res_p['exact'] else 'Monte Carlo'} null "
              f"({res_p['n_perm']} permutations)")
        print(f"    Pearson  r={res_p['r']:+.3f}  p(greater)={res_p['p_greater']:.4f}  "
              f"p(two-sided)={res_p['p_two_sided']:.4f}")
        print(f"    Spearman r={res_s['r']:+.3f}  p(greater)={res_s['p_greater']:.4f}  "
              f"p(two-sided)={res_s['p_two_sided']:.4f}")
        within = pairs[pairs["same_species"]] if not pairs.empty else pairs
        cross = pairs[~pairs["same_species"]] if not pairs.empty else pairs
        if len(within):
            print(f"    within-species BIDS distance: mean={within['distance'].mean():.3f} "
                  f"({len(within)} pair(s)) vs cross-species mean="
                  f"{cross['distance'].mean():.3f} ({len(cross)} pair(s))")
    return dict(species=species, D_func=D_func, D_phylo=D_phylo, C=C,
                bids_pairs=pairs, mantel_pearson=res_p, mantel_spearman=res_s)


# ═══════════════════════════════════════════════════════════════════════════
# Blomberg's K
# ═══════════════════════════════════════════════════════════════════════════

def blomberg_k(x, C):
    """Blomberg, Garland & Ives (2003) K. x: trait per species (tree
    order of C). K = 1 under Brownian motion."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    Ci = np.linalg.inv(C)
    one = np.ones(n)
    denom = one @ Ci @ one
    a = (one @ Ci @ x) / denom
    r = x - a
    mse0 = (r @ r) / (n - 1)
    mse = (r @ Ci @ r) / (n - 1)
    expected = (np.trace(C) - n / denom) / (n - 1)
    return float((mse0 / mse) / expected)


def blomberg_k_test(x, C, max_exact=8, n_perm=10000, seed=0):
    """K + permutation p (shuffle trait across tips), exact for n <= 8.
    One-sided: p = P(K_null >= K_obs)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    k_obs = blomberg_k(x, C)
    if n <= max_exact:
        null = np.array([blomberg_k(x[list(p)], C)
                         for p in itertools.permutations(range(n))])
        p = float(np.mean(null >= k_obs - 1e-12))
    else:
        rng = np.random.default_rng(seed)
        null = np.array([blomberg_k(rng.permutation(x), C) for _ in range(n_perm)])
        p = (1 + int(np.sum(null >= k_obs - 1e-12))) / (n_perm + 1)
    return k_obs, p


def _bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        prev = min(prev, p[idx] * n / rank)
        q[idx] = prev
    return q


def regionwise_phylo_signal(profiles_by_label, common, value_col="median_r",
                            zscore=True, verbose=True):
    """
    Per region: species trait = median across that species' BIDS of the
    (within-profile z-scored) region value. Blomberg's K + exact
    permutation p, BH-FDR across regions. Regions missing in any species
    are skipped (K needs the full tree).
    """
    by_species = {}
    for lab, (sp, p) in profiles_by_label.items():
        by_species.setdefault(sp, []).append(_vec(p, common, value_col, zscore))
    species = order_species(list(by_species))
    if len(species) < 4:
        if verbose:
            print(f"  [SKIP] region-wise K: only {len(species)} species")
        return pd.DataFrame()
    with np.errstate(invalid="ignore"):
        traits = np.vstack([np.nanmedian(np.vstack(by_species[s]), axis=0) for s in species])
    _, C = phylo_matrices(species)

    rows = []
    for j, region in enumerate(common):
        x = traits[:, j]
        if not np.isfinite(x).all() or np.nanstd(x) == 0:
            continue
        k, p = blomberg_k_test(x, C)
        rows.append(dict(region=region, K=k, p_perm=p, n_species=len(species)))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["q_fdr"] = _bh_fdr(df["p_perm"].values)
    df = df.sort_values("K", ascending=False).reset_index(drop=True)
    if verbose:
        n_sig = int((df["q_fdr"] < 0.05).sum())
        print(f"  [K] {len(df)} region(s) tested on {len(species)} species; "
              f"{n_sig} with q_FDR < 0.05; median K = {df['K'].median():.2f} "
              f"(K = 1 is the Brownian-motion expectation)")
    return df


# ═══════════════════════════════════════════════════════════════════════════
# PGLS
# ═══════════════════════════════════════════════════════════════════════════

def _lambda_cov(C, lam):
    V = C * lam
    np.fill_diagonal(V, np.diag(C))
    return V


def _gls_fit(y, X, V):
    n, k = X.shape
    Vi = np.linalg.inv(V)
    XtVi = X.T @ Vi
    beta = np.linalg.solve(XtVi @ X, XtVi @ y)
    r = y - X @ beta
    rss = float(r @ Vi @ r)
    sign, logdet = np.linalg.slogdet(V)
    sigma2_ml = rss / n
    ll = -0.5 * (n * np.log(2 * np.pi * sigma2_ml) + logdet + n)
    sigma2 = rss / (n - k)
    cov = sigma2 * np.linalg.inv(XtVi @ X)
    return beta, np.sqrt(np.diag(cov)), ll


def pgls(y, X, C, lambda_grid=np.linspace(0.0, 1.0, 101), fixed_lambda=None):
    """
    y ~ X under V = lambda-transformed C. lambda by ML on a grid unless
    fixed_lambda is given. t-tests on coefficients with n - k df.
    LRT for lambda vs 0 uses the boundary-corrected 0.5 * chi2(1).
    X must already include the intercept column.
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    n, k = X.shape
    if n - k < 1:
        return None
    if fixed_lambda is not None:
        lam = float(fixed_lambda)
        beta, se, ll = _gls_fit(y, X, _lambda_cov(C, lam))
        ll0 = _gls_fit(y, X, _lambda_cov(C, 0.0))[2]
    else:
        fits = [(lam, *_gls_fit(y, X, _lambda_cov(C, lam))) for lam in lambda_grid]
        lam, beta, se, ll = max(fits, key=lambda f: f[3])
        ll0 = fits[0][3]
    t = beta / se
    p = 2 * _stats.t.sf(np.abs(t), df=n - k)
    lr = max(0.0, 2 * (ll - ll0))
    p_lambda = 0.5 * _stats.chi2.sf(lr, 1) if lr > 0 else 1.0
    return dict(lam=lam, beta=beta, se=se, t=t, p=p, loglik=ll,
                lr_lambda=lr, p_lambda=p_lambda, n=n, k=k)


def pgls_battery(species_df, trait_col, covariate_cols, verbose=True):
    """
    species_df: one row per species with 'species', trait_col and
    covariate columns. Fits, for lambda = ML AND lambda = 1 (Brownian):
      - intercept-only (-> phylogenetic signal of the trait, lambda LRT)
      - trait ~ each covariate separately (n = 7 species cannot support
        more than one covariate per model: df = n - 2 = 5).
    Returns a long DataFrame.
    """
    species_df = species_df.dropna(subset=[trait_col]).copy()
    species = order_species(species_df["species"].tolist())
    species_df = species_df.set_index("species").loc[species].reset_index()
    if len(species) < 4:
        if verbose:
            print(f"  [SKIP] PGLS: only {len(species)} species")
        return pd.DataFrame()
    _, C = phylo_matrices(species)
    y = species_df[trait_col].values.astype(float)

    rows = []
    models = [("intercept_only", None)] + [(c, c) for c in covariate_cols
                                            if c in species_df.columns]
    for name, cov in models:
        if cov is None:
            X = np.ones((len(y), 1))
            ok = np.ones(len(y), dtype=bool)
        else:
            xc = pd.to_numeric(species_df[cov], errors="coerce").values
            ok = np.isfinite(xc)
            if ok.sum() < 4 or np.nanstd(xc[ok]) == 0:
                continue
            X = np.column_stack([np.ones(ok.sum()), xc[ok]])
        C_ok = C[np.ix_(ok, ok)]
        for lam_mode, fixed in (("ML", None), ("BM(lambda=1)", 1.0)):
            res = pgls(y[ok], X, C_ok, fixed_lambda=fixed)
            if res is None:
                continue
            slope_i = 1 if cov is not None else 0
            rows.append(dict(model=name, lambda_mode=lam_mode, lam=res["lam"],
                             n_species=res["n"],
                             coef=float(res["beta"][slope_i]),
                             se=float(res["se"][slope_i]),
                             t=float(res["t"][slope_i]),
                             p=float(res["p"][slope_i]),
                             lr_lambda_vs_0=res["lr_lambda"],
                             p_lambda_vs_0=res["p_lambda"]))
    df = pd.DataFrame(rows)
    if verbose and not df.empty:
        print(f"  [PGLS] trait={trait_col}, {len(species)} species "
              f"(lambda poorly identified at this n -- read the BM rows too)")
        print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    return df


def species_covariates_from_qc(df_qc_by_label, species_of_label, metrics):
    """Median of each QC metric per species, pooling that species' BIDS."""
    rows = []
    for label, df in (df_qc_by_label or {}).items():
        sp = species_of_label.get(label)
        if sp is None or df is None or df.empty:
            continue
        for m in metrics:
            if m in df.columns:
                vals = pd.to_numeric(df[m], errors="coerce").dropna()
                for v in vals:
                    rows.append(dict(species=sp, metric=m, value=float(v)))
    if not rows:
        return pd.DataFrame(columns=["species"])
    long = pd.DataFrame(rows)
    return (long.groupby(["species", "metric"])["value"].median()
            .unstack("metric").reset_index())


# ═══════════════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════════════

def plot_mantel(mres, title=None, output_path=None, dpi=300):
    """2x2: functional matrix | phylogenetic matrix | pairwise scatter |
    exact null."""
    if mres is None:
        print("  [plot] nothing to plot")
        return None
    sp = mres["species"]
    Df, Dp = mres["D_func"], mres["D_phylo"]
    mp, ms = mres["mantel_pearson"], mres["mantel_spearman"]
    pairs = mres["bids_pairs"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    for ax, M, ttl, cmap, fmt in ((axes[0, 0], Df, "Functional distance (1 - r)", "viridis", "{:.2f}"),
                                  (axes[0, 1], Dp, "MRCA age (Mya)", "magma", "{:.0f}")):
        im = ax.imshow(M, cmap=cmap)
        ax.set_xticks(range(len(sp))); ax.set_yticks(range(len(sp)))
        ax.set_xticklabels(sp, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(sp, fontsize=8)
        for i in range(len(sp)):
            for j in range(len(sp)):
                if i != j and np.isfinite(M[i, j]):
                    frac = (M[i, j] - np.nanmin(M)) / (np.nanmax(M) - np.nanmin(M) + 1e-12)
                    ax.text(j, i, fmt.format(M[i, j]), ha="center", va="center",
                            fontsize=7, color="black" if frac > 0.6 else "white")
        ax.set_title(ttl, fontsize=10, fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[1, 0]
    if not pairs.empty:
        cross = pairs[~pairs["same_species"]].copy()
        cross["mrca"] = [mrca_age(a, b) for a, b in zip(cross.species_a, cross.species_b)]
        ax.scatter(cross["mrca"], cross["distance"], s=14, color="#BBBBBB",
                   label="cross-species BIDS pairs", zorder=2)
        within = pairs[pairs["same_species"]]
        if len(within):
            ax.scatter(np.zeros(len(within)), within["distance"], s=60, marker="D",
                       color="#D55E00", edgecolors="black",
                       label="within-species BIDS pairs (not tested)", zorder=4)
            for k, (_, r) in enumerate(within.sort_values("distance").iterrows()):
                ax.annotate(f"{r['label_a']} vs {r['label_b']} ({r['distance']:.2f})",
                            (0, r["distance"]), textcoords="offset points",
                            xytext=(8, 10 * k), fontsize=6,
                            arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888"))
    iu = np.triu_indices(len(sp), 1)
    ax.scatter(Dp[iu], Df[iu], s=55, color="#0072B2", edgecolors="black",
               label="species-pair mean (tested)", zorder=3)
    ax.set_xlabel("MRCA age of the pair (Mya)", fontsize=9)
    ax.set_ylabel("Functional distance (1 - r)", fontsize=9)
    ax.set_title(f"Mantel Pearson r={mp['r']:+.3f}, p(>)={mp['p_greater']:.4f}\n"
                 f"Spearman r={ms['r']:+.3f}, p(>)={ms['p_greater']:.4f}",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=7)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax = axes[1, 1]
    ax.hist(mp["null"], bins=40, color="#AAAAAA", edgecolor="white")
    ax.axvline(mp["r"], color="#D55E00", lw=2.5,
               label=f"observed r={mp['r']:+.3f}")
    ax.set_xlabel("Mantel r under label permutation", fontsize=9)
    ax.set_ylabel("permutations", fontsize=9)
    ax.set_title(f"{'Exact' if mp['exact'] else 'Monte Carlo'} null "
                 f"({mp['n_perm']} permutations)", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    fig.suptitle(title or "Pairwise functional vs phylogenetic distance (Mantel)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig


def plot_regionwise_k(k_df, title=None, output_path=None, dpi=300, alpha=0.05):
    if k_df is None or k_df.empty:
        print("  [plot] nothing to plot")
        return None
    fig, ax = plt.subplots(figsize=(max(7, 0.5 * len(k_df)), 5))
    colors = ["#009E73" if q < alpha else "#BBBBBB" for q in k_df["q_fdr"]]
    ax.bar(range(len(k_df)), k_df["K"], color=colors, edgecolor="black", lw=0.4)
    for i, (k, p, q) in enumerate(zip(k_df["K"], k_df["p_perm"], k_df["q_fdr"])):
        ax.text(i, k, f"p={p:.3f}", ha="center", va="bottom", fontsize=6, rotation=90)
    ax.axhline(1.0, color="#D55E00", ls="--", lw=1.2, label="Brownian expectation (K = 1)")
    ax.set_xticks(range(len(k_df)))
    ax.set_xticklabels(k_df["region"], rotation=75, ha="right", fontsize=7)
    ax.set_ylabel("Blomberg's K", fontsize=10)
    ax.set_title(title or f"Phylogenetic signal per region (green: q_FDR < {alpha})",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {output_path}")
    return fig
