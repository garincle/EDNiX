"""
phylo_kmult.py -- Figure 4e: phylogenetic signal of the functional fingerprints with the
multivariate Blomberg K (K_mult, Adams 2014, Syst. Biol. 63:685-697) instead of a Mantel test.

Why K_mult (Clément 2026-11-25): the Spearman Mantel of the former panel reduces the 7 species x
30 values to 21 ranked distances against 21 MRCA ages that take 6 distinct values; with 7
species it needs rho >= ~0.68 for P < 0.05 and it is known to be weak in comparative analyses
(Harmon & Glor 2010; Guillot & Rousset 2013). K_mult uses the whole profile and the Brownian
covariance of the time tree (phylo.phylo_matrices): K = 1 is the similarity among relatives
expected under Brownian motion, K < 1 less. P = exact permutation of the species over the tips
(7! = 5,040).

Profiles. Group profile of a dataset = mean of its (normalised) run profiles, z-scored
(fig4_species_structure.group_profiles); species profile = mean of its datasets.
  raw      the species profiles themselves
  resLOSO  each dataset's group profile minus its projection on the common profile of the OTHER
           species (median of their species medians, z-scored). Leaving the species out avoids
           subtracting a component it helped estimate; what remains is what is specific to the
           species, beyond the pattern shared across species (the target of the run selection).

Pools: all runs, the selected runs, and N_DRAWS draws of as many runs per dataset at random.
selected vs random P = share of draws with K >= the selection; selected vs all runs: Delta K with
the same tip permutation applied to both pools (as phylo Mantel Delta rho).
Leave-one-species-out: K and P with each species removed (6 species, 720 permutations), all runs
and selected, so that no single species is shown to carry the signal on its own.
"""
import itertools

import numpy as np
import pandas as pd

N_DRAWS = 200
SEED = 0
KINDS = ("raw", "resLOSO")
KIND_LABEL = {"raw": "raw profiles", "resLOSO": "species-specific\n(shared removed)"}


# ── profiles ──────────────────────────────────────────────────────────────────
def _z(v):
    v = np.asarray(v, float)
    return (v - np.nanmean(v)) / np.nanstd(v)


def _common(groups, sp_of_ds, species):
    return _z(np.nanmedian([np.nanmedian([groups[d] for d in groups if sp_of_ds[d] == s], axis=0)
                            for s in species], axis=0))


def _resid(g, C):
    ok = np.isfinite(g) & np.isfinite(C)
    return g - np.cov(g[ok], C[ok])[0, 1] / np.var(C[ok], ddof=1) * C


def loso_residuals(groups, sp_of_ds, species):
    """{dataset: group profile minus its projection on the common profile of the other species}."""
    out = {}
    for s in species:
        Cs = _common(groups, sp_of_ds, [t for t in species if t != s])
        for d in groups:
            if sp_of_ds[d] == s:
                out[d] = _resid(groups[d], Cs)
    return out


def species_matrix(groups, sp_of_ds, species):
    """species x values: mean of the dataset profiles of each species (NaN -> 0)."""
    Y = np.vstack([np.nanmean([groups[d] for d in groups if sp_of_ds[d] == s], axis=0)
                   for s in species])
    return np.where(np.isfinite(Y), Y, 0.0)


def matrices(groups, sp_of_ds, species):
    groups = {d: g for d, g in groups.items() if sp_of_ds[d] in species}
    return dict(raw=species_matrix(groups, sp_of_ds, species),
                resLOSO=species_matrix(loso_residuals(groups, sp_of_ds, species), sp_of_ds,
                                       species))


# ── K_mult ────────────────────────────────────────────────────────────────────
class KTest:
    """K_mult for a fixed species set and tree, with every tip permutation precomputed."""

    def __init__(self, C):
        C = np.asarray(C, float)
        self.n = len(C)
        self.Ci = np.linalg.inv(C)
        one = np.ones(self.n)
        self.w = one @ self.Ci / (one @ self.Ci @ one)       # GLS weights of the phylogenetic mean
        self.k_exp = (np.trace(C) - self.n / (one @ self.Ci @ one)) / (self.n - 1)
        self.perms = np.array(list(itertools.permutations(range(self.n))))

    def k(self, Y):
        """K of Y (species x values, rows in the order of C), or of a stack (perm x species x values)."""
        Y = np.asarray(Y, float)
        a = np.einsum("i,...ij->...j", self.w, Y)
        E = Y - a[..., None, :]
        num = np.einsum("...ij,...ij->...", E, E)
        den = np.einsum("...ij,ik,...kj->...", E, self.Ci, E)
        return (num / den) / self.k_exp

    def null(self, Y):
        return self.k(np.asarray(Y, float)[self.perms])

    def test(self, Y):
        """(K, one-sided exact P)."""
        obs = float(self.k(Y))
        return obs, float(np.mean(self.null(Y) >= obs - 1e-12))

    def delta(self, Y_all, Y_sel):
        """K(selected) - K(all runs); null = the same tip permutation applied to both."""
        d = float(self.k(Y_sel) - self.k(Y_all))
        nd = self.null(Y_sel) - self.null(Y_all)
        return dict(delta=d, p=float(np.mean(nd >= d - 1e-12)), n_perm=len(nd))


def _tree(species):
    from Paper.analysis.fmri_fingerprints import phylo
    return phylo.phylo_matrices(list(species))[1]


def compute(vecs, ds, selected, sp_of_ds, n_draws=N_DRAWS, seed=SEED, species=None, C=None):
    """Every number of Figure 4e and of the leave-one-species-out check.
    vecs/ds/selected/sp_of_ds as selection_validation.run_vectors."""
    from Paper.analysis import fig4_species_structure as ss
    from Paper.analysis.fmri_fingerprints import phylo
    if species is None:
        species = phylo.order_species(sorted(set(sp_of_ds)))
    species = list(species)
    T = KTest(_tree(species) if C is None else C)
    sets = np.unique(ds)
    idx_all = {d: np.flatnonzero(ds == d) for d in sets}
    idx_sel = {d: np.flatnonzero((ds == d) & selected) for d in sets}
    idx_sel = {d: i for d, i in idx_sel.items() if len(i)}
    g_all = ss.group_profiles(vecs, ds, idx_all)
    g_sel = ss.group_profiles(vecs, ds, idx_sel)
    M_all, M_sel = matrices(g_all, sp_of_ds, species), matrices(g_sel, sp_of_ds, species)
    rng = np.random.default_rng(seed + 3)
    draws = []
    for _ in range(n_draws):
        g = ss.group_profiles(vecs, ds, {d: rng.choice(idx_all[d], len(idx_sel[d]), replace=False)
                                         for d in idx_sel})
        draws.append(matrices(g, sp_of_ds, species))
    out = dict(species=species, n_draws=n_draws, n_perm=len(T.perms))
    for k in KINDS:
        rk = [T.test(m[k]) for m in draws]
        rand = np.array([r[0] for r in rk])
        sel = T.test(M_sel[k])
        out[k] = dict(all=T.test(M_all[k]), selected=sel, random=rand,
                      random_p=np.array([r[1] for r in rk]),
                      p_vs_random=float((1 + np.sum(rand >= sel[0] - 1e-12)) / (n_draws + 1)),
                      delta=T.delta(M_all[k], M_sel[k]))
    out["loso_species"] = loso_species(g_all, g_sel, sp_of_ds, species)
    return out


def loso_species(g_all, g_sel, sp_of_ds, species):
    """K and exact P with each species left out, all runs and selected, raw and resLOSO."""
    rows = []
    for drop in species:
        keep = [s for s in species if s != drop]
        T = KTest(_tree(keep))
        for pool, g in (("all runs", g_all), ("selected", g_sel)):
            M = matrices(g, sp_of_ds, keep)
            for k in KINDS:
                K, p = T.test(M[k])
                rows.append(dict(dropped=drop, pool=pool, profiles=k, K=K, p=p,
                                 n_species=len(keep), n_perm=len(T.perms)))
    return pd.DataFrame(rows)


# ── tables / text ─────────────────────────────────────────────────────────────
def table(kres):
    rows = []
    for k in KINDS:
        r = kres[k]
        rows += [dict(profiles=k, pool="all runs", K=r["all"][0], p=r["all"][1]),
                 dict(profiles=k, pool="selected", K=r["selected"][0], p=r["selected"][1],
                      p_selected_vs_random=r["p_vs_random"],
                      delta_K_selected_minus_all=r["delta"]["delta"],
                      p_delta=r["delta"]["p"])]
        rows += [dict(profiles=k, pool="random", K=x, p=p_)
                 for x, p_ in zip(r["random"], r["random_p"])]
    return pd.DataFrame(rows)


def result_lines(kres):
    lines = [f"e  phylogenetic signal: multivariate Blomberg K (K_mult, Adams 2014) of the species "
             f"profiles on the time tree ({len(kres['species'])} species; K = 1 Brownian "
             f"expectation); exact P over the {kres['n_perm']} tip permutations, one-sided"]
    for k in KINDS:
        r = kres[k]
        lines.append(
            f"  {KIND_LABEL[k].replace(chr(10), ' ')}: all runs K = {r['all'][0]:.3f} "
            f"(P = {r['all'][1]:.4g}); random median K = {np.median(r['random']):.3f} (95th pct "
            f"{np.percentile(r['random'], 95):.3f}; median P of the draws "
            f"{np.median(r['random_p']):.3g}); selected K = {r['selected'][0]:.3f} "
            f"(P = {r['selected'][1]:.4g}); selected vs random P = {r['p_vs_random']:.4g} "
            f"({kres['n_draws']} draws); selected - all runs Delta K = "
            f"{r['delta']['delta']:+.3f} (P = {r['delta']['p']:.4g})")
    lo = kres["loso_species"]
    for k in KINDS:
        for pool in ("all runs", "selected"):
            t = lo[(lo.profiles == k) & (lo.pool == pool)]
            lines.append(f"  leave-one-species-out, {KIND_LABEL[k].replace(chr(10), ' ')}, "
                         f"{pool}: " + "; ".join(f"-{r.dropped} K = {r.K:.2f} (P = {r.p:.3g})"
                                                 for r in t.itertuples()))
    return lines


# ── drawing ───────────────────────────────────────────────────────────────────
def draw_k(ax, kres, fs=5.5):
    """Figure 4e: two groups (raw / species-specific), each all runs / random / selected."""
    from Paper.analysis.fmri_fingerprints import selection_validation as sv
    gap = 4.0
    tops, bots = [], []
    for g, k in enumerate(KINDS):
        r = kres[k]
        x0 = g * gap
        vp = ax.violinplot(r["random"], positions=[x0 + 1], widths=0.7, showextrema=False)
        for b_ in vp["bodies"]:
            b_.set_facecolor("#BDBDBD")
            b_.set_edgecolor("none")
            b_.set_alpha(1)
        ax.plot(x0, r["all"][0], marker="o", ms=3.4, color=sv.POOL_COLORS["random"], mec="white",
                mew=0.4, zorder=4)
        ax.plot(x0 + 2, r["selected"][0], marker="D", ms=3.6, color=sv.SEL_COLOR, mec="white",
                mew=0.4, zorder=4)
        tops.append(max(r["all"][0], r["selected"][0], np.max(r["random"])))
        bots.append(min(r["all"][0], r["selected"][0], np.min(r["random"])))
    top, bot = max(tops), min(bots)
    span = top - bot
    for g, k in enumerate(KINDS):
        r, x0, tk = kres[k], g * gap, tops[g]
        for x, y, p, col in ((x0, r["all"][0], r["all"][1], "#4D4D4D"),
                             (x0 + 2, r["selected"][0], r["selected"][1], sv.SEL_COLOR)):
            ax.text(x, y + 0.05 * span, sv._pfmt(p), ha="center", va="bottom",
                    fontsize=fs - 1.5, color=col)
        # median P of the random draws under their violin (no clash with the 'all runs' label)
        ax.text(x0 + 1, np.min(r["random"]) - 0.04 * span,
                "med.\n" + sv._pfmt(float(np.median(r["random_p"]))), ha="center", va="top",
                fontsize=fs - 1.5, color="#4D4D4D", linespacing=0.9)
        sv._bracket(ax, x0 + 1, x0 + 2, tk + 0.24 * span, 0.04 * span, sv._pfmt(r["p_vs_random"]),
                    sv.SEL_COLOR, fs - 0.3)
        sv._bracket(ax, x0, x0 + 2, tk + 0.42 * span, 0.04 * span, sv._pfmt(r["delta"]["p"]),
                    sv.SEL_COLOR, fs - 0.3)
        ax.text(x0 + 1, bot - 0.62 * span, KIND_LABEL[k], ha="center", va="top", fontsize=fs - 0.5)
    ax.set_ylim(bot - 0.30 * span, top + 0.68 * span)
    ax.set_xticks([g * gap + i for g in range(len(KINDS)) for i in range(3)])
    ax.set_xticklabels(["all", "rand.\n10", "sel.\n10"] * len(KINDS), fontsize=fs - 1.0)
    ax.set_xlim(-0.9, gap * (len(KINDS) - 1) + 2.9)
    ax.set_ylabel("phylogenetic signal, K$_{mult}$", fontsize=fs)
    ax.tick_params(labelsize=fs - 0.5, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)


def plot_loso_species(kres, out, width_mm=120, height_mm=48):
    """Supplementary panel: K with each species left out (all runs grey circles, selected red
    diamonds; filled = P < 0.05), the full-set values as dashed lines; raw | species-specific."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Paper.analysis.fmri_fingerprints import profile_stats as fa
    from Paper.analysis.fmri_fingerprints import selection_validation as sv
    lo, sps = kres["loso_species"], kres["species"]
    fs = 5.5
    with fa.paper_style():
        fig, axes = plt.subplots(1, 2, figsize=(width_mm / 25.4, height_mm / 25.4), sharey=True)
        for ax, k in zip(axes, KINDS):
            for pool, col, mk, dx in (("all runs", sv.POOL_COLORS["random"], "o", -0.15),
                                      ("selected", sv.SEL_COLOR, "D", 0.15)):
                t = lo[(lo.profiles == k) & (lo.pool == pool)].set_index("dropped").loc[sps]
                for i, r in enumerate(t.itertuples()):
                    ax.plot(i + dx, r.K, marker=mk, ms=3.2, mec=col, mew=0.7,
                            color=col if r.p < 0.05 else "white", zorder=3)
                full = kres[k]["all" if pool == "all runs" else "selected"][0]
                ax.axhline(full, color=col, lw=0.5, ls="--", zorder=1)
            ax.set_xticks(range(len(sps)))
            ax.set_xticklabels([f"-{s}" for s in sps], rotation=45, ha="right", fontsize=fs - 0.5)
            ax.set_title(KIND_LABEL[k].replace("\n", " "), fontsize=fs + 0.5, fontweight="bold",
                         loc="left", pad=1.5)
            ax.tick_params(labelsize=fs - 0.5, length=2)
            for s_ in ("top", "right"):
                ax.spines[s_].set_visible(False)
        axes[0].set_ylabel("K$_{mult}$, one species left out", fontsize=fs)
        from matplotlib.lines import Line2D
        axes[1].legend([Line2D([], [], color=sv.POOL_COLORS["random"], marker="o", ls="none", ms=3),
                        Line2D([], [], color=sv.SEL_COLOR, marker="D", ls="none", ms=3),
                        Line2D([], [], color="k", marker="o", ls="none", ms=3, mfc="white")],
                       ["all runs", "selected 10", "open: P ≥ 0.05"], fontsize=fs - 0.5,
                       frameon=False, loc="lower right", handletextpad=0.3)
        fa.save_paper(fig, out)
        plt.close(fig)
    return out
