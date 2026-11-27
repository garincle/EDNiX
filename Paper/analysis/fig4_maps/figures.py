"""Figures of the new Figure 4, the harmonisation figure and their supplements, drawn
at print size (183 mm) with the paper style (paper_style / profile_stats.paper_style).

  Figure 4
    a  phylogenetic signal of the all-runs 3D maps (P80/P20, retrosplenial seed): species
       distance (1 - r) vs MRCA age, raw and residual (LOSO common), Spearman Mantel; exact
       null of K_mult over the 7! species permutations
    b  common profile + one radar per dataset over its LOSO consensus (figure4_panels.
       fig4a_common_profile), all-runs maps by default (--profiles selected: kept maps)
    c  dataset x dataset similarity (r), raw profiles and residuals vs the LOSO common
  Harmonisation
    a  amplitude (RMS of the regional mean Fisher z) of every dataset, all runs vs selected,
       grey = random draws of as many runs (5-95 %)
    b  amplitude ratio within / between species, seed maps: null of random draws, all runs
       and selected
    c  same on the correlation matrices (when computed)
  Supplementary
    noise_bins   K per noise bin with its size control, identifiability, species distance
    bins_profiles  profile matrices and similarity matrices per group
    surfaces     mean map of every group (render_common_scale, needs wb_command)
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from Paper.analysis import paper_style as ps  # noqa: E402
from Paper.analysis.fig4_maps import stats as fs  # noqa: E402

MM = ps.MM
GREY = "#B8B8B8"
INK = "#222222"
_BOX = dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85)


def _style():
    from Paper.analysis.fmri_fingerprints import profile_stats as fa
    return fa.paper_style()


def _save(fig, out):
    from Paper.analysis.fmri_fingerprints import profile_stats as fa
    fa.save_paper(fig, out)
    return out + ".png"


def _ptxt(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"


# ── Figure 4 ────────────────────────────────────────────────────────────────────────────
def fig4a(ctx, ph, out, width_mm=183, height_mm=52):
    obs = ph["obs"]["all runs"]
    M = fs.species_profiles(ph["P_all"], ctx.sp_of_ds, ctx.sps, ph["Cg"], ph["Cs"])
    iu = np.triu_indices(len(ctx.sps), 1)
    age = np.asarray(ctx.Dp, float)[iu]
    with _style():
        fig, axes = plt.subplots(1, 3, figsize=(width_mm * MM, height_mm * MM),
                                 gridspec_kw=dict(wspace=0.45))
        for ax, key, D, title in ((axes[0], "M_raw", fs.dist_corr(M["raw"]), "Species profiles"),
                                  (axes[1], "M_resLOSO", fs.dist_corr(M["resLOSO"]),
                                   "Species-specific part (residual)")):
            y = D[iu]
            ax.scatter(age, y, s=7, color=INK, edgecolors="white", linewidths=0.3, zorder=3)
            ua = np.unique(age)
            ax.plot(ua, [y[age == a].mean() for a in ua], color=INK, lw=0.7, zorder=2)
            ax.set_xlabel("Divergence time (Mya)")
            ax.set_ylabel("Species distance (1 - r)")
            ax.set_title(title, loc="left", fontweight="bold")
            rho, p = obs[key][0], obs[key][1]
            ax.text(0.03, 0.97, f"Mantel rho = {rho:+.2f}\n{_ptxt(p)} (exact, 7! perm.)",
                    transform=ax.transAxes, va="top", fontsize=5.5, bbox=_BOX)
        ax = axes[2]
        for j, (key, color, lab) in enumerate((("K_raw", INK, "profiles"),
                                               ("K_resLOSO", "#C0392B", "residuals"))):
            k, p, null = obs[key]
            ax.hist(null, bins=60, color=color, alpha=0.25, density=True, lw=0)
            ax.axvline(k, color=color, lw=1.0)
            ax.text(0.02, 0.98 - 0.16 * j, f"{lab}: K = {k:.2f}, {_ptxt(p)}", color=color,
                    fontsize=5.2, va="top", transform=ax.transAxes, bbox=_BOX)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)          # room for the two lines of text
        ax.set_xlabel("Multivariate K (Adams 2014)")
        ax.set_ylabel("Density (species permutations)")
        ax.set_title("Phylogenetic signal, all runs", loc="left", fontweight="bold")
        return _save(fig, out)


def fig4b(ctx, pr, out):
    from Paper.analysis.fmri_fingerprints import figure4_panels as f4
    from Paper.analysis.fmri_fingerprints import fingerprints as ff
    c = ctx.c
    common = c["common"]
    vf1 = ff.group_vectors(c["ana_species_profiles"], common, c["cfg"]["VALUE_COL"],
                           c["cfg"]["ZSCORE"]) if c.get("ana_species_profiles") else {}
    order = (ff._region_order(ff.consensus(vf1, c["species_of_label"]), common) if vf1
             else np.arange(len(common)))
    f4.fig4a_common_profile(pr["vecs"], pr["dist"], c["species_of_label"], common, order,
                            c.get("region_mask"), out)
    return out + ".png"


def fig4c(ctx, pr, out, width_mm=183, height_mm=74):
    from Paper.analysis.fmri_fingerprints import figure4_panels as f4
    sol = ctx.c["species_of_label"]
    with _style():
        fig = plt.figure(figsize=(width_mm * MM, height_mm * MM))
        side = 46.0
        for k, (D, title) in enumerate(((pr["D"], "Profiles"),
                                        (pr["D_res"], "Residuals (LOSO common removed)"))):
            x0 = 20.0 + k * 92.0
            ax = fig.add_axes([x0 / width_mm, 4.0 / height_mm, side / width_mm, side / height_mm])
            cax = fig.add_axes([(x0 + side + 1.5) / width_mm, 4.0 / height_mm, 1.2 / width_mm,
                                side / height_mm])
            f4._matrix_r(ax, cax, D, sol)
            ax.set_title(title, loc="left", fontweight="bold", pad=24)
        return _save(fig, out)


def figure4(ctx, ph, pr, out_dir):
    pan = os.path.join(out_dir, "panels")
    os.makedirs(pan, exist_ok=True)
    pa = _try(fig4a, ctx, ph, os.path.join(pan, "a_phylogenetic_signal"))
    pb = _try(fig4b, ctx, pr, os.path.join(pan, f"b_profiles_{pr['pool']}"))
    pc = _try(fig4c, ctx, pr, os.path.join(pan, f"c_similarity_{pr['pool']}"))
    base = os.path.join(out_dir, "Figure4")
    ps.compose_images([pa, pb, pc], base, layout=[[0], [1], [2]], letters="abc")
    return base


# ── harmonisation ───────────────────────────────────────────────────────────────────────
def _amp_panel(ax, ctx, am):
    from Paper.analysis.fmri_fingerprints import supp_panels as fsp
    amp = am["amp"]
    labels = fsp._labels_sorted(list(amp.index), ctx.c["species_of_label"])
    x = np.arange(len(labels))
    for i, lab in enumerate(labels):
        r = amp.loc[lab]
        col = ps.species_color(r["species"])
        ax.plot([i, i], [r["random_p5"], r["random_p95"]], color=GREY, lw=4,
                solid_capstyle="round", zorder=1)
        ax.scatter(i, r["all"], s=10, facecolors="white", edgecolors=col, linewidths=0.8, zorder=3)
        ax.scatter(i, r["selected"], s=10, color=col, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels([fsp._lab(l) for l in labels], rotation=55, ha="right",
                       rotation_mode="anchor", fontsize=5)
    for t, lab in zip(ax.get_xticklabels(), labels):
        t.set_color(fsp._tcol(fsp._sp(lab, ctx.c["species_of_label"])))
    ax.set_ylabel("Amplitude (RMS regional mean z)")
    ax.scatter([], [], s=10, facecolors="white", edgecolors=INK, label="all runs")
    ax.scatter([], [], s=10, color=INK, label="selected")
    ax.plot([], [], color=GREY, lw=4, label="random runs, same n (5-95 %)")
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3, fontsize=5)


def _ratio_panel(ax, res, title):
    v = res["rat_sel"]
    ax.hist(v, bins=40, color=GREY, lw=0, density=True)
    a, s = res["ratio"]["all"], res["ratio"]["selected"]
    ax.axvline(a, color=INK, lw=0.8, ls=(0, (2.5, 1.5)))
    ax.axvline(s, color="#C0392B", lw=1.0)
    ytop = ax.get_ylim()[1]
    ax.text(0.02, 0.98, f"all runs (dashed) {a:.2f}", fontsize=5.2, va="top",
            transform=ax.transAxes, bbox=_BOX)
    ax.text(0.02, 0.86, f"selected {s:.2f}, {_ptxt(fs.p_low(v, s))} vs random runs",
            fontsize=5.2, va="top", color="#C0392B", transform=ax.transAxes, bbox=_BOX)
    ax.set_ylim(0, ytop * 1.3)
    ax.set_xlabel("Amplitude difference, same species / between species")
    ax.set_ylabel("Density (random runs)")
    ax.set_title(title, loc="left", fontweight="bold")


def harmonisation(ctx, am, mx, out_dir, width_mm=183, height_mm=58):
    base = os.path.join(out_dir, "Figure4_harmonisation")
    with _style():
        n = 3 if mx is not None else 2
        fig, axes = plt.subplots(1, n, figsize=(width_mm * MM, height_mm * MM),
                                 gridspec_kw=dict(wspace=0.4, width_ratios=[1.5] + [1] * (n - 1)))
        _amp_panel(axes[0], ctx, am)
        axes[0].set_title("Seed-map amplitude per dataset", loc="left", fontweight="bold",
                          pad=12)
        _ratio_panel(axes[1], am, "Seed maps")
        if mx is not None:
            _ratio_panel(axes[2], mx, "Correlation matrices")
        for ax, letter in zip(axes, "abc"):
            ax.text(-0.18, 1.12, letter, transform=ax.transAxes, fontsize=ps.PANEL_LABEL_SIZE,
                    fontweight="bold", va="bottom")
        _save(fig, base)
    return base


# ── supplementary ───────────────────────────────────────────────────────────────────────
def supp_noise_bins(ctx, ph, nb, out_dir, width_mm=183, height_mm=100):
    B = ctx.B
    groups = ["all", "selected"] + [f"bin{b + 1}" for b in range(B)]
    base = os.path.join(out_dir, "Supp_Figure4_noise_bins")
    with _style():
        fig, axes = plt.subplots(2, 3, figsize=(width_mm * MM, height_mm * MM),
                                 gridspec_kw=dict(wspace=0.4, hspace=0.75))
        ax = axes[0, 0]
        xb = np.arange(1, B + 1)
        for key, mk, col in (("K_raw", "o", INK), ("K_resLOSO", "s", "#C0392B")):
            pool = np.array([[x[key][0] for x in part] for part in nb["null"]]).ravel()
            ax.fill_between([0.7, B + 0.3], *np.percentile(pool, [5, 95]), color=col, alpha=0.12,
                            lw=0)
            ax.plot(xb, [o[key][0] for o in nb["obs"]], mk + "-", color=col, ms=3, label=key)
            ax.axhline(ph["obs"]["all runs"][key][0], color=col, lw=0.5, ls=":")
        ax.set_xticks(xb)
        ax.set_xlabel("Noise bin (1 = closest to the selection criterion)")
        ax.set_ylabel("Multivariate K")
        ax.legend(fontsize=5)
        ax.set_title("K per bin (band: random bins, 5-95 %)", loc="left", fontweight="bold")
        x = np.arange(len(groups))

        def band(ax, key):
            for i, g in enumerate(groups):
                nl = (nb["null_sel_ident"] if g == "selected" else
                      nb["null_ident"] if g.startswith("bin") else None)
                if not nl:
                    continue
                v = np.array([z_[key] for z_ in nl])
                ax.plot([i, i], np.percentile(v, [5, 95]), color=GREY, lw=4, zorder=1,
                        solid_capstyle="round")

        for ax, keys, ylab, title in (
                (axes[0, 1], (("r_same_raw", "o", "same species"), ("r_between_raw", "s", "between")),
                 "Mean r of dataset profiles", "Similarity"),
                (axes[0, 2], (("gap_raw", "o", "profiles"), ("gap_res", "s", "residuals")),
                 "r same - r between species", "Species identifiability"),
                (axes[1, 0], (("sibling_raw", "o", "profiles"), ("sibling_res", "s", "residuals")),
                 f"Sibling dataset nearest (of {nb['ident']['all']['n_sibling']})", "Sibling hits")):
            for key, mk, lab in keys:
                band(ax, key)
                ax.plot(x, [nb["ident"][g][key] for g in groups], mk + "-", color=INK, ms=3,
                        mfc=INK if mk == "o" else "white", label=lab, zorder=3)
            ax.set_xticks(x)
            ax.set_xticklabels(groups, rotation=45, ha="right")
            ax.set_ylabel(ylab)
            ax.legend(fontsize=5)
            ax.set_title(title, loc="left", fontweight="bold")
        ax = axes[1, 1]
        iu = np.triu_indices(len(ctx.sps), 1)
        age = np.asarray(ctx.Dp, float)[iu]
        ua = np.unique(age)
        blues = plt.get_cmap("Blues")
        for g in groups:
            kw = (dict(color="#555555", lw=1.2, ls="--") if g == "all" else
                  dict(color="#C0392B", lw=1.2) if g == "selected" else
                  dict(color=blues(0.35 + 0.6 * (int(g[3:]) - 1) / max(B - 1, 1)), lw=0.8))
            y = nb["ident"][g]["Dsp"][iu]
            ax.plot(ua, [y[age == a].mean() for a in ua], label=g, **kw)
        ax.set_xlabel("Divergence time (Mya)")
        ax.set_ylabel("Species distance (1 - r), mean per age")
        ax.legend(fontsize=4.5, ncol=2)
        ax.set_title("Distance vs divergence", loc="left", fontweight="bold")
        axes[1, 2].axis("off")
        _save(fig, base)
    return base


def supp_bins_profiles(ctx, ph, nb, out_dir):
    groups = ["all", "selected"] + [f"bin{b + 1}" for b in range(ctx.B)]
    P = {"all": ph["P_all"], "selected": ph["P_sel"],
         **{f"bin{b + 1}": nb["P_bins"][b] for b in range(ctx.B)}}
    order = sorted(ctx.sets, key=lambda d: (ctx.sps.index(ctx.sp_of_ds[d]), ctx.labels[d]))
    base = os.path.join(out_dir, "Supp_Figure4_bins_profiles")
    n = len(groups)
    with _style():
        fig, axes = plt.subplots(2, n, figsize=(183 * MM, 95 * MM),
                                 gridspec_kw=dict(wspace=0.08, hspace=0.25))
        edges = [i + 0.5 for i in range(len(order) - 1)
                 if ctx.sp_of_ds[order[i]] != ctx.sp_of_ds[order[i + 1]]]
        for j, g in enumerate(groups):
            M = np.vstack([P[g][d] for d in order])
            im0 = axes[0, j].imshow(M, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto",
                                    interpolation="nearest")
            R = np.array([[fs.corr(P[g][a], P[g][b]) for b in order] for a in order])
            im1 = axes[1, j].imshow(R, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
            for e in edges:
                axes[0, j].axhline(e, color="k", lw=0.3)
                axes[1, j].axhline(e, color="k", lw=0.3)
                axes[1, j].axvline(e, color="k", lw=0.3)
            axes[0, j].set_title(g)
            for a in axes[:, j]:
                a.set_xticks([])
                a.set_yticks(range(len(order)) if j == 0 else [])
            if j == 0:
                for a in axes[:, j]:
                    a.set_yticklabels([ctx.labels[d].replace("_", " ") for d in order], fontsize=4.5)
        fig.colorbar(im0, ax=axes[0, :], fraction=0.015, label="profile (z)")
        fig.colorbar(im1, ax=axes[1, :], fraction=0.015, label="r")
        _save(fig, base)
    return base


def supp_surfaces(ctx, ph, nb, out_dir):
    """Mean map of every group, nearest-neighbour resampled onto the oversampled grid of
    the all-runs map (as the pipeline does before rendering), rendered on the inflated
    surface of each species with one colour range per dataset (render_common_scale,
    the Fig. 5a code)."""
    import contextlib
    import io
    import nibabel as nib
    from matplotlib.image import imread
    from nilearn.image import resample_to_img
    from Statistics.Group_fMRI.run_selection import profiles as sri
    from Statistics.Group_fMRI.run_selection import surfaces as srf
    groups = ["all", "selected"] + [f"bin{b + 1}" for b in range(ctx.B)]
    rows = {"all": {d: None for d in ctx.sets}, "selected": ctx.idx_sel,
            **{f"bin{b + 1}": {d: ctx.bins[d][b] for d in ctx.sets} for b in range(ctx.B)}}
    sdir = os.path.join(out_dir, "surfaces")
    mp = {}
    for d in ctx.sets:
        lab = ctx.labels[d]
        target = nib.load(srf._oversampled(ctx.maps[lab]["all"]))
        mp[lab] = {}
        for g in groups:
            img = ctx.dm[d].mean_image(ctx.keys_of(d, rows[g][d]))
            if target.shape[:3] != img.shape[:3]:
                img = resample_to_img(img, target, interpolation="nearest")
            p = os.path.join(sdir, "maps", lab, f"{lab}_{g}_ttest-stat_fisher_cmap_oversampled.nii.gz")
            os.makedirs(os.path.dirname(p), exist_ok=True)
            nib.save(img, p)
            mp[lab][g] = p
    with contextlib.redirect_stdout(io.StringIO()):
        pngs, _ = srf.render_common_scale(sri, sdir, mp, ctx.c["species_of_label"],
                                          pct=ctx.c["cfg"].get("SURFACE_COMMON_PCT", 90.0))
    order = sorted(ctx.sets, key=lambda d: (ctx.sps.index(ctx.sp_of_ds[d]), ctx.labels[d]))
    base = os.path.join(out_dir, "Supp_Figure4_surfaces")
    with _style():
        fig, axes = plt.subplots(len(order), len(groups),
                                 figsize=(183 * MM, 183 * MM * len(order) / len(groups) * 0.62))
        axes = np.atleast_2d(axes)
        for i, d in enumerate(order):
            for j, g in enumerate(groups):
                ax = axes[i, j]
                ax.axis("off")
                p = pngs.get((ctx.labels[d], g))
                if p and os.path.isfile(p):
                    ax.imshow(imread(p))
                if i == 0:
                    ax.set_title(g)
                if j == 0:
                    ax.text(-0.05, 0.5, ctx.labels[d].replace("_", " "), transform=ax.transAxes,
                            ha="right", va="center", fontsize=5)
        _save(fig, base)
    return base


def _try(fn, *a, **k):
    try:
        return fn(*a, **k)
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"  [{fn.__name__}] [WARN] panel failed: {e}")
        return None
