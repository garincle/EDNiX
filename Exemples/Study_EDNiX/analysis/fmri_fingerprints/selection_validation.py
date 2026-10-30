"""
selection_validation.py -- Figure 4b-d: what the run selection buys, against keeping the
same number of runs at random (Clément, 2026-10-30).

Every comparison is selected runs vs 10 runs drawn at random per dataset, so that the
number of runs per dataset (25-195 before selection) never drives the result.

  b  generalisation (non-circular): the regional values are split at random into two
     halves; runs are selected on one half (leave-one-species-out consensus, top runs per
     dataset) and evaluated on the other half, which played no part in their choice.
     Three measures, selected (red) vs random draws (grey):
       cross-dataset r   median r between runs of different datasets
       r with other sp.  median over datasets of r(group profile, LOSO consensus)
       dataset-specific  G = median r same dataset - median r across datasets
  c  each dataset: r(group profile, LOSO consensus) on held-out values, random median ->
     selected, one line per dataset (species colour, dataset marker)
  d  noise removed, identity kept (summary of the former 4b histograms): median and IQR of
     the run-to-run r for the three pair classes, all runs / random / selected.

The selection in b-c is a simplified criterion (one round, no invariance weighting); d uses
the selection of the paper. The script versions of these tests are
analysis/fig4b_selection_improvement.py, fig4_selection_benefit.py, fig4_selection_heldout.py.
"""
import numpy as np
import pandas as pd

N_SPLITS = 50
N_DRAWS = 200
N_DRAWS_PAIRS = 1000
SEED = 0
CLASSES = ("same\ndataset", "same species,\nother dataset", "different\nspecies")
POOLS = ("all runs", "random", "selected")


# ── data ─────────────────────────────────────────────────────────────────────
def run_vectors(c):
    """Run profiles of the saved selection (same objects as Fig. 4): vecs (runs x values,
    normalised), ds (dataset index), selected (bool), sp_of_ds, labels."""
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints import profile_stats as fa
    cfg = c["cfg"]
    value_col, norm = cfg["VALUE_COL"], cfg["ZSCORE"]
    common, sol = c["common"], c["species_of_label"]
    labels = sorted(c["profiles_by_label"])
    vecs, ds, selected = [], [], []
    for i, lab in enumerate(labels):
        kept = set(c["selected_units"].get(lab, {}))
        for key, prof in c["profiles_by_label"][lab].items():
            v = fa._ordered_common_values(prof, common, value_col)
            if not np.isfinite(v).any():
                continue
            vecs.append(fa._normalize_values(v, norm))
            ds.append(i)
            selected.append(key in kept)
    return (np.vstack(vecs), np.array(ds), np.array(selected),
            np.array([sol.get(l, l) for l in labels]), labels)


def _z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def _corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() >= 3 else np.nan


def _corr_rows(m):
    return pd.DataFrame(m.T).corr(min_periods=3).to_numpy()


def _loso_targets(vecs, ds, sp_of_ds):
    """Per dataset: median over the OTHER species of the species medians of the dataset
    means (all runs)."""
    sets = np.unique(ds)
    means = {e: np.nanmean(vecs[ds == e], axis=0) for e in sets}
    spm = {s: np.nanmedian([means[e] for e in sets if sp_of_ds[e] == s], axis=0)
           for s in set(sp_of_ds[sets])}
    return {d: _z(np.nanmedian([spm[s] for s in spm if s != sp_of_ds[d]], axis=0))
            for d in sets}


def _pair_classes(r, ds, sp_of_ds, idx):
    d = ds[idx]
    iu = np.triu_indices(len(idx), k=1)
    rr = r[np.ix_(idx, idx)][iu]
    da, db = d[iu[0]], d[iu[1]]
    ok = np.isfinite(rr)
    rr, da, db = rr[ok], da[ok], db[ok]
    same_ds, same_sp = da == db, sp_of_ds[da] == sp_of_ds[db]
    return rr, same_ds, same_sp


def _class_q(r, ds, sp_of_ds, idx):
    """3 pair classes (same dataset / same species other dataset / different species) x
    (q25, median, q75) of the run-to-run r of the runs `idx`."""
    rr, same_ds, same_sp = _pair_classes(r, ds, sp_of_ds, idx)
    masks = (same_ds, same_sp & ~same_ds, ~same_sp)
    return np.array([np.percentile(rr[m], [25, 50, 75]) if m.any() else [np.nan] * 3
                     for m in masks])


def _heldout_eval(r, vecs, ds, sp_of_ds, idx, targets):
    rr, same_ds, _ = _pair_classes(r, ds, sp_of_ds, idx)
    e1 = float(np.median(rr[~same_ds]))
    g = float(np.median(rr[same_ds]) - e1)
    per = {e: _corr(_z(np.nanmean(vecs[idx[ds[idx] == e]], axis=0)), targets[e])
           for e in np.unique(ds[idx])}
    return e1, float(np.median(list(per.values()))), g, per


def heldout(vecs, ds, sp_of_ds, n_splits=N_SPLITS, n_draws=N_DRAWS, n_keep=10, seed=SEED):
    """Panel b-c numbers. Returns sel (E1, E2, G averaged over splits), rand (draws x 3,
    averaged over splits), p (one-sided, b), per_dataset DataFrame (c)."""
    rng = np.random.default_rng(seed)
    nreg = vecs.shape[1]
    per_ds = {d: np.flatnonzero(ds == d) for d in np.unique(ds)}
    sel_s, rand_s, sel_per, rand_per = [], [], [], []
    q_all, q_sel, q_rand = [], [], []
    for _ in range(n_splits):
        perm = rng.permutation(nreg)
        a, b = perm[: nreg // 2], perm[nreg // 2:]
        va, vb = vecs[:, a], vecs[:, b]
        ta, tb = _loso_targets(va, ds, sp_of_ds), _loso_targets(vb, ds, sp_of_ds)
        keep = []
        for d, idx in per_ds.items():
            sc = np.array([_corr(va[i], ta[d]) for i in idx])
            sc[~np.isfinite(sc)] = -np.inf
            keep.append(idx[np.argsort(-sc)[: min(n_keep, len(idx))]])
        rb = _corr_rows(vb)
        e1, e2, g, per = _heldout_eval(rb, vb, ds, sp_of_ds, np.concatenate(keep), tb)
        sel_s.append((e1, e2, g))
        sel_per.append(per)
        q_all.append(_class_q(rb, ds, sp_of_ds, np.arange(len(ds))))
        q_sel.append(_class_q(rb, ds, sp_of_ds, np.concatenate(keep)))
        draws, dper, dq = [], [], []
        for _ in range(n_draws):
            idx = np.concatenate([rng.choice(i, min(n_keep, len(i)), replace=False)
                                  for i in per_ds.values()])
            x1, x2, xg, p_ = _heldout_eval(rb, vb, ds, sp_of_ds, idx, tb)
            draws.append((x1, x2, xg))
            dper.append(p_)
            dq.append(_class_q(rb, ds, sp_of_ds, idx))
        rand_s.append(draws)
        q_rand.append(dq)
        rand_per.append({d: float(np.median([p_[d] for p_ in dper])) for d in per_ds})
    sel = np.mean(sel_s, axis=0)
    rand = np.mean(np.array(rand_s), axis=0)
    p = np.array([(1 + np.sum(rand[:, 0] >= sel[0])) / (n_draws + 1),
                  (1 + np.sum(rand[:, 1] >= sel[1])) / (n_draws + 1),
                  (1 + np.sum(rand[:, 2] <= sel[2])) / (n_draws + 1)])
    # d: pair classes on the test regions (averaged over splits); P = share of random draws
    # whose class median is at least that of the selection (one-sided)
    qa, qs = np.mean(q_all, axis=0), np.mean(q_sel, axis=0)
    qr = np.mean(np.array(q_rand), axis=0)                  # draws x 3 classes x 3 quantiles
    p_cls = np.array([(1 + np.sum(qr[:, k, 1] >= qs[k, 1])) / (n_draws + 1) for k in range(3)])
    rows = []
    for pool, arr in (("all runs", qa), ("random", np.median(qr, axis=0)), ("selected", qs)):
        for k, c_ in enumerate(CLASSES):
            rows.append(dict(pool=pool, pair_class=c_.replace("\n", " "), q25=arr[k, 0],
                             median=arr[k, 1], q75=arr[k, 2]))
    levels = pd.DataFrame(rows)
    per_dataset = pd.DataFrame(dict(
        dataset=list(per_ds),
        species=[sp_of_ds[d] for d in per_ds],
        r_random=[float(np.mean([x[d] for x in rand_per])) for d in per_ds],
        r_selected=[float(np.mean([x[d] for x in sel_per])) for d in per_ds]))
    # c: paired one-sided Wilcoxon signed-rank across datasets, selected vs random median
    from scipy import stats as _st
    diff = (per_dataset.r_selected - per_dataset.r_random).to_numpy()
    w = _st.wilcoxon(diff, alternative="greater") if np.any(diff != 0) else None
    p_c = float(w.pvalue) if w is not None else np.nan
    return dict(sel=sel, rand=rand, p=p, per_dataset=per_dataset, n_splits=n_splits,
                n_draws=n_draws, n_values=nreg, p_c=p_c, n_improved=int((diff > 0).sum()),
                levels=levels, p_classes=p_cls)


def pair_levels(vecs, ds, selected, sp_of_ds, n_draws=N_DRAWS_PAIRS, seed=SEED):
    """Panel d numbers: median and quartiles of the run-to-run r for the three pair
    classes, all runs / random subsets (median over draws of each quantile) / selected."""
    rng = np.random.default_rng(seed + 1)
    r = _corr_rows(vecs)
    per_ds = {d: np.flatnonzero(ds == d) for d in np.unique(ds)}
    n_sel = {d: int(selected[i].sum()) for d, i in per_ds.items()}

    def q(idx):
        rr, same_ds, same_sp = _pair_classes(r, ds, sp_of_ds, idx)
        masks = (same_ds, same_sp & ~same_ds, ~same_sp)
        return np.array([np.percentile(rr[m], [25, 50, 75]) if m.any() else [np.nan] * 3
                         for m in masks])                      # 3 classes x 3 quantiles

    out = {"all runs": q(np.arange(len(ds))), "selected": q(np.flatnonzero(selected))}
    draws = [q(np.concatenate([rng.choice(per_ds[d], n_sel[d], replace=False)
                               for d in per_ds if n_sel[d]])) for _ in range(n_draws)]
    out["random"] = np.nanmedian(np.array(draws), axis=0)
    rows = [dict(pool=p, pair_class=c.replace("\n", " "), q25=out[p][k, 0],
                 median=out[p][k, 1], q75=out[p][k, 2])
            for p in POOLS for k, c in enumerate(CLASSES)]
    return pd.DataFrame(rows)


def compute(c, n_splits=N_SPLITS, n_draws=N_DRAWS):
    vecs, ds, selected, sp_of_ds, labels = run_vectors(c)
    h = heldout(vecs, ds, sp_of_ds, n_splits=n_splits, n_draws=n_draws)
    h["per_dataset"].insert(0, "label", [labels[d] for d in h["per_dataset"]["dataset"]])
    return dict(heldout=h, levels=h["levels"])


MEASURES = ("between\ndatasets", "with other\nspecies", "dataset-\nspecific (G)")


def result_lines(v):
    h = v["heldout"]
    lines = [f"b  held-out values: {h['n_splits']} random splits of the {h['n_values']} "
             f"regional values; selection on one half (simplified criterion: one round, LOSO "
             f"consensus, 10 runs per dataset), evaluation on the other; {h['n_draws']} random "
             "draws of 10 runs per dataset per split"]
    for k, name in enumerate(MEASURES):
        lines.append(f"  {name.replace(chr(10), ' ')}: selected {h['sel'][k]:+.3f}, random "
                     f"median {np.median(h['rand'][:, k]):+.3f}; one-sided P = {h['p'][k]:.4g}"
                     + (" (minimum possible)" if h['p'][k] <= 1 / (h['n_draws'] + 1) + 1e-12
                        else ""))
    pdt = h["per_dataset"]
    up = int((pdt.r_selected > pdt.r_random).sum())
    lines.append(f"c  r with the LOSO consensus on test regions, per dataset: improved in "
                 f"{up}/{len(pdt)}; median {pdt.r_random.median():.3f} (random) -> "
                 f"{pdt.r_selected.median():.3f} (selected); paired Wilcoxon signed-rank "
                 f"one-sided P = {h['p_c']:.4g} (n = {len(pdt)} datasets)")
    lines += [f"  {r.label}: {r.r_random:.3f} -> {r.r_selected:.3f}" for r in pdt.itertuples()]
    lines.append("d  run-to-run r on test regions, median [IQR] (all runs / random 10 per "
                 "dataset / selected); P = share of random draws with a class median >= the "
                 "selection")
    for k, cl in enumerate(CLASSES):
        t = v["levels"][v["levels"].pair_class == cl.replace("\n", " ")]
        lines.append(f"  {cl.replace(chr(10), ' ')}: " + " / ".join(
            f"{r.median:.3f} [{r.q25:.3f}, {r.q75:.3f}]" for r in t.itertuples())
            + f"; P = {h['p_classes'][k]:.4g}")
    return lines


# ── drawing ──────────────────────────────────────────────────────────────────
SEL_COLOR = "#C0392B"
POOL_COLORS = {"all runs": "#BDBDBD", "random": "#7F7F7F", "selected": SEL_COLOR}


def _letter(fig, ax, letter, size, dx_mm=7.0):
    pos = ax.get_position()
    fig.text(pos.x0 - dx_mm / (fig.get_figwidth() * 25.4), 1.0, letter, fontsize=size,
             fontweight="bold", va="top", ha="left")


B_SHOWN = (0, 1)                      # b draws the two similarity measures; G stays in the tables
B_LABELS = ("run to run,\ndifferent datasets", "dataset profile vs\nother species")


def draw_b(ax, h, fs=5.5):
    """Held-out test. For each measure: grey violin = 10 random runs per dataset (one value
    per draw), red diamond = selected runs; both computed on regional values that played
    no part in the selection."""
    rand, sel, p = h["rand"], h["sel"], h["p"]
    pmin = 1 / (h["n_draws"] + 1)
    for x, k in enumerate(B_SHOWN):
        vp = ax.violinplot(rand[:, k], positions=[x], widths=0.6, showextrema=False)
        for b_ in vp["bodies"]:
            b_.set_facecolor("#BDBDBD")
            b_.set_edgecolor("none")
            b_.set_alpha(1)
        ax.plot(x, sel[k], marker="D", ms=3.6, color=SEL_COLOR, mec="white", mew=0.4, zorder=4)
        lab = f"P ≤ {pmin:.1g}" if p[k] <= pmin + 1e-12 else f"P = {p[k]:.2g}"
        ax.text(x, max(sel[k], np.max(rand[:, k])) + 0.04, lab, ha="center", va="bottom",
                fontsize=fs - 0.5)
    ax.set_xticks(range(len(B_SHOWN)))
    ax.set_xticklabels(B_LABELS, fontsize=fs - 0.5)
    ax.set_xlim(-0.55, len(B_SHOWN) - 0.45)
    ax.set_ylabel("r on test regions", fontsize=fs)
    ax.tick_params(labelsize=fs - 0.5, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend([Patch(fc="#BDBDBD", ec="none"),
               Line2D([], [], ls="", marker="D", ms=3.6, color=SEL_COLOR, mec="white")],
              ["random 10 runs", "selected 10 runs"], loc="lower right", frameon=False,
              fontsize=fs - 0.5, handlelength=1.0, handletextpad=0.4, borderaxespad=0.1)


def draw_b_schematic(ax, n=30, fs=5.5):
    """Strip of the regional values: dark = used to select the runs, white = used to test
    them (one of the random splits)."""
    rng = np.random.default_rng(3)
    half = set(rng.permutation(n)[: n // 2])
    for i in range(n):
        ax.add_patch(plt_rect(i, 0, 0.86, 1, "#4D4D4D" if i in half else "white"))
    ax.set_xlim(-0.5, n + 0.5)
    ax.set_ylim(-2.4, 1.2)
    ax.add_patch(plt_rect(0, -2.2, 0.86, 1, "#4D4D4D"))
    ax.text(1.4, -1.7, "select", fontsize=fs - 0.5, va="center", ha="left")
    ax.add_patch(plt_rect(n / 2, -2.2, 0.86, 1, "white"))
    ax.text(n / 2 + 1.4, -1.7, "test", fontsize=fs - 0.5, va="center", ha="left")
    ax.text(n / 2, 1.5, "regional values (random half split, x50)", fontsize=fs - 0.5,
            ha="center", va="bottom")
    ax.axis("off")


def plt_rect(x, y, w, h, fc):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, h, fc=fc, ec="#4D4D4D", lw=0.3)


def _pfmt(p, floor=None):
    if floor is not None and p <= floor + 1e-12:
        return f"P ≤ {floor:.1g}"
    return f"P = {p:.2g}" if p >= 0.01 else f"P = {p:.1g}"


def draw_c(ax, h, sol_color, marker_of, fs=5.5):
    """One line per dataset, random median -> selected (held-out values)."""
    pdt = h["per_dataset"]
    for r in pdt.itertuples():
        col = sol_color(r.species)
        ax.plot([0, 1], [r.r_random, r.r_selected], color=col, lw=0.8, alpha=0.9)
        for x, y in ((0, r.r_random), (1, r.r_selected)):
            ax.plot(x, y, marker=marker_of(r.label, r.species), ms=3.0, color=col,
                    mec="white", mew=0.3, zorder=3)
    ax.text(0.5, 1.0, f"paired Wilcoxon, n = {len(pdt)}: " + _pfmt(h["p_c"]),
            transform=ax.transAxes, ha="center", va="bottom", fontsize=fs - 0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["random\n10 runs", "selected\n10 runs"], fontsize=fs - 0.5)
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel("dataset profile vs other species,\nr on test regions", fontsize=fs)
    ax.tick_params(labelsize=fs - 0.5, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)


def draw_d(ax, levels, p_classes=None, n_draws=N_DRAWS, fs=5.5):
    """Median and IQR of run-to-run r per pair class, for the three pools."""
    off = {"all runs": -0.24, "random": 0.0, "selected": 0.24}
    for k, cl in enumerate(CLASSES):
        t = levels[levels.pair_class == cl.replace("\n", " ")].set_index("pool")
        for pool in POOLS:
            if pool not in t.index:
                continue
            x = k + off[pool]
            ax.plot([x, x], [t.loc[pool, "q25"], t.loc[pool, "q75"]], color=POOL_COLORS[pool],
                    lw=1.2, solid_capstyle="butt")
            ax.plot(x, t.loc[pool, "median"], marker="o", ms=3.0, color=POOL_COLORS[pool],
                    mec="white", mew=0.3, zorder=3)
    if p_classes is not None:
        for k, cl in enumerate(CLASSES):
            t = levels[levels.pair_class == cl.replace("\n", " ")].set_index("pool")
            ax.text(k + off["selected"], t.loc["selected", "q75"] + 0.03,
                    _pfmt(p_classes[k], 1 / (n_draws + 1)), ha="center", va="bottom",
                    fontsize=fs - 1.0, color=SEL_COLOR)
    ax.set_ylim(ax.get_ylim()[0], levels["q75"].max() + 0.22)
    ax.set_xticks(range(3))
    ax.set_xticklabels(CLASSES, fontsize=fs - 0.5)
    ax.set_xlim(-0.6, 2.6)
    ax.set_ylabel("run-to-run r on test regions", fontsize=fs)
    ax.tick_params(labelsize=fs - 0.5, length=2)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    from matplotlib.lines import Line2D
    ax.legend([Line2D([], [], color=POOL_COLORS[p], marker="o", ms=3, lw=1.2) for p in POOLS],
              ["all runs", "random 10", "selected 10"], loc="lower left",
              bbox_to_anchor=(0.0, 1.0), ncol=3, frameon=False, fontsize=fs - 0.5,
              handlelength=1.0, columnspacing=0.8, handletextpad=0.3, borderaxespad=0.1)
