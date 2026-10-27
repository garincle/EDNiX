"""
supp_panels.py -- print-size versions of the Fig 5 pipeline figures that the
paper uses as panels of Figures 4-5 and Supplementary Figures 7-11.

Same data and statistics as the exploration figures of profile_stats /
fingerprints (F1e, F1c, F3, F5, F9, F10, F16, F20), redrawn at their printed
width (Nature: 89 / 120 / 183 mm, 5-7 pt), without titles (the legend of the
article says what is shown; numbers stay in the .txt of each exploration figure),
in the colour code of the article (analysis/paper_style.py: one colour per
species, marker / line style per dataset, lobe colours of Figure 1) and with the
regions in the lobe order of Figure 1 (Atlas_nomenclature.atlas_order).

Every function writes <out>.png (600 dpi) and <out>.pdf (TrueType) through
fa.save_paper. Called by paper_figures.py (panels of the paper, from paper_stats.py) and
by sensitivity_analysis/fmri_selection/exploration_pipeline._run_figures (06_figures/Paper/).
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import Exemples.Study_EDNiX.analysis.fmri_fingerprints.profile_stats as fa
from Exemples.Study_EDNiX.analysis import paper_style as ps

MM = fa.MM


def _lab(label):
    return str(label).replace("_", " ")


def _sp(label, species_of_label=None):
    if species_of_label and label in species_of_label:
        return species_of_label[label]
    return ps.dataset_info(label)[0] or str(label).split("_")[0]


def _tcol(sp):
    """Species colour, darkened for text when it is too light to read."""
    from matplotlib.colors import to_rgb
    rgb = np.array(to_rgb(ps.species_color(sp)))
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return tuple(rgb * 0.62) if lum > 0.62 else tuple(rgb)


def _labels_sorted(labels, species_of_label=None):
    return ps.dataset_sort(list(labels), {l: _sp(l, species_of_label) for l in labels})


def _base(name):
    """'prefrontal_cortex_(pfc)|q80' -> ('prefrontal_cortex_(pfc)', 'q80')."""
    b, sep, ch = str(name).rpartition("|")
    return (b, ch) if sep else (str(name), "")


def atlas_entry_order(names):
    """Indices of profile entries in atlas order (Figure 1 lobes), all 'high'
    channels then all 'low' channels (two-channel profiles)."""
    from Exemples.Study_EDNiX.atlases import Atlas_nomenclature as AN
    block = [fa._channel_block(n) for n in names]
    keys = [AN.atlas_key(_base(n)[0]) for n in names]
    return sorted(range(len(names)), key=lambda i: (block[i], keys[i], i))


def short_entry(name):
    import Exemples.Study_EDNiX.analysis.fmri_fingerprints.fingerprints as ff
    return ff.short_region(name)


def _lobe_strip_x(ax, names, lw=2.5):
    ps_regions = [_base(n)[0] for n in names]
    import matplotlib.transforms as mt
    tr = mt.offset_copy(mt.blended_transform_factory(ax.transData, ax.transAxes),
                        fig=ax.figure, y=-(lw / 2 + 0.5), units="points")
    for j, r in enumerate(ps_regions):
        ax.plot([j - 0.5, j + 0.5], [0, 0], transform=tr, color=ps.lobe_color(r), lw=lw,
                solid_capstyle="butt", clip_on=False)
    ax.tick_params(axis="x", pad=lw + 2.5)


def _dataset_legend(fig_or_ax, labels, species_of_label=None, lines=False, **kw):
    h = []
    for l in _labels_sorted(labels, species_of_label):
        sp = _sp(l, species_of_label)
        m, ls = fa.dataset_style(l)
        h.append(Line2D([], [], color=ps.species_color(sp), ls=ls if lines else "none",
                        marker=None if lines else m, mec="black", mew=0.2, ms=3.5, lw=0.9,
                        label=_lab(l)))
    return fig_or_ax.legend(handles=h, frameon=False, **kw)


# ── F1e: convergence of the iterative selection ──────────────────────────────
def f1e_convergence(overlap_df, out, width_mm=60, height_mm=48):
    if overlap_df is None or overlap_df.empty:
        return None
    df = overlap_df.copy()
    df["frac_overlap"] = df["n_overlap_with_previous"] / df["n_select"]
    with fa.paper_style():
        fig, ax = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
        for lab in _labels_sorted(df["label"].unique()):
            g = df[df["label"] == lab].sort_values("round")
            m, ls = fa.dataset_style(lab)
            ax.plot(g["round"], g["frac_overlap"], ls=ls, marker=m, ms=2.5, mec="black",
                    mew=0.2, lw=0.8, color=ps.species_color(_sp(lab)), label=_lab(lab))
        ax.set_xlabel("Selection round")
        ax.set_ylabel("Kept runs unchanged\nfrom previous round")
        ax.set_ylim(-0.03, 1.03)
        ax.set_xticks(sorted(df["round"].unique()))
        fa.save_paper(fig, out)
    return out


# ── F1c: consensus by dataset ────────────────────────────────────────────────
def f1c_values(species_profiles, common, species_of_label, value_col, zscore):
    """Values drawn on F1c (Supplementary Figure 7b): {dataset: scaled profile over
    `common`} and the consensus (median across species of the per-species median).
    (None, None) when no dataset has a value."""
    per, by_sp = {}, {}
    for lab, prof in species_profiles.items():
        v = fa._ordered_common_values(prof, common, value_col)
        if not np.isfinite(v).any():
            continue
        v = fa._normalize_values(v, zscore)
        per[lab] = v
        by_sp.setdefault(species_of_label.get(lab, lab), []).append(v)
    if not per:
        return None, None
    with np.errstate(invalid="ignore"):
        cons = np.nanmedian(np.vstack([np.nanmedian(np.vstack(v), axis=0)
                                       for v in by_sp.values()]), axis=0)
    return per, cons


def f1c_consensus_by_dataset(species_profiles, common, species_of_label, value_col, zscore,
                             out, width_mm=95, height_mm=48):
    per, cons = f1c_values(species_profiles, common, species_of_label, value_col, zscore)
    if not per:
        return None
    order = atlas_entry_order(common)
    names = [common[i] for i in order]
    x = np.arange(len(order))
    with fa.paper_style():
        fig, ax = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
        for lab in _labels_sorted(per, species_of_label):
            _, ls = fa.dataset_style(lab)
            ax.plot(x, per[lab][order], ls=ls, lw=0.7, alpha=0.95,
                    color=ps.species_color(species_of_label.get(lab)), label=_lab(lab))
        ax.plot(x, cons[order], color="black", lw=1.4, zorder=5, label="consensus")
        ax.axhline(0, color="#999999", lw=0.4)
        nb = sum(1 for n in names if fa._channel_block(n) == 0)
        if 0 < nb < len(names):
            ax.axvline(nb - 0.5, color="#666666", lw=0.4, ls=":")
        ax.set_xticks(x)
        ax.set_xticklabels([short_entry(n) for n in names], rotation=90, fontsize=5)
        ax.set_xlim(-0.6, len(x) - 0.4)
        nl = (fa._norm_label(zscore) or "").strip()
        ax.set_ylabel(f"{nl} profile value".strip().replace("RMS-scaled", "RMS-scaled")
                      if nl else "Profile value")
        _lobe_strip_x(ax, names)
        ax.legend(frameon=False, fontsize=5, ncol=1, loc="upper left",
                  bbox_to_anchor=(1.0, 1.02), handlelength=2.2)
        fa.save_paper(fig, out)
    return out


# ── F5: region invariance ────────────────────────────────────────────────────
def f5_region_invariance(invariance, common, mask, out, width_mm=183, height_mm=38):
    inv = np.asarray(invariance, float)
    order = atlas_entry_order(common)
    names = [common[i] for i in order]
    used = (np.asarray(mask, bool)[order] if mask is not None else np.ones(len(order), bool))
    with fa.paper_style():
        fig, ax = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
        ax.bar(range(len(order)), inv[order], width=0.75,
               color=["#4D4D4D" if u else "#CCCCCC" for u in used], lw=0)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([short_entry(n) for n in names], rotation=90, fontsize=5)
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_ylabel("Invariance across species\n(1 \u2212 rescaled MAD)")
        _lobe_strip_x(ax, names)
        ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#4D4D4D", label="used by the selection"),
                           plt.Rectangle((0, 0), 1, 1, color="#CCCCCC", label="not used")],
                  frameon=False, loc="upper right", ncol=2, fontsize=5)
        fa.save_paper(fig, out)
    return out


# ── F3: score distributions ──────────────────────────────────────────────────
def f3_score_distributions(rankings, selections, out, width_mm=183, ncols=5, panel_h_mm=24):
    labels = _labels_sorted(rankings)
    if not labels:
        return None
    nrows = int(np.ceil(len(labels) / ncols))
    with fa.paper_style():
        fig, axes = plt.subplots(nrows, ncols, figsize=(width_mm * MM, (panel_h_mm * nrows + 6) * MM),
                                 squeeze=False)
        for i, lab in enumerate(labels):
            ax = axes[i // ncols][i % ncols]
            v = pd.to_numeric(rankings[lab]["pearson_r"], errors="coerce").dropna().to_numpy()
            if v.size == 0:
                ax.axis("off")
                continue
            sp = _sp(lab)
            n_sel = len(selections.get(lab, []))
            bins = np.linspace(v.min(), v.max(), 25)
            ax.hist(v, bins=bins, color=ps.species_color(sp), alpha=0.85, lw=0)
            if n_sel and n_sel <= v.size:
                cut = float(np.sort(v)[::-1][n_sel - 1])
                ax.hist(v[v >= cut], bins=bins, color="black", alpha=0.85, lw=0)
                ax.axvline(cut, color="black", lw=0.6, ls="--")
            ax.text(0.02, 0.97, _lab(lab), transform=ax.transAxes, ha="left", va="top",
                    fontsize=5.5, color=_tcol(sp))
            if i // ncols == nrows - 1 or i + ncols >= len(labels):
                ax.set_xlabel("Score (r with target)")
            if i % ncols == 0:
                ax.set_ylabel("Runs")
        for j in range(len(labels), nrows * ncols):
            axes[j // ncols][j % ncols].axis("off")
        fig.tight_layout(pad=0.3, w_pad=0.6, h_pad=0.4)
        fa.save_paper(fig, out)
    return out


# ── F9: run-to-run consistency ───────────────────────────────────────────────
def f9_run_distance_sd(dist_df, out, value_col="pearson_r", width_mm=70, height_mm=50):
    if dist_df is None or dist_df.empty:
        return None
    s = (dist_df.groupby(["label", "species"], as_index=False)[value_col]
         .agg(sd="std", n="count"))
    labels = _labels_sorted(s["label"], dict(zip(s["label"], s["species"])))
    s = s.set_index("label").loc[labels].reset_index()
    with fa.paper_style():
        fig, ax = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
        x = np.arange(len(s))
        ax.bar(x, s["sd"], color=[ps.species_color(sp) for sp in s["species"]], width=0.75,
               edgecolor="black", lw=0.3)
        for i, l in enumerate(s["label"]):
            if fa.dataset_style(l)[0] != "o":
                ax.bar(i, s["sd"].iloc[i], width=0.75, facecolor="none", edgecolor="white",
                       hatch="////", lw=0)
        ax.set_xticks(x)
        ax.set_xticklabels([_lab(l) for l in s["label"]], rotation=55, ha="right",
                           rotation_mode="anchor", fontsize=5)
        for t, sp in zip(ax.get_xticklabels(), s["species"]):
            t.set_color(_tcol(sp))
        ax.set_ylabel("SD of run r with reference")
        fa.save_paper(fig, out)
    return out


# ── F10: QC metrics vs distance (meta-analysis across datasets) ─────────────
QC_NAMES = {
    "func_TSNR": "tSNR", "func_TSNR_0": "tSNR", "func_avg_snr_gray": "SNR GM",
    "func_mean_fd": "mean FD", "func_mean_dvars": "DVARS", "func_gcor": "GCOR",
    "func_avg_snr_white": "SNR WM", "func_cnr": "CNR", "func_fwhm": "FWHM",
    "func_avg_outcount": "outlier count", "func_avg_enorm": "motion (enorm)",
}


def f10_qc_vs_distance(meta_df, out, alpha=0.05, max_metrics=20, width_mm=113, height_mm=50):
    if meta_df is None or meta_df.empty:
        return None
    pcol = "meta_q_fdr" if "meta_q_fdr" in meta_df.columns else "meta_p"
    sub = meta_df.reindex(meta_df["meta_r"].abs().sort_values(ascending=False).index)
    sub = sub.head(max_metrics).reset_index(drop=True)
    sig = (sub[pcol] < alpha).to_numpy()
    with fa.paper_style():
        fig, ax = plt.subplots(figsize=(width_mm * MM, height_mm * MM))
        x = np.arange(len(sub))
        ax.bar(x, sub["meta_r"], width=0.75, lw=0,
               color=["#4D4D4D" if g else "#CCCCCC" for g in sig])
        if "meta_ci_lo" in sub.columns and "meta_ci_hi" in sub.columns:
            ax.vlines(x, sub["meta_ci_lo"], sub["meta_ci_hi"], color="black", lw=0.5)
        for i, (r, g) in enumerate(zip(sub["meta_r"], sig)):
            if g:
                ax.text(i, r + (0.01 if r >= 0 else -0.01), fa.stars(sub[pcol].iloc[i]),
                        ha="center", va="bottom" if r >= 0 else "top", fontsize=5)
        ax.axhline(0, color="#888888", lw=0.4)
        ax.set_xticks(x)
        ax.set_xticklabels([QC_NAMES.get(m, str(m).replace("func_", "").replace("_", " "))
                            for m in sub["metric"]], rotation=55, ha="right",
                           rotation_mode="anchor", fontsize=5)
        ax.set_ylabel("r with distance to reference\n(meta-analysis across datasets)")
        fa.save_paper(fig, out)
    return out


# ── F16 a: dataset x dataset distance matrix (Figure 4c) ─────────────────────
def f16a_proximity_matrix(D, species_of_label, out, width_mm=63, annotate=True):
    labels = _labels_sorted(list(D.index), species_of_label)
    M = D.loc[labels, labels].to_numpy(float).copy()
    np.fill_diagonal(M, np.nan)
    n = len(labels)
    with fa.paper_style():
        fig = plt.figure(figsize=(width_mm * MM, width_mm * MM * 0.98))
        ax = fig.add_axes([0.30, 0.30, 0.56, 0.56])
        cax = fig.add_axes([0.88, 0.30, 0.025, 0.56])
        cmap = plt.get_cmap("viridis_r").with_extremes(bad="white")
        im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=np.nanmax(M),
                       interpolation="none")
        if annotate:
            for i in range(n):
                for j in range(n):
                    if i != j and np.isfinite(M[i, j]):
                        ax.text(j, i, f"{M[i, j]:.2f}".lstrip("0"), ha="center", va="center",
                                fontsize=3.6,
                                color="white" if M[i, j] > 0.55 * np.nanmax(M) else "black")
        sp = [species_of_label.get(l) for l in labels]
        for i in range(1, n):
            if sp[i] != sp[i - 1]:
                ax.axhline(i - 0.5, color="white", lw=0.8)
                ax.axvline(i - 0.5, color="white", lw=0.8)
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels([_lab(l) for l in labels], rotation=55, ha="right",
                           rotation_mode="anchor", fontsize=4.8)
        ax.set_yticklabels([_lab(l) for l in labels], fontsize=4.8)
        for t, s in zip(list(ax.get_xticklabels()) + list(ax.get_yticklabels()), sp + sp):
            t.set_color(_tcol(s))
        ax.tick_params(length=0, pad=1)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("1 - r", fontsize=5.5)
        cb.ax.tick_params(labelsize=5, length=1.5)
        cb.outline.set_linewidth(0.3)
        fa.save_paper(fig, out)
    return out


# ── F20: selection vs variability levels ─────────────────────────────────────
def f20_variability_levels(table, draws, species_of_label, out, width_mm=183, ncols=5,
                           panel_h_mm=34):
    from Exemples.Study_EDNiX.analysis.fmri_fingerprints.fingerprints import LEVELS
    labs = _labels_sorted(table["label"], species_of_label)
    if not labs:
        return None
    ncols = min(ncols, len(labs))
    nrows = int(np.ceil(len(labs) / ncols))
    x = np.arange(len(LEVELS))
    with fa.paper_style():
        fig, axes = plt.subplots(nrows, ncols, sharey=True, squeeze=False,
                                 figsize=(width_mm * MM, (panel_h_mm * nrows + 8) * MM))
        for k, lab in enumerate(labs):
            ax = axes[k // ncols, k % ncols]
            d = draws[lab]
            sp = species_of_label.get(lab)
            lo = [np.nanpercentile(d["null"][lv], 2.5) for lv in LEVELS]
            hi = [np.nanpercentile(d["null"][lv], 97.5) for lv in LEVELS]
            ax.vlines(x, lo, hi, color="#CCCCCC", lw=3.5, zorder=1)
            ax.plot(x, [d["obs_all"][lv] for lv in LEVELS], color="#555555", marker="s", ms=2,
                    lw=0.6, ls="--", zorder=3)
            ax.plot(x, [d["obs_sel"][lv] for lv in LEVELS], color=ps.species_color(sp),
                    marker="o", ms=3, mec="black", mew=0.2, lw=1.0, zorder=4)
            row = table[table["label"] == lab].iloc[0]
            ax.text(0.03, 0.03, _lab(lab), transform=ax.transAxes, fontsize=5.5,
                    color=_tcol(sp), va="bottom")
            p = row.get("p_gap_animal", np.nan)
            st = fa.stars(p)
            if st:
                ax.text(0.97, 0.97, st, transform=ax.transAxes, ha="right", va="top",
                        fontsize=6)
            ax.set_xticks(x)
            ax.set_xticklabels(["run", "session", "animal"], fontsize=5.5)
            if k % ncols == 0:
                ax.set_ylabel("Median run-to-run r")
        for k in range(len(labs), nrows * ncols):
            axes[k // ncols, k % ncols].axis("off")
        h = [Line2D([], [], color="#555555", marker="s", ms=2, ls="--", lw=0.6, label="all runs"),
             Line2D([], [], color="#444444", marker="o", ms=3, lw=1.0, label="selected"),
             Line2D([], [], color="#CCCCCC", lw=3.5, label="random subsets of the same size (95%)")]
        fig.legend(handles=h, loc="lower center", ncol=3, frameon=False,
                   bbox_to_anchor=(0.5, 0.0))
        fig.tight_layout(rect=(0, 0.07, 1, 1), pad=0.3, w_pad=0.4, h_pad=0.6)
        fa.save_paper(fig, out)
    return out


# ── F16 b: proximity of each dataset to every other one (Supp. Fig. 9) ─────
def f16b_proximity_rows(D, dist_df, species_of_label, test, out, width_mm=120, row_mm=5.2):
    """One row per dataset: dots = 1 - r to every other dataset (species colour of
    that dataset), star = the other dataset of the same species, black tick = the
    leave-one-species-out consensus, ring = nearest neighbour. Panel a of the
    exploration F16 is Figure 4c; it is not repeated here."""
    labels = _labels_sorted(list(D.index), species_of_label)
    n = len(labels)
    per = test["per_dataset"].set_index("label") if test else None
    cons = dist_df.set_index("label")["d_consensus"] if "d_consensus" in dist_df else None
    with fa.paper_style():
        fig, ax = plt.subplots(figsize=(width_mm * MM, (row_mm * n + 14) * MM))
        for i, lab in enumerate(labels):
            y = n - 1 - i
            for other in labels:
                if other == lab:
                    continue
                same = species_of_label.get(other) == species_of_label.get(lab)
                ax.scatter(D.loc[lab, other], y, s=18 if same else 6,
                           marker="*" if same else "o",
                           color=ps.species_color(species_of_label.get(other)),
                           edgecolor="black" if same else "none", lw=0.3,
                           zorder=3 if same else 2)
            if cons is not None and lab in cons.index and np.isfinite(cons[lab]):
                ax.scatter(cons[lab], y, marker="|", s=40, color="black", lw=1.0, zorder=4)
            if per is not None and lab in per.index:
                nn = per.loc[lab, "nearest_neighbour"]
                ax.scatter(D.loc[lab, nn], y, s=22, facecolor="none", edgecolor="#222222",
                           lw=0.5, zorder=5)
        ax.set_yticks(range(n))
        ax.set_yticklabels([_lab(l) for l in labels[::-1]])
        for t, l in zip(ax.get_yticklabels(), labels[::-1]):
            t.set_color(_tcol(species_of_label.get(l)))
        ax.set_ylim(-0.6, n - 0.4)
        ax.set_xlim(0, np.nanmax(D.to_numpy(float)) * 1.05)
        ax.set_xlabel("1 - r to the dataset of the row")
        ax.grid(axis="x", color="#EEEEEE", lw=0.4)
        h = [Line2D([], [], marker="*", color="#888888", mec="black", mew=0.3, ls="none", ms=5,
                    label="same species, other dataset"),
             Line2D([], [], marker="|", color="black", ls="none", ms=6, mew=1.0,
                    label="consensus of the other species"),
             Line2D([], [], marker="o", mfc="none", mec="#222222", mew=0.5, ls="none", ms=5,
                    label="nearest neighbour")]
        ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3,
                  frameon=False)
        if test:
            ax.text(1.0, 1.01, fa.stars(test.get("p_one_sided")), transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=6)
        fa.save_paper(fig, out)
    return out
