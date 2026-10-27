"""
Figure 2 -- cross-species cortical morphometry at atlas level 2 (lobes).

  a  surface area   b  volume   c  cortical thickness     (log10 axis, bilateral)

Columns = lobes in the order of the Figure 1 sunburst (lobe colour above each
column); x = datasets, grouped by species in phylogenetic order (tick label in the
species colour); grey violin = distribution of the dataset (one value per subject,
sessions averaged); points = subjects, species colour, one marker per dataset
(paper_style.dataset_marker). Clément, 2026-10-27: x = datasets, not species.

figure2_table: the statistics (paper_stats.figure2); plot_figure2: the drawing
(paper_figures.figure2). Former analysis/Figure_2_morphometry.py, functions unchanged.
Every other morphometry figure is in sensitivity_analysis/morphometry/.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

LEVEL = 2


ROW_HEIGHT_MM = 42      # per modality row: panels about as tall as wide (as lvl_2_combo)


def _log_ticks(ax):
    """log10 axis written as powers of ten (several decades) or 1-2-5 (less than ~1.5)."""
    lo, hi = ax.get_ylim()
    ticks = np.arange(np.floor(lo), np.ceil(hi) + 1)
    minor = np.log10(np.concatenate([np.arange(2, 10) * 10.0 ** t for t in ticks]))
    if hi - lo >= 1.5:
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"$10^{{{int(t)}}}$" for t in ticks])
    else:
        vals = np.concatenate([np.array([1, 2, 5]) * 10.0 ** t for t in ticks])
        vals = vals[(np.log10(vals) > lo) & (np.log10(vals) < hi)]
        ax.set_yticks(np.log10(vals))
        ax.set_yticklabels([f"{v:g}" for v in vals])
    ax.set_yticks(minor[(minor > lo) & (minor < hi)], minor=True)
    ax.set_ylim(lo, hi)


def plot_figure2(morpho, regions, out_path, letters="abc"):
    """morpho: {'surface'|'volume'|'thickness': DataFrame(species, bids_dir, subject,
    region, <value>)} already bilateral, one row per subject, at one atlas level.
    x = datasets grouped by species; grey violin = dataset; points = subjects in the
    species colour with the dataset marker."""
    regions = ps.region_order(regions)
    rows = [(m, c, l) for m, c, l in pdat.MORPHO_METRICS
            if morpho.get(m) is not None and not morpho[m].empty]
    datasets = pd.concat([morpho[m][['bids_dir', 'species']] for m, _, _ in rows]).drop_duplicates()
    sp_of = dict(zip(datasets['bids_dir'], datasets['species']))
    dsets = ps.dataset_sort(list(sp_of), sp_of)
    rng = np.random.default_rng(42)

    with plt.rc_context(ps.PAPER_RC):
        fig_h = (ROW_HEIGHT_MM * len(rows) + 34) * ps.MM
        fig, axes = plt.subplots(len(rows), len(regions), figsize=(ps.WIDTH_2COL, fig_h),
                                 sharey="row", squeeze=False)
        fig.subplots_adjust(left=0.075, right=0.995, top=0.95, bottom=0.14,
                            wspace=0.08, hspace=0.1)
        for i, (mod, col, label) in enumerate(rows):
            df = morpho[mod]
            for j, reg in enumerate(regions):
                ax = axes[i, j]
                rdf = df[df['region'] == reg]
                for x, bd in enumerate(dsets):
                    sp = sp_of[bd]
                    v = rdf.loc[rdf['bids_dir'] == bd, col].to_numpy(float)
                    v = np.log10(v[v > 0])
                    if v.size == 0:
                        continue
                    ps.violin(ax, v, x, ps.VIOLIN_GREY, width=0.8, alpha=0.45)
                    ps.dots(ax, v, x, ps.dataset_color(bd, sp), marker=ps.dataset_marker(bd, sp),
                            spread=0.16, size=3.0, rng=rng)
                ax.set_xlim(-0.6, len(dsets) - 0.4)
                ax.set_xticks(range(len(dsets)))
                if i == len(rows) - 1:
                    ax.set_xticklabels([ps.dataset_short(d, sp_of[d]) for d in dsets], rotation=55,
                                       ha="right", rotation_mode="anchor", fontsize=5)
                    for t, d in zip(ax.get_xticklabels(), dsets):
                        t.set_color(ps.text_color(sp_of[d]))
                else:
                    ax.set_xticklabels([])
                if j == 0:
                    ax.set_ylabel(label)
                if i == 0:
                    ps.lobe_strip(ax, reg, height=0.04)
                    ax.text(0.5, 1.06, ps.region_label(reg), transform=ax.transAxes,
                            ha="center", va="bottom")
            _log_ticks(axes[i, 0])
        fig.canvas.draw()
        # panel letter level with the top of each row, just left of its y-axis label
        r = fig.canvas.get_renderer()
        for i in range(len(rows)):
            ax0 = axes[i, 0]
            x_lab = ax0.yaxis.label.get_window_extent(r).x0 / fig.bbox.width
            fig.text(x_lab - 1.0 * ps.MM / fig.get_figwidth(), ax0.get_position().y1, letters[i],
                     fontsize=ps.PANEL_LABEL_SIZE, fontweight="bold", va="top", ha="right")
        h, l = ps.species_dataset_legend(dsets, sp_of)
        ps.legend_below(fig, axes, h, l, x=0.54)
        paths = ps.save(fig, out_path)
        plt.close(fig)
    return paths


def figure2_table(morpho, regions):
    """Values behind the figure: per-species median, IQR and n subjects."""
    out = []
    for mod, col, _ in pdat.MORPHO_METRICS:
        df = morpho.get(mod)
        if df is None or df.empty:
            continue
        g = df[df['region'].isin(regions)].groupby(['region', 'species'])[col]
        t = g.agg(n='count', median='median',
                  q25=lambda v: v.quantile(.25), q75=lambda v: v.quantile(.75)).reset_index()
        t.insert(0, 'modality', mod)
        out.append(t)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()
