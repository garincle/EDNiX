"""
Figure 2 -- cross-species cortical morphometry at atlas level 2 (lobes).

  a  surface area   b  volume   c  cortical thickness     (log10 axis, bilateral)

Columns = lobes in the order of the Figure 1 sunburst (lobe colour above each
column); x = species in phylogenetic order; violin = species distribution (one
value per subject, sessions averaged); points = subjects, marker = dataset.

Output: <PAPER_DIR>/Figure_2/Figure_2.pdf|png + Figure_2_values.xlsx
(this replaces figures/Figure2/lvl_2_combo_bilateral_log.png of Figure_2_anat.py;
every other morphometry figure is in sensitivity_analysis/morphometry/).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

LEVEL = 2
OUT = os.path.join(pdat.PAPER_DIR, 'Figure_2', 'Figure_2')


def plot_figure2(morpho, regions, out_path, letters="abc"):
    """morpho: {'surface'|'volume'|'thickness': DataFrame(species, bids_dir, subject,
    region, <value>)} already bilateral, one row per subject, at one atlas level."""
    regions = ps.region_order(regions)
    rows = [(m, c, l) for m, c, l in pdat.MORPHO_METRICS
            if morpho.get(m) is not None and not morpho[m].empty]
    species = ps.phylo_sort(pd.concat([morpho[m]['species'] for m, _, _ in rows]).unique())
    datasets = pd.concat([morpho[m][['bids_dir', 'species']] for m, _, _ in rows]).drop_duplicates()
    sp_of = dict(zip(datasets['bids_dir'], datasets['species']))
    rng = np.random.default_rng(42)

    with plt.rc_context(ps.PAPER_RC):
        fig_h = 33 * ps.MM * len(rows) + 14 * ps.MM
        fig, axes = plt.subplots(len(rows), len(regions), figsize=(ps.WIDTH_2COL, fig_h),
                                 sharey="row", squeeze=False)
        fig.subplots_adjust(left=0.065, right=0.995, top=0.955, bottom=0.17,
                            wspace=0.06, hspace=0.12)
        for i, (mod, col, label) in enumerate(rows):
            df = morpho[mod]
            for j, reg in enumerate(regions):
                ax = axes[i, j]
                rdf = df[df['region'] == reg]
                for x, sp in enumerate(species):
                    sdf = rdf[rdf['species'] == sp]
                    vals = sdf[col].to_numpy(float)
                    vals = vals[vals > 0]
                    if vals.size == 0:
                        continue
                    ps.violin(ax, np.log10(vals), x, ps.species_color(sp))
                    for bd, bdf in sdf.groupby('bids_dir'):
                        v = bdf[col].to_numpy(float)
                        ps.dots(ax, np.log10(v[v > 0]), x, ps.species_color(sp),
                                marker=ps.dataset_marker(bd, sp), rng=rng)
                ax.set_xlim(-0.6, len(species) - 0.4)
                ax.set_xticks(range(len(species)))
                if i == len(rows) - 1:
                    ax.set_xticklabels(species, rotation=45, ha="right",
                                       rotation_mode="anchor")
                else:
                    ax.set_xticklabels([])
                if j == 0:
                    ax.set_ylabel(label)
                if i == 0:
                    ps.lobe_strip(ax, reg, height=0.05)
                    ax.text(0.5, 1.07, ps.region_label(reg), transform=ax.transAxes,
                            ha="center", va="bottom")
            # log10 axis written as powers of ten
            ax0 = axes[i, 0]
            lo, hi = ax0.get_ylim()
            ticks = np.arange(np.floor(lo), np.ceil(hi) + 1)
            minor = np.log10(np.concatenate([np.arange(2, 10) * 10.0 ** t for t in ticks]))
            if hi - lo >= 1.5:      # several decades: powers of ten
                ax0.set_yticks(ticks)
                ax0.set_yticklabels([f"$10^{{{int(t)}}}$" for t in ticks])
                ax0.set_yticks(minor[(minor > lo) & (minor < hi)], minor=True)
            else:                   # less than ~1.5 decade (thickness): 1-2-5 in plain numbers
                vals = np.concatenate([np.array([1, 2, 5]) * 10.0 ** t for t in ticks])
                vals = vals[(np.log10(vals) > lo) & (np.log10(vals) < hi)]
                ax0.set_yticks(np.log10(vals))
                ax0.set_yticklabels([f"{v:g}" for v in vals])
                ax0.set_yticks(minor[(minor > lo) & (minor < hi)], minor=True)
            ax0.set_ylim(lo, hi)
            fig.canvas.draw()
            ps.row_label(fig, axes[i, 0], letters[i])

        h, l = ps.species_dataset_legend(list(sp_of), sp_of)
        fig.legend(h, l, loc="lower center", ncol=min(len(h), 7), bbox_to_anchor=(0.53, 0.0))
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


if __name__ == "__main__":
    regions = pdat.CORTICAL_REGIONS[LEVEL]
    morpho = pdat.load_morphometry(level=LEVEL)
    morpho = {m: pdat.select_regions(df, regions) if not df.empty else df
              for m, df in morpho.items()}
    plot_figure2(morpho, regions, OUT)
    figure2_table(morpho, regions).to_excel(OUT + '_values.xlsx', index=False)
    print(f"  [table] {OUT}_values.xlsx")
