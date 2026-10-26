"""
Supplementary Figures 6-8 -- regional outlier rates, a proxy of segmentation /
surface-reconstruction failures.   6 surface area   7 volume   8 cortical thickness

Atlas level 3 (cortical areas: same parcellation as Fig. 3 association/primary and
Figs 4-5), one value per subject (sessions averaged).

What is an outlier here (OUTLIER_MODEL = "allometric_residual", default)
  Inside each dataset, log10(region value) is regressed on log10(the subject's total
  value over the cortical regions of this level) (ordinary least squares); a region of
  a subject is an outlier when its residual deviates from the dataset median by more
  than 2.5 scaled MAD (|r - median| / (1.4826 MAD) > 2.5; Leys et al. 2013).
  Why: pooling the datasets of a species and testing raw values (previous version)
  counts as "outliers" what is biology or protocol -- infant vs adult macaques
  (Zhu_Garin), 20-90-year-old humans (Park), scanner differences between datasets.
  A region that is wrong relative to the rest of the same brain is what a failed
  segmentation produces; a whole brain that is small is not.
  OUTLIER_MODEL = "raw_pooled" reproduces the previous computation (raw values, all
  datasets of a species pooled) and goes to sensitivity_analysis/.

Figure: heat map, rows = datasets grouped by species (n subjects in the label),
columns = regions in the lobe order of Figure 1 (lobe colour above), cell = % of
subjects flagged. Grey = fewer than MIN_N subjects.
Output: <PAPER_DIR>/Supp_Fig_0{6,7,8}/... + outlier_rates_lvl3.xlsx (all cells, and
the per-modality summary used in the text).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

LEVEL = 3
THRESHOLD = 2.5
MIN_N = 5
OUTLIER_MODEL = "allometric_residual"     # or "raw_pooled"
SUPP_NUMBER = {"surface": 6, "volume": 7, "thickness": 8}


def outlier_table(morpho, regions, model=OUTLIER_MODEL, threshold=THRESHOLD, min_n=MIN_N):
    """One row per (modality, species, dataset, region): n, n_outliers, outlier_rate."""
    rows = []
    for mod, col, _ in pdat.MORPHO_METRICS:
        df = morpho.get(mod)
        if df is None or df.empty:
            continue
        df = df[df['region'].isin(regions) & (df[col] > 0)].copy()
        df['logv'] = np.log10(df[col])
        if model == "allometric_residual":
            agg = 'mean' if mod == 'thickness' else 'sum'
            tot = (df.groupby(['bids_dir', 'subject'])[col].agg(agg)
                   .rename('total').reset_index())
            df = df.merge(tot, on=['bids_dir', 'subject'])
            df['x'] = np.log10(df['total'])
            groups = df.groupby(['species', 'bids_dir', 'region'])
        else:
            groups = df.groupby(['species', 'region'])
        for key, g in groups:
            sp, reg = key[0], key[-1]
            bd = key[1] if len(key) == 3 else 'pooled'
            if len(g) < min_n:
                rows.append(dict(modality=mod, species=sp, bids_dir=bd, region=reg, n=len(g),
                                 n_outliers=np.nan, outlier_rate=np.nan))
                continue
            y = g['logv'].to_numpy()
            if model == "allometric_residual" and np.ptp(g['x']) > 0:
                b = np.polyfit(g['x'].to_numpy(), y, 1)
                y = y - np.polyval(b, g['x'].to_numpy())
            flag = pdat.robust_outliers(y, threshold)
            rows.append(dict(modality=mod, species=sp, bids_dir=bd, region=reg, n=len(g),
                             n_outliers=int(flag.sum()), outlier_rate=float(flag.mean())))
    return pd.DataFrame(rows)


def summary(tab):
    """Per modality: subject-weighted overall rate, median / max cell rate, % cells > 5 %."""
    t = tab.dropna(subset=['outlier_rate'])
    return (t.groupby('modality').apply(lambda g: pd.Series(dict(
        overall_rate=g['n_outliers'].sum() / g['n'].sum(),
        median_cell_rate=g['outlier_rate'].median(),
        max_cell_rate=g['outlier_rate'].max(),
        pct_cells_above_5pct=(g['outlier_rate'] > 0.05).mean() * 100,
        n_cells=len(g))), include_groups=False).reset_index())


def plot_outlier_heatmap(tab, mod, regions, out_path, vmax=0.20):
    t = tab[tab['modality'] == mod]
    if t.empty:
        return None
    regions = [r for r in ps.region_order(regions) if r in set(t['region'])]
    sp_of = dict(zip(t['bids_dir'], t['species']))
    dsets = ps.dataset_sort(list(sp_of), sp_of)
    M = (t.pivot_table(index='bids_dir', columns='region', values='outlier_rate', aggfunc='first')
         .reindex(index=dsets, columns=regions))
    N = t.groupby('bids_dir')['n'].max().reindex(dsets)
    cmap = plt.get_cmap('Reds').copy()
    cmap.set_bad('0.85')
    with plt.rc_context(ps.PAPER_RC):
        h = (4.2 * len(dsets) + 22) * ps.MM
        fig = plt.figure(figsize=(ps.WIDTH_15COL, h))
        ax = fig.add_axes([0.30, 18 * ps.MM / h, 0.60, 1 - 24 * ps.MM / h])
        im = ax.imshow(np.ma.masked_invalid(M.to_numpy(float)), cmap=cmap, vmin=0, vmax=vmax,
                       aspect='auto', interpolation='none')
        for (i, j), v in np.ndenumerate(M.to_numpy(float)):
            if np.isfinite(v) and v > 0:
                ax.text(j, i, f"{100 * v:.0f}", ha='center', va='center', fontsize=4.5,
                        color='white' if v > vmax * 0.6 else 'black')
        ax.set_xticks(range(len(regions)))
        ax.set_xticklabels([ps.region_label(r) for r in regions], rotation=45, ha='left',
                           rotation_mode='anchor')
        ax.xaxis.tick_top()
        ax.set_yticks(range(len(dsets)))
        ax.set_yticklabels([f"{sp_of[d]} {ps.dataset_short(d, sp_of[d])} (n={int(N[d])})"
                            for d in dsets])
        for tl, d in zip(ax.get_yticklabels(), dsets):
            tl.set_color(ps.species_color(sp_of[d]) if sp_of[d] not in ("Marmoset", "Macaque")
                         else "#9A7D0A")
        for j, r in enumerate(regions):
            ax.add_patch(plt.Rectangle((j - 0.5, -0.5 - 0.35), 1, 0.3, color=ps.lobe_color(r),
                                       clip_on=False, lw=0))
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        ax.tick_params(axis="x", pad=5)
        # separators between species
        sp_seq = [sp_of[d] for d in dsets]
        for i in range(1, len(dsets)):
            if sp_seq[i] != sp_seq[i - 1]:
                ax.axhline(i - 0.5, color='white', lw=1.2)
        cax = fig.add_axes([0.92, 18 * ps.MM / h, 0.015, 0.35])
        cb = fig.colorbar(im, cax=cax, extend='max')
        cb.set_label('Subjects flagged (%)')
        cb.set_ticks(np.linspace(0, vmax, 5))
        cb.set_ticklabels([f"{100 * x:.0f}" for x in np.linspace(0, vmax, 5)])
        cb.outline.set_linewidth(0.4)
        paths = ps.save(fig, out_path)
        plt.close(fig)
    return paths


if __name__ == "__main__":
    regions = pdat.CORTICAL_REGIONS[LEVEL]
    morpho = pdat.load_morphometry(level=LEVEL)
    morpho = {m: pdat.select_regions(df, regions) if not df.empty else df
              for m, df in morpho.items()}
    tab = outlier_table(morpho, regions)
    for mod, k in SUPP_NUMBER.items():
        plot_outlier_heatmap(tab, mod, regions,
                             os.path.join(pdat.PAPER_DIR, f'Supp_Fig_0{k}', f'Supp_Fig_0{k}'))
    out_xlsx = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_06', f'outlier_rates_lvl{LEVEL}.xlsx')
    with pd.ExcelWriter(out_xlsx) as w:
        summary(tab).to_excel(w, sheet_name='summary', index=False)
        tab.to_excel(w, sheet_name='cells', index=False)
    print(summary(tab).to_string(index=False))
    print(f"  [table] {out_xlsx}")
