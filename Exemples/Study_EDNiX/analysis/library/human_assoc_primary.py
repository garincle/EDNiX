"""
Supplementary Figure 7 -- association vs primary cortex in the human ageing dataset
(BIDS_Park): same analysis and drawing as Figure 3 g-i (macaque).

  a  surface area   b  cortical thickness   c  volume     (order of Figure 3 columns)

Per subject: % change of the summed association (PFC, posterior parietal) and summed
primary (striate, auditory, somatosensory, motor/premotor) values from the first
session; thin lines = subjects, thick = mean +/- SEM where >= 3 subjects are observed;
inset = paired t-test on the last-session % change (association - primary).

  load_supp6_data()   inputs (BIDS_Park morphometry and ages)
  supp6_stats()       every number (paper_stats.supp7)
  draw_supp6()        the drawing (paper_figures.supp7)

Former Supplementary/Supp_Fig_06_human_assoc_primary.py: plot_supp6 split into
supp6_stats + draw_supp6, same computations in the same order.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis.library import longitudinal as F3

# x axis of the trajectories. BIDS_Park is an accelerated longitudinal design (baseline
# ages ~20-90 y, a few years of follow-up): averaging "% change from baseline" at a
# given AGE mixes subjects followed for different durations, and the mean curve breaks
# into pieces. "years_from_baseline" puts every subject at 0 on its first session, so
# the mean is a mean over the same follow-up time. "age" = drawing of Figure 3 g-i.
X_AXIS = "years_from_baseline"


def _from_baseline(df):
    df = df.copy()
    df["age_years"] = df["age_years"] - df.groupby("subject")["age_years"].transform("min")
    return df


def load_supp6_data():
    """(morphometry, age table) of BIDS_Park, as Figure 3 loads them."""
    atlas_paths, _ = F3._build_atlas_paths(['Human'])
    morph = F3.load_morphometry({'Human': F3.HUMAN_BIDS_DIRS}, atlas_paths,
                                regions_of_interest=F3.REGIONS)
    age = F3.format_human_age(F3.HUMAN_BIDS_DIRS, age_prefix=F3.HUMAN_AGE_PREFIX)
    return morph, age


def supp6_stats(morph, age_df, species="Human", letters="abc", x_axis=X_AXIS):
    """One entry per panel (modality, Figure 3 column order): trajectories of
    association and primary % change (res) and the paired t-test on the last-session
    % change (test). Returns dict(panels, lines (Supp_Fig_07_stats.txt), table)."""
    panels, lines, rows = [], [], []
    for k, (mod, _, _, title) in enumerate(F3.FIG3_MODALITIES):
        res = F3.assoc_primary_trajectories(morph, age_df, F3.PRIMARY_REGIONS,
                                            F3.ASSOCIATION_REGIONS, species_label=species,
                                            data_key=mod)
        if res is None:
            panels.append(dict(letter=letters[k], modality=mod, title=title, res=None,
                               test=None))
            lines.append(f"{letters[k]}  {title}: no data")
            rows.append(dict(panel=letters[k], modality=mod))
            continue
        if x_axis == "years_from_baseline":
            res = (_from_baseline(res[0]), _from_baseline(res[1]), res[2])
        st = F3._assoc_primary_test(res[0], res[1])
        panels.append(dict(letter=letters[k], modality=mod, title=title, res=res, test=st))
        lines.append(f"{letters[k]}  {title}: n={st['n']} subjects; last-session % change "
                     f"association {st['assoc_pct']:.2f}% vs primary {st['primary_pct']:.2f}%; "
                     f"paired t({max(st['n'] - 1, 0)})={st['t']:.2f}, {F3._fmt_p(st['p'])}, "
                     f"dz={st['dz']:.2f}")
        rows.append(dict(panel=letters[k], modality=mod, n_subjects=st['n'],
                         assoc_pct=st['assoc_pct'], primary_pct=st['primary_pct'],
                         t=st['t'], p=st['p'], dz=st['dz']))
    return dict(panels=panels, lines=lines, table=pd.DataFrame(rows), x_axis=x_axis)


def stats_text(stats):
    return ("Supplementary Figure 7 -- association vs primary, human (BIDS_Park)\n"
            + "\n".join(stats["lines"]) + "\n")


def draw_supp6(stats, out_path):
    """Drawing of supp6_stats; same layout as the former plot_supp6."""
    x_axis = stats["x_axis"]
    with plt.rc_context(ps.PAPER_RC):
        fig, axes = plt.subplots(1, len(stats["panels"]), figsize=(ps.WIDTH_2COL, 48 * ps.MM))
        fig.subplots_adjust(left=0.07, right=0.995, top=0.93, bottom=0.2, wspace=0.35)
        for k, (ax, P) in enumerate(zip(axes, stats["panels"])):
            ps.panel_label(ax, P["letter"], x=-0.2, y=1.0)
            if P["res"] is None:
                ax.text(0.5, 0.5, 'regions not found', ha='center', va='center',
                        transform=ax.transAxes, color='0.5')
                continue
            F3.draw_assoc_vs_primary(ax, *P["res"], legend_fontsize=5.5,
                                     age_label="Years from baseline"
                                     if x_axis == "years_from_baseline" else "Age (years)")
            for ln in ax.lines:
                ln.set_linewidth(min(ln.get_linewidth(), 1.2))
            ax.set_ylabel(f"% {P['title'].lower()} change")
            if k:
                ax.get_legend().remove()
            st = P["test"]
            if np.isfinite(st['p']):
                ax.text(0.97, 0.03, f"paired {F3._fmt_p(st['p'])}", transform=ax.transAxes,
                        ha='right', va='bottom')
        paths = ps.save(fig, out_path)
        plt.close(fig)
    return paths


def plot_supp6(morph, age_df, out_path, species="Human", letters="abc", x_axis=X_AXIS):
    """supp6_stats + draw_supp6 + <out_path>_stats.txt (former entry point)."""
    stats = supp6_stats(morph, age_df, species=species, letters=letters, x_axis=x_axis)
    paths = draw_supp6(stats, out_path)
    with open(out_path + '_stats.txt', 'w') as fh:
        fh.write(stats_text(stats))
    return paths
