"""
Supplementary Figure 9 -- association vs primary cortex in the human ageing dataset
(BIDS_Park): same analysis and drawing as Figure 3 g-i (macaque).

  a  surface area   b  cortical thickness   c  volume     (order of Figure 3 columns)

Per subject: % change of the summed association (PFC, posterior parietal) and summed
primary (striate, auditory, somatosensory, motor/premotor) values from the first
session; thin lines = subjects, thick = mean +/- SEM where >= 3 subjects are observed;
inset = paired t-test on the last-session % change (association - primary).
Replaces Figure_3_longit/Human/fig3d/Fig3D_assoc_vs_primary_*_Human.png.
Output: <PAPER_DIR>/Supp_Fig_09/Supp_Fig_09.pdf|png + Supp_Fig_09_stats.txt
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis import Figure_3_longitudinal as F3

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_09', 'Supp_Fig_09')
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


def plot_supp9(morph, age_df, out_path, species="Human", letters="abc", x_axis=X_AXIS):
    lines = []
    with plt.rc_context(ps.PAPER_RC):
        fig, axes = plt.subplots(1, len(F3.FIG3_MODALITIES), figsize=(ps.WIDTH_2COL, 48 * ps.MM))
        fig.subplots_adjust(left=0.07, right=0.995, top=0.93, bottom=0.2, wspace=0.35)
        for k, (ax, (mod, _, _, title)) in enumerate(zip(axes, F3.FIG3_MODALITIES)):
            res = F3.assoc_primary_trajectories(morph, age_df, F3.PRIMARY_REGIONS,
                                                F3.ASSOCIATION_REGIONS, species_label=species,
                                                data_key=mod)
            ps.panel_label(ax, letters[k], x=-0.2, y=1.0)
            if res is None:
                ax.text(0.5, 0.5, 'regions not found', ha='center', va='center',
                        transform=ax.transAxes, color='0.5')
                lines.append(f"{letters[k]}  {title}: no data")
                continue
            if x_axis == "years_from_baseline":
                res = (_from_baseline(res[0]), _from_baseline(res[1]), res[2])
            F3.draw_assoc_vs_primary(ax, *res, legend_fontsize=5.5,
                                     age_label="Years from baseline"
                                     if x_axis == "years_from_baseline" else "Age (years)")
            for ln in ax.lines:
                ln.set_linewidth(min(ln.get_linewidth(), 1.2))
            ax.set_ylabel(f"% {title.lower()} change")
            if k:
                ax.get_legend().remove()
            st = F3._assoc_primary_test(res[0], res[1])
            if np.isfinite(st['p']):
                ax.text(0.97, 0.03, f"paired {F3._fmt_p(st['p'])}", transform=ax.transAxes,
                        ha='right', va='bottom')
            lines.append(f"{letters[k]}  {title}: n={st['n']} subjects; last-session % change "
                         f"association {st['assoc_pct']:.2f}% vs primary {st['primary_pct']:.2f}%; "
                         f"paired t({max(st['n'] - 1, 0)})={st['t']:.2f}, {F3._fmt_p(st['p'])}, "
                         f"dz={st['dz']:.2f}")
        paths = ps.save(fig, out_path)
        plt.close(fig)
    with open(out_path + '_stats.txt', 'w') as fh:
        fh.write("Supplementary Figure 9 -- association vs primary, human (BIDS_Park)\n")
        fh.write("\n".join(lines) + "\n")
    return paths


if __name__ == "__main__":
    atlas_paths, _ = F3._build_atlas_paths(['Human'])
    morph = F3.load_morphometry({'Human': F3.HUMAN_BIDS_DIRS}, atlas_paths,
                                regions_of_interest=F3.REGIONS)
    age = F3.format_human_age(F3.HUMAN_BIDS_DIRS, age_prefix=F3.HUMAN_AGE_PREFIX)
    plot_supp9(morph, age, OUT)
