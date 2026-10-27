"""
Supplementary Figure 8 -- selection scores and arousal state.

  a  distribution of the selection score (r of each run with its target) per dataset;
     black = runs kept by the final round, dashed = cut-off              (F3)
  b  correlation of each run with the selection criterion, by arousal state; within-
     animal contrasts for Mouse and Marmoset_Hori                         (6c / F8b)

Panels drawn at print size by the Fig 5 pipeline (Fig5/Launcher_Fig5.py,
RUN_PAPER_FIGURES=True) into <FIG5_ROOT>/06_figures/; this script lays them out.
Output: <PAPER_DIR>/Supp_Fig_08/Supp_Fig_08.pdf|png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_08', 'Supp_Fig_08')

if __name__ == "__main__":
    ps.compose_images([pdat.fig5_panel('F3_score_distributions'), pdat.fig5_panel('6c_arousal_state_F8b')], OUT, layout=[[0], [1]])
