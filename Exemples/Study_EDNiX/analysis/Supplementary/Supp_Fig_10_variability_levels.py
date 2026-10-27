"""
Supplementary Figure 10 -- does the selection erase inter-run, inter-session or
inter-animal variability? (F20) Median run-to-run profile r per level, all runs
(grey dashed), selected runs (species colour), random subsets of the same size (95 %).

Panels drawn at print size by the Fig 5 pipeline (Fig5/Launcher_Fig5.py,
RUN_PAPER_FIGURES=True) into <FIG5_ROOT>/06_figures/; this script lays them out.
Output: <PAPER_DIR>/Supp_Fig_10/Supp_Fig_10.pdf|png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_10', 'Supp_Fig_10')

if __name__ == "__main__":
    ps.compose_images([pdat.fig5_panel('F20_variability_levels')], OUT, layout=[[0]])
