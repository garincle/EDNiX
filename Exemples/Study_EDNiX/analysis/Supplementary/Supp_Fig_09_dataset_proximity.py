"""
Supplementary Figure 9 -- proximity of every dataset to every other dataset (F16 b;
F16 a is Figure 5c). Star = the other dataset of the same species, tick = consensus
of the other species, ring = nearest neighbour; stars = sibling permutation test.

Panels drawn at print size by the Fig 5 pipeline (Fig5/Launcher_Fig5.py,
RUN_PAPER_FIGURES=True) into <FIG5_ROOT>/06_figures/; this script lays them out.
Output: <PAPER_DIR>/Supp_Fig_09/Supp_Fig_09.pdf|png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_09', 'Supp_Fig_09')

if __name__ == "__main__":
    ps.compose_images([pdat.fig5_panel('F16b_proximity_rows')], OUT, layout=[[0]])
