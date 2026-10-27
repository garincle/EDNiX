"""
Figure 4 -- retrosplenial seed-based maps on the cortical surface of every dataset:
selected runs (centre) and excluded runs (sides), one colour scale per dataset (5a).

Drawn by the Fig 5 pipeline (fig5_paper.fig5a_surfaces, renders of part 1) into
<FIG5_ROOT>/06_figures/Fig5/5a_surface_maps_selected_vs_excluded; this script copies
it as the paper figure (single panel, no letter).
Output: <PAPER_DIR>/Figure_4/Figure_4.pdf|png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Figure_4', 'Figure_4')

if __name__ == "__main__":
    ps.compose_images([pdat.fig5_panel('5a_surface_maps_selected_vs_excluded')], OUT,
                      letters=" ", gutter_mm=0)
