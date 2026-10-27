"""
Figure 5 -- species identity of the retrosplenial connectivity fingerprint.

  a  each dataset's fingerprint (coloured) over the Human datasets (black dashed),
     most similar first; key radar = spoke names, ring = lobes of Figure 1  (F15b)
  b  run-to-run profile similarity, same species / other dataset vs different
     species, all runs and selected runs                                      (5b)
  c  dataset x dataset distance (1 - r) of the group profiles                 (F16 a)

Comes after Figure 4 (the seed maps the fingerprints are computed from).
Panels are drawn by the Fig 5 pipeline at print size (Fig5/Launcher_Fig5.py,
RUN_PAPER_FIGURES=True) into <FIG5_ROOT>/06_figures/; this script lays them out.
Numbers and tests: 06_figures/Fig5/5b_intra_vs_inter_species.txt,
supplementary/F15_fingerprint_vs_reference.txt, F16_dataset_proximity.txt.
Output: <PAPER_DIR>/Figure_5/Figure_5.pdf|png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Figure_5', 'Figure_5')

if __name__ == "__main__":
    panels = [pdat.fig5_panel('F15b_fingerprint_vs_human'),
              pdat.fig5_panel('5b_intra_vs_inter_species'),
              pdat.fig5_panel('F16a_proximity_matrix')]
    ps.compose_images(panels, OUT, layout=[[0], [1, 2]])
