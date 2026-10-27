"""
Supplementary Figure 7 -- behaviour of the run selection.

  a  convergence of the iterative selection: kept runs unchanged from the previous
     round, per dataset                                                  (F1e)
  b  regional profile of every dataset (all runs) and the cross-species consensus,
     regions in the lobe order of Figure 1                               (F1c)
  c  invariance of every region across species; dark = used by the selection score (F5)

Panels drawn at print size by the Fig 5 pipeline (Fig5/Launcher_Fig5.py,
RUN_PAPER_FIGURES=True) into <FIG5_ROOT>/06_figures/; this script lays them out.
Output: <PAPER_DIR>/Supp_Fig_07/Supp_Fig_07.pdf|png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_07', 'Supp_Fig_07')

if __name__ == "__main__":
    ps.compose_images([pdat.fig5_panel('F1e_convergence'), pdat.fig5_panel('F1c_consensus_by_dataset'), pdat.fig5_panel('F5_region_invariance')], OUT, layout=[[0, 1], [2]])
