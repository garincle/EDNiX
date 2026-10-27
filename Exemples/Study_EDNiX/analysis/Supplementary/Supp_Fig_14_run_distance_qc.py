"""
Supplementary Figure 14 -- run-level consistency and data quality.

  a  SD across runs of the correlation of each run with its dataset reference  (F9)
  b  QC metrics vs distance to the reference: meta-analytic r across datasets,
     dark = significant after Benjamini-Hochberg correction               (F10)

Panels drawn at print size by the Fig 5 pipeline (Fig5/Launcher_Fig5.py,
RUN_PAPER_FIGURES=True) into <FIG5_ROOT>/06_figures/; this script lays them out.
Output: <PAPER_DIR>/Supp_Fig_14/Supp_Fig_14.pdf|png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_14', 'Supp_Fig_14')

if __name__ == "__main__":
    ps.compose_images([pdat.fig5_panel('F9_run_distance_sd'), pdat.fig5_panel('F10_qc_vs_distance')], OUT, layout=[[0, 1]])
