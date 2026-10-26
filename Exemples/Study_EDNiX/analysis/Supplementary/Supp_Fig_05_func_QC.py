"""
Supplementary Figure 5 -- resting-state fMRI QC, one point per BOLD run.

  a  temporal SNR        b  grey-matter SNR of the mean functional image
  c  mean framewise displacement (Jenkinson, rotation radius = radius of the sphere
     of the subject's brain-mask volume)
  d  mean standardised DVARS     e  global correlation (GCOR)

x = datasets grouped by species (datasets of Figures 4-5).
Output: <PAPER_DIR>/Supp_Fig_05/Supp_Fig_05.pdf|png + values.xlsx
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_05', 'Supp_Fig_05')


def func_metrics(df):
    tsnr = next((c for c in ('func_TSNR', 'func_TSNR_0') if c in df.columns), 'func_TSNR')
    return [(tsnr, 'tSNR'),
            ('func_avg_snr_gray', 'SNR grey matter'),
            ('func_mean_fd', 'Mean FD (mm)'),
            ('func_mean_dvars', 'Mean DVARS (std.)'),
            ('func_gcor', 'GCOR')]


def plot_supp5(df_qc, out_path):
    return ps.plot_by_dataset(df_qc, func_metrics(df_qc), out_path)


if __name__ == "__main__":
    df = pdat.load_qc(pdat.FUNC_DATASETS, 'qc_func')
    plot_supp5(df, OUT)
    cols = ['species', 'bids_dir', 'subject', 'session', 'run'] + [c for c, _ in func_metrics(df)]
    df[[c for c in cols if c in df.columns]].to_excel(OUT + '_values.xlsx', index=False)
