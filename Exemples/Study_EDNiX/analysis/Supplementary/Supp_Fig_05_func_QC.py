"""
Supplementary Figure 5 -- resting-state fMRI QC, one point per BOLD run.

  first row (signal and registration of the mean functional image)
  a  grey-matter SNR        b  subcortical SNR (atlas level 1 'Subcortical areas')
  c  cortical contrast      d  normalised mutual information between the functional
                               image and the anatomical image
  second row (time series)
  e  temporal SNR           f  mean framewise displacement (Jenkinson, rotation radius =
                               radius of the sphere of the brain-mask volume)
  g  Jaccard overlap of the functional and anatomical brain masks
  h  mean standardised DVARS                i  global correlation (GCOR)

x = datasets grouped by species (datasets of Figures 4-5), one colour per dataset.
'func_avg_snr_subcortical' is written by EDNiX func QC since this version: runs processed
before it give an empty panel b until the QC step is re-run.
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
    return [('func_avg_snr_gray', 'SNR grey matter'),
            ('func_avg_snr_subcortical', 'SNR subcortical'),
            ('func_cortical_contrast', 'Cortical contrast'),
            ('func_nmi', 'NMI func-anat'),
            None,                               # second row starts here
            (tsnr, 'tSNR'),
            ('func_mean_fd', 'Mean FD (mm)'),
            ('func_jaccard', 'Jaccard func-anat mask'),
            ('func_mean_dvars', 'Mean DVARS (std.)'),
            ('func_gcor', 'GCOR')]


def plot_supp5(df_qc, out_path):
    return ps.plot_by_dataset(df_qc, func_metrics(df_qc), out_path, ncols=5)


if __name__ == "__main__":
    df = pdat.load_qc(pdat.FUNC_DATASETS, 'qc_func')
    plot_supp5(df, OUT)
    cols = ['species', 'bids_dir', 'subject', 'session', 'run'] + [
        m[0] for m in func_metrics(df) if m is not None]
    df[[c for c in cols if c in df.columns]].to_excel(OUT + '_values.xlsx', index=False)
