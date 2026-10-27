"""
Supplementary Figure 4 -- anatomical QC, one point per (subject, session).

  a  grey-matter SNR      b  white-matter SNR      c  GM/WM CNR
  d  correlation of the individual image with its species template
  e  smoothness (FWHM, AFNI 3dFWHMx)

x = datasets grouped by species (same datasets as Figure 2).
CNR: 'anat_cnr' (|mu_GM - mu_WM| / sqrt(sd_GM^2 + sd_WM^2), MRIQC definition) is written
by EDNiX anat QC since this version; subjects processed before it only have
'anat_cortical_contrast' (|mu_WM - mu_GM| / mean, a relative contrast, NOT a CNR) --
in that case the panel is drawn with that metric and the axis says so.
Output: <PAPER_DIR>/Supp_Fig_04/Supp_Fig_04.pdf|png + values.xlsx
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_04', 'Supp_Fig_04')


def anat_metrics(df):
    cnr = (('anat_cnr', 'GM/WM CNR') if 'anat_cnr' in df.columns and df['anat_cnr'].notna().any()
           else ('anat_cortical_contrast', 'GM/WM contrast (rel.)'))
    return [('anat_avg_snr_gray', 'SNR grey matter'),
            ('anat_avg_snr_white', 'SNR white matter'),
            cnr,
            ('anat_template_correlation', 'Correlation to template (r)'),
            ('anat_fwhm', 'FWHM (mm)')]


def plot_supp4(df_qc, out_path):
    return ps.plot_by_dataset(df_qc, anat_metrics(df_qc), out_path)


if __name__ == "__main__":
    df = pdat.load_qc(pdat.ANAT_DATASETS, 'qc_anat')
    plot_supp4(df, OUT)
    cols = ['species', 'bids_dir', 'subject', 'session'] + [c for c, _ in anat_metrics(df)]
    df[[c for c in cols if c in df.columns]].to_excel(OUT + '_values.xlsx', index=False)
