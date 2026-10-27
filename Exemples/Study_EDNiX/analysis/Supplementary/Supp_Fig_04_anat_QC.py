"""
Supplementary Figure 4 -- anatomical QC, one point per (subject, session).

  a  grey-matter SNR              b  subcortical SNR (atlas level 1 'Subcortical areas',
                                     mean of left and right)
  c  GM/WM CNR                    d  correlation of the individual image with its
                                     species template (template warped to native space)
  e  normalised mutual information between the image and the warped template

x = datasets grouped by species (same datasets as Figure 2), one colour per dataset.
CNR: 'anat_cnr' (|mu_GM - mu_WM| / sqrt(sd_GM^2 + sd_WM^2), MRIQC definition) is written
by EDNiX anat QC since PR #34; subjects processed before it only have
'anat_cortical_contrast' (|mu_WM - mu_GM| / mean, a relative contrast, NOT a CNR) --
in that case the panel is drawn with that metric and the axis says so.
NMI: 'anat_nmi_brain' (inside the brain, written since this version) when present,
otherwise 'anat_nmi' computed on the whole field of view (background included: it then
depends on how much of the field of view the brain fills, which differs by species).
Output: <PAPER_DIR>/Supp_Fig_04/Supp_Fig_04.pdf|png + values.xlsx
"""
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_style as ps
from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_04', 'Supp_Fig_04')
SUBCORTICAL = 'Subcortical areas'      # level-1 label of the QC atlas (anat QC regional_snr)


def add_subcortical_snr(df, prefix='anat'):
    """<prefix>_snr_subcortical = mean of the left and right SNR of the 'Subcortical
    areas' label (columns <prefix>_regional_snr_Subcortical areas_Left/_Right)."""
    cols = [f'{prefix}_regional_snr_{SUBCORTICAL}_{h}' for h in ('Left', 'Right')]
    cols = [c for c in cols if c in df.columns]
    if cols:
        df = df.copy()
        df[f'{prefix}_snr_subcortical'] = (df[cols].apply(pd.to_numeric, errors='coerce')
                                           .mean(axis=1, skipna=True))
    return df


def _first(df, options):
    for col, label in options:
        if col in df.columns and df[col].notna().any():
            return col, label
    return options[-1]


def anat_metrics(df):
    return [('anat_avg_snr_gray', 'SNR grey matter'),
            ('anat_snr_subcortical', 'SNR subcortical'),
            _first(df, [('anat_cnr', 'GM/WM CNR'),
                        ('anat_cortical_contrast', 'GM/WM contrast (rel.)')]),
            ('anat_template_correlation', 'Correlation to template (r)'),
            _first(df, [('anat_nmi_brain', 'NMI to template'),
                        ('anat_nmi', 'NMI to template (whole FOV)')])]


def plot_supp4(df_qc, out_path):
    df_qc = add_subcortical_snr(df_qc)
    return ps.plot_by_dataset(df_qc, anat_metrics(df_qc), out_path)


if __name__ == "__main__":
    df = add_subcortical_snr(pdat.load_qc(pdat.ANAT_DATASETS, 'qc_anat'))
    plot_supp4(df, OUT)
    cols = ['species', 'bids_dir', 'subject', 'session'] + [c for c, _ in anat_metrics(df)]
    df[[c for c in cols if c in df.columns]].to_excel(OUT + '_values.xlsx', index=False)
