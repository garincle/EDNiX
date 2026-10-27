"""
QC panels of the paper: Supplementary Figures 4 (anatomical) and 5 (functional).
anat_metrics / func_metrics name the panels, plot_supp4 / plot_supp5 draw them
(paper_figures), qc_summary gives the numbers (paper_stats). Former
Supplementary/Supp_Fig_04_anat_QC.py and Supp_Fig_05_func_QC.py, functions unchanged.

Supplementary Figure 5 -- anatomical QC, one point per (subject, session).

  a  grey-matter SNR              b  subcortical SNR (atlas level 1 'Subcortical areas',
                                     mean of left and right)
  c  GM/WM CNR                    d  correlation of the individual image with its
                                     species template (template warped to native space)
  e  normalised mutual information between the image and the warped template

x = datasets grouped by species (same datasets as Figure 2); colour = species, marker = dataset.
CNR: 'anat_cnr' (|mu_GM - mu_WM| / sqrt(sd_GM^2 + sd_WM^2), MRIQC definition) is written
by EDNiX anat QC since PR #34; subjects processed before it only have
'anat_cortical_contrast' (|mu_WM - mu_GM| / mean, a relative contrast, NOT a CNR) --
in that case the panel is drawn with that metric and the axis says so.
NMI: 'anat_nmi_brain' (inside the brain, written since this version) when present,
otherwise 'anat_nmi' computed on the whole field of view (background included: it then
depends on how much of the field of view the brain fills, which differs by species).

Supplementary Figure 6 -- resting-state fMRI QC, one point per BOLD run.

  2 rows x 4 (GCOR removed, Clément 2026-10-27)
  first row (signal and registration of the mean functional image)
  a  grey-matter SNR        b  subcortical SNR (atlas level 1 'Subcortical areas')
  c  cortical contrast      d  normalised mutual information between the functional
                               image and the anatomical image
  second row (time series)
  e  temporal SNR           f  mean framewise displacement (Jenkinson, rotation radius =
                               radius of the sphere of the brain-mask volume)
  g  Jaccard overlap of the functional and anatomical brain masks
  h  mean standardised DVARS

x = datasets grouped by species (datasets of Figures 4-5); colour = species, marker = dataset.
'func_avg_snr_subcortical' is written by EDNiX func QC (modalities/fMRI/_12_fMRI_QC.py)
only for runs whose QC step ran after it was added: until the QC step is re-run, panel b
shows the white-matter SNR ('func_avg_snr_white') and its axis says so.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from Exemples.Study_EDNiX.analysis import paper_style as ps

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


def func_metrics(df):
    tsnr = next((c for c in ('func_TSNR', 'func_TSNR_0') if c in df.columns), 'func_TSNR')
    return [('func_avg_snr_gray', 'SNR grey matter'),
            _first(df, [('func_avg_snr_subcortical', 'SNR subcortical'),
                        ('func_avg_snr_white', 'SNR white matter')]),
            ('func_cortical_contrast', 'Cortical contrast'),
            ('func_nmi', 'NMI func-anat'),
            (tsnr, 'tSNR'),
            ('func_mean_fd', 'Mean FD (mm)'),
            ('func_jaccard', 'Jaccard func-anat mask'),
            ('func_mean_dvars', 'Mean DVARS (std.)')]


def plot_supp5(df_qc, out_path):
    return ps.plot_by_dataset(df_qc, func_metrics(df_qc), out_path, ncols=4)


def qc_summary(df, metrics, keys=('species', 'bids_dir')):
    """Per dataset and metric: n, median, quartiles, min, max of the plotted values
    (finite values only, as on the figure)."""
    rows = []
    for m in metrics:
        if m is None or m[0] not in df.columns:
            continue
        col, label = m
        v = df[list(keys) + [col]].copy()
        v[col] = pd.to_numeric(v[col], errors='coerce')
        v = v[np.isfinite(v[col])]
        for k, g in v.groupby(list(keys), sort=False):
            x = g[col]
            rows.append(dict(zip(keys, k), metric=col, label=label, n=int(len(x)),
                             median=x.median(), q25=x.quantile(.25), q75=x.quantile(.75),
                             min=x.min(), max=x.max()))
    return pd.DataFrame(rows)
