"""
supp_freesurfer_comparison.py -- Supp. Fig. 14: subject-by-subject agreement between EDNiX and
standard FreeSurfer recon-all (Fischl 2012) on the same human T1w images.

Question: where a reference tool exists (adult human brain), does EDNiX's cross-species
design cost accuracy? Agreement, not superiority.

Data
  FreeSurfer  /scratch2/EDNiX/FreeSurfer_benchmark/BIDS_Park/sub-*_ses-*  (freesurfer_vs_ednix.py),
              subjects with scripts/recon-all.done only. Bilateral:
                surface area  sum of pial vertex areas over ?h.cortex.label
                thickness     mean of ?h.thickness over the same vertices (Fischl & Dale 2000:
                              symmetric nearest-point distance), mean of L and R
                volume        cortical ribbon voxels (ribbon.mgz labels 3 and 42) x voxel volume
  EDNiX       level-1 Isocortex, bilateral, same subject and session (paper_data.collect cache).
              EDNiX area = pial vertex areas, thickness = pial-to-white distance between
              corresponding vertices (wb_command -surface-to-surface-3d-distance), volume =
              atlas voxels.
  The thickness definitions differ, so the FreeSurfer surfaces are also measured with the
  EDNiX definition (thickness_vertex_mm, stats only): it separates a definition offset from
  a surface-placement offset. The FreeSurfer cortex label includes the medial temporal
  cortex that the EDNiX Isocortex excludes: a positive FreeSurfer offset on area and volume
  is expected; between-subject agreement (r, ICC(C,1)) is the main readout.

Statistics (per measure, n = subjects with both values)
  Pearson r [95% CI, Fisher z], Spearman rho, ICC(C,1) and ICC(A,1) [95% CI] (two-way,
  single measure; McGraw & Wong 1996), Bland-Altman on (EDNiX - FreeSurfer) / mean x 100:
  bias [95% CI, t] and 95% limits of agreement (Bland & Altman 1999), paired t-test of
  the raw difference.

Figure (183 mm wide): a-c EDNiX vs FreeSurfer per subject (dashed: identity; solid: OLS),
d-f Bland-Altman (solid: bias, shaded: its 95% CI, dashed: limits of agreement).

Outputs  <PAPER_DIR>/Supp_Fig_14/
  Supp_Fig_14.pdf / .png
  Supp_Fig_14_stats.xlsx   sheets paired (one row per subject), stats (one row per measure),
                           fs_status (recon-all outcome)
  Supp_Fig_14_stats.txt    the same, as sentences for Results/legend

    python supp_freesurfer_comparison.py
    python supp_freesurfer_comparison.py --refresh   # rebuild the EDNiX morphometry cache
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import nibabel as nib
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(0, MAIN_PATH)
from Exemples.Study_EDNiX.analysis import paper_data as pdat   # noqa: E402
from Exemples.Study_EDNiX.analysis import paper_style as ps    # noqa: E402

opj = os.path.join

SPECIES = 'Human'
BIDS = 'BIDS_Park'
FS_DIR = opj('/scratch2/EDNiX/FreeSurfer_benchmark', BIDS)
OUT_BASE = opj(pdat.PAPER_DIR, 'Supp_Fig_14', 'Supp_Fig_14')
STATUS_CSV = opj(pdat.PAPER_DIR, 'Validation', 'FreeSurfer', 'fs_status.csv')
REGION = 'Isocortex'
# (column, axis label, display divisor, display unit)
MEASURES = [('surface_area_mm2', 'Surface area', 100.0, 'cm$^2$'),
            ('thickness_mm', 'Thickness', 1.0, 'mm'),
            ('volume_mm3', 'Volume', 1000.0, 'cm$^3$')]
MODALITY = {'surface_area_mm2': 'surface', 'thickness_mm': 'thickness', 'volume_mm3': 'volume'}


# ----------------------------------------------------------------------------- FreeSurfer
def vertex_areas(coords, faces):
    a, b, c = coords[faces[:, 0]], coords[faces[:, 1]], coords[faces[:, 2]]
    fa = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    va = np.zeros(len(coords))
    for k in range(3):
        np.add.at(va, faces[:, k], fa / 3.0)
    return va


def fs_measures(sdir):
    """Bilateral area / thickness (FreeSurfer and EDNiX definitions) / volume of one subject."""
    area, thick, thick_v = 0.0, [], []
    for h in ('lh', 'rh'):
        pial, faces = nib.freesurfer.read_geometry(opj(sdir, 'surf', f'{h}.pial'))
        white, _ = nib.freesurfer.read_geometry(opj(sdir, 'surf', f'{h}.white'))
        th = nib.freesurfer.read_morph_data(opj(sdir, 'surf', f'{h}.thickness'))
        lab = opj(sdir, 'label', f'{h}.cortex.label')
        idx = nib.freesurfer.read_label(lab) if os.path.isfile(lab) else np.where(th > 0)[0]
        area += vertex_areas(pial, faces)[idx].sum()
        thick.append(th[idx].mean())
        thick_v.append(np.linalg.norm(pial[idx] - white[idx], axis=1).mean())
    rib = nib.load(opj(sdir, 'mri', 'ribbon.mgz'))
    r = np.asarray(rib.dataobj)
    vol = float(np.isin(r, (3, 42)).sum()) * float(np.prod(rib.header.get_zooms()[:3]))
    return dict(surface_area_mm2=area, thickness_mm=float(np.mean(thick)),
                thickness_vertex_mm=float(np.mean(thick_v)), volume_mm3=vol)


def collect_fs(fs_dir=FS_DIR):
    """One row per FreeSurfer subject with recon-all.done; unfinished ones are listed."""
    rows, skipped = [], []
    for sdir in sorted(glob.glob(opj(fs_dir, 'sub-*_ses-*'))):
        name = os.path.basename(sdir)
        m = re.match(r'sub-(.+)_ses-(.+)$', name)
        if not os.path.isfile(opj(sdir, 'scripts', 'recon-all.done')):
            skipped.append(name)
            continue
        try:
            rows.append(dict(subject=m.group(1), session=m.group(2), **fs_measures(sdir)))
        except (FileNotFoundError, OSError, ValueError) as e:
            skipped.append(f'{name} ({e})')
    if skipped:
        print(f'  [FreeSurfer] not used (unfinished or unreadable): {", ".join(skipped)}')
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- EDNiX
def collect_ednix(refresh=False):
    """Level-1 Isocortex, bilateral, one row per (subject, session) of BIDS."""
    raw = pdat.collect(pdat.ANAT_DATASETS, ('surface', 'volume', 'thickness'), 'morphometry', refresh)
    out = None
    for col, _, _, _ in MEASURES:
        df = raw.get(MODALITY[col]).copy()
        df['bids_dir'] = df['bids_dir'].map(pdat.bids_name)
        df = df[(df['bids_dir'] == BIDS) & (df['atlas_level'] == 1)]
        df = pdat.select_regions(pdat.bilateral(df, col), [REGION])
        df = df.groupby(['subject', 'session'], as_index=False)[col].mean()
        out = df if out is None else out.merge(df, on=['subject', 'session'], how='outer')
    return out


def _norm(s):
    s = str(s).replace('sub-', '').replace('ses-', '')
    return str(int(s)) if s.isdigit() else s


def paired_table(fs, ed):
    for d in (fs, ed):
        d['subject'] = d['subject'].map(_norm)
        d['session'] = d['session'].map(_norm)
    fs = fs.rename(columns={c: f'{c}_fs' for c in fs.columns if c not in ('subject', 'session')})
    ed = ed.rename(columns={c: f'{c}_ednix' for c in ed.columns if c not in ('subject', 'session')})
    m = fs.merge(ed, on=['subject', 'session'], how='left')
    miss = m.loc[m[[f'{c}_ednix' for c, *_ in MEASURES]].isna().all(axis=1), 'subject'].tolist()
    if miss:
        print(f'  [EDNiX] no level-1 {REGION} value for FreeSurfer subjects: {miss}')
    return m


# ----------------------------------------------------------------------------- statistics
def icc_two_way(x, y, alpha=0.05):
    """ICC(C,1) and ICC(A,1) with (1 - alpha) CIs, two raters (McGraw & Wong 1996,
    cases 3 and A; CI of ICC(A,1) with the Satterthwaite df)."""
    Y = np.column_stack([x, y])
    n, k = Y.shape
    gm = Y.mean()
    msr = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1)
    msc = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    sse = ((Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm) ** 2).sum()
    mse = sse / ((n - 1) * (k - 1))
    icc_c = (msr - mse) / (msr + (k - 1) * mse)
    icc_a = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
    q = 1 - alpha / 2
    f_obs = msr / mse
    df1, df2 = n - 1, (n - 1) * (k - 1)
    fl, fu = f_obs / stats.f.ppf(q, df1, df2), f_obs * stats.f.ppf(q, df2, df1)
    ci_c = ((fl - 1) / (fl + k - 1), (fu - 1) / (fu + k - 1))
    a = k * icc_a / (n * (1 - icc_a))
    b = 1 + k * icc_a * (n - 1) / (n * (1 - icc_a))
    v = (a * msc + b * mse) ** 2 / ((a * msc) ** 2 / (k - 1) + (b * mse) ** 2 / df2)
    fl_a, fu_a = stats.f.ppf(q, n - 1, v), stats.f.ppf(q, v, n - 1)
    lo_a = n * (msr - fl_a * mse) / (fl_a * (k * msc + (k * n - k - n) * mse) + n * msr)
    hi_a = n * (fu_a * msr - mse) / (k * msc + (k * n - k - n) * mse + n * fu_a * msr)
    return icc_c, ci_c, icc_a, (lo_a, hi_a)


def agreement(x, y):
    """Agreement statistics of EDNiX (y) vs FreeSurfer (x)."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x, float)[ok], np.asarray(y, float)[ok]
    n = len(x)
    if n < 3:
        return dict(n=n)
    r, p_r = stats.pearsonr(x, y)
    z, se = np.arctanh(np.clip(r, -0.999999, 0.999999)), 1 / np.sqrt(n - 3) if n > 3 else np.nan
    rho, p_rho = stats.spearmanr(x, y)
    ic, ci_c, ia, ci_a = icc_two_way(x, y)
    d = 100 * (y - x) / ((x + y) / 2)
    bias, sd = d.mean(), d.std(ddof=1)
    t = stats.t.ppf(0.975, n - 1)
    slope, intercept = np.polyfit(x, y, 1)
    t_pair = stats.ttest_rel(y, x)
    return dict(n=n, fs_median=np.median(x), ednix_median=np.median(y),
                pearson_r=r, r_ci_low=np.tanh(z - 1.96 * se), r_ci_high=np.tanh(z + 1.96 * se),
                pearson_p=p_r, spearman_rho=rho, spearman_p=p_rho,
                icc_c1=ic, icc_c1_ci_low=ci_c[0], icc_c1_ci_high=ci_c[1],
                icc_a1=ia, icc_a1_ci_low=ci_a[0], icc_a1_ci_high=ci_a[1],
                bias_pct=bias, bias_ci_low=bias - t * sd / np.sqrt(n),
                bias_ci_high=bias + t * sd / np.sqrt(n),
                loa_low=bias - 1.96 * sd, loa_high=bias + 1.96 * sd,
                paired_t=t_pair.statistic, paired_p=t_pair.pvalue,
                ols_slope=slope, ols_intercept=intercept)


def stats_table(m):
    rows = []
    for col, lab, _, _ in MEASURES:
        rows.append(dict(measure=lab, fs_column=f'{col}_fs', ednix_column=f'{col}_ednix',
                         **agreement(m[f'{col}_fs'].to_numpy(float), m[f'{col}_ednix'].to_numpy(float))))
    if 'thickness_vertex_mm_fs' in m:
        rows.append(dict(measure='Thickness (FreeSurfer surfaces, EDNiX definition)',
                         fs_column='thickness_vertex_mm_fs', ednix_column='thickness_mm_ednix',
                         **agreement(m['thickness_vertex_mm_fs'].to_numpy(float),
                                     m['thickness_mm_ednix'].to_numpy(float))))
    return pd.DataFrame(rows)


def p_text(p):
    return 'P < 0.001' if p < 0.001 else f'P = {p:.2g}'


def sentences(st, status):
    lines = [f'Supp. Fig. 14 -- EDNiX vs standard FreeSurfer recon-all, {SPECIES} ({BIDS}), '
             f'{REGION} (EDNiX level 1) vs FreeSurfer cortex label, bilateral, one session per '
             'subject (latest). Bias and limits of agreement (LoA) on (EDNiX - FreeSurfer) / mean, %.', '']
    if status is not None:
        ok = status['status'].astype(str).str.startswith(('ok', 'done')).sum()
        lines.append(f'recon-all completed without error for {ok}/{len(status)} subjects.')
        lines.append('')
    for _, s in st.iterrows():
        if s['n'] < 3:
            lines.append(f'{s["measure"]}: n = {int(s["n"])}, not enough pairs.')
            continue
        lines.append(
            f'{s["measure"]} (n = {int(s["n"])}): FreeSurfer median {s["fs_median"]:.4g}, EDNiX median '
            f'{s["ednix_median"]:.4g}; r = {s["pearson_r"]:.2f} [{s["r_ci_low"]:.2f}, {s["r_ci_high"]:.2f}], '
            f'{p_text(s["pearson_p"])}; Spearman rho = {s["spearman_rho"]:.2f}; '
            f'ICC(C,1) = {s["icc_c1"]:.2f} [{s["icc_c1_ci_low"]:.2f}, {s["icc_c1_ci_high"]:.2f}]; '
            f'ICC(A,1) = {s["icc_a1"]:.2f} [{s["icc_a1_ci_low"]:.2f}, {s["icc_a1_ci_high"]:.2f}]; '
            f'bias = {s["bias_pct"]:+.1f}% [{s["bias_ci_low"]:+.1f}, {s["bias_ci_high"]:+.1f}] '
            f'(paired t = {s["paired_t"]:.2f}, {p_text(s["paired_p"])}), '
            f'LoA {s["loa_low"]:+.1f} to {s["loa_high"]:+.1f}%; OLS slope {s["ols_slope"]:.2f}.')
    lines += ['', 'Statistics: Pearson r with Fisher-z 95% CI; ICC two-way, single measure, '
              'consistency (C,1) and absolute agreement (A,1), 95% CI (McGraw & Wong 1996); '
              'Bland-Altman bias with t-based 95% CI and 95% LoA (Bland & Altman 1999).']
    return lines


# ----------------------------------------------------------------------------- figure
def draw(m, st, out_base=OUT_BASE):
    plt.rcParams.update(ps.PAPER_RC)
    color = ps.species_color(SPECIES)
    marker = ps.dataset_marker(BIDS, SPECIES)
    fig, axes = plt.subplots(2, 3, figsize=(ps.WIDTH_2COL, 105 * ps.MM))
    for j, (col, lab, div, unit) in enumerate(MEASURES):
        x = m[f'{col}_fs'].to_numpy(float) / div
        y = m[f'{col}_ednix'].to_numpy(float) / div
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        s = st.iloc[j]
        ax = axes[0, j]
        if len(x) >= 3:
            lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
            pad = 0.08 * (hi - lo if hi > lo else abs(hi) or 1)
            lim = (lo - pad, hi + pad)
            ax.plot(lim, lim, color='0.55', lw=0.6, ls='--', zorder=1)
            xs = np.array(lim)
            ax.plot(xs, s['ols_slope'] * xs + s['ols_intercept'] / div, color=color, lw=0.8, zorder=2)
            ax.set_xlim(lim)
            ax.set_ylim(lim)
            ax.set_aspect('equal', adjustable='box')
            ax.text(0.03, 0.97,
                    f'n = {int(s["n"])}\nr = {s["pearson_r"]:.2f} [{s["r_ci_low"]:.2f}, {s["r_ci_high"]:.2f}]\n'
                    f'ICC(C,1) = {s["icc_c1"]:.2f}\nICC(A,1) = {s["icc_a1"]:.2f}',
                    transform=ax.transAxes, va='top', ha='left', fontsize=5)
        ax.scatter(x, y, s=9, marker=marker, facecolor=color, edgecolor='white', linewidth=0.3, zorder=3)
        ax.set_xlabel(f'FreeSurfer {lab.lower()} ({unit})')
        ax.set_ylabel(f'EDNiX {lab.lower()} ({unit})')
        ps.panel_label(ax, 'abc'[j], x=-0.28, y=1.02)

        ax = axes[1, j]
        mean = (x + y) / 2
        d = 100 * (y - x) / mean
        if len(x) >= 3:
            ax.axhspan(s['bias_ci_low'], s['bias_ci_high'], color=color, alpha=0.15, lw=0)
            ax.axhline(s['bias_pct'], color=color, lw=0.8)
            for v in (s['loa_low'], s['loa_high']):
                ax.axhline(v, color='0.4', lw=0.6, ls='--')
            ax.axhline(0, color='0.75', lw=0.5, ls=':')
            ax.text(1.0, s['bias_pct'], f' {s["bias_pct"]:+.1f}%', transform=ax.get_yaxis_transform(),
                    va='center', ha='left', fontsize=5, color=ps.text_color(SPECIES), clip_on=False)
            for v in (s['loa_low'], s['loa_high']):
                ax.text(1.0, v, f' {v:+.1f}%', transform=ax.get_yaxis_transform(), va='center',
                        ha='left', fontsize=5, color='0.35', clip_on=False)
            span = max(abs(s['loa_low']), abs(s['loa_high']), np.abs(d).max(), 1.0) * 1.2
            ax.set_ylim(min(-0.25 * span, s['loa_low'] - 0.15 * span), max(0.25 * span, s['loa_high'] + 0.15 * span))
        ax.scatter(mean, d, s=9, marker=marker, facecolor=color, edgecolor='white', linewidth=0.3, zorder=3)
        ax.set_xlabel(f'Mean {lab.lower()} ({unit})')
        ax.set_ylabel('EDNiX $-$ FreeSurfer (% of mean)')
        ps.panel_label(ax, 'def'[j], x=-0.28, y=1.02)
    fig.tight_layout(w_pad=2.5, h_pad=1.5)
    ps.save(fig, out_base)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--refresh', action='store_true', help='rebuild the EDNiX morphometry cache')
    a = ap.parse_args()
    fs = collect_fs()
    if fs.empty:
        sys.exit(f'No finished FreeSurfer subject in {FS_DIR} (scripts/recon-all.done).')
    m = paired_table(fs, collect_ednix(a.refresh))
    st = stats_table(m)
    status = pd.read_csv(STATUS_CSV) if os.path.isfile(STATUS_CSV) else None
    os.makedirs(os.path.dirname(OUT_BASE), exist_ok=True)
    with pd.ExcelWriter(f'{OUT_BASE}_stats.xlsx') as xw:
        m.to_excel(xw, sheet_name='paired', index=False)
        st.to_excel(xw, sheet_name='stats', index=False)
        if status is not None:
            status.to_excel(xw, sheet_name='fs_status', index=False)
    lines = sentences(st, status)
    with open(f'{OUT_BASE}_stats.txt', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    draw(m, st)
    print(f'-> {os.path.dirname(OUT_BASE)}')


if __name__ == '__main__':
    main()
