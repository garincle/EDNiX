"""
fs_ednix_matched_cortex.py -- FreeSurfer vs EDNiX on MATCHED cortex and MATCHED definitions
(diagnostic behind Supp. Fig. 14).

Supp. Fig. 14 compares EDNiX level-1 Isocortex (atlas-defined, excludes the periallocortex
and allocortex) with the FreeSurfer cortex label, and EDNiX's thickness definition with
FreeSurfer's. Region and definition then differ as well as the tools. Here both tools are
measured on the same region with the same code and the same definitions:

  region   whole cortical mantle of each tool's own surfaces: ?h.cortex.label of its own
           FreeSurfer subject (recon-all: SUBJECTS_DIR/sub-*_ses-*; EDNiX: the FreeSurfer
           subject it builds in sub-*/ses-*/anat/freesurfer/*). If a tool has no
           cortex.label, vertices with ?h.thickness > 0 (FreeSurfer sets the medial wall
           to 0) -- the mask used is written in the csv. No atlas is involved.
  area     pial vertex areas summed over the region, both hemispheres
  thickness_nearest  ?h.thickness (FreeSurfer nearest-point definition, Fischl & Dale 2000)
           averaged over the region, mean of L and R
  thickness_vertex   corresponding-vertex pial-white distance (EDNiX definition), same vertices
  volume   cortical grey-matter volume between the closed surfaces: V(pial) - V(white),
           each from the divergence theorem (sum over faces of v0 . (v1 x v2) / 6), both
           hemispheres. Identical computation for the two tools; includes the thin
           medial-wall slab in both.

Same scan for both tools (latest session of BIDS_Park, the Supp. Fig. 14 pairs).

Outputs  <PAPER_DIR>/Validation/FreeSurfer/matched_cortex/
  matched_cortex.png/.pdf  a-d recon-all vs EDNiX (area, thickness nearest, thickness
                           vertex, volume) with identity line, r, ICC(C,1), ICC(A,1), bias;
                           e thickness vs age (both tools); f Pearson r between all measures
  matched_cortex.csv       one row per subject
  matched_cortex_stats.txt agreement statistics (same functions as Supp. Fig. 14)
  QC_<subject>.png         white (blue) / pial (red) contours of each tool on its own T1
                           (3 coronal + 1 axial slices at the same relative positions) and
                           vertex-thickness histograms: shows which tool misplaces a surface

    python fs_ednix_matched_cortex.py
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import nibabel as nib
from scipy import ndimage, stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(0, MAIN_PATH)
from Exemples.Study_EDNiX.analysis import paper_data as pdat   # noqa: E402
from Exemples.Study_EDNiX.analysis import paper_style as ps    # noqa: E402
from Exemples.Study_EDNiX.analysis.validation import supp_freesurfer_comparison as S  # noqa: E402

opj = os.path.join

BIDS_DIR = '/scratch2/EDNiX/Human/BIDS_Park'
OUT_DIR = opj(pdat.PAPER_DIR, 'Validation', 'FreeSurfer', 'matched_cortex')
WHITE, PIAL = '#1F77B4', '#D62728'
# (column, label, display divisor, unit)
MEASURES = [('area', 'Cortical surface area (pial)', 100.0, 'cm$^2$'),
            ('thickness_nearest', 'Thickness, nearest-point', 1.0, 'mm'),
            ('thickness_vertex', 'Thickness, corresponding vertex', 1.0, 'mm'),
            ('volume', 'Cortical volume V(pial) $-$ V(white)', 1000.0, 'cm$^3$')]


# ----------------------------------------------------------------------------- inputs
def ednix_fs_subject(sub, ses):
    """EDNiX's FreeSurfer subject folder for one session (has surf/lh.white), or None."""
    hits = sorted(glob.glob(opj(BIDS_DIR, f'sub-{sub}', f'ses-{ses}', 'anat', 'freesurfer', '*',
                                'surf', 'lh.white')))
    return os.path.dirname(os.path.dirname(hits[0])) if hits else None


def t1_volume(sdir):
    # orig.mgz first for both tools: the subject's own T1 (EDNiX: the native T1 given to
    # FreeSurfer; its brain.mgz is a normalised helper image built from the warped seg-4FS,
    # useless to judge where the real grey matter is).
    for name in ('orig.mgz', 'orig_raw.mgz', 'T1.mgz', 'norm.mgz', 'brain.mgz', 'brain_raw.mgz'):
        f = opj(sdir, 'mri', name)
        if os.path.isfile(f):
            return f
    return None


def read_hemi(sdir, h):
    """white, pial (vertices, faces), thickness, region index and the mask rule used."""
    white, faces = nib.freesurfer.read_geometry(opj(sdir, 'surf', f'{h}.white'))
    pial, _ = nib.freesurfer.read_geometry(opj(sdir, 'surf', f'{h}.pial'))
    tf = opj(sdir, 'surf', f'{h}.thickness')
    th = nib.freesurfer.read_morph_data(tf) if os.path.isfile(tf) else None
    lab = opj(sdir, 'label', f'{h}.cortex.label')
    if os.path.isfile(lab):
        idx, rule = nib.freesurfer.read_label(lab), 'cortex.label'
    elif th is not None:
        idx, rule = np.where(th > 0)[0], 'thickness>0'
    else:
        idx, rule = np.arange(len(white)), 'all vertices'
    return dict(white=white, pial=pial, faces=faces, th=th, idx=idx, rule=rule)


def closed_volume(verts, faces):
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return abs(float(np.einsum('ij,ij->i', v0, np.cross(v1, v2)).sum()) / 6.0)


def measures(sdir):
    """Matched-region measures of one FreeSurfer subject (both hemispheres)."""
    hemis = {h: read_hemi(sdir, h) for h in ('lh', 'rh')}
    area, near, vert, vol = 0.0, [], [], 0.0
    for d in hemis.values():
        i = d['idx']
        area += S.vertex_areas(d['pial'], d['faces'])[i].sum()
        if d['th'] is not None:
            near.append(d['th'][i].mean())
        vert.append(np.linalg.norm(d['pial'][i] - d['white'][i], axis=1).mean())
        vol += closed_volume(d['pial'], d['faces']) - closed_volume(d['white'], d['faces'])
    return dict(area=area, thickness_nearest=float(np.mean(near)) if near else np.nan,
                thickness_vertex=float(np.mean(vert)), volume=vol,
                mask='/'.join(sorted({d['rule'] for d in hemis.values()}))), hemis


def human_ages():
    try:
        from Exemples.Study_EDNiX.analysis.library.longitudinal import format_human_age
        a = format_human_age([BIDS_DIR])
    except Exception as e:  # noqa: BLE001 -- ages are optional
        print(f'  [age] not available: {e}')
        return pd.DataFrame(columns=['subject', 'session', 'age'])
    a['subject'] = a['subject'].map(S._norm)
    a['session'] = a['session'].map(S._norm)
    return a


# ----------------------------------------------------------------------------- drawing
def plane_segments(verts, faces, axis, value):
    """Segments (n, 2, 3) where the mesh crosses the plane coord[axis] = value."""
    d = verts[:, axis] - value
    df = d[faces]
    cross = (df.min(1) < 0) & (df.max(1) > 0)
    segs = []
    for tri, dd in zip(faces[cross], df[cross]):
        pts = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if (dd[a] < 0) != (dd[b] < 0):
                t = dd[a] / (dd[a] - dd[b])
                pts.append(verts[tri[a]] + t * (verts[tri[b]] - verts[tri[a]]))
        if len(pts) == 2:
            segs.append(pts)
    return np.asarray(segs).reshape(-1, 2, 3)


def slice_image(img, axis, value, ranges, step=0.5):
    """Volume resampled on a plane of tkr RAS (axis fixed at value); returns image and
    the (h_axis, v_axis) extents in mm."""
    inv = np.linalg.inv(img.header.get_vox2ras_tkr())
    data = np.asarray(img.dataobj, dtype=np.float32)
    h_ax, v_ax = [a for a in range(3) if a != axis]
    hs = np.arange(ranges[h_ax][0], ranges[h_ax][1], step)
    vs = np.arange(ranges[v_ax][0], ranges[v_ax][1], step)
    H, V = np.meshgrid(hs, vs)
    ras = np.zeros((4, H.size))
    ras[h_ax], ras[v_ax], ras[axis], ras[3] = H.ravel(), V.ravel(), value, 1
    vox = inv @ ras
    out = ndimage.map_coordinates(data, vox[:3], order=1, mode='constant').reshape(H.shape)
    return out, (hs[0], hs[-1], vs[0], vs[-1]), h_ax, v_ax


def draw_slice(ax, img, hemis, axis, frac):
    allw = np.vstack([d['white'] for d in hemis.values()] + [d['pial'] for d in hemis.values()])
    lo, hi = allw.min(0) - 5, allw.max(0) + 5
    value = lo[axis] + 5 + frac * (hi[axis] - lo[axis] - 10)
    ranges = list(zip(lo, hi))
    im, ext, h_ax, v_ax = slice_image(img, axis, value, ranges)
    vmax = np.percentile(im[im > 0], 99) if (im > 0).any() else 1
    ax.imshow(im, origin='lower', extent=ext, cmap='gray', vmin=0, vmax=vmax, interpolation='nearest')
    for kind, col in (('white', WHITE), ('pial', PIAL)):
        for d in hemis.values():
            s = plane_segments(d[kind], d['faces'], axis, value)
            if len(s):
                ax.add_collection(LineCollection(s[:, :, [h_ax, v_ax]], colors=col, linewidths=0.35))
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'{"xyz"[axis]} = {value:.0f} mm', fontsize=5, pad=1)


def qc_page(fsid, tools, out):
    """tools: [(name, subject dir or None, hemis or None)]."""
    fig, axes = plt.subplots(2, 5, figsize=(ps.WIDTH_2COL, 80 * ps.MM),
                             gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 1.1]))
    hist = {}
    for i, (name, sdir, hemis) in enumerate(tools):
        if sdir is None or hemis is None or t1_volume(sdir) is None:
            for ax in axes[i, :4]:
                ax.axis('off')
            axes[i, 0].text(0, 0.5, f'{name}: no surfaces / T1', fontsize=6)
            continue
        img = nib.load(t1_volume(sdir))
        for j, (axis, frac) in enumerate(((1, 0.3), (1, 0.5), (1, 0.7), (2, 0.55))):
            draw_slice(axes[i, j], img, hemis, axis, frac)
        axes[i, 0].set_ylabel(f'{name}\n({os.path.basename(t1_volume(sdir))})', fontsize=6)
        th = [d['th'][d['idx']] for d in hemis.values() if d['th'] is not None]
        if th:
            hist[name] = np.concatenate(th)
    gs = axes[0, 4].get_gridspec()
    axes[0, 4].remove()
    axes[1, 4].remove()
    ax = fig.add_subplot(gs[:, 4])
    bins = np.linspace(0, 5, 101)
    for (name, v), col in zip(hist.items(), ('0.3', ps.species_color('Human'))):
        ax.hist(v, bins=bins, histtype='step', color=col, lw=0.8, density=True,
                label=f'{name}: median {np.median(v):.2f} mm')
    ax.set_xlabel('Vertex thickness (mm, nearest-point)')
    ax.set_ylabel('Density')
    ax.legend(fontsize=4.5, loc='upper right')
    fig.suptitle(f'{fsid}   white = blue, pial = red', fontsize=7, y=0.995)
    fig.tight_layout(pad=0.4, w_pad=0.3, h_pad=0.6)
    fig.savefig(out, dpi=300)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def scatter(ax, x, y, s, lab, unit, color):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) >= 2:
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        pad = 0.08 * (hi - lo if hi > lo else abs(hi) or 1)
        lim = (lo - pad, hi + pad)
        ax.plot(lim, lim, color='0.55', lw=0.6, ls='--')
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect('equal', adjustable='box')
    ax.scatter(x, y, s=9, color=color, edgecolor='white', linewidth=0.3, zorder=3)
    if s.get('n', 0) >= 3:
        ax.text(0.03, 0.97, f"n = {int(s['n'])}\nr = {s['pearson_r']:.2f} [{s['r_ci_low']:.2f}, {s['r_ci_high']:.2f}]"
                f"\nICC(C,1) = {s['icc_c1']:.2f}\nICC(A,1) = {s['icc_a1']:.2f}\nbias = {s['bias_pct']:+.1f}%",
                transform=ax.transAxes, va='top', ha='left', fontsize=4.8)
    ax.set_xlabel(f'recon-all ({unit})')
    ax.set_ylabel(f'EDNiX ({unit})')
    ax.set_title(lab, fontsize=6)


def main():
    plt.rcParams.update(ps.PAPER_RC)
    os.makedirs(OUT_DIR, exist_ok=True)
    rows = []
    for fs_dir in sorted(glob.glob(opj(S.FS_DIR, 'sub-*_ses-*'))):
        fsid = os.path.basename(fs_dir)
        if not os.path.isfile(opj(fs_dir, 'scripts', 'recon-all.done')):
            continue
        m = re.match(r'sub-(.+)_ses-(.+)$', fsid)
        sub, ses = m.group(1), m.group(2)
        ed_dir = ednix_fs_subject(sub, ses)
        row = dict(subject=S._norm(sub), session=S._norm(ses), ednix_fs_dir=ed_dir)
        fs_m, fs_h = measures(fs_dir)
        row.update({f'{k}_fs': v for k, v in fs_m.items()})
        ed_h = None
        if ed_dir is not None:
            ed_m, ed_h = measures(ed_dir)
            row.update({f'{k}_ednix': v for k, v in ed_m.items()})
        print(f"  {fsid}: EDNiX = {ed_dir}; masks recon-all {fs_m['mask']}, "
              f"EDNiX {row.get('mask_ednix', '-')}")
        rows.append(row)
        try:
            qc_page(fsid, [('recon-all', fs_dir, fs_h), ('EDNiX', ed_dir, ed_h)],
                    opj(OUT_DIR, f'QC_{fsid}.png'))
        except Exception as e:  # noqa: BLE001 -- one bad subject must not stop the others
            print(f'  [QC] {fsid}: {e}')
    if not rows:
        sys.exit(f'No finished recon-all subject in {S.FS_DIR}')
    df = pd.DataFrame(rows).merge(human_ages(), on=['subject', 'session'], how='left')
    df.to_csv(opj(OUT_DIR, 'matched_cortex.csv'), index=False)

    color = ps.species_color('Human')
    fig, axes = plt.subplots(2, 3, figsize=(ps.WIDTH_2COL, 125 * ps.MM))
    axes = axes.ravel()
    lines = ['recon-all vs EDNiX, matched cortex (each tool: whole cortex label of its own '
             'surfaces), identical definitions; bias = (EDNiX - recon-all) / mean, %.', '']
    for k, (col, lab, div, unit) in enumerate(MEASURES):
        x = df.get(f'{col}_fs', pd.Series(np.nan, index=df.index)).to_numpy(float)
        y = df.get(f'{col}_ednix', pd.Series(np.nan, index=df.index)).to_numpy(float)
        s = S.agreement(x, y)
        if s.get('n', 0) >= 3:
            lines.append(f"{lab} (n = {s['n']}): recon-all median {s['fs_median']:.4g}, EDNiX median "
                         f"{s['ednix_median']:.4g}; r = {s['pearson_r']:.2f} [{s['r_ci_low']:.2f}, "
                         f"{s['r_ci_high']:.2f}]; ICC(C,1) = {s['icc_c1']:.2f}; ICC(A,1) = {s['icc_a1']:.2f}; "
                         f"bias = {s['bias_pct']:+.1f}% (LoA {s['loa_low']:+.1f} to {s['loa_high']:+.1f}%)")
        scatter(axes[k], x / div, y / div, s, lab, unit, color)
        ps.panel_label(axes[k], 'abcd'[k], x=-0.25, y=1.02)

    ax = axes[4]
    if df['age'].notna().sum() >= 3:
        for tool, c in (('fs', '0.3'), ('ednix', color)):
            col = f'thickness_nearest_{tool}'
            if col not in df:
                continue
            ok = df['age'].notna() & df[col].notna()
            a, v = df.loc[ok, 'age'].to_numpy(float), df.loc[ok, col].to_numpy(float)
            if len(a) < 3:
                continue
            r, p = stats.pearsonr(a, v)
            ax.scatter(a, v, s=9, color=c, edgecolor='white', lw=0.3,
                       label=f"{'recon-all' if tool == 'fs' else 'EDNiX'}: r = {r:.2f}, P = {p:.2g}")
            sl, ic = np.polyfit(a, v, 1)
            xs = np.array([a.min(), a.max()])
            ax.plot(xs, sl * xs + ic, color=c, lw=0.7)
            lines.append(f"Thickness (nearest) vs age, {tool}: r = {r:.2f}, P = {p:.2g}, "
                         f"slope {sl:+.4f} mm/year")
        ax.legend(fontsize=4.5, loc='best')
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Thickness, nearest-point (mm)')
    ax.set_title('Plausibility: thinning with age', fontsize=6)
    ps.panel_label(ax, 'e', x=-0.25, y=1.02)

    ax = axes[5]
    cols = [c for c in df.columns if re.match(r'(area|thickness_nearest|thickness_vertex|volume)_(fs|ednix)$', c)]
    R = df[cols].astype(float).corr()
    im = ax.imshow(R.values, vmin=-1, vmax=1, cmap='RdBu_r')
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=60, ha='right', fontsize=4.5)
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(cols, fontsize=4.5)
    for i in range(len(cols)):
        for j in range(len(cols)):
            ax.text(j, i, f'{R.values[i, j]:.2f}', ha='center', va='center', fontsize=3.8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label('Pearson r', fontsize=5)
    ps.panel_label(ax, 'f', x=-0.45, y=1.02)
    fig.tight_layout(w_pad=1.5, h_pad=1.5)
    ps.save(fig, opj(OUT_DIR, 'matched_cortex'))
    plt.close(fig)
    with open(opj(OUT_DIR, 'matched_cortex_stats.txt'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print(df[['subject', 'session', 'age'] + cols].round(3).to_string(index=False))
    print(f'-> {OUT_DIR}')


if __name__ == '__main__':
    main()
