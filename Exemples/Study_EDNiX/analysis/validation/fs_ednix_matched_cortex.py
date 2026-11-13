"""
fs_ednix_matched_cortex.py -- FreeSurfer vs EDNiX on MATCHED cortex and MATCHED definitions
(diagnostic behind Supp. Fig. 7).

Supp. Fig. 7 compares EDNiX level-1 Isocortex (atlas-defined, excludes the periallocortex
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

Same scan for both tools (latest session of BIDS_Park, the Supp. Fig. 7 pairs).

Outputs  <PAPER_DIR>/Validation/FreeSurfer/matched_cortex/
  matched_cortex.png/.pdf  a-d recon-all vs EDNiX (area, thickness nearest, thickness
                           vertex, volume) with identity line, r, ICC(C,1), ICC(A,1), bias;
                           e thickness vs age (both tools); f Pearson r between all measures
  matched_cortex.csv       one row per subject
  matched_cortex_stats.txt agreement statistics (same functions as Supp. Fig. 7)
  QC_<subject>.png         both tools on the SAME image and the SAME slices: the recon-all
                           orig.mgz (the subject's own T1, not skull-stripped). EDNiX surfaces
                           are brought into it by a rigid mri_robust_register of the EDNiX T1
                           onto the recon-all brainmask (same scan), LTA cached in
                           registration/. Rows: recon-all white/pial, EDNiX white/pial, both
                           pial surfaces together; right: vertex-thickness histograms
  matched_cortex_stats.txt also gives thickness vs age partialling out surface area (EDNiX
                           thickness tracks brain size; an age effect could go through size)

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


# ----------------------------------------------------------------------------- common space
def lta_ras2ras(lta, src_img, dst_img):
    """4x4 scanner-RAS -> scanner-RAS matrix from a FreeSurfer LTA (vox2vox or ras2ras)."""
    lines = open(lta).read().splitlines()
    typ = next(int(l.split('=')[1].split()[0]) for l in lines if l.strip().startswith('type'))
    k = next(i for i, l in enumerate(lines) if l.split() == ['1', '4', '4'])
    M = np.array([[float(v) for v in lines[k + 1 + r].split()] for r in range(4)])
    if typ == 0:   # LINEAR_VOX_TO_VOX
        M = dst_img.affine @ M @ np.linalg.inv(src_img.affine)
    return M


def ednix_to_fs(fs_dir, ed_dir, fsid):
    """tkr-RAS(EDNiX orig.mgz) -> tkr-RAS(recon-all orig.mgz), by mri_robust_register of the
    EDNiX T1 onto the recon-all brainmask (same scan: rigid, 6 dof). Returns (4x4, note)."""
    from Exemples.Study_EDNiX.analysis.validation import freesurfer_vs_ednix as R
    mov, dst = opj(ed_dir, 'mri', 'orig.mgz'), opj(fs_dir, 'mri', 'brainmask.mgz')
    fs_orig = opj(fs_dir, 'mri', 'orig.mgz')
    if not all(os.path.isfile(f) for f in (mov, dst, fs_orig)):
        return None, 'no common space (orig/brainmask missing)'
    reg_dir = opj(OUT_DIR, 'registration')
    os.makedirs(reg_dir, exist_ok=True)
    lta = opj(reg_dir, f'{fsid}_ednix_to_reconall.lta')
    if not os.path.isfile(lta):
        lic_dirs, lic_env = R.license_bind()
        bind = ','.join(sorted({ed_dir, fs_dir, reg_dir, *lic_dirs}))
        cmd = (f'{lic_env}singularity exec --bind {bind} {R.FS_SIF} mri_robust_register '
               f'--mov {mov} --dst {dst} --lta {lta} --satit --iscale --maxit 20')
        rc = R.sh(cmd, opj(reg_dir, f'{fsid}.log'))
        if rc != 0 or not os.path.isfile(lta):
            return None, f'mri_robust_register failed (see {reg_dir}/{fsid}.log)'
    mov_img, dst_img, fs_img = nib.load(mov), nib.load(dst), nib.load(fs_orig)
    ras2ras = lta_ras2ras(lta, mov_img, dst_img)
    T = (np.linalg.inv(fs_img.affine @ np.linalg.inv(fs_img.header.get_vox2ras_tkr()))
         @ ras2ras @ mov_img.affine @ np.linalg.inv(mov_img.header.get_vox2ras_tkr()))
    return T, 'EDNiX registered to recon-all (mri_robust_register, rigid)'


def transform_hemis(hemis, T):
    out = {}
    for h, d in hemis.items():
        e = dict(d)
        for k in ('white', 'pial'):
            v = np.c_[d[k], np.ones(len(d[k]))] @ T.T
            e[k] = v[:, :3]
        out[h] = e
    return out


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
    """Volume resampled on a plane of tkr RAS (axis fixed at value)."""
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


def slice_planes(hemis):
    """The 4 planes (axis, value) and the display box, from one set of surfaces."""
    allv = np.vstack([d[k] for d in hemis.values() for k in ('white', 'pial')])
    lo, hi = allv.min(0) - 8, allv.max(0) + 8
    planes = [(axis, lo[axis] + 8 + f * (hi[axis] - lo[axis] - 16))
              for axis, f in ((1, 0.3), (1, 0.5), (1, 0.7), (2, 0.55))]
    return planes, list(zip(lo, hi))


def draw_slice(ax, img, layers, axis, value, ranges):
    """layers: [(hemis, white colour, pial colour, linewidth)]."""
    im, ext, h_ax, v_ax = slice_image(img, axis, value, ranges)
    vmax = np.percentile(im[im > 0], 99) if (im > 0).any() else 1
    ax.imshow(im, origin='lower', extent=ext, cmap='gray', vmin=0, vmax=vmax, interpolation='nearest')
    for hemis, cw, cp, lw in layers:
        for kind, col in (('white', cw), ('pial', cp)):
            for d in hemis.values():
                s = plane_segments(d[kind], d['faces'], axis, value)
                if len(s):
                    ax.add_collection(LineCollection(s[:, :, [h_ax, v_ax]], colors=col, linewidths=lw))
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'{"xyz"[axis]} = {value:.0f} mm', fontsize=5, pad=1)


def qc_page(fsid, fs_dir, fs_h, ed_dir, ed_h, out):
    """Same image (recon-all orig.mgz: the subject's own T1, not skull-stripped) and the
    same slices for both tools. Rows: recon-all surfaces, EDNiX surfaces (registered),
    both pial surfaces together. Falls back to EDNiX in its own space if no registration."""
    img = nib.load(opj(fs_dir, 'mri', 'orig.mgz'))
    planes, ranges = slice_planes(fs_h)
    T, note = (ednix_to_fs(fs_dir, ed_dir, fsid) if ed_h is not None else (None, 'no EDNiX surfaces'))
    ed_fs = transform_hemis(ed_h, T) if T is not None else None
    fig, axes = plt.subplots(3, 5, figsize=(ps.WIDTH_2COL, 118 * ps.MM),
                             gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 1.1]))
    rows = [('recon-all', [(fs_h, WHITE, PIAL, 0.35)]),
            ('EDNiX', [(ed_fs, WHITE, PIAL, 0.35)] if ed_fs else []),
            ('pial: recon-all red\nEDNiX yellow', [(fs_h, 'none', PIAL, 0.35)]
             + ([(ed_fs, 'none', '#FFD400', 0.35)] if ed_fs else []))]
    for i, (name, layers) in enumerate(rows):
        for j, (axis, value) in enumerate(planes):
            draw_slice(axes[i, j], img, layers, axis, value, ranges)
        axes[i, 0].set_ylabel(name, fontsize=5.5)
    if ed_fs is None:
        axes[1, 0].text(0.02, 0.5, note, transform=axes[1, 0].transAxes, fontsize=5, color='yellow')
    gs = axes[0, 4].get_gridspec()
    for r in range(3):
        axes[r, 4].remove()
    ax = fig.add_subplot(gs[:, 4])
    bins = np.linspace(0, 5, 101)
    for (name, hemis), col in zip((('recon-all', fs_h), ('EDNiX', ed_h)), ('0.3', ps.species_color('Human'))):
        if hemis is None:
            continue
        th = [d['th'][d['idx']] for d in hemis.values() if d['th'] is not None]
        if th:
            v = np.concatenate(th)
            ax.hist(v, bins=bins, histtype='step', color=col, lw=0.8, density=True,
                    label=f'{name}: median {np.median(v):.2f} mm')
    ax.set_xlabel('Vertex thickness (mm, nearest-point)')
    ax.set_ylabel('Density')
    ax.legend(fontsize=4.5, loc='upper right')
    fig.suptitle(f'{fsid}   recon-all orig.mgz, white = blue, pial = red   [{note}]', fontsize=6, y=0.997)
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
            qc_page(fsid, fs_dir, fs_h, ed_dir, ed_h, opj(OUT_DIR, f'QC_{fsid}.png'))
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
            acol = f'area_{tool}'
            ok2 = ok & df[acol].notna()
            if ok2.sum() >= 4:
                A = np.c_[np.ones(ok2.sum()), df.loc[ok2, acol].to_numpy(float)]
                res = [y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
                       for y in (df.loc[ok2, 'age'].to_numpy(float), df.loc[ok2, col].to_numpy(float))]
                rp = float(np.corrcoef(*res)[0, 1])
                n2 = int(ok2.sum())
                tp = rp * np.sqrt((n2 - 3) / max(1e-12, 1 - rp ** 2))
                pp = 2 * stats.t.sf(abs(tp), n2 - 3)
                lines.append(f"  partial r(thickness, age | area), {tool}: {rp:.2f}, P = {pp:.2g} (n = {n2})")
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
