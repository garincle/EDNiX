"""
Supplementary Figure 2 -- spatial registration in the fMRI module of EDNiX (modalities/fMRI).

  a  one run through its coordinate spaces: raw EPI -> motion correction -> distortion correction
     -> registration to the reference run -> functional AC-PC space (rigid + SyN restricted to the
     phase-encoding axis, to the anatomy) -> anatomy and atlases brought into the EPI ->
     denoised series -> template space
  b  the steps of the module (numbers = Skip_step of the functional launcher)
  c  how the transforms are composed: two interpolations of the time series in total

Images: one processed session of FUNC_BIDS (first subject / session with every file), else a
synthetic brain (MNI ICBM152, frames tagged 'synthetic').
Output: <PAPER_DIR>/Supp_Fig_02/Supp_Fig_02.pdf|png
Former Supplementary/Supp_Fig_02_*.py.
Run: python Exemples/Study_EDNiX/analysis/run_paper.py supp2, or python Exemples/Study_EDNiX/analysis/library/schematic_func.py [--synthetic] [--out DIR]
"""
import os
import sys
import argparse

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), *[os.pardir] * 4))
for _path in (_REPO, '/home/cgarin/PycharmProjects/EDNiX/'):
    if os.path.isdir(_path) and _path not in sys.path:
        sys.path.insert(0, _path)

from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis.library import pipeline_schematic as sch
from Exemples.Study_EDNiX.analysis.library.pipeline_schematic import Thumb, OPS, SPACES

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_02', 'Supp_Fig_02')

FUNC_BIDS = '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin'
REFERENCE = 'EDNiX'
ATLAS = 'EDNIxCSC'
ATLAS_LEVEL_INDEX = 1
T = 'T1w'
PATTERNS = {
    'bold':      'func/*_bold.nii.gz',
    'motion':    'func/raw/matrices/*_space-func_desc-motion_correction.1D',
    'n4':        'func/raw/preprocessing/*_space-func_desc-runMean_n4Bias.nii.gz',
    'unwarped':  'func/raw/preprocessing/*_space-func_desc-runMean_unwarped_reoriented.nii.gz',
    'inref':     'func/raw/preprocessing/*_space-func_desc-fMRI_runMean_inRef.nii.gz',
    'mean_acpc': 'func/acpc-func/preprocessing/all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS.nii.gz',
    'anat_func': f'func/acpc-func/preprocessing/anat_space-acpc-func_{T}.nii.gz',
    'atlas':     f'func/acpc-func/labels/*_seg-{ATLAS}_dseg.nii.gz',
    'residual':  'func/acpc-func/postprocessed_rs/*_space-acpc-func_desc-fMRI_residual.nii.gz',
    'mean_tmpl': f'func/templates/{REFERENCE}/postprocessed_rs/all_runs_space-template-func_desc-fMRI_Mean_Image_SS.nii.gz',
    'tmpl':      f'func/templates/{REFERENCE}/preprocessing/BASE_SS_fMRI.nii.gz',
}

FW, FH = 17.6, 15.8
RED, YEL, GRN = '#FF5A36', '#FFD400', '#3DDC97'


# ── plotted frames ───────────────────────────────────────────────────────────
def motion_plot(params, fd, thr):
    """6 rigid parameters (translations mm, rotations deg) + framewise displacement with threshold."""
    def draw(ax):
        sch.style_plot_axes(ax)
        ax.set_position(ax.get_position())
        t = np.arange(len(fd))
        cols = ['#0072B2', '#009E73', '#D55E00', '#56B4E9', '#E69F00', '#CC79A7']
        off = np.nanmax(np.abs(params)) * 1.2 + 1e-6
        for k in range(6):
            ax.plot(t, params[:, k] + off * 1.3, lw=0.35, color=cols[k])
        sc = off / max(np.nanmax(fd), thr * 1.5)
        ax.plot(t, fd * sc - off * 0.3, lw=0.4, color='#333333')
        ax.axhline(thr * sc - off * 0.3, lw=0.35, ls=(0, (2, 1)), color=RED)
        cens = np.where(fd > thr)[0]
        ax.scatter(cens, np.full(cens.size, -off * 0.62), s=0.6, color=RED, marker='|', lw=0.35)
        ax.set_xlim(0, len(fd) - 1)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['left'].set_visible(False)
        ax.text(0.03, 0.97, '6 rigid params', transform=ax.transAxes, fontsize=4.6, va='top', color='#555555')
        ax.text(0.03, 0.36, 'FD', transform=ax.transAxes, fontsize=4.6, va='top', color='#333333')
    return draw


def carpet_plot(carpet):
    def draw(ax):
        ax.imshow(carpet, cmap='gray', aspect='auto', vmin=-2, vmax=2, interpolation='nearest')
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.35)
    return draw


def synthetic_motion(n=300, seed=4):
    rng = np.random.default_rng(seed)
    drift = np.cumsum(rng.normal(0, 0.012, (n, 6)), axis=0)
    drift = sch.ndi.gaussian_filter1d(drift, 3, axis=0)
    for s in (90, 91, 210):
        drift[s:] += rng.normal(0, 0.12, 6)
    d = np.diff(drift, axis=0, prepend=drift[:1])
    fd = np.abs(d[:, :3]).sum(1) + 50 * np.deg2rad(np.abs(d[:, 3:])).sum(1)
    return drift, fd, 0.2


def fd_from(params, radius=50.0):
    """Framewise displacement (Power et al. 2012); params: 3 rotations (rad) + 3 translations (mm), AFNI-like
    columns of the EDNiX .1D (order detected from the magnitude)."""
    p = np.asarray(params, float)
    d = np.diff(p, axis=0, prepend=p[:1])
    rot_first = np.nanmax(np.abs(p[:, :3])) < np.nanmax(np.abs(p[:, 3:]))
    rot, tr = (d[:, :3], d[:, 3:]) if rot_first else (d[:, 3:], d[:, :3])
    if np.nanmax(np.abs(rot)) > 0.2:             # degrees
        rot = np.deg2rad(rot)
    return np.abs(tr).sum(1) + radius * np.abs(rot).sum(1)


# ── images ───────────────────────────────────────────────────────────────────
def thumbs_synthetic():
    sb = sch.SyntheticBrain()
    asp = FH / FW
    up = lambda a: sch.upsample(a, 3, order=1)[:sb.shape[0], :sb.shape[1]]
    upn = lambda a: sch.upsample(a, 3, order=0)[:sb.shape[0], :sb.shape[1]]
    pad3 = lambda a: np.pad(a, ((0, sb.shape[0] - a.shape[0]), (0, sb.shape[1] - a.shape[1])))
    epi_d = pad3(up(sb.epi(distort=7.0, dropout=0.65)))
    epi_u = pad3(up(sb.epi(distort=0.0, dropout=0.65)))
    rot, sh = 12.0, (5, -6)
    raw = np.abs(sb.noise(sb.rotate(epi_d, rot, sh), 0.03))
    by, bx = sch.fit_crop(epi_u, asp, sb.brain, pad=0.10)
    c = lambda a: sch.crop_pad(a, by, bx)
    ry, rx = sch.fit_crop(raw, asp, raw > 0.15, pad=0.10)
    cr = lambda a: sch.crop_pad(a, ry, rx)
    n4r = sb.rotate(epi_d, rot, sh)
    n4u = sb.rotate(epi_u, rot, sh)
    dist_edge = sb.smooth(n4r, 2) > 0.18
    other = sb.rotate(epi_u, 3.0, (2, -2))
    wm_edge = (sb.smooth(sb.wm, 1.0) > 0.5).astype(float)
    lab = np.where(sb.brain, sb.lobes(), 0)
    lab3 = np.where(pad3(upn((sch.block_mean(sb.brain.astype(float), 3) > 0.5).astype(float))) > 0,
                    pad3(upn(_mode3(lab))), 0)
    params, fd, thr = synthetic_motion()
    carpet = _synthetic_carpet()
    tm = (sb.smooth(sb.brain.astype(float), 1) > 0.5).astype(float)
    vmax = np.percentile(epi_u[sb.brain], 99.5)
    return dict(
        raw=Thumb(cr(raw), vmin=0, vmax=vmax, synthetic=True),
        motion=Thumb(plot=motion_plot(params, fd, thr), synthetic=True),
        sdc=Thumb(cr(n4u), vmin=0, vmax=vmax, contours=[(cr(dist_edge), RED, 0.5, (0, (1.5, 1)))], synthetic=True),
        inref=Thumb(c(epi_u), vmin=0, vmax=vmax, contours=[(c(sch.ndi.binary_fill_holes(sb.smooth(other, 2) > 0.18)), YEL, 0.45)],
                    synthetic=True),
        acpc=Thumb(c(epi_u), vmin=0, vmax=vmax, contours=[(c(wm_edge), YEL, 0.45)], synthetic=True),
        atlas=Thumb(c(epi_u), vmin=0, vmax=vmax, labels=c(lab3), label_alpha=0.6, synthetic=True),
        carpet=Thumb(plot=carpet_plot(carpet), synthetic=True),
        tmpl=Thumb(c(epi_u), vmin=0, vmax=vmax, contours=[(c(tm), YEL, 0.45)], synthetic=True),
    ), sb.source


def _mode3(lab):
    """Label map at 3 mm (label at the centre of each 3 x 3 block)."""
    H, Wd = (lab.shape[0] // 3) * 3, (lab.shape[1] // 3) * 3
    return lab[1:H:3, 1:Wd:3]


def _synthetic_carpet(nvox=160, n=300, seed=5):
    rng = np.random.default_rng(seed)
    nets = rng.normal(0, 1, (4, n))
    nets = sch.ndi.gaussian_filter1d(nets, 4, axis=1)
    w = np.zeros((nvox, 4))
    w[np.arange(nvox), np.repeat(np.arange(4), nvox // 4)] = 1
    x = w @ nets + rng.normal(0, 0.45, (nvox, n))
    x = (x - x.mean(1, keepdims=True)) / x.std(1, keepdims=True)
    return x


def thumbs_real(files):
    V = sch.Volume
    asp = FH / FW

    def sl(vol, support, contours=(), labels=None, pad=0.10, head=False):
        i = vol.sag_index(support > 0.5, offset_mm=-0.15 * np.ptp(np.where(support > 0.5)[0]) * vol.zooms[0])
        img = vol.sagittal(i)
        sup = img > np.percentile(img[img > 0], 25) if head else np.rot90((support > 0.5)[i])
        sy, sx = sch.fit_crop(img, asp * vol.zooms[1] / vol.zooms[2], sup, pad=pad)
        cp = lambda a: sch.crop_pad(a, sy, sx)
        cs = [(cp(np.rot90(m[i])), col, lw) + tuple(rest) for m, col, lw, *rest in contours]
        lab = cp(np.rot90(labels[i])) if labels is not None else None
        hi = np.percentile(img[img > 0], 99.5)
        return Thumb(cp(img), vmin=0, vmax=hi, contours=cs, labels=lab, aspect=vol.sag_aspect)

    bold = V(files['bold'], frame='mean')
    n4 = V(files['n4'])
    unw = V(files['unwarped'])
    inref = V(files['inref'])
    mean = V(files['mean_acpc'])
    anat = V(files['anat_func']).resample_like(mean, order=1)
    atlas = V(files['atlas'], frame=ATLAS_LEVEL_INDEX).resample_like(mean)
    mtpl = V(files['mean_tmpl'])
    tmpl = V(files['tmpl']).resample_like(mtpl, order=1)
    brain = atlas.data > 0
    thr = lambda v: v.data > np.percentile(v.data[v.data > 0], 35)
    a_wm = anat.data > np.percentile(anat.data[anat.data > 0], 70)
    params = np.loadtxt(files['motion'])
    fd = fd_from(params)
    tr = params[:, 3:] if np.nanmax(np.abs(params[:, :3])) < np.nanmax(np.abs(params[:, 3:])) else params[:, :3]
    ro = params[:, :3] if tr is params[:, 3:] else params[:, 3:]
    disp = np.column_stack([tr, ro])
    import nibabel as nib
    r4 = np.asanyarray(nib.as_closest_canonical(nib.load(files['residual'])).dataobj)
    vox = r4[brain] if r4.shape[:3] == brain.shape else r4[r4.std(-1) > 0]
    rng = np.random.default_rng(0)
    vox = vox[rng.choice(len(vox), min(300, len(vox)), replace=False)]
    vox = (vox - vox.mean(1, keepdims=True)) / (vox.std(1, keepdims=True) + 1e-9)
    return dict(
        raw=sl(bold, thr(bold), head=True),
        motion=Thumb(plot=motion_plot(disp, fd, 0.2)),
        sdc=sl(unw, thr(unw), contours=[(thr(n4.resample_like(unw, order=1)), RED, 0.5, (0, (1.5, 1)))], head=True),
        inref=sl(inref, thr(inref), head=True),
        acpc=sl(mean, brain, contours=[(a_wm, YEL, 0.45)]),
        atlas=sl(mean, brain, labels=atlas.data),
        carpet=Thumb(plot=carpet_plot(vox)),
        tmpl=sl(mtpl, thr(tmpl), contours=[(thr(tmpl), YEL, 0.45)]),
    ), files['mean_acpc']


def get_thumbs(force_synthetic=False):
    if not force_synthetic:
        files, ses = sch.find_session(FUNC_BIDS, PATTERNS)
        if files:
            try:
                th, _ = thumbs_real(files)
                print(f'  [Supp 2] images: {ses}')
                return th, ses
            except Exception as e:                                     # noqa: BLE001
                print(f'  [Supp 2] [WARN] real images failed ({e}) -- synthetic brain used')
        else:
            miss = sch.missing_report(FUNC_BIDS, PATTERNS) if os.path.isdir(FUNC_BIDS) else [FUNC_BIDS]
            print('  [Supp 2] [WARN] no complete session -- synthetic brain used. Missing: ' + '; '.join(miss))
    th, src = thumbs_synthetic()
    return th, 'synthetic: ' + src


# ── figure ───────────────────────────────────────────────────────────────────
def panel_a(cv, th, y0):
    cv.panel('a', 1.0, y0, 'Example run: image after each main step, by coordinate space')
    frames = [
        dict(thumb=th['raw'], title='Raw EPI', step=1, op='prepro', space='func',
             caption='despiked, slice timing'),
        dict(thumb=th['motion'], title='Motion', step=1, op='reg', space='func',
             caption='BOLDRigid, FD, censoring', link=dict(kind='flow')),
        dict(thumb=th['sdc'], title='Distortion', step=2, op='prepro', space='func',
             caption='field map / reverse PE', link=dict(kind='flow')),
        dict(thumb=th['inref'], title='Reference run', step=2, op='reg', space='func',
             caption='runs on one run', link=dict(kind='reg', label='Rigid')),
        dict(thumb=th['acpc'], title='To anatomy', step=5, op='reg', space='funcacpc',
             caption='anat. WM edges', link=dict(kind='reg', label='SyN')),
        dict(thumb=th['atlas'], title='Atlas in EPI', step=5, op='apply', space='funcacpc',
             caption='atlases, tissue masks', link=dict(kind='apply', label='T$^{-1}$')),
        dict(thumb=th['carpet'], title='Denoised', step=7, op='analysis', space='funcacpc',
             caption='voxels x time', link=dict(kind='flow')),
        dict(thumb=th['tmpl'], title='Template', step=9, op='apply', space='template',
             caption='template edges', link=dict(kind='apply', label='T')),
    ]
    boxes = sch.draw_strip(cv, frames, 3.0, y0 + 5.0, sch.W - 5.0, FW, FH)
    return boxes[0][1] + FH + 5.6


def panel_b(cv, y0):
    cv.panel('b', 1.0, y0, 'Processing steps of the functional module (numbers: EDNiX Skip_step)')
    y = y0 + 6.0
    colw, gap = 43.2, 2.0
    xs = [3.0 + i * (colw + gap) for i in range(4)]
    blocks = [
        ('Each run', 'functional space', [
            (1, 'Volume clean-up', ['first volumes removed (3dTcat),',
                                    'despiking (3dDespike), slice timing',
                                    '(3dTshift, wsinc9), outliers'], 'prepro', 'AFNI'),
            (1, 'Motion estimation', ['every volume to the run mean, BOLDRigid;',
                                      'framewise displacement scaled to the',
                                      'brain size; censoring of the volumes'], 'reg', 'ANTsPy motion_correction'),
            (2, 'Distortion correction', ['field map (FUGUE) or reverse phase',
                                          'encoding (TOPUP) on the run mean;',
                                          'kept as a displacement field'], 'prepro', 'FSL'),
            (2, 'Runs to reference run', ['corrected run means, BOLDRigid'], 'reg', 'ANTsPy'),
        ]),
        ('Masks and anatomy', 'functional AC-PC space', [
            (3, 'Functional brain mask', ['anatomical mask brought in, or a',
                                          'skull-stripping method of the EPI'], 'mask', 'Method_mask_func'),
            ('itk_1', 'Mask check', ['ITK-SNAP'], 'check', None),
            (5, 'EPI to anatomy', ['rigid (MI) of the mean EPI to the AC-PC',
                                   'anatomy, then SyN restricted to the',
                                   'phase-encoding axis'], 'reg', 'ANTsPy restrict_transformation'),
            (5, 'One-shot resampling', ['per volume: motion, distortion, run to',
                                        'reference, rigid: one interpolation;',
                                        'anatomy, atlases, masks brought in'], 'apply', 'ANTsPy apply_transforms'),
        ]),
        ('Signal', 'functional AC-PC space', [
            (6, 'ICA denoising (optional)', ['MELODIC components removed'], 'analysis', 'FSL melodic, fsl_regfilt'),
            (7, 'Nuisance regression', ['motion (12 / 24), WM and CSF',
                                        'components (CompCor), band-pass,',
                                        'censored volumes'], 'analysis', 'AFNI 3dTproject / 3dDeconvolve'),
            (10, 'Correlation matrices', ['regional time series of each atlas'], 'analysis', 'nilearn'),
            (11, 'Seed-based maps', ['one seed per atlas region'], 'analysis', 'nilearn'),
        ]),
        ('Other spaces', 'anatomy, template', [
            (8, 'EPI to anatomy space', ['SyN applied to the series'], 'apply', 'ANTsPy'),
            (9, 'EPI to template', ['SyN and anatomy-to-template chain', 'composed: one interpolation'],
             'apply', 'ANTsPy'),
            ('12-14', 'Quality control', ['tSNR, FD, DVARS, registration and',
                                          'seed-map specificity'], 'analysis', None),
            ('itk_3', 'Template check', ['EPI over the template'], 'check', 'ITK-SNAP'),
        ]),
    ]
    ymax = y
    for x, (name, sub, cards) in zip(xs, blocks):
        cv.block_header(x, y, colw, name, sub)
        yy = y + 6.4
        for step, title, lines, op, tool in cards:
            h = cv.card(x, yy, colw, step, title, lines, op, tool=tool)
            yy += h + 1.3
        ymax = max(ymax, yy)
    return ymax


def panel_c(cv, y0):
    """Composition of the transforms: what is estimated, what is resampled, and how many times."""
    cv.panel('c', 1.0, y0, 'Transform chain from each EPI volume to the template: two interpolations only')
    nodes = [('Volume j', 'func'), ('Run mean', 'func'), ('Undistorted', 'func'), ('Reference run', 'func'),
             ('Functional\nAC-PC grid', 'funcacpc'), ('Anatomy\n(AC-PC)', 'native'), ('Template', 'template')]
    edges = ['motion$_j$', 'SDC', 'to ref. run', 'rigid', 'SyN (PE)', 'to template']
    yi = y0 + 7.4                                      # inverse chain, above the nodes
    y = yi + 3.6
    nw, nh = 16.4, 7.0
    gap = (sch.W - 6.0 - len(nodes) * nw) / (len(nodes) - 1)
    xs = [3.0 + i * (nw + gap) for i in range(len(nodes))]
    for x, (lab, sp) in zip(xs, nodes):
        s = SPACES[sp]
        cv.rect(x, y, nw, nh, fc=s['fill'], ec=s['ink'], lw=0.45, r=1.4, z=3)
        cv.text(x + nw / 2, y + nh / 2, lab, size=sch.FS_TEXT - 0.2, ha='center', va='center', color=s['ink'],
                weight='bold', z=5)
    for i, e in enumerate(edges):
        x0, x1 = xs[i] + nw + 0.6, xs[i + 1] - 0.6
        cv.transition(x0, x1, y + nh / 2 + 0.8, kind='reg', label=e, label_dy=-0.9)
    # inverse: atlases into the EPI
    col = OPS['apply']['color']
    xa, xb = xs[6] + nw / 2, xs[4] + nw / 2
    cv.ax.plot([xa, xa, xb], [y - 0.4, yi, yi], color=col, lw=0.7, ls=(0, (2.2, 1.2)), zorder=4)
    cv.arrow(xb, yi, xb, y - 0.4, color=col, lw=0.7, head=3.2)
    cv.text(xb - 1.5, yi, 'inverse chain: atlases and tissue masks into the EPI (step 5)', size=sch.FS_SMALL,
            ha='right', va='center', color=col)
    # interpolation brackets
    yb = y + nh + 2.0
    for (i0, i1, txt, step) in [(0, 4, 'interpolation 1: each volume, into the functional AC-PC grid', 5),
                                (4, 6, 'interpolation 2: anatomy / template', '8-9')]:
        xa, xb = xs[i0] + 1.0 + (1.2 if i0 else 0), xs[i1] + nw - 1.0 - (1.2 if i1 == 4 else 0)
        cv.ax.plot([xa, xa, xb, xb], [yb - 0.8, yb, yb, yb - 0.8], color=col, lw=0.6, zorder=4)
        tw = 0.37 * sch.FS_SMALL * 0.3528 * 2 * len(txt) * 0.5
        xm = (xa + xb) / 2
        xr = cv.badge(xm - tw / 2 - 2.4, yb + 2.3, step, col, r=1.35, size=4.4)
        cv.text(xr + 1.0, yb + 2.3, txt, size=sch.FS_SMALL, va='center', color=col)
    return yb + 4.0


def make(out=OUT, force_synthetic=False):
    th, src = get_thumbs(force_synthetic)
    cv = sch.Canvas(220.0)
    ya = panel_a(cv, th, 1.0)
    yb = panel_b(cv, ya + 1.5)
    yc = panel_c(cv, yb + 1.5)
    yl = yc + 4.0
    cv.ops_legend(3.0, yl)
    yl2 = yl + 3.6
    cv.arrow_legend(3.0, yl2)
    cv.trim(yl2 + 3.0)
    print(f'  [Supp 2] image source: {src}')
    return cv.save(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--synthetic', action='store_true', help='synthetic images even if the data are found')
    ap.add_argument('--out', default=None, help='output folder (default <PAPER_DIR>/Supp_Fig_02)')
    a = ap.parse_args()
    make(os.path.join(a.out, 'Supp_Fig_02') if a.out else OUT, a.synthetic)


if __name__ == '__main__':
    main()
