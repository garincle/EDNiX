"""
Supplementary Figure 3 -- spatial registration in the FDG-PET module of EDNiX, with the CT path
(modalities/pet, modalities/ct; scenarios of Illustration/pet4.svg that use an anatomical MRI).

  a  registration of the motion-corrected PET session mean to the anatomy, by scenario:
       1v1   hybrid PET-MR: the AC-PC matrix of the T1w is reused, nothing is estimated
       1v2a  PET and MRI in different sessions: multi-start rigid registration to a synthetic
             FDG image of the T1w (coreg_target='synthetic')
       1v2b  PET/CT + MRI: PET -> CT (same scanner: identity) composed with a multi-start rigid
             CT -> T1w registration estimated by the CT module
     then the T1w -> template chain of the anatomical module (Supplementary Fig. 1)
  b  one session in the PET space: motion correction, atlases / masks / tissue fractions brought
     into the PET, parametric images computed there, export by one interpolation
  c  the steps of the CT and PET modules (numbers = Skip_step of the launchers)

Images: processed sessions of PET_BIDS (1v2a) and PETCT_BIDS (1v2b) when found, else a synthetic
brain (MNI ICBM152, frames tagged 'synthetic'); each frame falls back on its own.
Output: <PAPER_DIR>/Supp_Fig_03/Supp_Fig_03.pdf|png
Former Supplementary/Supp_Fig_03_*.py.
Run: python Exemples/Study_EDNiX/analysis/run_paper.py supp3, or python Exemples/Study_EDNiX/analysis/library/schematic_pet_ct.py [--synthetic] [--out DIR]
"""
import os
import sys
import glob
import argparse

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), *[os.pardir] * 4))
for _path in (_REPO, '/home/cgarin/PycharmProjects/EDNiX/'):
    if os.path.isdir(_path) and _path not in sys.path:
        sys.path.insert(0, _path)

from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis.library import pipeline_schematic as sch
from Exemples.Study_EDNiX.analysis.library.pipeline_schematic import Thumb, OPS, SPACES

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_03', 'Supp_Fig_03')

PET_BIDS = '/scratch2/EDNiX/Macaque/BIDS_Tremblay'       # 1v2a, coreg_target='synthetic'
PETCT_BIDS = '/scratch2/EDNiX/Human/BIDS_Merida'         # 1v2b (CERMEP, PET/CT + T1w)
REFERENCE = 'EDNiX'
ATLAS = 'EDNIxCSC'
ATLAS_LEVEL_INDEX = 1
T = 'T1w'
PET_PATTERNS = {
    'raw':      'pet/*_pet.nii.gz',
    'mean':     'pet/orig/preprocessing/*_desc-MoCo-between-mean_pet.nii.gz',
    'synth':    f'pet/orig/preprocessing/forMC/*_desc-synthetic_{T}.nii.gz',
    'pet_acpc': 'pet/orig/preprocessing/forMC/*_space-acpc_desc-mean_pet.nii.gz',
    'pet_tmpl': f'pet/orig/preprocessing/forMC/*_space-{REFERENCE}_desc-mean_pet.nii.gz',
    'atlas':    f'pet/orig/labels/*_space-orig_seg-{ATLAS}_dseg.nii.gz',
    'gm':       'pet/orig/masks/*_space-orig_label-GM_probseg.nii.gz',
}
PET_OPTIONAL = {
    'motion':   'pet/orig/preprocessing/forMC/*_desc-MoCo-within_motion.tsv',
    'cropbox':  'pet/orig/preprocessing/forMC/*_desc-cropbox_mask.nii.gz',
    'suv':      'pet/orig/postprocessing/*_desc-SUV-mean_pet.nii.gz',
}
CT_PATTERNS = {
    'ct':       'ct/orig/preprocessing/forMC/*_desc-clipped_ct.nii.gz',
    'mean':     'pet/orig/preprocessing/*_desc-MoCo-between-mean_pet.nii.gz',
    'pet_acpc': 'pet/orig/preprocessing/forMC/*_space-acpc_desc-mean_pet.nii.gz',
}
ANAT_IN_SUBJECT = f'ses-*/anat/preprocessing/*_space-acpc_desc-cropped_{T}.nii.gz'

FW, FH = 17.6, 15.8
PET_CMAP = 'afmhot'
CYAN, YEL, GRN = '#22D3EE', '#FFD400', '#3DDC97'


# ── synthetic images ─────────────────────────────────────────────────────────
class SynthPET:
    """All synthetic frames of the figure, from one SyntheticBrain."""

    def __init__(self):
        sb = self.sb = sch.SyntheticBrain()
        self.asp = FH / FW
        up = lambda a, k: np.pad(sch.upsample(a, k, order=1), 0)[:sb.shape[0], :sb.shape[1]]
        self.pad = lambda a: np.pad(a, ((0, max(0, sb.shape[0] - a.shape[0])), (0, max(0, sb.shape[1] - a.shape[1]))))[:sb.shape[0], :sb.shape[1]]
        fdg = self.pad(up(sb.fdg(), 2))
        self.rot, self.sh = 11.0, (7, -9)
        self.pet_space = lambda a, order=1: sb.rotate(a, self.rot, self.sh, order=order)
        rng = np.random.default_rng(11)
        self.fdg = fdg
        self.pet = np.clip(self.pet_space(fdg) + sch.ndi.gaussian_filter(rng.normal(0, 0.05, sb.shape), 1.2), 0, None)
        self.frame = np.clip(self.pet_space(0.35 * fdg + 0.25 * sb.smooth(sb.head.astype(float), 6))
                             + sch.ndi.gaussian_filter(rng.normal(0, 0.12, sb.shape), 1.0), 0, None)
        self.ct = self.pet_space(sb.ct(), order=1)
        self.ct[self.ct < -1000] = -1000
        self.head_pet = self.pet_space(sb.head.astype(float), order=0) > 0.5
        self.brain_pet = self.pet_space(sb.brain.astype(float), order=0) > 0.5
        self.t1 = np.abs(sb.noise(sb.t1_head(), 0.025))
        self.target = sb.synthetic_target()
        py, px = sch.fit_crop(self.pet, self.asp, self.head_pet, pad=0.06)
        self.cp = lambda a: sch.crop_pad(a, py, px)
        ay, ax = sch.fit_crop(self.t1, self.asp, sb.head, pad=0.05)
        self.ca = lambda a: sch.crop_pad(a, ay, ax)
        by, bx = sch.fit_crop(self.t1, self.asp, sb.brain, pad=0.10)
        self.cb = lambda a: sch.crop_pad(a, by, bx)

    def thumbs(self):
        sb, cp, ca, cb = self.sb, self.cp, self.ca, self.cb
        vmax = np.percentile(self.pet[self.brain_pet], 99.5)
        lab = np.where(self.brain_pet, self.pet_space(sb.lobes().astype(float), order=0), 0)
        lab2 = sch.upsample(lab[::2, ::2], 2, order=0)[:sb.shape[0], :sb.shape[1]]
        lab2 = self.pad(lab2)
        gm2 = self.pad(sch.upsample(sch.block_mean(self.pet_space(sb.gm), 2), 2, order=0))
        box = np.zeros(sb.shape)
        ys, xs = np.where(self.head_pet)
        box[ys.min() - 4:ys.max() + 4, xs.min() - 4:xs.max() + 4] = 1
        suv = self.pet * 2.6 / vmax
        tm = (sb.smooth(sb.brain.astype(float), 1) > 0.5).astype(float)
        return dict(
            pet=Thumb(cp(self.pet), cmap=PET_CMAP, vmin=0, vmax=vmax, synthetic=True),
            frame=Thumb(cp(self.frame), cmap=PET_CMAP, vmin=0, vmax=np.percentile(self.frame, 99.7), synthetic=True),
            motion=Thumb(plot=moco_plot(*synthetic_frames_motion()), synthetic=True),
            crop=Thumb(cp(self.pet), cmap=PET_CMAP, vmin=0, vmax=vmax, contours=[(cp(box), CYAN, 0.5, (0, (2, 1)))],
                       synthetic=True),
            ct=Thumb(cp(self.ct), cmap='gray', vmin=-100, vmax=250, synthetic=True,
                     contours=[(cp(self.head_pet), CYAN, 0.4, (0, (1.5, 1)))]),
            target=Thumb(cb(self.target), cmap=PET_CMAP, vmin=0, vmax=1.0, synthetic=True),
            t1=Thumb(ca(self.t1), vmin=0, vmax=1.0, synthetic=True),
            pet_tmpl=Thumb(cb(self.fdg), cmap=PET_CMAP, vmin=0, vmax=np.percentile(self.fdg[sb.brain], 99.5),
                           contours=[(cb(tm), CYAN, 0.45)], synthetic=True),
            atlas=Thumb(cp(self.pet), cmap='gray', vmin=0, vmax=vmax, labels=cp(lab2), label_alpha=0.6, synthetic=True),
            gm=Thumb(cp(gm2), cmap='viridis', vmin=0, vmax=1, synthetic=True),
            suv=Thumb(cp(np.where(self.brain_pet, suv, 0)), cmap='inferno', vmin=0, vmax=2.6, synthetic=True),
        )


def synthetic_frames_motion(n=26, seed=9):
    rng = np.random.default_rng(seed)
    t = np.cumsum(np.r_[np.full(6, 0.5), np.full(8, 2.0), np.full(12, 5.0)])
    p = np.cumsum(rng.normal(0, 0.12, (n, 6)), axis=0)
    p[17:, :3] += [0.8, -0.5, 0.6]
    d = np.abs(p[:, :3]).sum(1) * 0.6 + np.abs(p[:, 3:]).sum(1) * 0.4
    return t, p, d


def moco_plot(t, params, disp):
    """Rigid parameters of the frames of one run and mean displacement (sphere), against frame start."""
    def draw(ax):
        sch.style_plot_axes(ax)
        cols = ['#0072B2', '#009E73', '#D55E00', '#56B4E9', '#E69F00', '#CC79A7']
        off = np.nanmax(np.abs(params)) * 1.1 + 1e-6
        for k in range(params.shape[1]):
            ax.plot(t, params[:, k] + off * 1.4, lw=0.35, color=cols[k % 6], marker='o', ms=0.6)
        sc = off / max(np.nanmax(disp), 1e-6)
        ax.plot(t, disp * sc - off * 0.4, lw=0.45, color='#333333', marker='o', ms=0.7)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['left'].set_visible(False)
        ax.set_xlim(t.min() - 1, t.max() + 1)
        ax.text(0.03, 0.97, 'frames: rigid', transform=ax.transAxes, fontsize=4.6, va='top', color='#555555')
        ax.text(0.03, 0.38, 'displ.', transform=ax.transAxes, fontsize=4.6, va='top', color='#333333')
    return draw


# ── real images (each frame falls back on the synthetic one) ─────────────────
def thumbs_real(syn):
    V = sch.Volume
    asp = FH / FW
    th, used = {}, []

    def sl(vol, support, cmap='gray', contours=(), labels=None, pad=0.08, head=False, vmax=None, vmin=0):
        s = support > 0.5
        i = vol.sag_index(s, offset_mm=-0.15 * np.ptp(np.where(s)[0]) * vol.zooms[0])
        img = vol.sagittal(i)
        sup = img > np.percentile(img[np.isfinite(img)], 60) if head else np.rot90(s[i])
        sy, sx = sch.fit_crop(img, asp * vol.zooms[1] / vol.zooms[2], sup, pad=pad)
        cp = lambda a: sch.crop_pad(a, sy, sx)
        cs = [(cp(np.rot90(m[i])), col, lw) + tuple(rest) for m, col, lw, *rest in contours]
        lab = cp(np.rot90(labels[i])) if labels is not None else None
        hi = vmax if vmax is not None else np.percentile(img[img > 0], 99.5)
        return Thumb(cp(img), cmap=cmap, vmin=vmin, vmax=hi, contours=cs, labels=lab, aspect=vol.sag_aspect)

    files, ses = sch.find_session(PET_BIDS, PET_PATTERNS)
    if files:
        try:
            opt = {k: (sorted(glob.glob(os.path.join(ses, p))) or [None])[0] for k, p in PET_OPTIONAL.items()}
            mean = V(files['mean'])
            atlas = V(files['atlas'], frame=ATLAS_LEVEL_INDEX).resample_like(mean)
            brain = atlas.data > 0
            vmax = np.percentile(mean.data[brain], 99.5)
            th['pet'] = sl(mean, brain, cmap=PET_CMAP, head=True, vmax=vmax)
            raw = V(files['raw'], frame=0)
            th['frame'] = sl(raw, raw.data > np.percentile(raw.data, 90), cmap=PET_CMAP, head=True)
            th['atlas'] = sl(mean, brain, labels=atlas.data, vmax=vmax)
            gm = V(files['gm']).resample_like(mean, order=1)
            th['gm'] = sl(gm, brain, cmap='viridis', vmax=1.0)
            syn_t = V(files['synth'])
            th['target'] = sl(syn_t, syn_t.data > 0.02, cmap=PET_CMAP, vmax=1.0)
            pt = V(files['pet_tmpl'])
            th['pet_tmpl'] = sl(pt, pt.data > np.percentile(pt.data[pt.data > 0], 40), cmap=PET_CMAP)
            if opt['cropbox']:
                box = V(opt['cropbox']).resample_like(mean)
                th['crop'] = sl(mean, brain, cmap=PET_CMAP, head=True, vmax=vmax,
                                contours=[(box.data, CYAN, 0.5, (0, (2, 1)))])
            if opt['suv']:
                suv = V(opt['suv'])
                th['suv'] = sl(suv, V(files['atlas'], frame=ATLAS_LEVEL_INDEX).resample_like(suv).data,
                               cmap='inferno', vmax=np.percentile(suv.data[suv.data > 0], 99))
            if opt['motion']:
                import pandas as pd
                df = pd.read_csv(opt['motion'], sep='\t')
                t = df['start_s'].to_numpy() / 60 if df['start_s'].notna().all() else df['frame'].to_numpy()
                th['motion'] = Thumb(plot=moco_plot(t, df[['tx_mm', 'ty_mm', 'tz_mm', 'rx_deg', 'ry_deg', 'rz_deg']]
                                                    .to_numpy(), df['mean_displacement_mm'].to_numpy()))
            anat = sorted(glob.glob(os.path.join(os.path.dirname(ses) if os.path.basename(ses).startswith('ses-')
                                                 else ses, ANAT_IN_SUBJECT)))
            if anat:
                a = V(anat[0])
                th['t1'] = sl(a, a.data > np.percentile(a.data[a.data > 0], 30), head=True)
            used.append(ses)
        except Exception as e:                                     # noqa: BLE001
            print(f'  [Supp 3] [WARN] PET images of {ses} failed ({e})')
    else:
        print('  [Supp 3] [WARN] no complete PET session in ' + PET_BIDS + ' -- synthetic frames. Missing: '
              + '; '.join(sch.missing_report(PET_BIDS, PET_PATTERNS) if os.path.isdir(PET_BIDS) else [PET_BIDS]))
    files, ses = sch.find_session(PETCT_BIDS, CT_PATTERNS)
    if files:
        try:
            ct = V(files['ct'])
            head = ct.data > -500
            th['ct'] = sl(ct, head, head=False, vmin=-100, vmax=250, pad=0.05)
            used.append(ses)
        except Exception as e:                                     # noqa: BLE001
            print(f'  [Supp 3] [WARN] CT images of {ses} failed ({e})')
    else:
        print('  [Supp 3] [WARN] no complete PET/CT session in ' + PETCT_BIDS + ' -- synthetic CT frame.')
    out = syn.thumbs()
    out.update(th)
    return out, used


def get_thumbs(force_synthetic=False):
    syn = SynthPET()
    if force_synthetic:
        return syn.thumbs(), []
    return thumbs_real(syn)


# ── figure ───────────────────────────────────────────────────────────────────
def panel_a(cv, th, y0):
    cv.panel('a', 1.0, y0, 'Registration of the PET to the anatomy, by acquisition scenario')
    reg, app = OPS['reg']['color'], OPS['apply']['color']
    x_lab, x0 = 3.0, 27.0
    cols = [('pet', 'PET space'), ('ct', 'CT space'), ('native', 'Native AC-PC space'), ('native', None),
            ('template', 'Reference template')]
    n = len(cols)
    gap = (sch.W - 2.0 - x0 - n * FW) / (n - 1)
    xs = [x0 + i * (FW + gap) for i in range(n)]
    ytop = y0 + 5.5
    lane_h = [8.0, FH + 8.4, FH + 8.4]
    ylanes = [ytop + 7.5]
    for h in lane_h[:-1]:
        ylanes.append(ylanes[-1] + h + 1.4)
    ybot = ylanes[-1] + lane_h[-1]
    # columns (space bands): AC-PC spans the synthetic target and the T1w
    bands = [(0, 0, 'pet'), (1, 1, 'ct'), (2, 3, 'native'), (4, 4, 'template')]
    for i0, i1, sp in bands:
        cv.space_band(xs[i0] - 1.6, ytop, xs[i1] + FW - xs[i0] + 3.2, ybot - ytop + 1.0, sp)
    lanes = [
        ('1v1', 'hybrid PET-MR', 'same scanner session'),
        ('1v2a', 'PET and MRI', 'different sessions'),
        ('1v2b', 'PET/CT and MRI', 'CT links PET to MRI'),
    ]
    for (code, name, sub), y, h in zip(lanes, ylanes, lane_h):
        cv.text(x_lab, y + h / 2 - 2.2, code, size=sch.FS_TITLE + 0.5, weight='bold', va='center')
        cv.text(x_lab, y + h / 2 + 0.9, name, size=sch.FS_SMALL, va='center')
        cv.text(x_lab, y + h / 2 + 3.3, sub, size=sch.FS_SMALL, va='center', color=sch.SUBTLE, style='italic')
    # lane 1v1: pills (no image of a hybrid acquisition)
    y = ylanes[0]
    pill = dict(h=6.0)

    def node(i, yy, label, sp):
        s = SPACES[sp]
        cv.rect(xs[i], yy, FW, pill['h'], fc='white', ec=s['ink'], lw=0.5, r=1.2, z=3)
        cv.text(xs[i] + FW / 2, yy + pill['h'] / 2, label, size=sch.FS_SMALL, ha='center', va='center',
                color=s['ink'], weight='bold')
    ym = y + 1.0
    node(0, ym, 'PET mean', 'pet')
    node(3, ym, 'T1w', 'native')
    node(4, ym, 'Template', 'template')
    yc = ym + pill['h'] / 2
    cv.transition(xs[0] + FW + 0.4, xs[3] - 0.4, yc, kind='apply',
                  label='AC-PC matrix of the T1w reused, nothing estimated', label_dy=-0.9)
    cv.transition(xs[3] + FW + 0.4, xs[4] - 0.4, yc, kind='apply', label='T$_{anat}$', label_dy=-0.9)
    # lanes with images
    def frame(i, yy, key, title, caption=None):
        sch.draw_frame(cv, xs[i], yy, FW, FH, th[key], None, None, None, caption)
        if title:
            cv.text(xs[i] + FW / 2, yy - 0.8, title, size=sch.FS_SMALL, weight='bold', ha='center', va='bottom')
    # 1v2a
    y = ylanes[1] + 3.0
    frame(0, y, 'pet', 'PET mean', 'MoCo, cropped')
    frame(2, y, 'target', 'Synthetic FDG', 'GM 1, WM .35, CSF .05')
    frame(3, y, 't1', 'T1w', 'GM/WM/CSF masks')
    frame(4, y, 'pet_tmpl', 'PET in template', 'template edges')
    yc = y + FH / 2
    cv.transition(xs[0] + FW + 0.4, xs[2] - 0.4, yc, kind='reg', label='Rigid, multi-start', label_dy=-0.9)
    cv.text((xs[0] + FW + xs[2]) / 2, yc + 1.0, 'max. MI, brain\ncoverage ≥ 90%', size=sch.FS_SMALL,
            ha='center', va='top', color=reg)
    cv.arrow(xs[3] - 0.4, yc, xs[2] + FW + 0.4, yc, color=sch.FLOW, lw=0.6, head=3.0, ls=(0, (1.2, 1.0)))
    cv.text((xs[2] + FW + xs[3]) / 2, yc - 0.9, 'built\nfrom', size=sch.FS_SMALL, ha='center', va='bottom',
            color=sch.SUBTLE)
    cv.transition(xs[3] + FW + 0.4, xs[4] - 0.4, yc, kind='apply', label='T$_{anat}$', label_dy=-0.9)
    # 1v2b
    y = ylanes[2] + 3.0
    frame(0, y, 'pet', 'PET mean', None)
    frame(1, y, 'ct', 'CT', 'head mask > -500 HU')
    frame(3, y, 't1', 'T1w', None)
    frame(4, y, 'pet_tmpl', 'PET in template', None)
    yc = y + FH / 2
    cv.transition(xs[0] + FW + 0.4, xs[1] - 0.4, yc, kind='apply', label='identity', label_dy=-0.9)
    cv.text((xs[0] + FW + xs[1]) / 2, yc + 1.0, 'same scanner\n(else rigid)', size=sch.FS_SMALL, ha='center',
            va='top', color=sch.SUBTLE)
    cv.transition(xs[1] + FW + 0.4, xs[3] - 0.4, yc, kind='reg', label='Rigid, multi-start (CT module)',
                  label_dy=-0.9)
    cv.text((xs[1] + FW + xs[3]) / 2, yc + 1.0, 'MI on [-1000, 2000] HU, metric in the\nT1w brain + 10 mm (skull included)',
            size=sch.FS_SMALL, ha='center', va='top', color=reg)
    cv.transition(xs[3] + FW + 0.4, xs[4] - 0.4, yc, kind='apply', label='T$_{anat}$', label_dy=-0.9)
    cv.text(xs[4] + FW + 1.4, ybot + 2.0, 'T$_{anat}$: T1w -> template chain of the anatomical module '
            '(Supplementary Fig. 1); PET -> template = one composed transform',
            size=sch.FS_SMALL, ha='right', va='top', color=app)
    return ybot + 5.6


def panel_b(cv, th, y0):
    cv.panel('b', 1.0, y0, 'One PET session: everything is computed in the PET space')
    frames = [
        dict(thumb=th['frame'], title='Dynamic frames', step=1, op='prepro', space='pet', caption='reoriented'),
        dict(thumb=th['motion'], title='Motion', step=4, op='reg', space='pet', caption='frames, then runs',
             link=dict(kind='flow')),
        dict(thumb=th['crop'], title='Session mean', step=5, op='prepro', space='pet', caption='crop box (head)',
             link=dict(kind='flow')),
        dict(thumb=th['atlas'], title='Atlas in PET', step=6, op='apply', space='pet', caption='genericLabel',
             link=dict(kind='apply', label='T$^{-1}$')),
        dict(thumb=th['gm'], title='GM fraction', step=6, op='apply', space='pet', caption='for PVC',
             link=dict(kind='flow')),
        dict(thumb=th['suv'], title='SUV, ratio', step=7, op='analysis', space='pet', caption='kinetics, PVC',
             link=dict(kind='flow')),
        dict(thumb=th['pet_tmpl'], title='Export', step=10, op='apply', space='template',
             caption='one interpolation', link=dict(kind='apply', label='T')),
    ]
    boxes = sch.draw_strip(cv, frames, 3.0, y0 + 5.0, sch.W - 5.0, FW, FH)
    return boxes[0][1] + FH + 5.6


def panel_c(cv, y0):
    cv.panel('c', 1.0, y0, 'Steps of the CT and PET modules (numbers = Skip_step)')
    y = y0 + 6.0
    colw, gap = 43.2, 2.0
    xs = [3.0 + i * (colw + gap) for i in range(4)]
    blocks = [
        ('CT module', 'run before the PET (modalities/ct)', [
            (1, 'Reorientation', ['positions and header, as the anatomy'], 'prepro', None),
            (3, 'CT preparation', ['clipped to [-1000, 2000] HU; head mask',
                                   '> -500 HU (largest component, no',
                                   'bed); intracranial cavity filled'], 'mask', None),
            (3, 'CT to T1w (1v2b)', ['multi-start rigid, Mattes MI in the',
                                     'T1w brain + 10 mm; optional second',
                                     'pass in the 0-80 HU brain window'], 'reg', 'ANTsPy; registration record (JSON)'),
        ]),
        ('Each PET session', 'PET space (orig)', [
            (1, 'Reorientation', ['as the anatomy; frame timing, SIF'], 'prepro', None),
            (4, 'Within-run motion', ['frames rigid to the duration-weighted',
                                      'mean, two passes, in the crop box'], 'reg', 'ANTsPy'),
            (5, 'Between runs', ['run means rigid to their mean;',
                                 'session mean = registration image'], 'reg', 'ANTsPy'),
            ('2-3', 'PET study template', ['optional: common crop box'], 'prepro', None),
        ]),
        ('Registration', 'step 6, scenario of panel a', [
            (6, 'Multi-start rigid', ['starts: centres of mass, fields of view,',
                                      'AC-PC matrix, headers, rotation search,',
                                      'brain centroids; best MI kept if the',
                                      'brain is covered ≥ 90%'], 'reg', 'candidates.json'),
            (6, 'Synthetic FDG target', ['GM 1, WM 0.35, CSF 0.05 smoothed to',
                                         'the PET FWHM: the cortical ribbon,',
                                         'not the scalp, drives the alignment'], 'mask', 'coreg_target=\'synthetic\''),
            (6, 'Into the PET space', ['atlases (genericLabel), masks (linear,',
                                       '> 0.5), tissue fractions; boundary check'], 'apply', None),
        ]),
        ('PET-space analyses', 'then export', [
            (7, 'Parametric images', ['SUV, ratio to the brain mean, regional',
                                      'values of every atlas level'], 'analysis', None),
            (8, 'Kinetic models', ['SRTM, Logan (reference region)'], 'analysis', 'TPCCLIB'),
            ('pvc', 'Partial-volume correction', ['RBV / MG with the tissue fractions'], 'analysis', 'PETPVC or Python'),
            (10, 'Export', ['AC-PC and template spaces, composed',
                            'transforms, one interpolation'], 'apply', 'ANTsPy'),
        ]),
    ]
    ymax = y
    for x, (name, sub, cards) in zip(xs, blocks):
        cv.block_header(x, y, colw, name, sub, color=SPACES['ct']['ink'] if name == 'CT module' else sch.TEXT)
        yy = y + 6.4
        for step, title, lines, op, tool in cards:
            h = cv.card(x, yy, colw, step, title, lines, op, tool=tool)
            yy += h + 1.3
        ymax = max(ymax, yy)
    return ymax


def make(out=OUT, force_synthetic=False):
    th, used = get_thumbs(force_synthetic)
    cv = sch.Canvas(240.0)
    ya = panel_a(cv, th, 1.0)
    yb = panel_b(cv, th, ya + 1.0)
    yc = panel_c(cv, yb + 1.5)
    yl = yc + 3.0
    cv.ops_legend(3.0, yl)
    yl2 = yl + 3.6
    cv.arrow_legend(3.0, yl2)
    synthetic = any(getattr(t, 'synthetic', False) for t in th.values())
    cv.text(sch.W - 2.0, yl2, 'Images: ' + ('synthetic illustration' if not used else
                                            'example sessions' + (' + synthetic frames' if synthetic else '')),
            size=sch.FS_SMALL, ha='right', va='center', color=sch.SUBTLE)
    cv.trim(yl2 + 3.0)
    print('  [Supp 3] images: ' + (', '.join(used) if used else 'synthetic'))
    return cv.save(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--synthetic', action='store_true', help='synthetic images even if the data are found')
    ap.add_argument('--out', default=None, help='output folder (default <PAPER_DIR>/Supp_Fig_03)')
    a = ap.parse_args()
    make(os.path.join(a.out, 'Supp_Fig_03') if a.out else OUT, a.synthetic)


if __name__ == '__main__':
    main()
