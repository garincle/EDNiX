"""
Supplementary Figure 1 -- spatial registration in the anatomical MRI module of EDNiX
(modalities/anat: _loop1, studytemplate, _loop2, _loop3).

  a  one session through its coordinate spaces: raw T1w (scanner space) -> N4 -> brain mask
     pass 1 -> AC-PC alignment -> brain mask pass 2 -> registration to the template ->
     template atlas -> atlas brought back to the native AC-PC space
  b  the steps of the module (numbers = Skip_step of the anatomical launcher), by block
  c  the four registration designs (classical, study template, longitudinal, both)

Images: one processed session of ANAT_BIDS (first subject / session with every file), else a
synthetic brain (MNI ICBM152, frames tagged 'synthetic'). The schematic itself does not depend on
the data.
Output: <PAPER_DIR>/Supp_Fig_01/Supp_Fig_01.pdf|png
Run: python Exemples/Study_EDNiX/analysis/Supplementary/Supp_Fig_01_anat_pipeline.py [--synthetic] [--out DIR]
"""
import os
import sys
import glob
import argparse

import numpy as np

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), *[os.pardir] * 4))
for p in (_REPO, '/home/cgarin/PycharmProjects/EDNiX/'):
    if os.path.isdir(p) and p not in sys.path:
        sys.path.insert(0, p)

from Exemples.Study_EDNiX.analysis import paper_data as pdat
from Exemples.Study_EDNiX.analysis.Supplementary import pipeline_schematic as sch
from Exemples.Study_EDNiX.analysis.Supplementary.pipeline_schematic import Thumb, OPS, SPACES

OUT = os.path.join(pdat.PAPER_DIR, 'Supp_Fig_01', 'Supp_Fig_01')

# example session (first one with every file below)
ANAT_BIDS = '/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin'
SPECIES = 'Macaque'
REFERENCE = 'EDNiX'
ATLAS = 'EDNIxCSC'
ATLAS_LEVEL_INDEX = 1            # 4-D atlases: volume of level 2 (lobes)
T = 'T1w'
PATTERNS = {
    'raw':    f'anat/preprocessing/*_space-raw_desc-reorient_{T}.nii.gz',
    'n4':     f'anat/preprocessing/*_space-raw_desc-Bias_{T}.nii.gz',
    'mask1':  'anat/native/volumes/masks/*_desc-step1_mask.nii.gz',
    'acpc':   f'anat/preprocessing/*_space-acpc_desc-cropped_{T}.nii.gz',
    'mask2':  'anat/native/volumes/masks/*_space-acpc_mask.nii.gz',
    'ss':     f'anat/native/volumes/*_space-acpc_desc-SS_{T}.nii.gz',
    'warped': f'anat/templates/{REFERENCE}/volumes/native_space-{REFERENCE}_desc-SS_{T}.nii.gz',
    'atlas':  f'anat/native/volumes/labels/*_seg-{ATLAS}_dseg.nii.gz',
}

FW, FH = 17.6, 15.8            # frame size (mm)


def template_files():
    """Reference template of SPECIES in the atlas library (atlases/templatefeat.py layout)."""
    root = os.path.join(pdat.ATLAS_LIB, 'atlas')
    hits = sorted(glob.glob(os.path.join(root, '**', SPECIES, REFERENCE, 'volumes',
                                         f'{SPECIES}_space-acpc_desc-SS_{T}.nii.gz'), recursive=True))
    if not hits:
        return None
    vol = os.path.dirname(hits[0])
    out = dict(ss=hits[0], mask=os.path.join(vol, 'masks', f'{SPECIES}_mask.nii.gz'),
               atlas=os.path.join(vol, 'labels', f'{SPECIES}_seg-{ATLAS}_dseg.nii.gz'))
    return out if all(os.path.isfile(v) for v in out.values()) else None


# ── images ───────────────────────────────────────────────────────────────────
def thumbs_synthetic():
    sb = sch.SyntheticBrain()
    asp = FH / FW
    head = sb.t1_head()
    rot, shift = 26.0, (8, -12)
    raw = sb.noise(sb.rotate(head, rot, shift) * sb.rotate(sb.bias(0.55), 0), 0.035)
    raw = np.abs(raw)
    n4 = np.abs(sb.noise(sb.rotate(head, rot, shift), 0.03))
    m1 = sb.rotate(sb.smooth(sb.brain.astype(float), 6) > 0.32, rot, shift, order=0)
    sy, sx = sch.fit_crop(raw, asp, raw > 0.15, pad=0.04)
    c = lambda a: sch.crop_pad(a, sy, sx)
    hy, hx = sch.fit_crop(head, asp, head > 0.15, pad=0.04)
    ch = lambda a: sch.crop_pad(a, hy, hx)
    by, bx = sch.fit_crop(head, asp, sb.brain, pad=0.08)
    cb = lambda a: sch.crop_pad(a, by, bx)
    acpc = np.abs(sb.noise(head, 0.03))
    ss = sb.t1 * sb.brain
    warped = sb.smooth(ss, 1.0)
    tmpl_gm = sb.smooth(sb.wm, 1.5) > 0.5
    lab = sb.lobes()
    green, yellow = OPS['mask']['color'], '#FFD400'
    def line(deg, shift, sy_, sx_):
        (x0, y0), (x1, y1) = sb.acpc_line(deg, shift)
        return [((x0 - sx_.start, y0 - sy_.start), (x1 - sx_.start, y1 - sy_.start), yellow, 'AC-PC')]
    return dict(
        raw=Thumb(c(raw), vmin=0, vmax=1.15, synthetic=True, lines=line(rot, shift, sy, sx)),
        n4=Thumb(c(n4), vmin=0, vmax=1.0, synthetic=True),
        mask1=Thumb(c(n4), vmin=0, vmax=1.0, contours=[(c(m1), '#3DDC97', 0.6)], synthetic=True),
        acpc=Thumb(ch(acpc), vmin=0, vmax=1.0, synthetic=True, lines=line(0, (0, 0), hy, hx)),
        mask2=Thumb(cb(acpc), vmin=0, vmax=1.0, contours=[(cb(sb.brain), '#3DDC97', 0.6)], synthetic=True),
        warped=Thumb(cb(warped), vmin=0, vmax=1.0, contours=[(cb(tmpl_gm), yellow, 0.45)], synthetic=True),
        tatlas=Thumb(cb(ss), vmin=0, vmax=1.0, labels=cb(lab), synthetic=True),
        natlas=Thumb(cb(ss), vmin=0, vmax=1.0, labels=cb(lab), synthetic=True),
    ), sb.source


def thumbs_real(files, tf):
    V = sch.Volume
    asp = FH / FW
    raw, n4 = V(files['raw']), V(files['n4'])
    m1 = V(files['mask1']).resample_like(n4)
    acpc = V(files['acpc'])
    m2 = V(files['mask2']).resample_like(acpc)
    ss = V(files['ss'])
    atlas = V(files['atlas'], frame=ATLAS_LEVEL_INDEX).resample_like(ss)
    warped = V(files['warped'])
    tm = V(tf['mask']).resample_like(warped)
    tss = V(tf['ss'])
    tatlas = V(tf['atlas'], frame=ATLAS_LEVEL_INDEX).resample_like(tss)
    brain_w = np.ptp(np.where(m2.data > 0.5)[0]) * m2.zooms[0]
    off = -0.16 * brain_w                                      # parasagittal, left hemisphere

    def sl(vol, support, extra=(), labels=None, contours=(), pad=0.05, head=False):
        i = vol.sag_index(support > 0.5, off)
        img = vol.sagittal(i)
        sup = np.rot90((support > 0.5)[i]) if not head else img > np.percentile(img[img > 0], 20)
        sy, sx = sch.fit_crop(img, asp * vol.zooms[1] / vol.zooms[2], sup, pad=pad)
        cp = lambda a: sch.crop_pad(a, sy, sx)
        cs = [(cp(np.rot90(m[i])), col, lw) for m, col, lw in contours]
        lab = cp(np.rot90(labels[i])) if labels is not None else None
        hi = np.percentile(img[img > 0], 99.5)
        return Thumb(cp(img), vmin=0, vmax=hi, contours=cs, labels=lab, aspect=vol.sag_aspect)

    g, y = '#3DDC97', '#FFD400'
    tmask_w = tm.data > 0.5
    return dict(
        raw=sl(raw, m1.data, head=True),
        n4=sl(n4, m1.data, head=True),
        mask1=sl(n4, m1.data, head=True, contours=[(m1.data, g, 0.6)]),
        acpc=sl(acpc, m2.data, head=True),
        mask2=sl(acpc, m2.data, contours=[(m2.data, g, 0.6)], pad=0.12),
        warped=sl(warped, tmask_w, contours=[(tmask_w, y, 0.45)]),
        tatlas=sl(tss, V(tf['mask']).resample_like(tss).data, labels=tatlas.data),
        natlas=sl(ss, m2.resample_like(ss).data, labels=atlas.data),
    ), f'{files["ss"]}'


def get_thumbs(force_synthetic=False):
    if not force_synthetic:
        files, ses = sch.find_session(ANAT_BIDS, PATTERNS)
        tf = template_files()
        if files and tf:
            try:
                th, src = thumbs_real(files, tf)
                print(f'  [Supp 1] images: {ses}')
                return th, ses
            except Exception as e:                                     # noqa: BLE001
                print(f'  [Supp 1] [WARN] real images failed ({e}) -- synthetic brain used')
        else:
            miss = sch.missing_report(ANAT_BIDS, PATTERNS) if os.path.isdir(ANAT_BIDS) else [ANAT_BIDS]
            print('  [Supp 1] [WARN] no complete session / template -- synthetic brain used. Missing: '
                  + '; '.join(miss + ([] if tf else [f'template of {SPECIES} in {pdat.ATLAS_LIB}'])))
    th, src = thumbs_synthetic()
    return th, 'synthetic: ' + src


# ── figure ───────────────────────────────────────────────────────────────────
def panel_a(cv, th, y0):
    cv.panel('a', 1.0, y0, 'One session through its coordinate spaces')
    reg, app = OPS['reg']['color'], OPS['apply']['color']
    frames = [
        dict(thumb=th['raw'], title='Raw T1w', step=1, op='prepro', space='scanner',
             caption='reoriented, averaged'),
        dict(thumb=th['n4'], title='Bias field', step=1, op='prepro', space='scanner',
             caption='N4, head mask', link=dict(kind='flow')),
        dict(thumb=th['mask1'], title='Brain mask 1', step=2, op='mask', space='scanner',
             caption='coarse, species method', link=dict(kind='flow')),
        dict(thumb=th['acpc'], title='AC-PC', step=2, op='reg', space='native',
             caption='6 dof, MI, cropped',
             link=dict(kind='reg', label='Rigid')),
        dict(thumb=th['mask2'], title='Brain mask 2', step=5, op='mask', space='native',
             caption='fine mask, N4 in mask', link=dict(kind='flow')),
        dict(thumb=th['warped'], title='To template', step=6, op='reg', space='template',
             caption='template WM edges', link=dict(kind='reg', label='SyN')),
        dict(thumb=th['tatlas'], title='Template atlas', step=None, space='template',
             caption=f'{ATLAS} (level 2)'),
        dict(thumb=th['natlas'], title='Native atlas', step=8, op='apply', space='native',
             caption='genericLabel interp.', link=dict(kind='apply', label='T$^{-1}$')),
    ]
    boxes = sch.draw_strip(cv, frames, 3.0, y0 + 5.0, sch.W - 5.0, FW, FH)
    # inverse chain under the strip: from the template band back to native
    ybot = boxes[0][1] + FH + 5.6
    x5 = boxes[4][0] + FW / 2
    x7 = boxes[6][0] + FW / 2
    return ybot


def panel_b(cv, y0):
    cv.panel('b', 1.0, y0, 'Steps of the anatomical module (numbers = Skip_step)')
    y = y0 + 6.0
    colw, gap = 43.2, 2.0
    xs = [3.0 + i * (colw + gap) for i in range(4)]
    blocks = [
        ('Block 1', 'each session (_loop1)', [
            (1, 'Orientation, averaging, N4', ['header orientation, deoblique (3dWarp)',
                                               'runs: rigid motion-corrected mean',
                                               'N4 (or N3) bias field in a head mask'], 'prepro',
             'AFNI, ANTsPy, mri_robust_template'),
            (2, 'Brain mask, pass 1', ['coarse mask; method chosen per species',
                                       '(template registration, U-Net, 3dSkullStrip,',
                                       'BET, histogram, manual...)'], 'mask', 'brain_skullstrip_1'),
            (2, 'AC-PC alignment', ['translation then rigid (6 dof) to the',
                                    'template, Mattes MI; crop (3dAutobox)',
                                    'one matrix applied to every contrast'], 'reg', 'ANTsPy, AFNI'),
            ('itk_1', 'Mask check', ['edited mask re-runs step 2'], 'check', 'ITK-SNAP'),
        ]),
        ('Study template', 'optional, all subjects (studytemplate)', [
            (3, 'Study template', ['centre of mass, rigid, affine, then',
                                   'iterative non-linear (3dQwarp) to the',
                                   'mean of all subjects'], 'reg', 'AFNI anats_to_common (sammba-MRI)'),
            (4, 'Template brain mask', ['skull-stripping of the study template'], 'mask', 'template_skullstrip'),
            (7, 'Study to reference template', ['SyN registration (block stdyT);',
                                                'atlases and masks brought to the',
                                                'study template'], 'reg', 'ANTsPy'),
        ]),
        ('Block 2', 'each session (_loop2)', [
            (5, 'Brain mask, pass 2', ['fine mask from the AC-PC image,',
                                       'registered to the study / reference',
                                       'template; holes filled; adaptive N4'], 'mask', 'brain_skullstrip_2'),
            ('itk_2', 'Mask check', ['edited mask re-runs step 5'], 'check', 'ITK-SNAP'),
            (6, 'Registration to the template', ['translation, then SyNCC (default) to',
                                                 'reference / study template or to the',
                                                 'reference session (longitudinal)'], 'reg',
             'ANTsPy; forward and inverse QC mosaics'),
        ]),
        ('Block 3', 'each session (_loop3)', [
            (8, 'Back to native space', ['composed inverse transforms bring',
                                         'atlases (genericLabel) and brain /',
                                         'tissue masks to the AC-PC image'], 'apply', 'ANTsPy apply_transforms'),
            ('itk_3', 'Native mask check', ['edit re-runs steps 6 and 8'], 'check', 'ITK-SNAP'),
            (9, 'Masks for fMRI', ['grey matter, eroded white matter and',
                                   'ventricles for the functional module'], 'mask', None),
            ('10-15', 'Surfaces', ['FreeSurfer white / pial (small brains',
                                   'rescaled), flat maps, Workbench'], 'analysis', 'FreeSurfer, Connectome Workbench'),
            (16, 'Quality control', ['SNR, CNR, correlation to template'], 'analysis', None),
        ]),
    ]
    ymax = y
    for x, (name, sub, cards) in zip(xs, blocks):
        cv.block_header(x, y, colw, name, sub, color=SPACES['study']['ink'] if name == 'Study template' else sch.TEXT)
        yy = y + 6.4
        for step, title, lines, op, tool in cards:
            h = cv.card(x, yy, colw, step, title, lines, op, tool=tool, dashed=(name == 'Study template'))
            yy += h + 1.3
        ymax = max(ymax, yy)
    return ymax


def panel_c(cv, y0):
    """Registration designs (Illustration/coregistration_strategies.svg)."""
    cv.panel('c', 1.0, y0, 'Registration designs')
    designs = [
        ('Classical', [('Session', 'native'), ('Reference template', 'template')], ['6']),
        ('Study template', [('Session', 'native'), ('Study template', 'study'), ('Reference template', 'template')],
         ['6', '7']),
        ('Longitudinal', [('Session k', 'native'), ('Reference session', 'native'), ('Reference template', 'template')],
         ['6', '6']),
        ('Longitudinal + study template', [('Session k', 'native'), ('Reference session', 'native'),
                                           ('Study template', 'study'), ('Reference template', 'template')],
         ['6', '6', '7']),
    ]
    colw, gap = 43.2, 2.0
    y = y0 + 6.0
    nw, nh, step_y = 24.0, 4.8, 8.0
    ymax = y
    for k, (name, nodes, steps) in enumerate(designs):
        x = 3.0 + k * (colw + gap)
        cv.text(x, y, name, size=sch.FS_TEXT + 0.3, weight='bold')
        nx = x + 1.0
        ny = y + 3.6
        for i, (lab, sp) in enumerate(nodes):
            s = SPACES[sp]
            cv.rect(nx, ny + i * step_y, nw, nh, fc=s['fill'], ec=s['ink'], lw=0.4, r=1.2, z=3,
                    ls=(0, (2, 1.2)) if sp == 'study' else '-')
            cv.text(nx + nw / 2, ny + i * step_y + nh / 2, lab, size=sch.FS_TEXT - 0.3, ha='center', va='center',
                    color=s['ink'], weight='bold')
            if i > 0:
                ya, yb = ny + (i - 1) * step_y + nh + 0.2, ny + i * step_y - 0.2
                cv.arrow(nx + nw * 0.38, ya, nx + nw * 0.38, yb, color=OPS['reg']['color'], lw=0.9, head=3.4)
                cv.badge(nx + nw * 0.38 - 2.8, (ya + yb) / 2, steps[i - 1], OPS['reg']['color'], r=1.3, size=4.4)
        # inverse: atlases and masks of the reference template back to the session
        top, bot = ny + nh / 2, ny + (len(nodes) - 1) * step_y + nh / 2
        xr = cv.elbow(nx + nw + 0.2, bot, top, 3.2, OPS['apply']['color'])
        cv.badge(xr + 2.6, (top + bot) / 2 - 1.6, '8', OPS['apply']['color'], r=1.3, size=4.4)
        cv.text(xr + 1.2, (top + bot) / 2 + 0.6, 'atlases,\nmasks', size=sch.FS_SMALL,
                va='top', color=OPS['apply']['color'])
        ymax = max(ymax, ny + (len(nodes) - 1) * step_y + nh)
    return ymax


def make(out=OUT, force_synthetic=False):
    th, src = get_thumbs(force_synthetic)
    H = 200.0
    cv = sch.Canvas(H)
    ya = panel_a(cv, th, 1.0)
    yb = panel_b(cv, ya + 1.5)
    yc = panel_c(cv, yb + 1.5)
    # legend
    yl = yc + 4.0
    x = cv.ops_legend(3.0, yl)
    yl2 = yl + 3.6
    cv.arrow_legend(3.0, yl2)
    cv.text(sch.W - 2.0, yl2, 'Images: ' + ('synthetic illustration' if src.startswith('synthetic') else 'example session'),
            size=sch.FS_SMALL, ha='right', va='center', color=sch.SUBTLE)
    cv.trim(yl2 + 3.0)
    print(f'  [Supp 1] image source: {src}')
    return cv.save(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--synthetic', action='store_true', help='synthetic images even if the data are found')
    ap.add_argument('--out', default=None, help='output folder (default <PAPER_DIR>/Supp_Fig_01)')
    a = ap.parse_args()
    out = os.path.join(a.out, 'Supp_Fig_01') if a.out else OUT
    make(out, a.synthetic)


if __name__ == '__main__':
    main()
