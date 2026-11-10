"""
Human anat -> template: can the label overlap of EDNiX be improved beyond SyNCC? (BIDS_Park, recon-all aseg as
independent labels). Standalone (no EDNiX import): runs on any EDNiX checkout.

Variants (all from the same translation, fixed = EDNiX human template, moving = EDNiX native skull-stripped T1)
  A   SyNCC, as EDNiX step 6 (reference).
  B   SyNadaptive refinement: compliant masked SyN from A's affine (kept to complete the n = 10 of the first test).
  C   true refinement of A: compliant masked SyN between the template and the subject ALREADY warped by A
      (final transform = C o A). Tests "more flexibility" properly.
  D   multichannel SyN: CC on T1 + CC on GM and WM probability maps (subject posteriors vs template posteriors),
      from A's affine, same schedule as SyNCC. Tissue-guided registration (Avants et al., 2011, ANTs multivariate).
  E   A + posterior relabelling of the cortex / cerebral WM ribbon (CSF connected to the atlas ventricles becomes
      ventricle): a 3-class Gaussian mixture of the subject T1 with
      the template CSF/GM/WM priors warped by A as spatial prior (atlas x intensity likelihood, as in FreeSurfer,
      Fischl et al., 2002, and Atropos prior-based segmentation, Avants et al., 2011). Deep structures, ventricles and
      cerebellum keep the atlas labels; WM voxels far from the atlas cortex stay WM (no WM hypointensity -> cortex).
  DE  E applied after D.

Evaluation (native space, where EDNiX uses the atlases): template seg-4FS brought to native space by each variant,
Dice vs the recon-all aseg of the same scan (12 structures, mean of hemispheres), plus intensity checks.

Outputs (OUT_DIR): QC/<subject>_QC.png and a line in the log AS SOON AS a subject is finished; at the end
variants_subjects.csv, variants_structures.csv, variants_summary.json, variants.png/.pdf.
Registrations are cached (A and B are reused from the SyNadaptive calibration folder when present).

    cd ~/Downloads && nohup python calibrate_registration_variants.py > calib_variants.log 2>&1 &
    python calibrate_registration_variants.py --analyse-only
"""
import os
import sys
import glob
import json
import time
import argparse
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

opj = os.path.join
opi = os.path.isfile

# ----------------------------------------------------------------------------------------------------------- settings
MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
BIDS_DIR = '/scratch2/EDNiX/Human/BIDS_Park'
FS_SUBJECTS_DIR = '/scratch2/EDNiX/FreeSurfer_benchmark/BIDS_Park'      # recon-all: sub-<S>_ses-<T>
LIBRARY = opj(os.path.dirname(MAIN_PATH), 'Atlases_library')
SPECIES, REFERENCE = 'Human', 'EDNiX'
CACHE_DIR = '/scratch2/EDNiX/SyNadaptive_calibration/BIDS_Park'         # A, B, recon-all aseg of the first test
OUT_DIR = '/scratch2/EDNiX/registration_variants/BIDS_Park'
VARIANTS = ['A', 'B', 'C', 'D', 'E', 'DE']
N_JOBS = 5
THREADS_PER_JOB = 8
FS_ALIGN_MIN_CORR = 0.80
PRIOR_WEIGHT = 0.5          # exponent of the spatial prior in the posterior (E, and the tissue maps of D)
RIBBON_MAX_DIST = 3         # voxels: WM farther than this from the atlas cortex is never relabelled as cortex

STRUCTURES = {'cortex': (3, 42), 'cerebral_WM': (2, 41), 'lateral_ventricle': (4, 43), 'thalamus': (10, 49),
              'caudate': (11, 50), 'putamen': (12, 51), 'pallidum': (13, 52), 'hippocampus': (17, 53),
              'amygdala': (18, 54), 'cerebellum_cortex': (8, 47), 'cerebellum_WM': (7, 46), 'brainstem': (16,)}
REFINE = dict(syn_metric='CC', syn_sampling=4, grad_step=0.2, flow_sigma=2, total_sigma=0,
              reg_iterations=(100, 70, 50, 20))
SYNCC_LIKE = dict(syn_metric='CC', syn_sampling=4, grad_step=0.15, flow_sigma=3, total_sigma=0,
                  reg_iterations=(2100, 1200, 1200, 20))
COLORS = {'A': '#1F77B4', 'B': '#7F7F7F', 'C': '#FF7F0E', 'D': '#2CA02C', 'E': '#D62728', 'DE': '#9467BD'}


# ------------------------------------------------------------------------------------------------------------ inputs
# FreeSurfer codes of the template seg-4FS grouped into tissues, to derive priors when the library has none
PRIOR_CODES = {'CSF': (4, 5, 14, 15, 24, 43, 44, 72),
               'GM': (3, 42, 8, 47, 10, 11, 12, 13, 17, 18, 26, 28, 49, 50, 51, 52, 53, 54, 58, 60),
               'WM': (2, 41, 7, 46, 16, 77, 78, 85, 251, 252, 253, 254, 255)}


def template_files():
    """Template T1, seg-4FS and CSF/GM/WM priors. Priors: <volumes>/priors/ if present, otherwise derived from the
    seg-4FS (one-hot tissue groups, unlabelled brain voxels = CSF, Gaussian smoothing sigma 1 voxel, normalised) and
    cached in OUT_DIR/template_priors/."""
    hits = sorted(glob.glob(opj(LIBRARY, 'atlas', '**', SPECIES, REFERENCE, 'volumes',
                                f'{SPECIES}_space-acpc_desc-SS_T1w.nii.gz'), recursive=True),
                  key=lambda p: (p.count(os.sep), p))
    hits = [h for h in hits if f'{os.sep}freesurfer{os.sep}' not in h]
    if not hits:
        sys.exit(f'No {SPECIES}/{REFERENCE} template under {LIBRARY}/atlas')
    vol = os.path.dirname(hits[0])
    f = dict(t1=hits[0], seg=opj(vol, 'labels', f'{SPECIES}_seg-4FS_dseg.nii.gz'),
             **{k: opj(vol, 'priors', f'{SPECIES}_label-{k}_probseg.nii.gz') for k in ('CSF', 'GM', 'WM')})
    if not opi(f['seg']):
        sys.exit(f'Missing template segmentation: {f["seg"]}')
    if not all(opi(f[k]) for k in ('CSF', 'GM', 'WM')):
        import nibabel as nib
        from scipy import ndimage as ndi
        out = opj(OUT_DIR, 'template_priors')
        os.makedirs(out, exist_ok=True)
        seg_img = nib.load(f['seg'])
        seg = np.asanyarray(seg_img.dataobj).astype(int)
        brain = np.asanyarray(nib.load(f['t1']).dataobj) > 0
        hard = {k: np.isin(seg, codes).astype(np.float32) for k, codes in PRIOR_CODES.items()}
        hard['CSF'][brain & (seg == 0)] = 1
        soft = {k: ndi.gaussian_filter(v, 1.0) for k, v in hard.items()}
        tot = sum(soft.values()) + 1e-6
        for k in ('CSF', 'GM', 'WM'):
            f[k] = opj(out, f'{SPECIES}_label-{k}_probseg_from-seg4FS.nii.gz')
            nib.save(nib.Nifti1Image((soft[k] / tot * brain).astype(np.float32), seg_img.affine), f[k])
        print(f'INFO: no tissue priors in the library: derived from the seg-4FS -> {out}', flush=True)
    return f


def subjects():
    out = []
    for fs_dir in sorted(glob.glob(opj(FS_SUBJECTS_DIR, 'sub-*_ses-*'))):
        fsid = os.path.basename(fs_dir)
        if not opi(opj(fs_dir, 'scripts', 'recon-all.done')) or opi(opj(fs_dir, 'scripts', 'recon-all.error')):
            continue
        sub, ses = fsid[4:].split('_ses-')
        t1 = sorted(t for t in glob.glob(opj(BIDS_DIR, f'sub-{sub}', f'ses-{ses}', 'anat', 'native', 'volumes',
                                             '*_space-acpc_desc-SS_*T1w.nii.gz')) if 'desc-SS-step' not in t)
        if t1:
            out.append((fsid, t1[0], fs_dir))
        else:
            print(f'WARNING: {fsid}: no EDNiX space-acpc_desc-SS T1w, skipped', flush=True)
    return out


# ----------------------------------------------------------------------------------------------------------- metrics
def dice(a, b):
    s = a.sum() + b.sum()
    return float(2 * np.logical_and(a, b).sum() / s) if s else np.nan


def label_dice(ref, test):
    out = {}
    for name, codes in STRUCTURES.items():
        vals = [dice(ref == c, test == c) for c in codes if (ref == c).any() and (test == c).any()]
        out[name] = float(np.mean(vals)) if vals else np.nan
    return out


def kmeans_labels(values, k=3, iters=100):
    c = np.percentile(values, np.linspace(100 / (2 * k), 100 - 100 / (2 * k), k))
    for _ in range(iters):
        lab = np.abs(values[:, None] - c[None, :]).argmin(1)
        new = np.array([values[lab == i].mean() if (lab == i).any() else c[i] for i in range(k)])
        if np.allclose(new, c):
            break
        c = np.sort(new)
    return np.abs(values[:, None] - c[None, :]).argmin(1) + 1


def intensity_checks(fixed, moving, fwd, warp_for_folding):
    """Label-free checks in template space (as SyNadaptive): mask Dice, k-means tissue Dice, folding fraction."""
    import ants
    fm = fixed.numpy() > 0
    warped = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=fwd, interpolator='linear')
    wm_ = ants.apply_transforms(fixed=fixed, moving=ants.get_mask(moving), transformlist=fwd,
                                interpolator='nearestNeighbor').numpy() > 0.5
    lf = kmeans_labels(fixed.numpy()[fm].astype(float))
    lw = kmeans_labels(warped.numpy()[fm].astype(float))
    jac = ants.create_jacobian_determinant_image(fixed, warp_for_folding, do_log=False).numpy()
    return dict(mask_dice=dice(fm, wm_), tissue_dice=float(np.nanmean([dice(lf == c, lw == c) for c in (1, 2, 3)])),
                folding=float((jac[fm] <= 0).mean()))


# ------------------------------------------------------------------------------------------- tissue (prior x intensity)
def gmm_posteriors(img, priors, mask, w=PRIOR_WEIGHT, iters=15):
    """3-class Gaussian mixture of the intensities with spatial priors (CSF, GM, WM), posteriors within mask."""
    x = img[mask].astype(float)
    P = np.stack([np.clip(p[mask], 0, 1) for p in priors], 1) + 1e-6
    P /= P.sum(1, keepdims=True)
    Pw = P ** w
    post = P.copy()
    for _ in range(iters):
        nk = post.sum(0) + 1e-9
        mu = (post * x[:, None]).sum(0) / nk
        sd = np.sqrt((post * (x[:, None] - mu) ** 2).sum(0) / nk) + 1e-6
        lik = np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / sd
        post = Pw * lik
        post /= post.sum(1, keepdims=True) + 1e-300
    out = [np.zeros(img.shape, np.float32) for _ in range(3)]
    for k in range(3):
        out[k][mask] = post[:, k]
    return out


def relabel_ribbon(seg, posts, mask):
    """Cortex / cerebral WM / sulcal CSF of the atlas ribbon re-decided by the posteriors (hemisphere from the atlas)."""
    from scipy import ndimage as ndi
    L, R = np.isin(seg, (2, 3)), np.isin(seg, (41, 42))
    ribbon = L | R
    zone = mask & (ribbon | (ndi.binary_dilation(ribbon, iterations=1) & (seg == 0)))
    left = ndi.distance_transform_edt(~L) <= ndi.distance_transform_edt(~R)
    near_cortex = ndi.distance_transform_edt(~np.isin(seg, (3, 42))) <= RIBBON_MAX_DIST
    cls = np.argmax(np.stack(posts, 0), 0)                      # 0 CSF, 1 GM, 2 WM
    out = seg.copy()
    z = zone & (cls == 0)
    out[z] = 0
    # CSF connected to the atlas lateral ventricles = ventricle (enlarged ventricles of atrophied brains)
    vent = np.isin(seg, (4, 43))
    comp, _ = ndi.label(z | vent)
    keep = np.isin(comp, np.unique(comp[vent])) & z & (comp > 0)
    out[keep] = np.where(left[keep], 4, 43)
    z = zone & (cls == 1) & (near_cortex | np.isin(seg, (3, 42)))
    out[z] = np.where(left[z], 3, 42)
    z = zone & ((cls == 2) | ((cls == 1) & ~(near_cortex | np.isin(seg, (3, 42)))))
    out[z] = np.where(left[z], 2, 41)
    return out


# ---------------------------------------------------------------------------------------------------------- subject
def run_subject(args):
    fsid, t1_path, fs_dir, tpl = args
    for k in ('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[k] = str(THREADS_PER_JOB)
    import ants
    import nibabel as nib

    d, c = opj(OUT_DIR, fsid), opj(CACHE_DIR, fsid)
    os.makedirs(d, exist_ok=True)
    os.makedirs(c, exist_ok=True)
    res_json = opj(d, 'result.json')
    row = json.load(open(res_json)) if opi(res_json) else {}
    if row and all(f'{v}_label_mean' in row for v in VARIANTS):
        return row
    row.update(fsid=fsid, ednix_T1=t1_path)
    timing = row.setdefault('minutes', {})

    moving, fixed = ants.image_read(t1_path), ants.image_read(tpl['t1'])
    tpl_seg = ants.image_read(tpl['seg'])
    tpl_pri = [ants.image_read(tpl[k]) for k in ('CSF', 'GM', 'WM')]
    mmask = moving.numpy() > 0

    # recon-all aseg -> EDNiX native (rigid, same scan), cached
    aseg_f = opj(c, 'fs_aseg_in_ednix_native.nii.gz')
    if not opi(aseg_f):
        def mgz(p, out):
            im = nib.load(p)
            nib.save(nib.Nifti1Image(np.asanyarray(im.dataobj), im.affine), out)
            return ants.image_read(out)
        fsb = mgz(opj(fs_dir, 'mri', 'brainmask.mgz'), opj(c, 'fs_brainmask.nii.gz'))
        fsa = mgz(opj(fs_dir, 'mri', 'aseg.mgz'), opj(c, 'fs_aseg.nii.gz'))
        rig = ants.registration(fixed=moving, moving=fsb, type_of_transform='Rigid', aff_metric='mattes',
                                outprefix=opj(c, 'fs_to_ednix_'))
        ants.image_write(ants.apply_transforms(moving, fsa, rig['fwdtransforms'], interpolator='nearestNeighbor'),
                         aseg_f)
    if 'fs_align_corr' not in row:
        fsb_n = ants.apply_transforms(moving, ants.image_read(opj(c, 'fs_brainmask.nii.gz')),
                                      [opj(c, 'fs_to_ednix_0GenericAffine.mat')], interpolator='linear').numpy()
        m = mmask & (fsb_n > 0)
        row['fs_align_corr'] = float(np.corrcoef(moving.numpy()[m], fsb_n[m])[0, 1])
    ref = ants.image_read(aseg_f).numpy().astype(int)

    def files(p):
        return p + '0GenericAffine.mat', p + '1Warp.nii.gz', p + '1InverseWarp.nii.gz'

    # A (cached from the first test)
    t0 = time.time()
    tr = ants.registration(fixed=fixed, moving=moving, type_of_transform='Translation', aff_metric='mattes',
                           outprefix=opj(c, 'shift_'))
    aA, wA, iA = files(opj(c, 'A_'))
    if not opi(wA):
        ants.registration(fixed=fixed, moving=moving, type_of_transform='SyNCC',
                          initial_transform=tr['fwdtransforms'], outprefix=opj(c, 'A_'))
        timing['A'] = (time.time() - t0) / 60
    chains = {'A': dict(fwd=[wA, aA], inv=[aA, iA], inv_flags=[True, False], warp=wA)}

    # B (cached from the first test)
    if 'B' in VARIANTS:
        aB, wB, iB = files(opj(c, 'B_'))
        if not opi(wB):
            t0 = time.time()
            ants.registration(fixed=fixed, moving=moving, type_of_transform='SyNOnly', initial_transform=aA,
                              mask=ants.get_mask(fixed), moving_mask=ants.get_mask(moving), mask_all_stages=True,
                              outprefix=opj(c, 'B_'), **REFINE)
            timing['B'] = (time.time() - t0) / 60
        chains['B'] = dict(fwd=[wB, aB], inv=[aB, iB], inv_flags=[True, False], warp=wB)

    # C: refinement of the subject already warped by A
    if 'C' in VARIANTS:
        _, wC, iC = files(opj(d, 'C_'))
        if not opi(wC):
            t0 = time.time()
            movA = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=chains['A']['fwd'],
                                         interpolator='linear')
            ants.registration(fixed=fixed, moving=movA, type_of_transform='SyNOnly',
                              mask=ants.get_mask(fixed), moving_mask=ants.get_mask(movA), mask_all_stages=True,
                              outprefix=opj(d, 'C_'), **REFINE)
            timing['C'] = (time.time() - t0) / 60
        chains['C'] = dict(fwd=[wC, wA, aA], inv=[aA, iA, iC], inv_flags=[True, False, False], warp=wC)

    # tissue posteriors of the subject (priors warped by A) and of the template (its own priors)
    need_tissue = any(v in VARIANTS for v in ('D', 'E', 'DE'))
    if need_tissue:
        pri_n = [ants.apply_transforms(moving, p, chains['A']['inv'], whichtoinvert=chains['A']['inv_flags'],
                                       interpolator='linear').numpy() for p in tpl_pri]
        post_sub = gmm_posteriors(moving.numpy(), pri_n, mmask)
        post_tpl = gmm_posteriors(fixed.numpy(), [p.numpy() for p in tpl_pri], fixed.numpy() > 0)

    # D: multichannel SyN (T1 + GM + WM posteriors)
    if 'D' in VARIANTS or 'DE' in VARIANTS:
        aD, wD, iD = files(opj(d, 'D_'))
        if not opi(wD):
            t0 = time.time()
            sub_img = [moving.new_image_like(p) for p in post_sub]
            tpl_img = [fixed.new_image_like(p) for p in post_tpl]
            ants.registration(fixed=fixed, moving=moving, type_of_transform='SyNOnly', initial_transform=aA,
                              multivariate_extras=[['CC', tpl_img[1], sub_img[1], 1, 2],
                                                   ['CC', tpl_img[2], sub_img[2], 1, 2]],
                              outprefix=opj(d, 'D_'), **SYNCC_LIKE)
            timing['D'] = (time.time() - t0) / 60
        chains['D'] = dict(fwd=[wD, aD], inv=[aD, iD], inv_flags=[True, False], warp=wD)

    # evaluation of every variant
    segs = {}
    for v in VARIANTS:
        base = {'E': 'A', 'DE': 'D'}.get(v, v)
        ch = chains[base]
        seg_n = ants.apply_transforms(moving, tpl_seg, ch['inv'], whichtoinvert=ch['inv_flags'],
                                      interpolator='nearestNeighbor').numpy().astype(int)
        if v in ('E', 'DE'):
            if v == 'DE':   # posteriors with the priors warped by D
                pri_d = [ants.apply_transforms(moving, p, ch['inv'], whichtoinvert=ch['inv_flags'],
                                               interpolator='linear').numpy() for p in tpl_pri]
                seg_n = relabel_ribbon(seg_n, gmm_posteriors(moving.numpy(), pri_d, mmask), mmask)
            else:
                seg_n = relabel_ribbon(seg_n, post_sub, mmask)
        segs[v] = seg_n
        nib.save(nib.Nifti1Image(seg_n.astype(np.int16), nib.load(t1_path).affine),
                 opj(d, f'{v}_template_seg_in_native.nii.gz'))
        ld = label_dice(ref, seg_n)
        for k, val in ld.items():
            row[f'{v}_label_{k}'] = val
        row[f'{v}_label_mean'] = float(np.nanmean(list(ld.values())))
        if v not in ('E', 'DE') and f'{v}_tissue_dice' not in row:
            row.update({f'{v}_{k}': val for k, val in intensity_checks(fixed, moving, ch['fwd'], ch['warp']).items()})
    with open(res_json, 'w') as f:
        json.dump(row, f, indent=1)
    return row


# --------------------------------------------------------------------------------------------------- per-subject QC
QC_CLASSES = {'cortex': ((3, 42), '#D62728'), 'WM': ((2, 41), '#1F77B4'), 'ventricles': ((4, 43, 5, 44), '#FFD700'),
              'subcortical': ((10, 49, 11, 50, 12, 51, 13, 52), '#2CA02C'),
              'hippo/amyg': ((17, 53, 18, 54), '#9467BD'), 'cerebellum': ((7, 46, 8, 47), '#FF7F0E')}


def subject_qc(fsid):
    import nibabel as nib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = opj(OUT_DIR, fsid)
    r = json.load(open(opj(d, 'result.json')))

    def classes(a):
        out = np.zeros(a.shape, np.int16)
        for i, (codes, _) in enumerate(QC_CLASSES.values(), 1):
            out[np.isin(a, codes)] = i
        return out

    def sl(v, c):
        return [np.rot90(v[c[0], :, :]), np.rot90(v[:, c[1], :]), np.rot90(v[:, :, c[2]])]

    t1 = np.asanyarray(nib.load(r['ednix_T1']).dataobj).astype(float)
    ref = classes(np.asanyarray(nib.load(opj(CACHE_DIR, fsid, 'fs_aseg_in_ednix_native.nii.gz')).dataobj))
    nz = np.argwhere(t1 > 0)
    c = (nz.min(0) + nz.max(0)) // 2
    vmax = np.percentile(t1[t1 > 0], 99)
    cols = [col for _, col in QC_CLASSES.values()]
    rows = [('recon-all aseg', ref, None)] + [
        (v, classes(np.asanyarray(nib.load(opj(d, f'{v}_template_seg_in_native.nii.gz')).dataobj)), ref)
        for v in VARIANTS]
    fig, ax = plt.subplots(len(rows), 3, figsize=(10, 3.3 * len(rows)), facecolor='k')
    for i, (name, seg, filled) in enumerate(rows):
        for j, (img, s) in enumerate(zip(sl(t1, c), sl(seg, c))):
            a = ax[i, j]
            a.imshow(img, cmap='gray', vmin=0, vmax=vmax)
            if filled is not None:
                f = sl(filled, c)[j].astype(float)
                f[f == 0] = np.nan
                a.imshow(f, cmap=matplotlib.colors.ListedColormap(cols), vmin=0.5, vmax=len(cols) + 0.5,
                         alpha=0.25, interpolation='nearest')
            for k, col in enumerate(cols, 1):
                if (s == k).any():
                    a.contour(s == k, levels=[0.5], colors=[col], linewidths=0.6)
            a.set_axis_off()
        txt = (f"fs_align r={r['fs_align_corr']:.3f}" if filled is None else
               f"label Dice {r[f'{name}_label_mean']:.3f} | cortex {r[f'{name}_label_cortex']:.3f} | "
               f"WM {r[f'{name}_label_cerebral_WM']:.3f} | ventr. {r[f'{name}_label_lateral_ventricle']:.3f}")
        ax[i, 0].set_title(f'{name}: {txt}', color='w', loc='left', fontsize=8)
    fig.legend(handles=[matplotlib.patches.Patch(color=col, label=n) for n, (_, col) in QC_CLASSES.items()],
               loc='lower center', ncol=6, fontsize=7, facecolor='k', labelcolor='w')
    fig.suptitle(f'{fsid}: contours = labels of each variant, fill = recon-all aseg', color='w', fontsize=10)
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    os.makedirs(opj(OUT_DIR, 'QC'), exist_ok=True)
    out = opj(OUT_DIR, 'QC', f'{fsid}_QC.png')
    fig.savefig(out, dpi=100, facecolor='k')
    plt.close(fig)
    return out


def subject_line(r):
    return ' | '.join(f"{v} {r[f'{v}_label_mean']:.3f}" for v in VARIANTS) + \
        f" | cortex A {r['A_label_cortex']:.3f}" + ''.join(
            f" {v} {r[f'{v}_label_cortex']:.3f}" for v in VARIANTS if v != 'A')


# --------------------------------------------------------------------------------------------------------- analysis
def boot_ci(x, n=10000, seed=0):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return [np.nan, np.nan]
    m = np.random.default_rng(seed).choice(x, (n, len(x))).mean(1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def analyse(rows):
    import pandas as pd
    from scipy import stats
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    df = pd.DataFrame([{k: v for k, v in r.items() if k != 'minutes'} for r in rows]).sort_values('fsid')
    df.to_csv(opj(OUT_DIR, 'variants_subjects.csv'), index=False)
    df = df[df['fs_align_corr'] >= FS_ALIGN_MIN_CORR].reset_index(drop=True)
    summ, srows = dict(n_subjects=len(df), variants={}), []
    for v in VARIANTS:
        dl = df[f'{v}_label_mean'] - df['A_label_mean']
        p = stats.wilcoxon(df[f'{v}_label_mean'], df['A_label_mean']).pvalue if v != 'A' and len(df) >= 5 \
            and (dl != 0).any() else np.nan
        summ['variants'][v] = dict(label_mean=float(df[f'{v}_label_mean'].mean()), delta_vs_A=float(dl.mean()),
                                   delta_CI95=boot_ci(dl), n_better_than_A=int((dl > 0).sum()), wilcoxon_p=p)
        for s in STRUCTURES:
            a, b = df[f'A_label_{s}'], df[f'{v}_label_{s}']
            ok = a.notna() & b.notna()
            srows.append(dict(variant=v, structure=s, mean=b[ok].mean(), delta_vs_A=(b - a)[ok].mean(),
                              delta_CI95=boot_ci((b - a)[ok]), n_better=int(((b - a)[ok] > 0).sum())))
    sdf = pd.DataFrame(srows)
    sdf.to_csv(opj(OUT_DIR, 'variants_structures.csv'), index=False)
    best = max(summ['variants'], key=lambda v: summ['variants'][v]['label_mean'])
    summ['best'] = best
    with open(opj(OUT_DIR, 'variants_summary.json'), 'w') as f:
        json.dump(summ, f, indent=1, default=float)

    fig, ax = plt.subplots(1, 2, figsize=(14, 4.5), gridspec_kw=dict(width_ratios=[1, 2.2]))
    for i, v in enumerate(VARIANTS):
        y = df[f'{v}_label_mean']
        ax[0].scatter(np.full(len(y), i), y, color=COLORS.get(v, 'k'), s=14)
        ax[0].plot([i - 0.25, i + 0.25], [y.mean()] * 2, color='k')
    for k in range(len(df)):
        ax[0].plot(range(len(VARIANTS)), [df.loc[k, f'{v}_label_mean'] for v in VARIANTS], color='0.8', lw=0.5,
                   zorder=0)
    ax[0].set_xticks(range(len(VARIANTS)), VARIANTS)
    ax[0].set_ylabel('mean label Dice vs recon-all aseg')
    ax[0].set_title(f'a  per subject (n = {len(df)})', loc='left', fontsize=9)
    w = 0.8 / max(1, len(VARIANTS) - 1)
    for i, v in enumerate([v for v in VARIANTS if v != 'A']):
        s = sdf[sdf['variant'] == v]
        ax[1].bar(np.arange(len(STRUCTURES)) + (i - (len(VARIANTS) - 2) / 2) * w, s['delta_vs_A'], w,
                  color=COLORS.get(v, 'k'), label=v)
    ax[1].axhline(0, color='k', lw=0.6)
    ax[1].set_xticks(range(len(STRUCTURES)), list(STRUCTURES), rotation=45, ha='right', fontsize=8)
    ax[1].set_ylabel(r'$\Delta$ Dice vs A (SyNCC)')
    ax[1].legend(fontsize=7, ncol=len(VARIANTS))
    ax[1].set_title('b  per structure', loc='left', fontsize=9)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(opj(OUT_DIR, f'variants.{ext}'), dpi=200)
    plt.close(fig)
    print(json.dumps(summ, indent=1, default=float), flush=True)


# ------------------------------------------------------------------------------------------------------------- main
def main(analyse_only=False):
    os.makedirs(OUT_DIR, exist_ok=True)
    tpl = template_files()
    subs = subjects()
    print(f'template: {tpl["t1"]}\n{len(subs)} subject(s); variants {VARIANTS}', flush=True)
    rows = []
    if analyse_only:
        rows = [json.load(open(f)) for f in sorted(glob.glob(opj(OUT_DIR, '*', 'result.json')))]
    else:
        with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
            futs = {ex.submit(run_subject, (*s, tpl)): s[0] for s in subs}
            for f in as_completed(futs):
                try:
                    r = f.result()
                    rows.append(r)
                    print(f'{futs[f]}: done | label Dice {subject_line(r)}', flush=True)
                    try:
                        print(f'   QC: {subject_qc(futs[f])}', flush=True)
                    except Exception as e:
                        print(f'   QC figure failed: {e}', flush=True)
                except Exception:
                    print(f'{futs[f]}: FAILED\n{traceback.format_exc()}', flush=True)
    if len(rows) < 3:
        sys.exit(f'{len(rows)} subject(s): not enough for the group analysis.')
    analyse(rows)
    print(f'Outputs: {OUT_DIR}', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--analyse-only', action='store_true')
    main(ap.parse_args().analyse_only)
