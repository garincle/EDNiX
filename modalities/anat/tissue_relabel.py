"""
Tissue relabelling of the registered seg-4FS (cortex / cerebral WM / sulcal CSF).

The template seg-4FS warped to the subject places every label, but its GM/WM boundary is the template's boundary
smoothed by the deformation. Here, inside the cortical ribbon only, each voxel is re-assigned by a 3-class
(CSF / GM / WM) Gaussian mixture fitted on the subject's own intensities (T1w and/or T2w, full covariance), with the
warped labels as spatial prior (unified segmentation: Ashburner & Friston 2005; Atropos with priors: Avants et al.
2011, Tustison et al. 2014). Hemisphere and subcortical labels are kept from the registration.

Calibration (Human, BIDS_Park, n=5, recon-all aseg as reference): label Dice +0.025 vs SyNCC in 5/5 subjects, cortex
~0.58 -> ~0.76, cerebral WM ~0.75 -> ~0.85. Not validated against manual labels nor in non-human species yet ->
off by default ('tissue_relabel': True in the 'coreg' block of list_transfo turns it on).

Distances are expressed in units of the subject's cortical thickness, estimated from the warped cortex label, so the
same rule holds from mouse to human. Relabelling is skipped (registration labels kept) when GM and WM are not
separable in the data (Fisher distance below MIN_GM_WM_SEPARATION).
"""
import json
import os
import shutil

import numpy as np

opj, opi, opd, opb = os.path.join, os.path.isfile, os.path.dirname, os.path.basename

# FreeSurfer codes of seg-4FS
CORTEX, WM, VENT = (3, 42), (2, 41), (4, 43)
LEFT, RIGHT = (2, 3), (41, 42)
PRIOR_CODES = {'CSF': (4, 5, 14, 15, 24, 43, 44, 72),
               'GM': (3, 42, 8, 47, 10, 11, 12, 13, 17, 18, 26, 28, 49, 50, 51, 52, 53, 54, 58, 60),
               'WM': (2, 41, 7, 46, 16, 77, 78, 85, 251, 252, 253, 254, 255)}

DEFAULTS = dict(
    # thickness units; the values reproduce the human calibration (1 mm voxels, estimated thickness ~3 mm:
    # prior sigma 1 voxel, ribbon 3 voxels, ventricle growth 4 voxels)
    prior_weight=0.5,          # exponent of the spatial prior in the posterior
    prior_sigma=0.33,          # Gaussian smoothing of the warped labels -> priors
    ribbon_max_dist=1.0,       # WM farther than this from the warped cortex never becomes cortex
    vent_grow=1.33,            # max geodesic growth of the warped ventricles into CSF
    min_gm_wm_separation=1.0,  # Fisher distance between GM and WM classes below which nothing is changed
    iterations=20,
)


# ------------------------------------------------------------------------------------------------------------ core
def cortical_thickness_mm(seg, zooms):
    """Robust cortical thickness of the warped cortex label: 2 x the 90th percentile of the inner distance (mm)."""
    from scipy import ndimage as ndi
    ctx = np.isin(seg, CORTEX)
    if ctx.sum() < 100:
        return np.nan
    d = ndi.distance_transform_edt(ctx, sampling=zooms)
    return float(2 * np.percentile(d[ctx], 90))


def tissue_priors(seg, mask, sigma_vox):
    """CSF/GM/WM priors from the warped labels (unlabelled brain voxels = CSF), smoothed and normalised."""
    from scipy import ndimage as ndi
    hard = {k: np.isin(seg, c).astype(np.float32) for k, c in PRIOR_CODES.items()}
    hard['CSF'][mask & (seg == 0)] = 1
    soft = [ndi.gaussian_filter(hard[k], sigma_vox) for k in ('CSF', 'GM', 'WM')]
    tot = sum(soft) + 1e-6
    return [s / tot * mask for s in soft]


def gmm_posteriors(channels, priors, mask, w, iters):
    """Multichannel 3-class GMM (full covariance) with spatial priors. Returns posteriors and class stats."""
    X = np.stack([c[mask].astype(float) for c in channels], 1)
    X = (X - np.median(X, 0)) / (np.percentile(X, 99, 0) - np.percentile(X, 1, 0) + 1e-9)
    P = np.stack([np.clip(p[mask], 0, 1) for p in priors], 1) + 1e-6
    P /= P.sum(1, keepdims=True)
    logPw = w * np.log(P)
    post = P.copy()
    d = X.shape[1]
    mu = cov = None
    for _ in range(iters):
        nk = post.sum(0) + 1e-9
        mu = (post.T @ X) / nk[:, None]
        cov = np.empty((3, d, d))
        loglik = np.empty((len(X), 3))
        for k in range(3):
            r = X - mu[k]
            cov[k] = (post[:, k, None] * r).T @ r / nk[k] + 1e-6 * np.eye(d)
            inv = np.linalg.inv(cov[k])
            loglik[:, k] = -0.5 * np.einsum('ij,jk,ik->i', r, inv, r) - 0.5 * np.linalg.slogdet(cov[k])[1]
        lp = logPw + loglik
        lp -= lp.max(1, keepdims=True)
        post = np.exp(lp)
        post /= post.sum(1, keepdims=True)
    out = [np.zeros(mask.shape, np.float32) for _ in range(3)]
    for k in range(3):
        out[k][mask] = post[:, k]
    # Fisher distance GM vs WM along the line joining the two means
    delta = mu[2] - mu[1]
    pooled = 0.5 * (cov[1] + cov[2])
    sep = float(np.sqrt(delta @ np.linalg.solve(pooled, delta)))
    return out, dict(separation_gm_wm=sep, means=mu.tolist())


def relabel(seg, posts, mask, near_cortex_vox, vent_grow_vox):
    """Cortex / cerebral WM / sulcal CSF of the ribbon re-decided by the posteriors (hemisphere from the labels).

    Ribbon CSF becomes background (sulcal CSF is 0 in aseg), except CSF reached from the warped lateral ventricles
    by a geodesic growth of at most vent_grow_vox voxels away from the cortex (enlarged ventricles).
    """
    from scipy import ndimage as ndi
    L, R = np.isin(seg, LEFT), np.isin(seg, RIGHT)
    ribbon = L | R
    zone = mask & (ribbon | (ndi.binary_dilation(ribbon, iterations=1) & (seg == 0)))
    left = ndi.distance_transform_edt(~L) <= ndi.distance_transform_edt(~R)
    is_ctx = np.isin(seg, CORTEX)
    near_cortex = ndi.distance_transform_edt(~is_ctx) <= near_cortex_vox
    cls = np.argmax(np.stack(posts, 0), 0)                      # 0 CSF, 1 GM, 2 WM
    out = seg.copy()
    csf = zone & (cls == 0)
    out[csf] = 0
    vent = np.isin(seg, VENT)
    allowed = (csf & ~near_cortex) | vent
    grown = vent.copy()
    for _ in range(max(int(round(vent_grow_vox)), 0)):
        grown = ndi.binary_dilation(grown) & allowed
    keep = grown & csf
    out[keep] = np.where(left[keep], VENT[0], VENT[1])
    z = zone & (cls == 1) & (near_cortex | is_ctx)
    out[z] = np.where(left[z], CORTEX[0], CORTEX[1])
    z = zone & ((cls == 2) | ((cls == 1) & ~(near_cortex | is_ctx)))
    out[z] = np.where(left[z], WM[0], WM[1])
    return out


def wm_connectivity(seg):
    """Fraction of cerebral WM in its largest 26-connected component (thin WM blades broken -> lower)."""
    from scipy import ndimage as ndi
    wm = np.isin(seg, WM)
    if not wm.any():
        return np.nan
    lab, n = ndi.label(wm, structure=np.ones((3, 3, 3)))
    return float(np.bincount(lab.ravel())[1:].max() / wm.sum()) if n else np.nan


def boundary_gradient(img, seg, mask):
    """Mean intensity-gradient magnitude on the GM/WM boundary voxels (higher = boundary on the image edge)."""
    from scipy import ndimage as ndi
    g = np.sqrt(sum(ndi.sobel(ndi.gaussian_filter(img.astype(float), 0.5), a) ** 2 for a in range(3)))
    g /= np.percentile(g[mask], 99) + 1e-9
    wm = np.isin(seg, WM)
    border = wm & ndi.binary_dilation(np.isin(seg, CORTEX))
    return float(g[border].mean()) if border.any() else np.nan


def compute(seg, channels, mask, zooms, params=None):
    """Relabelled seg and a report (dict). channels: list of arrays on the seg grid (first = reference contrast)."""
    p = dict(DEFAULTS, **(params or {}))
    vox = float(np.mean(zooms))
    th = cortical_thickness_mm(seg, zooms)
    report = dict(params=p, voxel_mm=vox, cortical_thickness_mm=th, n_channels=len(channels), applied=False)
    if not np.isfinite(th) or not np.isin(seg, WM).any():
        report['reason'] = 'no cortex/WM label in the seg'
        return seg, report
    sigma = max(p['prior_sigma'] * th / vox, 0.5)
    priors = tissue_priors(seg, mask, sigma)
    posts, stats = gmm_posteriors(channels, priors, mask, p['prior_weight'], p['iterations'])
    report.update(stats)
    report.update(prior_sigma_vox=sigma, ribbon_max_dist_vox=p['ribbon_max_dist'] * th / vox,
                  vent_grow_vox=p['vent_grow'] * th / vox)
    if stats['separation_gm_wm'] < p['min_gm_wm_separation']:
        report['reason'] = (f"GM/WM not separable in the data (Fisher distance {stats['separation_gm_wm']:.2f} < "
                            f"{p['min_gm_wm_separation']}): registration labels kept")
        return seg, report
    out = relabel(seg, posts, mask, report['ribbon_max_dist_vox'], report['vent_grow_vox'])
    report['applied'] = True
    for name, codes in (('cortex', CORTEX), ('cerebral_WM', WM), ('lateral_ventricle', VENT)):
        a, b = np.isin(seg, codes).sum(), np.isin(out, codes).sum()
        report[f'volume_change_{name}'] = float(b / a - 1) if a else np.nan
    report['changed_fraction_ribbon'] = float((out != seg)[np.isin(seg, CORTEX + WM)].mean())
    report['wm_connectivity_before'] = wm_connectivity(seg)
    report['wm_connectivity_after'] = wm_connectivity(out)
    report['boundary_gradient_before'] = boundary_gradient(channels[0], seg, mask)
    report['boundary_gradient_after'] = boundary_gradient(channels[0], out, mask)
    return out, report


# ------------------------------------------------------------------------------------------------------------- QC
QC_CLASSES = {'cortex': (CORTEX, '#D62728'), 'WM': (WM, '#1F77B4'), 'ventricles': (VENT + (5, 44), '#FFD700'),
              'subcortical': ((10, 49, 11, 50, 12, 51, 13, 52), '#2CA02C'),
              'hippo/amyg': ((17, 53, 18, 54), '#9467BD'), 'cerebellum': ((7, 46, 8, 47), '#FF7F0E')}


def qc_figure(title, img, rows, out_png, zooms=(1.0, 1.0, 1.0)):
    """rows: [(name, seg or None)]; contours on the image, sagittal through a hemisphere, cropped to the brain.
    zooms: voxel sizes, so anisotropic voxels are displayed with their true proportions."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    nz = np.argwhere(img > 0)
    lo, hi = nz.min(0), nz.max(0) + 1
    c = (lo + hi) // 2
    c[0] = lo[0] + int(0.3 * (hi[0] - lo[0]))

    def sl(v):
        return [np.rot90(v[c[0], lo[1]:hi[1], lo[2]:hi[2]]), np.rot90(v[lo[0]:hi[0], c[1], lo[2]:hi[2]]),
                np.rot90(v[lo[0]:hi[0], lo[1]:hi[1], c[2]])]

    def classes(a):
        o = np.zeros(a.shape, np.int16)
        for i, (codes, _) in enumerate(QC_CLASSES.values(), 1):
            o[np.isin(a, codes)] = i
        return o

    cols = [col for _, col in QC_CLASSES.values()]
    vmax = np.percentile(img[img > 0], 99.5)
    z = [float(v) for v in zooms]
    ext = [(hi[k] - lo[k]) * z[k] for k in range(3)]                       # physical extent (mm)
    aspects = [z[2] / z[1], z[2] / z[0], z[1] / z[0]]                     # rows / columns of each displayed view
    fig, ax = plt.subplots(len(rows), 3, figsize=(13, 4.2 * len(rows)), facecolor='k', squeeze=False,
                           gridspec_kw=dict(width_ratios=[ext[1], ext[0], ext[0]]))
    for i, (name, seg) in enumerate(rows):
        ims = sl(img)
        ss = sl(classes(seg)) if seg is not None else [None] * 3
        for j in range(3):
            a = ax[i, j]
            a.imshow(ims[j], cmap='gray', vmin=0, vmax=vmax, interpolation='none', aspect=aspects[j])
            if ss[j] is not None:
                for k, col in enumerate(cols, 1):
                    if (ss[j] == k).any():
                        a.contour(ss[j] == k, levels=[0.5], colors=[col], linewidths=0.7)
            a.set_axis_off()
        ax[i, 0].set_title(name, color='w', loc='left', fontsize=9)
    fig.legend(handles=[matplotlib.patches.Patch(color=col, label=n) for n, (_, col) in QC_CLASSES.items()],
               loc='lower center', ncol=6, fontsize=8, facecolor='k', labelcolor='w')
    fig.suptitle(title, color='w', fontsize=11)
    fig.tight_layout(rect=(0, 0.015, 1, 0.985))
    os.makedirs(opd(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=110, facecolor='k')
    plt.close(fig)
    return out_png


# -------------------------------------------------------------------------------------------------------- pipeline
def find_channels(volumes_dir, ID, type_norm, targetsuffix='space-acpc_desc-SS'):
    """Reference image (type_norm) first, then the other of T1w/T2w when it exists on the same grid."""
    import nibabel as nib
    ref = opj(volumes_dir, f'{ID}_{targetsuffix}_{type_norm}.nii.gz')
    if not opi(ref):
        return []
    out = [ref]
    ref_img = nib.load(ref)
    for other in ('T1w', 'T2w'):
        if other == type_norm:
            continue
        f = opj(volumes_dir, f'{ID}_{targetsuffix}_{other}.nii.gz')
        if opi(f):
            im = nib.load(f)
            if im.shape[:3] == ref_img.shape[:3] and np.allclose(im.affine, ref_img.affine, atol=1e-3):
                out.append(f)
    return out


def run(ID, volumes_dir, labels_dir, bids_dir, type_norm, diary_file, params=None):
    """Relabel <labels_dir>/<ID>_seg-4FS_dseg.nii.gz in place. The registration-only seg, the report (json) and the
    QC are written to <labels_dir>/tissue_relabel/ and <bids_dir>/QC/tissue_relabel/."""
    import nibabel as nib
    from Tools import run_cmd
    seg_f = opj(labels_dir, f'{ID}_seg-4FS_dseg.nii.gz')
    chans = find_channels(volumes_dir, ID, type_norm)
    if not opi(seg_f) or not chans:
        run_cmd.msg(f'WARNING: tissue_relabel skipped, missing {seg_f if not opi(seg_f) else type_norm + " image"}',
                    diary_file, 'WARNING')
        return None
    work = opj(labels_dir, 'tissue_relabel')
    os.makedirs(work, exist_ok=True)
    orig_f = opj(work, f'{ID}_seg-4FS_desc-registration_dseg.nii.gz')
    shutil.copyfile(seg_f, orig_f)          # always from the fresh step-8 output
    seg_img = nib.load(seg_f)
    seg = np.asanyarray(seg_img.dataobj).astype(int)
    imgs = [np.asanyarray(nib.load(f).dataobj).astype(np.float32) for f in chans]
    if imgs[0].shape != seg.shape:
        run_cmd.msg(f'WARNING: tissue_relabel skipped, {chans[0]} and {seg_f} are not on the same grid',
                    diary_file, 'WARNING')
        return None
    mask = imgs[0] > 0
    out, rep = compute(seg, imgs, mask, seg_img.header.get_zooms()[:3], params)
    rep.update(seg=seg_f, registration_seg=orig_f, channels=chans)
    if rep['applied']:
        nib.save(nib.Nifti1Image(out.astype(np.int16), seg_img.affine, seg_img.header), seg_f)
        run_cmd.msg(f"INFO: tissue_relabel applied (GM/WM separation {rep['separation_gm_wm']:.2f}, "
                    f"{100 * rep['changed_fraction_ribbon']:.1f}% of the ribbon changed)", diary_file, 'OKGREEN')
    else:
        run_cmd.msg(f"WARNING: tissue_relabel not applied: {rep['reason']}", diary_file, 'WARNING')
    with open(opj(work, f'{ID}_tissue_relabel.json'), 'w') as f:
        json.dump(rep, f, indent=1, default=float)
    try:
        qc_figure(f'{ID}: tissue relabel ({"applied" if rep["applied"] else "NOT applied"})', imgs[0],
                  [(opb(chans[0]), None), ('registration (template seg-4FS warped)', seg),
                   ('tissue relabel', out)], opj(bids_dir, 'QC', 'tissue_relabel', f'{ID}_tissue_relabel.png'),
                  seg_img.header.get_zooms()[:3])
    except Exception as e:  # QC must never stop the pipeline
        run_cmd.msg(f'WARNING: tissue_relabel QC figure failed: {e}', diary_file, 'WARNING')
    return rep
