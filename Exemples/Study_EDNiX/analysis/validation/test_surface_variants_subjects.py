#!/usr/bin/env python3
"""
Surface-variant test on real EDNiX subjects (outputs of anat steps 1-8 already in the BIDS derivatives), same
variants as test_surface_variants_library.py (imported from the same folder). Stand-alone: reads the subject
files, writes everything under --out (default /scratch2/EDNiX/surface_tests/subjects). Nothing is written in
the datasets and the pipeline is not modified.

Per subject, the inputs are those of anat step 10 (_loop3.py):
  image : anat/native/volumes/<ID>_space-acpc_desc-SS[_res-iso]_T1w.nii.gz (the fake T1 EDNiX builds for T2
          datasets), else the T2w / acq-*_T2w of the same name
  aseg  : anat/native/volumes/labels/<ID>[_space-acpc]_seg-4FS[_res-iso]_dseg.nii.gz
  res-iso versions are preferred when present (EDNiX uses them when the native image is not isotropic).
  smallbrain rule from the brain volume with factor 20, as for subjects in _loop3.py.

Default variants: V0 (current EDNiX), V2m (legacy grid, real-contrast image), V2s (V2m + anti-aliased label
down-sampling), V6n (native grid, aseg smoothed sigma 2 voxels). V4n/V5n remain available with --variants.

Outputs:  <out>/summary.csv  (one row per species x subject x variant, thickness + fidelity + roughness)
          <out>/QC_<species>_<subject>.png

Run:  python3 test_surface_variants_subjects.py                         (default datasets, 2 subjects each)
      python3 test_surface_variants_subjects.py --species Rat Marmoset --n 3 --variants V0 V5n
      python3 test_surface_variants_subjects.py --list                  (subjects found, nothing run)
"""
import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_surface_variants_library as L   # noqa: E402

opj, opb, opd, opi = os.path.join, os.path.basename, os.path.dirname, os.path.isfile

# bids_dir of the anat launchers (Exemples/Study_EDNiX/Launcher_anat)
DATASETS = {
    'Mouse': '/scratch2/EDNiX/Mouse/BIDS_Grandjean2/',
    'Rat': '/scratch2/EDNiX/Rat/BIDS_Grandjean/',
    'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin/',
    'Marmoset': '/scratch2/EDNiX/Marmoset/BIDS_Hori/',
    'Macaque': '/scratch2/EDNiX/Macaque/BIDS_BenHamed/',
    'Dog': '/scratch2/EDNiX/Dog/BIDS_Boch_K9/',
}
SPECIES_DEFAULT = ['Mouse', 'Rat', 'Mouselemur', 'Marmoset', 'Macaque']
VARIANTS_DEFAULT = ['V0', 'V2m', 'V2s', 'V6n']
SUBJECT_FACTOR = 20          # smallbrain factor for subjects in _loop3.py
OUT_DEFAULT = '/scratch2/EDNiX/surface_tests/subjects'


def first_file(patterns):
    for p in patterns:
        hits = sorted(glob(p))
        if hits:
            return hits[0]
    return ''


def otsu(x, bins=256):
    h, e = np.histogram(x, bins)
    c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h)
    w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1)
    m1 = (np.sum(h * c) - np.cumsum(h * c)) / np.maximum(w1, 1)
    return c[np.argmax(w0 * w1 * (m0 - m1) ** 2)]


def fix_wm(img_file, aseg_file, out_file, rim_mm=0.4, max_frac=0.6):
    """Test correction of an aseg whose cortex label (3/42) swallows the underlying white matter (Marmoset
    template seg-4FS and every subject warped from it).
    1. N4 bias correction of the image inside the brain (a global threshold fails with the lateral bias of
       sub-02: first version relabelled 100 % of a hemisphere).
    2. Per hemisphere, Otsu threshold on cortex + WM label intensities clipped to their 1-99th percentiles
       (bright vessels / non-brain voxels no longer pull the threshold).
    3. Cortex voxels on the WM side of the threshold, deeper than rim_mm from the brain surface (a cortical
       rim is always kept) AND connected (26-n) to the WM label of the same hemisphere -> WM (3->2, 42->41).
    4. Guard: if more than max_frac of a hemisphere's cortex would be relabelled, that hemisphere is left
       unchanged and the report says so.
    Also writes <out_file>.png (WM label before = blue, added WM = yellow). Returns (out_file, report)."""
    import ants
    import nibabel as nib
    from scipy import ndimage
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    an = nib.load(aseg_file)
    A = np.round(np.asarray(an.dataobj)).astype(np.int32)
    brain = (A > 0) & (A != 24)
    img = ants.image_read(img_file)
    mask = img.new_image_like(brain.astype(np.float32))
    I = ants.n4_bias_field_correction(img, mask=mask).numpy().astype(np.float32)
    depth = ndimage.distance_transform_edt(brain, sampling=an.header.get_zooms()[:3])
    B = A.copy()
    rep = []
    for ctx_c, wm_c in ((3, 2), (42, 41)):
        ctx, wm = A == ctx_c, A == wm_c
        if not ctx.any() or not wm.any():
            continue
        v = I[ctx | wm]
        lo, hi = np.percentile(v, (1, 99))
        t = otsu(np.clip(v, lo, hi))
        bright_wm = np.median(I[wm]) > np.median(I[ctx])       # T1-like: WM brighter than cortex
        cand = ctx & (depth > rim_mm) & ((I > t) if bright_wm else (I < t))
        lab, _ = ndimage.label(cand | wm, structure=np.ones((3, 3, 3)))
        ids = np.unique(lab[wm])
        moved = cand & np.isin(lab, ids[ids > 0])
        frac = moved.sum() / ctx.sum()
        if frac > max_frac:
            rep.append(f'{ctx_c}->{wm_c}: REFUSED ({100 * frac:.0f} % of cortex > {100 * max_frac:.0f} %), unchanged')
            continue
        B[moved] = wm_c
        rep.append(f'{ctx_c}->{wm_c}: {int(moved.sum())} vox ({100 * frac:.0f} % of cortex), '
                   f'WM/cortex {wm.sum() / ctx.sum():.2f} -> {(B == wm_c).sum() / (B == ctx_c).sum():.2f}')
    nib.save(nib.Nifti1Image(B.astype(np.int16), an.affine, an.header), out_file)
    nz = np.argwhere(brain)
    lo3, hi3 = nz.min(0), nz.max(0)
    fig, ax = plt.subplots(1, 4, figsize=(16, 4.5))
    for k, frac in enumerate((.3, .5, .7, .55)):
        axis = 1 if k < 3 else 2
        idx = int(lo3[axis] + frac * (hi3[axis] - lo3[axis]))
        sl = [slice(None)] * 3
        sl[axis] = idx
        sl = tuple(sl)
        ax[k].imshow(I[sl].T, cmap='gray', origin='lower', vmin=np.percentile(I[brain], 1),
                     vmax=np.percentile(I[brain], 99))
        old, new = np.isin(A[sl], (2, 41)), np.isin(B[sl], (2, 41))
        ax[k].contour(np.isin(B[sl], (3, 42)).T, levels=[.5], colors='tab:red', linewidths=.5)
        ax[k].contour(old.T, levels=[.5], colors='tab:blue', linewidths=.6)
        ax[k].contour(new.T, levels=[.5], colors='yellow', linewidths=.6)
        ax[k].axis('off')
    fig.suptitle(opb(out_file) + ': cortex red, WM before blue, WM after yellow\n' + '; '.join(rep), fontsize=8)
    fig.tight_layout()
    fig.savefig(out_file.replace('.nii.gz', '.png'), dpi=110)
    plt.close(fig)
    return out_file, '; '.join(rep)


def find_subjects(bids_dir):
    """Subjects/sessions with the step-10 inputs. Returns [(tag, image, aseg)], tag = sub-X_ses-Y."""
    out = []
    for vol in sorted(glob(opj(bids_dir, '**', 'anat', 'native', 'volumes'), recursive=True)):
        ses_dir = opd(opd(opd(vol)))
        ses, sub = opb(ses_dir), opb(opd(ses_dir))
        if not (sub.startswith('sub-') and ses.startswith('ses-')):
            continue
        img = first_file([opj(vol, '*_space-acpc_desc-SS_res-iso_T1w.nii.gz'),
                          opj(vol, '*_space-acpc_desc-SS_T1w.nii.gz'),
                          opj(vol, '*_space-acpc_desc-SS_res-iso_*T1w.nii.gz'),
                          opj(vol, '*_space-acpc_desc-SS_*T1w.nii.gz'),
                          opj(vol, '*_space-acpc_desc-SS_res-iso_*T2w.nii.gz'),
                          opj(vol, '*_space-acpc_desc-SS_*T2w.nii.gz')])
        lab = opj(vol, 'labels')
        aseg = first_file([opj(lab, '*_seg-4FS_res-iso_dseg.nii.gz'),
                           opj(lab, '*_seg-4FS_dseg.nii.gz'),
                           opj(lab, '*_space-acpc_seg-4FS_dseg.nii.gz')])
        if img and aseg:
            out.append((f'{sub}_{ses}', img, aseg))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--main-path', default=L.MAIN_PATH_DEFAULT)
    ap.add_argument('--out', default=OUT_DEFAULT)
    ap.add_argument('--species', nargs='+', default=SPECIES_DEFAULT, choices=sorted(DATASETS))
    ap.add_argument('--variants', nargs='+', default=VARIANTS_DEFAULT, choices=L.VARIANTS)
    ap.add_argument('--n', type=int, default=2, help='subjects per species (random, fixed seed)')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--qc-only', action='store_true')
    ap.add_argument('--list', action='store_true', help='print the subjects found and stop')
    ap.add_argument('--fix-wm', action='store_true',
                    help='relabel bright cortex voxels connected to WM as WM (fix_wm) before FreeSurfer; '
                         'runs are tagged <subject>_fixwm, the original runs are kept')
    args = ap.parse_args()

    args.main_path = args.main_path.rstrip('/') + '/'
    lib = opj(opd(args.main_path), 'Atlases_library')
    os.makedirs(args.out, exist_ok=True)
    jobs, chosen = [], {}
    for sp in args.species:
        subs = find_subjects(DATASETS[sp])
        print(f'{sp}: {len(subs)} subjects with step-10 inputs in {DATASETS[sp]}')
        if not subs:
            continue
        subs = [subs[i] for i in np.random.default_rng(args.seed).permutation(len(subs))]   # reproducible pick
        chosen[sp] = subs[:args.n]
        if args.fix_wm:
            fixed = []
            for tag, img, aseg in chosen[sp]:
                d = opj(args.out, sp, tag + '_fixwm')
                os.makedirs(d, exist_ok=True)
                f, rep = fix_wm(img, aseg, opj(d, opb(aseg).replace('_dseg', '_desc-fixwm_dseg')))
                print(f'   fix_wm {tag}: {rep}')
                fixed.append((tag + '_fixwm', img, f))
            chosen[sp] = fixed
        for tag, img, aseg in chosen[sp]:
            print(f'   {tag}\n      {img}\n      {aseg}')
            for v in args.variants:
                jobs.append(dict(species=sp, variant=v, tag=tag, img=img, aseg=aseg, factor=SUBJECT_FACTOR,
                                 args=dict(main_path=args.main_path, lib=lib, out=args.out, force=args.force)))
    with open(opj(args.out, 'subjects_fixwm.json' if args.fix_wm else 'subjects.json'), 'w') as fh:
        json.dump(chosen, fh, indent=1)
    if args.list:
        return
    print(f'{len(jobs)} runs')

    rows = []
    if args.qc_only:
        for j in jobs:
            f = opj(args.out, j['species'], j['tag'], j['variant'], 'done.json')
            if opi(f):
                rows.append(json.loads(open(f).read()))
    else:
        lock = L.take_lock(args.out)
        with ProcessPoolExecutor(max(1, args.jobs)) as ex:
            futs = {ex.submit(L.run_one, j): j for j in jobs}
            for fu in as_completed(futs):
                r = fu.result()
                rows.append(r)
                print(f"[{len(rows)}/{len(jobs)}] {r['species']} {r.get('tag', '')} {r['variant']}: {r['status']} "
                      f"({r['seconds'] // 60} min) {r['note'][:120]}", flush=True)
        os.remove(lock)
    order = {v: i for i, v in enumerate(L.VARIANTS)}
    rows.sort(key=lambda r: (r['species'], r.get('tag', ''), order[r['variant']]))
    L.write_summary(args.out, rows)
    for sp, subs in chosen.items():
        for tag, _, _ in subs:
            try:
                L.qc_species(sp, args.out, rows, tag)
            except Exception as e:
                print(f'QC figure {sp} {tag} failed: {e}')
    print(f'\noutputs: {args.out}/summary.csv, {args.out}/QC_<species>_<subject>.png')


if __name__ == '__main__':
    main()
