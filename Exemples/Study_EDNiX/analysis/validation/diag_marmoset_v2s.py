#!/usr/bin/env python3
"""1) V2m vs V2s: label volumes in the FS inputs + thickness per hemisphere (find what V2s broke).
2) Marmoset aseg overlaid on the T1 (is the error in the aseg or in the surface step?).
Output: /scratch2/EDNiX/surface_tests/subjects/diag_marmoset_aseg.png + printed tables."""
import json
import os
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/scratch2/EDNiX/surface_tests/subjects'
subs = json.load(open(os.path.join(OUT, 'subjects.json')))

print('=== V2m vs V2s: voxel counts in mri/*.mgz and median thickness (real mm) per hemisphere ===')
for sp in ('Marmoset', 'Rat', 'Mouse'):
    for tag, _, _ in subs[sp]:
        for v in ('V2m', 'V2s'):
            d = os.path.join(OUT, sp, tag, v, 'freesurfer', sp)
            if not os.path.isdir(d):
                continue
            row = []
            f = nib.load(os.path.join(d, 'mri', 'filled.mgz')).get_fdata()
            row += [f'filled127={int((f == 127).sum())}', f'filled255={int((f == 255).sum())}']
            a = np.round(nib.load(os.path.join(d, 'mri', 'aseg.mgz')).get_fdata()).astype(int)
            row += [f'aseg{c}={int((a == c).sum())}' for c in (2, 3, 41, 42)]
            w = nib.load(os.path.join(d, 'mri', 'wm.mgz')).get_fdata()
            row += [f'wm110={int((w == 110).sum())}', f'wm250={int((w == 250).sum())}', f'wm_other={int(((w > 0) & (w != 110) & (w != 250)).sum())}']
            done = json.load(open(os.path.join(OUT, sp, tag, v, 'done.json')))
            k = float(done['real_mm_per_fs_mm'])
            for h in ('lh', 'rh'):
                t = os.path.join(d, 'surf', h + '.thickness')
                if os.path.isfile(t):
                    th = nib.freesurfer.read_morph_data(t) * k
                    row.append(f'{h}_thick={np.median(th):.2f} (p05 {np.percentile(th, 5):.2f})')
            print(f'{sp:9s} {tag:30s} {v}: ' + '  '.join(row))

print('\n=== Marmoset aseg overlay ===')
cols = {2: 'tab:blue', 41: 'tab:blue', 3: 'tab:red', 42: 'tab:red', 17: 'tab:green', 53: 'tab:green',
        18: 'tab:orange', 54: 'tab:orange'}
fig, ax = plt.subplots(len(subs['Marmoset']), 4, figsize=(16, 4.2 * len(subs['Marmoset'])), squeeze=False)
for r, (tag, img, aseg) in enumerate(subs['Marmoset']):
    I = nib.as_closest_canonical(nib.load(img)).get_fdata()
    A = np.round(nib.as_closest_canonical(nib.load(aseg)).get_fdata()).astype(int)
    print(tag, 'T1', I.shape, 'aseg', A.shape, 'labels', np.unique(A)[:40])
    nz = np.argwhere(A > 0)
    lo, hi = nz.min(0), nz.max(0)
    for c, (axis, frac) in enumerate([(1, .30), (1, .50), (1, .70), (2, .55)]):
        idx = int(lo[axis] + frac * (hi[axis] - lo[axis]))
        sl = [slice(None)] * 3
        sl[axis] = idx
        im, lab = I[tuple(sl)].T, A[tuple(sl)].T
        ax[r, c].imshow(im, cmap='gray', origin='lower', vmin=np.percentile(I[I > 0], 1), vmax=np.percentile(I[I > 0], 99))
        for code, col in cols.items():
            if (lab == code).any():
                ax[r, c].contour(lab == code, levels=[.5], colors=col, linewidths=.6)
        ax[r, c].set_title(f'{tag} {"coronal" if axis == 1 else "axial"} {int(frac * 100)} %', fontsize=8)
        ax[r, c].axis('off')
fig.suptitle('Marmoset aseg on T1: WM blue, cortex red, hippocampus green, amygdala orange', fontsize=10)
fig.tight_layout()
png = os.path.join(OUT, 'diag_marmoset_aseg.png')
fig.savefig(png, dpi=130)
print('figure:', png)
