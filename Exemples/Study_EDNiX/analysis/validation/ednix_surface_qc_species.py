"""
ednix_surface_qc_species.py -- EDNiX white (blue) and pial (red) surfaces on the subject's own
anatomy, several individuals per species, same slices within a species.

One figure per species: rows = individuals (latest session, spread over the species'
datasets of paper_data.ANAT_DATASETS), columns = 3 coronal + 1 axial + 1 sagittal slices.
Within a species every row uses the same slices, placed at the same fractions of each
brain's extent (30/50/70 % antero-posterior, 55 % inferior-superior; sagittal at 50 % of the
left hemisphere's own left-right extent, so it never falls in the interhemispheric fissure): positions are anatomical, not millimetres, so they also hold for the
"small brain" header scaling EDNiX applies to rodents. Row label: subject, dataset and
median vertex thickness (FreeSurfer nearest-point ?h.thickness over the cortex).

Anatomy: the whole-head (not skull-stripped) acpc T1 of EDNiX
(native/volumes/*desc-template*<T1>.nii.gz) shown in the geometry of the EDNiX FreeSurfer
subject (anat/freesurfer/<ID>/mri/orig.mgz, which carries the surfaces' coordinates).
orig.mgz is the skull-stripped acpc T1 converted by mri_convert; the voxel axis order /
flips between the NIfTI and orig.mgz are found by correlating the skull-stripped NIfTI
with orig.mgz over the 48 axis permutations x flips (r must be > 0.98, else the
skull-stripped orig.mgz is shown and the row says so).

Outputs  <PAPER_DIR>/Validation/EDNiX_surface_QC/EDNiX_surfaces_<species>.png

    python ednix_surface_qc_species.py              # 4 individuals per species
    python ednix_surface_qc_species.py --n 6
    python ednix_surface_qc_species.py --species Macaque Human
"""
import argparse
import glob
import itertools
import os
import sys

import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(0, MAIN_PATH)
from Exemples.Study_EDNiX.analysis import paper_data as pdat   # noqa: E402
from Exemples.Study_EDNiX.analysis import paper_style as ps    # noqa: E402
from Exemples.Study_EDNiX.analysis.validation import fs_ednix_matched_cortex as M  # noqa: E402

opj = os.path.join

OUT_DIR = opj(pdat.PAPER_DIR, 'Validation', 'EDNiX_surface_QC')
PLANES = [(1, 0.30), (1, 0.50), (1, 0.70), (2, 0.55), (0, 0.50)]   # (axis, fraction); sagittal = of lh
WHITE, PIAL = '#1F77B4', '#D62728'


# ----------------------------------------------------------------------------- inputs
def _ses_key(path):
    s = os.path.basename(path)[4:]
    return (0, int(s)) if s.isdigit() else (1, s)


def latest_sessions(bids_dir):
    """[(subject, session, EDNiX FreeSurfer subject dir)] -- latest session with surfaces."""
    out = []
    for sdir in sorted(glob.glob(opj(bids_dir, 'sub-*'))):
        for ses in sorted(glob.glob(opj(sdir, 'ses-*')), key=_ses_key, reverse=True):
            hits = sorted(glob.glob(opj(ses, 'anat', 'freesurfer', '*', 'surf', 'lh.pial')))
            if hits and os.path.isfile(opj(os.path.dirname(hits[0]), 'rh.pial')):
                out.append((os.path.basename(sdir)[4:], os.path.basename(ses)[4:],
                            os.path.dirname(os.path.dirname(hits[0]))))
                break
    return out


def pick_individuals(bids_dirs, n):
    """Up to n individuals, taken in turn from each dataset of the species."""
    pools = [[(pdat.bids_name(b),) + t for t in latest_sessions(b)] for b in bids_dirs]
    picked = []
    for k in range(max((len(p) for p in pools), default=0)):
        for p in pools:
            if k < len(p) and len(picked) < n:
                # spread over the dataset: 0, last, middle, ...
                idx = [0, len(p) - 1, len(p) // 2, len(p) // 4, 3 * len(p) // 4][k] if k < 5 else k
                cand = p[min(idx, len(p) - 1)]
                if cand not in picked:
                    picked.append(cand)
    return picked[:n]


def _match_axes(src, ref):
    """Axis permutation + flips turning array src into ref (same shape after permutation);
    returns (function, r)."""
    best = (None, -1.0)
    sub = tuple(slice(None, None, max(1, s // 64)) for s in ref.shape)
    refs = ref[sub].astype(np.float32).ravel()
    for perm in itertools.permutations(range(3)):
        if tuple(src.shape[p] for p in perm) != ref.shape:
            continue
        a = np.transpose(src, perm)
        for flips in itertools.product((False, True), repeat=3):
            b = a
            for ax, f in enumerate(flips):
                if f:
                    b = np.flip(b, ax)
            r = np.corrcoef(b[sub].astype(np.float32).ravel(), refs)[0, 1]
            if np.isfinite(r) and r > best[1]:
                best = ((perm, flips), r)
    if best[0] is None:
        return None, np.nan

    def fn(x, perm=best[0][0], flips=best[0][1]):
        y = np.transpose(x, perm)
        for ax, f in enumerate(flips):
            if f:
                y = np.flip(y, ax)
        return np.ascontiguousarray(y)
    return fn, best[1]


def anatomy(fs_subject):
    """(MGH image in the surfaces' geometry, note): whole-head T1 if it can be matched to
    orig.mgz, else the skull-stripped orig.mgz."""
    orig_f = opj(fs_subject, 'mri', 'orig.mgz')
    orig = nib.load(orig_f)
    ref = np.asarray(orig.dataobj)
    anat = os.path.dirname(os.path.dirname(fs_subject))           # .../ses-*/anat
    vol_dir = opj(anat, 'native', 'volumes')
    ss = [f for f in sorted(glob.glob(opj(vol_dir, '*desc-SS*.nii.gz')))
          if 'T2' not in os.path.basename(f) and 'step' not in os.path.basename(f)]
    head = [f for f in sorted(glob.glob(opj(vol_dir, '*desc-template*.nii.gz')))
            if 'T2' not in os.path.basename(f)]
    for f_ss in ss:
        img_ss = nib.load(f_ss)
        fn, r = _match_axes(np.asarray(img_ss.dataobj), ref)
        if fn is None or r < 0.98:
            continue
        for f_head in head:
            img_h = nib.load(f_head)
            if img_h.shape[:3] == img_ss.shape[:3] and np.allclose(img_h.affine, img_ss.affine, atol=1e-3):
                data = fn(np.asarray(img_h.dataobj)[..., 0] if img_h.ndim == 4 else np.asarray(img_h.dataobj))
                return (nib.MGHImage(data.astype(np.float32), orig.affine, orig.header),
                        f'{os.path.basename(f_head)} (r = {r:.3f})')
    return orig, 'orig.mgz, skull-stripped (whole-head T1 not matched)'


# ----------------------------------------------------------------------------- figure
def planes_for(hemis):
    """Slice positions at fixed fractions of the brain extent (sagittal: of the left
    hemisphere extent) and padded display ranges."""
    allv = np.vstack([d[k] for d in hemis.values() for k in ('white', 'pial')])
    lo, hi = allv.min(0), allv.max(0)
    lh = np.vstack([hemis['lh'][k] for k in ('white', 'pial')])
    lh_lo, lh_hi = lh.min(0), lh.max(0)
    pad = 0.08 * (hi - lo)
    planes = [(axis, (lh_lo[axis] + f * (lh_hi[axis] - lh_lo[axis])) if axis == 0
               else lo[axis] + f * (hi[axis] - lo[axis])) for axis, f in PLANES]
    return planes, list(zip(lo - pad, hi + pad))


def species_figure(species, rows, out):
    n = len(rows)
    fig, axes = plt.subplots(n, len(PLANES), figsize=(ps.WIDTH_2COL, (34 * n + 8) * ps.MM),
                             squeeze=False)
    for i, (dataset, sub, ses, fs_subject) in enumerate(rows):
        try:
            hemis = {h: M.read_hemi(fs_subject, h) for h in ('lh', 'rh')}
            img, note = anatomy(fs_subject)
            planes, ranges = planes_for(hemis)
            step = max(0.05, min(hi - lo for lo, hi in ranges) / 300)
            for j, (axis, value) in enumerate(planes):
                ax = axes[i, j]
                im, ext, h_ax, v_ax = M.slice_image(img, axis, value, ranges, step=step)
                vmax = np.percentile(im[im > 0], 99) if (im > 0).any() else 1
                ax.imshow(im, origin='lower', extent=ext, cmap='gray', vmin=0, vmax=vmax,
                          interpolation='nearest')
                for kind, col in (('white', WHITE), ('pial', PIAL)):
                    for d in hemis.values():
                        s = M.plane_segments(d[kind], d['faces'], axis, value)
                        if len(s):
                            ax.add_collection(M.LineCollection(s[:, :, [h_ax, v_ax]], colors=col,
                                                               linewidths=0.35))
                ax.set_xlim(ext[0], ext[1])
                ax.set_ylim(ext[2], ext[3])
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(f'{"sagittal (lh)" if axis == 0 else "coronal" if axis == 1 else "axial"} '
                                 f'{int(PLANES[j][1] * 100)} %', fontsize=5.5, pad=1)
            th = [d['th'][d['idx']] for d in hemis.values() if d['th'] is not None]
            med = f'{np.median(np.concatenate(th)):.2f} mm' if th else 'n/a'
            axes[i, 0].set_ylabel(f'sub-{sub} ses-{ses}\n{dataset.replace("BIDS_", "")}\n'
                                  f'thickness {med}', fontsize=5)
            axes[i, -1].text(1.02, 0.5, note, transform=axes[i, -1].transAxes, rotation=270,
                             va='center', ha='left', fontsize=3.5, color='0.4')
            print(f'  [{species}] sub-{sub} ses-{ses}: {note}')
        except Exception as e:  # noqa: BLE001 -- one bad subject must not stop the figure
            for ax in axes[i]:
                ax.axis('off')
            axes[i, 0].text(0, 0.5, f'sub-{sub} ses-{ses}: {e}', fontsize=5, transform=axes[i, 0].transAxes)
            print(f'  [{species}] sub-{sub} ses-{ses}: FAILED {e}')
    fig.suptitle(f'{species}: EDNiX white (blue) and pial (red) surfaces on the subject T1; '
                 f'same relative slices in every row', fontsize=6.5, y=0.998)
    fig.tight_layout(pad=0.3, w_pad=0.2, h_pad=0.4)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f'  [figure] {out}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--n', type=int, default=4, help='individuals per species')
    ap.add_argument('--species', nargs='*', default=None)
    a = ap.parse_args()
    plt.rcParams.update(ps.PAPER_RC)
    os.makedirs(OUT_DIR, exist_ok=True)
    for species in ps.phylo_sort(pdat.ANAT_DATASETS):
        if a.species and species not in a.species:
            continue
        rows = pick_individuals(pdat.ANAT_DATASETS[species], a.n)
        if not rows:
            print(f'  [{species}] no EDNiX FreeSurfer surfaces found')
            continue
        species_figure(species, rows, opj(OUT_DIR, f'EDNiX_surfaces_{species}.png'))
    print(f'-> {OUT_DIR}')


if __name__ == '__main__':
    main()
