"""
Cross-species test of the tissue relabelling (modalities/anat/tissue_relabel.py) on subjects already processed by
EDNiX (step 8 done). Read-only on the BIDS folders: everything is written to OUT_DIR.

For each subject: registration seg-4FS (native) vs relabelled seg, on the subject's own T1w/T2w.
Metrics without ground truth (all species):
  - separation_gm_wm: Fisher distance GM vs WM of the fitted mixture (< 1 -> relabelling refused);
  - boundary_gradient before/after: mean image-gradient magnitude on the GM/WM boundary (higher = on the edge);
  - wm_connectivity before/after: fraction of cerebral WM in its largest component (thin WM broken -> lower);
  - volume changes of cortex / WM / ventricles, fraction of the ribbon changed;
  - wm_cortex_ratio of the registration seg (sanity check of the template labels).
Human (BIDS_Park): label Dice vs the recon-all aseg of the SyNadaptive calibration when available.
A QC figure is written as soon as each subject is done (OUT_DIR/QC/<species>_<subject>.png), then summary.csv and
summary.png per species.

    cd ~/Downloads && python test_tissue_relabel_species.py
    cd ~/Downloads && python test_tissue_relabel_species.py --species Macaque Marmoset --n 3
(tissue_relabel.py is taken from the EDNiX checkout in MAIN_PATH if it has it, otherwise from this folder.)
"""
import argparse
import glob
import importlib.util
import json
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

opj, opb, opd, opi = os.path.join, os.path.basename, os.path.dirname, os.path.isfile

MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
DATASETS = {
    'Human': '/scratch2/EDNiX/Human/BIDS_Park/',
    'Macaque': '/scratch2/EDNiX/Macaque/BIDS_BenHamed/',
    'Dog': '/scratch2/EDNiX/Dog/BIDS_Boch_K9/',
    'Marmoset': '/scratch2/EDNiX/Marmoset/BIDS_Hori/',
    'Mouselemur': '/scratch2/EDNiX/Mouselemur/BIDS_Garin/',
    'Rat': '/scratch2/EDNiX/Rat/BIDS_Grandjean/',
    'Mouse': '/scratch2/EDNiX/Mouse/BIDS_Grandjean2/',
}
SPECIES_DEFAULT = ['Human', 'Macaque', 'Marmoset', 'Rat']
N_DEFAULT = 2
RECONALL_CACHE = '/scratch2/EDNiX/SyNadaptive_calibration/BIDS_Park'   # <fsid>/fs_aseg_in_ednix_native.nii.gz
OUT_DIR = '/scratch2/EDNiX/tissue_relabel_tests'
N_JOBS = 6

STRUCTURES = {'cortex': (3, 42), 'cerebral_WM': (2, 41), 'lateral_ventricle': (4, 43), 'thalamus': (10, 49),
              'caudate': (11, 50), 'putamen': (12, 51), 'pallidum': (13, 52), 'hippocampus': (17, 53),
              'amygdala': (18, 54), 'cerebellum_cortex': (8, 47), 'cerebellum_WM': (7, 46), 'brainstem': (16,)}


def load_module():
    sys.path.insert(0, MAIN_PATH)
    try:
        from modalities.anat import tissue_relabel
        return tissue_relabel
    except ImportError:
        sys.path.remove(MAIN_PATH)
        here = opj(opd(os.path.abspath(__file__)), 'tissue_relabel.py')
        if not opi(here):
            sys.exit(f'tissue_relabel.py found neither in {MAIN_PATH}modalities/anat nor next to this script')
        spec = importlib.util.spec_from_file_location('tissue_relabel', here)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod


T = load_module()


def first(patterns):
    for p in patterns:
        hits = sorted(h for h in glob.glob(p) if 'desc-SS-step' not in h and 'res-iso' not in h)
        if hits:
            return hits[0]
    return ''


def find_subjects(bids_dir):
    """[(tag, [reference image, other contrast?], registration seg)] on the native (non res-iso) grid."""
    import nibabel as nib
    out = []
    for vol in sorted(glob.glob(opj(bids_dir, '**', 'anat', 'native', 'volumes'), recursive=True)):
        ses_dir = opd(opd(opd(vol)))
        ses, sub = opb(ses_dir), opb(opd(ses_dir))
        if not (sub.startswith('sub-') and ses.startswith('ses-')):
            continue
        t1 = first([opj(vol, '*_space-acpc_desc-SS_T1w.nii.gz'), opj(vol, '*_space-acpc_desc-SS_*T1w.nii.gz')])
        t2 = first([opj(vol, '*_space-acpc_desc-SS_T2w.nii.gz'), opj(vol, '*_space-acpc_desc-SS_*T2w.nii.gz')])
        lab = opj(vol, 'labels')
        seg = first([opj(lab, 'tissue_relabel', '*_seg-4FS_desc-registration_dseg.nii.gz'),
                     opj(lab, '*_seg-4FS_dseg.nii.gz')])
        ref = t1 or t2
        if not (ref and seg):
            continue
        chans = [ref]
        other = t2 if ref == t1 else ''
        if other:
            a, b = nib.load(ref), nib.load(other)
            if a.shape[:3] == b.shape[:3] and np.allclose(a.affine, b.affine, atol=1e-3):
                chans.append(other)
        out.append((f'{sub}_{ses}', chans, seg))
    return out


def label_dice(ref, test):
    out = {}
    for name, codes in STRUCTURES.items():
        vals = [2 * ((ref == c) & (test == c)).sum() / ((ref == c).sum() + (test == c).sum())
                for c in codes if (ref == c).any() and (test == c).any()]
        out[name] = float(np.mean(vals)) if vals else np.nan
    return out


def run_subject(args):
    species, tag, chans, seg_f = args
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[k] = '2'
    import nibabel as nib
    seg_img = nib.load(seg_f)
    seg = np.asanyarray(seg_img.dataobj).astype(int)
    imgs = [np.asanyarray(nib.load(f).dataobj).astype(np.float32) for f in chans]
    if imgs[0].shape != seg.shape:
        raise ValueError(f'{chans[0]} {imgs[0].shape} and {seg_f} {seg.shape} are not on the same grid')
    mask = imgs[0] > 0
    zooms = seg_img.header.get_zooms()[:3]
    rows = {}
    for name, ch in (('T1T2', imgs), ('ref', imgs[:1])) if len(imgs) > 1 else (('ref', imgs),):
        out, rep = T.compute(seg, ch, mask, zooms)
        rows[name] = (out, rep)
    out, rep = rows['T1T2'] if 'T1T2' in rows else rows['ref']
    row = dict(species=species, subject=tag, seg=seg_f, channels=';'.join(opb(c) for c in chans),
               wm_cortex_ratio=float(np.isin(seg, T.WM).sum() / max(np.isin(seg, T.CORTEX).sum(), 1)),
               **{k: v for k, v in rep.items() if not isinstance(v, (dict, list))})
    if 'T1T2' in rows:
        row['applied_ref_only'] = rows['ref'][1]['applied']
        row['boundary_gradient_after_ref_only'] = rows['ref'][1].get('boundary_gradient_after', np.nan)
    d = opj(OUT_DIR, species, tag)
    os.makedirs(d, exist_ok=True)
    nib.save(nib.Nifti1Image(out.astype(np.int16), seg_img.affine), opj(d, f'{tag}_seg-4FS_desc-relabel_dseg.nii.gz'))
    qc_rows = [(opb(chans[0]), None)] + ([(opb(chans[1]), None)] if len(chans) > 1 else [])
    ref = None
    fa = opj(RECONALL_CACHE, tag, 'fs_aseg_in_ednix_native.nii.gz')
    if species == 'Human' and opi(fa):
        ref = np.asanyarray(nib.load(fa).dataobj).astype(int)
        if ref.shape == seg.shape:
            for nm, s in (('registration', seg), ('relabel', out)):
                ld = label_dice(ref, s)
                row[f'dice_{nm}_mean'] = float(np.nanmean(list(ld.values())))
                for k in ('cortex', 'cerebral_WM', 'lateral_ventricle'):
                    row[f'dice_{nm}_{k}'] = ld[k]
            qc_rows.append(('recon-all aseg', ref))
        else:
            ref = None

    def t(nm):
        if nm == 'registration' and ref is not None:
            return f"registration: Dice {row['dice_registration_mean']:.3f} | cortex {row['dice_registration_cortex']:.3f}"
        if nm == 'relabel' and ref is not None:
            return f"relabel: Dice {row['dice_relabel_mean']:.3f} | cortex {row['dice_relabel_cortex']:.3f}"
        return nm

    qc_rows += [(t('registration') + f" | grad {rep.get('boundary_gradient_before', np.nan):.3f}", seg),
                (t('relabel') + (f" | grad {rep.get('boundary_gradient_after', np.nan):.3f} | sep "
                                 f"{rep.get('separation_gm_wm', np.nan):.2f}" if rep['applied'] else
                                 f" NOT APPLIED: {rep.get('reason', '')}"), out)]
    qc = T.qc_figure(f'{species} {tag}: tissue relabel (contours on {opb(chans[0])})', imgs[0], qc_rows,
                     opj(OUT_DIR, 'QC', f'{species}_{tag}.png'), zooms)
    row['qc'] = qc
    with open(opj(d, 'result.json'), 'w') as f:
        json.dump(row, f, indent=1, default=float)
    return row


def line(r):
    if not r.get('applied'):
        return f"NOT applied ({r.get('reason', '')})"
    s = (f"sep {r['separation_gm_wm']:.2f} | thick {r['cortical_thickness_mm']:.2f} mm | ribbon changed "
         f"{100 * r['changed_fraction_ribbon']:.1f}% | grad {r['boundary_gradient_before']:.3f}->"
         f"{r['boundary_gradient_after']:.3f} | WM conn {r['wm_connectivity_before']:.3f}->"
         f"{r['wm_connectivity_after']:.3f} | WM/cortex {r['wm_cortex_ratio']:.2f}")
    if 'dice_relabel_mean' in r:
        s += f" | Dice {r['dice_registration_mean']:.3f}->{r['dice_relabel_mean']:.3f}"
    return s


def summarise(rows):
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    df = pd.DataFrame(rows)
    df.to_csv(opj(OUT_DIR, 'summary.csv'), index=False)
    ok = df[df['applied'] == True]  # noqa: E712
    if ok.empty:
        return
    metrics = [('boundary_gradient', 'GM/WM boundary gradient'), ('wm_connectivity', 'WM largest component')]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for a, (m, lab) in zip(ax, metrics):
        for i, (sp, g) in enumerate(ok.groupby('species', sort=False)):
            for _, r in g.iterrows():
                a.plot([i - 0.15, i + 0.15], [r[f'{m}_before'], r[f'{m}_after']], '-o', color='0.4', ms=3)
        a.set_xticks(range(ok['species'].nunique()), ok['species'].unique(), rotation=30)
        a.set_title(f'{lab} (left: registration, right: relabel)', fontsize=9)
    sp = ok.groupby('species', sort=False)
    ax[2].bar(range(len(sp)), sp['changed_fraction_ribbon'].mean() * 100, color='0.6')
    ax[2].set_xticks(range(len(sp)), list(sp.groups), rotation=30)
    ax[2].set_title('% of the ribbon relabelled (mean)', fontsize=9)
    fig.tight_layout()
    fig.savefig(opj(OUT_DIR, 'summary.png'), dpi=120)
    plt.close(fig)
    cols = ['separation_gm_wm', 'cortical_thickness_mm', 'changed_fraction_ribbon', 'boundary_gradient_before',
            'boundary_gradient_after', 'wm_connectivity_before', 'wm_connectivity_after', 'wm_cortex_ratio']
    cols += [c for c in ('dice_registration_mean', 'dice_relabel_mean') if c in df]
    print('\nper species (mean):\n' + df.groupby('species', sort=False)[cols].mean().round(3).to_string())


def main():
    global OUT_DIR
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--species', nargs='+', default=SPECIES_DEFAULT, choices=sorted(DATASETS))
    ap.add_argument('--n', type=int, default=N_DEFAULT, help='subjects per species (0 = all)')
    ap.add_argument('--out', default=OUT_DIR)
    a = ap.parse_args()
    OUT_DIR = a.out
    os.makedirs(opj(OUT_DIR, 'QC'), exist_ok=True)
    print(f'tissue_relabel module: {T.__file__}')
    jobs = []
    for sp in a.species:
        subs = find_subjects(DATASETS[sp])
        if sp == 'Human':   # subjects with a recon-all reference first
            subs.sort(key=lambda s: not opi(opj(RECONALL_CACHE, s[0], 'fs_aseg_in_ednix_native.nii.gz')))
        subs = subs[:a.n] if a.n else subs
        print(f'{sp}: {len(subs)} subject(s) from {DATASETS[sp]}')
        jobs += [(sp, *s) for s in subs]
    rows = []
    with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
        futs = {ex.submit(run_subject, j): j for j in jobs}
        for f in as_completed(futs):
            sp, tag = futs[f][:2]
            try:
                r = f.result()
                rows.append(r)
                print(f'{sp} {tag}: {line(r)}\n   QC: {r["qc"]}', flush=True)
            except Exception:
                print(f'{sp} {tag}: FAILED\n{traceback.format_exc()}', flush=True)
    if rows:
        summarise(rows)
    print(f'Outputs: {OUT_DIR}')


if __name__ == '__main__':
    main()
