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

Default variants: V0 (current EDNiX), V2m (legacy grid, real-contrast image), V4n (native grid, white on V2 image,
pial on constant-class image), V5n / V6n (V4n on the aseg smoothed by sigma 1 / 2 voxels).

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
VARIANTS_DEFAULT = ['V0', 'V2m', 'V4n', 'V5n', 'V6n']
SUBJECT_FACTOR = 20          # smallbrain factor for subjects in _loop3.py
OUT_DEFAULT = '/scratch2/EDNiX/surface_tests/subjects'


def first_file(patterns):
    for p in patterns:
        hits = sorted(glob(p))
        if hits:
            return hits[0]
    return ''


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
        for tag, img, aseg in chosen[sp]:
            print(f'   {tag}\n      {img}\n      {aseg}')
            for v in args.variants:
                jobs.append(dict(species=sp, variant=v, tag=tag, img=img, aseg=aseg, factor=SUBJECT_FACTOR,
                                 args=dict(main_path=args.main_path, lib=lib, out=args.out, force=args.force)))
    with open(opj(args.out, 'subjects.json'), 'w') as fh:
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
        with ProcessPoolExecutor(max(1, args.jobs)) as ex:
            futs = {ex.submit(L.run_one, j): j for j in jobs}
            for fu in as_completed(futs):
                r = fu.result()
                rows.append(r)
                print(f"[{len(rows)}/{len(jobs)}] {r['species']} {r.get('tag', '')} {r['variant']}: {r['status']} "
                      f"({r['seconds'] // 60} min) {r['note'][:120]}", flush=True)
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
