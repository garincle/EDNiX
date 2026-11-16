"""Dataset-level boundary check of the PET -> anatomy registrations already computed.

Step 6 of the PET pipeline already runs the boundary check (modalities/pet/boundary.py) after every registration to an
anatomy, session by session. This script re-reads those registrations without running the pipeline again and adds the
dataset-level summary: mean and 95 % confidence interval (Student t) over the sessions of the dilation and of the
translation along each axis. A systematic offset of the dataset shows as an interval that excludes 0, even when every
session is within the tolerance. Per-session results are written next to each registration record
(<session>_desc-registration_boundary.json / .png) and in <bids_dir>/QC/pet_boundary/.

    python Tools/pet_boundary_summary.py /path/to/BIDS_dataset [/path/to/other_BIDS ...] [--fwhm 4.5]

--fwhm: PET resolution (FWHM, mm) of the scanner + reconstruction; default: twice the largest PET voxel.
"""
import argparse
import glob
import os
import sys

import numpy as np

opj = os.path.join
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from modalities.pet import boundary  # noqa: E402


def report(rec, pet_fwhm, bids_dir):
    try:
        res = boundary.check_record(rec, pet_fwhm, bids_dir)
    except Exception as e:
        print(f'{rec}: not checked ({e})')
        return None
    sh = ', '.join(f'{k} {v:+.1f}' for k, v in res['shifts_mm'].items())
    print(f"{os.path.basename(rec)}\n    PET edge - MRI edge (mm): {sh}\n    clean sectors (no extracerebral uptake "
          f"beside): {res['clean_sectors']}; dilation {res['dilation_mm']} mm, translation {res['translation_mm']} mm "
          f"(FWHM {res['fwhm_mm']:.1f} mm, model {res['model']})\n    -> {res['verdict']}")
    return res


def group_summary(results):
    """Mean and 95 % confidence interval (Student t) over the sessions of the dilation and of the translation along
    each axis, where measured."""
    from scipy import stats
    rows = [('dilation', [r['dilation_mm'] for r in results if r['dilation_mm'] is not None])]
    rows += [('translation ' + a, [r['translation_mm'][a] for r in results if a in r['translation_mm']])
             for a in 'xyz']
    for name, v in rows:
        if len(v) < 3:
            print(f'    {name}: {len(v)} session(s) measured, no group estimate')
            continue
        m, se = float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(len(v)))
        h = float(stats.t.ppf(0.975, len(v) - 1)) * se
        print(f'    {name}: {m:+.2f} mm, 95 % CI [{m - h:+.2f}, {m + h:+.2f}] (n = {len(v)})')


def records_of(bids_dir):
    return sorted(glob.glob(opj(bids_dir, 'sub-*', '**', 'pet', 'orig', 'matrices', '*_desc-registration_pet.json'),
                            recursive=True))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('bids_dirs', nargs='+', help='BIDS dataset(s) already processed by the PET pipeline')
    ap.add_argument('--fwhm', type=float, default=None,
                    help='PET resolution (FWHM, mm); default: twice the largest PET voxel')
    args = ap.parse_args(argv)
    for bids_dir in args.bids_dirs:
        records = records_of(bids_dir)
        print(f'{bids_dir}: {len(records)} registration records')
        results = [r for r in (report(rec, args.fwhm, bids_dir) for rec in records) if r is not None]
        print(f'{bids_dir}: group summary')
        group_summary(results)


if __name__ == '__main__':
    main()
