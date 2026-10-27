"""
Why is the mean FD missing for some datasets? -- diagnostic, reads files only.

The functional QC (modalities/fMRI/_12_fMRI_QC.py) computes the Jenkinson FD from
  sub-*/[ses-*/]func/raw/matrices/<run>_space-func_desc-matrix_motion_params.1D
and writes sub-*/[ses-*/]func/QC/<run>_QC_values.json (EDNiX writes inside the BIDS tree,
not under derivatives/)
which must hold 12 values per volume (affine [R | t], written by step 1 since the
motion-matrix fix). When the file is missing or has another number of columns, the
QC writes mean_fd = NaN and the diary says "ERROR computing FD".

For every dataset of paper_data.FUNC_DATASETS this script counts, per run:
  qc_json      functional QC json found (*_QC_values.json under a func QC folder)
  mean_fd      mean_fd finite in that json
  motion_1D    *_desc-matrix_motion_params.1D found for that run
  n_columns    columns of that file (12 expected)
and prints, per dataset, the first 'ERROR computing FD' line found in the diaries.
Output: one line per dataset on screen + <SENSITIVITY_DIR>/fmri_qc/fd_inputs.xlsx
"""
import os
import sys
import glob
import json

import numpy as np
import pandas as pd

sys.path.insert(0, '/home/cgarin/PycharmProjects/EDNiX/')

from Exemples.Study_EDNiX.analysis import paper_data as pdat

OUT = os.path.join(pdat.SENSITIVITY_DIR, 'fmri_qc', 'fd_inputs.xlsx')


def _n_columns(path):
    try:
        a = np.atleast_2d(np.genfromtxt(path))
        return int(a.shape[1])
    except Exception as e:
        return f"unreadable ({e})"


def _finite(x):
    try:
        return bool(np.isfinite(float(x)))
    except (TypeError, ValueError):
        return False


def _func_files(bids_dir, sub_path):
    """Files sub-*/ses-*/func/<sub_path> and sub-*/func/<sub_path> (datasets without sessions)."""
    out = []
    for pattern in (('sub-*', 'ses-*', 'func'), ('sub-*', 'func')):
        out += glob.glob(os.path.join(bids_dir, *pattern, *sub_path))
    return sorted(set(out))


def _first_error(bids_dir, pattern="ERROR computing FD"):
    # diaries: func/<ID> session <S><date>.txt (Tools/diaryfile.create)
    for diary in _func_files(bids_dir, ('*session*.txt',)):
        try:
            with open(diary, errors='ignore') as fh:
                for line in fh:
                    if pattern in line:
                        return f"{os.path.relpath(diary, bids_dir)}: {line.strip()[:200]}"
        except OSError:
            continue
    return ""


def check_dataset(species, bids_dir):
    rows = []
    motion = {os.path.basename(p).split('_space-func')[0]: p for p in
              _func_files(bids_dir, ('raw', 'matrices', '*_space-func_desc-matrix_motion_params.1D'))}
    qcs = _func_files(bids_dir, ('QC', '*_QC_values.json'))
    for q in qcs:
        root = os.path.basename(q).replace('_QC_values.json', '')
        try:
            with open(q) as fh:
                fd = json.load(fh).get('mean_fd')
        except Exception:
            fd = None
        m = motion.get(root)
        rows.append(dict(species=species, bids_dir=pdat.bids_name(bids_dir), run=root,
                         qc_json=q, mean_fd=fd,
                         mean_fd_ok=_finite(fd),
                         motion_1D=m, n_columns=_n_columns(m) if m else None))
    if not qcs:
        rows.append(dict(species=species, bids_dir=pdat.bids_name(bids_dir), run=None,
                         qc_json=None, mean_fd=None, mean_fd_ok=False,
                         motion_1D=f"{len(motion)} files", n_columns=None))
    return rows


if __name__ == "__main__":
    all_rows = []
    for sp, bids_dirs in pdat.FUNC_DATASETS.items():
        for bd in bids_dirs:
            rows = check_dataset(sp, bd)
            all_rows += rows
            df = pd.DataFrame(rows)
            n = len(df[df['qc_json'].notna()])
            n_fd = int(df['mean_fd_ok'].sum())
            n_mot = int(df['motion_1D'].notna().sum()) if n else 0
            cols = sorted({str(c) for c in df['n_columns'].dropna()})
            print(f"{sp:11s} {pdat.bids_name(bd):18s} QC runs {n:4d} | mean_fd finite {n_fd:4d} | "
                  f"motion_params.1D {n_mot:4d} | columns {', '.join(cols) or '-'}")
            err = _first_error(bd)
            if err:
                print(f"{'':31s}first FD error: {err}")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    pd.DataFrame(all_rows).to_excel(OUT, index=False)
    print(f"  [table] {OUT}")
