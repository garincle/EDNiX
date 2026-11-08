"""
freesurfer_vs_ednix.py -- standard FreeSurfer recon-all on the raw T1w of the human Park
dataset, as the reference for the FreeSurfer vs EDNiX comparison (Supp. Fig. 14).

Human only (Clément, 2026-11-06): recon-all is not designed for non-human brains, so a
macaque run would test the header-scaling workaround, not EDNiX.

Dataset: /scratch2/EDNiX/Human/BIDS_Park, one session per subject = the latest (the EDNiX
longitudinal reference). recon-all -all, standard human settings, from the EDNiX container
(Singularity_library/Singularity/freesurfer_NHP.sif).

Restart-safe: a subject whose scripts/recon-all.done exists (and no recon-all.error) is
skipped. A subject interrupted mid-run (no .done) is deleted and restarted from scratch:
`recon-all -make all` could resume, but a file half-written at the crash would look up to
date. Subjects that failed (recon-all.error) are skipped unless --retry-failed.

Outputs: /scratch2/EDNiX/FreeSurfer_benchmark/BIDS_Park/   SUBJECTS_DIR (+ _work/<subject>/run.log)
         <PAPER_DIR>/Validation/FreeSurfer/fs_status.csv  recon-all outcome per subject

The comparison and the figure are in supp_freesurfer_comparison.py.

    python freesurfer_vs_ednix.py                  # 4 subjects in parallel
    python freesurfer_vs_ednix.py --jobs 6
    python freesurfer_vs_ednix.py --retry-failed   # also rerun subjects with recon-all.error
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(0, MAIN_PATH)
from Exemples.Study_EDNiX.analysis import paper_data as pdat  # noqa: E402

opj = os.path.join

SING_DIR = opj(MAIN_PATH, 'Singularity_library', 'Singularity')
FS_SIF = opj(SING_DIR, 'freesurfer_NHP.sif')

DATASETS = {'Human': '/scratch2/EDNiX/Human/BIDS_Park'}
FS_ROOT = '/scratch2/EDNiX/FreeSurfer_benchmark'
OUT_DIR = opj(pdat.PAPER_DIR, 'Validation', 'FreeSurfer')
HOSTLIB_DIR = opj(FS_ROOT, '_hostlibs')
LICENSE_CANDIDATES = [os.environ.get('FS_LICENSE', ''),
                      os.path.expanduser('~/license.txt'),
                      os.path.expanduser('~/.freesurfer/license.txt'),
                      opj(MAIN_PATH, 'license.txt'), opj(SING_DIR, 'license.txt'),
                      '/scratch2/EDNiX/license.txt']


def subjects_dir(bids):
    return opj(FS_ROOT, os.path.basename(bids.rstrip('/')))


# ----------------------------------------------------------------------------- inputs
def _ses_key(s):
    m = re.search(r'\d+', s)
    return (0, int(m.group()), s) if m else (1, 0, s)


def find_inputs(bids_dir):
    """[(subject, session, T1w path)] -- latest session of each subject, raw BIDS T1w only."""
    out = []
    for sdir in sorted(glob.glob(opj(bids_dir, 'sub-*'))):
        if not os.path.isdir(sdir):
            continue
        sub = os.path.basename(sdir)[4:]
        ses_dirs = sorted(glob.glob(opj(sdir, 'ses-*')), key=lambda p: _ses_key(os.path.basename(p)))
        cands = [(os.path.basename(d)[4:], opj(d, 'anat')) for d in ses_dirs] or [('1', opj(sdir, 'anat'))]
        for ses, adir in reversed(cands):
            t1 = [f for pat in ('*_T1w.nii*', '*_T1.nii*')
                  for f in sorted(glob.glob(opj(adir, pat)))
                  if 'space-' not in os.path.basename(f) and 'desc-' not in os.path.basename(f)]
            if t1:
                out.append((sub, ses, t1[0]))
                break
    return out


def license_bind():
    """(bind dirs, env prefix) for the FreeSurfer license; exits if none is usable."""
    for f in LICENSE_CANDIDATES:
        if f and os.path.isfile(f):
            return [os.path.dirname(f)], f'export SINGULARITYENV_FS_LICENSE="{f}"; '
    test = subprocess.run(['singularity', 'exec', FS_SIF, 'bash', '-c',
                           'test -f "${FS_LICENSE:-$FREESURFER_HOME/license.txt}"'])
    if test.returncode == 0:
        return [], ''
    sys.exit('No FreeSurfer license found (FS_LICENSE, ~/license.txt, the EDNiX folder or the '
             'container). Get one at https://surfer.nmr.mgh.harvard.edu/registration.html, '
             'save it as ~/license.txt and rerun.')


def _host_lib(name):
    """Path of a shared library on the host (ldconfig cache, then usual folders)."""
    try:
        out = subprocess.run(['/sbin/ldconfig', '-p'], capture_output=True, text=True).stdout
    except FileNotFoundError:
        out = ''
    for line in out.splitlines():
        if line.strip().startswith(name + ' ') and 'x86-64' in line:
            return line.split('=>')[-1].strip()
    for d in ('/usr/lib/x86_64-linux-gnu', '/lib/x86_64-linux-gnu', '/usr/lib64', '/lib64',
              opj(sys.prefix, 'lib'), os.path.expanduser('~/miniconda3/lib'),
              os.path.expanduser('~/anaconda3/lib')):
        f = opj(d, name)
        if os.path.isfile(f):
            return f
    return None


def container_libs():
    """The 4dfp tools of talairach_avi need libraries that freesurfer_NHP.sif does not ship
    (libquadmath.so.0). List what the container misses, copy it from the host into
    HOSTLIB_DIR, and return that folder (bound and put first in LD_LIBRARY_PATH)."""
    os.makedirs(HOSTLIB_DIR, exist_ok=True)
    probe = (f'export LD_LIBRARY_PATH={HOSTLIB_DIR}:$LD_LIBRARY_PATH; '
             'for b in $FREESURFER_HOME/bin/*4dfp* $FREESURFER_HOME/bin/mri_*; do '
             'ldd "$b" 2>/dev/null | grep "not found"; done | sort -u')
    missing = []
    for _ in range(3):
        r = subprocess.run(['singularity', 'exec', '--bind', HOSTLIB_DIR, FS_SIF, 'bash', '-c', probe],
                           capture_output=True, text=True)
        missing = sorted({l.split()[0] for l in r.stdout.splitlines() if 'not found' in l})
        if not missing:
            print(f'[libs] container libraries complete (extra folder: {HOSTLIB_DIR})')
            return HOSTLIB_DIR
        for name in missing:
            src = _host_lib(name)
            if src is None:
                sys.exit(f'{name} is missing in {FS_SIF} and not found on this machine. Install it '
                         f'(Debian: apt install libquadmath0 libgfortran5) or copy it into '
                         f'{HOSTLIB_DIR}, then rerun.')
            shutil.copy2(os.path.realpath(src), opj(HOSTLIB_DIR, name))
            print(f'[libs] {name} <- {src}')
    sys.exit(f'Container still misses libraries after copying them into {HOSTLIB_DIR}: {missing}')


def sh(cmd, log):
    with open(log, 'a') as fh:
        fh.write(f'\n$ {cmd}\n')
        fh.flush()
        r = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.STDOUT, executable='/bin/bash')
    return r.returncode


# ----------------------------------------------------------------------------- recon-all
def fs_state(sd, fsid):
    """'done', 'error' (recon-all.error, no .done), 'partial' (folder, no .done) or 'absent'."""
    scripts = opj(sd, fsid, 'scripts')
    if os.path.isfile(opj(scripts, 'recon-all.done')) and not os.path.isfile(opj(scripts, 'recon-all.error')):
        return 'done'
    if os.path.isfile(opj(scripts, 'recon-all.error')):
        return 'error'
    return 'partial' if os.path.isdir(opj(sd, fsid)) else 'absent'


def run_one(args):
    species, sub, ses, t1, sd, threads, lic_dirs, lic_env, libdir, retry_failed = args
    fsid = f'sub-{sub}_ses-{ses}'
    row = dict(species=species, subject=sub, session=ses, fsid=fsid)
    state = fs_state(sd, fsid)
    if state == 'done':
        return dict(row, status='done (cached)')
    if state == 'error' and not retry_failed:
        return dict(row, status='recon-all error (cached; --retry-failed to rerun)')
    work = opj(sd, '_work', fsid)
    os.makedirs(work, exist_ok=True)
    log = opj(work, 'run.log')
    bind = ','.join(sorted({os.path.dirname(t1), sd, SING_DIR, libdir, *lic_dirs}))
    try:
        with open(opj(work, 'scale_factor.txt'), 'w') as fh:
            fh.write('1.000000\n')
        if os.path.isdir(opj(sd, fsid)):
            shutil.rmtree(opj(sd, fsid))   # raises if NFS keeps open files (stale recon-all)
        cmd = (f'{lic_env}export SINGULARITYENV_SUBJECTS_DIR="{sd}"; '
               f'singularity exec --bind {bind} {FS_SIF} bash -c '
               f'"export LD_LIBRARY_PATH={libdir}:\\$LD_LIBRARY_PATH; '
               f'recon-all -s {fsid} -i {t1} -all -openmp {threads}"')
        print(f'[{species}] {fsid}: recon-all started ({state})', flush=True)
        rc = sh(cmd, log)
        status = 'ok' if rc == 0 else f'recon-all exit {rc}'
    except Exception as e:  # noqa: BLE001 -- the failure itself is recorded
        status = f'failed: {e}'
    print(f'[{species}] {fsid}: {status}', flush=True)
    return dict(row, status=status)


def run_all(jobs, retry_failed=False):
    lic_dirs, lic_env = license_bind()
    libdir = container_libs()
    threads = max(1, (os.cpu_count() or 4) // jobs)
    tasks = []
    for species, bids in DATASETS.items():
        sd = subjects_dir(bids)
        os.makedirs(sd, exist_ok=True)
        inputs = find_inputs(bids)
        states = [fs_state(sd, f'sub-{s}_ses-{e}') for s, e, _ in inputs]
        print(f'[{species}] {len(inputs)} subjects -> {sd}: '
              + ', '.join(f'{states.count(k)} {k}' for k in ('done', 'error', 'partial', 'absent')),
              flush=True)
        tasks += [(species, sub, ses, t1, sd, threads, lic_dirs, lic_env, libdir, retry_failed)
                  for sub, ses, t1 in inputs]
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        rows = list(ex.map(run_one, tasks))
    os.makedirs(OUT_DIR, exist_ok=True)
    st = pd.DataFrame(rows)
    st.to_csv(opj(OUT_DIR, 'fs_status.csv'), index=False)
    print(st.groupby(['species', 'status']).size().to_string())
    print(f'-> {opj(OUT_DIR, "fs_status.csv")}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--retry-failed', action='store_true')
    a = ap.parse_args()
    run_all(a.jobs, a.retry_failed)


if __name__ == '__main__':
    main()
