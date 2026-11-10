#!/usr/bin/env python3
"""
rebuild_all_ednix.py: rebuild EVERYTHING with one command and ONE confirmation.

    cd ~/PycharmProjects/EDNiX && python Exemples/Study_EDNiX/rebuild_all_ednix.py            # preview, YES, go
    python Exemples/Study_EDNiX/rebuild_all_ednix.py --plan                                 # preview only
    python Exemples/Study_EDNiX/rebuild_all_ednix.py --only Rat,Mouse --phases anat,func    # subset of launchers

In the terminal (interactive, a few minutes):
  0. shows the EDNiX commit, the library build plan (build_library_V2.py), the jobs of run_all_ednix.py and the dry
     run of the BIDS reset (manual masks and crops kept), then asks ONE 'YES' for all of it;
  then it detaches (keeps running after you close the terminal / ssh) and logs to run_reports/rebuild_<date>.log:
  1. builds Atlases_library_V2 (surface method V2h). STOPS if a species fails;
  2. verifies V2 against its manifest (sha256). STOPS if a file differs;
  3. swaps the libraries: Atlases_library -> Atlases_library_V1, Atlases_library_V2 -> Atlases_library;
  4. resets the BIDS folders of the selected launchers (what was previewed in step 0);
  5. runs run_all_ednix.py (all launchers as a dependency graph on all cores, then the QC reports).
Resumable: if Atlases_library is already the V2 (manifest), steps 1-3 are skipped. --skip-library uses the current
library as it is. Nothing is deleted before the library is built and verified.
"""
import argparse
import datetime
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAIN_PATH = HERE.parents[1]
sys.path.insert(0, str(MAIN_PATH))
PYTHON = sys.executable
REPORT_DIR = MAIN_PATH / 'run_reports'
LIB_SCRIPT = HERE / 'atlases' / 'build_library_V2.py'
RUN_ALL = HERE / 'run_all_ednix.py'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L = load(LIB_SCRIPT, 'build_library_V2')
R = load(RUN_ALL, 'run_all_ednix')
LIB = Path(L.LIBRARY_V1)            # <EDNiX>/Atlases_library (what the launchers use)
LIB_V2 = Path(L.LIBRARY_V2)         # <EDNiX>/Atlases_library_V2 (built here)


def say(msg):
    print(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def manifest_version(lib):
    f = Path(lib) / 'EDNiX_library_manifest.json'
    try:
        return json.loads(f.read_text()).get('library')
    except (OSError, ValueError):
        return None


def library_is_current():
    return manifest_version(LIB) == L.VERSION


def launchers(only, phases):
    return [x for x in R.discover(only) if x['kind'] in phases]


def preview(a):
    only, phases = a.only, a.phases
    git = L.git_info(str(MAIN_PATH))
    print(f"EDNiX {git['branch']} @ {git['commit'][:9]}" + (' (UNCOMMITTED CHANGES)' if git['uncommitted_changes'] else ''))
    print('\n=== 1-3. atlas library ===')
    if a.skip_library:
        print(f'skipped (--skip-library): the launchers use {LIB} as it is')
    elif library_is_current():
        print(f'{LIB} is already the {L.VERSION} library: build, verify and swap skipped')
    else:
        L.build(L.SPECIES, plan=True)
        print(f'-> builds {LIB_V2}, verifies it, then {LIB} -> {LIB}_V1 and {LIB_V2} -> {LIB}')
    print('\n=== 4. reset of the BIDS folders (dry run) ===')
    ls = launchers(only, phases)
    n_ses, n_files, n_func = (R.reset_all(ls, preview_only=True) if not a.no_reset else (0, 0, 0))
    if a.no_reset:
        print('skipped (--no-reset)')
    print('\n=== 5. run_all_ednix ===')
    work = REPORT_DIR / '_preview'
    work.mkdir(parents=True, exist_ok=True)
    jobs = R.make_jobs(ls, work)
    print(R.plan_text(jobs, R.dependencies(jobs)))
    shutil.rmtree(work, ignore_errors=True)
    return n_ses, n_files, n_func, len(jobs)


def build_and_swap():
    if library_is_current():
        say(f'{LIB} is already the {L.VERSION} library: steps 1-3 skipped')
        return
    say(f'1. building {LIB_V2} (surface method {L.SURFACE["method"]})')
    if subprocess.run([PYTHON, str(LIB_SCRIPT)], cwd=str(HERE / 'atlases')).returncode:
        sys.exit('STOP: build_library_V2.py failed (BIDS folders untouched).')
    man = json.loads((LIB_V2 / 'EDNiX_library_manifest.json').read_text())
    failed = [k for k, v in man.get('species', {}).items() if v.get('status') != 'ok']
    if failed or not man.get('species'):
        sys.exit(f'STOP: library species failed: {failed or "none built"} (see {LIB_V2}; BIDS folders untouched).')
    say('2. verifying the V2 library against its manifest')
    if subprocess.run([PYTHON, str(LIB_SCRIPT), '--verify'], cwd=str(HERE / 'atlases')).returncode:
        sys.exit('STOP: the V2 library differs from its manifest (BIDS folders untouched).')
    old = LIB.with_name(LIB.name + '_V1')
    if old.exists():
        old = LIB.with_name(LIB.name + '_V1_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    say(f'3. {LIB} -> {old}; {LIB_V2} -> {LIB}')
    LIB.rename(old)
    LIB_V2.rename(LIB)


def run(a):
    say(f'rebuild started (pid {os.getpid()})')
    if not a.skip_library:
        build_and_swap()
    ls = launchers(a.only, a.phases)
    if not a.no_reset:
        say('4. resetting the BIDS folders (confirmed at launch)')
        R.reset_all(ls, confirm=False)
    say('5. run_all_ednix')
    cmd = [PYTHON, str(RUN_ALL), '--only', ','.join(a.only), '--phases', ','.join(a.phases)]
    code = subprocess.run(cmd, cwd=str(MAIN_PATH)).returncode
    say(f'finished (run_all_ednix exit code {code}); reports in {REPORT_DIR}')
    sys.exit(code)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--only', default=','.join(R.ONLY), help='launcher name fragments, comma-separated')
    ap.add_argument('--phases', default=','.join(R.PHASES), help='anat,ct,pet,func,qc')
    ap.add_argument('--plan', action='store_true', help='preview only, ask nothing, change nothing')
    ap.add_argument('--skip-library', action='store_true', help='keep the current Atlases_library as it is')
    ap.add_argument('--no-reset', action='store_true', help='do not reset the BIDS folders')
    ap.add_argument('--confirmed', action='store_true', help=argparse.SUPPRESS)   # detached worker
    a = ap.parse_args()
    a.only = [x for x in a.only.split(',') if x]
    a.phases = [x for x in a.phases.split(',') if x]
    if a.confirmed:
        return run(a)
    n_ses, n_files, n_func, n_jobs = preview(a)
    if a.plan:
        return
    lib = 'keep the current library' if a.skip_library or library_is_current() else \
        f'build + verify the {L.VERSION} library and swap it in'
    rst = 'no reset' if a.no_reset else f'reset {n_ses} anat/PET/CT session(s) ({n_files} files) and {n_func} func'
    ans = input(f'\n🚨 {lib}, then {rst}, then run {n_jobs} job(s). Type YES: ')
    if ans.strip() != 'YES':
        print('Cancelled: nothing changed.')
        return
    REPORT_DIR.mkdir(exist_ok=True)
    log = REPORT_DIR / f"rebuild_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    args = [PYTHON, '-u', str(Path(__file__).resolve()), '--confirmed', '--only', ','.join(a.only),
            '--phases', ','.join(a.phases)] + (['--skip-library'] if a.skip_library else []) + \
        (['--no-reset'] if a.no_reset else [])
    with open(log, 'w') as fh:
        p = subprocess.Popen(args, stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                             cwd=str(MAIN_PATH), start_new_session=True)
    print(f'\nRunning in the background (pid {p.pid}); you can close the terminal.\n  tail -f {log}')


if __name__ == '__main__':
    main()
