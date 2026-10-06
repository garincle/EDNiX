#!/usr/bin/env python3
"""
ednix_run.py — EDNiX parallel orchestrator
===========================================
Runs all EDNiX launchers in parallel, anat first then func, maximising
machine usage by splitting large datasets across subject-level workers.

KEY IDEA
--------
Your launchers already use fmri_resume_groups / anat_resume_groups, which:
  - compute resume groups over ALL subjects in the BIDS dir
  - skip whatever is already finished

This orchestrator does NOT modify your launchers. For large datasets it
generates temporary worker copies that insert a subject filter AFTER the
resume_groups call but BEFORE the processing loop, so:
  - resume_groups still decides what needs running (skip logic intact)
  - each worker only processes its assigned slice of the remaining subjects
  - N workers run in parallel → full machine usage

The injected filter is a single line added right before
'for list_to_keep, Skip_step in groups:' :

    groups = _ednix_filter_groups(groups, WORKER_SUBJECTS)

where WORKER_SUBJECTS is the slice of subject IDs for that worker.

VALIDATION CONTROLS (run before launching)
------------------------------------------
  1. Launcher file exists and is readable
  2. BIDS dir exists and contains sub-* folders
  3. Python venv interpreter exists
  4. The launcher contains a resume_groups call and a processing loop
  5. The generated worker is valid Python (ast.parse)
  6. The injected filter helper is present in the worker
  7. Subject slices are non-overlapping and cover all subjects exactly once
  8. Dry-run prints the full plan for inspection before any launch

USAGE
-----
From terminal (recommended):
    nohup python ednix_run.py > run.log 2>&1 &
    tail -f run.log

From PyCharm console:
    from ednix_run import run, validate, plan
    validate()           # check everything is OK, launch nothing
    plan()               # show worker allocation, launch nothing
    run()                # full run: anat then func
    run(phase='anat')    # anat only
    run(only=['Rat'])    # subset
    run(dry_run=True)    # validate + plan + show commands, no launch
"""

import ast
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

# ===========================================================================
# CONFIG — edit paths if needed
# ===========================================================================
MAIN_PATH     = Path('/home/cgarin/PycharmProjects/EDNiX')
PYTHON        = str(MAIN_PATH / 'venvEDNiX/bin/python')
LAUNCHER_ROOT = MAIN_PATH / 'Exemples/Study_EDNiX'
REPORT_DIR    = MAIN_PATH / 'run_reports'
WORKER_ROOT   = MAIN_PATH / '_ednix_workers'   # temp worker scripts live here

# Machine resources
CPU_TOTAL        = 60
RAM_TOTAL_GB     = 206
RAM_RESERVED_GB  = 16    # keep free for OS / browser / PyCharm
CPU_RESERVED     = 4

# Per-worker resource estimates (used only for allocation maths)
WORKER_RAM_GB    = 7
WORKER_CORES     = 3

# Global concurrency budget. The orchestrator keeps the SUM of cores used by
# all running workers <= CORE_BUDGET. Set as (all cores - margin) so you keep
# a few cores free for navigating the machine.
CORE_MARGIN      = 5
def get_core_budget():
    try:
        return max(1, os.cpu_count() - CORE_MARGIN)
    except Exception:
        return CPU_TOTAL - CORE_MARGIN

# per-dataset core cost of one worker (override default if needed)
WORKER_CORES_BY_KEY = {}   # e.g. {'func_Human_ds004513': 8}
def worker_cores(ds):
    return WORKER_CORES_BY_KEY.get(ds['key'], WORKER_CORES)

# ===========================================================================
# DATASET REGISTRY
#   key       : unique label
#   launcher  : path relative to LAUNCHER_ROOT
#   bids      : BIDS directory (for subject counting)
#   split     : max subject-workers for this dataset (1 = no split)
# ===========================================================================
ANAT = [
    dict(key='anat_Rat',             launcher='Launcher_anat/Launcher_anat_Rat.py',             bids='/scratch2/EDNiX/Rat/BIDS_Grandjean/',        split=6),
    dict(key='anat_Mouse',           launcher='Launcher_anat/Launcher_anat_Mouse.py',           bids='/scratch2/EDNiX/Mouse/BIDS_Grandjean2/',     split=3),
    dict(key='anat_Marmoset',        launcher='Launcher_anat/Launcher_anat_marmoset_MBM.py',    bids='/scratch2/EDNiX/Marmoset/BIDS_Tian/',        split=2),
    dict(key='anat_Dog',             launcher='Launcher_anat/Launcher_anat_Dog.py',             bids='/scratch2/EDNiX/Dog/BIDS_Boch_K9/',          split=2),
    dict(key='anat_Lemur',           launcher='Launcher_anat/Launcher_anat_Lemur.py',           bids='/scratch2/EDNiX/Mouselemur/BIDS_Garin/',     split=2),
    dict(key='anat_Human_CERMEP',    launcher='Launcher_anat/Launcher_anat_Human_CERMEP.py',    bids='/scratch2/EDNiX/Human/BIDS_Merida/',         split=2),
    dict(key='anat_Human_ds004513',  launcher='Launcher_anat/Launcher_anat_Human_ds004513.py',  bids='/scratch2/EDNiX/Human/BIDS_Castrillon/',     split=2),
    dict(key='anat_Human_ds004856',  launcher='Launcher_anat/Launcher_anat_Human_ds004856.py',  bids='/scratch2/EDNiX/Human/BIDS_Park/',           split=1),
    dict(key='anat_Macaque_Zhu',     launcher='Launcher_anat/Launcher_anat_Macaque_Zhu.py',     bids='/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin/',    split=1),
    dict(key='anat_Macaque_Hamed',   launcher='Launcher_anat/Launcher_anat_Macaque_Hamed.py',   bids='/scratch2/EDNiX/Macaque/BIDS_BenHamed/',     split=1),
    dict(key='anat_Macaque_imagina', launcher='Launcher_anat/Launcher_anat_Macaque_imagina.py', bids='/scratch2/EDNiX/Macaque/BIDS_Tremblay',      split=1),
    dict(key='anat_Bat',             launcher='Launcher_anat/Launcher_anat_Bat.py',             bids='/scratch2/EDNiX/Bat/BIDS_Washington/',       split=1),
]

FUNC = [
    dict(key='func_Rat',             launcher='Launcher_func/Launcher_func_Rat.py',             bids='/scratch2/EDNiX/Rat/BIDS_Grandjean/',        split=6),
    dict(key='func_Mouse',           launcher='Launcher_func/Launcher_func_Mouse.py',           bids='/scratch2/EDNiX/Mouse/BIDS_Grandjean2/',     split=3),
    dict(key='func_Marmoset',        launcher='Launcher_func/Launcher_func_marmoset_MBM.py',    bids='/scratch2/EDNiX/Marmoset/BIDS_Tian/',        split=2),
    dict(key='func_Dog',             launcher='Launcher_func/Launcher_func_Dog.py',             bids='/scratch2/EDNiX/Dog/BIDS_Boch_K9/',          split=2),
    dict(key='func_Lemur',           launcher='Launcher_func/Launcher_func_Lemur.py',           bids='/scratch2/EDNiX/Mouselemur/BIDS_Garin/',     split=2),
    dict(key='func_Human_ds004513',  launcher='Launcher_func/Launcher_func_Human_ds004513.py',  bids='/scratch2/EDNiX/Human/BIDS_Castrillon/',     split=2),
    dict(key='func_Human_ds004856',  launcher='Launcher_func/Launcher_func_Human_ds004856.py',  bids='/scratch2/EDNiX/Human/BIDS_Park/',           split=1),
    dict(key='func_Macaque_Zhu',     launcher='Launcher_func/Launcher_func_Macaque_Zhu.py',     bids='/scratch2/EDNiX/Macaque/BIDS_Zhu_Garin/',    split=1),
    dict(key='func_Macaque_Hamed',   launcher='Launcher_func/Launcher_func_Macaque_Hamed.py',   bids='/scratch2/EDNiX/Macaque/BIDS_BenHamed/',     split=1),
    dict(key='func_Macaque_imagina', launcher='Launcher_func/Launcher_func_Macaque_imagina.py', bids='/scratch2/EDNiX/Macaque/BIDS_Tremblay/',     split=1),
    dict(key='func_Bat',             launcher='Launcher_func/Launcher_func_Bat.py',             bids='/scratch2/EDNiX/Bat/BIDS_Washington/',       split=1),
]

# ===========================================================================
# The filter helper injected into every worker.
# It keeps only the (ID, Session) groups whose ID is in the worker's slice.
# ===========================================================================
FILTER_HELPER = '''
# ===== injected by ednix_run.py =====
import os as _os
def _ednix_filter_groups(groups, keep_ids):
    """Keep only subjects whose ID is in keep_ids (set of strings)."""
    keep = set(str(x) for x in keep_ids)
    out = []
    for list_to_keep, Skip_step in groups:
        filtered = [(i, s) for (i, s) in list_to_keep if str(i) in keep]
        if filtered:
            out.append((filtered, Skip_step))
    return out
_EDNIX_WORKER_SUBJECTS = {worker_subjects!r}
# =====================================
'''

# ===========================================================================
# Logging
# ===========================================================================
_lock    = threading.Lock()
_run_log = []

def _log(msg, level='INFO', key=''):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    col = {'INFO':'', 'OK':'\033[92m', 'WARN':'\033[93m',
           'ERROR':'\033[91m', 'HEAD':'\033[96m', 'DIM':'\033[2m'}
    with _lock:
        print(f"[{ts}] [{level:5s}] {col.get(level,'')}{msg}\033[0m", flush=True)
        _run_log.append((ts, key, level, msg))

# ===========================================================================
# Helpers
# ===========================================================================
def get_subjects(bids_dir):
    p = Path(bids_dir)
    if not p.exists():
        return []
    return sorted(d.name.replace('sub-', '')
                  for d in p.iterdir()
                  if d.is_dir() and d.name.startswith('sub-'))

def get_free_ram_gb():
    try:
        mem = {}
        with open('/proc/meminfo') as f:
            for line in f:
                k, v = line.split(':')
                mem[k.strip()] = int(v.split()[0])
        return mem.get('MemAvailable', mem['MemFree']) / 1024 / 1024
    except Exception:
        return RAM_TOTAL_GB - RAM_RESERVED_GB

def get_load():
    try:
        with open('/proc/loadavg') as f:
            return float(f.read().split()[0])
    except Exception:
        return 0.0

# ===========================================================================
# VALIDATION
# ===========================================================================
def validate(verbose=True):
    """
    Run all pre-flight checks. Returns (ok:bool, problems:list[str]).
    Launches nothing.
    """
    problems = []

    # interpreter
    if not Path(PYTHON).exists():
        problems.append(f"Python interpreter not found: {PYTHON}")

    for ds in ANAT + FUNC:
        key = ds['key']
        lp  = LAUNCHER_ROOT / ds['launcher']

        # launcher exists
        if not lp.exists():
            problems.append(f"[{key}] launcher missing: {lp}")
            continue

        src = lp.read_text()

        # has resume_groups + processing loop
        has_rg   = ('fmri_resume_groups' in src or 'anat_resume_groups' in src)
        has_loop = bool(re.search(
            r"for\s+list_to_keep\s*,\s*Skip_step\s+in\s+groups\s*:", src))
        if not has_rg:
            problems.append(f"[{key}] no resume_groups call found")
        if not has_loop:
            problems.append(f"[{key}] no 'for list_to_keep, Skip_step in groups:' loop")

        # bids dir
        if not Path(ds['bids']).exists():
            problems.append(f"[{key}] BIDS dir missing: {ds['bids']}")
        else:
            subs = get_subjects(ds['bids'])
            if not subs:
                problems.append(f"[{key}] no sub-* folders in {ds['bids']}")

    if verbose:
        if problems:
            _log(f"VALIDATION FAILED — {len(problems)} problem(s):", 'ERROR')
            for p in problems:
                _log(f"  ✗ {p}", 'ERROR')
        else:
            _log("VALIDATION PASSED — all launchers OK", 'OK')

    return (len(problems) == 0), problems

# ===========================================================================
# WORKER GENERATION
# ===========================================================================
def make_workers(ds, dry_run=False):
    """
    Generate worker scripts for one dataset.
    Returns list of (worker_path, subject_slice) — or [(launcher, None)]
    if split==1.
    Validates each generated worker with ast.parse.
    """
    lp  = LAUNCHER_ROOT / ds['launcher']
    src = lp.read_text()
    subs = get_subjects(ds['bids'])
    n_split = min(ds['split'], len(subs)) if subs else 1

    # no split → run the launcher as-is (resume_groups handles skip)
    if n_split <= 1:
        return [(lp, None)]

    # find the processing loop
    loop_match = re.search(
        r"^(\s*)for\s+list_to_keep\s*,\s*Skip_step\s+in\s+groups\s*:",
        src, re.MULTILINE)
    if not loop_match:
        # can't split safely — fall back to whole launcher
        _log(f"[{ds['key']}] cannot locate groups loop — running unsplit",
             'WARN', ds['key'])
        return [(lp, None)]

    indent     = loop_match.group(1)
    loop_start = loop_match.start()

    # split subjects round-robin so slow subjects spread evenly across workers
    slices = [subs[i::n_split] for i in range(n_split)]
    slices = [sl for sl in slices if sl]

    # VALIDATION 7: non-overlapping + complete coverage
    flat = [s for sl in slices for s in sl]
    assert sorted(flat) == sorted(subs), \
        f"[{ds['key']}] subject slicing lost/duplicated subjects!"
    assert len(flat) == len(set(flat)), \
        f"[{ds['key']}] subject slices overlap!"

    WORKER_ROOT.mkdir(parents=True, exist_ok=True)
    outdir = WORKER_ROOT / ds['key']
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    workers = []
    for i, sl in enumerate(slices):
        # insert filter helper before the loop, and a filter call right after
        helper = FILTER_HELPER.format(worker_subjects=sl)
        filter_call = (f"\n{indent}groups = _ednix_filter_groups("
                       f"groups, _EDNIX_WORKER_SUBJECTS)\n")

        # inject sys.path BEFORE any import so 'modalities' is findable
        path_fix = f"import sys; sys.path.insert(0, '{MAIN_PATH}')\n"
        new_src = (
            path_fix
            + src[:loop_start]
            + helper
            + filter_call
            + src[loop_start:]
        )

        # VALIDATION 5: generated worker is valid Python
        try:
            ast.parse(new_src)
        except SyntaxError as e:
            _log(f"[{ds['key']}] worker {i+1} invalid Python: {e} — "
                 f"running unsplit", 'ERROR', ds['key'])
            return [(lp, None)]

        # VALIDATION 6: helper present
        if '_ednix_filter_groups' not in new_src:
            _log(f"[{ds['key']}] filter helper missing — running unsplit",
                 'ERROR', ds['key'])
            return [(lp, None)]

        wpath = outdir / f"worker_{i+1:02d}of{len(slices):02d}.py"
        if not dry_run:
            wpath.write_text(new_src)
        workers.append((wpath, sl))

    return workers

# ===========================================================================
# LAUNCH
# ===========================================================================
def launch_dataset(ds, results, dry_run=False):
    """Thread target: generate workers and run them in parallel."""
    key = ds['key']
    subs = get_subjects(ds['bids'])
    if not subs:
        _log(f"[{key}] no subjects — skipping", 'WARN', key)
        results[key] = {'ok': True, 'errors': [], 'note': 'no subjects'}
        return

    workers = make_workers(ds, dry_run=dry_run)
    _log(f"[{key}] {len(subs)} subjects → {len(workers)} worker(s)", 'HEAD', key)

    if dry_run:
        for wpath, sl in workers:
            n = len(sl) if sl else len(subs)
            _log(f"    [DRY] {Path(wpath).name}  ({n} subjects)", 'DIM', key)
        results[key] = {'ok': True, 'errors': [], 'note': 'dry-run'}
        return

    procs  = []
    errors = []
    for wpath, sl in workers:
        env = os.environ.copy()
        for v in ('OMP_NUM_THREADS', 'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',
                  'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            env[v] = str(WORKER_CORES)
        proc = subprocess.Popen(
            [PYTHON, str(wpath)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, env=env, start_new_session=True)
        procs.append({'name': Path(wpath).name, 'proc': proc,
                      'pid': proc.pid, 'start': time.time(), 'tail': []})
        _log(f"    PID {proc.pid:>7}  {Path(wpath).name}", 'INFO', key)

    while procs:
        still = []
        for p in procs:
            ret = p['proc'].poll()
            if ret is None:
                try:
                    line = p['proc'].stdout.readline()
                    if line.strip():
                        p['tail'].append(line.rstrip())
                        p['tail'] = p['tail'][-30:]
                        with _lock:
                            print(f"      [{key}/{p['name']}] {line.rstrip()}",
                                  flush=True)
                except Exception:
                    pass
                still.append(p)
            else:
                el = (time.time() - p['start']) / 60
                try:
                    rest = p['proc'].stdout.read()
                    if rest:
                        p['tail'].extend(rest.splitlines())
                        p['tail'] = p['tail'][-30:]
                except Exception:
                    pass
                if ret == 0:
                    _log(f"    Finished OK: {p['name']} ({el:.0f} min)", 'OK', key)
                else:
                    m = f"FAILED (exit {ret}): {p['name']} ({el:.0f} min)"
                    _log(f"    {m}", 'ERROR', key)
                    errors.append(m)
                    errors.append("      last output:\n" +
                                  '\n'.join(f"        {l}" for l in p['tail']))
        procs = still
        if procs:
            time.sleep(2)

    results[key] = {'ok': len(errors) == 0, 'errors': errors, 'note': ''}

# ===========================================================================
# PHASE RUNNER
# ===========================================================================
def run_phase(datasets, name, only=None, dry_run=False):
    """
    Core-budgeted scheduler. Builds a flat queue of ALL workers across all
    datasets, then runs them so the sum of cores in flight stays <= budget.
    Finished workers free their cores immediately and the next queued worker
    starts — so already-done subjects (skipped fast by resume_groups) never
    block the machine, and the machine stays full without ever exceeding it.
    """
    if only:
        datasets = [d for d in datasets
                    if any(o.lower() in d['key'].lower() for o in only)]
    if not datasets:
        _log(f"No datasets for phase {name}", 'WARN')
        return {}

    budget = get_core_budget()
    _log(f"\n{'='*68}", 'HEAD')
    _log(f"PHASE {name}: {len(datasets)} datasets | core budget = {budget} "
         f"(cpu {os.cpu_count()} - margin {CORE_MARGIN})", 'HEAD')
    _log(f"{'='*68}", 'HEAD')

    # Build flat worker queue: list of (ds, worker_path, subject_slice, cores)
    queue = []
    for d in datasets:
        subs = get_subjects(d['bids'])
        if not subs:
            _log(f"[{d['key']}] no subjects — skipping", 'WARN', d['key'])
            continue
        workers = make_workers(d, dry_run=dry_run)
        c = worker_cores(d)
        for wpath, sl in workers:
            queue.append({'ds': d, 'path': wpath, 'slice': sl, 'cores': c})

    _log(f"Queued {len(queue)} workers total", 'INFO')

    if dry_run:
        for w in queue:
            n = len(w['slice']) if w['slice'] else len(get_subjects(w['ds']['bids']))
            _log(f"    [DRY] {w['ds']['key']}/{Path(w['path']).name} "
                 f"({n} subj, {w['cores']} cores)", 'DIM')
        return {d['key']: {'ok': True, 'errors': [], 'note': 'dry-run'}
                for d in datasets}

    results = {d['key']: {'ok': True, 'errors': [], 'note': ''}
               for d in datasets}
    running = []          # list of dicts with proc + meta
    cores_in_flight = 0
    qi = 0

    def _start(w):
        nonlocal cores_in_flight
        env = os.environ.copy()
        for v in ('OMP_NUM_THREADS', 'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS',
                  'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            env[v] = str(w['cores'])
        proc = subprocess.Popen(
            [PYTHON, str(w['path'])],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, env=env, start_new_session=True)
        cores_in_flight += w['cores']
        _log(f"  ▶ START {w['ds']['key']}/{Path(w['path']).name}  "
             f"PID {proc.pid}  ({w['cores']} cores, "
             f"{cores_in_flight}/{budget} in flight)", 'INFO', w['ds']['key'])
        return {'proc': proc, 'w': w, 'start': time.time(), 'tail': []}

    # main scheduling loop
    while qi < len(queue) or running:
        # fill up to budget
        while qi < len(queue) and cores_in_flight + queue[qi]['cores'] <= budget:
            running.append(_start(queue[qi]))
            qi += 1
        # if nothing could start and nothing running (single worker > budget)
        if not running and qi < len(queue):
            # force-start one oversized worker
            running.append(_start(queue[qi]))
            qi += 1

        # poll running
        still = []
        for r in running:
            ret = r['proc'].poll()
            if ret is None:
                try:
                    line = r['proc'].stdout.readline()
                    if line.strip():
                        r['tail'].append(line.rstrip()); r['tail'] = r['tail'][-30:]
                        with _lock:
                            print(f"      [{r['w']['ds']['key']}/"
                                  f"{Path(r['w']['path']).name}] {line.rstrip()}",
                                  flush=True)
                except Exception:
                    pass
                still.append(r)
            else:
                cores_in_flight -= r['w']['cores']
                el = (time.time() - r['start']) / 60
                key = r['w']['ds']['key']
                try:
                    rest = r['proc'].stdout.read()
                    if rest:
                        r['tail'].extend(rest.splitlines()); r['tail'] = r['tail'][-30:]
                except Exception:
                    pass
                if ret == 0:
                    _log(f"  ✓ DONE  {key}/{Path(r['w']['path']).name} "
                         f"({el:.0f} min, freed {r['w']['cores']} cores)", 'OK', key)
                else:
                    m = f"FAILED (exit {ret}): {Path(r['w']['path']).name} ({el:.0f} min)"
                    _log(f"  ✗ {key}/{m}", 'ERROR', key)
                    results[key]['ok'] = False
                    results[key]['errors'].append(m)
                    results[key]['errors'].append(
                        "      last output:\n" +
                        '\n'.join(f"        {l}" for l in r['tail']))
        running = still
        if running or qi < len(queue):
            time.sleep(2)

    failed = [k for k, v in results.items() if not v['ok']]
    _log(f"\nPhase {name} done — {len(results)-len(failed)} OK / "
         f"{len(failed)} FAILED", 'HEAD')
    return results

# ===========================================================================
# PLAN (allocation preview)
# ===========================================================================
def plan(only=None):
    """Print the worker-allocation plan. Launches nothing."""
    free = get_free_ram_gb() - RAM_RESERVED_GB
    budget = get_core_budget()
    _log("ALLOCATION PLAN", 'HEAD')
    _log(f"Free RAM ~{free:.0f} GB | CPU {os.cpu_count()} | "
         f"core budget {budget} (margin {CORE_MARGIN}) | Load {get_load():.1f}",
         'INFO')
    _log("Workers run from a global queue capped at the core budget — "
         "totals below are QUEUE sizes, not simultaneous load.", 'DIM')
    for name, dss in [('ANAT', ANAT), ('FUNC', FUNC)]:
        sel = dss if not only else [d for d in dss
                                    if any(o.lower() in d['key'].lower() for o in only)]
        total_w = sum(min(d['split'], len(get_subjects(d['bids'])) or 1)
                      for d in sel)
        _log(f"\n[{name}] {len(sel)} datasets, ~{total_w} workers peak, "
             f"~{total_w*WORKER_RAM_GB} GB / ~{total_w*WORKER_CORES} cores", 'INFO')
        for d in sel:
            subs = get_subjects(d['bids'])
            nw   = min(d['split'], len(subs)) if subs else 1
            _log(f"  {d['key']:<26} {len(subs):>4} subjects → {nw} worker(s)", 'INFO')

# ===========================================================================
# REPORT
# ===========================================================================
def write_report(all_results, t0, run_id):
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    path = REPORT_DIR / f"run_report_{run_id}.txt"
    el = (time.time() - t0) / 3600
    L = ["="*72,
         "EDNiX PARALLEL RUN REPORT",
         f"Run ID  : {run_id}",
         f"End     : {datetime.now():%Y-%m-%d %H:%M:%S}",
         f"Duration: {el:.2f} h",
         "="*72, ""]

    L += ["DATASET RESULTS", "-"*72]
    for phase, res in all_results.items():
        L.append(f"\n  [{phase.upper()}]")
        for key, r in res.items():
            mark = '✓' if r['ok'] else '✗'
            note = f"  ({r['note']})" if r.get('note') else ''
            L.append(f"    {mark} {key:<28} {'OK' if r['ok'] else 'FAILED'}{note}")

    errs = [(p, k, r['errors']) for p, res in all_results.items()
            for k, r in res.items() if r['errors']]
    if errs:
        L += ["", "", "ERRORS", "-"*72]
        for p, k, e in errs:
            L.append(f"\n  [{p.upper()}] {k}:")
            L += [f"    {x}" for x in e]
    else:
        L += ["", "", "NO ERRORS DETECTED", "-"*72]

    # MAJOR_WARNING scan
    L += ["", "", "MAJOR_WARNING FILES", "-"*72]
    found = []
    for d in ANAT + FUNC:
        bp = Path(d['bids'])
        if bp.exists():
            for wf in bp.rglob('MAJOR_WARNING*.txt'):
                try:
                    found.append(f"  {wf}\n      {wf.read_text().strip()[:150]}")
                except Exception:
                    found.append(f"  {wf}")
    L += found if found else ["  None found."]

    L += ["", "", "FULL LOG", "-"*72]
    L += [f"  [{ts}] [{lv}] [{k}] {m}" for ts, k, lv, m in _run_log]

    path.write_text('\n'.join(L))
    _log(f"Report → {path}", 'OK')
    return path

# ===========================================================================
# MAIN
# ===========================================================================
def run(phase='all', only=None, dry_run=False, skip_validation=False):
    """
    Run the full pipeline.

    phase   : 'all' | 'anat' | 'func'
    only    : list of key fragments, e.g. ['Rat','Mouse']
    dry_run : validate + plan + show workers, launch nothing
    """
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S')
    t0 = time.time()

    _log("EDNiX PARALLEL ORCHESTRATOR", 'HEAD')
    _log(f"Run ID {run_id} | phase={phase} | only={only} | dry_run={dry_run}",
         'INFO')

    # --- pre-flight validation -------------------------------------------
    ok, problems = validate(verbose=True)
    if not ok and not skip_validation:
        _log("Aborting due to validation errors. "
             "Fix them or call run(skip_validation=True) to force.", 'ERROR')
        return {}

    load = get_load()
    if load > 40 and not dry_run:
        _log(f"WARNING: load is high ({load:.1f}). Consider waiting.", 'WARN')

    plan(only=only)

    all_results = {}
    try:
        if phase in ('all', 'anat'):
            all_results['anat'] = run_phase(ANAT, 'ANAT', only, dry_run)
            if not dry_run:
                _log("ANAT phase complete → starting FUNC", 'OK')
        if phase in ('all', 'func'):
            all_results['func'] = run_phase(FUNC, 'FUNC', only, dry_run)
    except KeyboardInterrupt:
        _log("Interrupted — writing partial report", 'WARN')
    finally:
        if not dry_run:
            write_report(all_results, t0, run_id)
            # cleanup temp workers
            if WORKER_ROOT.exists():
                shutil.rmtree(WORKER_ROOT, ignore_errors=True)
            _log(f"Done in {(time.time()-t0)/3600:.2f} h", 'OK')

    return all_results

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--ednix-phase', default='all', choices=['all','anat','func'])
    ap.add_argument('--only', default='')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--plan', action='store_true')
    a, _ = ap.parse_known_args()
    only = [x.strip() for x in a.only.split(',') if x.strip()] or None
    if a.validate:
        validate(); return
    if a.plan:
        plan(only=only); return
    run(phase=a.ednix_phase, only=only, dry_run=a.dry_run)

if __name__ == '__main__':
    main()
