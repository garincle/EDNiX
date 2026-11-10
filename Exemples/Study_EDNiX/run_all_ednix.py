#!/usr/bin/env python3
"""
run_all_ednix.py — launcher of launchers: (re)run every EDNiX launcher of the study in one go, on all the cores of the
machine, then build the QC reports.

    python Exemples/Study_EDNiX/run_all_ednix.py --plan          # what would run, nothing is launched
    python Exemples/Study_EDNiX/run_all_ednix.py                 # run (settings below)
    python Exemples/Study_EDNiX/run_all_ednix.py --only Rat,Mouse --phases anat,func
    nohup python Exemples/Study_EDNiX/run_all_ednix.py > run_all.log 2>&1 &

What it does
  1. Finds the launchers (Launcher_anat / _func / _pet / _ct, Archives excluded) and reads, without running them,
     their bids_dir, list_to_keep, list_to_remove and whether they build a study template.
  2. RESET (optional, RESET = True): previews the reset of every selected dataset (anat / func / PET / CT cleaners,
     manual masks and crops kept) and asks ONE 'YES' for all, in this terminal, before anything is launched.
  3. Runs the launchers as a dependency graph, not as phases: anat -> CT -> PET and anat -> func for the SAME BIDS
     folder; different datasets never wait for each other (the func of the rat starts as soon as the anat of the
     rat is done, even if the human anat is still running).
  4. Splits a launcher into subject workers (each worker = a copy of the launcher with list_to_keep set to its
     subjects; all the sessions of a subject stay together) unless it builds a study template (anat / PET / CT),
     which needs all the subjects in the same process.
  5. Keeps the sum of the threads of the running workers <= CORE_BUDGET (threads per worker set through
     OMP / ITK / OpenBLAS / MKL variables), starting the longest jobs first (fewer idle cores at the end).
  6. QC: Tools/QC_report.build on every BIDS folder (subject recaps, all-subject page, dataset page + Excel) and one
     Excel file across datasets (REPORT_DIR/<run>/EDNiX_study_summary.xlsx).
  One log per worker, a report with the duration and status of every job and the tail of the failed logs.

The workers run with stdin closed and EDNIX_SKIP_RESET=1: a launcher that would wait for a 'YES' or an ITK-SNAP check
fails instead of hanging; keep 'itk_*' in Skip_step for unattended runs.
"""
import ast
import os
import re
import sys
import glob
import time
import shutil
import argparse
import datetime
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAIN_PATH = HERE.parents[1]                       # the EDNiX repository
sys.path.insert(0, str(MAIN_PATH))

# =====================================================================================================================
# SETTINGS
# =====================================================================================================================
PYTHON = sys.executable                           # interpreter of the EDNiX environment
REPORT_DIR = MAIN_PATH / 'run_reports'
PHASES = ['anat', 'ct', 'pet', 'func', 'qc']      # 'qc' = QC reports at the end
ONLY = []                                         # launcher name fragments, e.g. ['Rat', 'Mouse']; [] = all
EXCLUDE = ['quicktest', '_custom', 'Check_pet_boundary']
RESET = False                                     # preview + ONE 'YES' before launching (see the cleaners)
RESET_ANAT_SCOPE = 'all'                          # 'all' or 'surfaces' (rebuild only the surfaces)
THREADS_PER_WORKER = 4                            # ANTs / AFNI / numpy threads of one worker
CORE_MARGIN = 4                                   # cores left free for the machine
MAX_SPLIT = 6                                     # max subject workers per launcher
WEIGHT = {'anat': 3.0, 'func': 2.0, 'pet': 1.0, 'ct': 0.5}   # relative cost per subject (scheduling order only)


def core_budget():
    return max(THREADS_PER_WORKER, (os.cpu_count() or 8) - CORE_MARGIN)


# =====================================================================================================================
# launcher discovery (static: the launchers are parsed, never imported)
# =====================================================================================================================
KIND_DIRS = {'anat': 'Launcher_anat', 'func': 'Launcher_func', 'pet': 'Launcher_pet', 'ct': 'Launcher_ct'}
MODALITY = {'anat': 'anat', 'func': 'func', 'pet': 'pet', 'ct': 'ct'}


def _top_assign(tree, name):
    """Last top-level `name = ...` node."""
    found = None
    for n in tree.body:
        if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name) \
                and n.targets[0].id == name:
            found = n
    return found


def _literal(tree, name, default=None):
    n = _top_assign(tree, name)
    if n is None:
        return default
    try:
        return ast.literal_eval(n.value)
    except Exception:
        strs = [c.value for c in ast.walk(n.value) if isinstance(c, ast.Constant) and isinstance(c.value, str)]
        return strs[0] if len(strs) == 1 else default


def parse_launcher(path, kind):
    src = Path(path).read_text(encoding='utf-8')
    tree = ast.parse(src)
    bids = _literal(tree, 'bids_dir')
    keep = _literal(tree, 'list_to_keep', []) or []
    remove = _literal(tree, 'list_to_remove', []) or []
    template = any(_literal(tree, v, False) is True for v in ('creat_study_template_pet', 'creat_study_template_ct')) \
        or (kind == 'anat' and _literal(tree, 'creat_study_template', False) is True)
    anchors = [n for n in (_top_assign(tree, 'list_to_keep'), _top_assign(tree, 'list_to_remove')) if n is not None]
    return dict(path=str(path), name=Path(path).stem, kind=kind, bids=bids, keep=[tuple(map(str, x)) for x in keep],
                remove=[tuple(map(str, x)) for x in remove], template=template,
                insert_line=max(n.end_lineno for n in anchors) if len(anchors) == 2 else None, src=src)


def discover(only=None, exclude=EXCLUDE):
    out = []
    for kind, d in KIND_DIRS.items():
        for p in sorted((HERE / d).glob('*.py')):
            if any(x in p.name for x in exclude):
                continue
            if only and not any(o.lower() in p.stem.lower() for o in only):
                continue
            try:
                info = parse_launcher(p, kind)
            except SyntaxError as e:
                print(f'SKIP {p.name}: {e}')
                continue
            if not info['bids']:
                print(f'SKIP {p.name}: no bids_dir found')
                continue
            out.append(info)
    return out


def subject_sessions(bids, modality):
    """{ID: [sessions]} of the BIDS folder that have a <modality> folder (raw data)."""
    out = {}
    for d in sorted(glob.glob(os.path.join(bids, 'sub-*', 'ses-*', modality))):
        ses = os.path.basename(os.path.dirname(d))[4:]
        sub = os.path.basename(os.path.dirname(os.path.dirname(d)))[4:]
        out.setdefault(sub, []).append(ses)
    return out


def selected_pairs(L):
    pairs = [(s, ses) for s, ss in subject_sessions(L['bids'], MODALITY[L['kind']]).items() for ses in ss]
    if L['keep']:
        pairs = [p for p in pairs if p in set(L['keep'])]
    return [p for p in pairs if p not in set(L['remove'])]


# =====================================================================================================================
# workers
# =====================================================================================================================
HELPER = '''
# ===== injected by run_all_ednix.py: this worker processes only its subjects =====
_EDNIX_WORKER_PAIRS = {pairs!r}
list_to_keep = [p for p in _EDNIX_WORKER_PAIRS
                if (not list_to_keep or p in {{(str(a), str(b)) for a, b in list_to_keep}})
                and p not in {{(str(a), str(b)) for a, b in list_to_remove}}]
if not list_to_keep:
    print('run_all_ednix worker: no subject left after list_to_keep / list_to_remove')
    sys.exit(0)
# ================================================================================
'''


def make_jobs(launchers, work_dir):
    jobs = []
    for L in launchers:
        pairs = selected_pairs(L)
        subs = sorted({s for s, _ in pairs})
        n = 1 if (L['template'] or L['insert_line'] is None) else max(1, min(MAX_SPLIT, len(subs)))
        cost = WEIGHT[L['kind']] * max(1, len(pairs))
        if n == 1:
            jobs.append(dict(launcher=L, script=L['path'], pairs=pairs, cost=cost, label=L['name']))
            continue
        lines = L['src'].split('\n')
        for i in range(n):
            mine = set(subs[i::n])
            wp = [p for p in pairs if p[0] in mine]
            src = '\n'.join(lines[:L['insert_line']]) + '\nimport sys\n' + HELPER.format(pairs=wp) + \
                  '\n'.join(lines[L['insert_line']:])
            ast.parse(src)                                   # never launch an invalid worker
            wpath = Path(work_dir) / f"{L['name']}__w{i + 1:02d}of{n:02d}.py"
            wpath.write_text(f'import sys; sys.path.insert(0, {str(MAIN_PATH)!r})\n' + src)
            jobs.append(dict(launcher=L, script=str(wpath), pairs=wp, cost=WEIGHT[L['kind']] * len(wp),
                             label=wpath.stem))
    return jobs


def dependencies(jobs):
    """job index -> set of job indices that must finish first (same BIDS folder: anat -> ct -> pet, anat -> func)."""
    before = {'anat': (), 'ct': ('anat',), 'pet': ('anat', 'ct'), 'func': ('anat',)}
    norm = lambda b: os.path.normpath(b)
    deps = {}
    for i, j in enumerate(jobs):
        L = j['launcher']
        deps[i] = {k for k, o in enumerate(jobs) if norm(o['launcher']['bids']) == norm(L['bids'])
                   and o['launcher']['kind'] in before[L['kind']]}
    return deps


def plan_text(jobs, deps):
    lines = [f'{len(jobs)} job(s), core budget {core_budget()} ({THREADS_PER_WORKER} threads / worker)']
    for i, j in enumerate(jobs):
        L = j['launcher']
        after = ', '.join(sorted({jobs[k]['launcher']['name'] for k in deps[i]})) or '-'
        lines.append(f"  [{i:02d}] {j['label']:<52} {L['kind']:<4} {len(j['pairs']):>4} session(s)  "
                     f"{'(study template: not split) ' if L['template'] else ''}after: {after}")
    return '\n'.join(lines)


def run_graph(jobs, deps, log_dir):
    budget = core_budget()
    env0 = os.environ.copy()
    for v in ('OMP_NUM_THREADS', 'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'ANTS_RANDOM_SEED'):
        env0.pop(v, None)
    env0.update(EDNIX_SKIP_RESET='1', MPLBACKEND='Agg', PYTHONUNBUFFERED='1',
                PYTHONPATH=os.pathsep.join([str(MAIN_PATH)] + ([env0['PYTHONPATH']] if env0.get('PYTHONPATH') else [])))
    state = {i: 'waiting' for i in range(len(jobs))}
    running, results = {}, {}
    order = sorted(range(len(jobs)), key=lambda i: -jobs[i]['cost'])       # longest first
    while any(s in ('waiting', 'running') for s in state.values()):
        # a job whose prerequisite failed is not started
        for i in order:
            if state[i] == 'waiting' and any(state[k] in ('failed', 'blocked') for k in deps[i]):
                state[i] = 'blocked'
                results[i] = dict(status='blocked', minutes=0, log='')
                print(f"  ⨯ BLOCK {jobs[i]['label']} (a prerequisite failed)", flush=True)
        used = len(running) * THREADS_PER_WORKER
        for i in order:
            if state[i] != 'waiting' or any(state[k] != 'done' for k in deps[i]):
                continue
            if used + THREADS_PER_WORKER > budget and running:
                break
            env = dict(env0)
            for v in ('OMP_NUM_THREADS', 'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS', 'OPENBLAS_NUM_THREADS',
                      'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
                env[v] = str(THREADS_PER_WORKER)
            log = Path(log_dir) / f"{jobs[i]['label']}.log"
            fh = open(log, 'w')
            proc = subprocess.Popen([PYTHON, jobs[i]['script']], stdout=fh, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, env=env, cwd=str(MAIN_PATH), start_new_session=True)
            running[i] = (proc, fh, time.time(), log)
            state[i] = 'running'
            used += THREADS_PER_WORKER
            print(f"  ▶ START {jobs[i]['label']}  pid {proc.pid}  ({used}/{budget} threads)", flush=True)
        for i, (proc, fh, t0, log) in list(running.items()):
            ret = proc.poll()
            if ret is None:
                continue
            fh.close()
            minutes = (time.time() - t0) / 60
            state[i] = 'done' if ret == 0 else 'failed'
            results[i] = dict(status=state[i], minutes=minutes, log=str(log), code=ret)
            print(f"  {'✓ DONE ' if ret == 0 else '✗ FAIL '} {jobs[i]['label']} ({minutes:.0f} min)", flush=True)
            del running[i]
        if running or any(s == 'waiting' for s in state.values()):
            time.sleep(5)
    return results


# =====================================================================================================================
# reset and QC
# =====================================================================================================================
def reset_all(launchers):
    """Preview the reset of every selected launcher, then ONE 'YES' for everything."""
    from Tools import data_reset as dr
    from modalities.anat import _000_Data_Reset as anat_r
    from modalities.pet import _000_Data_Reset as pet_r
    from modalities.fMRI import _000_Data_Reset as func_r
    plans, func_paths = [], []
    for L in launchers:
        for sub, ses in selected_pairs(L):
            dp = os.path.join(L['bids'], f'sub-{sub}', f'ses-{ses}')
            diary = os.path.join(dp, f"{sub}_session_{ses}_reset_{L['kind']}_diary.txt")
            if L['kind'] == 'anat':
                root, td, kept = anat_r.plan_anat(dp, RESET_ANAT_SCOPE)
            elif L['kind'] in ('pet', 'ct'):
                root, td, kept = pet_r.plan_emission(dp, L['kind'])
            else:
                func_paths.append((dp, diary))
                continue
            plans.append(dict(label=f"{L['name']}: sub-{sub} ses-{ses}", root=root, to_delete=td, kept=kept,
                              diary_file=diary))
    for pl in plans:
        dr.execute(pl['label'], pl['root'], pl['to_delete'], pl['kept'], pl['diary_file'], dry_run=True)
    for dp, diary in func_paths:
        func_r.reset_func(dp, diary, dry_run=True)
    n = sum(len(p['to_delete']) for p in plans)
    ans = input(f"\n🚨 Reset {len(plans)} anat/PET/CT session(s) ({n} files) and {len(func_paths)} func session(s) "
                f"listed above? Type 'YES': ")
    if ans.strip() != 'YES':
        print('Reset cancelled: nothing deleted, nothing launched.')
        return False
    for pl in plans:
        dr.execute(pl['label'], pl['root'], pl['to_delete'], pl['kept'], pl['diary_file'], dry_run=False,
                   skip_confirmation=True, verbose=False)
    for dp, diary in func_paths:
        func_r.reset_func(dp, diary, dry_run=False, skip_confirmation=True)
    return True


def build_qc(launchers, out_dir):
    from Tools import QC_report
    bids = sorted({os.path.normpath(L['bids']) for L in launchers if os.path.isdir(L['bids'])})
    for b in bids:
        try:
            QC_report.build(b)
        except Exception as e:
            print(f'QC report failed for {b}: {e}')
    if bids:
        print('Study summary:', QC_report.study_summary(bids, os.path.join(out_dir, 'EDNiX_study_summary.xlsx')))


def write_report(jobs, results, out_dir, t0):
    lines = [f'EDNiX run {datetime.datetime.now():%Y-%m-%d %H:%M}, {(time.time() - t0) / 3600:.2f} h', '']
    for i, j in enumerate(jobs):
        r = results.get(i, dict(status='not run', minutes=0, log=''))
        lines.append(f"{r['status']:<8} {r['minutes']:>7.1f} min  {j['label']}  {r['log']}")
    for i, r in results.items():
        if r['status'] == 'failed':
            tail = Path(r['log']).read_text(errors='replace').splitlines()[-25:]
            lines += ['', f"=== {jobs[i]['label']} (exit {r.get('code')})"] + ['    ' + l for l in tail]
    p = Path(out_dir) / 'report.txt'
    p.write_text('\n'.join(lines))
    print('\n'.join(lines[:len(jobs) + 2]))
    print(f'Report: {p}')
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--only', default=','.join(ONLY), help='launcher name fragments, comma-separated')
    ap.add_argument('--phases', default=','.join(PHASES), help='anat,ct,pet,func,qc')
    ap.add_argument('--plan', action='store_true', help='show the jobs and their order, launch nothing')
    ap.add_argument('--reset', action='store_true', default=RESET, help='reset (with ONE confirmation) first')
    a = ap.parse_args()
    only = [x for x in a.only.split(',') if x]
    phases = [x for x in a.phases.split(',') if x]
    launchers = [L for L in discover(only) if L['kind'] in phases]
    run_id = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir = REPORT_DIR / run_id
    work = out_dir / 'workers'
    work.mkdir(parents=True, exist_ok=True)
    jobs = make_jobs(launchers, work)
    deps = dependencies(jobs)
    print(plan_text(jobs, deps))
    if a.plan:
        shutil.rmtree(out_dir, ignore_errors=True)
        return
    if a.reset and not reset_all(launchers):
        return
    t0 = time.time()
    results = run_graph(jobs, deps, out_dir) if jobs else {}
    write_report(jobs, results, out_dir, t0)
    if 'qc' in phases:
        build_qc(discover(only), out_dir)


if __name__ == '__main__':
    main()
