#!/usr/bin/env python3
"""
End-to-end smoke test of EDNiX on real data, reduced for speed.

For every example launcher whose BIDS folder exists on this machine:
  1. copy ONE subject/session of the dataset into a work folder, keeping only the
     first --volumes volumes of each functional run;
  2. write a copy of the launcher pointing to that reduced dataset, restricted to that
     subject, skipping the interactive ITK-SNAP steps, MELODIC and the FreeSurfer /
     surface steps;
  3. run it with EDNIX_STRICT=1 (the first failing command stops the run).
Launchers sharing a BIDS folder use the same reduced copy, anatomical ones first
(the functional pipeline needs the anatomical outputs).

A report is written to <workdir>/smoke_report.md; the exit code is 1 if a run failed.

    python tests/e2e/run_smoke.py                    # every launcher with data here
    python tests/e2e/run_smoke.py --only anat
    python tests/e2e/run_smoke.py --launchers Exemples/Study_EDNiX/Launcher_anat/Launcher_anat_Mouse.py
    python tests/e2e/run_smoke.py --dry-run          # prepare data and launchers only
"""
import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'GUI'))
import ednix_gui_core as G  # noqa: E402

SKIP = {
    'preprocess_anat': ['itk_1', 'itk_2', 'itk_3', 9, 10, 11, 12, 'flat_map', 14, 15],
    'preprocess_fmri': ['itk_1', 'itk_2', 'itk_3', 6, 'Clean'],
}


def field(model, name):
    fs = [f for f in model.fields if f.name == name]
    return fs[-1] if fs else None


def pipeline(model):
    return model.calls[0][0] if model.calls else ''


def first_subject_session(bids):
    for sub in sorted(p for p in Path(bids).glob('sub-*') if p.is_dir()):
        sessions = sorted(p for p in sub.glob('ses-*') if p.is_dir())
        return sub.name[4:], (sessions[0].name[4:] if sessions else None)
    return None, None


def reduce_bids(src, dst, sub, ses, volumes):
    import nibabel as nib
    src, dst = Path(src), Path(dst)
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for f in src.iterdir():                       # dataset-level files (dataset_description.json, participants...)
        if f.is_file():
            shutil.copy2(f, dst / f.name)
    s_src = src / f'sub-{sub}' / (f'ses-{ses}' if ses else '')
    s_dst = dst / f'sub-{sub}' / (f'ses-{ses}' if ses else '')
    for f in s_src.rglob('*'):
        if f.is_dir():
            continue
        rel = f.relative_to(s_src)
        # raw data only: skip anything the pipeline wrote next to it
        if rel.parts[0] not in ('anat', 'func', 'fmap', 'pet', 'dwi') or len(rel.parts) > 2:
            continue
        out = s_dst / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        if rel.parts[0] == 'func' and f.name.endswith(('.nii', '.nii.gz')):
            img = nib.load(str(f))
            if img.ndim == 4 and img.shape[3] > volumes:
                img = img.slicer[..., :volumes]
            nib.save(img, str(out))
        else:
            shutil.copy2(f, out)
    return dst


def prepare(launcher, workdir, volumes):
    m = G.LauncherModel(launcher)
    fb = field(m, 'bids_dir')
    if fb is None or not Path(str(fb.value)).is_dir():
        return None, f'BIDS folder not on this machine ({fb.value if fb else "no bids_dir"})'
    bids = Path(str(fb.value))
    sub, ses = first_subject_session(bids)
    if sub is None:
        return None, f'no sub-* folder in {bids}'
    mini = Path(workdir) / 'data' / re.sub(r'[^\w.-]', '_', str(bids).strip('/'))
    values = {fb.key: str(mini) + ('/' if str(fb.value).endswith('/') else '')}
    for f in m.fields:
        if f.name == 'list_to_keep':
            values[f.key] = [(sub, ses if ses else '1')]
        elif f.name == 'list_to_remove':
            values[f.key] = []
        elif f.name == 'Skip_step' and isinstance(f.value, list):
            # run everything from the raw data, except interactive and FreeSurfer/surface steps
            # (the example launchers often skip early steps already done on their data)
            values[f.key] = list(SKIP.get(pipeline(m), []))
    src, _ = m.render(values)
    out = Path(workdir) / 'launchers' / (Path(launcher).stem + '_smoke.py')
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(src, encoding='utf-8')
    return dict(launcher=Path(launcher), smoke=out, bids=bids, mini=mini, sub=sub, ses=ses, pipeline=pipeline(m)), ''


def count_outputs(folder):
    return sum(1 for p in Path(folder).rglob('*') if p.is_file() and p.name.endswith(('.nii.gz', '.nii', '.gii', '.png')))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--launchers', nargs='*', help='launchers to test (default: every example launcher)')
    ap.add_argument('--only', default='all', choices=['all', 'anat', 'func'])
    ap.add_argument('--volumes', type=int, default=40, help='functional volumes kept per run')
    ap.add_argument('--workdir', default=str(Path(os.environ.get('TMPDIR', '/tmp')) / 'ednix_smoke'))
    ap.add_argument('--timeout-hours', type=float, default=12.0, help='per launcher')
    ap.add_argument('--dry-run', action='store_true', help='prepare data and launchers, do not run')
    a = ap.parse_args()

    launchers = [Path(p).resolve() for p in a.launchers] if a.launchers else \
        sorted((REPO / 'Exemples' / 'Study_EDNiX').glob('Launcher_*/*.py'))
    wanted = {'all': ('preprocess_anat', 'preprocess_fmri'), 'anat': ('preprocess_anat',), 'func': ('preprocess_fmri',)}[a.only]

    jobs, skipped = [], []
    for L in launchers:
        m = G.LauncherModel(L)
        if pipeline(m) not in wanted:
            continue
        job, why = prepare(L, a.workdir, a.volumes)
        (jobs.append(job) if job else skipped.append((L, why)))
    # one reduced copy per dataset, anatomical launchers first
    jobs.sort(key=lambda j: (str(j['bids']), j['pipeline'] != 'preprocess_anat'))
    made = set()
    for j in jobs:
        if j['bids'] not in made:
            reduce_bids(j['bids'], j['mini'], j['sub'], j['ses'], a.volumes)
            made.add(j['bids'])

    rows = []
    for j in jobs:
        if a.dry_run:
            rows.append((j, 'prepared', 0, 0.0, ''))
            continue
        log = Path(a.workdir) / 'logs' / (j['smoke'].stem + '.log')
        log.parent.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, EDNIX_STRICT='1', PYTHONPATH=str(REPO) + os.pathsep + os.environ.get('PYTHONPATH', ''))
        before, t0 = count_outputs(j['mini']), time.time()
        with open(log, 'w') as fh:
            try:
                rc = subprocess.run([sys.executable, '-u', str(j['smoke'])], cwd=REPO, env=env, stdout=fh,
                                    stderr=subprocess.STDOUT, timeout=a.timeout_hours * 3600).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
        new = count_outputs(j['mini']) - before
        status = 'passed' if rc == 0 and new > 0 else f'FAILED ({rc if rc else "no output"})'
        tail = ''
        if status != 'passed':
            lines = log.read_text(errors='replace').splitlines()
            tail = next((l for l in reversed(lines) if 'ERROR' in l or 'Error' in l), lines[-1] if lines else '')
        rows.append((j, status, new, (time.time() - t0) / 60, tail[:200]))
        print(f'{j["launcher"].name}: {status}', flush=True)

    report = Path(a.workdir) / 'smoke_report.md'
    with open(report, 'w') as fh:
        fh.write(f'# EDNiX smoke test, {datetime.datetime.now():%Y-%m-%d %H:%M}\n\n')
        fh.write(f'{a.volumes} functional volumes per run, one subject per dataset, EDNIX_STRICT=1.\n\n')
        fh.write('| launcher | subject | status | new images | minutes | last error |\n|---|---|---|---|---|---|\n')
        for j, st, n, mins, tail in rows:
            fh.write(f'| {j["launcher"].name} | {j["sub"]} {j["ses"] or ""} | {st} | {n} | {mins:.0f} | {tail} |\n')
        if skipped:
            fh.write('\nNot run:\n\n' + ''.join(f'- {L.name}: {why}\n' for L, why in skipped))
    print(f'report: {report}')
    sys.exit(1 if any(not r[1] in ('passed', 'prepared') for r in rows) else 0)


if __name__ == '__main__':
    main()
