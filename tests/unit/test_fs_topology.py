import os
from modalities.anat.freesurfer import FS_surf


def _setup(tmp_path):
    surf = tmp_path / 'sub' / 'surf'
    surf.mkdir(parents=True)
    (surf / 'lh.orig.nofix').write_text('mesh')
    return str(tmp_path)


def test_topology_retries_seeds_until_manifold(tmp_path, monkeypatch):
    fs = _setup(tmp_path)
    calls, euler = [], iter(['      2E=3F:  450956 != 450951 (5)', '      2E=3F:  450816 = 450816 (0)'])
    def do(cmd, diary):
        calls.append(cmd)
        return next(euler) if 'mris_euler_number' in cmd else ''
    monkeypatch.setattr(FS_surf.run_cmd, 'do', do)
    monkeypatch.setattr(FS_surf.run_cmd, 'run', lambda cmd, diary: ('', ''))
    monkeypatch.setattr(FS_surf.run_cmd, 'msg', lambda *a, **k: None)
    FS_surf._fix_topology(fs, 'sub', 0, None, 'S ', 'E ')
    fixes = [c for c in calls if 'mris_fix_topology' in c]
    assert len(fixes) == 2 and '-seed' not in fixes[0] and '-seed 1234 ' in fixes[1]
    assert os.path.isfile(os.path.join(fs, 'sub', 'surf', 'lh.orig'))


def test_topology_raises_when_no_seed_works(tmp_path, monkeypatch):
    fs = _setup(tmp_path)
    monkeypatch.setattr(FS_surf.run_cmd, 'do',
                        lambda cmd, diary: '2E=3F: 1 != 2 (1)' if 'mris_euler_number' in cmd else '')
    monkeypatch.setattr(FS_surf.run_cmd, 'run', lambda cmd, diary: ('', ''))
    monkeypatch.setattr(FS_surf.run_cmd, 'error', lambda m, d: m)
    try:
        FS_surf._fix_topology(fs, 'sub', 0, None, 'S ', 'E ')
    except RuntimeError as e:
        assert 'not manifold' in str(e)
    else:
        raise AssertionError('no error raised')


def test_anat_reference_and_loop1_skip(monkeypatch):
    from modalities.anat import _0_Pipeline_launcher as P
    default = ['x', '1', '/b/sub-x/ses-1']
    monkeypatch.delenv('EDNIX_ANAT_REF', raising=False)
    assert P._anat_ref('/b', default) == default
    monkeypatch.setenv('EDNIX_ANAT_REF', 'a|2')
    assert P._anat_ref('/b', default) == ['a', '2', '/b/sub-a/ses-2']
    stages = {'loop1'}
    items = [('a', '2', 'p1'), ('a', 1, 'p2'), ('b', '2', 'p3')]
    monkeypatch.setenv('EDNIX_ANAT_LOOP1_SKIP', 'a|2')
    assert P._loop1_items(stages, items) == items[1:]
    monkeypatch.delenv('EDNIX_ANAT_LOOP1_SKIP')
    assert P._loop1_items(stages, items) == items
    assert P._loop1_items({'loop2'}, items) == []


def test_topology_failure_is_a_surface_error(tmp_path, monkeypatch):
    fs = _setup(tmp_path)
    monkeypatch.setattr(FS_surf.run_cmd, 'do',
                        lambda cmd, diary: '2E=3F: 1 != 2 (1)' if 'mris_euler_number' in cmd else '')
    monkeypatch.setattr(FS_surf.run_cmd, 'run', lambda cmd, diary: ('', ''))
    monkeypatch.setattr(FS_surf.run_cmd, 'error', lambda m, d: m)
    import pytest
    with pytest.raises(FS_surf.SurfaceError):
        FS_surf._fix_topology(fs, 'sub', 0, None, 'S ', 'E ')
    with pytest.raises(FS_surf.SurfaceError):
        FS_surf._require(str(tmp_path / 'missing'), 'test', None)


def test_loop3_isolates_surface_steps():
    """steps 11 -> 15 sit in one try whose except catches FS_surf.SurfaceError only; step 16 is outside"""
    import ast, inspect
    from modalities.anat import _loop3
    run = next(n for n in ast.parse(inspect.getsource(_loop3)).body if isinstance(n, ast.FunctionDef) and n.name == 'run')
    tries = [n for n in run.body if isinstance(n, ast.Try)]
    assert len(tries) == 1
    t = tries[0]
    inside = {c.func.attr for c in ast.walk(ast.Module(body=t.body, type_ignores=[]))
              if isinstance(c, ast.Call) and isinstance(c.func, ast.Attribute)}
    assert {'Wcreate', 'Wconfig', 'Wsphere', 'Pcreate', 'FS_finalise', 'WB_prep'} <= inside
    assert 'anat_QC' not in inside
    assert ast.unparse(t.handlers[0].type) == 'FS_surf.SurfaceError'


def test_run_all_lists_surface_failures(tmp_path):
    import importlib.util
    from pathlib import Path
    repo = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location('run_all_ednix', repo / 'Exemples' / 'Study_EDNiX' / 'run_all_ednix.py')
    R = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(R)
    qc = tmp_path / 'sub-1' / 'ses-1' / 'anat' / 'QC_anat'
    qc.mkdir(parents=True)
    (qc / 'sub-1_ses-1_SURFACE_FAILED.txt').write_text('ERROR: lh.orig not manifold\n')
    (tmp_path / 'sub-2' / 'ses-1' / 'anat' / 'QC_anat').mkdir(parents=True)
    jobs = [dict(launcher=dict(bids=str(tmp_path))), dict(launcher=dict(bids=str(tmp_path) + '/'))]
    out = R.surface_failures(jobs)
    assert len(out) == 1 and out[0].endswith('ERROR: lh.orig not manifold')
