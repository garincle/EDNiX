"""rebuild_all_ednix: library built -> verified -> swapped -> reset -> run_all; stops before touching the BIDS
folders if the library build or its verification fails; resumes when the library is already V2."""
import importlib.util
import json
import os
import subprocess
from types import SimpleNamespace

import pytest

HERE = os.path.dirname(__file__)
SCRIPT = os.path.join(HERE, '..', '..', 'Exemples', 'Study_EDNiX', 'rebuild_all_ednix.py')


@pytest.fixture
def M(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('rebuild_all_ednix', SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    lib, v2 = tmp_path / 'Atlases_library', tmp_path / 'Atlases_library_V2'
    lib.mkdir()
    (lib / 'v1.txt').write_text('v1')
    monkeypatch.setattr(m, 'LIB', lib)
    monkeypatch.setattr(m, 'LIB_V2', v2)
    calls = []

    def fake_run(cmd, cwd=None):
        calls.append(cmd)
        if cmd[-1] == '--verify':
            return SimpleNamespace(returncode=m._verify_code)
        if str(cmd[1]).endswith('build_library_V2.py'):
            v2.mkdir(exist_ok=True)
            (v2 / 'EDNiX_library_manifest.json').write_text(json.dumps(
                dict(library='V2', species={'Rat/EDNiX': dict(status=m._species_status)})))
        return SimpleNamespace(returncode=0)

    m._verify_code, m._species_status = 0, 'ok'
    monkeypatch.setattr(subprocess, 'run', fake_run)
    monkeypatch.setattr(m.R, 'discover', lambda only: [])
    monkeypatch.setattr(m.R, 'reset_all', lambda ls, confirm=True, preview_only=False: calls.append('reset') or True)
    m._calls = calls
    return m


def _args(**k):
    return SimpleNamespace(only=[], phases=['anat'], skip_library=False, no_reset=False, **k)


def test_full_sequence(M):
    with pytest.raises(SystemExit) as e:
        M.run(_args())
    assert e.value.code == 0
    assert (M.LIB / 'EDNiX_library_manifest.json').exists()                 # V2 swapped in
    assert (M.LIB.with_name('Atlases_library_V1') / 'v1.txt').exists()      # V1 kept
    assert 'reset' in M._calls and str(M._calls[-1][1]).endswith('run_all_ednix.py')


@pytest.mark.parametrize('what', ['species', 'verify'])
def test_stops_before_bids(M, what):
    if what == 'species':
        M._species_status = 'failed'
    else:
        M._verify_code = 1
    with pytest.raises(SystemExit) as e:
        M.run(_args())
    assert 'STOP' in str(e.value.code)
    assert 'reset' not in M._calls and (M.LIB / 'v1.txt').exists()


def test_resume_when_library_is_v2(M):
    (M.LIB / 'EDNiX_library_manifest.json').write_text(json.dumps(dict(library='V2')))
    with pytest.raises(SystemExit):
        M.run(_args())
    assert not any('build_library_V2.py' in str(c[1]) for c in M._calls if isinstance(c, list))
    assert 'reset' in M._calls
