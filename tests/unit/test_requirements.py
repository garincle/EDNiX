"""requirements.txt and requirements-lock.txt cover every third-party import of the code (AST scan)."""
import importlib.util

from conftest import REPO

spec = importlib.util.spec_from_file_location('check_requirements', REPO / 'Tools' / 'check_requirements.py')
check_requirements = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_requirements)


def test_requirements_cover_imports():
    missing_req, missing_lock, _, _ = check_requirements.check(REPO)
    assert not missing_req, f'add to requirements.txt: {[(m, d, f[:3]) for m, d, f in missing_req]}'
    assert not missing_lock, f'regenerate requirements-lock.txt: {[(m, d) for m, d, _ in missing_lock]}'


def test_lock_is_fully_pinned():
    bad = []
    for line in (REPO / 'requirements-lock.txt').read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if line and not line.startswith('-') and '==' not in line:
            bad.append(line)
    assert not bad, f'unpinned lines in requirements-lock.txt: {bad}'


def test_lock_contains_requirements():
    req = check_requirements.read_requirements(REPO / 'requirements.txt')
    lock = check_requirements.read_requirements(REPO / 'requirements-lock.txt')
    assert req <= lock, f'in requirements.txt but not locked: {sorted(req - lock)}'


def test_freesurfer_license_is_bound_from_the_user(tmp_path, monkeypatch):
    """The released freesurfer_NHP.sif has no key: the user's license is passed via FS_LICENSE."""
    from Tools import Load_EDNiX_requirement as req
    lic = tmp_path / 'fs' / 'license.txt'
    lic.parent.mkdir()
    lic.write_text('user@example.org\n')
    monkeypatch.setenv('HOME', str(tmp_path / 'home'))
    monkeypatch.setenv('FS_LICENSE', str(lic))
    monkeypatch.delenv('SINGULARITYENV_FS_LICENSE', raising=False)
    main = tmp_path / 'EDNiX'
    main.mkdir()
    sing_fs = req.load_requirement(str(main), str(tmp_path / 'lib'), str(tmp_path / 'bids'), 'yes')[2]
    assert req.os.environ['SINGULARITYENV_FS_LICENSE'] == str(lic)
    assert str(lic.parent) in sing_fs.split('--bind ')[1].split(' ')[0].split(',')


def test_no_freesurfer_license_keeps_the_old_behaviour(tmp_path, monkeypatch, capsys):
    from Tools import Load_EDNiX_requirement as req
    monkeypatch.setenv('HOME', str(tmp_path / 'home'))
    monkeypatch.delenv('FS_LICENSE', raising=False)
    monkeypatch.delenv('SINGULARITYENV_FS_LICENSE', raising=False)
    main = tmp_path / 'EDNiX'
    main.mkdir()
    sing_fs = req.load_requirement(str(main), str(tmp_path / 'lib'), str(tmp_path / 'bids'), 'yes')[2]
    assert 'SINGULARITYENV_FS_LICENSE' not in req.os.environ
    assert sing_fs.startswith('singularity run --bind ')
    assert 'no FreeSurfer license found' in capsys.readouterr().out
