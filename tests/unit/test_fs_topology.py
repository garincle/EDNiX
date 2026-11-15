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
