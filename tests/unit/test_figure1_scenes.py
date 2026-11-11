"""Figure 1 scenes: installation in the atlas library (paths relative to the library) and panel a assembly."""
import importlib.util
import json
import os
from os.path import join as opj

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _load():
    spec = importlib.util.spec_from_file_location(
        'install_figure_scenes', opj(ROOT, 'Exemples', 'Study_EDNiX', 'atlases', 'install_figure_scenes.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


SCENE = """<?xml version="1.0" encoding="UTF-8"?>
<SceneFile Version="3">
 <Object Type="pathName" Name="a">{rel}</Object>
 <Object Type="pathName" Name="b">{abs_}</Object>
 <Object Type="pathName" Name="c">custom.dscalar.nii</Object>
 <Object Type="pathName" Name="d">gone/Human.r.inflated.surf.gii</Object>
 <ObjectArray Type="pathName" Name="e" Length="1"><Element Index="0">{rel}</Element></ObjectArray>
 <Object Type="string" Name="surfaceFileName">Human.l.inflated.surf.gii</Object>
</SceneFile>
"""


@pytest.fixture
def setup(tmp_path):
    old = tmp_path / 'old' / 'Atlases_library'
    new = tmp_path / 'EDNiX' / 'Atlases_library_V2'
    frag = 'atlas/mammals/primates/Human/EDNiX/surfaces/Native_resol'
    for lib in (old, new):
        (lib / frag).mkdir(parents=True)
        for f in ('Human.l.inflated.surf.gii', 'Human.r.inflated.surf.gii'):
            (lib / frag / f).write_text('x')
    src = tmp_path / 'Figures'
    src.mkdir()
    (src / 'custom.dscalar.nii').write_text('d')
    rel = os.path.relpath(old / frag / 'Human.l.inflated.surf.gii', src)
    txt = SCENE.format(rel=rel, abs_=str(old / frag / 'Human.l.inflated.surf.gii'))
    (src / 'Human.scene').write_text(txt)
    (src / 'Human_1.scene').write_text(txt)
    (src / 'human_1.png').write_bytes(b'png')
    (src / 'Unicorn.scene').write_text(txt)
    (new / 'EDNiX_library_manifest.json').write_text(json.dumps({'library': 'V2', 'files': {}}))
    return src, new, frag


def test_source_names():
    m = _load()
    assert m.re.fullmatch(r'([A-Za-z]+?)((?:_thickness)?(?:_1)?)', 'Mouselemur_1').groups() == ('Mouselemur', '_1')


def test_install_rewrites_relative(setup):
    src, lib, frag = setup
    m = _load()
    written = m.install(str(lib), str(src))
    assert set(written) == {'Human'}
    dst = lib / 'atlas/mammals/primates/Human/EDNiX/scenes'
    names = sorted(os.listdir(dst))
    assert names == ['Human_Figure1.scene', 'Human_Figure1_white.scene', 'Human_Figure1_white_reference.png', 'data']
    text = (dst / 'Human_Figure1_white.scene').read_text()
    paths = m._PATH_OBJECT.findall(text) + [('', e[1], '') for e in m._ELEMENT.findall(text)]
    for _, p, _ in paths:
        assert not os.path.isabs(p)
        assert os.path.isfile(os.path.normpath(opj(dst, p))), p
    assert '>Human.l.inflated.surf.gii</Object>' in text           # string objects untouched
    # the library can be renamed: relative paths still resolve
    moved = lib.parent / 'Atlases_library'
    lib.rename(moved)
    d2 = moved / 'atlas/mammals/primates/Human/EDNiX/scenes'
    for _, p, _ in m._PATH_OBJECT.findall((d2 / 'Human_Figure1.scene').read_text()):
        assert os.path.isfile(os.path.normpath(opj(d2, p)))
    man = json.loads((moved / 'EDNiX_library_manifest.json').read_text())
    assert any(k.endswith('Human_Figure1.scene') for k in man['files'])


def test_plan_writes_nothing(setup):
    src, lib, _ = setup
    _load().install(str(lib), str(src), plan=True)
    assert not (lib / 'atlas/mammals/primates/Human/EDNiX/scenes').exists()


def test_figure1_assembly_falls_back_on_reference(setup, tmp_path, monkeypatch):
    import matplotlib.image as mpimg
    src, lib, _ = setup
    img = np.ones((60, 80, 3))
    img[10:50, 10:70] = [0.8, 0.2, 0.2]
    mpimg.imsave(src / 'human_1.png', img)
    mpimg.imsave(src / 'human.png', 1 - img)
    _load().install(str(lib), str(src))
    from Exemples.Study_EDNiX.analysis.library import figure1_surfaces as f1s
    parts = tmp_path / 'parts'
    parts.mkdir()
    for n in ('cortical', 'subcortical', 'legend'):
        mpimg.imsave(parts / f'ednix_sunburst_{n}.png', img)
    out = tmp_path / 'Figure_1' / 'Figure_1'
    paths = f1s.figure1(str(parts), str(out), library=str(lib), sing_wb='true || ')
    assert all(os.path.isfile(p) for p in paths)
    assert os.path.isfile(tmp_path / 'Figure_1' / 'surfaces' / 'Figure_1a_black.png')
