"""Every module of the pipeline imports (PET, archives and legacy scripts excluded)."""
import importlib
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import REPO

EXCLUDE = ('/pet/', 'Archives', 'UNet_Model_training', '_T1_T2_quantif', 'old_orient_correct', 'Tools/niiRoi.py',
           'Tools/Clean_BIDS.py')
MODULES = sorted(
    str(p.relative_to(REPO))[:-3].replace('/', '.')
    for d in ('modalities', 'Tools', 'atlases', 'Statistics')
    for p in (REPO / d).rglob('*.py')
    if not any(x in str(p) for x in EXCLUDE) and '.ipynb_checkpoints' not in str(p)
)


@pytest.mark.slow
@pytest.mark.parametrize('module', MODULES)
def test_import(module):
    # a fresh interpreter per module: some modules run code at import time
    r = subprocess.run([sys.executable, '-W', 'ignore', '-c', f'import {module}'], cwd=REPO,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr.strip().splitlines()[-1] if r.stderr else 'import failed'


def test_every_python_file_compiles():
    bad = []
    for p in REPO.rglob('*.py'):
        if '.ipynb_checkpoints' in str(p):
            continue
        try:
            compile(p.read_text(encoding='utf-8'), str(p), 'exec')
        except SyntaxError as e:
            bad.append(f'{p}: line {e.lineno}')
    assert not bad, bad
