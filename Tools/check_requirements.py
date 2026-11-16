#!/usr/bin/env python3
"""
Compare the third-party modules imported by the EDNiX code with requirements.txt and requirements-lock.txt.

    python3 Tools/check_requirements.py            report; exit 1 if an import is not covered
    python3 Tools/check_requirements.py --verbose  also list the files importing each module

Every `import x` / `from x import y` (absolute imports, at any depth: lazy imports inside functions count) of
every .py file of the repository is read with the AST, without importing anything. Standard-library modules
(sys.stdlib_module_names) and the repository's own modules are ignored. Archives/ folders and the modules
listed in OPTIONAL are reported but do not fail the check.
"""
import argparse
import ast
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# import name -> pip distribution name, where they differ
IMPORT_TO_DIST = {
    'ants': 'antspyx',
    'bids': 'pybids',
    'skimage': 'scikit-image',
    'sklearn': 'scikit-learn',
    'yaml': 'pyyaml',
    'PIL': 'pillow',
    'cv2': 'opencv-python',
}

# Imported only by test code or by scripts that are not part of a run; not required in the venv.
OPTIONAL = {
    'pytest': 'tests only (pip install pytest pytest-xdist)',
}

SKIP_DIRS = ('.git', '.ipynb_checkpoints', '__pycache__', 'venv', '.venv')


def _norm(name):
    return re.sub(r'[-_.]+', '-', name).lower()


def local_modules(repo=REPO):
    names = set()
    for p in repo.rglob('*'):
        if any(s in p.parts for s in SKIP_DIRS):
            continue
        if p.is_dir():
            names.add(p.name)
        elif p.suffix == '.py':
            names.add(p.stem)
    return names


def scan_imports(repo=REPO):
    """{top-level module: set of files} for every third-party import of the repository."""
    std = set(sys.stdlib_module_names)
    local = local_modules(repo)
    uses = defaultdict(set)
    for p in repo.rglob('*.py'):
        if any(s in p.parts for s in SKIP_DIRS) or any(s.startswith('venv') for s in p.parts):
            continue
        try:
            tree = ast.parse(p.read_text(encoding='utf-8', errors='replace'))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                mods = [node.module]
            else:
                continue
            for m in mods:
                top = m.split('.')[0]
                if top in std or top in local or top == '__future__':
                    continue
                uses[top].add(str(p.relative_to(repo)))
    return uses


def read_requirements(path):
    """Normalised distribution names of a requirements file (comments, markers and versions dropped)."""
    dists = set()
    if not path.exists():
        return dists
    for line in path.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if not line or line.startswith('-'):
            continue
        name = re.split(r'[\s;<>=!~\[]', line, maxsplit=1)[0]
        dists.add(_norm(name))
    return dists


def check(repo=REPO):
    """(missing_in_requirements, missing_in_lock, archive_only, optional) -> lists of (module, dist, files)."""
    uses = scan_imports(repo)
    req = read_requirements(repo / 'requirements.txt')
    lock = read_requirements(repo / 'requirements-lock.txt')
    missing_req, missing_lock, archive_only, optional = [], [], [], []
    for mod in sorted(uses):
        dist = IMPORT_TO_DIST.get(mod, mod)
        files = sorted(uses[mod])
        if mod in OPTIONAL:
            optional.append((mod, dist, files))
            continue
        if all('Archives' in f or '/old/' in f for f in files):
            if _norm(dist) not in req:
                archive_only.append((mod, dist, files))
            continue
        if _norm(dist) not in req:
            missing_req.append((mod, dist, files))
        if lock and _norm(dist) not in lock:
            missing_lock.append((mod, dist, files))
    return missing_req, missing_lock, archive_only, optional


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    missing_req, missing_lock, archive_only, optional = check()

    def show(title, rows):
        if not rows:
            return
        print(title)
        for mod, dist, files in rows:
            print(f'  {mod:18s} (pip: {dist})  {len(files)} file(s)')
            if a.verbose:
                for f in files:
                    print(f'      {f}')

    show('NOT in requirements.txt:', missing_req)
    show('NOT in requirements-lock.txt:', missing_lock)
    show('Imported only by archived scripts (not required):', archive_only)
    show('Optional:', optional)
    if missing_req or missing_lock:
        sys.exit(1)
    print('OK: every third-party import of the pipeline is in requirements.txt and requirements-lock.txt.')


if __name__ == '__main__':
    main()
