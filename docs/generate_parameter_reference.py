#!/usr/bin/env python3
"""
Regenerate docs/reference/*.md from the pipeline code: parameters, defaults and
descriptions come from the signatures and docstrings of the pipeline entry points,
and the step identifiers from GUI/ednix_gui_core.py (the list shown in the interface).

    python docs/generate_parameter_reference.py
"""
import ast
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'GUI'))
import ednix_gui_core as G  # noqa: E402

ENTRY_POINTS = [
    ('modalities/anat/_0_Pipeline_launcher.py', 'preprocess_anat', 'anat_parameters.md', 'Anatomical pipeline'),
    ('modalities/fMRI/_0_Pipeline_launcher.py', 'preprocess_fmri', 'fmri_parameters.md', 'Functional pipeline'),
    ('modalities/pet/_0_Pipeline_launcher.py', 'preprocess_data', 'pet_parameters.md', 'PET pipeline (in development)'),
]


def _cell(text):
    return ' '.join(str(text).split()).replace('|', '\\|')


def reference(path, func):
    src = (REPO / path).read_text(encoding='utf-8')
    fn = next((n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == func), None)
    if fn is None:
        return None, False
    docs = G._parse_docstring(ast.get_docstring(fn) or '')
    args = fn.args.args
    defaults = [None] * (len(args) - len(fn.args.defaults)) + list(fn.args.defaults)
    rows = []
    for a, d in zip(args, defaults):
        default = 'required' if d is None else f'`{_cell(ast.get_source_segment(src, d))}`'
        rows.append((a.arg, default, docs.get(a.arg) or G.CURATED.get(a.arg, '')))
    return rows, bool(fn.args.kwarg)


def main():
    out_dir = REPO / 'docs' / 'reference'
    out_dir.mkdir(parents=True, exist_ok=True)
    for path, func, out, title in ENTRY_POINTS:
        rows, has_kwargs = reference(path, func)
        lines = [f'# {title}: parameters of `{func}`', '',
                 f'Generated from `{path}` by `docs/generate_parameter_reference.py`; do not edit by hand.', '']
        if rows is None:
            lines += [f'`{func}` is not active in this version of `{path}` (the pipeline is under development).']
            (out_dir / out).write_text('\n'.join(lines) + '\n', encoding='utf-8')
            print(f'{out_dir / out}: entry point not active')
            continue
        steps = G.STEPS.get(func)
        if steps:
            lines += ['## Steps (`Skip_step`)', '',
                      'Identifiers accepted by `Skip_step`: the listed steps are skipped.', '',
                      '| step | what it does |', '|---|---|']
            lines += [f'| `{s}` | {_cell(label)} |' for s, label in steps]
            lines += ['']
        lines += ['## Parameters', '', '| parameter | default | description |', '|---|---|---|']
        lines += [f'| `{n}` | {d} | {_cell(t) or "not documented"} |' for n, d, t in rows]
        if has_kwargs:
            lines += ['', f'`{func}` also accepts other keyword arguments, which it ignores: a misspelled '
                          'parameter name is therefore reported as a warning, not an error. '
                          'The interface and the unit tests check the example launchers for this.']
        (out_dir / out).write_text('\n'.join(lines) + '\n', encoding='utf-8')
        print(f'{out_dir / out}: {len(rows)} parameters')


if __name__ == '__main__':
    main()
