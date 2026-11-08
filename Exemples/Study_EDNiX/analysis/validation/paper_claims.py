"""
paper_claims.py -- link every number of the manuscript to the code line and the data that
produce it, and check that they agree.

Why not a Jupyter notebook
    A notebook mixes code, hidden execution state and outputs in one JSON file: it is hard
    to diff in git, can show numbers produced by cells run out of order, and cannot be run
    unattended after each re-analysis. Here the link "sentence -> number -> stats file ->
    code line -> datasets" is DATA (validation/paper_claims_registry.py, one entry per
    claim, reviewable line by line in a pull request) and this script CHECKS it against the
    current .docx, the current stats files and the current code. It runs in one command
    after run_paper.py, and from a notebook if wanted:
        from Exemples.Study_EDNiX.analysis.validation.paper_claims import check
        rep = check('EDNiX_PaperV3_Claude_v11.docx'); rep['claims']

What is checked, per claim of the registry
    1. TEXT    the claim's `prefix` is found in the manuscript (tracked insertions kept,
               deletions dropped); the number written right after it is read. Reported
               with its paragraph number, sentence number and section, so that it can be
               found with Ctrl+F (Word line numbers depend on the layout and cannot be
               computed from the file).
    2. VALUE   the number produced by the analysis is read from the stats file named by
               the claim (xlsx sheet + row filter + column, or regex in a txt file, or a
               small function for derived values) and compared with the written number at
               the precision written (0.021 -> +-0.0005; "P < 0.001" -> value < 0.001;
               "≈1" -> +-10% or the claim's tol).
    3. CODE    the function that writes this column is located in the source by its name
               and a literal that writes the value (AST, nothing is imported): file:line.
               If the literal is gone, the code changed after the registry was written.
    4. DATA    the datasets the text names (e.g. BIDS_Zhu_Garin) appear as literals in the
               loading code; the stats file is older than the last commit of the producing
               code -> STALE (re-run run_paper.py).
    Plus: METHODS checks (a statement of the Methods vs the code: constant, default
    argument or code snippet) and UNTRACED numbers (every number of the Results and
    legends that no claim covers), so that coverage is explicit.

Status of a claim
    Methods checks: OK, FAIL (the code does not do what the Methods say), NOT_IN_METHODS
    (the statement is absent from the manuscript), FILE_MISSING (script not in this checkout).
    OK, MISMATCH (written != computed), PENDING (placeholder [XX...] in the text; the
    computed value is reported so that it can be pasted), NOT_IN_TEXT (prefix not found:
    sentence rewritten), NO_SOURCE (stats file / sheet / row missing: run the analysis).

Usage
    python paper_claims.py MANUSCRIPT.docx                       # PAPER_DIR of paper_data
    python paper_claims.py MANUSCRIPT.docx --paper-dir /path/to/Paper --out /tmp/trace
Outputs <out>/paper_claims_report.xlsx (sheets claims, methods, untraced) and
<out>/paper_claims_report.md. Default <out> = <PAPER_DIR>/Traceability.
"""
import argparse
import ast
import datetime as _dt
import os
import re
import subprocess
import sys
import zipfile

import numpy as np
import pandas as pd

opj = os.path.join
HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)                                   # .../Study_EDNiX/analysis
REPO = os.path.dirname(os.path.dirname(os.path.dirname(ANALYSIS)))  # repository root
if REPO not in sys.path:
    sys.path.insert(0, REPO)

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'


# ══════════════════════════════════════════════════════════════════════════════
# 1. Manuscript text
# ══════════════════════════════════════════════════════════════════════════════
def norm(s):
    """Normalise spacing and signs so that the registry can be typed with ASCII:
    non-breaking / thin spaces -> space, unicode minus -> '-', runs of spaces -> one."""
    s = s.replace(' ', ' ').replace(' ', ' ').replace(' ', ' ')
    s = s.replace('−', '-')
    return re.sub(r'[ \t]+', ' ', s)


def read_docx(path):
    """Non-empty paragraphs as they read with tracked changes accepted: text of w:t
    (insertions included); w:delText (deletions) is a different tag and is skipped."""
    from lxml import etree
    with zipfile.ZipFile(path) as z:
        root = etree.fromstring(z.read('word/document.xml'))
    paras = []
    for p in root.iter(W + 'p'):
        t = ''.join(x.text or '' for x in p.iter(W + 't'))
        if t.strip():
            paras.append(norm(t))
    return paras


_ABBREV = r'(?<!Fig)(?<!Figs)(?<!e\.g)(?<!i\.e)(?<!al)(?<!vs)(?<!approx)(?<!Supp)(?<!ref)(?<!no)'
_SENT_END = re.compile(_ABBREV + r'[.;?!]\s+(?=[A-Z(\[])')


def sentence_index(text, pos):
    """1-based index of the sentence that contains character `pos`, and that sentence."""
    starts = [0] + [m.end() for m in _SENT_END.finditer(text)]
    k = max(i for i, s in enumerate(starts) if s <= pos)
    end = starts[k + 1] if k + 1 < len(starts) else len(text)
    return k + 1, text[starts[k]:end].strip()


def section_of(paras, i):
    """Nearest preceding short paragraph = heading (Results subsection, Methods...)."""
    for j in range(i, -1, -1):
        if len(paras[j]) < 90 and not re.search(r'\d{3,}', paras[j]):
            return paras[j].strip()
    return ''


# a number as written in the text: optional comparator, optional sign, thousands commas
NUM = re.compile(r'(?P<cmp>[<>≤≥≈~]\s*)?(?P<sgn>[-+])?\s?(?P<num>\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?)')
PLACEHOLDER = re.compile(r'\[XX[^\]]*\]')


_MASK = re.compile(r'ICC\([CA],\s?1\)|R²|r²|\{[^}]*\}|\[NEW-REF[^\]]*\]|\bFigs?\.\s*\d+[a-z]?'
                   r'|Supplementary (?:Figs?\.|Table)\s*\d+[a-z]?')


def read_written(text, start, k=0):
    """k-th number (or [XX] placeholder) after position `start` (0 = the first).
    ICC(C,1), citations and figure calls are masked so that their digits are not counted."""
    masked = _MASK.sub(lambda m: ' ' * len(m.group(0)), text)
    tok = re.compile(PLACEHOLDER.pattern + '|' + NUM.pattern)
    for j, m in enumerate(tok.finditer(masked, start)):
        if j < k:
            continue
        raw = text[m.start():m.end()]
        if raw.startswith('[XX'):
            return dict(kind='pending', raw=raw, pos=m.start())
        num = m.group('num')
        val = float(num.replace(',', ''))
        if m.group('sgn') == '-':
            val = -val
        dec = len(num.split('.')[1]) if '.' in num else 0
        return dict(kind='number', raw=raw.strip(), value=val, decimals=dec,
                    cmp=(m.group('cmp') or '').strip(), pos=m.start('num'), end=m.end())
    return None


# ══════════════════════════════════════════════════════════════════════════════
# 2. Computed values
# ══════════════════════════════════════════════════════════════════════════════
class Sources:
    """Stats files read once; path relative to PAPER_DIR."""

    def __init__(self, paper_dir):
        self.paper_dir = paper_dir
        self._x, self._t = {}, {}

    def path(self, rel):
        return opj(self.paper_dir, rel)

    def sheet(self, rel, sheet):
        key = (rel, sheet)
        if key not in self._x:
            self._x[key] = pd.read_excel(self.path(rel), sheet_name=sheet)
        return self._x[key]

    def text(self, rel):
        if rel not in self._t:
            with open(self.path(rel), encoding='utf-8') as fh:
                self._t[rel] = norm(fh.read())
        return self._t[rel]

    def cell(self, rel, sheet, where, column):
        """Value of `column` in the single row of `sheet` that matches every key of `where`."""
        df = self.sheet(rel, sheet)
        sel = df
        for k, v in where.items():
            sel = sel[sel[k].astype(str) == str(v)]
        if len(sel) != 1:
            raise LookupError(f'{len(sel)} rows match {where} in {rel}[{sheet}]')
        return float(sel[column].iloc[0])

    def regex(self, rel, pattern, group=1):
        m = re.search(pattern, self.text(rel))
        if not m:
            raise LookupError(f'/{pattern}/ not found in {rel}')
        return float(m.group(group).replace(',', ''))


def computed_value(claim, src):
    s = claim['source']
    if 'compute' in s:
        v = s['compute'](src)
    elif 'regex' in s:
        v = src.regex(s['file'], s['regex'], s.get('group', 1))
    else:
        v = src.cell(s['file'], s['sheet'], s.get('where', {}), s['column'])
    return claim.get('transform', lambda x: x)(v)


def compare(written, value, tol=None):
    """True if `value` rounds to what is written (comparators and ≈ handled)."""
    w, c = written['value'], written['cmp']
    if c in ('<',):
        return value < w
    if c in ('≤',):
        return value <= w + 1e-12
    if c in ('>',):
        return value > w
    if c in ('≥',):
        return value >= w - 1e-12
    if tol is None:
        tol = 0.1 * abs(w) if c in ('≈', '~') else 0.5 * 10 ** (-written['decimals'])
    return abs(value - w) <= tol + 1e-9 * max(1.0, abs(w))


# ══════════════════════════════════════════════════════════════════════════════
# 3. Code location (AST, nothing imported)
# ══════════════════════════════════════════════════════════════════════════════
_SRC = {}


def _source(rel):
    """(text, AST) of a repository file; ('', empty module) if the file is not in this
    checkout (e.g. a script still on an unmerged branch): reported as not found."""
    if rel not in _SRC:
        path = opj(REPO, rel)
        if not os.path.isfile(path):
            _SRC[rel] = ('', ast.parse(''))
        else:
            with open(path, encoding='utf-8') as fh:
                text = fh.read()
            _SRC[rel] = (text, ast.parse(text))
    return _SRC[rel]


def code_line(rel, func, key=None):
    """'file:line' of `key` inside function `func` of file `rel` (function line if key is
    None); None if the function or the key is not found."""
    text, tree = _source(rel)
    lines = text.splitlines()
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == func:
            if key is None:
                return f'{rel}:{node.lineno}'
            for ln in range(node.lineno, node.end_lineno + 1):
                if key in lines[ln - 1]:
                    return f'{rel}:{ln}'
            return None
    return None


def literals(rel):
    """Every string literal of a source file."""
    _, tree = _source(rel)
    return {n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)}


def last_commit_time(rel):
    """Unix time of the last commit of a file (its mtime outside git)."""
    if not os.path.isfile(opj(REPO, rel)):
        return 0
    try:
        out = subprocess.run(['git', '-C', REPO, 'log', '-1', '--format=%ct', '--', rel],
                             capture_output=True, text=True, timeout=20).stdout.strip()
        if out:
            return int(out)
    except (OSError, subprocess.SubprocessError):
        pass
    return int(os.path.getmtime(opj(REPO, rel)))


def module_constant(rel, name):
    _, tree = _source(rel)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == name for t in node.targets):
            return ast.literal_eval(node.value), node.lineno
    return None, None


def default_arg(rel, func, arg):
    _, tree = _source(rel)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == func:
            a = node.args
            pos = a.args[len(a.args) - len(a.defaults):]
            for p, d in list(zip(pos, a.defaults)) + list(zip(a.kwonlyargs, a.kw_defaults)):
                if p.arg == arg and d is not None:
                    try:
                        return ast.literal_eval(d), node.lineno
                    except ValueError:
                        return ast.unparse(d), node.lineno
    return None, None


# ══════════════════════════════════════════════════════════════════════════════
# 4. The check
# ══════════════════════════════════════════════════════════════════════════════
def locate(paras, prefix, within=None):
    """(paragraph index, position right after `prefix`) of the first match; `within` =
    another string that the paragraph must contain (to disambiguate)."""
    p = norm(prefix)
    for i, t in enumerate(paras):
        if within and norm(within) not in t:
            continue
        k = t.find(p)
        if k >= 0:
            return i, k + len(p)
    return None, None


def check_claims(paras, registry, src):
    rows, covered = [], set()
    for c in registry:
        r = dict(id=c['id'], figure=c.get('figure', ''), claim=c.get('what', ''))
        i, pos = locate(paras, c['prefix'], c.get('within'))
        written = read_written(paras[i], pos, c.get('k', 0)) if i is not None else None
        if written is None:
            r.update(status='NOT_IN_TEXT', prefix=c['prefix'])
        else:
            k, sent = sentence_index(paras[i], pos)
            r.update(paragraph=i + 1, sentence=k, section=section_of(paras, i - 1)[:60],
                     text=sent[:300], written=written['raw'])
            if written['kind'] == 'number':
                covered.add((i, written['pos']))
        # computed value
        try:
            val = computed_value(c, src)
            r['computed'] = val
        except (OSError, LookupError, KeyError, ValueError) as e:
            val = None
            r['computed'] = None
            r['source_error'] = str(e)[:200]
        s = c['source']
        r['source'] = s.get('file', '') + (f"[{s['sheet']}] {s.get('where', '')} -> {s['column']}"
                                           if 'sheet' in s else
                                           (f" /{s['regex']}/" if 'regex' in s else ' (computed)'))
        # code and data
        prod = c.get('producer')
        if prod:
            loc = code_line(*prod)
            if loc is None and not os.path.isfile(opj(REPO, prod[0])):
                loc = f'{prod[0]} -- FILE NOT IN THIS CHECKOUT (unmerged branch?)'
            r['code'] = loc or f'{prod[0]}::{prod[1]} -- key {prod[2]!r} NOT FOUND (code changed?)'
        data = c.get('data', [])
        if data:
            files = c.get('data_files') or ([prod[0]] if prod else [])
            lits = set().union(*(literals(f) for f in files)) if files else set()
            missing = [d for d in data if not any(d in x for x in lits)]
            r['data'] = ', '.join(data) + (f'  MISSING in {files}: {missing}' if missing else '  ok')
        if prod and 'file' in s and os.path.isfile(src.path(s['file'])):
            if os.path.getmtime(src.path(s['file'])) < last_commit_time(prod[0]):
                r['stale'] = 'stats older than the last commit of the producing code: re-run'
        # status
        if 'status' not in r:
            if written['kind'] == 'pending':
                r['status'] = 'PENDING'
            elif val is None:
                r['status'] = 'NO_SOURCE'
            else:
                r['status'] = 'OK' if compare(written, val, c.get('tol')) else 'MISMATCH'
        rows.append(r)
    return pd.DataFrame(rows), covered


def check_methods(paras, methods):
    rows = []
    for m in methods:
        i, _ = locate(paras, m['quote'])
        r = dict(id=m['id'], methods_text=m['quote'][:120], in_text='yes' if i is not None else 'NO',
                 paragraph=(i + 1) if i is not None else None, expected=repr(m['expected']))
        kind, rel = m['kind'], m['file']
        if kind == 'constant':
            v, ln = module_constant(rel, m['name'])
        elif kind == 'default':
            v, ln = default_arg(rel, m['func'], m['arg'])
        else:                                                    # 'contains'
            loc = code_line(rel, m['func'], m['snippet'])
            v, ln = (m['snippet'] if loc else None), (int(loc.rsplit(':', 1)[1]) if loc else None)
        r.update(code_value=repr(v), code=f'{rel}:{ln}' if ln else f'{rel} (not found)',
                 status=('FILE_MISSING' if not os.path.isfile(opj(REPO, rel)) else
                         'NOT_IN_METHODS' if i is None else 'OK' if v == m['expected'] else 'FAIL'),
                 note=m.get('note', ''))
        rows.append(r)
    return pd.DataFrame(rows)


_SKIP = re.compile(r'\{[^}]*\}|\[NEW-REF[^\]]*\]|\[XX[^\]]*\]|\bFigs?\.\s*\d+[a-z]?(?:[–,-]\s*\d*[a-z]?)*'
                   r'|Supplementary (?:Figs?\.|Table)\s*\d+(?:[–,-]\d+)*|\blevel[- ]\d|ses-\d+|sub-\d+'
                   r'|\b\d+(?:st|nd|rd|th)\b|BIDS_\w+|¹⁸F|\[\d+\]|\(\d{4}\)|step \d+|doi:\S+|Table \d')


def untraced(paras, covered):
    """Numbers of the Results and of the legends that no claim covers."""
    rows, in_results = [], False
    for i, t in enumerate(paras):
        if t.strip() == 'Results':
            in_results = True
        elif t.strip() in ('Discussion', 'Methods'):
            in_results = False
        legend = re.match(r'(Figure \d+|Supplementary Fig\. \d+) \|', t)
        if not (in_results or legend):
            continue
        masked = _SKIP.sub(lambda m: ' ' * len(m.group(0)), t)
        for m in NUM.finditer(masked):
            if (i, m.start('num')) in covered:
                continue
            k, _ = sentence_index(t, m.start())
            rows.append(dict(paragraph=i + 1, sentence=k, number=m.group(0).strip(),
                             section=(legend.group(1) if legend else section_of(paras, i - 1)[:60]),
                             context=t[max(0, m.start() - 60):m.end() + 30]))
    return pd.DataFrame(rows)


def check(docx, paper_dir=None, registry=None, methods=None):
    from Exemples.Study_EDNiX.analysis.validation import paper_claims_registry as reg
    if paper_dir is None:
        from Exemples.Study_EDNiX.analysis import paper_data as pdat
        paper_dir = pdat.PAPER_DIR
    paras = read_docx(docx)
    claims, covered = check_claims(paras, registry if registry is not None else reg.CLAIMS,
                                   Sources(paper_dir))
    meth = check_methods(paras, methods if methods is not None else reg.METHODS)
    unt = untraced(paras, covered)
    return dict(claims=claims, methods=meth, untraced=unt, docx=docx, paper_dir=paper_dir)


def write_report(rep, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    xlsx = opj(out_dir, 'paper_claims_report.xlsx')
    with pd.ExcelWriter(xlsx) as xw:
        for k in ('claims', 'methods', 'untraced'):
            rep[k].to_excel(xw, sheet_name=k, index=False)
    c, m, u = rep['claims'], rep['methods'], rep['untraced']
    counts = c['status'].value_counts().to_dict() if not c.empty else {}
    L = [f'# Paper claims -- {os.path.basename(rep["docx"])}',
         f'{_dt.datetime.now():%Y-%m-%d %H:%M}; stats from `{rep["paper_dir"]}`', '',
         f'**Claims**: {len(c)} registered; ' + ', '.join(f'{k} {v}' for k, v in counts.items()),
         f'**Methods vs code**: {int((m["status"] == "OK").sum()) if not m.empty else 0}/{len(m)} OK',
         f'**Untraced numbers** in Results and legends: {len(u)}', '',
         '| status | id | ¶ | sent. | written | computed | code | source |', '|---|---|---|---|---|---|---|---|']
    for r in c.itertuples():
        comp = '' if pd.isna(getattr(r, 'computed', np.nan)) else f'{r.computed:.6g}'
        L.append(f'| {r.status} | {r.id} | {getattr(r, "paragraph", "")} | {getattr(r, "sentence", "")} | '
                 f'{getattr(r, "written", "")} | {comp} | {getattr(r, "code", "")} | {r.source} |')
    L += ['', '## Methods vs code', '', '| status | id | Methods says | code | code value | note |',
          '|---|---|---|---|---|---|']
    for r in m.itertuples():
        L.append(f'| {r.status} | {r.id} | {r.methods_text} | {r.code} | {r.code_value} | {r.note} |')
    L += ['', '## Untraced numbers (first 200)', '', '| ¶ | sent. | number | context |', '|---|---|---|---|']
    for r in u.head(200).itertuples():
        L.append(f'| {r.paragraph} | {r.sentence} | {r.number} | {r.context.replace("|", "/")} |')
    md = opj(out_dir, 'paper_claims_report.md')
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L) + '\n')
    return xlsx, md


def main(argv=None):
    ap = argparse.ArgumentParser(description='Check every registered number of the manuscript '
                                             'against the stats files and the code.')
    ap.add_argument('docx', help='manuscript (.docx)')
    ap.add_argument('--paper-dir', default=None, help='folder of the stats (default paper_data.PAPER_DIR)')
    ap.add_argument('--out', default=None, help='report folder (default <paper-dir>/Traceability)')
    a = ap.parse_args(argv)
    rep = check(a.docx, a.paper_dir)
    out = a.out or opj(rep['paper_dir'], 'Traceability')
    xlsx, md = write_report(rep, out)
    print(open(md, encoding='utf-8').read().split('\n## Methods')[0])
    print(f'-> {xlsx}\n-> {md}')
    bad = rep['claims']['status'].isin(['MISMATCH', 'NOT_IN_TEXT']).sum()
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
