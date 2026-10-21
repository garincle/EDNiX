"""
Core of the EDNiX desktop GUI: launcher parsing/rewriting, atlas library
access, 3D surface meshes of atlas regions, background runs.
No GUI toolkit here: the web front-end talks to ednix_gui_server.py.
"""
import ast, os, re, sys, io, json, tokenize, difflib, subprocess, datetime, html
from pathlib import Path
import numpy as np

# =============================================================================
# Launcher model: parse a launcher script, expose its literal parameters
# =============================================================================
PIPELINE_FUNCS = {'preprocess_anat', 'preprocess_fmri', 'preprocess_pet', 'preprocess_ct', 'preprocess_data'}

STATIC_OPTIONS = {
    'which_on': ['all', 'max'],
    'MNIBcorrect_indiv': ['', 'N4', 'N3'],
    'normalize': ['Skip', 'zscore', 'psc'],
    'Align_img_to_template': ['Ants', 'No', ''],
    'doWARPonfunc': ['header', 'WARP'],
    'registration_type': ['origin', 'fast', 'V2'],
    'post_treatment_method': ['AFNI', 'Grandjean', 'fMRIPrep_style'],
    'fMRImasks': ['aseg', 'custom'],
    'transfo_message': ['do_it_for_me', 'do_as_I_said'],
    'resting_or_task': ['resting', 'task'],
    'orientation': ['', 'RAS', 'LPI', 'RPI', 'LAI', 'LIA', 'LSP', 'RAI'],
    'ICA_cleaning': ['Skip', 'Applied'],
    'n_for_ANTS': ['linear', 'lanczosWindowedSinc', 'hammingWindowedSinc', 'bSpline'],
    'type_of_transform': ['Rigid', 'Affine', 'SyN', 'SyNCC', 'BOLDRigid', 'BOLDAffine'],
    # PET / CT pipelines
    'crop_for_MoCo': ['auto', 'manual', 'none'],
    'export_resolution': ['pet', 'target'],
    'deoblique': ['WARP', 'header'],
    'MoCo_transform': ['Rigid', 'BOLDRigid'],
    'coreg_metric': ['mattes', 'meansquares', 'CC'],
    'type_of_transform_template': ['SyN', 'SyNCC', 'Affine'],
    'studyT_transform': ['SyN', 'Affine', 'Rigid'],
    'pvc_method': ['RBV', 'GTM', 'MG', 'LABBE', 'MTC', 'IY'],
    'T_PET': ['T1w', 'T2w'],
    'T_CT': ['T1w', 'T2w'],
}


def _branch_names(path, var):
    """String values compared to `var` (== 'X') in a module, e.g. skull-strip method names."""
    try:
        src = Path(path).read_text(encoding='utf-8')
    except OSError:
        return []
    return sorted(set(re.findall(var + r"\s*==\s*'([^']+)'", src)))


def discover_options(repo):
    repo = Path(repo)
    opts = dict(STATIC_OPTIONS)
    ss = _branch_names(repo / 'modalities/anat/skullstrip/Skullstrip_method.py', 'brain_skullstrip')
    for k in ('brain_skullstrip_1', 'brain_skullstrip_2', 'template_skullstrip'):
        opts[k] = sorted(set(ss + ['NoSkullStrip', 'Manual']))
    opts['Method_mask_func'] = _branch_names(repo / 'modalities/fMRI/Skullstrip_func.py', 'brain_skullstrip')
    return opts


class Field:
    def __init__(self, key, name, kind, start, end, value, section, help_txt, where):
        self.key, self.name, self.kind = key, name, kind
        self.start, self.end = start, end            # byte offsets in the source
        self.value, self.section, self.help, self.where = value, section, help_txt, where


def _line_offsets(src_bytes):
    offs, pos = [0], 0
    for line in src_bytes.split(b'\n'):
        pos += len(line) + 1
        offs.append(pos)
    return offs


def _span(node, offs):
    return offs[node.lineno - 1] + node.col_offset, offs[node.end_lineno - 1] + node.end_col_offset


def _comments(src):
    """{line: comment text} and set of lines that contain only a comment."""
    out, alone = {}, set()
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except (tokenize.TokenError, IndentationError):
        return out, alone
    code_lines = set()
    for t in toks:
        if t.type == tokenize.COMMENT:
            out[t.start[0]] = t.string.lstrip('#').strip()
        elif t.type not in (tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT,
                            tokenize.ENCODING, tokenize.ENDMARKER):
            code_lines.update(range(t.start[0], t.end[0] + 1))
    alone = {ln for ln in out if ln not in code_lines}
    return out, alone


def _is_header(raw_line):
    s = raw_line.strip()
    return s.startswith('###') or '----' in s or '====' in s


class LauncherModel:
    def __init__(self, path):
        self.path = Path(path)
        self.src = self.path.read_text(encoding='utf-8')
        self.tree = ast.parse(self.src)
        self.fields = []
        self.calls = []            # (func_name, module_attr_base, call_node)
        self.warnings = []
        self._parse()

    def _parse(self):
        src_b = self.src.encode('utf-8')
        offs = _line_offsets(src_b)
        lines = self.src.split('\n')
        comments, alone = _comments(self.src)
        section, pending_help = 'General', []
        body_by_line = {n.lineno: n for n in self.tree.body}
        last_line = 0

        def absorb_comments(upto):
            nonlocal section, pending_help
            for ln in range(last_line + 1, upto):
                if ln in alone:
                    if _is_header(lines[ln - 1]):
                        txt = comments[ln].strip('#-= ').strip()
                        if txt:
                            section, pending_help = txt[:70], []
                    else:
                        pending_help.append(comments[ln])

        for node in self.tree.body:
            absorb_comments(node.lineno)
            last_line = node.end_lineno
            if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
                name = node.targets[0].id
                help_txt = ' '.join(pending_help + ([comments[node.lineno]] if node.lineno in comments
                                                    and node.lineno not in alone else []))
                pending_help = []
                try:
                    val = ast.literal_eval(node.value)
                    s, e = _span(node.value, offs)
                    self.fields.append(Field(f'var:{name}:{node.lineno}', name, 'literal', s, e, val, section, help_txt, 'variable'))
                except (ValueError, SyntaxError, TypeError):
                    consts = [n for n in ast.walk(node.value) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
                    if len(consts) == 1:
                        s, e = _span(consts[0], offs)
                        self.fields.append(Field(f'var:{name}:{node.lineno}', name, 'string', s, e, consts[0].value,
                                                 section, help_txt, 'variable (text inside an expression)'))
            else:
                pending_help = []

        seen = {}
        for f in self.fields:
            seen.setdefault(f.name, []).append(f)
        for name, fs in seen.items():
            if len(fs) > 1:
                self.warnings.append(f'{name} is defined {len(fs)} times: only the last value is used.')

        for call in [n for n in ast.walk(self.tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)]:
            attr = call.func.attr
            base = call.func.value.id if isinstance(call.func.value, ast.Name) else ''
            if attr in PIPELINE_FUNCS:
                self.calls.append((attr, base, call))
                sec = f'Call to {base}.{attr}'
            elif attr == 'reset_all':
                sec = 'Data reset'
                for kw in call.keywords:
                    if kw.arg == 'dry_run' and isinstance(kw.value, ast.Constant) and kw.value.value is False:
                        self.warnings.append(f'line {call.lineno}: {base}.reset_all(dry_run=False) deletes results '
                                             'every time this launcher runs.')
            else:
                continue
            for kw in call.keywords:
                if kw.arg is None:
                    continue
                try:
                    val = ast.literal_eval(kw.value)
                except (ValueError, SyntaxError, TypeError):
                    continue
                s, e = _span(kw.value, offs)
                self.fields.append(Field(f'kw:{attr}:{kw.arg}:{call.lineno}', kw.arg, 'literal', s, e, val,
                                         sec, '', f'argument of {attr}()'))

    def render(self, new_values):
        """Source with the changed values substituted (new_values: {field.key: python value})."""
        src_b = bytearray(self.src.encode('utf-8'))
        changes = []
        for f in self.fields:
            if f.key in new_values and new_values[f.key] != f.value:
                changes.append((f.start, f.end, repr(new_values[f.key])))
        for s, e, txt in sorted(changes, reverse=True):
            src_b[s:e] = txt.encode('utf-8')
        return src_b.decode('utf-8'), len(changes)

    @staticmethod
    def check_calls(src, repo):
        """Keyword arguments that the called pipeline function does not declare (silently lost in **kwargs)."""
        repo = Path(repo)
        tree = ast.parse(src)
        issues = []
        for call in [n for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                     and n.func.attr in PIPELINE_FUNCS]:
            fname = call.func.attr
            params = None
            for mod in repo.glob('modalities/*/_0_Pipeline_launcher.py'):
                for n in ast.parse(mod.read_text(encoding='utf-8')).body:
                    if isinstance(n, ast.FunctionDef) and n.name == fname:
                        a = n.args
                        params = {x.arg for x in a.posonlyargs + a.args + a.kwonlyargs}
                        break
                if params is not None:
                    break
            if params is None:
                issues.append(f'{fname} not found in modalities/*/_0_Pipeline_launcher.py')
                continue
            for kw in call.keywords:
                if kw.arg and kw.arg not in params:
                    issues.append(f'line {call.lineno}: {fname}() has no parameter {kw.arg!r}; it would be silently ignored')
        return issues


# =============================================================================
# Atlas library model
# =============================================================================
ATLAS_RE = re.compile(r'^(?P<prefix>.+?)_seg-(?P<atlas>[^_]+)(?:_desc-(?P<desc>[^_]+))?_dseg\.nii(\.gz)?$')


def list_species(lib):
    lib = Path(lib)
    out = {}
    for dirpath, dirnames, _ in os.walk(lib / 'atlas'):
        if 'label_code' in dirnames and Path(dirpath).name != 'freesurfer':
            out.setdefault(Path(dirpath).name, Path(dirpath))
    return dict(sorted(out.items()))


def list_references(species_dir):
    return sorted(p.name for p in Path(species_dir).iterdir()
                  if p.is_dir() and (p / 'volumes' / 'labels').is_dir())


def list_atlases(ref_dir):
    out = {}
    for f in sorted((Path(ref_dir) / 'volumes' / 'labels').glob('*_dseg.nii*'), key=lambda p: p.name.replace('_', ' ')):
        m = ATLAS_RE.match(f.name)
        if m:
            key = m['atlas'] + (f" ({m['desc']})" if m['desc'] else '')
            out[key] = (m['atlas'], f)
    return out


def list_templates(ref_dir):
    vol = Path(ref_dir) / 'volumes'
    return sorted(f.name for f in vol.glob('*.nii*') if 'template' in f.name or 'SS' in f.name) or \
        sorted(f.name for f in vol.glob('*.nii*'))


# =============================================================================
# Launcher fields <-> JSON
# =============================================================================
def field_kind(v):
    if isinstance(v, bool):
        return 'bool'
    if isinstance(v, str):
        return 'str'
    if isinstance(v, (int, float)):
        return 'num'
    return 'list'


def launcher_to_json(model, options, docs=None):
    docs = docs or {}
    pipeline = model.calls[0][0] if model.calls else ''
    pdoc = docs.get(pipeline, {})
    out = []
    for f in model.fields:
        k = field_kind(f.value)
        out.append(dict(key=f.key, name=f.name, section=f.section, kind=k, where=f.where,
                        doc=pdoc.get(f.name) or CURATED.get(f.name, ''),
                        note=f.help,
                        value=f.value if k in ('bool', 'str') else repr(f.value),
                        options=options.get(f.name, [])))
    steps = [dict(id=s_id, label=lbl) for s_id, lbl in STEPS.get(pipeline, [])]
    return dict(path=str(model.path), pipeline=pipeline, fields=out, warnings=model.warnings, steps=steps)


def values_from_json(model, raw):
    """Convert the front-end values ({key: bool | str}) to Python values; returns (values, errors)."""
    vals, errors = {}, []
    for f in model.fields:
        if f.key not in raw:
            continue
        r = raw[f.key]
        k = field_kind(f.value)
        try:
            if k == 'bool':
                vals[f.key] = bool(r)
            elif k == 'str':
                t = str(r)
                vals[f.key] = ast.literal_eval(t) if t.strip()[:1] in ('[', '(', '{') else t
            else:
                vals[f.key] = ast.literal_eval(str(r))
        except (ValueError, SyntaxError) as e:
            errors.append(f'{f.name}: invalid value ({e})')
    return vals, errors


# =============================================================================
# 3D meshes of atlas regions
# =============================================================================
MESH_VERSION = 'v4'


def _b64(arr):
    import base64
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode('ascii')


def _surface(mask_f, level, step, offset, affine):
    from skimage.measure import marching_cubes
    v, f, _, _ = marching_cubes(mask_f, level, step_size=step, allow_degenerate=False)
    v = v + offset
    world = v @ affine[:3, :3].T + affine[:3, 3]
    return world.astype(np.float32), f.astype(np.uint32)


def build_meshes(atlas_path, level, template_path=None, cache_dir=None):
    """
    Surface mesh (world mm coordinates) of every label of one atlas level, plus a
    brain shell. Cached on disk (key: file, modification time, level, version).
    """
    import gzip, json, hashlib
    import nibabel as nib
    from scipy import ndimage
    atlas_path = Path(atlas_path)
    key_src = f'{atlas_path}|{atlas_path.stat().st_mtime}|{level}|{template_path}|{MESH_VERSION}'
    key = hashlib.sha1(key_src.encode()).hexdigest()[:20]
    if cache_dir:
        cache = Path(cache_dir) / f'mesh_{key}.json.gz'
        if cache.exists():
            with gzip.open(cache, 'rt') as fh:
                return json.load(fh)

    img = nib.load(str(atlas_path))
    data = img.dataobj[..., level] if len(img.shape) == 4 else img.dataobj
    lab = np.rint(np.asarray(data)).astype(np.int32)
    aff = img.affine
    vox_mm3 = float(abs(np.linalg.det(aff[:3, :3])))
    step = 1 if max(lab.shape) <= 200 else (2 if max(lab.shape) <= 400 else 3)

    ids, counts = np.unique(lab[lab > 0], return_counts=True)
    objs = ndimage.find_objects(lab)
    regions = []
    for lid, cnt in zip(ids, counts):
        sl = objs[lid - 1]
        if sl is None:
            continue
        sl = tuple(slice(max(s.start - 2, 0), min(s.stop + 2, n)) for s, n in zip(sl, lab.shape))
        sub = (lab[sl] == lid).astype(np.float32)
        centroid_vox = np.array(ndimage.center_of_mass(sub)) + [s.start for s in sl]
        sub = ndimage.gaussian_filter(sub, 0.6)
        try:
            v, f = _surface(sub, 0.5 * sub.max() if sub.max() < 1 else 0.5,
                            step if cnt > 200 else 1, np.array([s.start for s in sl]), aff)
        except (ValueError, RuntimeError):
            continue
        if len(f) == 0:
            continue
        regions.append(dict(id=int(lid), voxels=int(cnt), mm3=round(float(cnt) * vox_mm3, 4),
                            centroid=[round(float(c), 3) for c in (aff[:3, :3] @ centroid_vox + aff[:3, 3])],
                            nv=int(len(v)), nf=int(len(f)), v=_b64(v), f=_b64(f)))

    shell = None
    try:
        if template_path and Path(template_path).exists():
            t = nib.load(str(template_path))
            td = np.asarray(t.dataobj, dtype=np.float32)
            td = td[..., 0] if td.ndim == 4 else td
            pos = td[td > 0]
            if pos.size:
                from skimage.filters import threshold_otsu
                thr = threshold_otsu(pos[:: max(1, pos.size // 2_000_000)])   # tissue vs background
            else:
                thr = 0
            m = td > thr
            lab_cc, n_cc = ndimage.label(m)
            if n_cc > 1:   # keep the largest connected component (drops background specks)
                sizes = ndimage.sum(m, lab_cc, index=np.arange(1, n_cc + 1))
                m = lab_cc == (int(np.argmax(sizes)) + 1)
            m = ndimage.binary_fill_holes(m)
            taff = t.affine
        else:
            m = lab > 0
            taff = aff
        m = ndimage.gaussian_filter(m.astype(np.float32), 1.0)
        sstep = 2 if max(m.shape) <= 300 else 3
        v, f = _surface(m, 0.5, sstep, np.zeros(3), taff)
        shell = dict(nv=int(len(v)), nf=int(len(f)), v=_b64(v), f=_b64(f))
    except (ValueError, RuntimeError, MemoryError):
        shell = None

    out = dict(atlas=str(atlas_path), level=int(level), vox_mm3=vox_mm3, regions=regions, shell=shell)
    if cache_dir:
        Path(cache_dir).mkdir(parents=True, exist_ok=True)
        with gzip.open(cache, 'wt') as fh:
            json.dump(out, fh)
    return out


def atlas_levels(atlas_path):
    import nibabel as nib
    shp = nib.load(str(atlas_path)).shape
    return shp[3] if len(shp) == 4 else 1


# =============================================================================
# Background runs
# =============================================================================
class RunManager:
    """Launchers run detached (own session): they survive the GUI. State in <state_dir>/runs.json."""

    def __init__(self, state_dir):
        self.dir = Path(state_dir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.file = self.dir / 'runs.json'
        self.procs = {}
        self.runs = json.loads(self.file.read_text()) if self.file.exists() else []

    def _save(self):
        self.file.write_text(json.dumps(self.runs, indent=1))

    def start(self, launcher, repo, strict=True, python=None):
        import secrets
        rid = datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + '_' + secrets.token_hex(2)
        log = self.dir / f'run_{rid}_{Path(launcher).stem}.log'
        env = dict(os.environ)
        env['EDNIX_STRICT'] = '1' if strict else '0'
        env['PYTHONPATH'] = str(repo) + os.pathsep + env.get('PYTHONPATH', '')
        py = python or sys.executable
        fh = open(log, 'w')
        fh.write(f'[EDNiX GUI] {datetime.datetime.now().isoformat(timespec="seconds")}  {py} {launcher}\n')
        fh.flush()
        if os.name == 'posix':
            env['EDNIX_PY'], env['EDNIX_LAUNCHER'] = py, str(launcher)
            cmd = ['bash', '-c', '"$EDNIX_PY" -u "$EDNIX_LAUNCHER"; echo "[EDNIX_EXIT_CODE] $?"']
            p = subprocess.Popen(cmd, cwd=str(repo), env=env, stdout=fh, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, start_new_session=True)
        else:
            p = subprocess.Popen([py, '-u', str(launcher)], cwd=str(repo), env=env, stdout=fh,
                                 stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                 creationflags=getattr(subprocess, 'CREATE_NEW_PROCESS_GROUP', 0))
        self.procs[rid] = p
        self.runs.insert(0, dict(id=rid, launcher=str(launcher), log=str(log), pid=p.pid,
                                 started=datetime.datetime.now().isoformat(timespec='seconds'), strict=strict))
        self._save()
        return rid

    @staticmethod
    def _exit_code_from_log(log):
        try:
            with open(log, 'rb') as fh:
                fh.seek(0, 2)
                size = fh.tell()
                fh.seek(max(size - 4096, 0))
                tail = fh.read().decode('utf-8', errors='replace')
        except OSError:
            return None
        m = re.findall(r'\[EDNIX_EXIT_CODE\] (\d+)', tail)
        return int(m[-1]) if m else None

    def status(self, run):
        p = self.procs.get(run['id'])
        if p is not None:
            rc = p.poll()
            if rc is None:
                return 'running', None
        code = self._exit_code_from_log(run['log'])
        if code is not None:
            return ('finished' if code == 0 else 'failed'), code
        alive = False
        try:
            os.kill(run['pid'], 0)
            alive = Path(f'/proc/{run["pid"]}').exists() or os.name != 'posix'
        except (OSError, ProcessLookupError):
            alive = False
        return ('running', None) if alive else ('interrupted', None)

    STEP_RE = re.compile(r'Working on step\s*([\w.]+)\s*\(function:\s*([\w.]+)\)')

    def progress(self, run, nbytes=400_000):
        """Last pipeline step announced in the log ('##  Working on step N (function: X)  ##')."""
        try:
            with open(run['log'], 'rb') as fh:
                fh.seek(0, 2)
                size = fh.tell()
                fh.seek(max(size - nbytes, 0))
                tail = fh.read().decode('utf-8', errors='replace')
        except OSError:
            return None
        hits = self.STEP_RE.findall(tail)
        if not hits:
            return None
        step, fn = hits[-1]
        return dict(step=step, function=fn.lstrip('_'), steps_seen=len(hits))

    def listing(self):
        out = []
        for r in self.runs[:100]:
            st, code = self.status(r)
            out.append(dict(r, status=st, code=code, name=Path(r['launcher']).name, progress=self.progress(r)))
        return out

    def tail(self, rid, nbytes=60000):
        run = next((r for r in self.runs if r['id'] == rid), None)
        if run is None:
            return ''
        try:
            with open(run['log'], 'rb') as fh:
                fh.seek(0, 2)
                size = fh.tell()
                fh.seek(max(size - nbytes, 0))
                return fh.read().decode('utf-8', errors='replace')
        except OSError:
            return ''

    def stop(self, rid):
        import signal
        run = next((r for r in self.runs if r['id'] == rid), None)
        if run is None:
            return False
        try:
            if os.name == 'posix':
                os.killpg(run['pid'], signal.SIGTERM)
            else:
                subprocess.run(['taskkill', '/PID', str(run['pid']), '/T', '/F'], capture_output=True)
            return True
        except (OSError, ProcessLookupError):
            return False


# =============================================================================
# Label tables (self-contained: no dependency on the EDNiX repository modules)
# =============================================================================
def _hemi(name, lid):
    n = name.lower()
    if n.startswith(('l_', 'left', 'lh_', 'lh.', 'ctx-lh', 'ctx_lh')) or n.endswith(('_l', '_left')):
        return 'left'
    if n.startswith(('r_', 'right', 'rh_', 'rh.', 'ctx-rh', 'ctx_rh')) or n.endswith(('_r', '_right')):
        return 'right'
    return 'bilateral'


def parse_label_text(path):
    """
    {id: {name, rgb, hemi}} from a Workbench label file (name line, then 'id r g b a'),
    a FreeSurfer LUT / StatsLUT ('id name r g b a') or an ITK-SNAP label file
    ('id r g b a vis mesh "name"').
    """
    lines = [l.rstrip('\n').lstrip('\ufeff') for l in open(path, encoding='utf-8', errors='replace')]
    lines = [l for l in lines if l.strip() and not l.lstrip().startswith('#')]
    out = {}
    itk = re.compile(r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+[\d.]+\s+\d+\s+\d+\s+"(.*)"\s*$')
    lut = re.compile(r'^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)(?:\s+\d+)?\s*$')
    wb_num = re.compile(r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)(?:\s+\d+)?\s*$')
    if any(itk.match(l) for l in lines[:20]):
        for l in lines:
            m = itk.match(l)
            if m and int(m.group(1)) != 0:
                out[int(m.group(1))] = dict(name=m.group(5), rgb=[int(m.group(2)), int(m.group(3)), int(m.group(4))])
    elif sum(bool(wb_num.match(l)) for l in lines[:40]) >= max(1, len(lines[:40]) // 3):
        for i in range(len(lines) - 1):
            m = wb_num.match(lines[i + 1])
            if not wb_num.match(lines[i]) and m and int(m.group(1)) != 0:
                out[int(m.group(1))] = dict(name=lines[i].strip(), rgb=[int(m.group(2)), int(m.group(3)), int(m.group(4))])
    else:
        for l in lines:
            m = lut.match(l)
            if m and int(m.group(1)) != 0:
                out[int(m.group(1))] = dict(name=m.group(2), rgb=[int(m.group(3)), int(m.group(4)), int(m.group(5))])
    for lid, d in out.items():
        d['hemi'] = _hemi(d['name'], lid)
    return out


def label_file_for(species_dir, atlas_name):
    lc = Path(species_dir) / 'label_code'
    for fn in (f'{atlas_name}_label.txt', f'{atlas_name}_StatsLUT.txt', f'{atlas_name}_ITK.txt'):
        if (lc / fn).is_file():
            return lc / fn
    return None


def load_labels(species_dir, atlas_name):
    f = label_file_for(species_dir, atlas_name)
    if f is None:
        raise FileNotFoundError(f'no {atlas_name}_label.txt, _StatsLUT.txt or _ITK.txt in {Path(species_dir) / "label_code"}')
    return parse_label_text(f), f


# =============================================================================
# Cortical surfaces (Workbench native-resolution GIFTI)
# =============================================================================
SURFACE_KINDS = ('pial', 'midthickness', 'white', 'inflated', 'very_inflated')


def list_surfaces(ref_dir, species):
    d = Path(ref_dir) / 'surfaces' / 'Native_resol'
    if not d.is_dir():
        return dict(surfaces=[], atlases=[])
    kinds = [k for k in SURFACE_KINDS if (d / f'{species}.l.{k}.surf.gii').is_file() and (d / f'{species}.r.{k}.surf.gii').is_file()]
    atl = sorted({m.group(1) for f in d.glob(f'{species}.l.*.native.label.gii')
                  for m in [re.match(rf'{re.escape(species)}\.l\.(.+)\.native\.label\.gii$', f.name)] if m})
    return dict(surfaces=kinds, atlases=atl)


def build_surface(ref_dir, species, kind, atlas, species_dir=None):
    import nibabel as nib
    d = Path(ref_dir) / 'surfaces' / 'Native_resol'
    hemis, labels = [], {}
    for h in ('l', 'r'):
        surf = nib.load(str(d / f'{species}.{h}.{kind}.surf.gii'))
        coords = np.asarray(surf.agg_data('NIFTI_INTENT_POINTSET'), dtype=np.float32)
        faces = np.asarray(surf.agg_data('NIFTI_INTENT_TRIANGLE'), dtype=np.uint32)
        lab_img = nib.load(str(d / f'{species}.{h}.{atlas}.native.label.gii'))
        keys = np.asarray(lab_img.darrays[0].data).astype(np.int32).ravel()
        for key, lab in lab_img.labeltable.get_labels_as_dict().items():
            if key == 0 or not lab or lab.lower() in ('???', 'unknown', 'medial.wall', 'medial_wall'):
                continue
            L = next((x for x in lab_img.labeltable.labels if x.key == key), None)
            rgb = [int(round(c * 255)) for c in (L.red, L.green, L.blue)] if L is not None and L.red is not None else None
            labels.setdefault(int(key), dict(name=lab, rgb=rgb, hemi=_hemi(lab, key)))
        hemis.append(dict(hemi=h, nv=int(len(coords)), nf=int(len(faces)), v=_b64(coords), f=_b64(faces), lab=_b64(keys)))
    if species_dir is not None:   # names/colours from label_code when the label.gii table is empty
        try:
            txt, _ = load_labels(species_dir, atlas)
            for k, v in txt.items():
                labels.setdefault(k, v)
        except FileNotFoundError:
            pass
    present = set()
    for hm in hemis:
        import base64
        present.update(np.unique(np.frombuffer(base64.b64decode(hm['lab']), dtype=np.int32)).tolist())
    labels = {k: v for k, v in labels.items() if k in present}
    return dict(hemis=hemis, labels={str(k): v for k, v in labels.items()})


# =============================================================================
# Launcher documentation: docstrings of preprocess_* + curated notes + steps
# =============================================================================
STEPS = {
    'preprocess_anat': [
        (1, 'Orientation fix, averaging of repeated scans, bias-field correction'),
        (2, 'Cropping, denoising, AC-PC alignment and first skull-strip'),
        ('itk_1', 'Manual check of the first skull-strip in ITK-SNAP'),
        (5, 'Individual brain template (second skull-strip)'),
        ('itk_2', 'Manual check of the second skull-strip in ITK-SNAP'),
        (6, 'Registration to the study / reference template'),
        (8, 'Atlases and masks brought back to native space'),
        ('itk_3', 'Manual check of the native-space atlases in ITK-SNAP'),
        (9, 'FreeSurfer preparation (masks, T2w to T1w)'),
        (10, 'FreeSurfer input images (rescaling of small brains)'),
        (11, 'FreeSurfer white surface'),
        (12, 'FreeSurfer pial surface'),
        ('flat_map', 'Flat maps'),
        ('fs_extras', 'FreeSurfer parcellations Destrieux / DKT and FreeSurfer statistics (humans only; not used by EDNiX)'),
        (14, 'Atlases projected on the surfaces'),
        (15, 'FreeSurfer to Workbench conversion'),
        (16, 'Anatomical quality control (SNR, CNR...)'),
        ('Clean', 'Delete intermediate files'),
    ],
    'preprocess_fmri': [
        (1, 'Native-space preprocessing: slice timing, despiking, motion correction, FD censoring'),
        (2, 'Run-to-reference coregistration; distortion field from the field maps (applied at step 5)'),
        (3, 'Functional brain mask'),
        ('itk_1', 'Manual check of the functional masks in ITK-SNAP'),
        (5, 'Anatomy to functional registration; one-shot resampling of every volume (motion, distortion, registration)'),
        ('itk_2', 'Manual check of the functional data in ITK-SNAP'),
        (6, 'ICA (MELODIC) cleaning, interactive component selection'),
        (7, 'Nuisance regression, temporal filtering, smoothing'),
        (8, 'Functional data to anatomical space'),
        (9, 'Functional data to template space'),
        (10, 'Correlation matrices between atlas regions'),
        (11, 'Seed-based analysis'),
        (12, 'Functional quality control'),
        (13, 'Connectivity quality control (FC specificity)'),
        ('itk_3', 'Manual check of the final data in ITK-SNAP'),
        ('Clean', 'Delete intermediate files'),
    ],
    'preprocess_pet': [
        (1, 'Frame timing, SIF and mean image of the raw runs; reorientation of the header (as the anatomy)'),
        (2, 'PET study template (creat_study_template_pet)'),
        (3, 'Crop drawn on the study template in freeview (box used by the motion correction)'),
        (4, 'Motion correction within each dynamic run (two passes, rigid)'),
        (5, 'Motion correction between the runs of the session; mean image of the session'),
        (6, 'Registration to the anatomy / CT / template (scenarios of pet4.svg); atlases and masks in the PET space'),
        (7, 'Parametric images (SUV, SUVglu, ratio to the brain mean, time blocks) and regional values / TACs'),
        (8, 'Reference-region kinetic models: SRTM (bfmsrtm, imgbfbp, imgsrtm), Logan DVR (imgdv)'),
        ('pvc', 'Partial-volume correction (PETPVC, needs pvc_fwhm)'),
        (9, 'PET during the fMRI runs (simultaneous PET-MR, do_pet_to_func)'),
        (10, 'Export to the AC-PC and template spaces (+ skull-stripped), quality control'),
    ],
    'preprocess_ct': [
        (1, 'Reorientation of the header (as the PET of the same PET/CT acquisition)'),
        (2, 'CT study template (creat_study_template_ct, no anatomy)'),
        (3, 'Registration to the anatomy / CT template / study template, record for the PET pipeline, QC'),
    ],
}

CURATED = {
    'Skip_step': 'Pipeline steps NOT to run. Ticked steps are skipped (useful to resume a run or to skip manual ITK-SNAP checks).',
    'IhaveanANAT': 'True if every subject has an anatomical image; if False the functional data are registered directly to a template.',
    'aff_metric_ants': "Similarity metric of the functional-to-anatomical affine registration (ANTs), e.g. 'meansquares', 'mattes' (mutual information).",
    'aff_metric_ants_Transl': 'Similarity metric of the initial translation-only registration step (ANTs).',
    'blur': 'Spatial smoothing (FWHM, mm) applied during denoising; 0 disables it.',
    'dilate_mask': 'Number of voxels by which the functional brain mask is dilated.',
    'doWARPonfunc': "How oblique functional images are de-obliqued: 'header' edits the header only, 'WARP' resamples the data (3dWarp).",
    'n_for_ANTS': "Interpolation used by ANTs to resample functional images, e.g. 'lanczosWindowedSinc', 'linear'.",
    'type_of_transform': "ANTs transform of the functional-to-anatomical registration, e.g. 'BOLDAffine', 'SyN'.",
    'selected_atlases_SBA': "Seed atlases for the seed-based analysis: list of [atlas, level], e.g. [['EDNIxCSCLR', 3]]. Easiest from the Atlas tab.",
    'panda_files_SBA': 'Label IDs of the seed regions, one list per entry of selected_atlases_SBA, e.g. [[162, 1162]]. Easiest from the Atlas tab.',
    'selected_atlases_matrix': "Atlases and levels of the correlation matrices: 'all' or [[atlas names], [levels]].",
    'wanted_level_matrix': "Atlas level used for the matrices: 'all' or a level number.",
    'use_erode_V_func_masks': 'Use the eroded ventricle mask for nuisance regression (less grey-matter partial volume).',
    'use_erode_WM_func_masks': 'Use the eroded white-matter mask for nuisance regression (less grey-matter partial volume).',
    'endjson': 'Suffix of the JSON sidecar files.',
    'folderforTemplate_Anat': 'Folder of the anatomical template used when no individual anatomy exists (IhaveanANAT = False).',
    'BIDStype': 'Layout identifier of the BIDS dataset (see Tools/load_bids.py).',
    'selected_atlases': "Atlas / level pairs used to label the seed-based results, e.g. [['EDNIxCSCLR', 3]].",
    'dry_run': 'True only lists the files that would be deleted; False deletes them.',
    'matrix_fit': "Connectivity estimators (nilearn ConnectivityMeasure), e.g. 'correlation', 'partial correlation'.",
    'registration_fast': 'Not a parameter of preprocess_fmri: it is ignored. Use registration_type instead.',
    'MAIN_PATH': 'Folder of your EDNiX installation (the one that contains Atlases_library).',
    'motion_method': "Program that estimates the head motion: 'ants' (default), '3dvolreg' (AFNI) or 'mcflirt' (FSL); the last two are converted to ANTs matrices and checked.",
    'motion_mask': "Brain mask for the motion estimation: 'None' (whole image), 'previous' (the functional mask of an earlier execution: run once, then again from step 1) or a path.",
    'motion_mask_erode': 'Voxels removed from the border of the motion mask: at least the largest head displacement in voxels (default 3).',
    'extract_exterior_CSF': 'Not available: must be False (no step produces its mask, and the signals come from skull-stripped data).',
    'extract_WM': "Add the white-matter signal as a regressor. Used by the 'AFNI' method only ('Grandjean' and 'fMRIPrep_style' always use white matter).",
    'extract_Vc': "Add the ventricle signal as a regressor. Used by the 'AFNI' method only ('Grandjean' and 'fMRIPrep_style' always use the ventricles).",
    'acompcor_mode': "'fixed': 5 components per tissue; 'parallel': the components above chance; 'adaptive': variance criterion, usually limited by the cap. 'fMRIPrep_style' only.",
    'acompcor_n_fixed': "Number of aCompCor components per tissue in 'fixed' mode.",
    'acompcor_variance': "Variance to explain in 'adaptive' mode.",
    'acompcor_max': "Cap on the number of aCompCor components per tissue ('adaptive' and 'parallel' modes).",
    'doBALSAspace': 'Human, chimpanzee and macaque: also bring the subject to the BALSA template space (S1200, CY29, MY19), resample its surfaces to fs_LR, and for the macaque build the BALSA atlas in native space.',
    'extra_erode': 'Number of voxels by which the edge of the functional brain mask is eroded, removing dark edge voxels of the mean image; 0 disables it.',
    'list_to_keep': "Only these (subject, session) pairs are processed, e.g. [('01', '1')]. Empty list = all.",
    'list_to_remove': "These (subject, session) pairs are excluded, e.g. [('03', '2')].",
}


def _parse_docstring(doc):
    out = {}
    for m in re.finditer(r':param\s+(\w+):\s*(.*?)(?=\n\s*:(?:param|type|return)|\Z)', doc, re.S):
        out[m.group(1)] = ' '.join(m.group(2).split())
    lines, i = doc.split('\n'), 0
    while i < len(lines):
        m = re.match(r'^(\w[\w, ]*?)\s*:\s*(.*)$', lines[i])
        if m and not lines[i].startswith((' ', ':')):
            names, desc = [n.strip() for n in m.group(1).split(',')], []
            i += 1
            while i < len(lines) and (lines[i].startswith(' ') or not lines[i].strip()):
                if lines[i].strip():
                    desc.append(lines[i].strip())
                i += 1
            for n in names:
                if desc:
                    out.setdefault(n, ' '.join(desc))
            continue
        i += 1
    return out


def pipeline_docs(repo):
    docs = {}
    if not repo:
        return docs
    for mod in Path(repo).glob('modalities/*/_0_Pipeline_launcher.py'):
        try:
            tree = ast.parse(mod.read_text(encoding='utf-8'))
        except (SyntaxError, OSError):
            continue
        for n in tree.body:
            if isinstance(n, ast.FunctionDef) and n.name in PIPELINE_FUNCS:
                docs[n.name] = _parse_docstring(ast.get_docstring(n) or '')
    return docs


# =============================================================================
# Quality-control files written by the pipeline (<bids_dir>/QC)
# =============================================================================
IMG_EXT = ('.png', '.jpg', '.jpeg', '.gif', '.svg')
TAB_EXT = ('.csv', '.tsv', '.xlsx', '.xls')


# Per-subject reports of the pipeline (written next to the data, not in <BIDS>/QC):
#   sub-*/ses-*/func/QC/      <run>_QC_values.json, <run>_QC_summary.png (step 12), <run>_qc_report.png,
#                             all_runs_classification.csv (step 13)
#   sub-*/ses-*/anat/QC_anat/ <name>_QC_values.json, summary images (step 16)
#   sub-*/ses-*/pet/QC/       PET in the AC-PC / template spaces (step 10 of the PET pipeline), motion
#   sub-*/ses-*/ct/QC/        CT in the AC-PC / template spaces
#   sub-*/ses-*/{anat,func,pet,ct}/<ID> session <ses><date>.txt   processing logs (diary)
REPORT_DIRS = (('Functional reports', 'func', 'QC'), ('Anatomical reports', 'anat', 'QC_anat'),
               ('PET reports', 'pet', 'QC'), ('CT reports', 'ct', 'QC'))
VIRTUAL_TABLES = {'@func_metrics': 'Functional QC metrics (all runs)',
                  '@anat_metrics': 'Anatomical QC metrics (all subjects)',
                  '@specificity': 'Connectivity specificity (all runs)'}
# Pipeline messages are upper case (run_cmd: 'ERROR: ...', 'WARNING: ...'); FreeSurfer / AFNI write 'error:' /
# 'warning:' at the start of a line; Python ends a traceback with 'SomethingError: ...'. Case-sensitive on purpose:
# 'error' inside a sentence of a normal output is not an error.
ERROR_RE = re.compile(r'\bERROR\b|Traceback \(most recent call last\)|exited with code [1-9]|^\s*(?:\*\* )?error:|'
                      r'^\s*\w*(?:Error|Exception)\b(?::|$)|\bFATAL\b|Segmentation fault', re.MULTILINE)
WARN_RE = re.compile(r'\bWARNING\b|^\s*(?:\*\* )?warning:', re.MULTILINE)


def _subject_dirs(bids_dir, modality, sub_folder=''):
    root = Path(bids_dir)
    pats = [f'sub-*/ses-*/{modality}', f'sub-*/{modality}']
    out = []
    for pat in pats:
        for d in sorted(root.glob(pat)):
            d = d / sub_folder if sub_folder else d
            if d.is_dir():
                out.append(d)
    return out


def _ids(path, bids_dir):
    parts = Path(path).relative_to(bids_dir).parts
    sub = next((p[4:] for p in parts if p.startswith('sub-')), '')
    ses = next((p[4:] for p in parts if p.startswith('ses-')), '')
    return sub, ses


def qc_scan(bids_dir, max_images=5000):
    """Images and tables of the dataset QC folder, the per-subject reports, and the processing logs."""
    bids = Path(bids_dir)
    root = bids / 'QC'
    steps, tables, n = {}, [], 0
    if root.is_dir():
        for p in sorted(root.rglob('*')):
            if not p.is_file():
                continue
            rel = p.relative_to(bids).as_posix()
            inner = p.relative_to(root).as_posix()
            if p.suffix.lower() in IMG_EXT and n < max_images:
                steps.setdefault(inner.split('/')[0] if '/' in inner else '(top level)', []).append(rel)
                n += 1
            elif p.suffix.lower() in TAB_EXT:
                tables.append(dict(rel=rel, name=inner))
    found = dict(func_json=False, anat_json=False, classif=False)
    for label, modality, folder in REPORT_DIRS:
        for d in _subject_dirs(bids, modality, folder):
            for p in sorted(d.iterdir()):
                if p.suffix.lower() in IMG_EXT and n < max_images:
                    steps.setdefault(label, []).append(p.relative_to(bids).as_posix())
                    n += 1
                elif p.name.endswith('_QC_values.json'):
                    found['func_json' if modality == 'func' else 'anat_json'] = True
                elif p.name == 'all_runs_classification.csv':
                    found['classif'] = True
    virtual = [k for k, flag in (('@func_metrics', found['func_json']), ('@anat_metrics', found['anat_json']),
                                  ('@specificity', found['classif'])) if flag]
    tables = [dict(rel=k, name=VIRTUAL_TABLES[k]) for k in virtual] + tables
    logs = _log_files(bids)
    return dict(root=str(root), exists=root.is_dir() or bool(steps) or bool(logs), tables=tables,
                truncated=n >= max_images, steps=[dict(name=k, images=v) for k, v in steps.items()],
                n_logs=len(logs))


def qc_path(bids_dir, rel):
    """Resolve a file of the dataset and refuse anything outside the BIDS folder."""
    root = Path(bids_dir).resolve()
    p = (root / rel).resolve()
    if root not in p.parents:
        raise ValueError('file outside the BIDS folder')
    if not p.is_file():
        raise FileNotFoundError(rel)
    return p


def _flatten(d, prefix=''):
    out = {}
    for k, v in d.items():
        key = f'{prefix}{k}'
        if isinstance(v, dict):
            out.update(_flatten(v, key + '.'))
        elif isinstance(v, bool):
            out[key] = int(v)
        elif isinstance(v, (int, float)) or v is None:
            out[key] = v
        elif isinstance(v, str) and len(v) < 60:
            out[key] = v
    return out


def _virtual_table(bids_dir, rel):
    import pandas as pd
    bids, rows = Path(bids_dir), []
    if rel in ('@func_metrics', '@anat_metrics'):
        modality, folder = ('func', 'QC') if rel == '@func_metrics' else ('anat', 'QC_anat')
        for d in _subject_dirs(bids, modality, folder):
            for f in sorted(d.glob('*_QC_values.json')):
                try:
                    values = json.loads(f.read_text())
                except (OSError, ValueError):
                    continue
                sub, ses = _ids(f, bids)
                rows.append(dict(subject=sub, session=ses, run=f.name[:-len('_QC_values.json')], **_flatten(values)))
        return pd.DataFrame(rows)
    if rel == '@specificity':
        frames = []
        for d in _subject_dirs(bids, 'func', 'QC'):
            f = d / 'all_runs_classification.csv'
            if f.is_file():
                df = pd.read_csv(f)
                sub, ses = _ids(f, bids)
                df.insert(0, 'session', ses)
                df.insert(0, 'subject', sub)
                frames.append(df)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    raise ValueError(f'unknown table {rel}')


def qc_table(bids_dir, rel, max_rows=3000, z_thresh=3.5):
    """Table with robust z-scores (median / MAD) of numeric columns; |z| > z_thresh flags an outlier."""
    import pandas as pd
    if rel.startswith('@'):
        df = _virtual_table(bids_dir, rel)
    else:
        p = qc_path(bids_dir, rel)
        if p.suffix.lower() in ('.xlsx', '.xls'):
            df = pd.read_excel(p)
        else:
            df = pd.read_csv(p, sep='\t' if p.suffix.lower() == '.tsv' else None, engine='python')
    df = df.head(max_rows)
    flags = pd.DataFrame(False, index=df.index, columns=df.columns)
    for c in df.columns:
        if c in ('subject', 'session', 'run'):
            continue
        v = pd.to_numeric(df[c], errors='coerce')
        if v.notna().sum() < 5:
            continue
        med = v.median()
        mad = (v - med).abs().median() * 1.4826
        if mad > 0:
            flags[c] = ((v - med).abs() / mad) > z_thresh
    cells = df.astype(object).where(df.notna(), None).values.tolist()
    cells = [[(round(x, 4) if isinstance(x, float) else (x if isinstance(x, (int, str)) or x is None else str(x)))
              for x in row] for row in cells]
    return dict(columns=[str(c) for c in df.columns], rows=cells, flags=flags.values.tolist(),
                n_flagged_rows=int(flags.any(axis=1).sum()), z_thresh=z_thresh,
                link={'@func_metrics': 'Functional reports', '@specificity': 'Functional reports',
                      '@anat_metrics': 'Anatomical reports'}.get(rel))


def _log_files(bids):
    out = []
    for modality in ('anat', 'func', 'pet', 'ct'):
        for d in _subject_dirs(bids, modality):
            out += [p for p in sorted(d.glob('* session *.txt')) if 'launcher_parameters' not in p.name]
    return out


def classify_line(line):
    if ERROR_RE.search(line):
        return 'error'
    if WARN_RE.search(line):
        return 'warning'
    return ''


def qc_logs(bids_dir):
    """Errors and warnings of every processing log (diary) of the dataset, most errors first."""
    bids, out = Path(bids_dir), []
    for f in _log_files(bids):
        try:
            lines = f.read_text(errors='replace').splitlines()
        except OSError:
            continue
        kinds = [classify_line(l) for l in lines]
        errs = [l for l, k in zip(lines, kinds) if k == 'error']
        sub, ses = _ids(f, bids)
        out.append(dict(rel=f.relative_to(bids).as_posix(), subject=sub, session=ses,
                        modality=f.parent.name, errors=len(errs), warnings=kinds.count('warning'),
                        last_error=errs[-1][:300] if errs else '', size=len(lines)))
    out.sort(key=lambda r: (-r['errors'], -r['warnings'], r['rel']))
    return out


def qc_logtext(bids_dir, rel, only_issues=False, context=2, max_lines=20000):
    """Lines of a log with their kind ('error', 'warning', ''); with only_issues, errors / warnings and their context."""
    p = qc_path(bids_dir, rel)
    if not p.name.endswith('.txt'):
        raise ValueError('not a log file')
    lines = p.read_text(errors='replace').splitlines()[-max_lines:]
    kinds = [classify_line(l) for l in lines]
    keep = range(len(lines))
    if only_issues:
        idx = set()
        for i, k in enumerate(kinds):
            if k:
                idx.update(range(max(0, i - context), min(len(lines), i + context + 1)))
        keep = sorted(idx)
    return dict(lines=[[i + 1, lines[i], kinds[i]] for i in keep], errors=kinds.count('error'), warnings=kinds.count('warning'))
