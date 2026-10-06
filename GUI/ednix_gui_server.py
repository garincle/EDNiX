#!/usr/bin/env python3
"""
EDNiX desktop GUI - local server.

Serves the web interface (web/) on 127.0.0.1 and exposes a small JSON API
to EDNiX: launcher editing, background runs, atlas library and 3D meshes.
A random token (passed in the URL fragment, never in the query) protects the
API against other users of a shared machine.

Usage: python ednix_gui_server.py [--port 8765] [--no-browser]
"""
import argparse, json, os, sys, secrets, threading, webbrowser, socket, mimetypes, traceback
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ednix_gui_core as C  # noqa: E402

WEB = HERE / 'web'
STATE = Path(os.environ.get('EDNIX_GUI_STATE', Path.home() / '.ednix_gui'))
TOKEN = secrets.token_urlsafe(18)


def is_ednix(p):
    return (p / 'modalities').is_dir() and (p / 'Tools').is_dir() and (p / 'Statistics').is_dir()


def find_repo():
    cands = []
    if os.environ.get('EDNIX_REPO'):
        cands.append(Path(os.environ['EDNIX_REPO']).expanduser())
    cands += [HERE.parent, HERE.parent.parent, Path.cwd(), Path.home() / 'PycharmProjects' / 'EDNiX', Path.home() / 'EDNiX']
    return next((p.resolve() for p in cands if is_ednix(p)), None)


REPO = find_repo()
RUNS = C.RunManager(STATE)
OPTIONS = C.discover_options(REPO) if REPO else dict(C.STATIC_OPTIONS)
ATLAS_LIB = Path(os.environ.get('EDNIX_ATLAS_LIBRARY', (REPO / 'Atlases_library') if REPO else ''))
if (ATLAS_LIB / 'atlas').is_dir():
    OPTIONS['species'] = sorted(C.list_species(ATLAS_LIB))
MODELS = {}
HTTPD = None
DOCS = C.pipeline_docs(REPO)
PIPE_RE = C.re.compile(r'([\w.]*)\.(preprocess_anat|preprocess_fmri|preprocess_pet|preprocess_data)\s*\(')


def pipeline_of(path):
    """anat / func / pet from the pipeline function the launcher calls, else from its folder name."""
    try:
        m = PIPE_RE.search(Path(path).read_text(encoding='utf-8', errors='replace'))
    except OSError:
        m = None
    if m:
        fn, mod = m.group(2), m.group(1)
        if fn == 'preprocess_data':
            return 'pet' if 'pet' in mod else 'other'
        return {'preprocess_anat': 'anat', 'preprocess_fmri': 'func', 'preprocess_pet': 'pet'}[fn]
    folder = Path(path).parent.name.lower()
    for key, pipe in (('anat', 'anat'), ('func', 'func'), ('pet', 'pet')):
        if key in folder:
            return pipe
    return 'other'


def check_python():
    missing = []
    for mod in ('numpy', 'nibabel', 'scipy', 'skimage', 'pandas', 'ants', 'nilearn'):
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    return missing


class Handler(BaseHTTPRequestHandler):
    server_version = 'EDNiXGUI/1.0'

    def log_message(self, fmt, *args):
        pass

    # ---- helpers
    def _json(self, obj, code=200):
        body = json.dumps(obj).encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def _err(self, msg, code=400):
        self._json({'error': msg}, code)

    def _authorized(self, query=None):
        tok = self.headers.get('X-EDNiX-Token', '') or (query or {}).get('token', '')
        return secrets.compare_digest(tok, TOKEN)

    def _body(self):
        n = int(self.headers.get('Content-Length', '0') or 0)
        return json.loads(self.rfile.read(n).decode('utf-8')) if n else {}

    def _static(self, rel):
        rel = rel.lstrip('/') or 'index.html'
        target = (WEB / rel).resolve()
        if WEB.resolve() not in target.parents and target != WEB.resolve() / 'index.html':
            return self._err('not found', 404)
        if not target.is_file():
            return self._err('not found', 404)
        data = target.read_bytes()
        ctype = mimetypes.guess_type(str(target))[0] or 'application/octet-stream'
        if target.suffix == '.js':
            ctype = 'text/javascript'
        self.send_response(200)
        self.send_header('Content-Type', ctype + ('; charset=utf-8' if ctype.startswith('text') else ''))
        self.send_header('Content-Length', str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    # ---- routing
    def do_GET(self):
        u = urlparse(self.path)
        if not u.path.startswith('/api/'):
            return self._static(u.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        # images are loaded by <img> tags, which cannot send headers: the key may come in the query
        if not self._authorized(q if u.path == '/api/qc/file' else None):
            return self._err('missing or invalid access key: reopen the interface with the EDNiX launcher', 403)
        try:
            return self._api_get(u.path, q)
        except Exception as e:
            traceback.print_exc()
            return self._err(f'{type(e).__name__}: {e}', 500)

    def do_POST(self):
        u = urlparse(self.path)
        if not self._authorized():
            return self._err('missing or invalid access key', 403)
        try:
            return self._api_post(u.path, self._body())
        except Exception as e:
            traceback.print_exc()
            return self._err(f'{type(e).__name__}: {e}', 500)

    # ---- API
    def _species_dir(self, name):
        sp = C.list_species(ATLAS_LIB)
        if name not in sp:
            raise ValueError(f'unknown species: {name}')
        return sp[name]

    def _api_get(self, path, q):
        if path == '/api/info':
            launchers = []
            if REPO:
                for p in sorted((REPO / 'Exemples' / 'Study_EDNiX').glob('Launcher_*/*.py')):
                    launchers.append(dict(path=str(p), name=p.name, group=p.parent.name.replace('Launcher_', ''),
                                          pipeline=pipeline_of(p)))
            return self._json(dict(repo=str(REPO) if REPO else None, python=sys.executable,
                                   python_version=sys.version.split()[0], missing=check_python(),
                                   atlas_library=str(ATLAS_LIB), atlas_library_ok=(ATLAS_LIB / 'atlas').is_dir(),
                                   launchers=launchers, state_dir=str(STATE)))
        if path == '/api/atlas/species':
            return self._json(sorted(C.list_species(ATLAS_LIB)) if (ATLAS_LIB / 'atlas').is_dir() else [])
        if path == '/api/atlas/references':
            return self._json(C.list_references(self._species_dir(q['species'])))
        if path == '/api/atlas/list':
            ref_dir = self._species_dir(q['species']) / q['reference']
            out = []
            sp_dir = self._species_dir(q['species'])
            for key, (name, p) in C.list_atlases(ref_dir).items():
                out.append(dict(key=key, name=name, file=p.name, levels=C.atlas_levels(p),
                                has_labels=C.label_file_for(sp_dir, name) is not None))
            out.sort(key=lambda a: (not a['has_labels'], a['key'].lower()))
            return self._json(dict(atlases=out, templates=C.list_templates(ref_dir),
                                   surfaces=C.list_surfaces(ref_dir, q['species'])))
        if path == '/api/atlas/mesh':
            sp_dir = self._species_dir(q['species'])
            ref_dir = sp_dir / q['reference']
            atlases = C.list_atlases(ref_dir)
            name, p = atlases[q['atlas']]
            tpl = ref_dir / 'volumes' / q['template'] if q.get('template') else None
            mesh = C.build_meshes(p, int(q.get('level', 0)), tpl, STATE / 'cache')
            labels = {}
            try:
                labels, lfile = C.load_labels(sp_dir, name)
                mesh['label_file'] = str(lfile)
            except Exception as e:
                mesh['label_warning'] = f'label names not found: {e}'
            for reg in mesh['regions']:
                info = labels.get(reg['id'], {})
                reg['name'] = info.get('name', f'label {reg["id"]}')
                reg['hemi'] = info.get('hemi', 'bilateral')
                reg['rgb'] = info.get('rgb')
            missing = sum(1 for r in mesh['regions'] if r['id'] not in labels)
            if labels and missing:
                mesh['label_warning'] = f'{missing} label IDs of the image are not in {mesh.get("label_file")}'
            mesh['atlas_name'] = name
            return self._json(mesh)
        if path == '/api/atlas/surface':
            sp_dir = self._species_dir(q['species'])
            ref_dir = sp_dir / q['reference']
            out = C.build_surface(ref_dir, q['species'], q['surface'], q['atlas'], sp_dir)
            out['atlas_name'] = q['atlas']
            return self._json(out)
        if path == '/api/launcher':
            m = C.LauncherModel(q['path'])
            MODELS[str(m.path)] = m
            return self._json(C.launcher_to_json(m, OPTIONS, DOCS))
        if path == '/api/runs':
            return self._json(RUNS.listing())
        if path == '/api/qc/candidates':
            out = []
            for p in sorted((REPO / 'Exemples' / 'Study_EDNiX').glob('Launcher_*/*.py')) if REPO else []:
                try:
                    m = C.LauncherModel(p)
                except Exception:
                    continue
                for f in m.fields:
                    if f.name == 'bids_dir' and isinstance(f.value, str) and Path(f.value).is_dir():
                        out.append(dict(path=f.value.rstrip('/'), launcher=p.name,
                                        has_qc=(Path(f.value) / 'QC').is_dir()))
            seen, uniq = set(), []
            for o in out:
                if o['path'] not in seen:
                    seen.add(o['path'])
                    uniq.append(o)
            return self._json(uniq)
        if path == '/api/qc/scan':
            return self._json(C.qc_scan(q['bids']))
        if path == '/api/qc/table':
            return self._json(C.qc_table(q['bids'], q['rel']))
        if path == '/api/qc/logs':
            return self._json(C.qc_logs(q['bids']))
        if path == '/api/qc/logtext':
            return self._json(C.qc_logtext(q['bids'], q['rel'], only_issues=q.get('issues') == '1'))
        if path == '/api/qc/file':
            f = C.qc_path(q['bids'], q['rel'])
            if f.suffix.lower() not in C.IMG_EXT:
                return self._err('only images are served', 403)
            data = f.read_bytes()
            self.send_response(200)
            self.send_header('Content-Type', mimetypes.guess_type(str(f))[0] or 'application/octet-stream')
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'private, max-age=300')
            self.end_headers()
            self.wfile.write(data)
            return None
        if path == '/api/run/log':
            return self._json(dict(text=RUNS.tail(q['id'])))
        return self._err('unknown route', 404)

    def _api_post(self, path, b):
        if path in ('/api/launcher/check', '/api/launcher/save'):
            m = MODELS.get(b['path']) or C.LauncherModel(b['path'])
            vals, errors = C.values_from_json(m, b.get('values', {}))
            if errors:
                return self._json(dict(ok=False, errors=errors))
            src, n = m.render(vals)
            try:
                compile(src, b.get('out') or b['path'], 'exec')
            except SyntaxError as e:
                return self._json(dict(ok=False, errors=[f'syntax error: {e}']))
            issues = C.LauncherModel.check_calls(src, REPO) + m.warnings
            if path.endswith('check'):
                return self._json(dict(ok=True, changed=n, issues=issues))
            out = Path(b['out']).expanduser()
            if not out.suffix == '.py':
                return self._json(dict(ok=False, errors=['the file name must end with .py']))
            backup = None
            if out.exists():
                backup = out.with_name(out.name + '.bak_' + C.datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
                backup.write_bytes(out.read_bytes())
            out.write_text(src, encoding='utf-8')
            return self._json(dict(ok=True, changed=n, issues=issues, written=str(out),
                                   backup=str(backup) if backup else None))
        if path == '/api/run':
            launcher = Path(b['path']).expanduser()
            if not launcher.is_file():
                return self._err(f'launcher not found: {launcher}')
            rid = RUNS.start(launcher, REPO, strict=bool(b.get('strict', True)),
                             python=os.environ.get('EDNIX_PYTHON') or sys.executable)
            return self._json(dict(id=rid))
        if path == '/api/run/stop':
            return self._json(dict(ok=RUNS.stop(b['id'])))
        if path == '/api/shutdown':
            self._json(dict(ok=True))
            threading.Thread(target=HTTPD.shutdown, daemon=True).start()
            return None
        return self._err('unknown route', 404)


def free_port(preferred):
    for port in [preferred] + list(range(preferred + 1, preferred + 50)):
        with socket.socket() as s:
            try:
                s.bind(('127.0.0.1', port))
                return port
            except OSError:
                continue
    raise RuntimeError('no free port')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=8765)
    ap.add_argument('--no-browser', action='store_true')
    args = ap.parse_args()
    if REPO is None:
        print('EDNiX repository not found: put the GUI folder inside it or set EDNIX_REPO.')
    port = free_port(args.port)
    global HTTPD
    httpd = HTTPD = ThreadingHTTPServer(('127.0.0.1', port), Handler)
    url = f'http://127.0.0.1:{port}/#token={TOKEN}'
    STATE.mkdir(parents=True, exist_ok=True)
    (STATE / 'last_url.txt').write_text(url)
    print(f'EDNiX GUI: {url}\n(repository: {REPO})\nClose this window or press Ctrl+C to stop the interface; '
          f'runs already started keep going in the background.', flush=True)
    if not args.no_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
