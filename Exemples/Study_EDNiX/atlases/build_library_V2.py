"""
Build a reproducible V2 of the EDNiX atlas library (Atlases_library -> Atlases_library_V2).

What it does
  1. copies the V1 library to a versioned folder (V1 is never modified: everything below runs in the copy),
  2. rebuilds the template surfaces of every species (atlases/template.modif) with the parameters of SURFACE,
     in parallel (one process per species, N_JOBS at a time),
  3. writes the surface QC of each species (Tools/QC_report.surface_qc: white/pial/mid on 3 planes, thickness in real
     mm, aseg fidelity) to <V2>/QC_library/,
  4. installs the Workbench scenes of the paper Figure 1 into <Species>/EDNiX/scenes/ (install_figure_scenes.py,
     paths relative to the library),
  5. writes <V2>/EDNiX_library_manifest.json: EDNiX git commit (+ uncommitted changes flag), package versions, the
     exact parameters, per-species status / duration / QC metrics, and the sha256 of every file of the library.

`python build_library_V2.py --verify` re-hashes the library and reports any file that differs from the manifest.
`python build_library_V2.py --plan` lists the species found and what would be done, without running anything.

Surface method: 'V2h', chosen from the surface tests of PR #84 (5 library templates + 10 BIDS subjects): synthetic
brain.mgz with non-overlapping classes (CSF 10-35, GM 55-85 with the template texture, WM 110, sigma 0.5 voxel),
hippocampus/amygdala in the GM class, pial search capped at 3 real mm (modalities/anat/freesurfer/preFS.py). The
smallbrain rescale (factor 30) is kept and now also applies to the Marmoset template. 'legacy' rebuilds as in V1.
"""
import os
import sys
import json
import time
import shutil
import hashlib
import argparse
import datetime
import platform
import subprocess
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

opj = os.path.join
opd = os.path.dirname
opi = os.path.isfile

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, MAIN_PATH)

# ---------------------------------------------------------------------------------------------------------------------
# settings
# ---------------------------------------------------------------------------------------------------------------------
VERSION = 'V2'
LIBRARY_V1 = opj(opd(MAIN_PATH), 'Atlases_library')          # same rule as atlases/template.modif
LIBRARY_V2 = opj(opd(MAIN_PATH), 'Atlases_library_' + VERSION)
SPECIES = 'all'                                             # 'all' (every template found) or e.g. ['Rat', 'Marmoset']
REFERENCES = ['EDNiX']                                       # reference folder(s) of each species to rebuild
SMALLBRAIN = ['Mouse', 'Rat', 'Mouselemur', 'Bat', 'Marmoset']  # species rescaled for FreeSurfer (smallbrain.get);
#   Marmoset added in V2 (template < 10 cm3, like the subjects): without it FreeSurfer ran in high-resolution mode
# per-species smallbrain.get settings, e.g. {'Marmoset': dict(factor=20, volume_ref=300000)}; default: factor
# SURFACE['smallbrain_factor'] and smallbrain.volume_ref (107 cm3). The FreeSurfer grid is ~factor/10 times coarser than
# the template (x3 with 30, x2 with 20); with factor 20, volume_ref must be >= 1/8 of the brain voxel count so that a
# FreeSurfer voxel stays >= 1 FS mm (no high-resolution mode). Test with fix_library_wb_surfaces.py --rebuild --sandbox.
SMALLBRAIN_SETTINGS = {}
ATLAS_FOLLOWERS = [['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]]
SURFACE = dict(method='V2h',           # preFS.SURFACE_METHOD: 'V2h' (surface tests, PR #84) or 'legacy' (EDNiX v1)
               smallbrain_factor=30,   # smallbrain.get factor for the SMALLBRAIN species (V1 value: 30)
               sphere=True,            # spherical registration (needed by the library for fs_LR / BALSA)
               doBALSAspace=True,
               doflat=0,               # flat maps (slow)
               proj='ribbon')
N_JOBS = 4                    # species rebuilt in parallel
THREADS_PER_JOB = 4           # ITK / OpenMP threads per species
CHECKSUMS = True              # sha256 of every file in the manifest (a few minutes for a full library)


# ---------------------------------------------------------------------------------------------------------------------
def find_templates(library, species='all', references=('EDNiX',)):
    """[(species, reference, template_dir)] for every <template_dir>/<reference>/volumes/<sp>_space-acpc_desc-SS_T1w."""
    out = []
    for dirpath, dirnames, _ in os.walk(opj(library, 'atlas')):
        if os.path.basename(dirpath) == 'freesurfer':
            dirnames[:] = []
            continue
        sp = os.path.basename(dirpath)
        for ref in references:
            if opi(opj(dirpath, ref, 'volumes', f'{sp}_space-acpc_desc-SS_T1w.nii.gz')):
                if species == 'all' or sp in species:
                    out.append((sp, ref, dirpath))
    return sorted(out)


def git_info(path):
    def run(*a):
        try:
            return subprocess.run(['git', '-C', path] + list(a), capture_output=True, text=True,
                                  timeout=30).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ''
    return dict(commit=run('rev-parse', 'HEAD'), branch=run('rev-parse', '--abbrev-ref', 'HEAD'),
                uncommitted_changes=bool(run('status', '--porcelain', '--untracked-files=no')))


def versions():
    v = dict(python=platform.python_version(), platform=platform.platform())
    for mod in ('ants', 'numpy', 'nibabel', 'pandas'):
        try:
            v[mod] = __import__(mod).__version__
        except Exception:
            v[mod] = None
    return v


def sha256(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(block), b''):
            h.update(chunk)
    return h.hexdigest()


def checksums(root, skip=('EDNiX_library_manifest.json',)):
    out = {}
    for dirpath, _, files in os.walk(root):
        for f in files:
            p = opj(dirpath, f)
            rel = os.path.relpath(p, root)
            if rel in skip or os.path.islink(p):
                continue
            out[rel] = dict(size=os.path.getsize(p), sha256=sha256(p))
    return dict(sorted(out.items()))


def copy_library(src, dst):
    if os.path.exists(dst):
        print(f'INFO: {dst} exists: copy skipped (delete it to start again from V1).')
        return False
    if not os.path.isdir(src):
        raise FileNotFoundError(f'V1 library not found: {src}')
    print(f'INFO: copying {src} -> {dst}')
    shutil.copytree(src, dst, symlinks=True)
    return True


def build_species(sp, ref, template_dir, library, surface):
    """Rebuild the surfaces of one species template in `library`. Runs in a worker process."""
    os.environ.update({k: str(THREADS_PER_JOB) for k in ('ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS', 'OMP_NUM_THREADS',
                                                         'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')})
    os.environ.setdefault('MPLBACKEND', 'Agg')
    from atlases import template
    from modalities.anat.freesurfer import preFS
    if surface['method'] not in ('V2h', 'legacy'):
        raise NotImplementedError(f"surface method {surface['method']!r}: 'V2h' or 'legacy'")
    preFS.SURFACE_METHOD = surface['method']
    t0 = time.time()
    diary = opj(template_dir, ref, f'build_library_{VERSION}_diary.txt')
    open(diary, 'a').close()
    # start from empty surfaces, so no V1 file (e.g. an old sphere.reg) survives silently in V2
    shutil.rmtree(opj(template_dir, ref, 'freesurfer', sp), ignore_errors=True)
    from modalities.anat.freesurfer import smallbrain
    sb = SMALLBRAIN_SETTINGS.get(sp, {})
    smallbrain.volume_ref = sb.get('volume_ref', smallbrain.volume_ref)
    template.modif(surface['doBALSAspace'], SMALLBRAIN, sp, ref, ATLAS_FOLLOWERS, MAIN_PATH, template_dir, diary,
                   proj=surface['proj'], doflat=surface['doflat'], library_path=library,
                   smallbrain_factor=sb.get('factor', surface['smallbrain_factor']), sphere=surface['sphere'])
    return time.time() - t0


def species_qc(sp, ref, template_dir, qc_dir):
    from Tools import QC_report
    FS_dir = opj(template_dir, ref, 'freesurfer')
    k = QC_report._fs_scale(opj(template_dir, ref, 'preprocessing'), sp)
    return QC_report.surface_qc(FS_dir, sp, opj(qc_dir, f'{sp}_{ref}_surface_QC.png'), k, title=f'{sp} ({ref}) {VERSION}')


def build(species=SPECIES, plan=False):
    templates = find_templates(LIBRARY_V1 if not os.path.isdir(LIBRARY_V2) else LIBRARY_V2, species, REFERENCES)
    print(f'{len(templates)} template(s): ' + ', '.join(f'{s}/{r}' for s, r, _ in templates))
    if plan:
        print(f'V1: {LIBRARY_V1}\nV2: {LIBRARY_V2}\nsurface: {SURFACE}\njobs: {N_JOBS} x {THREADS_PER_JOB} threads')
        return None

    started = datetime.datetime.now().isoformat(timespec='seconds')
    copied = copy_library(LIBRARY_V1, LIBRARY_V2)
    templates = find_templates(LIBRARY_V2, species, REFERENCES)
    qc_dir = opj(LIBRARY_V2, 'QC_library')
    os.makedirs(qc_dir, exist_ok=True)

    results = {}
    with ProcessPoolExecutor(max_workers=N_JOBS) as ex:
        futs = {ex.submit(build_species, sp, ref, tdir, LIBRARY_V2, SURFACE): (sp, ref, tdir)
                for sp, ref, tdir in templates}
        for f in as_completed(futs):
            sp, ref, tdir = futs[f]
            key = f'{sp}/{ref}'
            try:
                results[key] = dict(status='ok', minutes=round(f.result() / 60, 1))
            except Exception as e:
                results[key] = dict(status='failed', error=f'{type(e).__name__}: {e}',
                                    traceback=traceback.format_exc()[-2000:])
            print(f'{key}: {results[key]["status"]}')
            if results[key]['status'] == 'ok':
                try:
                    results[key]['surface_QC'] = species_qc(sp, ref, tdir, qc_dir)
                except Exception as e:
                    results[key]['surface_QC'] = {'error': f'{type(e).__name__}: {e}'}

    # Figure 1 scenes (Workbench), paths rewritten relative to the library: before the checksums
    scenes = {}
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            'install_figure_scenes', opj(os.path.dirname(os.path.abspath(__file__)), 'install_figure_scenes.py'))
        ifs = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ifs)
        scenes = ifs.install(LIBRARY_V2, manifest=False)
    except Exception as e:
        print(f'WARNING: Figure 1 scenes not installed: {type(e).__name__}: {e}')

    manifest = dict(library=VERSION, built_from=LIBRARY_V1, copied_from_V1_in_this_run=copied,
                    started=started, finished=datetime.datetime.now().isoformat(timespec='seconds'),
                    ednix=git_info(MAIN_PATH), software=versions(),
                    settings=dict(SPECIES=species, REFERENCES=REFERENCES, SMALLBRAIN=SMALLBRAIN,
                                  SMALLBRAIN_SETTINGS=SMALLBRAIN_SETTINGS,
                                  ATLAS_FOLLOWERS=ATLAS_FOLLOWERS, SURFACE=SURFACE),
                    species=dict(sorted(results.items())),
                    figure_scenes={sp: len(v) for sp, v in sorted(scenes.items())})
    if CHECKSUMS:
        print('INFO: hashing the library (sha256)...')
        manifest['files'] = checksums(LIBRARY_V2)
    out = opj(LIBRARY_V2, 'EDNiX_library_manifest.json')
    with open(out, 'w') as f:
        json.dump(manifest, f, indent=1)
    n_fail = sum(r['status'] != 'ok' for r in results.values())
    print(f'\n{len(results) - n_fail}/{len(results)} species rebuilt. Manifest: {out}\nQC: {qc_dir}')
    if manifest['ednix']['uncommitted_changes']:
        print('WARNING: EDNiX had uncommitted changes: commit them for this library to be reproducible.')
    return manifest


def verify(library=LIBRARY_V2):
    """Compare the library with the checksums of its manifest. Returns the list of differences."""
    with open(opj(library, 'EDNiX_library_manifest.json')) as f:
        ref = json.load(f).get('files', {})
    now = checksums(library)
    diff = [f'missing: {k}' for k in ref if k not in now] + [f'added: {k}' for k in now if k not in ref] + \
           [f'changed: {k}' for k in ref if k in now and ref[k]['sha256'] != now[k]['sha256']]
    print('\n'.join(diff) if diff else f'OK: {len(ref)} files identical to the manifest.')
    return diff


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--species', nargs='+', help='only these species (default: SPECIES)')
    ap.add_argument('--plan', action='store_true', help='list what would be done, run nothing')
    ap.add_argument('--verify', action='store_true', help='check the V2 library against its manifest')
    a = ap.parse_args()
    if a.verify:
        sys.exit(1 if verify() else 0)
    build(a.species or SPECIES, plan=a.plan)
