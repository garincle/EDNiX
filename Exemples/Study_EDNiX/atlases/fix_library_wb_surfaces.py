#!/usr/bin/env python3
"""
fix_library_wb_surfaces.py: repair the surfaces of some library templates without rebuilding the whole library.

1. Redo only the FreeSurfer -> Workbench conversion (keeps the FreeSurfer surfaces; minutes per species):
    python Exemples/Study_EDNiX/atlases/fix_library_wb_surfaces.py                         # Marmoset, V2 library
    python Exemples/Study_EDNiX/atlases/fix_library_wb_surfaces.py --species Marmoset Rat --library <path>
   Why: small brains were brought back from FreeSurfer with mri_surf2surf --reg-inv <rescale.dat>, whose translation
   is a per-orientation guess ('pos' in smallbrain.get): the Marmoset V2 surfaces were shifted from the volume.
   FS_2_WB_short now applies the exact inverse (scaling by 1/s, no translation).

2. Full rebuild of a species (FreeSurfer + Workbench) with other smallbrain settings, in a SANDBOX copy (the library
   is not touched): template folder copied, the rest of the library symlinked.
    python Exemples/Study_EDNiX/atlases/fix_library_wb_surfaces.py --rebuild --factor 20 --volume-ref 300000 \\
        --sandbox /scratch2/EDNiX/surface_tests/marmoset_grid
   -> <sandbox>/QC_library/<Species>_EDNiX_surface_QC.png to compare with <library>/QC_library/.
   Without --sandbox, --rebuild replaces the species in --library (and updates its manifest).
"""
import argparse
import datetime
import json
import os
import shutil
import sys
from os.path import join as opj

HERE = os.path.dirname(os.path.abspath(__file__))
MAIN_PATH = os.path.dirname(os.path.dirname(HERE)) + '/'
sys.path.insert(0, MAIN_PATH)
sys.path.insert(0, HERE)

import build_library_V2 as L  # noqa: E402


def fix(species, library):
    """FreeSurfer -> Workbench conversion only. Returns {species/ref: template_dir}."""
    from atlases import template
    done = {}
    for sp, ref, tdir in L.find_templates(library, species, L.REFERENCES):
        diary = opj(tdir, ref, 'fix_wb_surfaces_diary.txt')
        open(diary, 'a').close()
        print(f'[fix] {sp}/{ref}: FreeSurfer -> Workbench conversion ({tdir})', flush=True)
        template.modif(L.SURFACE['doBALSAspace'], L.SMALLBRAIN, sp, ref, L.ATLAS_FOLLOWERS, L.MAIN_PATH, tdir, diary,
                       proj=L.SURFACE['proj'], doflat=L.SURFACE['doflat'], library_path=library,
                       smallbrain_factor=L.SURFACE['smallbrain_factor'], sphere=L.SURFACE['sphere'], wb_only=True)
        _qc(sp, ref, tdir, library)
        done[f'{sp}/{ref}'] = tdir
    return done


def make_sandbox(library, species, sandbox):
    """<sandbox>: copy of the template folders of `species`, symlinks to everything else of the library."""
    os.makedirs(sandbox, exist_ok=True)
    for name in os.listdir(library):
        if name in ('atlas', 'QC_library', 'EDNiX_library_manifest.json'):
            continue
        dst = opj(sandbox, name)
        if not os.path.lexists(dst):
            os.symlink(opj(library, name), dst)
    for sp, ref, tdir in L.find_templates(library, species, L.REFERENCES):
        dst = opj(sandbox, os.path.relpath(tdir, library))
        if not os.path.isdir(dst):
            print(f'[sandbox] copying {tdir} -> {dst}', flush=True)
            shutil.copytree(tdir, dst, symlinks=True)
    return sandbox


def rebuild(species, library, factor=None, volume_ref=None):
    """Full FreeSurfer + Workbench rebuild (build_library_V2.build_species) with optional smallbrain settings."""
    done = {}
    for sp, ref, tdir in L.find_templates(library, species, L.REFERENCES):
        sb = dict(L.SMALLBRAIN_SETTINGS.get(sp, {}))
        if factor:
            sb['factor'] = factor
        if volume_ref:
            sb['volume_ref'] = volume_ref
        L.SMALLBRAIN_SETTINGS[sp] = sb
        print(f'[rebuild] {sp}/{ref} in {library}: smallbrain {sb or "defaults"}', flush=True)
        minutes = L.build_species(sp, ref, tdir, library, L.SURFACE) / 60
        print(f'[rebuild] {sp}/{ref}: {minutes:.0f} min', flush=True)
        _qc(sp, ref, tdir, library)
        done[f'{sp}/{ref}'] = tdir
    return done


def _qc(sp, ref, tdir, library):
    try:
        qc = L.species_qc(sp, ref, tdir, opj(library, 'QC_library'))
        print(f'[qc] {sp}: {qc}')
    except Exception as e:
        print(f'[qc] {sp}: surface QC failed: {type(e).__name__}: {e}')


def update_manifest(library, done, what):
    path = opj(library, 'EDNiX_library_manifest.json')
    if not os.path.isfile(path):
        print(f'[fix] no manifest in {library}: nothing to update')
        return
    with open(path) as f:
        man = json.load(f)
    if 'files' in man:
        files = {k: v for k, v in man['files'].items()
                 if not any(k.startswith(os.path.relpath(t, library) + os.sep) for t in done.values())}
        for t in done.values():
            files.update({os.path.relpath(opj(t, k), library): v for k, v in L.checksums(t).items()})
        man['files'] = dict(sorted(files.items()))
    man.setdefault('fixes', []).append(dict(date=datetime.datetime.now().isoformat(timespec='seconds'), what=what,
                                            species=sorted(done), ednix=L.git_info(L.MAIN_PATH)))
    with open(path, 'w') as f:
        json.dump(man, f, indent=1)
    print(f'[fix] manifest updated: {path}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--species', nargs='+', default=['Marmoset'])
    ap.add_argument('--library', default=L.LIBRARY_V2 if os.path.isdir(L.LIBRARY_V2) else L.LIBRARY_V1,
                    help='default: Atlases_library_V2 while it exists (before the swap), else Atlases_library')
    ap.add_argument('--rebuild', action='store_true', help='full FreeSurfer + Workbench rebuild of the species')
    ap.add_argument('--factor', type=int, help='--rebuild: smallbrain factor (library default 30)')
    ap.add_argument('--volume-ref', type=float, help='--rebuild: smallbrain volume_ref in mm3 (default 107000)')
    ap.add_argument('--sandbox', help='--rebuild: work in this copy, the library is not touched')
    a = ap.parse_args()
    if not a.rebuild:
        done = fix(a.species, a.library)
        what = 'FreeSurfer -> Workbench conversion redone (exact 1/s unscale)'
        target = a.library
    else:
        target = make_sandbox(a.library, a.species, a.sandbox) if a.sandbox else a.library
        done = rebuild(a.species, target, a.factor, a.volume_ref)
        what = f'species rebuilt (smallbrain factor={a.factor}, volume_ref={a.volume_ref})'
    if not done:
        sys.exit(f'no template found for {a.species} in {target}')
    if not a.sandbox:
        update_manifest(target, done, what)
    print(f'[done] QC: {opj(target, "QC_library")}')
