#!/usr/bin/env python3
"""
fix_library_wb_surfaces.py: redo only the FreeSurfer -> Workbench conversion of library templates, keeping their
FreeSurfer surfaces (template.modif(..., wb_only=True)), then update the library manifest.

    python Exemples/Study_EDNiX/atlases/fix_library_wb_surfaces.py                       # Marmoset, Atlases_library
    python Exemples/Study_EDNiX/atlases/fix_library_wb_surfaces.py --species Marmoset Rat --library <path>

Why: small brains (smallbrain.get) were brought back from FreeSurfer with mri_surf2surf --reg-inv <rescale.dat>,
whose translation is a per-orientation guess ('pos' in smallbrain.get): the Marmoset V2 surfaces were shifted
relative to the volume. FS_2_WB_short now applies the exact inverse (scaling by 1/s, no translation). Minutes per
species instead of a full rebuild.
"""
import argparse
import datetime
import json
import os
import sys
from os.path import join as opj

HERE = os.path.dirname(os.path.abspath(__file__))
MAIN_PATH = os.path.dirname(os.path.dirname(HERE)) + '/'
sys.path.insert(0, MAIN_PATH)
sys.path.insert(0, HERE)

import build_library_V2 as L  # noqa: E402


def fix(species, library):
    from atlases import template
    done = {}
    for sp, ref, tdir in L.find_templates(library, species, L.REFERENCES):
        diary = opj(tdir, ref, 'fix_wb_surfaces_diary.txt')
        open(diary, 'a').close()
        print(f'[fix] {sp}/{ref}: FreeSurfer -> Workbench conversion ({tdir})', flush=True)
        template.modif(L.SURFACE['doBALSAspace'], L.SMALLBRAIN, sp, ref, L.ATLAS_FOLLOWERS, L.MAIN_PATH, tdir, diary,
                       proj=L.SURFACE['proj'], doflat=L.SURFACE['doflat'], library_path=library,
                       smallbrain_factor=L.SURFACE['smallbrain_factor'], sphere=L.SURFACE['sphere'], wb_only=True)
        try:
            L.species_qc(sp, ref, tdir, opj(library, 'QC_library'))
        except Exception as e:
            print(f'[fix] {sp}: surface QC failed: {type(e).__name__}: {e}')
        done[f'{sp}/{ref}'] = tdir
    return done


def update_manifest(library, done):
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
    man.setdefault('fixes', []).append(dict(date=datetime.datetime.now().isoformat(timespec='seconds'),
                                            what='FreeSurfer -> Workbench conversion redone (exact 1/s unscale)',
                                            species=sorted(done), ednix=L.git_info(L.MAIN_PATH)))
    with open(path, 'w') as f:
        json.dump(man, f, indent=1)
    print(f'[fix] manifest updated: {path}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--species', nargs='+', default=['Marmoset'])
    ap.add_argument('--library', default=L.LIBRARY_V1, help='default: Atlases_library (the V2 once swapped)')
    a = ap.parse_args()
    done = fix(a.species, a.library)
    if not done:
        sys.exit(f'no template found for {a.species} in {a.library}')
    update_manifest(a.library, done)
