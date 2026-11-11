#!/usr/bin/env python3
"""
install_figure_scenes.py: put the Workbench scenes of the paper Figure 1 (panel a) into the atlas library.

    python Exemples/Study_EDNiX/atlases/install_figure_scenes.py                          # into Atlases_library
    python Exemples/Study_EDNiX/atlases/install_figure_scenes.py --library .../Atlases_library_V2
    python Exemples/Study_EDNiX/atlases/install_figure_scenes.py --plan                   # report only

Source (SOURCE_DIR): <Species>.scene (black background), <Species>_1.scene (white background), optional
<Species>_thickness.scene, and the PNGs exported by hand from them (reference images). Names are case-insensitive
(dog.png, Dog_1.png, macaque.png...).

Destination, per species: <library>/atlas/mammals/<group>/<Species>/EDNiX/scenes/
    <Species>_Figure1.scene              <- <Species>.scene
    <Species>_Figure1_white.scene        <- <Species>_1.scene
    <Species>_Figure1_thickness.scene    <- <Species>_thickness.scene
    <Species>_Figure1[_white]_reference.png
    data/                                files a scene uses that are not in the library (copied)

Every file path of a scene (Object/ObjectArray of Type "pathName") is rewritten RELATIVE to the new scene file:
the library stays valid when it is renamed (rebuild_all_ednix swaps Atlases_library_V2 -> Atlases_library). A path
that points into any Atlases_library*/atlas/... is re-anchored on <library>/atlas/...; a path outside the library is
copied into scenes/data/; a missing file is looked up by its name in <Species>/EDNiX/. Paths that stay unresolved are
listed (the scene is still written, Workbench will report them).

build_library_V2.py calls install() before writing its manifest; run standalone, the files added here are appended to
the manifest of the library (sha256), so build_library_V2.py --verify stays clean.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from os.path import join as opj

SOURCE_DIR = '/home/cgarin/Documents/3_Articles/Figures/EDNiX'
LIBRARY = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library'
REFERENCE = 'EDNiX'
SPECIES = {'human': 'Human', 'chimpanzee': 'Chimpanzee', 'macaque': 'Macaque', 'marmoset': 'Marmoset',
           'mouselemur': 'Mouselemur', 'pig': 'Pig', 'rat': 'Rat', 'mouse': 'Mouse', 'dog': 'Dog', 'cat': 'Cat'}
VARIANTS = {'': 'Figure1', '_1': 'Figure1_white', '_thickness': 'Figure1_thickness',
            '_thickness_1': 'Figure1_thickness_white'}
MANIFEST = 'EDNiX_library_manifest.json'

_PATH_OBJECT = re.compile(r'(<Object\b[^>]*\bType="pathName"[^>]*>)([^<]*)(</Object>)')
_PATH_ARRAY = re.compile(r'(<ObjectArray\b[^>]*\bType="pathName"[^>]*>)(.*?)(</ObjectArray>)', re.S)
_ELEMENT = re.compile(r'(<Element\b[^>]*>)([^<]*)(</Element>)')
_LIBRARY_PART = re.compile(r'/Atlases_library[^/]*/(atlas/.*)$')


def species_dirs(library, reference=REFERENCE):
    """{species: <library>/atlas/.../<Species>/<reference>} for every species folder of the library."""
    out = {}
    for dirpath, dirnames, _ in os.walk(opj(library, 'atlas')):
        if os.path.basename(dirpath) in ('freesurfer', reference):
            dirnames[:] = []
            continue
        if os.path.isdir(opj(dirpath, reference)):
            out[os.path.basename(dirpath)] = opj(dirpath, reference)
    return out


def source_files(src=SOURCE_DIR):
    """{species: {variant: {'scene': path, 'png': path}}} from the file names of `src`."""
    out = {}
    for f in sorted(os.listdir(src)):
        stem, ext = os.path.splitext(f)
        if ext.lower() not in ('.scene', '.png'):
            continue
        m = re.fullmatch(r'([A-Za-z]+?)((?:_thickness)?(?:_1)?)', stem, re.I)
        if not m or m.group(1).lower() not in SPECIES:
            print(f'[skip] {f}: unknown species or variant')
            continue
        sp, var = SPECIES[m.group(1).lower()], m.group(2).lower()
        out.setdefault(sp, {}).setdefault(VARIANTS[var], {})['scene' if ext.lower() == '.scene' else 'png'] = opj(src, f)
    return out


def _find_by_name(name, root):
    for dirpath, _, files in os.walk(root):
        if name in files:
            return opj(dirpath, name)
    return None


def rewrite_paths(text, src_dir, dst_dir, library, species_dir, copy_data=True):
    """Scene text with every pathName rewritten relative to dst_dir. Returns (text, report rows)."""
    report = []

    def new_path(old):
        old_s = old.strip()
        if not old_s:
            return old
        full = os.path.normpath(old_s if os.path.isabs(old_s) else opj(src_dir, old_s))
        m = _LIBRARY_PART.search(full)
        if m:
            target, how = opj(library, m.group(1)), 'library'
            if not os.path.exists(target):
                alt = _find_by_name(os.path.basename(target), species_dir)
                target, how = (alt, 'library (found by name)') if alt else (target, 'MISSING in library')
        elif os.path.isfile(full):
            target, how = opj(dst_dir, 'data', os.path.basename(full)), 'copied to scenes/data'
            if copy_data:
                os.makedirs(os.path.dirname(target), exist_ok=True)
                shutil.copy2(full, target)
        else:
            alt = _find_by_name(os.path.basename(full), species_dir)
            if alt is None:
                report.append((old_s, old_s, 'UNRESOLVED (kept as is)'))
                return old
            target, how = alt, 'library (found by name)'
        rel = os.path.relpath(target, dst_dir)
        report.append((old_s, rel, how))
        return rel

    text = _PATH_OBJECT.sub(lambda m: m.group(1) + new_path(m.group(2)) + m.group(3), text)
    text = _PATH_ARRAY.sub(lambda a: a.group(1) + _ELEMENT.sub(
        lambda e: e.group(1) + new_path(e.group(2)) + e.group(3), a.group(2)) + a.group(3), text)
    return text, report


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def update_manifest(library, written):
    """Append the files written here to the library manifest (if any), so --verify stays clean."""
    man_path = opj(library, MANIFEST)
    if not os.path.isfile(man_path):
        return False
    with open(man_path) as f:
        man = json.load(f)
    if 'files' not in man:
        return False
    for p in written:
        man['files'][os.path.relpath(p, library)] = dict(size=os.path.getsize(p), sha256=_sha256(p))
    man['files'] = dict(sorted(man['files'].items()))
    man.setdefault('figure_scenes', {})['source'] = SOURCE_DIR
    with open(man_path, 'w') as f:
        json.dump(man, f, indent=1)
    return True


def install(library=LIBRARY, src=SOURCE_DIR, plan=False, manifest=True):
    """Copy the Figure 1 scenes + reference PNGs into the library. Returns {species: [written files]}."""
    if not os.path.isdir(src):
        print(f'[scenes] {src} not found: nothing installed')
        return {}
    dirs, files = species_dirs(library), source_files(src)
    written, unresolved = {}, 0
    for sp, variants in sorted(files.items()):
        if sp not in dirs:
            print(f'[scenes] {sp}: no {REFERENCE} folder in {library}: skipped')
            continue
        dst = opj(dirs[sp], 'scenes')
        print(f'\n[scenes] {sp} -> {dst}')
        for var, fs in sorted(variants.items()):
            if 'scene' in fs:
                with open(fs['scene'], encoding='utf-8') as f:
                    text, rep = rewrite_paths(f.read(), os.path.dirname(fs['scene']), dst, library, dirs[sp],
                                              copy_data=not plan)
                for old, new, how in rep:
                    unresolved += how.startswith(('UNRESOLVED', 'MISSING'))
                    print(f'   {os.path.basename(fs["scene"])}: {how:24s} {old} -> {new}')
                out = opj(dst, f'{sp}_{var}.scene')
                if not plan:
                    os.makedirs(dst, exist_ok=True)
                    with open(out, 'w', encoding='utf-8') as f:
                        f.write(text)
                    written.setdefault(sp, []).append(out)
                    written[sp] += [opj(dst, 'data', os.path.basename(n)) for _, n, h in rep if h.startswith('copied')]
            if 'png' in fs:
                out = opj(dst, f'{sp}_{var}_reference.png')
                if not plan:
                    os.makedirs(dst, exist_ok=True)
                    shutil.copy2(fs['png'], out)
                    written.setdefault(sp, []).append(out)
                print(f'   {os.path.basename(fs["png"])} -> {os.path.basename(out)}')
    print(f'\n[scenes] {sum(map(len, written.values()))} file(s) written for {len(written)} species; '
          f'{unresolved} unresolved path(s)' + (' (plan: nothing written)' if plan else ''))
    if written and manifest and update_manifest(library, sorted({p for v in written.values() for p in v})):
        print(f'[scenes] manifest updated: {opj(library, MANIFEST)}')
    return written


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--library', default=LIBRARY)
    ap.add_argument('--source', default=SOURCE_DIR)
    ap.add_argument('--plan', action='store_true', help='report what would be done, write nothing')
    a = ap.parse_args()
    install(a.library, a.source, plan=a.plan)
    sys.exit(0)
