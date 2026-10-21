"""
Transform chains and registration records of the PET and CT pipelines.

Chain
-----
A chain is a list of [file, invert] pairs in the order of ants.apply_transforms (transformlist,
whichtoinvert). A pair of chains {'to': ..., 'from': ...} links a modality space M (PET or CT,
'orig') with a target space T:
    'to'   resamples an image of M into T   (apply_transforms(fixed=T image, moving=M image))
    'from' resamples an image of T into M   (apply_transforms(fixed=M image, moving=T image))
Points of the fixed space go through the transforms in list order, so, for M -> A -> T,
    to(M->T)   = to(A->T) + to(M->A)
    from(T->M) = from(M->A) + from(A->T)

Registration record
-------------------
JSON file in orig/matrices, written by step 6 (PET) or step 3 (CT), read by the next steps:
    {"modality": "pet", "scenario": "1v2a", "reference": <session mean image>,
     "targets": {"acpc":  {"image": ..., "brain": ..., "mask": ..., "to": [...], "from": [...]},
                 "EDNiX": {...}, "studyTemplate": {...}}}
"""
import os
import json
import ants

opj = os.path.join
opi = os.path.isfile

ANTS_AFF = '0GenericAffine.mat'
ANTS_WARP = '1Warp.nii.gz'
ANTS_INVWARP = '1InverseWarp.nii.gz'


def affine(mat):
    """Chains of one affine matrix of an ANTs registration (fixed = target, moving = modality)."""
    return {'to': [[mat, False]], 'from': [[mat, True]]}


def registration(outprefix):
    """Chains of an ANTs registration written with `outprefix` (affine, + warp if non-linear)."""
    aff = outprefix + ANTS_AFF
    if not opi(aff):
        raise FileNotFoundError(f'registration matrix not found: {aff}')
    if opi(outprefix + ANTS_WARP):
        return {'to': [[outprefix + ANTS_WARP, False], [aff, False]],
                'from': [[aff, True], [outprefix + ANTS_INVWARP, False]]}
    return affine(aff)


def from_backtonative(package):
    """Chains native (AC-PC) <-> template of one package of modalities.anat.loop3.backtonative.get."""
    return {'to': [[f, bool(w)] for f, w in zip(package[9], package[10])],
            'from': [[f, bool(w)] for f, w in zip(package[7], package[8])]}


def compose(first, then):
    """Chains M -> T from M -> A (first) and A -> T (then)."""
    return {'to': then['to'] + first['to'], 'from': first['from'] + then['from']}


def check(chain):
    missing = [f for f, _ in chain['to'] + chain['from'] if not opi(f)]
    if missing:
        raise FileNotFoundError('transform file(s) not found: ' + ', '.join(missing))
    return chain


def apply(fixed, moving, chain, direction, interpolator='linear', imagetype=None, out=None):
    """Resample `moving` (path or ANTsImage) on `fixed` with chain[direction] ('to' or 'from')."""
    fx = ants.image_read(fixed) if isinstance(fixed, str) else fixed
    mv = ants.image_read(moving) if isinstance(moving, str) else moving
    if imagetype is None:
        imagetype = 3 if mv.dimension == 4 else 0
    tl = [f for f, _ in chain[direction]]
    inv = [bool(w) for _, w in chain[direction]]
    res = ants.apply_transforms(fixed=fx, moving=mv, transformlist=tl, whichtoinvert=inv,
                                interpolator=interpolator, imagetype=imagetype)
    if out:
        ants.image_write(res, out, ri=False)
    return res


def write_record(path, record):
    with open(path, 'w') as f:
        json.dump(record, f, indent=2)
    return path


def read_record(path):
    if not opi(path):
        raise FileNotFoundError(f'registration record not found: {path} (run the registration step first)')
    with open(path) as f:
        return json.load(f)
