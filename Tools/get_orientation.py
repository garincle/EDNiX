import subprocess
import numpy as np
import ants
import re
import nibabel as nib

spgo = subprocess.getoutput

# purpose for FS : LIA
FS ={'R':'L','L':'R',
     'A':'P','P':'A',
     'S':'I','I':'S'}


def fromITK_to_newFS(orig,desired):
    transfo = []

    for index,item in enumerate(orig):
        if FS[item] in desired:
            for i, j in enumerate(desired):
                if FS[item]==j:
                    transfo.append(i+1)
        else:
            for i, j in enumerate(desired):
                if item==j:
                    transfo.append((i+1)*-1)
    transfo = ' '.join(str(x) for x in transfo)
    return transfo

#reorient_img = fromITK_to_newFS(orient_itk,'LIA')
'''
def use_FS(img: object, sing_fs: object) -> object:

    cmd = (sing_fs + 'mri_info --orientation ' + img)
    orient_img = spgo(cmd).split('\n')[-1]

    orient_raw = ''
    reorient   = ''

    for index, orient in enumerate(orient_FS):
        if orient_img   == orient:
            orient_raw = orient_FS[index]
            reorient = ' -r ' + reorient_img[index]+ ' '

    fwdFS_cmd = ' --in_orientation ' + orient_raw + reorient
    bckFS_cmd = ' --in_orientation LIA' + reorient

    return [orient_raw,reorient,fwdFS_cmd,bckFS_cmd]
'''
def use_ants(img):
    img_hd = ants.image_read(img)
    orient_img = img_hd.orientation
    orientFS = ''.join([FS[i] for i in orient_img])
    reorient = fromITK_to_newFS(orient_img, 'LIA')

    fwdFS_cmd = ' --in_orientation ' + orientFS + ' -r ' +  reorient + ' '
    bckFS_cmd = ' --in_orientation LIA -r ' + reorient + ' '

    return [orient_img,reorient,fwdFS_cmd,bckFS_cmd]


def use_afni(img,sing_afni):
    cmd = sing_afni + '3dinfo -orient ' + img
    orient_img = spgo(cmd).split('\n')[-1]
    orientFS = ''.join([FS[i] for i in orient_img])
    reorient = fromITK_to_newFS(orient_img, 'LIA')

    fwdFS_cmd = ' --in_orientation ' + orientFS + ' -r ' + reorient + ' '
    bckFS_cmd = ' --in_orientation LIA ' + reorient + ' '

    return [orient_img, reorient, fwdFS_cmd, bckFS_cmd]


def get_orientation_nibabel(nifti_path):
    """Get 3-letter orientation code using NiBabel."""
    img = nib.load(nifti_path)
    aff = img.affine

    # Extract orientation from affine matrix
    ornt = nib.orientations.io_orientation(aff)
    codes = nib.orientations.ornt2axcodes(ornt)
    orient_code = ''.join(codes)

    # Validate (should already be valid from NiBabel)
    if not re.fullmatch(r'^[RLAPSI]{3}$', orient_code):
        raise ValueError(f"Invalid orientation: {orient_code}")
    return orient_code

# Directions (scanner x, y, z; z = into the bore, y = up) of the anatomical
# Anterior and Superior axes.
# Human positions: labels written by the scanner (DICOM PatientPosition).
_HUMAN_AS = {'HFS': ((0, 1, 0), (0, 0, 1)),
             'FFS': ((0, 1, 0), (0, 0, -1)),
             'HFP': ((0, -1, 0), (0, 0, 1)),
             'FFP': ((0, -1, 0), (0, 0, -1))}
# Animal positions, with the semantics implied by the existing lookup tables
# (all physically consistent entries of the previous version are reproduced):
#   AHF: head first, prone (sphinx)      AFF: feet first, prone
#   AFS: head first, supine              AHS: feet first, supine
# NB: the AFS / AHS names are counter-intuitive but kept for backward compatibility.
_ANIMAL_AS = {'AHF': ((0, 0, 1), (0, 1, 0)),
              'AFF': ((0, 0, -1), (0, 1, 0)),
              'AFS': ((0, 0, 1), (0, -1, 0)),
              'AHS': ((0, 0, -1), (0, -1, 0))}


def _frame(a_dir, s_dir):
    a = np.asarray(a_dir, float)
    s = np.asarray(s_dir, float)
    r = np.cross(a, s)                      # RAS right-handed: R = A x S
    return {'R': r, 'L': -r, 'A': a, 'P': -a, 'S': s, 'I': -s}


def getreal(humanPosition, animalPosition, orient):
    """
    Translate the orientation code written by the scanner (which assumes a human
    in `humanPosition`) into the true anatomical orientation of an animal placed
    in `animalPosition`.

    The mapping is derived geometrically (scanner label direction -> animal
    anatomical axis with the same direction), so it is always a proper rotation.
    The previous lookup tables contained five reflections (AFF/FFS, AFF/FFP,
    AFS/FFS, AFS/FFP, AHS/FFP) that silently mirrored left and right.

    humanPosition : 'HFS', 'FFS', 'HFP' or 'FFP'
    animalPosition: 'AHF', 'AFF', 'AFS', 'AHS' (see _ANIMAL_AS) or 'humanlike'
    orient        : 3-letter orientation code, e.g. 'RAI'
    """
    if animalPosition == 'humanlike':
        return orient
    if humanPosition not in _HUMAN_AS:
        raise ValueError(f"humanPosition '{humanPosition}' not supported with animalPosition "
                         f"'{animalPosition}' (expected one of {sorted(_HUMAN_AS)}).")
    if animalPosition not in _ANIMAL_AS:
        raise ValueError(f"animalPosition '{animalPosition}' not recognised "
                         f"(expected one of {sorted(_ANIMAL_AS)} or 'humanlike').")
    scanner = _frame(*_HUMAN_AS[humanPosition])
    animal = _frame(*_ANIMAL_AS[animalPosition])
    mapping = {lab: next(k for k, v in animal.items() if np.allclose(v, vec))
               for lab, vec in scanner.items()}
    return ''.join(mapping[c] for c in orient)
