import subprocess
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

def getreal(humanPosition, animalPosition, orient):

    '''
    from https://dicom.innolitics.com/ciods/ct-image/general-series/00185100:
    patient position should be :
    HFP:  Head First-Prone
    HFS:  Head First-Supine
    HFDR: Head First-Decubitus Right
    HFDL: Head First-Decubitus Left
    FFDR: Feet First-Decubitus Right
    FFDL: Feet First-Decubitus Left
    FFP:  Feet First-Prone
    FFS:  Feet First-Supine
    LFP:  Left First-Prone
    LFS:  Left First-Supine
    RFP:  Right First-Prone
    RFS:  Right First-Supine
    AFDR: Anterior First-Decubitus Right
    AFDL: Anterior First-Decubitus Left
    PFDR: Posterior First-Decubitus Right
    PFDL: Posterior First-Decubitus Left

    Hence following the same logic : "AHF" stands for "Animal Head First",  "AFF" stands for "Animal Feet First",
    "humanlike" means no change to be done.
    '''

    OPPOSITES = {
        ('AHF', 'HFS'): {'R': 'L', 'L': 'R', 'A': 'S', 'P': 'I', 'S': 'A', 'I': 'P'},
        ('AHF', 'FFS'): {'R': 'R', 'L': 'L', 'A': 'S', 'P': 'I', 'S': 'P', 'I': 'A'},
        ('AHF', 'HFP'): {'R': 'R', 'L': 'L', 'A': 'I', 'P': 'S', 'S': 'A', 'I': 'P'},
        ('AHF', 'FFP'): {'R': 'L', 'L': 'R', 'A': 'I', 'P': 'S', 'S': 'P', 'I': 'A'},

        ('AFF', 'HFS'): {'R': 'R', 'L': 'L', 'A': 'S', 'P': 'I', 'S': 'P', 'I': 'A'},
        ('AFF', 'FFS'): {'R': 'L', 'L': 'R', 'A': 'S', 'P': 'I', 'S': 'A', 'I': 'P'},
        ('AFF', 'HFP'): {'R': 'L', 'L': 'R', 'A': 'I', 'P': 'S', 'S': 'P', 'I': 'A'},
        ('AFF', 'FFP'): {'R': 'R', 'L': 'L', 'A': 'I', 'P': 'S', 'S': 'A', 'I': 'P'},

        ('AFS', 'HFS'): {'R': 'R', 'L': 'L', 'A': 'I', 'P': 'S', 'S': 'A', 'I': 'P'},
        ('AFS', 'FFS'): {'R': 'L', 'L': 'R', 'A': 'I', 'P': 'S', 'S': 'A', 'I': 'P'},
        ('AFS', 'HFP'): {'R': 'L', 'L': 'R', 'A': 'S', 'P': 'I', 'S': 'A', 'I': 'P'},
        ('AFS', 'FFP'): {'R': 'R', 'L': 'L', 'A': 'S', 'P': 'I', 'S': 'A', 'I': 'P'},

        ('AHS', 'HFS'): {'R': 'L', 'L': 'R', 'A': 'I', 'P': 'S', 'S': 'P', 'I': 'A'},
        ('AHS', 'FFS'): {'R': 'R', 'L': 'L', 'A': 'I', 'P': 'S', 'S': 'A', 'I': 'P'},
        ('AHS', 'HFP'): {'R': 'R', 'L': 'L', 'A': 'S', 'P': 'I', 'S': 'P', 'I': 'A'},
        ('AHS', 'FFP'): {'R': 'L', 'L': 'R', 'A': 'S', 'P': 'I', 'S': 'A', 'I': 'P'},
    }

    default = {'R': 'R', 'L': 'L', 'A': 'A', 'P': 'P', 'S': 'S', 'I': 'I'}

    opposite = OPPOSITES.get(humanPosition, animalPosition, default)

    neworient = ''.join(opposite[c] for c in orient)

    return neworient