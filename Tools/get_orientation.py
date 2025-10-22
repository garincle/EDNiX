import subprocess
import ants
import re
import nibabel as nib

spgo = subprocess.getoutput

# purpose for FS : LIA
orient_FS =  ['RAS','LAS','RSA','LSA',
              'RIA','LIA','RSP','LSP',
              'RPS','LPS','RPI','LPI',
              'RIP','LIP','RAI','LAI']

orient_itk = ['LPI','RPI','LIP','RIP',
              'LSP','RSP','LIA','RIA',
              'LAI','RAI','LAS','RAS',
              'LSA','RSA','LPS','RPS']

orient_ants = [[-1,-2,3], [1,-2,3], [-1,3,-2], [1,3,-2],
               [-1,-3,-2],[1,-3,-2],[-1,3,2],  [1,3,2],
               [-1,2,3],  [1,2,3],  [-1,2,-3], [1,2,-3],
               [-1,-3,2], [1,-3,2], [-1,-2,-3],[1,-2,-3]]

reorient_img = ['-1 3 -2',  '1 3 -2',  '-1 -2 3',  '1 -2 3',
                '-1 2 3',   '1 2 3',   '-1 -2 -3', '1 -2 -3',
                '-1 -3 -2', '1 -3 -2', '-1 -3 2',  '1 -3 2',
                '-1 2 -3',  '1 2 -3',  '-1 3 2',   '1 3 2']

def use_FS(img,sing_fs):

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

def use_ants(img,sing_fs):
    img_hd = ants.image_header_info(img)
    matrix= img_hd['direction']
    orient_img = [0,0,0]
    
    orient_raw = ''
    reorient   = ''

    for i in range(len(matrix)):
        for index, item in enumerate(matrix[i]):
            if item != 0:
                if item < 0:
                    orient_img[i]=(index+1)*-1
                else:
                    orient_img[i]=index+1
    
    for index, orient in enumerate(orient_ants):
        if orient_img == orient:
            orient_raw = orient_FS[index]
            reorient = ' -r ' + reorient_img[index] + ' '
            deobl=0
            break
        else :
            orient_raw,reorient,_,_ = use_FS(img,sing_fs)
            deobl=1

    fwdFS_cmd = ' --in_orientation ' + orient_raw + reorient
    bckFS_cmd = ' --in_orientation LIA' + reorient

    return [orient_raw,reorient,fwdFS_cmd,bckFS_cmd,deobl]


def use_afni(img,sing_afni):
    cmd = sing_afni + '3dinfo -orient ' + img
    orient_img = spgo(cmd).split('\n')[-1]

    for index, orient in enumerate(orient_itk):
        if orient_img == orient:
            orient_raw = orient_FS[index]
            reorient = ' -r ' + reorient_img[index] + ' '

    fwdFS_cmd = ' --in_orientation ' + orient_raw + reorient
    bckFS_cmd = ' --in_orientation LIA' + reorient

    return [orient_raw, reorient, fwdFS_cmd, bckFS_cmd]

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

    neworient = orient
    pos    = ['A','P','S','I']
    newpos = ['A','P','S','I']

    if animalPosition == 'AHF':
        if humanPosition == 'HFS':
            if orient[0] == 'L': neworient[0] = 'R'
            else: neworient[0] = 'L'
            newpos = ['S','I','A','P']

        elif humanPosition == 'FFS':
            neworient[0] = orient[0]
            newpos = ['S', 'I', 'P', 'A']

        elif humanPosition == 'HFP':
            neworient[0] = orient[0]
            newpos = ['I', 'S', 'A', 'P']

        elif humanPosition == 'FFP':
            if orient[0] == 'L': neworient[0] = 'R'
            else: neworient[0] = 'L'
            newpos = ['I', 'S', 'P', 'A']

        for z in [1, 2]:
            for i, j in enumerate(pos):
                if orient[z] == j:
                    neworient[z] = newpos[i]

    elif animalPosition == 'AFF':
        if humanPosition == 'HFS':
            neworient[0] = orient[0]
            newpos = ['S', 'I', 'P', 'A']

        elif humanPosition == 'FFS':
            if orient[0] == 'L': neworient[0] = 'R'
            else: neworient[0] = 'L'
            newpos = ['S', 'I', 'A', 'P']

        elif humanPosition == 'HFP':
            if orient[0] == 'L': neworient[0] = 'R'
            else: neworient[0] = 'L'
            newpos = ['I', 'S', 'P', 'A']

        elif humanPosition == 'FFP':
            neworient[0] = orient[0]
            newpos = ['I', 'S', 'A', 'P']

        for z in [1, 2]:
            for i, j in enumerate(pos):
                if orient[z] == j:
                    neworient[z] = newpos[i]


    elif animalPosition == 'AFS':
        if humanPosition == 'HFS':
            neworient[0] = orient[0]
            newpos = ['I', 'S', 'A', 'P']

        elif humanPosition == 'FFS':
            if orient[0] == 'L':neworient[0] = 'R'
            else:neworient[0] = 'L'
            newpos = ['I', 'S', 'P', 'A']

        elif humanPosition == 'HFP':
            if orient[0] == 'L':neworient[0] = 'R'
            else:neworient[0] = 'L'
            newpos = ['S', 'I', 'A', 'P']

        elif humanPosition == 'FFP':
            neworient[0] = orient[0]
            newpos = ['S', 'I', 'P', 'A']

        for z in [1, 2]:
            for i, j in enumerate(pos):
                if orient[z] == j:
                    neworient[z] = newpos[i]

    elif animalPosition == 'AHS':
        if humanPosition == 'HFS':
            if orient[0] == 'L':neworient[0] = 'R'
            else:neworient[0] = 'L'
            newpos = ['I', 'S', 'P', 'A']

        elif humanPosition == 'FFS':
            neworient[0] = orient[0]
            newpos = ['I', 'S', 'A', 'P']

        elif humanPosition == 'HFP':
            neworient[0] = orient[0]
            newpos = ['S', 'I', 'P', 'A']

        elif humanPosition == 'FFP':
            if orient[0] == 'L':neworient[0] = 'R'
            else:neworient[0] = 'L'
            newpos = ['S', 'I', 'A', 'P']

        for z in [1, 2]:
            for i, j in enumerate(pos):
                if orient[z] == j:
                    neworient[z] = newpos[i]

    elif animalPosition == 'humanlike':
        neworient = orient


    return neworient