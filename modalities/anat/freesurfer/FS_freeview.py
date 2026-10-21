#import
import os
import subprocess

opj = os.path.join
ope = os.path.exists
opi = os.path.isfile
spco = subprocess.check_output

def FSview(FS_dir, subID, surface,sing_fs,FS_refs):
    # from https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/TopologicalDefect_freeview
    
    Aseg_lut    = opj(FS_refs, 'FreeSurferAllLut.txt')

    if surface == 'white':
        if opi(opj(FS_dir,subID,'mri','segmentation.nii.gz')) == True:
            command = (sing_fs + 'freeview -v ' + opj(FS_dir,subID,'mri','orig.mgz') + ' ' +
            opj(FS_dir,subID,'mri','wm.mgz') + ':colormap=heat:opacity=0.3:visible=0 ' +
            opj(FS_dir,subID,'mri','aseg.mgz') + ':visible=0:colormap=lut:lut=' + Aseg_lut + ':opacity=0.4 ' +
            opj(FS_dir,subID,'mri','segmentation.nii.gz') + ':colormap=GE_COLOR:colorscale=0,4:opacity=0.2:visible=1' +
            ' -f ' + opj(FS_dir,subID,'surf','lh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','rh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','lh.orig.nofix') + ':edgecolor=green:visible=0 ' +
            opj(FS_dir,subID,'surf','rh.orig.nofix') + ':edgecolor=green:visible=0')
        else:
            command = (sing_fs + 'freeview -v ' + opj(FS_dir,subID,'mri','orig.mgz') + ' ' +
            opj(FS_dir,subID,'mri','wm.mgz') + ':colormap=heat:opacity=0.3:visible=0 ' +
            opj(FS_dir,subID,'mri','aseg.mgz') + ':visible=0:colormap=lut:lut=' + Aseg_lut + ':opacity=0.4 ' +
            ' -f ' + opj(FS_dir,subID,'surf','lh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','rh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','lh.orig.nofix') + ':edgecolor=green:visible=0 ' +
            opj(FS_dir,subID,'surf','rh.orig.nofix') + ':edgecolor=green:visible=0')

        spco([command], shell=True)
            
    elif surface == 'pial':
        if ope(opj(FS_dir,subID,'mri','segmentation.nii.gz')) == True:
            command = (sing_fs + 'freeview -v ' + opj(FS_dir,subID,'mri','orig.mgz') + ' ' +
            opj(FS_dir,subID,'mri','segmentation.nii.gz') + ':colormap=GE_COLOR:opacity=0.2:visible=0 ' +
            '-f ' + opj(FS_dir,subID,'surf','lh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','rh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','lh.pial') + ':edgecolor=blue ' +
            opj(FS_dir,subID,'surf','rh.pial') + ':edgecolor=blue')
        else:
            command = (sing_fs + 'freeview -v ' + opj(FS_dir,subID,'mri','orig.mgz') + ' ' +
            '-f ' + opj(FS_dir,subID,'surf','lh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','rh.white') + ':edgecolor=yellow ' +
            opj(FS_dir,subID,'surf','lh.pial') + ':edgecolor=blue ' +
            opj(FS_dir,subID,'surf','rh.pial') + ':edgecolor=blue')

        spco([command], shell=True)



# colour maps of check_coreg (suggested by the PET engineer)
colormap = {'Swhite': 'yellow', 'Spial': 'blue', 'pet': 'heat', 'ct': 'grayscale', 'seg': 'lut'}


def check_coreg(img, refimg, imgtype, sing_fs, surfaces=(), transform_volume=True):
    """
    Open freeview to check (and, with transform_volume, correct by hand) the registration of
    `img` (PET or CT) on the anatomical image `refimg`.

    surfaces : white / pial GIFTI surfaces in the space of `refimg` (optional), drawn as outlines.
    With transform_volume=True, freeview opens its 'Transform Volume' tool on `img`: after a manual
    correction the user saves the moved image under its own name ('Save Volume' with resampling),
    and the calling code re-estimates the registration from it.
    """
    for f in (img, refimg):
        if not opi(f):
            raise FileNotFoundError(f'check_coreg: file not found: {f}')
    cmd = (sing_fs + 'freeview -v ' + refimg + ':colormap=grayscale:opacity=1:visible=1 ' +
           img + ':colormap=' + colormap.get(imgtype, 'heat') + ':opacity=0.5:visible=1')
    if transform_volume:
        cmd += ' -transform-volume'
    outlines = [s for s in surfaces if opi(s)]
    if outlines:
        cmd += ' -f ' + ' '.join(s + ':edgecolor=' + (colormap['Spial'] if 'pial' in os.path.basename(s)
                                                       else colormap['Swhite']) for s in outlines)
    spco([cmd], shell=True)
