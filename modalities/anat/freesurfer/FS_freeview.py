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

