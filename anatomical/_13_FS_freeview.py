#import
import os
import subprocess

#Path to the excels files and data structure
opj = os.path.join
spco = subprocess.check_output

def FS_Freeview(FS_dir, monkey, surface, Aseg_lut):
    # from https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/TopologicalDefect_freeview
    
    if surface == 'white':
        command = 'freeview -v ' + opj(FS_dir,monkey,'mri','orig.mgz') + ' ' + \
        opj(FS_dir,monkey,'mri','wm.mgz') + ':colormap=heat:opacity=0.3:visible=0 ' + \
        opj(FS_dir,monkey,'mri','aseg.mgz') + ':visible=0:colormap=lut:lut=' + Aseg_lut + ':opacity=0.4 ' + \
        ' -f ' + opj(FS_dir,monkey,'surf','lh.white') + ':edgecolor=yellow ' + \
        opj(FS_dir,monkey,'surf','rh.white') + ':edgecolor=yellow ' + \
        opj(FS_dir,monkey,'surf','lh.orig.nofix') + ':edgecolor=green:visible=0 ' + \
        opj(FS_dir,monkey,'surf','rh.orig.nofix') + ':edgecolor=green:visible=0'
        spco([command], shell=True)
    elif surface == 'pial':
        command = 'freeview -v ' + opj(FS_dir,monkey,'mri','orig.mgz') + ' ' + \
        '-f ' + opj(FS_dir,monkey,'surf','lh.white') + ':edgecolor=yellow ' + \
        opj(FS_dir,monkey,'surf','rh.white') + ':edgecolor=yellow ' + \
        opj(FS_dir,monkey,'surf','lh.pial') + ':edgecolor=blue ' + \
        opj(FS_dir,monkey,'surf','rh.pial') + ':edgecolor=blue'
        spco([command], shell=True)

