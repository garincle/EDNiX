#import
import os
import subprocess

opj = os.path.join
ope = os.path.exists
opi = os.path.isfile
opb = os.path.basename
spco = subprocess.check_output


BALSAsamp= {'resol':['native','native','native_LR',''],
            '164k':['fsaverage_164k','164k_fs_LR','fsaverage_LR_164k','atlasroi'],
            '32k':['fsaverage_LR_32k','32k_fs_LR','fsaverage_LR_32k','atlasroi']}

NATIVEsamp= {'resol':['Native_resol','','native_LR',''],
            '32k':['fsaverage_LR_32k','32k_fs_LR','fsaverage_LR_32k','atlasroi']}

colormap = {'Swhite':'yellow','Spial':'blue','Vwhite':'heat','seg':'GE_COLOR','pet':'GE_COLOR','ct':'grayscale'}


def check_seg(FS_dir, subID, surface, sing_fs, FS_refs):
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

    def check_coreg(img, refimg, space, S_folder,type, sing_fs):

        cmd = (sing_fs + 'freeview -v ' + refimg + ':colormap=grayscale:opacity=1:visible=1 ' +
               img + ':colormap=' + colormap[type] + ':opacity=0.5:visible=1 -transform-volume')

        if space == 'acpc':
            preference = ['32k', 'resol']
            name = NATIVEsamp
        elif space == 'template':
            preference = ['32k', '164k', 'resol']
            name = BALSAsamp

        # to see the surfaces may help .....
        SURF = ''
        fname = opb(refimg).split('_')[0]

        if S_folder:
            for i in preference:
                if ope(opj(S_folder, BALSAsamp[i][0])):
                    SURF = BALSAsamp[i]
                    break

        if SURF:
            lwhite = opj(S_folder, SURF[0],
                         '.'.join(filter(None, [fname, 'l', 'white', SURF[1], 'surf', 'gii'])))
            rwhite = opj(S_folder, SURF[0],
                         '.'.join(filter(None, [fname, 'r', 'white', SURF[1], 'surf', 'gii'])))
            lpial = opj(S_folder, SURF[0],
                        '.'.join(filter(None, [fname, 'l', 'pial', SURF[1], 'surf', 'gii'])))
            rpial = opj(S_folder, SURF[0],
                        '.'.join(filter(None, [fname, 'r', 'pial', SURF[1], 'surf', 'gii'])))

            print(lwhite)
            cmd = cmd + (' -f ' + lwhite + ':edgecolor=' + colormap['Swhite'] +
                         ' ' + rwhite + ':edgecolor=' + colormap['Swhite'] +
                         ' ' + lpial + ':edgecolor=' + colormap['Spial'] +
                         ' ' + rpial + ':edgecolor=' + colormap['Spial'])

        spco([cmd], shell=True)

