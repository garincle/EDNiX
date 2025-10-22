#import
import os
import glob

opj = os.path.join
ope = os.path.exists
opb = os.path.basename
opi = os.path.isfile

from Tools import run_cmd
from anatomical.connectomeWB import WB_label

surface = ['Left','Right']
Hmin    = ['l','r']
DESGCS = 'h.destrieux.simple.2009-07-29.gcs'
DKTGCS = 'h.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs'

def FS_finalise(FS_dir, subject, label, CTAB,subC,change_hd,diary_name,sing_fs,export_fs,path_code_label,specie,FS_refs,sing_afni,sing_wb):

    # 1) Create the ribbon.mgz
    nl = 'Creation of the ribbon.nii.gz labelling volume'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    cmd = export_fs + 'mris_volmask --aseg_name aseg --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon ' + subject
    run_cmd.do(cmd,diary_name)

    os.remove(opj(FS_dir,subject,'mri', 'lh.ribbon.mgz'))
    os.remove(opj(FS_dir, subject, 'mri', 'rh.ribbon.mgz'))

    list_annot = glob.glob(opj(FS_dir, subject,'label','*' + subject + '*.annot'))
    for i in list_annot:
        os.remove(i)

    for H in range(2):

        # 2) Midthickness surface
        nl = 'Creation of the Midthickness surface'
        run_cmd.msg(nl, diary_name, 'OKGREEN')

        cmd = sing_fs + 'mris_expand -thickness ' + opj(FS_dir,subject,'surf', Hmin[H] + 'h.white') + ' 0.5 ' + opj(FS_dir,subject,'surf',Hmin[H] + 'h.mid')
        run_cmd.run(cmd,diary_name)

        if specie =='Human':
            cmd = (export_fs + 'mris_anatomical_stats  -th3 -mgz -cortex ' +
                               opj(FS_dir, subject, 'label', Hmin[H] + 'h.cortex.label') +
                               ' -f ' + opj(FS_dir, subject, 'stats', Hmin[H] + 'h.aparc.stats') +
                               ' -b -a ' + opj(FS_dir, subject, 'label', Hmin[H] + 'h.aparc.annot') +
                               ' -c ' + opj(FS_dir, subject, 'label', 'aparc.annot.ctab') +
                               ' ' + subject + ' ' + Hmin[H] + 'h')
            run_cmd.do(cmd,diary_name)

            cmd = (export_fs + 'mris_ca_label -l ' + opj(FS_dir, subject, 'label', Hmin[H] + 'h.cortex.label') +
                   ' -aseg ' + opj(FS_dir, subject, 'mri', 'aseg.mgz') +
                   ' ' + subject + ' ' + Hmin[H] + 'h' +
                   ' ' + opj(FS_dir, subject, 'surf', Hmin[H] + 'h.sphere.reg') +
                   ' ' + opj(FS_refs, Hmin[H] + DESGCS) +
                   ' ' + Hmin[H] + 'h.aparc.a2009s.annot')
            run_cmd.do(cmd,diary_name)

            cmd = (export_fs + 'mris_anatomical_stats -th3 -mgz -cortex ' +
                   opj(FS_dir, subject, 'label', Hmin[H] + 'h.cortex.label') +
                   ' -f ' + opj(FS_dir, subject, 'stats', Hmin[H] + 'h.aparc.a2009s.stats') +
                   ' -b -a ' + opj(FS_dir, subject, 'label', Hmin[H] + 'h.aparc.a2009s.annot') +
                   ' -c ' + opj(FS_dir, subject, 'label','aparc.annot.a2009s.ctab') +
                   ' ' + subject  +' ' + Hmin[H] + 'h')
            run_cmd.do(cmd,diary_name)

            cmd = (export_fs + 'mris_ca_label -l ' + opj(FS_dir, subject, 'label', Hmin[H] + 'h.cortex.label') +
                  ' -aseg ' + opj(FS_dir, subject, 'mri', 'aseg.mgz') +
                  ' ' + subject + ' ' + Hmin[H] + 'h' +
                  ' ' + opj(FS_dir, subject, 'surf', Hmin[H] + 'h.sphere.reg') +
                  ' ' + opj(FS_refs, Hmin[H] + DKTGCS) +
                  ' ' + Hmin[H] + 'h.aparc.DKTatlas.annot')
            run_cmd.do(cmd,diary_name)

            cmd = (export_fs + 'mris_anatomical_stats -th3 -mgz -cortex ' +
                   opj(FS_dir, subject, 'label', Hmin[H] + 'h.cortex.label') +
                   ' -f ' + opj(FS_dir, subject, 'stats', Hmin[H] + 'h.aparc.DKTatlas.stats') +
                   ' -b -a ' + opj(FS_dir, subject, 'label', Hmin[H] + 'h.aparc.DKTatlas.annot') +
                   ' -c ' + opj(FS_dir, subject, 'label','aparc.annot.DKTatlas.ctab') +
                   ' ' + subject  + ' ' + Hmin[H] + 'h')
            run_cmd.do(cmd,diary_name)

        nl = surface[H] + ' surface done!'
        run_cmd.msg(nl, diary_name, 'OKGREEN')
        

    # 3) Add the atlases (remaining questions: What is the most appropriate surface to use ? what is the best projection fraction ?)

    for i,j,k in zip(label,CTAB,subC):
        if opi(i)==True:
            if k == 1:
                WB_label.surfFS(i, j, change_hd, diary_name, sing_fs, path_code_label, export_fs,sing_afni,sing_wb)


    '''
    if specie == 'Human':
        
        cmd = (export_fs + 'mri_relabel_hypointensities ' + opj(FS_dir, subject,'mri','aseg.mgz') +
               ' ' + opj(FS_dir, subject, 'surf') + ' ' + opj(FS_dir, subject,'mri','aseg.hypos.mgz'))
        nl = spgo(cmd)
        print(nl)
        

        cmd = (export_fs + 'mri_aparc2aseg --s ' + subject + ' --volmask --aseg aseg.hypos' +
               ' --relabel ' + opj(FS_dir, subject, 'mri', 'brain.mgz') +
               ' ' + opj(FS_refs, 'RB_all_2008-03-26.gca') +
               ' ' + opj(FS_dir, subject, 'mri', 'aseg.auto_noCCseg.label_intensities.txt'))
        nl = spgo(cmd)
        print(nl)

        cmd = (export_fs + 'mri_aparc2aseg --s ' + subject + ' --volmask --annot aparc.a2009s --aseg ' + opj(FS_dir, subject,'mri','aseg.hypos.mgz') +
               ' --relabel ' + opj(FS_dir, subject, 'mri', 'brain.mgz') +
               ' ' + opj(FS_refs, 'RB_all_2008-03-26.gca') +
               ' ' + opj(FS_dir, subject, 'mri', 'aseg.auto_noCCseg.label_intensities.txt'))
        nl = spgo(cmd)
        print(nl)

        cmd = (export_fs + 'mri_aparc2aseg --s ' + subject + ' --volmask --annot aparc.DKTatlas --aseg ' + opj(FS_dir, subject,'mri','aseg.hypos.mgz') +
                    ' --relabel ' + opj(FS_dir, subject, 'mri', 'brain.mgz') +
                    ' ' + opj(FS_refs, 'RB_all_2008-03-26.gca') +
                    ' ' + opj(FS_dir, subject, 'mri', 'aseg.auto_noCCseg.label_intensities.txt'))
        nl = spgo(cmd)
        print(nl)
        cmd = (export_fs + 'apas2aseg --i aparc+aseg.mgz --o aseg.mgz ')
        nl = spgo(cmd)
        print(nl)
        cmd = (export_fs + 'mri_segstats --seg ' + opj(FS_dir, subject, 'mri', 'aseg.mgz') +
               ' --sum ' + opj(FS_dir, subject, 'stats', 'aseg.stats') +
               ' --pv ' + opj(FS_dir, subject, 'mri', 'brain.mgz') +
               ' --empty --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray' +
               ' --in ' + opj(FS_dir, subject, 'mri', 'brain.mgz') +
               ' --in-intensity-name brain --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol'
               ' --totalgray --euler --ctab ' + opj(path_code_label, 'ASegStatsLUT.txt') + ' --subject ' + subject)
        nl = spgo(cmd)
        print(nl)
        cmd = (export_fs + 'mri_aparc2aseg --s ' + subject + '--labelwm --hypo-as-wm --rip-unknown --volmask' +
               ' --o ' + opj(FS_dir, subject, 'mri', 'wmparc.mgz') + ' --ctxseg aparc+aseg.mgz ')
        nl = spgo(cmd)
        print(nl)

        cmd = (export_fs + 'mri_segstats --seg ' + opj(FS_dir, subject, 'mri', 'wmparc.mgz') +
               ' --sum ' + opj(FS_dir, subject, 'stats', 'wmparc.stats') +
               ' --pv ' + opj(FS_dir, subject, 'mri', 'brain.mgz') +
               ' --brain-vol-from-seg --excludeid 0' +
               ' --in ' + opj(FS_dir, subject, 'mri', 'brain.mgz') +
               ' --in-intensity-name brain --in-intensity-units MR --etiv --surf-wm-vol'
               ' --ctab ' + opj(path_code_label, 'WMParcStatsLUT.txt') + ' --subject ' + subject)
        nl = spgo(cmd)
        print(nl)
        # cmd = (export_fs + 'mri_label2label --srcsubject fsaverage --srclabel fsaverage/label/?h.BA*.label' +
        #       ' --trgsubject ' + subject + ' --trglabel ' + ? + 'h.BA*.label --hemi ' + ? + 'h --regmethod surface ')
    '''


    

    
