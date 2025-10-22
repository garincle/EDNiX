#import
import os
import numpy as np
import shutil
import ants
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname

opi = os.path.isfile


from Tools import run_cmd,make_iso_img
from anatomical.freesurfer import preFS, smallbrain

Hmin    = ['l', 'r']
surface = ['Left', 'Right']
CORTEX  = ['CORTEX_LEFT', 'CORTEX_RIGHT']

def vol(atlas,anat,label_dir,spec,transfo,w2i,label,diary_name,sing_wb):

    nl = 'Make volume ' + atlas + ' readable for WB'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    volume = opj(label_dir, opb(anat).split('_')[0] + '_seg-' + opb(label).split('_')[0] + '_dseg.nii.gz')

    if opi(atlas) == True:

        header = ants.image_header_info(atlas)
        if len(header['dimensions']) > 3:
            dim4 = 3
        else:
            dim4 = 0

        if transfo == '':
            shutil.copyfile(atlas,volume)
        else :
            img = ants.image_read(atlas)
            new = ants.image_read(anat)
            moved = ants.apply_transforms(fixed=new, moving=img,
                                          transformlist=transfo,
                                          interpolator='genericLabel',
                                          whichtoinvert=w2i,
                                          imagetype=dim4)
            ants.image_write(moved,volume,ri=False)

        cmd = (sing_wb + 'wb_command -volume-label-import ' + volume + ' ' + label +
               ' ' + volume + ' -drop-unused-labels')
        run_cmd.wb(cmd,diary_name)

        dictionary = {"Sources": [atlas,
                                  anat],
                      "Description": 'coregistration (Antspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(volume.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        if not spec == '':
            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + spec +' INVALID ' + volume)
            run_cmd.wb(cmd,diary_name)

    

def surfFS(volume,ctab,change_hd,diary_name,sing_fs,path_code_label,export_fs,sing_afni,sing_wb):

    nl = 'Create annot label files for Freesurfer '
    run_cmd.msg(nl, diary_name,'OKGREEN')


    # the labelled volume is supposed to be in animal/native/volumes/label/
    # the Freesurfer folder in ../../freesurfer/animal
    data_path  = opd(opd(opd(opd(volume))))
    FS_dir     = opj(data_path,'freesurfer')
    dir_prepro = opj(data_path,'preprocessing')

    volname = opb(volume).split('_')
    animal = volname[0]
    name =''
    for i, j in enumerate(volname[:-1]):
        if 'seg' in j:
            name = volname[i].split('-')[-1]


    if change_hd == 1:
        _, spacing, _, new_size, _ = smallbrain.check(opj(dir_prepro, animal + '_rescale.json'))

        ref      = opj(dir_prepro, animal + '_desc-norm_T1w_resamp-' + str(spacing[0]) +'.nii.gz')
        run_cmd.msg(ref, diary_name, 'OKGREEN')
    else:
        ref = opj(dir_prepro, animal + '_desc-norm_T1w.nii.gz')
        spacing  = ''
        new_size = ''

    header = ants.image_header_info(volume)
    iso = len(np.unique(np.array(header['spacing'][0:3]) * 1000))

    if iso>1:
        img_out = '_'.join(volname[:-1] + ['res-iso'] + [volname[-1]])
        newvolume = opj(opd(volume), img_out)
        make_iso_img.make_iso(volume, newvolume, diary_name, sing_afni, 'seg', ' -overwrite')
        cmd = (sing_wb + 'wb_command -volume-label-import' + ' ' + newvolume +
               ' ' + opj(path_code_label, name + '_label.txt') + ' ' + newvolume + ' -drop-unused-labels')
        run_cmd.run(cmd, diary_name)
        dictionary = {"Sources": volume,
                      "Description": 'resampling to isometric voxels.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(newvolume.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)
        volume=newvolume


    if header['nDimensions'] > 3:
        img4D = ants.image_read(volume)
        img4D_split = ants.ndimage_to_list(img4D)
        for i in range(int(header['dimensions'][3])):
            ants.image_write(img4D_split[i], opj(opd(volume), animal + '_seg-' + name + str(i + 1) + '_dseg.nii.gz'),
                             ri=False)

            
            preFS.toannot(opj(opd(volume), animal + '_seg-' + name + str(i + 1) + '_dseg.nii.gz'),
                          ref, change_hd, spacing, new_size, FS_dir, diary_name, sing_fs)

            for H in range(2):
                if ctab == 'txt':
                    ext = '_StatsLUT.txt'
                elif ctab == 'ctab':
                    ext = '_' + Hmin[H] + '.ctab'
                cmd = (export_fs + 'mris_sample_parc -ct ' + opj(path_code_label, name) + ext +
                       ' -surf mid -projfrac 0.01 ' + animal + ' ' + Hmin[H] + 'h ' + animal + '_' + name + str(i+1) + '.mgz ' +
                       Hmin[H] + 'h.' + animal + '_' + name + str(i+1) + '.annot')
                run_cmd.do(cmd,diary_name)

            os.remove(opj(opd(volume), animal + '_seg-' + name +  str(i+1) + '_dseg.nii.gz'))

    else :
        preFS.toannot(volume, ref, change_hd, spacing, new_size, FS_dir, diary_name, sing_fs)

    for H in range(2):
        if ctab == 'txt':
            ext = '_StatsLUT.txt'
        elif ctab == 'ctab':
            ext = '_' + Hmin[H] +  '.ctab'

        cmd = (export_fs + 'mris_sample_parc -ct ' + opj(path_code_label,name) + ext +
               ' -surf mid -projfrac 0.01 ' + animal + ' ' + Hmin[H] + 'h ' + animal + '_' + name + '.mgz ' +
               Hmin[H] + 'h.' + animal + '_' + name + '.annot')
        run_cmd.do(cmd,diary_name)

        nl = surface[H] + ' surface done!'
        run_cmd.msg(nl, diary_name,'OKGREEN')

    

def surfWB(cmd_mris_convert,FS_dir, animal,label,WB_dir,surftype,h,diary_name,sing_fs,sing_wb):

    nl = 'create label surfaces for WB'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    newgii = opj(WB_dir, animal + '.' + Hmin[h] + '.' + label + surftype + '.label.gii')

    cmd = (sing_fs + cmd_mris_convert + ' --annot ' + opj(FS_dir, animal, 'label',Hmin[h] + 'h.' + animal + '_' + label + '.annot') +
           ' ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.white') +
           ' ' + newgii)
    run_cmd.wb(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -set-structure ' + newgii + ' ' + CORTEX[h])
    run_cmd.wb(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -set-map-names ' + newgii + ' -map 1 ' + animal + '_' + Hmin[h] + '_' + label)
    run_cmd.wb(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -gifti-label-add-prefix ' + newgii + ' ' + Hmin[h] + '_ ' + newgii)
    run_cmd.wb(cmd,diary_name)

    nl = surface[h] + ' surface done!'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    

def mergelabel(animal,WB_dir,surftype,h,diary_name,label,nb,sing_wb):

    nl = 'Merge ' + label + 'label surfaces for WB'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    newlabel = []
    
    for z in range(nb):
        newlabel.append('-label ' + opj(WB_dir, animal + '.' + Hmin[h] + '.' + label + str(z+1) + surftype + '.label.gii'))
    print(newlabel)
    to_merge = ' '.join(newlabel)
    
    cmd = (sing_wb + 'wb_command -label-merge ' +
            opj(WB_dir,animal + '.' + Hmin[h] + '.' + label + surftype + '.label.gii') +
            ' ' + to_merge)
    run_cmd.wb(cmd,diary_name)
    
    for z in range(nb):
        os.remove(opj(WB_dir, animal + '.' + Hmin[h] + '.' + label + str(z+1) + surftype + '.label.gii'))


def cifti(animal,label,WB_dir,surftype,roi,spec,diary_name,sing_wb):

    nl = 'Create label surfaces for WB'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    newcifti = opj(WB_dir, animal + '.' + label + '.dlabel.nii')

    cmd = (sing_wb + 'wb_command -cifti-create-label ' + newcifti +
           ' -left-label '  + opj(WB_dir, animal + '.l.' + label + surftype + '.label.gii') +
           ' -roi-left '    + opj(WB_dir, animal + '.l.' + roi + '.shape.gii') +
           ' -right-label ' + opj(WB_dir, animal + '.r.' + label + surftype + '.label.gii') +
           ' -roi-right '   + opj(WB_dir, animal + '.r.' + roi + '.shape.gii'))
    run_cmd.wb(cmd,diary_name)


    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + spec +' INVALID ' + newcifti)
    run_cmd.wb(cmd,diary_name)


def vo2surfWB(animal,list_atlas,vol_dir,WB_dir,surftype,hemi_type,proj,roi,spec,diary_name,sing_wb):

    nl = 'Create label surfaces for WB from volume'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    for n in range(len(list_atlas)):
        g_name = list_atlas[0][n]
        nb = list_atlas[2][n]
        if list_atlas[3][n] ==1:

            vol = opj(vol_dir,animal + '_seg-' + g_name + '_dseg.nii.gz')
            
            for h in hemi_type :
                newgii = opj(WB_dir,animal + '.' + h + '.' + g_name + surftype + '.label.gii')
                if proj == 'ribbon':
                    cmd = (sing_wb + 'wb_command -volume-label-to-surface-mapping ' + vol + 
                    ' ' + opj(WB_dir,animal + '.' + h + '.midthickness' + surftype + '.surf.gii') + 
                    ' ' + newgii + 
                    ' -ribbon-constrained ' + opj(WB_dir,animal + '.' + h+ '.white' + surftype + '.surf.gii') + 
                    ' ' + opj(WB_dir,animal + '.' + h + '.pial' + surftype + '.surf.gii'))
                elif proj == 'vox':
                    cmd = (sing_wb + 'wb_command -volume-label-to-surface-mapping ' + vol + 
                    ' ' + opj(WB_dir,animal + '.' + h + '.midthickness' + surftype + '.surf.gii') + 
                    ' ' + newgii )
                run_cmd.wb(cmd,diary_name)

                for z in range(nb):
                    cmd = (sing_wb + 'wb_command -set-map-names ' + newgii
                     + ' -map ' + str(z +1)+ ' ' + g_name + '-level-' + str(z +1))
                    run_cmd.wb(cmd,diary_name)

                    if list_atlas[1][n] =='txt':
                        cmd = (sing_wb + 'wb_command -gifti-label-add-prefix ' + newgii + ' ' + h + '_ ' + newgii)
                        run_cmd.wb(cmd,diary_name)
            
        cifti(animal,g_name,WB_dir,surftype,roi,spec,diary_name,sing_wb)

                    

