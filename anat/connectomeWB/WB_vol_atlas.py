#import
import os
import shutil
import numpy as np
import json
import ants

opj = os.path.join
opi = os.path.isfile

from Tools import run_cmd
from Tools import getpath


def MY19(cmd_mris_convert,animal,FS_dir,path_BALSA,dir_path,Ref_file,path_label,Hmin,CORTEX,diary_name,sing_wb,sing_fs):

    BALSAname = 'MY19'
    lalbenames = opj(path_label,BALSAname + '_label.txt')

    nl = 'Create MY19 volumetric atlas'
    run_cmd.msg(nl,diary_name,'HEADER')

    (path_anat, matrices_dir, _, prepro_dir, dir_native, volumes_dir, labels_dir, masks_dir, wb_template_dir,
     wb_template_vol, wb_template_labels, wb_template_masks,
     wb_balsa_dir, wb_balsa_vol, wb_balsa_labels, wb_balsa_masks) = getpath.anat(dir_path, BALSAname, BALSAname, '', '', 'template')

    dir_native_resol, dir_native_32, dir_balsa_resol, dir_balsa_32, dir_balsa_64 = getpath.surf(path_anat, BALSAname,
                                                                                                BALSAname)


    Yerkes_label = opj(path_BALSA,'surfaces','fsaverage_LR_32k',BALSAname + '_Parcellations_v2.32k_fs_LR.dlabel.nii')

    list_Yerkes = ['MW', 'LV00', 'FV91', 'PHT00', 'M129', 'M132', 'B05', 'BB47', 'UD86',
                   'SP78+','LK02', 'FOA00', 'V6','MOD', 'KMA09', 'LV00_FOA00_PHT00', 'PFC']

    nl = '# 1- Adapt the template dlabel.nii file for every spec files'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    shutil.copyfile(Yerkes_label, opj(dir_balsa_32, animal + '.' + BALSAname + '.32_fs_LR.dlabel.nii'))
    nl = '## 1a Create a link for the surfaces that have the same resolution'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal + '_fsaverage_LR_32k.wb.spec') +
           ' INVALID ' + opj(dir_balsa_32, animal + '.' + BALSAname +  '.32_fs_LR.dlabel.nii'))
    run_cmd.wb(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_balsa_32, animal + '_fsaverage_LR_32k.wb.spec') +
           ' INVALID ' + opj(dir_balsa_32, animal + '.' + BALSAname +  '.32_fs_LR.dlabel.nii'))
    run_cmd.wb(cmd,diary_name)

    nl = '# 1b - Change the surface resolution'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    for s in range(len(Hmin)):
        cmd = (sing_wb + 'wb_command -cifti-separate ' + opj(dir_balsa_32, animal + '.' + BALSAname +  '.32_fs_LR.dlabel.nii') +
               ' COLUMN -label ' + CORTEX[s] +
               ' ' + opj(dir_balsa_32, BALSAname +  '_Parcellations_v2.' + Hmin[s] + '.32k_fs_LR.label.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -label-resample ' + opj(dir_balsa_32, BALSAname +  '_Parcellations_v2.' + Hmin[s] + '.32k_fs_LR.label.gii') +
               ' ' + opj(dir_balsa_32, animal + '.' + Hmin[s] + '.sphere.32k_fs_LR.surf.gii') +
               ' ' + opj(dir_balsa_resol, animal + '.' + Hmin[s] + '.sphere.reg.reg_LR.native.surf.gii') +
               ' BARYCENTRIC ' +
               opj(dir_balsa_resol, animal + '.' + Hmin[s] + '.' + BALSAname +  '_fs_LR.label.gii') + ' -largest')
        run_cmd.wb(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -cifti-create-label ' + opj(dir_balsa_resol, animal + '.' + BALSAname +  '.dlabel.nii') +
           ' -left-label '  + opj(dir_balsa_resol, animal + '.l.' + BALSAname +  '_fs_LR.label.gii') +
           ' -roi-left '    + opj(dir_balsa_resol, animal + '.l.roi.shape.gii') +
           ' -right-label ' + opj(dir_balsa_resol, animal + '.r.' + BALSAname +  '_fs_LR.label.gii') +
           ' -roi-right '   + opj(dir_balsa_resol, animal + '.r.roi.shape.gii'))
    run_cmd.wb(cmd,diary_name)

    shutil.copyfile(opj(dir_balsa_resol, animal + '.' + BALSAname +  '.dlabel.nii'),
                    opj(dir_native_resol, animal + '.' + BALSAname +  '.dlabel.nii'))

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_balsa_resol, animal + '_native_LR.wb.spec') +
           ' INVALID ' + opj(dir_balsa_resol, animal + '.' + BALSAname +  '.dlabel.nii'))
    run_cmd.wb(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') +
           ' INVALID ' + opj(dir_native_resol, animal + '.' + BALSAname +  '.dlabel.nii'))
    run_cmd.wb(cmd,diary_name)

    nl = '# 2 - create an annot version for Freesurfer'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    #   To do so it has to be the LAST map in the cifti file so we reorder it

    for i in range(len(list_Yerkes)):
        #old_list = list_Yerkes.copy()
        #old_list.remove(list_Yerkes[i])
        #new_list = old_list + [list_Yerkes[i]]

        olist = np.arange(1,len(list_Yerkes)+1)
        olist = olist[olist!=i+1]
        nlist = np.append(olist,i+1)

        with open(opj(dir_native_resol, 'new_list.txt'), 'wb') as f:
            np.savetxt(f, nlist , fmt='%.0f')
        cmd = (sing_wb + 'wb_command -cifti-reorder ' + opj(dir_native_resol, animal + '.' + BALSAname +  '.dlabel.nii') + ' ROW ' +
               opj(dir_native_resol, 'new_list.txt') + ' ' + opj(dir_native_resol,'new.dlabel.nii'))
        run_cmd.wb(cmd,diary_name)

        for s in range(len(Hmin)):
            cmd = (sing_wb + 'wb_command -cifti-separate ' + opj(dir_native_resol,'new.dlabel.nii') +
                   ' COLUMN -label ' + CORTEX[s] +
                   ' ' + opj(dir_native_resol, 'new.' + Hmin[s] + '.label.gii'))
            run_cmd.wb(cmd, diary_name)

            cmd = (sing_fs + cmd_mris_convert + ' --annot ' +
                   opj(dir_native_resol, 'new.' + Hmin[s] + '.label.gii') + ' ' +
                   opj(FS_dir, animal, 'surf', Hmin[s] + 'h.white') + ' ' +
                   opj(FS_dir, animal, 'label', Hmin[s] + 'h.' + animal + '_' + BALSAname +  '-' + list_Yerkes[i] + '.annot'))
            run_cmd.wb(cmd,diary_name)

        os.remove(opj(dir_native_resol, 'new_list.txt'))
        os.remove(opj(dir_native_resol, 'new.dlabel.nii'))
        os.remove(opj(dir_native_resol, 'new.l.label.gii'))
        os.remove(opj(dir_native_resol, 'new.r.label.gii'))

    nl = '# 3 - convert the available surface Atlases into volumes Atlases'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    #img = ants.image_read(Ref_file)

    for s in range(len(Hmin)):
        cmd = (sing_wb + 'wb_command -cifti-separate ' + opj(dir_native_resol, animal + '.' + BALSAname +  '.dlabel.nii') +
               ' COLUMN -label ' + CORTEX[s] +
               ' ' + opj(dir_native_resol, 'new.' + Hmin[s] + '.label.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -label-to-volume-mapping ' + opj(dir_native_resol, 'new.' + Hmin[s] + '.label.gii') + ' ' +
               opj(dir_native_resol, animal + '.' + Hmin[s] + '.white.surf.gii') + ' -ribbon-constrained ' +
               opj(dir_native_resol, animal + '.' + Hmin[s] + '.white.surf.gii') + ' ' +
               opj(dir_native_resol, animal + '.' + Hmin[s] + '.pial.surf.gii') + ' ' +
               Ref_file +
               ' ' + opj(volumes_dir, animal + '_' + Hmin[s] + '_' + BALSAname + '.nii.gz'))
        run_cmd.wb(cmd,diary_name)

    img1  = ants.image_read(opj(volumes_dir, animal + '_l_' + BALSAname + '.nii.gz'))
    tmp_1 = ants.ndimage_to_list(img1)
    img2  = ants.image_read(opj(volumes_dir, animal + '_r_' + BALSAname + '.nii.gz'))
    tmp_2 = ants.ndimage_to_list(img2)
    img3  = list()
    for i in range(len(list_Yerkes)):
        #test = tmp_2[i] +2000
        #test[test<2000.1]=0
        test = tmp_1[i] + tmp_2[i]
        img3.append(test)

    new_img = ants.list_to_ndimage(img1,img3)
    ants.image_write(new_img,opj(labels_dir, animal + '_seg-' + BALSAname + '_dseg.nii.gz'),ri=False)
    
    cmd = (sing_wb + 'wb_command -volume-label-import ' + opj(labels_dir, animal + '_seg-' + BALSAname + '_dseg.nii.gz') +
           ' ' + lalbenames +
           ' ' + opj(labels_dir, animal + '_seg-' + BALSAname + '_dseg.nii.gz') + ' -drop-unused-labels')
    run_cmd.wb(cmd,diary_name)

    dictionary = {"Sources": [opj(dir_native_resol, animal + '.' + BALSAname + '.dlabel.nii'),
                              Ref_file],
                  "Description": 'Get a label volume (wb_command and Antspy)', }
    json_object = json.dumps(dictionary, indent=2)

    with open(opj(labels_dir, animal + '_seg-' + BALSAname + '_dseg.json'),"w") as outfile:
        outfile.write(json_object)

    os.remove(opj(volumes_dir, animal + '_l_' + BALSAname + '.nii.gz'))
    os.remove(opj(volumes_dir, animal + '_r_' + BALSAname + '.nii.gz'))
    os.remove(opj(dir_native_resol, 'new.l.label.gii'))
    os.remove(opj(dir_native_resol, 'new.r.label.gii'))

    nl = 'volumetric atlas of yerkes: done!'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    
