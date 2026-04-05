#import
import os
import shutil
import numpy as np
import json
import ants

opj = os.path.join

from Tools import run_cmd

Hmin        = ['l','r']
CORTEX      = ['CORTEX_LEFT','CORTEX_RIGHT']


def get(BALSAname, path_BALSA):

        if BALSAname == 'MY19':
            origFile = opj(path_BALSA, 'surfaces', 'fsaverage_LR_32k',BALSAname + '_Parcellations_v2.32k_fs_LR.dlabel.nii')
            parcellation_list = ['MW', 'LV00', 'FV91', 'PHT00', 'M129', 'M132', 'B05', 'BB47', 'UD86',
                                 'SP78+', 'LK02', 'FOA00', 'V6', 'MOD', 'KMA09', 'LV00_FOA00_PHT00', 'PFC']
            surfsamp = '32k'

        else:
            print('sorry not ready yet')

        return (origFile,parcellation_list, surfsamp)


def set(origFile, animal, BALSAname, parcellation_list, labelnames, dir_native_resol, dir_native_32, dir_balsa_resol,
        dir_balsa_32, dir_balsa_64, volumes_dir, labels_dir, FS_dir, Ref_file, surfsamp, Hmin,
        CORTEX, diary_name, cmd_mris_convert, sing_wb, sing_fs):
    if surfsamp == '32k':
        dir_orig = dir_balsa_32
        dir_native_samp = dir_native_32
    elif surfsamp == '164k':
        dir_orig = dir_balsa_64
        # dir_native_samp = dir_native_64

    label_suffix = '_'.join([surfsamp, 'fs', 'LR'])
    labelsurf_end = '.'.join([label_suffix, 'label', 'gii'])

    typefile = '.'.join(origFile.split('.')[-2:])

    if typefile == 'dlabel.nii':
        cifti(origFile, animal, BALSAname, dir_orig, dir_native_samp, surfsamp, Hmin, CORTEX, diary_name, sing_wb)

    gifti(opj(dir_orig, '.'.join([BALSAname, Hmin[0], labelsurf_end])),
          opj(dir_orig, '.'.join([BALSAname, Hmin[1], labelsurf_end])),
          animal, BALSAname, dir_orig, dir_balsa_resol, dir_native_resol, surfsamp, Hmin, diary_name, sing_wb)

    FSannot(animal, BALSAname, parcellation_list, dir_native_resol, FS_dir, Hmin, CORTEX, cmd_mris_convert, diary_name,
            sing_wb, sing_fs)

    volumes(animal, BALSAname, parcellation_list, labelnames, dir_native_resol, volumes_dir, labels_dir, Ref_file, Hmin,
            CORTEX, diary_name, sing_wb)


def cifti(origFile, animal, BALSAname, dir_orig, dir_native_samp, surfsamp, Hmin, CORTEX, diary_name, sing_wb):
    specsuffix = '_'.join(['fsaverage', 'LR', surfsamp])
    spec_end = '.'.join([specsuffix, 'wb', 'spec'])
    label_suffix = '_'.join([surfsamp, 'fs', 'LR'])
    labelcifti_end = '.'.join([label_suffix, 'dlabel', 'nii'])
    labelsurf_end = '.'.join([label_suffix, 'label', 'gii'])

    nl = '# 1- Add the template cifti file for every spec files'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    shutil.copyfile(origFile, opj(dir_orig, '.'.join([animal, BALSAname, labelcifti_end])))

    nl = '## 1a Create a link for the surfaces that have the same resolution'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_samp, '_'.join([animal, spec_end])) +
           ' INVALID ' + opj(dir_orig, '.'.join([animal, BALSAname, labelcifti_end])))
    run_cmd.wb(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_orig, '_'.join([animal, spec_end])) +
           ' INVALID ' + opj(dir_orig, '.'.join([animal, BALSAname, labelcifti_end])))
    run_cmd.wb(cmd, diary_name)

    nl = '# 1b - Change the surface resolution'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    for s in range(len(Hmin)):
        cmd = (sing_wb + 'wb_command -cifti-separate ' + opj(dir_orig, '.'.join([animal, BALSAname, labelcifti_end])) +
               ' COLUMN -label ' + CORTEX[s] + ' ' + opj(dir_orig, '.'.join([animal, BALSAname, Hmin[s], labelsurf_end])))
        run_cmd.wb(cmd, diary_name)


def gifti(origFileL, origFileR, animal, BALSAname, dir_orig, dir_balsa_resol, dir_native_resol, surfsamp, Hmin,
          diary_name, sing_wb):

    specsuffix    = '_'.join(['native', 'LR'])
    spec_end      = '.'.join([specsuffix, 'wb', 'spec'])
    label_suffix  = '_'.join([surfsamp, 'fs', 'LR'])
    labelsurf_end = '.'.join([label_suffix, 'label', 'gii'])
    surf_end      = '.'.join([label_suffix, 'surf', 'gii'])

    if not opj(dir_orig, '.'.join([animal, BALSAname, Hmin[0], labelsurf_end])):
        shutil.copyfile(origFileL, opj(dir_orig, '.'.join([animal, BALSAname, Hmin[0], labelsurf_end])))
    if not opj(dir_orig, '.'.join([animal, BALSAname, Hmin[1], labelsurf_end])):
        shutil.copyfile(origFileR, opj(dir_orig, '.'.join([animal, BALSAname, Hmin[1], labelsurf_end])))

    for s in range(len(Hmin)):
        cmd = (sing_wb + 'wb_command -label-resample ' + opj(dir_orig, '.'.join([animal, BALSAname, Hmin[s], labelsurf_end])) +
               ' ' + opj(dir_orig, '.'.join([animal, Hmin[s], 'sphere', surf_end])) +
               ' ' + opj(dir_balsa_resol,
                         '.'.join([animal, Hmin[s], 'sphere', 'reg', 'reg_LR.native', 'surf', 'gii'])) +
               ' BARYCENTRIC ' +
               opj(dir_balsa_resol, '.'.join([animal, Hmin[s], BALSAname, 'fs_LR', 'label', 'gii'])) + ' -largest')
        run_cmd.wb(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -cifti-create-label ' + opj(dir_balsa_resol,
                                                             '.'.join([animal, BALSAname, 'dlabel', 'nii'])) +
           ' -left-label ' + opj(dir_balsa_resol, '.'.join([animal, 'l', BALSAname, 'fs_LR', 'label', 'gii'])) +
           ' -roi-left ' + opj(dir_balsa_resol, '.'.join([animal,'l', 'roi', 'shape', 'gii'])) +
           ' -right-label ' + opj(dir_balsa_resol, '.'.join([animal, 'r', BALSAname, 'fs_LR', 'label', 'gii'])) +
           ' -roi-right ' + opj(dir_balsa_resol, '.'.join([animal,'r', 'roi', 'shape', 'gii'])))
    run_cmd.wb(cmd, diary_name)

    shutil.copyfile(opj(dir_balsa_resol, '.'.join([animal, BALSAname, 'dlabel', 'nii'])),
                    opj(dir_native_resol, '.'.join([animal, BALSAname, 'dlabel', 'nii'])))

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_balsa_resol, '_'.join([animal, spec_end])) +
           ' INVALID ' + opj(dir_balsa_resol, '.'.join([animal, BALSAname, 'dlabel', 'nii'])))
    run_cmd.wb(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, '_'.join([animal, spec_end])) +
           ' INVALID ' + opj(dir_native_resol, '.'.join([animal, BALSAname, 'dlabel', 'nii'])))
    run_cmd.wb(cmd, diary_name)


def FSannot(animal, BALSAname, parcellation_list, dir_native_resol, FS_dir, Hmin, CORTEX, cmd_mris_convert, diary_name,
            sing_wb, sing_fs):

    nl = '# 2 - create an annot version for Freesurfer'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    #   To do so it has to be the LAST map in the cifti file so we reorder it
    for i in range(len(parcellation_list)):

        olist = np.arange(1, len(parcellation_list) + 1)
        olist = olist[olist != i + 1]
        nlist = np.append(olist, i + 1)

        with open(opj(dir_native_resol, 'new_list.txt'), 'wb') as f:
            np.savetxt(f, nlist, fmt='%.0f')
        cmd = (sing_wb + 'wb_command -cifti-reorder ' + opj(dir_native_resol,
                                                            '.'.join([animal, BALSAname, 'dlabel', 'nii'])) + ' ROW ' +
               opj(dir_native_resol, 'new_list.txt') + ' ' + opj(dir_native_resol, 'new.dlabel.nii'))
        run_cmd.wb(cmd, diary_name)

        for s in range(len(Hmin)):
            cmd = (sing_wb + 'wb_command -cifti-separate ' + opj(dir_native_resol, 'new.dlabel.nii') +
                   ' COLUMN -label ' + CORTEX[s] +
                   ' ' + opj(dir_native_resol, 'new.' + Hmin[s] + '.label.gii'))
            run_cmd.wb(cmd, diary_name)

            cmd = (sing_fs + cmd_mris_convert + ' --annot ' +
                   opj(dir_native_resol, 'new.' + Hmin[s] + '.label.gii') + ' ' +
                   opj(FS_dir, animal, 'surf', Hmin[s] + 'h.white') + ' ' +
                   opj(FS_dir, animal, 'label',
                       Hmin[s] + 'h.' + animal + '_' + BALSAname + '-' + parcellation_list[i] + '.annot'))
            run_cmd.run(cmd, diary_name)

        os.remove(opj(dir_native_resol, 'new_list.txt'))
        os.remove(opj(dir_native_resol, 'new.dlabel.nii'))
        os.remove(opj(dir_native_resol, 'new.l.label.gii'))
        os.remove(opj(dir_native_resol, 'new.r.label.gii'))


def volumes(animal, BALSAname, parcellation_list, labelnames, dir_native_resol, volumes_dir, labels_dir, Ref_file, Hmin,
            CORTEX, diary_name, sing_wb):
    nl = '# 3 - convert the available surface Atlases into volumes Atlases'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    for s in range(len(Hmin)):
        cmd = (sing_wb + 'wb_command -cifti-separate ' + opj(dir_native_resol,
                                                             '.'.join([animal, BALSAname, 'dlabel', 'nii'])) +
               ' COLUMN -label ' + CORTEX[s] +
               ' ' + opj(dir_native_resol, '.'.join(['new', Hmin[s], 'label', 'gii'])))
        run_cmd.wb(cmd, diary_name)

        cmd = (sing_wb + 'wb_command -label-to-volume-mapping ' + opj(dir_native_resol,
                                                                      '.'.join(['new', Hmin[s], 'label', 'gii']))
               + ' ' +
               opj(dir_native_resol, '.'.join([animal, Hmin[s], 'white', 'surf', 'gii'])) + ' -ribbon-constrained ' +
               opj(dir_native_resol, '.'.join([animal, Hmin[s], 'white.surf.gii'])) + ' ' +
               opj(dir_native_resol, '.'.join([animal, Hmin[s], 'pial', 'surf', 'gii'])) + ' ' +
               Ref_file +
               ' ' + opj(volumes_dir, '_'.join([animal, Hmin[s], BALSAname + '.nii.gz'])))
        run_cmd.wb(cmd, diary_name)

    img1 = ants.image_read(opj(volumes_dir, '_'.join([animal, 'l', BALSAname + '.nii.gz'])))
    tmp_1 = ants.ndimage_to_list(img1)
    img2 = ants.image_read(opj(volumes_dir, '_'.join([animal, 'r', BALSAname + '.nii.gz'])))
    tmp_2 = ants.ndimage_to_list(img2)
    img3 = list()
    for i in range(len(parcellation_list)):
        test = tmp_1[i] + tmp_2[i]
        img3.append(test)

    new_img = ants.list_to_ndimage(img1, img3)
    ants.image_write(new_img, opj(labels_dir, '_'.join([animal, 'seg-' + BALSAname, 'dseg.nii.gz'])), ri=False)

    cmd = (sing_wb + 'wb_command -volume-label-import ' + opj(labels_dir,
                                                              '_'.join([animal, 'seg-' + BALSAname, 'dseg.nii.gz'])) +
           ' ' + labelnames +
           ' ' + opj(labels_dir, '_'.join([animal, 'seg-' + BALSAname, 'dseg.nii.gz'])) + ' -drop-unused-labels')
    run_cmd.wb(cmd, diary_name)

    dictionary = {"Sources": [opj(dir_native_resol, '.'.join([animal, BALSAname, 'dlabel', 'nii'])),
                              Ref_file],
                  "Description": 'Get a label volume (wb_command and Antspy)', }
    json_object = json.dumps(dictionary, indent=2)

    with open(opj(labels_dir, '_'.join([animal, 'seg-' + BALSAname, 'dseg.json'])), "w") as outfile:
        outfile.write(json_object)

    os.remove(opj(volumes_dir, '_'.join([animal, 'l', BALSAname + '.nii.gz'])))
    os.remove(opj(volumes_dir, '_'.join([animal, 'r', BALSAname + '.nii.gz'])))

    nl = 'volumetric atlas of yerkes: done!'
    run_cmd.msg(nl, diary_name, 'OKGREEN')