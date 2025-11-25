# import
import os
import glob
import ants
import json

opj  = os.path.join
opi = os.path.isfile

from Tools import run_cmd
from Tools import getpath

def create(dir_path,Ref_file,label_FS_all,animal,BALSAname,version,Hmin,RibbonValue,diary_name,sing_wb):

    nl = 'Creation of the ' + version + ' ribbon volume'
    run_cmd.msg(nl,diary_name,'HEADER')

    (path_anat, _, _, _, _, _, labels_dir, masks_dir, _,_, _, _,_, _,
     wb_balsa_labels, _) = getpath.anat(dir_path, BALSAname, BALSAname, '', '', 'template')

    dir_native_resol, _, dir_balsa_resol, _, _ = getpath.surf(path_anat, BALSAname,BALSAname)

    if version == 'native':
        dir_surf = dir_native_resol
        dir_out = labels_dir
        Surf_in  = 'white'
        Surf_out = ['pial', 'midthickness']
        vol_out  = ['ribbon', 'infra']
    elif version == 'template':
        dir_surf = dir_balsa_resol
        dir_out = wb_balsa_labels
        Surf_in  = 'white.native'
        Surf_out = ['pial.native']
        vol_out  = ['ribbon']
    else:
        print('Sorry, "version" should be set to "native" or "template" only !')
        return

    for h in range(2):
        cmd = (sing_wb + 'wb_command -create-signed-distance-volume ' +
               opj(dir_surf,animal + '.' + Hmin[h] + '.' + Surf_in + '.surf.gii') +
               ' ' + Ref_file +
               ' ' + opj(dir_surf, animal + '.' + Hmin[h] + '.' + Surf_in + '.nii.gz'))
        run_cmd.run(cmd,diary_name)

        img_w1 = ants.image_read(opj(dir_surf,animal + '.' + Hmin[h] + '.' + Surf_in + '.nii.gz'))
        img_w = ants.image_clone(img_w1)
        img_w[img_w<0]=0
        img_w[img_w>0]=1

        for s in range(len(Surf_out)):

            cmd = (sing_wb + 'wb_command -create-signed-distance-volume ' +
                   opj(dir_surf,animal + '.' + Hmin[h] + '.' + Surf_out[s] + '.surf.gii') +
                   ' ' + Ref_file +
                   ' ' + opj(dir_surf, animal + '.' + Hmin[h] + '.' + Surf_out[s] + '.nii.gz'))
            run_cmd.run(cmd,diary_name)

            img_u = ants.image_read(opj(dir_surf, animal + '.' + Hmin[h] + '.' + Surf_out[s] + '.nii.gz'))
            img_u[img_u > 0] = 0
            img_u[img_u < 0] = 1
            img_u = (img_u*img_w)*RibbonValue[h][1]

            img_r1 = ants.image_clone(img_w1)
            img_r1[img_r1>0]=0
            img_r1[img_r1<0] = RibbonValue[h][0]

            img_u = img_u + img_r1

            ants.image_write(img_u, opj(dir_surf, animal + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz'),
                             ri=False)

    nl = '# 1- Create the ribbon volume'
    run_cmd.msg(nl,diary_name,'OKGREEN')


    img1 = ants.image_read(opj(dir_surf, animal + '.l.ribbon.nii.gz'))
    img2 = ants.image_read(opj(dir_surf,animal + '.r.ribbon.nii.gz'))
    img1= img1+img2
    ants.image_write(img1, opj(dir_out, animal + '_seg-ribbon_dseg.nii.gz'),ri=False)

    dictionary = {"Sources": [opj(dir_surf, animal + '.lh.' + Surf_in + '.surf.gii'),
                              opj(dir_surf, animal + '.rh.' + Surf_in + '.surf.gii'),
                              opj(dir_surf, animal + '.lh.' + Surf_out[0] + '.surf.gii'),
                              opj(dir_surf, animal + '.rh.' + Surf_out[0] + '.surf.gii'),
                              Ref_file],
                  "Description": 'Get a white matter-cortical segmentation volume according to the surfaces', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_out, animal + '_seg-ribbon_dseg.json'), "w") as outfile:
        outfile.write(json_object)

    cmd = (sing_wb + 'wb_command -volume-label-import ' + opj(dir_out, animal + '_seg-ribbon_dseg.nii.gz') + label_FS_all +
           opj(dir_out, animal + '_seg-ribbon_dseg.nii.gz') + ' -drop-unused-labels')
    run_cmd.wb(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_surf, animal + '_native_LR.wb.spec') +
           ' INVALID ' + opj(dir_out, animal + '_seg-ribbon_dseg.nii.gz'))
    run_cmd.run(cmd,diary_name)


    if version == 'native':

        nl = '# 2- Create the cortical only version of the ribbon volume'
        run_cmd.msg(nl,diary_name,'OKGREEN')

        img_c = ants.mask_image(img1, img1, [3, 42])

        nl = '# 3- Create the infra-granular layer version of the ribbon volume'
        run_cmd.msg(nl,diary_name,'OKGREEN')

        img1 = ants.image_read(opj(dir_surf, animal + '.l.infra.nii.gz'))
        img2 = ants.image_read(opj(dir_surf, animal + '.r.infra.nii.gz'))
        img1 = img1 + img2
        img = ants.mask_image(img1,img1,[3,42])

        ants.image_write(img, opj(masks_dir, animal + '_desc-infra-layer_mask.nii.gz'), ri=False)

        dictionary = {"Sources": [opj(dir_surf, animal + '.lh.' + Surf_in + '.surf.gii'),
                                  opj(dir_surf, animal + '.rh.' + Surf_in + '.surf.gii'),
                                  opj(dir_surf, animal + '.lh.' + Surf_out[1] + '.surf.gii'),
                                  opj(dir_surf, animal + '.rh.' + Surf_out[1] + '.surf.gii'),
                                  Ref_file],
                      "Description": 'Get a "pseudo infra-granular layer" masks according tot the surfaces', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, animal + '_desc-infra-layer_mask.json'), "w") as outfile:
            outfile.write(json_object)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_surf, animal + '_native_LR.wb.spec') +
               ' INVALID ' + opj(masks_dir, animal + '_desc-infra-layer_mask.nii.gz'))
        run_cmd.run(cmd,diary_name)

        nl = '# 4- Create the supra-granular layer version of the ribbon volume'
        run_cmd.msg(nl,diary_name,'OKGREEN')

        img = img_c-img
        ants.image_write(img, opj(masks_dir, animal + '_desc-supra-layer_mask.nii.gz'), ri=False)

        dictionary = {"Sources": [opj(labels_dir, animal + '_seg-ribbon_dseg.nii.gz'),
                                  opj(masks_dir, animal + '_desc-infra-layer_mask.nii.gz')],
                      "Description": 'Get a "pseudo infra-granular layer" masks according tot the surfaces', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(masks_dir, animal + '_desc-supra-layer_mask.json'), "w") as outfile:
            outfile.write(json_object)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_surf, animal + '_native_LR.wb.spec') +
               ' INVALID ' + opj(masks_dir, animal + '_desc-supra-layer_mask.nii.gz'))
        run_cmd.run(cmd,diary_name)

    nii = glob.glob(opj(dir_surf, animal + '.*.nii.gz'))
    for i in nii:
        os.remove(i)

    nl = 'creation of the ribbon files: done!'
    run_cmd.msg(nl,diary_name,'OKGREEN')





