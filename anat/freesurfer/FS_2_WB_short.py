# import
import os
import subprocess
import glob
import shutil
import numpy as np
import ants
import json

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile
spco = subprocess.check_output

from Tools import run_cmd

from anat.connectomeWB import WB_ribbon, WB_label, WB_vol_atlas, WB_sphere_reg

# Constant
SurfaceSmoothingFWHM = 2
Hmin = ['l', 'r']
Hcap = ['L', 'R']
Hsurf = ['left', 'right']
CORTEX = ['CORTEX_LEFT', 'CORTEX_RIGHT']
RibbonValue = [[2, 3], [41, 42]]
cras_name = ['c_ras', 'c_ras_r']


def WB_prep(cmd_mris_convert, FS_dir, FS_refs, animal, species, Ref_file, datfile, reference, list_atlas, balsa_folder,
            BALSAname,
            path_label_code, template_dir, proj, diary_name, sing_fs, sing_wb, export_fs):
    label_FS_all = ' ' + opj(FS_refs, 'FreeSurferAllLut.txt') + ' '

    conv_voxmm = 0.8
    if species == 'Human':
        conv_voxmm = 0.75
    elif species == 'Macaque':
        conv_voxmm = 0.75
    elif species == 'Chimpanzee':
        conv_voxmm = 0.8

    vol_dir = opj(template_dir, reference,'volumes')
    labels_dir = opj(vol_dir, 'labels')
    masks_dir = opj(vol_dir, 'masks')
    dir_native_resol = opj(template_dir, reference, 'surfaces', 'Native_resol')

    ################################################################################################################
    # Get the transformation matrices to convert from Freesurfer to wb
    nl = 'calculate the transform matrix'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    cmd = sing_fs + 'mri_convert -odt float ' + Ref_file + ' ' + opj(vol_dir, 'tmp.nii.gz')
    run_cmd.run(cmd, diary_name)

    cras = spco(sing_fs + 'mri_info --center ' + opj(vol_dir, 'tmp.nii.gz'), shell=True).decode('ascii').split('\n')
    os.remove(opj(vol_dir, 'tmp.nii.gz'))
    mm = cras[0].split(' ')
    ras = [[[1, 0, 0, float(mm[0])], [0, 1, 0, float(mm[1])], [0, 0, 1, float(mm[2])]],
           [[-1, 0, 0, float(mm[0])], [0, -1, 0, float(mm[1])], [0, 0, 1, float(mm[2])]]]
    cras_name = ['c_ras', 'c_ras_r']
    for h in range(2):
        print(opj(vol_dir, cras_name[h] + '.mat'))
        with open(opj(vol_dir, cras_name[h] + '.txt'), 'wb') as f:
            np.savetxt(f, ras[h], fmt='%.10f')
        shutil.move(opj(vol_dir, cras_name[h] + '.txt'),
                    opj(vol_dir, cras_name[h] + '.mat'))

    #####################################################################################################################################################################################################################
    #  2 Proceed with the surfaces

    nl = 'Creation of the surface files\n'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    for h in range(2):

        for FS_surf, WB_surf, S_type in zip(['white', 'pial', 'mid'],
                                            ['white', 'pial', 'midthickness'],
                                            ['GRAY_WHITE', 'PIAL', 'MIDTHICKNESS']):
            if opi(datfile) == True:
                cmd = (export_fs + 'mri_surf2surf --s ' + animal + ' --sval-xyz ' + FS_surf + ' --reg-inv ' + datfile +
                       ' --tval-xyz ' + opj(FS_dir, animal, 'mri', 'orig_raw.mgz') + ' --hemi ' + Hmin[
                           h] + 'h --tval ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.' + FS_surf + '_temp'))
                run_cmd.do(cmd, diary_name)
                surfname = opj(FS_dir, animal, 'surf', Hmin[h] + 'h.' + FS_surf + '_temp')
            else:
                surfname = opj(FS_dir, animal, 'surf', Hmin[h] + 'h.' + FS_surf)

            cmd = (sing_fs + cmd_mris_convert + ' --to-scanner ' + surfname + ' ' + opj(dir_native_resol, animal + '.' + Hmin[
                h] + '.' + WB_surf + '.surf.gii'))
            run_cmd.run(cmd, diary_name)

            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + WB_surf + '.surf.gii') +
                   ' ' + CORTEX[h] + ' -surface-type ANATOMICAL -surface-secondary-type ' + S_type)
            run_cmd.run(cmd, diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                                   animal + '_native_LR.wb.spec') + ' ' + CORTEX[h] +
                   ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + WB_surf + '.surf.gii'))
            run_cmd.run(cmd, diary_name)

        # FS flat map (only for Native spec file)
        if opi(opj(FS_dir, animal, 'surf', Hmin[h] + 'h.full.flatten.patch.3d')) == True:
            cmd = (sing_fs + cmd_mris_convert + ' -p ' + opj(FS_dir, animal, 'surf',
                                                             Hmin[h] + 'h.full.flatten.patch.3d') +
                   ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.flat.surf.gii'))
            run_cmd.run(cmd, diary_name)
            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol,
                                                                animal + '.' + Hmin[h] + '.flat.surf.gii') +
                   ' ' + CORTEX[h] + ' -surface-type FLAT -surface-secondary-type GRAY_WHITE')
            run_cmd.run(cmd, diary_name)
            cmd = (sing_wb + 'wb_command -surface-apply-affine ' + opj(dir_native_resol,
                                                                       animal + '.' + Hmin[h] + '.flat.surf.gii') +
                   ' ' + opj(vol_dir, cras_name[h] + '.mat') + ' ' + opj(dir_native_resol,
                                                                         animal + '.' + Hmin[
                                                                             h] + '.flat.surf.gii'))
            run_cmd.run(cmd, diary_name)
            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                                   animal + '_native_LR.wb.spec') + ' ' +
                   CORTEX[h] +
                   ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.flat.surf.gii'))
            run_cmd.run(cmd, diary_name)

        # from midthickness generation of inflated and very inflated surface
        cmd = (sing_wb + 'wb_command -file-information ' + opj(dir_native_resol,
                                                               animal + '.' + Hmin[h] + '.midthickness.surf.gii'))
        #      ' | grep "Number of Vertices:" | cut -f2 -d: | tr -d "[:space:]"'
        out, _ = run_cmd.run(cmd, diary_name)
        result = out.decode("utf-8").split('\n')
        nb = 1
        for item in result:
            if item.find("Number of Vertices:") != -1:
                nb = int(item.split(' ')[-1])
        NativeInflationScale = nb * conv_voxmm / 32492  # 32492 (so 32k) for human or MY19 low res or  74k for F99

        cmd = sing_wb + 'wb_command -surface-generate-inflated ' + \
              opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') + \
              ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.inflated.surf.gii') + \
              ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.very_inflated.surf.gii') + \
              ' -iterations-scale ' + str(NativeInflationScale)
        run_cmd.run(cmd, diary_name)

        cmd = sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') + ' ' + \
              CORTEX[h] + \
              ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.inflated.surf.gii')
        run_cmd.run(cmd, diary_name)

        cmd = sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') + ' ' + \
              CORTEX[h] + \
              ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.very_inflated.surf.gii')
        run_cmd.run(cmd, diary_name)

        # sulc thickness and curvature for native brain

        for SHAPE_FS, SHAPE_WB in zip(['sulc', 'curv'],
                                      ['sulc', 'curvature']):
            cmd = (sing_fs + cmd_mris_convert + ' -c ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.' + SHAPE_FS) +
                   ' ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.white') +
                   ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii'))
            run_cmd.run(cmd, diary_name)

            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol, animal + '.' + Hmin[
                h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' ' + CORTEX[h])
            run_cmd.run(cmd, diary_name)
            cmd = (sing_wb + 'wb_command -metric-math "var * -1" ' + opj(dir_native_resol, animal + '.' + Hmin[
                h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' -var var ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii'))
            run_cmd.do(cmd, diary_name)
            cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol, animal + '.' + Hmin[
                h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' -map 1 ' + animal + '_' + Hmin[h] + '_' + SHAPE_WB)
            run_cmd.run(cmd, diary_name)
            cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_native_resol, animal + '.' + Hmin[
                h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true')
            run_cmd.run(cmd, diary_name)

        # extract thickness
        cmd = (sing_fs + cmd_mris_convert + ' -c ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.thickness') +
               ' ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.white') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol,
                                                            animal + '.' + Hmin[h] + '.thickness.shape.gii') + ' ' +
               CORTEX[h])
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-math "var * -1" ' + opj(dir_native_resol,
                                                                     animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' -var var ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol,
                                                            animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' -map 1 ' + animal + '_' + Hmin[h] + '_Thickness')
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_native_resol,
                                                             animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true')
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-math "abs(thickness)" ' + opj(dir_native_resol, animal + '.' + Hmin[
            h] + '.thickness.shape.gii') +
               ' -var thickness ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_native_resol,
                                                             animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false')
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-math "thickness > 0" ' + opj(dir_native_resol,
                                                                          animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' -var thickness ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd, diary_name)

        cmd = (sing_wb + 'wb_command -metric-fill-holes ' + opj(dir_native_resol,
                                                                animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-remove-islands ' + opj(dir_native_resol,
                                                                    animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol,
                                                            animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' -map 1 ' + animal + '_' + Hmin[h] + '_ROI')
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-dilate ' + opj(dir_native_resol,
                                                            animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' 10 ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') + ' -nearest')
        run_cmd.run(cmd, diary_name)
        cmd = (sing_wb + 'wb_command -metric-dilate ' + opj(dir_native_resol,
                                                            animal + '.' + Hmin[h] + '.curvature.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' 10 ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.curvature.shape.gii') + ' -nearest')
        run_cmd.run(cmd, diary_name)

    nl = 'creation of the surfaces and labels files: done!\n'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    ###################################################################################################################

    ##   3 Create CIFTI once both hemispheres have been proceed

    nl = 'Creation of the cifti files\n'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    # in wb_native_dir,'Native',

    cmd = (sing_wb + 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol, animal + '.sulc.dscalar.nii') +
           ' -left-metric ' + opj(dir_native_resol, animal + '.l.sulc.shape.gii') +
           ' -right-metric ' + opj(dir_native_resol, animal + '.r.sulc.shape.gii'))
    run_cmd.run(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol,
                                                        animal + '.sulc.dscalar.nii') + ' -map 1 ' + animal + '_Sulc')
    run_cmd.run(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -cifti-palette ' + opj(dir_native_resol, animal + '.sulc.dscalar.nii') +
           ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal + '.sulc.dscalar.nii') +
           ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true')
    run_cmd.run(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') +
           ' INVALID ' + opj(dir_native_resol, animal + '.sulc.dscalar.nii'))
    run_cmd.run(cmd, diary_name)

    for surf, palette in zip(['curvature', 'thickness'],
                             [
                                 ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true',
                                 ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false']):
        cmd = (sing_wb + 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol,
                                                                        animal + '.' + surf + '.dscalar.nii') +
               ' -left-metric ' + opj(dir_native_resol, animal + '.l.' + surf + '.shape.gii') +
               ' -roi-left ' + opj(dir_native_resol, animal + '.l.roi.shape.gii') +
               ' -right-metric ' + opj(dir_native_resol, animal + '.r.' + surf + '.shape.gii') +
               ' -roi-right ' + opj(dir_native_resol, animal + '.r.roi.shape.gii'))
        run_cmd.run(cmd, diary_name)

        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol,
                                                           animal + '.' + surf + '.dscalar.nii') + ' -map 1 ' + animal + '_' + surf)
        run_cmd.run(cmd, diary_name)

        cmd = (sing_wb + 'wb_command -cifti-palette ' + opj(dir_native_resol, animal + '.' + surf + '.dscalar.nii') +
               ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal + '.' + surf + '.dscalar.nii') + palette)
        run_cmd.run(cmd, diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') +
               ' INVALID ' + opj(dir_native_resol, animal + '.' + surf + '.dscalar.nii'))
        run_cmd.run(cmd, diary_name)

    # parcellations
    WB_label.vol2surfWB(0,species, list_atlas, labels_dir, dir_native_resol, 'native', Hmin,
                        proj, 'roi', opj(dir_native_resol, animal + '_native_LR.wb.spec'), diary_name,
                        sing_wb)

    nl = 'creation of the cifti files: done!'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    #                 2.0 ribbon              #########################################################"
    #     for native space:  + supra / infra masks

    img = WB_ribbon.create(Ref_file, label_FS_all, animal, dir_native_resol, labels_dir, 'white',
                           ['pial', 'midthickness'], ['ribbon', 'infra'], Hmin, RibbonValue,
                           diary_name, sing_wb)

    WB_ribbon.layers(img, animal, dir_native_resol, labels_dir, masks_dir, 'white', ['pial', 'midthickness'], Ref_file,
                     diary_name, sing_wb)

    if not BALSAname == '':
        dir_native_32 = opj(template_dir, reference, 'surfaces', 'fsaverage_LR_32k')
        dir_balsa_resol = opj(template_dir, BALSAname, 'surfaces', 'native')
        dir_balsa_32 = opj(template_dir, BALSAname, 'surfaces', 'fsaverage_LR_32k')
        dir_balsa_64 = opj(template_dir, BALSAname, 'surfaces', 'fsaverage_164k')
        wb_balsa_vol = opj(template_dir, BALSAname, 'volumes')
        wb_balsa_masks = opj(wb_balsa_vol, 'masks')

        Ref_norm = opj(wb_balsa_vol, '_'.join([species, 'space-' + BALSAname, 'desc-SS', 'T1w.nii.gz']))
        ref_transfo1 = opj(balsa_folder, 'standard_' + BALSAname + '.nii.gz')
        ref_transfo2 = opj(balsa_folder, 'affine_' + BALSAname + '.nii.gz')

        # warping of the surfaces
        WB_sphere_reg.warp(animal, ['white', 'pial', 'midthickness'], dir_native_resol, dir_balsa_resol,
                           list_atlas, ref_transfo1,
                           ref_transfo2, conv_voxmm, cmd_mris_convert, FS_dir, 'WB', diary_name, sing_wb, sing_fs)

        for h in range(2):
            labels = sorted(glob.glob(opj(dir_native_resol,'.'.join([animal,Hmin[h],'*label','gii']))))
            print(labels)
            for l in range(len(labels)):
                shutil.copyfile(labels[l],opj(dir_balsa_resol,opb(labels[l])))
                print(labels[l])
                print(opj(dir_balsa_resol,opb(labels[l])))

        # sphere and registration
        WB_sphere_reg.sphere(cmd_mris_convert, animal, FS_dir, dir_balsa_resol, diary_name, sing_fs, sing_wb)
        # sphere to sphere registration
        WB_sphere_reg.reg(animal, dir_balsa_resol, dir_balsa_64, diary_name, sing_wb)
        # resampling in high resolution : 164k
        WB_sphere_reg.surf_resamp164(animal, dir_native_resol, dir_balsa_resol, dir_balsa_64, conv_voxmm,
                                 FS_refs, list_atlas, diary_name, sing_wb)
        # resampling in low resolution : 32k
        WB_sphere_reg.surf_resamp32(animal, dir_native_resol, dir_balsa_resol, dir_balsa_32, conv_voxmm,
                                FS_refs, list_atlas, diary_name, sing_wb)
        # now back to native space but with the low resolution (32k)
        WB_sphere_reg.surf_resamp_native(animal, dir_native_resol, dir_balsa_resol, dir_native_32, dir_balsa_32,
                                         conv_voxmm,
                                         diary_name, sing_wb)

        WB_sphere_reg.cifti(animal, dir_balsa_resol, dir_balsa_32, dir_balsa_64, diary_name, sing_wb)
        #     for template space
        WB_ribbon.create(Ref_norm, label_FS_all, animal,
                      dir_balsa_resol, wb_balsa_masks, 'white.native', ['pial.native'], ['ribbon'], Hmin, RibbonValue,
                      diary_name, sing_wb)

        # extras ###########################################################################################

        if species == 'Macaque':

            # parcellation from template surfaces to native space volume
            nl = 'Create ' + BALSAname + ' volumetric atlas'
            run_cmd.msg(nl, diary_name,'OKGREEN')

            labelfile = opj(path_label_code, BALSAname + '_label.txt')

            origFile, parcellation_list, ssamp = WB_vol_atlas.get(BALSAname,balsa_folder)

            WB_vol_atlas.set(origFile,animal, BALSAname, parcellation_list, labelfile, dir_native_resol,
                dir_native_32,
                dir_balsa_resol, dir_balsa_32, dir_balsa_64, vol_dir, labels_dir, FS_dir, Ref_file,
                ssamp, Hmin,
                CORTEX, diary_name, cmd_mris_convert, sing_wb, sing_fs)

    nl = 'Done,congratulations!'
    run_cmd.msg(nl, diary_name, 'OKGREEN')
