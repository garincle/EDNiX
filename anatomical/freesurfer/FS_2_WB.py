#import
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

from Tools import run_cmd, get_orientation
from Tools import getpath

from anatomical.connectomeWB import WB_nuclei, WB_ribbon, WB_label, WB_vol_atlas, WB_sphere_reg
from anatomical import norm2template

# Constant
SurfaceSmoothingFWHM  = 2
Hmin        = ['l','r']
Hcap        = ['L','R']
Hsurf       = ['left','right']
CORTEX      = ['CORTEX_LEFT','CORTEX_RIGHT']
RibbonValue = [[2,3], [41, 42]]

def WB_prep(cmd_mris_convert,FS_dir,FS_refs,dir_path,animal,species,spacing,resamp,datfile,reference,list_atlas,balsa_folder,BALSAname,balsa_brainT1,
            path_label_code,targetsuffix,listTimage,do_Nu,diary_name,sing_fs,sing_wb,sing_afni,export_fs, list_transfo,type_norm,iso):

    #  spec file used by Workbench %
    #  inspired from FreeSurfer2CaretConvertAndRegisterNonlinear.sh within the HCP pipeline and the NHP pipeline
    # 2023.10.10 Simon Clavagnier made modifications to fit python
    # 2025.03.28 Simon Clavagnier made modifications to fit any species
    # 2025.06 : improved flexibility

    nl = 'Creation of the ' + animal + ' spec file template'
    run_cmd.msg(nl, diary_name,'HEADER')

    ###################################################################################################################
    # organization
    #    0.0 settings

    (path_anat, matrices_dir, _, prepro_dir, dir_native, volumes_dir, labels_dir, masks_dir, wb_template_dir, wb_template_vol, wb_template_labels, wb_template_masks,
    wb_balsa_dir, wb_balsa_vol, wb_balsa_labels, wb_balsa_masks) = getpath.anat(dir_path,reference,BALSAname, '_','_','template')

    dir_native_resol, dir_native_32, dir_balsa_resol, dir_balsa_32, dir_balsa_64 = getpath.surf(path_anat,reference,BALSAname)
    if iso == 0:
        Ref_file = opj(volumes_dir, '_'.join([animal,targetsuffix,'res-iso',type_norm]) + '.nii.gz')
    elif iso == 1:
        Ref_file = opj(volumes_dir, '_'.join([animal, targetsuffix, type_norm]) + '.nii.gz')
    else :
        print('sorry I cannot find if your image is isotropic of not')
    print(iso)
    print(Ref_file)
    print(ope(Ref_file))
    _, _, _, bckFS_cmd, _ = get_orientation.use_ants(Ref_file, sing_fs)

    conv_voxmm = 0.8
    towarp     = 0

    if ope(dir_native_resol) == False:
        os.makedirs(dir_native_resol)

    if not BALSAname == '':
        towarp = 1
        if species == 'Human':
            conv_voxmm = 0.75
            atlasfile = 'Atlas_ROIs.1.60.nii.gz'
        elif species == 'Macaque':
            conv_voxmm = 0.75
            atlasfile = 'Atlas_ROIs.0.5.nii.gz'
        elif species == 'Chimpanzee':
            conv_voxmm = 0.8
            atlasfile = 'Atlas_ROIs.0.8.nii.gz'

        balsa_label     = opj(path_label_code, BALSAname + '_label.txt')

        if 'T2' in type_norm:
            suffix = 'T2w'
        else:
            suffix = 'T1w'
        Ref_norm        = opj(wb_balsa_vol, '_'.join([animal,'space-' + BALSAname,'desc-SS',suffix]) + '.nii.gz')
        ref_transfo1    = opj(matrices_dir, 'standard_' + BALSAname + '.nii.gz')
        ref_transfo2    = opj(matrices_dir, 'affine_' + BALSAname + '.nii.gz')

        if not ope(dir_native_32):
            os.makedirs(dir_native_32)
        if not ope(wb_balsa_vol):
            os.makedirs(wb_balsa_vol)
        if not ope(wb_balsa_labels):
            os.makedirs(wb_balsa_labels)
        if not ope(wb_balsa_masks):
            os.makedirs(wb_balsa_masks)
        if not ope(dir_balsa_resol):
            os.makedirs(dir_balsa_resol)
        if not ope(dir_balsa_64):
            os.makedirs(dir_balsa_64)
        if not ope(dir_balsa_32):
            os.makedirs(dir_balsa_32)

        if not opi(Ref_norm):
            if not ope(matrices_dir):
                os.makedirs(matrices_dir)
            refnb = 0
            for i, j in enumerate(list_transfo):
                if list_transfo[i]["name"] == 'coreg':
                    refnb = i

            # important variables : check if it is ok .........................................................................
            targetname  = 'native'
            transfoname = [targetname, 'to', BALSAname]
            transfoT    = 'Shift'
            transfoS    = 'Final'
            suffix      = '_desc-SS_T1w'
            transfonameR = list(reversed(transfoname))
            transfonameT = '_'.join(transfoname + [transfoT])
            transfonameS = '_'.join(transfoname + [transfoS])


            norm2template.norm(animal,Ref_file,'',wb_balsa_vol,BALSAname,balsa_brainT1,'',
                               matrices_dir,list_transfo[refnb]["type_of_transform"],transfonameT,transfonameS,
                             list_transfo[refnb]["affmetricT"], list_transfo[refnb]["affmetric"], list_transfo[refnb]["interpol"],diary_name,sing_wb,suffix,1)


    label_FS_all = ' ' + opj(FS_refs,'FreeSurferAllLut.txt') + ' '
    label_FS_Sc  = ' ' + opj(FS_refs,'FreeSurferSubcorticalLabelTableLut.txt') + ' '
    if opi(opj(dir_native_resol, animal + '_native_LR.wb.spec')):
        os.remove(opj(dir_native_resol, animal + '_native_LR.wb.spec'))

    ################################################################################################################
    # Get the transformation matrices to convert from Freesurfer to wb
    nl = 'calculate the transform matrix'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    cmd = sing_fs + 'mri_convert -odt float ' + Ref_file + ' ' + opj(volumes_dir,'tmp.nii.gz')
    run_cmd.run(cmd,diary_name)

    cras = spco(sing_fs + 'mri_info --cras ' + opj(volumes_dir,'tmp.nii.gz'), shell=True).decode('ascii').split('\n')
    os.remove(opj(volumes_dir,'tmp.nii.gz'))
    mm = cras[0].split(' ')
    ras = [[[1, 0, 0, float(mm[0])], [0, 1, 0, float(mm[1])], [0, 0, 1, float(mm[2])]],
          [[-1, 0, 0, float(mm[0])], [0, -1, 0, float(mm[1])], [0, 0, 1, float(mm[2])]]]
    cras_name = ['c_ras','c_ras_r']
    for h in range(2):
        with open(opj(volumes_dir, cras_name[h] + '.txt'), 'wb') as f:
            np.savetxt(f, ras[h], fmt='%.4f')
        shutil.move(opj(volumes_dir,cras_name[h] + '.txt'), 
                    opj(volumes_dir,cras_name[h] + '.mat'))


    ################################################################################################################

    #  1.0 Copy Files
    if towarp == 1:
        WB_sphere_reg.copyref(FS_refs, dir_balsa_64, dir_balsa_32, wb_balsa_labels, atlasfile, animal, diary_name, sing_wb)

    # 1.1a get the volumes
    # 1.1b get subcortical structures...

    cmd = sing_fs + 'mri_convert ' + bckFS_cmd  + opj(FS_dir,animal,'mri','aseg.mgz') + ' ' + opj(labels_dir,animal + '_seg-wmparc_dseg.nii.gz')
    run_cmd.run(cmd,diary_name)

    if opi(opj(FS_dir,animal,'mri','segmentation.mgz')) == True:
        cmd = sing_fs + 'mri_convert ' + bckFS_cmd  + opj(FS_dir,animal,'mri','segmentation.mgz') + \
        ' ' + opj(labels_dir,animal + '_label-3-FS_dseg.nii.gz')
        run_cmd.run(cmd,diary_name)

    if opi(datfile) == True:
        run_cmd.msg(datfile, diary_name, 'OKGREEN')
        aseg      = ants.image_read(opj(labels_dir,animal + '_seg-wmparc_dseg.nii.gz'))
        aseg_new  = ants.image_clone(aseg)
        brain_img = ants.image_read(Ref_file)

        S = brain_img.spacing
        X = aseg.origin
        new_orig = np.zeros(3)
        new_orig[0] = X[0] / spacing[0] * S[0]
        new_orig[1] = X[1] / spacing[1] * S[1]
        new_orig[2] = X[2] / spacing[2] * S[2]

        ants.set_origin(aseg_new, [new_orig[0], new_orig[1], new_orig[2]])
        if resamp ==1:
            aseg_new = ants.resample_image(aseg_new, [spacing[0], spacing[1], spacing[2]], False, 1)
            aseg_new[aseg_new < 0] = 0
        ants.set_spacing(aseg_new, S)
        ants.image_write(aseg_new, opj(labels_dir, animal + '_seg-wmparc_dseg.nii.gz'))
        cmd = (sing_afni + '3dresample -overwrite -master ' + Ref_file +
               ' -prefix ' + opj(labels_dir,animal + '_seg-wmparc_dseg.nii.gz') +
               ' -input ' + opj(labels_dir,animal + '_seg-wmparc_dseg.nii.gz'))
        run_cmd.run(cmd,diary_name)

        if opi(opj(FS_dir,animal,'mri','segmentation.mgz')) == True:
            mseg = ants.image_read(opj(labels_dir,animal + '_label-3-FS_dseg.nii.gz'))
            mseg_new = ants.image_clone(mseg)
            ants.set_origin(mseg_new, [new_orig[0], new_orig[1], new_orig[2]])
            mseg_new = ants.resample_image(mseg_new, [spacing[0], spacing[1], spacing[2]], False, 1)
            aseg_new[mseg_new < 0] = 0
            ants.set_spacing(mseg_new, S)
            ants.image_write(mseg_new, opj(labels_dir,animal + '_label-3-FS_dseg.nii.gz'))
            cmd = (sing_afni + '3dresample -overwrite -master ' + Ref_file +
                   ' -prefix ' + opj(labels_dir,animal + '_label-3-FS_dseg.nii.gz') +
                   ' -input ' + opj(labels_dir,animal + '_label-3-FS_dseg.nii.gz'))
            run_cmd.run(cmd,diary_name)

            dictionary = {"Sources": opj(FS_dir,animal,'mri','segmentation.mgz'),
                          "Description": 'manually improved tissue classification (Freeview)', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(labels_dir,animal + '_label-3-FS_dseg.json'), "w") as outfile:
                outfile.write(json_object)


    cmd = (sing_wb + 'wb_command -volume-label-import ' + opj(labels_dir,animal + '_seg-wmparc_dseg.nii.gz') +
           label_FS_all + opj(labels_dir,animal + '_seg-wmparc_dseg.nii.gz') + ' -drop-unused-labels')
    run_cmd.wb(cmd,diary_name)
    
    cmd = (sing_wb + 'wb_command -volume-label-import ' + opj(labels_dir,animal + '_seg-wmparc_dseg.nii.gz') +
           label_FS_Sc + opj(labels_dir,animal + '_seg-sc_dseg.nii.gz') + ' -discard-others')
    run_cmd.wb(cmd,diary_name)

    if ope(opj(FS_dir,animal,'mri','segmentation.nii.gz')) == True:
        cmd = sing_fs + 'mri_convert ' + bckFS_cmd  + opj(FS_dir,animal,'mri','segmentation.nii.gz') + \
        ' ' + opj(labels_dir,animal + '-label-3-FS_dseg.nii.gz')
        run_cmd.run(cmd,diary_name)


    #####################################################################################################################################################################################################################
    #  2 Proceed with the surfaces

    nl = 'Creation of the surface files\n'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    for h in range(2):

        labels = glob.glob(opj(FS_dir, animal, 'label', Hmin[h] + 'h.' + animal + '_*.annot'))

        for FS_surf, WB_surf, S_type in zip(['white', 'pial', 'mid'],
                                            ['white', 'pial', 'midthickness'],
                                            ['GRAY_WHITE', 'PIAL', 'MIDTHICKNESS']):  
            
            if opi(datfile) == True:
                cmd =(export_fs + 'mri_surf2surf --s ' + animal  + ' --sval-xyz ' + FS_surf + ' --reg-inv ' + datfile +
                      ' --tval-xyz ' + opj(FS_dir,animal,'mri','orig.mgz') + ' --hemi ' + Hmin[h] + 'h --tval ' + opj(FS_dir,animal,'surf', Hmin[h] +'h.' + FS_surf + '_temp'))
                run_cmd.do(cmd,diary_name)
                surfname = opj(FS_dir,animal,'surf', Hmin[h] +'h.'+ FS_surf + '_temp')
            else:
                surfname = opj(FS_dir, animal, 'surf', Hmin[h] + 'h.'+ FS_surf)
            
            cmd = (sing_fs + cmd_mris_convert + ' '  + surfname + ' ' + opj(dir_native_resol,animal + '.' + Hmin[h] + '.' + WB_surf + '.surf.gii'))
            run_cmd.run(cmd,diary_name)
            
            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol,animal + '.' + Hmin[h] + '.' + WB_surf + '.surf.gii') +
                   ' ' + CORTEX[h] + ' -surface-type ANATOMICAL -surface-secondary-type ' + S_type)
            run_cmd.run(cmd,diary_name)
            
            cmd = (sing_wb + 'wb_command -surface-apply-affine ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + WB_surf + '.surf.gii') + ' ' + opj(volumes_dir,'c_ras.mat') + ' ' +
                   opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + WB_surf + '.surf.gii'))
            run_cmd.run(cmd,diary_name)
            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') + ' ' + CORTEX[h] +
                   ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + WB_surf + '.surf.gii'))
            run_cmd.run(cmd,diary_name)

        # FS flat map (only for Native spec file)

        if opi(opj(FS_dir, animal, 'surf', Hmin[h] + 'h.full.flatten.patch.3d')) == True:
            cmd = (sing_fs + cmd_mris_convert + ' -p ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.full.flatten.patch.3d') +
                    ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.flat.surf.gii'))
            run_cmd.run(cmd,diary_name)
            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.flat.surf.gii') +
                    ' ' + CORTEX[h] + ' -surface-type FLAT -surface-secondary-type GRAY_WHITE')
            run_cmd.run(cmd,diary_name)
            cmd = (sing_wb + 'wb_command -surface-apply-affine ' + opj(dir_native_resol,
                                                                 animal + '.' + Hmin[h] + '.flat.surf.gii') +
                    ' ' + opj(volumes_dir, cras_name[h] + '.mat') + ' ' + opj(dir_native_resol,
                                                                                              animal + '.' + Hmin[
                                                                                                  h] + '.flat.surf.gii'))
            run_cmd.run(cmd,diary_name)
            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') + ' ' +
                       CORTEX[h] +
                       ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.flat.surf.gii'))
            run_cmd.run(cmd,diary_name)


        # from midthickness generation of inflated and very inflated surface
        cmd = (sing_wb + 'wb_command -file-information ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii'))
        #      ' | grep "Number of Vertices:" | cut -f2 -d: | tr -d "[:space:]"'
        out,_=run_cmd.run(cmd,diary_name)
        result = out.decode("utf-8").split('\n')
        nb=1
        for item in result:
            if item.find("Number of Vertices:") != -1:
                nb = int(item.split(' ')[-1])
        NativeInflationScale = nb * conv_voxmm / 32492  # 32492 (so 32k) for human or MY19 low res or  74k for F99

        cmd = sing_wb + 'wb_command -surface-generate-inflated ' + \
              opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') + \
              ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.inflated.surf.gii') + \
              ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.very_inflated.surf.gii') + \
              ' -iterations-scale ' + str(NativeInflationScale)
        run_cmd.run(cmd,diary_name)

        cmd = sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') + ' ' + CORTEX[h] + \
            ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.inflated.surf.gii')
        run_cmd.run(cmd,diary_name)

        cmd = sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') + ' ' + CORTEX[h] + \
            ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.very_inflated.surf.gii')
        run_cmd.run(cmd,diary_name)


        # sulc thickness and curvature for native brain

        for SHAPE_FS,SHAPE_WB in zip(['sulc', 'curv'],
                                     ['sulc', 'curvature']):

            cmd = (sing_fs + cmd_mris_convert + ' -c ' + opj(FS_dir,animal,'surf', Hmin[h] + 'h.' + SHAPE_FS) +
                ' ' + opj(FS_dir,animal,'surf',Hmin[h] + 'h.white') +
                ' ' + opj(dir_native_resol,animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii'))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol,animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
               ' ' + CORTEX[h])
            run_cmd.run(cmd,diary_name)
            cmd = (sing_wb + 'wb_command -metric-math "var * -1" ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
               ' -var var ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii'))
            run_cmd.do(cmd,diary_name)
            cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
               ' -map 1 ' + animal + '_' + Hmin[h] + '_' + SHAPE_WB)
            run_cmd.run(cmd,diary_name)
            cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
               ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true')
            run_cmd.run(cmd,diary_name)

        # extract thickness
        cmd = (sing_fs + cmd_mris_convert + ' -c ' + opj(FS_dir,animal,'surf', Hmin[h] + 'h.thickness') +
               ' ' + opj(FS_dir,animal,'surf', Hmin[h] + 'h.white') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') + ' ' + CORTEX[h])
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-math "var * -1" ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' -var var ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' -map 1 ' + animal + '_' + Hmin[h] + '_Thickness')
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true')
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-math "abs(thickness)" ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' -var thickness ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false')
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-math "thickness > 0" ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' -var thickness ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-fill-holes ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-remove-islands ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' -map 1 ' + animal + '_' + Hmin[h] + '_ROI')
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-dilate ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' 10 ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') + ' -nearest')
        run_cmd.run(cmd,diary_name)
        cmd = (sing_wb + 'wb_command -metric-dilate ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.curvature.shape.gii') +
               ' ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
               ' 10 ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.curvature.shape.gii') + ' -nearest')
        run_cmd.run(cmd,diary_name)

        # parcellations (.annot)

        for j in labels:
            name = opb(j).split('.')[1].split('_')[-1]
            print(name)
            WB_label.surfWB(cmd_mris_convert, FS_dir, animal, name, opj(dir_native_resol), '.native', h, diary_name,
                            sing_fs, sing_wb)

        for a in range(len(list_atlas[0])):
            if list_atlas[3][a] == 1:
                if list_atlas[2][a] > 1:
                    WB_label.mergelabel(animal, opj(dir_native_resol), '.native', h, diary_name, list_atlas[0][a], list_atlas[2][a], sing_wb)

    if towarp == 1:
        # warping of the surfaces
        WB_sphere_reg.warp(animal, ['white', 'pial', 'midthickness'], dir_native_resol, dir_balsa_resol, list_atlas, ref_transfo1,
                           ref_transfo2, conv_voxmm, cmd_mris_convert, FS_dir, diary_name, sing_wb, sing_fs)
        # sphere and registration
        WB_sphere_reg.sphere(cmd_mris_convert, animal, FS_dir, dir_balsa_resol, diary_name, sing_fs, sing_wb)
        # sphere to sphere registration
        WB_sphere_reg.reg(animal, dir_balsa_resol, dir_balsa_64, diary_name, sing_wb)
        # resampling in high resolution : 164k
        WB_sphere_reg.surf_resamp164(animal, dir_native_resol, dir_balsa_resol, dir_balsa_64, conv_voxmm, FS_refs, list_atlas, diary_name, sing_wb)
        # resampling in low resolution : 32k
        WB_sphere_reg.surf_resamp32(animal, dir_native_resol, dir_balsa_resol, dir_balsa_32, conv_voxmm, FS_refs, list_atlas, diary_name, sing_wb)
        # now back to native space but with the low resolution (32k)
        WB_sphere_reg.surf_resamp_native(animal, dir_native_resol, dir_balsa_resol, dir_native_32, dir_balsa_32, conv_voxmm,
                                         diary_name, sing_wb)

    nl = 'creation of the surfaces and labels files: done!\n'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    ###################################################################################################################

    ##   3 Create CIFTI once both hemispheres have been proceed

    nl = 'Creation of the cifti files\n'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    # in wb_native_dir,'Native',
    
    cmd = (sing_wb + 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol,animal + '.sulc.dscalar.nii') +
           ' -left-metric ' + opj(dir_native_resol,animal + '.l.sulc.shape.gii') +
           ' -right-metric ' + opj(dir_native_resol,animal + '.r.sulc.shape.gii'))
    run_cmd.run(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol,animal + '.sulc.dscalar.nii') + ' -map 1 ' + animal + '_Sulc')
    run_cmd.run(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -cifti-palette ' + opj(dir_native_resol,animal + '.sulc.dscalar.nii') +
           ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol,animal + '.sulc.dscalar.nii') +
           ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true')
    run_cmd.run(cmd,diary_name)

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') +
           ' INVALID ' + opj(dir_native_resol, animal + '.sulc.dscalar.nii'))
    run_cmd.run(cmd,diary_name)

    for surf,palette in zip(['curvature','thickness'],
                            [' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true',
                             ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false']):

        cmd = (sing_wb + 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol,animal + '.' + surf + '.dscalar.nii') +
               ' -left-metric '  + opj(dir_native_resol,animal + '.l.' + surf + '.shape.gii') +
               ' -roi-left '  + opj(dir_native_resol,animal + '.l.roi.shape.gii') +
               ' -right-metric ' + opj(dir_native_resol,animal + '.r.' + surf + '.shape.gii') +
               ' -roi-right ' + opj(dir_native_resol,animal + '.r.roi.shape.gii'))
        run_cmd.run(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_native_resol,animal + '.' + surf + '.dscalar.nii') + ' -map 1 ' + animal + '_' + surf)
        run_cmd.run(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -cifti-palette ' + opj(dir_native_resol,animal + '.' + surf + '.dscalar.nii') +
               ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol,animal + '.' + surf + '.dscalar.nii') + palette)
        run_cmd.run(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') +
               ' INVALID ' + opj(dir_native_resol, animal + '.' + surf + '.dscalar.nii'))
        run_cmd.run(cmd,diary_name)

    labels = glob.glob(opj(dir_native_resol, animal + '.l.*.label.gii'))

    for j in labels:

        name = opb(j).split('.')[2]
        WB_label.cifti(animal, name, opj(dir_native_resol), '.native', 'roi',
                       opj(dir_native_resol, animal + '_native_LR.wb.spec'), diary_name, sing_wb)
        if iso == 1:
            atlasfile = opj(labels_dir, animal + '_seg-' + name + '_res-iso_dseg.nii.gz')
        else:
            atlasfile = opj(labels_dir, animal + '_seg-' + name + '_dseg.nii.gz')
        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') +
               ' INVALID ' + atlasfile)
        run_cmd.run(cmd,diary_name)

    for i in listTimage:
        if iso == 0:
            imfile = opj(volumes_dir, '_'.join([animal, targetsuffix, 'res-iso', i]) + '.nii.gz')
        else:
            imfile = opj(volumes_dir, '_'.join([animal, targetsuffix, i]) + '.nii.gz')

    cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal + '_native_LR.wb.spec') +
           ' INVALID ' + imfile)
    run_cmd.run(cmd,diary_name)


    if towarp ==1 :
        WB_sphere_reg.cifti(animal, dir_balsa_resol, dir_balsa_32, dir_balsa_64, diary_name, sing_wb)


    nl ='creation of the cifti files: done!'
    run_cmd.msg(nl, diary_name, 'OKGREEN')

    #                 2.0 ribbon              #########################################################"
    #     for native space:  + supra / infra masks
    WB_ribbon.create(dir_path, Ref_file, label_FS_all, animal, BALSAname, 'native', Hmin, RibbonValue, diary_name, sing_wb)

    if towarp == 1:
        #     for template space
        WB_ribbon.create(dir_path, Ref_norm, label_FS_all, animal, BALSAname, 'template', Hmin, RibbonValue, diary_name, sing_wb)
    

    # extras ###########################################################################################

    if species == 'Macaque':
        # parcellation from template surfaces (resolution = 32k) to native space volume
        WB_vol_atlas.MY19(cmd_mris_convert, animal, FS_dir, balsa_folder, dir_path, Ref_file, path_label_code,
                          Hmin, CORTEX, diary_name, sing_wb, sing_fs)

    # get the gray nuclei surfaces
    if do_Nu[0] ==1:
        WB_nuclei.native(cmd_mris_convert, FS_dir, dir_native_resol, animal, species, do_Nu[1], diary_name, sing_fs, sing_wb)
        if towarp == 1:
            WB_nuclei.template(dir_native_resol, dir_balsa_resol, animal, ref_transfo1, ref_transfo2, diary_name, sing_wb)

    nl = 'Done,congratulations!'
    run_cmd.msg(nl, diary_name,'OKGREEN')


        
