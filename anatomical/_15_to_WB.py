# import
import os
import subprocess
import shutil
import numpy as np

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

# Constant
conv_voxmm = 0.75  # 0.75 if FS to wb (1mm to 0.7mm volume)
SurfaceSmoothingFWHM = 2

Hmin = ['l', 'r']
Hcap = ['L', 'R']
CORTEX = ['CORTEX_LEFT', 'CORTEX_RIGHT']
RibbonValue = [[2, 3], [41, 42]]

from fonctions.extract_filename import extract_filename

def WB_prep(FS_dir, dir_native, animal_folder, Ref_file, species, list_atlases_2,s_bind,afni_sif,fsl_sif,fs_sif,wb_sif):

    cmd = 'singularity run' + s_bind + fs_sif + 'mri_info --orientation ' + Ref_file
    orient_raw = spgo(cmd)

    # purpouse for FS : LIA
    if orient_raw   == 'RAS': reorient = ' -r -1 3 -2 '
    elif orient_raw == 'LAS': reorient = ' -r 1 3 -2 '
    elif orient_raw == 'RSA': reorient = ' -r -1 -2 3'
    elif orient_raw == 'LSA': reorient = ' -r 1 -2 3 '
    elif orient_raw == 'RIA': reorient = ' -r -1 2 3 '
    elif orient_raw == 'LIA': reorient = ' -r 1 2 3 '
    elif orient_raw == 'RSP': reorient = ' -r -1 -2 -3 '
    elif orient_raw == 'LSP': reorient = ' -r 1 -2 -3 '
    elif orient_raw == 'RPS': reorient = ' -r -1 -3 -2 '
    elif orient_raw == 'LPS': reorient = ' -r 1 -3 -2 '

    bckFS_cmd = ' --in_orientation LIA' + reorient

    #  spec file used by Workbench %
    #  inspired from FreeSurfer2CaretConvertAndRegisterNonlinear.sh within the HCP pipeline  / NHP pipeline
    # 2023.02 Simon Clavagnier made modifications to fit python

    # organization
    #    0.0 settings

    ##### specific to macaque ###### to change for each species!!!!!! XXXX creat a roulette to change XXX skip the step if not surface atlas
    if species == "Human": ####XXX where are they?
        FS_mesh = opj(FS_dir, 'standard_mesh_atlases_macaque')
        label_FS_all = ' ' + opj(FS_mesh, 'FreeSurferAllLut.txt') + ' '
        label_FS_Sc = ' ' + opj(FS_mesh, 'FreeSurferSubcorticalLabelTableLut.txt') + ' '
    else:
        FS_mesh = opj(FS_dir, 'standard_mesh_atlases_macaque')
        label_FS_all = ' ' + opj(FS_mesh, 'FreeSurferAllLut.txt') + ' '
        label_FS_Sc = ' ' + opj(FS_mesh, 'FreeSurferSubcorticalLabelTableLut.txt') + ' '

    ################################ animal_folder space image anat ################################

    # Path the animal_folder we are working on
    wb_native_dir = opj(dir_native, '02_Wb')
    dir_native_resol = opj(wb_native_dir, 'surfaces', 'Native_resol')

    #  1.0 start copying
    if ope(opj(wb_native_dir, 'surfaces')) == False:
        os.makedirs(opj(dir_native_resol))

        # transform matrix

    cmd = 'singularity run' + s_bind + fs_sif + 'mri_convert -odt float ' + Ref_file + ' ' + opj(wb_native_dir,'dummy.nii.gz')
    spco(cmd,shell=True)

    cras = str(spco('mri_info --cras ' + opj(wb_native_dir,'dummy.nii.gz'), shell=True))
    os.remove(opj(wb_native_dir,'dummy.nii.gz'))
    mm = cras[2:-3].split(" ")
    c_ras = [[1, 0, 0, float(mm[0])], [0, 1, 0, float(mm[1])], [0, 0, 1, float(mm[2])]]
    c_ras_r = [[-1, 0, 0, float(mm[0])], [0, -1, 0, float(mm[1])], [0, 0, 1, float(mm[2])]]

    with open(opj(wb_native_dir, 'volumes', 'c_ras.txt'), 'wb') as f:
        np.savetxt(f, c_ras, fmt='%.2f')
    shutil.move(opj(wb_native_dir, 'volumes', 'c_ras.txt'), opj(wb_native_dir, 'volumes', 'c_ras.mat'))
    
    with open(opj(wb_native_dir, 'volumes', 'c_ras_r.txt'), 'wb') as f:
        np.savetxt(f, c_ras_r, fmt='%.2f')
    shutil.move(opj(wb_native_dir, 'volumes', 'c_ras_r.txt'), opj(wb_native_dir, 'volumes', 'c_ras_r.mat'))


    #####################################################################################################################################################################################################################
    #                    1.2 Proceed with the surfaces                                                  #####################################################################################################################################################################################################################
    for h in range(2):
        # white surface
        command = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + ' ' + CORTEX[h] + \
                  ' -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -surface-apply-affine ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + ' ' + opj(wb_native_dir, 'volumes', 'c_ras.mat') + ' ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii')
        spco([command], shell=True)

        # pial surface
        command = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.pial') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + ' ' + CORTEX[h] + \
                  ' -surface-type ANATOMICAL -surface-secondary-type PIAL' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -surface-apply-affine ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + ' ' + opj(wb_native_dir, 'volumes', 'c_ras.mat') + ' ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii')
        spco([command], shell=True)

        # midthickness surface + inflated and very inflated
        command = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.mid') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + ' ' + CORTEX[h] + \
                  ' -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -surface-apply-affine ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + ' ' + opj(wb_native_dir, 'volumes', 'c_ras.mat') + ' ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii')
        spco([command], shell=True)

        # from midthickness generation of inflated and very inflated surface in native space
        command = 'singularity run' + s_bind + wb_sif + 'wb_command -file-information ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' | grep "Number of Vertices:" | cut -f2 -d: | tr -d "[:space:]"'
        NativeVerts = spco([command], shell=True)
        NativeInflationScale = int(NativeVerts) * conv_voxmm / 32492  # 32492 (so 32k) for human or MY19 low res or  74k for F99

        command = 'singularity run' + s_bind + wb_sif + 'wb_command -surface-generate-inflated ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.inflated.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.very_inflated.surf.gii') + \
                  ' -iterations-scale ' + str(NativeInflationScale)
        spco([command], shell=True)

        command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.inflated.surf.gii')
        spco([command], shell=True)

        command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.very_inflated.surf.gii')
        spco([command], shell=True)

        # flat map
        if h == 0 and opi(opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.full.flatten.patch.3d')) == True :
            command = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 -p ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.full.flatten.patch.3d') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + ' ' + \
                      CORTEX[h] + \
                      ' -surface-type FLAT -surface-secondary-type GRAY_WHITE' + \
                      ';singularity run' + s_bind + wb_sif + 'wb_command -surface-apply-affine ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + ' ' + opj(
                    wb_native_dir, 'volumes', 'c_ras.mat') + \
                      ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + \
                      ';singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                                                                 animal_folder + '_native_LR.wb.spec') + ' ' + \
                      CORTEX[h] + \
                      ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii')
            nl = spgo(command)

        elif h == 1 and opi(opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.full.flatten.patch.3d')) == True :
            command = 'singularity run' + s_bind + wb_sif + 'mris_convert_v6 -p ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.full.flatten.patch.3d') + \
                      ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + \
                      ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + ' ' + \
                      CORTEX[h] + \
                      ' -surface-type FLAT -surface-secondary-type GRAY_WHITE' + \
                      ';singularity run' + s_bind + wb_sif + 'wb_command -surface-apply-affine ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + ' ' + opj(wb_native_dir, 'volumes', 'c_ras_r.mat') + \
                      ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + \
                      ';singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,animal_folder + '_native_LR.wb.spec') + ' ' + \
                      CORTEX[h] + \
                      ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii')
            nl = spgo(command)
                
        ##### apply transfo aux fichier informations (that color the vertex) XXX
        # sulc thickness and curvature for native brain
        command = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf',Hmin[h] + 'h.sulc') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + ' ' + CORTEX[h] + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-math "var * -1" ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ' -var var ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Sulc' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true'
        spco([command], shell=True)

        command = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.curv') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' ' + CORTEX[h] + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-math "var * -1" ' + ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' -var var ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Curvature' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true'
        spco([command], shell=True)

        command = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.thickness') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' ' + CORTEX[h] + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-math "var * -1" ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' -var var ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Thickness' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-math "abs(thickness)" ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' -var thickness ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-math "thickness > 0" ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' -var thickness ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii')
        spco([command], shell=True)

        command = 'singularity run' + s_bind + wb_sif + 'wb_command -metric-fill-holes ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-remove-islands ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_ROI' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-dilate ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' 10 ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' -nearest' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -metric-dilate ' + opj(dir_native_resol,animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' 10 ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' -nearest'
        spco([command], shell=True)



        for atlas in list_atlases_2:
            cmd = 'singularity run' + s_bind + fs_sif + 'mris_convert_v6 --annot ' + opj(FS_dir, animal_folder, 'label', Hmin[h] + 'h.' + animal_folder + '_' + opb(atlas) + '.annot') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_native_resol, opb(atlas) + '.' + Hmin[h] + '.native.label.gii') + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-structure ' + opj(dir_native_resol, opb(atlas) + '.' + Hmin[h] + '.native.label.gii') + ' ' + CORTEX[h] + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol, opb(atlas) + '.' + Hmin[h] + '.native.label.gii') + ' -map 1 atlas' + \
                  ';singularity run' + s_bind + wb_sif + 'wb_command -gifti-label-add-prefix ' + opj(dir_native_resol, opb(atlas) + '.' + Hmin[h] + '.native.label.gii') + \
                  ' ' + Hmin[h] + '_ ' + opj(dir_native_resol, opb(atlas) + '.' + Hmin[h] + '.native.label.gii')
            spco([cmd], shell=True)

    print('creation of the spec files: done!')


    ##        1.3 Create CIFTI once both hemispheres have been proceed           #################################################################
    # in wb_native_dir,'Native',

    command = 'singularity run' + s_bind + wb_sif + 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' -left-metric ' + opj(dir_native_resol, animal_folder + '.l.sulc.shape.gii') + \
              ' -right-metric ' + opj(dir_native_resol, animal_folder + '.r.sulc.shape.gii') + \
              ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.sulc.dscalar.nii') + ' -map 1 ' + animal_folder + '_Sulc' + \
              ';singularity run' + s_bind + wb_sif + 'wb_command -cifti-palette ' + opj(dir_native_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'singularity run' + s_bind + wb_sif + 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' -left-metric ' + opj(dir_native_resol, animal_folder + '.l.curvature.shape.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-metric ' + opj(dir_native_resol, animal_folder + '.r.curvature.shape.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.curvature.dscalar.nii') + ' -map 1 ' + animal_folder + '_Curvature' + \
              ';singularity run' + s_bind + wb_sif + 'wb_command -cifti-palette ' + opj(dir_native_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'singularity run' + s_bind + wb_sif + 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' -left-metric ' + opj(dir_native_resol, animal_folder + '.l.thickness.shape.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-metric ' + opj(dir_native_resol, animal_folder + '.r.thickness.shape.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.thickness.dscalar.nii') + ' -map 1 ' + animal_folder + '_thickness' + \
              ';wb_command -cifti-palette ' + opj(dir_native_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false'
    spco([command], shell=True)

    for atlas in list_atlases_2:
        cmd = 'singularity run' + s_bind + wb_sif + 'wb_command -cifti-create-label ' + opj(dir_native_resol, opb(atlas) + '.dlabel.nii') + \
              ' -left-label ' + opj(dir_native_resol, opb(atlas) + '.l.native.label.gii') + ' -roi-left ' + opj(dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_native_resol, opb(atlas) + '.r.native.label.gii') + ' -roi-right ' + opj(dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';singularity run' + s_bind + wb_sif + 'wb_command -set-map-names ' + opj(dir_native_resol, opb(atlas) + '.dlabel.nii') + ' -map 1 atlas' + \
              ';singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' +opj(dir_native_resol, opb(atlas) + '.dlabel.nii')
        spco([cmd], shell=True)

    command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.sulc.dscalar.nii')
    spco([command], shell=True)
    command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.curvature.dscalar.nii')
    spco([command], shell=True)
    command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.thickness.dscalar.nii')
    spco([command], shell=True)

    command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + Ref_file
    spco([command], shell=True)

    print('creation of the cifti files: done!')

    #                 2.0 ribbon              #########################################################"

    #     for native space:  + supra / infra masks

    Surf_in = 'white'
    Surf_out = ['pial', 'midthickness']
    vol_out = ['ribbon', 'infra']

    for h in range(0, 2):
        command = 'singularity run' + s_bind + wb_sif + 'wb_command -create-signed-distance-volume ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
            h] + '.' + Surf_in + '.surf.gii') + \
                  ' ' + Ref_file + ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '.nii.gz')
        spco([command], shell=True)
        command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                    animal_folder + '.' + Hmin[h] + '.' + Surf_in + '.nii.gz') + ' -thr 0 -bin -mul 255 ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz')
        spco([command], shell=True)
        command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                    animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz') + ' -bin ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz')
        spco([command], shell=True)

        for s in range(0, 2):
            command = 'singularity run' + s_bind + wb_sif + 'wb_command -create-signed-distance-volume ' + opj(dir_native_resol,
                                                                         animal_folder + '.' + Hmin[h] + '.' + Surf_out[
                                                                             s] + '.surf.gii') + \
                      ' ' + Ref_file + ' ' + opj(dir_native_resol,
                                                 animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[
                s] + '.nii.gz') + ' -uthr 0 -abs -bin -mul 255 ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz') + ' -bin ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz') + ' -mas ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz') + ' -mul 255 ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz') + ' -bin -mul ' + \
                      str(RibbonValue[h][1]) + ' ' + opj(dir_native_resol,
                                                         animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
                h] + '.' + Surf_in + '.nii.gz') + ' -uthr 0 -abs -bin -mul 255 ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz') + ' -bin ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz') + ' -mul ' + \
                      str(RibbonValue[h][0]) + ' ' + opj(dir_native_resol,
                                                         animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_mask.nii.gz')
            spco([command], shell=True)
            command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz') + ' -add ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_mask.nii.gz') + \
                      ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz')
            spco([command], shell=True)

            os.remove(opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz'))
            os.remove(opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_mask.nii.gz'))
            os.remove(opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '.nii.gz'))
            os.remove(opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz'))

        os.remove(opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '.nii.gz'))
        os.remove(opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz'))

    command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol, animal_folder + '.l.ribbon.nii.gz') + ' -add ' + opj(dir_native_resol,
                                                                                                animal_folder + '.r.ribbon.nii.gz') + \
              ' ' + opj(wb_native_dir, 'volumes', 'masks', 'ribbon.nii.gz')
    spco([command], shell=True)

    os.remove(opj(dir_native_resol, animal_folder + '.l.ribbon.nii.gz'))
    os.remove(opj(dir_native_resol, animal_folder + '.r.ribbon.nii.gz'))

    command = 'singularity run' + s_bind + wb_sif + 'wb_command -volume-label-import ' + opj(wb_native_dir, 'volumes', 'masks',
                                                       'ribbon.nii.gz') + label_FS_all + \
              opj(wb_native_dir, 'volumes', 'masks', 'ribbon.nii.gz') + ' -drop-unused-labels'
    spco([command], shell=True)

    command = 'singularity run' + s_bind + fsl_sif + 'fslmaths ' + opj(dir_native_resol, animal_folder + '.l.infra.nii.gz') + ' -add ' + opj(dir_native_resol,
                                                                                               animal_folder + '.r.infra.nii.gz') + \
              ' ' + opj(dir_native_resol, animal_folder + '.infra_1.nii.gz')
    spco([command], shell=True)

    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -prefix ' + opj(wb_native_dir, 'volumes', 'masks', animal_folder + '.infra.nii.gz') + ' -a ' + \
              opj(dir_native_resol, animal_folder + '.infra_1.nii.gz') + ' -expr "amongst(a,3,42)"'
    spco([command], shell=True)

    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -prefix ' + opj(dir_native_resol, 'ribbon_cort.nii.gz') + ' -a ' + opj(wb_native_dir, 'volumes',
                                                                                             'masks',
                                                                                             'ribbon.nii.gz') + ' -expr "amongst(a,3,42)"'
    spco([command], shell=True)
    command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + opj(dir_native_resol, 'ribbon_cort.nii.gz') + ' -b ' + opj(wb_native_dir, 'volumes',
                                                                                        'masks',
                                                                                        animal_folder + '.infra.nii.gz') + \
              ' -expr "a-b" -prefix ' + opj(wb_native_dir, 'volumes', 'masks', animal_folder + '.supra.nii.gz')
    spco([command], shell=True)

    os.remove(opj(dir_native_resol, animal_folder + '.l.infra.nii.gz'))
    os.remove(opj(dir_native_resol, animal_folder + '.r.infra.nii.gz'))
    os.remove(opj(dir_native_resol, animal_folder + '.infra_1.nii.gz'))
    os.remove(opj(dir_native_resol, 'ribbon_cort.nii.gz'))

    command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(wb_native_dir,
                                                                                                       'volumes',
                                                                                                       'masks',
                                                                                                       'ribbon.nii.gz')
    spco([command], shell=True)
    command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(wb_native_dir,
                                                                                                       'volumes',
                                                                                                       'masks',
                                                                                                       animal_folder + '.infra.nii.gz')
    spco([command], shell=True)
    command = 'singularity run' + s_bind + wb_sif + 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(wb_native_dir,
                                                                                                       'volumes',
                                                                                                       'masks',
                                                                                                       animal_folder + '.supra.nii.gz')
    spco([command], shell=True)
