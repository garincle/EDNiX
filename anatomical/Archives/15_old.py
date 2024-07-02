# import
import os
import subprocess
import glob
import shutil
import numpy as np

# import numpy as np
# import nibabel as nib

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output

# Constant
conv_voxmm = 0.75  # 0.75 if FS to wb (1mm to 0.7mm volume)
SurfaceSmoothingFWHM = 2

Hmin = ['l', 'r']
Hcap = ['L', 'R']
CORTEX = ['CORTEX_LEFT', 'CORTEX_RIGHT']
RibbonValue = [[2, 3], [41, 42]]


def WB_prep(FS_dir, path_ATLAS, data_path, animal_folder, Ref_file):
    cmd = 'mri_info --orientation ' + Ref_file
    orient_raw = spgo(cmd)

    # purpouse for FS : LIA
    if orient_raw == 'RAS':
        reorient = ' -r -1 3 -2 '
    elif orient_raw == 'LAS':
        reorient = ' -r 1 3 -2 '
    elif orient_raw == 'RSA':
        reorient = ' -r -1 -2 3'
    elif orient_raw == 'LSA':
        reorient = ' -r 1 -2 3 '
    elif orient_raw == 'RIA':
        reorient = ' -r -1 2 3 '
    elif orient_raw == 'LIA':
        reorient = ' -r 1 2 3 '
    elif orient_raw == 'RSP':
        reorient = ' -r -1 -2 -3 '
    elif orient_raw == 'LSP':
        reorient = ' -r 1 -2 -3 '
    elif orient_raw == 'RPS':
        reorient = ' -r -1 -3 -2 '
    elif orient_raw == 'LPS':
        reorient = ' -r 1 -3 -2 '

    fwdFS_cmd = ' --in_orientation ' + orient_raw + reorient

    #  spec file used by Workbench %
    #  inspired from FreeSurfer2CaretConvertAndRegisterNonlinear.sh within the HCP pipeline  / NHP pipeline
    # 2023.02 Simon Clavagnier made modifications to fit python

    # organization
    #    0.0 settings

    ##### specific to macaque ###### to change for each species!!!!!! XXXX creat a roulette to change XXX skip the step if not surface atlas

    FS_mesh = opj(FS_dir, 'standard_mesh_atlases_macaque')
    label_FS_all = ' ' + opj(FS_mesh, 'FreeSurferAllLut.txt') + ' '
    label_FS_Sc = ' ' + opj(FS_mesh, 'FreeSurferSubcorticalLabelTableLut.txt') + ' '
    REF_yerkes = opj(path_ATLAS, 'Yerkes', 'macaque19')
    Yerkes_label = opj(REF_yerkes, 'MacaqueYerkes19_v1.2_976nz', 'parcellations',
                       'Yerkes19_Parcellations_v2.32k_fs_LR.dlabel.nii')
    # gray             = opj(REF_yerkes,'atlas','Gray.nii.gz')

    wb_native_dir = opj(data_path, 'native', '02_Wb')
    dir_native_resol = opj(wb_native_dir, 'surfaces', 'Native_resol')
    dir_native_32 = opj(wb_native_dir, 'surfaces', 'fsaverage_LR_32')
    Ref_file = opj(wb_native_dir, 'volumes', animal_folder + '_brain.nii.gz')

    wb_atlas_dir = opj(data_path, 'atlas_BALSA')
    dir_atlas_resol = opj(wb_atlas_dir, 'surfaces', 'native')
    dir_atlas_32 = opj(wb_atlas_dir, 'surfaces', 'fsaverage_LR_32')
    dir_atlas_64 = opj(wb_atlas_dir, 'surfaces', 'fsaverage_164')

    Ref_norm = opj(wb_atlas_dir, 'volumes', animal_folder + '_to_Y19_SyN.nii.gz')
    Ref_warp = opj(data_path, 'matrices', animal_folder + '_to_Y19_SyN_1Warp.nii.gz')
    Ref_warp2 = opj(data_path, 'matrices', animal_folder + '_to_Y19_SyN_0GenericAffine.mat')

    ##### convert ANTs tranfo into WB ones
    command = 'ConvertTransformFile 3 ' + \
              ' ' + opj(dir_transfo, monkey + '_to_fascicularis_SyN_0GenericAffine.mat') + \
              ' ' + opj(dir_transfo, monkey + '_to_fascicularis_for_surface.mat') + \
              ' --homogeneousMatrix --RAS'
    spco([command], shell=True)

    command = 'wb_command -convert-affine -from-world ' + \
              opj(dir_transfo, monkey + '_to_fascicularis_for_surface.mat') + ' -inverse' + \
              ' -to-world ' + opj(dir_transfo, 'affine2NHP_fasci.nii.gz')
    spco([command], shell=True)

    command = 'wb_command -convert-warpfield -from-itk ' + \
              opj(dir_transfo, monkey + '_to_fascicularis_SyN_1InverseWarp.nii.gz') + \
              ' -to-world ' + opj(dir_transfo, 'standard2NHP_fasci.nii.gz')
    spco([command], shell=True)

    ###transfo surface specific to WB (ANTs to suraces)
    ref_transfo1 = opj(data_path, 'matrices', 'standard2NHP.nii.gz')
    ref_transfo2 = opj(data_path, 'matrices', 'affine2NHP.nii.gz')

    #  1.0 start copying
    if ope(opj(wb_native_dir, 'surfaces')) == False:
        os.makedirs(opj(dir_native_resol))
        os.makedirs(opj(dir_native_32))

    if ope(opj(wb_atlas_dir, 'surfaces')) == False:
        os.makedirs(dir_atlas_resol)
        os.makedirs(dir_atlas_64)
        os.makedirs(dir_atlas_32)

    # transform matrix creat a .matin wich WB can understandd the transofarmation (creat empty affin matrix)
    command = 'mri_info ' + opj(FS_dir, animal_folder, 'mri', 'orig.mgz') + ' | grep "c_r" | cut -d "=" -f 5 | sed s/" "/""/g'
    mX = spco([command], shell=True)
    command = 'mri_info ' + opj(FS_dir, animal_folder, 'mri', 'orig.mgz') + ' | grep "c_a" | cut -d "=" -f 5 | sed s/" "/""/g'
    mY = spco([command], shell=True)
    command = 'mri_info ' + opj(FS_dir, animal_folder, 'mri', 'orig.mgz') + ' | grep "c_s" | cut -d "=" -f 5 | sed s/" "/""/g'
    mZ = spco([command], shell=True)

    ###### adapt funtion of the species !!! XXX
    c_ras = [[1, 0, 0, float(mX)], [0, 1, 0, float(mY)], [0, 0, 1, float(mZ)]]
    print(c_ras)

    ## save it
    with open(opj(wb_native_dir, 'volumes', 'c_ras.txt'), 'wb') as f:
        np.savetxt(f, c_ras, fmt='%.2f')
    shutil.move(opj(wb_native_dir, 'volumes', 'c_ras.txt'), opj(wb_native_dir, 'volumes', 'c_ras.mat'))

    # Copy References Files for atlas file donc if atlas file extists!XXX
    for h in range(0, 2):
        #  High res
        shutil.copy(opj(FS_mesh, 'fs_' + Hcap[h], 'fsaverage.' + Hcap[h] + '.sphere.164k_fs_' + Hcap[h] + '.surf.gii'), \
                    opj(dir_atlas_64, animal_folder + '.sphere.164k_fs_' + Hmin[h] + '.surf.gii'))
        shutil.copy(opj(FS_mesh, 'fs_' + Hcap[h], 'fs_' + Hcap[h] + '-to-fs_LR_fsaverage.' + Hcap[h] + '_LR.spherical_std.164k_fs_' + Hcap[h] + '.surf.gii'), \
                    opj(dir_atlas_64, animal_folder + '.def_sphere.164k_fs_' + Hmin[h] + '.surf.gii'))
        shutil.copy(opj(FS_mesh, 'fsaverage.' + Hcap[h] + '_LR.spherical_std.164k_fs_LR.surf.gii'), \
                    opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii'))

        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' ' + \
                  CORTEX[h] + ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        shutil.copy(opj(FS_mesh, Hcap[h] + '.atlasroi.164k_fs_LR.shape.gii'),
                    opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.atlasroi.164k_fs_LR.shape.gii'))
        shutil.copy(opj(FS_mesh, Hcap[h] + '.refsulc.164k_fs_LR.shape.gii'),
                    opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.refsulc.164k_fs_LR.shape.gii'))

        # low res
        shutil.copy(opj(FS_mesh, Hcap[h] + '.sphere.32k_fs_LR.surf.gii'),
                    opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii'))
        shutil.copy(opj(FS_mesh, Hcap[h] + '.atlasroi.32k_fs_LR.shape.gii'),
                    opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.atlasroi.32k_fs_LR.shape.gii'))
        shutil.copy(opj(FS_mesh, 'Atlas_ROIs.0.5.nii.gz'),
                    opj(wb_atlas_dir, 'volumes', 'labels', 'Atlas_ROIs.0.5.nii.gz'))

    # 1.1a get the volumes
    # 1.1b get subcortical structures... ## tranform asegmodified on surface but not reaaly usefull anymore!! XXX if remove take the orignal nii
    command = 'mri_convert ' + fwdFS_cmd + ' ' + opj(FS_dir, animal_folder, 'mri', 'aseg.mgz') + ' ' + opj(wb_native_dir, 'volumes', 'labels', 'wmparc.nii.gz')
    spco([command], shell=True)

    command = 'mri_convert ' + fwdFS_cmd + ' ' +  opj(FS_dir, animal_folder, 'mri', 'segmentation.nii.gz') + \
              ' ' + opj(wb_native_dir, 'volumes', 'labels', 'segmentation.nii.gz')
    spco([command], shell=True)

    ##### put the modified indiv in the common space (so inverse of template to indiv)
    command = 'antsApplyTransforms -d 3 -i ' + opj(wb_native_dir, 'volumes', 'labels', 'wmparc.nii.gz') + \
              ' -r ' + Ref_norm + \
              ' -o ' + opj(wb_atlas_dir, 'volumes', 'labels', 'Atlas_wmparc.nii.gz') + \
              ' -t ' + Ref_warp2 + \
              ' -t ' + Ref_warp + \
              ' -n NearestNeighbor'
    spco([command], shell=True)


    ### import label from FS to WB (only for ASEG)
    command = 'wb_command -volume-label-import ' + opj(wb_native_dir, 'volumes', 'labels', 'wmparc.nii.gz') + label_FS_all + opj(wb_native_dir, 'volumes', 'labels', 'wmparc.nii.gz') + ' -drop-unused-labels'
    spco([command], shell=True)

    command = 'wb_command -volume-label-import ' + opj(wb_atlas_dir, 'volumes', 'labels', 'Atlas_wmparc.nii.gz') + label_FS_all + opj(wb_atlas_dir, 'volumes', 'labels', 'Atlas_wmparc.nii.gz') + ' -drop-unused-labels'
    spco([command], shell=True)

    command = 'wb_command -volume-label-import ' + opj(wb_native_dir, 'volumes', 'labels', 'wmparc.nii.gz') + label_FS_Sc + opj(wb_native_dir, 'volumes', 'labels', 'parc_sc.nii.gz') + ' -discard-others'
    spco([command], shell=True)

    command = 'wb_command -volume-label-import ' + opj(wb_atlas_dir, 'volumes', 'labels', 'Atlas_wmparc.nii.gz') + label_FS_Sc + opj(wb_atlas_dir, 'volumes', 'labels', 'Atlas_parc_sc.nii.gz') + ' -discard-others'
    spco([command], shell=True)

    #####################################################################################################################################################################################################################
    #                    1.2 Proceed with the surfaces                                                  #####################################################################################################################################################################################################################
    for h in range(0, 2):
        # white surface
        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + ' ' + CORTEX[h] + \
                  ' -surface-type ANATOMICAL -surface-secondary-type GRAY_WHITE' + \
                  ';wb_command -surface-apply-affine ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + ' ' + opj(wb_native_dir, 'volumes', 'c_ras.mat') + ' ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii')
        spco([command], shell=True)

        # normalisation surfacique to Yerkes
        command = 'wb_command -surface-apply-warpfield ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + \
                  ' ' + ref_transfo1 + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white.transient.surf.gii') + \
                  ';wb_command -surface-apply-affine ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white.transient.surf.gii') + \
                  ' ' + ref_transfo2 + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white.native.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[h] + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white.native.surf.gii')
        spco([command], shell=True)

        # pial surface
        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.pial') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + ' ' + CORTEX[h] + \
                  ' -surface-type ANATOMICAL -surface-secondary-type PIAL' + \
                  ';wb_command -surface-apply-affine ' + opj(dir_native_resol,
                                                             animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + ' ' + opj(
            wb_native_dir, 'volumes', 'c_ras.mat') + ' ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii')
        spco([command], shell=True)
        # normalisation to Yerkes
        command = 'wb_command -surface-apply-warpfield ' + opj(dir_native_resol,
                                                               animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + \
                  ' ' + ref_transfo1 + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial.transient.surf.gii') + \
                  ';wb_command -surface-apply-affine ' + opj(dir_atlas_resol,
                                                             animal_folder + '.' + Hmin[h] + '.pial.transient.surf.gii') + \
                  ' ' + ref_transfo2 + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial.native.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
                      h] + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial.native.surf.gii')
        spco([command], shell=True)

        # midthickness surface + inflated and very inflated

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf',
                                                                                            Hmin[h] + 'h.mid') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + ' ' + CORTEX[
                      h] + \
                  ' -surface-type ANATOMICAL -surface-secondary-type MIDTHICKNESS' + \
                  ';wb_command -surface-apply-affine ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
            h] + '.midthickness.surf.gii') + ' ' + opj(wb_native_dir, 'volumes', 'c_ras.mat') + ' ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii')
        spco([command], shell=True)
        # normalisation to Yerkes
        command = 'wb_command -surface-apply-warpfield ' + opj(dir_native_resol,
                                                               animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + ref_transfo1 + ' ' + opj(dir_atlas_resol,
                                                 animal_folder + '.' + Hmin[h] + '.midthickness.transient.surf.gii') + \
                  ';wb_command -surface-apply-affine ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.midthickness.transient.surf.gii') + \
                  ' ' + ref_transfo2 + ' ' + opj(dir_atlas_resol,
                                                 animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
                      h] + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii')
        spco([command], shell=True)







        # from midthickness generation of inflated and very inflated surface in native space
        command = 'wb_command -file-information ' + opj(dir_native_resol,
                                                        animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' | grep "Number of Vertices:" | cut -f2 -d: | tr -d "[:space:]"'
        NativeVerts = spco([command], shell=True)
        NativeInflationScale = int(
            NativeVerts) * conv_voxmm / 32492  # 32492 (so 32k) for human or MY19 low res or  74k for F99

        command = 'wb_command -surface-generate-inflated ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.inflated.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.very_inflated.surf.gii') + \
                  ' -iterations-scale ' + str(NativeInflationScale)
        spco([command], shell=True)

        command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
            h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.inflated.surf.gii')
        spco([command], shell=True)
        command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
            h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.very_inflated.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -file-information ' + opj(dir_atlas_resol,
                                                        animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ' | grep "Number of Vertices:" | cut -f2 -d: | tr -d "[:space:]"'
        NativeVerts = spco([command], shell=True)
        NativeInflationScale = int(
            NativeVerts) * conv_voxmm / 32492  # 32492 (so 32k) for human or MY19 low res or  74k for F99

        command = 'wb_command -surface-generate-inflated ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.inflated.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.very_inflated.surf.gii') + \
                  ' -iterations-scale ' + str(NativeInflationScale)
        spco([command], shell=True)

        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
            h] + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.inflated.surf.gii')
        spco([command], shell=True)
        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
            h] + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.very_inflated.surf.gii')
        spco([command], shell=True)

        # sphere and registration
        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf',
                                                                                            Hmin[h] + 'h.sphere') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.native.surf.gii') + \
                  ';wb_command -set-structure ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.sphere.native.surf.gii') + ' ' + \
                  CORTEX[h] + \
                  ' -surface-type SPHERICAL -surface-secondary-type GRAY_WHITE' + \
                  ';/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 ' + opj(FS_dir, animal_folder, 'surf',
                                                                                             Hmin[h] + 'h.sphere.reg') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.native.surf.gii') + \
                  ';wb_command -set-structure ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.sphere.reg.native.surf.gii') + ' ' + \
                  CORTEX[h] + \
                  ' -surface-type SPHERICAL -surface-secondary-type GRAY_WHITE'
        spco([command], shell=True)

        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
            h] + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.native.surf.gii')
        spco([command], shell=True)
        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + CORTEX[
            h] + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.native.surf.gii')
        spco([command], shell=True)

        ##### if flat map done, otherwise don't use it!!!
        # FS flat map (only for Native spec file)
        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 -p ' + opj(FS_dir, animal_folder, 'surf',
                                                                                               Hmin[
                                                                                                   h] + 'h.full.flatten.patch.3d') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + ' ' + CORTEX[h] + \
                  ' -surface-type FLAT -surface-secondary-type GRAY_WHITE' + \
                  ';wb_command -surface-apply-affine ' + opj(dir_native_resol,
                                                             animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + ' ' + opj(
            wb_native_dir, 'volumes', 'c_ras.mat') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.flat.surf.gii')
        spco([command], shell=True)

        ##### apply transfo aux fichier informations ( that color the vertex)
        # sulc thickness and curvature for native brain
        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf',
                                                                                               Hmin[h] + 'h.sulc') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + ' ' + CORTEX[h] + \
                  ';wb_command -metric-math "var * -1" ' + opj(dir_native_resol,
                                                               animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ' -var var ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
            h] + '.sulc.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Sulc' + \
                  ';wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true'
        spco([command], shell=True)

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.curv') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' ' + CORTEX[h] + \
                  ';wb_command -metric-math "var * -1" ' + ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' -var var ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Curvature' + \
                  ';wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true'
        spco([command], shell=True)

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.thickness') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' ' + CORTEX[h] + \
                  ';wb_command -metric-math "var * -1" ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' -var var ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Thickness' + \
                  ';wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true' + \
                  ';wb_command -metric-math "abs(thickness)" ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' -var thickness ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';wb_command -metric-palette ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false' + \
                  ';wb_command -metric-math "thickness > 0" ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' -var thickness ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-fill-holes ' + opj(dir_native_resol,
                                                         animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';wb_command -metric-remove-islands ' + opj(dir_native_resol,
                                                              animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_ROI' + \
                  ';wb_command -metric-dilate ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' 10 ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' -nearest' + \
                  ';wb_command -metric-dilate ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' 10 ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' -nearest'
        spco([command], shell=True)



        # parcellations (.annot) transfo atlas surfacic from FS to WB
        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 --annot ' + opj(FS_dir, animal_folder, 'label', Hmin[h] + 'h.' + animal_folder + '_D99.annot') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + ' ' + CORTEX[h] + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_D99' + \
                  ';wb_command -gifti-label-add-prefix ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + \
                  ' ' + Hmin[h] + '_ ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.D99.native.label.gii')
        spco([command], shell=True)

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 --annot ' + opj(FS_dir, animal_folder,
                                                                                                    'label', Hmin[
                                                                                                        h] + 'h.' + animal_folder + '_CIVM.annot') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_native_resol,
                                                                                     animal_folder + '.' + Hmin[
                                                                                         h] + '.CIVM.native.label.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.CIVM.native.label.gii') + ' ' + CORTEX[
                      h] + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
            h] + '.CIVM.native.label.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_CIVM' + \
                  ';wb_command -gifti-label-add-prefix ' + opj(dir_native_resol,
                                                               animal_folder + '.' + Hmin[h] + '.CIVM.native.label.gii') + \
                  ' ' + Hmin[h] + '_ ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.CIVM.native.label.gii')
        spco([command], shell=True)

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 --annot ' + opj(FS_dir, animal_folder,
                                                                                                    'label', Hmin[
                                                                                                        h] + 'h.' + animal_folder + '_Inia19.annot') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_native_resol,
                                                                                     animal_folder + '.' + Hmin[
                                                                                         h] + '.Inia19.native.label.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.Inia19.native.label.gii') + ' ' + \
                  CORTEX[h] + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
            h] + '.Inia19.native.label.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Inia19' + \
                  ';wb_command -gifti-label-add-prefix ' + opj(dir_native_resol,
                                                               animal_folder + '.' + Hmin[h] + '.Inia19.native.label.gii') + \
                  ' ' + Hmin[h] + '_ ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.Inia19.native.label.gii')
        spco([command], shell=True)

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 --annot ' + opj(FS_dir, animal_folder,
                                                                                                    'label', Hmin[
                                                                                                        h] + 'h.' + animal_folder + '_atlas.annot') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_native_resol,
                                                                                     animal_folder + '.' + Hmin[
                                                                                         h] + '.atlas.native.label.gii') + \
                  ';wb_command -set-structure ' + opj(dir_native_resol,
                                                      animal_folder + '.' + Hmin[h] + '.atlas.native.label.gii') + ' ' + \
                  CORTEX[h] + \
                  ';wb_command -set-map-names ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
            h] + '.atlas.native.label.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_atlas' + \
                  ';wb_command -gifti-label-add-prefix ' + opj(dir_native_resol,
                                                               animal_folder + '.' + Hmin[h] + '.atlas.native.label.gii') + \
                  ' ' + Hmin[h] + '_ ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.atlas.native.label.gii')
        spco([command], shell=True)

        #  same thing for Atlas dir

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf',
                                                                                               Hmin[
                                                                                                   h] + 'h.sulc') + ' ' + \
                  opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.sulc.shape.gii') + \
                  ';wb_command -set-structure ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + ' ' + CORTEX[h] + \
                  ';wb_command -metric-math "var * -1" ' + opj(dir_atlas_resol,
                                                               animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ' -var var ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Sulc' + \
                  ';wb_command -metric-palette ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true'
        spco([command], shell=True)

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf',
                                                                                               Hmin[h] + 'h.curv') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_atlas_resol,
                                                                                     animal_folder + '.' + Hmin[
                                                                                         h] + '.curvature.shape.gii') + \
                  ';wb_command -set-structure ' + opj(dir_atlas_resol, animal_folder + '.l.curvature.shape.gii') + ' ' + \
                  CORTEX[h] + \
                  ';wb_command -metric-math "var * -1" ' + opj(dir_atlas_resol,
                                                               animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' -var var ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Curvature' + \
                  ';wb_command -metric-palette ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true'
        spco([command], shell=True)

        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 -c ' + opj(FS_dir, animal_folder, 'surf',
                                                                                               Hmin[
                                                                                                   h] + 'h.thickness') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_atlas_resol,
                                                                                     animal_folder + '.' + Hmin[
                                                                                         h] + '.thickness.shape.gii') + \
                  ';wb_command -set-structure ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' ' + CORTEX[
                      h] + \
                  ';wb_command -metric-math "var * -1" ' + opj(dir_atlas_resol,
                                                               animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' -var var ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.thickness.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_Thickness' + \
                  ';wb_command -metric-palette ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true' + \
                  ';wb_command -metric-math "abs(thickness)" ' + opj(dir_atlas_resol,
                                                                     animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' -var thickness ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ';wb_command -metric-palette ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false' + \
                  ';wb_command -metric-math "thickness > 0" ' + opj(dir_atlas_resol,
                                                                    animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ' -var thickness ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-fill-holes ' + opj(dir_atlas_resol,
                                                         animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + ' ' + opj(dir_atlas_resol,
                                                                                                    animal_folder + '.' + Hmin[
                                                                                                        h] + '.roi.shape.gii') + \
                  ';wb_command -metric-remove-islands ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.midthickness.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + ' ' + opj(dir_atlas_resol,
                                                                                                    animal_folder + '.' + Hmin[
                                                                                                        h] + '.roi.shape.gii') + \
                  ';wb_command -set-map-names ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.roi.shape.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_ROI' + \
                  ';wb_command -metric-dilate ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ' 10 ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' -nearest' + \
                  ';wb_command -metric-dilate ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ' 10 ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' -nearest'
        spco([command], shell=True)


        # parcellations (.annot) loop for each atlas!!! XXX


        command = 'wb_command -volume-label-import' + \
                  ' ' + opj(dir_prepro, 'atlas_tmp.nii.gz') + \
                  ' ' + opj(REF_fasci, 'volumes', 'labels', 'fascicularis_atlas_label.txt') + \
                  ' ' + opj(labels_dir, monkey + '_atlas.nii.gz') + ' -drop-unused-labels'
        spco([command], shell=True)



        command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 --annot ' + opj(FS_dir, animal_folder,
                                                                                                    'label', Hmin[
                                                                                                        h] + 'h.' + animal_folder + '_D99.annot') + \
                  ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[h] + 'h.white') + ' ' + opj(dir_atlas_resol,
                                                                                     animal_folder + '.' + Hmin[
                                                                                         h] + '.D99.native.label.gii') + \
                  ';wb_command -set-structure ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + ' ' + CORTEX[
                      h] + \
                  ';wb_command -set-map-names ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.D99.native.label.gii') + ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_D99' + \
                  ';wb_command -gifti-label-add-prefix ' + opj(dir_atlas_resol,
                                                               animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + \
                  ' ' + Hmin[h] + '_ ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.D99.native.label.gii')
        spco([command], shell=True)





        # to ref reference surface if exists!!!
        command = 'wb_command -surface-sphere-project-unproject ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.sphere.reg.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.sphere.164k_fs_' + Hmin[h] + '.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.def_sphere.164k_fs_' + Hmin[h] + '.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii')
        spco([command], shell=True)

        # Calculate and visualise distorsion.......................................
        command = 'wb_command -surface-vertex-areas ' + opj(dir_atlas_resol,
                                                            animal_folder + '.' + Hmin[h] + '.sphere.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.native.shape.gii') + \
                  ';wb_command -surface-vertex-areas ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.shape.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-math "ln(spherereg / sphere) / ln(2)" ' + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.ArealDistorsion_FS.shape.gii') + \
                  ' -var sphere ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.native.shape.gii') + \
                  ' -var spherereg ' + opj(dir_atlas_resol,
                                           animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.shape.gii')
        spco([command], shell=True)

        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.native.shape.gii'))
        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.shape.gii'))

        command = 'wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                     animal_folder + '.' + Hmin[h] + '.ArealDistorsion_FS.shape.gii') + \
                  ' -map 1 ' + animal_folder + '_' + Hmin[h] + '_ArealDistorsion_FS' + \
                  ';wb_command -metric-palette ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.ArealDistorsion_FS.shape.gii') + \
                  ' MODE_AUTO_SCALE -palette-name ROY_BIG_BL -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE -1 1'
        spco([command], shell=True)

        command = 'wb_command -surface-distortion ' + opj(dir_atlas_resol,
                                                          animal_folder + '.' + Hmin[h] + '.sphere.native.surf.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol,
                            animal_folder + '.' + Hmin[h] + '.EdgeDistorsion_FS.native.shape.gii') + ' -edge-method' + \
                  ';wb_command -surface-distortion ' + opj(dir_atlas_resol,
                                                           animal_folder + '.' + Hmin[h] + '.sphere.native.surf.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol,
                            animal_folder + '.' + Hmin[h] + '.Strain_FS.native.shape.gii') + ' -local-affine-method' + \
                  ';wb_command -metric-merge ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.StrainJ_FS.native.shape.gii') + ' -metric ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.Strain_FS.native.shape.gii') + ' -column 1' + \
                  ';wb_command -metric-merge ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.StrainR_FS.native.shape.gii') + ' -metric ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.Strain_FS.native.shape.gii') + ' -column 2' + \
                  ';wb_command -metric-math "ln(var) / ln(2)" ' + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.StrainJ_FS.native.shape.gii') + \
                  ' -var var ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.StrainJ_FS.native.shape.gii') + \
                  ';wb_command -metric-math "ln(var) / ln(2)" ' + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.StrainR_FS.native.shape.gii') + \
                  ' -var var ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.StrainR_FS.native.shape.gii')
        spco([command], shell=True)

        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.Strain_FS.native.shape.gii'))

        # resampling in high resolution : 164k

        command = 'wb_command -surface-resample ' + opj(dir_atlas_resol,
                                                        animal_folder + '.' + Hmin[h] + '.white.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + ' BARYCENTRIC ' + \
                  opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.white.164k_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.white.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -surface-resample ' + opj(dir_atlas_resol,
                                                        animal_folder + '.' + Hmin[h] + '.pial.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.pial.164k_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.pial.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -surface-resample ' + opj(dir_atlas_resol,
                                                        animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.l.midthickness.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        HighResInflationScale = int(164) * conv_voxmm / 32

        command = 'wb_command -surface-generate-inflated ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[
            h] + '.midthickness.164k_fs_LR.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.inflated.164k_fs_LR.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.very_inflated.164k_fs_LR.surf.gii') + \
                  ' -iterations-scale ' + str(HighResInflationScale) + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.inflated.164k_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.very_inflated.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        shutil.copy(opj(FS_mesh, 'colin.cerebral.' + Hcap[h] + '.flat.164k_fs_LR.surf.gii'),
                    opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.flat.164k.surf.gii'))

        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.flat.164k.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.curvature.164k_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii') + \
                  ' -current-roi ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';wb_command -metric-mask ' + opj(dir_atlas_64,
                                                    animal_folder + '.' + Hmin[h] + '.curvature.164k_fs_LR.shape.gii') + ' ' + \
                  opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.atlasroi.164k_fs_LR.shape.gii') + ' ' + \
                  opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.curvature.164k_fs_LR.shape.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sulc.164k_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.thickness.164k_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii') + \
                  ' -current-roi ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';wb_command -metric-mask ' + opj(dir_atlas_64,
                                                    animal_folder + '.' + Hmin[h] + '.thickness.164k_fs_LR.shape.gii') + ' ' + \
                  opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.atlasroi.164k_fs_LR.shape.gii') + ' ' + \
                  opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.thickness.164k_fs_LR.shape.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.ArealDistorsion_FS.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64,
                                           animal_folder + '.' + Hmin[h] + '.ArealDistorsion_FS.164k_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.EdgeDistorsion_FS.native.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64,
                                           animal_folder + '.' + Hmin[h] + '.EdgeDistorsion_FS.164k_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.StrainJ_FS.native.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.StrainJ_FS.164k_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.StrainR_FS.native.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.StrainR_FS.164k_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii')
        spco([command], shell=True)

        # labels loop
        command = 'wb_command -label-resample ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_64, animal_folder + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_atlas_64,
                                        animal_folder + '.' + Hmin[h] + '.D99.164k_fs_LR.label.gii') + ' -largest'
        spco([command], shell=True)



        # resampling in low resolution : 32k
        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command  -surface-resample  ' + opj(dir_atlas_resol,
                                                          animal_folder + '.' + Hmin[h] + '.white.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.white.32_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.white.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command  -surface-resample  ' + opj(dir_atlas_resol,
                                                          animal_folder + '.' + Hmin[h] + '.pial.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.pial.32_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.pial.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command  -surface-resample  ' + opj(dir_atlas_resol,
                                                          animal_folder + '.' + Hmin[h] + '.midthickness.native.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii')
        spco([command], shell=True)

        LowResInflationScale = int(32) * conv_voxmm / 32

        command = 'wb_command -surface-generate-inflated ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[
            h] + '.midthickness.32_fs_LR.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.inflated.32_fs_LR.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.very_inflated.32_fs_LR.surf.gii') + \
                  ' -iterations-scale ' + str(LowResInflationScale) + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.inflated.32_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.very_inflated.32_fs_LR.surf.gii')
        spco([command], shell=True)

        shutil.copy(opj(FS_mesh, 'colin.cerebral.' + Hcap[h] + '.flat.32k_fs_LR.surf.gii'),
                    opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.flat.32k.surf.gii'))
        command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.flat.32k.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.curvature.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.curvature.32_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii') + \
                  ' -current-roi ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';wb_command -metric-mask ' + opj(dir_atlas_32,
                                                    animal_folder + '.' + Hmin[h] + '.curvature.32_fs_LR.shape.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.atlasroi.32k_fs_LR.shape.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.curvature.32_fs_LR.shape.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.sulc.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sulc.32_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.thickness.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.thickness.32_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii') + \
                  ' -current-roi ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.roi.shape.gii') + \
                  ';wb_command -metric-mask ' + opj(dir_atlas_32,
                                                    animal_folder + '.' + Hmin[h] + '.thickness.32_fs_LR.shape.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.atlasroi.32k_fs_LR.shape.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.thickness.32_fs_LR.shape.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.ArealDistorsion_FS.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_32,
                                           animal_folder + '.' + Hmin[h] + '.ArealDistorsion_FS.32_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.EdgeDistorsion_FS.native.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_32,
                                           animal_folder + '.' + Hmin[h] + '.EdgeDistorsion_FS.32_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.StrainJ_FS.native.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.StrainJ_FS.32_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -metric-resample ' + opj(dir_atlas_resol,
                                                       animal_folder + '.' + Hmin[h] + '.StrainR_FS.native.shape.gii') + ' ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.StrainR_FS.32_fs_LR.shape.gii') + \
                  ' -area-surfs ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii')
        spco([command], shell=True)

        # labels # labels loop
        command = 'wb_command -label-resample ' + opj(dir_atlas_resol,
                                                      animal_folder + '.' + Hmin[h] + '.D99.native.label.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.D99.32_fs_LR.label.gii') + ' -largest'
        spco([command], shell=True)


        # now back to native space but with the low resolution (32k)
        command = 'wb_command -surface-resample  ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.white.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.white.32_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_native_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.white.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -surface-resample  ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.pial.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.pial.32_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_native_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.pial.32_fs_LR.surf.gii')
        spco([command], shell=True)

        command = 'wb_command -surface-resample  ' + opj(dir_native_resol,
                                                         animal_folder + '.' + Hmin[h] + '.midthickness.surf.gii') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') + \
                  ' ' + opj(dir_atlas_32, animal_folder + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') + \
                  ' BARYCENTRIC ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii') + \
                  ';wb_command -add-to-spec-file ' + opj(dir_native_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.midthickness.32_fs_LR.surf.gii')
        spco([command], shell=True)

        LowResInflationScale = int(32) * conv_voxmm / 32
        command = 'wb_command -surface-generate-inflated ' + opj(dir_native_32, animal_folder + '.' + Hmin[
            h] + '.midthickness.32_fs_LR.surf.gii') + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.inflated.32_fs_LR.surf.gii') + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.very_inflated.32_fs_LR.surf.gii') + \
                  ' -iterations-scale ' + str(LowResInflationScale)
        spco([command], shell=True)

        command = 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.inflated.32_fs_LR.surf.gii')
        spco([command], shell=True)
        command = 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.very_inflated.32_fs_LR.surf.gii')
        spco([command], shell=True)
        command = 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' ' + \
                  CORTEX[h] + \
                  ' ' + opj(dir_native_32, animal_folder + '.' + Hmin[h] + '.very_inflated.32_fs_LR.surf.gii')
        spco([command], shell=True)

    print('creation of the spec files: done!')





    ##        1.3 Create CIFTI once both hemispheres have been proceed           #################################################################
    # in wb_native_dir,'Native',

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' -left-metric ' + opj(dir_native_resol, animal_folder + '.l.sulc.shape.gii') + \
              ' -right-metric ' + opj(dir_native_resol, animal_folder + '.r.sulc.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.sulc.dscalar.nii') + ' -map 1 ' + animal_folder + '_Sulc' + \
              ';wb_command -cifti-palette ' + opj(dir_native_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' -left-metric ' + opj(dir_native_resol, animal_folder + '.l.curvature.shape.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-metric ' + opj(dir_native_resol, animal_folder + '.r.curvature.shape.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.curvature.dscalar.nii') + ' -map 1 ' + animal_folder + '_Curvature' + \
              ';wb_command -cifti-palette ' + opj(dir_native_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_native_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' -left-metric ' + opj(dir_native_resol, animal_folder + '.l.thickness.shape.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-metric ' + opj(dir_native_resol, animal_folder + '.r.thickness.shape.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.thickness.dscalar.nii') + ' -map 1 ' + animal_folder + '_thickness' + \
              ';wb_command -cifti-palette ' + opj(dir_native_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_native_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_native_resol, animal_folder + '.D99.dlabel.nii') + \
              ' -left-label ' + opj(dir_native_resol, animal_folder + '.l.D99.native.label.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_native_resol, animal_folder + '.r.D99.native.label.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.D99.dlabel.nii') + ' -map 1 ' + animal_folder + '_D99'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_native_resol, animal_folder + '.CIVM.dlabel.nii') + \
              ' -left-label ' + opj(dir_native_resol, animal_folder + '.l.CIVM.native.label.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_native_resol, animal_folder + '.r.CIVM.native.label.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.CIVM.dlabel.nii') + ' -map 1 ' + animal_folder + '_CIVM'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_native_resol, animal_folder + '.Inia19.dlabel.nii') + \
              ' -left-label ' + opj(dir_native_resol, animal_folder + '.l.Inia19.native.label.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_native_resol, animal_folder + '.r.Inia19.native.label.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.Inia19.dlabel.nii') + ' -map 1 ' + animal_folder + '_Inia19'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_native_resol, animal_folder + '.atlas.dlabel.nii') + \
              ' -left-label ' + opj(dir_native_resol, animal_folder + '.l.atlas.native.label.gii') + ' -roi-left ' + opj(
        dir_native_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_native_resol, animal_folder + '.r.atlas.native.label.gii') + ' -roi-right ' + opj(
        dir_native_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_native_resol,
                                                  animal_folder + '.atlas.dlabel.nii') + ' -map 1 ' + animal_folder + '_atlas'
    spco([command], shell=True)

    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.sulc.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.curvature.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.thickness.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.D99.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.Inia19.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.CIVM.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.atlas.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + Ref_file
    spco([command], shell=True)



    #       dir_atlas_resol                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_resol, animal_folder + '.l.sulc.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_resol, animal_folder + '.r.sulc.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.sulc.dscalar.nii') + ' -map 1 ' + animal_folder + '_Sulc' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_resol, animal_folder + '.sulc.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_resol, animal_folder + '.l.curvature.shape.gii') + ' -roi-left ' + opj(
        dir_atlas_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_resol, animal_folder + '.r.curvature.shape.gii') + ' -roi-right ' + opj(
        dir_atlas_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.curvature.dscalar.nii') + ' -map 1 ' + animal_folder + '_Curvature' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_resol, animal_folder + '.curvature.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_resol, animal_folder + '.l.thickness.shape.gii') + ' -roi-left ' + opj(
        dir_atlas_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_resol, animal_folder + '.r.thickness.shape.gii') + ' -roi-right ' + opj(
        dir_atlas_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.thickness.dscalar.nii') + ' -map 1 ' + animal_folder + '_thickness' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_resol, animal_folder + '.thickness.dscalar.nii') + \
              ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_resol,
                                                             animal_folder + '.ArealDistorsion_FS.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_resol, animal_folder + '.l.ArealDistorsion_FS.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_resol, animal_folder + '.r.ArealDistorsion_FS.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.ArealDistorsion_FS.dscalar.nii') + ' -map 1 ' + animal_folder + '_ArealDistorsion_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_resol, animal_folder + '.ArealDistorsion_FS.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_resol, animal_folder + '.ArealDistorsion_FS.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_resol,
                                                             animal_folder + '.EdgeDistorsion_FS.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_resol, animal_folder + '.l.EdgeDistorsion_FS.native.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_resol, animal_folder + '.r.EdgeDistorsion_FS.native.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.EdgeDistorsion_FS.dscalar.nii') + ' -map 1 ' + animal_folder + '_EdgeDistorsion_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_resol, animal_folder + '.EdgeDistorsion_FS.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_resol, animal_folder + '.EdgeDistorsion_FS.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_resol, animal_folder + '.StrainJ_FS.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_resol, animal_folder + '.l.StrainJ_FS.native.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_resol, animal_folder + '.r.StrainJ_FS.native.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.StrainJ_FS.dscalar.nii') + ' -map 1 ' + animal_folder + '_StrainJ_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_resol, animal_folder + '.StrainJ_FS.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_resol, animal_folder + '.StrainJ_FS.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_resol, animal_folder + '.StrainR_FS.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_resol, animal_folder + '.l.StrainR_FS.native.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_resol, animal_folder + '.r.StrainR_FS.native.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.StrainR_FS.dscalar.nii') + ' -map 1 ' + animal_folder + '_StrainR_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_resol, animal_folder + '.StrainR_FS.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_resol, animal_folder + '.StrainR_FS.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_resol, animal_folder + '.D99.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_resol, animal_folder + '.l.D99.native.label.gii') + ' -roi-left ' + opj(
        dir_atlas_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_resol, animal_folder + '.r.D99.native.label.gii') + ' -roi-right ' + opj(
        dir_atlas_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.D99.dlabel.nii') + ' -map 1 ' + animal_folder + '_D99'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_resol, animal_folder + '.CIVM.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_resol, animal_folder + '.l.CIVM.native.label.gii') + ' -roi-left ' + opj(
        dir_atlas_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_resol, animal_folder + '.r.CIVM.native.label.gii') + ' -roi-right ' + opj(
        dir_atlas_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.CIVM.dlabel.nii') + ' -map 1 ' + animal_folder + '_CIVM'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_resol, animal_folder + '.Inia19.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_resol, animal_folder + '.l.Inia19.native.label.gii') + ' -roi-left ' + opj(
        dir_atlas_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_resol, animal_folder + '.r.Inia19.native.label.gii') + ' -roi-right ' + opj(
        dir_atlas_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.Inia19.dlabel.nii') + ' -map 1 ' + animal_folder + '_Inia19'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_resol, animal_folder + '.atlas.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_resol, animal_folder + '.l.atlas.native.label.gii') + ' -roi-left ' + opj(
        dir_atlas_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_resol, animal_folder + '.r.atlas.native.label.gii') + ' -roi-right ' + opj(
        dir_atlas_resol, animal_folder + '.r.roi.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                  animal_folder + '.atlas.dlabel.nii') + ' -map 1 ' + animal_folder + '_atlas'
    spco([command], shell=True)

    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.sulc.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.curvature.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.thickness.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.ArealDistorsion_FS.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.EdgeDistorsion_FS.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.StrainJ_FS.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.StrainR_FS.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.D99.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.Inia19.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.CIVM.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.atlas.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + Ref_norm
    spco([command], shell=True)

    # wb_atlas_dir,'32k' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_32, animal_folder + '.sulc.32_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_32, animal_folder + '.l.sulc.32_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_32, animal_folder + '.r.sulc.32_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.sulc.32_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_Sulc' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_32, animal_folder + '.sulc.32_fs_LR.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_32, animal_folder + '.sulc.32_fs_LR.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_32, animal_folder + '.curvature.32_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_32, animal_folder + '.l.curvature.32_fs_LR.shape.gii') + ' -roi-left ' + opj(
        dir_atlas_32, animal_folder + '.l.atlasroi.32k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_32, animal_folder + '.r.curvature.32_fs_LR.shape.gii') + ' -roi-right ' + opj(
        dir_atlas_32, animal_folder + '.r.atlasroi.32k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.curvature.32_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_Curvature' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_32, animal_folder + '.curvature.32_fs_LR.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_32, animal_folder + '.curvature.32_fs_LR.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_32, animal_folder + '.thickness.32_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_32, animal_folder + '.l.thickness.32_fs_LR.shape.gii') + ' -roi-left ' + opj(
        dir_atlas_32, animal_folder + '.l.atlasroi.32k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_32, animal_folder + '.r.thickness.32_fs_LR.shape.gii') + ' -roi-right ' + opj(
        dir_atlas_32, animal_folder + '.r.atlasroi.32k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.thickness.32_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_thickness' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_32, animal_folder + '.thickness.32_fs_LR.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_32, animal_folder + '.thickness.32_fs_LR.dscalar.nii') + \
              ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_32,
                                                             animal_folder + '.ArealDistorsion_FS.32_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_32, animal_folder + '.l.ArealDistorsion_FS.32_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_32, animal_folder + '.r.ArealDistorsion_FS.32_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.ArealDistorsion_FS.32_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_ArialDistortion_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_32, animal_folder + '.ArealDistorsion_FS.32_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_32, animal_folder + '.ArealDistorsion_FS.32_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_32,
                                                             animal_folder + '.EdgeDistorsion_FS.32_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_32, animal_folder + '.l.EdgeDistorsion_FS.32_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_32, animal_folder + '.r.EdgeDistorsion_FS.32_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.EdgeDistorsion_FS.32_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_EdgeDistortion_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_32, animal_folder + '.EdgeDistorsion_FS.32_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_32, animal_folder + '.EdgeDistorsion_FS.32_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_32,
                                                             animal_folder + '.StrainJ_FS.32_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_32, animal_folder + '.l.StrainJ_FS.32_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_32, animal_folder + '.r.StrainJ_FS.32_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.StrainJ_FS.32_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_StrainJ' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_32, animal_folder + '.StrainJ_FS.32_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_32, animal_folder + '.StrainJ_FS.32_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_32,
                                                             animal_folder + '.StrainR_FS.32_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_32, animal_folder + '.l.StrainR_FS.32_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_32, animal_folder + '.r.StrainR_FS.32_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.StrainR_FS.32_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_StrainR' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_32, animal_folder + '.StrainR_FS.32_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_32, animal_folder + '.StrainR_FS.32_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_32, animal_folder + '.D99.32_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_32, animal_folder + '.l.D99.32_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_32, animal_folder + '.l.atlasroi.32k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_32, animal_folder + '.r.D99.32_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_32, animal_folder + '.r.atlasroi.32k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.D99.32_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_D99'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_32, animal_folder + '.CIVM.32_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_32, animal_folder + '.l.CIVM.32_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_32, animal_folder + '.l.atlasroi.32k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_32, animal_folder + '.r.CIVM.32_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_32, animal_folder + '.r.atlasroi.32k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.CIVM.32_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_CIVM'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_32, animal_folder + '.Inia19.32_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_32, animal_folder + '.l.Inia19.32_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_32, animal_folder + '.l.atlasroi.32k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_32, animal_folder + '.r.Inia19.32_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_32, animal_folder + '.r.atlasroi.32k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.Inia19.32_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_Inia19'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_32, animal_folder + '.atlas.32_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_32, animal_folder + '.l.atlas.32_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_32, animal_folder + '.l.atlasroi.32k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_32, animal_folder + '.r.atlas.32_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_32, animal_folder + '.r.atlasroi.32k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_32,
                                                  animal_folder + '.atlas.32_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_atlas'
    spco([command], shell=True)

    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.sulc.32_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.curvature.32_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.thickness.32_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + \
              opj(dir_atlas_32, animal_folder + '.ArealDistorsion_FS.32_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + \
              opj(dir_atlas_32, animal_folder + '.EdgeDistorsion_FS.32_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.StrainJ_FS.32_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.StrainR_FS.32_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.D99.32_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.Inia19.32_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.CIVM.32_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + opj(
        dir_atlas_32, animal_folder + '.atlas.32_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32,
                                                    animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + Ref_norm
    spco([command], shell=True)





    #         wb_atlas_dir,'164k'                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_64, animal_folder + '.sulc.164_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_64, animal_folder + '.l.sulc.164k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_64, animal_folder + '.r.sulc.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.sulc.164_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_Sulc' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_64, animal_folder + '.sulc.164_fs_LR.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_64, animal_folder + '.sulc.164_fs_LR.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_64,
                                                             animal_folder + '.curvature.164_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_64, animal_folder + '.l.curvature.164k_fs_LR.shape.gii') + ' -roi-left ' + opj(
        dir_atlas_64, animal_folder + '.l.atlasroi.164k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_64,
                                      animal_folder + '.r.curvature.164k_fs_LR.shape.gii') + ' -roi-right ' + opj(dir_atlas_64,
                                                                                                           animal_folder + '.r.atlasroi.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.curvature.164_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_Curvature' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_64, animal_folder + '.curvature.164_fs_LR.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_64, animal_folder + '.curvature.164_fs_LR.dscalar.nii') + \
              ' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_64,
                                                             animal_folder + '.thickness.164_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_64, animal_folder + '.l.thickness.164k_fs_LR.shape.gii') + ' -roi-left ' + opj(
        dir_atlas_64, animal_folder + '.l.atlasroi.164k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_64,
                                      animal_folder + '.r.thickness.164k_fs_LR.shape.gii') + ' -roi-right ' + opj(dir_atlas_64,
                                                                                                           animal_folder + '.r.atlasroi.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.thickness.164_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_thickness' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_64, animal_folder + '.thickness.164_fs_LR.dscalar.nii') + \
              ' MODE_AUTO_SCALE_PERCENTAGE ' + opj(dir_atlas_64, animal_folder + '.thickness.164_fs_LR.dscalar.nii') + \
              ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_64,
                                                             animal_folder + '.ArealDistorsion_FS.164_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_64, animal_folder + '.l.ArealDistorsion_FS.164k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_64, animal_folder + '.r.ArealDistorsion_FS.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.ArealDistorsion_FS.164_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_ArialDistortion_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_64, animal_folder + '.ArealDistorsion_FS.164_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_64, animal_folder + '.ArealDistorsion_FS.164_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_64,
                                                             animal_folder + '.EdgeDistorsion_FS.164_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_64, animal_folder + '.l.EdgeDistorsion_FS.164k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_64, animal_folder + '.r.EdgeDistorsion_FS.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.EdgeDistorsion_FS.164_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_EdgeDistortion_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_64, animal_folder + '.EdgeDistorsion_FS.164_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_64, animal_folder + '.EdgeDistorsion_FS.164_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_64,
                                                             animal_folder + '.StrainJ_FS.164_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_64, animal_folder + '.l.StrainJ_FS.164k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_64, animal_folder + '.r.StrainJ_FS.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.StrainJ_FS.164_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_StrainJ_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_64, animal_folder + '.StrainJ_FS.164_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_64, animal_folder + '.StrainJ_FS.164_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-dense-scalar ' + opj(dir_atlas_64,
                                                             animal_folder + '.StrainR_FS.164_fs_LR.dscalar.nii') + \
              ' -left-metric ' + opj(dir_atlas_64, animal_folder + '.l.StrainR_FS.164k_fs_LR.shape.gii') + \
              ' -right-metric ' + opj(dir_atlas_64, animal_folder + '.r.StrainR_FS.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.StrainR_FS.164_fs_LR.dscalar.nii') + ' -map 1 ' + animal_folder + '_StrainR_FS' + \
              ';wb_command -cifti-palette ' + opj(dir_atlas_64, animal_folder + '.StrainR_FS.164_fs_LR.dscalar.nii') + \
              ' MODE_USER_SCALE ' + opj(dir_atlas_64, animal_folder + '.StrainR_FS.164_fs_LR.dscalar.nii') + \
              ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_64, animal_folder + '.D99.164_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_64, animal_folder + '.l.D99.164k_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_64, animal_folder + '.l.atlasroi.164k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_64, animal_folder + '.r.D99.164k_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_64, animal_folder + '.r.atlasroi.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.D99.164_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_D99'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_64, animal_folder + '.CIVM.164_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_64, animal_folder + '.l.CIVM.164k_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_64, animal_folder + '.l.atlasroi.164k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_64, animal_folder + '.r.CIVM.164k_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_64, animal_folder + '.r.atlasroi.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.CIVM.164_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_CIVM'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_64, animal_folder + '.Inia19.164_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_64, animal_folder + '.l.Inia19.164k_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_64, animal_folder + '.l.atlasroi.164k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_64, animal_folder + '.r.Inia19.164k_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_64, animal_folder + '.r.atlasroi.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.Inia19.164_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_Inia19'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_64, animal_folder + '.atlas.164_fs_LR.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_64, animal_folder + '.l.atlas.164k_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_64, animal_folder + '.l.atlasroi.164k_fs_LR.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_64, animal_folder + '.r.atlas.164k_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_64, animal_folder + '.r.atlasroi.164k_fs_LR.shape.gii') + \
              ';wb_command -set-map-names ' + opj(dir_atlas_64,
                                                  animal_folder + '.atlas.164_fs_LR.dlabel.nii') + ' -map 1 ' + animal_folder + '_atlas'
    spco([command], shell=True)

    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.sulc.164_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.curvature.164_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.thickness.164_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + \
              opj(dir_atlas_64, animal_folder + '.ArealDistorsion_FS.164_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + \
              opj(dir_atlas_64, animal_folder + '.EdgeDistorsion_FS.164_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.StrainJ_FS.164_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.StrainR_FS.164_fs_LR.dscalar.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.D99.164_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.Inia19.164_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.CIVM.164_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + opj(
        dir_atlas_64, animal_folder + '.atlas.164_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_64,
                                                    animal_folder + '_fsaverage_LR_164.wb.spec') + ' INVALID ' + Ref_norm
    spco([command], shell=True)

    print('creation of the cifti files: done!')






    #                 2.0 ribbon              #########################################################"

    #     for native space:  + supra / infra masks

    Surf_in = 'white'
    Surf_out = ['pial', 'midthickness']
    vol_out = ['ribbon', 'infra']

    for h in range(0, 2):
        command = 'wb_command -create-signed-distance-volume ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
            h] + '.' + Surf_in + '.surf.gii') + \
                  ' ' + Ref_file + ' ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_native_resol,
                                    animal_folder + '.' + Hmin[h] + '.' + Surf_in + '.nii.gz') + ' -thr 0 -bin -mul 255 ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_native_resol,
                                    animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz') + ' -bin ' + \
                  opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz')
        spco([command], shell=True)

        for s in range(0, 2):
            command = 'wb_command -create-signed-distance-volume ' + opj(dir_native_resol,
                                                                         animal_folder + '.' + Hmin[h] + '.' + Surf_out[
                                                                             s] + '.surf.gii') + \
                      ' ' + Ref_file + ' ' + opj(dir_native_resol,
                                                 animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[
                s] + '.nii.gz') + ' -uthr 0 -abs -bin -mul 255 ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz') + ' -bin ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_out[s] + '_uthr0.nii.gz') + ' -mas ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_thr0.nii.gz') + ' -mul 255 ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz') + ' -bin -mul ' + \
                      str(RibbonValue[h][1]) + ' ' + opj(dir_native_resol,
                                                         animal_folder + '.' + Hmin[h] + '.' + vol_out[s] + '.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol, animal_folder + '.' + Hmin[
                h] + '.' + Surf_in + '.nii.gz') + ' -uthr 0 -abs -bin -mul 255 ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz') + ' -bin ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol,
                                        animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_uthr0.nii.gz') + ' -mul ' + \
                      str(RibbonValue[h][0]) + ' ' + opj(dir_native_resol,
                                                         animal_folder + '.' + Hmin[h] + '.' + Surf_in + '_mask.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol,
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

    command = 'fslmaths ' + opj(dir_native_resol, animal_folder + '.l.ribbon.nii.gz') + ' -add ' + opj(dir_native_resol,
                                                                                                animal_folder + '.r.ribbon.nii.gz') + \
              ' ' + opj(wb_native_dir, 'volumes', 'masks', 'ribbon.nii.gz')
    spco([command], shell=True)

    os.remove(opj(dir_native_resol, animal_folder + '.l.ribbon.nii.gz'))
    os.remove(opj(dir_native_resol, animal_folder + '.r.ribbon.nii.gz'))

    command = 'wb_command -volume-label-import ' + opj(wb_native_dir, 'volumes', 'masks',
                                                       'ribbon.nii.gz') + label_FS_all + \
              opj(wb_native_dir, 'volumes', 'masks', 'ribbon.nii.gz') + ' -drop-unused-labels'
    spco([command], shell=True)

    command = 'fslmaths ' + opj(dir_native_resol, animal_folder + '.l.infra.nii.gz') + ' -add ' + opj(dir_native_resol,
                                                                                               animal_folder + '.r.infra.nii.gz') + \
              ' ' + opj(dir_native_resol, animal_folder + '.infra_1.nii.gz')
    spco([command], shell=True)

    command = '3dcalc -prefix ' + opj(wb_native_dir, 'volumes', 'masks', animal_folder + '.infra.nii.gz') + ' -a ' + \
              opj(dir_native_resol, animal_folder + '.infra_1.nii.gz') + ' -expr "amongst(a,3,42)"'
    spco([command], shell=True)

    command = '3dcalc -prefix ' + opj(dir_native_resol, 'ribbon_cort.nii.gz') + ' -a ' + opj(wb_native_dir, 'volumes',
                                                                                             'masks',
                                                                                             'ribbon.nii.gz') + ' -expr "amongst(a,3,42)"'
    spco([command], shell=True)
    command = '3dcalc -a ' + opj(dir_native_resol, 'ribbon_cort.nii.gz') + ' -b ' + opj(wb_native_dir, 'volumes',
                                                                                        'masks',
                                                                                        animal_folder + '.infra.nii.gz') + \
              ' -expr "a-b" -prefix ' + opj(wb_native_dir, 'volumes', 'masks', animal_folder + '.supra.nii.gz')
    spco([command], shell=True)

    os.remove(opj(dir_native_resol, animal_folder + '.l.infra.nii.gz'))
    os.remove(opj(dir_native_resol, animal_folder + '.r.infra.nii.gz'))
    os.remove(opj(dir_native_resol, animal_folder + '.infra_1.nii.gz'))
    os.remove(opj(dir_native_resol, 'ribbon_cort.nii.gz'))

    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(wb_native_dir,
                                                                                                       'volumes',
                                                                                                       'masks',
                                                                                                       'ribbon.nii.gz')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(wb_native_dir,
                                                                                                       'volumes',
                                                                                                       'masks',
                                                                                                       animal_folder + '.infra.nii.gz')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(wb_native_dir,
                                                                                                       'volumes',
                                                                                                       'masks',
                                                                                                       animal_folder + '.supra.nii.gz')
    spco([command], shell=True)



    #  for atlas space
    for h in range(0, 2):
        command = 'wb_command -create-signed-distance-volume ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[
            h] + '.white.native.surf.gii') + \
                  ' ' + Ref_norm + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white.native.nii.gz')
        spco([command], shell=True)
        command = 'wb_command -create-signed-distance-volume ' + opj(dir_atlas_resol,
                                                                     animal_folder + '.' + Hmin[h] + '.pial.native.surf.gii') + \
                  ' ' + Ref_norm + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol,
                                    animal_folder + '.' + Hmin[h] + '.white.native.nii.gz') + ' -thr 0 -bin -mul 255 ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_thr0.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_thr0.native.nii.gz') + \
                  ' -bin ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_thr0.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol,
                                    animal_folder + '.' + Hmin[h] + '.pial.native.nii.gz') + ' -uthr 0 -abs -bin -mul 255 ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial_uthr0.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial_uthr0.native.nii.gz') + ' -bin ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial_uthr0.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial_uthr0.native.nii.gz') + ' -mas ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_thr0.native.nii.gz') + \
                  ' -mul 255 ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.ribbon.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.ribbon.nii.gz') + ' -bin -mul  ' + str(
            RibbonValue[h][1]) + \
                  '  ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.ribbon.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol,
                                    animal_folder + '.' + Hmin[h] + '.white.native.nii.gz') + ' -uthr 0 -abs -bin -mul 255 ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_uthr0.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_uthr0.native.nii.gz') + ' -bin ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_uthr0.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol,
                                    animal_folder + '.' + Hmin[h] + '.white_uthr0.native.nii.gz') + ' -mul  ' + str(
            RibbonValue[h][1]) + \
                  '  ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_mask.native.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.ribbon.nii.gz') + ' -add ' + \
                  opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_mask.native.nii.gz') + \
                  ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.ribbon.nii.gz')
        spco([command], shell=True)

        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white.native.nii.gz'))
        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_thr0.native.nii.gz'))
        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_uthr0.native.nii.gz'))
        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.white_mask.native.nii.gz'))
        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial.native.nii.gz'))
        os.remove(opj(dir_atlas_resol, animal_folder + '.' + Hmin[h] + '.pial_uthr0.native.nii.gz'))

    command = 'fslmaths ' + opj(dir_atlas_resol, animal_folder + '.l.ribbon.nii.gz') + ' -add ' + opj(dir_atlas_resol,
                                                                                               animal_folder + '.r.ribbon.nii.gz') + \
              ' ' + opj(wb_atlas_dir, 'volumes', 'masks', 'ribbon.nii.gz')
    spco([command], shell=True)

    os.remove(opj(dir_atlas_resol, animal_folder + '.l.ribbon.nii.gz'))
    os.remove(opj(dir_atlas_resol, animal_folder + '.r.ribbon.nii.gz'))

    command = 'wb_command -volume-label-import ' + opj(wb_atlas_dir, 'volumes', 'masks',
                                                       'ribbon.nii.gz') + label_FS_all + \
              opj(wb_atlas_dir, 'volumes', 'masks', 'ribbon.nii.gz') + ' -drop-unused-labels'
    spco([command], shell=True)

    print('creation of the ribbon files: done!')




    #### atlas surfacing in volumique and in indiv space, only for macaque
    #      3.0 Yerkes parcellation from atlas space 32k to native space native k ##################################################################################################

    shutil.copy(Yerkes_label, opj(dir_atlas_32, animal_folder + '.Yerkes19.32_fs_LR.dlabel.nii'))

    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + \
              opj(dir_atlas_32, animal_folder + '.Yerkes19.32_fs_LR.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal_folder + '_fsaverage_LR_32.wb.spec') + ' INVALID ' + \
              opj(dir_atlas_32, animal_folder + '.Yerkes19.32_fs_LR.dlabel.nii')
    spco([command], shell=True)

    command = 'wb_command -cifti-separate ' + opj(dir_atlas_32, animal_folder + '.Yerkes19.32_fs_LR.dlabel.nii') + \
              ' COLUMN -label CORTEX_LEFT ' + opj(dir_atlas_32, 'Yerkes19_Parcellations_v2.l.32k_fs_LR.label.gii')
    spco([command], shell=True)
    command = 'wb_command -cifti-separate ' + opj(dir_atlas_32, animal_folder + '.Yerkes19.32_fs_LR.dlabel.nii') + \
              ' COLUMN -label CORTEX_RIGHT ' + opj(dir_atlas_32, 'Yerkes19_Parcellations_v2.r.32k_fs_LR.label.gii')
    spco([command], shell=True)

    command = 'wb_command -label-resample ' + opj(dir_atlas_32, 'Yerkes19_Parcellations_v2.l.32k_fs_LR.label.gii') + \
              ' ' + opj(dir_atlas_32, animal_folder + '.l.sphere.32k_fs_LR.surf.gii') + \
              ' ' + opj(dir_atlas_resol, animal_folder + '.l.sphere.reg.reg_LR.native.surf.gii') + \
              ' BARYCENTRIC ' + opj(dir_atlas_resol, animal_folder + '.l.Yerkes19_fs_LR.label.gii') + ' -largest'
    spco([command], shell=True)

    command = 'wb_command -label-resample ' + opj(dir_atlas_32, 'Yerkes19_Parcellations_v2.r.32k_fs_LR.label.gii') + \
              ' ' + opj(dir_atlas_32, animal_folder + '.r.sphere.32k_fs_LR.surf.gii') + \
              ' ' + opj(dir_atlas_resol, animal_folder + '.r.sphere.reg.reg_LR.native.surf.gii') + \
              ' BARYCENTRIC ' + opj(dir_atlas_resol, animal_folder + '.r.Yerkes19_fs_LR.label.gii') + ' -largest'
    spco([command], shell=True)

    command = 'wb_command -cifti-create-label ' + opj(dir_atlas_resol, animal_folder + '.Yerkes19.dlabel.nii') + \
              ' -left-label ' + opj(dir_atlas_resol, animal_folder + '.l.Yerkes19_fs_LR.label.gii') + ' -roi-left ' + opj(
        dir_atlas_resol, animal_folder + '.l.roi.shape.gii') + \
              ' -right-label ' + opj(dir_atlas_resol, animal_folder + '.r.Yerkes19_fs_LR.label.gii') + ' -roi-right ' + opj(
        dir_atlas_resol, animal_folder + '.r.roi.shape.gii')
    spco([command], shell=True)

    shutil.copy(opj(dir_atlas_resol, animal_folder + '.Yerkes19.dlabel.nii'),
                opj(dir_native_resol, animal_folder + '.Yerkes19.dlabel.nii'))
    command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(
        dir_atlas_resol, animal_folder + '.Yerkes19.dlabel.nii')
    spco([command], shell=True)
    command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol,
                                                    animal_folder + '_native_LR.wb.spec') + ' INVALID ' + opj(dir_native_resol,
                                                                                                       animal_folder + '.Yerkes19.dlabel.nii')
    spco([command], shell=True)






    # get markov 129 to .annot for FS
    #        To do so it should be the LAST map in the cifti file so we reorder it
    markov_order = [1, 2, 3, 4, 17, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 5]
    list_Yerkes = ['MW', 'LV00', 'FV91', 'PHT00', 'PFC', 'M132', 'B05', 'BB47', 'UD86', 'SP78+', 'LK02', 'FOA00', 'V6',
                   'MOD', 'KMA09', 'LV00_FOA00_PHT00', 'M129']
    with open(opj(REF_yerkes, 'markov_map.txt'), 'wb') as f:
        np.savetxt(f, markov_order, fmt='%.0f')

    command = 'wb_command -cifti-reorder ' + opj(dir_native_resol, animal_folder + '.Yerkes19.dlabel.nii') + ' ROW ' + \
              opj(REF_yerkes, 'markov_map.txt') + ' ' + opj(dir_native_resol, animal_folder + '.Yerkes19.reord.dlabel.nii')
    spco([command], shell=True)
    os.remove(opj(REF_yerkes, 'markov_map.txt'))

    command = 'wb_command -cifti-separate ' + opj(dir_native_resol, animal_folder + '.Yerkes19.reord.dlabel.nii') + \
              ' COLUMN -label CORTEX_LEFT ' + opj(dir_native_resol, animal_folder + '.l.Yerkes19.label.gii') + \
              ';wb_command -label-export-table ' + opj(dir_native_resol, animal_folder + '.l.Yerkes19.label.gii') + \
              ' ' + opj(REF_yerkes, 'yerkes_label.txt')
    spco([command], shell=True)

    command = 'wb_command -cifti-separate ' + opj(dir_native_resol, animal_folder + '.Yerkes19.reord.dlabel.nii') + \
              ' COLUMN -label CORTEX_RIGHT ' + opj(dir_native_resol, animal_folder + '.r.Yerkes19.label.gii')
    spco([command], shell=True)

    command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 --annot ' + \
              opj(dir_native_resol, animal_folder + '.l.Yerkes19.label.gii') + ' ' + \
              opj(FS_dir, animal_folder, 'surf', 'lh.white') + ' ' + \
              opj(FS_dir, animal_folder, 'label', 'lh.' + animal_folder + '_Y19_markov129.annot')
    spco([command], shell=True)

    command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 --annot ' + \
              opj(dir_native_resol, animal_folder + '.r.Yerkes19.label.gii') + ' ' + \
              opj(FS_dir, animal_folder, 'surf', 'rh.white') + ' ' + \
              opj(FS_dir, animal_folder, 'label', 'rh.' + animal_folder + '_Y19_markov129.annot')
    spco([command], shell=True)

    # build the 17 available Atlas as label volumes

    command = 'wb_command -label-to-volume-mapping ' + opj(dir_native_resol, animal_folder + '.l.Yerkes19.label.gii') + ' ' + \
              opj(dir_native_resol, animal_folder + '.l.white.surf.gii') + ' -ribbon-constrained ' + \
              opj(dir_native_resol, animal_folder + '.l.white.surf.gii') + ' ' + \
              opj(dir_native_resol, animal_folder + '.l.pial.surf.gii') + ' ' + \
              Ref_file + ' ' + opj(wb_native_dir, 'volumes', animal_folder + '_l_Yerkes19.nii.gz')
    spco([command], shell=True)

    command = 'wb_command -label-to-volume-mapping ' + opj(dir_native_resol, animal_folder + '.r.Yerkes19.label.gii') + ' ' + \
              opj(dir_native_resol, animal_folder + '.r.white.surf.gii') + ' -ribbon-constrained ' + \
              opj(dir_native_resol, animal_folder + '.r.white.surf.gii') + ' ' + \
              opj(dir_native_resol, animal_folder + '.r.pial.surf.gii') + ' ' + \
              Ref_file + ' ' + opj(wb_native_dir, 'volumes', animal_folder + '_r_Yerkes19.nii.gz')
    spco([command], shell=True)

    os.makedirs(opj(dir_native_resol, 'temp'))

    command = '3dTsplit4D -prefix ' + opj(dir_native_resol, 'temp', 'LEFT.nii.gz') + ' ' + opj(wb_native_dir, 'volumes',
                                                                                               animal_folder + '_l_Yerkes19.nii.gz')
    spco([command], shell=True)
    command = '3dTsplit4D -prefix ' + opj(dir_native_resol, 'temp', 'RIGHT.nii.gz') + ' ' + opj(wb_native_dir,
                                                                                                'volumes',
                                                                                                animal_folder + '_r_Yerkes19.nii.gz')
    spco([command], shell=True)

    for i in range(0, len(list_Yerkes)):
        if i < 10:
            command = 'fslmaths ' + opj(dir_native_resol, 'temp', 'RIGHT.0' + str(i) + '.nii.gz') + ' -add 2000 -mas ' + \
                      opj(dir_native_resol, 'temp', 'RIGHT.0' + str(i) + '.nii.gz') + ' ' + \
                      opj(dir_native_resol, 'temp', 'RIGHT.0' + str(i) + '.nii.gz')
            spco([command], shell=True)
            command = 'fslmaths ' + opj(dir_native_resol, 'temp', 'LEFT.0' + str(i) + '.nii.gz') + ' -add ' + \
                      opj(dir_native_resol, 'temp', 'RIGHT.0' + str(i) + '.nii.gz') + ' ' + \
                      opj(wb_native_dir, 'volumes', 'labels', animal_folder + '_' + list_Yerkes[i] + '.nii.gz')
            spco([command], shell=True)
    else:
        command = 'fslmaths ' + opj(dir_native_resol, 'temp', 'RIGHT.' + str(i) + '.nii.gz') + ' -add 2000 -mas ' + \
                  opj(dir_native_resol, 'temp', 'RIGHT.' + str(i) + '.nii.gz') + ' ' + \
                  opj(dir_native_resol, 'temp', 'RIGHT.' + str(i) + '.nii.gz')
        spco([command], shell=True)
        command = 'fslmaths ' + opj(dir_native_resol, 'temp', 'LEFT.' + str(i) + '.nii.gz') + ' -add ' + \
                  opj(dir_native_resol, 'temp', 'RIGHT.' + str(i) + '.nii.gz') + ' ' + \
                  opj(wb_native_dir, 'volumes', 'labels', animal_folder + '_' + list_Yerkes[i] + '.nii.gz')
        spco([command], shell=True)

    shutil.rmtree(opj(dir_native_resol, 'temp'))
    os.remove(opj(wb_native_dir, 'volumes', animal_folder + '_l_Yerkes19.nii.gz'))
    os.remove(opj(wb_native_dir, 'volumes', animal_folder + '_r_Yerkes19.nii.gz'))

    print('volumetric atlas of yerkes: done!')

   ##specific to an atlas
    #           Nuclei ##################################################################################################################################"

    nuclei_name = ['caudate', 'putamen', 'ventral_striatum', 'amygdala', 'hippocampus', 'thalamus', 'substancia_nigra',
                   'int_pallidum', 'ext_pallidum']
    nuclei_nb = [['2,12', '4,14,16', '6', '24', '22', '20,96,98,100,102,104,106,108', '18', '10', '8'],
                 ['1,11', '3,13,15', '5', '23', '21', '19,95,97,99,101,103,105,107', '17', '9', '7']]
    nuclei_sm = ['3', '2', '2', '4', '2', '3', '3', '2', '2']
    nuclei_struct = [['CAUDATE_LEFT', 'PUTAMEN_LEFT', 'DIENCEPHALON_VENTRAL_LEFT', 'AMYGDALA_LEFT', 'HIPPOCAMPUS_LEFT',
                      'THALAMUS_LEFT', 'INVALID', 'PALLIDUM_LEFT', 'OTHER'],
                     ['CAUDATE_RIGHT', 'PUTAMEN_RIGHT', 'DIENCEPHALON_VENTRAL_RIGHT', 'AMYGDALA_RIGHT',
                      'HIPPOCAMPUS_RIGHT', 'THALAMUS_RIGHT', 'INVALID', 'PALLIDUM_RIGHT', 'OTHER']]
    for i in range(0, 2):
        for j in range(0, len(nuclei_name)):
            command = 'mri_convert ' + opj(FS_dir, animal_folder, 'mri', 'atlas.mgz') + ' ' + opj(FS_dir, animal_folder, 'mri', 'atlas.nii.gz')
            spco([command], shell=True)
            command = '3dcalc -a ' + opj(FS_dir, animal_folder, 'mri', 'atlas.nii.gz') + ' -expr "amongst(a,' + nuclei_nb[i][j] + ')" -prefix ' + \
                      opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz')
            spco([command], shell=True)
            command = 'ImageMath 3 ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + ' FillHoles ' + \
                      opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + ' 2'
            spco([command], shell=True)
            command = 'ImageMath 3 ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + \
                      ' GetLargestComponent ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz')
            spco([command], shell=True)
            command = 'ImageMath 3 ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + \
                      ' G ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + ' 0.5'  # if not perfect try with 1 or 0.99
            spco([command], shell=True)
            command = 'ThresholdImage 3 ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + \
                      ' ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + ' 0.5 1 1 0'
            spco([command], shell=True)

            command = 'mri_tessellate ' + opj(FS_dir, animal_folder, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') + \
                      ' 1 ' + opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j])
            spco([command], shell=True)
            command = 'mris_extract_main_component ' + opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j]) + \
                      ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j])
            spco([command], shell=True)
            command = 'mris_euler_number ' + opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j])
            spco([command], shell=True)
            command = 'mris_remove_intersection ' + opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j]) + \
                      ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j])
            spco([command], shell=True)
            command = 'mris_smooth -nw -n ' + nuclei_sm[j] + ' ' + opj(FS_dir, animal_folder, 'surf',
                                                                       Hmin[i] + 'h_' + nuclei_name[j]) + \
                      ' ' + opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j] + '_smoothed')
            spco([command], shell=True)

            command = '/home/cgarin/Documents/0000_CODE/2023/EasyMRI_brain_CG/mris_convert_v6 ' + \
                      opj(FS_dir, animal_folder, 'surf', Hmin[i] + 'h_' + nuclei_name[j] + '_smoothed') + ' ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii')
            spco([command], shell=True)
            command = 'wb_command -set-structure ' + opj(dir_native_resol, animal_folder + '.' + Hmin[i] + '.' + nuclei_name[
                j] + '_smoothed.surf.gii') + \
                      ' ' + nuclei_struct[i][j] + ' -surface-type ANATOMICAL -surface-secondary-type INVALID'
            spco([command], shell=True)
            command = 'wb_command -surface-apply-affine ' + opj(dir_native_resol,
                                                                animal_folder + '.' + Hmin[i] + '.' + nuclei_name[
                                                                    j] + '_smoothed.surf.gii') + \
                      ' ' + opj(wb_native_dir, 'volumes', 'c_ras.mat') + ' ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii')
            spco([command], shell=True)
            command = 'wb_command -add-to-spec-file ' + opj(dir_native_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                      nuclei_struct[i][j] + ' ' + \
                      opj(dir_native_resol, animal_folder + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii')
            spco([command], shell=True)

            # for Yerkes
            command = 'wb_command -surface-apply-warpfield ' + opj(dir_native_resol,
                                                                   animal_folder + '.' + Hmin[i] + '.' + nuclei_name[
                                                                       j] + '_smoothed.surf.gii') + \
                      ' ' + ref_transfo1 + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[i] + '.' + nuclei_name[
                j] + '_smoothed.transient.surf.gii')
            spco([command], shell=True)
            command = 'wb_command -surface-apply-affine ' + \
                      opj(dir_atlas_resol,
                          animal_folder + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.transient.surf.gii') + \
                      ' ' + ref_transfo2 + ' ' + opj(dir_atlas_resol, animal_folder + '.' + Hmin[i] + '.' + nuclei_name[
                j] + '_smoothed.native.surf.gii')
            spco([command], shell=True)
            command = 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal_folder + '_native_LR.wb.spec') + ' ' + \
                      nuclei_struct[i][j] + ' ' + \
                      opj(dir_atlas_resol, animal_folder + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.native.surf.gii')
            spco([command], shell=True)

