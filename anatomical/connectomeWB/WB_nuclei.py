#import
import os
import ants

opj = os.path.join
opd = os.path.dirname

from Tools import run_cmd

# Constant
Hmin        = ['l','r']
Hcap        = ['L','R']
Hsurf       = ['left','right']

nuclei_name = ['caudate', 'putamen', 'ventral_striatum', 'amygdala', 'hippocampus', 'thalamus', 'ext_pallidum']

nuclei_struct = [['CAUDATE_LEFT', 'PUTAMEN_LEFT', 'DIENCEPHALON_VENTRAL_LEFT', 'AMYGDALA_LEFT', 
                  'HIPPOCAMPUS_LEFT','THALAMUS_LEFT', 'PALLIDUM_LEFT'],
                 ['CAUDATE_RIGHT', 'PUTAMEN_RIGHT', 'DIENCEPHALON_VENTRAL_RIGHT', 'AMYGDALA_RIGHT',
                  'HIPPOCAMPUS_RIGHT', 'THALAMUS_RIGHT', 'PALLIDUM_RIGHT']]
nuclei_sm = ['3', '2', '2', '4', '2', '3', '2']

def native(cmd_mris_convert,FS_dir,dir_surf,animal,species,atlas,diary_name,sing_fs,sing_wb):

    nl = 'Creation of the nuclei surfaces for the native brain'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    dir_vol = opj(opd(opd(dir_surf)),'volumes')

    if atlas =='BAL2024':
        nuclei_nb = [[[2,12], [4,14,16], [6], [24], [22], [20,96,98,100,102,104,106,108], [8]],
                     [[1,11], [3,13,15], [5], [23], [21], [19,95,97,99,101,103,105,107], [7]]]

    atlas_file_mgz = opj(FS_dir, animal, 'mri', animal + '_' + atlas + '.mgz')
    atlas_file_nii = opj(FS_dir, animal, 'mri', animal + '_' + atlas + '.nii.gz')

    cmd = (sing_fs + 'mri_convert ' + atlas_file_mgz + ' ' + atlas_file_nii)
    run_cmd.run(cmd,diary_name)

    for i in range(2):
        for j in range(len(nuclei_name)):

            atlas_img = ants.image_read(atlas_file_nii)
            img = ants.mask_image(atlas_img, atlas_img, nuclei_nb[i][j])
            mask_img = ants.image_clone(img)
            mask_img[mask_img > 0] = 1
            IMG = ants.morphology(mask_img, operation='close',radius=2)
            IMG = ants.iMath(IMG, operation='GetLargestComponent')
            IMG = ants.smooth_image(IMG, 0.5)
            IMG = ants.threshold_image(IMG, 0.5, 1, 1, 0, True)
            ants.image_write(IMG, opj(FS_dir, animal, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz'), ri=False)

            cmd = (sing_fs + 'mri_tessellate ' + opj(FS_dir, animal, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz') +
                   ' 1 ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j]))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_fs + 'mris_extract_main_component ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j]) +
                   ' ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j]))
            run_cmd.run(cmd,diary_name)

            cmd = sing_fs + 'mris_euler_number ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j])
            run_cmd.run(cmd,diary_name)

            cmd = (sing_fs + 'mris_remove_intersection ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j]) +
                   ' ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j]))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_fs + 'mris_smooth -nw -n ' + nuclei_sm[j] +
                   ' ' + opj(FS_dir, animal, 'surf',Hmin[i] + 'h_' + nuclei_name[j]) +
                   ' ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j] + '_smoothed'))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_fs + cmd_mris_convert + ' ' + opj(FS_dir, animal, 'surf', Hmin[i] + 'h_' + nuclei_name[j] + '_smoothed') +
                   ' ' + opj(dir_surf, animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii'))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_surf, animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii') +
                   ' ' + nuclei_struct[i][j] + ' -surface-type ANATOMICAL -surface-secondary-type INVALID')
            run_cmd.run(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -surface-apply-affine ' + opj(dir_surf,animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii') +
                   ' ' + opj(dir_vol, 'c_ras.mat') +
                   ' ' + opj(dir_surf, animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii'))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_surf, animal + '_native_LR.wb.spec') + ' ' + nuclei_struct[i][j] +
                   ' ' + opj(dir_surf, animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii'))
            run_cmd.run(cmd,diary_name)

            os.remove(opj(FS_dir, animal, 'mri', Hmin[i] + 'h_' + nuclei_name[j] + '.nii.gz'))


def template(dir_surf1,dir_surf2, animal, ref_transfo1, ref_transfo2, diary_name,sing_wb):

    nl = 'Creation of the nuclei surfaces for the template'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    for i in range(2):
        for j in range(len(nuclei_name)):
            cmd = (sing_wb + 'wb_command -surface-apply-warpfield ' +
                   opj(dir_surf1,animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.surf.gii') +
                   ' ' + ref_transfo1 +
                   ' ' + opj(dir_surf2, animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.transient.surf.gii'))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -surface-apply-affine ' +
                  opj(dir_surf2,animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.transient.surf.gii') +
                  ' ' + ref_transfo2 +
                   ' ' + opj(dir_surf2, animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.native.surf.gii'))
            run_cmd.run(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_surf2, animal + '_native_LR.wb.spec') +
                   ' ' + nuclei_struct[i][j] +
                   ' ' + opj(dir_surf2, animal + '.' + Hmin[i] + '.' + nuclei_name[j] + '_smoothed.native.surf.gii'))
            run_cmd.run(cmd,diary_name)

