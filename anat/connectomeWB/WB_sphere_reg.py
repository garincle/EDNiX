#import
import os
import glob
import shutil

from Tools import run_cmd
from anat.connectomeWB import WB_label

opj = os.path.join
opb = os.path.basename
opi = os.path.isfile

# Constant
SurfaceSmoothingFWHM  = 2
Hmin        = ['l','r']
Hcap        = ['L','R']
Hsurf       = ['left','right']
CORTEX      = ['CORTEX_LEFT','CORTEX_RIGHT']

def copyref(FS_mesh,dir_atlas_64,dir_atlas_32,wb_balsa_labels,atlasfile,animal,diary_name,sing_wb):

    nl = '# Copy of the std sphere References Files\n'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    for h in range(2):

        #  High res
        if not opi(opj(dir_atlas_64, animal + '.sphere.164k_fs_' + Hmin[h] + '.surf.gii')):
            shutil.copyfile(opj(FS_mesh, 'fs_' + Hcap[h], 'fsaverage.' + Hcap[h] + '.sphere.164k_fs_' + Hcap[h] + '.surf.gii'),
                        opj(dir_atlas_64, animal + '.sphere.164k_fs_' + Hmin[h] + '.surf.gii'))
        if not opi(opj(dir_atlas_64, animal + '.def_sphere.164k_fs_' + Hmin[h] + '.surf.gii')):
            shutil.copyfile(opj(FS_mesh, 'fs_' + Hcap[h],'fs_' + Hcap[h] + '-to-fs_LR_fsaverage.' + Hcap[h] + '_LR.spherical_std.164k_fs_' + Hcap[h] + '.surf.gii'),
                        opj(dir_atlas_64, animal + '.def_sphere.164k_fs_' + Hmin[h] + '.surf.gii'))
        if not opi(opj(dir_atlas_64, animal + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii')):
            shutil.copyfile(opj(FS_mesh, 'fsaverage.' + Hcap[h] + '_LR.spherical_std.164k_fs_LR.surf.gii'),
                        opj(dir_atlas_64, animal + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii'))

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal + '_fsaverage_LR_164k.wb.spec') +
               ' ' + CORTEX[h] +
               ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)
        if not opi(opj(dir_atlas_64, animal + '.' + Hmin[h] + '.atlasroi.164k_fs_LR.shape.gii')):
            shutil.copyfile(opj(FS_mesh, Hcap[h] + '.atlasroi.164k_fs_LR.shape.gii'),
                            opj(dir_atlas_64, animal + '.' + Hmin[h] + '.atlasroi.164k_fs_LR.shape.gii'))
        if not opi(opj(dir_atlas_64, animal + '.' + Hmin[h] + '.refsulc.164k_fs_LR.shape.gii')):
            shutil.copyfile(opj(FS_mesh, Hcap[h] + '.refsulc.164k_fs_LR.shape.gii'),
                            opj(dir_atlas_64, animal + '.' + Hmin[h] + '.refsulc.164k_fs_LR.shape.gii'))

        # low res
        if not opi(opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii')):
            shutil.copyfile(opj(FS_mesh, Hcap[h] + '.sphere.32k_fs_LR.surf.gii'),
                            opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii'))
        if not opi(opj(dir_atlas_32, animal + '.' + Hmin[h] + '.atlasroi.32k_fs_LR.shape.gii')):
            shutil.copyfile(opj(FS_mesh, Hcap[h] + '.atlasroi.32k_fs_LR.shape.gii'),
                            opj(dir_atlas_32, animal + '.' + Hmin[h] + '.atlasroi.32k_fs_LR.shape.gii'))
        if not opi(opj(wb_balsa_labels, atlasfile)):
            shutil.copyfile(opj(FS_mesh, atlasfile),opj(wb_balsa_labels, atlasfile))


def warp(animal,WB_surf,dir_native_resol,dir_atlas_resol,list_atlas,ref_transfo1,ref_transfo2,conv_voxmm,cmd_mris_convert,FS_dir,diary_name,sing_wb,sing_fs):

    nl ='# transform surfaces into atlas space\n'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    for h in range(2):
        labels = glob.glob(opj(FS_dir, animal, 'label', Hmin[h] + 'h.' + animal + '_*.annot'))
        for FSs in range(len(WB_surf)):
            cmd = (sing_wb + 'wb_command -surface-apply-warpfield ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.' + WB_surf[FSs] + '.surf.gii') +
               ' ' + ref_transfo1 + ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + WB_surf[FSs] + '.transient.surf.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -surface-apply-affine ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + WB_surf[FSs] + '.transient.surf.gii') +
                   ' ' + ref_transfo2 + ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + WB_surf[FSs] + '.native.surf.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal + '_native_LR.wb.spec') + ' ' + CORTEX[h] +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + WB_surf[FSs] + '.native.surf.gii'))
            run_cmd.wb(cmd,diary_name)


        nl = 'Creation of the ' + animal + ' spec file template step 2'
        run_cmd.msg(nl,diary_name,'OKGREEN')

        # from midthickness generation of inflated and very inflated surface
        cmd = (sing_wb + 'wb_command -file-information ' + opj(dir_atlas_resol,
                                                    animal + '.' + Hmin[h] + '.midthickness.native.surf.gii'))
        #       ' | grep "Number of Vertices:" | cut -f2 -d: | tr -d "[:space:]"')
        out,_=run_cmd.wb(cmd,diary_name)
        result = out.decode("utf-8").split('\n')
        for item in result:
            if item.find("Number of Vertices:") != -1:
                nb = int(item.split(' ')[-1])

        NativeInflationScale = nb * conv_voxmm / 32492  # 32492 (so 32k) for human or MY19 low res or  74k for F99

        cmd = (sing_wb + 'wb_command -surface-generate-inflated ' +
               opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.midthickness.native.surf.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.inflated.surf.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.very_inflated.surf.gii') +
               ' -iterations-scale ' + str(NativeInflationScale))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal + '_native_LR.wb.spec') +
               ' ' + CORTEX[h] +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.inflated.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal + '_native_LR.wb.spec') +
               ' ' + CORTEX[h] +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.very_inflated.surf.gii'))
        run_cmd.wb(cmd,diary_name)


        # sulc thickness and curvature for native brain

        for SHAPE_FS, SHAPE_WB in zip(['sulc', 'curv'],
                                      ['sulc', 'curvature']):
            cmd = (sing_fs + cmd_mris_convert + ' -c ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.' + SHAPE_FS) +
                   ' ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.white') +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii'))
            run_cmd.wb(cmd,diary_name)


            cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_atlas_resol,
                                                      animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' ' + CORTEX[h])
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -metric-math "var * -1" ' + opj(dir_atlas_resol,
                                                               animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' -var var ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii'))
            run_cmd.do(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_atlas_resol,
                                                      animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' -map 1 ' + animal + '_' + Hmin[h] + '_' + SHAPE_WB)
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_atlas_resol,
                                                       animal + '.' + Hmin[h] + '.' + SHAPE_WB + '.shape.gii') +
                   ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true')
            run_cmd.wb(cmd,diary_name)


        cmd = (sing_fs + cmd_mris_convert + ' -c ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.thickness') +
               ' ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.white') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' ' + CORTEX[h])
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-math "var * -1" ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.thickness.shape.gii') +
              ' -var var ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' -map 1 ' + animal + '_' + Hmin[h] + '_Thickness')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
              ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-neg true -disp-zero true')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-math "abs(thickness)" ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' -var thickness ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' MODE_AUTO_SCALE_PERCENTAGE -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-math "thickness > 0" ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' -var thickness ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii'))
        run_cmd.do(cmd,diary_name)


        cmd = (sing_wb + 'wb_command -metric-fill-holes ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.midthickness.native.surf.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-remove-islands ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.midthickness.native.surf.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
        run_cmd.wb(cmd, diary_name)

        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.roi.shape.gii') +
               ' -map 1 ' + animal + '_' + Hmin[h] + '_ROI')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-dilate ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.midthickness.native.surf.gii') +
               ' 10 ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.thickness.shape.gii') + ' -nearest')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-dilate ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.curvature.shape.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.midthickness.native.surf.gii') +
               ' 10 ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.curvature.shape.gii') + ' -nearest')
        run_cmd.wb(cmd,diary_name)


        # parcellations (.annot)

        for j in labels:
            name = opb(j).split('.')[1].split('_')[-1]
            WB_label.surfWB(cmd_mris_convert, FS_dir, animal, name, opj(dir_atlas_resol), '.native', h, diary_name, sing_fs, sing_wb)

        for a in range(len(list_atlas[0])):
            if list_atlas[3][a] == 1:
                if list_atlas[2][a] > 1:
                    WB_label.mergelabel(animal, opj(dir_atlas_resol), '.native', h, diary_name, list_atlas[0][a], list_atlas[2][a], sing_wb)


def sphere(cmd_mris_convert,animal,FS_dir,dir_atlas_resol,diary_name,sing_fs,sing_wb):

    nl ='sphere and registration\n'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    for h in range(2):
        # native sphere ###############################################################################################

        nl = ' building the ' + Hsurf[h] + ' native sphere surface '
        run_cmd.msg(nl,diary_name, 'OKGREEN')

        cmd = (sing_fs + cmd_mris_convert + ' ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.sphere') +
        ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.native.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.sphere.native.surf.gii') +
               ' ' + CORTEX[h] + ' -surface-type SPHERICAL -surface-secondary-type GRAY_WHITE')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal + '_native_LR.wb.spec') +
               ' ' + CORTEX[h] +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.native.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        # native sphere for registration #############################################################################

        nl = ' building the ' + Hsurf[h] + ' native reg sphere surface '
        run_cmd.msg(nl,diary_name, 'OKGREEN')

        cmd = (sing_fs + cmd_mris_convert + ' ' + opj(FS_dir, animal, 'surf', Hmin[h] + 'h.sphere.reg') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.native.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -set-structure ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.native.surf.gii') +
               ' ' + CORTEX[h] + ' -surface-type SPHERICAL -surface-secondary-type GRAY_WHITE')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_resol, animal + '_native_LR.wb.spec') +
               ' ' + CORTEX[h] +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.native.surf.gii'))
        run_cmd.wb(cmd, diary_name)


def reg(animal,dir_atlas_resol,dir_atlas_64,diary_name,sing_wb):

    nl ='# sphere to sphere registration\n'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    # to ref surface
    for h in range(2):

        cmd = (sing_wb + 'wb_command -surface-sphere-project-unproject ' +
               opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.native.surf.gii') +
              ' ' + opj(dir_atlas_64, animal + '.sphere.164k_fs_' + Hmin[h] + '.surf.gii') +
              ' ' + opj(dir_atlas_64, animal + '.def_sphere.164k_fs_' + Hmin[h] + '.surf.gii') +
              ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        # Calculate and visualise distortion.......................................
        cmd = (sing_wb + 'wb_command -surface-vertex-areas ' +
               opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.sphere.native.surf.gii') +
              ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.native.shape.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -surface-vertex-areas ' +
               opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.shape.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-math "ln(spherereg / sphere) / ln(2)" ' +
              ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.ArealDistortion_FS.shape.gii') +
              ' -var sphere ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.native.shape.gii') +
              ' -var spherereg ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.shape.gii'))
        run_cmd.do(cmd,diary_name)

        os.remove(opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.native.shape.gii'))
        os.remove(opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.shape.gii'))

        cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.ArealDistortion_FS.shape.gii') +
              ' -map 1 ' + animal + '_' + Hmin[h] + '_ArealDistortion_FS')
        run_cmd.wb(cmd,diary_name)
        
        cmd = (sing_wb + 'wb_command -metric-palette ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.ArealDistortion_FS.shape.gii') +
               ' MODE_AUTO_SCALE -palette-name ROY_BIG_BL -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE -1 1')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -surface-distortion ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.sphere.native.surf.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
               ' ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.EdgeDistortion_FS.shape.gii') +
               ' -edge-method')
        run_cmd.wb(cmd,diary_name)
        
        cmd = (sing_wb + 'wb_command -surface-distortion ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.sphere.native.surf.gii') +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
               ' ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.Strain_FS.shape.gii') +
               ' -local-affine-method')
        run_cmd.wb(cmd,diary_name)
        
        cmd =(sing_wb + 'wb_command -metric-merge ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.StrainJ_FS.shape.gii') +
              ' -metric ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.Strain_FS.shape.gii') +
              ' -column 1')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-math "ln(var) / ln(2)" ' +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.StrainJ_FS.shape.gii') +
               ' -var var ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.StrainJ_FS.shape.gii'))
        run_cmd.do(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-merge ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.StrainR_FS.shape.gii') +
        ' -metric ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.Strain_FS.shape.gii') +
        ' -column 2')
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -metric-math "ln(var) / ln(2)" ' +
               ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.StrainR_FS.shape.gii') +
               ' -var var ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.StrainR_FS.shape.gii'))
        run_cmd.do(cmd,diary_name)

        os.remove(opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.Strain_FS.shape.gii'))

        
def surf_resamp164(animal,dir_native_resol,dir_atlas_resol,dir_atlas_64,conv_voxmm,FS_mesh,list_atlas,diary_name,sing_wb):

    nl ='# resampling the surfaces in high resolution : 164k\n'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    for h in range(2):

        for surf in ['white','pial','midthickness']:

            cmd = (sing_wb + 'wb_command -surface-resample ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.' + surf + '.native.surf.gii') +
                  ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                  ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') +
                  ' BARYCENTRIC ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.' + surf + '.164k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal + '_fsaverage_LR_164k.wb.spec') +
                   ' ' + CORTEX[h] +' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.' + surf + '.164k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)


        # inflation ##############################################################################################"

        HighResInflationScale = int(164) * conv_voxmm / 32

        cmd = (sing_wb + 'wb_command -surface-generate-inflated ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii') +
              ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.inflated.164k_fs_LR.surf.gii') +
              ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.very_inflated.164k_fs_LR.surf.gii') +
              ' -iterations-scale ' + str(HighResInflationScale))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal + '_fsaverage_LR_164k.wb.spec') +
               ' ' + CORTEX[h] + ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.inflated.164k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal + '_fsaverage_LR_164k.wb.spec') +
               ' ' +  CORTEX[h] + ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.very_inflated.164k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        # flat maps ##############################################################################################
        if opi(opj(FS_mesh, 'colin.cerebral.' + Hcap[h] + '.flat.164k_fs_LR.surf.gii'))==True:
            if not opi(opj(dir_atlas_64, animal + '.' + Hmin[h] + '.flat.164k.surf.gii')):
                shutil.copyfile(opj(FS_mesh, 'colin.cerebral.' + Hcap[h] + '.flat.164k_fs_LR.surf.gii'),
                                opj(dir_atlas_64, animal + '.' + Hmin[h] + '.flat.164k.surf.gii'))
            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_64, animal + '_fsaverage_LR_164k.wb.spec') +
                   ' ' + CORTEX[h] +
                   ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.flat.164k.surf.gii'))
            run_cmd.wb(cmd,diary_name)

        for surf in ['curvature', 'thickness']:

            cmd = (sing_wb + 'wb_command -metric-resample ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.' + surf + '.shape.gii') +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                   ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') +
                  ' ADAP_BARY_AREA ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.' + surf + '.164k_fs_LR.shape.gii') +
                  ' -area-surfs ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
                  ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii') +
                  ' -current-roi ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -metric-mask ' + opj(dir_atlas_64,animal + '.' + Hmin[h] + '.' + surf + '.164k_fs_LR.shape.gii') +
                   ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.atlasroi.164k_fs_LR.shape.gii') +
                   ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.' + surf + '.164k_fs_LR.shape.gii'))
            run_cmd.wb(cmd,diary_name)

        for surf in ['sulc', 'ArealDistortion_FS','EdgeDistortion_FS','StrainJ_FS','StrainR_FS']:

            cmd = (sing_wb + 'wb_command -metric-resample ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.' + surf + '.shape.gii') +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                   ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') +
                   ' ADAP_BARY_AREA ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.' + surf + '.164k_fs_LR.shape.gii') +
                   ' -area-surfs ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
                   ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.midthickness.164k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)

        # labels
        for j, k in zip(list_atlas[0], list_atlas[3]):
            if k == 1:
                cmd = (sing_wb + 'wb_command -label-resample ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + j + '.native.label.gii') +
                      ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                      ' ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.sphere.164k_fs_LR.surf.gii') +
                      ' BARYCENTRIC ' + opj(dir_atlas_64, animal + '.' + Hmin[h] + '.' + j + '.164k_fs_LR.label.gii') + ' -largest')
                run_cmd.wb(cmd,diary_name)


def surf_resamp32(animal,dir_native_resol,dir_atlas_resol,dir_atlas_32,conv_voxmm,FS_mesh,list_atlas,diary_name,sing_wb):

    nl ='# resampling the surfaces in low resolution : 32k\n'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    for h in range(2):

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal + '_fsaverage_LR_32k.wb.spec') +
               ' ' + CORTEX[h] + ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        for surf in ['white','pial','midthickness']:

            cmd = (sing_wb + 'wb_command -surface-resample ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.' + surf + '.native.surf.gii') +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') +
                   ' BARYCENTRIC ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.' + surf + '.32k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd =(sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal + '_fsaverage_LR_32k.wb.spec') +
                  ' ' + CORTEX[h] + ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.' + surf + '.32k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)

        # inflation ##############################################################################################

        LowResInflationScale = int(32) * conv_voxmm / 32

        cmd = (sing_wb + 'wb_command -surface-generate-inflated ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.midthickness.32k_fs_LR.surf.gii') +
               ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.inflated.32k_fs_LR.surf.gii') +
               ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.very_inflated.32k_fs_LR.surf.gii') +
               ' -iterations-scale ' + str(LowResInflationScale))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal + '_fsaverage_LR_32k.wb.spec') +
               ' ' + CORTEX[h] + ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.inflated.32k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd=(sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal + '_fsaverage_LR_32k.wb.spec') +
             ' ' + CORTEX[h] + ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.very_inflated.32k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        # flat surfaces ###########################################################################################
        if opi(opj(FS_mesh, 'colin.cerebral.' + Hcap[h] + '.flat.32k_fs_LR.surf.gii'))==True:
            if not opi(opj(dir_atlas_32, animal + '.' + Hmin[h] + '.flat.32k.surf.gii')):
                shutil.copyfile(opj(FS_mesh, 'colin.cerebral.' + Hcap[h] + '.flat.32k_fs_LR.surf.gii'),
                                opj(dir_atlas_32, animal + '.' + Hmin[h] + '.flat.32k.surf.gii'))
            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_atlas_32, animal + '_fsaverage_LR_32k.wb.spec') +
                   ' ' + CORTEX[h] + ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.flat.32k.surf.gii'))
            run_cmd.wb(cmd,diary_name)


        for surf in ['curvature', 'thickness']:

            cmd = (sing_wb + 'wb_command -metric-resample ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.' + surf + '.shape.gii') +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') +
                   ' ADAP_BARY_AREA ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.' + surf + '.32k_fs_LR.shape.gii') +
                   ' -area-surfs ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.midthickness.32k_fs_LR.surf.gii') +
                   ' -current-roi ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.roi.shape.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -metric-mask ' + opj(dir_atlas_32,animal + '.' + Hmin[h] + '.' + surf + '.32k_fs_LR.shape.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.atlasroi.32k_fs_LR.shape.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.' + surf + '.32_fs_LR.shape.gii'))
            run_cmd.wb(cmd,diary_name)

        for surf in ['sulc', 'ArealDistortion_FS', 'EdgeDistortion_FS', 'StrainJ_FS', 'StrainR_FS']:

            cmd = (sing_wb + 'wb_command -metric-resample ' + opj(dir_atlas_resol,animal + '.' + Hmin[h] + '.' + surf + '.shape.gii') +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') +
                   ' ADAP_BARY_AREA ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.' + surf + '.32k_fs_LR.shape.gii') +
                   ' -area-surfs ' + opj(dir_native_resol, animal + '.' + Hmin[h] + '.midthickness.surf.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.midthickness.32k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)

        # labels
        for j,k in zip(list_atlas[0],list_atlas[3]):
            if k==1:
                cmd = (sing_wb + 'wb_command -label-resample ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.' + j + '.native.label.gii') +
                      ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                      ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') +
                      ' BARYCENTRIC ' +
                       opj(dir_atlas_32, animal + '.' + Hmin[h] + '.' + j + '.32k_fs_LR.label.gii') + ' -largest')
                run_cmd.wb(cmd,diary_name)


def surf_resamp_native(animal, dir_native_resol, dir_atlas_resol,dir_native_32,dir_atlas_32, conv_voxmm, diary_name,sing_wb):

    nl ='# resampling the 32k surfaces in native space\n'
    run_cmd.msg(nl,diary_name, 'OKGREEN')

    for h in range(2):

        for surf in ['white', 'pial', 'midthickness']:

            cmd = (sing_wb + 'wb_command -surface-resample ' + opj(dir_native_resol,animal + '.' + Hmin[h] + '.' + surf + '.surf.gii') +
                   ' ' + opj(dir_atlas_resol, animal + '.' + Hmin[h] + '.sphere.reg.reg_LR.native.surf.gii') +
                   ' ' + opj(dir_atlas_32, animal + '.' + Hmin[h] + '.sphere.32k_fs_LR.surf.gii') +
                   ' BARYCENTRIC ' + opj(dir_native_32, animal + '.' + Hmin[h] + '.' + surf + '.32k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal + '_fsaverage_LR_32k.wb.spec') +
                   ' ' + CORTEX[h] + ' ' + opj(dir_native_32, animal + '.' + Hmin[h] + '.' + surf + '.32k_fs_LR.surf.gii'))
            run_cmd.wb(cmd,diary_name)

        # inflation ###############################################################################################
        LowResInflationScale = int(32) * conv_voxmm / 32
        cmd = (sing_wb + 'wb_command -surface-generate-inflated ' + opj(dir_native_32, animal + '.' + Hmin[h] + '.midthickness.32k_fs_LR.surf.gii') +
              ' ' + opj(dir_native_32, animal + '.' + Hmin[h] + '.inflated.32k_fs_LR.surf.gii') +
              ' ' + opj(dir_native_32, animal + '.' + Hmin[h] + '.very_inflated.32k_fs_LR.surf.gii') +
              ' -iterations-scale ' + str(LowResInflationScale))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal + '_fsaverage_LR_32k.wb.spec') +
               ' ' + CORTEX[h] + ' ' + opj(dir_native_32, animal + '.' + Hmin[h] + '.inflated.32k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)

        cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_native_32, animal + '_fsaverage_LR_32k.wb.spec') +
               ' ' + CORTEX[h] + ' ' + opj(dir_native_32, animal + '.' + Hmin[h] + '.very_inflated.32k_fs_LR.surf.gii'))
        run_cmd.wb(cmd,diary_name)
            

def cifti(animal ,dir_atlas_resol,dir_atlas_32,dir_atlas_64,diary_name,sing_wb):

    nl='# create the differents cifti files\n'
    run_cmd.msg(nl,diary_name,'OKGREEN')

    scale = [' MODE_AUTO_SCALE_PERCENTAGE ',' MODE_USER_SCALE ']

    palette = [' -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true',
               ' -pos-percent 4 96 -interpolate true -palette-name videen_style -disp-pos true -disp-neg false -disp-zero false',
               ' -pos-user 0 1 -neg-user 0 -1 -interpolate true -palette-name ROY-BIG-BL -disp-pos true -disp-neg true -disp-zero false']


    for dir_surf,suffix,roi,specf,gii in zip([dir_atlas_resol,dir_atlas_32,dir_atlas_64],
                                             ['','.32k_fs_LR','.164k_fs_LR'],
                                             ['roi','atlasroi','atlasroi'],
                                             ['_native_LR','_fsaverage_LR_32k','_fsaverage_LR_164k'],
                                             ['.native','.32k_fs_LR','.164k_fs_LR']):

        for surf ,color,mode_scale in zip(['curvature', 'thickness'],[palette[0],palette[1]],[scale[0],scale[0]]):

            cmd = (sing_wb + 'wb_command -cifti-create-dense-scalar ' + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii') +
                   ' -left-metric ' +  opj(dir_surf, animal + '.l.' + surf + suffix + '.shape.gii') +
                   ' -roi-left ' +     opj(dir_surf, animal + '.l.' + roi + suffix + '.shape.gii') +
                   ' -right-metric ' + opj(dir_surf, animal + '.r.' + surf + suffix + '.shape.gii') +
                   ' -roi-right ' +    opj(dir_surf, animal + '.r.' + roi + suffix + '.shape.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_surf,animal + '.' + surf + suffix + '.dscalar.nii') +
                   ' -map 1 ' + animal + '_' +  surf)
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -cifti-palette ' + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii') +
                   mode_scale + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii') + color)
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_surf,animal + specf + '.wb.spec') +
                   ' INVALID ' + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii'))
            run_cmd.wb(cmd,diary_name)


        for surf ,color,mode_scale in zip(['sulc', 'ArealDistortion_FS', 'EdgeDistortion_FS', 'StrainJ_FS', 'StrainR_FS'],
                                          [palette[0],palette[2],palette[2],palette[2],palette[2]],
                                          [scale[0],scale[1],scale[1],scale[1],scale[1]]):

            cmd = (sing_wb + 'wb_command -cifti-create-dense-scalar ' + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii') +
                   ' -left-metric ' + opj(dir_surf, animal + '.l.' + surf + suffix + '.shape.gii') +
                   ' -right-metric ' + opj(dir_surf, animal + '.r.' + surf + suffix + '.shape.gii'))
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -set-map-names ' + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii') +
                   ' -map 1 ' + animal + '_' + surf)
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -cifti-palette ' + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii') +
                   mode_scale + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii') + color)
            run_cmd.wb(cmd,diary_name)

            cmd = (sing_wb + 'wb_command -add-to-spec-file ' + opj(dir_surf, animal + specf + '.wb.spec') +
                   ' INVALID ' + opj(dir_surf, animal + '.' + surf + suffix + '.dscalar.nii'))
            run_cmd.wb(cmd,diary_name)

        labels = glob.glob(opj(dir_surf, animal + '.l.*.label.gii'))

        for j in labels:
            name = opb(j).split('.')[2]
            run_cmd.msg(name,diary_name,'OKGREEN')
            WB_label.cifti(animal,
                           name,
                           opj(dir_surf),
                           gii,
                           roi + suffix,
                           opj(dir_surf, animal + specf + '.wb.spec'),
                           diary_name,
                           sing_wb)
