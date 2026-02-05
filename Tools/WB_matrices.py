import os
import ants
import numpy as np

opj = os.path.join

from tools import run_cmd

def ants2wb(MTX,dir_transfo,REF_name,diary_name,sing_wb):
    affine = ants.read_transform(MTX)
    PTX = affine.apply_to_point((0, 0, 0))
    test = [None] * 4, [None] * 4, [None] * 4, [None] * 4

    test[0][:-1] = affine.parameters[0:3]
    test[0][2] = test[0][2] * -1
    test[0][3] = PTX[0] * -1

    test[1][:-1] = affine.parameters[3:6]
    test[1][2] = test[1][2] * -1
    test[1][3] = PTX[1] * -1

    test[2][:-1] = affine.parameters[6:9]
    test[2][0] = test[2][0] * -1
    test[2][1] = test[2][1] * -1
    test[2][3] = PTX[2]

    test[3][:] = 0, 0, 0, 1

    np.savetxt(opj(dir_transfo, 'anat2' + REF_name + '_for_surface.mat'), test, fmt='%10.8f', delimiter=' ')

    cmd = (sing_wb + 'wb_command -convert-affine -from-world ' +
           opj(dir_transfo,'anat2' + REF_name + '_for_surface.mat') +
           ' -inverse -to-world ' + opj(dir_transfo, 'affine_' + REF_name + '.nii.gz'))
    run_cmd.run(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -convert-warpfield -from-itk ' +
           opj(dir_transfo,'anat2' + REF_name + '_SyN_1InverseWarp.nii.gz') +
           ' -to-world ' + opj(dir_transfo, 'standard_' + REF_name + '.nii.gz'))
    run_cmd.run(cmd, diary_name)