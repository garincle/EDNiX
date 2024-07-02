import os
import subprocess
import glob
import shutil
import sys

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput
from fonctions.extract_filename import extract_filename


def to_anat_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    nb_run, RS, n_for_ANTS, TR, overwrite):
    #########################################################################################################
    ################################### registration to anat space ##########################################
    #########################################################################################################

    #### transfo to brain template space (centred) are
    #    mvt_shft_INV_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1InverseWarp.nii.gz') + \
    #' -t ' + str([opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_0GenericAffine.mat'), 1])

    ############################### ############################### ############################### 
    ############################### apply transfo to anat space to each volume of each func image #
    ############################### ############################### ############################### 

    mvt_shft_INV_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1InverseWarp.nii.gz') + \
    ' -t ' + str([opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_0GenericAffine.mat'), 1])

    mvt_shft_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1Warp.nii.gz') + \
    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_0GenericAffine.mat')

    ## test on mean img (to see spatially that is works)
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.nii.gz') + \
    mvt_shft_ANTs
    spco([command], shell=True)

    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])
        for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro2,'tmp')):
            os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro2,'tmp'), f))

        ## apply on pre-processed imgs
        command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','dummy.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')
        spco([command], shell=True)

        #####transfo originally linearXXXX (line 157) should be linear too here????
        WW = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro2,'tmp','*.nii.gz')))
        for j in range(0, len(WW)):
            if j < 10:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.00' + str(j) + '.nii.gz') + \
                mvt_shft_ANTs
            elif j < 100 and j >= 10:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.0' + str(j) + '.nii.gz') + \
                mvt_shft_ANTs
            else:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.' + str(j) + '.nii.gz') + \
                mvt_shft_ANTs
            
            spco([command], shell=True)

        list_name = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.*.nii.gz')))
        name = ''
        for m in range(0,len(list_name)):
            name = name + ' ' + list_name[m]

        command = '3dTcat' + overwrite + ' ' + name + ' -tr ' + str(TR) + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz')
        spco([command], shell=True)