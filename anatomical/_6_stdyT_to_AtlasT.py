##############################################################################
####### Co-register your study template to the atlas template ################
##############################################################################
#Path to the excels files and data structure

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
ops = os.path.splitext

#################################################
########    creat brain image of animal  ########
#################################################


def stdyT_to_AtlasT(list_atlases, Aseg_ref, Aseg_refLR, BASE_SS, dir_out, n_for_ANTS, study_template_atlas_forlder, Atemplate_to_Stemplate, overwrite):


    stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + Atemplate_to_Stemplate, 'study_template.nii.gz')

    #######################################################################
    ################coregistration temaplte to template####################
    #######################################################################

    if not os.path.exists(dir_out):
        os.mkdir(dir_out)


    def flatten_comprehension(matrix):
        return [item for row in matrix for item in row]

    ####################################################################################
    list_atlases2 = []
    list_atlases2.extend(list_atlases)
    list_atlases2.extend([Aseg_ref])
    list_atlases2.extend([Aseg_refLR])

    command = '3dAllineate' + ' -overwrite' + ' -warp shift_rotate -cmass -source ' + BASE_SS + \
    ' -nomask' + \
    ' -prefix ' + opj(dir_out, 'template_in_stdy_template.nii.gz') + \
    ' -base ' + stdy_template + ' -1Dmatrix_save ' + opj(dir_out,'Align_Center_shft.1D')
    spco(command, shell=True)


    ######Coregistration!!!!
    command = 'antsRegistration -d 3 --float 0 --verbose 1 -n ' + n_for_ANTS +\
        ' -o [' + opj(dir_out,'NMT_to_anat_SyN_final_') + ',' + opj(dir_out,'NMT_to_anat_SyN_final.nii.gz') + ']' + \
        ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m MI[' + stdy_template + ',' + opj(dir_out, 'template_in_stdy_template.nii.gz') + ',1,32,Regular,0.2]' + \
        ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
        ' -m CC[' + stdy_template + ',' + opj(dir_out, 'template_in_stdy_template.nii.gz') + ',1,4,Regular,0.2]'
    spco([command], shell=True)


    #command = '3dQwarp -overwrite -iwarp' + ' -base ' + BASE_SS + ' -prefix ' + opj(dir_out,'NMT_to_anat_SyN_final.nii.gz') + ' -source ' + opj('/media/cgarin/Clement_1/1_Macaques/1_PFC_study/10_2023_PP_U/Study_template/studytemplate2_T2/template2_Al.nii.gz')
    #spco(command, shell=True)

    ####################################################################################
    ########################## seg into the native space ###################

    for atlas in list_atlases2:

        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + opj(dir_out,'Align_Center_shft.1D') + \
        ' -prefix ' + opj(dir_out, opb(ops(ops(atlas)[0])[0]) + '_Allin.nii.gz') + \
        ' -master ' + stdy_template + \
        ' -input  ' + atlas
        spco([command], shell=True)


        command = 'antsApplyTransforms -d 3 -i ' + opj(dir_out, opb(ops(ops(atlas)[0])[0]) + '_Allin.nii.gz') + \
            ' -r ' + stdy_template + \
            ' -o ' + opj(dir_out, opb(atlas)) + \
            ' -t ' + opj(dir_out,'NMT_to_anat_SyN_final_1Warp.nii.gz') + \
            ' -t ' + opj(dir_out,'NMT_to_anat_SyN_final_0GenericAffine.mat') + \
            ' -n NearestNeighbor'
        spco([command], shell=True)

    '''
    for atlas in list_atlases:
        spco(['3dNwarpApply', '-nwarp', opj(dir_out,'NMT_to_anat_SyN_final_WARPINV.nii.gz'), '-interp', 'NN',
        '-source', atlas, '-master', opj('/media/cgarin/Clement_1/1_Macaques/1_PFC_study/10_2023_PP_U/Study_template/studytemplate2_T2/template2_Al.nii.gz'), '-prefix', opj(dir_out, opb(atlas)), '-overwrite'])


    for atlas in list_atlases:
        spco(['3dcalc', '-overwrite', '-a', 
            opj(dir_out, opb(atlas)), 
            '-b', stdy_template,
            '-prefix', opj(dir_out, opb(atlas)), '-expr', 'step(b)*a'])

    '''