##############################################################################
####### Co-register your study template to the atlas template ################
##############################################################################
#Path to the excels files and data structure

import os
import subprocess
import ants

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


def stdyT_to_AtlasT(list_atlases, Aseg_ref, Aseg_refLR, BASE_SS, dir_out, n_for_ANTS, study_template_atlas_forlder, Atemplate_to_Stemplate, overwrite,
                    s_bind,afni_sif):


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

    command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + ' -overwrite' + ' -warp shift_rotate -cmass -source ' + BASE_SS + \
    ' -nomask' + \
    ' -prefix ' + opj(dir_out, 'template_in_stdy_template.nii.gz') + \
    ' -base ' + stdy_template + ' -1Dmatrix_save ' + opj(dir_out,'Align_Center_shft.1D')
    spco(command, shell=True)

    ######Coregistration!!!!

    REF = ants.image_read(stdy_template)
    IMG = ants.image_read(opj(dir_out, 'template_in_stdy_template.nii.gz'))

    mTx = ants.registration(fixed=REF, moving=IMG,
                            type_of_transform='SyNCC',
                            outprefix=opj(dir_out,'NMT_to_anat_SyN_final_'))
    TRANS = ants.apply_transforms(fixed=REF, moving=IMG,
                                  transformlist=mTx['fwdtransforms'], interpolator=n_for_ANTS)
    ants.image_write(TRANS, opj(dir_out,'NMT_to_anat_SyN_final.nii.gz'), ri=False)


    ####################################################################################
    ########################## seg into the native space ###################

    for atlas in list_atlases2:

        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + opj(dir_out,'Align_Center_shft.1D') + \
        ' -prefix ' + opj(dir_out, opb(ops(ops(atlas)[0])[0]) + '_Allin.nii.gz') + \
        ' -master ' + stdy_template + \
        ' -input  ' + atlas
        spco([command], shell=True)

        IMG = ants.image_read(opj(dir_out, opb(ops(ops(atlas)[0])[0]) + '_Allin.nii.gz'))

        TRANS = ants.apply_transforms(fixed=REF, moving=IMG,
                                      transformlist=mTx['fwdtransforms'], interpolator='nearestNeighbor')

        ants.image_write(TRANS, opj(dir_out, opb(atlas)), ri=False)



