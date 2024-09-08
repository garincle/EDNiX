import os
import subprocess
from nilearn import plotting
import ants
import anatomical.Skullstrip_method
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def brainT_to_T(BASE_SS, creat_study_template, dir_prepro, ID, Session, listTimage, n_for_ANTS, dir_transfo, type_norm, BASE_SS_coregistr, Ref_file, volumes_dir, transfo_concat_inv,w2inv_inv, study_template_atlas_forlder, otheranat, bids_dir, type_of_transform, template_skullstrip, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, BASE_SS_mask, check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif, overwrite):

    #### if you want to use a sty template let's skullstrip this to improve the co-registration
    if creat_study_template== True:
        warp_adj = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'warped_3_adjusted_mean.nii.gz')
        stdy_template = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template.nii.gz')

        step_skullstrip = 3
        stdy_template_mask = anatomical.Skullstrip_method.Skullstrip_method(step_skullstrip, template_skullstrip, study_template_atlas_forlder, masking_img, brain_skullstrip_1, brain_skullstrip_2, masks_dir, volumes_dir, dir_prepro, type_norm, dir_transfo, BASE_SS_coregistr, BASE_SS_mask,
        otheranat, ID, Session, check_visualy_final_mask, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif)

        ##extract brain
        command = 'singularity run' + s_bind + afni_sif + '3dcalc -overwrite -a ' + stdy_template_mask + ' -b ' + warp_adj + ' -expr "(a*b)" -prefix ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz')
        spco(command, shell=True)

        ##align template to BASE_SS (atlas template)
        command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -warp shift_rotate -cmass -overwrite -base ' + BASE_SS + \
        ' -nomask -onepass' + \
        ' -master ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz') + \
        ' -prefix ' + stdy_template + \
        ' -source ' + opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm, 'study_template_not_align.nii.gz') + ' -1Dmatrix_save ' + \
        opj(study_template_atlas_forlder,ID + '_allign_template_to_stdy_template.1D')
        spco(command, shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dAllineate -overwrite -interp NN -1Dmatrix_apply ' + opj(study_template_atlas_forlder,ID + '_allign_template_to_stdy_template.1D') + \
        ' -prefix ' + stdy_template_mask + \
        ' -master ' + stdy_template + \
        ' -input  ' + stdy_template_mask
        spco([command], shell=True)

    #################################################################### coregistration with ANTs ####################################################################
    ######Coregistration!!!! (template to anat)
    IMG = ants.image_read(Ref_file)
    REF = ants.image_read(BASE_SS_coregistr)

    print(bcolors.OKGREEN + 'INFO: Co-registration ready' + bcolors.ENDC)
    print(bcolors.OKGREEN + 'INFO: WE WILL COREGISTER: ' + Ref_file + '(the anat img)' + ' to ' + BASE_SS_coregistr + '(the template)' + bcolors.ENDC)
    print(bcolors.OKGREEN + 'INFO: type_of_transform=' + type_of_transform + bcolors.ENDC)

    mTx  = ants.registration(fixed=IMG,moving=REF,
                              type_of_transform=type_of_transform,
                              outprefix=opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_'),
                              grad_step=0.1,
                              flow_sigma=3,
                              total_sigma=0,
                              aff_sampling=32,
                              aff_random_sampling_rate=0.2,
                              syn_sampling=32,
                              aff_iterations=(1000, 500, 250, 100),
                              aff_shrink_factors=(8, 4, 2, 1),
                              aff_smoothing_sigmas=(3, 2, 1, 0),
                              reg_iterations=(1000, 500, 250,100),
                              reg_smoothing_sigmas=(3, 2, 1, 0),
                              reg_shrink_factors=(8, 4, 2, 1),
                              verbose=True)

    TRANS = ants.apply_transforms(fixed=IMG, moving=REF,
                                      transformlist=mTx['fwdtransforms'], interpolator=n_for_ANTS)
    ants.image_write(TRANS,opj(dir_prepro,'template_to_' + type_norm + '_SyN_final.nii.gz'), ri=False)

    ######### apply the transformation to the original images to see the reverse transfo paramters works
    for Timage in listTimage:
        img_ref = opj(volumes_dir,ID + Timage + '_brain_step_1.nii.gz')
        IMG = ants.image_read(img_ref)
        TRANS = ants.apply_transforms(fixed=REF, moving=IMG,
                                      transformlist=transfo_concat_inv,
                                      interpolator=n_for_ANTS,whichtoinvert=w2inv_inv)
        ants.image_write(TRANS, opj(dir_prepro, Timage + '_to_template_SyN_final.nii.gz'), ri=False)

        ####plot the QC
        try:
            display = plotting.plot_anat(BASE_SS_coregistr, threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.add_contours(opj(dir_prepro, Timage + '_to_template_SyN_final.nii.gz'),
                                 linewidths=.2, colors=['red'])
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + 'to_template.png')
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(opj(dir_prepro,'template_to_' + Timage + '_SyN_final_test_matrix.nii.gz'), threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.savefig(bids_dir + '/QC/' + ID + '_' + str(Session) + '_' + Timage + 'to_template.png')
            # Don't forget to close the display
            display.close()