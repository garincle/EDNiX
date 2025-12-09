#################################################
########    create brain image of animal  ########
#################################################
#co-register linear to the new sty template...
##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################

import os
import json
import ants

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from Tools import QC_plot
import anat.skullstrip.Skullstrip_method

def create_indiv_template_brain(dir_prepro, list_transfo, ID, Session, listTimage, volumes_dir,
                                masking_img, brain_skullstrip_2,masks_dir, BASE_SS_coregistr, BASE_SS_mask,
                                check_visualy_final_mask,bids_dir,MNIBcorrect_indiv,sing_afni,sing_fsl,sing_fs, sing_itk,
                                sing_synstrip,Unetpath,diary_file,preftool):
    dir_prepro, list_transfo, ID,
    Session, listTimage, volumes_dir, masking_img,
    brain_skullstrip_2, masks_dir, BASE_SS_coregistr,
    BASE_SS_mask, check_visualy_final_mask,
    bids_dir, MNIBcorrect_indiv, sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip,
    diary_file, preftool


    nl = 'Run anat._5_create_template_brain.create_indiv_template_brain'
    run_cmd.msg(nl, diary_file,'HEADER')

    refnb = 0
    for i, j in enumerate(list_transfo):
        if list_transfo[i]["name"] == 'SS2':
            refnb = i
    nl = list_transfo[refnb]["type_of_transform"]
    run_cmd.msg(nl, diary_file, 'ENDC')

    end_maskname   = '_'.join([ID, 'final', 'mask', '2.nii.gz'])
    output4mask    = opj(masks_dir, ID + '_desc-step2_mask.nii.gz')

    input4msk      = opj(volumes_dir, ID + '_space-acpc_desc-SS-step1_' + masking_img + '.nii.gz')
    input_template = opj(volumes_dir, ID + '_space-acpc_desc-template_')
    if not MNIBcorrect_indiv == 'N3':
        N4_img         = opj(dir_prepro,  ID + '_space-acpc_desc-n3Bias_')
    else:
        N4_img         = opj(dir_prepro,  ID + '_space-acpc_desc-n4Bias_')
    SS_img         = opj(volumes_dir, ID + '_space-acpc_desc-SS-step2_')



    anat.skullstrip.Skullstrip_method.Skullstrip_method(brain_skullstrip_2, end_maskname, input4msk, output4mask,
                                                        masking_img, dir_prepro, masks_dir,
                                                        BASE_SS_coregistr, BASE_SS_mask,
                                                        list_transfo[refnb]["type_of_transform"], ID, list_transfo[refnb]["affmetric"],
                                                        check_visualy_final_mask,
                                                        sing_afni, sing_fsl, sing_fs, sing_itk,
                                                        sing_synstrip, Unetpath,
                                                        diary_file, preftool)

    nl = 'Run anat._5_create_template_brain.create_indiv_template_brain (section after skulls stripping)'
    run_cmd.msg(nl, diary_file,'HEADER')

    command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output4mask +
               ' -input ' + output4mask + ' -fill_holes')
    run_cmd.run(command, diary_file)

    #################################### signal correction
    # if T1w and T2w

    for Timage in listTimage:

        BRAIN = ants.image_read(input_template + Timage  + '.nii.gz')
        MSK = ants.image_read(output4mask)
        N4 = ants.n4_bias_field_correction(BRAIN, mask=MSK,
                                           shrink_factor=4,
                                           convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07},
                                           spline_param=200)
        ants.image_write(N4, N4_img + Timage + '.nii.gz', ri=False)

        dictionary = {"Sources": [input_template + Timage  + '.nii.gz',
                                  output4mask],
                      "Description": 'Non uniformity bias field correction (N4 from ANTSpy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(N4_img + Timage +'.json', "w") as outfile:
            outfile.write(json_object)

        # Check if all the data is zero
        if N4.sum() == 0:
            nl = "WARNING: N4BiasFieldCorrection failed, we can continue it is not the end of the world =)"
            run_cmd.msg(nl, diary_file, 'WARNING')

            command = (sing_afni + '3dcalc -overwrite -a ' + input_template + Timage  + '.nii.gz' +
                       ' -expr "a" -prefix ' + N4_img + Timage + '.nii.gz')
            run_cmd.do(command, diary_file)

            dictionary = {"Sources": input_template + Timage  + '.nii.gz',
                          "Description": 'Copy.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(N4_img + Timage +'.json', "w") as outfile:
                outfile.write(json_object)

        ### masking !!!!
        ### create brain of the subject template

        IMG = N4*MSK

        ants.image_write(IMG, SS_img + Timage + '.nii.gz', ri=False)

        dictionary = {"Sources": [N4_img + Timage + '.nii.gz',
                                  output4mask],
                      "Description": 'Skull strip.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(SS_img + Timage + '.json', "w") as outfile:
            outfile.write(json_object)

        ## plot several QC of the first skullstripped brain (check  that it is ok)

        QC_plot.mosaic(N4_img + Timage + '.nii.gz',
                       SS_img + Timage + '.nii.gz',
                       opj(bids_dir, 'QC','skullstrip_step2', ID + '_' + str(Session) + '_' + Timage + '_brain.png'))








