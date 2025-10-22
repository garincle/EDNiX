import os
import numpy as np
import ants
import json
import nibabel as nib
from nilearn.image import resample_to_img

opj = os.path.join
opb = os.path.basename
opi = os.path.isfile
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from Tools import QC_plot
from Tools import getpath

from anatomical.loop1 import acpcalign
import anatomical.skullstrip.Skullstrip_method

def clean_anat(Align_img_to_template, cost3dAllineate, bids_dir, listTimage, type_of_transform, ID, aff_metric_ants, Session, otheranat,
               type_norm, data_path,masking_img, brain_skullstrip_1, BASE_SS_coregistr,
               BASE_SS_mask, BASE_SS, IgotbothT1T2, check_visualy_each_img, check_visualy_final_mask,overwrite,n_for_ANTS,
               sing_afni,sing_fsl,sing_fs, sing_itk, sing_synstrip,Unetpath,diary_file,preftool):

    _, _, dir_transfo, dir_prepro, _, volumes_dir, _, masks_dir = getpath.anat(data_path,
                                                                     '', '', False, False, 'native')

    nl = 'Run anatomical._2_clean_anat.clean_anat'
    run_cmd.msg(nl, diary_file,'HEADER')
    anat_input0 = opj(dir_prepro, ID + '_space-raw_desc-reorient_')
    anat_input1 = opj(dir_prepro, ID + '_space-raw_desc-n4Bias_')
    anat_input2 = opj(dir_prepro,  ID + '_brain4AlignCenter_')
    anat_input3 = opj(dir_prepro,  ID + '_space-acpc_desc-64-orig-3dAllineate_')
    anat_input4 = opj(dir_prepro,  ID + '_space-acpc_desc-64_')

    anat_input5 = opj(dir_prepro,  ID + '_space-acpc_test_QC_')
    anat_input6 = opj(dir_prepro,  ID + '_space-acpc_desc-cropped_')
    anat_input7 = opj(volumes_dir, ID + '_space-acpc_desc-template_')

    end_maskname = '_'.join([ID, 'final', 'mask.nii.gz'])
    output4mask  = opj(masks_dir, ID + '_desc-step1_mask.nii.gz')
    maskRSTnorm  = opj(masks_dir, ID + '_desc-step1+norm_mask.nii.gz')

    anatomical.skullstrip.Skullstrip_method.Skullstrip_method(brain_skullstrip_1, end_maskname,
                                                              anat_input1 + masking_img + '.nii.gz',
                                                              output4mask,masking_img, dir_prepro,masks_dir,
                                                              BASE_SS_coregistr, BASE_SS_mask,type_of_transform, ID, aff_metric_ants,
                                                              check_visualy_final_mask, sing_afni, sing_fsl, sing_fs, sing_itk,
                                                              sing_synstrip, Unetpath, diary_file,preftool)

    nl = 'Run anatomical._2_clean_anat.clean_anat (section after skulls stripping)'
    run_cmd.msg(nl, diary_file,'HEADER')
    
    #### Apply masking to other anat imagesother
    for Timage in listTimage:

        command = (sing_afni + '3dcalc -overwrite -a ' + anat_input1 + Timage + '.nii.gz' + ' -b ' + output4mask +
                   ' -expr "a*b" -prefix ' + anat_input2 + Timage + '.nii.gz')
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": [anat_input1 + Timage + '.nii.gz',
                                  output4mask],
                      "Description": 'Skull stripping.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anat_input2 + Timage + '.json', "w") as outfile:
            outfile.write(json_object)

        #  QC
        QC_plot.mosaic(anat_input1 + Timage + '.nii.gz',
                       anat_input2 + Timage + '.nii.gz',
                       opj(bids_dir, 'QC','skullstrip_step1', '_'.join([ID, str(Session), Timage, 'skullstriped.png'])))


        if not Align_img_to_template == 'Ants':
            # Normalization of the images size (Brain)
            command = (sing_afni + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' +
                       anat_input2 + Timage + '.nii.gz' + ' ' + anat_input2 + Timage + '.nii.gz' + ' -overwrite')
            run_cmd.run(command, diary_file)

    if not Align_img_to_template == 'Ants':
        command = (sing_afni + '3dresample -master ' + anat_input2 + type_norm + '.nii.gz' +
                   ' -prefix ' + maskRSTnorm +
                   ' -input ' +  output4mask + ' -overwrite')
        run_cmd.run(command, diary_file)

    ####################################################################################
    ########################## transfo rigid to atlas template (BASE_SS)   #############
    ####################################################################################


    if Align_img_to_template == '3dAllineate':

        acpcalign.allineate(data_path,
                            anat_input2 + type_norm + '.nii.gz',
                            anat_input3 + type_norm + '.nii.gz',
                            maskRSTnorm,
                            BASE_SS,
                            cost3dAllineate,sing_afni,overwrite,diary_file)

    elif Align_img_to_template == '@Align_Centers':

        acpcalign.aligncenter(data_path,
                            anat_input2 + type_norm + '.nii.gz',
                            anat_input3 + type_norm + '.nii.gz',
                            BASE_SS,
                            sing_afni, overwrite, diary_file)

    elif Align_img_to_template == 'No':


        acpcalign.none(anat_input2 + type_norm + '.nii.gz',
                       anat_input3 + type_norm + '.nii.gz',
                       overwrite,diary_file,sing_afni)

    elif Align_img_to_template == 'Ants':
        acpcalign.ants(data_path, anat_input2 + type_norm + '.nii.gz',
                       anat_input4 + type_norm + '.nii.gz',
                       output4mask,
                       BASE_SS, n_for_ANTS, diary_file)

    else:
        nl = 'ERROR: Align_img_to_template need to be define as string (3dAllineate, @Align_Centers, No)'
        raise Exception(run_cmd.error(nl, diary_file))

    if not Align_img_to_template == 'Ants':
        # Get voxel sizes
        img = nib.load(anat_input2 + type_norm + '.nii.gz')
        delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]

        command = (sing_afni + '3dresample' + overwrite +
                   ' -prefix ' + anat_input3 + type_norm + '.nii.gz' +
                   ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' +
                   ' -input ' + anat_input3 + type_norm + '.nii.gz')
        run_cmd.run(command, diary_file)

    if IgotbothT1T2 == True:
        listTimage2 = list(listTimage)
        listTimage2.append('T1wdividedbyT2w')
    else:
        listTimage2 = list(listTimage)

    ##### Apply transfo to other anat images

    for Timage in listTimage2:

        if Timage == 'T1wdividedbyT2w':
            input_Allin = anat_input0 + Timage + '.nii.gz'
        else:
            input_Allin = anat_input2 + Timage + '.nii.gz'

        if not Align_img_to_template == 'ants':

            command = (sing_afni + '3dAllineate -overwrite -interp NN -1Dmatrix_apply ' + anat_input2 + type_norm + '.1D' +
                       ' -prefix ' + anat_input4 + Timage + '.nii.gz' +
                       ' -master ' + anat_input3 + type_norm + '.nii.gz' +
                       ' -input ' + input_Allin)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": [input_Allin,
                                      anat_input3 + type_norm + '.nii.gz'],
                          "Description": 'Rigid normalization.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(anat_input3 + type_norm + '.json', "w") as outfile:
                outfile.write(json_object)

            if Timage in [str(type_norm), str(otheranat)]:
                    command = (sing_afni + '3dAllineate -overwrite -interp NN -1Dmatrix_apply ' + anat_input2 + type_norm + '.1D' +
                               ' -prefix ' + anat_input5 + Timage + '.nii.gz' +
                               ' -master ' + anat_input4 + Timage + '.nii.gz' +
                               ' -input '  + anat_input0 + Timage + '.nii.gz')
                    run_cmd.run(command, diary_file)

                    dictionary = {"Sources": [anat_input0 + Timage + '.nii.gz',
                                              anat_input4 + type_norm + '.nii.gz'],
                                  "Description": 'Rigid normalization.', }
                    json_object = json.dumps(dictionary, indent=2)
                    with open(anat_input5 + Timage + '.json', "w") as outfile:
                        outfile.write(json_object)

                    command = (sing_afni + '3dAutobox -input ' + anat_input4 + Timage + '.nii.gz' +
                               ' -prefix ' + anat_input5 + Timage + '.nii.gz' + ' -noclust -overwrite')
                    run_cmd.run(command, diary_file)

            command = (sing_afni + '3dAutobox' + overwrite + ' -input ' + anat_input4 + Timage + '.nii.gz' +
                       ' -prefix ' + anat_input6 + Timage + '.nii.gz' + ' -noclust -overwrite')
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": anat_input4 + Timage + '.nii.gz',
                          "Description": 'Crop.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(anat_input6 + Timage + '.json', "w") as outfile:
                outfile.write(json_object)

        else:
            print('apply_transforms')
            IMG = ants.image_read(input_Allin)
            REF = ants.image_read(anat_input4 + type_norm + '.nii.gz')
            new_img = ants.apply_transforms(fixed=REF, moving=IMG,
                                            transformlist=opj(dir_transfo, 'acpc_0GenericAffine.mat'),
                                            interpolator=n_for_ANTS, whichtoinvert=[False])
            ants.image_write(new_img, anat_input4 + Timage + '.nii.gz')

            dictionary = {"Sources": [input_Allin,
                                      anat_input4 + type_norm + '.nii.gz'],
                          "Description": 'Rigid normalization.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(anat_input4 + Timage + '.json', "w") as outfile:
                outfile.write(json_object)
            print('3dAutobox')
            command = (sing_afni + '3dAutobox -input ' + anat_input4 + Timage + '.nii.gz' +
                       ' -prefix ' + anat_input6 + Timage + '.nii.gz' + ' -noclust -overwrite')
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": anat_input4 + Timage + '.nii.gz',
                          "Description": 'Crop.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(anat_input6 + Timage + '.json', "w") as outfile:
                outfile.write(json_object)


            ####### optional crop !! save as "ID_acpc_cropped.nii.gz" #######





        ####### optional manual acpc center #######
        if check_visualy_each_img == True:

            # Check the result and do manually the finest correction save the file as "ID_acpc_tmp.nii.gz"
            command = (sing_fs + 'freeview -viewport x -v ' + BASE_SS + ' ' + anat_input6 + type_norm + '.nii.gz' +
                       ':opacity=0.6:visible=1 -ras 0 0 0 -cc -transform-volume -subtitle "SET-THE-VOLUME-TO-THE-ACPC-PLANE!"')
            run_cmd.do(command, diary_file)
            if opi(opj(dir_prepro, ID + '_acpc_tmp.nii.gz')):
                acpcalign.ants(data_path, anat_input2 + type_norm + '.nii.gz',
                               anat_input6 + type_norm + '.nii.gz','',
                               opj(dir_prepro, ID + '_acpc_tmp.nii.gz'), n_for_ANTS, diary_file)

                for Timage in [str(type_norm), str(otheranat)]:

                    IMG = ants.image_read(anat_input2 + Timage + '.nii.gz')
                    REF = ants.image_read(opj(dir_prepro, ID + '_acpc_tmp.nii.gz'))
                    new_img =  ants.apply_transforms(fixed=REF, moving=IMG,
                                                     transformlist=opj(dir_transfo, 'acpc_0GenericAffine.mat'),
                                                     interpolator=n_for_ANTS,whichtoinvert=[False])
                    ants.image_write(new_img,anat_input6 + Timage + '.nii.gz')






    ######################################################################
    ####### B0 correction ####### ==> create INDIV template !!!!!!!!!!!####
    ######################################################################

    for Timage in listTimage:

        IMG = ants.image_read(anat_input6 + Timage + '.nii.gz')
        IMG = ants.denoise_image(IMG,r=3,noise_model='Gaussian')
        ants.image_write(IMG, anat_input7 + Timage + '.nii.gz', ri=False)

        dictionary = {"Sources": anat_input6 + Timage + '.nii.gz',
                      "Description": 'Denoising (ANTspy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anat_input7 + Timage + '.json', "w") as outfile:
            outfile.write(json_object)

        # QC
        QC_plot.mosaic(anat_input7 + Timage + '.nii.gz',
                       BASE_SS_coregistr,
                       opj(bids_dir, 'QC','align_rigid_to_template',
                           ID + '_' + str(Session) + '_' + Timage + '_align_rigid_to_template.png'))


        ####### check visually indiv template #######
        if check_visualy_each_img == True:
            # now look at the coordinate of the brain and ajust bet2 according to them
            command = (sing_fs + 'freeview -v ' + anat_input7 + Timage + '.nii.gz')
            run_cmd.do(command, diary_file)






