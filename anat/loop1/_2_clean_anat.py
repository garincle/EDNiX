import os
import ants
import json
import nibabel as nib
opj = os.path.join
opb = os.path.basename
opi = os.path.isfile
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from Tools import QC_plot
from Tools import getpath
from Tools import check_nii
from Tools import run_cmd, get_orientation

from anat.loop1 import acpcalign
from anat.loop2 import _mask2_Data_QC
from anat import norm2template
import anat.skullstrip.Skullstrip_method

def clean_anat(Align_img_to_template, bids_dir, listTimage, list_transfo, ID, Session,
               type_norm, data_path,masking_img, brain_skullstrip_1, BASE_SS_coregistr,
               BASE_SS_mask, BASE_SS, IgotbothT1T2, check_visualy_each_img, check_visualy_final_mask,overwrite, anat_ref_path,
               sing_afni,sing_fsl,sing_fs, sing_itk, sing_wb,sing_synstrip,Unetpath, Skip_step, diary_file,preftool):

    _,  dir_transfo,_, dir_prepro, _, volumes_dir, _, masks_dir = getpath.anat(data_path,
                                                                     '', '', False, False, 'native')

    _,  dir_transfo_ref,_, dir_prepro_ref, _, volumes_dir_ref, _, masks_dir_ref = getpath.anat(anat_ref_path[2],
                                                                     '', '', False, False, 'native')

    nl = 'Run anat._2_clean_anat.clean_anat'
    run_cmd.msg(nl, diary_file,'HEADER')

    anat_input0 = opj(dir_prepro, ID + '_space-raw_desc-reorient_')
    anat_input1 = opj(dir_prepro, ID + '_space-raw_desc-n4Bias_')

    anat_input2 = opj(dir_prepro,  ID + '_brain4AlignCenter_')
    anat_input4 = opj(dir_prepro,  ID + '_space-acpc_desc-64_')

    anat_input6 = opj(dir_prepro,  ID + '_space-acpc_desc-cropped_')
    anat_input7 = opj(volumes_dir, ID + '_space-acpc_desc-template_')
    anat_input8 = opj(volumes_dir, ID + '_space-acpc_desc-SS-step1_')

    end_maskname = '_'.join([ID, 'final', 'mask.nii.gz'])
    output4mask  = opj(masks_dir, ID + '_desc-step1_mask.nii.gz')
    msk_img      = opj(masks_dir, ID + '_space-acpc_desc-step1_mask.nii.gz')


    BASE_SS_anat_res = opj(bids_dir, 'BASE_SS_anat_res.nii.gz')
    # Load the image directly
    img = nib.load(opj(dir_prepro_ref, anat_ref_path[0] + '_space-raw_desc-reorient_' + type_norm + '.nii.gz'))
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]
    # Usage
    orient_meanimg = get_orientation.get_orientation_nibabel(opj(dir_prepro_ref, anat_ref_path[0] + '_space-raw_desc-reorient_' + type_norm + '.nii.gz'))

    if os.path.exists(BASE_SS_anat_res):
        img_BASE_SS_anat_res = nib.load(BASE_SS_anat_res)
        nl = "We found already a template supposed to be in anat resolution, let's check that the resolution of your reference anat image (img 0 in the list of anat) has not change!!"
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        res_img = img.header.get_zooms()[:3]
        res_base = img_BASE_SS_anat_res.header.get_zooms()[:3]

        if all(round(abs(a - b), 10) == 0 for a, b in zip(res_img, res_base)):
            nl = "OK, Same resolution ! (" + str(res_img) + ")"
            run_cmd.msg(nl, diary_file, 'OKGREEN')
        else:
            nl = "WARNING !!!!!!  Different resolutions:"
            run_cmd.msg(nl, diary_file, 'WARNING')
            nl = " Image:" + str(res_img)
            run_cmd.msg(nl, diary_file, 'WARNING')
            nl = " BASE_SS_anat_res:"+ str(res_base)
            run_cmd.msg(nl, diary_file, 'WARNING')

            command = (sing_afni + '3dresample' + overwrite +
                       ' -orient ' + orient_meanimg +
                       ' -prefix ' + BASE_SS_anat_res +
                       ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                       ' -rmode Cu -input ' + BASE_SS)
            run_cmd.run(command, diary_file)

    else:
        print('No BASE_SS_anat_res found, we create it now')
        command = (sing_afni + '3dresample' + overwrite +
                   ' -orient ' + orient_meanimg +
                   ' -prefix ' + BASE_SS_anat_res +
                   ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                   ' -rmode Cu -input ' + BASE_SS)
        run_cmd.run(command, diary_file)

    refnb1 = 0
    for i, j in enumerate(list_transfo):
        if list_transfo[i]["name"] == 'SS1':
            refnb1 = i

    anat.skullstrip.Skullstrip_method.Skullstrip_method(brain_skullstrip_1, end_maskname,
                                                        anat_input1 + masking_img + '.nii.gz',
                                                        output4mask, masking_img, dir_prepro, masks_dir,
                                                        BASE_SS_coregistr, BASE_SS_mask, list_transfo[refnb1]["type_of_transform"], ID, list_transfo[refnb1]["affmetric"],
                                                        check_visualy_final_mask, sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip, Unetpath, diary_file, preftool)

    nl = 'Run anat._2_clean_anat.clean_anat (section after skulls stripping)'
    run_cmd.msg(nl, diary_file,'HEADER')
    
    #### Apply masking to other anat imagesother

    command = (sing_afni + '3dcalc -overwrite -a ' + anat_input1 + type_norm + '.nii.gz' + ' -b ' + output4mask +
               ' -expr "a*b" -prefix ' + anat_input2 + type_norm + '.nii.gz')
    run_cmd.do(command, diary_file)

    dictionary = {"Sources": [anat_input1 + type_norm + '.nii.gz',
                              output4mask],
                  "Description": 'Skull stripping.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anat_input2 + type_norm + '.json', "w") as outfile:
        outfile.write(json_object)

    if 'itk_1' in Skip_step:
        run_cmd.msg('INFO: skip step ' + str('itk_1'), diary_file, 'OKGREEN')
    else:
        output4mask = output4mask
        end_maskname = '_'.join([ID, 'final', 'mask', '.nii.gz'])
        input4msk = anat_input1 + type_norm + '.nii.gz'
        _mask2_Data_QC._itk_check_masks(output4mask, input4msk, end_maskname, masks_dir, sing_itk, diary_file)

    #  QC
    QC_plot.mosaic(anat_input1 + type_norm + '.nii.gz',
                   anat_input2 + type_norm + '.nii.gz',
                   opj(bids_dir, 'QC','skullstrip_step1', '_'.join([ID, str(Session), type_norm, 'skullstriped.png'])))


    ####################################################################################
    ########################## transfo rigid to atlas template (BASE_SS)   #############
    ####################################################################################
    refnb2 = 0

    for i, j in enumerate(list_transfo):
        if list_transfo[i]["name"] == 'align':
            refnb2 = i

    if Align_img_to_template == 'No':
        acpcalign.afni_empty(dir_transfo,diary_file)

    elif Align_img_to_template == 'Ants':
        nl = 'acpc alignement'
        run_cmd.msg(nl, diary_file,'OKGREEN')

        norm2template.norm(ID, anat_input2 + type_norm + '.nii.gz', output4mask, dir_prepro,
                           'acpc', BASE_SS_anat_res, '',
                           dir_transfo, list_transfo[refnb2]["type_of_transform"], 'Translation', 'acpc',
                           list_transfo[refnb2]["affmetricT"], list_transfo[refnb2]["affmetric"], list_transfo[refnb2]["interpol"],
                           diary_file, sing_wb, '_desc-64_' + type_norm, 0)

    else:
        nl = 'ERROR: Align_img_to_template need to be define as string ("No" or "Ants")'
        raise Exception(run_cmd.error(nl, diary_file))


    command = (sing_afni + '3dAutobox -input ' + anat_input4 + type_norm + '.nii.gz' +
               ' -npad 2 -prefix ' + anat_input6 + type_norm + '.nii.gz' + ' -noclust -overwrite')
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": anat_input4 + type_norm + '.nii.gz',
                  "Description": 'Crop.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anat_input6 + type_norm + '.json', "w") as outfile:
        outfile.write(json_object)


    if check_visualy_each_img == True:
        ####### optional manual acpc center #######
        # Check the result and do manually the finest correction save the file as "ID_acpc_tmp.nii.gz"
        command = (sing_fs + 'freeview -viewport x -v ' + BASE_SS + ' ' + anat_input6 + type_norm + '.nii.gz' +
                   ':opacity=0.6:visible=1 -ras 0 0 0 -cc -transform-volume -subtitle "SET-THE-VOLUME-TO-THE-ACPC-PLANE!"')
        run_cmd.do(command, diary_file)

        if opi(opj(dir_prepro, ID + '_acpc_tmp.nii.gz')):
            norm2template.norm(ID, anat_input2 + type_norm + '.nii.gz', '', dir_prepro,
                               'acpc', opj(dir_prepro, ID + '_acpc_tmp.nii.gz'), '',
                               dir_transfo, list_transfo[refnb2]["type_of_transform"], 'shift', 'acpc',
                               list_transfo[refnb2]["affmetricT"], list_transfo[refnb2]["affmetric"], list_transfo[refnb2]["interpol"],
                               diary_file, sing_wb, '_desc-cropped_' + type_norm, 0)


    ##### Apply transfo to other anat images
    IMG = ants.image_read(output4mask)
    REF = ants.image_read(anat_input6 + type_norm + '.nii.gz')
    msk = ants.apply_transforms(fixed=REF, moving=IMG,
                                    transformlist=opj(dir_transfo, 'acpc_0GenericAffine.mat'),
                                    interpolator='genericLabel', whichtoinvert=[False])
    ants.image_write(msk, msk_img)

    dictionary = {"Sources": [output4mask,
                              anat_input6 + type_norm + '.nii.gz'],
                  "Description": 'Rigid normalization.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(msk_img.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    if IgotbothT1T2 == True:
        listTimage2 = list(listTimage)
        listTimage2.append('T1wdividedbyT2w')
    else:
        listTimage2 = list(listTimage)

    for Timage in listTimage2:
        if Timage == 'T1wdividedbyT2w':
            input_Allin = anat_input0 + Timage + '.nii.gz'
        else:
            input_Allin = anat_input1 + Timage + '.nii.gz'

        if not Align_img_to_template == 'No':
            IMG = ants.image_read(input_Allin)
            REF = ants.image_read(anat_input6 + type_norm + '.nii.gz')
            new_img = ants.apply_transforms(fixed=REF, moving=IMG,
                                            transformlist=opj(dir_transfo, 'acpc_0GenericAffine.mat'),
                                            interpolator=list_transfo[refnb1]["interpol"], whichtoinvert=[False])
            ants.image_write(new_img, anat_input6 + Timage + '.nii.gz')
            check_nii.resamp(anat_input6 + Timage + '.nii.gz',anat_input6 + type_norm + '.nii.gz','T1w','','',diary_file, sing_wb)


            dictionary = {"Sources": [input_Allin,
                                      anat_input6 + type_norm + '.nii.gz'],
                          "Description": 'Rigid normalization.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(anat_input6 + Timage + '.json', "w") as outfile:
                outfile.write(json_object)

        else:
            #### ensure that both image are now in the same space, with the same header
            check_nii.resamp(msk_img, anat_input6 + type_norm + '.nii.gz', 'msk', '', '',
                             diary_file,
                             '')
        ######################################################################
        ###     B0 correction ####### ==> create INDIV template !!!!!      ###
        ######################################################################

        IMG = ants.image_read(anat_input6 + Timage + '.nii.gz')
        IMG = ants.denoise_image(IMG, r=3, noise_model='Gaussian')

        ants.image_write(IMG, anat_input7 + Timage + '.nii.gz', ri=False)

        dictionary = {"Sources": anat_input6 + Timage + '.nii.gz',
                      "Description": 'Denoising (ANTspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anat_input7 + Timage + '.json', "w") as outfile:
            outfile.write(json_object)


        ######################################################################
        ###     Skull stripping  !!                                        ###
        ######################################################################
        new_img = IMG * msk
        ants.image_write(new_img, anat_input8 + Timage + '.nii.gz')

        dictionary = {"Sources": [anat_input7 + Timage + '.nii.gz',
                                  msk_img],
                      "Description": 'Skull stripping.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anat_input8 + Timage + '.json', "w") as outfile:
            outfile.write(json_object)

        ######################################################################
        ###     QC               !!                                        ###
        ######################################################################
        QC_plot.mosaic(anat_input8 + Timage + '.nii.gz',
                       BASE_SS_coregistr,
                       opj(bids_dir, 'QC','align_rigid_to_template',
                           ID + '_' + str(Session) + '_' + Timage + '_align_rigid_to_template.png'))


        ####### check visually indiv template #######
        if check_visualy_each_img == True:
            # now look at the coordinate of the brain and ajust bet2 according to them
            command = (sing_fs + 'freeview -v ' + anat_input8 + Timage + '.nii.gz')
            run_cmd.do(command, diary_file)



