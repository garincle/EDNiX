###################################
###      Skullstrip method      ###
###################################

import os
import shutil
import nibabel as nib
import ants
import nilearn
import json
import numpy as np

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd

from anat.skullstrip import afniskullstrip
from anat.skullstrip import epimask
from anat.skullstrip import sammbamask
from anat.skullstrip import antsmask
from anat.skullstrip import qwarpmask
from anat.skullstrip import noSstrip
from anat.skullstrip import betmask
from anat.skullstrip import manualmask
from anat.skullstrip import Unetmask
from anat.skullstrip import sstripmask
from Tools import check_nii
############################################
### a revoir le naming avec logique !!!! ###
############################################

SSmethodslist= ['3dSkullStrip','3dSkullStrip_noDilate','3dSkullStrip_macaque','3dSkullStrip_dog_macaque','3dSkullStrip_Rat','3dSkullStrip_rat','3dSkullStrip_rat_no_dil',
                'Custum_Macaque','Custum_Baboon','Custum_mouse','Custum_bat','Custum_dog','custum_rat','CustumNilearn_','CustumNilearnExcludeZeros_',
                'Custum_Chimp',
                'Vol_sammba_','sammba_rat','sammba_dog','sammba_mouse',
                '_bet',
                'bet2_ANTS','Custum_ANTS_NL','Custum_ANTS_Garin',
                'CustumThreshold_',
                'synthstrip',
                'muSkullStrip_Human','Custum_Macaque2',
                'QWARP','Custum_QWARP','Custum_QWARPT2','Custum_QWARPT2_dil',
                'NoSkullStrip','Manual']


def Skullstrip_method(brain_skullstrip, end_maskname, input_for_msk,output_for_mask,masking_img, dir_prepro,masks_dir,
                      BASE_SS_coregistr, BASE_SS_mask,type_of_transform, ID, aff_metric_ants, check_visualy_final_mask,
                      sing_afni, sing_fsl, sing_fs, sing_itk, sing_synstrip,Unetpath,diary_file,preftool):

    nl = 'Run anat.Skullstrip_method.Skullstrip_method'
    run_cmd.msg(nl, diary_file,'HEADER')

    nl = ('INFO: If you can not find a good solution for Skullstriping due to bad image quality, you can always modify it by hands and save it as: ' +
          opj(masks_dir,end_maskname))
    run_cmd.msg(nl, diary_file,'OKGREEN')

    nl = 'WARNING: Note that any manual segmentation (saved with the correct name) it will AUTOMATICALLY be selected for Skullstriping'
    run_cmd.msg(nl, diary_file,'WARNING')


    if brain_skullstrip == '':
        nl ='No step_skullstrip ?'
        run_cmd.msg(nl, diary_file,'ENDC')
        nl = 'Nothing will be done.'
        run_cmd.msg(nl, diary_file,'ENDC')
        raise Exception(run_cmd.error(nl + ': '+ nl, diary_file))


    nl = "INFO: brain_skullstrip method is " + brain_skullstrip
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    nl = 'INFO: looking for manual segmentation named:' + opj(masks_dir,end_maskname) + '...'
    run_cmd.msg(nl, diary_file,'OKGREEN')

    ###########################################################################################################################################################################
    ######################################################### Brain Skullstrip Library (do your own if you need to!!!) ########################################################
    ###########################################################################################################################################################################

    # This Brain Skullstrip Library will be evolving and can/should be personalized, don't hesitate to propose other efficient solution on GIT!!!

    if ope(opj(masks_dir,end_maskname)):
        nl = 'WARNING: We found an already existing final mask !!! No Skullstrip will be performed!'
        run_cmd.msg(nl, diary_file, 'WARNING')
        shutil.copyfile(opj(masks_dir,end_maskname), output_for_mask)

        #### ensure that both image are now in the same space, with the same header
        check_nii.resamp(output_for_mask, input_for_msk, 'msk', '', '',
                         diary_file,
                         '')
        dictionary = {"Sources": opj(masks_dir,end_maskname),
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(output_for_mask[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

        nl = 'INFO: please delete' + opj(masks_dir, end_maskname) + ' if you want retry to create a skulstripp image'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

    else:
        nl = 'WARNING: Manual segmentation not found'
        run_cmd.msg(nl, diary_file, 'WARNING')
        nl = 'INFO: continuing with the selected skullstriping method'
        run_cmd.msg(nl, diary_file, 'OKGREEN')


        if brain_skullstrip == '3dSkullStrip':

            afniskullstrip.run(input_for_msk, output_for_mask,' -blur_fwhm 2 -orig_vol -mask_vol -use_skull',
                           ' -fill_holes -dilate_input 2',
                           sing_afni, diary_file)

        elif brain_skullstrip == '3dSkullStrip_noDilate':
            afniskullstrip.run(input_for_msk, output_for_mask,' -blur_fwhm 2 -orig_vol -mask_vol -use_skull',
                           ' -fill_holes',
                           sing_afni, diary_file)

        #### species specific ######################################## MACAQUE ####################################
        elif brain_skullstrip == '3dSkullStrip_macaque':

            afniskullstrip.run(input_for_msk,output_for_mask,' -blur_fwhm 1 -orig_vol -mask_vol -monkey',
                           ' -fill_holes -dilate_input 4',
                           sing_afni,diary_file)

        elif brain_skullstrip == 'Custum_Macaque':

            epimask.run(input_for_msk, output_for_mask, [0.2, 0.90], False, ' -fill_holes -dilate_input 8 -2',
                        sing_afni, sing_fs, diary_file)

        #### species specific ######################################## MARMOSET ####################################

        elif brain_skullstrip == '3dSkullStrip_marmoset':
            afniskullstrip.run(input_for_msk, output_for_mask,' -blur_fwhm 1 -shrink_fac 0.1 -fac 100 -norm_vol -mask_vol -marmoset',
                           ' -fill_holes',
                           sing_afni, diary_file)

        elif brain_skullstrip == '3dSkullStrip_dog_macaque':

            afniskullstrip.run(input_for_msk, output_for_mask,' -blur_fwhm 2 -orig_vol -mask_vol -monkey',
                           ' -fill_holes -dilate_input 4 -4',
                           sing_afni, diary_file)

        #################################### BABOON ####################################

        elif brain_skullstrip == 'Custum_Baboon':

            epimask.run(input_for_msk,output_for_mask,[0.85,0.95],False,' -fill_holes -dilate_input 15',
                        sing_afni,sing_fs,diary_file)

        #################################### CHIMPANZEE ####################################

        elif brain_skullstrip == 'Custum_Chimp':

            loadimg       = nib.load(input_for_msk).get_fdata()
            loadimgsort85 = np.percentile(np.abs(loadimg)[np.abs(loadimg)>0], 40)
            mask_imag     = nilearn.image.threshold_img(input_for_msk, loadimgsort85, cluster_threshold=10)
            mask_imag.to_filename(output_for_mask)

            command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output_for_mask +
                       ' -input ' + output_for_mask + ' -fill_holes -dilate_input 3')
            run_cmd.run(command, diary_file)

            command = (sing_afni + '3dresample -master ' + input_for_msk +
                       ' -prefix ' + output_for_mask + ' -input ' + output_for_mask +
                       ' -overwrite -bound_type SLAB')
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": input_for_msk,
                          "Description": 'Brain mask (image threshold from Nilearn).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(output_for_mask[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        #################################### Mouse ####################################

        elif brain_skullstrip == 'Custum_mouse':

            epimask.run(input_for_msk, output_for_mask, [0.2, 0.60], False, ' -fill_holes -dilate_input 3 -1',
                        sing_afni, sing_fs, diary_file)

        elif brain_skullstrip == 'sammba_mouse':

            sammbamask.run(input_for_msk, output_for_mask, [360, 0.2, 0.8, 500, 2, 10], sing_fs, diary_file)

        #################################### BAT ####################################

        elif brain_skullstrip == 'Custum_bat':

            epimask.run(input_for_msk, output_for_mask, [0.7, 0.80], False, ' -fill_holes -dilate_input 8 -2',
                        sing_afni, sing_fs, diary_file)

        #################################### DOG ####################################

        elif brain_skullstrip == 'Custum_dog':

            epimask.run(input_for_msk, output_for_mask, [0.85, 0.95], False, ' -fill_holes -dilate_input 23',
                        sing_afni, sing_fs, diary_file)

        elif brain_skullstrip == 'sammba_dog':

            sammbamask.run(input_for_msk, output_for_mask, [1650, 0.85, 0.7, 500, 2, 10], sing_fs, diary_file)

            command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output_for_mask +
                       ' -input ' + output_for_mask + ' -fill_holes -dilate_input -5 7')
            run_cmd.run(command, diary_file)

        #################################### RAT ####################################
        ### specific for GD study

        elif brain_skullstrip == '3dSkullStrip_Rat':
            if ID in ['301502', '302101', '302105', '302106', '301603', '300908', '301500', '301501', '301503',
                      '301504', '301505', '301508', '301509', '302107', '302108', '300600']:

                afniskullstrip.run(input_for_msk, output_for_mask,' -orig_vol -mask_vol -rat',
                               ' -fill_holes -dilate_input 15',
                               sing_afni, diary_file)

            else:
                afniskullstrip.run(input_for_msk, output_for_mask,' -orig_vol -mask_vol -surface_coil -rat',
                               ' -fill_holes -dilate_input 4',
                               sing_afni, diary_file)

        elif brain_skullstrip == '3dSkullStrip_rat':

            afniskullstrip.run(input_for_msk, output_for_mask,' -blur_fwhm 1 -orig_vol -mask_vol -rat',
                           ' -fill_holes -dilate_input 2',
                           sing_afni, diary_file)

        elif brain_skullstrip == '3dSkullStrip_rat_no_dil':

            afniskullstrip.run(input_for_msk, output_for_mask,' -blur_fwhm 1 -orig_vol -mask_vol -rat',
                           ' -fill_holes',
                           sing_afni, diary_file)

        elif brain_skullstrip == 'sammba_rat':

            sammbamask.run(input_for_msk,output_for_mask,[2500,0.2,0.85,500,2,10],sing_fs,diary_file)


        elif brain_skullstrip == 'custum_rat':

            epimask.run(input_for_msk, output_for_mask, [0.8, 0.85], True, ' -fill_holes -dilate_input 3',
                        sing_afni, sing_fs, diary_file)

        #################################### General function ####################################

        elif brain_skullstrip.startswith('Vol_sammba_'):

            volume = int(brain_skullstrip.split('_')[2])
            sammbamask.run(input_for_msk, output_for_mask, [volume, 0.2, 0.8, 500, 2, 10], sing_fs, diary_file)

        elif brain_skullstrip.startswith('_bet'):
            # Extract the last two digits to use as the -f value and Create the approximate brain mask using bet2
            betmask.do(input_for_msk, output_for_mask, brain_skullstrip, sing_fsl, diary_file)


        elif brain_skullstrip == 'bet2_ANTS':

            tmp_mask = antsmask.create(input_for_msk,dir_prepro,masking_img,
                         type_of_transform,aff_metric_ants,BASE_SS_coregistr,BASE_SS_mask,
                         'nonlinear',1,'',sing_fsl,diary_file)

            antsmask.clean1(input_for_msk,output_for_mask,tmp_mask,BASE_SS_coregistr,BASE_SS_mask)


        elif brain_skullstrip == 'Custum_ANTS_NL':

            tmp_mask = antsmask.create(input_for_msk, dir_prepro, masking_img,
                         type_of_transform, aff_metric_ants, BASE_SS_coregistr, BASE_SS_mask,
                         'translation', 0, '',sing_fsl, diary_file)
            antsmask.clean1(input_for_msk, output_for_mask, tmp_mask, BASE_SS_coregistr, BASE_SS_mask)


        elif brain_skullstrip == 'Custum_ANTS_Garin':

            tmp_mask = antsmask.create(input_for_msk, dir_prepro, masking_img,
                                       type_of_transform, aff_metric_ants, BASE_SS_coregistr, BASE_SS_mask,
                                       'nonlinear', 0, '',sing_fsl, diary_file)

            antsmask.clean2(input_for_msk, output_for_mask, tmp_mask, BASE_SS_coregistr, BASE_SS_mask, ' -fill_holes', sing_afni, diary_file)


        elif brain_skullstrip.startswith('CustumANTSGarin_'):

            tmp_mask = antsmask.create(input_for_msk, dir_prepro, masking_img,
                                       type_of_transform, aff_metric_ants, BASE_SS_coregistr, BASE_SS_mask,
                                       'nonlinear', 0, '',sing_fsl, diary_file)

            # Extract the percentile from the string
            dilate_b = int(brain_skullstrip.split('_')[1])
            option2 = ' -fill_holes -dilate_input ' + str(dilate_b) + ' -' + str(dilate_b)

            antsmask.clean2(input_for_msk, output_for_mask, tmp_mask, BASE_SS_coregistr, BASE_SS_mask, option2,
                            sing_afni, diary_file)

        elif brain_skullstrip.startswith('CustumNilearn_'):

            # Extract the cutoff values from the string
            _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
            lower_cutoff = float(lower_cutoff)
            upper_cutoff = float(upper_cutoff)

            epimask.run(input_for_msk, output_for_mask, [lower_cutoff, upper_cutoff], False,
                        ' -fill_holes -dilate_input 1',sing_afni, sing_fs, diary_file)

            tmp_mask1 = ants.image_read(output_for_mask)
            antsmask.clean3(input_for_msk, output_for_mask, tmp_mask1,sing_afni, diary_file)


        elif brain_skullstrip.startswith('CustumNilearnExcludeZeros_'):

            # Extract the cutoff values from the string
            _, lower_cutoff, upper_cutoff = brain_skullstrip.split('_')
            lower_cutoff = float(lower_cutoff)
            upper_cutoff = float(upper_cutoff)

            epimask.run(input_for_msk, output_for_mask, [lower_cutoff, upper_cutoff], True,
                        ' -fill_holes -dilate_input 1',sing_afni, sing_fs, diary_file)

            tmp_mask1 = ants.image_read(output_for_mask)
            antsmask.clean3(input_for_msk, output_for_mask, tmp_mask1,sing_afni, diary_file)


        elif brain_skullstrip.startswith('CustumThreshold_'):
            # Handle 'Custum' type skullstripping dynamically for percentile
            percentile = int(brain_skullstrip.split('_')[1]) # Extract the percentile from the string

            loadimg     = nib.load(input_for_msk).get_fdata() # Load the image and calculate the threshold at the given percentile
            loadimgsort = np.percentile(np.abs(loadimg)[np.abs(loadimg) > 0], percentile)
            # Threshold the image using nilearn
            mask_imag = nilearn.image.threshold_img(input_for_msk, threshold=loadimgsort, cluster_threshold=10)
            mask_imag.to_filename(output_for_mask)
            # Use AFNI to process the mask
            command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output_for_mask +
                       ' -input ' + output_for_mask + ' -fill_holes -dilate_input 1')
            run_cmd.run(command, diary_file)

            tmp_mask1 = ants.image_read(output_for_mask)
            antsmask.clean3(input_for_msk, output_for_mask, tmp_mask1,sing_afni, diary_file)

        elif brain_skullstrip =='synthstrip':
            sstripmask.do(input_for_msk, output_for_mask, sing_synstrip, diary_file)

        elif brain_skullstrip == 'muSkullStrip_Human':
            Unetmask.do(input_for_msk,output_for_mask,'Site-Human-epoch_08.model',Unetpath,' -fill_holes',sing_afni,diary_file)

        elif brain_skullstrip == 'muSkullStrip_macaque':
            Unetmask.do(input_for_msk, output_for_mask, 'model-20_macaque-epoch', Unetpath, ' -fill_holes',
                        sing_afni, diary_file)

        elif brain_skullstrip == 'muSkullStrip_cross_species':
            Unetmask.do(input_for_msk, output_for_mask, 'model-02-_cross_species-epoch', Unetpath, ' -fill_holes',
                        sing_afni, diary_file)

        elif brain_skullstrip == 'muSkullStrip_cross_species_dil':
            Unetmask.do(input_for_msk, output_for_mask, 'model-02-_cross_species-epoch', Unetpath, ' -fill_holes -dilate_result 20 -10',
                        sing_afni, diary_file)

        elif brain_skullstrip == 'Custum_Macaque2':
            import re
            # Use a regular expression to find the session number (ses-XX)
            match = re.search(r'/ses-(\d+)/', masks_dir)
            if match:
                Session = match.group(1)
                nl = 'Session number: ' + Session
            else:
                nl = 'No session number found.'
            run_cmd.msg(nl, diary_file,'OKGREEN')

            if ID == 'Pickle' and Session == 5:
                Unetmask.do(input_for_msk, output_for_mask, 'model-20_macaque-epoch', Unetpath,
                            ' -fill_holes -dilate_input 10',
                            sing_afni, diary_file)

            elif ID == 'Trinity' and Session == 6:
                Unetmask.do(input_for_msk, output_for_mask, 'model-20_macaque-epoch', Unetpath,
                            ' -fill_holes -dilate_input 10',
                            sing_afni, diary_file)

            else:
                Unetmask.do(input_for_msk, output_for_mask, 'model-20_macaque-epoch', Unetpath,
                            ' -fill_holes -dilate_input 3',
                            sing_afni, diary_file)

        elif brain_skullstrip == 'QWARP':
            qwarpmask.run(input_for_msk, output_for_mask, dir_prepro, masks_dir, masking_img, BASE_SS_coregistr, BASE_SS_mask,
                          '',5,' -fill_holes',sing_afni, diary_file)

        elif brain_skullstrip =='Custum_QWARP':
            qwarpmask.run(input_for_msk, output_for_mask, dir_prepro, masks_dir, masking_img, BASE_SS_coregistr, BASE_SS_mask,
                          ' -ls', 5, ' -fill_holes', sing_afni, diary_file)

        elif brain_skullstrip =='Custum_QWARPT2':
            qwarpmask.run(input_for_msk, output_for_mask, dir_prepro, masks_dir, masking_img, BASE_SS_coregistr, BASE_SS_mask,
                          ' -lpa', 3, ' -fill_holes', sing_afni, diary_file)

        elif brain_skullstrip =='Custum_QWARPT2_dil':
            qwarpmask.run(input_for_msk, output_for_mask, dir_prepro, masks_dir, masking_img, BASE_SS_coregistr, BASE_SS_mask,
                          ' -lpa', 3, ' -fill_holes -dilate_input 2', sing_afni, diary_file)

        elif brain_skullstrip == 'NoSkullStrip':
            noSstrip.run(input_for_msk,output_for_mask,sing_afni,diary_file)

        elif brain_skullstrip == 'Manual':
            manualmask.do(input_for_msk, output_for_mask, sing_afni, preftool, sing_itk, sing_fs, diary_file)

        else:
            nl = "ERROR: brain_skullstrip not recognized, check that brain_skullstrip_1 or brain_skullstrip_2 are correctly written!!"
            raise Exception(run_cmd.error(nl, diary_file))

        if check_visualy_final_mask == True:
            manualmask.do(input_for_msk, output_for_mask, sing_afni, preftool, sing_itk, sing_fs, diary_file)

    return(output_for_mask)





