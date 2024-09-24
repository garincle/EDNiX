import os
import subprocess
from nilearn.image import resample_to_img
from fonctions.extract_filename import extract_filename
from nilearn.image.image import mean_img
import fonctions.Skullstrip_func
import ants

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def Refimg_to_meanfMRI(anat_func_same_space, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                       dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                       costAllin, anat_subject, doMaskingfMRI, Method_mask_func, lower_cutoff, upper_cutoff, overwrite, type_of_transform, aff_metric_ants,
                       s_bind,afni_sif,fs_sif, fsl_sif, itk_sif):

    ##### creat new variable for template space (we will need to store and downsample template img to func resolution)

    if ope(dir_fMRI_Refth_RS_prepro3) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro3)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro3,'tmp'))

    #####################################################################################################################################
    #####################################################################################################################################
    ############################ Step 2 creat mean of all func (Mean image) into native anat space ######################################
    #####################################################################################################################################
    #####################################################################################################################################
    ### creat a list of the image to be corrected
    ###add the ref image to MEAN_im_list (otherwise will be forgoten in the following loop)
    root_RS_ref = extract_filename(RS[REF_int])

    MEAN_im_list = ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS_ref + '_xdtr_mean_deob.nii.gz') #image "reference" (to be created)
    MEAN_im_list_1 = [opj(dir_fMRI_Refth_RS_prepro1,root_RS_ref + '_xdtr_mean_deob.nii.gz')] #image "reference" (to be created)


    for r in range(0, int(nb_run)):
        root_RS = extract_filename(RS[r])

        if not RS[r] == RS[REF_int]: # do not process ref...
            #add the coregitered mean image to the list
            MEAN_im_list = MEAN_im_list + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz') #add images in the same space
            MEAN_im_list_1.append(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz'))

    ################################################################# creat a mean image to use for the anat to func and recenter
    ###### average all func data and clean the image #####
    command = 'singularity run' + s_bind + afni_sif + '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.nii.gz') + MEAN_im_list
    spco([command], shell=True)

    #################################### production of Mean image ####################################
    #command = '3dMean' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + MEAN_im_list
    #spco([command], shell=True)

    mean_haxby = mean_img(MEAN_im_list_1)
    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'))
    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.nii.gz'))


    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    ADI = spgo(command).split('\n')
    delta_x = str(abs(round(float(ADI[-1]), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    ADJ = spgo(command).split('\n')
    delta_y= str(abs(round(float(ADJ[-1]), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    ADK = spgo(command).split('\n')
    delta_z = str(abs(round(float(ADK[-1]), 10)))

    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    orient_meanimg = spgo(command).split('\n')[2]

    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + BASE_SS_coregistr +  \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + ' -expr "a"'
    spco([command], shell=True)

    #command = '3drefit ' + overwrite + ' -orient ' + orient_meanimg + ' ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz')
    #spco([command], shell=True)

    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
    ' -orient ' + orient_meanimg + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + \
    ' -rmode Cu -input  ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz')
    spco([command], shell=True)
    '''

    caca2 = resample_to_img(opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'), interpolation='nearest')
    caca2.to_filename(opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'))
    '''
    #creat anat space dir
    if ope(dir_fMRI_Refth_RS_prepro2) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro2)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro2,'tmp'))

    ############################### ############################### ############################### 
    ##################### resample the masks for signal extraction ################################
    ############################### ############################### ############################### 
    ##### first you need to re-create a dilate ref anat image (works better for co-registration anat to fMRI)
    #dilate MORE the "maskDilat"
    if dilate_mask >0:
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + brainmask + ' -fill_holes -dilate_input ' + str(dilate_mask)
        spco(command, shell=True)
    else:
        command = 'singularity run' + s_bind + afni_sif + '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + brainmask + ' -fill_holes'
        spco(command, shell=True)

    # B0 correction of the anat   < ====== Why ?
    IMG = ants.image_read(anat_subject)
    IMG = ants.denoise_image(IMG,r=3,noise_model='Gaussian')
    ants.image_write(IMG,opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz'),ri=False)

    if anat_func_same_space == True:
        command = 'singularity run' + s_bind + afni_sif + 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' | tail -n +3 >' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D')
        spco([command], shell=True)

        mvt_shft = opj(dir_prepro,ID + '_brain_for_Align_Center_inv.1D')
        command = 'singularity run' + s_bind + afni_sif + 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' -I | tail -n +3 > ' + mvt_shft
        spco([command], shell=True)

        #don't work for two different 1d matrices... so lets do it speparatly....
        for input1, output2 in zip([opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz'), brainmask,
                                    opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), V_mask, W_mask, G_mask],
            [opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'),
             opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'),
            opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):
            print(input1)
            print(ope(input1))
            if ope(input1):
                ##### apply the recenter fmri

                command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + output2 + ' ' + input1 + ' -overwrite'
                spco([command], shell=True)

                command = 'singularity run' + s_bind + afni_sif + '3dAllineate -final NN' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
                ' -prefix ' + output2 + \
                ' -master ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz') + \
                ' -input  ' + output2
                spco([command], shell=True)
                print('to check !!')

                #' -master ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + \
                ### do not work with macaque!!!!!!!!!!!!!!!!!!!!!!

                #command = '3drefit -atrcopy ' + opj(dir_prepro, ID + '_mprage_reorient' + TfMRI + '.nii.gz') + ' IJK_TO_DICOM_REAL ' + output2
                #spco([command], shell=True)

                caca2 = resample_to_img(output2, opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'), interpolation='nearest')
                caca2.to_filename(output2)
            else:
                print('WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, then no extraction of WM or Ventricules or GM will be possible... pls check that!!!!')
        ####skullstrip the anat
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + ' -b ' + opj(
            dir_fMRI_Refth_RS_prepro2, 'orig_anat_for_fMRI.nii.gz') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)
    else:
        for input1, output2 in zip([opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz'), brainmask, opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), V_mask, W_mask, G_mask],
            [opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'), 
            opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):
            if ope(input1):
                if input1 == opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz'):

                    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                    ' -prefix ' +opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
                    ' -master ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + \
                    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz')
                    spco([command], shell=True)
                    ####skullstrip the anat
                    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + \
                              ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + ' -expr "a*b"'
                    spco([command], shell=True)
                    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + \
                    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz')
                    spco([command], shell=True)

                else:
                    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                    ' -prefix ' + output2 + \
                    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                    ' -input  ' + input1
                    spco([command], shell=True)
            else:
                print('WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, then no extraction of WM or Ventricules or GM will be possible... pls check that!!!!')


    ############################### ############################### ############################### 
    ############################ put anat IN Mean image space ##################################
    ############################### ############################### ###############################

    if doMaskingfMRI == True:
        if ope(opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz')):
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz') + \
                ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
                spco([command], shell=True)
        else:
            if anat_func_same_space == True:
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz')) + \
                          ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
                spco([command], shell=True)

                print('you are using the mask from the anat img!!!')
            else:
                fonctions.Skullstrip_func.Skullstrip_func(Method_mask_func, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, overwrite, costAllin, lower_cutoff, upper_cutoff, type_of_transform, aff_metric_ants, s_bind, afni_sif, fsl_sif, fs_sif, itk_sif)

    else:
        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + \
                  ' -a  ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
                  ' -expr "step(a)" ' + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz')
        spco([command], shell=True)
