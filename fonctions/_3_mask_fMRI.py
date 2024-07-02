import os
import subprocess
import glob
import shutil
import sys
import nilearn
from nilearn.image import resample_to_img

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

from fonctions.extract_filename import extract_filename

def Refimg_to_meanfMRI(SED, anat_func_same_space, BASE_SS_coregistr, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, REF_int, ID, dir_prepro, n_for_ANTS, brainmask, V_mask, W_mask, G_mask, dilate_mask,
                       list_atlases, labels_dir, costAllin, anat_subject, IhaveanANAT, doMaskingfMRI, do_anat_to_func, Method_mask_func, lower_cutoff, upper_cutoff, overwrite):

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
    command = '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.nii.gz') + MEAN_im_list
    spco([command], shell=True)

    #################################### production of Mean image ####################################
    #command = '3dMean' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + MEAN_im_list
    #spco([command], shell=True)
    import nilearn
    from nilearn.image.image import mean_img
    mean_haxby = mean_img(MEAN_im_list_1)
    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz'))
    mean_haxby.to_filename(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_test.nii.gz'))


    command = '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_x = str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_y= str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_z = str(abs(round(float(spgo([command])[-8:]), 10)))

    command = '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    orient_meanimg = str(spgo([command]))

    command = '3dcalc' + overwrite + ' -a ' + BASE_SS_coregistr +  \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + ' -expr "a"'
    spco([command], shell=True)

    #command = '3drefit ' + overwrite + ' -orient ' + orient_meanimg + ' ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz')
    #spco([command], shell=True)

    command = '3dresample' + overwrite + \
    ' -orient ' + orient_meanimg + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
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
        command = '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + brainmask + ' -fill_holes -dilate_input ' + str(dilate_mask)
        spco(command, shell=True)
    else:
        command = '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + brainmask + ' -fill_holes'
        spco(command, shell=True)

    # B0 correction of the anat
    command = 'DenoiseImage' + ' -d 3 -i ' + anat_subject + ' -n Gaussian -r 3 ' + \
        ' -o ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz')
    spco([command], shell=True)


    if anat_func_same_space == True:
        command = 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' > ' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D')
        spco([command], shell=True)

        mvt_shft = opj(dir_prepro,ID + '_brain_for_Align_Center_inv.1D')
        command = 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' -I > ' + mvt_shft
        spco([command], shell=True)

        #don't work for two different 1d matrices... so lets do it speparatly....
        for input1, output2 in zip([opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz'), brainmask,
                                    opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), V_mask, W_mask, G_mask],
            [opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'),
             opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'),
            opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):


            ##### apply the recenter fmri

            command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + output2 + ' ' + input1 + ' -overwrite'
            spco([command], shell=True)

            command = '3dAllineate -final NN' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
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

        ####skullstrip the anat
        command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz') + ' -b ' + opj(
            dir_fMRI_Refth_RS_prepro2, 'orig_anat_for_fMRI.nii.gz') + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)

    else:
        for input1, output2 in zip([opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz'), brainmask, opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), V_mask, W_mask, G_mask],
            [opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'), 
            opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):

            if input1 == opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz'):

                command = '3dresample' + overwrite + \
                ' -prefix ' +opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
                ' -master ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + \
                ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz')
                spco([command], shell=True)

                ####skullstrip the anat
                command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + \
                          ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + ' -expr "a*b"'
                spco([command], shell=True)

                command = '3dresample' + overwrite + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz') + \
                ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz')
                spco([command], shell=True)


            else:
                command = '3dresample' + overwrite + \
                ' -prefix ' + output2 + \
                ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                ' -input  ' + input1
                spco([command], shell=True)

    ############################### ############################### ############################### 
    ############################ put anat IN Mean image space ##################################
    ############################### ############################### ###############################

    if doMaskingfMRI == True:
        if ope(opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz')):
                command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'manual_mask.nii.gz') + \
                ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
                spco([command], shell=True)
        else:
            if anat_func_same_space == True:
                command = '3dcalc' + overwrite + ' -a ' + opj(opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz')) + \
                          ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz') + ' -expr "a"'
                spco([command], shell=True)
                print('you are using the mask from the anat img!!!')
            else:
                if Method_mask_func=="3dAllineate":
                    ##### mask the func img
                    command = '3dAllineate' + overwrite + ' -cmass -EPI -final NN -float -twobest 5 -fineblur 0 -nomask -base ' + \
                    opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask.nii.gz') + \
                    ' -source ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' -' + costAllin + ' -1Dmatrix_save ' + \
                    opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask.1D') + \
                    ' -master ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz')
                    spco(command, shell=True)

                    command = 'cat_matvec ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask.1D') + \
                    ' -I > ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask_INV.1D')
                    spco([command], shell=True)

                    command = '3dAllineate' + overwrite + ' -final NN -1Dmatrix_apply ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask_INV.1D') + \
                    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + \
                    ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
                    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz')
                    spco([command], shell=True)

                elif Method_mask_func == "nilearn":
                    from nilearn.masking import compute_epi_mask
                    # convert to float
                    command = 'mri_convert -odt float ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')[:-7] + '_float.nii.gz'
                    spco([command], shell=True)
                    mask_img = compute_epi_mask(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')[:-7] + '_float.nii.gz', lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
                                                connected=True, opening=1,
                                                exclude_zeros=False, ensure_finite=True)
                    mask_img.to_filename(opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_nilearn.nii.gz'))
                    command = '3dmask_tool -overwrite -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat_Allineate_in_func.nii.gz') + \
                              ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_nilearn.nii.gz') + ' -fill_holes -dilate_input 1'
                    spco(command, shell=True)