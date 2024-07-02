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

def Refimg_to_meanfMRI(anat_func_same_space, volumes_dir, BASE_SS_coregistr, TfMRI, correction_direction, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, 
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, recordings, REF_int, ID, TR, dir_prepro, n_for_ANTS, brainmask_D, brainmask, V_mask, W_mask, G_mask, dilate_mask, list_atlases, labels_dir, costAllin, overwrite):

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
    root_RS_ref, extension_RS_ref = os.path.splitext(RS[REF_int])

    MEAN_im_list = ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS_ref + '_xdtr_mean_deob.nii.gz') #image "reference" (to be created)


    for r in range(0, int(nb_run)):
        root_RS, extension_RSr = os.path.splitext(RS[r])

        if not RS[r] == RS[REF_int]: # do not process ref...
            #add the coregitered mean image to the list
            MEAN_im_list = MEAN_im_list + ' ' + opj(dir_fMRI_Refth_RS_prepro1, + '_xdtrf_2ref.nii.gz') #add images in the same space


    ################################################################# creat a mean image to use for the anat to func and recenter
    ###### average all func data and clean the image #####
    command = '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.nii.gz') + MEAN_im_list
    spco([command], shell=True)

    #################################### production of Mean image ####################################
    command = '3dMean' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + MEAN_im_list
    spco([command], shell=True)

    command = '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_x = spgo([command])
    command = '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_y= spgo([command])
    command = '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    delta_z = spgo([command])

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
    ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
    ' -rmode Cu -input  ' + BASE_SS_coregistr
    spco([command], shell=True)

    #recenter func img
    os.chdir(dir_fMRI_Refth_RS_prepro1)
    command = '@Align_Centers -overwrite -base ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + ' -dset ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' -prefix ' + ' Mean_Image_RcT.nii.gz'
    spco([command], shell=True)
    os.chdir('/home/cgarin/')

    recenter_fmri = opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.1D')

    ############################### ############################### ############################### 
    ##################### if anat and func are in the samme space #################################
    ############################### ############################### ############################### 
    #### center each func images
    for r in range(0, int(nb_run)):
        root_RS, extension_RSr = os.path.splitext(RS[r])
        command = '3dAllineate' + overwrite + ' -final NN -warp shift_only -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref_RcT.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz')
        spco([command], shell=True)

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
        ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + ' -fill_holes -dilate_input ' + str(dilate_mask)
        spco(command, shell=True)
    else:
        command = '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + ' -fill_holes'
        spco(command, shell=True)

    if anat_func_same_space == True:

        command = 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' > ' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D')
        spco([command], shell=True)

        mvt_shft = opj(dir_prepro,ID + '_brain_for_Align_Center_inv.1D')
        command = 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' -I > ' + mvt_shft
        spco([command], shell=True)

        #don't work for two different 1d matrices... so lets do it speparatly....
        #concate_1D = opj(dir_fMRI_Refth_RS_prepro2, 'concate_1D.1D')
        #command = 'cat ' + mvt_shft + ' ' + recenter_fmri + ' > ' + concate_1D
        #spco([command], shell=True)

        ##### apply the recenter fmri
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz')
        spco([command], shell=True)

        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz')
        spco([command], shell=True)

        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz')
        spco([command], shell=True)

        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz')
        spco([command], shell=True)

        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz')
        spco([command], shell=True)

        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz') + \
        ' -input  ' +  opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')
        spco([command], shell=True)

        ##### apply the recenter fmri
        command = '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz')
        spco([command], shell=True)

        ##### apply the recenter fmri
        command = '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' +opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz')
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz')
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + brainmask
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + V_mask
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + W_mask
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + recenter_fmri + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + G_mask
        spco([command], shell=True)

    else:
        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + \
        ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
        ' -input  ' + opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
        spco([command], shell=True)

        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz')
        spco([command], shell=True)

        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + \
        ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
        ' -input  ' + brainmask
        spco([command], shell=True)

        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
        ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
        ' -input  ' + V_mask
        spco([command], shell=True)

        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + \
        ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
        ' -input  ' + W_mask
        spco([command], shell=True)

        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz') + \
        ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
        ' -input  ' + G_mask
        spco([command], shell=True)


    ############################### ############################### ############################### 
    ############################ put anat IN Mean image space ##################################
    ############################### ############################### ############################### 
    ####skullstrip the indiv pre-template (and clean it)
    # B0 correction
    command = 'DenoiseImage' + ' -d 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz') + ' -n Gaussian -r 3 ' + \
        ' -o ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz')
    spco([command], shell=True)

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + \
    ' -master ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz')
    spco([command], shell=True)

    ####skullstrip the anat
    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    if anat_func_same_space == True:

        print('lets use the anat mask')

    else:
        ##### mask the func img 
        command = '3dAllineate' + overwrite + ' -cmass -EPI -final NN -float -twobest 5 -fineblur 0 -base ' + \
        opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask.nii.gz') + \
        ' -source ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + ' -' + costAllin + ' -1Dmatrix_save ' + \
        opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask.1D')
        spco(command, shell=True)

        command = 'cat_matvec ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask.1D') + \
        ' -I > ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask_INV.1D')
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_for_mask_INV.1D') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz')
        spco([command], shell=True)

    if anat_func_same_space == True:
    for input1, output2 in zip([opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')], [opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_fMRI.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'), opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')]):
        
        command = '3dWarp' + overwrite + ' -card2oblique ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI) + ' -prefix ' + output2 + ' ' + input1
        spco([command], shell=True)


    #### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    command = '3dcalc' + overwrite + ' -a ' +  opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################
    #average indiv anat =>coreg_ref
    # mean of all func SS => opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz')


    ###### a remplacer par afni?? ou juste la dernière étape?? XXXXX
    if 'x' in correction_direction:
        restrict = '1x0.1x0.1'
    elif 'y' in correction_direction:
        restrict = '0.1x1x0.1' 
    elif 'x' in correction_direction:
        restrict = '0.1x0.1x1'

    command = 'antsRegistration -d 3 --float 0 --verbose 1 -w [0.01,0.99] -n ' + n_for_ANTS + \
    ' -o [' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_') + ',' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped.nii.gz') + ']' + \
    ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
    ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS.nii.gz') + ',1,32,Regular,0.2]' + \
    ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
    ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS.nii.gz') + ',1,4,Regular,0.2]'  + \
    ' --restrict-deformation ' + restrict
    spco([command], shell=True)

    mvt_shft_INV_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1InverseWarp.nii.gz') + \
    ' -t ' + str([opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_0GenericAffine.mat'), 1])

    mvt_shft_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1Warp.nii.gz') + \
    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_0GenericAffine.mat')

    ########## ########## help for antsRegistration ########## ########## ########## ##########
    # -f -s -c X1 x X2 x X3 x X4 from large mvt to small si 4 prend toujours 4
    ## erreur accepté a la fin de la coregistration
    # -w [0.01,0.99] threshold
    # function pour -u pour normaliser les images entre elles
    #  -n Linear ==> transfo comme NN
    ########## ########## ########## ########## ########## ########## ########## ##########

    # mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz') + \
    ' --interpolation nearestNeighbor'  + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)

    # dilate mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
    ' --interpolation nearestNeighbor' + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)

    # white matter mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz') + \
    ' --interpolation nearestNeighbor' + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)

    # ventricule mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + \
    ' --interpolation nearestNeighbor' + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)

    # Gray mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'Gmask.nii.gz') + \
    ' --interpolation nearestNeighbor' + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)

    #### apply to all atlases
    for atlas in list_atlases:
        ## in anat space resample to func 

        if anat_func_same_space == True:
            command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' -overwrite'
            spco([command], shell=True)
            command = '3dAllineate' + overwrite + ' -final NN -overwrite -1Dmatrix_apply ' + mvt_shft + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
            ' -input  ' +  opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
            spco([command], shell=True)

            ##### apply the recenter fmri
            command = '3dAllineate -final NN' + overwrite + ' -overwrite -1Dmatrix_apply ' + recenter_fmri + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
            ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
            ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
            spco([command], shell=True)

            ## in anat space resample to func 
            command = '3dWarp' + overwrite + ' -card2oblique ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI) + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas))
            spco([command], shell=True)

        else:
            command = '3dresample' + overwrite + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
            ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
            ' -input  ' + opj(labels_dir, TfMRI + opb(atlas))
            spco([command], shell=True)

        ## in func space resample to func 
        command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' -n ' + n_for_ANTS + \
        ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
        ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, opb(atlas)) + \
        ' --interpolation nearestNeighbor' + \
        mvt_shft_INV_ANTs
        spco([command], shell=True)

    ## in func space resample to func 
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI.nii.gz') + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)


    command = '3dinfo -di ' + opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
    delta_x1 = spgo([command])
    command = '3dinfo -dj ' + opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
    delta_y1 = spgo([command])
    command = '3dinfo -dk ' + opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
    delta_z1 = spgo([command])

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz') + \
    ' -dxyz ' + delta_x1 + ' ' + delta_y1 + ' ' + delta_z1 + ' ' + \
    ' -rmode Cu -input  ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz')
    spco([command], shell=True)

    #### creat a nice anat in func space
    if anat_func_same_space == True:
        anatstd = opj(dir_fMRI_Refth_RS_prepro2,'orig_anat_for_plot.nii.gz')

    else:
        anatstd = opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')

    ## in func space resample to func 
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + anatstd + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.nii.gz') + \
    mvt_shft_INV_ANTs
    spco([command], shell=True)


    ###finaly mask the func with mask
    for i in range(0, int(nb_run)):
        root_RS, extension_RSr = os.path.splitext(RS[i])
        ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        
        command = '3dcalc' + overwrite + ' -a ' +   opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref_RcT.nii.gz') + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, '_xdtrf_2ref_RcT_masked.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)
