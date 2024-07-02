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

def Refimg_to_meanfMRI(anat_func_same_space, type_norm, Ref_file, volumes_dir, BASE_SS_coregistr, TfMRI, masks_dir, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, 
    dir_fMRI_Refth_RS_prepro3, RS, nb_run, recordings, REF_int, ID, TR, dir_prepro, n_for_ANTS, brainmask_D, brainmask, V_mask, W_mask, dilate_mask, list_atlases, overwrite):

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

    ############################### ############################### ############################### 
    ##################### if anat and func are in the samme space #################################
    ############################### ############################### ############################### 
    ##recenter mean Mean_Image (will be the co-registration reference)
    if anat_func_same_space == True:
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz')
        spco([command], shell=True) 


    ############################### ############################### ############################### 
    ##################### if they are NOT in the samme space ######################################
    ############################### ############################### ############################### 
    else:     
        #dim = nib.load(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')).header.get_zooms()
        #dim2 = ''
        #for x in dim:
        #    dim2 = dim2 + ' ' + str(round(x, 1))
        command = '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        delta_x = str(round(float(spgo([command])[-8:]), 2))
        command = '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        delta_y= str(round(float(spgo([command])[-8:]), 2))
        command = '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        delta_z = str(round(float(spgo([command])[-8:]), 2))

        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
        ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
        ' -rmode Cu -input  ' + BASE_SS_coregistr
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -cmass -EPI -final cubic -float -twobest 5 -fineblur 0 -base ' + \
        opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
        ' -source ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' -nmi -1Dmatrix_save ' + \
        opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.1D')
        spco(command, shell=True)
        recenter_fmri = opj(dir_fMRI_Refth_RS_prepro1,'_brain_for_Align_Center.1D')

    ############################### ############################### ############################### 
    ##################### if anat and func are in the samme space #################################
    ############################### ############################### ############################### 
    #### center each func images
    for r in range(0, int(nb_run)):
        root_RS, extension_RSr = os.path.splitext(RS[r])
        if anat_func_same_space == True:
            ###recenter to an acceptable center...
            command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref_RcT.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz') + ' -overwrite'
            spco([command], shell=True)
            command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref_RcT.nii.gz') + \
            ' -input  ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref_RcT.nii.gz')
            spco([command], shell=True) 

    ############################### ############################### ############################### 
    ##################### if they are NOT in the samme space ######################################
    ############################### ############################### ############################### 
        else:
            command = '3dAllineate' + overwrite + ' -interp NN -overwrite -1Dmatrix_apply ' + recenter_fmri + \
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
        ' -input ' + brainmask_D + ' -fill_holes -dilate_input ' + str(dilate_mask)
        spco(command, shell=True)
    else:
        command = '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + \
        ' -input ' + brainmask_D + ' -fill_holes'
        spco(command, shell=True)

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
    ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz')
    spco([command], shell=True)

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz') + \
    ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
    ' -input  ' + brainmask
    spco([command], shell=True)

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + \
    ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
    ' -input  ' + V_mask
    spco([command], shell=True)

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz') + \
    ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
    ' -input  ' + W_mask
    spco([command], shell=True)

    ############################### ############################### ############################### 
    ############################ put anat IN Mean image space ##################################
    ############################### ############################### ############################### 
    ####skullstrip the indiv pre-template (and clean it)
    # B0 correction
    command = 'DenoiseImage' + ' -d 3 -i ' +  opj(dir_prepro, ID + '_acpc' + TfMRI + '.nii.gz') + ' -n Gaussian -r 3 ' + \
        ' -o ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz')
    spco([command], shell=True)

    command = '3dresample' + overwrite + ' -master ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + ' -prefix ' +  \
    opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz')
    spco([command], shell=True)

    ####skullstrip the anat
    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilatanat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise_maskeDil.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    ####resample the anat as reference for realignment
    command = '3dresample' + overwrite + ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + ' -prefix ' +  \
    opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ' -idnput ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_NU_denoise_maskeDil.nii.gz')
    spco([command], shell=True)

    ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    command = '3dcalc' + overwrite + ' -a ' +  opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################
    #average indiv anat =>coreg_ref
    # mean of all func SS => opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz')

    command = 'antsRegistration -d 3 --float 0 --verbose 1 -w [0.01,0.99] -n ' + n_for_ANTS + \
    ' -o [' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_') + ',' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped.nii.gz') + ']' + \
    ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
    ' -m Mattes[' + opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + ',' +opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS.nii.gz') + ',1,4,Regular,0.2]'
    spco([command], shell=True)

    mvt_shft_INV_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1InverseWarp.nii.gz')
    mvt_shft_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1Warp.nii.gz')

    ########## ########## help for antsRegistration ########## ########## ########## ##########
    # -f -s -c X1 x X2 x X3 x X4 from large mvt to small si 4 prend toujours 4
    ## erreur accepté a la fin de la coregistration
    # -w [0.01,0.99] threshold
    # function pour -u pour normaliser les images entre elles
    #  -n Linear ==> transfo comme NN
    ########## ########## ########## ########## ########## ########## ########## ##########

    mvt_shft = opj(dir_prepro,ID + '_brain_for_Align_Center_inv.1D')
    if os.path.isfile(mvt_shft) == False and os.path.isfile(opj(dir_prepro,ID + '_brain_for_Align_Center.1D')) == True:
    print('dd')
    command = 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
    ' -I > ' + mvt_shft
    spco([command], shell=True)

    # mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz') + \
    ' --interpolation nearestNeighbor'  + \
    mvt_shft_ANTs
    spco([command], shell=True)

    command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
    ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
    ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz')
    spco([command], shell=True) 

    # dilate mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
    ' --interpolation nearestNeighbor' + \
    mvt_shft_ANTs
    spco([command], shell=True)

    # white matter mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz') + \
    ' --interpolation nearestNeighbor' + \
    mvt_shft_ANTs
    spco([command], shell=True)

    # ventricule mask
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + \
    ' --interpolation nearestNeighbor' + \
    mvt_shft_ANTs
    spco([command], shell=True)

    #### apply to all atlases
    for atlas in list_atlases:
        ## in anat space resample to func 
        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_RcT.nii.gz') + \
        ' -input  ' + opj(labels_dir, TfMRI + opb(atlas))
        spco([command], shell=True)

        ## in func space resample to func 
        command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro2, opb(atlas)) + ' -n ' + n_for_ANTS + \
        ' -r ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS.nii.gz') + \
        ' -o ' + opj(dir_fMRI_Refth_RS_prepro1, opb(atlas)) + \
        ' --interpolation nearestNeighbor' + \
        mvt_shft_ANTs
        spco([command], shell=True)

    for i in range(0, int(nb_run)):
        root_RS, extension_RSr = os.path.splitext(RS[i])
        ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
        
        command = '3dcalc' + overwrite + ' -a ' +   opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref_RcT.nii.gz') + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, '_xdtrf_2ref_RcT_masked.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)