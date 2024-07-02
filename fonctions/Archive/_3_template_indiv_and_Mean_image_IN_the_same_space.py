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

def Refimg_to_meanfMRI(anat_func_same_space, type_norm, Ref_file, volumes_dir, TfMRI, masks_dir, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, nb_run, recordings, REF_int, ID, TR, dir_prepro, n_for_ANTS, overwrite):

    #####################################################################################################################################
    #####################################################################################################################################
    ############################ Step 2 creat mean of all func (Mean image) into native anat space ######################################
    #####################################################################################################################################
    #####################################################################################################################################
    ###add the ref image to MEAN_im_list (otherwise will be forgoten in the following loop)
    root_RS_ref, extension_RS_ref = os.path.splitext(RS[REF_int])

    MEAN_im_list = ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS_ref + '_xdtr_mean_deob.nii.gz') #image "reference" (to be created)


    for r in range(0, int(nb_run)):
        root_RS, extension_RSr = os.path.splitext(RS[r])

        if not RS[r] == RS[REF_int]: # do not process ref...
            #add the coregitered mean image to the list
            MEAN_im_list = MEAN_im_list + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_') #add images in the same space


    #################################################################creat a mean image to use for the anat to func
    ###### average all func data and clean the image #####
    command = '3dTcat' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'residual_motion.nii.gz') + MEAN_im_list
    spco([command], shell=True)

    #################################### production of Mean image ####################################
    command = '3dMean' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + MEAN_im_list
    spco([command], shell=True)

    ''' XXX
    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_filter.nii.gz') + ' -expr "astep(a,1)"'
    spco([command], shell=True)
    command = '3dedge3' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_edge_tmp.nii.gz') + ' -input ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    spco([command], shell=True)
    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_edge_tmp.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_edge.nii.gz') + ' -expr "step(a)"'
    spco([command], shell=True)
    '''

    #creat anat space dir
    if ope(dir_fMRI_Refth_RS_prepro2) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro2)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro2,'tmp'))

    ############################### ############################### ############################### 
    ############################ put anat IN Mean image space ##################################
    ############################### ############################### ############################### 
    ##### first you need to re-create a dilate ref anat image (works better for co-registration anat to fMRI)
    #dilate MORE the "maskDilat"
    command = '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
    ' -input ' + brainmask_D + ' -fill_holes -dilate_input 2'
    spco(command, shell=True)

    mvt_shft = opj(dir_prepro,ID + '_brain_for_Align_Center_inv.1D')
    if os.path.isfile(mvt_shft) == False and os.path.isfile(opj(dir_prepro,ID + '_brain_for_Align_Center.1D')) == True:
        print('dd')
        command = 'cat_matvec ' + opj(dir_prepro,ID + '_brain_for_Align_Center.1D') + \
        ' -I > ' + mvt_shft
        spco([command], shell=True)


    ############################### ############################### ############################### 
    ##################### if they are in the samme space ##########################################
    ############################### ############################### ############################### 

    if anat_func_same_space == 'YES':

        ####then

        print('then lets use the anatomical mask!!!')

        #####XXX might not work... to test

        # mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + ' ' + brainmask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz')
        spco([command], shell=True)  

        # dilate mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_orig.nii.gz') + ' ' + brainmask_D + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_orig.nii.gz')
        spco([command], shell=True)  

        # white matter mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' ' + W_mask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz')
        spco([command], shell=True)  

        # ventricule mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' ' + V_mask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz')
        spco([command], shell=True) 

        ####skullstrip the indiv pre-template (and clean it)
        # B0 correction
        command = 'DenoiseImage' + ' -d 3 -i ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '.nii.gz') + ' -n Gaussian -r 3 ' + \
            ' -o ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '_NU_denoise.nii.gz')
        spco([command], shell=True)

        command = '3dresample' + overwrite + ' -master ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '_NU_denoise.nii.gz') + ' -prefix ' +  \
        opj(dir_prepro,ID + '_mask_ref_resempled.nii.gz') + ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz')
        spco([command], shell=True)

        command = '3dcalc' + overwrite + ' -a ' + opj(dir_prepro,ID + '_mask_ref_resempled.nii.gz') + ' -b ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '_NU_denoise.nii.gz') + \
        ' -prefix ' +  opj(dir_prepro, ID + '_' + TfMRI + '_mprage_reorient_NU_SS.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)

        command = '3dresample' + overwrite + ' -master ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' -prefix ' +  \
        opj(dir_fMRI_Refth_RS_prepro2,'anat_in_func.nii.gz') + ' -input ' + opj(dir_prepro, ID + '_' + TfMRI + '_mprage_reorient_NU_SS.nii.gz')
        spco([command], shell=True)

        #############  ############# ############# ############ ############
        ############# put masks in original anat (before linear mvt) #######
        ############# ############# ############# ############# ############

    ############################### ############################### ############################### 
    ##################### if they are NOT in the samme space ######################################
    ############################### ############################### ############################### 

    else:
        print('harder but we can do it !!! Hopefully...=(')

        # mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + ' ' + brainmask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + \
        ' -master ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz')
        spco([command], shell=True)  

        # dilate mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_orig.nii.gz') + ' ' + brainmask_D + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_pre.nii.gz') + \
        ' -master ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_orig.nii.gz')
        spco([command], shell=True)

        # white matter mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + ' ' + W_mask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + \
        ' -master ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz')
        spco([command], shell=True)  

        # ventricule mask
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + ' ' + V_mask + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -interp NN -1Dmatrix_apply ' + mvt_shft + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
        ' -master ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz')
        spco([command], shell=True) 

        ####skullstrip the indiv pre-template (and clean it)
        # B0 correction
        command = 'DenoiseImage' + ' -d 3 -i ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '.nii.gz') + ' -n Gaussian -r 3 ' + \
            ' -o ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '_NU_denoise.nii.gz')
        spco([command], shell=True)

        command = '3dresample' + overwrite + ' -master ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '_NU_denoise.nii.gz') + ' -prefix ' +  \
        opj(dir_prepro,ID + '_mask_ref_resempled.nii.gz') + ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_pre.nii.gz')
        spco([command], shell=True)

        command = '3dcalc' + overwrite + ' -a ' + opj(dir_prepro,ID + '_mask_ref_resempled.nii.gz') + ' -b ' + opj(dir_prepro, ID + '_brain_for_Align_Center' + TfMRI + '_NU_denoise.nii.gz') + \
        ' -prefix ' +  opj(dir_prepro, ID + '_' + TfMRI + '_mprage_reorient_NU_SS.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)

        #############  ############# ############# ###############
        ############# put masks in the template indiv anat #######
        ############# ############# ############# ################ 

        command = '3dAllineate' + overwrite + ' -cmass -EPI -final cubic -float -twobest 5 -fineblur 0 -base ' + \
        opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_in_func_Allin_AFNI.nii.gz') + \
        ' -source ' + opj(dir_prepro, ID + '_' + TfMRI + '_mprage_reorient_NU_SS.nii.gz') + ' -nmi -1Dmatrix_save ' + \
        opj(dir_prepro, ID + 'anat_in_func_Allin_AFNI.1D')
        spco(command, shell=True)

        ##### apply that to the other masks
        command = '3dAllineate' + overwrite + ' -interp NN -overwrite -1Dmatrix_apply ' + opj(dir_prepro, ID + 'anat_in_func_Allin_AFNI.1D') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_in_func_Allin_AFNI.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat_pre.nii.gz')
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -interp NN -overwrite -1Dmatrix_apply ' + opj(dir_prepro, ID + 'anat_in_func_Allin_AFNI.1D') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_in_func_Allin_AFNI.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz')
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -interp NN -overwrite -1Dmatrix_apply ' + opj(dir_prepro, ID + 'anat_in_func_Allin_AFNI.1D') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_in_func_Allin_AFNI.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz')
        spco([command], shell=True)

        command = '3dAllineate' + overwrite + ' -interp NN -overwrite -1Dmatrix_apply ' + opj(dir_prepro, ID + 'anat_in_func_Allin_AFNI.1D') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_in_func_Allin_AFNI.nii.gz') + \
        ' -input  ' + opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz')
        spco([command], shell=True)


    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################

    ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')
    command = '3dcalc' + overwrite + ' -a ' +  opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_fMRImage_SS.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    #keep only func space
    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_filter.nii.gz') + ' -b ' + opj(dir_fMRI_Refth_RS_prepro2,'anat_in_func_Allin_AFNI.nii.gz') + \
    ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_template_forFMRI.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    #average indiv anat =>coreg_ref
    # mean of all func SS => opj(dir_fMRI_Refth_RS_prepro1,'Mean_fMRImage_SS.nii.gz')

    command = 'antsRegistration -d 3 --float 0 --verbose 1 -w [0.01,0.99] -n Linear' + \
    ' -o [' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_') + ',' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped.nii.gz') + ']' + \
    ' -t Affine[0.1] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
    ' -m MI[' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_fMRImage_SS.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_template_forFMRI.nii.gz') + ',1,32,Regular,0.2]' + \
    ' -t Syn[0.1,3,0] -f 8x4x2x1 -s 3x2x1x0vox -c [1000x500x250x100,1e-6,10]' + \
    ' -m Mattes[' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_fMRImage_SS.nii.gz') + ',' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_fMRImage_SS.nii.gz') + ',1,4,Regular,0.2]'
    spco([command], shell=True)

    ########## ########## help for antsRegistration ########## ########## ########## ##########
    # -f -s -c X1 x X2 x X3 x X4 from large mvt to small si 4 prend toujours 4
    ## erreur accepté a la fin de la coregistration
    # -w [0.01,0.99] threshold
    # function pour -u pour normaliser les images entre elles
    #  -n Linear ==> transfo comme NN
    ########## ########## ########## ########## ########## ########## ########## ##########

    #XXX why?
    ###extract edge of the func
    command = '3dedge3' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_edge_tmp.nii.gz') + ' -input ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped.nii.gz')
    spco([command], shell=True)

    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_edge_tmp.nii.gz') + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_edge.nii.gz') + ' -expr "step(a)"'  
    spco([command], shell=True)

    command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_edge.nii.gz') +  ' -b ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_edge_brain.nii.gz') + ' -expr "a*b"'
    spco([command], shell=True)

    os.remove(opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_edge_tmp.nii.gz'))



    #########################################################################################################
    ################################### registration to anat space ##########################################
    #########################################################################################################

    ############################### ############################### ############################### 
    ############################### apply transfo to anat space to each volume of each func image #
    ############################### ############################### ############################### 

    MEAN = ''
    for i in range(0, int(nb_run)):

        root_RS, extension_RSr = os.path.splitext(RS[i])
        """
        if i == 0 and recordings == 'new':
            command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','dummy.nii.gz') + ' ' +  opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf.nii.gz')

        elif i == 0 and recordings == 'old':
            command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','dummy.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_deob.nii.gz')

        elif i == 0 and recordings == 'very_old':
            command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','dummy.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr.nii.gz')

        else:
        """
        command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','dummy.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz')
        spco([command], shell=True)

        #####transfo originally linearXXXX (line 157) should be linear too here????
        WW = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro2,'tmp','*.nii.gz')))
        for j in range(0, len(WW)):
            if j < 10:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n ' + 'Linear ' + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_template_forFMRI.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.00' + str(j) + '.nii.gz') + \
                ' -t ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_0GenericAffine.mat') + \
                ' -t ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_1Warp.nii.gz')
            elif j < 100 and j >= 10:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n ' + 'Linear ' + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_template_forFMRI.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.0' + str(j) + '.nii.gz') + \
                ' -t ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_0GenericAffine.mat') + \
                ' -t ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_1Warp.nii.gz')
            else:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n ' + 'Linear ' + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro2, TfMRI + '_template_forFMRI.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.0' + str(j) + '.nii.gz') + \
                ' -t ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_0GenericAffine.mat') + \
                ' -t ' + opj(dir_fMRI_Refth_RS_prepro2,'Mean_Image_unwarped_1Warp.nii.gz')
            
            spco([command], shell=True)
        list_name = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro2,'tmp','warp1.*.nii.gz')))
        name = ''
        for m in range(0,len(list_name)):
            name = name + ' ' + list_name[m]
        command = '3dTcat' + overwrite + ' ' + name + ' -tr ' + str(TR) + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfw.nii.gz')
        spco([command], shell=True)

        '''
        if extract_exterior_CSF == True
            #dilate MORE the "maskDilat"
            command = '3dmask_tool' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_erod_forextline.nii.gz') + \
            ' -input ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_erod.nii.gz') + ' -fill_holes -dilate_input -1'
            spco(command, shell=True)

            command = '3dcalc' + overwrite + ' -a ' + opj(mask_erod_forextline,'mask_erod.nii.gz') +  ' -b ' + opj(dir_fMRI_Refth_RS_prepro2,'mask_erod.nii.gz') + \
            ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2,'exterior_ligne.nii.gz') + ' -expr "b-a"'
            spco([command], shell=True)
        '''
        
        ####keep only signal in mask dilate
        command = '3dcalc -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfw.nii.gz') +  ' -b ' + opj(dir_fMRI_Refth_RS_prepro2,'maskDilat.nii.gz') + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfw_masked.nii.gz') + ' -expr "a*b"'
        spco([command], shell=True)
        
        for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro2,'tmp')):
            os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro2,'tmp'), f))
        
        command = '3dTstat' + overwrite + ' -mean -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfw_mean.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfw.nii.gz')
        spco([command], shell=True)
        MEAN = MEAN + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfw_mean.nii.gz')

    command = '3dMean' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + MEAN
    spco([command], shell=True)

    if check_visualy_each_img == 'YES':
        command = 'itksnap -g ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz') + ' -s ' + opj(dir_fMRI_Refth_RS_prepro2, '_mask_ref_resempled_func.nii.gz')
        spco([command], shell=True)