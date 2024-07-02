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

####################################################################################
########################## Step 3 normalisation to template atlas space ############
####################################################################################
def to_common_template_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, BASE_SS_coregistr, 
                nb_run, RS, transfo_concat_inv, TR, n_for_ANTS, list_atlases, TfMRI, BASE_SS_mask, GM_mask, GM_mask_studyT, creat_sutdy_template, 
                anat_func_same_space, orientation, path_anat, ID, overwrite):

    mvt_shft_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1Warp.nii.gz') + \
    ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_0GenericAffine.mat')

    ##### creat new variable  for template space
    if ope(dir_fMRI_Refth_RS_prepro3) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro3)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro3,'tmp'))

        #########################################################################################################
        ################################### registration to anat space ##########################################
        #########################################################################################################

        #### transfo to brain template space (centred) are
        #    mvt_shft_ANTs = ' -t ' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_1InverseWarp.nii.gz') + \
        #' -t ' + str([opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_unwarped_0GenericAffine.mat'), 1])

        ############################### ############################### ############################### 
        ############################### apply transfo to anat space to each volume of each func image #
        ############################### ############################### ###############################

    if anat_func_same_space == True:

        for input2, output2 in zip([opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS.nii.gz')], [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz')]):
        
            command = '3dcalc' + overwrite + ' -a ' +   input2 + \
            ' -prefix ' +  output2 + ' -expr "a"'
            spco([command], shell=True)

            command = '3dinfo -orient ' + opj(path_anat, ID + 'template_indiv' + TfMRI + '.nii.gz') + ' ' +  output2
            orientation_orig = spgo([command], shell=True)

            if str(nx)=="b'1\\n1\\n'":
                print('same oblique')

                # reorient the fiedls according to the json file
                command = '3dWarp' + overwrite + ' -deoblique -wsinc5 -prefix ' +  output2 + \
                ' ' +  output2
                spco([command], shell=True)

            else:
                orientation_orig = '3dinfo -same_obl ' +opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz')

                print('oblique different')
                # reorient the fiedls according to the json file
                command = '3drefit' + overwrite + ' -orient ' + orientation_orig + ' ' +  output2
                spco([command], shell=True)
                command = '3dWarp' + overwrite + ' -deoblique -wsinc5 -prefix ' +  output2 + \
                ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
                spco([command], shell=True)
                command = '3drefit' + overwrite + ' -orient ' + orientation + ' ' +  output2
                spco([command], shell=True)


            ##### apply the recenter fmri
            command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' +  output2 + ' ' +  output2 + ' -overwrite'
            spco([command], shell=True)
            command = '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
            ' -prefix ' + output2 + \
            ' -input  ' + output2
            spco([command], shell=True)
            command = '3dAutobox' + overwrite + ' -input ' + output2 + ' -prefix ' + output2 + ' -noclust  -overwrite'
            spco(command, shell=True)


    else: 
        command = '3dcalc' + overwrite + ' -a ' +   opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS.nii.gz') + \
        ' -prefix ' +  opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz') + ' -expr "a"'
        spco([command], shell=True)

    ## test on mean img (to see spatially that is works)
    command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz') + ' -n ' + n_for_ANTS + \
    ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
    ' -o ' + opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz') + \
    mvt_shft_ANTs + \
    transfo_concat_inv
    spco([command], shell=True)



    ##### go for BOLD img preTTT
    for i in range(0, int(nb_run)):
        root_RS, extension_RSr = os.path.splitext(RS[i])
        for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro3,'tmp')):
            os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro3,'tmp'), f))

        output2 =  opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual_correct.nii.gz')
        input2 = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')

        command = '3dcalc' + overwrite + ' -a ' +   input2 + \
        ' -prefix ' +  output2 + ' -expr "a"'
        spco([command], shell=True)

        command = '3dinfo -same_obl ' + opj(path_anat, ID + 'template_indiv' + TfMRI + '.nii.gz') + ' ' +  output2
        nx = spco([command], shell=True)

        if str(nx)=="b'1\\n1\\n'":
            print('same oblique')
            # reorient the fiedls according to the json file
            command = '3dWarp' + overwrite + ' -deoblique -wsinc5 -prefix ' +  output2 + \
            ' ' +  output2
            spco([command], shell=True)

        else:
            print('oblique different')
            command = '3drefit' + overwrite + ' -orient ' + orientation_orig + ' ' +  output2
            spco([command], shell=True)
            # reorient the fiedls according to the json file
            command = '3dWarp' + overwrite + ' -deoblique -wsinc5 -prefix ' +  output2 + \
            ' ' +  opj(dir_fMRI_Refth_RS_prepro1,root_RS + imgI)
            spco([command], shell=True)
            command = '3drefit' + overwrite + ' -orient ' + orientation + ' ' +  output2
            spco([command], shell=True)

        ##### apply the recenter fmri
        command = '3dZeropad -I 200 -S 200 -A 200 -P 200 -L 200 -R 200 -S 200 -prefix ' +  output2 + ' ' +  output2 + ' -overwrite'
        spco([command], shell=True)
        command = '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
        ' -prefix ' + output2 + \
        ' -input  ' + output2
        spco([command], shell=True)
        command = '3dAutobox' + overwrite + ' -input ' + output2 + ' -prefix ' + output2 + ' -noclust  -overwrite'
        spco(command, shell=True)

        ## apply on pre-processed imgs
        command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'tmp','dummy.nii.gz') + ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual_correct.nii.gz')
        spco([command], shell=True)

        #####transfo originally linearXXXX (line 157) should be linear too here????
        WW = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro3,'tmp','*.nii.gz')))
        for j in range(0, len(WW)):
            if j < 10:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n ' + n_for_ANTS + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.00' + str(j) + '.nii.gz') + \
                mvt_shft_ANTs + \
                transfo_concat_inv
            elif j < 100 and j >= 10:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n ' + n_for_ANTS + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.0' + str(j) + '.nii.gz') + \
                mvt_shft_ANTs + \
                transfo_concat_inv
            else:
                command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n ' + n_for_ANTS + \
                ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
                ' -o ' + opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.' + str(j) + '.nii.gz') + \
                mvt_shft_ANTs + \
                transfo_concat_inv
            
            spco([command], shell=True)

        list_name = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.*.nii.gz')))
        name = ''
        for m in range(0,len(list_name)):
            name = name + ' ' + list_name[m]

        command = '3dTcat' + overwrite + ' ' + name + ' -tr ' + str(TR) + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
        spco([command], shell=True)

    #### apply to all atlases
    for atlas in list_atlases:
        ## in anat space resample to func 
        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, opb(atlas)) + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
        ' -input  ' + atlas
        spco([command], shell=True)

    command = '3dresample' + overwrite + \
    ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'mask_brain.nii.gz') + \
    ' -master ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
    ' -input  ' + BASE_SS_mask
    spco([command], shell=True)

    if creat_sutdy_template== True:
        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
        ' -input  ' + GM_mask_studyT
        spco([command], shell=True)
    else:
        command = '3dresample' + overwrite + \
        ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz') + \
        ' -master ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
        ' -input  ' + GM_mask
        spco([command], shell=True)


