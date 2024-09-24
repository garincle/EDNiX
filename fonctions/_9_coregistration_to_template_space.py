import os
import subprocess
import glob
from fonctions.extract_filename import extract_filename
import ants

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
def to_common_template_space(deoblique, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3,
                nb_run, RS, transfo_concat_Anat, w2inv_Anat,do_anat_to_func, n_for_ANTS, list_atlases, TfMRI, BASE_SS_mask, GM_mask, GM_mask_studyT, creat_study_template,
                anat_func_same_space, orientation, dir_prepro, ID, REF_int, IhaveanANAT, overwrite,s_bind,afni_sif):

    if IhaveanANAT == False:
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz') + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz') + ' -expr "a"'
            spco([command], shell=True)
    else:
        if do_anat_to_func == True:
            mvt_shft_ANTs = []
            w2inv_fwd = []
            for elem1, elem2 in zip([  # opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift_0GenericAffine.mat'),
                opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz'),
                opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat')], [False, False]):
                if ope(elem1):
                    mvt_shft_ANTs.append(elem1)
                    w2inv_fwd.append(elem2)
        elif do_anat_to_func == False and anat_func_same_space == True:
            mvt_shft_ANTs = []
            w2inv_fwd = []
        else : print('ERROR: If Anat and Func are not in the same space you need to perform that trasnformation (do_anat_to_func = True)')

    ##### creat new variable  for template space
    if ope(dir_fMRI_Refth_RS_prepro3) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro3)
        os.makedirs(opj(dir_fMRI_Refth_RS_prepro3,'tmp'))

        ###############################################################################################
        ################################### registration to anat space ################################
        ###############################################################################################
        ############################### ############################### ###############################
        ############################### apply transfo to anat space to Mean_Image image for test ######
        ############################### ############################### ###############################

    for input2, output2, output3 in zip([opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI.nii.gz')],
                               [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_test_pre.nii.gz')],
                                [opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz'), opj(dir_fMRI_Refth_RS_prepro3, 'Ref_anat_in_fMRI_test_in_template.nii.gz')]):

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + input2 + \
        ' -prefix ' + output2 + ' -expr "a"'
        spco([command], shell=True)

        if anat_func_same_space == True:
            r = REF_int
            root_RS = extract_filename(RS[r])
            command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz') + ' ' +  output2
            orientation_orig = spgo(command).split('\n')[-1]
            print(orientation_orig)
            ##### apply the recenter fmri
            #######!!!!!! I don't know why if deobblique1=WARP no gridset if exeption2 require gridset!!!!!
            if deoblique == 'header':
                command = 'singularity run' + s_bind + afni_sif + '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2
                spco([command], shell=True)

            elif deoblique == 'WARP' or deoblique == 'WARP_without_3drefit':
                # reorient the fiedls according to the json file
                command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + \
                ' ' +  output2
                spco([command], shell=True)

            elif deoblique == 'WARP_Gridset':  # do nothing
                print('WARP_Gridset')
                # reorient the fiedls according to the json file
                command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + ' -gridset ' + output2 + \
                          ' ' + output2
                spco([command], shell=True)

            elif deoblique == 'header_WO_deob':
                print('header_WO_deob')
                command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' + output2
                spco([command], shell=True)

            elif deoblique == 'no_deoblique':  # do nothing
                print('no_deoblique')

            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 50 -S 50 -A 50 -P 50 -L 50 -R 50 -S 50 -prefix ' +  output2 + ' ' +  output2 + ' -overwrite'
            spco([command], shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
            ' -prefix ' + output2 + \
            ' -input  ' + output2
            #' -master ' + opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
            spco([command], shell=True)
            '''
            command = '3dinfo -di ' + input2
            delta_x = str(abs(round(float(spgo([command])[-8:]), 10)))
            command = '3dinfo -dj ' + input2
            delta_y = str(abs(round(float(spgo([command])[-8:]), 10)))
            command = '3dinfo -dk ' + input2
            delta_z = str(abs(round(float(spgo([command])[-8:]), 10)))

            command = '3dresample' + overwrite + \
                      ' -prefix ' + output2 + \
                      ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                      ' -input  ' + output2
            spco([command], shell=True)
        '''
        else:
            print('nothing to do')

        ## test on mean img (to see spatially that is works)
        print("starting est on mean image and anat test to standart space")
        MEAN = ants.image_read(output2)
        REF  = ants.image_read(opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'))
        print(mvt_shft_ANTs + transfo_concat_Anat)
        print(w2inv_Anat + w2inv_fwd)
        TRANS = ants.apply_transforms(fixed=REF, moving=MEAN,
                                      transformlist=mvt_shft_ANTs + transfo_concat_Anat,
                                      interpolator='nearestNeighbor',
                                      whichtoinvert=w2inv_fwd + w2inv_Anat)
        ants.image_write(TRANS, output3, ri=False)
        print(str(output3) + ' done!')

            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
            ## ## ## ## ## ## ## ## ## ## ## ## ## ##  Work on all FUNC ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    for i in range(0, int(nb_run)):
        ##### go for BOLD img preTTT
        root_RS = extract_filename(RS[i])
        for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro3, 'tmp')):
            if ope(os.path.join(opj(dir_fMRI_Refth_RS_prepro3, 'tmp'), f)):
                os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro3, 'tmp'), f))

        output3 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_anat_sfht.nii.gz')
        output2 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_anat_reorient.nii.gz')
        input2  = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')

        if anat_func_same_space == True:
            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + input2 + \
                      ' -prefix ' + output2 + ' -expr "a"'
            spco([command], shell=True)

            command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -same_obl ' + opj(dir_prepro, ID + 'template_indiv' + TfMRI + '.nii.gz') + ' ' +  output2
            nx = spgo(command).split('\n')[-1]

            r = REF_int
            root_RS = extract_filename(RS[r])
            command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz') + ' ' +  output2
            orientation_orig = spgo(command).split('\n')[-1]
            print(orientation_orig)

            ##### apply the recenter fmri
            #######!!!!!! I don't know why if deobblique1=WARP no gridset if exeption2 require gridset!!!!!
            if deoblique == 'header':
                command = 'singularity run' + s_bind + afni_sif + '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2
                spco([command], shell=True)

            elif deoblique == 'WARP' or deoblique == 'WARP_without_3drefit':
                # reorient the fiedls according to the json file
                command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + \
                          ' ' + output2
                spco([command], shell=True)

            elif deoblique == 'WARP_Gridset':  # do nothing
                print('WARP_Gridset')
                # reorient the fiedls according to the json file
                command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + ' -gridset ' + output2 + \
                          ' ' + output2
                spco([command], shell=True)

            elif deoblique == 'header_WO_deob':
                print('header_WO_deob')
                command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' + output2
                spco([command], shell=True)

            elif deoblique == 'no_deoblique':  # do nothing
                print('no_deoblique')

            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I 50 -S 50 -A 50 -P 50 -L 50 -R 50 -S 50 -prefix ' + output3 + ' ' + output2 + ' -overwrite'
            spco([command], shell=True)

            command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' + opj(
                dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
                      ' -prefix ' + output3 + \
                      ' -input  ' + output3
                      # ' -master ' + opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
            spco([command], shell=True)
            ## apply on pre-processed imgs
            FUNC = ants.image_read(output3)

        else:
            ## apply on pre-processed imgs
            FUNC = ants.image_read(input2)
        #####transfo

        TRANS = ants.apply_transforms(fixed=REF, moving=FUNC,
                                      transformlist=transfo_concat_Anat + mvt_shft_ANTs,
                                      interpolator='nearestNeighbor',
                                      whichtoinvert=w2inv_Anat + w2inv_fwd,imagetype=3)
        ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz'),
                         ri=False)


    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
    delta_x = str(abs(round(float(spgo(command).split('\n')[-1]), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
    delta_y = str(abs(round(float(spgo(command).split('\n')[-1]), 10)))
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
    delta_z = str(abs(round(float(spgo(command).split('\n')[-1]), 10)))
    ## in anat space resample to func
    command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
    orient_meanimg = spgo(command).split('\n')[-1]

    #### apply to all atlases
    if len(list_atlases) > 0:
        for atlas in list_atlases:
            command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                      ' -orient ' + orient_meanimg + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, opb(atlas)) + \
                      ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                      ' -input  ' + atlas
            spco([command], shell=True)

        command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                  ' -orient ' + orient_meanimg + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'mask_brain.nii.gz') + \
                  ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                  ' -input  ' + BASE_SS_mask
        spco([command], shell=True)
    else:
        print('WARNING: list_atlases is empty!')

    if creat_study_template== True:
        command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                  ' -orient ' + orient_meanimg + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz') + \
                  ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                  ' -input  ' + GM_mask_studyT
        spco([command], shell=True)
    else:
        command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                  ' -orient ' + orient_meanimg + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz') + \
                  ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                  ' -input  ' + GM_mask
        spco([command], shell=True)


