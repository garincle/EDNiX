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
from fonctions.extract_filename import extract_filename


####################################################################################
########################## Step 3 normalisation to template atlas space ############
####################################################################################
def to_common_template_space(Session, deoblique_exeption1, deoblique_exeption2, deoblique, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, BASE_SS_coregistr,
                nb_run, RS, transfo_concat_inv, TR, n_for_ANTS, list_atlases, TfMRI, BASE_SS_mask, GM_mask, GM_mask_studyT, creat_sutdy_template, 
                anat_func_same_space, orientation, path_anat, ID, REF_int, dir_prepro, IhaveanANAT, overwrite):


    if IhaveanANAT == False:
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz') + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz') + ' -expr "a"'
            spco([command], shell=True)

    else:
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
            ############################### apply transfo to anat space to Mean_Image image for test ######
            ############################### ############################### ###############################
        for input2, output2 in zip([opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')],
                                   [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz')]):
            command = '3dcalc' + overwrite + ' -a ' + input2 + \
                      ' -prefix ' + output2 + ' -expr "a"'
            spco([command], shell=True)

            if anat_func_same_space == True:
                debolique_spe = ID + 'ses-' + str(Session)
                if debolique_spe in deoblique_exeption1:
                    deoblique_1 = 'exeption1'
                    print(deoblique_1)
                elif debolique_spe in deoblique_exeption2:
                    deoblique_1 = 'exeption2'
                else:
                    deoblique_1 = deoblique
                    print(deoblique_1)

                command = '3dinfo -same_obl ' + opj(path_anat, ID + 'template_indiv' + TfMRI + '.nii.gz') + ' ' +  output2
                nx = spco([command], shell=True)

                r = REF_int
                root_RS = extract_filename(RS[r])
                command = '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz') + ' ' +  output2
                orientation_orig = spgo([command])
                print(orientation_orig)
                ##### apply the recenter fmri
                #######!!!!!! I don't know why if deobblique1=WARP no gridset if exeption2 require gridset!!!!!
                if deoblique_1 == 'header':
                    command = '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2
                    spco([command], shell=True)

                elif deoblique_1 == 'WARP':
                    # reorient the fiedls according to the json file
                    command = '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + \
                    ' ' +  output2
                    spco([command], shell=True)

                elif deoblique_1 == 'exeption1':  # do nothing
                    print('exeption1: nothing to do hear')

                elif deoblique_1 == 'exeption2':  # re-alineate
                    # reorient the fiedls according to the json file
                    command = '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + ' -gridset ' + output2 + \
                              ' ' + output2
                    spco([command], shell=True)

                command = '3dZeropad -I 50 -S 50 -A 50 -P 50 -L 50 -R 50 -S 50 -prefix ' +  output2 + ' ' +  output2 + ' -overwrite'
                spco([command], shell=True)

                command = '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
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
                print("nothing to do")


            ## test on mean img (to see spatially that is works)
            command = 'antsApplyTransforms -d 3 -e 3 -i ' + opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz') + ' -n ' + n_for_ANTS + \
            ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
            ' -o ' + opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz') + \
            transfo_concat_inv  + \
            mvt_shft_ANTs
            spco([command], shell=True)


            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
            ## ## ## ## ## ## ## ## ## ## ## ## ## ##  Work on all FUNC ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


            for i in range(0, int(nb_run)):
                ##### go for BOLD img preTTT
                root_RS = extract_filename(RS[i])
                for f in os.listdir(opj(dir_fMRI_Refth_RS_prepro3, 'tmp')):
                    os.remove(os.path.join(opj(dir_fMRI_Refth_RS_prepro3, 'tmp'), f))

                output3 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_anat_sfht.nii.gz')
                output2 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_anat_reorient.nii.gz')
                input2 = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')

                if anat_func_same_space == True:
                    command = '3dcalc' + overwrite + ' -a ' + input2 + \
                              ' -prefix ' + output2 + ' -expr "a"'
                    spco([command], shell=True)

                    command = '3dinfo -same_obl ' + opj(path_anat, ID + 'template_indiv' + TfMRI + '.nii.gz') + ' ' +  output2
                    nx = spco([command], shell=True)

                    r = REF_int
                    root_RS = extract_filename(RS[r])
                    command = '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_mean_preWARP' + '.nii.gz') + ' ' +  output2
                    orientation_orig = spgo([command])
                    print(orientation_orig)

                    ##### apply the recenter fmri
                    #######!!!!!! I don't know why if deobblique1=WARP no gridset if exeption2 require gridset!!!!!
                    if deoblique_1 == 'header':
                        command = '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2
                        spco([command], shell=True)

                    elif deoblique_1 == 'WARP':
                        # reorient the fiedls according to the json file
                        command = '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + \
                                  ' ' + output2
                        spco([command], shell=True)

                    elif deoblique_1 == 'exeption1':  # do nothing
                        print('exeption1: nothing to do hear')

                    elif deoblique_1 == 'exeption2':  # re-alineate
                        # reorient the fiedls according to the json file
                        command = '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + ' -gridset ' + output2 + \
                                  ' ' + output2
                        spco([command], shell=True)

                    command = '3dZeropad -I 50 -S 50 -A 50 -P 50 -L 50 -R 50 -S 50 -prefix ' + output3 + ' ' + output2 + ' -overwrite'
                    spco([command], shell=True)

                    command = '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' + opj(
                        dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
                              ' -prefix ' + output3 + \
                              ' -input  ' + output3
                              # ' -master ' + opj(dir_prepro, ID + '_acpc_cropped' + TfMRI + '.nii.gz')
                    spco([command], shell=True)

                    ## apply on pre-processed imgs
                    command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, 'tmp', 'dummy.nii.gz') + ' ' + \
                              output3
                    spco([command], shell=True)

                    '''
                    except:
                        ##### apply the recenter fmri
                        command = '3dAllineate' + overwrite + ' -overwrite -1Dmatrix_apply ' + opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
                        ' -prefix ' + output3 + \
                        ' -input  ' + output2 + \
                        ' -master ' + opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz')
                        spco([command], shell=True)
                    
                        ## apply on pre-processed imgs
                        command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, 'tmp', 'dummy.nii.gz') + ' ' + \
                                  output3
                        spco([command], shell=True)
                    '''
                else:
                    ## apply on pre-processed imgs
                    command = '3dTsplit4D' + overwrite + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, 'tmp',
                                                                           'dummy.nii.gz') + ' ' + input2
                    spco([command], shell=True)

                #####transfo
                WW = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro3,'tmp','*.nii.gz')))
                for j in range(0, len(WW)):
                    if j < 10:
                        command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                        ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
                        ' -o ' + opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.00' + str(j) + '.nii.gz') + \
                        transfo_concat_inv + \
                        mvt_shft_ANTs
                    elif j < 100 and j >= 10:
                        command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                        ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
                        ' -o ' + opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.0' + str(j) + '.nii.gz') + \
                        transfo_concat_inv + \
                        mvt_shft_ANTs
                    else:
                        command = 'antsApplyTransforms -d 3 -e 3 -i ' + WW[j] + ' -n NearestNeighbor' + \
                        ' -r ' + opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz') + \
                        ' -o ' + opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.' + str(j) + '.nii.gz') + \
                        transfo_concat_inv + \
                        mvt_shft_ANTs

                    spco([command], shell=True)

                list_name = sorted(glob.glob(opj(dir_fMRI_Refth_RS_prepro3,'tmp','warp1.*.nii.gz')))
                name = ''
                for m in range(0,len(list_name)):
                    name = name + ' ' + list_name[m]

                command = '3dTcat' + overwrite + ' ' + name + ' -tr ' + str(TR) + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
                spco([command], shell=True)


    command = '3dinfo -di ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
    delta_x = str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dj ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
    delta_y = str(abs(round(float(spgo([command])[-8:]), 10)))
    command = '3dinfo -dk ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')

    delta_z = str(abs(round(float(spgo([command])[-8:]), 10)))
    ## in anat space resample to func
    command = '3dinfo -orient ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
    orient_meanimg = str(spgo([command]))

    #### apply to all atlases
    for atlas in list_atlases:
        command = '3dresample' + overwrite + \
                  ' -orient ' + orient_meanimg + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, opb(atlas)) + \
                  ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                  ' -input  ' + atlas
        spco([command], shell=True)

    command = '3dresample' + overwrite + \
              ' -orient ' + orient_meanimg + \
              ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'mask_brain.nii.gz') + \
              ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
              ' -input  ' + BASE_SS_mask
    spco([command], shell=True)

    if creat_sutdy_template== True:
        command = '3dresample' + overwrite + \
                  ' -orient ' + orient_meanimg + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz') + \
                  ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                  ' -input  ' + GM_mask_studyT
        spco([command], shell=True)
    else:
        command = '3dresample' + overwrite + \
                  ' -orient ' + orient_meanimg + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz') + \
                  ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                  ' -input  ' + GM_mask
        spco([command], shell=True)


