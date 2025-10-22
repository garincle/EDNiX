import os
import ants
import json
import nibabel as nib
import re
import numpy as np

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd,get_orientation
from fonctions import plot_QC_func
from fonctions.extract_filename import extract_filename

####################################################################################
########################## Step 3 normalisation to template atlas space ############
####################################################################################

def to_common_template_space(deoblique, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3,
                             nb_run, RS, transfo_concat_Anat, w2inv_Anat,do_anat_to_func, list_atlases,
                             BASE_SS_mask, GM_mask, GM_mask_studyT, creat_study_template,anat_func_same_space,
                             orientation, template,template_labeldir, IhaveanANAT, overwrite,sing_afni,diary_file):

    nl = '##  Working on step ' + str(9) + '(function: _9_coregistration_to_template_space).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    # Extract ID
    sub_path = opn(dir_fMRI_Refth_RS_prepro3).split(os.sep)
    ID = [segment.split('-')[1] for segment in sub_path if segment.startswith('sub-')][0]

    ref_img   = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')
    anat_func = opj(dir_fMRI_Refth_RS_prepro1, ('_').join([ID, '_space-func', TfMRI + '.nii.gz']))

    if IhaveanANAT == False:
        for i in range(int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = (sing_afni + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2,root_RS + '_residual_in_anat.nii.gz') +
                       ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,root_RS + '_residual_in_template.nii.gz') + ' -expr "a"')
            run_cmd.do(command, diary_file)

            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz'),
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.json'), "w") as outfile:
                outfile.write(json_object)
    else:
        if do_anat_to_func == True:
            mvt_shft_ANTs = []
            w2inv_fwd     = []

            for elem1, elem2 in zip([  # opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift_0GenericAffine.mat'),
                ref_img.replace('.nii.gz', '_unwarped_1Warp.nii.gz'),
                ref_img.replace('.nii.gz', '_0GenericAffine.mat')], [False, False]):
                if ope(elem1):
                    mvt_shft_ANTs.append(elem1)
                    w2inv_fwd.append(elem2)
        elif do_anat_to_func == False and anat_func_same_space == True:
            mvt_shft_ANTs = []
            w2inv_fwd     = []
        else :
            nl = 'ERROR: If Anat and Func are not in the same space you need to perform that transformation (do_anat_to_func = True)'
            raise Exception(run_cmd.error(nl, diary_file))


    ##### create new variable for template space
    if ope(dir_fMRI_Refth_RS_prepro3) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro3)

        ############################### ############################### ###############################
        ############################### apply transfo to anat space to Mean_Image image for test ######
        ############################### ############################### ###############################

    for input2, output2, output3 in zip([ref_img, anat_func],
                               [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution_pre.nii.gz')],
                                [opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz'), opj(dir_fMRI_Refth_RS_prepro3, 'Ref_anat_in_fMRI_anat_resolution_test_in_template.nii.gz')]):

        command = (sing_afni + '3dcalc' + overwrite + ' -a ' + input2 + ' -prefix ' + output2 + ' -expr "a"')
        run_cmd.do(command, diary_file)

        ##### apply the recenter fmri
        if deoblique == 'header':
            command = (sing_afni + '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2)
            run_cmd.run(command, diary_file)
            desc = 'change header orientation + deoblique (3drefit, AFNI).'

        elif deoblique == 'deob_WO_orient':
            command = (sing_afni + '3drefit -deoblique ' + overwrite + ' ' + output2)
            run_cmd.run(command, diary_file)
            desc= 'change header orientation + deoblique (3drefit, AFNI).'

        elif deoblique == 'WARP' or deoblique == 'WARP_without_3drefit' or deoblique == 'WARPbaboon':
            # reorient the fields according to the json file
            command = (sing_afni + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 +
                       ' ' + output2)
            run_cmd.run(command, diary_file)
            desc= 'reorientation + deoblique (3dWarp, AFNI).'

        elif deoblique == 'WARP_Gridset':
            # reorient the fiedls according to the json file
            command = (sing_afni + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 +
                       ' -gridset ' + output2 + ' ' + output2)
            run_cmd.run(command, diary_file)
            desc ='reorientation + deoblique (3dWarp, AFNI).'

        elif deoblique == 'header_WO_deob':
            command = (sing_afni + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' + output2)
            run_cmd.run(command, diary_file)
            desc ='change header orientation (3drefit, AFNI).'

        elif deoblique == 'no_deoblique':  # do nothing
            nl = 'nothing done here'
            run_cmd.msg(nl, diary_file, 'OKGREEN')
            desc = 'Copy.'

        else:
            nl = 'ERROR: wrong deoblique name'
            desc = 'Copy.'
            dictionary = {"Sources": input2,
                          "Description": desc},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)
            raise Exception(run_cmd.error(nl, diary_file))


        dictionary = {"Sources": input2,
                     "Description": desc},
        json_object = json.dumps(dictionary, indent=2)
        with open(output2[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

        '''
        if anat_func_same_space == True:
            
            # --- Step 1: Load and parse the affine matrix ---
            translation_file = opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D')
            affine_params = np.loadtxt(translation_file)  # Load 12-parameter affine

            # Reshape into 3x4 matrix (rotation + translation)
            affine_matrix = affine_params.reshape(3, 4)
            translations_mm = affine_matrix[:, 3]  # Extract [tx, ty, tz] in mm

            # --- Step 2: Calculate worst-case displacement at image edges ---
            ref_img = nib.load(output2)
            voxel_size = ref_img.header.get_zooms()[:3]  # (dx, dy, dz) in mm
            img_shape = ref_img.shape[:3]  # (nx, ny, nz)

            # Convert translations to voxel shifts (account for rotation)
            max_shift_voxels = np.ceil(np.abs(translations_mm) / voxel_size).astype(int)

            # Add safety margin (5 voxels) + edge protection (10% of FOV)
            safety_margin = 5
            edge_margin = (0.1 * np.array(img_shape)).astype(int)  # Extra 10% of image size
            padding_voxels = max_shift_voxels + safety_margin + edge_margin

            direction = ['I', 'S', 'A', 'P', 'L', 'R']
            cmd = []
            for i in range(len(direction)):
                if direction[i] in ['I','S']:
                    vox = padding_voxels[0]
                elif direction[i] in ['A','P']:
                    vox = padding_voxels[1]
                elif direction[i] in ['L', 'R']:
                    vox = padding_voxels[2]

                cmd.append('-' + direction[i])
                cmd.append(str(vox))
            pad = ' '.join(cmd)
            print(f"Applied padding (voxels): X={padding_voxels[0]}, Y={padding_voxels[1]}, Z={padding_voxels[2]}")

            command = (sing_afni + '3dZeropad ' + pad + ' -prefix ' + output2 + ' ' + output2 + ' -overwrite')
            run_cmd.run(command, diary_file)

            command = (sing_afni + '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' +
                       opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') +
                       ' -prefix ' + output2 + ' -input ' + output2)

            run_cmd.run(command, diary_file)
        '''
        ## test on mean img (to see spatially that is works)
        nl = "starting func to template space on mean image and anat test transformation"
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        MEAN = ants.image_read(output2)
        REF  = ants.image_read(opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'))
        TRANS = ants.apply_transforms(fixed=REF, moving=MEAN,
                                      transformlist=transfo_concat_Anat + mvt_shft_ANTs,
                                      interpolator='nearestNeighbor',
                                      whichtoinvert=w2inv_Anat + w2inv_fwd)
        ants.image_write(TRANS, output3, ri=False)
        dictionary = {"Sources": [output2,
                                  opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz')],
                      "Description": ' Non linear normalization (ANTspy).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(output3[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

        nl = str(output3) + ' done!'
        run_cmd.msg(nl, diary_file, 'OKGREEN')

        # Freeing memory
        del MEAN
        del TRANS

    ### plot the QC
    bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))

    if not ope(opj(bids_dir, 'QC','meanIMG_in_template')):
        os.mkdir(opj(bids_dir, 'QC','meanIMG_in_template'))



    plot_QC_func.plot_qc(opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'),
            opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz'),
            opj(bids_dir, 'QC','meanIMG_in_template', ID + 'meanIMG_in_template.png'))

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ##                                          Work on all FUNC                                                ## ##
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    for i in range(int(nb_run)):
        ##### go for BOLD img preTTT

        root_RS = extract_filename(RS[i])
        residual_in_template = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz')
        output3 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_anat_sfht.nii.gz')
        output2 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_anat_reorient.nii.gz')
        input2  = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')

        if ope(input2) == False:
            residual_in_template = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template_failed.nii.gz')
            output3 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_3dDeconvolve_failed_in_anat_sfht.nii.gz')
            output2 = opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_3dDeconvolve_failed_in_anat_reorient.nii.gz')
            input2  = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_3dDeconvolve_failed.nii.gz')

            if ope(input2) == False:
                nl = 'ERROR: corrected func image in func space not found'
                raise Exception(run_cmd.error(nl, diary_file))

        command = (sing_afni + '3dcalc' + overwrite + ' -a ' + input2 + ' -prefix ' + output2 + ' -expr "a"')
        run_cmd.do(command, diary_file)

        ##### apply the recenter fmri
        #######!!!!!! I don't know why if deoblique1=WARP no gridset if exception2 require gridset!!!!!

        if deoblique == 'header':
            command = (sing_afni + '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2)
            run_cmd.run(command, diary_file)
            desc ='change header orientation + deoblique (3drefit, AFNI).'

        elif deoblique == 'deob_WO_orient':
            command = sing_afni + '3drefit -deoblique ' + overwrite + ' ' + output2
            run_cmd.run(command, diary_file)
            desc = 'change header orientation + deoblique (3drefit, AFNI).'

        elif deoblique == 'header_WO_deob':
            command = sing_afni + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' + output2
            run_cmd.run(command, diary_file)
            desc = 'change header orientation (3drefit, AFNI).'

        elif deoblique == 'WARP' or deoblique == 'WARP_without_3drefit' or deoblique == 'WARPbaboon':
            # reorient the fields according to the json file
            command = (sing_afni + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + ' ' + output2)
            run_cmd.run(command, diary_file)
            desc = 'reorientation + deoblique (3dWarp, AFNI).'

        elif deoblique == 'WARP_Gridset':  # do nothing
            # reorient the fiedls according to the json file
            command = (sing_afni + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' +
                       output2 + ' -gridset ' + output2 + ' ' + output2)
            run_cmd.run(command, diary_file)
            desc = 'reorientation + deoblique (3dWarp, AFNI).'

        elif deoblique == 'no_deoblique':  # do nothing
            nl = 'nothing done here'
            run_cmd.msg(nl, diary_file, 'OKGREEN')
            desc = 'Copy.'

        else:
            nl = 'ERROR: wrong deoblique name'
            dictionary = {"Sources": input2,
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)
            raise Exception(run_cmd.error(nl,diary_file))

        dictionary = {"Sources": input2,
                      "Description": desc},
        json_object = json.dumps(dictionary, indent=2)
        with open(output2[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

        '''
        if anat_func_same_space == True:
    
            command = (sing_afni + '3dZeropad ' + pad + ' -prefix ' + output3 + ' ' + output2 + ' -overwrite')
            run_cmd.run(command, diary_file)

            ### delete zeropad and add master mean image in 3D Allineate?
            command = (sing_afni + '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' +
                       opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') +
                       ' -prefix ' + output3 + ' -input ' + output3)
            run_cmd.run(command, diary_file)

            ## apply on pre-processed imgs
            FUNC = ants.image_read(output3)
            img_name = output2
        else:
            ## apply on pre-processed imgs
            FUNC = ants.image_read(input2)
            img_name = input2
        '''

        #  transfo
        TRANS = ants.apply_transforms(fixed=REF, moving=FUNC,
                                      transformlist=transfo_concat_Anat + mvt_shft_ANTs,
                                      interpolator='nearestNeighbor',
                                      whichtoinvert=w2inv_Anat + w2inv_fwd, imagetype=3)
        ants.image_write(TRANS, residual_in_template, ri=False)
        dictionary = {"Sources": [img_name,
                                  opj(dir_fMRI_Refth_RS_prepro3, 'BASE_SS_fMRI.nii.gz')],
                      "Description": ' Non linear normalization (ANTspy).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.json'), "w") as outfile:
            outfile.write(json_object)

    # Load the image directly
    img = nib.load(residual_in_template)
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]

    # Usage
    orient_meanimg = get_orientation.get_orientation_nibabel(residual_in_template)
    run_cmd.msg('Orientation: ' + orient_meanimg, diary_file, 'ENDC')

    #### apply to every atlas
    if len(list_atlases) > 0:
        for atlas in list_atlases:
            atlasfile = '_seg-' + atlas + '_dseg.nii.gz'
            command = (sing_afni + '3dresample' + overwrite + ' -orient ' + orient_meanimg +
                       ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, ID + atlasfile) +
                       ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' +
                       ' -input ' + opj(template_labeldir,template + atlasfile))
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": [atlas,
                                      residual_in_template],
                          "Description": ' Resampling (3dresample, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro3, atlasfile.replace('.nii.gz','.json')), "w") as outfile:
                outfile.write(json_object)
    else:
        nl = 'WARNING: list_atlases is empty!'
        run_cmd.msg(nl, diary_file, 'WARNING')

    #### apply to every mask (brain and Gray matter)

    command = (sing_afni + '3dresample' + overwrite + ' -orient ' + orient_meanimg +
               ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, 'mask_brain.nii.gz') +
               ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' +
               ' -input ' + BASE_SS_mask)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [BASE_SS_mask,
                              residual_in_template],
                  "Description": ' Resampling (3dresample, AFNI).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro3, 'mask_brain.json'), "w") as outfile:
        outfile.write(json_object)

    if creat_study_template== True:
        MASK = GM_mask_studyT
    else:
        MASK = GM_mask

    command = (sing_afni + '3dresample' + overwrite + ' -orient ' + orient_meanimg +
               ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, 'Gmask.nii.gz') +
               ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + ' -input  ' + MASK)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [MASK,
                              residual_in_template],
                  "Description": ' Resampling (3dresample, AFNI).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro3, 'Gmask.json'), "w") as outfile:
        outfile.write(json_object)



