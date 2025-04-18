import os
import subprocess
from fonctions.extract_filename import extract_filename
import ants
import datetime
import json
from fonctions.plot_QC_func import plot_qc
import nibabel as nib
import re
import numpy as np
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


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
                             nb_run, RS, transfo_concat_Anat, w2inv_Anat,do_anat_to_func, list_atlases,
                             BASE_SS_mask, GM_mask, GM_mask_studyT, creat_study_template,anat_func_same_space,
                             orientation, REF_int, IhaveanANAT, overwrite,s_bind,afni_sif,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(9) + '(function: _9_coregistration_to_template_space).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    if IhaveanANAT == False:
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz') + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, root_RS + '_residual_in_template.nii.gz') + ' -expr "a"'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
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
                opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz'),
                opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat')], [False, False]):
                if ope(elem1):
                    mvt_shft_ANTs.append(elem1)
                    w2inv_fwd.append(elem2)
        elif do_anat_to_func == False and anat_func_same_space == True:
            mvt_shft_ANTs = []
            w2inv_fwd     = []
        else :
            nl = 'ERROR: If Anat and Func are not in the same space you need to perform that transformation (do_anat_to_func = True)'
            diary.write(f'\n{nl}')
            raise Exception(bcolors.FAIL + nl + bcolors.ENDC)


    ##### create new variable for template space
    if ope(dir_fMRI_Refth_RS_prepro3) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro3)

        ############################### ############################### ###############################
        ############################### apply transfo to anat space to Mean_Image image for test ######
        ############################### ############################### ###############################

    for input2, output2, output3 in zip([opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.nii.gz')],
                               [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_RcT_SS_pre.nii.gz'), opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution_pre.nii.gz')],
                                [opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz'), opj(dir_fMRI_Refth_RS_prepro3, 'Ref_anat_in_fMRI_anat_resolution_test_in_template.nii.gz')]):

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + input2 + \
        ' -prefix ' + output2 + ' -expr "a"'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        ##### apply the recenter fmri
        if deoblique == 'header':
            command = 'singularity run' + s_bind + afni_sif + '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'change header orientation + deoblique (3drefit, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'WARP' or deoblique == 'WARP_without_3drefit':
            # reorient the fields according to the json file
            command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + \
            ' ' +  output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'reorientation + deoblique (3dWarp, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'WARP_Gridset':  # do nothing
            # reorient the fiedls according to the json file
            command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + ' -gridset ' + output2 + \
                      ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'reorientation + deoblique (3dWarp, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'header_WO_deob':
            command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'change header orientation (3drefit, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'no_deoblique':  # do nothing
            nl = 'nothing done here'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            dictionary = {"Sources": input2,
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        else:
            nl = 'ERROR: wrong deoblique name'
            diary.write(f'\n{nl}')
            dictionary = {"Sources": input2,
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)
            raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

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

            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I ' +  str(padding_voxels[0]) + \
                      ' -S ' + str(padding_voxels[0]) + \
                      ' -A ' + str(padding_voxels[1]) + \
                      ' -P ' + str(padding_voxels[1]) + \
                      ' -L ' + str(padding_voxels[2]) + \
                      ' -R ' + str(padding_voxels[2]) + \
                      ' -prefix ' + output2 + ' ' + output2 + ' -overwrite'
            nl = spgo(command)
            print(nl)
            diary.write(f'\n{nl}')
            print(nl)


            command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' + \
                      opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
                      ' -prefix ' + output2 + ' -input  ' + output2

            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

        ## test on mean img (to see spatially that is works)
        nl = "starting func to template space on mean image and anat test transformation"
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')

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
        print(bcolors.OKGREEN + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')
        # Freeing memory
        del MEAN
        del TRANS

    ### plot the QC
    bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))

    if not os.path.exists(opj(bids_dir, 'QC','meanIMG_in_template')):
        os.mkdir(opj(bids_dir, 'QC','meanIMG_in_template'))

    # Extract ID
    sub_path = os.path.normpath(dir_fMRI_Refth_RS_prepro3).split(os.sep)
    ID = [segment.split('-')[1] for segment in sub_path if segment.startswith('sub-')][0]

    plot_qc(opj(dir_fMRI_Refth_RS_prepro3,'BASE_SS_fMRI.nii.gz'),
            opj(dir_fMRI_Refth_RS_prepro3, 'Mean_Image_RcT_SS_in_template.nii.gz'),
            opj(bids_dir, 'QC','meanIMG_in_template', ID + 'meanIMG_in_template.png'))

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ##                                          Work on all FUNC                                                ## ##
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    for i in range(0, int(nb_run)):
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
                diary.write(f'\n{nl}')
                raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + input2 + \
                  ' -prefix ' + output2 + ' -expr "a"'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        ##### apply the recenter fmri
        #######!!!!!! I don't know why if deoblique1=WARP no gridset if exception2 require gridset!!!!!
        if deoblique == 'header':
            command = 'singularity run' + s_bind + afni_sif + '3drefit -deoblique ' + overwrite + ' -orient ' + orientation + ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'change header orientation + deoblique (3drefit, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'deob_WO_orient':
            command = 'singularity run' + s_bind + afni_sif + '3drefit -deoblique ' + overwrite + ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'change header orientation + deoblique (3drefit, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'header_WO_deob':
            command = 'singularity run' + s_bind + afni_sif + '3drefit ' + overwrite + ' -orient ' + orientation + ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'change header orientation (3drefit, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'WARP' or deoblique == 'WARP_without_3drefit':
            # reorient the fields according to the json file
            command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix ' + output2 + \
                      ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'reorientation + deoblique (3dWarp, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'WARP_Gridset':  # do nothing
            # reorient the fiedls according to the json file
            command = 'singularity run' + s_bind + afni_sif + '3dWarp' + overwrite + ' -deoblique -NN -prefix '\
                      + output2 + ' -gridset ' + output2 + ' ' + output2
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": input2,
                          "Description": 'reorientation + deoblique (3dWarp, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        elif deoblique == 'no_deoblique':  # do nothing
            nl = 'nothing done here'
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            dictionary = {"Sources": input2,
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)
        else:
            nl = 'ERROR: wrong deoblique name'
            diary.write(f'\n{nl}')
            dictionary = {"Sources": input2,
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

            raise Exception(bcolors.FAIL + nl + bcolors.ENDC)


        if anat_func_same_space == True:
            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -I ' + str(padding_voxels[0]) + \
                      ' -S ' + str(padding_voxels[0]) + \
                      ' -A ' + str(padding_voxels[1]) + \
                      ' -P ' + str(padding_voxels[1]) + \
                      ' -L ' + str(padding_voxels[2]) + \
                      ' -R ' + str(padding_voxels[2]) + \
                      ' -prefix ' + output3 + ' ' + output2 + ' -overwrite'
            nl = spgo(command)
            print(nl)
            diary.write(f'\n{nl}')
            print(nl)

            print(f"Applied padding (voxels): X={padding_voxels[0]}, Y={padding_voxels[1]}, Z={padding_voxels[2]}")

            ### delet zeropad and add master mean image in 3D Allineate?
            command = 'singularity run' + s_bind + afni_sif + '3dAllineate' + overwrite + ' -overwrite -interp NN -1Dmatrix_apply ' + \
                      opj(dir_fMRI_Refth_RS_prepro2, '_brain_for_Align_Center.1D') + \
                      ' -prefix ' + output3 + \
                      ' -input  ' + output3
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            ## apply on pre-processed imgs
            FUNC = ants.image_read(output3)
            img_name = output2
        else:
            ## apply on pre-processed imgs
            FUNC = ants.image_read(input2)
            img_name = input2

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

    def get_orientation_nibabel(nifti_path):
        """Get 3-letter orientation code using NiBabel."""
        img = nib.load(nifti_path)
        aff = img.affine

        # Extract orientation from affine matrix
        ornt = nib.orientations.io_orientation(aff)
        codes = nib.orientations.ornt2axcodes(ornt)
        orient_code = ''.join(codes)

        # Validate (should already be valid from NiBabel)
        if not re.fullmatch(r'^[RLAPSI]{3}$', orient_code):
            raise ValueError(f"Invalid orientation: {orient_code}")
        return orient_code

    # Usage
    orient_meanimg = get_orientation_nibabel(residual_in_template)
    print(f"Orientation: {orient_meanimg}")

    #### apply to every atlas
    if len(list_atlases) > 0:
        for atlas in list_atlases:
            command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
                      ' -orient ' + orient_meanimg + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, opb(atlas)) + \
                      ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
                      ' -input  ' + atlas
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": [atlas,
                                      residual_in_template],
                          "Description": ' Resampling (3dresample, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro3, opb(atlas)[:-7] + '.json'), "w") as outfile:
                outfile.write(json_object)
    else:
        nl = 'WARNING: list_atlases is empty!'
        print(bcolors.WARNING + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')

    #### apply to every mask (brain and Gray matter)

    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
              ' -orient ' + orient_meanimg + \
              ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3, 'mask_brain.nii.gz') + \
              ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
              ' -input  ' + BASE_SS_mask
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
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

    command = 'singularity run' + s_bind + afni_sif + '3dresample' + overwrite + \
              ' -orient ' + orient_meanimg + \
              ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro3,'Gmask.nii.gz') + \
              ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' + \
              ' -input  ' + MASK
    nl = spgo(command)
    diary.write(f'\n{nl}')
    print(nl)
    dictionary = {"Sources": [MASK,
                              residual_in_template],
                  "Description": ' Resampling (3dresample, AFNI).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro3, 'Gmask.json'), "w") as outfile:
        outfile.write(json_object)

    diary.write(f'\n')
    diary.close()


