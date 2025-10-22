import os
import ants
import math
import json

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from Tools.extract_filename import extract_filename

def create(input_for_msk,process_dir,masking_img,
        type_of_transform,aff_metric_ants,BASE_SS_coregistr,BASE_SS_mask,
        Tversion,BET,fvalue,
        sing_fsl,diary_file):

    nl = 'INFO: type_of_transform is ' + str(type_of_transform)
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    hd_IMG   = ants.image_header_info(input_for_msk)
    IMG      = ants.image_read(input_for_msk)
    REF_BET  = ants.image_read(BASE_SS_coregistr)
    REF_MASK = ants.image_read(BASE_SS_mask)

    if BET == 1:
        if fvalue == '' :
            fvalue = 0.40
        output = opj(process_dir, extract_filename(input_for_msk) + '_bet' + masking_img)
        command = (sing_fsl + 'bet2 ' + input_for_msk +
                   ' ' + output +
                   ' -f ' + str(fvalue) + ' -c' +
                   ' ' + str(math.ceil(int(hd_IMG['dimensions'][0]) / 2)) +
                   ' ' + str(math.ceil(int(hd_IMG['dimensions'][1]) / 2)) +
                   ' ' + str(math.ceil(int(hd_IMG['dimensions'][2]) / 2)) + ' -m')
        run_cmd.run(command, diary_file)
        refIMG = ants.image_read(output + '.nii.gz')
    else:
        refIMG = IMG



    nameT = opj(process_dir, extract_filename(input_for_msk) + 'template_to_' + masking_img + '_SyN_')

    mTx = ants.registration(fixed=refIMG, moving=REF_BET, type_of_transform='Translation',
                      outprefix=nameT + 'shift_')

    # MEAN_tr = ants.apply_transforms(fixed=refIMG, moving=REF_BET, transformlist=nameT + 'shift_0GenericAffine.mat',
    #                                 interpolator='hammingWindowedSinc')
    # ants.image_write(MEAN_tr, nameT + 'shift.nii.gz', ri=False)

    if Tversion == 'translation':
        transfo_concat = [nameT + 'shift_0GenericAffine.mat']
        run_cmd.msg(mTx, diary_file, 'ENDC')

    elif Tversion == 'nonlinear':
        mTx = ants.registration(fixed=refIMG, moving=REF_BET, outprefix=nameT,
                          type_of_transform=type_of_transform,
                          initial_transform=nameT + 'shift_0GenericAffine.mat',
                          aff_metric=aff_metric_ants,
                          grad_step=0.1,
                          flow_sigma=3,
                          total_sigma=0,
                          aff_sampling=32,
                          aff_random_sampling_rate=0.2,
                          syn_sampling=32,
                          aff_iterations=(1000, 500, 250, 100),
                          aff_shrink_factors=(8, 4, 2, 1),
                          aff_smoothing_sigmas=(3, 2, 1, 0),
                          reg_iterations=(1000, 500, 250, 100),
                          reg_smoothing_sigmas=(3, 2, 1, 0),
                          reg_shrink_factors=(8, 4, 2, 1),
                          verbose=True)

        transfo_concat = [nameT + '1Warp.nii.gz', nameT + '0GenericAffine.mat']
        run_cmd.msg(mTx, diary_file, 'ENDC')

    tmp_mask1 = ants.apply_transforms(fixed=refIMG, moving=REF_MASK,
                                      transformlist=transfo_concat, interpolator='genericLabel')

    return tmp_mask1


def clean1(input_for_msk,output_for_mask,tmp_mask1,BASE_SS_coregistr,BASE_SS_mask):

    IMG = ants.image_read(input_for_msk)

    tmp_mask1 = ants.threshold_image(tmp_mask1, 0.5, 1, 1, 0, True)
    tmp_mask1 = ants.morphology(tmp_mask1, operation='dilate', radius=2, mtype='binary', shape='ball')
    tmp_mask1 = ants.iMath(tmp_mask1, operation='GetLargestComponent')
    seg_tmp = ants.atropos(a=IMG, m='[0.1,1x1x1]', c='[3,0]', i='kmeans[3]', x=tmp_mask1)

    #  Clean up:
    seg_tmp2 = ants.iMath(seg_tmp['segmentation'], 'Pad', 10)
    W = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
    W = ants.iMath(W, operation='GetLargestComponent')
    W = W * 3
    G = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
    G = ants.iMath(G, operation='GetLargestComponent')
    TMP1 = ants.iMath(G, 'FillHoles', 2)
    G = G * TMP1
    C = ants.threshold_image(seg_tmp2, 1, 1, 1, 0, True)
    TMP2 = ants.morphology(C, operation='erode', radius=10, mtype='binary', shape='ball')
    G[G == 0] = TMP2[G == 0]
    G = G * 2
    seg_tmp2 = W
    seg_tmp2[W == 0] = G[W == 0]

    #  clean the Brainmask
    tmp_mask3 = ants.threshold_image(seg_tmp2, 3, 3, 1, 0, True)
    TMP3 = ants.threshold_image(seg_tmp2, 2, 2, 1, 0, True)
    tmp_mask3[tmp_mask3 == 0] = TMP3[tmp_mask3 == 0]
    tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=2, mtype='binary', shape='ball')
    tmp_mask3 = ants.iMath(tmp_mask3, operation='GetLargestComponent')
    tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=4, mtype='binary', shape='ball')
    tmp_mask3 = ants.iMath(tmp_mask3, 'FillHoles', 2)
    tmp_mask3 = ants.morphology(tmp_mask3, operation='dilate', radius=5, mtype='binary', shape='ball')
    tmp_mask3 = ants.morphology(tmp_mask3, operation='erode', radius=5, mtype='binary', shape='ball')
    tmp_mask3 = ants.iMath(tmp_mask3, 'Pad', -10)

    ants.image_write(tmp_mask3, output_for_mask, ri=False)

    dictionary = {"Sources": [input_for_msk,
                              BASE_SS_coregistr,
                              BASE_SS_mask],
                  "Description": 'Brain mask (bet2 from FSL and atropos from ANTspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)

def clean2(input_for_msk,output_for_mask,tmp_mask1,BASE_SS_coregistr,BASE_SS_mask,option2,sing_afni,diary_file):

    spacing = tmp_mask1.spacing  # This will give you the voxel size in x, y, z (e.g., (1.0, 1.0, 1.2) mm)

    # Use voxel size as sigma for Gaussian smoothing
    sigma = spacing  # Set sigma to voxel size for each dimension

    # Apply Gaussian smoothing (sigma controls the amount of smoothing)
    smoothed_mask = ants.smooth_image(tmp_mask1, sigma=sigma)  # Adjust sigma for more or less smoothing

    # Threshold to return to binary mask
    binary_smoothed_mask = ants.threshold_image(smoothed_mask, low_thresh=0.5, high_thresh=1)
    binary_smoothed_mask = ants.iMath(binary_smoothed_mask, operation='GetLargestComponent')
    ants.image_write(binary_smoothed_mask, output_for_mask, ri=False)

    command = (sing_afni + '3dmask_tool -overwrite -prefix ' + output_for_mask +
               ' -input ' + output_for_mask + option2)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [input_for_msk,
                              BASE_SS_coregistr,
                              BASE_SS_mask],
                  "Description": 'Brain mask (Atropos from ANTspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)

def clean3(input_for_msk,output_for_mask,tmp_mask1,sing_afni,diary_file):

    spacing = tmp_mask1.spacing  # This will give you the voxel size in x, y, z (e.g., (1.0, 1.0, 1.2) mm)

    # Use voxel size as sigma for Gaussian smoothing
    sigma = spacing  # Set sigma to voxel size for each dimension

    # Apply Gaussian smoothing (sigma controls the amount of smoothing)
    smoothed_mask = ants.smooth_image(tmp_mask1, sigma=sigma)  # Adjust sigma for more or less smoothing

    # Threshold to return to binary mask
    binary_smoothed_mask = ants.threshold_image(smoothed_mask, low_thresh=0.5, high_thresh=1)
    binary_smoothed_mask = ants.iMath(binary_smoothed_mask, operation='GetLargestComponent')
    ants.image_write(binary_smoothed_mask, output_for_mask, ri=False)

    # Resample the mask image
    command = (sing_afni + '3dresample -master ' + input_for_msk + ' -prefix ' + output_for_mask +
               ' -input ' + output_for_mask + ' -overwrite -bound_type SLAB')
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": input_for_msk,
                  "Description": 'Brain mask (compute_epi_mask from Nilearn).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(output_for_mask[:-7] + '.json', "w") as outfile:
        outfile.write(json_object)