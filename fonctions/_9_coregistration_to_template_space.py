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
opi = os.path.isfile

from Tools import run_cmd,get_orientation
from fonctions import plot_QC_func
from fonctions.extract_filename import extract_filename

####################################################################################
########################## Step 3 normalisation to template atlas space ############
####################################################################################

def to_common_template_space(dir_prepro_template_process, bids_dir, ID, dir_prepro_template_labels, n_for_ANTS,
                            dir_prepro_orig_postprocessed, dir_prepro_acpc_postprocessed, dir_prepro_template_postprocessed,
                             nb_run, RS, do_anat_to_func, list_atlases, info, dir_prepro_orig_process, species,
                             template_dir_labels, template_dir_masks ,anat_func_same_space, dir_prepro_acpc_process,
                             dir_prepro_template_masks, IhaveanANAT, use_erode_WM_func_masks, overwrite,sing_afni,diary_file):

    nl = '##  Working on step ' + str(9) + '(function: _9_coregistration_to_template_space).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    Mean_Image_acpc = opj(dir_prepro_orig_process, 'all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS.nii.gz')
    Mean_Image_template = opj(dir_prepro_template_postprocessed, 'all_runs_space-template-func_desc-fMRI_Mean_Image_SS.nii.gz')
    Mean_Image_unwarped = opj(dir_prepro_acpc_process, 'all_runs_space-anat_desc-fMRI_Mean_Image_unwarped.nii.gz')
    Template_res_func = opj(dir_prepro_template_process, 'BASE_SS_fMRI.nii.gz')
    if IhaveanANAT == False:
        for i in range(int(nb_run)):
            root_RS = extract_filename(RS[i])
            if ope(opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')) == False:
                residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
                residual_anat = opj(dir_prepro_acpc_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')
                residual_template = opj(dir_prepro_template_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')
            else:
                residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
                residual_anat = opj(dir_prepro_acpc_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')
                residual_template = opj(dir_prepro_template_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')

            command = (sing_afni + '3dcalc' + overwrite + ' -a ' + residual_anat +
                       ' -prefix ' + residual_template + ' -expr "a"')
            run_cmd.do(command, diary_file)

            dictionary = {"Sources": residual_anat,
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(residual_template.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)
    else:
        if do_anat_to_func == True:
            mvt_shft_ANTs = []
            w2inv_fwd = []
            for elem1, elem2 in zip([Mean_Image_unwarped.replace('.nii.gz', '_1Warp.nii.gz'),
                                     Mean_Image_unwarped.replace('.nii.gz', '_0GenericAffine.mat')],
                                    [True, False]):
                if opi(elem1):
                    mvt_shft_ANTs.append(elem1)
                    w2inv_fwd.append(elem2)
        elif do_anat_to_func == False and anat_func_same_space == True:
            mvt_shft_ANTs = []
            w2inv_fwd     = []
        else :
            nl = 'ERROR: If Anat and Func are not in the same space you need to perform that transformation (do_anat_to_func = True)'
            raise Exception(run_cmd.error(nl, diary_file))

    ############################### ############################### ###############################
    ############################### apply transfo to anat space to Mean_Image image for test ######
    ############################### ############################### ###############################

    ## test on mean img (to see spatially that is works)
    nl = "starting func to template space on mean image and anat test transformation"
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    MEAN = ants.image_read(Mean_Image_acpc)
    REF  = ants.image_read(Template_res_func)

    TRANS = ants.apply_transforms(fixed=REF, moving=MEAN,
                                  transformlist=info[0][9] + mvt_shft_ANTs,
                                  interpolator=n_for_ANTS,
                                  whichtoinvert=info[0][10] + w2inv_fwd)
    ants.image_write(TRANS, Mean_Image_template, ri=False)

    dictionary = {"Sources": [Mean_Image_acpc,
                              Template_res_func],
                  "Description": ' Non linear normalization (ANTspy).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(Mean_Image_template.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    nl = str(Mean_Image_template) + ' done!'
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    # Freeing memory
    del MEAN
    del TRANS

    if not ope(opj(bids_dir, 'QC','meanIMG_in_template')):
        os.mkdir(opj(bids_dir, 'QC','meanIMG_in_template'))

    plot_QC_func.plot_qc(Template_res_func,
            Mean_Image_template,
            opj(bids_dir, 'QC','meanIMG_in_template', ID + 'meanIMG_in_template.png'))

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ##                                          Work on all FUNC                                                ## ##
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    if IhaveanANAT == True:
        for i in range(int(nb_run)):
            ##### go for BOLD img preTTT

            for i in range(int(nb_run)):
                root_RS = extract_filename(RS[i])
                if ope(opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')) == False:
                    residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
                    residual_anat = opj(dir_prepro_acpc_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')
                    residual_template = opj(dir_prepro_template_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')
                else:
                    residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
                    residual_anat = opj(dir_prepro_acpc_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')
                    residual_template = opj(dir_prepro_template_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')

            FUNC = ants.image_read(residual)
            REF = ants.image_read(Template_res_func)

            TRANS = ants.apply_transforms(fixed=REF, moving=FUNC,
                                          transformlist=info[0][9] + mvt_shft_ANTs,
                                          interpolator=n_for_ANTS,
                                          whichtoinvert=info[0][10] + w2inv_fwd,
                                          imagetype=3)
            ants.image_write(TRANS, residual_template, ri=False)

            dictionary = {"Sources": [residual_template,
                                      Template_res_func],
                          "Description": ' Non linear normalization (ANTspy).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(residual_template.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            nl = str(residual_template) + ' done!'
            run_cmd.msg(nl, diary_file, 'OKGREEN')
            # Freeing memory
            del FUNC
            del REF

    # Load the image directly
    img = nib.load(residual_template)
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]

    # Usage
    orient_meanimg = get_orientation.get_orientation_nibabel(residual_template)
    run_cmd.msg('Orientation: ' + orient_meanimg, diary_file, 'ENDC')

    #### apply to every atlas
    if len(list_atlases[0]) > 0:
        for atlas in list_atlases[0]:
            atlasfile = ID + '_seg-' + atlas + '_dseg.nii.gz'
            command = (sing_afni + '3dresample' + overwrite + ' -orient ' + orient_meanimg +
                       ' -prefix ' + opj(dir_prepro_template_labels, atlasfile) +
                       ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                       ' -input ' + opj(template_dir_labels, species + '_seg-' + atlas + '_dseg.nii.gz'))
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": [opj(template_dir_labels, species + '_seg-' + atlas + '_dseg.nii.gz'),
                                      residual_template],
                          "Description": ' Resampling (3dresample, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro_template_labels, atlasfile).replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)
    else:
        nl = 'WARNING: list_atlases is empty!'
        run_cmd.msg(nl, diary_file, 'WARNING')
    if use_erode_WM_func_masks == True:
        descmask = '_desc-erod-'
    else:
        descmask = '_desc-'

    ## MASKS in func space:
    for input1, output2 in zip([opj(template_dir_masks, species + '_mask.nii.gz'),
                                opj(template_dir_masks, species + descmask + 'Vent_mask.nii.gz'),
                                opj(template_dir_masks, species + descmask + 'White_mask.nii.gz'),
                                opj(template_dir_masks, species + '_desc-Gray_mask.nii.gz')],
                               [opj(dir_prepro_template_masks, 'mask_ref.nii.gz'),
                                opj(dir_prepro_template_masks, 'Vmask.nii.gz'),
                                opj(dir_prepro_template_masks, 'Wmask.nii.gz'),
                                opj(dir_prepro_template_masks, 'Gmask.nii.gz')]):
        command = (sing_afni + '3dresample' + overwrite +
                   ' -prefix ' + output2 +
                   ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                   ' -input ' + input1)
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": [input1,
                                  residual_template],
                      "Description": ' Resampling (3dresample, AFNI).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(output2.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)


