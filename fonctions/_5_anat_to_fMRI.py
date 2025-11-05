import os
import subprocess
from nilearn.image import resample_to_img
from fonctions.extract_filename import extract_filename
import ants
import json
from fonctions import plot_QC_func
import nibabel as nib
from fonctions import transfo_for_func
opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd, check_nii
from fonctions import _2b_fix_orient
def Refimg_to_meanfMRI(SED, anat_func_same_space, TfMRI, dir_prepro_raw_process, RS, nb_run, ID, bids_dir, dir_prepro_raw_masks, REF_int, dir_prepro_raw_matrices, recordings,
                       n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl, list_atlases, labels_dir, anat_subject, dir_transfo, IhaveanANAT, do_anat_to_func, TR_val,
                       type_of_transform, registration_fast, dir_prepro_acpc_masks, dir_prepro_acpc_process, dir_prepro_orig_masks, dir_prepro_acpc_labels,
                       dir_prepro_orig_labels, BASE_atlas_folder, dir_prepro_orig_process, doWARPonfunc, template_dir_labels, species, template_dir_masks,
                       overwrite, sing_afni,diary_file):

    nl = '##  Working on step ' + str(5) + '(function: _5_anat_to_fMRI).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    vox = 200
    direction = ['I', 'S', 'A', 'P', 'L', 'R']
    cmd = []
    for i in range(len(direction)):
        cmd.append('-' + direction[i])
        cmd.append(str(vox))
    pad = ' '.join(cmd)

    if IhaveanANAT == True:
        orig_mask = dir_prepro_acpc_masks
    else:
        orig_mask = template_dir_masks

    Mean_Image = opj(dir_prepro_raw_process, 'all_runs_space-func_desc-fMRI_Mean_Image.nii.gz')
    Mean_Image_SS = opj(dir_prepro_raw_process, 'all_runs_space-func_desc-fMRI_Mean_Image_SS.nii.gz')
    Prepro_fMRI_mask = opj(dir_prepro_raw_masks, ID + '_fMRI_mask.nii.gz')

    Mean_Image_res_anat = opj(dir_prepro_orig_process, 'all_runs_space-acpc-func_desc-fMRI_Mean_Image_res-anat.nii.gz')
    maskDilat_funcspace = opj(dir_prepro_orig_masks, ID + '_space-acpc-func_desc-fMRI_mask_dilated.nii.gz')
    Mean_Image_acpc = opj(dir_prepro_orig_process, 'all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS.nii.gz')
    anat_in_acpc_func = opj(dir_prepro_orig_process, ('_').join(['anat_space-acpc-func', TfMRI + '.nii.gz']))

    Mean_Image_unwarped = opj(dir_prepro_acpc_process, 'all_runs_space-anat_desc-fMRI_Mean_Image_unwarped.nii.gz')
    maskDilatfunc = opj(orig_mask, ID + '_space-acpc_mask_dilated_res_func.nii.gz')
    anat_res_func = opj(dir_prepro_acpc_process, ('_').join(['anat_space-acpc_res-func', TfMRI + '.nii.gz']))

    ####################################################################################
    ########################## skullstrip func images  #################################
    ####################################################################################
    # Load the image directly
    meanImg = nib.load(Mean_Image)
    # Skull strip
    nl = f"Applied skull stripping: {Mean_Image_SS} masked by {Prepro_fMRI_mask}"
    run_cmd.msg(nl, diary_file, 'OKGREEN')
    mask_img   = nib.load(Prepro_fMRI_mask)
    masked_img = nib.Nifti1Image(meanImg.get_fdata() * mask_img.get_fdata(), meanImg.affine, meanImg.header)
    nib.save(masked_img, Mean_Image_SS)

    # modify the json file
    json_file = Mean_Image_SS.replace('.nii.gz', '.json')
    dictionary = {
        "ADD_Sources": [Mean_Image_SS, Prepro_fMRI_mask],
        "ADD_Description": "Skull stripping (Nilearn/Nibabel equivalent of AFNI 3dcalc)."}
    json_object = json.dumps(dictionary, indent=2)
    with open(Mean_Image_SS.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################
    ## remove previous transfo files:
    for mat in ['_0GenericAffine.mat','_1Warp.nii.gz','_1InverseWarp.nii.gz']:
        if opi(Mean_Image_unwarped.replace('.nii.gz',mat)):
            os.remove(Mean_Image_unwarped.replace('.nii.gz',mat))

    MEAN        = ants.image_read(Mean_Image_SS)
    ANAT        = ants.image_read(anat_res_func)
    mask        = ants.image_read(maskDilatfunc)
    moving_mask = ants.image_read(Prepro_fMRI_mask)

    if 'i' in SED:      restrict = (1,0.1,0.1)
    elif 'j' in SED:    restrict = (0.1,1,0.1)
    elif 'k' in SED:    restrict = (0.1,0.1,1)
    else:               restrict = (1, 1, 1)

    if anat_func_same_space == True:
        matrix = opj(dir_transfo, 'acpc_0GenericAffine.mat')
        MEAN_tr = ants.apply_transforms(fixed=ANAT, moving=MEAN, whichtoinvert=[False],
                                        transformlist=matrix, interpolator='nearestNeighbor')
    else:
        matrix = Mean_Image_acpc.replace('.nii.gz', '_0GenericAffine.mat')
        ants.registration(fixed=ANAT, moving=MEAN, type_of_transform='Translation',
                          aff_metric=aff_metric_ants_Transl,
                          outprefix=Mean_Image_acpc.replace('.nii.gz', '_'))
        MEAN_tr = ants.apply_transforms(fixed=ANAT, moving=MEAN,
                                        transformlist=matrix, interpolator='nearestNeighbor')
    ants.image_write(MEAN_tr, Mean_Image_acpc, ri=False)
    dictionary = {"Sources": [Mean_Image_SS,
                              anat_res_func],
                  "Description": 'Co-registration (translation, ANTspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(Mean_Image_acpc.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    MEAN = ants.image_read(Mean_Image_acpc)
    if do_anat_to_func == False:
        nl = 'No anat to func step required'
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        command = (sing_afni + '3dcalc -a ' + Mean_Image_acpc + ' -prefix ' + Mean_Image_unwarped + ' -expr "a"' + overwrite)
        run_cmd.do(command, diary_file)
    else:
        nl = 'Anat to func step required'
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        if registration_fast == False:
            mTx2 = ants.registration(fixed=ANAT, moving=MEAN,
                                     type_of_transform=type_of_transform,
                                     outprefix=Mean_Image_unwarped.replace('.nii.gz','_'),
                                     grad_step=0.1,flow_sigma=3,total_sigma=0,aff_sampling=32,
                                     aff_random_sampling_rate=0.2,
                                     syn_sampling=32,
                                     aff_iterations=(1000, 500, 250, 100),
                                     aff_shrink_factors=(8, 4, 2, 1),
                                     aff_smoothing_sigmas=(3, 2, 1, 0),
                                     reg_iterations=(1000, 500, 250, 100),
                                     reg_smoothing_sigmas=(3, 2, 1, 0),
                                     reg_shrink_factors=(8, 4, 2, 1),
                                     verbose=True,
                                     mask=mask,
                                     moving_mask=moving_mask,
                                     aff_metric=aff_metric_ants,
                                     restrict_transformation=restrict)
    
        if registration_fast == True:
            nl = "registration_fast selected"
            run_cmd.msg(nl, diary_file, 'OKGREEN')
    
            mTx2 = ants.registration(fixed=ANAT, moving=MEAN,
                                     type_of_transform=type_of_transform,
                                     outprefix=Mean_Image_unwarped.replace('.nii.gz','_'),
                                     grad_step=0.1,  # Step size for Rigid transform
                                     flow_sigma=3, total_sigma=0,  # Flow sigma
                                     aff_sampling=32,  # Sampling rate
                                     aff_random_sampling_rate=0.25,  # Random sampling rate
                                     aff_iterations=(100, 50, 25, 10),  # Number of iterations
                                     aff_shrink_factors=(8, 4, 2, 1),  # Shrink factors
                                     aff_smoothing_sigmas=(3, 2, 1, 0),  # Smoothing sigmas
                                     aff_metric_params=(ANAT, MEAN, 1, 32, 'Regular', 0.25),  # Metric parameters
                                     convergence=(1e-6, 10),  # Convergence criteria
                                     verbose=True,
                                     mask=mask,
                                     moving_mask=moving_mask,
                                     aff_metric=aff_metric_ants,
                                     restrict_transformation=restrict)

        transfo_concat = []
        for elem1 in (Mean_Image_unwarped.replace('.nii.gz','_1Warp.nii.gz'),
                      Mean_Image_unwarped.replace('.nii.gz','_0GenericAffine.mat')):
            if ope(elem1):
                transfo_concat.append(elem1)
        run_cmd.msg(str(transfo_concat), diary_file, 'ENDC')
    
        MEAN_trAff = ants.apply_transforms(fixed=ANAT, moving=MEAN, transformlist=transfo_concat,
                                           interpolator='nearestNeighbor')
        ants.image_write(MEAN_trAff, Mean_Image_unwarped, ri=False)
    
        dictionary = {"Sources": [Mean_Image_acpc,
                                  anat_res_func,
                                  matrix],
                      "Description": 'Co-registration (Non linear, ANTspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(Mean_Image_unwarped.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)

    if do_anat_to_func == True:
        mvt_shft_INV_ANTs = []
        w2inv_inv = []
        for elem1, elem2 in zip([Mean_Image_unwarped.replace('.nii.gz','_0GenericAffine.mat'),
                                 Mean_Image_unwarped.replace('.nii.gz','_1InverseWarp.nii.gz')],
                                [True, False]):
            if opi(elem1):
                mvt_shft_INV_ANTs.append(elem1)
                w2inv_inv.append(elem2)

    meanImg = nib.load(Mean_Image_acpc)
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in meanImg.header.get_zooms()[:3]]

    if do_anat_to_func == True:
        ## ANATOMICAL IMAGE in acpc-func space:
        BRAIN = ants.image_read(anat_subject)
        TRANS = ants.apply_transforms(fixed=MEAN, moving=BRAIN,
                                      transformlist=mvt_shft_INV_ANTs,
                                      interpolator=n_for_ANTS,
                                      whichtoinvert=w2inv_inv)
        ants.image_write(TRANS, anat_in_acpc_func, ri=False)
        dictionary = {"Sources": [anat_subject,
                                  Mean_Image_SS],
                      "Description": 'Normalization (ANTspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(anat_in_acpc_func.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)
    else:
        # resample to func
        command = (sing_afni + '3dresample' + overwrite +
                   ' -prefix ' + anat_in_acpc_func +
                   ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                   ' -input ' + anat_subject)
        print(command)
        dictionary = {"Sources": [anat_subject,
                                  Mean_Image_acpc],
                      "Description": 'Resampling (3dresample, AFNI)', "Command": command, }
        json_object = json.dumps(dictionary, indent=3)
        with open(anat_in_acpc_func.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

    ## MASKS in func space:
    for input1, output2 in zip([opj(orig_mask,'mask_ref.nii.gz'),
                                maskDilatfunc,
                                opj(orig_mask,'Vmask.nii.gz'),
                                opj(orig_mask,'Wmask.nii.gz'),
                                opj(orig_mask,'Gmask.nii.gz')],
                               [opj(dir_prepro_orig_masks,'mask_ref.nii.gz'),
                                maskDilat_funcspace,
                                opj(dir_prepro_orig_masks,'Vmask.nii.gz'),
                                opj(dir_prepro_orig_masks,'Wmask.nii.gz'),
                                opj(dir_prepro_orig_masks,'Gmask.nii.gz')]):
        if ope(input1):
            if do_anat_to_func == True:
                # mask
                MEAN = ants.image_read(Mean_Image_acpc)
                IMG = ants.image_read(input1)
                TRANS = ants.apply_transforms(fixed=MEAN, moving=IMG,
                                              transformlist=mvt_shft_INV_ANTs,
                                              interpolator='nearestNeighbor',
                                              whichtoinvert=w2inv_inv)
                ants.image_write(TRANS,output2,ri=False)
            elif do_anat_to_func == False and anat_func_same_space == True:
                nl = 'No anat to func step required'
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                # resample to func
                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + output2 +
                           ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                           ' -input ' + input1)
                run_cmd.run(command, diary_file)
            else:
                nl = 'ERROR: If Anat and Func are not in the same space you need to perform that transformation (do_anat_to_func = True)'
                run_cmd.msg(nl, diary_file, 'FAIL')

            command = (sing_afni + '3dmask_tool' + overwrite + ' -prefix ' + output2 +
                       ' -input ' + output2 + ' -fill_holes')
            run_cmd.do(command, diary_file)

            command = (sing_afni + '3dclust' + overwrite + ' -NN1 10 -prefix ' + output2 + ' ' + output2)
            dictionary = {"Sources": [Mean_Image_SS,
                                      input1],
                          "Description": ['1. Normalization (nearestNeighbor,ANTspy).',
                                          '2. fill holes (3dmask_tool, AFNI',
                                          '3. make sure there is enough voxels (3dclust, AFNI)'], "Command": command,}
            json_object = json.dumps(dictionary, indent=3)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)
            run_cmd.run(command, diary_file)
        else:
            nl = 'WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, ' + \
                 ' then no extraction of WM or Ventricles or GM will be possible... pls check that!'
            run_cmd.msg(nl, diary_file, 'WARNING')
    print("Work on list atlas" + str(list_atlases))
    #### each available atlas in anat space:
    if len(list_atlases[0]) > 0:
        for atlas in list_atlases[0]:
            if IhaveanANAT == True:
                atlasfile = ID + '_seg-' + str(atlas) + '_dseg.nii.gz'
                input1    = opj(labels_dir, atlasfile)
            else:
                input1    = opj(template_dir_labels, species + '_seg-' + atlas + '_dseg.nii.gz')
            if not ope(input1):
                nl = 'WARNING: atlas file ' + str(input1) + ' not found!'
                print(nl)
            else:
                # resample to func
                command = (sing_afni + '3dresample' + overwrite +
                           ' -prefix ' + opj(dir_prepro_acpc_labels, atlasfile) +
                           ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                           ' -input ' + input1)
                print(command)
                dictionary = {"Sources": [input1,
                                          Mean_Image_acpc],
                              "Description": 'Resampling (3dresample, AFNI)', "Command": command, }
                json_object = json.dumps(dictionary, indent=3)
                with open(opj(dir_prepro_acpc_labels, atlasfile.replace('.nii.gz','.json')), "w") as outfile:
                    outfile.write(json_object)
                run_cmd.run(command, diary_file)

                if do_anat_to_func == True:
                    ## in func space
                    ATLAS = ants.image_read(opj(dir_prepro_acpc_labels, atlasfile))
                    # Detect image dimension (3D vs 4D)
                    img_dim = len(ATLAS.shape)
                    MEAN = ants.image_read(Mean_Image_acpc)
                    # Decide image type: 0 = scalar, 3 = time series
                    image_type = 3 if img_dim == 4 else 0
                    TRANS = ants.apply_transforms(
                        fixed=MEAN,
                        moving=ATLAS,
                        transformlist=mvt_shft_INV_ANTs,
                        interpolator='genericLabel',
                        whichtoinvert=w2inv_inv,
                        imagetype=image_type)
                    ants.image_write(TRANS, opj(dir_prepro_orig_labels, atlasfile), ri=False)

                    dictionary = {"Sources": [opj(dir_prepro_acpc_labels, atlasfile),
                                              Mean_Image_acpc],
                                  "Description": 'Normalization (nearestNeighbor, AFNI)'}
                    json_object = json.dumps(dictionary, indent=2)
                    with open(opj(dir_prepro_orig_labels, atlasfile.replace('.nii.gz','.json')), "w") as outfile:
                        outfile.write(json_object)
                else:
                    # resample to func
                    command = (sing_afni + '3dresample' + overwrite +
                               ' -prefix ' + opj(dir_prepro_orig_labels, atlasfile) +
                               ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z +
                               ' -input ' + input1)
                    print(command)
                    dictionary = {"Sources": [input1,
                                              Mean_Image_acpc],
                                  "Description": 'Resampling (3dresample, AFNI)', "Command": command, }
                    json_object = json.dumps(dictionary, indent=3)
                    with open(opj(dir_prepro_orig_labels, atlasfile.replace('.nii.gz', '.json')), "w") as outfile:
                        outfile.write(json_object)
                    run_cmd.run(command, diary_file)
    else:
        nl = 'WARNING: list_atlases is empty!'
        run_cmd.msg(nl, diary_file, 'WARNING')

    # Load the image directly
    img = nib.load(anat_subject)  # <== How it is compatible with "IhaveanANAT == False"
    # Get voxel sizes
    delta_x1, delta_y1, delta_z1 = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]

    command = (sing_afni + '3dresample' + overwrite +
               ' -prefix ' + Mean_Image_res_anat +
               ' -dxyz ' + delta_x1 + ' ' + delta_y1 + ' ' + delta_z1 +
               ' -rmode Cu -input ' + Mean_Image_acpc)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [Mean_Image_acpc,
                              anat_subject],
                  "Description": 'resampling (3dresample, AFNI)',
                  "Command": command, }
    json_object = json.dumps(dictionary, indent=3)
    with open(Mean_Image_res_anat.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)
    run_cmd.run(command, diary_file)

    ### finally mask the func with mask
    for i in range(int(nb_run)):
        root_RS     = extract_filename(RS[i])

        #### create the QC folder
        if not ope(opj(bids_dir,'QC')):
            os.mkdir(opj(bids_dir,'QC'))
        if not ope(opj(bids_dir,'QC','mask_to_fMRI_orig')):
            os.mkdir(opj(bids_dir,'QC','mask_to_fMRI_orig'))

        ####plot the QC
        plot_QC_func.plot_qc(Mean_Image_acpc,
                             opj(dir_prepro_orig_masks,'Gmask.nii.gz'),
                             opj(bids_dir, 'QC', 'mask_to_fMRI_orig', root_RS + '_fMRI_in_anat.png'))


    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###
    ###                                          co-registration of each run to the norm                                                ###
    ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ### ### ### #### ###
    root_RS_ref = extract_filename(RS[REF_int])
    for r in range(0, int(nb_run)):
        root_RS = extract_filename(RS[r])

        fMRI_BASE = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_BASE.nii.gz')
        fMRI_BASE_Mean = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_BASE_Mean.nii.gz')

        fMRI_runMean_inRef_acpc = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_runMean_inRef.nii.gz')
        fMRI_run_inRef_acpc = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_run_inRef.nii.gz')
        fMRI_run_inRef_acpc_SS = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_SS.nii.gz')

        imagetype = [3, 2]
        for image_to_send_acpc, image_to_creat_acpc, imagetype in zip([fMRI_BASE, fMRI_BASE_Mean], [fMRI_run_inRef_acpc, fMRI_runMean_inRef_acpc], imagetype):
            additional_transformations = []
            if root_RS == root_RS_ref and recordings == 'very_old' and doWARPonfunc == True:  # do not process ref not corrected...
                # pas de transfo run to MEAN
                print('No transformation for the reference run')
            else:
                fMRI_run_inRef_WARP = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-fMRI_run_inRef1Warp.nii.gz')
                fMRI_run_inRef_mat = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-fMRI_run_inRef0GenericAffine.mat')

                # AJOUTER DANS L'ORDRE INVERSE
                additional_transformations.append(fMRI_run_inRef_WARP)  # 3. Warp (appliqué en premier)
                additional_transformations.append(fMRI_run_inRef_mat)  # 2. Affine

            if anat_func_same_space == True:
                matrix = opj(dir_transfo, 'acpc_0GenericAffine.mat')
            else:
                matrix = Mean_Image_acpc.replace('.nii.gz', '_0GenericAffine.mat')

            additional_transformations.append(matrix)

            print(additional_transformations)
            if imagetype == 3:
                transfo_for_func.apply_motion_correction_and_transforms(
                    input_4d=image_to_send_acpc,
                    reference_image=Mean_Image_acpc,
                    output_4d=opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_run_inRef'),
                    mat_files_pattern=opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-motion_correction_*.mat'),
                    additional_transformations=additional_transformations,
                    TR=TR_val,
                    sing_afni=sing_afni,
                    interpolator=n_for_ANTS,
                    metadata_sources=[image_to_send_acpc, Mean_Image_acpc],
                    description='Individual motion correction + ACPC normalization')

                nl = f"Applied skull stripping: {fMRI_run_inRef_acpc} masked by {maskDilat_funcspace}"
                command = (sing_afni + '3dcalc -overwrite -a ' + fMRI_run_inRef_acpc + ' -b ' + maskDilat_funcspace +
                           ' -expr "a*b" -prefix ' + fMRI_run_inRef_acpc_SS)
                run_cmd.do(command, diary_file)

                dictionary = {
                    "ADD_Sources": [fMRI_run_inRef_acpc_SS, Prepro_fMRI_mask],
                    "ADD_Description": "Skull stripping (Nilearn/Nibabel equivalent of AFNI 3dcalc)."}
                json_object = json.dumps(dictionary, indent=2)
                with open(fMRI_run_inRef_acpc_SS.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)

            elif imagetype == 2:
                FUNCACPC = ants.image_read(Mean_Image_acpc)
                FUNC = ants.image_read(image_to_send_acpc)
                MEAN_tr = ants.apply_transforms(fixed=FUNCACPC, moving=FUNC,
                                                transformlist=additional_transformations, interpolator=n_for_ANTS, imagetype=0)
                ants.image_write(MEAN_tr, image_to_creat_acpc, ri=False)
                dictionary = {"Sources": [image_to_send_acpc,
                                          Mean_Image_acpc],
                              "Description": 'Co-registration (translation, ANTspy).', }
                json_object = json.dumps(dictionary, indent=2)
                with open(image_to_creat_acpc.replace('.nii.gz', '.json'), "w") as outfile:
                    outfile.write(json_object)

                if not ope(opj(bids_dir, 'QC')):
                    os.mkdir(opj(bids_dir, 'QC'))
                if not ope(opj(bids_dir, 'QC', 'fMRI_runMean_in_REF')):
                    os.mkdir(opj(bids_dir, 'QC', 'fMRI_runMean_in_REF'))

                ####plot the QC
                plot_QC_func.plot_qc(Mean_Image_acpc,
                                     image_to_creat_acpc,
                                     opj(bids_dir, 'QC', 'fMRI_MEAN_in_REF_FUNC_ACPC', root_RS + 'fMRI_in_REF_FUNC_ACPC.png'))






