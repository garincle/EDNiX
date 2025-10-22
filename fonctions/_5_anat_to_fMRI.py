import os
import subprocess
from nilearn.image import resample_to_img
from fonctions.extract_filename import extract_filename
import ants
import json
from fonctions import plot_QC_func
import nibabel as nib

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd

def Refimg_to_meanfMRI(REF_int, SED, anat_func_same_space, TfMRI, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, RS, nb_run, ID, dir_prepro,
                       n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl, list_atlases, labels_dir, anat_subject, IhaveanANAT, do_anat_to_func, type_of_transform, registration_fast,
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

    mean_img  = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')
    mean_img2 = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift.nii.gz')
    mean_img3 = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped.nii.gz')
    mean_img4 = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_res-anat.nii.gz')

    maskDilat = opj(dir_fMRI_Refth_RS_prepro2, 'maskDilat.nii.gz')
    maskAlli  = opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat_Allineate_in_func.nii.gz')
    anat_acpc = opj(dir_fMRI_Refth_RS_prepro2, ('_').join([ID, 'res-func', TfMRI + '.nii.gz']))
    anat_func = opj(dir_fMRI_Refth_RS_prepro1, ('_').join([ID, '_space-func', TfMRI + '.nii.gz']))
    anat4fig = opj(dir_fMRI_Refth_RS_prepro2, 'orig_anat_for_plot.nii.gz')


    if IhaveanANAT == True and anat_func_same_space == True:
            refanat = opj(dir_prepro, ID + '_space-raw_desc-n4Bias_' + TfMRI + '.nii.gz')
    else:
        refanat = anat_subject

    # Load the image directly
    meanImg = nib.load(mean_img)
    
    # Get voxel sizes
    delta_x, delta_y, delta_z = [str(round(abs(x), 10)) for x in meanImg.header.get_zooms()[:3]]

    # Skull strip
    nl = f"Applied skull stripping: {mean_img} masked by {maskAlli}"
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    mask_img   = nib.load(maskAlli)
    masked_img = nib.Nifti1Image(mean_img.get_fdata() * mask_img.get_fdata(), meanImg.affine, meanImg.header)
    nib.save(masked_img, mean_img)

    # modify the json file
    json_file = mean_img.replace('.nii.gz', '.json')
    dictionary = {
        "ADD_Sources": [mean_img, maskAlli],
        "ADD_Description": "Skull stripping (Nilearn/Nibabel equivalent of AFNI 3dcalc)."}
    
    # Load existing JSON data
    with open(json_file, "r") as infile:
        existing_data = json.load(infile)
    # Update the existing data with the new dictionary
    existing_data.update(dictionary)
    # Save the updated content back to the file
    with open(json_file, "w") as outfile:
        json.dump(existing_data, outfile, indent=2)

    ####################################################################################
    ########################## use template and transfo to anat (average indiv anat)  ##
    ####################################################################################
    ## remove previous transfo files:
    for mat in ['_0GenericAffine.mat','_1Warp.nii.gz','_1InverseWarp.nii.gz']:
        if opi(mean_img3.replace('.nii.gz',mat)):
            os.remove(mean_img3.replace('.nii.gz',mat))

    if do_anat_to_func == False:
        nl = 'No anat to func step required'
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        command = (sing_afni + '3dcalc -a ' + mean_img + ' -prefix ' + mean_img3 + ' -expr "a"' + overwrite)
        run_cmd.do(command, diary_file)

    else:
        if 'i' in SED:      restrict = (1,0.1,0.1)
        elif 'j' in SED:    restrict = (0.1,1,0.1)
        elif 'k' in SED:    restrict = (0.1,0.1,1)
        elif 'None' in SED: restrict = (1, 1, 1)

        MEAN        = ants.image_read(mean_img)
        ANAT        = ants.image_read(anat_acpc)
        mask        = ants.image_read(maskDilat)
        moving_mask = ants.image_read(maskAlli)

        if anat_func_same_space == True:
            matrix = opj(dir_transfo, 'acpc_0GenericAffine.mat')
            MEAN_tr = ants.apply_transforms(fixed=ANAT, moving=MEAN,
                                            transformlist=matrix,interpolator=n_for_ANTS)
        else:
            matrix = mean_img2.replace('.nii.gz', '_0GenericAffine.mat')
            ants.registration(fixed=ANAT, moving=MEAN, type_of_transform='Translation',
                              aff_metric=aff_metric_ants_Transl,
                              outprefix=mean_img2.replace('.nii.gz', '_'))
            MEAN_tr = ants.apply_transforms(fixed=ANAT, moving=MEAN,
                                            transformlist=matrix,interpolator=n_for_ANTS)

        ants.image_write(MEAN_tr, mean_img2, ri=False)
        dictionary = {"Sources": [mean_img,
                                  anat_acpc],
                      "Description": 'Co-registration (translation, ANTspy).', },
        json_object = json.dumps(dictionary, indent=2)
        with open(mean_img2.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        if registration_fast == False:

            mTx2 = ants.registration(fixed=ANAT, moving=MEAN,
                                     type_of_transform=type_of_transform,
                                     initial_transform=matrix,
                                     outprefix=mean_img3.replace('.nii.gz','_'),
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
                                     initial_transform=mtx1['fwdtransforms'],
                                     outprefix=mean_img3.replace('.nii.gz','_'),
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
        for elem1 in (mean_img3.replace('.nii.gz','_1Warp.nii.gz'),
                      mean_img3.replace('.nii.gz','_0GenericAffine.mat')):
            if ope(elem1):
                transfo_concat.append(elem1)
        run_cmd.msg(transfo_concat, diary_file, 'ENDC')

        MEAN_trAff = ants.apply_transforms(fixed=ANAT, moving=MEAN, transformlist=transfo_concat,
                                           interpolator='nearestNeighbor')
        ants.image_write(MEAN_trAff, mean_img3, ri=False)

        dictionary = {"Sources": [mean_img,
                                  anat_acpc,
                                  mtx1['fwdtransforms']],
                      "Description": 'Co-registration (Non linear, ANTspy).', },
        json_object = json.dumps(dictionary, indent=2)
        with open(mean_img3.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)

    if do_anat_to_func == True:
        mvt_shft_INV_ANTs = []
        w2inv_inv = []
        for elem1, elem2 in zip([mean_img3.replace('.nii.gz','_0GenericAffine.mat'),
                                 mean_img3.replace('.nii.gz','_1InverseWarp.nii.gz')],
                                [True, False]):
            if opi(elem1):
                mvt_shft_INV_ANTs.append(elem1)
                w2inv_inv.append(elem2)
    elif do_anat_to_func == False and anat_func_same_space == True:
        mvt_shft_INV_ANTs = []
        w2inv_inv         = []
    else:
        nl = 'ERROR: If Anat and Func are not in the same space you need to perform that transformation (do_anat_to_func = True)'
        run_cmd.msg(nl, diary_file, 'FAIL')

    # doesn't work for two different 1d matrices... so let's do it separately....

    ## ANATOMICAL IMAGE in func space:
    BRAIN = ants.image_read(anat_subject)
    TRANS = ants.apply_transforms(fixed=MEAN, moving=BRAIN,
                                  transformlist=mvt_shft_INV_ANTs,
                                  interpolator=n_for_ANTS,
                                  whichtoinvert=w2inv_inv)
    ants.image_write(TRANS, anat_func, ri=False)
    dictionary = {"Sources": [anat_subject,
                              mean_img],
                  "Description": 'Normalization (ANTspy).', },
    json_object = json.dumps(dictionary, indent=2)
    with open(anat_func.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    ## MASKS in func space:
    for input1, output2 in zip([opj(dir_fMRI_Refth_RS_prepro2,'mask_ref.nii.gz'),
                                maskDilat,
                                opj(dir_fMRI_Refth_RS_prepro2,'Vmask.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro2,'Wmask.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro2,'Gmask.nii.gz')],
                               [opj(dir_fMRI_Refth_RS_prepro1,'mask_ref.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro1,'Wmask.nii.gz'),
                                opj(dir_fMRI_Refth_RS_prepro1,'Gmask.nii.gz')]):
        if ope(input1):
            # mask
            MEAN = ants.image_read(mean_img)
            IMG = ants.image_read(input1)
            TRANS = ants.apply_transforms(fixed=MEAN, moving=IMG,
                                          transformlist=mvt_shft_INV_ANTs,
                                          interpolator='nearestNeighbor',
                                          whichtoinvert=w2inv_inv)
            ants.image_write(TRANS,output2,ri=False)

            command = (sing_afni + '3dmask_tool' + overwrite + ' -prefix ' + output2 +
                       ' -input ' + output2 + ' -fill_holes')
            run_cmd.do(command, diary_file)

            command = (sing_afni + '3dclust -NN1 10 -prefix ' + output2 + ' ' + output2)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": [mean_img,
                                      input1],
                          "Description": ['1. Normalization (nearestNeighbor,ANTspy).',
                                          '2. fill holes (3dmask_tool, AFNI',
                                          '3. make sure there is enough voxels (3dclust, AFNI)'], }
            json_object = json.dumps(dictionary, indent=2)
            with open(output2[:-7] + '.json', "w") as outfile:
                outfile.write(json_object)

        else:
            nl = 'WARNING:' + str(input1) + ' not found!!! this may be because you have not provided an aseg file, ' + \
                 ' then no extraction of WM or Ventricles or GM will be possible... pls check that!'
            run_cmd.msg(nl, diary_file, 'WARNING')


    #### each available atlas in func space:
    if len(list_atlases) > 0:
        for atlas in list_atlases:
            if IhaveanANAT == True:
                atlasfile = ID + '_seg-' + atlas + '_dseg.nii.gz'
                input1    = opj(labels_dir, atlasfile)
            else:
                atlasfile = ID + '_seg-' + opb(atlas).split('.')[0] + '_dseg.nii.gz'
                input1    = atlas

            # resample to func
            command = (sing_afni + '3dresample' + overwrite +
                       ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro2, atlasfile) +
                       ' -dxyz ' + delta_x + ' ' + delta_y + ' ' + delta_z + ' ' +
                       ' -input  ' + input1)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": [input1,
                                      mean_img],
                          "Description": 'Resampling (3dresample, AFNI)'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro2, atlasfile.replace('.nii.gz','.json')), "w") as outfile:
                outfile.write(json_object)

            ## in func space
            ATLAS = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2, atlasfile))
            TRANS = ants.apply_transforms(fixed=MEAN, moving=ATLAS,
                                          transformlist=mvt_shft_INV_ANTs,
                                          interpolator='genericLabel',
                                          whichtoinvert=w2inv_inv)
            ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, atlasfile), ri=False)

            dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro2, atlasfile),
                                      mean_img],
                          "Description": 'Normalization (nearestNeighbor, AFNI)'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, atlasfile.replace('.nii.gz','.json')), "w") as outfile:
                outfile.write(json_object)

    else:
        nl = 'WARNING: list_atlases is empty!'
        run_cmd.msg(nl, diary_file, 'WARNING')

    # Load the image directly
    img = nib.load(anat_subject)  # <== How it is compatible with "IhaveanANAT == False"
    # Get voxel sizes
    delta_x1, delta_y1, delta_z1 = [str(round(abs(x), 10)) for x in img.header.get_zooms()[:3]]

    command = (sing_afni + '3dresample' + overwrite +
               ' -prefix ' + mean_img4 +
               ' -dxyz ' + delta_x1 + ' ' + delta_y1 + ' ' + delta_z1 + ' ' +
               ' -rmode Cu -input ' + mean_img)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [mean_img,
                              anat_subject],
                  "Description": 'resampling (3dresample, AFNI)'},
    json_object = json.dumps(dictionary, indent=2)
    with open(mean_img4.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)

    #### create a nice anat in func space
    if anat_func_same_space == True:

        matrix    = opj(dir_transfo,'acpc_0GenericAffine.mat')
        movANAT   = ants.image_read(anat_subject)
        refANA    = ants.image_read(refanat)

        TRANS = ants.apply_transforms(fixed=refANA, moving=movANAT,
                                      transformlist=mvt_shft_INV_ANTs,
                                      interpolator='hammingWindowedSinc',
                                      whichtoinvert=w2inv_inv)
        ants.image_write(TRANS,anat4fig)
        dictionary = {"Sources": [anat_subject,
                                  refanat],
                      "Description": 'Normalization (ANTSpy).'},

    else:
        command = (sing_afni + '3dcalc' + overwrite + ' -a ' + refanat + ' -prefix ' + anat4fig + ' -expr "a"')
        run_cmd.do(command, diary_file)
        dictionary = {"Sources": refanat,
                      "Description": 'Copy.'},

    json_object = json.dumps(dictionary, indent=2)
    with open(anat4fig.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    '''
    ## in func space resample to func
    ANAT = ants.image_read(anat_subject)
    MEAN_RES = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz'))
    TRANS = ants.apply_transforms(fixed=MEAN_RES, moving=ANAT,
                                  transformlist=mvt_shft_INV_ANTs,
                                  interpolator='nearestNeighbor',
                                  whichtoinvert=w2inv_inv)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.nii.gz'), ri=False)
    dictionary = {"Sources": [anatstd,
                              opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image_RcT_SS_anat_resolution.nii.gz')],
                  "Description": 'Normalization (nearestNeighbore, ANTSpy).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro1, 'Ref_anat_in_fMRI_anat_resolution.json'), "w") as outfile:
        outfile.write(json_object)
    '''

    ### finally mask the func with mask
    for i in range(int(nb_run)):
        root_RS     = extract_filename(RS[i])
        root_RS_ref = extract_filename(RS[REF_int])

        ### take opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz')

        command = (sing_afni + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') +
                   ' -b ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref.nii.gz') +
                   ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') +
                   ' -expr "a*b"')
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1,root_RS + '_xdtrf_2ref.nii.gz'),
                                  maskDilat],
                      "Description": 'Skull stripping (3dcalc, AFNI).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.json'), "w") as outfile:
            outfile.write(json_object)


        if not root_RS == root_RS_ref:  # do not process ref...
            #### send mask the bold "not in norm", essentially for QC
            mvt_shft_ANTs_func_to_norm = [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_0GenericAffine.mat'),
                                          opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_warp_1InverseWarp.nii.gz')]
            w2inv_fwd_func_to_norm = [True, False]

            mask = ants.image_read(maskDilat)
            MEAN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz'))
            TRANS = ants.apply_transforms(fixed=MEAN, moving=mask,
                                          transformlist=mvt_shft_ANTs_func_to_norm,
                                          interpolator='nearestNeighbor',
                                          whichtoinvert=w2inv_fwd_func_to_norm)
            ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.nii.gz'), ri=False)

            dictionary = {"Sources": [maskDilat,
                                      opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz')],
                          "Description": 'Normalization (ANTspy).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.json'), "w") as outfile:
                outfile.write(json_object)

        else:
            command = (sing_afni + '3dcopy ' + maskDilat +
                       ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.nii.gz') + overwrite)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": maskDilat,
                          "Description": 'Copy.'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.json'), "w") as outfile:
                outfile.write(json_object)


        bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))
        if not ope(opj(bids_dir,'QC')):
            os.mkdir(opj(bids_dir,'QC'))

        bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))
        if not ope(opj(bids_dir,'QC','mask_to_fMRI_orig')):
            os.mkdir(opj(bids_dir,'QC','mask_to_fMRI_orig'))

        ####plot the QC
        plot_QC_func.plot_qc(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_mean_deob.nii.gz'),
                             opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_mask_final_in_fMRI_orig.nii.gz'),
                             opj(bids_dir, 'QC', 'mask_to_fMRI_orig', root_RS + '_mask_final_in_fMRI_orig.png'))

