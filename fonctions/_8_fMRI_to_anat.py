import os
import ants
import json

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from fonctions.extract_filename import extract_filename
from fonctions import plot_QC_func


def to_anat_space(dir_prepro_acpc_process, dir_prepro_orig_process, bids_dir, ID, TfMRI,
    nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space, dir_prepro_acpc_postprocessed, dir_prepro_orig_postprocessed, diary_file):

    #########################################################################################################
    ################################### registration to anat space ##########################################
    #########################################################################################################

    nl = '##  Working on step ' + str(8) + '(function: _8_fMRI_to_anat).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    anat_res_func = opj(dir_prepro_acpc_process, ('_').join(['anat_space-acpc_res-func', TfMRI + '.nii.gz']))
    Mean_Image_unwarped = opj(dir_prepro_acpc_process, 'all_runs_space-anat_desc-fMRI_Mean_Image_unwarped.nii.gz')
    Mean_Image_in_anat = opj(dir_prepro_acpc_process, 'all_runs_space-anat_desc-fMRI_Mean_Image_apply.nii.gz')
    Mean_Image_acpc = opj(dir_prepro_orig_process, 'all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS.nii.gz')

    ############################### ############################### ############################### 
    ############################### TEST apply transfo to MEAN inamge in acpc                ######
    ############################### ############################### ############################### 
    if do_anat_to_func == True:
        mvt_shft_ANTs = []
        w2inv_fwd     = []
        for elem1, elem2 in zip([Mean_Image_unwarped.replace('.nii.gz','_unwarped_1Warp.nii.gz'),
                                 Mean_Image_unwarped.replace('.nii.gz','_0GenericAffine.mat')], [False, False]):
            if ope(elem1):
                mvt_shft_ANTs.append(elem1)
                w2inv_fwd.append(elem2)

    elif do_anat_to_func == False and anat_func_same_space == True:
        mvt_shft_ANTs = []
        w2inv_fwd     = []
    else:
        nl = 'ERROR: If Anat and Func are not in the same space you need to perform that transformation (do_anat_to_func = True)'
        raise Exception(run_cmd.error(nl, diary_file))

    ## test on mean img (to see spatially it works)
    MEAN  = ants.image_read(Mean_Image_acpc)
    BRAIN = ants.image_read(anat_res_func)
    TRANS = ants.apply_transforms(fixed=BRAIN, moving=MEAN,
                                  transformlist=mvt_shft_ANTs,
                                  interpolator=n_for_ANTS,
                                  whichtoinvert=w2inv_fwd)
    ants.image_write(TRANS, Mean_Image_in_anat, ri=False)
    dictionary = {"Sources": [Mean_Image_unwarped,
                              anat_res_func],
                  "Description": ' Non linear normalization (ANTspy).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(Mean_Image_in_anat.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    # QC plot
    if not ope(opj(bids_dir,'QC','meanIMG_in_anat')):
        os.mkdir(opj(bids_dir,'QC','meanIMG_in_anat'))

    plot_QC_func.plot_qc(anat_res_func,
            Mean_Image_in_anat,
            opj(bids_dir,'QC','meanIMG_in_anat', ID + '_Mean_Image_RcT_SS_in_anat.png'))

    ############################### ############################### ############################### 
    ###############################  apply transfo to func inamge in acpc                ##########
    ############################### ############################### ############################### 

    for i in range(int(nb_run)):
        root_RS = extract_filename(RS[i])
        if ope(opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')) == False:
            residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
            residual_anat = opj(dir_prepro_acpc_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')
        else:
            residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
            residual_anat = opj(dir_prepro_acpc_postprocessed, root_RS + '_space-acpc-anat_desc-fMRI_residual.nii.gz')

        FUNC  = ants.image_read(residual)

        TRANS = ants.apply_transforms(fixed=BRAIN, moving=FUNC,
                                      transformlist=mvt_shft_ANTs,
                                      interpolator='bSpline',
                                      whichtoinvert=w2inv_fwd,
                                      imagetype=3)
        ants.image_write(TRANS, residual_anat, ri=False)
        dictionary = {"Sources": [residual,
                                  anat_res_func],
                      "Description": ' Non linear normalization (ANTspy).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(residual_anat[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)



