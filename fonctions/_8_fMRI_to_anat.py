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


def to_anat_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space,diary_file):


    #########################################################################################################
    ################################### registration to anat space ##########################################
    #########################################################################################################

    nl = '##  Working on step ' + str(8) + '(function: _8_fMRI_to_anat).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    # Extract ID
    sub_path = opn(dir_fMRI_Refth_RS_prepro2).split(os.sep)
    ID = [segment.split('-')[1] for segment in sub_path if segment.startswith('sub-')][0]

    anat_acpc = opj(dir_fMRI_Refth_RS_prepro2, ('_').join([ID, 'res-func', TfMRI + '.nii.gz']))
    ref_img = opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz')

    ############################### ############################### ############################### 
    ############################### apply transfo to anat space to each volume of each func image #
    ############################### ############################### ############################### 
    if do_anat_to_func == True:

        mvt_shft_ANTs = []
        w2inv_fwd     = []

        for elem1, elem2 in zip([ref_img.replace('.nii.gz','_unwarped_1Warp.nii.gz'),
                                 ref_img.replace('.nii.gz','_0GenericAffine.mat')], [False, False]):
            if ope(elem1):
                mvt_shft_ANTs.append(elem1)
                w2inv_fwd.append(elem2)

    elif do_anat_to_func == False and anat_func_same_space == True:

        mvt_shft_ANTs = []
        w2inv_fwd     = []

    else:
        nl = 'ERROR: If Anat and Func are not in the same space you need to perform that transformation (do_anat_to_func = True)'
        raise Exception(run_cmd.error(nl, diary_file))

    ## test on mean img (to see spatially that is works)
    MEAN  = ants.image_read(ref_img)
    BRAIN = ants.image_read(anat_acpc)
    TRANS = ants.apply_transforms(fixed=BRAIN, moving=MEAN,
                                  transformlist=mvt_shft_ANTs,
                                  interpolator=n_for_ANTS,
                                  whichtoinvert=w2inv_fwd)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.nii.gz'), ri=False)
    dictionary = {"Sources": [ref_img,
                              anat_acpc],
                  "Description": ' Non linear normalization (ANTspy).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.json'), "w") as outfile:
        outfile.write(json_object)


    bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))
    if not ope(opj(bids_dir,'QC','meanIMG_in_anat')):
        os.mkdir(opj(bids_dir,'QC','meanIMG_in_anat'))



    plot_QC_func.plot_qc(anat_acpc,
            opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.nii.gz'),
            opj(bids_dir,'QC','meanIMG_in_anat', ID + '_Mean_Image_RcT_SS_in_anat.png'))

    for i in range(int(nb_run)):
        root_RS = extract_filename(RS[i])
        if ope(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')) == False:
            FUNC_path = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_3dDeconvolve_failed.nii.gz')
            OUTPUT    = opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_3dDeconvolve_failed_in_anat.nii.gz')
        else:
            FUNC_path = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')
            OUTPUT    = opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz')

        FUNC  = ants.image_read(FUNC_path)

        TRANS = ants.apply_transforms(fixed=BRAIN, moving=FUNC,
                                      transformlist=mvt_shft_ANTs,
                                      interpolator='nearestNeighbor',
                                      whichtoinvert=w2inv_fwd,imagetype=3)
        ants.image_write(TRANS, OUTPUT, ri=False)
        dictionary = {"Sources": [FUNC_path,
                                  anat_acpc],
                      "Description": ' Non linear normalization (ANTspy).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(OUTPUT[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)



