import os
from fonctions.extract_filename import extract_filename
import ants
import datetime
import json
from fonctions import plot_QC_func
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

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists


def to_anat_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space,diary_file):


    #########################################################################################################
    ################################### registration to anat space ##########################################
    #########################################################################################################

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(8) + '(function: _8_fMRI_to_anat).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    ############################### ############################### ############################### 
    ############################### apply transfo to anat space to each volume of each func image #
    ############################### ############################### ############################### 
    if do_anat_to_func == True:

        mvt_shft_ANTs = []
        w2inv_fwd     = []

        for elem1, elem2 in zip([opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz'),
                                 opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat')], [False, False]):
            if ope(elem1):
                mvt_shft_ANTs.append(elem1)
                w2inv_fwd.append(elem2)

    elif do_anat_to_func == False and anat_func_same_space == True:

        mvt_shft_ANTs = []
        w2inv_fwd     = []

    else:
        nl = 'ERROR: If Anat and Func are not in the same space you need to perform that trasnformation (do_anat_to_func = True)'
        diary.write(f'\n{nl}')
        raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

    ## test on mean img (to see spatially that is works)
    MEAN  = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'))
    BRAIN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'))
    TRANS = ants.apply_transforms(fixed=BRAIN, moving=MEAN,
                                  transformlist=mvt_shft_ANTs,
                                  interpolator=n_for_ANTS,
                                  whichtoinvert=w2inv_fwd)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.nii.gz'), ri=False)
    dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'),
                              opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz')],
                  "Description": ' Non linear normalization (ANTspy).'},
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.json'), "w") as outfile:
        outfile.write(json_object)



    bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))
    if not os.path.exists(opj(bids_dir,'QC','meanIMG_in_anat')):
        os.mkdir(opj(bids_dir,'QC','meanIMG_in_anat'))

    # Extract ID
    sub_path = os.path.normpath(dir_fMRI_Refth_RS_prepro2).split(os.sep)
    ID = [segment.split('-')[1] for segment in sub_path if segment.startswith('sub-')][0]

    plot_QC_func.plot_qc(opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'),
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
                                  opj(dir_fMRI_Refth_RS_prepro2, 'anat_rsp_in_func.nii.gz')],
                      "Description": ' Non linear normalization (ANTspy).'},
        json_object = json.dumps(dictionary, indent=2)
        with open(OUTPUT[:-7] + '.json', "w") as outfile:
            outfile.write(json_object)

    diary.write(f'\n')
    diary.close()


