import os
import subprocess
import glob
from fonctions.extract_filename import extract_filename
import ants
from nilearn import plotting
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
spco = subprocess.check_output
spgo = subprocess.getoutput

def to_anat_space(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
    nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space):
    #########################################################################################################
    ################################### registration to anat space ##########################################
    #########################################################################################################

    ############################### ############################### ############################### 
    ############################### apply transfo to anat space to each volume of each func image #
    ############################### ############################### ############################### 
    if do_anat_to_func == True:
        mvt_shft_ANTs = []
        w2inv_fwd = []
        for elem1, elem2 in zip([  # opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_shift_0GenericAffine.mat'),
            opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_1Warp.nii.gz'),
            opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image_unwarped_0GenericAffine.mat')], [False, False]):
            if ope(elem1):
                mvt_shft_ANTs.append(elem1)
                w2inv_fwd.append(elem2)
    elif do_anat_to_func == False and anat_func_same_space == True:
        mvt_shft_ANTs = []
        w2inv_fwd = []
    else:
        raise Exception(bcolors.FAIL + 'ERROR: If Anat and Func are not in the same space you need to perform that trasnformation (do_anat_to_func = True)' + bcolors.ENDC)

    ## test on mean img (to see spatially that is works)
    MEAN  = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, 'Mean_Image.nii.gz'))
    BRAIN = ants.image_read(opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'))
    TRANS = ants.apply_transforms(fixed=BRAIN, moving=MEAN,
                                  transformlist=mvt_shft_ANTs,
                                  interpolator=n_for_ANTS,
                                  whichtoinvert=w2inv_fwd)
    ants.image_write(TRANS, opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.nii.gz'), ri=False)



    bids_dir = opd(opd(opd(opd(opd(dir_fMRI_Refth_RS_prepro1)))))
    if not os.path.exists(bids_dir + '/QC/meanIMG_in_anat/'): os.mkdir(bids_dir + '/QC/meanIMG_in_anat/')

    ####plot the QC
    try:
        display = plotting.plot_anat(opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'),
                                     threshold='auto',
                                     display_mode='mosaic', dim=4)
        display.add_contours(opj(dir_fMRI_Refth_RS_prepro2, 'Mean_Image_RcT_SS_in_anat.nii.gz'),
                             linewidths=.2, colors=['red'])
        display.savefig(bids_dir + '/QC/meanIMG_in_anat/Mean_Image_RcT_SS_in_anat.png')
        # Don't forget to close the display
        display.close()
    except:
        display = plotting.plot_anat(opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz'),
                                     threshold='auto',
                                     display_mode='mosaic', dim=4)
        display.savefig(bids_dir + '/QC/meanIMG_in_anat/Mean_Image_RcT_SS_in_anat.png')
        # Don't forget to close the display
        display.close()

    for i in range(int(nb_run)):
        root_RS = extract_filename(RS[i])
        if ope(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')) == False:
            FUNC_path = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_3dDeconvolve_failed.nii.gz')
            OUTPUT = opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_3dDeconvolve_failed_in_anat.nii.gz')
        else:
            FUNC_path = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')
            OUTPUT = opj(dir_fMRI_Refth_RS_prepro2, root_RS + '_residual_in_anat.nii.gz')

        FUNC = ants.image_read(FUNC_path)
        TRANS = ants.apply_transforms(fixed=BRAIN, moving=FUNC,
                                      transformlist=mvt_shft_ANTs,
                                      interpolator='nearestNeighbor',
                                      whichtoinvert=w2inv_fwd,imagetype=3)
        ants.image_write(TRANS, OUTPUT, ri=False)

