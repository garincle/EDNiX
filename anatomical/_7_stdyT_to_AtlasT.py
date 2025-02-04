##############################################################################
####### Co-register your study template to the atlas template ################
##############################################################################
#Path to the excels files and data structure

import os
import subprocess
import ants
import datetime
import json

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
ops = os.path.splitext

spco = subprocess.check_output
spgo = subprocess.getoutput



#################################################
########    create brain image of animal  ########
#################################################


def stdyT_to_AtlasT(list_atlases, Aseg_ref, Aseg_refLR, BASE_SS, dir_out, n_for_ANTS, aff_metric_ants, study_template_atlas_folder, Atemplate_to_Stemplate, type_of_transform_stdyT, overwrite,
                    s_bind,afni_sif,diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anatomical._7_stdyT_to_AtlasT.stdyT_to_AtlasT'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    stdy_template = opj(study_template_atlas_folder, 'studytemplate2_' + Atemplate_to_Stemplate, 'study_template.nii.gz')

    #######################################################################
    ################coregistration temaplte to template####################
    #######################################################################
    if not os.path.exists(dir_out):
        os.mkdir(dir_out)

    def flatten_comprehension(matrix):
        return [item for row in matrix for item in row]

    ####################################################################################
    list_atlases2 = []
    if list_atlases:
        list_atlases2.extend(list_atlases)
    if Aseg_ref:
        list_atlases2.extend([Aseg_ref])
    if Aseg_refLR:
        list_atlases2.extend([Aseg_refLR])

    ###### Coregistration!!!!
    print(stdy_template)
    print(BASE_SS)
    REF = ants.image_read(stdy_template)
    IMG = ants.image_read(BASE_SS)

    mtx1 = ants.registration(fixed=IMG, moving=REF, type_of_transform='Translation',
                             outprefix=opj(dir_out,'NMT_to_anat_SyN_final_shift_'))
    MEAN_tr = ants.apply_transforms(fixed=IMG, moving=REF, transformlist=mtx1['fwdtransforms'],
                                    interpolator=n_for_ANTS)
    ants.image_write(MEAN_tr, opj(dir_out,'NMT_to_anat_SyN_final_shift.nii.gz'), ri=False)

    dictionary = {"Sources": [BASE_SS,
                              stdy_template],
                  "Description": 'Co-registration (translation,ANTspy).', }
    json_object = json.dumps(dictionary, indent=2)

    with open(opj(dir_out, 'NMT_to_anat_SyN_final_shift.json'), "w") as outfile:
        outfile.write(json_object)

    mTx  = ants.registration(fixed=REF,moving=IMG,
                              type_of_transform=type_of_transform_stdyT,
                              outprefix=opj(dir_out,'NMT_to_anat_SyN_final_'),
                              aff_metric=aff_metric_ants,
                              initial_transform=mtx1['fwdtransforms'],
                              grad_step=0.1,
                              flow_sigma=3,
                              total_sigma=0,
                              aff_sampling=32,
                              aff_random_sampling_rate=0.2,
                              syn_sampling=32,
                              aff_iterations=(1000, 500, 250, 100),
                              aff_shrink_factors=(8, 4, 2, 1),
                              aff_smoothing_sigmas=(3, 2, 1, 0),
                              reg_iterations=(1000, 500, 250,100),
                              reg_smoothing_sigmas=(3, 2, 1, 0),
                              reg_shrink_factors=(8, 4, 2, 1),
                              verbose=True)

    TRANS = ants.apply_transforms(fixed=REF, moving=IMG,
                                      transformlist=mTx['fwdtransforms'], interpolator=n_for_ANTS)
    ants.image_write(TRANS, opj(dir_out,'NMT_to_anat_SyN_final.nii.gz'), ri=False)
    dictionary = {"Sources": [BASE_SS,
                              stdy_template],
                  "Description": 'Co-registration (non linear,ANTspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_out,'NMT_to_anat_SyN_final.json'), "w") as outfile:
        outfile.write(json_object)


    ####################################################################################
    ########################## seg into the native space ###################
    for atlas in list_atlases2:
        print(bcolors.OKGREEN + 'INFO: Working in sending ' + atlas + ' in anat space')
        if ope(atlas):

            IMG = ants.image_read(atlas)
            TRANS = ants.apply_transforms(fixed=REF, moving=IMG,
                                          transformlist=mTx['fwdtransforms'], interpolator='nearestNeighbor')
            ants.image_write(TRANS, opj(dir_out, opb(atlas)), ri=False)

            dictionary = {"Sources": [atlas,
                                      stdy_template],
                          "Description": 'Co-registration (non linear,ANTspy).', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_out, opb(atlas[:-7])) + '.json', "w") as outfile:
                outfile.write(json_object)

        else:
            nl = 'WARNING: ' + str(atlas) + ' not found in list_atlases, we can continue but it might restrict several outcome of the script'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()



