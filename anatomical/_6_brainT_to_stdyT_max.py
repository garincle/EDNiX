import os
import subprocess
from nilearn import plotting
import ants
import json
import datetime

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

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

def brainT_to_T_max(aff_metric_ants, creat_study_template, dir_prepro, ID, Session, listTimage, n_for_ANTS,
                    dir_transfo, type_norm, BASE_SS_coregistr, Ref_file, volumes_dir, transfo_concat_inv,w2inv_inv,
                    study_template_atlas_forlder, otheranat, bids_dir, type_of_transform,
                    which_on, all_data_path_max, IgotbothT1T2, all_data_path,
                    s_bind, afni_sif,diary_file):

    ct = datetime.datetime.now()
    nl = 'Run anatomical._6_brainT_to_stdyT_max.brainT_to_T_max'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    #################################################################### coregistration with ANTs ####################################################################
    ######Coregistration!!!! (template to anat)
    IMG = ants.image_read(Ref_file)
    REF = ants.image_read(BASE_SS_coregistr)

    nl = 'INFO: Co-registration ready'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')
    nl = 'INFO: WE WILL COREGISTER: ' + Ref_file + '(the anat img)' + ' to ' + BASE_SS_coregistr + '(the template)'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')
    nl = 'INFO: type_of_transform=' + type_of_transform
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    mtx1 = ants.registration(fixed=IMG, moving=REF, type_of_transform='Translation',
                             outprefix=opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_shift_'))

    MEAN_tr = ants.apply_transforms(fixed=IMG, moving=REF, transformlist=mtx1['fwdtransforms'],
                                    interpolator=n_for_ANTS)
    ants.image_write(MEAN_tr, opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_shift.nii.gz'), ri=False)

    dictionary = {"Sources": [Ref_file,
                              BASE_SS_coregistr],
                  "Description": 'Co-registration (translation).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_prepro, 'template_to_' + type_norm + '_SyN_final_max_shift.json'), "w") as outfile:
        outfile.write(json_object)


    mTx  = ants.registration(fixed=IMG,moving=REF,
                              outprefix=opj(dir_transfo,'template_to_' + type_norm + '_SyN_final_max_'),
                            type_of_transform=type_of_transform,
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

    TRANS = ants.apply_transforms(fixed=IMG, moving=REF,
                                  transformlist=mTx['fwdtransforms'], interpolator=n_for_ANTS)
    ants.image_write(TRANS, opj(dir_prepro, 'template_to_' + type_norm + '_SyN_final_max.nii.gz'), ri=False)
    dictionary = {"Sources": [Ref_file,
                              BASE_SS_coregistr],
                  "Description": 'Co-registration (non linear,ANTspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_prepro, 'template_to_' + type_norm + '_SyN_final_max.json'), "w") as outfile:
        outfile.write(json_object)


    ######### apply the transformation to the original images to see the reverse transfo paramters works
    for Timage in listTimage:
        img_ref = opj(volumes_dir,ID + Timage + '_brain_step_1.nii.gz')
        IMG = ants.image_read(img_ref)
        TRANS = ants.apply_transforms(fixed=REF, moving=IMG,
                                      transformlist=transfo_concat_inv,
                                      interpolator=n_for_ANTS, whichtoinvert=w2inv_inv)
        ants.image_write(TRANS, opj(dir_prepro, Timage + '_to_template_SyN_final_max.nii.gz'), ri=False)

        dictionary = {"Sources": [opj(volumes_dir, ID + Timage + '_brain_step_1.nii.gz'),
                                  BASE_SS_coregistr],
                      "Description": 'Co-registration (non linear,ANTspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, Timage + '_to_template_SyN_final_max.json'), "w") as outfile:
            outfile.write(json_object)

        if not ope(opj(bids_dir + 'QC','anat_to_template')): os.mkdir(opj(bids_dir + 'QC','anat_to_template'))
        if not ope(opj(bids_dir + 'QC','template_to_anat')): os.mkdir(opj(bids_dir + 'QC','template_to_anat'))

        ####plot the QC
        try:
            display = plotting.plot_anat(BASE_SS_coregistr, threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.add_contours(opj(dir_prepro, Timage + '_to_template_SyN_final_max.nii.gz'),
                                 linewidths=.2, colors=['red'])
            display.savefig(opj(bids_dir + 'QC','anat_to_template' + ID + '_' + str(Session) + '_' + Timage + 'to_template_max.png'))
            # Don't forget to close the display
            display.close()
        except:
            display = plotting.plot_anat(opj(dir_prepro,'template_to_' + Timage + '_SyN_final_test_matrix_max.nii.gz'), threshold='auto',
                                         display_mode='mosaic', dim=4)
            display.savefig(opj(bids_dir + 'QC','template_to_anat' + ID + '_' + str(Session) + '_' + Timage + 'to_template_max.png'))
            # Don't forget to close the display
            display.close()


    if creat_study_template == True and IgotbothT1T2 == True:
        ###creat brain of the subject template
        # average the available T1 files
        if which_on == 'max':  # all or max
            all_data_path_temp = all_data_path_max

        elif which_on == 'all':
            all_data_path_temp = all_data_path

        second_template = []
        for data_path in all_data_path_temp:
            path_anat  = opj(data_path, 'anat')
            dir_native = opj(path_anat, 'native')
            dir_prepro = opj(dir_native, '01_preprocess')
            second_template.append(opj(dir_prepro, Timage + '_to_template_test_matrix_max.nii.gz'))

        # Check if all paths exist
        all_paths_exist = all(ope(path) for path in second_template)

        if all_paths_exist:
            TEMP = ' '.join(second_template)
            if not ope(opj(study_template_atlas_forlder, 'studytemplate2_' + otheranat)): os.mkdir(
                opj(study_template_atlas_forlder, 'studytemplate2_' + otheranat))

            command = 'singularity run' + s_bind + afni_sif + '3dMean -prefix ' + opj(study_template_atlas_forlder, 'studytemplate2_' + otheranat, 'study_template.nii.gz') + ' ' + TEMP
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            dictionary = {"Sources": TEMP,
                          "Description": 'Mean.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(study_template_atlas_forlder, 'studytemplate2_' + otheranat, 'study_template.json'), "w") as outfile:
                outfile.write(json_object)

        diary.write(f'\n')
        diary.close()