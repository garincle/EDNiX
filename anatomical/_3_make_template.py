##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################
import os
import subprocess
import anatomical.anat_to_common_EMB
import json
import datetime

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

############# use #############
#######Name: nilearn ##########
#######Version: 0.9.0 #########
###############################

def make_template(which_on, all_ID_max, all_Session_max, all_data_path_max, all_ID, all_Session, all_data_path, type_norm, study_template_atlas_folder,
                  s_bind, afni_sif,diary_file):

    ct = datetime.datetime.now()
    nl = 'Run anatomical._3_make_template.make_template'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    if which_on == 'max': # all or max
        all_ID_temp = all_ID_max
        all_Session_temp =  [all_Session_max[all_ID.index(ID)] for ID in all_ID_max]
        all_data_path_temp = all_data_path_max

    elif which_on == 'all':
        all_ID_temp = all_ID
        all_Session_temp = all_Session
        all_data_path_temp = all_data_path

    else:
        nl = 'ERROR: with which_on name, need to be all or max'
        diary.write(f'\n{nl}')
        raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

    template_list = []
    for ID, Session, data_path in zip(all_ID_temp, all_Session_temp, all_data_path_temp):

        # The anatomy
        path_anat     = opj(data_path,'anat')
        dir_native    = opj(path_anat,'native')
        dir_prepro    = opj(dir_native,'01_preprocess')

        ############ load the image for the template
        # create a large image for co-registration
        if not ope(opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz')):
            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -overwrite -I 20 -S 20 -A 20 -P 20 -L 20 -R 20 -S 20 -prefix ' + \
                      opj(all_data_path[0], 'anat', 'native', '01_preprocess', all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + \
                      ' ' + opj(all_data_path[0], 'anat', 'native', '01_preprocess', all_ID[0] + '_acpc_cropped' + type_norm + '.nii.gz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": opj(all_data_path[0], 'anat', 'native', '01_preprocess', all_ID[0] + '_acpc_cropped' + type_norm + '.nii.gz'),
                          "Description": 'Size image normalization.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(all_data_path[0], 'anat', 'native', '01_preprocess', all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.json'), "w") as outfile:
                outfile.write(json_object)


        # resample the other anat to this large img
        command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + ' -overwrite' + \
        ' -input ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz') + ' -prefix ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + ' -bound_type SLAB'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz'),
                      "Description": 'Size image normalization.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_prepro, ID + '_acpc_cropped' + type_norm + 'Zp.json'), "w") as outfile:
            outfile.write(json_object)

        template_list.append(opj(dir_prepro, ID + '_acpc_cropped' + type_norm + 'Zp.nii.gz'))

    if not ope(study_template_atlas_folder):
        os.mkdir(study_template_atlas_folder)
    templatedir2 = opj(study_template_atlas_folder, 'studytemplate2_' + type_norm)
    if not ope(templatedir2):
        os.mkdir(templatedir2)

    diary.write(f'\n')
    diary.close()

    anatomical.anat_to_common_EMB.anats_to_common(
        s_bind,afni_sif,diary_file,
        template_list,
        templatedir2,
        blur_radius_coarse=11,
        convergence=0,
        registration_kind='nonlinear',
        nonlinear_levels=[4,7],
        nonlinear_minimal_patches=[13,9])

