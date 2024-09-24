##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################
import os
import subprocess
import anatomical.anat_to_common_EMB
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
#######Name: nilearn
#######Version: 0.9.0
################################


def make_template(which_on, all_ID_max, max_session, all_data_path_max, all_ID, all_Session, all_data_path, type_norm, study_template_atlas_forlder,
                  s_bind, afni_sif):

    if which_on == 'max': # all or max
        all_ID_temp = all_ID_max
        all_Session_temp = max_session
        all_data_path_temp = all_data_path_max

    elif which_on == 'all':
        all_ID_temp = all_ID
        all_Session_temp = all_Session
        all_data_path_temp = all_data_path

    else:
        raise Exception(bcolors.FAIL + 'ERROR: with which_on name, need to be all or max' + bcolors.ENDC)

    template_list = []
    for ID, Session, data_path in zip(all_ID_temp, all_Session_temp, all_data_path_temp):

        # The anatomy
        path_anat    = opj(data_path,'anat/')
        dir_transfo  = opj(path_anat,'matrices')

        dir_native    = opj(path_anat,'native')
        dir_prepro    = opj(dir_native,'01_preprocess')
        wb_native_dir = opj(dir_native,'02_Wb')
        volumes_dir   = opj(wb_native_dir,'volumes')
        labels_dir    = opj(volumes_dir,'labels')
        masks_dir     = opj(volumes_dir,'masks')

        ############load the image for the template
        #creat a large image for co-registration
        if not ope(opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz')):
            command = 'singularity run' + s_bind + afni_sif + '3dZeropad -overwrite -I 20 -S 20 -A 20 -P 20 -L 20 -R 20 -S 20 -prefix ' + opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + \
                      ' ' + opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + '.nii.gz')
            spco(command, shell=True)


        #resemple the other anat to this large img
        command = 'singularity run' + s_bind + afni_sif + '3dresample -overwrite -master ' + opj(opj(all_data_path[0], 'anat', 'native', '01_preprocess'), all_ID[0] + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + ' -overwrite' + \
        ' -input ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + '.nii.gz') + ' -prefix ' + opj(dir_prepro, ID + '_acpc_cropped' + type_norm + 'Zp.nii.gz') + ' -bound_type SLAB'
        spco(command, shell=True)
        template_list.append(opj(dir_prepro, ID + '_acpc_cropped' + type_norm + 'Zp.nii.gz'))

            ###########################template

    if not os.path.exists(study_template_atlas_forlder): os.mkdir(study_template_atlas_forlder)
    templatedir2 = opj(study_template_atlas_forlder, 'studytemplate2_' + type_norm)
    if not os.path.exists(templatedir2): os.mkdir(templatedir2)

    anatomical.anat_to_common_EMB.anats_to_common(
        s_bind,
        afni_sif,
        template_list,
        templatedir2,
        blur_radius_coarse=11,
        convergence=0,
        registration_kind='nonlinear',
        nonlinear_levels=[4,7],
        nonlinear_minimal_patches=[13,9])

