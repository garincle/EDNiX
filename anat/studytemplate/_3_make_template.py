##############################################################################
####### CREATE THE STUDY TEMPLATE (IF YOU WANT ON) ###########################
##############################################################################

import os
import subprocess
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

spco = subprocess.check_output
spgo = subprocess.getoutput

from Tools import run_cmd
from Tools import getpath

from anat.studytemplate import anat_to_common_EMB

############# use #############
#######Name: nilearn ##########
#######Version: 0.9.0 #########
###############################

def make_template(which_on, all_ID_max, all_Session_max, all_data_path_max, all_ID, all_Session, all_data_path, type_norm, study_template_atlas_folder,
                  sing_afni,diary_file):

    nl = 'Run anat._3_make_template.make_template'
    run_cmd.msg(nl, diary_file,'HEADER')

    fname        = '_space-acpc_desc-SS-step1_'


    if which_on == 'max': # all or max
        all_ID_temp        = all_ID_max
        all_Session_temp   =  [all_Session_max[all_ID.index(ID)] for ID in all_ID_max]
        all_data_path_temp = all_data_path_max
    elif which_on == 'all':
        all_ID_temp        = all_ID
        all_Session_temp   = all_Session
        all_data_path_temp = all_data_path
    else:
        nl = 'ERROR: with which_on name, need to be "all" or "max"'
        raise Exception(run_cmd.error(nl, diary_file))

    for IDp, Sessionp in zip(all_ID_temp, all_Session_temp):
        run_cmd.printcolor(f'INFO: Subject selected for study template creation: {IDp} {Sessionp}', 'OKGREEN')

    (_, _, _, _, _, volumes_dir_ref, _, _) = getpath.anat(all_data_path_temp[0],'','',0,0,'native')
    dummy_img = opj(volumes_dir_ref, all_ID_temp[0] + fname + type_norm + '.nii.gz')

    zpad_suffix  = 'Zp'
    dist          = 0
    # create a large image for co-registration

    ref_img   = dummy_img.replace('.nii.gz', zpad_suffix + '.nii.gz')

    side = 'I', 'S', 'A', 'P', 'L', 'R', 'S'
    mat = []
    for i in side:
        mat.append(' '.join(['-' + i, str(dist)]))
    zpside = ' '.join(mat)

    command = (sing_afni + '3dZeropad -overwrite ' + zpside + ' -prefix ' + ref_img + ' ' + dummy_img)
    print(command)
    run_cmd.do(command, diary_file)

    dictionary = {"Sources": dummy_img,
                  "Description": 'Size image normalization.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(ref_img.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

    (studyacpc_dir,studyprepro_dir,studytransfo_dir, studyvolume_dir, studylabels_dir, studymasks_dir,
     _,_,_,_,_,_,_,_) = getpath.stytemplate(study_template_atlas_folder,'','')

    template_list = []
    for ID, Session, data_path in zip(all_ID_temp, all_Session_temp, all_data_path_temp):

        (_, _, _, _, _, volumes_dir, _, _) = getpath.anat(data_path, '', '', 0, 0, 'native')

        # resample the other anat to this large img
        command = (sing_afni + '3dresample -overwrite -master ' + ref_img +
                   ' -input '  + opj(volumes_dir, ID + fname + type_norm + '.nii.gz') +
                   ' -prefix ' + opj(volumes_dir, ID + fname + zpad_suffix + '.nii.gz') + ' -bound_type SLAB')
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": opj(volumes_dir, ID + fname + type_norm + '.nii.gz'),
                      "Description": 'Size image normalization.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(volumes_dir, ID + fname + zpad_suffix +'.json'), "w") as outfile:
            outfile.write(json_object)
        template_list.append(opj(volumes_dir, ID + fname + zpad_suffix + '.nii.gz'))

    if not ope(study_template_atlas_folder):
        os.mkdir(study_template_atlas_folder)

    if not ope(opd(studyprepro_dir)):
        os.mkdir(opd(studyprepro_dir))
    if not ope(studyprepro_dir):
        os.mkdir(studyprepro_dir)

    anat_to_common_EMB.anats_to_common(sing_afni, diary_file,
                                       template_list,
                                       studyprepro_dir,
                                       blur_radius_coarse=11,
                                       convergence=0,
                                       registration_kind='nonlinear',
                                       nonlinear_levels=[4,7],
                                       nonlinear_minimal_patches=[13,9])



