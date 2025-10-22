import os
import numpy as np
import ants
import json


opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from Tools import getpath
from Tools import create1Dmatrix

def allineate(data_path,anatinput,anatoutput,maskinput,BASE_SS,cost3dAllineate,sing_afni,overwrite,diary_file):

    _, _, _, dir_prepro, _, _, _, masks_dir, = getpath.anat(data_path,'','',False,False,'native')

    command = (sing_afni + '3dAllineate' + overwrite + ' -warp shift_rotate -cmass -base ' + BASE_SS +
               ' -cost ' + cost3dAllineate + ' -interp NN' +
               ' -source_mask ' + maskinput +
               ' -prefix '      + anatoutput +
               ' -source '      + anatinput +
               ' -1Dmatrix_save ' + anatinput.replace('.nii.gz','.1D'))
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [anatinput,
                              BASE_SS,
                              maskinput],
                  "Description": 'rigid co-registration (AFNI, 3dAllineate).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatoutput.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)

def aligncenter(data_path,anatinput,anatoutput,BASE_SS,sing_afni,overwrite,diary_file):

    _, _, _, dir_prepro, _, _, _, masks_dir, = getpath.anat(data_path,'','',False,False,'native')

    current_working_directory = os.getcwd()

    os.chdir(dir_prepro)

    command = (sing_afni + '@Align_Centers -base ' + BASE_SS +
               ' -dset ' + anatinput + ' -cm -prefix ' + opb(anatinput) + overwrite)
    run_cmd.run(command, diary_file)

    os.chdir(str(current_working_directory))

    command = (sing_afni + '3dcopy ' + anatinput + ' ' + anatoutput + overwrite)
    run_cmd.run(command, diary_file)

    dictionary = {"Sources": [anatinput,
                              BASE_SS],
                  "Description": 'rigid co-registration (AFNI, @Align_Centers).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatoutput.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)


def none(anatinput,anatoutput,overwrite,diary_file,sing_afni):

    command = (sing_afni + '3dcopy ' + anatinput + ' ' + anatoutput + overwrite)
    run_cmd.run(command, diary_file)

    filename = anatoutput.replace('.nii.gz','.1D')
    matrix = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0]])

    create1Dmatrix(filename, matrix,diary_file)
    nl = f"INFO: .1D matrix file saved as {filename}"
    run_cmd.msg(nl, diary_file, 'OKGREEN')

    dictionary = {"Sources": anatinput,
                  "Description": 'Copy.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(anatoutput.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)







