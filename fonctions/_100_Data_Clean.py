import shutil
import subprocess
import os
from fonctions.extract_filename import extract_filename
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
# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname
ope = os.path.exists

#################################################################################################
####Seed base analysis
#################################################################################################
def clean(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, RS, nb_run):
    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])

        for direction_results in [dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3]:
            if direction_results == dir_fMRI_Refth_RS_prepro1:

                list_to_remove = [opj(direction_results, 'Mean_Image_RcT_SS_pre.nii.gz'),
                opj(direction_results, 'Mean_Image_test.nii.gz'),
                opj(direction_results, root_RS + '.nii.gz'),
                opj(direction_results, root_RS + '_xd.nii.gz'),
                opj(direction_results, root_RS + '_xdt.nii.gz'),
                opj(direction_results, root_RS + '_xdt_mean.nii.gz'),
                opj(direction_results, root_RS + '_xdtr.nii.gz'),
                opj(direction_results, root_RS + '_xdtr_deob.nii.gz'),
                opj(direction_results, root_RS + '_xdtrf_2ref.nii.gz'),
                opj(direction_results, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz'),
                opj(direction_results, root_RS + '_xdtrf_mean_preWARP.nii.gz'),
                opj(direction_results, root_RS + '_xdtr_mean.nii.gz'),
                opj(direction_results, root_RS + '_xdtr_mean_deob.nii.gz'),
                opj(direction_results, root_RS + '_xdtr_mean_deob_ref_fudge.nii.gz'),
                opj(direction_results, root_RS + '_xdtr_mean_preWARP.nii.gz')]

                for remove_data in list_to_remove:
                    if ope(remove_data):
                        os.remove(remove_data)
                        print(bcolors.OKGREEN + 'INFO: ' + remove_data + ' removed !' + bcolors.ENDC)
                    else:
                        print(bcolors.OKGREEN + 'INFO: ' +remove_data + ' not found' + bcolors.ENDC)

                if ope(opj(direction_results, 'tmp')):
                    shutil.rmtree(opj(direction_results, 'tmp'))
                    print(bcolors.OKGREEN + 'INFO: ' + opj(direction_results, 'tmp')+ ' removed !' + bcolors.ENDC)
                else:
                    print(bcolors.OKGREEN + 'INFO: ' + opj(direction_results, 'tmp') + ' not found' + bcolors.ENDC)

            if direction_results == dir_fMRI_Refth_RS_prepro2:
                if ope(opj(direction_results, 'tmp')):
                    shutil.rmtree(opj(direction_results, 'tmp'))
                    print(bcolors.OKGREEN + 'INFO: ' + opj(direction_results, 'tmp')+ ' removed !' + bcolors.ENDC)
                else:
                    print(bcolors.OKGREEN + 'INFO: ' + opj(direction_results, 'tmp') + ' not found' + bcolors.ENDC)


            if direction_results == dir_fMRI_Refth_RS_prepro3:

                list_to_remove = [opj(direction_results, 'Mean_Image_RcT_SS_pre.nii.gz'),
                                  opj(direction_results, root_RS + '_residual_in_anat_sfht.nii.gz')]

                for remove_data in list_to_remove:
                    if ope(remove_data):
                        os.remove(remove_data)
                        print(bcolors.OKGREEN + 'INFO: ' + remove_data + ' removed !' + bcolors.ENDC)
                    else:
                        print(bcolors.OKGREEN + 'INFO: ' + remove_data + ' not found' + bcolors.ENDC)

                if ope(opj(direction_results, 'tmp')):
                    shutil.rmtree(opj(direction_results, 'tmp'))
                    print(bcolors.OKGREEN + 'INFO: ' + opj(direction_results, 'tmp') + ' removed !' + bcolors.ENDC)
                else:
                    print(bcolors.OKGREEN + 'INFO: ' + opj(direction_results, 'tmp') + ' not found' + bcolors.ENDC)
