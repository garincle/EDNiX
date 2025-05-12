import os
from fonctions.extract_filename import extract_filename
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
# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

#################################################################################################
####Seed base analysis
#################################################################################################
def clean(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3, RS, nb_run,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(100) + '(function: _100_Data_Clean).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    ext_file = ['nii.gz','json']

    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])

        for direction_results in [dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3]:

            if direction_results == dir_fMRI_Refth_RS_prepro1:
                for ext in ext_file:
                    list_to_remove = [opj(direction_results, 'Mean_Image_RcT_SS_pre' + ext),
                                      opj(direction_results, 'Mean_Image_test' + ext),
                                      opj(direction_results, root_RS + ext),
                                      opj(direction_results, root_RS + '_xd' + ext),
                                      opj(direction_results, root_RS + '_xdt' + ext),
                                      opj(direction_results, root_RS + '_xdt_mean' + ext),
                                      opj(direction_results, root_RS + '_xdtr' + ext),
                                      opj(direction_results, root_RS + '_xdtr_deob' + ext),
                                      opj(direction_results, root_RS + '_xdtrf_2ref' + ext),
                                      opj(direction_results, root_RS + '_xdtrf_mean_preWARP' + ext),
                                      opj(direction_results, root_RS + '_xdtr_mean' + ext),
                                      opj(direction_results, root_RS + '_xdtr_mean_deob' + ext),
                                      opj(direction_results, root_RS + '_xdtr_mean_deob_ref_fudge' + ext),
                                      opj(direction_results, root_RS + '_xdtr_mean_preWARP' + ext)]

                    for remove_data in list_to_remove:
                        if ope(remove_data):
                            os.remove(remove_data)
                            nl = 'INFO: ' + remove_data + ' removed !'
                            print(bcolors.OKGREEN + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')
                        else:
                            nl = 'INFO: ' +remove_data + ' not found'
                            print(bcolors.WARNING + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')


            if direction_results == dir_fMRI_Refth_RS_prepro3:
                for ext in ext_file:
                    list_to_remove = [opj(direction_results, 'Mean_Image_RcT_SS_pre' + ext),
                                      opj(direction_results, root_RS + '_residual_in_anat_sfht' + ext)]

                    for remove_data in list_to_remove:
                        if ope(remove_data):
                            os.remove(remove_data)
                            nl = 'INFO: ' + remove_data + ' removed !'
                            print(bcolors.OKGREEN + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')
                        else:
                            nl = 'INFO: ' + remove_data + ' not found'
                            print(bcolors.WARNING + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')


    diary.write(f'\n')
    diary.close()