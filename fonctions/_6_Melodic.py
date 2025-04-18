import os
import subprocess
import glob
import datetime
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


opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput


def Melodic_correct(dir_RS_ICA_native_PreTT, dir_RS_ICA_native, dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2,
                    nb_ICA_run, nb_run, RS, ICA_cleaning, MAIN_PATH,s_bind,fsl_sif,itk_sif,TR, diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(6) + '(function: _6_Melodic).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    #   6.0 first pass : 'individual runs' ICA analysis
    if not ope(dir_RS_ICA_native_PreTT):
        os.mkdir(dir_RS_ICA_native_PreTT)
    if not ope(dir_RS_ICA_native):
        os.mkdir(dir_RS_ICA_native)

    if ICA_cleaning == 'MELODIC':
        # Normalized: use fsl_regfilt to remove noises.
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = 'singularity run' + s_bind + fsl_sif + 'melodic -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + \
            ' -o ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica') + \
            ' --tr=' + str(TR) + \
            ' --bgimage=' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
            ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native,'report.html') + \
            ' --mask=' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66'

            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

            if ope(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt'))==True:
                with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica'), 'melodic_select.txt') as readtxt:
                    readtxt_f = readtxt.read()
                nl = 'found txt file of the previous ICA, remove it if you want to create a new one!'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')


            else:
                command = 'gnome-terminal -- singularity run' + s_bind + itk_sif + 'itksnap -g ' +  opj(dir_fMRI_Refth_RS_prepro2,'anat_rsp_in_func.nii.gz') + \
                ' -o ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_IC.nii.gz')
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                f = []
                for nb in list(range(1, nb_ICA_run+1)):
                    while True:
                        Your_choice = input('Is ' + str(nb) + ' an artefact ?')
                        if Your_choice == '' or not Your_choice in ['y','n']:
                            print('Please answer with y or n !')
                        else:
                            break
                    if Your_choice=="y":
                        f.append(nb-1)
                        print(bcolors.OKBLUE + "to be removed" + bcolors.ENDC)
                    elif Your_choice=="n":
                        print(bcolors.OKBLUE + "let's keep it" + bcolors.ENDC)
                f = ','.join(str(e) for e in f)
                # open text file
                with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt'), "w") as file:
                    file.write(f)
                with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt')) as readtxt:
                    readtxt_f = readtxt.read()

            nl = 'processing fsl regflirt !! Close itksnap please !'
            print(bcolors.OKBLUE + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            #-a switch on aggressive filtering (full instead of partial regression)
            command = 'singularity run' + s_bind + fsl_sif + 'fsl_regfilt -a -i ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz') + \
            ' -d ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica','melodic_mix') + \
            ' -o ' + opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz') + ' -f ' + str(readtxt_f)
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)

        DATA = sorted(glob.glob(opj(dir_RS_ICA_native, root_RS + '*_norm_final_clean.nii.gz')))
        name  = ''
        DATA2 = ''
        for m in range(0,len(DATA)):
            name  = name  + ',' + DATA[m]
            DATA2 = DATA2 + ' ' + DATA[m]

        command = 'singularity run' + s_bind + fsl_sif + 'melodic -i ' + name[1:] + \
        ' -o ' + opj(dir_RS_ICA_native_PreTT, root_RS + '_residual.ica') + \
        ' --tr=' + str(TR) + \
        ' --bgimage=' + opj(dir_fMRI_Refth_RS_prepro1,'Mean_Image.nii.gz') + \
        ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native_PreTT,'report.html') + \
        ' --mask=' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        diary.write(f'\n')
        diary.close()

    elif ICA_cleaning == 'Skip':
        print(bcolors.OKGREEN + "as requested do not perform ICA cleaning" + bcolors.ENDC)

    else:
        nl = "ERROR: ICA_cleaning must be str, and one of Skip, MELODIC or ICA-ARMOA"
        diary.write(f'\n{nl}')
        raise Exception(bcolors.FAIL + nl + bcolors.ENDC)
