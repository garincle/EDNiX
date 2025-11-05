import os
import glob

opj = os.path.join
ope = os.path.exists

from Tools import run_cmd
from fonctions.extract_filename import extract_filename

def Melodic_correct(dir_RS_ICA_native_PreTT, dir_RS_ICA_native, dir_prepro_orig_process, dir_prepro_orig_masks, ID,
                    nb_ICA_run, nb_run, RS, ICA_cleaning, sing_fsl,sing_itk,TR, TfMRI, diary_file):

    nl = '##  Working on step ' + str(6) + '(function: _6_Melodic).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')
    
    #   6.0 first pass : 'individual runs' ICA analysis
    if not ope(dir_RS_ICA_native_PreTT):
        os.mkdir(dir_RS_ICA_native_PreTT)
    if not ope(dir_RS_ICA_native):
        os.mkdir(dir_RS_ICA_native)
        
    maskDilat_funcspace = opj(dir_prepro_orig_masks, ID + '_space-acpc-func_desc-fMRI_mask_dilated.nii.gz')
    Mean_Image_acpc = opj(dir_prepro_orig_process, 'all_runs_space-acpc-func_desc-fMRI_Mean_Image_SS.nii.gz')
    anat_in_acpc_func = opj(dir_prepro_orig_process, ('_').join(['anat_space-acpc-func', TfMRI + '.nii.gz']))
    
    if ICA_cleaning == 'MELODIC':
        # Normalized: use fsl_regfilt to remove noises.
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])

            fMRI_run_inRef_acpc_SS = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_SS.nii.gz')

            command = (sing_fsl + 'melodic -i ' + fMRI_run_inRef_acpc_SS +
                       ' -o ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica') +
                       ' --tr=' + str(TR) +
                       ' --bgimage=' + Mean_Image_acpc +
                       ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native,'report.html') +
                       ' --mask=' + maskDilat_funcspace + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66')

            run_cmd.run(command, diary_file)

            if ope(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt'))==True:
                with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica'), 'melodic_select.txt') as readtxt:
                    readtxt_f = readtxt.read()
                nl = 'found txt file of the previous ICA, remove it if you want to create a new one!'
                run_cmd.msg(nl, diary_file, 'WARNING')
            else:
                command = (sing_itk + 'itksnap -g ' +  anat_in_acpc_func +
                           ' -o ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_IC.nii.gz'))
                run_cmd.run(command, diary_file)
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
                        nl = "to be removed"
                    elif Your_choice=="n":
                        nl = "let's keep it"
                    run_cmd.msg(nl, diary_file, 'OKBLUE')

                f = ','.join(str(e) for e in f)
                # open text file
                with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt'), "w") as file:
                    file.write(f)
                with open(opj(dir_RS_ICA_native, root_RS + '_residual.ica', 'melodic_select.txt')) as readtxt:
                    readtxt_f = readtxt.read()

            nl = 'processing fsl regflirt !! Close itksnap please !'
            run_cmd.msg(nl, diary_file, 'OKBLUE')

            #-a switch on aggressive filtering (full instead of partial regression)
            command = (sing_fsl + 'fsl_regfilt -a -i ' + fMRI_run_inRef_acpc_SS +
                       ' -d ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica','melodic_mix') +
                       ' -o ' + opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz') + ' -f ' + str(readtxt_f))
            run_cmd.run(command, diary_file)

        DATA = sorted(glob.glob(opj(dir_RS_ICA_native, root_RS + '*_norm_final_clean.nii.gz')))
        name  = ''
        DATA2 = ''
        for m in range(0,len(DATA)):
            name  = name  + ',' + DATA[m]
            DATA2 = DATA2 + ' ' + DATA[m]

        command = (sing_fsl + 'melodic -i ' + name[1:] +
                   ' -o ' + opj(dir_RS_ICA_native_PreTT, root_RS + '_residual.ica') +
                   ' --tr=' + str(TR) +
                   ' --bgimage=' + Mean_Image_acpc +
                   ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native_PreTT,'report.html') +
                   ' --mask=' + maskDilat_funcspace + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66')
        run_cmd.run(command, diary_file)

    elif ICA_cleaning == 'Skip':
        nl = "as requested do not perform ICA cleaning"
        run_cmd.msg(nl, diary_file, 'OKGREEN')
    else:
        nl = "ERROR: ICA_cleaning must be str, and one of Skip, MELODIC or ICA-ARMOA"
        raise Exception(run_cmd.error(nl, diary_file))
