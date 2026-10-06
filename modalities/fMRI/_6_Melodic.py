import os
import glob

opj = os.path.join
ope = os.path.exists

from Tools import run_cmd
from modalities.fMRI.extract_filename import extract_filename


def Melodic_correct(dir_RS_ICA_native_PreTT, dir_RS_ICA_native, dir_prepro_orig_process, dir_prepro_orig_masks, ID,
                    nb_ICA_run, nb_run, RS, ICA_cleaning, sing_fsl, sing_itk, TR, TfMRI, diary_file):
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

            fMRI_run_inRef_acpc_SS = opj(dir_prepro_orig_process,
                                         root_RS + '_space-acpc-func_desc-fMRI_run_inRef_SS.nii.gz')

            command = (sing_fsl + 'melodic -i ' + fMRI_run_inRef_acpc_SS +
                       ' -o ' + opj(dir_RS_ICA_native, root_RS + '_residual.ica') +
                       ' --tr=' + str(TR) +
                       ' --bgimage=' + Mean_Image_acpc +
                       ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native,
                                                                                          'report.html') +
                       ' --mask=' + maskDilat_funcspace + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66')

            run_cmd.run(command, diary_file)

            # fsl_regfilt -f expects 1-BASED component numbers (as displayed by MELODIC and
            # ITK-SNAP). The previous version stored nb - 1, so the component just before
            # each selected one was regressed out. Selections are now stored 1-based in
            # melodic_select_1based.txt; an old 0-based melodic_select.txt is converted.
            ica_dir = opj(dir_RS_ICA_native, root_RS + '_residual.ica')
            sel_new = opj(ica_dir, 'melodic_select_1based.txt')
            sel_old = opj(ica_dir, 'melodic_select.txt')
            if ope(sel_new):
                with open(sel_new) as readtxt:
                    readtxt_f = readtxt.read().strip()
                nl = 'found ' + sel_new + ' from a previous run, remove it to select components again'
                run_cmd.msg(nl, diary_file, 'WARNING')
            elif ope(sel_old):
                with open(sel_old) as readtxt:
                    old_sel = [int(v) for v in readtxt.read().strip().split(',') if v.strip() != '']
                readtxt_f = ','.join(str(v + 1) for v in old_sel)
                with open(sel_new, 'w') as fh:
                    fh.write(readtxt_f)
                nl = ('WARNING: 0-based selection ' + sel_old + ' converted to 1-based (' + readtxt_f +
                      '); check it against the MELODIC report.')
                run_cmd.msg(nl, diary_file, 'WARNING')
            else:
                command = (sing_itk + 'itksnap -g ' + anat_in_acpc_func +
                           ' -o ' + opj(ica_dir, 'melodic_IC.nii.gz'))
                run_cmd.run(command, diary_file)
                f = []
                for nb in list(range(1, nb_ICA_run + 1)):
                    while True:
                        Your_choice = input('Is ' + str(nb) + ' an artefact ?')
                        if Your_choice == '' or not Your_choice in ['y', 'n']:
                            print('Please answer with y or n !')
                        else:
                            break
                    if Your_choice == "y":
                        f.append(nb)
                        nl = "to be removed"
                    elif Your_choice == "n":
                        nl = "let's keep it"
                    run_cmd.msg(nl, diary_file, 'OKBLUE')
                readtxt_f = ','.join(str(e) for e in f)
                with open(sel_new, "w") as file:
                    file.write(readtxt_f)
            nl = 'processing fsl regflirt !! Close itksnap please !'
            run_cmd.msg(nl, diary_file, 'OKBLUE')

            # -a switch on aggressive filtering (full instead of partial regression)
            if readtxt_f == '':
                nl = 'No component selected for ' + root_RS + ': data copied unchanged.'
                run_cmd.msg(nl, diary_file, 'WARNING')
                command = (sing_fsl + 'imcp ' + fMRI_run_inRef_acpc_SS + ' ' +
                           opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz'))
            else:
                command = (sing_fsl + 'fsl_regfilt -a -i ' + fMRI_run_inRef_acpc_SS +
                           ' -d ' + opj(ica_dir, 'melodic_mix') +
                           ' -o ' + opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz') + ' -f ' + readtxt_f)
            run_cmd.run(command, diary_file)

        DATA = sorted(glob.glob(opj(dir_RS_ICA_native, '*_norm_final_clean.nii.gz')))
        name = ''
        DATA2 = ''
        for m in range(0, len(DATA)):
            name = name + ',' + DATA[m]
            DATA2 = DATA2 + ' ' + DATA[m]

        command = (sing_fsl + 'melodic -i ' + name[1:] +
                   ' -o ' + opj(dir_RS_ICA_native_PreTT, root_RS + '_residual.ica') +
                   ' --tr=' + str(TR) +
                   ' --bgimage=' + Mean_Image_acpc +
                   ' --nobet --bgthreshold=10 -a symm -v --report --guireport=' + opj(dir_RS_ICA_native_PreTT,
                                                                                      'report.html') +
                   ' --mask=' + maskDilat_funcspace + ' --Oall -d ' + str(nb_ICA_run) + ' --mmthresh=0.66')
        run_cmd.run(command, diary_file)

    elif ICA_cleaning == 'Skip':
        nl = "as requested do not perform ICA cleaning"
        run_cmd.msg(nl, diary_file, 'OKGREEN')
    else:
        nl = "ERROR: ICA_cleaning must be str, and one of Skip, MELODIC or ICA-ARMOA"
        raise Exception(run_cmd.error(nl, diary_file))
