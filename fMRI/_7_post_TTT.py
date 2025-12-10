import os
import subprocess
import json
import ants
import pandas as pd

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

spgo = subprocess.getoutput

from Tools import run_cmd
from fMRI.extract_filename import extract_filename

def signal_regression(dir_prepro_orig_process, dir_RS_ICA_native, dir_prepro_orig_masks, dir_prepro_raw_process,
    nb_run, RS, blur, TR, ICA_cleaning, extract_exterior_CSF, extract_WM, normalize, ID, post_treatment_method,
    do_not_correct_signal, band, extract_Vc, extract_GS, dir_prepro_orig_postprocessed, dir_prepro_raw_matrices, overwrite, sing_afni, sing_fsl, diary_file):

    nl = '##  Working on step ' + str(7) + '(function: _7_post_TTT).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    for i in range(int(nb_run)):
        root_RS = extract_filename(RS[i])

        if ICA_cleaning == 'Skip':
            input = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_SS.nii.gz')
        else:
            input = opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz')

        for extract_type, img_name ,suffix in zip([extract_exterior_CSF, extract_WM, extract_GS],
                                                  ['exterior_ligne.nii.gz', 'Wmask.nii.gz', 'maskDilat.nii.gz'],
                                                  ['_NonB','_Wc','-GS']):
            if extract_type == True:
                command = (sing_afni + '3dmaskSVD' + overwrite + ' -polort 2 -vnorm -mask ' +
                           opj(dir_prepro_orig_masks, img_name) +
                           ' ' + input)
                val, _ = run_cmd.get(command, diary_file)
                list1D = val.decode("utf-8").split('\n')
                Bp1d = open(opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.1D'), 'w')
                for i in range(len(list1D)):
                    line = ' '.join(list1D[i].split()) + '\n'
                    Bp1d.write(line)
                Bp1d.close()

        if extract_Vc == True:
                msk = ants.image_read(opj(dir_prepro_orig_masks,'Vmask.nii.gz'))
                countVmask = int(msk.sum())

                if countVmask>10:
                    command = (sing_afni + '3dmaskSVD' + overwrite + ' -polort 2 -vnorm -mask ' + opj(
                        dir_prepro_orig_masks, 'Vmask.nii.gz') +
                               ' ' + input)
                    val, _ = run_cmd.get(command, diary_file)
                    list1D = val.decode("utf-8").split('\n')
                    Bp1d = open(opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_Vc.1D'), 'w')
                    for i in range(len(list1D)):
                        line = ' '.join(list1D[i].split()) + '\n'
                        Bp1d.write(line)
                    Bp1d.close()

        for option_type, suffix,descript in zip([' -cvarinv',' -cvar', ' -tsnr',' -stdev'],
                                                ['_tsnr1','_cvar','_tsnr2','_stdev'],
                                                ['fabs(mean)/stdev (with detrend)',
                                                 'stdev/fabs(mean) (with detrend)',
                                                 'fabs(mean)/stdev NOT DETRENDED',
                                                 'standard deviation with detrend']):
            command = (sing_afni + '3dTstat' + overwrite + option_type +
                       ' -prefix ' + opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.nii.gz') +
                       ' ' + input)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": input,
                          "Description": descript + ' (3dTstat, AFNI)'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.json'), "w") as outfile:
                outfile.write(json_object)

        ############################### ############################### ###############################
        ##                                Corrections of the signal                                  ##
        ############################### ############################### ###############################

        # 5.0 Regress out most of the noise from the data: bandpass filter, motion correction white mater noise and cbf noise , plus drift and derivatives
        # after filtering : blur within the mask and normalise the data.


        # create bandpass regressors (instead of using 3dBandpass, say)
        hd = ants.image_header_info(input)

        command = (sing_afni + '1dBport' + overwrite + ' -nodata ' + str(int(hd['dimensions'][3])) +
                   ' ' + str(TR) + ' -band ' + band + ' -invert -nozero')
        Bp, _ = run_cmd.get(command, diary_file)
        list1D = Bp.decode("utf-8").split('\n')
        Bp1d = open(opj(dir_prepro_orig_postprocessed, root_RS + 'bandpass_rall.1D'), 'w')
        for i in range(len(list1D)):
            line = ' '.join(list1D[i].split()) + '\n'
            Bp1d.write(line)
        Bp1d.close()

        if ope(dir_RS_ICA_native) == False:
            os.makedirs(dir_RS_ICA_native)

    if do_not_correct_signal == False:
        for i in range(int(nb_run)):
            try:
                root_RS = extract_filename(RS[i])
                original_dir = os.getcwd()
                os.chdir(dir_prepro_orig_process)
                
                if ICA_cleaning == 'Skip':
                    input = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_run_inRef.nii.gz')
                else:
                    input = opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz')

                Mean = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_Mean.nii.gz')
                Sdev = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_Sdev.nii.gz')
                maskDilat_funcspace = opj(dir_prepro_orig_masks, ID + '_space-acpc-func_desc-fMRI_mask_dilated.nii.gz')
                bandpass = opj(dir_prepro_orig_postprocessed, root_RS + 'bandpass_rall.1D')
                censore1D = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-censor.1D')
                demean = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-demean.1D')
                deriv = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-deriv.1D')
                residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')

                if post_treatment_method == 'AFNI':
                    command = (sing_afni + '3dDeconvolve -input ' + input +
                               ' -mask ' + maskDilat_funcspace +
                               ' -ortvec ' + bandpass + ' bandpass_rall' +
                               ' -ortvec ' + demean + ' mot_demean' +
                               ' -ortvec ' + deriv + ' mot_deriv' +
                               ' -censor ' + censore1D +
                               ' -polort A -float' +
                               ' -num_stimts 0' + overwrite +
                               ' -fout -tout' +
                               ' -x1D ' + opj(dir_prepro_orig_postprocessed, root_RS + 'X.xmat.1D') +
                               ' -xjpeg ' + opj(dir_prepro_orig_postprocessed, root_RS + 'X.jpg') +
                               ' -x1D_uncensored ' + opj(dir_prepro_orig_postprocessed, root_RS + 'X.nocensor.xmat.1D') +
                               ' -fitts ' + opj(dir_prepro_orig_postprocessed, root_RS + 'Xfittssubj') +
                               ' -errts ' + opj(dir_prepro_orig_postprocessed, root_RS + 'errts') +
                               ' -x1D_stop' +
                               ' -bucket ' + opj(dir_prepro_orig_postprocessed, root_RS + 'statssubj'))

                    for extract_type, suffix in zip([extract_exterior_CSF, extract_WM, extract_GS, extract_Vc],
                                                    ['_NonB', '_Wc', '-GS','-Vc']):
                        if extract_type == True:
                            command = command + ' -ortvec ' + opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.1D') + ' residual_norm' + suffix + ' '

                    if ope(opj(dir_prepro_orig_postprocessed, root_RS + '_confounds_correct.tsv')):
                        confounds_df = pd.read_csv(opj(dir_prepro_orig_postprocessed, root_RS + '_confounds_correct.tsv'), sep='\t')

                        for column in confounds_df.columns:
                            column_file = opj(dir_prepro_orig_postprocessed, root_RS + f"_{column}.1D")
                            confounds_df[[column]].to_csv(column_file, sep=' ', index=False, header=False)

                            # Append the column file to the command
                            command += ' -ortvec ' +  column_file + ' residual_norm_' + column

                    nl = 'INFO: 3dDeconvolve command is ' + command
                    run_cmd.msg(nl, diary_file, 'OKGREEN')
                    #run_cmd.run(command, diary_file)
                    subprocess.run(command,
                        shell=True, check=True)

                    command = (sing_afni + '3dTproject -polort 0' + overwrite + ' -input ' +
                               input +
                               ' -censor ' + censore1D +
                               ' -cenmode ZERO -ort ' + opj(dir_prepro_orig_postprocessed, root_RS + 'X.nocensor.xmat.1D') +
                               ' -prefix ' + residual)
                    if blur > 0:
                        command += ' -blur ' + str(blur)
                    run_cmd.run(command, diary_file)

                    # Normalization options
                    zscore_command = (sing_afni + '3dTstat -mean -overwrite' +
                                      ' -prefix ' + Mean +
                                      ' ' + residual +
                                      ' && 3dTstat -stdev -prefix ' + Sdev +
                                      ' ' + residual +
                                      ' && 3dcalc -overwrite -a ' + residual +
                                      ' -b ' + Mean +
                                      ' -c ' + Sdev +
                                      ' -expr "(a-b)/c" -prefix ' + residual)

                    psc_command = (sing_afni + '3dTstat -mean -overwrite -prefix ' + Mean +
                                   ' ' + residual +
                                   ' && 3dcalc -overwrite -a ' + residual +
                                   ' -b ' + Mean +
                                   ' -expr "((a - b) / b) * 100" -prefix ' + residual)

                    # Choose normalization method
                    if normalize == 'zscore':
                        run_cmd.do(zscore_command, diary_file)
                    elif normalize == 'psc':
                        run_cmd.do(psc_command, diary_file)
                    else:
                        run_cmd.msg('no normalization',diary_file, 'OKCYAN')
                elif post_treatment_method == 'Grandjean':
                    # --- Définition des chemins principaux ---
                    regressors_file = opj(dir_prepro_orig_postprocessed, 'regressors.1D')
                    matrix_motion_correction = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-motion_correction.1D')
                    wm_signal = opj(dir_prepro_orig_postprocessed, 'wm.1D')
                    csf_signal = opj(dir_prepro_orig_postprocessed, 'csf.1D')

                    wm_mask = opj(dir_prepro_orig_masks, 'Wmask.nii.gz')
                    csf_mask = opj(dir_prepro_orig_masks, 'Vmask.nii.gz')

                    if ope(wm_mask) == True and ope(csf_mask) == True:
                        run_cmd.msg("WM and CSF masks found", diary_file, 'OKCYAN')
                        # --- Extraction des signaux de référence (WM et CSF) ---
                        run_cmd.msg("Extracting WM signal", diary_file, 'OKCYAN')
                        cmd = f"{sing_fsl}fslmeants -i {input} -o {wm_signal} -m {wm_mask}"
                        run_cmd.run(cmd, diary_file)

                        run_cmd.msg("Extracting CSF signal", diary_file, 'OKCYAN')
                        cmd = f"{sing_fsl}fslmeants -i {input} -o {csf_signal} -m {csf_mask}"
                        run_cmd.run(cmd, diary_file)

                        # --- Combinaison des régressors (motion + WM + CSF) ---
                        run_cmd.msg("Combining regressors (motion + WM + CSF)", diary_file, 'OKCYAN')

                        subprocess.run(f"{sing_afni}1dcat {matrix_motion_correction} {wm_signal} {csf_signal} > {regressors_file}",
                            shell=True, check=True)

                    elif ope(wm_mask) == True and ope(csf_mask) == False:
                        run_cmd.msg("Only WM mask found", diary_file, 'OKCYAN')
                        # --- Extraction des signaux de référence (WM) ---
                        run_cmd.msg("Extracting WM signal", diary_file, 'OKCYAN')
                        cmd = f"{sing_fsl}fslmeants -i {input} -o {wm_signal} -m {wm_mask}"
                        run_cmd.run(cmd, diary_file)

                        # --- Combinaison des régressors (motion + WM) ---
                        run_cmd.msg("Combining regressors (motion + WM)", diary_file, 'OKCYAN')

                        subprocess.run(f"{sing_afni}1dcat {matrix_motion_correction} {wm_signal} > {regressors_file}",
                            shell=True, check=True)

                    elif ope(wm_mask) == False and ope(csf_mask) == True:
                        run_cmd.msg("Only CSF mask found", diary_file, 'OKCYAN')
                        # --- Extraction des signaux de référence (CSF) ---
                        run_cmd.msg("Extracting CSF signal", diary_file, 'OKCYAN')
                        cmd = f"{sing_fsl}fslmeants -i {input} -o {csf_signal} -m {csf_mask}"
                        run_cmd.run(cmd, diary_file)

                        # --- Combinaison des régressors (motion + CSF) ---
                        run_cmd.msg("Combining regressors (motion + CSF)", diary_file, 'OKCYAN')

                        subprocess.run(f"{sing_afni}1dcat {matrix_motion_correction} {csf_signal} > {regressors_file}",
                            shell=True, check=True)
                    else:
                        run_cmd.msg("No WM or CSF masks found, using only motion regressors", diary_file, 'OKCYAN')
                        # --- Utilisation uniquement des régressors de mouvement ---
                        subprocess.run(f"{sing_afni}1dcat {matrix_motion_correction} > {regressors_file}",
                            shell=True, check=True)
                    # --- Application du bandpass + régression des confounds ---
                    fwhm_str = str(blur) if blur > 0 else "0"
                    prefiltered = opj(dir_prepro_orig_postprocessed, root_RS + '_prefiltered_func_data_tempfilt.nii.gz')

                    run_cmd.msg("Running 3dTproject (bandpass + regressors removal)", diary_file, 'OKCYAN')
                    cmd = (
                        f"{sing_afni}3dTproject "
                        f"-blur {fwhm_str} "
                        f"-passband {band} "
                        f"-ort {regressors_file} "
                        f"-prefix {prefiltered} "
                        f"-input {input} "
                        f"{overwrite}")
                    run_cmd.run(cmd, diary_file)

                    # --- Restauration du niveau de base (Grandjean) ---
                    run_cmd.msg("Restoring mean intensity with fslmaths", diary_file, 'OKCYAN')
                    cmd = f"{sing_fsl}fslmaths {input} -Tmean {Mean}"
                    run_cmd.run(cmd, diary_file)
                    cmd = f"{sing_fsl}fslmaths {prefiltered} -add {Mean} {residual}"
                    run_cmd.run(cmd, diary_file)
                else:
                    print("post treatement method name: " + str(post_treatment_method) + " do not exists")
            except Exception as e:
                print('3dDeconvolve failed for run: ' + str(RS[i]) + ' with error: ' + str(e))
                root_RS = extract_filename(RS[i])
                Deconvolve_failed = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-3dDeconvolve_failed.nii.gz')
                command = (sing_afni + '3dcalc' + overwrite +
                           ' -a ' + input +
                           ' -prefix ' + Deconvolve_failed +
                           ' -expr "a"')
                run_cmd.do(command, diary_file)

                dictionary = {"Sources": input,
                              "Description": 'Copy'},
                json_object = json.dumps(dictionary, indent=2)
                with open(Deconvolve_failed.replace('.nii.gz','.json'), "w") as outfile:
                    outfile.write(json_object)
                os.chdir(original_dir)

            if not ope(residual):
                print('path to residual do not exists: ope(residual)=' + str(ope(residual)))
                root_RS = extract_filename(RS[i])
                command = (sing_afni + '3dcalc' + overwrite +
                           ' -a ' + input +
                           ' -prefix ' + Deconvolve_failed + ' -expr "a"')
                run_cmd.do(command, diary_file)

                dictionary = {"Sources": input,
                              "Description": 'Copy'},
                json_object = json.dumps(dictionary, indent=2)
                with open(Deconvolve_failed.replace('.nii.gz','.json'), "w") as outfile:
                    outfile.write(json_object)
                os.chdir(original_dir)
    else:
        for i in range(int(nb_run)):
            root_RS = extract_filename(RS[i])
            if ICA_cleaning == 'Skip':
                input = opj(dir_prepro_orig_process, root_RS + '_space-acpc-func_desc-fMRI_run_inRef.nii.gz')
            else:
                input = opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz')
            residual = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')

            command = (sing_afni + '3dcalc' + overwrite +
                       ' -a ' + input +
                       ' -prefix ' + residual + ' -expr "a"')
            run_cmd.do(command, diary_file)
            dictionary = {"Sources": input,
                          "Description": 'Copy'},
            json_object = json.dumps(dictionary, indent=2)
            with open(residual.replace('.nii.gz','.json'), "w") as outfile:
                outfile.write(json_object)

