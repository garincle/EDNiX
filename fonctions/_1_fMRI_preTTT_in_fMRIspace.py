import os
import ants
import json
import numpy as np
import nibabel as nib
import nilearn

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,diaryfile, check_nii
from fonctions.extract_filename import extract_filename

##### XXX add ICA or DL to visualize pre-processing effect

def preprocess_data(dir_prepro_orig_process, RS, list_RS, nb_run, T1_eq, TR, Slice_timing_info,
                    overwrite, sing_afni, diary_file,diary_WARNING):

    nl = '##  Working on step ' + str(1) + '(function: _1_fMRI_preTTT_in_fMRIspace).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    for i in range(0, int(nb_run)):

        nl = 'work on ' + str(dir_prepro_orig_process) + ' run ' + str(i +1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        root = extract_filename(RS[i])
        raw_func = list_RS[i]
        base_fMRI_targeted = opj(dir_prepro_orig_process, root + '_space-func_desc-vol_rmv_from_txt.nii.gz')
        base_fMRI = opj(dir_prepro_orig_process, root + '_space-func_desc-vol_rmv.nii.gz')
        fMRI_despike = opj(dir_prepro_orig_process, root + '_space-func_desc-despiked.nii.gz')
        fMRI_SliceT = opj(dir_prepro_orig_process, root + '_space-func_desc-despiked.nii.gz')
        fMRI_runMean = opj(dir_prepro_orig_process, root + '_space-func_desc-runMean.nii.gz')
        fMRI_stc = opj(dir_prepro_orig_process, 'stc.txt')
        fMRI_outcount = opj(dir_prepro_orig_process, root + '_space-func_desc-outcount.r$run.1D')
        file_volreg = opj(dir_prepro_orig_process, root + '_volreg_dfile.1D')
        fMRI_volreg = opj(dir_prepro_orig_process,'_space-func_desc-volreg.nii.gz')
        matrix_volreg = opj(dir_prepro_orig_process, root + 'volreg.aff12.1D')
        censore1D = opj(dir_prepro_orig_process, root + '_space-func_desc-censor.1D')
        censoretxt = opj(dir_prepro_orig_process, root + '_space-func_desc-censor.txt')
        motion_enorm = opj(dir_prepro_orig_process, root + '_space-func_desc-motion_enorm.1D')
        deriv = opj(dir_prepro_orig_process, root + '_space-func_desc-deriv.1D')
        fMRI_runMean_align = opj(dir_prepro_orig_process, root + '_space-func_desc-runMean_align.nii.gz')
        fMRI_runMean_n4Bias = opj(dir_prepro_orig_process, root + '_space-func_desc-runMean_n4Bias.nii.gz')

        # Clean bad volumes
        if ope(opj(opd(list_RS[i]), root, '.txt')) == True:
            # Open the file in read mode
            with open(opj(opd(list_RS[i]), root, '.txt'), 'r') as file:
                # Read the first line and convert it to an integer
                cut_low = int(file.readline().strip())

                # Read the second line and convert it to an integer
                cut_high = int(file.readline().strip())

            command = (sing_afni + '3dTcat -prefix ' + base_fMRI_targeted +
                       ' ' + raw_func + '[' + str(cut_low) + '-' + str(cut_high - 1) + ']' + overwrite)
            run_cmd.do(command, diary_file)

            dictionary = {"Sources": [base_fMRI_targeted,
                                      opj(opd(raw_func), root, '.txt')],
                          "Description": 'Remove volumes.',
                          "Command": command,}
            json_object = json.dumps(dictionary, indent=2)
            with open(base_fMRI_targeted.replace('.nii.gz','json'), "w") as outfile:
                outfile.write(json_object)
            raw_func = base_fMRI_targeted

        img = nib.load(raw_func)
        nb_vol = img.shape[-1]  # For 4D data, last dimension is time

        command = (sing_afni + '3dTcat -prefix ' + base_fMRI +
                   ' ' + raw_func + '[' + str(T1_eq) + '-' + str(nb_vol - 1) + ']' + overwrite)
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": base_fMRI,
                      "Description": 'Remove first volumes.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(base_fMRI.replace('.nii.gz','json'), "w") as outfile:
            outfile.write(json_object)

        # Despiking
        command = (sing_afni + '3dDespike -NEW -nomask' + overwrite + ' -prefix ' +
                   fMRI_despike +
                   ' ' + base_fMRI)
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": base_fMRI,
                      "Description": 'Despiking.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_despike.replace('.nii.gz','json'), "w") as outfile:
            outfile.write(json_object)

        # slice-timing correction -heptic!!!!!!
        if Slice_timing_info == 'Auto':
            if opi(fMRI_stc):
                command = (sing_afni + '3dTshift -wsinc9 ' + overwrite +
                           ' -TR ' + str(TR) + ' -tpattern @' + fMRI_stc +
                           ' -prefix ' + fMRI_SliceT +
                           ' ' + fMRI_despike)
                run_cmd.run(command, diary_file)

                dictionary = {"Sources": fMRI_despike,
                              "Description": 'Slice timing correction.',
                              "Command": command,}
                json_object = json.dumps(dictionary, indent=3)
                with open(fMRI_SliceT.replace('.nii.gz','json'), "w") as outfile:
                    outfile.write(json_object)

            else:
                cmd = ('export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";' + sing_afni + '3dinfo -slice_timing ' + base_fMRI)
                nl,_= run_cmd.get(cmd,diary_file)
                info = nl.decode("utf-8").split('\n')[0]
                STC = info.split('|')

                STC = list(map(float, STC))
                if np.sum(STC) > 0:
                    # means that AFNI has access to the slice timing in the nifti header

                    command = (sing_afni + '3dTshift -wsinc9 ' + overwrite +
                               ' -prefix ' + fMRI_SliceT +
                               ' ' + fMRI_despike)
                    run_cmd.run(command, diary_file)

                    dictionary = {"Sources": fMRI_despike,
                                  "Description": 'Slice timing correction.',
                                  "Command": command,}
                    json_object = json.dumps(dictionary, indent=3)
                    with open(fMRI_SliceT.replace('.nii.gz','json'), "w") as outfile:
                        outfile.write(json_object)
                else:
                    nl = "WARNING: Slice Timing not found, this will be particularly DANGEROUS, you SHOULD PROVIDE MANUALLY ONE!"
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    diary_WARNING_file = diaryfile.create(diary_WARNING,nl)

                    command = (sing_afni + '3dcalc -a ' + fMRI_despike +
                               ' -prefix ' + fMRI_SliceT +
                               ' -expr "a"' + overwrite)
                    run_cmd.do(command, diary_file)

                    dictionary = {"Sources": fMRI_despike,
                                  "Description": 'copy.',
                                  "Command": command,}
                    json_object = json.dumps(dictionary, indent=3)

                    with open(fMRI_SliceT.replace('.nii.gz','json'), "w") as outfile:
                        outfile.write(json_object)

        elif isinstance(Slice_timing_info, list) == False and Slice_timing_info.split(' ')[0] == '-tpattern':
            command = (sing_afni+ '3dTshift -wsinc9 ' + overwrite +
                       ' -TR ' + str(TR) + ' ' + Slice_timing_info +
                       ' -prefix ' + fMRI_SliceT +
                       ' ' + fMRI_despike)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": fMRI_despike,
                          "Description": 'Slice timing correction.',
                          "Command": command,}
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_SliceT.replace('.nii.gz','json'), "w") as outfile:
                outfile.write(json_object)

        elif isinstance(Slice_timing_info, list) == True:

            command = (sing_afni + '3dTshift -wsinc9 ' + overwrite +
                       ' -TR ' + str(TR) + ' -tpattern @' + fMRI_stc +
                       ' -prefix ' + fMRI_SliceT +
                       ' ' + fMRI_despike)
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": fMRI_despike,
                          "Description": 'Slice timing correction.',
                          "Command": command,}
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_SliceT.replace('.nii.gz','json'), "w") as outfile:
                outfile.write(json_object)
        else:
            nl = 'ERROR : please check Slice_timing_info, this variable is not define as it should'
            raise ValueError(run_cmd.error(nl, diary_file))

        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   fMRI_runMean +
                   ' ' + fMRI_SliceT)
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": fMRI_SliceT,
                      "Description": 'Mean image.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_runMean.replace('.nii.gz','json'), "w") as outfile:
            outfile.write(json_object)

        # outlier fraction for each volume
        command = (sing_afni +  '3dToutcount' + overwrite + ' -automask -fraction -polort 4 -legendre ' +
                   fMRI_SliceT + ' > ' +
                   fMRI_outcount)
        run_cmd.run(command, diary_file)

        # realignment intra-run (volreg)
        # register each volume to the base image
        command = (sing_afni + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' +
                   fMRI_runMean +
                   ' -1Dfile ' + file_volreg +
                   ' -prefix ' + fMRI_volreg +
                   ' -cubic' +
                   ' -twodup ' +
                   ' -1Dmatrix_save ' + matrix_volreg +
                   ' ' + fMRI_SliceT)
        run_cmd.run(command, diary_file)

        check_nii.keep_header(fMRI_volreg, fMRI_SliceT)

        dictionary = {"Sources": [fMRI_SliceT,
                                  fMRI_runMean],
                      "Description": 'Rigid realignment (3dVolreg from AFNI).',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_volreg.replace('.nii.gz','json'), "w") as outfile:
            outfile.write(json_object)

        # censoring # see ex 10 in 1d_tool
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + file_volreg +
                   ' -derivative -censor_prev_TR -collapse_cols euclidean_norm' +
                   ' -moderate_mask -1.2 1.2 -show_censor_count' +
                   ' -write_censor ' + censore1D +
                   ' -write_CENSORTR ' + censoretxt)
        run_cmd.run(command, diary_file)

        # compute motion magnitude time series: the Euclidean norm
        # (sqrt(sum squares)) of the motion parameter derivatives
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + file_volreg +
                   ' -set_nruns 1 -derivative -collapse_cols euclidean_norm' +
                   ' -write ' + motion_enorm)
        run_cmd.run(command, diary_file)

        # writing regressors # get the first derivative
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + file_volreg +
                   ' -derivative -write ' + deriv)
        run_cmd.run(command, diary_file)

        # writing regressors get demean
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_prepro_orig_process, root + '_dfile.1D') +
                   ' -demean -write ' + opj(dir_prepro_orig_process, root + '_xdtr_demean.1D'))
        run_cmd.run(command, diary_file)

        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   fMRI_runMean_align +
                   ' ' + fMRI_volreg)
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": fMRI_volreg,
                      "Description": 'Mean image.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_runMean_align.replace('.nii.gz','json'), "w") as outfile:
            outfile.write(json_object)

        # BiasFieldCorrection
        IMG = ants.image_read(fMRI_runMean_align)
        N4 = ants.n4_bias_field_correction(IMG,
                                           shrink_factor=4,
                                           convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07},
                                           spline_param=200)
        ants.image_write(N4, fMRI_runMean_n4Bias, ri=False)
        dictionary = {"Sources": opj(dir_prepro_orig_process, root + '_xdtr_mean_preWARP.nii.gz'),
                      "Description": 'Bias field correction (N4).',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_runMean_n4Bias.replace('.nii.gz','json'), "w") as outfile:
            outfile.write(json_object)
