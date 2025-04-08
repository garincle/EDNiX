import os
import subprocess
import datetime
import json
import pandas as pd
import nibabel as nib
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
from fonctions.extract_filename import extract_filename

def signal_regression(dir_fMRI_Refth_RS_prepro1, dir_RS_ICA_native,
    nb_run, RS, blur, TR, ICA_cleaning, extract_exterior_CSF, extract_WM, normalize,
    do_not_correct_signal, band, extract_Vc, extract_GS, overwrite, s_bind,afni_sif,diary_file):

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(7) + '(function: _7_post_TTT).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])

        if ICA_cleaning == 'Skip':
            input = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrf_2ref_RcT_masked.nii.gz')
        else:
            input = opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz')


        command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + input + \
                  ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + ' -expr "a"'
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)
        dictionary = {"Sources": input,
                      "Description": 'Copy.'},
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.json'), "w") as outfile:
            outfile.write(json_object)


    # 4.2 Get the SVD values from the masks
    for i in range(0, int(nb_run)):
        root_RS = extract_filename(RS[i])

        for extract_type, img_name ,suffix in zip([extract_exterior_CSF, extract_WM, extract_GS],
                                                  ['exterior_ligne.nii.gz', 'Wmask.nii.gz', 'maskDilat.nii.gz'],
                                                  ['_NonB','_Wc','-GS']):
            if extract_type == True:
                command = 'singularity run' + s_bind + afni_sif + '3dmaskSVD' + overwrite + ' -polort 2 -vnorm -mask ' + \
                          opj(dir_fMRI_Refth_RS_prepro1, img_name) + \
                          ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                          ' > ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS' + suffix + '.1D')
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

        if extract_Vc == True:
                command = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dROIstats -nomeanout -nzvoxels -mask ' + \
                          opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + ' ' +  opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz')
                count = spgo(command).split('\n')[-1].split('\t')
                # check if ok....
                if count[-1] == '0[?]':
                    countVmask = 0
                else:
                    countVmask = int(count[-1])

                if countVmask>10:
                    command = 'singularity run' + s_bind + afni_sif + '3dmaskSVD' + overwrite + ' -polort 2 -vnorm -mask ' + opj(dir_fMRI_Refth_RS_prepro1,'Vmask.nii.gz') + \
                    ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                    ' > ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS_Vc.1D')
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)


        for option_type, suffix,descript in zip([' -cvarinv',' -cvar', ' -tsnr',' -stdev'],
                                                ['_tsnr1','_cvar','_tsnr2','_stdev'],
                                                ['fabs(mean)/stdev (with detrend)',
                                                 'stdev/fabs(mean) (with detrend)',
                                                 'fabs(mean)/stdev NOT DETRENDED',
                                                 'standard deviation with detrend']):

            command = 'singularity run' + s_bind + afni_sif + '3dTstat' + overwrite + option_type + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS' + suffix + '.nii.gz') + \
                      ' ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz')
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz'),
                          "Description": descript + ' (3dTstat, AFNI)'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS' + suffix + '.json'), "w") as outfile:
                outfile.write(json_object)

        ############################### ############################### ###############################
        ##                                Corrections of the signal                                  ##
        ############################### ############################### ###############################

        # 5.0 Regress out most of the noise from the data: bandpass filter, motion correction white mater noise and cbf noise , plus drift and derivatives
        # after filtering : blur within the mask and normalise the data.

        def get_number_of_trs(nifti_path):
            """Safely get number of TRs using nibabel"""
            try:
                img = nib.load(nifti_path)
                return str(img.shape[-1])  # Returns TRs as string
            except Exception as e:
                raise ValueError(f"nibabel failed to read {nifti_path}: {str(e)}")
        # Usage
        nifti_path = opj(dir_fMRI_Refth_RS_prepro1, f"{root_RS}_xdtrfwS.nii.gz")
        NumberofTR = get_number_of_trs(nifti_path)

        # create bandpass regressors (instead of using 3dBandpass, say)
        command = 'singularity run' + s_bind + afni_sif + '1dBport' + overwrite + ' -nodata ' + NumberofTR + ' ' + str(TR) + ' -band ' + band + \
                  ' -invert -nozero > ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'bandpass_rall.1D')
        nl = spgo(command)
        diary.write(f'\n{nl}')
        print(nl)

        if ope(dir_RS_ICA_native) == False:
            os.makedirs(dir_RS_ICA_native)

    if do_not_correct_signal == False:
        for i in range(0, int(nb_run)):
            try:
                root_RS = extract_filename(RS[i])
                original_dir = os.getcwd()
                os.chdir(dir_fMRI_Refth_RS_prepro1)

                command = 'singularity run' + s_bind + afni_sif + '3dDeconvolve -input ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                ' -mask ' + opj(dir_fMRI_Refth_RS_prepro1,'maskDilat.nii.gz') + \
                ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1,root_RS + 'bandpass_rall.1D') + ' bandpass_rall' + \
                ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_demean.1D') + ' mot_demean' + \
                ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_deriv.1D') + ' mot_deriv' + \
                ' -censor ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_censor.1D') + \
                ' -polort A -float ' +                                                  \
                ' -num_stimts 0 ' + overwrite +                                                      \
                ' -fout -tout '  +                                                          \
                ' -x1D ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.xmat.1D ') +                                                       \
                ' -xjpeg ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.jpg')  +                                                    \
                ' -x1D_uncensored ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.nocensor.xmat.1D') +                                 \
                ' -fitts ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'Xfittssubj') +                                                 \
                ' -errts ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'errts') +                                                         \
                ' -x1D_stop ' +                                                          \
                ' -bucket ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'statssubj')

                for extract_type, suffix in zip([extract_exterior_CSF, extract_WM, extract_GS, extract_Vc],
                                                ['_NonB', '_Wc', '-GS','-Vc']):
                    if extract_type == True:
                        command = command + ' -ortvec ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS' + suffix + '.1D') + ' residual_norm' + suffix + ' '

                if ope(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_confounds_correct.tsv')):
                    confounds_df = pd.read_csv(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_confounds_correct.tsv'), sep='\t')

                    for column in confounds_df.columns:
                        column_file = opj(dir_fMRI_Refth_RS_prepro1, root_RS + f"_{column}.1D")
                        confounds_df[[column]].to_csv(column_file, sep=' ', index=False, header=False)

                        # Append the column file to the command
                        command += f' -ortvec {column_file} residual_norm_{column} '


                nl = 'INFO: 3dDeconvolve command is ' + command
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                command = 'singularity run' + s_bind + afni_sif + '3dTproject -polort 0' + overwrite + ' -input ' + \
                opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                ' -censor ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_censor.1D') + \
                ' -cenmode ZERO -ort ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'X.nocensor.xmat.1D') + \
                ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')
                if blur>0:
                    command = command + ' -blur ' + str(blur)
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)

                # Normalization options
                zscore_command = 'singularity run ' + s_bind + ' ' + afni_sif + ' ' + \
                                 '3dTstat -mean -overwrite -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'mean.nii.gz') + ' ' + \
                                 opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + ' && ' + \
                                 '3dTstat -stdev -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'stdev.nii.gz') + ' ' + \
                                 opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + ' && ' + \
                                 '3dcalc -overwrite -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + \
                                 ' -b ' + opj(dir_fMRI_Refth_RS_prepro1, 'mean.nii.gz') + \
                                 ' -c ' + opj(dir_fMRI_Refth_RS_prepro1, 'stdev.nii.gz') + \
                                 ' -expr "(a-b)/c" -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')

                psc_command = 'singularity run ' + s_bind + ' ' + afni_sif + ' ' + \
                              '3dTstat -mean -overwrite -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, 'mean_func.nii.gz') + ' ' + \
                              opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + ' && ' + \
                              '3dcalc -overwrite -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + \
                              ' -b ' + opj(dir_fMRI_Refth_RS_prepro1, 'mean_func.nii.gz') + \
                              ' -expr "((a - b) / b) * 100" -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')

                # Choose normalization method
                if normalize == 'zscore':
                    nl = spgo(zscore_command)
                    diary.write(f'\n{nl}')
                    print(nl)
                elif normalize == 'psc':
                    nl = spgo(psc_command)
                    diary.write(f'\n{nl}')
                    print(nl)

                else:
                    print('no normalization')

                dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz'),
                                          opj(dir_fMRI_Refth_RS_prepro1, 'maskDilat.nii.gz'),
                                          opj(dir_fMRI_Refth_RS_prepro1, root_RS + 'bandpass_rall.1D'),
                                          opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_demean.1D'),
                                          opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_deriv.1D'),
                                          opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtr_censor.1D')],
                              "Description": 'Filtering and blurring: ' + str(blur) + ' (3dDeconvolve and 3dTproject, AFNI)'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.json'), "w") as outfile:
                    outfile.write(json_object)

            except:
                root_RS = extract_filename(RS[i])
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + \
                          ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                          ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_3dDeconvolve_failed.nii.gz') + ' -expr "a"'
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz'),
                              "Description": 'Copy'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_3dDeconvolve_failed.json'), "w") as outfile:
                    outfile.write(json_object)
                os.chdir(original_dir)

            if not ope(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')):
                root_RS = extract_filename(RS[i])
                command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + \
                          ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                          ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_3dDeconvolve_failed.nii.gz') + ' -expr "a"'
                nl = spgo(command)
                diary.write(f'\n{nl}')
                print(nl)
                dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz'),
                              "Description": 'Copy'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_3dDeconvolve_failed.json'), "w") as outfile:
                    outfile.write(json_object)
                os.chdir(original_dir)
    else:
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + overwrite + ' -a ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz') + \
                      ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz') + ' -expr "a"'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_xdtrfwS.nii.gz'),
                          "Description": 'Copy'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.json'), "w") as outfile:
                outfile.write(json_object)

    diary.write(f'\n')
    diary.close()
