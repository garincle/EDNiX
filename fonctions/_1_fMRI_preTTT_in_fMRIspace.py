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

from Tools import run_cmd,diaryfile
from fonctions.extract_filename import extract_filename

##### XXX add ICA or DL to visualize pre-processing effect

def preprocess_data(dir_fMRI_Refth_RS_prepro1, RS, list_RS, nb_run, T1_eq, TR, Slice_timing_info,
                    overwrite, sing_afni, diary_file,diary_WARNING):

    if ope(dir_fMRI_Refth_RS_prepro1) == False:
        os.makedirs(dir_fMRI_Refth_RS_prepro1)


    nl = '##  Working on step ' + str(1) + '(function: _1_fMRI_preTTT_in_fMRIspace).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    for i in range(0, int(nb_run)):
        nl = 'work on ' + str(dir_fMRI_Refth_RS_prepro1) + ' run ' + str(i +1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        root = extract_filename(RS[i])

        # copy func imag
        command = (sing_afni + '3dcalc -a ' + list_RS[i] + ' -prefix ' +
                   opj(dir_fMRI_Refth_RS_prepro1, root + '.nii.gz') + ' -expr "a"' + overwrite)
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": list_RS[i],
                      "Description": 'Copy.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '.json'), "w") as outfile:
            outfile.write(json_object)

        base_fMRI = opj(dir_fMRI_Refth_RS_prepro1, root + '.nii.gz')

        # Clean bad volumes
        if ope(opj(opd(list_RS[i]), root, '.txt')) == True:
            # Open the file in read mode
            with open(opj(opd(list_RS[i]), root, '.txt'), 'r') as file:
                # Read the first line and convert it to an integer
                cut_low = int(file.readline().strip())

                # Read the second line and convert it to an integer
                cut_high = int(file.readline().strip())

            command = (sing_afni + '3dTcat -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,root + '_x0.nii.gz') +
                       ' ' + base_fMRI + '[' + str(cut_low) + '-' + str(cut_high - 1) + ']' + overwrite)
            run_cmd.do(command, diary_file)

            base_fMRI = opj(dir_fMRI_Refth_RS_prepro1, root + '_x0.nii.gz')
            dictionary = {"Sources": [base_fMRI,
                                      opj(opd(list_RS[i]), root, '.txt')],
                          "Description": 'Remove volumes.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_x0.json'), "w") as outfile:
                outfile.write(json_object)

        img = nib.load(base_fMRI)
        nb_vol = img.shape[-1]  # For 4D data, last dimension is time

        command = (sing_afni + '3dTcat -prefix ' + opj(dir_fMRI_Refth_RS_prepro1,root + '_x.nii.gz') +
                   ' ' + base_fMRI + '[' + str(T1_eq) + '-' + str(nb_vol - 1) + ']' + overwrite)
        run_cmd.do(command, diary_file)

        dictionary = {"Sources": base_fMRI,
                      "Description": 'Remove first volumes.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_x.json'), "w") as outfile:
            outfile.write(json_object)

        # Despiking
        command = (sing_afni + '3dDespike -NEW -nomask' + overwrite + ' -prefix ' +
                   opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz') +
                   ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_x.nii.gz'))
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_x.nii.gz'),
                      "Description": 'Despiking.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.json'), "w") as outfile:
            outfile.write(json_object)

        # slice-timing correction -heptic!!!!!!
        if Slice_timing_info == 'Auto':
            if opi(opj(dir_fMRI_Refth_RS_prepro1, 'stc.txt')):
                command = (sing_afni + '3dTshift -wsinc9 ' + overwrite +
                           ' -TR ' + str(TR) + ' -tpattern @' + opj(dir_fMRI_Refth_RS_prepro1, 'stc.txt') +
                           ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') +
                           ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'))
                run_cmd.run(command, diary_file)

                dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                              "Description": 'Slice timing correction.', }
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
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
                               ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') +
                               ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'))
                    run_cmd.run(command, diary_file)

                    dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                                  "Description": 'Slice timing correction.', }
                    json_object = json.dumps(dictionary, indent=2)
                    with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                        outfile.write(json_object)
                else:
                    nl = "WARNING: Slice Timing not found, this will be particularly DANGEROUS, you SHOULD PROVIDE MANUALLY ONE!"
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    diary_WARNING_file = diaryfile.create(diary_WARNING,nl)

                    command = (sing_afni + '3dcalc -a ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz') +
                               ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz') +
                               ' -expr "a"' + overwrite)
                    run_cmd.do(command, diary_file)

                    dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                                  "Description": 'copy.', }
                    json_object = json.dumps(dictionary, indent=2)

                    with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                        outfile.write(json_object)

        elif isinstance(Slice_timing_info, list) == False and Slice_timing_info.split(' ')[0] == '-tpattern':
            command = (sing_afni+ '3dTshift -wsinc9 ' + overwrite +
                       ' -TR ' + str(TR) + ' ' + Slice_timing_info +
                       ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') +
                       ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'))
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                          "Description": 'Slice timing correction.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                outfile.write(json_object)

        elif isinstance(Slice_timing_info, list) == True:

            command = (sing_afni + '3dTshift -wsinc9 ' + overwrite +
                       ' -TR ' + str(TR) + ' -tpattern @' + opj(dir_fMRI_Refth_RS_prepro1, 'stc.txt') +
                       ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') +
                       ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'))
            run_cmd.run(command, diary_file)

            dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xd.nii.gz'),
                          "Description": 'Slice timing correction.', }
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.json'), "w") as outfile:
                outfile.write(json_object)
        else:
            nl = 'ERROR : please check Slice_timing_info, this variable is not define as it should'
            raise ValueError(run_cmd.error(nl, diary_file))

        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz') +
                   ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'))
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'),
                      "Description": 'Mean image.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.json'), "w") as outfile:
            outfile.write(json_object)

        # outlier fraction for each volume
        command = (sing_afni +  '3dToutcount' + overwrite + ' -automask -fraction -polort 4 -legendre ' +
                   opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz ') + ' > ' +
                   opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_outcount.r$run.1D'))
        run_cmd.run(command, diary_file)

        # realignment intra-run (volreg)
        # register each volume to the base image
        command = (sing_afni + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' +
                   opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz') +
                   ' -1Dfile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') +
                   ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz ') +
                   ' -cubic' +
                   ' -twodup ' +
                   ' -1Dmatrix_save ' + opj(dir_fMRI_Refth_RS_prepro1, root + '.aff12.1D') +
                   ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'))
        run_cmd.run(command, diary_file)

        #with stat values
        extracted_data = nib.load(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz')).get_fdata()
        labeled_img2 = nilearn.image.new_img_like(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'), extracted_data, copy_header=True )
        labeled_img2.to_filename(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz'))

        dictionary = {"Sources": [opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt.nii.gz'),
                                  opj(dir_fMRI_Refth_RS_prepro1, root + '_xdt_mean.nii.gz')],
                      "Description": 'Rigid realignment (3dVolreg from AFNI).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.json'), "w") as outfile:
            outfile.write(json_object)

        # censoring # see ex 10 in 1d_tool
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') +
                   ' -derivative -censor_prev_TR -collapse_cols euclidean_norm' +
                   ' -moderate_mask -1.2 1.2 -show_censor_count' +
                   ' -write_censor ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_censor.1D') +
                   ' -write_CENSORTR ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_CENSORTR.txt'))
        run_cmd.run(command, diary_file)

        # compute motion magnitude time series: the Euclidean norm
        # (sqrt(sum squares)) of the motion parameter derivatives
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') +
                   ' -set_nruns 1 -derivative -collapse_cols euclidean_norm' +
                   ' -write ' + opj(dir_fMRI_Refth_RS_prepro1, root + 'motion_enorm.1D'))
        run_cmd.run(command, diary_file)

        # writing regressors # get the first derivative
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') +
                   ' -derivative -write ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_deriv.1D'))
        run_cmd.run(command, diary_file)

        # writing regressors get demean
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_dfile.1D') +
                   ' -demean -write ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_demean.1D'))
        run_cmd.run(command, diary_file)

        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz') +
                   ' ' + opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz'))
        run_cmd.run(command, diary_file)

        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr.nii.gz'),
                      "Description": 'Mean image.', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.json'), "w") as outfile:
            outfile.write(json_object)

        # BiasFieldCorrection
        IMG = ants.image_read(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz'))
        N4 = ants.n4_bias_field_correction(IMG,
                                           shrink_factor=4,
                                           convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07},
                                           spline_param=200)
        ants.image_write(N4, opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean.nii.gz'), ri=False)
        dictionary = {"Sources": opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean_preWARP.nii.gz'),
                      "Description": 'Bias field correction (N4).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(opj(dir_fMRI_Refth_RS_prepro1, root + '_xdtr_mean.json'), "w") as outfile:
            outfile.write(json_object)
