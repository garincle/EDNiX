import os
import ants
import json
import numpy as np
import nibabel as nib
import subprocess
import glob
from pathlib import Path
opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile
spgo = subprocess.getoutput
from Tools import run_cmd,diaryfile, check_nii
from fonctions.extract_filename import extract_filename
from fonctions import _2b_fix_orient

##### XXX add ICA or DL to visualize pre-processing effect

def preprocess_data(dir_prepro_raw_process, RS, list_RS, nb_run, T1_eq, TR, Slice_timing_info, dir_prepro_raw_matrices,
                    overwrite, sing_afni, diary_file,animalPosition, humanPosition, orientation, doWARPonfunc, diary_WARNING):

    nl = '##  Working on step ' + str(1) + '(function: _1_fMRI_preTTT_in_fMRIspace).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    for i in range(0, int(nb_run)):
        nl = 'work on ' + str(dir_prepro_raw_process) + ' run ' + str(i +1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        root = extract_filename(RS[i])
        raw_func = list_RS[i]
        base_fMRI_targeted = opj(dir_prepro_raw_process, root + '_space-func_desc-vol_rmv_from_txt.nii.gz')
        base_fMRI = opj(dir_prepro_raw_process, root + '_space-func_desc-vol_rmv.nii.gz')
        fMRI_despike = opj(dir_prepro_raw_process, root + '_space-func_desc-despiked.nii.gz')
        fMRI_SliceT = opj(dir_prepro_raw_process, root + '_space-func_desc-SliceTfixed.nii.gz')
        fMRI_runMean = opj(dir_prepro_raw_process, root + '_space-func_desc-runMean.nii.gz')
        fMRI_stc = opj(dir_prepro_raw_process, 'stc.txt')
        fMRI_outcount = opj(dir_prepro_raw_process, root + '_space-func_desc-outcount.r$run.1D')

        outpuprefix_motion = opj(dir_prepro_raw_matrices, root + '_space-func_desc-motion_correction')
        mat_files_pattern = opj(dir_prepro_raw_matrices, root + '_space-func_desc-motion_correction_*.mat')
        fMRI_run_motion_corrected = opj(dir_prepro_raw_process, root + '_space-func_desc-motion_corrected.nii.gz')
        fMRI_run_motion_corrected_orient = opj(dir_prepro_raw_process, root + '_space-func_desc-motion_corrected_orient.nii.gz')
        fMRI_run_motion_corrected_orientMEAN = opj(dir_prepro_raw_process, root + '_space-func_desc-motion_corrected_orientMEAN.nii.gz')
        matrix_motion_correction = opj(dir_prepro_raw_matrices, root + '_space-func_desc-motion_correction.1D')

        censore1D = opj(dir_prepro_raw_process, root + '_space-func_desc-censor.1D')
        censoretxt = opj(dir_prepro_raw_process, root + '_space-func_desc-censor.txt')
        demean = opj(dir_prepro_raw_process, root + '_space-func_desc-demean.1D')
        deriv = opj(dir_prepro_raw_process, root + '_space-func_desc-deriv.1D')
        motion_enorm = opj(dir_prepro_raw_process, root + '_space-func_desc-motion_enorm.1D')

        fMRI_runMean_align = opj(dir_prepro_raw_process, root + '_space-func_desc-runMean_align.nii.gz')
        fMRI_runMean_n4Bias = opj(dir_prepro_raw_process, root + '_space-func_desc-runMean_n4Bias.nii.gz')

        fMRI_BASE = opj(dir_prepro_raw_process, root + '_space-func_desc-fMRI_BASE.nii.gz')
        fMRI_BASE_Mean = opj(dir_prepro_raw_process, root + '_space-func_desc-fMRI_BASE_Mean.nii.gz')

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
            dictionary = {"Sources": [base_fMRI_targeted,
                                      opj(opd(raw_func), root, '.txt')],
                          "Description": 'Remove volumes.',
                          "Command": command,}
            json_object = json.dumps(dictionary, indent=3)
            with open(base_fMRI_targeted.replace('.nii.gz','.json'), "w") as outfile:
                outfile.write(json_object)
            raw_func = base_fMRI_targeted
            run_cmd.do(command, diary_file)

        img = nib.load(raw_func)
        nb_vol = img.shape[-1]  # For 4D data, last dimension is time

        command = (sing_afni + '3dTcat -prefix ' + base_fMRI +
                   ' ' + raw_func + '[' + str(T1_eq) + '-' + str(nb_vol - 1) + ']' + overwrite)
        dictionary = {"Sources": base_fMRI,
                      "Description": 'Remove first volumes.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(base_fMRI.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.do(command, diary_file)

        # Despiking
        command = (sing_afni + '3dDespike -NEW -nomask' + overwrite + ' -prefix ' +
                   fMRI_despike +
                   ' ' + base_fMRI)
        dictionary = {"Sources": base_fMRI,
                      "Description": 'Despiking.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_despike.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

        # slice-timing correction -heptic!!!!!!
        if Slice_timing_info == 'Auto':
            if opi(fMRI_stc):
                command = (sing_afni + '3dTshift -wsinc9' + overwrite +
                           ' -TR ' + str(TR) + ' -tpattern @' + fMRI_stc +
                           ' -prefix ' + fMRI_SliceT +
                           ' ' + fMRI_despike)
                dictionary = {"Sources": fMRI_despike,
                              "Description": 'Slice timing correction.',
                              "Command": command,}
                json_object = json.dumps(dictionary, indent=3)
                with open(fMRI_SliceT.replace('.nii.gz','.json'), "w") as outfile:
                    outfile.write(json_object)
                run_cmd.run(command, diary_file)

            else:
                cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";' + sing_afni + '3dinfo -slice_timing ' + \
                      list_RS[0]
                nl = spgo(cmd).split('\n')
                STC = nl[-1].split('|')
                STC = list(map(float, STC))
                if np.sum(STC) > 0:
                    nl = "INFO: SliceTiming = " + str(STC)
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    command = (sing_afni + '3dTshift -wsinc9' + overwrite +
                               ' -prefix ' + fMRI_SliceT +
                               ' ' + fMRI_despike)
                    dictionary = {"Sources": fMRI_despike,
                                  "Description": 'Slice timing correction.',
                                  "Command": command,}
                    json_object = json.dumps(dictionary, indent=3)
                    with open(fMRI_SliceT.replace('.nii.gz','.json'), "w") as outfile:
                        outfile.write(json_object)
                    run_cmd.run(command, diary_file)

                else:
                    nl = "WARNING: Slice Timing not found, this will be particularly DANGEROUS, you SHOULD PROVIDE MANUALLY ONE!"
                    run_cmd.msg(nl, diary_file, 'WARNING')

                    diary_WARNING_file = diaryfile.create(diary_WARNING,nl)

                    command = (sing_afni + '3dcalc -a ' + fMRI_despike +
                               ' -prefix ' + fMRI_SliceT +
                               ' -expr "a"' + overwrite)
                    dictionary = {"Sources": fMRI_despike,
                                  "Description": 'copy.',
                                  "Command": command,}
                    json_object = json.dumps(dictionary, indent=3)
                    with open(fMRI_SliceT.replace('.nii.gz','.json'), "w") as outfile:
                        outfile.write(json_object)
                    run_cmd.do(command, diary_file)

        elif isinstance(Slice_timing_info, list) == False and Slice_timing_info.split(' ')[0] == '-tpattern':
            command = (sing_afni+ '3dTshift -wsinc9' + overwrite +
                       ' -TR ' + str(TR) + ' ' + Slice_timing_info +
                       ' -prefix ' + fMRI_SliceT +
                       ' ' + fMRI_despike)
            dictionary = {"Sources": fMRI_despike,
                          "Description": 'Slice timing correction.',
                          "Command": command,}
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_SliceT.replace('.nii.gz','.json'), "w") as outfile:
                outfile.write(json_object)
            run_cmd.run(command, diary_file)

        elif isinstance(Slice_timing_info, list) == True:
            command = (sing_afni + '3dTshift -wsinc9' + overwrite +
                       ' -TR ' + str(TR) + ' -tpattern @' + fMRI_stc +
                       ' -prefix ' + fMRI_SliceT +
                       ' ' + fMRI_despike)
            dictionary = {"Sources": fMRI_despike,
                          "Description": 'Slice timing correction.',
                          "Command": command,}
            json_object = json.dumps(dictionary, indent=3)
            with open(fMRI_SliceT.replace('.nii.gz','.json'), "w") as outfile:
                outfile.write(json_object)
            run_cmd.run(command, diary_file)
        else:
            nl = 'ERROR : please check Slice_timing_info, this variable is not define as it should'
            raise ValueError(run_cmd.error(nl, diary_file))

        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   fMRI_runMean +
                   ' ' + fMRI_SliceT)

        dictionary = {"Sources": fMRI_SliceT,
                      "Description": 'Mean image.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_runMean.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

        # outlier fraction for each volume
        command = (sing_afni +  '3dToutcount' + overwrite + ' -automask -fraction -polort 4 -legendre ' +
                   fMRI_SliceT + ' > ' +
                   fMRI_outcount)
        run_cmd.run(command, diary_file)

        '''
        # realignment intra-run (volreg)
        # register each volume to the base image
        command = (sing_afni + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' +
                   fMRI_runMean +
                   ' -1Dfile ' + file_motion_correction +
                   ' -prefix ' + fMRI_run_motion_corrected +
                   ' -cubic' +
                   ' -twodup' +
                   ' -1Dmatrix_save ' + matrix_motion_correction +
                   ' ' + fMRI_SliceT)
        dictionary = {"Sources": [fMRI_SliceT,
                                  fMRI_runMean],
                      "Description": 'Rigid realignment (3dVolreg from AFNI).',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_run_motion_corrected.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

        check_nii.keep_header(fMRI_run_motion_corrected, fMRI_SliceT)
        '''

        ### 2.0 Start fix_orient
        _2b_fix_orient.fix_orient(fMRI_BASE, fMRI_SliceT, list_RS,
                                  animalPosition, humanPosition, orientation, doWARPonfunc, sing_afni, diary_file)

        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   fMRI_BASE_Mean +
                   ' ' + fMRI_BASE)
        dictionary = {"Sources": fMRI_BASE,
                      "Description": 'Mean image.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_BASE_Mean.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

        # Realignment intra-run avec ANTs motion_correction
        motion_result = ants.motion_correction(
            image=ants.image_read(fMRI_BASE),
            fixed=ants.image_read(fMRI_BASE_Mean),  # Image de base
            verbose=True,
            type_of_transform='BOLDRigid',
            interpolator='bSpline',
            outprefix=outpuprefix_motion)

        def extract_motion_params_antspy(matrices_dir, output_1D):
            """
            Extrait les paramètres de mouvement (rotations + translations)
            à partir des fichiers *_0GenericAffine.mat produits par ants.motion_correction().
            Sortie compatible 3dDeconvolve (.1D).
            """
            mat_files = sorted(Path(matrices_dir).glob('*_0GenericAffine.mat'))

            if not mat_files:
                raise FileNotFoundError(f"Aucun fichier *_0GenericAffine.mat trouvé dans {matrices_dir}")

            motion_params = []

            for mat in mat_files:
                t = ants.read_transform(str(mat))
                params = t.parameters  # [Rx, Ry, Rz, Tx, Ty, Tz]

                if len(params) < 6:
                    raise ValueError(f"Le fichier {mat} ne contient pas 6 paramètres attendus (trouvé {len(params)}).")

                Rx, Ry, Rz, Tx, Ty, Tz = params[:6]

                # ANTs donne les rotations en radians → conversion degrés pour compatibilité AFNI
                roll, pitch, yaw = np.degrees([Rx, Ry, Rz])

                motion_params.append([roll, pitch, yaw, Tx, Ty, Tz])
            motion_params = np.array(motion_params)
            np.savetxt(output_1D, motion_params, fmt="%.6f", delimiter=" ")

            print(f"[INFO] Paramètres de mouvement sauvegardés dans {output_1D}")
            print(f"[INFO] Format: roll pitch yaw dS dL dP (deg, mm)")

        extract_motion_params_antspy(
            matrices_dir=dir_prepro_raw_matrices,
            output_1D=matrix_motion_correction)

        # Sauvegarder l'image realignée
        motion_result['motion_corrected'].to_filename(fMRI_run_motion_corrected_orient)
        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   fMRI_run_motion_corrected_orientMEAN +
                   ' ' + fMRI_run_motion_corrected_orient)
        dictionary = {"Sources": fMRI_run_motion_corrected_orient,
                      "Description": 'Mean image.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_run_motion_corrected_orientMEAN.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

        # Realignment intra-run avec ANTs motion_correction
        motion_result = ants.motion_correction(
            image=ants.image_read(fMRI_SliceT),
            fixed=ants.image_read(fMRI_runMean),  # Image de base
            verbose=True,
            type_of_transform='BOLDRigid',
            interpolator='bSpline')

        # Sauvegarder l'image realignée
        motion_result['motion_corrected'].to_filename(fMRI_run_motion_corrected)

        """Créer une transformation composite de tous les mouvements"""
        mat_files = sorted(glob.glob(mat_files_pattern))
        print(f"Combining {len(mat_files)} motion transformation files")
        if not mat_files:
            raise ValueError("Aucun fichier .mat trouvé")
        # Commencer avec la transformation identité

        # Créer le JSON de métadonnées
        command = (f"ants.motion_correction with Rigid transform, "
                   f"fixed base image, linear interpolation")

        dictionary = {
            "Sources": [fMRI_SliceT, fMRI_runMean],
            "Description": 'Rigid realignment (ANTs motion_correction).',
            "Command": command,}

        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_run_motion_corrected.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        # censoring # see ex 10 in 1d_tool
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + matrix_motion_correction +
                   ' -derivative -censor_prev_TR -collapse_cols euclidean_norm' +
                   ' -moderate_mask -1.2 1.2 -show_censor_count' +
                   ' -write_censor ' + censore1D +
                   ' -write_CENSORTR ' + censoretxt)
        run_cmd.run(command, diary_file)

        # compute motion magnitude time series: the Euclidean norm
        # (sqrt(sum squares)) of the motion parameter derivatives
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + matrix_motion_correction +
                   ' -set_nruns 1 -derivative -collapse_cols euclidean_norm' +
                   ' -write ' + motion_enorm)
        run_cmd.run(command, diary_file)

        # writing regressors # get the first derivative
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + matrix_motion_correction +
                   ' -derivative -write ' + deriv)
        run_cmd.run(command, diary_file)

        # writing regressors get demean
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + opj(dir_prepro_raw_process, root + '_dfile.1D') +
                   ' -demean -write ' + demean)
        run_cmd.run(command, diary_file)

        command = (sing_afni + '3dTstat' + overwrite + ' -mean -prefix ' +
                   fMRI_runMean_align +
                   ' ' + fMRI_run_motion_corrected)
        dictionary = {"Sources": fMRI_run_motion_corrected,
                      "Description": 'Mean image.',
                      "Command": command,}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_runMean_align.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
        run_cmd.run(command, diary_file)

        # BiasFieldCorrection
        IMG = ants.image_read(fMRI_runMean_align)
        N4 = ants.n4_bias_field_correction(IMG,
                                           shrink_factor=4,
                                           convergence={'iters': [50, 50, 50, 50], 'tol': 1e-07},
                                           spline_param=200)
        ants.image_write(N4, fMRI_runMean_n4Bias, ri=False)
        dictionary = {"Sources": opj(dir_prepro_raw_process, root + '_xdtr_mean_preWARP.nii.gz'),
                      "Description": 'Bias field correction (N4).',}
        json_object = json.dumps(dictionary, indent=3)
        with open(fMRI_runMean_n4Bias.replace('.nii.gz','.json'), "w") as outfile:
            outfile.write(json_object)
