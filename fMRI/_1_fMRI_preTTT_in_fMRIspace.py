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
from fMRI.extract_filename import extract_filename
from fMRI import _2b_fix_orient

##### XXX add ICA or DL to visualize pre-processing effect

def preprocess_data(dir_prepro_raw_process, RS, list_RS, nb_run, T1_eq, TR, Slice_timing_info, dir_prepro_raw_matrices, n_for_ANTS,
                    overwrite, sing_afni, diary_file,animalPosition, humanPosition, orientation, doWARPonfunc, diary_WARNING):

    nl = '##  Working on step ' + str(1) + '(function: _1_fMRI_preTTT_in_fMRIspace).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    for i in range(0, int(nb_run)):
        nl = 'work on ' + str(dir_prepro_raw_process) + ' run ' + str(i +1)
        run_cmd.msg(nl, diary_file, 'OKGREEN')
        root_RS = extract_filename(RS[i])
        raw_func = list_RS[i]
        base_fMRI_targeted = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-vol_rmv_from_txt.nii.gz')
        base_fMRI = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-vol_rmv.nii.gz')
        fMRI_despike = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-despiked.nii.gz')
        fMRI_SliceT = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-SliceTfixed.nii.gz')
        fMRI_runMean = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean.nii.gz')
        fMRI_stc = opj(dir_prepro_raw_process, 'stc.txt')
        fMRI_outcount = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-outcount_run' + str(i) + '.1D')

        outpuprefix_motion_folder = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-motion_correction')
        outpuprefix_motion = opj(outpuprefix_motion_folder, 'img')
        mat_files_pattern = opj(outpuprefix_motion_folder, 'img*.mat')
        fMRI_run_motion_corrected = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-motion_corrected.nii.gz')
        file_motion_correction = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-motion_correction.1D')
        matrix_motion_params = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-matrix_motion_params.1D')

        censore1D = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-censor.1D')
        censoretxt = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-censor.txt')
        demean = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-demean.1D')
        deriv = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-deriv.1D')
        motion_enorm = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-motion_enorm.1D')

        fMRI_runMean_align = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_align.nii.gz')
        fMRI_runMean_n4Bias = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-runMean_n4Bias.nii.gz')

        fMRI_BASE = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_BASE.nii.gz')
        fMRI_BASE_Mean = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-fMRI_BASE_Mean.nii.gz')

        # Clean bad volumes
        if ope(opj(opd(list_RS[i]), root_RS, '.txt')) == True:
            # Open the file in read mode
            with open(opj(opd(list_RS[i]), root_RS, '.txt'), 'r') as file:
                # Read the first line and convert it to an integer
                cut_low = int(file.readline().strip())

                # Read the second line and convert it to an integer
                cut_high = int(file.readline().strip())

            command = (sing_afni + '3dTcat -prefix ' + base_fMRI_targeted +
                       ' ' + raw_func + '[' + str(cut_low) + '-' + str(cut_high - 1) + ']' + overwrite)
            dictionary = {"Sources": [base_fMRI_targeted,
                                      opj(opd(raw_func), root_RS, '.txt')],
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
                   fMRI_SliceT + ' > ' + fMRI_outcount)
        subprocess.run(command, shell=True, check=True)

        ### 2.0 Start fix_orient
        print(animalPosition)
        _2b_fix_orient.fix_orient(fMRI_BASE, fMRI_SliceT, list_RS,
                                  animalPosition, humanPosition, orientation, doWARPonfunc, sing_afni, diary_file)
        '''
        # realignment intra-run (volreg)
        # register each volume to the base image
        command = (sing_afni + '3dvolreg' + overwrite + ' -verbose -zpad 1 -base ' +
                   fMRI_runMean +
                   ' -1Dfile ' + file_motion_correction +
                   ' -prefix ' + fMRI_run_motion_corrected +
                   ' -cubic' +
                   ' -twodup' +
                   ' -1Dmatrix_save ' + matrix_volreg +
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
        if not ope(outpuprefix_motion_folder):
            os.makedirs(outpuprefix_motion_folder)

        # Realignment intra-run avec ANTs motion_correction
        motion_result = ants.motion_correction(
            image=ants.image_read(fMRI_SliceT),
            fixed=ants.image_read(fMRI_runMean),  # Image de base
            verbose=True,
            type_of_transform='BOLDRigid',
            interpolator=n_for_ANTS,
            outprefix=outpuprefix_motion)

        # Sauvegarder l'image realignée
        motion_result['motion_corrected'].to_filename(fMRI_run_motion_corrected)
        check_nii.keep_header(fMRI_run_motion_corrected, fMRI_SliceT)
        def _ants_affine_to_distance(affine):
            """
            Convert ANTs affine parameters to motion distances and rotations
            Following ANTsPy convention for parameter ordering
            """
            dx, dy, dz = affine[9:]
            rot_x = np.arcsin(affine[6])
            cos_rot_x = np.cos(rot_x)
            rot_y = np.arctan2(affine[7] / cos_rot_x, affine[8] / cos_rot_x)
            rot_z = np.arctan2(affine[3] / cos_rot_x, affine[0] / cos_rot_x)

            deg = np.degrees
            return dx, dy, dz, deg(rot_x), deg(rot_y), deg(rot_z)

        def extract_motion_params_antspy(matrices_dir, output_1D, output_matrix=None):
            """
            Extract motion parameters from ANTs transformation files

            Parameters:
            -----------
            matrices_dir : str
                Directory containing *_0GenericAffine.mat files
            output_1D : str
                Output file for motion parameters (roll pitch yaw Tx Ty Tz)
            output_matrix : str, optional
                Output file for transformation matrices in AFNI format

            Returns:
            --------
            motion_params : np.array
                Motion parameters (N_volumes x 6)
            """
            mat_files = sorted(Path(matrices_dir).glob('*_0GenericAffine.mat'))
            if not mat_files:
                raise FileNotFoundError(f"No *_0GenericAffine.mat files found in {matrices_dir}")

            motion_params = []
            all_matrices = []

            for mat in mat_files:
                # Read the transformation
                transform = ants.read_transform(str(mat))

                # Get affine parameters
                affine = transform.parameters

                # Convert to motion parameters
                dx, dy, dz, rot_x, rot_y, rot_z = _ants_affine_to_distance(affine)

                # Store motion parameters: roll, pitch, yaw, Tx, Ty, Tz
                motion_params.append([dx, dy, dz, rot_x, rot_y, rot_z])

                # Extract 4x4 transformation matrix for matrix output
                if output_matrix is not None:
                    ants_aff = affine.reshape(4, 3).T
                    afni_aff12 = ants_aff.reshape(-1)
                    # AFNI format: flatten row-major and store all 12 values (excluding last row)
                    all_matrices.append(afni_aff12)

            # Save motion parameters
            motion_params = np.array(motion_params)
            np.savetxt(output_1D, motion_params, fmt="%.6f", delimiter=" ")
            print(f"[INFO] {len(motion_params)} volumes processed")
            print(f"[INFO] Motion parameters saved: roll pitch yaw Tx Ty Tz (deg, mm)")

            # Save transformation matrices if requested
            if output_matrix is not None:
                all_matrices = np.array(all_matrices)
                np.savetxt(output_matrix, all_matrices, fmt="%.6f", delimiter=" ")
                print(f"[INFO] Transformation matrices saved to: {output_matrix}")

            return motion_params

        # Extract motion parameters and transformation matrices from ANTs
        extract_motion_params_antspy(
            matrices_dir=outpuprefix_motion_folder,
            output_1D=file_motion_correction,
            output_matrix=matrix_motion_params)

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
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + file_motion_correction +
                   ' -derivative -censor_prev_TR -collapse_cols euclidean_norm' +
                   ' -moderate_mask -1.2 1.2 -show_censor_count' +
                   ' -write_censor ' + censore1D +
                   ' -write_CENSORTR ' + censoretxt)
        run_cmd.run(command, diary_file)

        # compute motion magnitude time series: the Euclidean norm
        # (sqrt(sum squares)) of the motion parameter derivatives
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + file_motion_correction +
                   ' -set_nruns 1 -derivative -collapse_cols euclidean_norm' +
                   ' -write ' + motion_enorm)
        run_cmd.run(command, diary_file)

        # writing regressors # get the first derivative
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + file_motion_correction +
                   ' -derivative -write ' + deriv)
        run_cmd.run(command, diary_file)

        # writing regressors get demean
        command = (sing_afni + '1d_tool.py' + overwrite + ' -infile ' + file_motion_correction +
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
        def adaptive_n4_correction(input_image, output_path, mask=None,
                                   auto_tune=True, species_info=None):
            """
            Adaptive N4 bias field correction for multi-species data

            Parameters:
            -----------
            input_image : str or ANTsImage
                Input image path or ANTs image object
            output_path : str
                Output file path
            mask : str or ANTsImage, optional
                Brain mask to constrain correction
            auto_tune : bool
                Automatically adjust parameters based on brain size
            species_info : dict, optional
                Dictionary with 'species' name for logging
            """
            # Read image
            if isinstance(input_image, str):
                IMG = ants.image_read(input_image)
            else:
                IMG = input_image

            # Get image properties
            voxel_volume = np.prod(IMG.spacing)
            image_shape = IMG.shape

            # Create or read mask
            if mask is None:
                # Adaptive threshold based on image intensity distribution
                img_array = IMG.numpy()
                # Use Otsu-like approach: threshold at valley between background and brain
                hist, bins = np.histogram(img_array[img_array > 0], bins=100)
                # Find a reasonable threshold (e.g., 10% of mean of non-zero voxels)
                threshold = np.mean(img_array[img_array > 0]) * 0.15
                mask_img = ants.get_mask(IMG, low_thresh=threshold)
            elif isinstance(mask, str):
                mask_img = ants.image_read(mask)
            else:
                mask_img = mask

            # Calculate brain volume (in voxels and mm³)
            brain_voxels = np.sum(mask_img.numpy() > 0)
            brain_volume_mm3 = brain_voxels * voxel_volume

            print(f"[INFO] Brain volume: {brain_volume_mm3:.2f} mm³ ({brain_voxels} voxels)")
            print(f"[INFO] Voxel dimensions: {IMG.spacing}")

            # Adaptive parameter selection based on brain size
            if auto_tune:
                # Adjust spline_param based on brain size
                # Smaller brains need finer control (lower spline_param)
                # Larger brains can use coarser control (higher spline_param)

                # Reference: human brain ~1,200,000 mm³, mouse ~500 mm³, rat ~2,000 mm³
                if brain_volume_mm3 < 1000:  # Small animals (mouse, small rat)
                    spline_param = 50
                    shrink_factor = 2
                    convergence_iters = [50, 50, 30, 20]
                    species_category = "small"
                elif brain_volume_mm3 < 5000:  # Medium animals (rat, marmoset)
                    spline_param = 100
                    shrink_factor = 3
                    convergence_iters = [50, 50, 40, 30]
                    species_category = "medium"
                elif brain_volume_mm3 < 100000:  # Large animals (cat, macaque)
                    spline_param = 150
                    shrink_factor = 3
                    convergence_iters = [50, 50, 50, 40]
                    species_category = "large"
                else:  # Very large (human, chimpanzee)
                    spline_param = 200
                    shrink_factor = 4
                    convergence_iters = [50, 50, 50, 50]
                    species_category = "very_large"

                print(f"[INFO] Auto-tuned for {species_category} brain")
                print(f"[INFO] spline_param={spline_param}, shrink_factor={shrink_factor}")
            else:
                # Default conservative parameters
                spline_param = 100
                shrink_factor = 3
                convergence_iters = [50, 50, 40, 30]
                species_category = "default"

            # Apply N4 correction
            N4 = ants.n4_bias_field_correction(
                IMG,
                mask=mask_img,
                shrink_factor=shrink_factor,
                convergence={'iters': convergence_iters, 'tol': 1e-07},
                spline_param=spline_param,
                verbose=True
            )

            # Save corrected image
            ants.image_write(N4, output_path, ri=False)

            # Save bias field for QC
            bias_field_path = output_path.replace('.nii.gz', '_biasfield.nii.gz')
            bias_field = IMG / (N4 + 1e-10)  # Avoid division by zero
            ants.image_write(bias_field, bias_field_path, ri=False)

            # Create comprehensive metadata
            metadata = {
                "Sources": input_image if isinstance(input_image, str) else "ANTsImage",
                "Description": "Bias field correction (N4) with adaptive parameters",
                "BrainVolume_mm3": float(brain_volume_mm3),
                "BrainVolume_voxels": int(brain_voxels),
                "VoxelSpacing": list(IMG.spacing),
                "ImageShape": list(image_shape),
                "SpeciesCategory": species_category,
                "Parameters": {
                    "shrink_factor": shrink_factor,
                    "convergence_iters": convergence_iters,
                    "convergence_tol": 1e-07,
                    "spline_param": spline_param,
                    "mask_used": True,
                    "auto_tuned": auto_tune
                },
                "BiasField": bias_field_path
            }

            if species_info:
                metadata.update(species_info)

            json_object = json.dumps(metadata, indent=3)
            with open(output_path.replace('.nii.gz', '.json'), "w") as outfile:
                outfile.write(json_object)

            return N4, metadata

        # Option 1: Auto-adaptive (recommended for mixed species)
        N4, metadata = adaptive_n4_correction(
            fMRI_runMean_align,
            fMRI_runMean_n4Bias,
            mask=None,  # Will auto-generate
            auto_tune=True,
            species_info=None)

        ### 2.0 Start fix_orient
        _2b_fix_orient.fix_orient(fMRI_BASE_Mean, fMRI_runMean_n4Bias, list_RS,
                                  animalPosition, humanPosition, orientation, doWARPonfunc, sing_afni, diary_file)
