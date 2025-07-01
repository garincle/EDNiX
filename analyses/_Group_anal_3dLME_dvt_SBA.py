import nilearn
from nilearn import plotting
import glob
import subprocess
import os
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm
import Tools.Load_EDNiX_requirement
import Launcher_analysis.Load_BIDS_data_for_analysis
#################################################################################################
#### LOADER YUNG LEMUR
#################################################################################################
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

#################################################################################################
#### Seed base analysis
#################################################################################################
def _3dLME_dev_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, folder_atlases, xcell_extrernal_data, panda_files, selected_atlases,
              lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha ,all_ID, all_Session, all_data_path, max_sessionlist, endfmri, mean_imgs,
                ntimepoint_treshold):

    s_path, afni_sif, fsl_sif, fs_sif, itk_sif, wb_sif, strip_sif, s_bind = Tools.Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, bids_dir, FS_dir)
    output_results1 = opj(bids_dir, 'Results')
    if not os.path.exists(output_results1): os.mkdir(output_results1)

    # If oversampling is enabled, use the base template, else use a predefined atlas
    if oversample_map == True:
        studytemplatebrain = templatehigh
    else:
        studytemplatebrain = templatelow

    # Concatenate all the images in `mean_imgs`
    mean_imgs_rs = nilearn.image.concat_imgs(mean_imgs, ensure_ndim=None, memory=None, memory_level=0,
                                             auto_resample=True, verbose=0)
    mask_img = compute_epi_mask(mean_imgs_rs,
                                lower_cutoff=lower_cutoff, upper_cutoff=upper_cutoff,
                                connected=True, opening=1,
                                exclude_zeros=True, ensure_finite=True)

    mask_img.to_filename(opj(output_results1, 'mask_mean_func.nii.gz'))

    # Apply mask using AFNI 3dmask_tool
    command = f"singularity run {s_bind} {afni_sif} 3dmask_tool -overwrite -prefix {opj(output_results1, 'mask_mean_func.nii.gz')} " \
              f"-input {opj(output_results1, 'mask_mean_func.nii.gz')} -fill_holes"
    nl = spgo(command)
    print(nl)

    # Resample to match the mask function
    command = f"singularity run {s_bind} {afni_sif} 3dresample -master {opj(output_results1, 'mask_mean_func.nii.gz')} -prefix {opj(output_results1, 'mask_mean_func_orig.nii.gz')} " \
              f"-input {mask_func} -overwrite -bound_type SLAB"
    nl = spgo(command)
    print(nl)

    # Apply the mask to the original data
    command = f"singularity run {s_bind} {afni_sif} 3dcalc -a {opj(output_results1, 'mask_mean_func.nii.gz')} -b {opj(output_results1, 'mask_mean_func_orig.nii.gz')} " \
              f"-expr 'a*b' -prefix {opj(output_results1, 'mask_mean_func_overlapp.nii.gz')} -overwrite"
    nl = spgo(command)
    print(nl)

    # Iterate through each atlas and corresponding region
    for panda_file, atlas in zip(panda_files, selected_atlases):
        output_results = opj(output_results1, 'Grp_SBA_3dLME_network_dvt')
        if not os.path.exists(opj(output_results1)):
            os.mkdir(opj(output_results1))
        if not os.path.exists(output_results):
            os.mkdir(output_results)

        # Define function to format seed names
        def format_seed_name(seed_name):
            replace_chars = [' ', '(', ')', ',', '/', ':', ';', '.', '-']
            for char in replace_chars:
                seed_name = seed_name.replace(char, '_')
            formatted_name = ''.join(char for char in seed_name if char.isalnum() or char == '_')
            formatted_name = formatted_name.strip('_')
            return formatted_name

        def run_command(command):
            """Helper function to run shell commands."""
            print(f"Running command: {command}")
            return subprocess.getoutput(command)

        # Loop through each region in the pandas dataframe
        for column, row in panda_file.T.items():
            Seed_name = row['region']
            Seed_name = format_seed_name(Seed_name)

            output_folder = opj(output_results, Seed_name)
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)

            # List to store individual subject Fisher maps for averaging later
            mean_per_subject = []
            # Loop over all subjects and collect their Fisher maps for each region
            all_images = []
            all_Run_updated = []
            all_Session_updated = []
            all_ID_updated = []
            for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
                dir_fMRI_Refth_RS = opj(data_path, 'func')
                dir_fMRI_Refth_RS_prepro = opj(dir_fMRI_Refth_RS, '01_prepro')
                dir_fMRI_Refth_RS_prepro3 = opj(dir_fMRI_Refth_RS_prepro, '03_atlas_space')

                # Check the func runs
                list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))
                RS = [os.path.basename(i) for i in list_RS]

                if len(list_RS) == 0:
                    nl = 'ERROR : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endfmri) and here it is ' + str(
                        opj(dir_fMRI_Refth_RS, endfmri)) + ' I would check how you define "endfmri"'
                    raise ValueError(nl)
                list_RS_list = list_RS.copy()
                list_pop_index = []

                for imageF in list_RS_list:
                    # Load the fMRI NIfTI image
                    fmri_image = nib.load(imageF)
                    # Get the shape of the image (x, y, z, t)
                    image_shape = fmri_image.shape
                    # Check the number of time points (4th dimension)
                    if len(image_shape) == 4:
                        ntimepoint = image_shape[3]  # The 4th dimension represents time

                        if int(ntimepoint) < ntimepoint_treshold:
                            index_of_imageF = list_RS.index(imageF)
                            list_RS.pop(index_of_imageF)
                            list_pop_index.append(index_of_imageF)

                nb_run = len(list_RS)
                # Setup for distortion correction using Fieldmaps
                for i in range(0, int(nb_run)):
                    root_RS = extract_filename(RS[i])
                    input_results = opj(dir_fMRI_Refth_RS_prepro3, '10_Results', 'SBA', Seed_name)
                    if ope(opj(input_results, root_RS + '_correlations_fish.nii.gz')):
                        all_images.append(opj(input_results, root_RS + '_correlations_fish.nii.gz'))
                        all_Run_updated.append(str(i))
                        all_Session_updated.append(Session)
                        all_ID_updated.append(ID)

                        if len(all_images) == 1:
                            Resample_master = opj(input_results, root_RS + '_correlations_fish.nii.gz')
                        if len(all_images) > 1:
                            # Check if image matches reference grid
                            check_cmd = f"singularity run {s_bind} {afni_sif} 3dMatch -quiet -source {Resample_master} -input {opj(input_results, root_RS + '_correlations_fish.nii.gz')}"
                            result = subprocess.run(check_cmd, shell=True, capture_output=True, text=True)

                            if result.returncode != 0:  # Grid mismatch
                                all_images = all_images[:-1]
                                resampled_path = opj(input_results, root_RS + '_correlations_fish.nii.gz').replace('.nii.gz', '_resampled.nii.gz')
                                resample_cmd = f"singularity run {s_bind} {afni_sif} 3dresample -overwrite -master {Resample_master} -input {opj(input_results, root_RS + '_correlations_fish.nii.gz')} -prefix {resampled_path}"
                                subprocess.run(resample_cmd, shell=True, check=True)
                                all_images.append(opj(input_results, root_RS + '_correlations_fish_resampled.nii.gz'))

            # Design matrix for GLM
            # Create a DataFrame with unique subject-run pairs
            panda_design_matrix = pd.DataFrame({'Subj': all_ID_updated, 'Sess': all_Session_updated, 'run': all_Run_updated, 'InputFile': all_images})

            ### add extra variable
            panda_design_matrix = panda_design_matrix.merge(
                xcell_extrernal_data[['Subj', 'Sess', 'maturity', 'age']],
                on=['Subj', 'Sess'],
                how='left'  # Use 'inner' if you only want matching rows
            )
            panda_design_matrix = panda_design_matrix[['Subj', 'Sess', 'run', 'maturity', 'age', 'InputFile']]

            # Ensure no duplicates and sort for consistency
            panda_design_matrix = panda_design_matrix.drop_duplicates().sort_values(by=['Subj', 'run'])

            # Save design matrix
            base_filename = 'design_matrix.txt'
            design_matrix_txt = opj(output_folder, base_filename)

            if os.path.exists(design_matrix_txt):
                os.remove(design_matrix_txt)
            panda_design_matrix.to_csv(design_matrix_txt, index=False, sep='\t')


            stat_maps = opj(output_folder, '3dLME_glt.nii.gz')
            if os.path.exists(stat_maps):
                os.remove(stat_maps)
            if os.path.exists(opj(output_folder, 'resid.nii.gz')):
                os.remove(opj(output_folder, 'resid.nii.gz'))

            # Set model for fixed and random effects
            # Set model for fixed and random effects
            os.chdir(output_results)
            command = f"singularity run {s_bind} {afni_sif} 3dLME" + \
                      ' -prefix ' + output_folder + '/3dLME_glt.nii.gz' + \
                      ' -jobs' + ' 20' + ' -mask ' + opj(output_results1, 'mask_mean_func_overlapp.nii.gz') + \
                      ' -model' + ' "maturity*age"' + \
                      ' -qVars' + ' "age"' + ' -ranEff' + ' "~1|Subj"' + ' -num_glt' + ' 5' + \
                      ' -gltLabel 1 "ageNM" -gltCode 1 "age :" -gltLabel 2 "1MNM" -gltCode 2 "maturity : ' + \
                      '1*M -1*NM age :" -gltLabel 3 "1M" -gltCode 3 "maturity : 1*M age :" -gltLabel 4 ' + \
                      '"1NM" -gltCode 4 "maturity : 1*NM age :" -gltLabel 5 "groupeffect" -gltCode 5 "maturity : 0.5*M +0.5*NM"' + ' -dataTable @' + \
                      output_folder + '/design_matrix.txt' + ' -resid ' + output_folder + '/resid.nii.gz'
            spco(command, shell=True)

            # Further analysis with 3dFWHMx, ClustSim, etc. (remaining code...)
            #command = f"singularity run {s_bind} {afni_sif} 3dFWHMx -detrend {subprocess.list2cmdline(all_images)} " \
            #          f" -input {output_results + 'resid.nii.gz'}  -unif"
            #spco(command, shell=True)

            command = f"singularity run {s_bind} {afni_sif} 3dClustSim -mask " \
            f" {opj(output_results1, 'mask_mean_func_overlapp.nii.gz')} -LOTS -prefix {output_folder + '/Clust_'}"
            spco(command, shell=True)

            stat_maps = output_folder + '/3dLME_glt.nii.gz'

            for i, gltlabel in zip([5, 7, 9, 11, 13], ['age', '1MNM', '1M', '1NM', 'groupeffect']):

                img_glt = output_folder + '/' + Seed_name + '_' + str(gltlabel) + '.nii.gz'
                output_z = nib.load(stat_maps).get_fdata()[:, :, :, 0, i]
                output_z = nilearn.image.new_img_like(stat_maps, output_z, copy_header=True)
                output_z.to_filename(img_glt)

                cluster_file = output_folder + '/Clust_.NN1_2sided.1D'

                def extract_cluster_threshold(filename, pthr, alpha):
                    # Read the file
                    with open(filename, 'r') as f:
                        lines = f.readlines()

                    # Extract alpha values from the header
                    header_line = lines[6]  # The alpha values are located in the 4th line
                    alpha_values = [float(a) for a in header_line.split("|")[1].split()]

                    # Check if the requested alpha exists in the file
                    if alpha not in alpha_values:
                        print(f"Invalid alpha={alpha}. Available alphas: {alpha_values}")
                        return

                    # Find the column index for the specified alpha
                    col_idx = alpha_values.index(alpha) + 1  # Adjust for 0-based indexing

                    # Load the numeric data (skip comments)
                    data = np.loadtxt(filename, comments="#")

                    # Find the row index for the specified pthr
                    row_idx = np.where(np.isclose(data[:, 0], pthr, atol=1e-6))[0]
                    if len(row_idx) == 0:
                        print(f"No matching pthr={pthr} found in {filename}.")
                        return

                    # Extract the cluster size from the data
                    cluster_size = data[row_idx[0], col_idx]

                    # Print the result
                    print(f"Cluster size threshold for pthr={pthr}, alpha={alpha}: {cluster_size:.1f} voxels")
                    return cluster_size

                # Example usage:
                cluster_size = extract_cluster_threshold(cluster_file, alpha, 0.05)
                # Apply thresholding
                loadimg = nib.load(img_glt).get_fdata()
                loadimgsort99 = np.percentile(np.abs(loadimg)[np.abs(loadimg) > 0], 99)

                z_score = norm.ppf(1 - alpha / 2)
                print(f"The z-score corresponding to an alpha level of {alpha} is {z_score:.2f}")

                # Threshold the image
                mask_imag = nilearn.image.threshold_img(img_glt, z_score, cluster_threshold=float(cluster_size))
                mask_imag.to_filename(opj(output_folder, Seed_name + '_' + str(gltlabel) + '_thresholded.nii.gz'))

                # Visualization
                display = plotting.plot_stat_map(opj(output_folder, Seed_name + '_' + str(gltlabel) + '_thresholded.nii.gz'), dim=0, threshold=z_score, vmax=loadimgsort99,
                                                 colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=10)
                display.savefig(opj(output_folder, Seed_name + '_' + str(gltlabel) + '_thresholded_stat_mosaic.jpg'))
                display.close()
                plt.close('all')