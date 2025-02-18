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
def _3dMEMA_EDNiX(bids_dir, BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY, cut_coordsZ, panda_files, selected_atlases,
              lower_cutoff, upper_cutoff, s_bind, afni_sif, alpha ,all_ID, all_Session, all_data_path, max_sessionlist, endfmri, mean_imgs, ntimepoint_treshold):

    from fonctions.extract_filename import extract_filename

    output_results1 = opj(bids_dir, 'Results')
    if not os.path.exists(output_results1): os.mkdir(output_results1)

    # If oversampling is enabled, use the base template, else use a predefined atlas
    if oversample_map == True:
        studytemplatebrain = BASE_SS
    else:
        studytemplatebrain = opj(folder_atlases, 'BASE_SS_fMRI.nii.gz')

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
        output_results = opj(output_results1, 'Grp_SBA_3dMEMA')
        if not os.path.exists(opj(output_results1, 'Grp_SBA_3dMEMA')):
            os.mkdir(opj(output_results1, 'Grp_SBA_3dMEMA'))
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

            # List to store individual runs (not averaged) for each subject
            all_runs_data = []
            subject_run_labels = []

            # Loop over all subjects and collect their Fisher maps for each region
            for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
                dir_fMRI_Refth_RS = opj(data_path, 'func')
                dir_fMRI_Refth_RS_prepro = opj(dir_fMRI_Refth_RS, '01_prepro')
                dir_fMRI_Refth_RS_prepro3 = opj(dir_fMRI_Refth_RS_prepro, '03_atlas_space')

                # Check the func runs
                list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))
                RS = [os.path.basename(i) for i in list_RS]

                if len(list_RS) == 0:
                    raise ValueError(f'ERROR: No func image found at {opj(dir_fMRI_Refth_RS, endfmri)}')

                list_pop_index = []
                for imageF in list_RS.copy():
                    fmri_image = nib.load(imageF)
                    image_shape = fmri_image.shape

                    if len(image_shape) == 4 and image_shape[3] < ntimepoint_treshold:
                        list_RS.remove(imageF)
                        list_pop_index.append(list_RS.index(imageF))

                nb_run = len(list_RS)
                for i in range(nb_run):
                    root_RS = extract_filename(RS[i])
                    input_results = opj(dir_fMRI_Refth_RS_prepro3, '10_Results', 'SBA', Seed_name)

                    if ope(opj(input_results, root_RS + '_correlations_fish.nii.gz')):
                        all_runs_data.append(opj(input_results, root_RS + '_correlations_fish.nii.gz'))
                        subject_run_labels.append(f"sub-{ID}_run-{i + 1}")

            if len(all_runs_data) > 1:
                os.chdir(output_folder)

                # Prepare and run `3dMEMA`
                mema_output = opj(output_folder, Seed_name + '_mema-stat_fisher.nii.gz')

                # Building `3dMEMA` command properly
                mema_input_list = []
                for subject, run_data in zip(subject_run_labels, all_runs_data):
                    mema_input_list.append(f"{subject} {run_data}")

                mema_input_str = " \\\n    ".join(mema_input_list)

                mema_command = f"""
                singularity run {s_bind} {afni_sif} 3dMEMA \\
                    -prefix {mema_output} \\
                    -jobs 8 \\
                    -set {Seed_name} \\
                    {mema_input_str} \\
                    -max_zeros 4 \\
                    -model_outliers \\
                    -residual_Z
                """

                nl = spgo(mema_command)
                print(nl)

                command = f"singularity run {s_bind} {afni_sif} 3dcalc -overwrite -a {opj(output_folder, 'TTnew+orig.HEAD[1]')} -expr a -prefix {opj(output_folder, Seed_name + 'ttest-stat_fisher_zmap.nii.gz')}"
                nl = spgo(command)
                print(nl)

                # Extract cluster size information
                cluster_size_command = f'singularity run {s_bind} {afni_sif} 1d_tool.py -infile {output_folder}/TTnew.CSimA.NN1_2sided.1D -csim_show_clustsize -verb 0 -csim_pthr {str(alpha)}'
                cluster_size_output = run_command(cluster_size_command)
                # Extract the cluster size from the output
                try:
                    # Split the output into lines and find the line containing the cluster size
                    for line in cluster_size_output.splitlines():
                        if line.strip().isdigit():  # Check if the line is a number
                            cluster_size = int(line.strip())
                            print(f"Extracted Cluster Size: {cluster_size}")
                            break
                    else:
                        raise ValueError("Cluster size not found in the output.")
                except ValueError as e:
                    print(f"Error extracting cluster size: {e}")
                    cluster_size = 0  # Set a default value or handle the error as needed

                print(str(cluster_size))

                # Thresholding and visualization
                z_map = opj(output_folder, Seed_name + 'ttest-stat_fisher_zmap.nii.gz')
                loadimg = nib.load(z_map).get_fdata()
                loadimgsort99 = np.percentile(np.abs(loadimg)[np.abs(loadimg) > 0], 99)

                # Calculate the z-score for a two-tailed test
                z_score = norm.ppf(1 - alpha / 2)
                print(f"The z-score corresponding to an alpha level of {alpha} is {z_score:.2f}")

                # Threshold the image
                mask_imag = nilearn.image.threshold_img(z_map, z_score, cluster_threshold=int(cluster_size))  # Cluster threshold could be adjusted
                mask_imag.to_filename(opj(output_folder, Seed_name + 'thresholded_ttest-stat_fisher.nii.gz'))

                # Visualization
                display = plotting.plot_stat_map(opj(output_folder, Seed_name+ 'thresholded_ttest-stat_fisher.nii.gz'), dim=0, threshold=z_score, vmax=loadimgsort99,
                                                 colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=10)
                display.savefig(opj(output_folder, Seed_name+ 'thresholded_stat_mosaic.jpg'))
                display.close()
                plt.close('all')
