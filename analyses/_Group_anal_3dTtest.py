import nilearn
from nilearn import plotting
import glob
import subprocess
import os
import numpy as np
import nibabel as nib
from nilearn.masking import compute_epi_mask
import matplotlib.pyplot as plt
from scipy.stats import norm
import Tools.Load_EDNiX_requirement
from fonctions.extract_filename import extract_filename
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
def _3dttest_EDNiX(bids_dir, templatehigh, templatelow, oversample_map, mask_func, cut_coords, panda_files, selected_atlases,
              lower_cutoff, upper_cutoff, MAIN_PATH, FS_dir, alpha ,all_ID, all_Session, all_data_path, endfmri, mean_imgs, ntimepoint_treshold):

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
        output_results = opj(output_results1, 'Grp_SBA_3dTTEST')
        if not os.path.exists(opj(output_results1, 'Grp_SBA_3dTTEST')):
            os.mkdir(opj(output_results1, 'Grp_SBA_3dTTEST'))
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

            unique_IDs = np.unique(all_ID)  # Liste des ID uniques

            for ID in unique_IDs:
                # Trouver les index où l'ID correspond
                indices = [i for i, x in enumerate(all_ID) if x == ID]

                subject_results = []

                for idx in indices:
                    Session = all_Session[idx]
                    data_path = all_data_path[idx]

                    dir_fMRI_Refth_RS = opj(data_path, 'func')
                    dir_fMRI_Refth_RS_prepro = opj(dir_fMRI_Refth_RS, '01_prepro')
                    dir_fMRI_Refth_RS_prepro3 = opj(dir_fMRI_Refth_RS_prepro, '03_atlas_space')

                    list_RS = sorted(glob.glob(opj(dir_fMRI_Refth_RS, endfmri)))
                    list_RS_list = list_RS.copy()

                    for imageF in list_RS_list:
                        fmri_image = nib.load(imageF)
                        image_shape = fmri_image.shape

                        if len(image_shape) == 4 and image_shape[3] >= ntimepoint_treshold:
                            root_RS = extract_filename(os.path.basename(imageF))
                            input_results = opj(dir_fMRI_Refth_RS_prepro3, '10_Results', 'SBA', Seed_name)

                            if ope(opj(input_results, root_RS + '_correlations_fish.nii.gz')):
                                subject_results.append(opj(input_results, root_RS + '_correlations_fish.nii.gz'))

                if subject_results:
                    output_file = opj(output_folder, f"{ID}_avg_fisher_map.nii.gz")
                    if ope(output_file):
                        os.remove(output_file)

                    command = f"singularity run {s_bind} {afni_sif} 3dMean -prefix {output_file} {subprocess.list2cmdline(subject_results)}"
                    nl = spgo(command)
                    print(nl)

                # **AVERAGE the Fisher maps across subjects** for each region
                if ope(opj(output_folder, str(ID) + '_avg_fisher_map.nii.gz')):
                    if ope(opj(output_folder, str(ID) + '_avg_fisher_map.nii.gz')):
                        mean_per_subject.append(opj(output_folder, str(ID) + '_avg_fisher_map.nii.gz'))

            if len(mean_per_subject)>1:
                # Find and remove all files starting with 'TTnew'
                for file in glob.glob(os.path.join(output_folder, "TTnew*")):
                    os.remove(file)
                    print(f"Removed: {file}")

                os.chdir(output_folder)
                print(subprocess.list2cmdline(mean_per_subject))
                if len(mean_per_subject)>14:
                    command = f"singularity run {s_bind} {afni_sif} 3dttest++ -setA {subprocess.list2cmdline(mean_per_subject)} " \
                              f"-toz -Clustsim -mask {opj(output_results1, 'mask_mean_func_overlapp.nii.gz')} "
                    nl = spgo(command)
                    print(nl)

                    command = f"singularity run {s_bind} {afni_sif} 3dcalc -overwrite -a {opj(output_folder, 'TTnew+orig.HEAD[1]')} -expr a -prefix {opj(output_folder, Seed_name + 'ttest-stat_fisher_zmap.nii.gz')}"
                    nl = spgo(command)
                    print(nl)

                    # Extract cluster size information
                    cluster_size_command = f'singularity run {s_bind} {afni_sif} 1d_tool.py -infile {output_folder}/TTnew.CSimA.NN1_2sided.1D -csim_show_clustsize -verb 0 -csim_pthr 0.05 -csim_alpha {str(alpha)}'
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

                else:
                    command = f"singularity run {s_bind} {afni_sif} 3dttest++ -setA {subprocess.list2cmdline(mean_per_subject)} " \
                              f"-toz -mask {opj(output_results1, 'mask_mean_func_overlapp.nii.gz')} "
                    nl = spgo(command)
                    print(nl)
                    cluster_size = 10

                    command = f"singularity run {s_bind} {afni_sif} 3dcalc -overwrite -a {opj(output_folder, 'TTnew+orig.HEAD[1]')} -expr a -prefix {opj(output_folder, Seed_name + 'ttest-stat_fisher_zmap.nii.gz')}"
                    nl = spgo(command)
                    print(nl)

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
                                                 colorbar=True, bg_img=studytemplatebrain, display_mode='mosaic', cut_coords=cut_coords)
                display.savefig(opj(output_folder, Seed_name+ 'thresholded_stat_mosaic.jpg'))
                display.close()
                plt.close('all')
