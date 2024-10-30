import nilearn
import os
import subprocess
import re
import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname
ope = os.path.exists
spgo = subprocess.getoutput
#################################################################################################
####Seed base analysis
#################################################################################################
def anat_QC(type_norm, labels_dir, dir_prepro, ID, listTimage, masks_dir, s_bind, afni_sif):

    for Timage in listTimage:
        direction = opj(labels_dir, type_norm)
        atlas_filename = opj(direction, 'atlaslvl1_LR.nii.gz')
        lines = []
        anat_filename = opj(dir_prepro, ID + '_acpc_test_QC' + Timage + '.nii.gz')
        brain_mask = opj(masks_dir,'brain_mask_in_anat_DC.nii.gz')
        output_results =  opj(dir_prepro, 'QC_anat')
        if not os.path.exists(output_results): os.mkdir(output_results)
        
        if ope(anat_filename):
            if ope(atlas_filename) == False:
                print(bcolors.WARNING + 'WARNING: no altlas lvl 1 LR found, this is a requirement for some of QC analysis')
            else:
                # Load the NIfTI images
                img1 = nib.load(atlas_filename)
                img2 = nib.load(anat_filename)

                # Extract headers
                header1 = img1.header
                header2 = img2.header

                # Compare specific fields
                fields_to_compare = ['dim', 'pixdim']
                tolerance = 0.0000006
                differences = {}
                for field in fields_to_compare:
                    value1 = header1[field][1:4]
                    value2 = header2[field][1:4]
                    if not np.allclose(value1, value2, atol=tolerance):
                        differences[field] = (value1, value2)

                # Print differences
                if differences:
                    print(bcolors.OKGREEN + "Differences found in the following fields:" + bcolors.ENDC)
                    for field, values in differences.items():
                        print(bcolors.OKGREEN + f"{field}:" + bcolors.ENDC)
                        print(bcolors.OKGREEN + f"  Image 1: {values[0]}" + bcolors.ENDC)
                        print(bcolors.OKGREEN + f"  Image 2: {values[1]}" + bcolors.ENDC)
                        caca = nilearn.image.resample_to_img(atlas_filename, anat_filename, interpolation='nearest')
                        caca.to_filename(atlas_filename)
                        extracted_data = nib.load(atlas_filename).get_fdata()
                        labeled_img2 = nilearn.image.new_img_like(anat_filename, extracted_data, copy_header=True)
                        labeled_img2.to_filename(atlas_filename)
                else:
                    print(bcolors.OKGREEN + "No differences found in the specified fields." + bcolors.ENDC)

                # Atlas labels as described by you
                labels = {2: '3rd ventricles', 3: '4th ventricles', 7: 'Cerebellum', 5: 'Cerebellum White',
                          6: 'Cortical White matter', 1: 'CSF', 4: 'Lateral ventricles', 8: 'Isocortex',
                          9: 'Allocortex', 10: 'Periallocortex', 11: 'Subcortical areas', 12: 'Diencephalon',
                          13: 'Brain stem'}

                # Load the atlas and fMRI data
                atlas_img = nib.load(atlas_filename)
                atlas_data = atlas_img.get_fdata()
                anat_img = nib.load(anat_filename)
                anat_data = anat_img.get_fdata()

                # Function to calculate SNR for each region and hemisphere
                def compute_snr_single_frame(atlas_data, anat_data, labels, hemisphere_offset=1000):
                    """
                    Computes SNR for each region in a single-frame anatomical image.

                    :param atlas_data: 3D array with labeled regions
                    :param anat_data: 3D anatomical image data array
                    :param labels: Dictionary of region labels and names
                    :param hemisphere_offset: Offset for right hemisphere label (not needed for single-frame gray matter)
                    :return: SNR values per region and the calculated noise value
                    """
                    snr_values = {}

                    for label, region_name in labels.items():
                        # Left hemisphere
                        left_mask = atlas_data == label
                        left_signal = anat_data[left_mask]

                        # Right hemisphere
                        right_mask = atlas_data == label + hemisphere_offset
                        right_signal = anat_data[right_mask]

                        # Compute signal for each hemisphere
                        left_signal_avg = np.mean(left_signal) if left_signal.size > 0 else np.nan
                        right_signal_avg = np.mean(right_signal) if right_signal.size > 0 else np.nan

                        # Compute noise as the standard deviation of the lowest 10% values across the brain
                        brain_values = anat_data[anat_data > 0].flatten()  # Exclude zero values
                        noise = np.std(brain_values[brain_values <= np.percentile(brain_values, 10)]) if brain_values.size > 0 else np.nan

                        # SNR = signal / noise
                        snr_values[region_name] = {
                            'Left': left_signal_avg / noise if noise != 0 else np.nan,
                            'Right': right_signal_avg / noise if noise != 0 else np.nan
                        }

                    return snr_values, noise

                # Example usage
                atlas_data = nib.load(atlas_filename).get_fdata()
                anat_data = nib.load(anat_filename).get_fdata()
                snr_results, noise = compute_snr_single_frame(atlas_data, anat_data, labels)

                # Prepare data for CSV and plotting
                data = {'Region': list(snr_results.keys())}
                data['Left_SNR'] = [values['Left'] for values in snr_results.values()]
                data['Right_SNR'] = [values['Right'] for values in snr_results.values()]

                # Save SNR values to CSV
                df = pd.DataFrame(data)
                df.to_csv(opj(output_results, f'{Timage}_snr_regions_single_frame.csv'), index=False)

                # Plot SNR values per region
                plt.figure(figsize=(10, 6))
                plt.bar(df['Region'], df['Left_SNR'], label='Left Hemisphere', alpha=0.6)
                plt.bar(df['Region'], df['Right_SNR'], label='Right Hemisphere', alpha=0.6, bottom=df['Left_SNR'])
                plt.xlabel('Region')
                plt.ylabel('SNR')
                plt.title('SNR of Brain Regions in Single Frame')
                plt.xticks(rotation=45, ha='right')
                plt.legend()
                plt.tight_layout()
                plt.savefig(opj(output_results, f'{Timage}_snr_single_frame_plot.png'))
                plt.close()

                # Function to compute average SNR across all regions for gray and white matter (no left and right separation)
                def compute_average_snr_gray_white(snr_results, gray_regions, white_regions):
                    avg_gray_snr = []
                    avg_white_snr = []
                    time_points = len(next(iter(snr_results.values()))['Left'])  # Length of time points

                    for t in range(time_points):
                        # Average gray matter SNR
                        gray_snr = []
                        for region in gray_regions:
                            gray_snr.append(np.mean([snr_results[region]['Left'][t], snr_results[region]['Right'][t]]))
                        avg_gray_snr.append(np.mean(gray_snr))

                        # Average white matter SNR
                        white_snr = []
                        for region in white_regions:
                            white_snr.append(np.mean([snr_results[region]['Left'][t], snr_results[region]['Right'][t]]))
                        avg_white_snr.append(np.mean(white_snr))

                    return avg_gray_snr, avg_white_snr

                # Define which regions belong to gray and white matter
                gray_regions = ['Isocortex', 'Allocortex', 'Periallocortex', 'Subcortical areas', 'Diencephalon']
                white_regions = ['Cortical White matter']

                # Compute the average SNR for gray and white matter
                avg_gray_snr, avg_white_snr = compute_average_snr_gray_white(snr_results, gray_regions, white_regions)

                # Add gray and white matter SNR values to the DataFrame
                df['Gray_Matter'] = avg_gray_snr
                df['White_Matter'] = avg_white_snr

                # Extract the average SNR across all time points for each region, including gray and white
                average_snr_results = {}
                # Compute the average SNR for each region
                for region, hemispheres in snr_results.items():
                    avg_left_snr = np.mean(hemispheres['Left'])
                    avg_right_snr = np.mean(hemispheres['Right'])
                    average_snr_results[region] = {'Left': avg_left_snr, 'Right': avg_right_snr}

                # Add the average gray and white matter to the dictionary
                average_snr_results['Gray_Matter'] = np.mean(avg_gray_snr)
                average_snr_results['White_Matter'] = np.mean(avg_white_snr)

                # Print the average SNR for each region
                print("Average SNR for each region:")
                for region, values in average_snr_results.items():
                    if isinstance(values, dict):
                        print(f"{region}: Left = {values['Left']:.2f}, Right = {values['Right']:.2f}")
                    else:
                        print(f"{region}: Average SNR = {values:.2f}")

                # Build the `line_snr` variable, including the average SNR for gray and white matter, and for each brain region.
                line_snr = [
                    f"  avg snr Gray_Matter: {str(average_snr_results['Gray_Matter'])}",
                    f"  avg snr White_Matter: {str(average_snr_results['White_Matter'])}"]

                # Add average SNR for each specific region (left and right hemispheres)
                for region, values in average_snr_results.items():
                    if region not in ['Gray_Matter', 'White_Matter']:  # Skip already handled gray/white matter
                        if isinstance(values, dict):  # If region has left and right hemisphere data
                            line_snr.append(f"  avg snr {region} Left: {str(values['Left'])}")
                            line_snr.append(f"  avg snr {region} Right: {str(values['Right'])}")
                        else:  # If region is a single value (like gray/white matter)
                            line_snr.append(f"  avg snr {region}: {str(values)}")

                cnr_val = np.abs(average_snr_results['Gray_Matter'] - average_snr_results['White_Matter']) / noise
                cortical_contrast = (average_snr_results['White_Matter'] - average_snr_results['Gray_Matter']) / ((average_snr_results['White_Matter'] + average_snr_results['Gray_Matter']) / 2)

                line_QC_func_atlaslvl1 = [f"  cortical_contrast: {cortical_contrast}",
                                f"  cnr: {cnr_val}",
                                f"  average_snr_Gray_Matter: {str(average_snr_results['Gray_Matter'])}",
                                f"  average_snr_White_Matter: {str(average_snr_results['White_Matter'])}",]

                lines.append(line_QC_func_atlaslvl1)
                lines.append(line_snr)

            #################### QC that doesn't require atlaslvl1 ####################

            # FWHM Calculation
            def fwhm(anat_file, mask_file):
                """Calculate the FWHM of the input image using AFNI's 3dFWHMx.
                - Uses AFNI 3dFWHMx. More details here:
                    https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dFWHMx.html

                :type anat_file: str
                :param anat_file: The filepath to the anatomical image NIFTI file.
                :type mask_file: str
                :param mask_file: The filepath to the binary head mask NIFTI file.
                :type out_vox: bool
                :param out_vox: (default: False) Output the FWHM as number of voxels
                                instead of mm (the default).
                :rtype: tuple
                :return: A tuple of the FWHM values (x, y, z, and combined).
                """

                command = f'singularity run {s_bind} {afni_sif} 3dFWHMx -overwrite -combined ' \
                          f'-mask {mask_file} ' \
                          f'-input {anat_file}'
                fwhm_string_list = spgo([command])
                print(fwhm_string_list)

                try:
                    # Use regex to find the line with the four numeric values
                    match = re.search(r'(\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+)', fwhm_string_list)
                    if match:
                        retcode = str(match.group(1))
                        print(retcode)
                    vals = np.array(retcode.split(), dtype=np.float64)

                except Exception as e:
                    err = "\n\n[!] Something went wrong with AFNI's 3dFWHMx. Error " \
                          "details: %s\n\n" % e
                    raise Exception(err)

                return list(vals)[3]

            fwhm_val = fwhm(anat_filename,
                            atlas_filename)
            print(fwhm_val)

            line_QC_func = [f"  fwhm_val: {fwhm_val}"]

            def compute_snr_brain_mask(gray_mask_path, anat_data_path, output_path, root_name):
                """
                Computes SNR within a specified gray matter mask, plots the SNR over time, and saves the average SNR.

                :param gray_mask_path: Path to binary gray matter mask NIfTI file
                :param anat_data_path: Path to 4D fMRI data NIfTI file
                :param output_path: Directory path to save results
                :param root_name: Root name for output files
                """
                # Load the data
                gray_mask_img = nib.load(gray_mask_path)
                fmri_img = nib.load(anat_data_path)
                gray_mask = gray_mask_img.get_fdata().astype(bool)
                anat_data = fmri_img.get_fdata()

                time_points = anat_data.shape[-1]
                snr_values = []

                for t in range(time_points):
                    # Signal in gray matter mask at time point t
                    mask_signal = anat_data[gray_mask, t]

                    # Compute noise from bottom 10% of histogram values in the whole brain, excluding zeros
                    brain_values = anat_data[..., t].flatten()
                    non_zero_values = brain_values[brain_values > 0]

                    if len(non_zero_values) > 0:
                        noise_threshold = np.percentile(non_zero_values, 10)
                        noise_values = non_zero_values[non_zero_values <= noise_threshold]
                        noise = np.std(noise_values)
                    else:
                        raise ValueError("10% of noise is just 0; noise calculation cannot be completed.")

                    # Calculate SNR for the time point
                    snr = np.mean(mask_signal) / noise
                    snr_values.append(snr)

                # Calculate the average SNR
                avg_snr = np.mean(snr_values)
                print(f"Average SNR in gray matter mask: {avg_snr}")
                return avg_snr

            # Example usage
            compute_snr_brain_mask_val = compute_snr_brain_mask(
                brain_mask,
                anat_filename,
                output_results, Timage)
            line_QC_func.append(f"  snr_brain_mask {compute_snr_brain_mask_val}")

            lines.append(line_QC_func)
            flattened_list = [item for sublist in lines for item in sublist]

            print(bcolors.OKGREEN + 'QC will look like ' + str(flattened_list) + bcolors.ENDC)

            with open(output_results + '/' + Timage + 'QC_result.txt', 'w') as f:
                for line in flattened_list:
                    f.write(line)
                    f.write('\n')
        else:
            print(bcolors.WARNING + 'WARNING: ' + str(anat_filename) + ' not found!!' + bcolors.ENDC)