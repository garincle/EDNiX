import nilearn
import os
import subprocess
import re
import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
import datetime

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
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput
#################################################################################################
####Seed base analysis
#################################################################################################
def anat_QC(type_norm, labels_dir, dir_prepro, ID, listTimage, masks_dir, s_bind, afni_sif,diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anatomical._16_anat_QC_SNR.anat_QC'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    for Timage in listTimage:
        direction          = opj(labels_dir)
        atlas_filename_rsp = opj(direction, type_norm + 'atlaslvl1_LRrspQC' + Timage + '.nii.gz')
        atlas_filename     = opj(direction, type_norm + 'atlaslvl1_LR.nii.gz')
        lines              = []
        anat_filename      = opj(dir_prepro, ID + '_acpc_test_QC_' + Timage + '.nii.gz')
        brain_mask         = opj(masks_dir, Timage + 'brain_mask_final_QCrsp.nii.gz')
        output_results     =  opj(dir_prepro, 'QC_anat')

        if not ope(output_results): os.mkdir(output_results)
        
        if ope(anat_filename):
            if ope(atlas_filename) == False:
                nl = 'WARNING: no atlas lvl 1 LR found, this is a requirement for some of QC analysis'
                print(bcolors.WARNING +  + bcolors.ENDC)
                diary.write(f'\n{nl}')
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
                    nl = "Differences found in the following fields:"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    for field, values in differences.items():
                        nl = field + ':'
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        nl = "  Image 1: " + str(values[0])
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        nl = f"  Image 2: {values[1]}"
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        caca = nilearn.image.resample_to_img(atlas_filename, anat_filename, interpolation='nearest')
                        caca.to_filename(atlas_filename)
                        extracted_data = nib.load(atlas_filename).get_fdata()
                        extracted_data = np.rint(extracted_data).astype(np.int32)
                        labeled_img2 = nilearn.image.new_img_like(anat_filename, extracted_data, copy_header=True)
                        labeled_img2.to_filename(atlas_filename_rsp)
                        atlas_filename = opj(direction, type_norm + 'atlaslvl1_LRrspQC' + Timage + '.nii.gz')
                else:
                    nl = "No differences found in the specified fields."
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                # Atlas labels as described by you
                labels = {2: '3rd ventricles', 3: '4th ventricles', 7: 'Cerebellum', 5: 'Cerebellum White',
                          6: 'Cortical White matter', 1: 'CSF', 4: 'Lateral ventricles', 8: 'Isocortex',
                          9: 'Allocortex', 10: 'Periallocortex', 11: 'Subcortical areas', 12: 'Diencephalon',
                          13: 'Brain stem'}

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
                        if noise == 0:
                            print(bcolors.WARNING + 'std of noise values = 0, try with 20% of the lowest values in the img' + bcolors.ENDC)
                            noise = np.std(brain_values[brain_values <= np.percentile(brain_values,20)]) if brain_values.size > 0 else np.nan
                            if noise == 0:
                                print(bcolors.WARNING + 'std of noise values = 0, try with 30% of the lowest values in the img' + bcolors.ENDC)
                                noise = np.std(brain_values[brain_values <= np.percentile(brain_values,30)]) if brain_values.size > 0 else np.nan
                                if noise == 0:
                                    print(bcolors.WARNING + 'std of noise values = 0, try with 40% of the lowest values in the img' + bcolors.ENDC)
                                    noise = np.std(brain_values[brain_values <= np.percentile(brain_values,40)]) if brain_values.size > 0 else np.nan
                                    if noise == 0:
                                        print(bcolors.WARNING + 'std of noise values = 0, try with 50% of the lowest values in the img' + bcolors.ENDC)
                                        noise = np.std(brain_values[brain_values <= np.percentile(brain_values, 50)]) if brain_values.size > 0 else np.nan
                                    else:
                                        raise ValueError(
                                            bcolors.FAIL + "more than 50% background is the same value; noise calculation cannot be completed." + bcolors.ENDC)

                        # SNR = signal / noise
                        snr_values[region_name] = {
                            'Left': left_signal_avg / noise if noise != 0 else np.nan,
                            'Right': right_signal_avg / noise if noise != 0 else np.nan
                        }
                    return snr_values, noise

                # Example usage
                atlas_data = nib.load(atlas_filename).get_fdata()
                atlas_data = np.rint(atlas_data).astype(np.int32)
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
                    """
                    Computes the average SNR for gray and white matter from single-frame SNR results.

                    :param snr_results: Dictionary with SNR values for each region
                    :param gray_regions: List of region names for gray matter
                    :param white_regions: List of region names for white matter
                    :return: Average SNR values for gray and white matter
                    """
                    gray_snr = []
                    white_snr = []

                    # Calculate average SNR for gray matter regions
                    for region in gray_regions:
                        if region in snr_results:
                            region_snr = np.mean([snr_results[region]['Left'], snr_results[region]['Right']])
                            gray_snr.append(region_snr)

                    # Calculate average SNR for white matter regions
                    for region in white_regions:
                        if region in snr_results:
                            region_snr = np.mean([snr_results[region]['Left'], snr_results[region]['Right']])
                            white_snr.append(region_snr)

                    # Return the average SNR across all gray and white matter regions
                    avg_gray_snr = np.mean(gray_snr) if gray_snr else np.nan
                    avg_white_snr = np.mean(white_snr) if white_snr else np.nan

                    return avg_gray_snr, avg_white_snr

                # Define which regions belong to gray and white matter
                gray_regions = ['Isocortex', 'Allocortex', 'Periallocortex', 'Subcortical areas', 'Diencephalon']
                white_regions = ['Cortical White matter']

                # Compute the average SNR for gray and white matter
                avg_gray_snr, avg_white_snr = compute_average_snr_gray_white(snr_results, gray_regions, white_regions)

                # Create a dictionary to store the average SNR for each region
                average_snr_results = {}
                for region, hemispheres in snr_results.items():
                    avg_left_snr = hemispheres['Left']
                    avg_right_snr = hemispheres['Right']
                    average_snr_results[region] = {'Left': avg_left_snr, 'Right': avg_right_snr}

                # Add the average gray and white matter SNR to the dictionary
                average_snr_results['Gray_Matter'] = avg_gray_snr
                average_snr_results['White_Matter'] = avg_white_snr

                # Print the average SNR for each region
                nl  = "Average SNR for each region:"
                print(nl)
                diary.write(f'\n{nl}')

                for region, values in average_snr_results.items():
                    if isinstance(values, dict):
                        nl = f"{region}: Left = {values['Left']:.2f}, Right = {values['Right']:.2f}"
                        print(nl)
                        diary.write(f'\n{nl}')
                    else:
                        nl=f"{region}: Average SNR = {values:.2f}"
                        print(nl)
                        diary.write(f'\n{nl}')

                # Build the `line_snr` variable, including the average SNR for gray and white matter, and each brain region
                line_snr = [
                    f"  avg snr Gray_Matter: {average_snr_results['Gray_Matter']:.2f}",
                    f"  avg snr White_Matter: {average_snr_results['White_Matter']:.2f}"
                ]

                # Add average SNR for each specific region (left and right hemispheres)
                for region, values in average_snr_results.items():
                    if region not in ['Gray_Matter', 'White_Matter']:
                        line_snr.append(f"  avg snr {region} Left: {values['Left']:.2f}")
                        line_snr.append(f"  avg snr {region} Right: {values['Right']:.2f}")


                cnr_val = np.abs(average_snr_results['Gray_Matter'] - average_snr_results['White_Matter']) / noise
                cortical_contrast = (average_snr_results['White_Matter'] - average_snr_results['Gray_Matter']) / ((average_snr_results['White_Matter'] + average_snr_results['Gray_Matter']) / 2)

                line_QC_func_atlaslvl1 = [f"  cortical_contrast: {cortical_contrast}",
                                f"  cnr: {cnr_val}",
                                f"  average_snr_Gray_Matter: {str(average_snr_results['Gray_Matter'])}",
                                f"  average_snr_White_Matter: {str(average_snr_results['White_Matter'])}",]

                lines.append(line_QC_func_atlaslvl1)
                lines.append(line_snr)

            #################### QC that doesn't require atlaslvl1 ####################
            line_QC_func = []
            try:
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

                    original_dir = os.getcwd()
                    os.chdir(dir_prepro)

                    command = f'singularity run {s_bind} {afni_sif} 3dFWHMx -overwrite ' \
                              f'-mask {mask_file} ' \
                              f'-input {anat_file}'
                    fwhm_string_list = spgo([command])
                    print(fwhm_string_list)
                    os.chdir(original_dir)

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
                line_QC_func.append(f"  fwhm_val: {fwhm_val}")
            except:
                print(bcolors.WARNING + 'FWHM calculation failed (the fix is complicated, but it is not a big deal)' + bcolors.ENDC)

            def compute_snr_brain_mask(gray_mask_path, anat_data_path, output_path, root_name):
                """
                Computes SNR within a specified gray matter mask for single-frame data, plots the SNR histogram,
                and saves the average SNR to a text file.

                :param gray_mask_path: Path to binary gray matter mask NIfTI file
                :param anat_data_path: Path to single-frame anatomical data NIfTI file
                :param output_path: Directory path to save results
                :param root_name: Root name for output files
                """
                # Load the data
                gray_mask_img = nib.load(gray_mask_path)
                anat_img = nib.load(anat_data_path)
                gray_mask = gray_mask_img.get_fdata()
                gray_mask = np.rint(gray_mask).astype(np.int32)
                anat_data = anat_img.get_fdata()

                # Check if anat_data has a fourth dimension (multi-frame data)
                if anat_data.ndim == 4:
                    # Select the first frame along the fourth axis
                    anat_data = anat_data[..., 0]
                else:
                    # Use anat_data as-is if it’s already 3D
                    anat_data = anat_data

                # Signal within the gray matter mask
                mask_signal = anat_data[gray_mask > 0]

                # Compute noise from the bottom 10% of histogram values in the whole brain, excluding zeros
                brain_values = anat_data.flatten()
                non_zero_values = brain_values[brain_values > 0]

                if len(non_zero_values) > 0:
                    noise_threshold = np.percentile(non_zero_values, 10)
                    noise_values = non_zero_values[non_zero_values <= noise_threshold]
                    noise = np.std(noise_values)
                    if noise == 0:
                        print(bcolors.WARNING + 'std of noise values = 0, try with 20% of the lowest values in the img' + bcolors.ENDC)
                        noise_threshold = np.percentile(non_zero_values, 20)
                        noise_values = non_zero_values[non_zero_values <= noise_threshold]
                        noise = np.std(noise_values)
                        if noise == 0:
                            print(bcolors.WARNING + 'std of noise values = 0, try with 30% of the lowest values in the img' + bcolors.ENDC)
                            noise_threshold = np.percentile(non_zero_values, 30)
                            noise_values = non_zero_values[non_zero_values <= noise_threshold]
                            noise = np.std(noise_values)
                            if noise == 0:
                                print(bcolors.WARNING + 'std of noise values = 0, try with 40% of the lowest values in the img' + bcolors.ENDC)
                                noise_threshold = np.percentile(non_zero_values, 40)
                                noise_values = non_zero_values[non_zero_values <= noise_threshold]
                                noise = np.std(noise_values)
                                if noise == 0:
                                    print(bcolors.WARNING + 'std of noise values = 0, try with 50% of the lowest values in the img' + bcolors.ENDC)
                                    noise_threshold = np.percentile(non_zero_values, 50)
                                    noise_values = non_zero_values[non_zero_values <= noise_threshold]
                                    noise = np.std(noise_values)
                                else:
                                    raise ValueError(bcolors.FAIL + "more than 50% background is the same value; noise calculation cannot be completed." + bcolors.ENDC)
                else:
                    raise ValueError(bcolors.FAIL + "10% of noise is just 0; noise calculation cannot be completed." + bcolors.ENDC)

                # Calculate SNR for the single frame
                snr = np.mean(mask_signal) / noise

                # Plot histogram of SNR values within the gray matter mask
                plt.figure(figsize=(8, 5))
                plt.hist(mask_signal, bins=30, color='skyblue', edgecolor='black')
                plt.axvline(snr, color='red', linestyle='--', label=f'Average SNR: {snr:.2f}')
                plt.xlabel('Signal Intensity')
                plt.ylabel('Frequency')
                plt.title('SNR Histogram within Gray Matter Mask')
                plt.legend()
                plt.tight_layout()
                plt.savefig(f"{output_path}/{root_name}_snr_histogram.png")
                plt.close()

                # Print and return the average SNR
                print(f"Average SNR in gray matter mask: {snr}")
                return snr

            # Example usage
            compute_snr_brain_mask_val = compute_snr_brain_mask(
                brain_mask,
                anat_filename,
                output_results, Timage)
            line_QC_func.append(f"  snr_brain_mask {compute_snr_brain_mask_val}")

            lines.append(line_QC_func)
            flattened_list = [item for sublist in lines for item in sublist]

            nl = 'QC will look like ' + str(flattened_list)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            with open(opj(output_results, Timage + 'QC_result.txt'), 'w') as f:
                for line in flattened_list:
                    f.write(line)
                    f.write('\n')
        else:
            nl = 'WARNING: ' + str(anat_filename) + ' not found!!'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()