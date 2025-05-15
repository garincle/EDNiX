import nilearn
import os
import subprocess
import re
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import ants
import datetime
import json
from matplotlib.gridspec import GridSpec
from nilearn import plotting
from sklearn.preprocessing import StandardScaler

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
def save_qc_values(output_results, root_RS, qc_values):
    """Save QC values to JSON file with proper serialization."""

    def numpy_to_python(obj):
        if isinstance(obj, (np.floating, np.integer)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (list, tuple)):
            return [numpy_to_python(x) for x in obj]
        elif isinstance(obj, dict):
            return {k: numpy_to_python(v) for k, v in obj.items()}
        return obj

    try:
        # Convert all numpy types to native Python types
        py_qc_values = numpy_to_python(qc_values)

        # Save to JSON file
        with open(opj(output_results, root_RS + '_QC_values.json'), 'w') as f:
            json.dump(py_qc_values, f, indent=2)

    except Exception as e:
        print(f"ERROR saving QC values: {str(e)}")
        # Try saving with simplified values
        try:
            simple_values = {}
            for k, v in qc_values.items():
                if isinstance(v, (np.floating, np.integer)):
                    simple_values[k] = float(v)
                elif isinstance(v, (list, tuple, np.ndarray)):
                    simple_values[k] = [float(x) if isinstance(x, (np.floating, np.integer)) else str(x)
                                        for x in v]
                else:
                    simple_values[k] = str(v)

            with open(opj(output_results, root_RS + '_QC_values_simple.json'), 'w') as f:
                json.dump(simple_values, f, indent=2)
        except Exception as e2:
            print(f"ERROR saving simplified QC values: {str(e2)}")


def create_qc_figure(output_results, Timage, snr_results, Ref_file,
                     template_in_anat, average_snr_results, noise,
                     anat_filename, BASE_SS_coregistr, fwhm_val, nmi,
                     compute_snr_brain_mask_val, cnr_val, cortical_contrast):
    """
    Create a comprehensive QC figure with multiple subplots including:
    - SNR plots
    - Image correlation
    - Coverage overlay
    - Comprehensive metrics table
    """
    try:
        # Create figure with constrained layout
        fig = plt.figure(figsize=(20, 18), constrained_layout=True)
        plt.style.use('seaborn' if 'seaborn' in plt.style.available else 'ggplot')

        # Set font sizes
        plt.rcParams.update({
            'font.size': 10,
            'axes.titlesize': 12,
            'axes.labelsize': 10,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 9
        })

        # Create grid layout (3 rows, 2 columns)
        gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 1.2])

        # Subplot A: Regional SNR
        ax1 = fig.add_subplot(gs[0, 0])
        regions = list(snr_results.keys())
        left_snr = [snr_results[region]['Left'] for region in regions]
        right_snr = [snr_results[region]['Right'] for region in regions]

        x = np.arange(len(regions))
        width = 0.35
        ax1.bar(x - width / 2, left_snr, width, label='Left', alpha=0.7)
        ax1.bar(x + width / 2, right_snr, width, label='Right', alpha=0.7)

        ax1.set_xticks(x)
        ax1.set_xticklabels(regions, rotation=45, ha='right')
        ax1.set_ylabel('SNR')
        ax1.set_title('A: Regional SNR')
        ax1.legend()
        ax1.grid(True, linestyle=':', alpha=0.5)

        # Subplot B: Normalized Intensity Distribution
        ax3 = fig.add_subplot(gs[0, 1])
        try:
            anat_img = nib.load(Ref_file)
            anat_data = anat_img.get_fdata()
            template_img = nib.load(template_in_anat)
            template_data = template_img.get_fdata()

            # Resample template to match anat dimensions if needed
            if anat_data.shape != template_data.shape:
                template_img_resampled = nilearn.image.resample_to_img(
                    template_img, anat_img, interpolation='continuous')
                template_data = template_img_resampled.get_fdata()

            # Normalize intensities to [0,1] range
            mask = (template_data != 0) & (anat_data != 0) & np.isfinite(template_data) & np.isfinite(anat_data)
            anat_vals = anat_data[mask]
            template_vals = template_data[mask]
            anat_vals = StandardScaler().fit_transform(anat_vals.reshape(-1, 1)).flatten()
            template_vals = StandardScaler().fit_transform(template_vals.reshape(-1, 1)).flatten()
            # Subsample for plotting
            anat_vals_sub = anat_vals[::10]
            template_vals_sub = template_vals[::10]

            ax3.hist(anat_vals_sub, bins=50, alpha=0.5, label='Anatomical', density=True, color='blue')
            ax3.hist(template_vals_sub, bins=50, alpha=0.5, label='Template', density=True, color='orange')
            ax3.set_xlabel('Normalized Intensity')
            ax3.set_ylabel('Density')
            ax3.set_title('B: Normalized Intensity Distribution')
            ax3.legend()
            ax3.grid(True, linestyle=':', alpha=0.5)
        except Exception as e:
            print(f"Error creating intensity distribution plot: {str(e)}")
            ax3.axis('off')
            ax3.set_title('B: Intensity Distribution (Failed)')

        # Subplot C: Image Correlation
        ax4 = fig.add_subplot(gs[1, 0])
        try:
            if 'anat_vals' in locals() and 'template_vals' in locals():
                # Calculate correlation
                cc = np.corrcoef(anat_vals, template_vals)[0, 1]

                ax4.scatter(template_vals_sub, anat_vals_sub, alpha=0.3, s=5, color='green')
                ax4.set_xlabel('Template Intensity')
                ax4.set_ylabel('Anatomical Intensity')
                ax4.set_title(f'C: Intensity Correlation (r={cc:.2f})')
                ax4.grid(True, linestyle=':', alpha=0.5)
        except Exception as e:
            print(f"Error creating correlation plot: {str(e)}")
            ax4.axis('off')
            ax4.set_title('C: Intensity Correlation (Failed)')

        # Subplot D: Coverage Overlay with proper cuts
        ax5 = fig.add_subplot(gs[1, 1])
        try:
            # Find center of mass for better cuts
            anat_img = nib.load(anat_filename)
            anat_data = anat_img.get_fdata()

            display = plotting.plot_anat(anat_filename, axes=ax5,
                                         title='D: Anat-Template registration',
                                         draw_cross=False)
            display.add_contours(BASE_SS_coregistr, colors='r', linewidths=0.5)
        except Exception as e:
            print(f"Error creating coverage plot: {str(e)}")
            ax5.axis('off')
            ax5.set_title('D: Coverage (Failed)')

        # Subplot E: Comprehensive Metrics Table
        ax6 = fig.add_subplot(gs[2, :])  # Span both columns
        ax6.axis('off')

        # Safe formatting function
        def safe_format(value, format_spec=".2f", default="N/A", unit=""):
            if value is None:
                return f"{default}{unit}"
            try:
                if isinstance(value, (np.ndarray, list)):
                    if len(value) == 0:
                        return f"{default}{unit}"
                    value = np.nanmean(value)
                float_val = float(value)
                return f"{float_val:{format_spec}}{unit}" if not np.isnan(float_val) else f"{default}{unit}"
            except (TypeError, ValueError):
                return f"{default}{unit}"

        # Create simplified metrics table without CNR
        metrics_data = [
            ["Metric", "Value", "Description"],
            ["FWHM", safe_format(fwhm_val, '.2f', unit=' mm'), "Image smoothness"],
            ["NMI", safe_format(nmi, '.3f'), "Registration quality"],
            ["Template Correlation", safe_format(cc if 'cc' in locals() else None, '.2f'), "Intensity correlation"],
            ["Global SNR", safe_format(compute_snr_brain_mask_val, '.1f'), "Whole brain SNR"],
            ["Gray Matter SNR", safe_format(average_snr_results['Gray_Matter'], '.1f'), "Average in GM"],
            ["White Matter SNR", safe_format(average_snr_results['White_Matter'], '.1f'), "Average in WM"],
            ["Noise Estimate", safe_format(noise, '.2f'), "Background noise level"]
        ]

        # Create table with adjusted column widths
        table = ax6.table(
            cellText=metrics_data,
            loc='center',
            cellLoc='left',
            colWidths=[0.25, 0.15, 0.6]  # Adjusted column widths
        )

        # Style table
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)  # Adjusted scaling

        # Highlight important metrics
        for (i, j), cell in table.get_celld().items():
            if i == 0:  # Header row
                cell.set_text_props(weight='bold', fontsize=11)
            if j == 0 and i > 0 and metrics_data[i][0] in ['FWHM', 'NMI', 'Global SNR']:
                cell.set_text_props(weight='bold')

        # Add title to the table
        ax6.set_title('E: Key QC Metrics', y=1.05, fontsize=12)

        # Save figure
        plt.savefig(
            opj(output_results, f'{Timage}_QC_summary.png'),
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()

        # Save QC values including regional SNR
        qc_values = {
            'fwhm': fwhm_val,
            'nmi': nmi,
            'template_correlation': cc if 'cc' in locals() else None,
            'global_snr': compute_snr_brain_mask_val,
            'avg_snr_gray': average_snr_results['Gray_Matter'],
            'avg_snr_white': average_snr_results['White_Matter'],
            'cortical_contrast': cortical_contrast,
            'regional_snr': snr_results,
            'noise_estimate': noise
        }
        save_qc_values(output_results, Timage, qc_values)

        # Save comprehensive QC results to text file
        with open(opj(output_results, Timage + '_QC_result.txt'), 'w') as f:
            f.write("=== COMPREHENSIVE QC REPORT ===\n\n")
            f.write("=== IMAGE QUALITY METRICS ===\n")
            f.write(f"FWHM (image smoothness): {safe_format(fwhm_val, '.2f', unit=' mm')}\n")
            f.write(f"Global SNR (brain mask): {safe_format(compute_snr_brain_mask_val, '.1f')}\n")
            f.write(f"Noise estimate: {safe_format(noise, '.2f')}\n\n")

            f.write("=== REGISTRATION QUALITY ===\n")
            f.write(f"NMI (template registration): {safe_format(nmi, '.3f')}\n")
            f.write(f"Template correlation: {safe_format(cc if 'cc' in locals() else None, '.2f')}\n\n")

            f.write("=== TISSUE CONTRAST METRICS ===\n")
            f.write(f"Gray Matter SNR: {safe_format(average_snr_results['Gray_Matter'], '.1f')}\n")
            f.write(f"White Matter SNR: {safe_format(average_snr_results['White_Matter'], '.1f')}\n\n")

            f.write("=== REGIONAL SNR VALUES ===\n")
            for region, values in snr_results.items():
                f.write(f"{region}:\n")
                f.write(f"  Left = {safe_format(values['Left'], '.1f')}\n")
                f.write(f"  Right = {safe_format(values['Right'], '.1f')}\n")

            f.write("\n=== END OF REPORT ===\n")

    except Exception as e:
        print(f"ERROR in create_qc_figure: {str(e)}")
        raise
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
        noise = np.std(brain_values[brain_values <= np.percentile(brain_values,
                                                                  10)]) if brain_values.size > 0 else np.nan
        if noise == 0:
            print(
                bcolors.WARNING + 'std of noise values = 0, try with 20% of the lowest values in the img' + bcolors.ENDC)
            noise = np.std(brain_values[brain_values <= np.percentile(brain_values,
                                                                      20)]) if brain_values.size > 0 else np.nan
            if noise == 0:
                print(
                    bcolors.WARNING + 'std of noise values = 0, try with 30% of the lowest values in the img' + bcolors.ENDC)
                noise = np.std(brain_values[brain_values <= np.percentile(brain_values,
                                                                          30)]) if brain_values.size > 0 else np.nan
                if noise == 0:
                    print(
                        bcolors.WARNING + 'std of noise values = 0, try with 40% of the lowest values in the img' + bcolors.ENDC)
                    noise = np.std(brain_values[brain_values <= np.percentile(brain_values,
                                                                              40)]) if brain_values.size > 0 else np.nan
                    if noise == 0:
                        print(
                            bcolors.WARNING + 'std of noise values = 0, try with 50% of the lowest values in the img' + bcolors.ENDC)
                        noise = np.std(brain_values[brain_values <= np.percentile(brain_values,
                                                                                  50)]) if brain_values.size > 0 else np.nan
                    else:
                        raise ValueError(
                            bcolors.FAIL + "more than 50% background is the same value; noise calculation cannot be completed." + bcolors.ENDC)

        # SNR = signal / noise
        snr_values[region_name] = {
            'Left': left_signal_avg / noise if noise != 0 else np.nan,
            'Right': right_signal_avg / noise if noise != 0 else np.nan
        }
    return snr_values, noise


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


def compute_snr_brain_mask(gray_mask_path, anat_data_path, output_path, root_name):
    """
    Computes SNR within a specified gray matter mask for single-frame data, plots the SNR histogram,
    and saves the average SNR to a text file.

    :param gray_mask_path: Path to binary gray matter mask NIfTI file
    :param anat_data_path: Path to single-frame anatomical data NIfTI file
    :param output_path: Directory path to save results
    :param root_name: Root name for output files
    :return: SNR value
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
        # Use anat_data as-is if it's already 3D
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
            print(
                bcolors.WARNING + 'std of noise values = 0, try with 20% of the lowest values in the img' + bcolors.ENDC)
            noise_threshold = np.percentile(non_zero_values, 20)
            noise_values = non_zero_values[non_zero_values <= noise_threshold]
            noise = np.std(noise_values)
            if noise == 0:
                print(
                    bcolors.WARNING + 'std of noise values = 0, try with 30% of the lowest values in the img' + bcolors.ENDC)
                noise_threshold = np.percentile(non_zero_values, 30)
                noise_values = non_zero_values[non_zero_values <= noise_threshold]
                noise = np.std(noise_values)
                if noise == 0:
                    print(
                        bcolors.WARNING + 'std of noise values = 0, try with 40% of the lowest values in the img' + bcolors.ENDC)
                    noise_threshold = np.percentile(non_zero_values, 40)
                    noise_values = non_zero_values[non_zero_values <= noise_threshold]
                    noise = np.std(noise_values)
                    if noise == 0:
                        print(
                            bcolors.WARNING + 'std of noise values = 0, try with 50% of the lowest values in the img' + bcolors.ENDC)
                        noise_threshold = np.percentile(non_zero_values, 50)
                        noise_values = non_zero_values[non_zero_values <= noise_threshold]
                        noise = np.std(noise_values)
                    else:
                        raise ValueError(
                            bcolors.FAIL + "more than 50% background is the same value; noise calculation cannot be completed." + bcolors.ENDC)
    else:
        raise ValueError(
            bcolors.FAIL + "10% of noise is just 0; noise calculation cannot be completed." + bcolors.ENDC)

    # Calculate SNR
    signal_avg = np.mean(mask_signal) if mask_signal.size > 0 else np.nan
    snr = signal_avg / noise if noise != 0 else np.nan

    return snr


def fwhm(anat_file, mask_file, dir_prepro, s_bind, afni_sif):
    """Calculate the FWHM of the input image using AFNI's 3dFWHMx.

    :type anat_file: str
    :param anat_file: The filepath to the anatomical image NIFTI file.
    :type mask_file: str
    :param mask_file: The filepath to the binary head mask NIFTI file.
    :type dir_prepro: str
    :param dir_prepro: Directory for preprocessing
    :type s_bind: str
    :param s_bind: Singularity bind paths
    :type afni_sif: str
    :param afni_sif: Path to AFNI singularity container
    :rtype: float
    :return: The combined FWHM value.
    """
    original_dir = os.getcwd()
    os.chdir(dir_prepro)

    try:
        command = f'singularity run' + s_bind + afni_sif + '3dFWHMx -overwrite ' \
                  f'-mask {mask_file} ' \
                  f'-input {anat_file}'
        fwhm_string_list = spgo(command)
        print(fwhm_string_list)

        # Use regex to find the line with the four numeric values
        match = re.search(r'(\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+)', fwhm_string_list)
        if match:
            retcode = str(match.group(1))
            print(retcode)
            vals = np.array(retcode.split(), dtype=np.float64)
            return vals[3]  # Return the combined FWHM value
        else:
            raise ValueError("Could not parse FWHM output from AFNI")

    except Exception as e:
        print(f"WARNING: FWHM calculation failed - {str(e)}")
        return None
    finally:
        os.chdir(original_dir)


def anat_QC(BASE_SS_coregistr, Ref_file, type_norm, labels_dir, dir_prepro, ID, listTimage, masks_dir,
            s_bind, afni_sif, diary_file):
    ct = datetime.datetime.now()
    nl = 'Run anatomical._16_anat_QC_SNR.anat_QC'
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n{nl}')

    for Timage in listTimage:
        direction = opj(labels_dir)
        atlas_filename_rsp = opj(direction, type_norm + 'atlaslvl1_LRrspQC' + Timage + '.nii.gz')
        atlas_filename_orig = opj(direction, type_norm + 'atlaslvl1_LR.nii.gz')
        atlas_filename = opj(direction, type_norm + 'atlaslvl1_LR_rsp_anat.nii.gz')
        lines = []
        anat_filename = opj(dir_prepro, ID + '_acpc_test_QC_' + Timage + '.nii.gz')
        brain_mask = opj(masks_dir, Timage + 'brain_mask_final_QCrsp.nii.gz')
        output_results = opj(dir_prepro, 'QC_anat')
        template_in_anat = opj(dir_prepro,'template_in_anat_DC.nii.gz')


        if not ope(output_results):
            os.mkdir(output_results)

        if ope(anat_filename):
            if not ope(atlas_filename_orig):
                nl = 'WARNING: no atlas lvl 1 LR found, this is a requirement for some of QC analysis'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            else:
                # Atlas labels as described by you
                labels = {2: '3rd ventricles', 3: '4th ventricles', 7: 'Cerebellum', 5: 'Cerebellum White',
                          6: 'Cortical White matter', 1: 'CSF', 4: 'Lateral ventricles', 8: 'Isocortex',
                          9: 'Allocortex', 10: 'Periallocortex', 11: 'Subcortical areas', 12: 'Diencephalon',
                          13: 'Brain stem'}

                # Example usage
                atlas_data = nib.load(atlas_filename).get_fdata()
                atlas_data = np.rint(atlas_data).astype(np.int32)
                anat_data = nib.load(anat_filename).get_fdata()
                snr_results, noise = compute_snr_single_frame(atlas_data, anat_data, labels)

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

                cnr_val = np.abs(average_snr_results['Gray_Matter'] - average_snr_results['White_Matter']) / noise
                cortical_contrast = (average_snr_results['White_Matter'] - average_snr_results['Gray_Matter']) / (
                        (average_snr_results['White_Matter'] + average_snr_results['Gray_Matter']) / 2)

                #################### QC that doesn't require atlaslvl1 ####################
                line_QC_func = []
                try:
                    fwhm_val = fwhm(anat_filename, atlas_filename, dir_prepro, s_bind, afni_sif)
                    print(fwhm_val)
                    line_QC_func.append(f"  fwhm_val: {fwhm_val}")
                except Exception as e:
                    print(bcolors.WARNING + f'FWHM calculation failed: {str(e)}' + bcolors.ENDC)
                    fwhm_val = None

                # Compute SNR with brain mask
                compute_snr_brain_mask_val = compute_snr_brain_mask(
                    brain_mask,
                    anat_filename,
                    output_results, Timage)
                line_QC_func.append(f"  snr_brain_mask: {compute_snr_brain_mask_val}")

                ## NMI index
                fixed = ants.image_read(Ref_file)
                moving = ants.image_read(template_in_anat)
                nmi = ants.image_mutual_information(fixed, moving)
                line_QC_func.append(f"  NMI_index: {nmi}")

                # Create the comprehensive QC figure
                create_qc_figure(output_results, Timage, snr_results, Ref_file,
                                template_in_anat, average_snr_results, noise,
                                anat_filename, BASE_SS_coregistr, fwhm_val, nmi,
                                compute_snr_brain_mask_val, cnr_val, cortical_contrast)

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