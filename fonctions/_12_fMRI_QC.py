import nilearn
import os
from nilearn import image
from nilearn.input_data import NiftiLabelsMasker
import shutil
import subprocess
from shutil import copyfile
import math
import scipy
from nitime.lazy import scipy_linalg as linalg
import nitime.utils as utils
from fonctions.extract_filename import extract_filename
import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import json


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


opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
opa = os.path.abspath

spco = subprocess.check_output
spgo = subprocess.getoutput

#################################################################################################
####Seed base analysis
#################################################################################################
def fMRI_QC(correction_direction, dir_fMRI_Refth_RS_prepro1, RS, nb_run,
            s_bind, afni_sif,diary_file):

    # Several functions were adapted from
    # https://github.com/preprocessed-connectomes-project/quality-assessment-protocol
    # See http://preprocessed-connectomes-project.org/quality-assessment-protocol/#installing-the-qap-package for details on the reasoning
    # and https://hal.science/hal-02326351/file/2018_Neuron_PRIME_DE_NeuroResource_Milham.pdf

    # rajouter estimation nb de composant à partir de mélodic

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(12) + '(function: _12_fMRI_QC).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    dir_path      = dir_fMRI_Refth_RS_prepro1
    output_results = opj(dir_path, '10_Results', 'fMRI_QC_SNR')

    if not ope(opj(dir_path, '10_Results')):
        os.mkdir(opj(dir_path, '10_Results'))

    if ope(output_results):
        shutil.rmtree(output_results)
        os.mkdir(output_results)
    else:
        os.mkdir(output_results)

    atlas_filename = opj(dir_path, 'atlaslvl1_LR.nii.gz')

    for i in range(0, int(nb_run)):
        lines = []
        root_RS = extract_filename(RS[i])
        func_filename = opj(dir_path, root_RS + '_xdtrf_2ref.nii.gz')

        if ope(func_filename):
            if ope(atlas_filename) == False:
                nl = 'WARNING: no altlas lvl 1 LR found, this is a requirement for some of QC analysis'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            else:
                # Load the NIfTI images
                img1 = nib.load(atlas_filename)
                img2 = nib.load(func_filename)

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

                        nl = f"{field}:"
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        nl = f"  Image 1: {values[0]}"
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')
                        nl = f"  Image 2: {values[1]}"
                        print(bcolors.OKGREEN + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')

                        dummy = nilearn.image.resample_to_img(atlas_filename, func_filename, interpolation='nearest')
                        dummy.to_filename(atlas_filename)
                        extracted_data = nib.load(atlas_filename).get_fdata()
                        extracted_data = np.rint(extracted_data).astype(np.int32)
                        labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                        labeled_img2.to_filename(atlas_filename)
                else:
                    nl =  "No differences found in the specified fields."
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                # Atlas labels as described by you
                labels = {2: '3rd ventricles', 3: '4th ventricles', 7: 'Cerebellum', 5: 'Cerebellum White',
                          6: 'Cortical White matter', 1: 'CSF', 4: 'Lateral ventricles', 8: 'Isocortex',
                          9: 'Allocortex', 10: 'Periallocortex', 11: 'Subcortical areas', 12: 'Diencephalon',
                          13: 'Brain stem'}

                # Load the atlas and fMRI data

                atlas_img  = nib.load(atlas_filename)
                atlas_data = atlas_img.get_fdata()
                atlas_data = np.rint(atlas_data).astype(np.int32)

                fmri_img   = nib.load(func_filename)
                fmri_data  = fmri_img.get_fdata()

                # Function to calculate SNR for each region and hemisphere
                def compute_snr(atlas_data, fmri_data, labels, hemisphere_offset=1000):
                    time_points = fmri_data.shape[-1]
                    snr_values = {}

                    for label, region_name in labels.items():
                        snr_values[region_name] = {'Left': [], 'Right': []}

                        # Left hemisphere
                        left_mask   = atlas_data == label
                        left_signal = fmri_data[left_mask]  # Shape should be (voxels, time_points)

                        if left_signal.ndim == 1:
                            left_signal = left_signal[:, np.newaxis]  # Add a new axis if it's 1D

                        # Right hemisphere
                        right_mask   = atlas_data == label + hemisphere_offset
                        right_signal = fmri_data[right_mask]

                        if right_signal.ndim == 1:
                            right_signal = right_signal[:, np.newaxis]  # Add a new axis if it's 1D

                        for t in range(time_points):
                            # Signal at time point t (mean of region)
                            left_signal_t  = np.mean(left_signal[:, t])
                            right_signal_t = np.mean(right_signal[:, t])

                            # Compute noise from bottom 10% of the histogram for the whole brain at t, excluding zero values
                            brain_values    = fmri_data[:, :, :, t].flatten()
                            non_zero_values = brain_values[brain_values > 0]  # Exclude zeros
                            if len(non_zero_values) > 0:
                                noise_threshold = np.percentile(non_zero_values, 10)
                                noise_values    = non_zero_values[non_zero_values <= noise_threshold]
                                noise           = np.std(noise_values)
                            else:
                                raise ValueError(bcolors.FAIL + '10% of noise is just 0... so calculation of noise cannot be calculated correctly')

                            # SNR = signal / noise
                            left_snr  = left_signal_t / noise
                            right_snr = right_signal_t / noise

                            snr_values[region_name]['Left'].append(left_snr)
                            snr_values[region_name]['Right'].append(right_snr)

                    return snr_values, noise

                # Compute SNR for each region
                snr_results, noise = compute_snr(atlas_data, fmri_data, labels)

                # Create DataFrame and save to CSV
                data = {'Time': np.arange(fmri_data.shape[-1])}
                for region, hemispheres in snr_results.items():
                    data[f'{region}_Left']  = hemispheres['Left']
                    data[f'{region}_Right'] = hemispheres['Right']

                df = pd.DataFrame(data)
                df.to_csv(opj(output_results, root_RS + 'snr_regions.csv'), index=False)

                # Plot the results
                plt.figure(figsize=(10, 6))
                for region in labels.values():
                    plt.plot(df['Time'], df[f'{region}_Left'], label=f'{region} Left', linestyle='--')
                    plt.plot(df['Time'], df[f'{region}_Right'], label=f'{region} Right')

                plt.xlabel('Time')
                plt.ylabel('SNR')
                plt.title('SNR of Brain Regions Over Time')
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                plt.savefig(opj(output_results, root_RS + 'snr_plot.png'))
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
                            gray_snr.append(np.mean([snr_results[region]['Left'][t],
                                                     snr_results[region]['Right'][t]]))
                        avg_gray_snr.append(np.mean(gray_snr))

                        # Average white matter SNR
                        white_snr = []
                        for region in white_regions:
                            white_snr.append(np.mean([snr_results[region]['Left'][t],
                                                      snr_results[region]['Right'][t]]))
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

                # Plot Gray and White Matter SNR over time
                plt.figure(figsize=(10, 6))
                plt.plot(df['Time'], df['Gray_Matter'], label='Gray Matter', color='blue')
                plt.plot(df['Time'], df['White_Matter'], label='White Matter', color='green')

                plt.xlabel('Time')
                plt.ylabel('Average SNR')
                plt.title('Average SNR for Gray and White Matter Over Time')
                plt.legend()
                plt.tight_layout()

                # Save the plot to file
                plt.savefig(opj(output_results, root_RS + 'snr_gray_white_plot.png'))
                plt.close()

                # Save to a new CSV file with the averages of gray and white matter SNR
                output_avg_csv = opj(output_results,root_RS + 'snr_gray_white.csv')
                df[['Time', 'Gray_Matter', 'White_Matter']].to_csv(output_avg_csv, index=False)

                # Extract the average SNR across all time points for each region, including gray and white
                average_snr_results = {}

                # Compute the average SNR for each region
                for region, hemispheres in snr_results.items():
                    avg_left_snr  = np.mean(hemispheres['Left'])
                    avg_right_snr = np.mean(hemispheres['Right'])
                    average_snr_results[region] = {'Left': avg_left_snr, 'Right': avg_right_snr}

                # Add the average gray and white matter to the dictionary
                average_snr_results['Gray_Matter'] = np.mean(avg_gray_snr)
                average_snr_results['White_Matter'] = np.mean(avg_white_snr)

                # Print the average SNR for each region
                nl = "Average SNR for each region:"
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

                for region, values in average_snr_results.items():

                    if isinstance(values, dict):
                        nl = f"{region}: Left = {values['Left']:.2f}, Right = {values['Right']:.2f}"
                    else:
                        nl = f"{region}: Average SNR = {values:.2f}"

                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

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
                cortical_contrast = ((average_snr_results['White_Matter'] - average_snr_results['Gray_Matter'])
                                     / ((average_snr_results['White_Matter'] + average_snr_results['Gray_Matter']) / 2))

                line_QC_func_atlaslvl1 = [f"  cortical_contrast: {cortical_contrast}",
                                          f"  cnr: {cnr_val}",
                                          f"  average_snr_Gray_Matter: {str(average_snr_results['Gray_Matter'])}",
                                          f"  average_snr_White_Matter: {str(average_snr_results['White_Matter'])}",]

                ##### extract signal of Tsnr values in the test regions

                # Step 1: Extract the time series from the image

                for imageQC, QCexplain in zip([opj(dir_path, root_RS + '_xdtrfwS_stdev.nii.gz'),
                                               opj(dir_path, root_RS + '_xdtrfwS_tsnr1.nii.gz'),
                                               opj(dir_path, root_RS + '_xdtrfwS_tsnr2.nii.gz')],
                                              ['stdev', 'TSNRcvarinv', 'TSNR']):

                    nifti_img = image.load_img(imageQC)

                    masker = NiftiLabelsMasker(labels_img= atlas_filename,detrend=False,smoothing_fwhm=None,
                                               standardize=False,low_pass=None,high_pass=None,t_r=None,
                                               memory=None, verbose=5)


                    time_series = masker.fit_transform(nifti_img)

                    # Step 5: The result is a 2D array (time points x regions)

                    nl = "Time series shape (time points x regions): " + str(time_series.shape)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    mean_signal_per_region = np.mean(time_series, axis=1)
                    nl = "Mean signal per region: " + str(mean_signal_per_region)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    if QCexplain == 'stdev':
                        stdev = mean_signal_per_region
                    elif QCexplain == 'TSNRcvarinv':
                        TSNRcvarinv = mean_signal_per_region
                    elif QCexplain == 'TSNR':
                        TSNR = mean_signal_per_region

                line_TSNR = [f"  stdev signal: {stdev}", f"  TSNR signal: {TSNR}", f"  TSNR cvarinv signal: {TSNRcvarinv}"]

                lines.append(line_TSNR)
                lines.append(line_QC_func_atlaslvl1)
                lines.append(line_snr)

            #################### QC that doesn't require atlaslvl1 ####################

            line_QC_func  = []
            line_QC_func2 = []

            try:
                anat_file = opj(dir_path, root_RS + '_residual.nii.gz')
                mask_file = opj(dir_path,'maskDilat.nii.gz')

                command = 'singularity run' + s_bind + afni_sif + '3dFWHMx -overwrite -combined' + \
                              ' -mask ' + mask_file + ' -input ' + anat_file + \
                              ' -acf ' + anat_file[:-7] + '.acf.txt' + \
                              '> ' + anat_file[:-7] + '_FWHMx.txt'
                nl = spgo([command])
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

                df = pd.read_csv(anat_file[:-7] + '.acf.txt',sep='\s+',
                                 names=["a","b","c", "combined_estimation"])

                fwhm_val = np.float64(df["combined_estimation"][1])

                if not np.isnan(fwhm_val):
                    nl = str(fwhm_val)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    line_QC_func.append(f"  fwhm_val: {nl}")

                else:
                    nl = "[!] Something went wrong with AFNI's 3dFWHMx."
                    print(bcolors.FAIL + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

            except:
                nl = 'FWHM calculation failed (the fix is complicated, but it is not a big deal)'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            def load(func_file, mask_file):
                """
                Load the functional timeseries data from a NIFTI file into Nibabel data
                format, check/validate the data, and remove voxels with zero variance.
                """
                try:
                    func_img = nib.load(func_file)
                    mask_img = nib.load(mask_file)
                except:
                    raise Exception(bcolors.FAIL + 'ERROR' + bcolors.ENDC)

                mask = mask_img.get_fdata()
                mask_var_filtered = mask
                func = func_img.get_fdata().astype(float)

                # Calculate variance across time for each voxel
                variance_map = func.var(axis=-1)

                # Update the mask to exclude voxels with variance below the threshold
                mask_var_filtered[variance_map < 1e-6] = 0

                func = func[mask_var_filtered.nonzero()].T  # will have ntpts x nvoxs
                return func


            def fd_jenkinson(in_file, rmax=80., out_file=None, out_array=True):
                """Calculate Jenkinson's Mean Framewise Displacement (aka RMSD) and save
                the Mean FD values to a file.

                - Method to calculate Framewise Displacement (FD) calculations
                  (Jenkinson et al., 2002).
                - Implementation written by @ Krsna, May 2013.
                - Jenkinson FD from 3dvolreg's *.affmat12.1D file from -1Dmatrix_save
                  option input: subject ID, rest_number, name of 6 parameter motion
                  correction file (an output of 3dvolreg) output: FD_J.1D file
                - in_file should have one 3dvolreg affine matrix in one row - NOT the
                  motion parameters.

                :type in_file: str
                :param in_file: Filepath to the coordinate transformation output vector
                                of AFNI's 3dvolreg (generated by running 3dvolreg with
                                the -1Dmatrix_save option).
                :type rmax: float
                :param rmax: (default: 80.0) The default radius of a sphere that
                             represents the brain.
                :type out_file: str
                :param out_file: (default: None) The filepath to where the output file
                                 should be written.
                :type out_array: bool
                :param out_array: (default: False) Flag to return the data in a Python
                                  NumPy array instead of an output file.
                :rtype: str
                :return: (if out_array=False) The filepath to the output file containing
                         the Mean FD values.
                :rtype: NumPy array
                :return: (if out_array=True) An array of the output Mean FD values.
                """

                if out_file is None:
                    fname, ext = ops(opb(in_file))
                    out_file   = opa('%s_fdfile%s' % (fname, ext))

                # if in_file (coordinate_transformation) is actually the rel_mean output
                # of the MCFLIRT command, forward that file
                if 'rel.rms' in in_file:
                    copyfile(in_file, out_file)
                    return out_file

                try:
                    pm_ = np.genfromtxt(in_file)
                except:
                    raise Exception(bcolors.FAIL + 'ERROR: ' + bcolors.ENDC)

                original_shape = pm_.shape
                pm = np.zeros((pm_.shape[0], pm_.shape[1] + 4))
                pm[:, :original_shape[1]] = pm_
                pm[:, original_shape[1]:] = [0.0, 0.0, 0.0, 1.0]

                # rigid body transformation matrix
                T_rb_prev = np.matrix(np.eye(4))

                flag = 0
                X = [0]  # First timepoint
                for i in range(0, pm.shape[0]):
                    # making use of the fact that the order of aff12 matrix is "row-by-row"
                    T_rb = np.matrix(pm[i].reshape(4, 4))

                    if flag == 0:
                        flag = 1
                    else:
                        M = np.dot(T_rb, T_rb_prev.I) - np.eye(4)
                        A = M[0:3, 0:3]
                        b = M[0:3, 3]

                        FD_J = math.sqrt(
                            (rmax * rmax / 5) * np.trace(np.dot(A.T, A)) + np.dot(b.T, b))
                        X.append(FD_J)

                    T_rb_prev = T_rb

                try:
                    X = np.array(X).reshape(-1)
                    np.savetxt(out_file, X)
                except:
                    raise Exception(bcolors.FAIL + 'ERROR: ' + bcolors.ENDC)

                if out_array:
                    return np.array(X).reshape(-1)
                else:
                    return out_file

            fd_jenkinson_array = fd_jenkinson(opj(dir_path, root_RS + '.aff12.1D'), rmax=80.,
                                              out_file=opj(dir_path, root_RS + '.aff12_fdfile.1D'),
                                              out_array=True)

            try:
                # Define timepoints (assuming each value in fd_jenkinson_val represents one timepoint)
                timepoints = np.arange(len(fd_jenkinson_array))

                # Plot the Framewise Displacement values
                plt.figure(figsize=(10, 6))
                plt.plot(timepoints, fd_jenkinson_array, label='FD Jenkinson', linestyle='-', color='b')

                # Customize the plot similar to your reference example
                plt.xlabel('Timepoint')
                plt.ylabel('Framewise Displacement (FD)')
                plt.title('Jenkinson Framewise Displacement Over Time')
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.grid()
                # Optimize layout and save the plot if desired
                plt.tight_layout()
                plt.savefig(opj(output_results, root_RS + '_fd_jenkinson_plot.png'))
                plt.close()
                mean_fd = np.mean(fd_jenkinson_array)
                line_QC_func2.append(f"  mean_fd: {mean_fd}")
            except:
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                print(bcolors.WARNING + 'mean_fd calculation failed (the fix is complicated, but it is not a big deal)' + bcolors.ENDC)

            try:
                def global_correlation(func_reorient, func_mask):
                    """Calculate the global correlation (GCOR) of the functional timeseries.

                    - From "Correcting Brain-Wide Correlation Differences in Resting-State
                      fMRI", Ziad S. Saad et al. More info here:
                        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3749702

                    :type func_reorient: str
                    :param func_reorient: Filepath to the deobliqued, reoriented functional
                                          timeseries NIFTI file.
                    :type func_mask: str
                    :param func_mask: Filepath to the functional brain mask NIFTI file.
                    :rtype: float
                    :return: The global correlation (GCOR) value.
                    """


                    zero_variance_func = load(func_reorient, func_mask)

                    list_of_ts = zero_variance_func.transpose()

                    # get array of z-scored values of each voxel in each volume of the
                    # timeseries
                    demeaned_normed = []

                    for ts in list_of_ts:
                        demeaned_normed.append(scipy.stats.mstats.zscore(ts))

                    demeaned_normed = np.asarray(demeaned_normed)

                    # make an average of the normalized timeseries, into one averaged
                    # timeseries, a vector of N volumes
                    volume_list = demeaned_normed.transpose()

                    avg_ts = []

                    for voxel in volume_list:
                        avg_ts.append(voxel.mean())

                    avg_ts = np.asarray(avg_ts)

                    # calculate the global correlation
                    gcor = (avg_ts.transpose().dot(avg_ts)) / len(avg_ts)

                    return gcor

                gcor_val =  global_correlation(opj(dir_path, root_RS + '_xdtr_deob.nii.gz'),
                                               opj(dir_path, root_RS + '_mask_final_in_fMRI_orig.nii.gz'))
                nl = str(gcor_val)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                line_QC_func2.append(f"  gcor: {gcor_val}")
            except:
                nl = 'gcor calculation failed (the fix is complicated, but it is not a big deal)'
                diary.write(f'\n{nl}')
                print(bcolors.WARNING + nl + bcolors.ENDC)
                raise ValueError()

            try:
                def robust_stdev(func):
                    """Compute robust estimation of standard deviation.

                    :type func: Nibabel data
                    :param func: The functional timeseries data.
                    :rtype: float
                    :return: The standard deviation value.
                    """

                    lower_qs = np.percentile(func, 25, axis=0)
                    upper_qs = np.percentile(func, 75, axis=0)
                    stdev = (upper_qs - lower_qs) / 1.349
                    return stdev

                def ar_nitime(x, order=1, center=False):
                    """Derive a model of the noise present in the functional timeseries for
                    the calculation of the standardized DVARS.

                    - Borrowed from nipy.algorithms.AR_est_YW. aka "from nitime import
                      algorithms as alg".

                    :type x: Nibabel data
                    :param x: The vector of one voxel's timeseries.
                    :type order: int
                    :param order: (default: 1) Which lag of the autocorrelation of the
                                  timeseries to use in the calculation.
                    :type center: bool
                    :param center: (default: False) Whether to center the timeseries (to
                                   demean it).mean_dvar
                    :rtype: float
                    :return: The modeled noise value for the current voxel's timeseries.
                    """

                    if center:
                        x = x.copy()
                        x = x - x.mean()
                    r_m = utils.autocorr(x)[:order + 1]
                    Tm = linalg.toeplitz(r_m[:order])
                    y = r_m[1:]
                    ak = linalg.solve(Tm, y)
                    return ak[0]

                def ar1(func, method=ar_nitime):
                    """Apply the 'ar_nitime' function across the centered functional
                    timeseries.

                    :type func: Nibabel data
                    :param func: The functional timeseries data.
                    :type method: Python function
                    :param method: (default: ar_nitime) The algorithm to use to calculate AR1.
                    :rtype: NumPy array
                    :return: The vector of AR1 values.
                    """
                    func_centered = func - func.mean(0)
                    ar_vals = np.apply_along_axis(method, 0, func_centered)
                    return ar_vals

                def calc_dvars(func_file, mask_file, output_all=False):
                    """Calculate the standardized DVARS metric.

                    :type func_file: str
                    :param func_file: The filepath to the NIFTI file containing the functional
                                       timeseries.
                    :type mask_file: str
                    :param mask_file: The filepath to the NIFTI file containing the binary
                                      functional brain mask.
                    :type output_all: bool
                    :param output_all: (default: False) Whether to output all versions of
                                       DVARS measure (non-standardized, standardized and
                                       voxelwise standardized).
                    :rtype: NumPy array
                    :return: The output DVARS values vector.
                    """

                    # load data
                    func = load(func_file, mask_file)

                    # Robust standard deviation
                    func_sd = robust_stdev(func)

                    # AR1
                    func_ar1 = ar1(func)

                    # Predicted standard deviation of temporal derivative
                    func_sd_pd = np.sqrt(2 * (1 - func_ar1)) * func_sd
                    diff_sd_mean = func_sd_pd.mean()

                    # Compute temporal difference time series
                    func_deriv = np.diff(func, axis=0)

                    # DVARS
                    # (no standardization)
                    dvars_plain = func_deriv.std(1, ddof=1)  # TODO: Why are we not ^2 this & getting the sqrt?
                    # standardization
                    dvars_stdz = dvars_plain / diff_sd_mean
                    # voxelwise standardization
                    diff_vx_stdz = func_deriv / func_sd_pd
                    dvars_vx_stdz = diff_vx_stdz.std(1, ddof=1)

                    if output_all:
                        try:
                            out = np.vstack((dvars_stdz, dvars_plain, dvars_vx_stdz))
                        except:
                            raise Exception(bcolors.FAIL + 'ERROR: ' + bcolors.ENDC)
                    else:
                        try:
                            out = dvars_stdz.reshape(len(dvars_stdz), 1)
                        except:
                            raise Exception(bcolors.FAIL + 'ERROR: ' + bcolors.ENDC)

                    return out

                calc_dvars_array = calc_dvars(opj(dir_path, root_RS + '_xdtr_deob.nii.gz'),
                                              opj(dir_path, root_RS + '_mask_final_in_fMRI_orig.nii.gz'))

                timepoints = np.arange(len(calc_dvars_array))
                plt.figure(figsize=(10, 6))
                plt.plot(timepoints, calc_dvars_array, label='Dvars', linestyle='-', color='b')

                # Customize the plot similar to your reference example
                plt.xlabel('Timepoint')
                plt.ylabel('dvars')
                plt.title('Dvars Over Time')
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.grid()

                # Optimize layout and save the plot if desired
                plt.tight_layout()
                plt.savefig(opj(output_results,root_RS + '_dvars_plot.png'))
                plt.close()
                mean_dvars_val = np.mean(calc_dvars_array)
                line_QC_func2.append(f"  mean_dvar: {mean_dvars_val}")
            except:
                nl = 'mean_dvar calculation failed (the fix is complicated, but it is not a big deal)'
                diary.write(f'\n{nl}')
                raise ValueError(bcolors.WARNING + nl + bcolors.ENDC)


            def plot_motion_parameters(motion_file, output_results):
                # Load the motion parameters
                motion_params = np.genfromtxt(motion_file)
                # Check if motion_params has the right shape (n_frames x 6)
                if motion_params.shape[1] != 6:
                    raise ValueError(
                        "Expected motion parameters file with 6 columns (3 translations, 3 rotations).")
                # Create time vector based on the number of frames
                time = np.arange(motion_params.shape[0])
                # Create a figure and axis
                plt.figure(figsize=(12, 8))
                # Plot each motion parameter
                labels = ['Translational X', 'Translational Y', 'Translational Z',
                          'Rotational Pitch', 'Rotational Roll', 'Rotational Yaw']
                for i in range(6):
                    plt.plot(time, motion_params[:, i], label=labels[i])
                # Customize the plot
                plt.xlabel('Time (frames)')
                plt.ylabel('Motion Parameters')
                plt.title('Motion Parameters Over Time')
                plt.axhline(y=0, color='k', linestyle='--', lw=0.5)  # Add a horizontal line at y=0
                plt.legend()
                plt.tight_layout()

                # Save the plot
                plt.savefig(opj(output_results, root_RS + '_motion_parameters_plot.png'))
                plt.close()

            # Example usage
            plot_motion_parameters(opj(dir_path, root_RS + '_dfile.1D'), output_results)

            def ghost_direction(epi_data_path, mask_data_path, direction="y", ref_file=None,
                                out_file=None):
                """Calculate the Ghost to Signal Ratio of EPI images.

                - GSR from Giannelli 2010. More details here:
                    https://www.ncbi.nlm.nih.gov/pubmed/21081879
                - This should be used for EPI images where the phase encoding direction
                  is known.

                :type epi_data: Nibabel data
                :param epi_data: The mean of the functional timeseries.
                :type mask_data: Nibabel data
                :param mask_data: The functional brain binary mask data.
                :type direction: str
                :param direction: (default: 'y') The phase encoding direction of the EPI
                                  image.
                :type ref_file: str
                :param ref_file: (default: None) If you are saving the Nyquist ghost mask,
                                  this is the filepath of the reference file to use to
                                 populate the header of the ghost mask NIFTI file.
                :type out_file: str
                :param out_file: (default: None) If you are saving the Nyquist ghost mask,
                                  this is the filepath to the ghost mask NIFTI file.
                :rtype: float
                :return: The ghost-to-signal ratio (GSR) value.
                """

                # first we need to make a nyquist ghost mask, we do this by circle
                # shifting the original mask by N/2 and then removing the intersection
                # with the original mask
                epi_img  = nib.load(epi_data_path)
                epi_data = epi_img.get_fdata()

                mask_img = nib.load(mask_data_path)
                mask_data = mask_img.get_fdata()

                n2_mask_data = np.zeros_like(mask_data)

                # rotate by n/2
                if direction == "x":
                    n2 = int(np.floor(mask_data.shape[0] / 2))
                    n2_mask_data[:n2, :, :] = mask_data[n2:(n2 * 2), :, :]
                    n2_mask_data[n2:(n2 * 2), :, :] = mask_data[:n2, :, :]
                elif direction == "y":
                    n2 = int(np.floor(mask_data.shape[1] / 2))
                    n2_mask_data[:, :n2, :] = mask_data[:, n2:(n2 * 2), :]
                    n2_mask_data[:, n2:(n2 * 2), :] = mask_data[:, :n2, :]
                elif direction == "z":
                    n2 = int(np.floor(mask_data.shape[2] / 2))
                    n2_mask_data[:, :, :n2] = mask_data[:, :, n2:(n2 * 2)]
                    n2_mask_data[:, :, n2:(n2 * 2)] = mask_data[:, :, :n2]
                else:
                    raise Exception("Unknown direction %s, should be x, y, or z" \
                                    % direction)

                # now remove the intersection with the original mask
                n2_mask_data = n2_mask_data * (1 - mask_data)

                # now create a non-ghost background region, that contains 2s
                n2_mask_data = n2_mask_data + 2 * (1 - n2_mask_data - mask_data)

                # Save mask
                if ref_file is not None and out_file is not None:
                    ref = nib.load(ref_file)
                    out = nib.Nifti1Image(n2_mask_data, ref.affine, ref.header)
                    out.to_filename(out_file)

                # now we calculate the Ghost to signal ratio, but here we define signal
                # as the entire foreground image
                gsr = (epi_data[n2_mask_data == 1].mean() - epi_data[n2_mask_data == 2].mean()) / epi_data[
                    n2_mask_data == 0].mean()

                return gsr

            if correction_direction in ['y', 'y-', 'z', 'z-', 'x', 'x-']:
                if correction_direction in ['y', 'y-']:
                    direct_aqc = 'y'
                elif correction_direction in ['z', 'z-']:
                    direct_aqc = 'z'
                elif correction_direction in ['x', 'x-']:
                    direct_aqc = 'x'

                ghost_direction_val = ghost_direction(opj(dir_path,root_RS + '_xdtr_deob.nii.gz'),
                                                      opj(dir_path, root_RS + '_mask_final_in_fMRI_orig.nii.gz'),
                                                      direction=direct_aqc,
                                                      ref_file=opj(dir_path,root_RS + '_xdtr_deob.nii.gz'),
                                                      out_file=opj(output_results,root_RS + '_ghost_mask.nii.gz'))

                dictionary = {"Sources": [opj(dir_path, root_RS + '_xdtr_deob.nii.gz'),
                                          opj(dir_path, root_RS + '_mask_final_in_fMRI_orig.nii.gz')],
                              "Description": 'ghost_direction.'},
                json_object = json.dumps(dictionary, indent=2)
                with open(opj(output_results,root_RS + '_ghost_mask.json'), "w") as outfile:
                    outfile.write(json_object)


                line_QC_func.append(f"  ghost_acq_direction: {ghost_direction_val}")
                nl = str(line_QC_func)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                nl = str(ghost_direction_val)
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')

            def ghost_all(epi_data, mask_data):
                """Call the 'ghost_direction' function on all possible phase encoding
                directions.

                :type epi_data: Nibabel data
                :param epi_data: The mean of the functional timeseries.
                :type mask_data: Nibabel data
                :param mask_data: The functional brain binary mask data.
                :rtype: tuple
                :return: The ghost-to-signal ratios (GSR) of each phase encoding direction.
                """

                directions = ["x", "y"]
                gsrs = [ghost_direction(epi_data, mask_data, d) for d in directions]

                return tuple(gsrs + [None])

            ghost_all_val = ghost_all(opj(dir_path,root_RS + '_xdtr_deob.nii.gz'),
                                    opj(dir_path, root_RS + '_mask_final_in_fMRI_orig.nii.gz'))

            line_QC_func.append(f"  ghost_all_direction: {ghost_all_val}")
            nl = str(line_QC_func)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            nl = str(ghost_all_val)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            def compute_snr_brain_mask(gray_mask_path, fmri_data_path, output_path, root_name):
                """
                Computes SNR within a specified gray matter mask, plots the SNR over time, and saves the average SNR.

                :param gray_mask_path: Path to binary gray matter mask NIfTI file
                :param fmri_data_path: Path to 4D fMRI data NIfTI file
                :param output_path: Directory path to save results
                :param root_name: Root name for output files
                """
                # Load the data
                gray_mask_img = nib.load(gray_mask_path)
                fmri_img = nib.load(fmri_data_path)
                gray_mask = gray_mask_img.get_fdata()

                # Ensure the mask is binary (0 or 1)
                gray_mask = np.rint(gray_mask).astype(np.int32)
                fmri_data = fmri_img.get_fdata()

                time_points = fmri_data.shape[-1]
                snr_values = []

                # Iterate over time points to compute SNR for each
                for t in range(time_points):
                    # Extract the 3D fMRI volume at time point t
                    fmri_volume = fmri_data[..., t]

                    # Signal within gray matter mask at this time point
                    mask_signal = fmri_volume[gray_mask == 1]

                    # Compute noise from bottom percentile values in the whole brain, excluding zeros
                    brain_values = fmri_volume[fmri_volume > 0].flatten()  # Exclude zeros directly
                    noise = None
                    for threshold in [10, 20, 30, 40, 50]:  # Try different percentiles
                        if len(brain_values) == 0:
                            raise ValueError("No non-zero values found in the brain data; check your fMRI data.")

                        noise_threshold = np.percentile(brain_values, threshold)
                        noise_values    = brain_values[brain_values <= noise_threshold]
                        noise           = np.std(noise_values)

                        if noise > 0:
                            break  # Valid noise found, break out of loop
                        else:
                            print(f"Warning: Noise std is zero at {threshold}% threshold, trying higher threshold.")

                    if noise == 0:
                        raise ValueError("Unable to calculate noise; more than 50% of background values are identical.")

                    # Calculate SNR for the time point
                    snr = np.mean(mask_signal) / noise
                    snr_values.append(snr)

                # Save SNR values
                snr_output_path = opj(output_path, root_name + '_snr_values.npy')
                np.save(snr_output_path, snr_values)
                print(f"SNR values saved to {snr_output_path}")

                # Plot SNR over time
                plt.figure(figsize=(10, 6))
                plt.plot(snr_values, label='SNR over Time')
                plt.xlabel('Time Points')
                plt.ylabel('SNR')
                plt.title('SNR within Gray Matter Mask Over Time')
                plt.legend()
                plt.grid(True)
                plt.tight_layout()

                # Save the plot
                plot_output_path = opj(output_path, root_name + '_snr_plot.png')
                plt.savefig(plot_output_path)
                print(f"SNR plot saved to {plot_output_path}")

                # Return the mean SNR across time points
                mean_snr = np.mean(snr_values)
                print(f"Average SNR: {mean_snr}")
                return mean_snr

            # Example usage
            compute_snr_brain_mask_val = compute_snr_brain_mask(
                opj(dir_path, root_RS + '_mask_final_in_fMRI_orig.nii.gz'),
                opj(dir_path, root_RS + '_xdtr_deob.nii.gz'),
                output_results, root_RS)

            line_QC_func.append(f"  snr_brain_mask {compute_snr_brain_mask_val}")

            ####### motion metrics #####

            # Step 1: Load the motion metrics
            motion_enorm = np.loadtxt(opj(dir_path, root_RS + 'motion_enorm.1D'))   # Euclidean norm values
            derivatives  = np.loadtxt(opj(dir_path, root_RS + '_xdtr_deriv.1D'))    # Derivatives (velocity)
            censor_1d    = np.loadtxt(opj(dir_path, root_RS + '_xdtr_censor.1D'))   # Censored time points
            outcount     = np.loadtxt(opj(dir_path, root_RS + '_xdt_outcount.r.1D'),
                                      skiprows=2)                                   # Censored time points

            # Step 2: Calculate Average Euclidean Norm
            avg_enorm = np.mean(motion_enorm)

            # Step 3: Calculate the Censor Fraction

            # The censor file has 1s for good time points and 0s for censored ones
            censor_fraction = 1 - np.mean(censor_1d)  # Censored fraction is the inverse of the mean (1 means kept)

            # Step 4: Calculate Average Absolute Velocity (from the derivatives)
            # The derivatives file has six columns (3 translations, 3 rotations), we average their absolute values
            avg_outcount = np.mean(outcount)
            avg_velocity = np.mean(np.abs(derivatives), axis=0).mean()

            nl = f"  avg motion velocity: {avg_velocity}"
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            nl = f"  avg motion magnitude: {avg_enorm}"
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            nl = f"  censor_fraction: {censor_fraction}"
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')
            nl = f"  outlier proportions: {avg_outcount}"
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            line_motion = [f"  avg motion velocity: {avg_velocity}", f"  avg motion magnitude: {avg_enorm}", f"  censor_fraction: {censor_fraction}", f"  outlier proportions: {avg_outcount}"]

            lines.append(line_motion)
            lines.append(line_QC_func)
            lines.append(line_QC_func2)
            flattened_list = [item for sublist in lines for item in sublist]

            nl = 'QC will look like ' + str(flattened_list)
            print(bcolors.OKGREEN + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

            with open(opj(output_results, root_RS + 'QC_result.txt'), 'w') as f:
                for line in flattened_list:
                    f.write(line)
                    f.write('\n')
        else:
            nl = 'WARNING: ' + str(func_filename) + ' not found!!'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()