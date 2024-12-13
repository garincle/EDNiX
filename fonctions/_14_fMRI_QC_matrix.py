import nilearn
import subprocess
import os
from nilearn import image
from nilearn.input_data import NiftiLabelsMasker
import numpy as np
import nibabel as nib
from nilearn.connectome import ConnectivityMeasure
import pandas as pd
import matplotlib.pyplot as plt
from nilearn import plotting
import pingouin as pg
from pingouin import pairwise_ttests
from pingouin import ttest
import seaborn as sns
import shutil
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from scipy.spatial.distance import pdist, squareform
from fonctions.extract_filename import extract_filename
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

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath

opd = os.path.dirname
ope = os.path.exists

spgo = subprocess.getoutput

#################################################################################################
####Seed base analysis
#################################################################################################
def fMRI_QC_matrix(ID, Session, segmentation_name_list,
                   dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3,
                   specific_roi_tresh, unspecific_ROI_thresh, RS, nb_run, s_bind,afni_sif,diary_file):

    # This is a function to estimate functional connectivity specificity. See Grandjean 2020 for details on the reasoning

    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  Working on step ' + str(14) + '(function: _14_fMRI_QC_matrix).  ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    def specific_FC(specific_roi, unspecific_ROI):
        if (specific_roi >= specific_roi_tresh) and (unspecific_ROI < unspecific_ROI_thresh):
            cat = 'Specific'
        elif (specific_roi >= specific_roi_tresh) and (unspecific_ROI >= unspecific_ROI_thresh):
            cat = 'Unspecific'
        elif (abs(specific_roi) < specific_roi_tresh) and (abs(unspecific_ROI) < unspecific_ROI_thresh):
            cat = 'No'
        else:
            cat = 'Spurious'
        return cat

    for direction in [dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3]:

        out_results = opj(direction, '10_Results')

        if not ope(out_results):
            os.mkdir(out_results)

        if not ope(out_results):
            os.mkdir(out_results)
        out_results_V = opj(out_results, 'fMRI_QC_matrix')

        if ope(out_results_V):
            shutil.rmtree(out_results_V)

        os.mkdir(out_results_V)

        atlas                = 'atlaslvl3_LR.nii.gz'
        atlas_filenamelvl3LR = opj(direction, 'atlaslvl3_LR.nii.gz')
        atlas_filenameQC     = opj(direction, atlas[:-7] + '_forfMRIQC.nii.gz')

        if not ope(atlas_filenamelvl3LR):
            nl = 'ERROR: no atlas lvl 3 LR found, Step 14 will failed!'
            diary.write(f'\n{nl}')
            raise Exception(bcolors.FAIL + nl + bcolors.ENDC)

        else:
            #### Build the atlas for LR somato ant anterior cingulate
            string_build_atlas = str('')
            for new_label, old_label in zip([3, 3, 1, 2], [54, 1054, 58, 1058]):
                string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(
                    int(old_label)) + ')))+'
            string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"

            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas_filenamelvl3LR + ' -expr ' + string_build_atlas2 + ' -prefix ' + opj(
                direction, atlas[:-7] + '_forfMRIQC.nii.gz') + ' -overwrite'
            nl = spgo(command)
            diary.write(f'\n{nl}')
            print(nl)
            dictionary = {"Sources": atlas_filenamelvl3LR,
                          "Description": string_build_atlas2 + ' (3dcalc, AFNI).'},
            json_object = json.dumps(dictionary, indent=2)
            with open(opj(
                direction, atlas[:-7] + '_forfMRIQC.json'), "w") as outfile:
                outfile.write(json_object)

            for atlas_filename in [atlas_filenameQC, atlas_filenamelvl3LR]:
                for i in range(0, int(nb_run)):
                    root_RS = extract_filename(RS[i])
                    if direction == dir_fMRI_Refth_RS_prepro1:
                        func_filename = opj(direction, root_RS + '_residual.nii.gz')
                    elif direction == dir_fMRI_Refth_RS_prepro2:
                        func_filename = opj(direction, root_RS + '_residual_in_anat.nii.gz')
                    elif direction == dir_fMRI_Refth_RS_prepro3:
                        func_filename = opj(direction, root_RS + '_residual_in_template.nii.gz')

                    if ope(func_filename):
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
                                labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                                labeled_img2.to_filename(atlas_filename)
                        else:
                            nl = "No differences found in the specified fields."
                            print(bcolors.OKGREEN + nl + bcolors.ENDC)
                            diary.write(f'\n{nl}')
                    else:
                        nl = 'WARNING: ' + str(func_filename) + ' not found!!'
                        print(bcolors.WARNING + nl + bcolors.ENDC)
                        diary.write(f'\n{nl}')


            ##########################################################################
            ############     WORK in QC 1: Adapted Grandean's method:     ############
            ############     correlations LR somatoSenory compare to ACC  ############
            ##########################################################################

            for i in range(0, int(nb_run)):
                root_RS = extract_filename(RS[i])
                if direction == dir_fMRI_Refth_RS_prepro1:
                    func_filename = opj(direction, root_RS + '_residual.nii.gz')
                if direction == dir_fMRI_Refth_RS_prepro2:
                    func_filename = opj(direction, root_RS + '_residual_in_anat.nii.gz')
                if direction == dir_fMRI_Refth_RS_prepro3:
                    func_filename = opj(direction, root_RS + '_residual_in_template.nii.gz')

                if ope(func_filename):

                    output_results = opj(direction, '10_Results','fMRI_QC_matrix')

                    NAD_masker = NiftiLabelsMasker(labels_img=atlas_filenameQC,
                                                    detrend=False,
                                                    smoothing_fwhm=None,
                                                    standardize=False,
                                                    low_pass=None,
                                                    high_pass=None,
                                                    t_r=None,
                                                    memory=None, verbose=5)

                    time_series = NAD_masker.fit_transform(func_filename)
                    correlation_measure = ConnectivityMeasure(kind="correlation")
                    correlation_matrix  = correlation_measure.fit_transform([time_series])[0]

                    # Mask the main diagonal for visualization:
                    np.fill_diagonal(correlation_matrix, 0)

                    # matrices are ordered for block-like representation
                    nl = str(correlation_matrix.shape)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = str(correlation_matrix)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    specific_roi = correlation_matrix[0][1]
                    unspecific_ROI = correlation_matrix[0][2]

                    cat = specific_FC(specific_roi, unspecific_ROI)

                    line_specificity = ['QC result', 'specific_roi:', str(specific_roi), 'unspecific_ROI:', str(unspecific_ROI), 'Result:', str(cat)]
                    nl = str(line_specificity)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    ##########################################################################
                    ############                      WORK on QC 2:               ############
                    ############     correlation matrix for the entire lvl 3      ############
                    ##########################################################################

                    panda_file = segmentation_name_list[2]

                    panda_file_right = panda_file.copy()
                    panda_file_right['label'] += 1000
                    # Combine both DataFrames
                    panda_file_combined = pd.concat([panda_file, panda_file_right], ignore_index=True)

                    # Create combined region labels with side information
                    labels_with_side = ['L-' + region for region in panda_file['region']] + ['R-' + region for
                                                                                             region in
                                                                                             panda_file['region']]
                    # Combine the original regions twice (once for left, once for right)

                    # Create a new pandas DataFrame with the labels_with_side
                    # Create a new pandas DataFrame with the labels_with_side and original region names
                    panda_file_combined['region'] = labels_with_side

                    # Load your brain atlas data
                    atlas_img = nib.load(atlas_filenamelvl3LR)
                    atlas_data = atlas_img.get_fdata()

                    # Flatten the atlas data for easy comparison
                    atlas_flat = np.unique(atlas_data.flatten().astype(int))
                    atlas_flat[atlas_flat != 0]

                    # Check if labels in the segmentation DataFrame are in the atlas data
                    filtered_labels = panda_file_combined[panda_file_combined['label'].isin(atlas_flat)]
                    filtered_labels_list = list(filtered_labels['label'])
                    atlas_filtered_list = list(filtered_labels['region'])

                    # Create new lists to store the filtered values
                    final_labels_list = []
                    final_atlas_list = []

                    # Loop through the filtered_labels_list and filter both the labels and corresponding regions in atlas_filtered_list
                    for idx, num in enumerate(filtered_labels_list):
                        if int(num) < 1000:
                            # Check if num + 1000 exists in the list
                            if (int(num) + 1000) in filtered_labels_list:
                                final_labels_list.append(num)
                                final_atlas_list.append(atlas_filtered_list[idx])  # Add the corresponding region
                        else:
                            # Check if num - 1000 exists in the list
                            if (int(num) - 1000) in filtered_labels_list:
                                final_labels_list.append(num)
                                final_atlas_list.append(atlas_filtered_list[idx])  # Add the corresponding region

                    # Output the final filtered lists
                    nl = "Filtered labels list: " + str(final_labels_list)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = "Filtered atlas regions list:" + str(final_atlas_list)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    ### create a new list of label with zero if not in the new label list (for 3dcalc)
                    list_old_label = []
                    for lab in atlas_flat:
                        if lab in np.array(final_labels_list):
                            list_old_label.append(lab)
                        else:
                            list_old_label.append(0)
                    nl = str(list_old_label)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    string_build_atlas = str('')
                    for new_label, old_label in zip(list_old_label, list(atlas_flat)):
                        string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(
                            int(old_label)) + ')))+'
                    string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"

                    command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas_filename + ' -expr ' + string_build_atlas2 + ' -prefix ' + atlas_filename[:-7] + '_filtered.nii.gz' + ' -overwrite'
                    nl = spgo(command)
                    diary.write(f'\n{nl}')
                    print(nl)
                    dictionary = {"Sources": atlas_filename,
                                  "Description": string_build_atlas2 + ' (3dcalc, AFNI).'},
                    json_object = json.dumps(dictionary, indent=2)
                    with open(atlas_filename[:-7] + '_filtered.json', "w") as outfile:
                        outfile.write(json_object)

                    ##########################################################################
                    NAD_masker = NiftiLabelsMasker(
                        labels_img=atlas_filename[:-7] + '_filtered.nii.gz',
                        detrend=False,
                        smoothing_fwhm=None,
                        standardize=False,
                        low_pass=None,
                        high_pass=None,
                        t_r=None,
                        memory=None, verbose=5)
                    try:
                        time_series = NAD_masker.fit_transform(func_filename)
                    except:
                        dummy = nilearn.image.resample_to_img(
                            atlas_filename[:-7] + '_filtered.nii.gz', func_filename,
                            interpolation='nearest')
                        dummy.to_filename(atlas_filename[:-7] + '_filtered.nii.gz')
                        extracted_data = nib.load( atlas_filename[:-7] + '_filtered.nii.gz').get_fdata()
                        labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                        labeled_img2.to_filename(atlas_filename[:-7] + '_filtered.nii.gz')

                        NAD_masker = NiftiLabelsMasker(
                            labels_img=atlas_filename[:-7] + '_filtered.nii.gz',
                            detrend=False,
                            smoothing_fwhm=None,
                            standardize=False,
                            low_pass=None,
                            high_pass=None,
                            t_r=None,
                            memory=None, verbose=5)

                        time_series = NAD_masker.fit_transform(func_filename)

                    correlation_measure = ConnectivityMeasure(kind="correlation")
                    correlation_matrix = correlation_measure.fit_transform([time_series])[0]

                    # Plot the correlation matrix

                    # Make a large figure
                    # Mask the main diagonal for visualization:
                    np.fill_diagonal(correlation_matrix, 0)
                    # The labels we have start with the background (0), hence we skip the
                    # first label
                    # matrices are ordered for block-like representation
                    nl = str(correlation_matrix.shape)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = str(len(panda_file['region']))
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = str(correlation_matrix)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    plotting.plot_matrix(
                        correlation_matrix,
                        figure=(10, 8),
                        labels=final_atlas_list,
                        vmax=0.8,
                        vmin=-0.8)

                    plt.savefig(output_results + '/' + root_RS + 'lvl3LR_correlation_matrix.png')

                    # Save the connectivity matrix
                    matrix_filename = opj(output_results, root_RS + 'lvl3LR_correlation_matrix.npy')
                    np.save(matrix_filename, correlation_matrix)

                    # Convert correlation matrix to a DataFrame
                    corr_df = pd.DataFrame(correlation_matrix, index=final_atlas_list, columns=final_atlas_list)

                    connection = []
                    value = []
                    # Iterate through the matrix to populate the list
                    for i in range(len(final_atlas_list)):
                        for j in range(i + 1, len(final_atlas_list)):
                            connection.append(f"{final_atlas_list[i]} to {final_atlas_list[j]}")
                            value.append(correlation_matrix[i, j])

                    # Create DataFrame from the list
                    flattened_df = pd.DataFrame([value], columns=connection)
                    flattened_df['ID'] = ID
                    flattened_df['Session'] = Session

                    # Display the flattened DataFrame
                    nl = str(flattened_df)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    # Optionally, save the flattened DataFrame to a CSV file
                    flattened_df.to_csv(opj(output_results,root_RS + 'lvl3LR_correlation__correlation_matrix.csv'),
                                        index=False)

                    ##########################################################################
                    ###############################     R vs L Intra     #####################
                    ##########################################################################
                    ##########################################################################

                    # Extract the labels
                    left_labels = panda_file[panda_file['label'].isin(final_labels_list)]['label'].values
                    right_labels = panda_file_right[panda_file_right['label'].isin(final_labels_list)]['label'].values

                    # Initialize lists to store connectivity values
                    left_connectivity = []
                    # Iterate through each pair of left hemisphere labels
                    for i, left_label_A in enumerate(left_labels):
                        left_index_A = np.where(left_labels == left_label_A)[0][0]
                        for j, left_label_B in enumerate(left_labels):
                            left_index_B = np.where(left_labels == left_label_B)[0][0]
                            # Ensure that correlations between the same region are excluded
                            if i < j:  # This condition ensures that only unique pairs are considered
                                left_connectivity.append(correlation_matrix[left_index_A, left_index_B])
                    left_connectivity = np.array(left_connectivity)

                    right_connectivity = []
                    # Iterate through each pair of left hemisphere labels
                    for i, left_label_A in enumerate(right_labels):
                        left_index_A = np.where(right_labels == left_label_A)[0][0] + len(left_labels)
                        for j, left_label_B in enumerate(right_labels):
                            left_index_B = np.where(right_labels == left_label_B)[0][0] + len(left_labels)
                            # Ensure that correlations between the same region are excluded
                            if i < j:  # This condition ensures that only unique pairs are considered
                                right_connectivity.append(correlation_matrix[left_index_A, left_index_B])
                    right_connectivity = np.array(right_connectivity)

                    nl = str(len(right_connectivity))
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = str(len(left_connectivity))
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    # Ensure the arrays are of the same length
                    assert left_connectivity.shape == right_connectivity.shape, "Arrays must be of the same shape"

                    # Number of regions
                    num_regions = left_connectivity.shape[0]

                    # Prepare the data
                    data = {'Correlation': np.concatenate([left_connectivity, right_connectivity]),
                        'Hemisphere': ['Left'] * num_regions + ['Right'] * num_regions,
                        'Region': list(range(1, num_regions + 1)) * 2}

                    # Create the DataFrame
                    Correl_pd = pd.DataFrame(data)

                    nl = str(Correl_pd)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    # Now you can use this DataFrame with pairwise_ttests
                    posthocs = pairwise_ttests(dv='Correlation', within='Hemisphere', return_desc=True, interaction=True,
                                               within_first=False, subject='Region', data=Correl_pd, parametric=False)
                    nl =  str(posthocs)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    nl = f"T-statistic: {float(posthocs['W-val'])}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"P-value: {float(posthocs['p-unc'])}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    line_withinH_result = [f"W-statistic: {float(posthocs['W-val'])}", f"P-value: {float(posthocs['p-unc'])}"]
                    nl =  str(line_withinH_result)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')


                    ################## Figure

                    sns.set(style="ticks")
                    fig, ax1 = plt.subplots(1, 1, figsize=(6, 8))
                    ax1.set_ylim([-0.4, 1])
                    ax = pg.plot_paired(dv='Correlation', within='Hemisphere', ax=ax1, boxplot_in_front=True,
                                        pointplot_kwargs={'scale': 0.5, 'markers': '.'},
                                        subject='Region', data=Correl_pd,
                                        boxplot_kwargs={"palette": "Set1", 'linewidth': 2, 'width': 0.5, 'fliersize': 2.5})
                    ax.grid(False)
                    ax.set(xlabel=None, ylabel=None)
                    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=25)
                    ax.locator_params(axis="y", nbins=5)

                    ylabels = ['{:,.1f}'.format(x) for x in ax.get_yticks()]
                    ylabels_list = [float(i) for i in ylabels]
                    ax.set_yticks(ylabels_list)
                    ax.set_yticklabels(ylabels_list, fontsize=25)

                    for _, s in ax.spines.items():
                        s.set_linewidth(1.5)
                        s.set_color('black')

                    fig = ax.get_figure()
                    file_name = root_RS + 'paired_ttest_results_intraH.png'
                    file_path = opj(output_results, file_name)
                    nl = f"Saving figure to: {file_path}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"DPI: {400}, Bbox_inches: {'tight'}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    fig.savefig(file_path, dpi=400, bbox_inches='tight')

                    ##########################################################################
                    ###############################     LR vs Intra     ######################
                    ##########################################################################
                    ##########################################################################

                    # Initialize lists to store connectivity values and corresponding labels

                    right_to_left_connectivity = []
                    right_to_left_labels = []

                    # Iterate through each right hemisphere label
                    for right_label in right_labels:
                        # Find the corresponding left hemisphere label
                        left_label = right_label - 1000

                        # Find the indices of the right and left hemisphere labels
                        right_index = np.where(right_labels == right_label)[0][0] + len(left_labels)
                        left_index  = np.where(left_labels == left_label)[0][0]

                        # Extract the connectivity value between the right and left regions
                        connectivity_value = correlation_matrix[right_index, left_index]

                        # Append the connectivity value and corresponding labels to the lists
                        right_to_left_connectivity.append(connectivity_value)
                        right_to_left_labels.append(right_label)

                    # Convert the lists to numpy arrays
                    right_to_left_connectivity = np.array(right_to_left_connectivity)
                    right_to_left_labels       = np.array(right_to_left_labels)

                    # Convert the lists to numpy arrays
                    LR_connectivity = np.array(right_to_left_connectivity)    # < seems a repeatition of the previous line
                    LR_labels       = np.array(right_to_left_labels)

                    nl = str(right_to_left_connectivity)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = str(right_to_left_labels)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    # Combine left and right connectivity values for plotting
                    combined_connectivity = np.concatenate([np.concatenate((left_connectivity, right_connectivity), axis=None), LR_connectivity])

                    # Create labels for the bar plot
                    labels = ['IntraH Connectivity'] * len(np.concatenate((left_connectivity, right_connectivity), axis=None)) + ['LR Connectivity'] * len(
                        LR_connectivity)


                    # Prepare the data
                    data = {
                        'Correlation': combined_connectivity,
                        'Hemisphere': labels,
                        'Region': list(range(1, len(combined_connectivity) + 1))
                    }

                    # Create the DataFrame
                    Correl_pd = pd.DataFrame(data)

                    # define samples
                    group1 = Correl_pd[Correl_pd['Hemisphere'] == 'IntraH Connectivity']
                    group2 = Correl_pd[Correl_pd['Hemisphere'] == 'LR Connectivity']

                    # perform independent two sample t-test
                    posthocs = ttest(group1['Correlation'], group2['Correlation'])

                    nl = f"T-statistic: {float(posthocs['T'])}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"P-value: {float(posthocs['p-val'])}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    lines_withinH_vs_LR = [f"T-statistic: {float(posthocs['T'])}", f"P-value: {float(posthocs['p-val'])}"]
                    nl = str(lines_withinH_vs_LR)
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    ################## Figure

                    sns.set(style="ticks")

                    fig, ax1 = plt.subplots(1, 1, figsize=(6, 8))
                    ax1.set_ylim([-0.4, 1])
                    ax = sns.boxplot(data=Correl_pd, x="Hemisphere", y="Correlation", hue="Hemisphere")
                    ax.grid(False)
                    ax.set(xlabel=None, ylabel=None)

                    fig = ax.get_figure()
                    fig.savefig(opj(output_results, root_RS + 'ttest_results_LR_vs_intraH.png'), dpi=400, bbox_inches='tight')

                    ##########################################################################
                    ###############################     histogram     ########################
                    ##########################################################################
                    ##########################################################################

                    # Flatten the connectivity matrix and plot histogram
                    mask_corre = np.tri(correlation_matrix.shape[0], k=-1)
                    # Convert the binary mask to a boolean mask
                    boolean_mask = mask_corre.astype(bool)

                    # Extracting the masked (True) values
                    correlation_values = correlation_matrix[boolean_mask]

                    plt.figure(figsize=(10, 6))
                    plt.hist(correlation_values, bins=50, color='skyblue', edgecolor='black')
                    plt.title('Distribution of Correlation Values')
                    plt.xlabel('Correlation Coefficient')
                    plt.ylabel('Frequency')
                    plt.savefig(opj(output_results, root_RS + 'correlation_histogram.png'))

                    ##########################################################################
                    ##################     Network characteristics     ######################
                    ##########################################################################
                    ##########################################################################
                    # Eigenvalue analysis
                    def eigenvalue_analysis(matrix):
                        eigenvalues, _ = np.linalg.eig(matrix)
                        return np.sort(np.abs(eigenvalues))[::-1]

                    # Clustering stability
                    def clustering_stability(matrix, n_clusters, n_iter=10):
                        labels_list = []
                        for _ in range(n_iter):
                            kmeans = KMeans(n_clusters=n_clusters, random_state=None).fit(matrix)
                            labels_list.append(kmeans.labels_)
                        return np.mean([np.mean(labels == labels_list[0]) for labels in labels_list])

                    # Dunn Index calculation
                    def dunn_index(matrix, labels):
                        distances = squareform(pdist(matrix))
                        unique_labels = np.unique(labels)
                        intra_dists = [np.mean(distances[labels == label]) for label in unique_labels]
                        inter_dists = [np.min(distances[np.ix_(labels == label1, labels == label2)])
                                       for i, label1 in enumerate(unique_labels) for label2 in unique_labels[i + 1:]]
                        return np.min(inter_dists) / np.max(intra_dists)

                    # Spectral clustering and evaluation
                    def spectral_clustering(matrix, n_clusters):
                        clustering = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors',
                                                        assign_labels='kmeans', random_state=42)
                        labels = clustering.fit_predict(matrix)
                        silhouette = silhouette_score(matrix, labels)
                        davies_bouldin = davies_bouldin_score(matrix, labels)
                        calinski_harabasz = calinski_harabasz_score(matrix, labels)
                        dunn = dunn_index(matrix, labels)
                        return silhouette, davies_bouldin, calinski_harabasz, dunn
                    # Spectral clustering and evaluation
                    def spectral_clustering(matrix, n_clusters):
                        clustering = SpectralClustering(n_clusters=n_clusters, affinity='nearest_neighbors',
                                                        assign_labels='kmeans', random_state=42)
                        labels = clustering.fit_predict(matrix)
                        silhouette = silhouette_score(matrix, labels)
                        davies_bouldin = davies_bouldin_score(matrix, labels)
                        calinski_harabasz = calinski_harabasz_score(matrix, labels)
                        dunn = dunn_index(matrix, labels)
                        return silhouette, davies_bouldin, calinski_harabasz, dunn

                    # Perform analyses and print results
                    n_clusters = 7

                    # Eigenvalue analysis
                    # Provides insights into the structure of the matrix.
                    eigenvalues = eigenvalue_analysis(correlation_matrix)
                    nl = f"  Top 5 Eigenvalues: {eigenvalues[:5]}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    # Clustering stability
                    # Measures how consistent the clustering results are when clustering multiple times with different
                    # initial conditions. Here we choose 7 as in Yeo et al.
                    stability = clustering_stability(correlation_matrix, n_clusters)
                    nl = f"  Clustering Stability: {stability}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    # Spectral clustering evaluation
                    silhouette, davies_bouldin, calinski_harabasz, dunn = spectral_clustering(correlation_matrix, n_clusters)

                    nl = f"  Silhouette Score: {silhouette}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"  Davies-Bouldin Index: {davies_bouldin}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"  Calinski-Harabasz Index: {calinski_harabasz}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    #Dunn Index: Helps in identifying compact and well-separated clusters.
                    #The Dunn Index evaluates clustering quality by considering the ratio between the minimum inter-cluster distance and the maximum intra-cluster distance.
                    nl = f"  Dunn Index: {dunn}"
                    diary.write(f'\n{nl}')
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)

                    line_network = [f"  Top 5 Eigenvalues: {eigenvalues[:5]}", f"  Clustering Stability: {stability}", f"  Silhouette Score: {silhouette}",
                                        f"  Davies-Bouldin Index: {davies_bouldin}", f"  Calinski-Harabasz Index: {calinski_harabasz}", f"  Dunn Index: {dunn}"]
                    # Function to calculate distribution statistics
                    def calculate_distribution_stats(matrix):
                        flattened = matrix.flatten()
                        median = np.median(flattened)
                        variance = np.var(flattened)
                        min_val = np.min(flattened)
                        max_val = np.max(flattened)
                        return median, variance, min_val, max_val

                    # Function to calculate density of correlation coefficients in different ranges
                    def calculate_density(matrix):
                        flattened = matrix.flatten()
                        density_neg = np.sum((flattened >= -1) & (flattened < -0.05)) / len(flattened)
                        density_zero = np.sum((flattened >= -0.05) & (flattened <= 0.05)) / len(flattened)
                        density_pos = np.sum((flattened > 0.05) & (flattened <= 1)) / len(flattened)
                        return density_neg, density_zero, density_pos

                    # Distribution statistics
                    median, variance, min_val, max_val = calculate_distribution_stats(correlation_matrix)
                    nl = f"  Median: {median}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"  Variance: {variance}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"  Min: {min_val}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"  Max: {max_val}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    # Density analysis
                    density_neg, density_zero, density_pos = calculate_density(correlation_matrix)
                    nl = f"  Density (Negative Correlations): {density_neg}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"  Density (Around Zero): {density_zero}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')
                    nl = f"  Density (Positive Correlations): {density_pos}"
                    print(bcolors.OKGREEN + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

                    line_descriptive = [f"  Median: {median}", f"  Variance: {variance}", f"  Min: {min_val}",
                                        f"  Max: {max_val}", f"  Density (Negative Correlations): {density_neg}", f"  Density (Around Zero): {density_zero}", f"  Density (Positive Correlations): {density_pos}"]
                    lines = line_descriptive + line_withinH_result + lines_withinH_vs_LR + line_specificity

                    with open(output_results + '/'+ root_RS + 'QC_result.txt', 'w') as f:
                        for line in lines:
                            f.write(line)
                            f.write('\n')
                else:
                    nl = 'WARNING: ' + str(func_filename) + ' not found!!'
                    print(bcolors.WARNING + nl + bcolors.ENDC)
                    diary.write(f'\n{nl}')

    diary.write(f'\n')
    diary.close()

