import nilearn
import subprocess
import pandas as pd
import os
from nilearn.input_data import NiftiLabelsMasker
import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting
import nibabel as nib
from nilearn.connectome import ConnectivityMeasure
from fonctions.extract_filename import extract_filename

# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
spco = subprocess.check_output
opd = os.path.dirname



#################################################################################################
####Seed base analysis
#################################################################################################
def correl_matrix(dir_fMRI_Refth_RS_prepro1, RS, nb_run, selected_atlases_matrix, segmentation_name_list, ID, Session, bids_dir,s_bind,afni_sif):

    out_results = opj(bids_dir, 'Results')
    if not os.path.exists(out_results): os.mkdir(out_results)
    out_results_V = opj(out_results, 'fMRI_matrix')
    if not os.path.exists(out_results_V): os.mkdir(out_results_V)

    for panda_file, atlas in zip(segmentation_name_list, selected_atlases_matrix):
        for i in range(0, int(nb_run)):
            root_RS = extract_filename(RS[i])
            func_filename = opj(dir_fMRI_Refth_RS_prepro1, root_RS + '_residual.nii.gz')
            atlas = extract_filename(atlas)
            atlas_filename = opj(dir_fMRI_Refth_RS_prepro1, atlas + '.nii.gz')
            output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results')
            if not os.path.exists(output_results): os.mkdir(output_results)
            output_results = opj(dir_fMRI_Refth_RS_prepro1, '10_Results/correl_matrix')
            if not os.path.exists(output_results): os.mkdir(output_results)
            print(panda_file['region'])

            # Load the NIfTI images
            img1 = nib.load(atlas_filename)
            img2 = nib.load(func_filename)

            # Extract headers
            header1 = img1.header
            header2 = img2.header

            # Compare specific fields
            fields_to_compare = ['dim', 'pixdim']
            tolerance = 0.00006
            differences = {}
            for field in fields_to_compare:
                value1 = header1[field][:4]
                value2 = header2[field][:4]
                if not np.allclose(value1, value2, atol=tolerance):
                    differences[field] = (value1, value2)

            # Print differences
            if differences:
                print("Differences found in the following fields:")
                for field, values in differences.items():
                    print(f"{field}:")
                    print(f"  Image 1: {values[0]}")
                    print(f"  Image 2: {values[1]}")
                    caca = nilearn.image.resample_to_img(atlas_filename, func_filename, interpolation='nearest')
                    caca.to_filename(atlas_filename)
                    extracted_data = nib.load(atlas_filename).get_fdata()
                    labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                    labeled_img2.to_filename(atlas_filename)
            else:
                print("No differences found in the specified fields.")

            # Load your brain atlas data
            atlas_img = nib.load(atlas_filename)
            atlas_data = atlas_img.get_fdata()

            # Flatten the atlas data for easy comparison
            atlas_flat = np.unique(atlas_data.flatten().astype(int))
            atlas_flat[atlas_flat != 0]

            # Check if labels in the segmentation DataFrame are in the atlas data
            filtered_labels = panda_file[panda_file['label'].isin(atlas_flat)]
            atlas_filtered_list = list(filtered_labels['region'])
            print("label that you provided")
            print(panda_file)

            print("label found in the atlas after extraction of the signal :" + str(atlas))
            print(atlas_flat)

            print("if we combine them (analyse similarity we ware able to keep:")
            print(filtered_labels)
            print("corresponding to the labels")
            print(atlas_filtered_list)

            print("if they are different from the one you provided, you should investigate why!: region too small for fMRI? error in the labeling? error in the coregistration?, this is not always a big dill but you should know what happend!!")



            ### creat a new list of label with zero if not in the new label list (for 3dcalc)
            filtered_labels_list = list(filtered_labels['label'])
            list_old_label = []
            for lab in atlas_flat:
                if lab in np.array(filtered_labels_list):
                    list_old_label.append(lab)
                else:
                    list_old_label.append(0)

            string_build_atlas = str('')
            for new_label, old_label in zip(list_old_label, list(atlas_flat)):
                string_build_atlas = string_build_atlas + '(' + str(int(new_label)) + '*(equals(a,' + str(
                    int(old_label)) + ')))+'
            string_build_atlas2 = "'" + string_build_atlas[:-1] + "'"

            command = 'singularity run' + s_bind + afni_sif + '3dcalc' + ' -a ' + atlas_filename + ' -expr ' + string_build_atlas2 + ' -prefix ' + opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered.nii.gz') + ' -overwrite'
            spco(command, shell=True)

            ##########################################################################
            NAD_masker = NiftiLabelsMasker(labels_img=opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered.nii.gz'),
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
                caca = nilearn.image.resample_to_img(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz'), func_filename, interpolation='nearest')
                caca.to_filename(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz'))
                extracted_data = nib.load(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz')).get_fdata()
                labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                labeled_img2.to_filename(opj(dir_fMRI_Refth_RS_prepro1, atlas + '_filtered.nii.gz'))

                NAD_masker = NiftiLabelsMasker(
                    labels_img=opj(dir_fMRI_Refth_RS_prepro1, atlas + '_run_' + str(i) + '_filtered.nii.gz'),
                    detrend=False,
                    smoothing_fwhm=None,
                    standardize=False,
                    low_pass=None,
                    high_pass=None,
                    t_r=None,
                    memory=None, verbose=5)

                time_series = NAD_masker.fit_transform(func_filename)

            correlation_measure = ConnectivityMeasure(kind="partial correlation")
            correlation_matrix = correlation_measure.fit_transform([time_series])[0]

            # Plot the correlation matrix

            # Make a large figure
            # Mask the main diagonal for visualization:
            np.fill_diagonal(correlation_matrix, 0)
            # The labels we have start with the background (0), hence we skip the
            # first label
            # matrices are ordered for block-like representation
            print(correlation_matrix.shape)
            print(len(panda_file['region']))
            print(correlation_matrix)

            plotting.plot_matrix(
                correlation_matrix,
                figure=(10, 8),
                labels=atlas_filtered_list,
                vmax=0.8,
                vmin=-0.8,
            )

            plt.savefig(output_results + '/' + atlas + '_run_' + str(i) + 'matrix.png')

            # Convert correlation matrix to a DataFrame
            corr_df = pd.DataFrame(correlation_matrix, index=atlas_filtered_list, columns=atlas_filtered_list)

            connection = []
            value  = []
            # Iterate through the matrix to populate the list
            for i in range(len(atlas_filtered_list)):
                for j in range(i + 1, len(atlas_filtered_list)):
                    connection.append(f"{atlas_filtered_list[i]} to {atlas_filtered_list[j]}")
                    value.append(correlation_matrix[i, j])

            print(connection)
            print(value)
            print(len(connection))
            print(len(value))
            # Create DataFrame from the list
            flattened_df = pd.DataFrame([value], columns=connection)
            flattened_df['ID'] = ID
            flattened_df['Session'] = Session

            # Display the flattened DataFrame
            print(flattened_df)

            # Optionally, save the flattened DataFrame to a CSV file
            flattened_df.to_csv(output_results + '/' + atlas + '_run_' + str(i) + '_flattened_correlation_matrix.csv', index=False)

