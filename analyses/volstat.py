#import
import os
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import Tools.Load_EDNiX_requirement

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
import warnings
from os.path import join as opj

warnings.filterwarnings('ignore')
##########################################
########### Subject loader################
##########################################
#https://bids-standard.github.io/pybids/reports/index.html

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
spco = subprocess.check_output
spgo = subprocess.getoutput

#################################################################################################
####first parameters
#################################################################################################

##################################################################find the values that are droped because too small
def extractVol(MAIN_PATH, FS_dir, allinfo_study_c, regressor_list, all_ID, all_Session, all_data_path, type_norm, segmentation_name_list,
               segmentation_ID_list, atlas_names_Seg_list, list_atlases, bids_dir):

    s_path, afni_sif, fsl_sif, fs_sif, itk_sif, wb_sif, strip_sif, s_bind = Tools.Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, bids_dir, FS_dir)
    for segmentation, IDSEG, atlas_ref, atlas_names_Seg in zip(atlas_names_Seg_list, segmentation_ID_list, list_atlases, segmentation_name_list):

        atlasval = nib.load(atlas_ref).get_fdata()
        label_atlas = np.unique(atlasval)
        listexlidedvalues = np.array([])
        #### Anat_warp
        ####################################################################################
        ########################## Start the pipeline !!!!!!!!!!!!!!!!!!!!!!   #############
        ####################################################################################

        for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):

            # The anatomy
            path_anat    = opj(data_path ,'anat/')
            dir_transfo  = opj(path_anat ,'matrices')

            dir_native    = opj(path_anat ,'native')
            dir_prepro    = opj(dir_native ,'01_preprocess')
            wb_native_dir = opj(dir_native ,'02_Wb')
            volumes_dir   = opj(wb_native_dir ,'volumes')
            labels_dir    = opj(volumes_dir ,'labels')
            masks_dir     = opj(volumes_dir ,'masks')

            label_img = opj(labels_dir ,type_norm + segmentation)
            txtfile = opj(labels_dir ,type_norm + segmentation[:-7] + '.txt')

            command = 'singularity run' + s_bind + afni_sif + '3dhistog -int ' + label_img
            table = spgo(command)
            # Remove the first line (AFNI version info)
            table = table.split('\n', 1)[1]  # Split on first newline and keep everything after it

            with open(txtfile, 'w') as myfile:
                myfile.write(table)

            data = pd.read_table(txtfile, delim_whitespace=True)
            data2 = data.rename(index=str, columns={"#Magnitude": "label", "Freq": "count"})
            data2.drop(['Cum_Freq'], axis="columns")
            dim = np.prod(nib.load(label_img).header.get_zooms())
            data2['volume'] = data2['count'] * dim

            volume = data2['volume'].where(data2['volume' ] >0).dropna()
            # list of the index in the volume file
            list_volume_index = volume.index.values.tolist()
            # transform string into int()
            list_volume_index = np.array([int(i) for i in list_volume_index])
            # concat previous list with list of the label in the atlas
            arraycombine = np.array(np.concatenate((list_volume_index, label_atlas), axis=0))
            # find values that are not repeated
            uniqueValues, occurCount = np.unique(arraycombine, return_counts=True)
            uniquevalues_array = label_atlas[np.where(occurCount == 1)]
            # list them for each sub
            listexlidedvalues = np.append(listexlidedvalues, uniquevalues_array)

        # find all the values that are remvoved from the coregistration
        uniqueValues = np.unique(listexlidedvalues)
        len(uniqueValues)
        # transform string into int()
        uniqueValues = [int(i) for i in uniqueValues]

        ##################################################################create a panda array with all the volume values
        all_data_volumes = pd.DataFrame()
        #### Anat_warp

        for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):

            # The anatomy
            path_anat    = opj(data_path ,'anat/')
            dir_native    = opj(path_anat ,'native')
            wb_native_dir = opj(dir_native ,'02_Wb')
            volumes_dir   = opj(wb_native_dir ,'volumes')
            labels_dir    = opj(volumes_dir ,'labels')

            label_img = opj(labels_dir ,type_norm + segmentation)
            txtfile = opj(labels_dir ,type_norm + segmentation[:-7] + '.txt')

            data = pd.read_table(txtfile, delim_whitespace=True)
            data2 = data.rename(index=str, columns={"#Magnitude": "label", "Freq": "count"})
            data2.drop(['Cum_Freq'], axis="columns")
            dim = np.prod(nib.load(label_img).header.get_zooms())
            data2['volume'] = data2['count'] * dim

            volume = data2['volume'].where(data2['volume' ] >0).dropna()

            # index of each volume
            list_volume_index = volume.index.values.tolist()
            list_volume_index = np.array([int(i) for i in list_volume_index])

            # volume to list
            list_volume = volume.values.tolist()

            print(len(list_volume_index))
            ####rajout

            # extract their index in the atlas and in the file
            list_volume_index_indexesall = []
            label_atlas_index_indexesall = []

            for val in uniqueValues:
                ttt = np.argwhere(list_volume_index == val).flatten()
                list_volume_index_indexesall.extend(ttt.tolist())
                ddd = np.argwhere(label_atlas == val).flatten()
                label_atlas_index_indexesall.extend(ddd.tolist())

            # Convert to numpy arrays of integer type
            list_volume_index_indexesall = np.array(list_volume_index_indexesall, dtype=int)
            label_atlas_index_indexesall = np.array(label_atlas_index_indexesall, dtype=int)

            # Filter the atlas names
            filtered_df = atlas_names_Seg[atlas_names_Seg['label'].isin(label_atlas)].reset_index(drop=True)

            # Keep the names of regions NOT in the atlas (opposite logic?)
            list_name_atlas_final = filtered_df[~filtered_df.label.isin(label_atlas_index_indexesall)]

            # Remove regions from sub-index
            list_volume_final = np.delete(np.array(list_volume), list_volume_index_indexesall)[1:]
            list_index_final = np.delete(np.array(list_volume_index), list_volume_index_indexesall)[1:]

            final_pd_orig = pd.DataFrame({'values': list_volume_final, 'label': list_index_final})
            data_volumes = pd.merge(final_pd_orig, list_name_atlas_final, on='label').transpose()
            data_volumes.columns = data_volumes.loc['region']
            data_volumes = data_volumes.drop(index=['region', 'label'])
            IDsession = pd.DataFrame([[ID, Session]], columns=['subject', 'session'], index=['values'])

            ########concatenate all volume files
            data_volumes2 = data_volumes.join(IDsession)
            data_volumes3 = data_volumes2.set_index(['session', 'subject'])
            all_data_volumes = pd.concat([all_data_volumes, data_volumes3])

        ####################################################################contat that the excel array
        allinfo_study = allinfo_study_c.set_index(['session', 'subject'])
        result = pd.merge(allinfo_study, all_data_volumes, on=['subject', 'session'])

        sns.set(style="whitegrid")
        plt.switch_backend('Agg')

        def plot_longitudinal_results(result, regressor_list, list_name_atlas_final, IDSEG, bids_dir):
            # --- Handle subject column ---
            if isinstance(result.index, pd.MultiIndex) or 'subject' in result.index.names:
                result = result.reset_index()

            if 'subject' not in result.columns:
                raise ValueError("The dataframe must contain a 'subject' column.")

            result['subject'] = result['subject'].astype(str)
            print("Available columns after subject handling:", result.columns.tolist())

            # --- Convert region columns to numeric ---
            roi_list = list_name_atlas_final['region'].astype(str).tolist()
            for col in roi_list:
                if col in result.columns:
                    result[col] = pd.to_numeric(result[col], errors='coerce')

            # --- Loop over regressors ---
            for regressor in regressor_list:
                out_path = opj(bids_dir, 'Results', 'Volumes', IDSEG, regressor)
                os.makedirs(out_path, exist_ok=True)

                result.to_excel(opj(out_path, f'{IDSEG}.xlsx'), index=False)

                # --- Loop over ROIs ---
                for region in roi_list:
                    if region not in result.columns:
                        print(f"Skipping {region} - not found in results")
                        continue

                    try:
                        # Create safe names for the model
                        safe_region = f"region_{region}".replace(' ', '_').replace('-', '_').replace('(', '').replace(
                            ')', '').replace('.', '')
                        safe_regressor = f"reg_{regressor}".replace(' ', '_').replace('-', '_').replace('(',
                                                                                                        '').replace(')',
                                                                                                                    '').replace(
                            '.', '')

                        # Make a working copy with safe column names
                        df = result[['subject', regressor, region]].copy().dropna()
                        df.columns = ['subject', safe_regressor, safe_region]

                        if df.empty:
                            print(f"No data for {region} vs {regressor}")
                            continue

                        df = df.sort_values(['subject', safe_regressor])
                        subjects = df['subject'].unique()

                        plt.figure(figsize=(10, 6))

                        # Plot subject trajectories
                        colors = plt.cm.rainbow(np.linspace(0, 1, len(subjects)))
                        for subj, color in zip(subjects, colors):
                            subj_df = df[df['subject'] == subj]
                            plt.plot(subj_df[safe_regressor], subj_df[safe_region],
                                     'o-', color=color, alpha=0.6, markersize=4, linewidth=1)

                        # Prepare data for LME
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore', ConvergenceWarning)

                            # Standardize the regressor
                            df['_reg_z'] = (df[safe_regressor] - df[safe_regressor].mean()) / df[safe_regressor].std()

                            # Fit the mixed effects model
                            model = smf.mixedlm(f"{safe_region} ~ _reg_z", data=df, groups=df["subject"])
                            result_mixed = model.fit(reml=False, method="lbfgs", maxiter=200)

                            coef = result_mixed.params.get("_reg_z", np.nan)
                            intercept = result_mixed.params.get("Intercept", np.nan)
                            pval = result_mixed.pvalues.get("_reg_z", np.nan)

                            # Create grid for predictions
                            x_grid = np.linspace(df[safe_regressor].min(), df[safe_regressor].max(), 100)
                            x_grid_z = (x_grid - df[safe_regressor].mean()) / df[safe_regressor].std()

                            # Get fixed effects parameters and covariance matrix
                            params = result_mixed.params
                            cov = result_mixed.cov_params()

                            # We only want the fixed effects part (Intercept and _reg_z)
                            fe_params = params[['Intercept', '_reg_z']]
                            fe_cov = cov.loc[['Intercept', '_reg_z'], ['Intercept', '_reg_z']]

                            # Create design matrix for predictions
                            X = np.column_stack([np.ones_like(x_grid_z), x_grid_z])

                            # Calculate predictions (fixed effects only)
                            pred_vals = X @ fe_params

                            # Calculate standard errors for predictions
                            pred_var = np.sum(X @ fe_cov * X, axis=1)
                            pred_se = np.sqrt(pred_var)

                            # Calculate confidence intervals
                            ci_lower_pred = pred_vals - 1.96 * pred_se
                            ci_upper_pred = pred_vals + 1.96 * pred_se

                            # Plot the population trend and CI
                            plt.plot(x_grid, pred_vals, color='black', linewidth=2, label='Population Trend')
                            plt.fill_between(x_grid, ci_lower_pred, ci_upper_pred,
                                             color='gray', alpha=0.3, label='95% CI')

                            stats_text = (f"Slope (z): {coef:.2f}\n"
                                          f"p-value: {pval:.2f}\n"
                                          f"intercept: {intercept:.2f}\n")

                        # Final formatting
                        plt.gca().text(0.05, 0.95, stats_text,
                                       transform=plt.gca().transAxes,
                                       verticalalignment='top',
                                       bbox=dict(facecolor='white', alpha=0.8))
                        plt.xlabel(regressor, fontsize=12)
                        plt.ylabel(f'{region} Volume', fontsize=12)
                        plt.title(f'Longitudinal: {region} vs {regressor}', fontsize=14)
                        plt.grid(True, alpha=0.3)
                        plt.tight_layout()
                        plt.legend(loc='lower right')

                        plt.savefig(opj(out_path, f'{region}_longitudinal.jpg'),
                                    dpi=300, bbox_inches='tight')
                        plt.close()

                    except Exception as e:
                        print(f"Error plotting {region} vs {regressor}: {str(e)}")
                        import traceback
                        traceback.print_exc()
                        plt.close()

        # Usage:
        plot_longitudinal_results(result, regressor_list, list_name_atlas_final, IDSEG, bids_dir)