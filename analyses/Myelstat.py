import sys
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
spco = subprocess.check_output
spgo = subprocess.getoutput

MAIN_PATH   = opj('/','srv','projects','easymribrain')
sys.path.append(opj(MAIN_PATH + 'Code','EasyMRI_brain-master'))

### singularity set up
s_bind   = ' --bind ' + opj('/','scratch','in_Process/') + ',' + MAIN_PATH
afni_sif = ' ' + opj(MAIN_PATH,'code','singularity','afni_make_build_AFNI_23.1.10.sif') + ' '

def extractMyel(otheranat, max_sessionlist, Pd_allinfo_study, regressor_list,
                         all_ID, all_Session, all_data_path, type_norm, segmentation_name_list,
                         segmentation_ID_list, atlas_names_Seg_list, atlas_ref_list, bids_dir):

    ###########################################################################################################################################################
    ############################################################## start the proces ###########################################################################
    ###########################################################################################################################################################
    ######### define other usefull paramater automatically (do no touch)#########



    ####################################################################################
    ########################## Start the pipeline !!!!!!!!!!!!!!!!!!!!!!   #############
    ####################################################################################

    for segmentation, IDSEG, atlas_ref, atlas_names_Seg in zip(atlas_names_Seg_list, segmentation_ID_list, atlas_ref_list, segmentation_name_list):

        subjects_datas = []
        all_data_volumes = pd.DataFrame()
        for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, max_sessionlist):
            # The anatomy
            path_anat    = opj(data_path,'anat/')
            dir_native    = opj(path_anat,'native')
            dir_prepro    = opj(dir_native,'01_preprocess')
            wb_native_dir = opj(dir_native, '02_Wb')
            wb_T1T2_dir = opj(dir_native, '03_T1T2')


            volumes_dir = opj(wb_native_dir, 'volumes')
            masks_dir = opj(volumes_dir, 'masks')
            labels_dir = opj(volumes_dir, 'labels')

            data_raw = opj(wb_T1T2_dir, ID + '_' + type_norm + '_' + otheranat + '_brain_RSPL.nii.gz')
            label_img =  opj(wb_T1T2_dir, segmentation)
            '''
            label_imgrsp = opj(wb_T1T2_dir, 'RSP_' + segmentation)
            spco(['3dresample', '-master', data_raw, '-prefix',
                  label_imgrsp, '-input',
                  label_img, '-overwrite'])
            '''

            value = spco(['export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";singularity run' + s_bind + afni_sif + '3dROIstats', '-mask', label_img, data_raw + '<0..5>'])
            with open(opj(wb_T1T2_dir, ID + '_' + type_norm + '_' + otheranat +  '_' + IDSEG + '_extracted.csv'), 'wb') as myfile:
                myfile.write(value)
            data_pd = pd.read_table(opj(wb_T1T2_dir, ID + '_' + type_norm + '_' + otheranat +  '_' + IDSEG + '_extracted.csv'), delim_whitespace=True)
            data_pd.columns = [col.replace('Mean_', '') for col in data_pd.columns]
            data_pd = data_pd.drop(['File'], axis="columns")
            data_pd = data_pd.drop(['Sub-brick'], axis="columns")

            mapping_dict = dict(zip(atlas_names_Seg['label'], atlas_names_Seg['region']))

            # Convert the column mapping DataFrame to a dictionary
            column_mapping = dict(zip(atlas_names_Seg['label'], atlas_names_Seg['region']))
            string_mapping_dict = {str(key): value for key, value in mapping_dict.items()}

            # Rename columns using the dictionary
            data_pd.rename(columns=string_mapping_dict, inplace=True)
            # Rename the index to 'region'


            data_volumes = data_pd.transpose()
            data_volumes.columns = ['T1/T2']
            data_volumes.index.name = 'region'
            data_volumes = data_volumes.transpose()

            # Add new columns 'Oliver' and 'Session'
            data_volumes['Session'] = Session
            data_volumes['ID'] = ID



            ########concatenate all volume files
            data_volumes3 = data_volumes.set_index(['Session', 'ID'])
            all_data_volumes = pd.concat([all_data_volumes, data_volumes3])



        ####################################################################contat that the excel array
        ####################################################################contat that the excel array
        allinfo_study = Pd_allinfo_study.set_index(['Session', 'ID'])
        result = pd.merge(allinfo_study, all_data_volumes, on=['ID', 'Session'])

        for regressor in regressor_list:
            out_results = opj(bids_dir, 'Results')
            if not os.path.exists(out_results): os.mkdir(out_results)
            out_results_V = opj(out_results, 'T1T2/')
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            out_results_V = out_results_V + '/' + IDSEG
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            out_results_V = out_results_V + '/' + regressor
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            result.to_excel(out_results_V + '/' + IDSEG + '.xlsx')
            result.reset_index(inplace=True)


            if IDSEG == 'Aseg_mac':
                result['White_Matter'] = result['Cerebral_White_Matter_R'] + result['Cerebral_White_Matter_L']

                result['SC'] = result['Thalamus_Proper_L'] + result['Caudate_L'] + result['Putamen_L'] + result['Pallidum_L'] + \
                               result['Accumbens_area_L'] + result['VentralDC_L'] + result['Thalamus_Proper_R'] + result['Caudate_R'] + \
                               result['Putamen_R'] + result['Pallidum_R'] + result['Accumbens_area_R'] + result['VentralDC_R']
                result['Cortex'] = result['Cerebral_Cortex_R'] + result['Cerebral_Cortex_L']
                result['Hippocampus'] = result['Hippocampus_L'] + result['Hippocampus_R']
                result['Amygdala'] = result['Amygdala_L'] + result['Amygdala_R']
                result['Cerebellum_White_Matter'] = result['Cerebellum_White_Matter_L'] + result['Cerebellum_White_Matter_R']
                result['Cerebellum_Cortex'] = result['Cerebellum_Cortex_L'] + result['Cerebellum_Cortex_R']


            for n, region in enumerate(list(atlas_names_Seg['region'])):
                x =regressor
                y = str(region)
                subject = 'ID'

                data = result[[x, y, subject]].dropna(axis=0)

                for key, grp in data.groupby(subject):
                    plt.plot(grp[x], grp[y], 'o-', label=key)
                plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
                plt.tight_layout()
                plt.subplots_adjust(bottom=0.13)
                plt.xlabel(regressor)
                plt.savefig(out_results_V + '/' + region + '.jpg', dpi=300)
                plt.close('all')
