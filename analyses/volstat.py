#import
import os
import subprocess
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt

##########################################
########### Subject loader################
##########################################
#https://bids-standard.github.io/pybids/reports/index.html
from sammba import io_conversions, registration
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
def extractVol(Pd_allinfo_study, regressor_list, all_ID, all_Session, all_data_path, type_norm, segmentation_name_list, segmentation_ID_list, atlas_names_Seg_list, atlas_ref_list, bids_dir):

    for segmentation, IDSEG, atlas_ref, atlas_names_Seg in zip(atlas_names_Seg_list, segmentation_ID_list, atlas_ref_list, segmentation_name_list):

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

            table = spco(['3dhistog', '-int', label_img])
            with open(txtfile, 'wb') as myfile:
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
                ttt = np.argwhere(list_volume_index==val)
                list_volume_index_indexesall = np.append(list_volume_index_indexesall, ttt)
                ddd = np.argwhere(label_atlas==val)
                label_atlas_index_indexesall = np.append(label_atlas_index_indexesall, ddd)

            # keep the names of the regions included in the atlas
            list_name_atlas_final = atlas_names_Seg[~atlas_names_Seg.label.isin(label_atlas_index_indexesall)]

            ####indiv
            # remvove them from the sub index
            list_volume_final = np.delete(np.array(list_volume), list_volume_index_indexesall)[1:]
            list_index_final = np.delete(np.array(list_volume_index), list_volume_index_indexesall)[1:]
            final_pd_orig = pd.DataFrame({'values': list_volume_final, 'label': list_index_final})
            data_volumes = pd.merge(final_pd_orig, list_name_atlas_final, on='label').transpose()
            data_volumes.columns = data_volumes.loc['region']
            data_volumes = data_volumes.drop(index=['region', 'label'])
            IDsession = pd.DataFrame([[ID, Session]], columns=['ID', 'Session'], index=['values'])

            ########concatenate all volume files
            data_volumes2 = data_volumes.join(IDsession)
            data_volumes3 = data_volumes2.set_index(['Session', 'ID'])
            all_data_volumes = pd.concat([all_data_volumes, data_volumes3])

        ####################################################################contat that the excel array
        allinfo_study = Pd_allinfo_study.set_index(['Session', 'ID'])
        result = pd.merge(allinfo_study, all_data_volumes, on=['ID', 'Session'])

        for regressor in regressor_list:
            out_results = opj(bids_dir, 'Results')
            if not os.path.exists(out_results): os.mkdir(out_results)
            out_results_V = opj(out_results, 'Volumes')
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            out_results_V = out_results_V + '/' + IDSEG
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            out_results_V = out_results_V + '/' + regressor
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            result.to_excel(out_results_V + '/' + IDSEG + '.xlsx')
            result.reset_index(inplace=True)

            if IDSEG == 'Aseg_mac':
                out_results_Aseg_mac = out_results_V
                result['White_Matter'] = result['Cerebral_White_Matter_R'] + result['Cerebral_White_Matter_L']
                result['SC'] = result['Thalamus_Proper_L'] + result['Caudate_L'] + result['Putamen_L'] + result['Pallidum_L'] + \
                               result['Accumbens_area_L'] + result['VentralDC_L'] + result['Thalamus_Proper_R'] + result['Caudate_R'] + \
                               result['Putamen_R'] + result['Pallidum_R'] + result['Accumbens_area_R'] + result['VentralDC_R']
                result['Cortex'] = result['Cerebral_Cortex_R'] + result['Cerebral_Cortex_L']
                result['Hippocampus'] = result['Hippocampus_L'] + result['Hippocampus_R']
                result['Amygdala'] = result['Amygdala_L'] + result['Amygdala_R']
                result['Cerebellum_White_Matter'] = result['Cerebellum_White_Matter_L'] + result['Cerebellum_White_Matter_R']
                result['Cerebellum_Cortex'] = result['Cerebellum_Cortex_L'] + result['Cerebellum_Cortex_R']


            for n, region in enumerate(list(list_name_atlas_final['region'])):
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

            ############################################################################
            ################### whole brain Session 1###################################
            ############################################################################
            file_path = out_results_Aseg_mac + '/Aseg_mac.xlsx'
            result2 = result.copy()
            if os.path.exists(file_path):

                result3 = pd.read_excel(file_path)
                result3['ID'] = result3['ID'].fillna(method='ffill')

                result2['Whole_Brain'] = result3['Cerebral_White_Matter_L'] + result3['Cerebral_Cortex_L'] + result3['Lateral_Ventricle_L'] + result3['Inf_Lat_Vent_L'] + \
                     result3['Cerebellum_White_Matter_L'] + \
                     result3['Cerebellum_Cortex_L'] + result3['Thalamus_Proper_L'] + result3['Caudate_L'] + result3['Putamen_L'] + result3['Pallidum_L'] + result3['3rd_Ventricle'] + \
                     result3['4th_Ventricle'] + result3['Brain_Stem'] + result3['Hippocampus_L'] + result3['Amygdala_L'] + result3['Accumbens_area_L'] + \
                     result3['VentralDC_L'] + result3['choroid_plexus_L'] + result3['Cerebral_White_Matter_R'] + result3['Cerebral_Cortex_R'] + \
                     result3['Lateral_Ventricle_R'] + result3['Inf_Lat_Vent_R'] + result3['Cerebellum_White_Matter_R'] + result3['Cerebellum_Cortex_R'] + \
                     result3['Thalamus_Proper_R'] + \
                     result3['Caudate_R'] + result3['Putamen_R'] + result3['Pallidum_R'] + result3['Hippocampus_R'] + result3['Amygdala_R'] + result3['Accumbens_area_R'] + \
                     result3['VentralDC_R'] + result3['choroid_plexus_R'] + result3['WM_hypointensities'] + result3['Optic_Chiasm']


                out_results_V2 = out_results_V + '/' + 'divided_by_whole_brain'
                if not os.path.exists(out_results_V2): os.mkdir(out_results_V2)

                # Divide the content of each column in the list by the value
                for n, region in enumerate(list(list_name_atlas_final['region'])):
                    result2[str(region)] = result2[str(region)] / result2['Whole_Brain']

                result2.to_excel(out_results_V2 + '/' + IDSEG + '.xlsx')

                for n, region in enumerate(list(list_name_atlas_final['region'])):

                    x = regressor
                    y = str(region)
                    subject = 'ID'

                    data = result2[[x, y, subject]].dropna(axis=0)

                    for key, grp in data.groupby(subject):
                        plt.plot(grp[x], grp[y], 'o-', label=key)
                    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
                    plt.tight_layout()
                    plt.subplots_adjust(bottom=0.13)
                    plt.xlabel(regressor)
                    plt.savefig(out_results_V2 + '/' + region + '.jpg', dpi=300)
                    plt.close('all')


            ############################################################################
            ################### whole brain Session 1###################################
            ############################################################################

            file_path = out_results_Aseg_mac + '/Aseg_mac.xlsx'
            result2 = result.copy()
            result2['ID'] = result2['ID'].fillna(method='ffill')
            if os.path.exists(file_path):
                result3 = pd.read_excel(file_path)
                result3['ID'] = result3['ID'].fillna(method='ffill')

                out_results_V2 = out_results_V + '/' + 'divided_by_whole_brain_session1'
                if not os.path.exists(out_results_V2): os.mkdir(out_results_V2)

                for ID in list(pd.array(all_ID).unique()):
                    Whole_Brain = result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebral_White_Matter_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebral_Cortex_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Lateral_Ventricle_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Inf_Lat_Vent_L'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebellum_White_Matter_L'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebellum_Cortex_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Thalamus_Proper_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Caudate_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Putamen_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Pallidum_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), '3rd_Ventricle'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), '4th_Ventricle'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Brain_Stem'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Hippocampus_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Amygdala_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Accumbens_area_L'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'VentralDC_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'choroid_plexus_L'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebral_White_Matter_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebral_Cortex_R'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Lateral_Ventricle_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Inf_Lat_Vent_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebellum_White_Matter_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cerebellum_Cortex_R'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Thalamus_Proper_R'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Caudate_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Putamen_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Pallidum_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Hippocampus_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Amygdala_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Accumbens_area_R'] + \
                     result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'VentralDC_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'choroid_plexus_R'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'WM_hypointensities'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Optic_Chiasm']

                    # Divide the content of each column in the list by the value
                    for n, region in enumerate(list(list_name_atlas_final['region'])):
                        print(result2.loc[(result2['ID'] == ID), str(region)])

                        result2.loc[(result2['ID'] == ID), str(region)] = result2.loc[(result2['ID'] == ID), str(region)] / float(Whole_Brain)

                        print(result2.loc[(result2['ID'] == ID), str(region)] / float(Whole_Brain))
                        print(float(Whole_Brain))

                result2.to_excel(out_results_V2 + '/' + IDSEG + '.xlsx')

                for n, region in enumerate(list(list_name_atlas_final['region'])):

                    x = regressor
                    y = str(region)
                    subject = 'ID'

                    data = result2[[x, y, subject]].dropna(axis=0)

                    for key, grp in data.groupby(subject):
                        plt.plot(grp[x], grp[y], 'o-', label=key)
                    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
                    plt.tight_layout()
                    plt.subplots_adjust(bottom=0.13)
                    plt.xlabel(regressor)
                    plt.savefig(out_results_V2 + '/' + region + '.jpg', dpi=300)
                    plt.close('all')


            ############################################################################
            ################### whole brain Session 1###################################
            ############################################################################

            out_results_V2 = out_results_V + '/' + 'divided_by_regions_session1'
            if not os.path.exists(out_results_V2): os.mkdir(out_results_V2)
            result2 = result.copy()
            result2['ID'] = result2['ID'].fillna(method='ffill')

            for ID in list(pd.array(all_ID).unique()):
                # Divide the content of each column in the list by the value
                for n, region in enumerate(list(list_name_atlas_final['region'])):
                    print(result2.loc[(result2['ID'] == ID), str(region)])
                    print(float(Whole_Brain))
                    print(float(result2.loc[result2['Session'] == 1 & (result2['ID'] == ID), str(region)]))
                    result2.loc[(result2['ID'] == ID), str(region)] = 100*((result2.loc[(result2['ID'] == ID), str(region)] / float(result2.loc[result2['Session'] == 1 & (result2['ID'] == ID), str(region)])))
                    print(result2.loc[(result2['ID'] == ID), str(region)])

            result2.to_excel(out_results_V2 + '/' + IDSEG + '.xlsx')

            for n, region in enumerate(list(list_name_atlas_final['region'])):

                x = regressor
                y = str(region)
                subject = 'ID'

                data = result2[[x, y, subject]].dropna(axis=0)

                for key, grp in data.groupby(subject):
                    plt.plot(grp[x], grp[y], 'o-', label=key)
                plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
                plt.tight_layout()
                plt.subplots_adjust(bottom=0.13)
                plt.xlabel(regressor)
                plt.savefig(out_results_V2 + '/' + region + '.jpg', dpi=300)
                plt.close('all')