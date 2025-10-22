import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

######################################################################### ADD VARIANCE!!!!!!!! ###################################

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
def extractThick(s_bind, wb_sif, Pd_allinfo_study, regressor_list, all_ID, all_Session, all_data_path, segmentation_name_list, segmentation_ID_list, atlas_names_Seg_list, atlas_ref_list, bids_dir, overwrite_option):

    for segmentation, IDSEG, atlas_ref, atlas_names_Seg in zip(atlas_names_Seg_list, segmentation_ID_list, atlas_ref_list, segmentation_name_list):
        all_data_volumes = pd.DataFrame()
        for ID, Session, data_path in zip(all_ID, all_Session, all_data_path):
            animal_folder = 'sub-' + ID + '_ses-' + str(Session)
            # The anatomy
            path_anat    = opj(data_path ,'anat/')
            dir_transfo  = opj(path_anat ,'matrices')

            dir_native    = opj(path_anat ,'native')

            dir_prepro    = opj(dir_native ,'01_preprocess')
            wb_native_dir = opj(dir_native ,'02_Wb')
            volumes_dir   = opj(wb_native_dir ,'volumes')
            dir_native_resol = opj(wb_native_dir, 'surfaces', 'Native_resol')
            labels_dir    = opj(volumes_dir ,'labels')
            masks_dir     = opj(volumes_dir ,'masks')

            Folder_ROIs = opj(dir_native_resol, IDSEG)

            # creat path
            if ope(Folder_ROIs) == False:
                os.makedirs(Folder_ROIs)

            H_SIDE = ['LEFT','RIGHT']
            h = ['l','r']


            ### creat empty dataframe for everybody

            final_pd_orig = pd.DataFrame([[ID, Session]], columns=['ID', 'Session'])
            final_pd_orig.set_index(['Session', 'ID'])

            atlas_regions = []
            print(atlas_names_Seg)

            indiv_pd = pd.DataFrame([[ID, Session]], columns=['ID', 'Session'])
            final_pd_orig.set_index(['Session', 'ID'])

            if ope(Folder_ROIs + '/thickness.xlsx') == False or ope(Folder_ROIs + '/thickness.xlsx') == True and overwrite_option == True:
                for index, row in atlas_names_Seg.iterrows():
                    key = row['region']
                    for i in range(0, 2):

                        cmd = 'singularity run' + s_bind + wb_sif + 'wb_command -cifti-label-to-roi ' + opj(dir_native_resol, segmentation + '.dlabel.nii') + \
                              ' ' + opj(Folder_ROIs, str(key) + '.' + h[i] +  '_rois.dscalar.nii') + ' -map 1 -name ' + h[i] +  '_' + str(key)
                        spco(cmd, shell=True)

                        print(str(key) + '.' + h[i] +  '_rois.func.gii')
                        cmd = 'singularity run' + s_bind + wb_sif + 'wb_command -cifti-separate ' + opj(Folder_ROIs, str(key) + '.' + h[i] +  '_rois.dscalar.nii') + ' COLUMN -metric CORTEX_' + H_SIDE[i] + ' ' + \
                        opj(Folder_ROIs, str(key) + '.' + h[i] +  '_rois.func.gii')
                        spco(cmd, shell=True)
                        print(str(key) + '.' + h[i] + '_rois.func.gii22')

                        DATA = spco('singularity run' + s_bind + wb_sif + 'wb_command -metric-stats ' + opj(dir_native_resol, animal_folder  + '.' + h[i] + '.thickness.shape.gii') + \
                                    ' -reduce MEAN -roi ' + opj(Folder_ROIs, str(key) + '.' + h[i] + '_rois.func.gii'), shell=True)

                        #DATA = spco('wb_command -metric-vertex-sum ' + opj(dir_native_resol, animal_folder  + '.' + h[i] + '.thickness.shape.gii') + \
                        #           ' -roi ' + opj(Folder_ROIs, str(key) + '.' + h[i] + '_rois.func.gii'), shell=True)

                        # Decode the byte string to a regular string
                        DATA = DATA.decode('utf-8')
                        # Convert the regular string to a float
                        DATA = float(DATA)
                        print(str(DATA))

                        '''
                        Vertex = spco('wb_command -metric-vertex-sum ' + opj(Folder_ROIs, str(key) + '.' + h[i] +  '_rois.func.gii') + \
                                      ' -roi ' + opj(Folder_ROIs, str(key) + '.' + h[i] +  '_rois.func.gii'), shell=True)
                        Vertex = Vertex.decode('utf-8')
                        # Convert the regular string to a float
                        Vertex = float(Vertex)
                        print(str(Vertex))
                        # Values to add to specific columns
                        indiv_pd[str(key) + '_' + str(h[i])] = float(DATA/Vertex)
                        '''

                        indiv_pd[str(key) + '_' + str(h[i])] = float(DATA)
                        #print('thickness of ' + str(DATA/Vertex))
                        print('thickness of ' + str(DATA))

                        atlas_regions.append(str(key) + '_' + str(h[i]))
                        indiv_pd.to_excel(Folder_ROIs + '/thickness.xlsx')

            else:
                indiv_pd = pd.read_excel(Folder_ROIs + '/thickness.xlsx')
                atlas_regions = indiv_pd.columns.tolist()[3:]

            all_data_volumes = pd.concat([all_data_volumes, indiv_pd])

        ####################################################################contat that the excel array
        allinfo_study = Pd_allinfo_study.set_index(['Session', 'ID'])
        result = pd.merge(allinfo_study, all_data_volumes, on=['ID', 'Session'])



        for regressor in regressor_list:
            out_results = opj(bids_dir, 'Results')
            if not os.path.exists(out_results): os.mkdir(out_results)
            out_results_V = opj(out_results, 'Thickness')
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            out_results_V = out_results_V + '/' + IDSEG
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            out_results_V = out_results_V + '/' + regressor
            if not os.path.exists(out_results_V): os.mkdir(out_results_V)
            result.to_excel(out_results_V + '/' + IDSEG + '.xlsx')
            if IDSEG == 'atlaslvl1':
                out_results_Aseg_mac = out_results_V

            for n, region in enumerate(atlas_regions):
                x = regressor
                y = str(region)
                subject = 'ID'

                data = result[[x, y, subject]].dropna(axis=0)

                for key, grp in data.groupby(subject):
                    plt.plot(grp[x], grp[y], 'o-', label=key)
                plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
                plt.tight_layout()
                plt.subplots_adjust(bottom=0.13)
                plt.xlabel(regressor)
                plt.savefig(out_results_V + '/' + str(region) + '.jpg', dpi=300)
                plt.close('all')

            ############################################################################
            ################### whole brain Session 1###################################
            ############################################################################
            file_path = out_results_Aseg_mac + '/atlaslvl1.xlsx'
            result2 = result.copy()
            if os.path.exists(file_path):

                result3 = pd.read_excel(file_path)
                result3['ID'] = result3['ID'].fillna(method='ffill')

                result2['Cortex'] = result3['Cortex_l'] + result3['Cortex_r']

                out_results_V2 = out_results_V + '/' + 'divided_by_whole_Cortex'
                if not os.path.exists(out_results_V2): os.mkdir(out_results_V2)

                # Divide the content of each column in the list by the value
                for n, region in enumerate(atlas_regions):
                    result2[str(region)] = result2[str(region)] / result2['Cortex']

                result2.to_excel(out_results_V2 + '/' + IDSEG + '.xlsx')

                for n, region in enumerate(atlas_regions):

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

            file_path = out_results_Aseg_mac + '/atlaslvl1.xlsx'
            result2 = result.copy()
            result2['ID'] = result2['ID'].fillna(method='ffill')

            if os.path.exists(file_path):
                result3 = pd.read_excel(file_path)
                result3['ID'] = result3['ID'].fillna(method='ffill')

                out_results_V2 = out_results_V + '/' + 'divided_by_whole_Cortex_session1'
                if not os.path.exists(out_results_V2): os.mkdir(out_results_V2)

                for ID in list(pd.array(all_ID).unique()):
                    Whole_Brain = result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cortex_l'] + result3.loc[result3['Session'] == 1 & (result3['ID'] == ID), 'Cortex_r']
                    # Divide the content of each column in the list by the value
                    for n, region in enumerate(atlas_regions):
                        print(result2.loc[(result2['ID'] == ID), str(region)])

                        result2.loc[(result2['ID'] == ID), str(region)] = result2.loc[(result2['ID'] == ID), str(region)] / float(Whole_Brain)

                        print(result2.loc[(result2['ID'] == ID), str(region)] / float(Whole_Brain))
                        print(float(Whole_Brain))

                result2.to_excel(out_results_V2 + '/' + IDSEG + '.xlsx')

                for n, region in enumerate(atlas_regions):

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
                for n, region in enumerate(atlas_regions):
                    print(result2.loc[(result2['ID'] == ID), str(region)])
                    print(float(Whole_Brain))
                    print(float(result2.loc[result2['Session'] == 1 & (result2['ID'] == ID), str(region)]))
                    result2.loc[(result2['ID'] == ID), str(region)] = 100*((result2.loc[(result2['ID'] == ID), str(region)] / float(result2.loc[result2['Session'] == 1 & (result2['ID'] == ID), str(region)]))-1)
                    print(result2.loc[(result2['ID'] == ID), str(region)])
            result2.to_excel(out_results_V2 + '/' + IDSEG + '.xlsx')

            for n, region in enumerate(atlas_regions):

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


'''

            ###stat et j'ai une valeur par / region L and R
    dictionnary = ancienne valeur et nouvelles

    if ref_atlas == 'fascicularis':
        REF_fasci   = opj(path_ATLAS,'fascicularis','Wb','Native')
        surf_left   = opj(REF_fasci,"fsaverage_LR_32","fascicularis.l.pial.32_fs_LR.surf.gii")
        surf_right  = opj(REF_fasci,"fsaverage_LR_32","fascicularis.r.pial.32_fs_LR.surf.gii")
        ref_ref     = opj(path_ATLAS,'fascicularis','Wb','atlas','fsaverage_LR_32')
        curv_left   = opj(ref_ref,'fascicularis.l.curvature.32_fs_LR.shape.gii') #CHARM_1_in_NMT_v2.0_sym.nii.gz.dlabel.nii
        curv_right  = opj(ref_ref,'fascicularis.l.curvature.32_fs_LR.shape.gii')
        atlas_surf  = opj(REF_fasci ,'Native_resol')

    left                  = nb.load(surf_left)
    coords_l, triangles_l = left.agg_data(('pointset', 'triangle'))
    x1 ,y1, z1            = coords_l.T
    i1, j1, k1            = triangles_l.T
    curvi_l               = nb.load(curv_left)
    color_l               = curvi_l.agg_data()

    right                 = nb.load(surf_right)
    coords_r, triangles_r = right.agg_data(('pointset', 'triangle'))
    x2 ,y2, z2            = coords_r.T
    i2, j2, k2            = triangles_r.T
    curvi_r               = nb.load(curv_right)
    color_r               = curvi_r.agg_data()

    # Assign the values to the vertices based on the labels
    for label, value in zip(np.unique(label_data), values_list):
        vertex_values[label_data == label] = value

    Figures = go.Mesh3d(x=x1,
                        y=y1,
                        z=z1,
                        colorbar_title='z',
                        colorscale=[[0, 'black'],
                                    [0.5, 'gray'],
                                    [1, 'white']],
                        intensity=color_l.T,
                        opacity=0.50,
                        i=i1,
                        j=j1,
                        k=k1,
                        name='y',
                        showscale=False)
'''