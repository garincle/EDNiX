#########################################################################################################################################################
#########################################################################################################################################################
##################################################### Myeline 3D GLM ##################################################################################
#########################################################################################################################################################
#########################################################################################################################################################

#####################################################################fix necessary?????????????????


MAIN_PATH   = opj('/','srv','projects','easymribrain')
### singularity set up
s_bind   = ' --bind ' + opj('/','scratch','in_Process/') + ',' + MAIN_PATH
afni_sif = ' ' + opj(MAIN_PATH,'code','singularity','afni_make_build_AFNI_23.1.10.sif') + ' '

mask = '/media/cgarin/Clement_11/1_Macaques/1_PFC_study/8_Study_template/studytemplate2/template2_mask.nii.gz'
reg = 'brain_mask'

#filter for selecting monkeys
filter1 = allinfo_study_b["ID"].isin(["Trinity"]) & allinfo_study_b["Name_BOLD"].isin(["mbep2d_boldHiRes_AP_mbep2d_boldHiRes_AP_31834207.nii"])
filter2 = allinfo_study_b["ID"].isin(["Roshan"]) & allinfo_study_b["Name_BOLD"].isin(["mbep2d_boldHiRes_AP_mbep2d_boldHiRes_AP_59110276.nii"])

allinfo_study = allinfo_study_b[~filter1]
allinfo_study = allinfo_study[~filter2]

#creat animal session comumn
animal_session = []
for ID, DOS in zip(allinfo_study.ID, allinfo_study.DOS):
    timestampStr = DOS.strftime("%Y_%m_%d")
    animal_session.append(ID + timestampStr)


studytemplatebrain = opj('/media/cgarin/Clement_11/1_Macaques/1_PFC_study/8_Study_template/studytemplate2/template2.nii.gz')

list_of_img_to_analyze = []
for ind, DOS in zip(allinfo_study.index, allinfo_study.DOS):
    timestampStr = DOS.strftime("%Y_%m_%d")
    animal_folder = allinfo_study['ID'][ind] + timestampStr
    if animal_folder=='Unity2020_03_14':
        list_of_img_to_analyze.append(opj(output_dir, animal_folder + 'T2fix'  + '/myeline_map.nii.gz'))
    elif animal_folder=='Trinity2020_06_09':
        list_of_img_to_analyze.append(opj(output_dir, animal_folder + 'T2fix'  + '/myeline_map.nii.gz'))
    elif animal_folder=='Trinity2020_12_17':
        list_of_img_to_analyze.append(opj(output_dir, 'Trinity2020_12_18' + '/myeline_map.nii.gz'))
    elif animal_folder=='Roshan2021_09_30':
        print('Roshan2021_09_30 movement on T2')
    else:
        list_of_img_to_analyze.append(opj(output_dir, animal_folder + '/myeline_map.nii.gz'))




template_filename = 
file_results = '/media/cgarin/Clement_11/1_Macaques/1_PFC_study/4_Results_GP/2_Myeline_' + reg + '/'
mask
panda_disign_matrix
list_of_img_to_analyze

for regressor in regressor_list:

    panda_disign_matrix = allinfo_study[['ID', regressor, 'Sexe']]
    panda_disign_matrix.rename(columns={'ID':'Subj'}, inplace=True)
    panda_disign_matrix.rename(columns={'Sexe':'Subj'}, inplace=True)
    panda_disign_matrix['InputFile'] = list_of_img_to_analyze

    filename_disign_matrix = 'disign_matrix.txt'
    disign_matrix_txt = file_results + filename_disign_matrix
    ### remove the disign matrix in case it exists
    if os.path.exists(disign_matrix_txt):
        os.remove(disign_matrix_txt)

    #creat the path
    if not os.path.exists(file_results): os.mkdir(file_results)

    if not os.path.exists(disign_matrix_txt):
        panda_disign_matrix.to_csv(disign_matrix_txt, index=None, sep='\t', mode='a')

    if os.path.exists(file_results + '3dLME_glt.nii.gz'):
        os.remove(file_results + '3dLME_glt.nii.gz')
        os.remove(file_results + 'resid.nii.gz')


    ## lauch 3dLME
    command = ('singularity run' + s_bind + afni_sif + '3dLME' + ' -prefix ' + file_results + '3dLME_glt.nii.gz' + ' -jobs' + ' 20' + \
               ' -mask ' + mask + ' -model' + ' "Sexe*' + regressor + '"' + ' -qVars' + \
               ' "' + regressor + '"' + ' -ranEff' + ' "~1+' + regressor + '"' + ' -num_glt' + ' 4' + \
               ' -gltLabel 1 "' + regressor + '" -gltCode 1 "' + regressor + ' :"' + \
               ' -gltLabel 2 "diffMF" -gltCode 2 "gender : 1*M -1*F ' + regressor + ' :"' + \
               ' -gltLabel 3 "1M" -gltCode 3 "gender : 1*M ' + regressor + ' :"' + \
               ' -gltLabel 4 "1F" -gltCode 4 "gender : 1*F ' + regressor + ' :"' + \
               ' -dataTable' + ' @' + file_results + 'disign_matrix.txt' + ' -resid ' + file_results + 'resid.nii.gz')
    spco(command, shell=True)

    # run cluster analysis  3dClusterize !!! change
    spco(['singularity run' + s_bind + afni_sif + '3dClustSim', '-mask', mask, '-prefix', file_results + 'Clust_'])
    3dLME = file_results + '3dLME_glt.nii.gz'

    for i, gltlabel in zip([5, 7, 9, 11], [regressor, '1MNM', '1M', '1NM']):
        if i ==5:

            spco(['singularity run' + s_bind + afni_sif + '1d_tool.py', '-infile', file_results + 'Clust_.NN1_2sided.1D[2]{5}', '-write', file_results + 'clustertresh.txt', '-overwrite'])
            print(open(file_results + 'clustertresh.txt', "r").read())
            cluster_threshold = float(open(file_results + 'clustertresh.txt', "r").read()[:4])

        output_z = nib.load(3dLME).get_data()[:,:,:,0,i]
        output_z = nilearn.image.new_img_like(3dLME, output_z, copy_header=True)
        output_z.to_filename(file_results + 'Seed' + str(gltlabel) + '_3dLME.nii.gz')

        thresholded_map1, min_n = map_threshold(output_z, alpha=0.001, height_control="fpr", cluster_threshold=cluster_threshold)
        thresholded_map1.to_filename(file_results + str(gltlabel) + '3dLME_cluster_OY.nii.gz')

        display = plotting.plot_stat_map(thresholded_map1, dim=0,
                                colorbar=True, bg_img=template_filename, display_mode='y', cut_coords=[30, 25, 20, 15, 10, 5, 0, -5, -10, -15, -20], cmap='jet')
        display.savefig(file_results + 'Seed' + str(gltlabel) + '_3dLME_cluster_y_.jpg')
        display.close()

        thresholded_map1 = file_results + str(gltlabel) + '3dLME_cluster.nii.gz'

        display = plotting.plot_stat_map(thresholded_map1, dim=0,
                                colorbar=True, bg_img=template_filename, display_mode='z', cut_coords=[-10, -7, -4, -1, 2, 5, 8, 11, 14, 17, 20], cmap='jet')
        display.savefig(file_results + 'Seed' + str(gltlabel) + '_3dLME_cluster_z_.jpg')
        display.close()
