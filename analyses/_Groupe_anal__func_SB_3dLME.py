import matplotlib.pyplot as plt

MAIN_PATH   = opj('/','srv','projects','easymribrain')
### singularity set up
s_bind   = ' --bind ' + opj('/','scratch','in_Process/') + ',' + MAIN_PATH
afni_sif = ' ' + opj(MAIN_PATH,'code','singularity','afni_make_build_AFNI_23.1.10.sif') + ' '

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

#################################################################################################
####Seed base analysis
#################################################################################################

def LME_EDNiX(BASE_SS, oversample_map, mask_func, folder_atlases, cut_coordsX, cut_coordsY,
              cut_coordsZ, bids_dir, images_dir, min_size):

    ## lauch 3dLME
    command = ('singularity run' + s_bind + afni_sif + '3dLME' + ' -prefix ' + file_results + '3dLME_glt.nii.gz' + ' -jobs' + ' 20' + \
               ' -mask ' + mask_func + ' -model' + ' "Sexe*' + regressor + '"' + ' -qVars' + \
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
                                colorbar=True, bg_img=template_filename,
                                display_mode='mosaic',
                                cut_coords=(len(cut_coordsY), len(cut_coordsY), len(cut_coordsY)),
                                cmap='jet')
        display.savefig(file_results + 'Seed' + str(gltlabel) + '_3dLME_cluster.jpg')
        display.close()
        plt.close('all')