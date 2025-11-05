import os
import pandas as pd

opj = os.path.join

from Tools import Load_subject_with_BIDS,Read_atlas

def setup(MAIN_PATH,species,reference,selected_atlases_matrix,wanted_level,selected_atlases,panda_files, template_dir):

    # set up links with available atlases

    Lut_dir = opj(MAIN_PATH, "Atlases_library", "LUT_files")

    if reference == 'EDNiX':
        available = ['EDNIxCSC', 'EDNIxCSCLR']
        lvl_nb = 4
        atlasfoldername = "atlas"
        atlasconfigname = "atlas_config_V3.json"

    elif reference == 'other':
        print()

    config_file = opj(MAIN_PATH, "Atlases_library", atlasfoldername, atlasconfigname)
    path_vars   = {'atlas_dir': template_dir, 'species': species, 'Lut_dir': Lut_dir}
    config      = Read_atlas.load_config(Load_subject_with_BIDS.linux_path(config_file), path_vars)
    atlas_dfs   = Read_atlas.extract_atlas_definitions(config)

    segmentation_name_list = []
    if selected_atlases_matrix == 'all':
        selected = []
        level    = []
        for i in available:
            for j in range(lvl_nb):
                selected.append(i)
                level.append(j)
        selected_atlases_matrix = [selected,level]
    else:
        selected_atlases_matrix = selected_atlases_matrix
    if wanted_level == 'all':
        for i in range(len(selected_atlases_matrix[0])):
            segmentation_name_list.append(atlas_dfs[selected_atlases_matrix[0][i] + '[' + str(selected_atlases_matrix[1][i]) + ']'])
    else:
        for i in range(len(selected_atlases_matrix[0])):
            segmentation_name_list.append(
                atlas_dfs[selected_atlases_matrix[0][i] + '[' + wanted_level + ']'])

    if selected_atlases == 'default':
        # Using NEW VERSION format (single atlas)
        selected_atlases = ['EDNIxCSCLR',3]
    else :
        selected_atlases = selected_atlases

    if panda_files == 'default':
        panda_files = [pd.DataFrame({'region':['retrosplenial'],'label':[162]})]  # Using NEW VERSION format (single DataFrame)
    else:
        panda_files = panda_files


    return Lut_dir,selected_atlases_matrix,wanted_level,segmentation_name_list,selected_atlases,panda_files