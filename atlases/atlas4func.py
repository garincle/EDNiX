import os
import pandas as pd

opj = os.path.join

from Tools import Load_subject_with_BIDS,Read_atlas

def setup(MAIN_PATH,species,reference,selected_atlases_matrix,wanted_level,selected_atlases,panda_files):

    # set up links with available atlases

    Lut_dir = opj(MAIN_PATH, "Atlas_library", "LUT_files")

    if reference == 'EDNIx':
        available = ['EDNIxCSC', 'EDNIxCSCLR']
        lvl_nb = 4
        atlasfoldername = "Atlases_V3"
        atlasconfigname = "atlas_config_V3"

    elif reference == 'other':
        print()

    atlas_dir   = opj(MAIN_PATH, "Atlas_library", atlasfoldername, species)
    config_file = opj(MAIN_PATH, "Atlas_library", atlasfoldername, atlasconfigname)
    path_vars   = {'atlas_dir': atlas_dir, 'species': species, 'Lut_dir': Lut_dir}
    config      = Read_atlas.load_config(Load_subject_with_BIDS.linux_path(config_file), path_vars)
    atlas_dfs   = Read_atlas.extract_atlas_definitions(config)

    if selected_atlases_matrix == 'all':
        selected_atlases_matrix = available
    else:
        selected_atlases_matrix = available[selected_atlases_matrix]
    if wanted_level == 'all':
        wanted_level = range(lvl_nb)

    segmentation_name_list = []

    for i in selected_atlases_matrix:
        for j in wanted_level:
            segmentation_name_list.append(atlas_dfs[i + '[' + str(j) + ']'])

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