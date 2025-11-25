import os
opn = os.path.normpath
opj = os.path.join
opb = os.path.basename
opd = os.path.dirname

from Tools import Load_subject_with_BIDS
from atlases import templatefeat

def get(MAIN_PATH,bids_dir,allinfo_study_c, species,list_to_keep,list_to_remove,reference,type_norm,MNIBcorrect_indiv, masking_img, atlas_followers):

    # reference templates
    if reference == '':
        reference = 'EDNiX'

    reftemplate_path = opj(opd(MAIN_PATH), "Atlases_library")
    fs_tools = opj(reftemplate_path, 'freesurfer_refs')
    path_BALSA = opj(reftemplate_path, 'BALSA')

    (FS_refs, template_dir, reference,
     balsa_folder, BALSAname, balsa_brainT1,
     BASE_folder,BASE_atlas_folder, BASE_template, BASE_SS, BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask,
     CSF, GM, WM, Aseg_ref,
     list_atlas, path_label_code) = templatefeat.get(species, reftemplate_path, fs_tools, path_BALSA, reference,
                                                     '', '',type_norm,'anat', atlas_followers)


    all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max = Load_subject_with_BIDS.load_data_bids(
        allinfo_study_c, bids_dir, list_to_keep, list_to_remove)

    '''
    # Load and process config.
    config_file_path = opj(reftemplate_path, "atlas_config_V3.json")
    config = Read_atlas.load_config(Load_subject_with_BIDS.linux_path(config_file_path),
                                    {'atlas_dir': template_dir, 'Lut_dir': path_label_code})

    # BASE_SS         = config["paths"]["BASE_SS"]
    # BASE_mask       = config["paths"]["BASE_mask"]
    # GM_mask         = config["paths"]["GM_mask"]
    # Aseg_ref        = config["paths"]["Aseg_ref"]
    # Aseg_refLR      = config["paths"]["Aseg_refLR"]
    FreeSlabel_ctab = config["lookup_tables"]["FreeSlabel_ctab_list"]

    atlas_dfs = Read_atlas.extract_atlas_definitions(config)
    (lvl1, lvl1LR, lvl2, lvl2LR,
     lvl3, lvl3LR, lvl4, lvl4LR) = (
        atlas_dfs['lvl1'], atlas_dfs['lvl1LR'],
        atlas_dfs['lvl2'], atlas_dfs['lvl2LR'],
        atlas_dfs['lvl3'], atlas_dfs['lvl3LR'],
        atlas_dfs['lvl4'], atlas_dfs['lvl4LR'])
    # Get all atlas file paths
    (lvl1_file, lvl1LR_file, lvl2_file, lvl2LR_file,
     lvl3_file, lvl3LR_file, lvl4_file, lvl4LR_file) = (
        config["atlas_definitions"]["lvl1"]["atlas_file"],
        config["atlas_definitions"]["lvl1LR"]["atlas_file"],
        config["atlas_definitions"]["lvl2"]["atlas_file"],
        config["atlas_definitions"]["lvl2LR"]["atlas_file"],
        config["atlas_definitions"]["lvl3"]["atlas_file"],
        config["atlas_definitions"]["lvl3LR"]["atlas_file"],
        config["atlas_definitions"]["lvl4"]["atlas_file"],
        config["atlas_definitions"]["lvl4LR"]["atlas_file"])

    # Create combined lists
    list_atlases = [lvl1_file, lvl2_file, lvl3_file, lvl4_file,
                    lvl1LR_file, lvl2LR_file, lvl3LR_file, lvl4LR_file]
    #### for 14 (surfaces) ####
    list_atlases_2 = [lvl1_file, lvl2_file, lvl3_file, lvl4_file,
                      lvl1LR_file, lvl2LR_file, lvl3LR_file, lvl4LR_file]

    ### file to standardize space

    FreeSlabel_ctab_list = [FreeSlabel_ctab[0], FreeSlabel_ctab[0], FreeSlabel_ctab[0], FreeSlabel_ctab[0],
                            FreeSlabel_ctab[1], FreeSlabel_ctab[1], FreeSlabel_ctab[1], FreeSlabel_ctab[1]]
    Lut_file = opj(FreeSlabel_ctab[0])
    '''
    # by default values

    if MNIBcorrect_indiv == '':
        MNIBcorrect_indiv = 'N4'  # "N4" or "N3"

    if masking_img == '':
        masking_img = type_norm # "N4" or "N3"

    return (FS_refs, template_dir, reference,
     balsa_folder, BALSAname, balsa_brainT1,
     BASE_atlas_folder, BASE_template, BASE_SS, BASE_mask, BASE_Gmask, BASE_Wmask, BASE_Vmask,
     CSF, GM, WM, Aseg_ref,
     list_atlas, path_label_code,all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
            fs_tools,reftemplate_path,MNIBcorrect_indiv, masking_img)