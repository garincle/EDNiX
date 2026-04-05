import os
import pandas as pd

opj = os.path.join

from Tools import Load_subject_with_BIDS,Read_atlas
from Plotting.ednix_bids_tools import parse_label_file, get_atlas_label_path, find_species_path

import os
import pandas as pd

opj = os.path.join

from Tools import Load_subject_with_BIDS
from Plotting.ednix_bids_tools import parse_label_file, get_atlas_label_path, find_species_path


def setup(MAIN_PATH, species, reference, selected_atlases_matrix, wanted_level,
          selected_atlases, default_labels, template_dir):

    Lut_dir          = opj(MAIN_PATH, "Atlases_library")
    species_fragment = find_species_path(Lut_dir, species)

    if reference == 'EDNiX':
        available = ['EDNIxCSC', 'EDNIxCSCLR']
        lvl_nb    = 4

    # ── selected_atlases_matrix ───────────────────────────────────────────────
    if selected_atlases_matrix == 'all':
        selected, level = [], []
        for i in available:
            for j in range(lvl_nb):
                selected.append(i)
                level.append(j)
        selected_atlases_matrix = [selected, level]

    # ── Charge un label_df par atlas unique ───────────────────────────────────
    # EDNIxCSC  → 138 labels bilatéraux  (1–182)
    # EDNIxCSCLR → 276 labels L+R        (1–182 gauche, 1001–1182 droite)
    label_dfs = {}
    for atlas_name in set(selected_atlases_matrix[0]):
        lpath = get_atlas_label_path(Lut_dir, species_fragment, atlas_name,
                                     prefer_statslut=False)
        df = parse_label_file(lpath)
        label_dfs[atlas_name] = (
            df[['region_name', 'label_id']]
            .rename(columns={'region_name': 'region', 'label_id': 'label'})
            .copy()
        )

    # ── segmentation_name_list — DataFrame correct par entrée atlas/niveau ────
    # correl_matrix: for ndim, panda_file in enumerate(segmentation_name_list)
    #                atlas = [selected_atlases_matrix[0][ndim],
    #                          selected_atlases_matrix[1][ndim]]
    #                atlas_data[:,:,:, atlas[1]]  ← niveau via sous-brique 4D
    segmentation_name_list = []
    filtered_selected = []
    filtered_level    = []

    for atlas_name, lvl in zip(selected_atlases_matrix[0], selected_atlases_matrix[1]):
        if wanted_level != 'all' and lvl != wanted_level:
            continue
        segmentation_name_list.append(label_dfs[atlas_name].copy())
        filtered_selected.append(atlas_name)
        filtered_level.append(lvl)

    selected_atlases_matrix = [filtered_selected, filtered_level]

    # ── selected_atlases ─────────────────────────────────────────────────────
    if selected_atlases == 'default':
        selected_atlases = ['EDNIxCSCLR', 3]

    # ── panda_files — labels résolus depuis le bon atlas ─────────────────────
    if default_labels == 'default':
        default_labels = [162, 1162]

    # Utilise le label_df de l'atlas de selected_atlases
    panda_atlas = selected_atlases[0] if isinstance(selected_atlases, list) else 'EDNIxCSCLR'
    panda_df    = label_dfs.get(panda_atlas)
    if panda_df is None:
        # Si cet atlas n'est pas dans selected_atlases_matrix, le charger quand même
        lpath    = get_atlas_label_path(Lut_dir, species_fragment, panda_atlas,
                                        prefer_statslut=False)
        raw_df   = parse_label_file(lpath)
        panda_df = raw_df[['region_name', 'label_id']].rename(
            columns={'region_name': 'region', 'label_id': 'label'}
        )

    subset = panda_df[panda_df['label'].isin(default_labels)].copy()
    if subset.empty:
        raise ValueError(
            f"None of the label IDs {default_labels} were found for atlas '{panda_atlas}'"
        )

    panda_files = [pd.DataFrame({
        'region': subset['region'].tolist(),
        'label':  subset['label'].tolist(),
    })]

    return Lut_dir, selected_atlases_matrix, wanted_level, segmentation_name_list, selected_atlases, panda_files