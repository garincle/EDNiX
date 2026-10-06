import os
import pandas as pd

opj = os.path.join

from Tools.atlas_labels import parse_label_file, get_atlas_label_path, find_species_path


def setup(MAIN_PATH, species, reference, selected_atlases_matrix, wanted_level,
          selected_atlases, default_labels, template_dir):
    Lut_dir = opj(MAIN_PATH, "Atlases_library")
    species_fragment = find_species_path(Lut_dir, species)

    if reference == 'EDNiX':
        available = ['EDNIxCSC', 'EDNIxCSCLR']
        lvl_nb = 4
    elif selected_atlases_matrix == 'all':
        raise ValueError(f"selected_atlases_matrix='all' is only defined for reference 'EDNiX' "
                         f"(got '{reference}'); pass an explicit [atlas_names, levels] list.")

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
    filtered_level = []

    for atlas_name, lvl in zip(selected_atlases_matrix[0], selected_atlases_matrix[1]):
        if wanted_level != 'all' and lvl != wanted_level:
            continue
        segmentation_name_list.append(label_dfs[atlas_name].copy())
        filtered_selected.append(atlas_name)
        filtered_level.append(lvl)

    selected_atlases_matrix = [filtered_selected, filtered_level]

    # ── selected_atlases ─────────────────────────────────────────────────────
    if selected_atlases == 'default':
        selected_atlases = [['EDNIxCSCLR', 3]]  # Can be expanded: [['EDNIxCSCLR', 3], ['EDNIxCSCLR', 1]]
    # If selected_atlases is a single list like ['EDNIxCSCLR', 3], convert to nested list
    elif isinstance(selected_atlases, list) and len(selected_atlases) > 0 and not isinstance(selected_atlases[0], list):
        selected_atlases = [selected_atlases]

    # ── default_labels ─────────────────────────────────────────────────────
    if default_labels == 'default':
        default_labels = [[162, 1162]]  # Can be expanded: [[162, 1162], [1, 3]]
    # If default_labels is a single list like [162, 1162], convert to nested list
    elif isinstance(default_labels, list) and len(default_labels) > 0 and not isinstance(default_labels[0], list):
        default_labels = [default_labels]

    # ── panda_files — labels résolus depuis le bon atlas ─────────────────────
    # Loop over each atlas configuration and its corresponding labels
    panda_files = []

    # Handle case where selected_atlases and default_labels might have different lengths
    # If default_labels has only one entry, use it for all atlas configs
    if len(default_labels) == 1 and len(selected_atlases) > 1:
        default_labels = default_labels * len(selected_atlases)
    # If selected_atlases has only one entry, use it for all label configs
    elif len(selected_atlases) == 1 and len(default_labels) > 1:
        selected_atlases = selected_atlases * len(default_labels)
    # If they don't match and neither is length 1, raise error
    elif len(selected_atlases) != len(default_labels):
        raise ValueError(
            f"Length mismatch: selected_atlases ({len(selected_atlases)}) and default_labels ({len(default_labels)}) must have same length or one must be length 1")

    for atlas_config, label_list in zip(selected_atlases, default_labels):
        # atlas_config is like ['EDNIxCSCLR', 3]
        panda_atlas = atlas_config[0]
        panda_level = atlas_config[1] if len(atlas_config) > 1 else None

        # Get the label dataframe for this atlas
        panda_df = label_dfs.get(panda_atlas)
        if panda_df is None:
            # Si cet atlas n'est pas dans selected_atlases_matrix, le charger quand même
            lpath = get_atlas_label_path(Lut_dir, species_fragment, panda_atlas,
                                         prefer_statslut=False)
            raw_df = parse_label_file(lpath)
            panda_df = raw_df[['region_name', 'label_id']].rename(
                columns={'region_name': 'region', 'label_id': 'label'}
            )

        # Get the labels for this configuration
        subset = panda_df[panda_df['label'].isin(label_list)].copy()
        if subset.empty:
            raise ValueError(
                f"None of the label IDs {label_list} were found in atlas '{panda_atlas}'"
            )

        # Create DataFrame for this configuration
        panda_files.append(pd.DataFrame({
            'region': subset['region'].tolist(),
            'label': subset['label'].tolist(),
        }))

    return Lut_dir, selected_atlases_matrix, wanted_level, segmentation_name_list, selected_atlases, panda_files