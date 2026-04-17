"""
EDNiX BIDS Tools v2  —  Usage examples
========================================
This file shows how to use every major entry point.
Copy/adapt the relevant block for your analysis.
"""

import os, sys
sys.path.insert(0, '//')
from Plotting.ednix_bids_tools import (
    # Part 0 – label parsing
    parse_label_file, get_atlas_label_path,
    # Part 1 – path extractors
    extract_bold_paths, extract_surface_paths, extract_thickness_paths,
    extract_volume_paths, extract_qc_paths,
    # Part 2 – data processors
    process_surfaces, process_thickness, process_volumes,
    process_volumes_and_save, process_qc,
    extract_and_process_surfaces_wb, extract_and_process_thickness_wb,
    # Part 3 – multi-species collector
    collect_multi_species,
    # Part 4 – export
    export_to_excel, export_to_csv, export_summary_stats,
    # Part 5 – plots
    plot_morphometry_intra_bids, plot_morphometry_inter_species,
    plot_multi_bids_comparison, plot_qc_dashboard,
    # pipeline wrapper
    run_full_pipeline,
)

opj = os.path.join

# ─────────────────────────────────────────────────────────────────────────────
# 0.  Label files
# ─────────────────────────────────────────────────────────────────────────────

# Parse a label file directly
label_df = parse_label_file(
    '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/primates/Mouselemur/label_code/EDNIxCSC_label.txt')
print(label_df.head(10))

# Auto-locate label file from the library
ATLAS_LIB = (
    '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/'
)
label_path_ml = get_atlas_label_path(
    ATLAS_LIB,
    'atlas/mammals/primates/Mouselemur',
    'EDNIxCSC',
)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Single-BIDS, single-species – manual workflow
# ─────────────────────────────────────────────────────────────────────────────

BIDS_ML = '/scratch2/EDNiX/Mouse_lemur/BIDS_Garin'

# Subjects to exclude
REMOVE = []
# --- Compute surfaces live with wb_command ---
SING_WB = 'vglrun singularity run --bind /srv/projects/,/srv/projects/easymribrain,/scratch2/,/scratch/ /home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/connectome_workbench_1.5.0-freesurfer-update.sif'  # or 'singularity exec /path/to/container.sif '


# --- Compute thickness live with wb_command ---
df_thick_live = extract_and_process_thickness_wb(
    SING_WB, BIDS_ML, list_to_remove=REMOVE,
    regions_to_process=False,
    overwrite=True,
)

df_surf_live = extract_and_process_surfaces_wb(
    SING_WB, BIDS_ML, list_to_remove=REMOVE,
    regions_to_process=False,
    should_average_hemispheres=False,
    overwrite=True,
)
# --- Paths ---
bold_paths   = extract_bold_paths(BIDS_ML, list_to_remove=REMOVE)
surf_paths   = extract_surface_paths(BIDS_ML, list_to_remove=REMOVE)
thick_paths  = extract_thickness_paths(BIDS_ML, list_to_remove=REMOVE)
vol_paths    = extract_volume_paths(BIDS_ML, 'EDNIxCSC', list_to_remove=REMOVE)
qc_paths     = extract_qc_paths(BIDS_ML, list_to_remove=REMOVE)

# --- DataFrames ---
df_surf  = process_surfaces(surf_paths)
df_thick = process_thickness(thick_paths)
df_qc    = process_qc(qc_paths)

# Volume processing needs the label file
df_vol   = process_volumes(vol_paths, label_path_ml)

# Save volumes per subject and get the combined DataFrame
df_vol   = process_volumes_and_save(
    BIDS_ML, 'EDNIxCSC', label_path_ml,
    list_to_remove=REMOVE,
    overwrite=False)

'''
species_bids_dict = {'Rat': '/scratch2/EDNiX/Rat/BIDS_Gd',
 'Mouse': '/scratch2/EDNiX/Mouse/BIDS_Gd',
 'Dog': '/scratch2/EDNiX/Dog/BIDS_k9',
 'Marmoset': '/scratch2/EDNiX/Marmoset/BIDS_NIH_MBM',
 'Macaque': '/scratch2/EDNiX/Macaque/BIDS_BenHamed',
 'Macaque': '/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin',
 'Macaque': '/scratch2/EDNiX/Macaque/imagina',
 'Human': '/scratch2/EDNiX/Human/BIDS_ds004513-raw-data',
 'Human': '/scratch2/EDNiX/Human/CERMEP_MXFDG',
 'Human': '/scratch2/EDNiX/Human/ds004856',
 'Bat': '/scratch2/EDNiX/Bat/BIDS_bat',
 'Mouse_lemur': '/scratch2/EDNiX/Mouse_lemur/BIDS_Garin'}


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Multi-BIDS, multi-species – recommended workflow
# ─────────────────────────────────────────────────────────────────────────────

species_config = {
    'Macaque': {
        'bids_dirs': [
            '/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin',
            '/scratch2/EDNiX/Macaque/BIDS_Saleem',
        ],
        'list_to_keep':   [],
        'list_to_remove': [('BadSub', '1')],
    },
    'Mouse_lemur': {
        'bids_dirs': [
            '/scratch2/EDNiX/Mouse_lemur/BIDS_Garin',
        ],
        'list_to_remove': [('285AB', '01')],
    },
    'Dog': {
        'bids_dirs': [
            '/scratch2/EDNiX/Dog/Auburn_bids',
            '/scratch2/EDNiX/Dog/BIDS_Dora',
            '/scratch2/EDNiX/Dog/BIDS_k9',
        ],
        'list_to_remove': [('19', '1'), ('30', '1')],
    },
}

# Provide label paths per species (one approach)
atlas_label_paths = {
    'Macaque':     get_atlas_label_path(ATLAS_LIB, 'atlas/mammals/primates/Macaque',     'EDNIxCSC'),
    'Mouse_lemur': get_atlas_label_path(ATLAS_LIB, 'atlas/mammals/primates/Mouselemur',  'EDNIxCSC'),
    'Dog':         get_atlas_label_path(ATLAS_LIB, 'atlas/mammals/carnivorans/Dog',       'EDNIxCSC'),
}

# OR let the tool resolve them automatically from the library
species_atlas_fragments = {
    'Macaque':     'atlas/mammals/primates/Macaque',
    'Mouse_lemur': 'atlas/mammals/primates/Mouselemur',
    'Dog':         'atlas/mammals/carnivorans/Dog',
}

data = collect_multi_species(
    species_config,
    regions_of_interest=['Isocortex', 'Frontal_lobe', 'Somatosensory_cortex',
                         'Hippocampal_formation', 'Amygdala'],
    extract=('surface', 'volume', 'thickness', 'qc'),
    atlas_name='EDNIxCSC',
    atlas_label_paths=atlas_label_paths,  # or use atlas_library_root + species_atlas_fragments
)

# Access individual DataFrames
df_surface_all  = data['surface']   # all species combined
df_volume_all   = data['volume']
df_thickness_all= data['thickness']
df_qc_all       = data['qc']

# Export
OUT = '/scratch2/EDNiX/results/multispecies_analysis'
export_to_excel(data, opj(OUT, 'ednix_results.xlsx'))
export_to_csv(data,   opj(OUT, 'csv'))
export_summary_stats(data, opj(OUT, 'ednix_summary_stats.xlsx'),
                     groupby=('species', 'bids_dir'))


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Plots
# ─────────────────────────────────────────────────────────────────────────────

REGIONS = ['Isocortex', 'Frontal_lobe', 'Somatosensory_cortex',
           'Hippocampal_formation', 'Amygdala']
FIG_DIR = opj(OUT, 'figures')

# Intra-BIDS (single species, left hemisphere)
plot_morphometry_intra_bids(
    df_surface_all[df_surface_all['species'] == 'Macaque'],
    metric_col='surface_area_mm2',
    metric_label='Surface area (mm²)',
    regions=REGIONS,
    output_path=opj(FIG_DIR, 'surface_left_Macaque_intra.png'),
    species='Macaque',
    hemisphere='left',
)

# Inter-species comparison
plot_morphometry_inter_species(
    df_surface_all,
    metric_col='surface_area_mm2',
    metric_label='Surface area (mm²)',
    regions=REGIONS,
    output_path=opj(FIG_DIR, 'surface_left_interspecies.png'),
    hemisphere='left',
    show_stats=True,
)

# Multi-BIDS within species (batch comparison)
plot_multi_bids_comparison(
    df_surface_all,
    metric_col='surface_area_mm2',
    metric_label='Surface area (mm²)',
    regions=REGIONS,
    output_path=opj(FIG_DIR, 'surface_left_Dog_multibids.png'),
    species='Dog',
    hemisphere='left',
)

# Volume inter-species
plot_morphometry_inter_species(
    df_volume_all,
    metric_col='volume_mm3',
    metric_label='Volume (mm³)',
    regions=REGIONS,
    output_path=opj(FIG_DIR, 'volume_left_interspecies.png'),
    hemisphere='left',
    show_stats=True,
)

# Thickness inter-species
plot_morphometry_inter_species(
    df_thickness_all,
    metric_col='thickness_mm',
    metric_label='Cortical thickness (mm)',
    regions=['Isocortex', 'Frontal_lobe', 'Somatosensory_cortex'],
    output_path=opj(FIG_DIR, 'thickness_left_interspecies.png'),
    hemisphere='left',
)

# QC dashboard
plot_qc_dashboard(df_qc_all, opj(FIG_DIR, 'QC_dashboard.png'))


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Or run everything in one call
# ─────────────────────────────────────────────────────────────────────────────

results = run_full_pipeline(
    species_config=species_config,
    output_dir=OUT,
    regions_of_interest=['Isocortex', 'Frontal_lobe', 'Somatosensory_cortex',
                         'Hippocampal_formation', 'Amygdala'],
    plot_regions=['Isocortex', 'Frontal_lobe', 'Somatosensory_cortex'],
    hemispheres=('left', 'right'),
    extract=('surface', 'volume', 'thickness', 'qc'),
    atlas_name='EDNIxCSC',
    atlas_label_paths=atlas_label_paths,
)

print("Generated plots:")
for p in results['plots']:
    print(' ', p)
'''