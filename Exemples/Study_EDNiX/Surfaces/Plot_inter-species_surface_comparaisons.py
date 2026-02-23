import os
from os.path import join as opj
from Plotting import Plot_intespecies_surfaces
from QC_group import Extract_surfaces_from_BIDS
# Usage examples:
from Tools import Load_EDNiX_requirement

# Setup paths
MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
reftemplate_path = opj(os.path.dirname(MAIN_PATH), "Atlases_library")

bids_list = ['/scratch2/EDNiX/Rat/BIDS_Gd/',
              '/scratch2/EDNiX/Dog/BIDS_k9/',
              '/scratch2/EDNiX/Marmoset/BIDS_NIH_MBM/',
              '/scratch2/EDNiX/Human/BIDS_ds004513-raw-data/',
              '/scratch2/EDNiX/Human/ds004856/',
              '/scratch2/EDNiX/Macaque/BIDS_Cdt_Garin/',
              '/scratch2/EDNiX/Macaque/BIDS_BenHamed/',
              '/scratch2/EDNiX/Mouse_lemur/BIDS_Garin/',
              '/scratch2/EDNiX/Mouse/BIDS_Gd/']

def build_bids_dirs(bids_paths):
    bids_dirs = {}
    for p in bids_paths:
        norm = os.path.normpath(p)
        parts = norm.split(os.sep)

        # Expect structure: /scratch2/EDNiX/<SPECIES>/<BIDS_DIR>/
        # Find species as element after "EDNiX"
        if "EDNiX" in parts:
            idx = parts.index("EDNiX")
            if idx + 1 < len(parts):
                species = parts[idx + 1]
                bids_dirs[species] = p  # preserve trailing slash
    return bids_dirs

bids_dirs = build_bids_dirs(bids_list)

print(bids_dirs)

specific_regions_list = [[
        "Isocortex",
        "Periallocortex"],
        ["Frontal_lobe",
        "Parietal_lobe",
        "Occipital_lobe",
        "Temporal_lobe"],
        ["Somatosensory_cortex",
        "Visual_striate_cortex",
        "Visual_pre_and_extra_striate_cortex",
        "Auditory_cortex",
        "Motor_and_premotor_cortex"]]

def flatten(xss):
    return [x for xs in xss for x in xs]
REGIONS = flatten(specific_regions_list)

for bids_dir in bids_list:
    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')

    for specific_regions, map_number, over1, over2,  in zip(specific_regions_list, [1,2,3], [True,False,False], [False,False,False]):
        # OPTION 3: Process with averaged hemispheres
        averaged_data = Extract_surfaces_from_BIDS.process_all_subjects(sing_wb, bids_dir,
                                             regions_to_process=specific_regions,
                                             should_average_hemispheres=True,
                                             output_file_name = 'all_surface_data',
                                             overwrite=over2,
                                             overwrite_indiv=over1,
                                             map_number=map_number)

# Output directory for plots
output_dir = "/scratch2/surface_analysis"
os.makedirs(output_dir, exist_ok=True)

# Collect all data
print("Collecting surface area data from all species...")
print(f"Regions of interest: {REGIONS}")
df = Plot_intespecies_surfaces.collect_surface_data(bids_dirs, REGIONS)

if df.empty:
    print("No data found! Check your file paths and region names.")
else:
    print(f"Collected data for {len(df)} region measurements")
    print(f"Species distribution: {df['species'].value_counts().to_dict()}")
    print(f"Region distribution: {df['region'].value_counts().to_dict()}")

    # Save the collected data
    data_file = opj(output_dir, 'all_surface_data.csv')
    df.to_csv(data_file, index=False)
    print(f"Saved raw data to: {data_file}")

    # Create individual violin plots
    print("\nCreating individual violin plots...")
    Plot_intespecies_surfaces.create_violin_plots(df, REGIONS, output_dir)

    # Create combined violin plot
    #print("\nCreating combined violin plot...")
    #create_combined_violin_plot(df, REGIONS, output_dir)

    # Print summary statistics
    print("\nSummary Statistics:")
    summary = df.groupby(['species', 'region'])['surface_area_mm2'].agg([
        'count', 'mean', 'std', 'min', 'max'
    ]).round(3)
    print(summary)