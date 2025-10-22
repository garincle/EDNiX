import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
from os.path import join as opj


def collect_surface_data(bids_dirs_dict, regions_of_interest):
    """
    Collect surface area data from multiple BIDS directories for different species

    Parameters:
    bids_dirs_dict: dict with {species_name: bids_directory_path}
    regions_of_interest: list of region names to extract
    """
    all_data = []

    for species, bids_dir in bids_dirs_dict.items():
        print(f"Processing {species} from {bids_dir}")

        # Find all surface.xlsx files
        surface_files = glob.glob(opj(bids_dir, "sub-*", "ses-*", "anat", "native", "surfaces", "Native_resol", "surface.xlsx"))

        for surface_file in surface_files:
            try:
                # Extract subject and session info from path
                path_parts = surface_file.split('/')
                subject_id = [p for p in path_parts if p.startswith('sub-')][0]
                session_id = [p for p in path_parts if p.startswith('ses-')][0]

                # Read the surface data
                df = pd.read_excel(surface_file, index_col=0)

                # Extract data for regions of interest
                for region in regions_of_interest:
                    # Look for matching regions (handles l_ and r_ prefixes)
                    matching_regions = [col for col in df.index if region in col]

                    for matched_region in matching_regions:
                        surface_area = df.loc[matched_region, 'Surface_Area']

                        all_data.append({
                            'species': species,
                            'subject': f"{subject_id}_{session_id}",
                            'region': matched_region,
                            'surface_area_mm2': surface_area,
                            'hemisphere': 'left' if matched_region.startswith('l_') else 'right'
                        })

            except Exception as e:
                print(f"Error processing {surface_file}: {e}")
                continue

    return pd.DataFrame(all_data)


def create_violin_plots(df, regions_to_plot, output_dir):
    """Create violin plots for specified regions with quartile lines and swarm plots"""

    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("Set2")

    # Create one plot per region
    for region in regions_to_plot:
        # Filter data for this region
        region_data = df[df['region'].str.contains(region)]

        if region_data.empty:
            print(f"No data found for region: {region}")
            continue

        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 6))

        # Violin plot with quartile lines
        sns.violinplot(data=region_data, x='species', y='surface_area_mm2',
                       inner='quartile', palette='Set2', ax=ax)

        # Add swarm plot for individual data points
        sns.swarmplot(data=region_data, x='species', y='surface_area_mm2',
                      color='black', alpha=0.5, size=3, ax=ax)

        # Customize the plot
        ax.set_title(f'Surface Area Distribution: {region}', fontsize=14, fontweight='bold')
        ax.set_ylabel('Surface Area (mm²)', fontsize=12)
        ax.set_xlabel('Species', fontsize=12)

        # Add some statistics to the plot
        species_stats = region_data.groupby('species')['surface_area_mm2'].agg(['mean', 'std', 'count'])
        stats_text = "\n".join([f"{species}: n={row['count']}, mean={row['mean']:.2f}±{row['std']:.2f}mm²"
                                for species, row in species_stats.iterrows()])

        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)

        plt.tight_layout()

        # Save the plot
        output_file = opj(output_dir, f'surface_area_{region.replace(" ", "_").replace("/", "_")}_violin.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved plot: {output_file}")

        #plt.show()
        plt.close()


def create_combined_violin_plot(df, regions_to_plot, output_dir):
    """Create a combined violin plot with both regions"""

    # Filter data for the regions of interest
    plot_data = pd.concat([df[df['region'].str.contains(region)] for region in regions_to_plot])

    if plot_data.empty:
        print("No data found for the specified regions")
        return

    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Plot 1: First region
    region1_data = plot_data[plot_data['region'].str.contains(regions_to_plot[0])]
    if not region1_data.empty:
        # Violin plot with quartile lines
        sns.violinplot(data=region1_data, x='species', y='surface_area_mm2',
                       inner='quartile', palette='Set2', ax=ax1)
        # Add swarm plot
        sns.swarmplot(data=region1_data, x='species', y='surface_area_mm2',
                      color='black', alpha=0.5, size=3, ax=ax1)
        ax1.set_title(regions_to_plot[0], fontsize=14, fontweight='bold')
        ax1.set_ylabel('Surface Area (mm²)', fontsize=12)

    # Plot 2: Second region
    region2_data = plot_data[plot_data['region'].str.contains(regions_to_plot[1])]
    if not region2_data.empty:
        # Violin plot with quartile lines
        sns.violinplot(data=region2_data, x='species', y='surface_area_mm2',
                       inner='quartile', palette='Set2', ax=ax2)
        # Add swarm plot
        sns.swarmplot(data=region2_data, x='species', y='surface_area_mm2',
                      color='black', alpha=0.5, size=3, ax=ax2)
        ax2.set_title(regions_to_plot[1], fontsize=14, fontweight='bold')
        ax2.set_ylabel('Surface Area (mm²)', fontsize=12)

    plt.tight_layout()

    # Save the combined plot
    output_file = opj(output_dir, 'surface_area_combined_violin.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved combined plot: {output_file}")

    #plt.show()
    plt.close()