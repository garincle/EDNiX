import os
import pandas as pd
import subprocess
import glob
from os.path import join as opj


def spco(cmd, shell=False):
    """Execute shell command"""
    result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    return result.stdout if result.returncode == 0 else None


def get_region_names_from_segmentation(sing_wb, segmentation_file, map_number=1):
    """Extract region names from segmentation file"""
    temp_dir = "/tmp/surface_extraction"
    os.makedirs(temp_dir, exist_ok=True)
    temp_file = opj(temp_dir, "label_table.txt")

    cmd = f'{sing_wb}wb_command -cifti-label-export-table "{segmentation_file}" {map_number} "{temp_file}"'
    spco(cmd, shell=True)

    region_names = []
    if os.path.exists(temp_file):
        with open(temp_file, 'r') as f:
            for line in f.read().strip().split('\n'):
                line = line.strip()
                if line and (line.startswith('l_') or line.startswith('r_')):
                    region_names.append(line)
        os.remove(temp_file)

    return region_names


def extract_region_surface_area(sing_wb, segmentation_file, region_name, subject_id, session_dir, map_number=1):
    """Extract surface area for a single region"""
    Folder_ROIs = opj(session_dir, 'ROIs')
    dir_native_resol = opj(session_dir, 'anat', 'native', 'surfaces', 'Native_resol')
    os.makedirs(Folder_ROIs, exist_ok=True)

    # Determine hemisphere from region name
    if region_name.startswith('l_'):
        hemisphere, H_SIDE = 'l', 'LEFT'
    elif region_name.startswith('r_'):
        hemisphere, H_SIDE = 'r', 'RIGHT'
    else:
        return None

    # Check if surface file exists
    surface_file = opj(dir_native_resol, f"{subject_id}.{hemisphere}.midthickness.surf.gii")
    if not os.path.exists(surface_file):
        return None

    # Create ROI
    roi_file = opj(Folder_ROIs, f"{region_name}_rois.dscalar.nii")
    cmd = f'{sing_wb}wb_command -cifti-label-to-roi "{segmentation_file}" "{roi_file}" -map {map_number} -name "{region_name}"'
    if spco(cmd, shell=True) is None:
        return None

    # Separate by hemisphere
    shape_file = opj(Folder_ROIs, f"{region_name}_rois.shape.gii")
    cmd = f'{sing_wb}wb_command -cifti-separate "{roi_file}" COLUMN -metric CORTEX_{H_SIDE} "{shape_file}"'
    if spco(cmd, shell=True) is None:
        return None

    # Calculate surface areas
    area_file = opj(Folder_ROIs, f"{region_name}_midthickness_shape.gii")
    cmd = f'{sing_wb}wb_command -surface-vertex-areas "{surface_file}" "{area_file}"'
    if spco(cmd, shell=True) is None:
        return None

    # Get surface area
    cmd = f'{sing_wb}wb_command -metric-stats "{area_file}" -reduce SUM -roi "{shape_file}"'
    DATA = spco(cmd, shell=True)

    try:
        return float(DATA.strip()) if DATA else None
    except ValueError:
        return None


def get_regions_to_process(sing_wb, subject_dir, subject_id, regions_to_process=None, map_number=1):
    """Get list of regions to process for a subject"""
    dir_native_resol = opj(subject_dir, 'anat', 'native', 'surfaces', 'Native_resol')
    segmentation_files = glob.glob(opj(dir_native_resol, f"{subject_id}.EDNIxCSC*.dlabel.nii"))

    if not segmentation_files:
        print(f"No segmentation files found for {subject_id}")
        return None, None

    segmentation_file = next((seg for seg in segmentation_files if "LR" not in seg), segmentation_files[0])
    all_regions = get_region_names_from_segmentation(sing_wb, segmentation_file, map_number)

    if regions_to_process is None:
        return all_regions, segmentation_file

    regions_to_extract = []
    for region_pattern in regions_to_process:
        if region_pattern.startswith('l_') or region_pattern.startswith('r_'):
            if region_pattern in all_regions:
                regions_to_extract.append(region_pattern)
        else:
            matching_regions = [r for r in all_regions if region_pattern in r]
            regions_to_extract.extend(matching_regions)

    return regions_to_extract, segmentation_file


def extract_regions_surface_areas(sing_wb, segmentation_file, regions_to_extract, subject_id, subject_dir,
                                  map_number=1):
    """Extract surface areas for specified regions"""
    results = {}
    print(f"Processing {len(regions_to_extract)} regions for {subject_id} (map {map_number})")

    for region_name in regions_to_extract:
        surface_area = extract_region_surface_area(sing_wb, segmentation_file, region_name,
                                                   subject_id, subject_dir, map_number)
        if surface_area is not None:
            results[region_name] = surface_area
            print(f"✓ {region_name}: {surface_area:.6f}")
        else:
            print(f"✗ {region_name}: skipped")

    return results


def average_hemispheres(results):
    """Average left and right hemisphere values"""
    if not results:
        return results

    averaged_results = {}
    processed_regions = set()

    for region_name in results:
        if region_name in processed_regions:
            continue

        base_name = region_name[2:]  # Remove l_ or r_ prefix
        opposite_hemisphere = f"r_{base_name}" if region_name.startswith('l_') else f"l_{base_name}"

        if opposite_hemisphere in results:
            avg_value = (results[region_name] + results[opposite_hemisphere]) / 2
            averaged_results[base_name] = avg_value
            processed_regions.update([region_name, opposite_hemisphere])
            print(f"✓ {base_name} (avg): {avg_value:.6f}")
        else:
            averaged_results[region_name] = results[region_name]
            processed_regions.add(region_name)
            print(f"✓ {base_name}: WARNING Only one hemisphere available")

    return averaged_results


def save_surface_results(results, subject_dir, overwrite=True):
    """Save surface area results to Excel file"""
    if not results:
        return None

    output_file = opj(subject_dir, "anat", "native", "surfaces", "Native_resol", "surface.xlsx")
    indiv_pd = pd.DataFrame.from_dict(results, orient='index', columns=['Surface_Area'])
    indiv_pd.index.name = 'Region'

    if os.path.exists(output_file) and not overwrite:
        existing_df = pd.read_excel(output_file, index_col=0)
        for region, value in results.items():
            existing_df.loc[region] = value
        existing_df.to_excel(output_file)
        print(f"Updated results: {output_file}")
    else:
        indiv_pd.to_excel(output_file)
        print(f"Saved results to: {output_file}")

    return output_file


def process_subject_regions(sing_wb, subject_dir, subject_id, regions_to_process=None,
                            should_average_hemispheres=False, map_number=1, overwrite=True):
    """Main function to process regions for a subject"""
    # Step 1: Get regions to process
    regions_to_extract, segmentation_file = get_regions_to_process(
        sing_wb, subject_dir, subject_id, regions_to_process, map_number
    )

    if not regions_to_extract:
        return None

    # Step 2: Extract surface areas
    results = extract_regions_surface_areas(
        sing_wb, segmentation_file, regions_to_extract, subject_id, subject_dir, map_number
    )

    # Step 3: Average hemispheres if requested
    if should_average_hemispheres:
        results = average_hemispheres(results)

    # Step 4: Save results
    save_surface_results(results, subject_dir, overwrite)

    return results


def process_all_subjects(sing_wb, bids_root, regions_to_process=None, should_average_hemispheres=False,
                         output_file_name='all_surface_data', overwrite=False, overwrite_indiv=True, map_number=1):
    """Process all subjects in the BIDS directory"""
    all_subjects_data = {}
    subject_dirs = glob.glob(opj(bids_root, "sub-*"))

    print(f"Processing {len(subject_dirs)} subjects (map {map_number})")

    for subject_dir in subject_dirs:
        subject_id = os.path.basename(subject_dir).replace('sub-', '')

        for session_dir in glob.glob(opj(subject_dir, "ses-*")):
            native_surfaces_dir = opj(session_dir, 'anat', 'native', 'surfaces', 'Native_resol')

            if os.path.exists(native_surfaces_dir):
                print(f"\nProcessing {subject_id} - {os.path.basename(session_dir)}")
                subject_data = process_subject_regions(
                    sing_wb, session_dir, subject_id, regions_to_process,
                    should_average_hemispheres, map_number, overwrite=overwrite_indiv
                )

                if subject_data:
                    session_id = os.path.basename(session_dir)
                    all_subjects_data[f"{subject_id}_{session_id}"] = subject_data

    # Save combined results
    if all_subjects_data:
        combined_df = pd.DataFrame(all_subjects_data).T
        combined_df.index.name = 'Subject_Session'
        combined_output = opj(bids_root, output_file_name + '.xlsx')

        if os.path.exists(combined_output) and not overwrite:
            existing_df = pd.read_excel(combined_output, index_col=0)
            updated_df = pd.concat([existing_df, combined_df], ignore_index=False)
            updated_df.to_excel(combined_output)
        else:
            combined_df.to_excel(combined_output)

    return None