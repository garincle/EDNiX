import os
import glob
import shutil
from os.path import join as opj
import subprocess

def spco(cmd, shell=False):
    """Execute shell command"""
    result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    return result.stdout if result.returncode == 0 else None

def plot_statistical_map_on_surfaces(sing_wb, statistical_map, subject_dir, subject_id,
                                     output_folder, template_scene, palette='-palette-JET256-pos-max-30'):
    """
    Plot statistical map on surfaces for a single subject

    Parameters:
    sing_wb: singularity workbench command prefix
    statistical_map: path to statistical nifti file
    subject_dir: subject directory (e.g., /path/to/BIDS/sub-XXX/ses-XX)
    subject_id: subject identifier
    output_folder: where to save the results
    template_scene: path to template scene file
    palette: color palette for workbench
    """

    # Get surface directory
    atlas_surf = opj(subject_dir, 'anat', 'native', 'surfaces', 'Native_resol')

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Base name for output files
    base_name = os.path.basename(statistical_map).replace('.nii.gz', '').replace('.nii', '')
    output_base = opj(output_folder, f"{subject_id}_{base_name}")

    Hmin = ['l', 'r']

    # Process each hemisphere
    for h in Hmin:
        # Map volume to surface
        surface_file = opj(output_folder, f"{subject_id}_{h}_{base_name}.func.gii")
        command = f'{sing_wb}wb_command -volume-to-surface-mapping "{statistical_map}" ' \
                  f'"{opj(atlas_surf, subject_id)}.{h}.midthickness.surf.gii" ' \
                  f'"{surface_file}" -enclosing'
        spco(command, shell=True)

        # Set map names
        command = f'{sing_wb}wb_command -set-map-names "{surface_file}" -map 1 "statistical_map"'
        spco(command, shell=True)

    # Create dense scalar file
    dscalar_file = opj(output_folder, f"{subject_id}_{base_name}.dscalar.nii")
    command = f'{sing_wb}wb_command -cifti-create-dense-scalar "{dscalar_file}" ' \
              f'-left-metric "{opj(output_folder, f"{subject_id}_r_{base_name}.func.gii")}" ' \
              f'-right-metric "{opj(output_folder, f"{subject_id}_l_{base_name}.func.gii")}"'
    spco(command, shell=True)

    # Set map names for dscalar
    command = f'{sing_wb}wb_command -set-map-names "{dscalar_file}" -map 1 "statistical_map"'
    spco(command, shell=True)

    # Apply palette
    command = f'{sing_wb}wb_command -cifti-palette "{dscalar_file}" MODE_USER_SCALE "{dscalar_file}" {palette}'
    spco(command, shell=True)

    # Create and modify scene file
    with open(template_scene, 'r') as f:
        scenetxt = f.read()

    # Replace file names in scene
    scenetxt = scenetxt.replace('Drew_run1_l_coord_norm.func.gii', f"{subject_id}_l_{base_name}.func.gii")
    scenetxt = scenetxt.replace('Drew_run1_r_coord_norm.func.gii', f"{subject_id}_r_{base_name}.func.gii")

    # Save modified scene
    scene_file = opj(output_folder, f"{subject_id}_{base_name}.scene")
    with open(scene_file, 'w') as f:
        f.write(scenetxt)

    # Generate final image
    output_image = opj(output_folder, f"{subject_id}_{base_name}.png")
    command = f'{sing_wb}wb_command -show-scene "{scene_file}" fmri_stat "{output_image}" 500 100 -use-window-size'
    spco(command, shell=True)

    print(f"✓ Surface plot created: {output_image}")
    return output_image


def create_surface_qc_summary(sing_wb, bids_root, output_dir, template_scene, scene_ID_name, scene_name="Exemple1",
                              statistical_map_pattern=None, max_subjects=200):
    """
    Create QC summary of surface plots for all subjects in BIDS directory

    Parameters:
    sing_wb: singularity workbench command prefix
    bids_root: root BIDS directory
    output_dir: where to save QC summary
    template_scene: path to template scene file
    statistical_map_pattern: pattern to find statistical maps (optional)
    max_subjects: maximum number of subjects to process
    """

   # Find all subject directories
    subject_dirs = glob.glob(opj(bids_root, "sub-*", "ses-*"))

    print(f"Found {len(subject_dirs)} subject-sessions")

    all_plots = []
    processed_count = 0

    for subject_dir in subject_dirs[:max_subjects]:
        subject_id = os.path.basename(os.path.dirname(subject_dir)).replace('sub-', '')
        session_id = os.path.basename(subject_dir).replace('ses-', '')

        print(f"Processing {subject_id} - {session_id}")

        # Check if surface files exist
        surf_dir = opj(bids_root, f"sub-{subject_id}", f"ses-{session_id}", 'anat', 'native', 'surfaces', 'Native_resol')
        print(surf_dir)
        left_surf = opj(surf_dir, f"{subject_id}.l.midthickness.surf.gii")
        right_surf = opj(surf_dir, f"{subject_id}.r.midthickness.surf.gii")
        print(left_surf)
        if not os.path.exists(left_surf) or not os.path.exists(right_surf):
            print(f"  ✗ Missing surface files for {subject_id}")
            continue

        # Create subject-specific output directory
        os.makedirs(output_dir, exist_ok=True)
        uniq_ID =  f"sub-{subject_id}_ses-{session_id}"
        # If statistical map pattern provided, find and plot statistical maps
        if statistical_map_pattern:
            stat_maps = glob.glob(opj(subject_dir, statistical_map_pattern))
            for stat_map in stat_maps:
                try:
                    plot_path = plot_statistical_map_on_surfaces(
                        sing_wb, stat_map, subject_dir, subject_id,
                        output_dir, template_scene
                    )
                    all_plots.append(plot_path)
                except Exception as e:
                    print(f"  ✗ Error processing statistical map: {e}")
        else:
            # Just create surface visualization without statistical map
            try:
                # Create scene and plot
                with open(template_scene, 'r') as f:
                    scenetxt = f.read()

                scenetxt = scenetxt.replace(scene_ID_name, subject_id)
                scenetxt = scenetxt.replace(scene_ID_name, subject_id)

                scene_file = opj(surf_dir, f"{uniq_ID}_surface.scene")
                with open(scene_file, 'w') as f:
                    f.write(scenetxt)

                output_image = opj(output_dir, f"{uniq_ID}_surface.png")
                command = f'{sing_wb}wb_command -show-scene "{scene_file}" "{scene_name}" "{output_image}" 500 100 -use-window-size'
                spco(command, shell=True)
                print(command)
                all_plots.append(output_image)
                print(f"  ✓ Surface QC created: {output_image}")

            except Exception as e:
                print(f"  ✗ Error creating surface QC: {e}")

        processed_count += 1

    # Create summary HTML page
    create_qc_html_summary(all_plots, output_dir)

    print(f"\nQC summary complete! Processed {processed_count} subjects")
    print(f"Results saved to: {output_dir}")


def create_qc_html_summary(plot_paths, output_dir):
    """Create HTML summary page with all QC plots"""

    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Surface QC Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            .plot-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 20px; }
            .plot-item { border: 1px solid #ddd; padding: 10px; text-align: center; }
            .plot-item img { max-width: 100%; height: auto; }
            .subject-id { font-weight: bold; margin-bottom: 10px; }
        </style>
    </head>
    <body>
        <h1>Surface Quality Control Summary</h1>
        <div class="plot-grid">
    """

    for plot_path in plot_paths:
        plot_name = os.path.basename(plot_path)
        subject_id = plot_name.split('_')[0]
        html_content += f"""
            <div class="plot-item">
                <div class="subject-id">{subject_id}</div>
                <img src="{plot_name}" alt="{subject_id}">
                <div>{plot_name}</div>
            </div>
        """

    html_content += """
        </div>
    </body>
    </html>
    """

    # Save HTML file
    with open(opj(output_dir, 'qc_summary.html'), 'w') as f:
        f.write(html_content)

    # Copy all plot images to output directory if they're not already there
    for plot_path in plot_paths:
        if not plot_path.startswith(output_dir):
            dest_path = opj(output_dir, os.path.basename(plot_path))
            shutil.copy2(plot_path, dest_path)