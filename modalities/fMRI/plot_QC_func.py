import os
from nilearn import plotting
from os.path import join as opj
import matplotlib.pyplot as plt

def plot_qc(brain_image, mask_image, output_path):
    """
    Generate QC plots for the provided fMRI data and mask.
    """
    # Ensure directories exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    try:
        print(f"Attempting to plot brain image: {brain_image}")
        display = plotting.plot_anat(brain_image, threshold='auto', display_mode='mosaic', dim=4)
        print(f"Adding contours from mask: {mask_image}")
        display.add_contours(mask_image, linewidths=0.2, colors=['red'])
        print(f"Saving QC plot with contours to: {output_path}")
        display.savefig(output_path)
        display.close()
    except Exception as e:
        print(f"Failed to plot with contours. Error: {e}")
        try:
            print(f"Plotting brain image without contours: {brain_image}")
            display = plotting.plot_anat(brain_image, threshold='auto', display_mode='mosaic', dim=4)
            display.savefig(output_path)
            display.close()
        except Exception as e_inner:
            print(f"Failed to plot even the brain image. Error: {e_inner}")
    finally:
        plt.close('all')
