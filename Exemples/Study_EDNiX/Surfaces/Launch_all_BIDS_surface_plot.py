import os
from os.path import join as opj

# Load your requirements first
from Tools import Load_EDNiX_requirement
from Plotting import Plot_BIDS_surface_for_QC


MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
reftemplate_path = opj(os.path.dirname(MAIN_PATH), "Atlases_library")
bids_dir = "/srv/projects/easymribrain/scratch/EDNiX/Rat/BIDS_Gd/"

sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
    MAIN_PATH, reftemplate_path, bids_dir, 'yes')

# Function 2: Create QC summary for all subjects
Plot_BIDS_surface_for_QC.create_surface_qc_summary(
    sing_wb=sing_wb,
    bids_root=bids_dir,
    output_dir=bids_dir + "/QC/Surface",
    template_scene=bids_dir + "/sub-301105/ses-1/anat/native/surfaces/Native_resol/Exemple1.scene",
    scene_ID_name="301105",
    scene_name="Exemple1")

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
reftemplate_path = opj(os.path.dirname(MAIN_PATH), "Atlases_library")
bids_dir = "/srv/projects/easymribrain/scratch/EDNiX/Mouse_lemur/BIDS_Garin/"

sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
    MAIN_PATH, reftemplate_path, bids_dir, 'yes')

# Function 2: Create QC summary for all subjects
Plot_BIDS_surface_for_QC.create_surface_qc_summary(
    sing_wb=sing_wb,
    bids_root=bids_dir,
    output_dir=bids_dir + "/QC/Surface",
    template_scene=bids_dir + "/sub-967HACA/ses-01/anat/native/surfaces/Native_resol/Exemple1.scene",
    scene_ID_name="967HACA",
    scene_name="Exemple1")

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
reftemplate_path = opj(os.path.dirname(MAIN_PATH), "Atlases_library")
bids_dir = '/srv/projects/easymribrain/scratch/EDNiX/Mouse/BIDS_Gd'

sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
    MAIN_PATH, reftemplate_path, bids_dir, 'yes')

# Function 2: Create QC summary for all subjects
Plot_BIDS_surface_for_QC.create_surface_qc_summary(
    sing_wb=sing_wb,
    bids_root=bids_dir,
    output_dir=bids_dir + "/QC/Surface",
    template_scene=bids_dir + "/sub-jgrAesMEDISOc22r/ses-2/anat/native/surfaces/Native_resol/Exemple1.scene",
    scene_ID_name="jgrAesMEDISOc22r",
    scene_name="Exemple1")