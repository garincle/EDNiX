import subprocess
import ants
import nibabel as nib
import numpy as np
import os
from scipy.ndimage import binary_dilation, generate_binary_structure
import glob
# Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput
ops = os.path.splitext

# Define Singularity paths and bindings
def load_singularity():
    MAIN_PATH = r'/srv/projects/easymribrain/code/EDNiX/' # Update this path
    atlas_folder = '/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2'
    s_bind = f' --bind {atlas_folder},{MAIN_PATH}'
    s_path = opj(MAIN_PATH, 'Tool_library', 'Singularity')

    return {
        'afni_sif': opj(s_path, 'afni_ub24_latest.sif'),
        'fs_sif': opj(s_path, 'freesurfer_NHP.sif'),
        'bind': s_bind}

# Load Singularity containers
sing = load_singularity()

# Helper function to run commands in Singularity
def run_singularity(container, command):
    full_cmd = f'singularity run {sing["bind"]} {container} {command}'
    return spco(full_cmd, shell=True)


###
atlas_folder_orig = '/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/Macaque/'
BASE_SS = atlas_folder_orig + 'template.nii.gz'
list_atlases_Cat = glob.glob('/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/Macaque/*.nii.gz')
template_for_dog_atlas = '/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/Baboon/Transfo/Haiko89_Asymmetric.Template_n89.nii.gz'

dir_out = '/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/Baboon/Transfo/'
if not os.path.exists(dir_out): os.mkdir(dir_out)
fixed = ants.image_read(BASE_SS)
moving = ants.image_read(template_for_dog_atlas)

mtx1 = ants.registration(fixed=fixed, moving=moving, type_of_transform='Translation',
                         outprefix=opj(dir_out,'template_to__SyN_final_shift_'))
MEAN_tr = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=mtx1['fwdtransforms'],
                                interpolator='hammingWindowedSinc')
ants.image_write(MEAN_tr, opj(dir_out + 'template_shft_MacinBaboon.nii.gz'), ri=False)

# run dog into cat space transformation (work better that way)
mTx7 = ants.registration(fixed=fixed, moving=moving,
                         type_of_transform='SyNRA',
                         initial_transform=mtx1['fwdtransforms'],
                         grad_step=0.1,
                         flow_sigma=3,
                         total_sigma=0,
                         aff_sampling=32,
                         aff_random_sampling_rate=0.2,
                         syn_sampling=32,
                         aff_iterations=(1000, 500, 250, 100),
                         aff_shrink_factors=(8, 4, 2, 1),
                         aff_smoothing_sigmas=(3, 2, 1, 0),
                         reg_iterations=(1000, 500, 250, 100),
                         reg_smoothing_sigmas=(3, 2, 1, 0),
                         reg_shrink_factors=(8, 4, 2, 1),
                         outprefix=opj(dir_out + 'template_MacinBaboon_SyN_'),
                         verbose=True)

new_cat = ants.apply_transforms(fixed=moving, moving=fixed, transformlist=mTx7['invtransforms'],
                                interpolator='hammingWindowedSinc')
ants.image_write(new_cat, dir_out + 'template_MacinBaboon.nii.gz', ri=False)


for atlas in list_atlases_Cat:
    fixed = ants.image_read(template_for_dog_atlas)  # dog brain
    moving = ants.image_read(atlas)  # cat brain
    new_cat = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=['/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/Baboon/Transfo/template_MacinBaboon_SyN_0GenericAffine.mat', '/srv/projects/easymribrain/code/EDNiX/Atlas_library/Atlases_V2/Baboon/Transfo/template_MacinBaboon_SyN_1InverseWarp.nii.gz'],
                                    whichtoinvert=[True, False], interpolator='nearestNeighbor')
    ants.image_write(new_cat, opj(dir_out, opb(atlas)), ri=False)
    command = f'3dcalc -a {template_for_dog_atlas} -b {opj(dir_out, opb(atlas))} -expr "step(a)*b" -prefix {opj(dir_out, opb(atlas))} -overwrite'
    run_singularity(sing['afni_sif'], command)