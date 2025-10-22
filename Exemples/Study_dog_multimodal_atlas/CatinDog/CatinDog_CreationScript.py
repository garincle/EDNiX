import subprocess
import ants
import nibabel as nib
import numpy as np
import os
from scipy.ndimage import binary_dilation, generate_binary_structure

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
    atlas_folder = '/srv/projects/easymribrain/data/Atlas/'
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

n_for_ANTS = 'HammingWindowedSinc'
###
atlas_folder_orig = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/0_Atlas_modify/'
BASE_SS = atlas_folder_orig + 'Atlas/Cat/templatecc.nii.gz'
list_atlases_Cat = [
    atlas_folder_orig + 'Atlas/Cat_cortical/atlascconly.nii.gz',
    atlas_folder_orig + 'Atlas/Cat/atlas.nii.gz']
template_for_dog_atlas = atlas_folder_orig + 'Atlas/Dog/templatecc.nii.gz'

dir_out = atlas_folder_orig + 'Atlas/CatinDog/transformation/'
if not os.path.exists(dir_out): os.mkdir(dir_out)

# FreeSurfer commands to correct for big signal difference
command = 'mri_nu_correct.mni --i ' + atlas_folder_orig + 'Atlas/Dog/template.nii.gz --o ' + \
    atlas_folder_orig + 'Atlas/Dog/templatecc.nii.gz --distance 60 --proto-iters 150 --stop 0.000'
run_singularity(sing['fs_sif'], command)

command = 'mri_nu_correct.mni --i ' + atlas_folder_orig + 'Atlas/Cat/template.nii.gz --o ' + \
    atlas_folder_orig + 'Atlas/Cat/templatecc.nii.gz --distance 60 --proto-iters 150 --stop 0.000'
run_singularity(sing['fs_sif'], command)

# run dog into cat space transformation (work better that way)
fixed = ants.image_read(BASE_SS)
moving = ants.image_read(template_for_dog_atlas)
mTx7 = ants.registration(fixed=fixed, moving=moving,
                         type_of_transform='SyNCC',
                         initial_transform=None,
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
                         outprefix=opj(dir_out + 'Dog2Cat_SyN_'))

new_cat = ants.apply_transforms(fixed=moving, moving=fixed, transformlist=mTx7['invtransforms'],
                                interpolator='hammingWindowedSinc')
ants.image_write(new_cat, dir_out + 'template_CatinDog.nii.gz', ri=False)

for atlas in list_atlases_Cat:
    fixed = ants.image_read(template_for_dog_atlas)  # dog brain
    moving = ants.image_read(atlas)  # cat brain
    new_cat = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=mTx7['invtransforms'],
                                    interpolator='nearestNeighbor')
    ants.image_write(new_cat, opj(dir_out, opb(atlas)), ri=False)
    command = f'3dcalc -a {template_for_dog_atlas} -b {opj(dir_out, opb(atlas))} -expr "step(a)*b" -prefix {opj(dir_out, opb(atlas))} -overwrite'
    run_singularity(sing['afni_sif'], command)

##### test with Cat cortical atlas (subcortical regions are not great in cat)
atlas = atlas_folder_orig + 'Atlas/Cat_cortical/atlascconly_correct.nii.gz'
command = f'3dcalc -a {template_for_dog_atlas} -b {opj(dir_out, "atlascconly.nii.gz")} -expr "step(a)*b" -prefix {opj(dir_out, opb(atlas))} -overwrite'
run_singularity(sing['afni_sif'], command)

command = f'3dcalc -a {opj(dir_out, opb(atlas))} -expr "a" -prefix {opj(dir_out, opb(atlas))} -overwrite'
run_singularity(sing['afni_sif'], command)

for ndil, radius in enumerate([1, 1, 1]):  # dilatation incrémentale
    print(f"\n==== DILATATION radius = {radius} voxels ====")

    if ope(os.path.join(dir_out, f'atlas_dilated_r{ndil-1}.nii.gz')):
        print(f'base for dil {ndil-1} is atlas_dilated_r{ndil-1}.nii.gz')
        atlas_nii = nib.load(os.path.join(dir_out, f'atlas_dilated_r{ndil-1}.nii.gz'))
        atlas_data = atlas_nii.get_fdata()
        affine = atlas_nii.affine
        header = atlas_nii.header
        labels = np.unique(atlas_data)
        labels = labels[labels != 0]  # On ignore le fond

    else:
        print(f'base for dil {ndil} is ' + str(opj(dir_out, opb(atlas))))
        atlas_nii = nib.load(opj(dir_out, opb(atlas)))
        atlas_data = atlas_nii.get_fdata()
        affine = atlas_nii.affine
        header = atlas_nii.header
        labels = np.unique(atlas_data)
        labels = labels[labels != 0]  # On ignore le fond

    dilated_data = np.zeros_like(atlas_data)

    for label in labels:
        print(f"  → Label {int(label)}")
        binary_mask = (atlas_data == label)

        # Structure de dilatation 3D connectée
        structure = generate_binary_structure(3, 1)
        dilated_mask = binary_dilation(binary_mask, structure=structure, iterations=radius)

        # Évite recouvrement
        dilated_data[(dilated_mask) & (dilated_data == 0)] = label

    out_path = os.path.join(dir_out, f'atlas_dilated_r{ndil}.nii.gz')
    nib.save(nib.Nifti1Image(dilated_data.astype(np.int16), affine, header), out_path)
    print(f"✔️ Sauvegardé : {out_path}")

atlasdog = atlas_folder_orig + 'Atlas/Dog/atlas.nii.gz'
command = f'3dresample -overwrite -master {atlasdog} -rmode NN -input {opj(dir_out, opb("atlas_dilated_r2.nii.gz"))} -prefix {opj(dir_out, opb(atlas))}'
run_singularity(sing['afni_sif'], command)

##### cortical mask of dog
command = f'3dcalc -a /srv/projects/easymribrain/data/Atlas/13_Atlas_project/New_atlas_Dual/Dog/atlaslvl1.nii.gz -expr "equals(a,8)" -prefix {opj(dir_out, "dog_cortical_mask.nii.gz")} -overwrite'
run_singularity(sing['afni_sif'], command)

command = f'3dresample -overwrite -master {atlasdog} -rmode NN -input {opj(dir_out, opb(atlas))} -prefix {opj(dir_out, "Atlas_final_respl.nii.gz")}'
run_singularity(sing['afni_sif'], command)

##### overlap cat - dog  cortex
command = f'3dcalc -a {opj(dir_out, "Atlas_final_respl.nii.gz")} -b {opj(dir_out, "dog_cortical_mask.nii.gz")} -short -expr "a*b" -prefix {opj(dir_out, "new_cortical_atlas_dog.nii.gz")} -overwrite'
run_singularity(sing['afni_sif'], command)

##### combine dog new cortex (BA) with the rest of the brain of the old atlas
command = f'3dcalc -a {opj(dir_out, "new_cortical_atlas_dog.nii.gz")} -b {atlasdog} -short -expr "ifelse(a,a,step(b)*(b+10000))" -prefix {opj(dir_out, "new_atlas_dog_BA_full.nii.gz")} -overwrite'
run_singularity(sing['afni_sif'], command)



import nibabel as nib
import numpy as np
from scipy.ndimage import binary_fill_holes, gaussian_filter
from skimage.morphology import remove_small_objects, binary_closing, ball


def process_atlas(atlas_path, output_path, min_size=10, fill=True, smooth=True, sigma=1.0, closing_radius=1):
    """
    Process a NIfTI atlas by:
    1. Removing regions smaller than min_size voxels
    2. Filling holes in each 3D region
    3. Smoothing each region's boundaries

    Parameters:
    - atlas_path: Path to input NIfTI file
    - output_path: Path to save processed NIfTI file
    - min_size: Minimum voxel size for regions to keep (default 10)
    - fill: Whether to fill holes in regions (default True)
    - smooth: Whether to smooth region boundaries (default True)
    - sigma: Standard deviation for Gaussian smoothing (default 1.0)
    - closing_radius: Radius for morphological closing to fill small gaps (default 1)
    """

    # Load the atlas
    img = nib.load(atlas_path)
    data = img.get_fdata()
    affine = img.affine
    header = img.header

    # Get all unique region labels (excluding 0 which is typically background)
    unique_labels = np.unique(data)
    unique_labels = unique_labels[unique_labels != 0]

    # Initialize output array
    processed_data = np.zeros_like(data)

    for label_val in unique_labels:
        # Create binary mask for current region
        region_mask = (data == label_val)

        # Remove small islands
        cleaned_mask = remove_small_objects(region_mask, min_size=min_size)

        if not np.any(cleaned_mask):
            continue

        # Fill holes in 3D
        if fill:
            # Create spherical footprint for 3D morphological operations
            footprint = ball(closing_radius)
            # Perform morphological closing to fill small gaps
            cleaned_mask = binary_closing(cleaned_mask, footprint=footprint)
            # Fill remaining holes
            cleaned_mask = binary_fill_holes(cleaned_mask)

        # Smooth the region
        if smooth:
            # Convert to float for smoothing
            region_float = cleaned_mask.astype(float)

            # Apply 3D Gaussian smoothing
            smoothed_region = gaussian_filter(region_float, sigma=sigma, mode='nearest')

            # Threshold to maintain region shape
            smoothed_region = (smoothed_region > 0.5).astype(int)

            # Add to output with original label value
            processed_data[smoothed_region > 0] = label_val
        else:
            processed_data[cleaned_mask] = int(label_val)

    # Save the processed atlas
    processed_img = nib.Nifti1Image(processed_data.astype(int), affine, header)
    nib.save(processed_img, output_path)
    print(f"Processed atlas saved to {output_path}")


# Example usage
process_atlas(opj(dir_out, "new_atlas_dog_BA_full.nii.gz"), opj(dir_out, "atlas_smoothed.nii.gz"), min_size=10, fill=True, smooth=False, sigma=1.0, closing_radius=1)

######## Dog in cat
dir_out2 = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/Atlases_V2/DoginCat/'
moving = ants.image_read(opj(dir_out, "atlas_smoothed.nii.gz"))  # dog brain
fixed = ants.image_read(dir_out2 + 'template.nii.gz')  # cat brain

new_cat = ants.apply_transforms(fixed=fixed, moving=moving, transformlist=[
    dir_out + 'Dog2Cat_SyN_1Warp.nii.gz',
    dir_out + 'Dog2Cat_SyN_0GenericAffine.mat'], interpolator='nearestNeighbor')
ants.image_write(new_cat, opj(dir_out2, opb(atlas)), ri=False)

command = f'3dresample -overwrite -master {BASE_SS} -rmode NN -input {opj(dir_out2, opb(atlas))} -prefix {opj(dir_out2, opb(atlas))}'
run_singularity(sing['afni_sif'], command)

command = f'3dcalc -overwrite -a {opj(dir_out2, opb(atlas))} -b {BASE_SS} -prefix {opj(dir_out2, opb(atlas))} -expr "step(b)*a"'
run_singularity(sing['afni_sif'], command)
