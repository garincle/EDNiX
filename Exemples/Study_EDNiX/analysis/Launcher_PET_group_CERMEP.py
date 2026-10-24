"""
Group analysis of the [18F]FDG PET of CERMEP-IDB-MRXFDG (37 healthy adults, Merida et al. 2021) after
Launcher_pet/Launcher_fdg_CERMEP.py (steps 1-10, export to the EDNiX template space).

 1. voxel-wise: effect of age on FDG SUVR (cerebellum reference), sex as covariate, cluster FWE by
    permutation (Statistics/Group_PET/pet_group.py);
 2. surfaces of the result (Statistics/Group_PET/pet_surface.py, as the SBA figures);
 3. regional: same model on the regional SUV tables of the PET space (Statistics/Group_PET/pet_roi_group.py).

Expected FDG decrease with age in the medial / lateral prefrontal and anterior cingulate cortex and the insula
(Kalpouzos et al. 2009 Neurobiol Aging 30:112; Yoshizawa et al. 2014 NeuroImage Clin 5:123).
Columns 'age' and 'sex' of participants.tsv: rename COVARIATES if the dataset uses other names.
"""
import os
import sys

opj = os.path.join
MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(1, MAIN_PATH)

from Tools import Load_EDNiX_requirement
from Tools.atlas_labels import parse_label_file, get_atlas_label_path
from Statistics.Group_PET import pet_group, pet_surface, pet_roi_group

bids_dir = '/scratch2/EDNiX/Human/BIDS_Merida/'
species = 'Human'
template = 'EDNiX'                                  # `reference` of the PET launcher (pet/templates/EDNiX)
ATLAS_LIB = opj(MAIN_PATH, 'Atlases_library')
species_fragment = 'atlas/mammals/primates/Human'
atlas_ednix = opj(ATLAS_LIB, species_fragment, 'EDNiX')
COVARIATES = ['age', 'sex']

label_tables = {a: parse_label_file(get_atlas_label_path(ATLAS_LIB, species_fragment, a))
                for a in ('EDNIxCSCLR', 'EDNIxCSC')}
cerebellum = {'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': 'default', 'name': 'Cer'}

# ---------------------------------------------------------------------------------------- 1 voxel-wise
res = pet_group.run_group_analysis(
    bids_dir, template, 'regression', desc='SUV', covariates=COVARIATES, test_covariate='age',
    normalization=cerebellum, label_tables=label_tables,
    tissue_mask=opj(atlas_ednix, 'volumes', 'masks', species + '_desc-Whole_Brain_Gray_mask.nii.gz'),
    rel_threshold=0.8, fwhm=8.0, voxel_p=0.001, alpha=0.05, cluster_method='permutation', n_perm=5000,
    background=opj(atlas_ednix, 'volumes', species + '_space-acpc_desc-SS_T1w.nii.gz'),
    peak_atlas='EDNIxCSCLR', output_suffix='age')

# ---------------------------------------------------------------------------------------- 2 surfaces
sing_wb = Load_EDNiX_requirement.load_requirement(MAIN_PATH, ATLAS_LIB, bids_dir, 'yes')[4]
pet_surface.render_pet_surfaces(
    res, species, ATLAS_LIB, species_fragment, sing_wb,
    scene_template=opj(atlas_ednix, 'surfaces', 'Native_resol', species + '_thickness.scene'),
    scene_label='New Scene 1', projection='ribbon', surface_type='inflated')

# ---------------------------------------------------------------------------------------- 3 regions
pet_roi_group.run_roi_analysis(
    bids_dir, 'EDNIxCSCLR', 0, 'SUV', 'regression', covariates=COVARIATES, test_covariate='age',
    normalization=[7, 1007], min_voxels=10, output_suffix='age')
