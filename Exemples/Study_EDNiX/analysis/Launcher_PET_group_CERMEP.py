"""
Group analysis of the [18F]FDG PET of CERMEP-IDB-MRXFDG (37 healthy adults, Merida et al. 2021) after
Launcher_pet/Launcher_fdg_CERMEP.py (steps 1-10, export to the EDNiX template space).

 1. voxel-wise: effect of age on the FDG brain ratio (bratio), sex as covariate, cluster FWE by permutation
    (Statistics/Group_PET/pet_group.py);
 2. surfaces of the result (Statistics/Group_PET/pet_surface.py, as the SBA figures);
 3. regional: same model on the regional bratio tables of the PET space (Statistics/Group_PET/pet_roi_group.py).

bratio (step 7 of the PET pipeline) = activity / mean activity over the brain mask of the PET space: the global
proportional scaling of SPM, already applied, so normalization='none' here (a second division by the brain mean
of the analysis mask would only rescale every image by a near-constant factor). Same measure as the cross-species
profiles (PET_interspecies_profiles.py). Caveat: a global decrease with age is removed by the ratio and a
widespread decrease can appear as relative increases in the regions that are spared (Borghammer et al. 2008
NeuroImage 39:1608; Yakushev et al. 2009 NeuroImage 44:43). Check: the same model on desc='SUV' with
normalization={'type': 'atlas', 'atlas': 'EDNIxCSCLR', 'level': 0, 'labels': 'default', 'name': 'Cer'}
(cerebellum reference); the brain-mean SUV of every subject is in inputs.tsv of a run with normalization='global'.

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
DESC = 'bratio'                                     # static run: desc-bratio (desc-bratio-mean for a dynamic run)
TISSUE_MASK = opj(atlas_ednix, 'volumes', 'masks', species + '_desc-Whole_Brain_Gray_mask.nii.gz')
BACKGROUND = opj(atlas_ednix, 'volumes', species + '_space-acpc_desc-SS_T1w.nii.gz')
SCENE = opj(atlas_ednix, 'surfaces', 'Native_resol', species + '_thickness.scene')
FWHM= 4
# ---------------------------------------------------------------------------------------- 1 group mean
mean = pet_group.run_group_mean(
    bids_dir, template, DESC, tissue_mask=TISSUE_MASK, rel_threshold=0, fwhm=FWHM, background=BACKGROUND)
sing_wb = Load_EDNiX_requirement.load_requirement(MAIN_PATH, ATLAS_LIB, bids_dir, 'yes')[4]

# ---------------------------------------------------------------------------------------- 2 surfaces of the mean
pet_surface.render_pet_surfaces(
    mean, species, ATLAS_LIB, species_fragment, sing_wb, scene_template=SCENE, scene_label='New Scene 1',
    maps=['mean_all', 'cv_all'], projection='ribbon', surface_type='inflated')

# ---------------------------------------------------------------------------------------- 1 voxel-wise
res = pet_group.run_group_analysis(
    bids_dir, template, 'regression', desc=DESC, covariates=COVARIATES, test_covariate='age',
    normalization='none', label_tables=label_tables,
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
    bids_dir, 'EDNIxCSCLR', 0, DESC, 'regression', covariates=COVARIATES, test_covariate='age',
    normalization='none', min_voxels=10, output_suffix='age')
