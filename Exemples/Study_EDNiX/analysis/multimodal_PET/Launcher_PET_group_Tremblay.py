"""
Group analysis of the [18F]FDG PET of the macaques of BIDS_Tremblay after
Launcher_pet/Launcher_fdg_Macaque_imagina.py (steps 1-10, export to the EDNiX template space).

 1. group mean of the brain ratio (pet_group.run_group_mean: mean, SD, CV maps; any number of animals);
 2. surfaces of the mean (pet_surface.render_pet_surfaces, as the SBA figures; ribbon-constrained projection);
 3. regional means (pet_roi_group.run_roi_mean: regions of >= 10 PET voxels in every session);
 4. statistics (pet_group.run_group_analysis / pet_roi_group.run_roi_analysis) only with RUN_GLM = True and at
    least MIN_SUBJECTS animals: with 2-3 animals a GLM has no inferential value.

Dynamic acquisition: step 7 averages the frames of SUV_window and names the results desc-<x>-mean
(desc-bratio-mean); desc='bratio' finds them. File names without a ses- entity (sub-X_trc-FDG_...) take the
session of their ses-<label> folder. Sessions of an animal are averaged (one value per animal).

Measure: brain ratio (bratio = activity / mean activity over the brain mask of the PET space), the measure of the
human CERMEP launcher and of the cross-species profiles; it is already normalised: normalization='none'.
Compare with the human images only when SUV_window of the macaque PET launcher matches the human acquisition
(CERMEP: one 10-min frame 50 min after injection); SUV_window=None averages all the frames, perfusion included.
"""
import os
import sys

opj = os.path.join
MAIN_PATH = '/home/cgarin/PycharmProjects/EDNiX/'
sys.path.insert(1, MAIN_PATH)

from Tools import Load_EDNiX_requirement
from Tools.atlas_labels import parse_label_file, get_atlas_label_path
from Statistics.Group_PET import pet_group, pet_surface, pet_roi_group

bids_dir = '/scratch2/EDNiX/Macaque/BIDS_Tremblay/'
species = 'Macaque'
template = 'EDNiX'                                  # `reference` of the PET launcher (pet/templates/EDNiX)
ATLAS_LIB = opj(MAIN_PATH, 'Atlases_library')
species_fragment = 'atlas/mammals/primates/Macaque'
atlas_ednix = opj(ATLAS_LIB, species_fragment, 'EDNiX')
DESC = 'bratio'
TISSUE_MASK = opj(atlas_ednix, 'volumes', 'masks', species + '_desc-Whole_Brain_Gray_mask.nii.gz')
BACKGROUND = opj(atlas_ednix, 'volumes', species + '_space-acpc_desc-SS_T1w.nii.gz')
SCENE = opj(atlas_ednix, 'surfaces', 'Native_resol', species + '_thickness.scene')
FWHM = None                                         # None: 2 x the voxel of the export grid (8 mm = human value)

RUN_GLM = False                                     # statistics (section 4)
DESIGN = 'regression'                               # 'one-sample', 'two-sample', 'regression', 'paired'
COVARIATES = ['age', 'sex']                         # columns of participants.tsv; None or []: no covariate
TEST_COVARIATE = 'age'                              # regression only
GROUP_ASSIGNMENT = None                             # two-sample without participants.tsv: {'A': [IDs], 'B': [IDs]}
MIN_SUBJECTS = 8

sing_wb = Load_EDNiX_requirement.load_requirement(MAIN_PATH, ATLAS_LIB, bids_dir, 'yes')[4]

# ---------------------------------------------------------------------------------------- 1 group mean
mean = pet_group.run_group_mean(
    bids_dir, template, DESC, tissue_mask=TISSUE_MASK, rel_threshold=0, fwhm=FWHM, background=BACKGROUND)

# ---------------------------------------------------------------------------------------- 2 surfaces of the mean
pet_surface.render_pet_surfaces(
    mean, species, ATLAS_LIB, species_fragment, sing_wb, scene_template=SCENE, scene_label='New Scene 1',
    maps=['mean_all', 'cv_all'], projection='ribbon', surface_type='inflated')

# ---------------------------------------------------------------------------------------- 3 regional means
pet_roi_group.run_roi_mean(bids_dir, 'EDNIxCSCLR', 0, DESC, normalization='none', min_voxels=10)

# ---------------------------------------------------------------------------------------- 4 statistics
if not RUN_GLM:
    print('RUN_GLM = False: group mean only')
elif mean['n'] < MIN_SUBJECTS:
    print(f'{mean["n"]} animal(s): statistics skipped (at least {MIN_SUBJECTS} needed)')
else:
    model = dict(covariates=COVARIATES, group_assignment=GROUP_ASSIGNMENT,
                 test_covariate=TEST_COVARIATE if DESIGN == 'regression' else None)
    label_tables = {a: parse_label_file(get_atlas_label_path(ATLAS_LIB, species_fragment, a))
                    for a in ('EDNIxCSCLR', 'EDNIxCSC')}
    res = pet_group.run_group_analysis(
        bids_dir, template, DESIGN, desc=DESC, **model, normalization='none', label_tables=label_tables,
        aggregate='subject', tissue_mask=TISSUE_MASK, rel_threshold=0.8, fwhm=FWHM, voxel_p=0.001, alpha=0.05,
        cluster_method='permutation', n_perm=5000, background=BACKGROUND, peak_atlas='EDNIxCSCLR',
        output_suffix=DESIGN)
    pet_surface.render_pet_surfaces(
        res, species, ATLAS_LIB, species_fragment, sing_wb, scene_template=SCENE, scene_label='New Scene 1',
        projection='ribbon', surface_type='inflated')
    pet_roi_group.run_roi_analysis(
        bids_dir, 'EDNIxCSCLR', 0, DESC, DESIGN, **model, normalization='none', min_voxels=10,
        output_suffix=DESIGN)
