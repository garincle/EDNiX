"""
Group analysis of the [18F]FDG PET of the macaques of BIDS_Tremblay after
Launcher_pet/Launcher_fdg_Macaque_imagina.py (steps 1-10, export to the EDNiX template space).

Dynamic acquisition: step 7 averages the frames of SUV_window and names the results desc-<x>-mean
(desc-bratio-mean, desc-SUV-mean); desc='bratio' finds them (and desc-bratio of static runs, if any).
File names without a ses- entity (sub-X_trc-FDG_...) take the session of their ses-<label> folder.

Measure: brain ratio (bratio = activity / mean activity over the brain mask of the PET space), the measure of
the human CERMEP launcher (Launcher_PET_group_CERMEP.py) and of the cross-species profiles: normalization='none'.

Several sessions per animal are averaged (aggregate='subject': one observation per animal, as the GLM assumes
independent observations); to compare conditions within animals use design='paired' with paired_col='session'.

Design: DESIGN / COVARIATES / GROUP_ASSIGNMENT below. COVARIATES = None (or []): no covariate, and with no group
column to read (one-sample, paired on sessions, two-sample with GROUP_ASSIGNMENT) participants.tsv is not read.
With fewer than MIN_SUBJECTS animals the script stops after listing the images found.
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
DESIGN = 'regression'                               # 'one-sample', 'two-sample', 'regression', 'paired'
COVARIATES = ['age', 'sex']                         # columns of participants.tsv; None or []: no covariate
TEST_COVARIATE = 'age'                              # regression only
GROUP_ASSIGNMENT = None                             # two-sample without participants.tsv: {'A': [IDs], 'B': [IDs]}
MIN_SUBJECTS = 8                                    # permutation p-values need enough relabellings
SUFFIX = 'age'

# ---------------------------------------------------------------------------------------- 0 inputs
images = pet_group.collect_images(bids_dir, template, DESC)
print(images[['subject', 'session', 'desc']].to_string(index=False))
n_animals = images.subject.nunique()
if n_animals < MIN_SUBJECTS:
    raise SystemExit(f'{n_animals} animal(s) with desc-{DESC} images: at least {MIN_SUBJECTS} are needed for a group '
                     f'analysis (run steps 1-10 of the PET pipeline on the other animals)')
model = dict(covariates=COVARIATES, group_assignment=GROUP_ASSIGNMENT,
             test_covariate=TEST_COVARIATE if DESIGN == 'regression' else None)

label_tables = {a: parse_label_file(get_atlas_label_path(ATLAS_LIB, species_fragment, a))
                for a in ('EDNIxCSCLR', 'EDNIxCSC')}

# ---------------------------------------------------------------------------------------- 1 voxel-wise
res = pet_group.run_group_analysis(
    bids_dir, template, DESIGN, desc=DESC, **model,
    normalization='none', label_tables=label_tables, aggregate='subject',
    # fwhm=None: 2 x the voxel of the export grid (export_resolution='pet'); 8 mm, the human value, would cover a
    # macaque gyrus and its neighbours (brain ~15 times smaller in volume)
    tissue_mask=opj(atlas_ednix, 'volumes', 'masks', species + '_desc-Whole_Brain_Gray_mask.nii.gz'),
    rel_threshold=0.8, fwhm=None, voxel_p=0.001, alpha=0.05, cluster_method='permutation', n_perm=5000,
    background=opj(atlas_ednix, 'volumes', species + '_space-acpc_desc-SS_T1w.nii.gz'),
    peak_atlas='EDNIxCSCLR', output_suffix=SUFFIX)

# ---------------------------------------------------------------------------------------- 2 surfaces
sing_wb = Load_EDNiX_requirement.load_requirement(MAIN_PATH, ATLAS_LIB, bids_dir, 'yes')[4]
pet_surface.render_pet_surfaces(
    res, species, ATLAS_LIB, species_fragment, sing_wb,
    scene_template=opj(atlas_ednix, 'surfaces', 'Native_resol', species + '_thickness.scene'),
    scene_label='New Scene 1', projection='ribbon', surface_type='inflated')

# ---------------------------------------------------------------------------------------- 3 regions
pet_roi_group.run_roi_analysis(
    bids_dir, 'EDNIxCSCLR', 0, DESC, DESIGN, **model,
    normalization='none', min_voxels=10, output_suffix=SUFFIX)
