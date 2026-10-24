"""
FDG-PET regional profiles, human vs macaque: fingerprints, bar plots, statistics.

Inputs: the ROI tables written by step 7 of the PET pipeline (Launcher_pet/Launcher_fdg_CERMEP.py,
Launcher_pet/Launcher_fdg_Macaque_imagina.py) with the EDNiX atlas (shared label table: regions are homologous by
name). Static human images give desc 'bratio', dynamic macaque images 'bratio-mean' (window SUV_window): each
species takes the first desc of DESC that exists. Method and outputs: Plotting/Plot_interspecies_regional_profiles.py.
"""
import os
import sys

opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from Tools import Load_subject_with_BIDS
from Plotting import Plot_interspecies_regional_profiles as PIRP

sources = {
    'Human': Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Human/BIDS_Merida/')),
    'Macaque': Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Macaque/BIDS_Tremblay/')),
}
out_dir = Load_subject_with_BIDS.linux_path(opj('/scratch2/EDNiX/Comparative/PET_FDG_profiles/'))

config = dict(
    ATLAS='EDNIxCSCLR',
    LEVELS=[0, 1, 2, 3],          # every level written by step 7 (atlas_followers [..., [4, 4], ...])
    DESC=('bratio-PVC', 'bratio', 'bratio-mean-PVC', 'bratio-mean'),   # PVC first when step 'pvc' was run
    HEMISPHERES='merge',          # 'separate' for L / R fingerprints
    NORMALISATION=None,           # bratio is already a ratio to the whole-brain mean; 'zscore' = profile shape only
    REF_SPECIES='Human',
    MIN_VOXELS=10,                # in every subject and hemisphere label (rule of the former pipeline)
    FWHM_MM=None,                 # {'Human': ..., 'Macaque': ...} reconstructed FWHM (mm) -> size_in_fwhm column
    MIN_SIZE_FWHM=None,           # e.g. 2 to drop regions smaller than 2 FWHM in any species (needs FWHM_MM)
    TEST='permutation', N_PERM=10000, N_BOOT=5000, ALPHA=0.05,
)

PIRP.run(sources, out_dir, config)
