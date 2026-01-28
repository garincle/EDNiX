from anat.studytemplate import anat_to_common_EMB
import glob
import os
opj = os.path.join

MAIN_PATH = '/home/cgarin/Documents/Snake/'
bids_dir = '/home/cgarin/Documents/Snake/'
reftemplate_path = opj(bids_dir, 'template', 'anat_to_common_EMB')
diary_file = opj(MAIN_PATH, 'Create_template_anat_to_common_EMB.txt')
sing_afni = 'singularity run --bind /home/cgarin/Documents/Snake/,/home/cgarin/PycharmProjects/EDNiX/Atlases_library,/home/cgarin/PycharmProjects/EDNiX/Tool_library/Singularity /home/cgarin/PycharmProjects/EDNiX/Tool_library/Singularity/afni_ub24_latest.sif '
os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
anat_filenames = sorted(glob.glob(bids_dir + 'sub-*/ses-*/anat/*_acq-3D*_T2w.nii.gz'))


if not os.path.exists(reftemplate_path):
    os.makedirs(reftemplate_path)

anat_to_common_EMB.anats_to_common(sing_afni,diary_file, anat_filenames, reftemplate_path,
                    registration_kind='affine',
                    nonlinear_levels=[1, 2, 3],
                    nonlinear_minimal_patches=[75],
                    nonlinear_weight_file=None,
                    convergence=0.005, blur_radius_coarse=1.1,
                    caching=False, verbose=1,
                    unifize_kwargs=None, brain_masking_unifize_kwargs=None)