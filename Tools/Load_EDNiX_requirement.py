import os
opj = os.path.join

def load_requirement(MAIN_PATH, bids_dir, FS_dir):
    s_bind = ' --bind ' + bids_dir + ',' + FS_dir + ',' + MAIN_PATH
    s_path = opj(MAIN_PATH, 'Tool_library', 'Singularity')
    afni_sif = ' ' + opj(s_path, 'afni_ub24_latest.sif') + ' '
    fsl_sif = ' ' + opj(s_path, 'fsl_6.0.5.1-cuda9.1.sif') + ' '
    fs_sif = ' ' + opj(s_path, 'freesurfer_NHP.sif') + ' '
    itk_sif = ' ' + opj(s_path, 'itksnap_5.0.9.sif') + ' '
    wb_sif = ' ' + opj(s_path, 'connectome_workbench_1.5.0-freesurfer-update.sif') + ' '
    strip_sif = ' ' + opj(s_path, 'synthstrip.1.5.sif') + ' '

    return s_path, afni_sif, fsl_sif, fs_sif, itk_sif, wb_sif, strip_sif, s_bind

