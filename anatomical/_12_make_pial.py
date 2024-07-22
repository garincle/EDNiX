####################################################################################
########################## Pial Surface construction ##############################
####################################################################################
import os
import subprocess
import ants

#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
spco = subprocess.check_output
spgo = subprocess.getoutput

def make_pial(FS_dir, animal_folder, type_norm, otheranat, Hmin, Ref_file, do_surfacewith, overwrite,s_bind,fs_sif):

    export_FS = 'export SINGULARITYENV_SUBJECTS_DIR="' + FS_dir + '"'

    IMG = ants.image_read(Ref_file.replace('.nii.gz', '_norm_' + type_norm + '.nii.gz'))
    IMG = ants.smooth_image(IMG, 0.5)
    ants.image_write(IMG,opj(FS_dir,animal_folder,'mri','brain_norm_T1.nii.gz'),ri=False)

    if do_surfacewith == 'T1andT2':

        if 'T1' in type_norm:

            for H in range(2):
                command = export_FS + ';singularity run' + s_bind + fs_sif + 'mris_make_surfaces -white NOWRITE -aseg aseg -orig white -noaparc -mgz -T1 brain -T2 ' + opj(FS_dir,animal_folder, 'mri', otheranat + 'brain') + ' ' + animal_folder + ' ' + Hmin[H] + 'h'
                spco([command], shell=True)
        else:
            for H in range(2):
                command = export_FS + ';singularity run' + s_bind + fs_sif + 'mris_make_surfaces -white NOWRITE -aseg aseg -orig white -noaparc -mgz -T1 ' + 'T1brain' + ' -T2 ' + opj(FS_dir,animal_folder, 'mri', 'brain ') + animal_folder + ' ' + Hmin[H] + 'h'
                spco([command], shell=True)

    elif do_surfacewith == 'T1':
        for H in range(2):
            command = export_FS + ';singularity run' + s_bind + fs_sif + 'mris_make_surfaces -white NOWRITE -aseg aseg -orig white -noaparc -mgz -T1 brain ' + animal_folder + ' ' + Hmin[H] + 'h'
            spco([command], shell=True)
