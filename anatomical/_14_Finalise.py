# import
import os
import subprocess

opj = os.path.join
spco = subprocess.check_output
spgo = subprocess.getoutput
opb = os.path.basename
surface = ['Left', 'Right']
Hmin = ['l', 'r']


def FS_finalise(FS_dir, animal_folder, FreeSlabel_ctab_list, list_atlases_2, labels_dir, type_norm,
                Ref_file,s_bind,fs_sif):

    export_FS = 'export SINGULARITYENV_SUBJECTS_DIR="' + FS_dir + '"'

    cmd = 'singularity run' + s_bind + fs_sif + 'mri_info --orientation ' + Ref_file
    orient_raw = spgo(cmd)

    # purpouse for FS : LIA
    if orient_raw   == 'RAS': reorient = ' -r -1 3 -2 '
    elif orient_raw == 'LAS': reorient = ' -r 1 3 -2 '
    elif orient_raw == 'RSA': reorient = ' -r -1 -2 3'
    elif orient_raw == 'LSA': reorient = ' -r 1 -2 3 '
    elif orient_raw == 'RIA': reorient = ' -r -1 2 3 '
    elif orient_raw == 'LIA': reorient = ' -r 1 2 3 '
    elif orient_raw == 'RSP': reorient = ' -r -1 -2 -3 '
    elif orient_raw == 'LSP': reorient = ' -r 1 -2 -3 '
    elif orient_raw == 'RPS': reorient = ' -r -1 -3 -2 '
    elif orient_raw == 'LPS': reorient = ' -r 1 -3 -2 '

    fwdFS_cmd = ' --in_orientation ' + orient_raw + reorient

    # 1) Create the ribbon.mgz
    command = export_FS + ';singularity run' + s_bind + fs_sif + 'mris_volmask --aseg_name aseg --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon ' + animal_folder
    spco([command], shell=True)

    for H in range(2):
        # 2) Midthickness surface
        command = 'singularity run' + s_bind + fs_sif + 'mris_expand -thickness ' + opj(FS_dir, animal_folder, 'surf', Hmin[H] + 'h.white') + ' 0.5 ' + opj(FS_dir, animal_folder, 'surf', Hmin[H] + 'h.mid')
        spco([command], shell=True)

        # 3) Add the available atlases
        for atlas, FreeSlabel_ctab in zip(list_atlases_2, FreeSlabel_ctab_list):
            # convert to mgz
            command = 'singularity run' + s_bind + fs_sif + 'mri_convert ' + fwdFS_cmd + ' ' + opj(labels_dir, type_norm + opb(atlas))  + ' ' + opj(FS_dir, animal_folder, 'mri', opb(atlas) + '.mgz')
            spco([command], shell=True)

            # create label surface
            command = export_FS + ';singularity run' + s_bind + fs_sif + 'mris_sample_parc -ct ' + FreeSlabel_ctab + ' -surf mid -projfrac 0.01 ' + animal_folder + ' ' + Hmin[H] + 'h ' + opb(atlas) + '.mgz ' + \
            Hmin[H] + 'h.' + animal_folder + '_' + opb(atlas) + '.annot'
            spco([command], shell=True)

        print(surface[H] + ' surface done!')