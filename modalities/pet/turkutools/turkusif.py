# import
import os


opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd

def weight(file, diary_name, sing_turku):

    if not opi(file):
        raise FileNotFoundError(f"Input file not found: {file}")

    cmd = 'gzip -d ' + file
    run_cmd.run(cmd, diary_name)

    cmd = (sing_turku + 'imgweigh ' + file.replace('.nii.gz', '.nii') + ' ' + file.replace('.nii.gz', '.sif'))
    run_cmd.run(cmd, diary_name)

    cmd = 'gzip -f ' + file.replace('.nii.gz', '.nii')
    run_cmd.run(cmd, diary_name)