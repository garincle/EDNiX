##########################################################################

import pandas as pd
import glob
import os
import subprocess
import ants
spgo = subprocess.getoutput
spco = subprocess.check_output
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

output_atlas = '/srv/projects/easymribrain/data/Atlas/Atlases_V2/'
input_atlas = '/srv/projects/easymribrain/data/Atlas/Atlases_V2/'
orig_atlases = '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/0_Atlas_modify/Atlas/'
atlas_exel = '/srv/projects/easymribrain/data/Atlas/Classiff/Classif_EDNiX.xlsx'

MAIN_PATH = r'/home/cgarin/PycharmProjects/EDNiX/'  # Update this path
atlas_folder = '/srv/projects/easymribrain/data/Atlas/'
s_bind = f' --bind {atlas_folder},{MAIN_PATH} '
s_path = opj(MAIN_PATH, 'Singularity_library', 'Singularity')
afni_sif = opj(s_path, 'afni_ub24_latest.sif ')
species = [
    'Pig',
    'Mouse_lemur',
    'Macaque',
    'CatinDog',
    'DoginCat']

transformations = [
    '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/WIP/mammals/ungulates/Pig/atlas/ednix/modified/matrices/to_ac_0GenericAffine.mat',
    '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/WIP/mammals/primates/Mouselemur/atlas/ednix/modified/matrices/from_MIRcen2018_0GenericAffine.mat',
    '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/WIP/mammals/primates/Macaque/atlas/ednix/modified/matrices/to_ac_0GenericAffine.mat',
    '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/WIP/mammals/carnivorans/Dog/atlas/ednix/modified/matrices/to_ac1_0GenericAffine.mat',
    '/srv/projects/easymribrain/data/Atlas/13_Atlas_project/WIP/mammals/carnivorans/Cat/atlas/ednix/modified/matrices/to_ac1_0GenericAffine.mat']

ref_brains  = [
    '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/ungulates/Pig/EDNiX/volumes/Pig_space-acpc_desc-SS_T1w.nii.gz',
    '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/primates/Mouselemur/EDNiX/volumes/Mouselemur_space-acpc_desc-SS_T1w.nii.gz',
    '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/primates/Macaque/EDNiX/volumes/Macaque_space-acpc_desc-SS_T1w.nii.gz',
    '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/carnivorans/Dog/EDNiX/volumes/Dog_space-acpc_desc-SS_T1w.nii.gz',
    '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/carnivorans/Cat/EDNiX/volumes/Cat_space-acpc_desc-SS_T1w.nii.gz']

### all species
for animal, transfo, ref_brain in zip(species, transformations, ref_brains):
    for lr in ['_LR', '']:

        new = ants.image_read(ref_brain)
        for i in range(4):
            orig = ants.image_read(opj(input_atlas,animal, 'atlaslvl' + str(i + 1) + lr + '.nii.gz'))
            moved = ants.apply_transforms(fixed=new, moving=orig,
                                          transformlist=transfo,
                                          interpolator='genericLabel',
                                          whichtoinvert=[True])
            ants.image_write(moved, opj(input_atlas,animal, animal + '_seg-atlaslvl' + str(i + 1) + '-' + lr[1:] + '_dseg.nii.gz'), ri=False)

species = [
    'Pig',
    'Mouse_lemur',
    'Marmoset',
    'Macaque',
    'Chimpanzee',
    'CatinDog',
    'DoginCat',
    'Mouse',
    'RatWHS',
    'Human']


for animal in species:
    print(animal)
    for lr in ['_LR', '']:

        lr_suffix = lr[1:]  # 'LR' or ''
        print(lr_suffix)
        new_atlas = ' '.join([
            opj(input_atlas, animal, animal + '_seg-atlaslvl1-' + lr_suffix + '_dseg.nii.gz'),
            opj(input_atlas, animal, animal + '_seg-atlaslvl2-' + lr_suffix + '_dseg.nii.gz'),
            opj(input_atlas, animal, animal + '_seg-atlaslvl3-' + lr_suffix + '_dseg.nii.gz'),
            opj(input_atlas, animal, animal + '_seg-atlaslvl4-' + lr_suffix + '_dseg.nii.gz'),
        ])

        if animal not in ['Pig','Mouse_lemur','Macaque','CatinDog','DoginCat','RatWHS']:
            new_atlas = ' '.join([
                opj(input_atlas,animal, 'atlaslvl1' + lr + '.nii.gz'),
                opj(input_atlas,animal, 'atlaslvl2' + lr + '.nii.gz'),
                opj(input_atlas, animal, 'atlaslvl3' + lr + '.nii.gz'),
                opj(input_atlas,animal, 'atlaslvl4' + lr + '.nii.gz')])

        out_path = opj(input_atlas, animal, animal + '_seg-EDNIxCSC' + lr_suffix + '_dseg.nii.gz')
        EDNIx_label = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/primates/Macaque/label_code/EDNIxCSCLR_label.txt'
        cmd = 'singularity run' + s_bind + afni_sif + '3dTcat -overwrite -prefix ' + out_path + ' ' + new_atlas
        nl = spgo(cmd)
        print(nl)


##################
species = ['Human']
ref_brain = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/mammals/primates/Human/EDNiX/volumes/Human_space-acpc_desc-SS_T1w.nii.gz'
#change space for human
### all species
for animal in species:
    for lr in ['_LR', '']:
        lr_suffix = lr[1:]  # 'LR' or ''

        new = ants.image_read(ref_brain)
        orig = ants.image_read(opj(input_atlas, animal, animal + '_seg-EDNIxCSC' + lr_suffix + '_dseg.nii.gz'))
        moved = ants.apply_transforms(fixed=new, moving=orig,
                                      transformlist=[
                                        '/srv/projects/easymribrain/data/Atlas/Atlases_V2/Human/template_to_Collin_SyN_0GenericAffine.mat',
                                        '/srv/projects/easymribrain/data/Atlas/Atlases_V2/Human/template_to_Collin_SyN_1InverseWarp.nii.gz'

                                          ],
                                      interpolator='genericLabel',
                                      whichtoinvert=[True, False],
                                      imagetype=3)
        ants.image_write(moved, opj(input_atlas,animal, animal + '_seg-EDNIxCSC' + lr_suffix + '_dseg_collin2.nii.gz'), ri=False)




