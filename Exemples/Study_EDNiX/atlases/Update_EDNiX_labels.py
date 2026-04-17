# import
import os
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
import csv
import fnmatch
import ants
import subprocess

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile
spgo = subprocess.getoutput

from Tools import run_cmd,diaryfile
from atlases import label_translate
from anat.connectomeWB import WB_label

def run(LUTfile,species_list,Atlaspath,wantSurf,sing_wb):

    label_translate.EDNiXprepare(LUTfile,species_list)

    fname = opb(LUTfile)
    name = fname.split('_')[0]
    newname = name[:-2]

    for i in range(len(species_list)):

        # search for the relevant folder
        for dirpath, dirnames, filenames in os.walk(opj(Atlaspath)):
            folder_name = fnmatch.filter(dirnames, species_list[i])
            if folder_name:
                if not opb(dirpath) == 'freesurfer':
                    path_ref = opj(dirpath, folder_name[0])

        if 'path_ref' not in locals() or not path_ref:
            raise ValueError("path_ref is not defined — atlas path or species or MAIN_PATH is likely incorrect. " +
                             'specie = ' + str(species_list[i]) + ' path_ATLAS = ' + str(Atlaspath) + ' dirpath = ' + str(
                dirpath) + ' dirnames = ' + str(dirnames))

        path_label_code = opj(path_ref, 'label_code')


        labels_dir       = opj(path_ref, 'EDNiX', 'volumes', 'labels')
        dir_native_resol = opj(path_ref, 'EDNiX', 'surfaces', 'Native_resol')

        diary_file = diaryfile.create(labels_dir, 'update EDNiX ')


        for f in list_atlas[0]:
            volfile = opj(path_ref, 'EDNiX', 'volumes', 'labels', '_'.join([species_list[i], 'seg-' + f, 'dseg.nii.gz']))
            if opi(volfile):
                cmd = (sing_wb + 'wb_command -volume-label-import ' + volfile + ' ' + opj(path_label_code, f + '_label.txt') +
                       ' ' + volfile + ' -drop-unused-labels')
                print(cmd)
                spgo(cmd)

        if wantSurf ==1:
            surffile = opj(dir_native_resol,'.'.join([species_list[i], 'l', 'white', 'surf', 'gii']))

            if opi(surffile):
                WB_label.vol2surfWB(0,species_list[i],
                                    list_atlas,
                                    labels_dir, dir_native_resol, 'native', ['l','r'],
                                    'vox_dil', 'roi', opj(dir_native_resol, species_list[i] + '_native_LR.wb.spec'), diary_file,
                                    sing_wb)

        # for BALSA
        if species_list[i] in BALSA:

            dir_native_32   = opj(path_ref, 'EDNiX', 'surfaces', 'fsaverage_LR_32k')
            dir_balsavol    = opj(path_ref, BALSA[species_list[i]], 'volumes')
            dir_balsa       = opj(path_ref, BALSA[species_list[i]], 'surfaces')

            Ref_norm = opj(dir_balsavol, '_'.join([species_list[i], 'space-' + BALSA[species_list[i]], 'desc-SS', 'T1w.nii.gz']))

            matrices_dir = opj(opd(opd(Atlaspath)), 'BALSA', BALSA[species_list[i]], 'modified')

            BALSAsamp={'resol':['native','native','native_LR',''],
                       '164k':['fsaverage_164k','164k_fs_LR','fsaverage_LR_164k','atlasroi'],
                       '32k':['fsaverage_LR_32k','32k_fs_LR','fsaverage_LR_32k','atlasroi']}

            spherereg = [opj(dir_balsa,BALSAsamp['resol'][0],
                             '.'.join([species_list[i], Hmin[0], 'sphere','reg','reg_LR',BALSAsamp['resol'][1],'surf','gii'])),
                         opj(dir_balsa,BALSAsamp['resol'][0],
                             '.'.join([species_list[i], Hmin[1], 'sphere','reg','reg_LR',BALSAsamp['resol'][1],'surf','gii']))]

            ref = ants.image_read(Ref_norm)

            for j in list_atlas[0]:
                volfile  = opj(labels_dir,'_'.join([species_list[i], 'seg-' + j, 'dseg.nii.gz']))
                labelref = opj(dir_balsavol,'labels','_'.join([species_list[i], 'space-' + BALSA[species_list[i]],'seg-' + j, 'dseg.nii.gz']))
                img = ants.image_read(volfile)
                print("img.dimension=" + str(img.dimension))
                if img.dimension > 3:
                    dim4 = 3
                else:
                    dim4 = 0

                moved = ants.apply_transforms(fixed=ref, moving=img,
                                      transformlist=[opj(matrices_dir,'_'.join(['EDNiX','to',BALSA[species_list[i]],'SyN',antsaff])),
                                                     opj(matrices_dir,'_'.join(['EDNiX','to',BALSA[species_list[i]],'SyN',antswarp]))],
                                      imagetype=dim4,
                                      interpolator=labelinterp)
                ants.image_write(moved,labelref)
                cmd = (sing_wb + 'wb_command -volume-label-import ' + labelref + ' ' + opj(path_label_code, j + '_label.txt') +
                       ' ' + labelref+ ' -drop-unused-labels')
                run_cmd.wb(cmd, diary_file)

            if wantSurf == 1:
                print(('1'))
                print(BALSAsamp['resol'][0])
                WB_label.vol2surfWB(1,species_list[i],
                                    list_atlas,
                                    opj(dir_balsavol,'labels'), opj(dir_balsa,BALSAsamp['resol'][0]), BALSAsamp['resol'][1], ['l', 'r'],
                                    'vox_dil', 'roi', opj(dir_balsa,BALSAsamp['resol'][0], '_'.join([species_list[i], BALSAsamp['resol'][2] + '.wb.spec'])),
                                    diary_file,
                                    sing_wb)

                keys = list(BALSAsamp.keys())
                for B in range(1, 3):  # Pour les indices 1 et 2
                    key = keys[B]  # key sera '164k' (index 1) puis '32k' (index 2)
                    value = BALSAsamp[key]
                    dir = opj(dir_balsa, value[0])  # Premier élément de la liste
                    print(dir)
                    for h in range(2):
                        cmd = (sing_wb + 'wb_command -label-resample ' +
                               opj(dir_balsa,BALSAsamp['resol'][0], '.'.join([species_list[i],Hmin[h],j,BALSAsamp['resol'][1],'label','gii'])) +
                                   ' ' + spherereg[h] +
                                   ' ' + opj(dir, '.'.join([species_list[i], Hmin[h], 'sphere', value[1], 'surf','gii'])) +
                                   ' BARYCENTRIC ' + opj(dir, '.'.join([species_list[i],Hmin[h],j,value[1],'label','gii'])) + ' -largest')
                        run_cmd.wb(cmd, diary_file)

                    WB_label.cifti(species_list[i], j, dir,value[1],'.'.join([value[3],value[1]]),
                                                opj(dir, '_'.join([species_list[i], value[2] + '.wb.spec'])),
                                                diary_file,sing_wb)

Hmin = ['l','r']
list_atlas=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]]
#list_atlas=[['CIVM', 'D99v2', 'INIA19', 'CHARM', 'SARM'], ['txt', 'txt', 'ctab', 'txt', 'txt'], [1, 1, 1, 6, 6], [1, 1, 1, 1, 0]]
EDNIX_removelist = []

# for antspy
antsaff     = '0GenericAffine.mat'
antswarp    = '1Warp.nii.gz'
antsinvwarp = '1InverseWarp.nii.gz'

# interpolation
T1interp    = 'hammingWindowedSinc'
Boldinterp  = 'lanczosWindowedSinc'
labelinterp = 'genericLabel'
petinterp   = 'linear'

BALSA    = {"Human": 'S1200',"Chimpanzee": 'CY29',"Macaque":'MY19'}
BALSAsamp= {'resol':['native','native','native_LR',''],
                       '164k':['fsaverage_164k','164k_fs_LR','fsaverage_LR_164k','atlasroi'],
                       '32k':['fsaveçrage_LR_32k','32k_fs_LR','fsaverage_LR_32k','atlasroi']}


Atlaspath = opj('/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/')
LUTfile = '/home/cgarin/PycharmProjects/EDNiX/Atlases_library/atlas/ednix_lut/EDNIxCSCLR_StatsLUT.txt'
species_list = ['Cat', 'Dog']
sing_wb = 'vglrun singularity exec --bind /home/cgarin/PycharmProjects/EDNiX/Atlases_library/,/scratch2/,/scratch/ /home/cgarin/PycharmProjects/EDNiX/Singularity_library/workbench_2.1.0.sif '
#sing_wb = 'vglrun singularity run --bind /srv/projects/,/srv/projects/easymribrain,/scratch2/,/scratch/ /home/cgarin/PycharmProjects/EDNiX/Singularity_library/Singularity/connectome_workbench_1.5.0-freesurfer-update.sif '
run(LUTfile,species_list,Atlaspath,1,sing_wb)