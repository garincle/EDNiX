import os
import json
from nilearn.image import math_img
import ants
import numpy as np

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from Tools import check_nii

def check(T1wfile, T2wfile,sing_afni,diary_name):

    test = []
    cmd = (sing_afni + '3dinfo -same_obl ' + T1wfile + ' ' + T2wfile)
    var,_ = run_cmd.get(cmd, diary_name)
    list1D = var.decode("utf-8").split('\n')[:2]
    tt = sum(map(int,list1D))
    if tt == 2:
        test.append(True)
    else:
        test.append(False)
        nl = 'Beware : the two images have different oblique'
        run_cmd.msg(nl, diary_name, 'WARNING')

    hd1 = ants.image_header_info(T1wfile)
    hd2 = ants.image_header_info(T2wfile)

    if hd1['dimensions'] != hd2['dimensions']:
        nl = 'Beware : the two images do not have the same size'
        run_cmd.msg(nl, diary_name, 'WARNING')
        nl = 'if it is the only issue the myelin process can still be done'
        run_cmd.msg(nl, diary_name, 'WARNING')


    test.append(hd1['spacing'] == hd2['spacing'])
    ori = all(np.ceil(hd1['origin']) == np.ceil(hd2['origin']))
    if ori == False:
        dist1 = np.sqrt(
            np.sqrt(np.power((hd1['origin'][0] - hd2['origin'][0]), 2) +
                    np.power((hd1['origin'][1] - hd2['origin'][1]), 2) +
                    np.power((hd1['origin'][2] - hd2['origin'][2]), 2)))
        dist2 = np.sqrt(
            np.sqrt(np.power(hd1['spacing'][0], 2) +
                    np.power(hd1['spacing'][1], 2) +
                    np.power(hd1['spacing'][2], 2)))

        if dist2 - dist1 > 0:
            ori = True
        else:
            nl = 'Beware : the two images do not have the same origin'
            run_cmd.msg(nl, diary_name, 'WARNING')

    test.append(ori)
    for i in range(3):
        res = all(np.around(hd1['direction'][i], 3) == np.around(hd2['direction'][i], 3))
        if res == False:
            print('beware that the geometry seems to differ between the two files')
            print(hd1['direction'][i])
            print(hd2['direction'][i])
        test.append(res)


    if not all(test):
        nl = 'Sorry the two images cannot be combined for estimating the myelin distribution'
        run_cmd.msg(nl, diary_name, 'FAIL')
        return False
    else :
        return True


def divT1T2v2(T1, T2, type_norm, anatname,newfile,diary_name, sing_wb, sing_afni):

    if T1 == type_norm:
        refimg = T1
        padfile = T2
        T1wfile = anatname + refimg + '.nii.gz'
        T2wfile = anatname + 'formyeline_' + padfile + '.nii.gz'

    else:
        refimg = T2
        padfile = T1
        T2wfile = anatname + refimg + '.nii.gz'
        T1wfile = anatname + 'formyeline_' + padfile + '.nii.gz'


    nl = 'Combining T1w and T2w images V2'
    run_cmd.msg(nl, diary_name,'HEADER')

    command = (sing_afni + '3dcalc -overwrite -a ' + anatname + padfile + '.nii.gz'  +
               ' -prefix ' + anatname + 'formyeline_' + padfile + '.nii.gz' +
               ' -expr "a"')
    run_cmd.do(command, diary_name)

    check_nii.resamp(anatname + 'formyeline_' + padfile + '.nii.gz', anatname + refimg + '.nii.gz', 'msk', '', '', diary_name,
                     sing_wb)

    cmd = (sing_wb + 'wb_command -volume-math "clamp((T1w/T2w),0,100)" ' + newfile +
           ' -var T1w ' + T1wfile +
           ' -var T2w ' + T2wfile + ' -fixnan 0')
    run_cmd.do(cmd, diary_name)

    dictionary = {"Sources": [T1wfile,
                              T2wfile],
                  "Description": 'T1w Divided By T2w', }
    json_object = json.dumps(dictionary, indent=2)
    with open(newfile.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)


def newversion(T1wfile,T2wfile,ref,list_transfo,newfile,dir_transfo,diary_name,sing_afni,sing_wb):

    if 'T1' in ref:
        reffile = T1wfile
        padfile = T2wfile
        name = 'T2'

    if 'T2' in ref:
        reffile = T2wfile
        padfile = T1wfile
        name = 'T1'


    '''
        command = (sing_afni + '3dcalc -overwrite -a ' + padfile +
               ' -prefix ' + padfile.replace('_' + name + '.nii.gz', '_desc-backup_' + name + '.nii.gz') +
               ' -expr "a"')
    print(command)
    #run_cmd.run(command, diary_name)
    import subprocess
    subprocess.run(command, shell=True, check=True)
    dictionary = {"Sources": padfile,
                  "Description": 'backup of the ' + padfile + ' image before rotation.', }
    json_object = json.dumps(dictionary, indent=2)
    with open(padfile.replace('_' + name + '.nii.gz', '_desc-backup_' + name + '.json'), "w") as outfile:
        outfile.write(json_object)
        
    cmd = ('3dWarp -deoblique -wsinc5 -prefix ' + opj(main_path, 'T1w.nii.gz') + ' ' + file1)
    run_cmd.run(cmd, diary_name)

    cmd = ('3dWarp -deoblique -wsinc5 -prefix ' + opj(main_path, 'T2w.nii.gz') + ' ' + file2)
    run_cmd.run(cmd, diary_name)
    '''
    img1 = ants.image_read(reffile)
    img2 = ants.image_read(padfile)
    refnb = 0
    for i, j in enumerate(list_transfo):
        if list_transfo[i]["name"] == 'between':
            refnb = i

    mTx0 = ants.registration(fixed=img2, moving=img1,
                            type_of_transform='Translation',aff_metric=list_transfo[refnb]["affmetricT"])

    mTx = ants.registration(fixed=img2, moving=img1,initial_transform=mTx0['fwdtransforms'],
                            type_of_transform=list_transfo[refnb]["type_of_transform"],
                            aff_metric=list_transfo[refnb]["affmetric"],
                            outprefix=opj(dir_transfo, 'from_' + ref + '_to_' + name + '_'))

    #new_img = mTx['warpedfixout']
    #ants.image_write(new_img, opj(main_path, 'test2.nii.gz'))

    moved = ants.apply_transforms(fixed=img1, moving=img2,
                                  transformlist=mTx['invtransforms'],
                                  interpolator=list_transfo[refnb]["interpol"],
                                  whichtoinvert=[True])
    ants.image_write(moved, padfile)

    cmd = (sing_afni + '3dZeropad -overwrite -master ' + reffile + ' -prefix ' + padfile + ' ' + padfile)
    run_cmd.run(cmd, diary_name)

    cmd = (sing_wb + 'wb_command -volume-math "clamp((T1w/T2w),0,100)" ' + newfile +
           ' -var T1w ' + T1wfile +
           ' -var T2w ' + T2wfile + ' -fixnan 0')
    run_cmd.do(cmd, diary_name)

    dictionary = {"Sources": [T1wfile,
                              T2wfile],
                  "Description": 'T1w Divided By T2w', }
    json_object = json.dumps(dictionary, indent=2)
    with open(newfile.replace('.nii.gz', '.json'), "w") as outfile:
        outfile.write(json_object)