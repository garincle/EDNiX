import os
import ants
import json
import datetime
import time
import numpy as np

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,filename
from pet import roi_signal

def getframes(funcfile,T1_eq,funcjson,petjson):

    f = open(petjson)
    info_pet = json.load(f)
    FS = info_pet["FrameTimesStart"]
    startpet = datetime.datetime.strptime(info_pet["AcquisitionTime"], '%H:%M:%S.%f')
    f.close()
    addtime = []
    for i in FS:
        test = startpet + datetime.timedelta(seconds=i)
        addtime.append(test.strftime('%H:%M:%S'))

    run_hd = ants.image_header_info(funcfile)
    duration = T1_eq + run_hd["dimensions"][3] * run_hd["spacing"][3]
    f = open(funcjson)
    info_run = json.load(f)
    startrun = datetime.datetime.strptime(info_run["AcquisitionTime"], '%H:%M:%S.%f') + datetime.timedelta(
        seconds=T1_eq * run_hd["spacing"][3])
    test2 = startrun.strftime('%H:%M:%S.%f')
    stoprun = datetime.datetime.strptime(info_run["AcquisitionTime"], '%H:%M:%S.%f') + datetime.timedelta(
        seconds=duration)
    test4 = stoprun.strftime('%H:%M:%S.%f')
    #print([test2, test4, str(stoprun - startrun)])

    test3 = []
    test5 = []
    for i, j in enumerate(addtime):
        a = time.mktime(time.strptime(j, '%H:%M:%S'))
        b = time.mktime(time.strptime(test2, '%H:%M:%S.%f'))
        c = time.mktime(time.strptime(test4, '%H:%M:%S.%f'))
        test3.append(a - b)
        test5.append(a - c)

    mindiffstart = np.min(np.abs(test3))
    mindiffstop  = np.min(np.abs(test5))

    for i, j in enumerate(addtime):
        if np.abs(test3[i]) == mindiffstart:
            if test3[i] < 0:
                startframe = i
                d = time.mktime(time.strptime(j, '%H:%M:%S'))
            else:
                startframe = i - 1
                d = time.mktime(time.strptime(addtime[i - 1], '%H:%M:%S'))

        if np.abs(test5[i]) == mindiffstop:
            if test5[i] < 0:
                stopframe = i + 1
                e = time.mktime(time.strptime(addtime[i + 1], '%H:%M:%S'))
            else:
                stopframe = i
                e = time.mktime(time.strptime(j, '%H:%M:%S'))

    # print([str(stoprun3 - startrun3), str(e - d)])
    return startframe,stopframe


def cut(anadir,funcfile, petfile,dseg,startframe, stopframe, atlas_label, diary_name, sing_afni, sing_wb):

    Palette = (' MODE_USER_SCALE -thresholding THRESHOLD_TYPE_OFF THRESHOLD_TEST_SHOW_INSIDE 0 2.5' +
               ' -pos-user 0 2.5 -interpolate true -palette-name ROY-BIG -disp-pos true -disp-neg false -disp-zero false')

    # check up
    for fneeded in [funcfile, petfile,dseg]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    listn = []
    for t in ['task','run']:
        _, n = filename.get_desc(opb(funcfile), t)
        listn.append(n.split('-')[-1])

    desc1   = '-'.join(['desc'] + listn)
    desc2   = '-'.join(['desc','mean'] + listn)
    desc3   = '-'.join(['desc', 'mean'] + listn + ['bratio'])
    endfile = 'pet.nii.gz'

    fname,_ = filename.get_desc(opb(petfile),'desc')


    step1    = opj(anadir,'_'.join([fname, desc1, endfile]))
    step2    = opj(anadir,'_'.join([fname, desc2, endfile]))
    finalimg = opj(anadir,'_'.join([fname, desc3, endfile]))
    
    if not ope(anadir):
        os.makedirs(anadir)
    nl = 'extract and mean the pet signal during the functional run'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    cmd = (sing_afni + '3dTcat -prefix ' + step1 + ' ' + petfile + '[' + str(startframe) + '-' +str(stopframe) + ']')
    run_cmd.run(cmd, diary_name)

    cmd = (sing_afni + '3dTstat -mean -prefix ' + step2 +' ' + step1)
    run_cmd.run(cmd, diary_name)

    roi_signal.read(step2, dseg, atlas_label, 0, 0,'', diary_name, sing_afni)

    img     = ants.image_read(step2)
    meanS   = img[img > 0].mean()
    new_img = ants.image_clone(img)
    new_img = new_img / meanS
    ants.image_write(new_img, finalimg)

    cmd = (sing_wb + 'wb_command -volume-palette ' + finalimg + Palette)
    run_cmd.run(cmd, diary_name)

    dictionary = {"Sources": [funcfile,
                              petfile],
                  "Description": 'cutting of the pet file according to the timing ogf the fMRI run file and normalisation of the signal by the brain mean (ANTspy)', }
    json_object = json.dumps(dictionary, indent=2)
    with open(finalimg.replace('.nii.gz','.json'), "w") as outfile:
        outfile.write(json_object)


