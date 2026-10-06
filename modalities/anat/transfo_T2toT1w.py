#import
import os
import json
import ants

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname

from Tools import run_cmd

def prepa(file,mask,diary_name):

    nl = 'transformation of T2w image into a close version of a T1w image'
    run_cmd.msg(nl, diary_name,'OKGREEN')

    img = ants.image_read(file)
    msk = ants.image_read(mask)

    dir_T2 = opd(file)
    name = '_'.join(opb(file).split('_')[0:-1])

    new_img = img*msk
    R = 110 / (new_img.max() - new_img.min())
    T1w = (new_img*(-1) + new_img.max()) * R * msk
    ants.image_write(T1w, opj(dir_T2, name + '_T1w.nii.gz'), ri=False)

    dictionary = {"Sources": [file,
                             mask],
                  "Description": 'transformation of T2w image into a close version of a T1w image (Antspy).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(opj(dir_T2, name + '_T1w.json'), "w") as outfile:
        outfile.write(json_object)
