#import
import os
import sys
import json
import ants
import numpy as np

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname

sys.path.insert(1, '../tools')
import run_cmd
sys.path.insert(1, 'pet')
import paramImg

def doit(img,steps,stops,duration,ratio,diary_name,sing_afni):

    Dir_pet = opd(img)
    N = opb(img).split('_')
    E=-1
    for i in range(len(N)):
        a = N[i].split('-')
        if a[0] == 'desc':
            E = i
    Name =  '_'.join(N[0:E])

    hd = ants.image_header_info(img)
    nb_vol = int(hd['dimensions'][3])

    tt = np.zeros(len(steps))
    for j in range(len(steps)):
        a = stops-steps[j]
        diff=np.min(a[a>=0])
        for index, stop in enumerate(stops):
            diff0 = stop-steps[j]
            if diff0==diff:
                tt[j] = index +1
    tt = tt.astype(int)
    print(tt)
    for j in range(len(tt)):

        nl = 'Create a mean image per ' + str(duration) + 'min long blocks\n'
        run_cmd.msg(nl, diary_name)

        new_file = opj(Dir_pet, Name + '_desc-mean+' + str(int(steps[j] / 60)) + '_pet.nii.gz')
        if j < len(tt)-1:
            cutting = '[' + str(tt[j]) + '-' +str(tt[j+1]-1) + ']'
        else:
            cutting = '[' + str(tt[j]) + '-' + str(nb_vol-1) + ']'

        cmd = (sing_afni + '3dTstat -mean -prefix ' + new_file + ' ' + img+cutting)
        run_cmd.run(cmd,diary_name)

        IMG = ants.image_read(new_file)
        new_IMG = ants.image_clone(IMG)
        msk = ants.image_read(opj(Dir_pet, Name + '_desc-coreg_mask.nii.gz'))

        msk[msk>0]=1
        mean_brain = ants.mask_image(new_IMG, msk, 1)
        new_IMG= new_IMG/mean_brain.mean()

        ants.image_write(new_IMG,new_file, ri=False)

        dictionary = {"Sources": [img + cutting,
                                  str(ratio)],
                      "Description": 'Mean image (3dTstat, AFNI and ANtspy).', }
        json_object = json.dumps(dictionary, indent=2)
        with open(new_file.replace('.nii.gz', '.json'), "w") as outfile:
            outfile.write(json_object)

        nl = ' mean value of ' + opb(new_file) + ': ' + str(new_IMG.mean())
        run_cmd.msg(nl, diary_name)










