#import
import os
import json
import ants
import numpy as np

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile

from Tools import run_cmd,filename

def setTime(file,duration):

    steps, frame_stop= [], []

    if opi(file.replace('.nii.gz', '.sif')):
        f = open(file.replace('.nii.gz', '.sif'), "r")
        matrix = f.readlines()
        matrix = matrix[1:]
        starts = np.zeros(len(matrix))
        stops = np.zeros(len(matrix))
        Bq = np.zeros(len(matrix))
        for q in range(len(matrix)):
            starts[q] = int(matrix[q].split(' ')[0])
            stops[q] = int(matrix[q].split(' ')[1])
            Bq[q] = int(matrix[q].split(' ')[2])

        steps.append(np.arange(starts[0], stops[-1], duration * 60))  # blocks
        frame_stop.append(stops)

    else:
        print(file.replace('.nii.gz', '.sif') + " not found")

    return steps, frame_stop

def doit(img,maskfile,steps,stops,duration,ratio,diary_name,sing_afni):

    Dir_pet = opd(img)
    fname,suffix = filename.get_desc(opb(img), 'desc')
    info = suffix.split('-')[1:]
    if len(info)>1:
        infodesc = '-'.join(info)
    else:
        infodesc = info[0]

    for fneeded in [img,maskfile]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

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
    files=[]
    for j in range(len(tt)):

        nl = 'Create a mean image per ' + str(duration) + 'min long blocks\n'
        run_cmd.msg(nl, diary_name,'OKGREEN')

        new_file = opj(Dir_pet, '_'.join([fname,'-'.join(['desc',infodesc,'period',str(int(steps[j] / 60))]),'pet.nii.gz']))
        files.append(new_file)
        if j < len(tt)-1:
            cutting = '[' + str(tt[j]) + '-' +str(tt[j+1]-1) + ']'
        else:
            cutting = '[' + str(tt[j]) + '-' + str(nb_vol-1) + ']'

        cmd = (sing_afni + '3dTstat -mean -prefix ' + new_file + ' ' + img+cutting)
        run_cmd.run(cmd,diary_name)

        IMG = ants.image_read(new_file)
        new_IMG = ants.image_clone(IMG)
        msk = ants.image_read(maskfile)

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
        run_cmd.msg(nl, diary_name,'OKGREEN')

    return files










