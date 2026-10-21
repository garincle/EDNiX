# import
import os
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,filename

def run(img,dftfile2_step,DV_time,diary_name,sing_turku):

    nl = ('# Computation of parametric image with Distribution Volume Ratio (DVR) '
          'applying the multiple time graphical analysis (MTGA) approach for reversible tracer uptake\n')
    run_cmd.msg(nl, diary_name,'OKGREEN')
    # http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgdv.html

    fname, _ = filename.get_desc(opb(img), 'desc')
    res_dir = opd(img)

    DVRimg  = opj(res_dir, '_'.join([fname, 'desc-DVR', 'pet.nii']))
    Plotimg = opj(res_dir, '_'.join([fname, 'desc-DVR-plot', 'pet.nii']))

    # check up
    for fneeded in [img, dftfile2_step]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    # execute

    # turku functions do not accept compressed gz files
    cmd = 'gzip -d ' + img
    run_cmd.run(cmd, diary_name)

    cmd = (sing_turku + 'imgdv -n=' + Plotimg.split('.')[0] +
           ' ' + dftfile2_step +
           ' ' + img.replace('.nii.gz', '.nii') +
           ' ' + str(DV_time) + ' ' + DVRimg)
    run_cmd.run(cmd, diary_name)

    for i in [img.replace('.nii.gz', '.nii'), DVRimg, Plotimg]:
        cmd = 'gzip -f ' + i
        run_cmd.run(cmd, diary_name)

    dictionary = {"Sources": [img,
                              dftfile2_step],
                  "Description": 'Computation of parametric image with Distribution Volume Ratio (DVR) '
                                 'applying the multiple time graphical analysis (MTGA) approach for '
                                 'reversible tracer uptake (turku,tpcclib: imgdv).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(DVRimg.replace('.nii','.json'), "w") as outfile:
        outfile.write(json_object)