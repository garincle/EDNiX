# import
import os
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,filename

def run(img,dftfile2_step,diary_name,sing_turku):


    nl = ('# Computation of parametric image with binding potential (BP) using simplified reference tissue model '
          '(SRTM) with basis function approach\n')
    run_cmd.msg(nl, diary_name,'OKGREEN')
    # http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgbfbp.html

    fname, _ = filename.get_desc(opb(img), 'desc')
    res_dir = opd(img)
    BPimg    = opj(res_dir, '_'.join([fname, 'desc-BP', 'pet.nii']))
    BPerrimg = opj(res_dir, '_'.join([fname, 'desc-BP-err', 'pet.nii']))

    # check up
    for fneeded in [img, img.replace('.nii.gz', '.sif'), dftfile2_step]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    # execute

    # turku functions do not accept compressed gz files
    cmd = 'gzip -d ' + img
    run_cmd.run(cmd, diary_name)

    # main function
    cmd = (sing_turku + 'imgbfbp -err=' + BPerrimg +
           ' ' + img.replace('.nii.gz', '.nii') +
           ' ' + dftfile2_step +
           ' ' + BPimg +
           ' ' + img.replace('.nii.gz', '.sif'))
    run_cmd.run(cmd, diary_name)

    for i in [img.replace('.nii.gz', '.nii'), BPimg, BPerrimg]:
        cmd = 'gzip -f ' + i
        run_cmd.run(cmd, diary_name)


    dictionary = {"Sources": [img,
                              dftfile2_step,
                              img.replace('.nii.gz', '.sif')],
                  "Description": 'Computation of parametric image with binding potential (BP) '
                                 'using simplified reference tissue model (SRTM) '
                                 'with basis function approach (turku,tpcclib: imgbfbp).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(BPimg.replace('.nii', '.json'), "w") as outfile:
        outfile.write(json_object)

