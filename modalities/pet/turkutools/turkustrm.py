# import
import os
import json
import pandas as pd
import matplotlib.pyplot as plt

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd,filename

def bfm(dftfile1,siffile,list_rois,isotope,REF_name,dopict,diary_name,sing_turku):

    nl = '# Computation of parametric image with binding potential (BP) using simplified reference tissue model (SRTM) with basis function approach\n'
    run_cmd.msg(nl, diary_name,'OKGREEN')
    # http://www.turkupetcentre.net/petanalysis/tpcclib/doc/imgbfbp.html

    res_dir = opd(dftfile1)

    dftfilefit = dftfile1.replace('ref.dft','fit.dft')
    fileBP     = dftfile1.replace('ref.dft', 'bp.txt')

    # check up
    for fneeded in [dftfile1, siffile]:
        if not opi(fneeded):
            raise FileNotFoundError(f"Input file not found: {fneeded}")

    # execute
    cmd = (sing_turku + 'dftweigh ' + dftfile1 + ' ' + siffile)
    run_cmd.run(cmd, diary_name)

    cmd = (sing_turku + 'bfmsrtm -fit=' + dftfilefit + ' -i=' + isotope +
           ' ' + dftfile1 + ' ' + REF_name + ' ' + fileBP)
    run_cmd.run(cmd, diary_name)

    if dopict ==1:
        pict_dir = opj(res_dir, 'pictures')
        if ope(pict_dir) == False:
            os.makedirs(pict_dir)

        # figures to see the time evolution of the pet signal .............................................................

        file0 = pd.read_csv(dftfile1, sep='\t', header=None, skiprows=4)

        x = file0[0]
        y1 = file0[len(list_rois) + 2]
        fig, ax = plt.subplots()

        for i in range(16):  # see the first 16 rois
            y = file0[2 + i]
            ax.plot(x, y, 'o-', linewidth=1, label=list_rois[i])

        ax.plot(x, y1, 'ro-', linewidth=3, label='reference: ' + REF_name)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        # fig.legend(loc='upper right')

        plt.savefig(opj(pict_dir, 'Before_fit.png'), format='png')
        plt.savefig(opj(pict_dir, 'Before_fit.svg'), format='svg')
        plt.close()

        file1 = pd.read_csv(dftfilefit, sep='\t', header=None, skiprows=4,
                            skipfooter=2, engine='python')
        x = file1[0]
        y1 = file1[len(list_rois) + 2]
        for i in range(len(list_rois)):
            fig, ax = plt.subplots()
            y = file1[2 + i]
            ax.plot(x, y, 'o-', linewidth=2)
            ax.plot(x, y1, 'o-', linewidth=2, )
            ax.set_title(list_rois[i])
            plt.savefig(opj(pict_dir, list_rois[i] + '_fit.png'), format='png')
            plt.close()

def img(img,dftfile2_step,diary_name,sing_turku):

    nl = ('# Computation of parametric image with binding potential (BP) using simplified reference tissue model '
          '(SRTM) with linear approach\n')
    run_cmd.msg(nl, diary_name,'OKGREEN')
    # http://www.turkupetcentre.net/petanalysis/tpcclib/doc/bfmsrtm.html

    fname, _ = filename.get_desc(opb(img), 'desc')
    strm_img = opj(opd(img),'_'.join([fname, 'desc-STRM','pet.nii']))

    # turku functions do not accept compressed gz files
    cmd = 'gzip -d ' + img
    run_cmd.run(cmd, diary_name)

    cmd = (sing_turku + 'imgsrtm' + ' ' + img.replace('.nii.gz', '.nii') +
           ' ' + dftfile2_step +
           ' ' + strm_img)
    run_cmd.run(cmd, diary_name)

    for g in [img.replace('.nii.gz', '.nii'),strm_img]:
        cmd = 'gzip -f ' + g
        run_cmd.run(cmd, diary_name)

    dictionary = {"Sources": [img,
                              dftfile2_step],
                  "Description": 'Computation of parametric image with binding potential (BP) using '
                                 'simplified reference tissue model (SRTM) with linear approach '
                                 '(turku,tpcclib: imgsrtm).', }
    json_object = json.dumps(dictionary, indent=2)
    with open(strm_img.replace('.nii','.json'), "w") as outfile:
        outfile.write(json_object)