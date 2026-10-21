import os
import numpy as np
import pandas as pd
import ants

opj = os.path.join
opb = os.path.basename

from tools import run_cmd

def get(img,dest,atlas,atlas_ctab,option,suffix,sing_afni,diary_name):
    imgnii = ants.image_read(atlas)
    if len(imgnii.shape) == 4:
        nb_vol = imgnii.shape[3]
        img_unmerged = ants.ndimage_to_list(imgnii)
    else:
        nb_vol = 1

    for V in range(nb_vol):
        if nb_vol > 1:
            roifile = opj(dest, opb(img.replace('.nii.gz', '_'.join(['-' + str(V), suffix,'ROIs.txt']))))
        else:
            roifile = opj(dest, opb(img.replace('.nii.gz', '_'.join(['', suffix, 'ROIs.txt']))))

        cmd = (sing_afni + '3dROIstats ' + option + '-quiet -mask ' + atlas + '[' + str(V) + '] ' + img)
        Bp, _ = run_cmd.get(cmd, diary_name)

        list1D = Bp.decode("utf-8").split('\n')
        Bp1d = open(roifile, 'w')
        for i in range(len(list1D)):
            line = '\t'.join(list1D[i].split()) + '\n'
            Bp1d.write(line)
        Bp1d.close()

        df = pd.read_csv(roifile, sep='\t', header=None,dtype=np.float64)
        # remove empty columns
        df.dropna(how='all', axis=1, inplace=True)

        ROIS = pd.read_csv(atlas_ctab, sep='\t')
        list_labels = list(ROIS['ROIs'])
        list_rois   = list(ROIS['labels'])

        #  3dROIstats (AFNI's function) read labels in ascending order only and the ctab are not necessarily organized as such
        order = np.arange(len(list_labels))
        reorderlabel, neworder = zip(*sorted(zip(list_labels, order)))

        reorderlist_rois = list(map(list_rois.__getitem__, neworder))
        _, reverseorder = zip(*sorted(zip(neworder, order)))

        if nb_vol > 1:
            rr = np.setdiff1d(list_labels, img_unmerged[V].unique()[1:])
        else:
            rr = np.setdiff1d(list_labels, img.unique()[1:])

        # add missing regions (those too small that have been missing during the normalization such as habenula nucleus in low resolution)
        if len(rr) > 0:
            for z in rr:
                e = reorderlabel.index(z)
                nl = 'beware : the column nb :' + str(e) + ' is missing in the segmented file'
                run_cmd.msg(nl, diary_name)
                df.insert(loc=e, column='new_' + str(e), value='Inf')

        df.columns = list(order)
        df = df[list(reverseorder)]

        df.to_csv(roifile, sep='\t', header=None, index=None)

    


