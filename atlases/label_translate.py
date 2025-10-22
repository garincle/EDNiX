# import
import os
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
import csv

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

def FS2_WB(label_file,REF):
    # for Workbench
    dir_label = opd(label_file)
    LABEL2 = pd.read_csv(label_file, sep=' ',
                         names=["Index", "Names", "C1", "C2", "C3", 'alpha'])

    with open(opj(dir_label, REF + '_label.txt'), "w") as f:
        for i in range(len(LABEL2['Index'])):
            line1 = LABEL2['Names'][i]
            line2 = ' '.join([str(LABEL2['Index'][i]),
                              str(LABEL2['C1'][i]),
                              str(LABEL2['C2'][i]),
                              str(LABEL2['C3'][i]),
                              '255'])
            if not LABEL2['Index'][i] == 0:
                f.write(line1 + '\n' + line2 + '\n')


def FS2_ITK(label_file, REF):
    dir_label = opd(label_file)

    LABEL = open(label_file, "r")
    lines = LABEL.readlines()

    with open(opj(dir_label, REF + '_ITK.txt'), "a") as f:

        # header
        f.write('################################################\n')
        f.write('# ITK-SnAP Label Description File\n')
        f.write('# File format: \n')
        f.write('# IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL\n')
        f.write('# Fields: \n')
        f.write('#    IDX:   Zero-based index \n')
        f.write('#    -R-:   Red color component (0..255)\n')
        f.write('#    -G-:   Green color component (0..255)\n')
        f.write('#    -B-:   Blue color component (0..255)\n')
        f.write('#    -A-:   Label transparency (0.00 .. 1.00)\n')
        f.write('#    VIS:   Label visibility (0 or 1)\n')
        f.write('#    IDX:   Label mesh visibility (0 or 1)\n')
        f.write('#  LABEL:   Label description \n')
        f.write('################################################\n')
        f.write('    0     0    0    0        0  0  0    "Clear Label"\n')

        for i in range(1, len(lines)):
            N = lines[i].split(' ')
            space1 = ' '
            for j in range(5 - len(N[0])):
                space1 = space1 + ' '
            space2 = ' '
            for j in range(4 - len(N[2])):
                space2 = space2 + ' '
            space3 = ' '
            for j in range(4 - len(N[3])):
                space3 = space3 + ' '
            space4 = ' '
            for j in range(8 - len(N[4])):
                space4 = space4 + ' '

            f.write('    ' + N[0] + space1 + N[2] + space2 + N[3] + space3 + N[4] + space4 + '1  1  1    "' + N[1] + '"\n')


def ITK2_FS(label_file,REF,LR,first):
    dir_label = opd(label_file)

    with open(label_file) as f:
        lines = f.readlines()

    if lines[1] == '# ITK-SnAP Label Description File\n':
        print('ITK_snap file')
    else:
        print('got it wrong')

    label_tmp = open(label_file.replace('.txt', '_tmp.txt'), "a")
    test = lines[3][2:].split()
    label_names = ' '.join(test[0:6]) + ' ' + test[7]
    label_tmp.write(f'{label_names}\n')
    for v in range(15, len(lines)):
        EE = lines[v].split('"')[1].split()
        RR = '_'.join(EE)
        TT = lines[v].split()
        aa = ' '.join(TT[0:6])
        label_tmp.write(f'{aa} {RR}\n')
    label_tmp.close()

    # For Freesurfer
    lines1 = list()
    lines2 = list()

    LABEL = pd.read_csv(label_file.replace('.txt', '_tmp.txt'), sep=" ")

    for i in range(0, len(LABEL['IDX'])):
        index = str(LABEL['IDX'][i])
        name = str(LABEL['LABEL'][i])
        if LR == 1:
            test = name.split('_')
            test.remove(test[first])
            name_lr = '_'.join(test)
        else:
            name_lr = name

        c_1 = str(LABEL['-R-'][i])
        c_2 = str(LABEL['-G-'][i])
        c_3 = str(LABEL['-B-'][i])

        LUT = ' '.join([index, name, c_1, c_2, c_3, '0\n'])
        lines1.append(LUT)

        if LR == 1:
            LUT_lr = ' '.join([index, name_lr, c_1, c_2, c_3, '0\n'])
            lines2.append(LUT_lr)

    if LR == 0:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f:
            for i in range(len(lines1)):
                f.write(lines1[i])
    elif LR == 1:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f,open(opj(dir_label, REF + '_l_ctab'), "w") as f_left, open(opj(dir_label, REF + '_r_ctab'), "w") as f_right:
            for i in range(len(lines1)):
                f.write(lines1[i])
                side = str(LABEL['LABEL'][i]).split('_')[first]
                if side.lower() == 'left' or side.lower() == 'l':
                    f_left.write(lines2[i])
                elif side.lower() == 'right' or side.lower() == 'r':
                    f_right.write(lines2[i])

    os.remove(label_file.replace('.txt', '_tmp.txt'))


def MRIcron2_FS(label_file,REF,LR,first,seg):

    dir_label = opd(label_file)
    lines1 = list()
    lines2 = list()

    # For Freesurfer
    LABEL = pd.read_csv(label_file, sep=" ", names=["Index", "Names", "unknown"])
    for i in range(len(LABEL['Index'])):
        index = str(LABEL['Index'][i])
        name  = LABEL['Names'][i]
        if LR == 1:
            test = name.split(seg)
            test.remove(test[first])
            name_lr = seg.join(test)
        else:
            name_lr = name

        c_1 = str(np.random.randint(256, size=1)).strip('[]')
        c_2 = str(np.random.randint(256, size=1)).strip('[]')
        c_3 = str(np.random.randint(256, size=1)).strip('[]')

        LUT = ' '.join([index, name, c_1, c_2, c_3, '0\n'])
        lines1.append(LUT)

        if LR == 1:
            LUT_lr = ' '.join([index, name_lr, c_1, c_2, c_3, '0\n'])
            lines2.append(LUT_lr)

    if LR == 0:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f:
            for i in range(len(lines1)):
                f.write(lines1[i])
    elif LR == 1:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f,open(opj(dir_label, REF + '_l_ctab'), "w") as f_left,open(opj(dir_label, REF + '_r_ctab'), "w") as f_right:
            for i in range(len(lines1)):
                f.write(lines1[i])
                side = str(LABEL['Names'][i]).split(seg)[first]
                if side.lower() == 'left' or side.lower() == 'l':
                    f_left.write(lines2[i])
                elif side.lower() == 'right' or side.lower() == 'r':
                    f_right.write(lines2[i])



def xml2_FS(label_file,REF,LR,first):

    dir_label = opd(label_file)
    lines1 = list()
    lines2 = list()

    LABEL = pd.read_xml(label_file, xpath="//label")
    for i in range(len(LABEL['index'])):

        index = str(LABEL['index'][i])

        NN = LABEL['label'][i].split(',')
        dummy = ' '.join(NN)
        NN = dummy.split(' ')
        name = '_'.join(NN)

        if LR == 1:
            test = name.split('_')
            test.remove(test[first])
            name_lr = '_'.join(test)
        else:
            name_lr = name

        c_1 = str(np.random.randint(256, size=1)).strip('[]')
        c_2 = str(np.random.randint(256, size=1)).strip('[]')
        c_3 = str(np.random.randint(256, size=1)).strip('[]')

        LUT = ' '.join([index, name, c_1, c_2, c_3, '0\n'])
        lines1.append(LUT)

        if LR == 1:
            LUT_lr = ' '.join([index, name_lr, c_1, c_2, c_3, '0\n'])
            lines2.append(LUT_lr)

    if LR == 0:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f:
            for i in range(len(lines1)):
                f.write(lines1[i])
    elif LR == 1:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f,open(opj(dir_label, REF + '_l_ctab'), "w") as f_left,open(opj(dir_label, REF + '_r_ctab'), "w") as f_right:
            for i in range(len(lines1)):
                f.write(lines1[i])

                NN = LABEL['label'][i].split(',')
                dummy = ' '.join(NN)
                NN = dummy.split(' ')
                name = '_'.join(NN)
                side = name.split('_')[first]
                if side.lower() == 'left' or side.lower() == 'l':
                    f_left.write(lines2[i])
                elif side.lower() == 'right' or side.lower() == 'r':
                    f_right.write(lines2[i])



def csv2_FS(label_file,REF,header,index_pos,name_pos,color_pos,color_type,left_pos,LR,first,seg,toadd):
    # For Freesurfer
    dir_label = opd(label_file)
    lines1 =list()
    if LR ==1:
        lines2 = list()

    LABEL = pd.read_csv(opj(dir_label,label_file), sep=",")
    for i in range(len(LABEL[header[0]])):
        index = str(LABEL[header[index_pos]][i])

        # check space:
        check = LABEL[header[name_pos]][i].split(' ')
        if len(check)>1:
            check.remove('')
            name = ' '.join(check)
        else: name = LABEL[header[name_pos]][i]

        if toadd == 1:
            name_lr = name
            side = LABEL[header[left_pos]][i].split(seg)[first]
            name = '_'.join([side,name])

        else:
            if LR==1:
                test = name[i].split(seg)
                test.remove(test[first])
                name_lr =seg.join(test)
            else:
                name_lr = name

        if color_pos == '':
            c_1 = str(np.random.randint(256, size=1)).strip('[]')
            c_2 = str(np.random.randint(256, size=1)).strip('[]')
            c_3 = str(np.random.randint(256, size=1)).strip('[]')
        else:
            if color_type == 'RGB':
                c_1 = str(LABEL[header[color_pos]][i])
                c_2 = str(LABEL[header[color_pos+1]][i])
                c_3 = str(LABEL[header[color_pos+2]][i])
            elif color_type == 'hex':
                RGB = mcolors.to_rgb(LABEL[header[color_pos]][i])
                c_1 = str(int(RGB[0] * 255))
                c_2 = str(int(RGB[1] * 255))
                c_3 = str(int(RGB[2] * 255))

        LUT = ' '.join([index,name,c_1,c_2,c_3,'0\n'])
        lines1.append(LUT)

        if LR==1:
            LUT_lr = ' '.join([index, name_lr, c_1, c_2, c_3, '0\n'])
            lines2.append(LUT_lr)

    if LR==0:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f:
            for i in range(len(lines1)):
                f.write(lines1[i])
    elif LR == 1:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f,open(opj(dir_label, REF + '_l_ctab'), "w") as f_left,open(opj(dir_label, REF + '_r_ctab'), "w") as f_right:
            for i in range(len(lines1)):
                f.write(lines1[i])
                side = LABEL[header[left_pos]][i].split(seg)[first]
                if side.lower() =='left' or side.lower() =='l':
                    f_left.write(lines2[i])
                elif side.lower() =='right' or side.lower() =='r':
                    f_right.write(lines2[i])


def BIDS2_FS(label_file, REF, LR, first,seg,toadd):
    # tsv file in BIDS standard or nonstandard
    # the header should look like :
    #  index name abbreviation color | index name abbreviation color mapping
    dir_label = opd(label_file)
    lines =list()
    LABEL = pd.read_csv(label_file, sep="\t")
    for i in range(len(LABEL['index'])):
        index = str(LABEL['index'][i])
        if toadd == 1:
            side = LABEL['abbreviation'][i].split(seg)[first]
            name = '_'.join([side, LABEL['name'][i]])
        else:
            name = LABEL['name'][i]

        RGB = mcolors.to_rgb(LABEL['color'][i])
        c_1 = str(int(RGB[0] * 255))
        c_2 = str(int(RGB[1] * 255))
        c_3 = str(int(RGB[2] * 255))

    LUT = ' '.join([index, name, c_1, c_2, c_3, '0\n'])
    lines.append(LUT)

    if LR == 0:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f:
            for i in range(len(lines)):
                f.write(lines[i])
    elif LR == 1:
        with open(opj(dir_label, REF + '_StatsLUT.txt'), "w") as f,open(opj(dir_label, REF + '_l_ctab'), "w") as f_left,open(opj(dir_label, REF + '_r_ctab'), "w") as f_right:
            for i in range(len(lines)):
                f.write(lines[i])
                side = LABEL['abbreviation'][i].split(seg)[first]
                if side.lower() == 'left' or side.lower() == 'l':
                    f_left.write(lines[i])
                elif side.lower() == 'right' or side.lower() == 'r':
                    f_right.write(lines[i])


def FS2_BIDS(label_folder,REF):

    LABEL2 = pd.read_csv(opj(label_folder, REF + '_StatsLUT.txt'),
                         sep=' ',
                         names=["Index", "Names", "C1", "C2", "C3", 'alpha'])
    tsvFileName = opj(label_folder, REF + '.tsv')
    tsv_hd = ['index','name','abbreviation','color']
    tsv_row = list()

    for i in range(len(LABEL2['Index'])):
        color_label = (LABEL2['C1'][i], LABEL2['C2'][i], LABEL2['C3'][i])
        color = mcolors.to_hex(tuple(v / 255. for v in color_label))
        tsv_row.append([LABEL2['Index'][i],LABEL2['Names'][i],LABEL2['Names'][i],color])

    with open(tsvFileName, "w") as tsvfile:
        tsv_writer = csv.writer(tsvfile, delimiter='\t')
        tsv_writer.writerow(tsv_hd)
        tsv_writer.writerows(tsv_row)
        tsvfile.close()




