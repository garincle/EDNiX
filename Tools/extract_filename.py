import os

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile


def extract_filename(file_path):
    # Use os.path.basename to get the base name of the file
    base_name = opb(file_path)
    # Check if the file has a ".gz" extension and remove it if present
    if base_name.endswith(".nii.gz"):
        base_name = base_name[:-7]
    # Return the desired part of the file path
    elif base_name.endswith(".nii"):
        base_name = base_name[:-4]
    # Return the desired part of the file path
    return base_name

def change(dir,file,old,new):
    if opi(file):
        #dir = opd(file)
        oldname = opb(file).split('.')
        listN = oldname[0].split('_')
        n=0
        for i,j in enumerate(listN):
            if old in j:
                listN[i]='-'.join([old,new])
                newname = opj(dir, '.'.join(['_'.join(listN)] + oldname[1:]))
                n=1

        if n==0:
            NN = '_'.join(listN[0:-1] + ['-'.join([old,new])] + [listN[-1]])
            newname = opj(dir,'.'.join([NN] + oldname[1:]))
        return newname
    else:
        print('sorry, impossible to find the requested file: ' + file)

def move(dir,file,old,new):
    if opi(file):
        #dir = opd(file)
        oldname = opb(file).split('.')
        listN = oldname[0].split('_')
        n=0
        for i,j in enumerate(listN):
            if old in j:
                listN[i]=new
                newname = opj(dir, '.'.join(['_'.join(listN)] + oldname[1:]))
                n=1

        if n==0:
            NN = '_'.join(listN[0:-1] + ['-'.join([old,new])] + [listN[-1]])
            newname = opj(dir,'.'.join([NN] + oldname[1:]))
        return newname
    else:
        print('sorry, impossible to find the requested file: ' + file)

def get_subname(s):
    return s.split('_')[0]

def get_runname(s):
    return '_'.join(s.split('_')[:3])

def get_typename(s):
    return '_'.join(s.split('_')[:2])

def get_endname(s):
    # [4] is supposed to be the place for desc-.....

    return '_'.join(s.split('_')[4:])

def get_desc(s,spec):

    N = opb(s).split('_')
    E = [i for i, elem in enumerate(N) if spec in elem]
    if E==[]:
        E=[-1]
    fname    = '_'.join(N[:E[0]])
    specname = N[E[0]]

    return fname,specname