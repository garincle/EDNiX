from atlases import template
from pathlib import Path
import fnmatch
import os

opj = os.path.join
opi = os.path.isfile
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists


atlas_followers=[['EDNIxCSCLR', 'EDNIxCSC'], ['ctab', 'txt'], [4, 4], [1, 1]]
MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
Hmin    = ['l','r']
smallbrainspecieslist = ['Mouse','Rat','Mouselemur','Bat']
for species in ['Mouse','Rat','Mouselemur', 'Marmoset']:
    reference = 'EDNiX' # search for the relevant folder
    for dirpath, dirnames, filenames in os.walk(os.path.join('/home/cgarin/PycharmProjects/EDNiX/Atlases_library/', 'atlas')):
        folder_name = fnmatch.filter(dirnames, species)
        if folder_name:
            if not opb(dirpath) == 'freesurfer':
                path_ref = os.path.join(dirpath, folder_name[0])
    bids_dir = path_ref
    diary_file = path_ref + '/diary.txt'
    p = Path(diary_file); p.parent.mkdir(parents=True, exist_ok=True); p.touch()
    template.modif(smallbrainspecieslist, species,reference,atlas_followers,MAIN_PATH,bids_dir,diary_file,proj='ribbon',doflat=0)