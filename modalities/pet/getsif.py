import os

from modalities.pet import sifFile

def get(all_ID, all_Session, all_data_path):
    for ID, Session, data_path, in zip(all_ID, all_Session, all_data_path):

        sifFile.check(file, dir_dcm, diary_file, endfile, sing_turku)