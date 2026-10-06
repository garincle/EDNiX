import os
import sys

opj = os.path.join

def get():
    # where is the EDNiX toolbox ?
    #MAIN_PATH = opj('srv','projects','easymribrain','code','EDNiX')


    MAIN_PATH = opj('media','clavagni','LaCie','travail','simon','current_works','pipeline','EDNiX')

    sys.path.insert(1, opj(MAIN_PATH))
    return MAIN_PATH

