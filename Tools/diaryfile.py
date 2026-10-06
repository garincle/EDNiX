#import
import os
import datetime

opi = os.path.isfile

def create(file,msg):

    date_file = datetime.date.today()
    name = file + str(date_file) + '.txt'
    ct = datetime.datetime.now()
    if not opi(name):
        diary = open(name, "w")
        diary.write(msg)
    else:
        diary = open(name, "a")
    diary.write(f'\n{ct}')
    diary.write(f'\n')
    diary.close()

    return name
