#import
import json
import time

def check(json1,json2):

    f2 = open(json1)
    info_map = json.load(f2)
    time_map = time.mktime(time.strptime(info_map["AcquisitionTime"], '%H:%M:%S.%f'))

    f = open(json2)
    info_bold = json.load(f)
    time_bold = time.mktime(time.strptime(info_bold["AcquisitionTime"], '%H:%M:%S.%f'))

    tdelta = time_bold - time_map

    return tdelta
