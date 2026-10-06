#import
import pydicom as pdcm


def extract(file):
    ds1 = pdcm.dcmread(file)
    
    Iso_code = ds1[0x54,0x16][0][0x54,0x300][0][0x8,0x104].value
    ISO      = Iso_code.split('^')
    isotope  = ISO[0] + '-' + ISO[1]
    
    Day_scan =  ds1[0x8,0x20].value
    acq_scan =  ds1[0x8,0x31].value

    START = Day_scan[6:8] + '/' + Day_scan[4:6] + '/' + Day_scan[0:4] + \
        ' ' + acq_scan[0:2] + ':' + acq_scan[2:4] + ':' + acq_scan[4:6]
    
    return isotope,START
