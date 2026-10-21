#import
import pydicom as pdcm
import os

opi = os.path.isfile

# info extracted
# from https://dicom.nema.org/medical/dicom/2025b/output/chtml/part16/sect_CID_18.html

def extract(file):

    if not opi(file):
        raise FileNotFoundError(f"Input file not found: {file}")

    ds1 = pdcm.dcmread(file)
    
    Iso_code = ds1[0x54,0x16][0][0x54,0x300][0][0x8,0x104].value
    ISO      = Iso_code.split('^')
    isotope  = ISO[0] + '-' + ISO[1]

    # beware that sometimes the Code within the dicom is wrong
    # (use of C-111A1 for Fluor when the isotope is something else)
    # so check halflife (around 110 min for F-18 and 20 min for C-11)
    halflife = ds1[0x54,0x16][0][0x18,0x1075].value
    # for an unknown reason the code value of C-11 (supposed to be 'C-105A1')

    # correction:
    if float(halflife)/60 <80:
        isotope ='C-11'

    Day_scan =  ds1[0x8,0x20].value
    acq_scan =  ds1[0x8,0x31].value

    START = ('/'.join([Day_scan[6:8], Day_scan[4:6], Day_scan[0:4]]) +
             ' ' + ':'.join([acq_scan[0:2], acq_scan[2:4], acq_scan[4:6]]))

    return isotope,START
