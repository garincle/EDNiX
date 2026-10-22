"""Isotope and scan start of a PET DICOM file (used to build the SIF when the BIDS sidecar lacks them)."""
import os

opi = os.path.isfile

from modalities.pet import petutils as pu

# DICOM tags (https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.8.9.html)
RADIOPHARMACEUTICAL_SEQ = (0x0054, 0x0016)
RADIONUCLIDE_CODE_SEQ = (0x0054, 0x0300)
CODE_MEANING = (0x0008, 0x0104)
HALF_LIFE = (0x0018, 0x1075)


def extract(file):
    """
    (isotope 'F-18', scan start 'dd/mm/yyyy hh:mm:ss').
    The isotope is identified from the radionuclide half-life (0018,1075), which is reliable, and only
    then from the code meaning (e.g. '^18^Fluorine'): the PET engineer noticed wrong radionuclide codes
    in some files (C-11 data coded as F-18). The start is the acquisition (else series) date and time.
    """
    import pydicom as pdcm
    if not opi(file):
        raise FileNotFoundError(f"Input file not found: {file}")
    ds = pdcm.dcmread(file, stop_before_pixels=True)
    isotope = None
    if RADIOPHARMACEUTICAL_SEQ in ds:
        item = ds[RADIOPHARMACEUTICAL_SEQ][0]
        if HALF_LIFE in item:
            isotope = pu.isotope_from_half_life(item[HALF_LIFE].value)
        if isotope is None and RADIONUCLIDE_CODE_SEQ in item:
            code = item[RADIONUCLIDE_CODE_SEQ][0]
            if CODE_MEANING in code:
                isotope = pu.normalise_isotope(code[CODE_MEANING].value)
    day = str(getattr(ds, 'AcquisitionDate', '') or getattr(ds, 'SeriesDate', '') or getattr(ds, 'StudyDate', ''))
    tim = str(getattr(ds, 'AcquisitionTime', '') or getattr(ds, 'SeriesTime', '') or getattr(ds, 'StudyTime', ''))
    start = None
    if len(day) >= 8 and len(tim) >= 6:
        start = '/'.join([day[6:8], day[4:6], day[0:4]]) + ' ' + ':'.join([tim[0:2], tim[2:4], tim[4:6]])
    return isotope, start
