"""
Scan Information Files (SIF, Turku PET Centre format) of the dynamic PET images.

NIfTI images carry neither the frame times nor the isotope: the Turku tools (tpcclib) read them from
a SIF with the same base name (<image>.sif). Format (https://www.turkupetcentre.net/petanalysis/sif.html):
    title line : scan start 'dd/mm/yyyy hh:mm:ss', number of frames, number of columns (4),
                 SIF version (1), study id, isotope (e.g. F-18)
    frame lines: frame start (s), frame end (s), prompts, randoms
Here the counts are estimated from the image by imgweigh (tpcclib), which is what the weights of the
kinetic fits are computed from.
"""
import os
import glob
import shutil
import datetime

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
opi = os.path.isfile
ope = os.path.exists

from Tools import run_cmd
from modalities.pet.turkutools import turkusif


def start_string(pet_info, dicom_start=None):
    """Scan start for the SIF title line: DICOM date/time when available, else BIDS TimeZero + ScanStart
    (BIDS has no date: 01/01/1970 is written, which tpcclib only uses for display)."""
    if dicom_start:
        return dicom_start
    t0 = pet_info.get('TimeZero') or pet_info.get('AcquisitionTime')
    try:
        t = datetime.datetime.strptime(str(t0).split('.')[0], '%H:%M:%S')
    except (TypeError, ValueError):
        t = datetime.datetime(1970, 1, 1)
    t = t + datetime.timedelta(seconds=float(pet_info.get('ScanStart') or 0.0))
    return '01/01/1970 ' + t.strftime('%H:%M:%S')


def write(sif_path, start, timing, isotope, study='EDNiX'):
    """SIF with the frame times of `timing` (counts set to 0, to be estimated by imgweigh)."""
    n = len(timing['start'])
    with open(sif_path, 'w') as f:
        f.write(f'{start} {n} 4 1 {str(study)[:8]} {isotope}\n')
        for s, e in zip(timing['start'], timing['end']):
            f.write(f'{s:g} {e:g} 0 0\n')
    return sif_path


def read(file):
    """(start 'dd/mm/yyyy hh:mm:ss', isotope, frame starts, frame ends) of a SIF."""
    with open(file) as f:
        lines = [l.split() for l in f if l.strip()]
    title = lines[0]
    start = ' '.join(title[0:2])
    isotope = title[-1]
    starts = [float(l[0]) for l in lines[1:]]
    ends = [float(l[1]) for l in lines[1:]]
    return start, isotope, starts, ends


def create(img, sif_path, pet_info, timing, isotope, diary_name, sing_turku, dicom_start=None, study='EDNiX'):
    """
    SIF of a dynamic image: frame times from the BIDS sidecar, counts estimated by imgweigh (weights of
    the kinetic fits). Without the Turku container (sing_turku is None) only the frame times are written.
    """
    nl = 'create the sif file ' + opb(sif_path)
    run_cmd.msg(nl, diary_name, 'OKGREEN')
    if not opi(img):
        raise FileNotFoundError(f"Input file not found: {img}")
    if timing is None:
        raise ValueError(f'{img}: a SIF needs the frame timing (FrameTimesStart / FrameDuration)')
    write(sif_path, start_string(pet_info, dicom_start), timing, isotope, study)
    if sing_turku is not None:
        turkusif.weight(img, sif_path, diary_name, sing_turku)
    return sif_path


def dicom_info(dir_dcm, endfile):
    """(isotope, start) from the first DICOM file of dir_dcm (pattern *.<endfile>), or (None, None)."""
    from modalities.pet import dcmFile
    if not dir_dcm or not ope(dir_dcm):
        return None, None
    files = sorted(glob.glob(opj(dir_dcm, '*' + (('.' + endfile) if endfile else ''))))
    files = [f for f in files if opi(f)]
    if not files:
        return None, None
    return dcmFile.extract(files[0])


def copy(sif, img):
    """Copy a SIF next to `img` with the matching name (<img without .nii.gz>.sif)."""
    if sif and opi(sif):
        dst = img[:-7] + '.sif' if img.endswith('.nii.gz') else img + '.sif'
        if os.path.abspath(sif) != os.path.abspath(dst):
            shutil.copyfile(sif, dst)
        return dst
    return ''
