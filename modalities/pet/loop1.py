"""
Block 1 of the PET pipeline for one subject / session (step 1):
frame timing, mean image and SIF of every raw run (pet/raw), then reorientation (pet/orig/preprocessing).
"""
import os
import json

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists
opi = os.path.isfile

from Tools import run_cmd
from modalities.pet import petutils as pu
from modalities.pet import reorient, sifFile


def position(pos, i):
    """Position of run i: launcher lists ([''] or one value per run) or a string."""
    if isinstance(pos, (list, tuple)):
        if not pos:
            return ''
        return pos[i] if i < len(pos) else pos[0]
    return pos or ''


def timing_file(paths, raw_img):
    return pu.out_name(paths['raw'], pu.base_name(raw_img), desc='frametiming', ext='.json')


def sif_file(paths, raw_img):
    return opj(paths['raw'], pu.strip_ext(raw_img) + '.sif')


def run(list_img, paths, humanPosition, animalPosition, orientation, check_ref, header_ref, Need_check_orient,
        PET_header_wrong, deoblique, isotope, srtm, dcm, Skip_step, diary_file, sing_afni, sing_fs, sing_turku,
        follow=None):
    """
    list_img: raw PET runs of the session. check_ref: image shown under the PET for the orientation check;
    header_ref: raw image used when PET_header_wrong=True; dcm: (folder, extension) of the DICOM files.
    follow: world transform of the reorientation of the anatomy / CT of the same scanner space (reorient.doit).
    Returns the reoriented images (also when step 1 is skipped, if they exist).
    """
    for d in (paths['path'], paths['raw'], paths['orig'], paths['orig_process']):
        os.makedirs(d, exist_ok=True)
    reoriented = []
    for i, img in enumerate(list_img):
        info = pu.raw_sidecar(img)
        nfr = pu.n_frames(img)
        timing = pu.frame_timing(info, nfr)
        out = pu.out_name(paths['orig_process'], pu.base_name(img), desc='reorient')
        reoriented.append(out)
        if 1 in Skip_step:
            continue
        run_cmd.msg('INFO: ' + opb(img) + f': {nfr} frame(s)', diary_file, 'OKGREEN')

        # frame timing (pet2b: raw/<name>_desc-frametiming_pet.json)
        if timing is not None:
            with open(timing_file(paths, img), 'w') as f:
                json.dump(pu.timing_to_json(timing), f, indent=2)

        # mean image of the raw run
        if nfr > 1:
            pu.write_mean(img, pu.out_name(paths['raw'], pu.base_name(img), desc='mean'), timing, None,
                          'Duration-weighted mean of the raw frames.')

        # SIF (frame times and counts, Turku tools), only for the kinetic models
        sif = ''
        if srtm and nfr > 1:
            dcm_iso, dcm_start = sifFile.dicom_info(*dcm) if dcm and dcm[0] else (None, None)
            if dcm_iso and isotope != dcm_iso:
                run_cmd.msg(f'WARNING: isotope {isotope} (launcher / sidecar) but {dcm_iso} in the DICOM half-life',
                            diary_file, 'WARNING')
            sif = sifFile.create(img, sif_file(paths, img), info, timing, isotope, diary_file, sing_turku,
                                 dicom_start=dcm_start)

        hp = position(humanPosition, i) or info.get('PatientPosition', '')
        ap = position(animalPosition, i)
        if follow is None and not reorient.no_reorientation(ap) and not hp and not orientation:
            raise ValueError(f'{opb(img)}: animalPosition={ap!r} needs the scanner position (humanPosition or '
                             'PatientPosition in the sidecar), or give the orientation code')
        reorient.doit(paths['orig_process'], img, hp, ap, orientation, check_ref, Need_check_orient,
                      PET_header_wrong, timing, diary_file, sing_afni, sing_fs, suffix='pet', deoblique=deoblique,
                      sif=sif, header_ref=header_ref, follow=follow)
    return reoriented
