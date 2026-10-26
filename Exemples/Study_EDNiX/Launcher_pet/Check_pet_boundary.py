import glob
import os
import sys
opj = os.path.join

MAIN_PATH = opj('/home/cgarin/PycharmProjects/EDNiX/')
sys.path.insert(1, opj(MAIN_PATH))

from modalities.pet import boundary

######################################################################################################################
# Boundary check of the PET -> T1w registrations already done (step 6 of Launcher_fdg_Macaque_imagina.py: PET -> T1w;
# Launcher_fdg_CERMEP.py: PET -> CT -> T1w), without running the pipeline again. For every session: is the dim outer
# cortex partial volume (PET boundary on the MRI boundary) or a registration error (PET brain larger / smaller /
# shifted)? Results next to each registration record (<session>_desc-registration_boundary.json / .png) and in
# <BIDS>/QC/pet_boundary/. Step 6 now runs it by itself.
######################################################################################################################
# dataset: PET resolution (FWHM, mm) of the scanner + reconstruction; None: twice the largest PET voxel
datasets = {opj('/scratch2/EDNiX/Macaque/BIDS_Tremblay/'): None,
            opj('/scratch2/EDNiX/Human/BIDS_Merida/'): None}


def report(rec, pet_fwhm, bids_dir):
    try:
        res = boundary.check_record(rec, pet_fwhm, bids_dir)
    except Exception as e:
        print(f'{rec}: not checked ({e})')
        return
    sh = ', '.join(f'{k} {v:+.1f}' for k, v in res['shifts_mm'].items())
    print(f"{os.path.basename(rec)}\n    PET edge - MRI edge (mm): {sh}\n    clean sectors (no extracerebral uptake "
          f"beside): {res['clean_sectors']}; dilation {res['dilation_mm']} mm, translation {res['translation_mm']} mm "
          f"(FWHM {res['fwhm_mm']:.1f} mm, model {res['model']})\n    -> {res['verdict']}")


def records_of(bids_dir):
    return sorted(glob.glob(opj(bids_dir, 'sub-*', '**', 'pet', 'orig', 'matrices', '*_desc-registration_pet.json'),
                            recursive=True))


for bids_dir, pet_fwhm in datasets.items():
    records = records_of(bids_dir)
    print(f'{bids_dir}: {len(records)} registration records')
    for rec in records:
        report(rec, pet_fwhm, bids_dir)

