"""
Reset of the PET (or CT) outputs of a subject / session before relaunching EDNiX.pet / EDNiX.ct, keeping the raw
BIDS data and what was drawn by hand. Same behaviour as modalities/fMRI/_000_Data_Reset.py (dry run first, one 'YES'
for a whole batch).

KEPT
  - raw BIDS files directly in pet/ (ct/): names starting with 'sub-' (images, sidecars, blood / recording tables),
  - the crops drawn in freeview (steps 3-6, crop_for_MoCo='manual' or 'auto' without study template):
      orig/preprocessing/forMC/<session>_desc-crop_pet.nii.gz (+ .json)      (MANUAL_PATTERNS)
  - the masks drawn by hand: raw/masks/ (whole folder) and any *_desc-manual_* file (as func/raw/masks/)
  - keep_registration=True also keeps orig/matrices/ (registrations, including a manual correction made with
    do_manual=True). Only valid if step 1 (reorientation) is not run again: the matrices are defined on the
    reoriented image.
REMOVED
  - every other file under pet/ (ct/): raw/, orig/, acpc/, templates/, QC/, diaries, launcher_parameters and
    Prepro_Settings.json (rewritten at every run), then the empty folders.

The PET study template and its manual crop (<BIDS>/sty_template/pet) and <BIDS>/QC are not touched.
"""
import os

from Tools import data_reset as dr

opj = os.path.join
opb = os.path.basename

MANUAL_PATTERNS = ('*_desc-crop_*', '*_desc-manual_*')


def plan_emission(data_path, modality='pet', keep_registration=False):
    """(root, to_delete, kept) for one subject / session; modality 'pet' or 'ct'."""
    if modality not in ('pet', 'ct'):
        raise ValueError(f"modality={modality!r}: expected 'pet' or 'ct'")
    root = opj(data_path, modality)
    matrices = opj(root, 'orig', 'matrices') + os.sep
    manual_masks = opj(root, 'raw', 'masks') + os.sep

    def keep_file(p):
        if dr.matches(opb(p), MANUAL_PATTERNS) or p.startswith(manual_masks):
            return True
        return keep_registration and p.startswith(matrices)

    to_delete, kept = dr.collect(root, dr.is_raw_bids, keep_file)
    return root, to_delete, kept


def reset_pet(data_path, diary_file, dry_run=False, skip_confirmation=False, keep_registration=False,
              modality='pet'):
    """Reset one subject / session (data_path = the folder that contains pet/ or ct/)."""
    root, to_delete, kept = plan_emission(data_path, modality, keep_registration)
    return dr.execute(f'{modality}._000_Data_Reset.reset (keep_registration={keep_registration})', root, to_delete,
                      kept, diary_file, dry_run=dry_run, skip_confirmation=skip_confirmation)


def reset_all(allinfo_study_c, bids_dir, list_to_keep=None, list_to_remove=None, dry_run=True,
              keep_registration=False, modality='pet'):
    """
    Reset every subject / session selected as in the launchers (pass the allinfo_study_c built with
    modalities=['pet'] (or ['ct']) and the same list_to_keep / list_to_remove as the preprocess_pet call).
    dry_run=True (default) only previews; dry_run=False previews every session then asks ONE 'YES'.
    """
    plans = []
    for ID, Session, data_path in dr.subjects(allinfo_study_c, bids_dir, list_to_keep, list_to_remove):
        root, to_delete, kept = plan_emission(data_path, modality, keep_registration)
        plans.append(dict(label=f'{modality} reset: sub-{ID} ses-{Session}', root=root, to_delete=to_delete,
                          kept=kept,
                          diary_file=opj(data_path, f'{ID}_session_{Session}_reset_{modality}_diary.txt')))
    return dr.run_batch(plans, dry_run=dry_run)
