"""
Reset of the CT outputs before relaunching EDNiX.ct: same rules as the PET reset (modalities/pet/_000_Data_Reset.py),
applied to the ct/ folder.
"""
from modalities.pet import _000_Data_Reset as _pet_reset


def reset_ct(data_path, diary_file, dry_run=False, skip_confirmation=False, keep_registration=False):
    return _pet_reset.reset_pet(data_path, diary_file, dry_run=dry_run, skip_confirmation=skip_confirmation,
                                keep_registration=keep_registration, modality='ct')


def reset_all(allinfo_study_c, bids_dir, list_to_keep=None, list_to_remove=None, dry_run=True,
              keep_registration=False):
    return _pet_reset.reset_all(allinfo_study_c, bids_dir, list_to_keep, list_to_remove, dry_run=dry_run,
                                keep_registration=keep_registration, modality='ct')
