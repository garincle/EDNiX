import os
opb = os.path.basename
def extract_filename(file_path):
    # Use os.path.basename to get the base name of the file
    base_name = opb(file_path)
    # Check if the file has a ".gz" extension and remove it if present
    if base_name.endswith(".nii.gz"):
        base_name = base_name[:-7]
    # Return the desired part of the file path
    elif base_name.endswith(".nii"):
        base_name = base_name[:-4]
    # Return the desired part of the file path
    return base_name


def get_desc(s, spec):
    """
    Split a BIDS-like file name at its first `<spec>-...` entity.

    get_desc('sub-01_ses-1_trc-FDG_desc-reorient_pet.nii.gz', 'desc')
        -> ('sub-01_ses-1_trc-FDG', 'desc-reorient')
    When the entity is absent, the name without its suffix (last '_' field) is returned
    with an empty string: ('sub-01_ses-1_trc-FDG', '') for 'sub-01_ses-1_trc-FDG_pet.nii.gz'.
    (Suggested by the PET engineer; the entity is matched on its key, not on a substring.)
    """
    parts = extract_filename(s).split('_')
    for i, elem in enumerate(parts):
        if elem.split('-')[0] == spec:
            return '_'.join(parts[:i]), elem
    return '_'.join(parts[:-1]), ''
