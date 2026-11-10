"""
Reset of the anatomical outputs of a subject / session before relaunching EDNiX.anat, keeping every manual file.
Same behaviour as modalities/fMRI/_000_Data_Reset.py (dry run first, one 'YES' for a whole batch).

KEPT
  - raw BIDS files directly in anat/ (names starting with 'sub-'),
  - masks and segmentations corrected by hand, wherever they are under anat/ (MANUAL_PATTERNS):
      <ID>_final_mask.nii.gz          itk_1 (read back by step 2)
      <ID>_final_mask_2.nii.gz        itk_2 (read back by step 5), and its _pre_itk2_backup
      <ID>_final_mask_orig.nii.gz     older name, kept by Tools/Clean_BIDS.py
      <ID>_backtonative_mask.nii.gz   itk_3 (read back by step 8)
      *_desc-manual*                  masks / aseg drawn with the 'manual' methods (skullstrip/manualmask.py,
                                      loop3/manualaseg.py)
      segmentation.mgz / .nii.gz      tissue classification corrected in freeview (read by FS_2_WB)
      *.scene                         Workbench scenes made by hand (e.g. Exemple1.scene for the surface QC)
REMOVED
  - every other file written under anat/ (preprocessing/, matrices/, native/, freesurfer/, templates/, QC_anat/,
    diaries and launcher_parameters files), then the empty folders.

scope='surfaces' removes only the surface outputs (anat/freesurfer/ and anat/native/surfaces/, manual files still
kept): use it to rebuild the surfaces (steps 10-15) with new parameters without redoing the registration.

The dataset-level folders (<BIDS>/sty_template, <BIDS>/QC) are not touched.
"""
import os

from Tools import data_reset as dr

opj = os.path.join
opb = os.path.basename

MANUAL_PATTERNS = ('*_final_mask*', '*_backtonative_mask*', '*_desc-manual*', 'segmentation.mgz',
                   'segmentation.nii.gz', '*.scene')
SURFACE_DIRS = (opj('freesurfer',), opj('native', 'surfaces'))


def plan_anat(data_path, scope='all'):
    """(root, to_delete, kept) for one subject / session. scope: 'all' or 'surfaces'."""
    if scope not in ('all', 'surfaces'):
        raise ValueError(f"scope={scope!r}: expected 'all' or 'surfaces'")
    root = opj(data_path, 'anat')

    def keep_file(p):
        if dr.matches(opb(p), MANUAL_PATTERNS):
            return True
        if scope == 'surfaces':
            rel = os.path.relpath(p, root)
            return not any(rel.startswith(d + os.sep) for d in SURFACE_DIRS)
        return False

    def keep_loose(name):
        return True if scope == 'surfaces' else dr.is_raw_bids(name)

    to_delete, kept = dr.collect(root, keep_loose, keep_file)
    return root, to_delete, kept


def reset_anat(data_path, diary_file, dry_run=False, skip_confirmation=False, scope='all'):
    """Reset one subject / session (data_path = the folder that contains anat/)."""
    root, to_delete, kept = plan_anat(data_path, scope)
    return dr.execute('anat._000_Data_Reset.reset_anat (scope=' + scope + ')', root, to_delete, kept, diary_file,
                      dry_run=dry_run, skip_confirmation=skip_confirmation)


def reset_all(allinfo_study_c, bids_dir, list_to_keep=None, list_to_remove=None, dry_run=True, scope='all'):
    """
    reset_anat on every subject / session selected as in the launchers (pass the same allinfo_study_c,
    list_to_keep and list_to_remove as your preprocess_anat call). dry_run=True (default) only previews;
    dry_run=False previews every session then asks ONE 'YES' for the batch.
    """
    plans = []
    for ID, Session, data_path in dr.subjects(allinfo_study_c, bids_dir, list_to_keep, list_to_remove):
        root, to_delete, kept = plan_anat(data_path, scope)
        plans.append(dict(label=f'anat reset ({scope}): sub-{ID} ses-{Session}', root=root, to_delete=to_delete,
                          kept=kept, diary_file=opj(data_path, f'{ID}_session_{Session}_reset_anat_diary.txt')))
    return dr.run_batch(plans, dry_run=dry_run)
