"""
Step 1 of the CT pipeline: orientation of the header, exactly as the PET (modalities/pet/reorient.py,
shared code), so that a PET and a CT of the same PET/CT acquisition stay in the same space. The PET
engineer's copy of the PET module is replaced by this call.
"""
from modalities.pet import reorient as _pet_reorient


def doit(new_dir, img, humanPosition, animalPosition, orientation, Ref_native_raw, Need_check, header_wrong,
         diary_name, sing_afni, sing_fs, deoblique='WARP', header_ref='', follow=None):
    """See modalities/pet/reorient.doit; follow: world transform of the anatomy of the same scanner space."""
    return _pet_reorient.doit(new_dir, img, humanPosition, animalPosition, orientation, Ref_native_raw, Need_check,
                              header_wrong, None, diary_name, sing_afni, sing_fs, suffix='ct', deoblique=deoblique,
                              header_ref=header_ref, follow=follow)
