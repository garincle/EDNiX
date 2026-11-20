import pytest

from modalities.anat.skullstrip import Unetmask


def test_missing_unet_model_stops_with_its_path(tmp_path):
    """A U-Net method whose model file is absent stops before running muSkullStrip.py."""
    unetpath = tmp_path / 'Singularity'
    (unetpath / 'NHP-BrainExtraction' / 'UNet_Model' / 'models').mkdir(parents=True)
    diary = tmp_path / 'diary.txt'
    with pytest.raises(Exception, match='model-legacy-ep9'):
        Unetmask.do(str(tmp_path / 'sub-01_T1w.nii.gz'), str(tmp_path / 'mask.nii.gz'),
                    'model-legacy-ep9', str(unetpath), ' -fill_holes', '', str(diary))
    assert 'U-Net model not found' in diary.read_text()


def test_v2_methods_are_listed():
    from modalities.anat.skullstrip.Skullstrip_method import SSmethodslist
    assert {'muSkullStrip_cross_speciesV2', 'muSkullStrip_cross_speciesV2_dil'} <= set(SSmethodslist)
