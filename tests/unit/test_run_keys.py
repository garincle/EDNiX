"""Run id rule shared by profiles, reliability and group_maps (run_keys.py)."""
import pytest

from Statistics.Group_fMRI.run_selection import run_keys as rk

CASES = [
    # unchanged ids (every paper dataset without acq-/task- variants)
    ("sub-01_ses-1_task-rest_run-3_bold.nii.gz", "3"),
    ("sub-01_ses-1_task-rest_run-01_bold.nii.gz", "01"),
    ("sub-01_task-rest_run-RL-4_bold.nii.gz", "RL-4"),
    ("sub-01_task-rest_run-LR-2_space-template_desc-fMRI_residual.nii.gz", "LR-2"),
    ("sub-01_ses-1_task-rest_bold.nii.gz", "NA"),
    # Babayan: acq-AP / acq-PA with the same run number must differ
    ("sub-01_ses-1_task-rest_acq-AP_run-01_bold.nii.gz", "acqAP-01"),
    ("sub-01_ses-1_task-rest_acq-PA_run-01_bold.nii.gz", "acqPA-01"),
    # Hori: restawake / restanest with the same run number must differ
    ("sub-01_ses-1_task-restawake_run-1_bold.nii.gz", "restawake-1"),
    ("sub-01_ses-1_task-restanest_run-1_bold.nii.gz", "restanest-1"),
    ("sub-01_ses-1_task-restawake_acq-AP_run-2_bold.nii.gz", "restawake-acqAP-2"),
    # variant label without a run token
    ("sub-01_ses-1_task-rest_acq-AP_bold.nii.gz", "acqAP"),
]


@pytest.mark.parametrize("name, expected", CASES)
def test_run_id(name, expected):
    assert rk.run_id(name) == expected
    assert rk.run_id("/some/dir/" + name) == expected


def test_babayan_hori_keys_are_distinct():
    pairs = [("acq-AP", "acq-PA"), ("task-restawake", "task-restanest")]
    for a, b in pairs:
        fa = f"sub-01_ses-1_{a if 'task' in a else 'task-rest_' + a}_run-01_bold.nii.gz"
        fb = f"sub-01_ses-1_{b if 'task' in b else 'task-rest_' + b}_run-01_bold.nii.gz"
        assert rk.run_id(fa) != rk.run_id(fb)


def test_with_variants_keeps_legacy_core():
    assert rk.with_variants("sub-1_task-rest_run-1_bold.nii.gz", "01") == "01"
    assert rk.with_variants("sub-1_task-rest_acq-AP_run-1_bold.nii.gz", "01") == "acqAP-01"
    assert rk.with_variants("sub-1_task-rest_acq-AP_bold.nii.gz", None) == "acqAP"


def test_profiles_and_reliability_use_the_same_rule():
    profiles = pytest.importorskip("Statistics.Group_fMRI.run_selection.profiles")
    reliability = pytest.importorskip("Statistics.Group_fMRI.run_selection.reliability")
    for name, expected in CASES:
        assert profiles._extract_run_id(name) == expected
        assert reliability.run_token(name) == expected
