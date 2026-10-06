# import
import os
import subprocess
import glob
import json
import nibabel as nib
import numpy as np
import pandas as pd
from modalities.anat import set_launcher
from atlases import atlas4func
import inspect

opj = os.path.join
opb = os.path.basename
opr = os.path.relpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
opi = os.path.isfile

spgo = subprocess.getoutput

from Tools import run_cmd, getpath, diaryfile, Load_EDNiX_requirement, check_nii
from modalities.anat.loop3 import backtonative
from modalities.fMRI import chooseanat, _1_fMRI_preTTT_in_fMRIspace, _2_coregistration_to_norm, _3_mask_fMRI, _4_check_mask, \
    _6_Melodic, _8_fMRI_to_anat
from modalities.fMRI import _7_post_TTT, _11_Seed_base_many_regionsatlas, _5_anat_to_fMRI, _200_Data_QC
from modalities.fMRI import _9_coregistration_to_template_space, _10_Correl_matrix, _12_fMRI_QC, _14_fMRI_QC_matrix, _100_Data_Clean, \
    itk_2
from modalities.fMRI.extract_filename import extract_filename

EES_MIN, EES_MAX = 5e-5, 5e-3     # s: effective echo spacing of an EPI (50 us - 5 ms)
TRT_MIN, TRT_MAX = 2e-3, 0.3       # s: total readout time of an EPI (2 - 300 ms)


class SliceTimingError(Exception):
    """SliceTiming impossible for this TR (kept apart: the reading of the sidecar below has a bare except)."""


def check_slice_timing(slice_timing, tr):
    """
    (kind, message) about the unit of a SliceTiming list (seconds, offset of each slice within its volume):
      ('error', ...)   a slice time reaches the TR or more: impossible (milliseconds?);
      ('warning', ...) the slices span less than 1 % of the TR: slice-timing correction would be negligible
                       (values 1000 times too small); the message gives the list x 1000 when it fits the TR;
      (None, '')       plausible, or all slices simultaneous.
    """
    st = np.asarray(slice_timing, dtype=float)
    if st.size == 0 or not np.all(np.isfinite(st)):
        return 'error', 'SliceTiming is empty or not numeric.'
    top, span = float(st.max()), float(st.max() - st.min())
    if top >= tr * 1.001:
        return 'error', (f'SliceTiming goes up to {top:g} s, beyond the TR ({tr:g} s): a slice cannot be acquired after the '
                         f'end of its volume. The values are probably in milliseconds: give Slice_timing_info in seconds '
                         f'in the launcher.')
    if 0 < span < 0.01 * tr:
        msg = (f'WARNING: the slices of SliceTiming span {span:g} s, {100 * span / tr:.3f} % of the TR ({tr:g} s): slice-timing '
               f'correction shifts the slices by less than that, i.e. it is negligible. ')
        if 0.3 * tr <= span * 1000 < tr:
            msg += (f'Multiplied by 1000 the slices span {span * 1000:g} s ({100 * span * 1000 / tr:.0f} % of the TR, as a '
                    f'continuous acquisition does): the sidecar is probably in the wrong unit. If so, give in the launcher '
                    f'Slice_timing_info={[round(float(v) * 1000, 6) for v in st]}')
        else:
            msg += 'Check the unit of SliceTiming and give Slice_timing_info (seconds) in the launcher.'
        return 'warning', msg
    return None, ''


def resolve_TRT(TRT, info_RS):
    """Total readout time (s) for TOPUP, as text, or 'None' (with a reason) when absent or implausible."""
    if str(TRT) != 'Auto':
        return str(TRT), f'TRT = {TRT} s (given in the launcher)', 'OKGREEN'
    v = (info_RS or {}).get('TotalReadoutTime')
    if v is None:
        return 'None', 'INFO: Total Readout Time not found in header', 'WARNING'
    v = float(v)
    if not TRT_MIN <= v <= TRT_MAX:
        return 'None', (f'WARNING: TotalReadoutTime = {v:g} s of the sidecar is outside the range of an EPI readout '
                        f'({TRT_MIN * 1000:g}-{TRT_MAX * 1000:g} ms): not used, so TOPUP is skipped. Wrong unit in the '
                        f'sidecar? Give TRT (seconds) in the launcher.'), 'WARNING'
    return str(v), f'INFO: Total Readout Time = {v} s', 'OKGREEN'


def resolve_dwell(DwellT, info_RS):
    """
    Effective echo spacing (seconds per phase-encoding line), the value FSL FUGUE expects
    for --dwell. Returns (value as text or 'None', log message, message level).

    'Auto' reads, in this order:
      1. EffectiveEchoSpacing (BIDS);
      2. TotalReadoutTime / (ReconMatrixPE - 1) (BIDS);
      3. DwellTime, only if its value is plausible for an echo spacing (0.05-5 ms): in BIDS,
         DwellTime is the receiver dwell time per readout sample (microseconds), a different
         quantity that would make the FUGUE correction a hundred times too small.
    The previous version read DwellTime first, then divided a 'TotalReadOutTimeEPI' by the
    number of slices instead of the number of phase-encoding lines.
    """
    if str(DwellT) != 'Auto':
        return str(DwellT), f'DwellT = {DwellT} s (given in the launcher)', 'OKGREEN'
    info = info_RS or {}
    unusable = ("WARNING: {what} = {v:g} s is outside the range of an effective echo spacing "
                f"({EES_MIN * 1000:g}-{EES_MAX * 1000:g} ms): not used, so no distortion correction. Wrong unit in the "
                "sidecar? Give DwellT (seconds) in the launcher.")
    if info.get('EffectiveEchoSpacing') is not None:
        v = float(info['EffectiveEchoSpacing'])
        if not EES_MIN <= v <= EES_MAX:
            return 'None', unusable.format(what='EffectiveEchoSpacing', v=v), 'WARNING'
        return '%.16f' % v, f'DwellT = {v} s (EffectiveEchoSpacing of the JSON sidecar)', 'OKGREEN'
    if info.get('TotalReadoutTime') is not None and info.get('ReconMatrixPE'):
        v = float(info['TotalReadoutTime']) / (float(info['ReconMatrixPE']) - 1.0)
        if not EES_MIN <= v <= EES_MAX:
            return 'None', unusable.format(what='TotalReadoutTime / (ReconMatrixPE - 1)', v=v), 'WARNING'
        return '%.16f' % v, (f'DwellT = {v:.6g} s (TotalReadoutTime / (ReconMatrixPE - 1) of the JSON sidecar)'), 'OKGREEN'
    if info.get('DwellTime') is not None:
        v = float(info['DwellTime'])
        if 5e-5 <= v <= 5e-3:
            return '%.16f' % v, (f'WARNING: DwellT = {v} s taken from DwellTime, read as an effective echo spacing '
                                 '(non-standard sidecar): CHECK IT.'), 'WARNING'
        return 'None', (f'WARNING: DwellTime = {v} s is a receiver sample dwell time, not an echo spacing: no '
                        'distortion correction. Give DwellT (effective echo spacing, s) in the launcher.'), 'WARNING'
    return 'None', ('WARNING: effective echo spacing not found (EffectiveEchoSpacing, TotalReadoutTime and '
                    'ReconMatrixPE): no distortion correction. Give DwellT in the launcher.'), 'WARNING'


TR_MIN, TR_MAX = 0.05, 30.0       # s: plausible repetition times of BOLD acquisitions


def _check_plausible_tr(tr, source):
    """Refuse a TR that cannot be a BOLD repetition time, suggesting the unit error that explains it."""
    if TR_MIN <= tr <= TR_MAX:
        return
    hint = ''
    for factor, what in ((1000.0, 'x 1000'), (1e-3, '/ 1000')):
        if TR_MIN <= tr * factor <= TR_MAX:
            hint = f' ({what} would give {tr * factor:g} s: wrong unit in the sidecar?)'
    raise ValueError(f'TR = {tr:g} s ({source}) is not a plausible repetition time for BOLD fMRI '
                     f'(expected {TR_MIN:g} to {TR_MAX:g} s){hint}. Give the TR of the acquisition protocol, in seconds, '
                     f"in the launcher (for example TR='1.0') instead of 'Auto'.")


def resolve_TR(TR, info_RS):
    """
    Repetition time in seconds, as a float.

    TR     : 'Auto' or a value given in the launcher (number or text such as '2.0').
    info_RS: the JSON sidecar of the run (dict), or None.
    Returns (TR in s, log message, message level); raises ValueError when no plausible TR can be found.

    A TR given in the launcher is used as it is (a warning if it is outside 0.05-30 s).
    'Auto' reads RepetitionTime, else the spacing of VolumeTiming (BIDS, seconds, regular only). If both
    are missing, the TR is estimated from SliceTiming (seconds in
    BIDS) as the number of distinct slice times multiplied by their spacing: this assumes no gap
    between volumes and the right unit, so the estimate is refused when it is not plausible
    (a sidecar whose SliceTiming is not in seconds gives a TR wrong by a factor 1000), and a
    warning asks to check it. Version 5 multiplied the estimate by 1000 a priori; that is wrong for
    BIDS sidecars but compensated a SliceTiming 1000 times too small, which it then turned into
    a TR 1000 times too small: no unit is assumed any more.
    """
    if str(TR) != 'Auto':
        try:
            tr = float(TR)
        except (TypeError, ValueError):
            raise ValueError(f'TR = {TR!r} is not a number (seconds).')
        if tr <= 0:
            raise ValueError(f'TR = {tr:g} s: the repetition time must be positive.')
        if not TR_MIN <= tr <= TR_MAX:
            return tr, (f'WARNING: TR = {tr:g} s given in the launcher is outside the usual range of BOLD fMRI '
                        f'({TR_MIN:g}-{TR_MAX:g} s): check the unit (seconds).'), 'WARNING'
        return tr, f'TR = {tr} s (given in the launcher)', 'OKGREEN'
    if info_RS and info_RS.get('RepetitionTime') is not None:
        tr = float(info_RS['RepetitionTime'])
        _check_plausible_tr(tr, 'RepetitionTime of the JSON sidecar')
        return tr, f'TR = {tr} s (RepetitionTime of the JSON sidecar)', 'OKGREEN'
    vt = (info_RS or {}).get('VolumeTiming')           # BIDS: onset of each volume (s), instead of RepetitionTime
    if vt is not None and len(vt) > 2:
        steps = np.diff(np.asarray(vt, dtype=float))
        tr = float(np.median(steps))
        if not np.allclose(steps, tr, rtol=1e-3, atol=1e-6):
            raise ValueError("VolumeTiming is not regular (sparse or clustered acquisition): there is no single TR. "
                             "Give the TR of the acquisition protocol, in seconds, in the launcher.")
        _check_plausible_tr(tr, 'spacing of VolumeTiming in the JSON sidecar')
        return tr, (f'TR = {tr} s (spacing of VolumeTiming in the JSON sidecar, which has no RepetitionTime)'), 'OKGREEN'
    st = (info_RS or {}).get('SliceTiming')
    if st:
        times = np.unique(np.round(np.asarray(st, dtype=float), 9))
        if len(times) > 1:
            tr = float(len(times) * np.median(np.diff(times)))
            _check_plausible_tr(tr, 'estimated from the SliceTiming of the JSON sidecar, which has no RepetitionTime')
            return tr, (f'WARNING: RepetitionTime missing from the JSON sidecar; TR estimated from SliceTiming: '
                        f'{tr:.4f} s (assumes no gap between volumes and SliceTiming in seconds). CHECK IT, or give '
                        f'TR in the launcher.'), 'WARNING'
    raise ValueError("TR is 'Auto' but the JSON sidecar has neither RepetitionTime nor usable SliceTiming: "
                     "give the TR in seconds in the launcher (for example TR='2.0').")


def preprocess_fmri(Skip_step, MAIN_PATH, bids_dir,
                    species, allinfo_study_c, endfmri, endjson, endmap, resting_or_task,
                    animalPosition, humanPosition, orientation,
                    Slice_timing_info,
                    TfMRI, type_norm, creat_study_template,
                    anat_func_same_space, coregistration_longitudinal,
                    Method_mask_func, extra_erode=0, do_anat_to_func=True, folderforTemplate_Anat='', IhaveanANAT=True,
                    ntimepoint_treshold=100, REF_int=0, T1_eq=5, correction_direction='Auto', overwrite_option=True,
                    DwellT='Auto', SED='Auto', TR='Auto', TRT='Auto',
                    nb_ICA_run=20, ICA_cleaning='Skip',
                    costAllin='lpa',
                    doWARPonfunc=False, registration_type='origin', type_of_transform='BOLDAffine', n_for_ANTS='Lanczos', aff_metric_ants='meansquares', aff_metric_ants_Transl='mattes', dilate_mask=0,
                    list_to_keep=[], list_to_remove=[], atlas_followers=[[], [], [], []],
                    reference='EDNiX', post_treatment_method='Grandjean',
                    band='0.01 0.1', blur=0, do_not_correct_signal = False, extract_WM=True, extract_Vc = False, extract_GS = False,
                    use_erode_WM_func_masks = True, use_erode_V_func_masks=True, normalize='Skip',
                    acompcor_mode='fixed', acompcor_n_fixed=5, acompcor_variance=0.50, acompcor_max=15,
                    motion_mask='None', motion_mask_erode=3, motion_method='ants',
                    selected_atlases_matrix='all', wanted_level_matrix='all', matrix_fit=['covariance', 'correlation', 'partial correlation'],
                    selected_atlases_SBA='default', panda_files_SBA='default',
                    SBAspace=['func', 'anat', 'atlas'], erod_seed=True, smoothSBA=False,
                    intra_thresh=0.240, delta_thresh=0.05,
                    oversample_map=False, use_cortical_mask_func=False, n_cut=10, threshold_val=10,doBALSAspace=False, **kwargs):

    """
    Preprocess functional MRI data within a BIDS-compliant study, integrating anat
    preprocessing, motion and ICA-based denoising, temporal and spatial correction, and
    optional surface-based analyses (SBA). Supports both resting-state and task fMRI, human
    and non-human species, and longitudinal designs.

    The pipeline is modular: steps can be skipped or selectively executed, with full control
    over registration, masking, smoothing, temporal filtering, and atlas transformations.

    Parameters
    ----------
    Skip_step : list
        Steps to skip in the pipeline (integers or descriptive labels).
    MAIN_PATH : str
        Base path for output and intermediate results.
    bids_dir : str
        Path to BIDS-formatted dataset.
    species : str
        Species label ('Human', 'Rat', 'Mouse', etc.).
    allinfo_study_c : pandas.DataFrame
        DataFrame containing BIDS subject/session information.
    endfmri, endjson, endmap : str
        File suffixes for fMRI, JSON, and derived maps.
    resting_or_task : str
        Type of fMRI acquisition ('resting' or 'task').
    animalPosition, humanPosition : list
        Orientation information for animals or humans during acquisition.
    orientation : str
        Target orientation for preprocessing (e.g., 'LPI').
    Slice_timing_info : dict or DataFrame
        Slice acquisition timing information.
    TfMRI : str
        Suffix of the anat use to register functional MRI data ('T2w', etc.).
    type_norm : str
        Anatomical type_norm contrast used originally for normalization ('T1w' or 'T2w', should be the same as in type_norm anat).
    creat_study_template : bool
        Whether  study-specific anat template has been generated in the anat regitration (should be the same as in creat_study_template anat)
    anat_func_same_space : bool
        Flag to indicate if anat and functional images share EXACTLY the same space.
    coregistration_longitudinal : bool
        Whether to apply longitudinal coregistration across sessions. (should be the same as in coregistration_longitudinal anat)
    Method_mask_func : str
        Method for generating functional masks.
    do_anat_to_func : bool, optional
        Perform anat-to-functional registration (default: True).
    folderforTemplate_Anat : str, optional
        Path to anat templates.
    IhaveanANAT : bool, optional
        Whether anat images are available (default: True).
    ntimepoint_treshold : int, optional
        Minimum number of timepoints required to include a run.
    REF_int : int, optional
        Reference volume index for registration.
    T1_eq : int, optional
        Number of T1-equilibration volumes to remove.
    correction_direction : str, optional
        Phase-encoding (unwarping) direction for FSL FUGUE: 'x', 'x-', 'y', 'y-', 'z' or 'z-';
        'Auto' reads PhaseEncodingDirection from the JSON sidecar. Also restricts the non-linear
        func-to-anat registration (do_anat_to_func) to that axis.
    overwrite_option : bool, optional
        Overwrite existing files (default: True).
    DwellT : str or float, optional
        Effective echo spacing in seconds, used by FUGUE; 'Auto' reads EffectiveEchoSpacing, or
        TotalReadoutTime / (ReconMatrixPE - 1), from the JSON sidecar.
    SED : str, optional
        Slice-encoding direction; 'Auto' reads SliceEncodingDirection from the JSON sidecar
        (logged; the registration restriction uses the phase-encoding direction).
    TR : str or float, optional
        Repetition time in seconds; 'Auto' reads RepetitionTime from the JSON sidecar (or estimates
        it from SliceTiming, with a warning).
    TRT : str or float, optional
        Total readout time in seconds, used by TOPUP; 'Auto' reads TotalReadoutTime from the JSON sidecar.
    nb_ICA_run : int, optional
        Number of ICA components for denoising.
    ICA_cleaning : str, optional
        ICA cleaning method ('Skip', 'MELODIC', etc.).
    costAllin : str, optional
        Registration cost function for 3dAllineate alignment ('lpa', etc.). ## will be removed
    doWARPonfunc : str, optional
        How oblique functional images are de-obliqued: 'WARP' resamples them (AFNI 3dWarp -deoblique),
        'header' keeps the data and only rewrites the header; other values (default False) leave
        oblique images as they are.
    registration_type : str, optional
        Parameter set of the func-to-anat registration when do_anat_to_func is True: 'origin'
        (default), 'fast' (fewer iterations) or 'V2'. The transform itself is type_of_transform
        (a non-linear one, e.g. 'SyN', is needed for a distortion-like correction).
    type_of_transform : str, optional
        Type of transformation for alignment ('BOLDAffine', 'Rigid', 'SyN', etc.).
    n_for_ANTS : str, optional
        Interpolation method for ANTs registration ('Lanczos', 'Linear', etc.).
    aff_metric_ants, aff_metric_ants_Transl : str, optional
        Affine registration metrics for ANTs.
    dilate_mask : int, optional
        Number of voxels to dilate functional masks.
    list_to_keep, list_to_remove : list, optional
        Subject/session tuples to include or exclude.
    atlas_followers : list of lists, optional
        Multi-level atlas selection for processing.
    reference : str, optional
        Reference dataset/session ('EDNiX' by default).
    post_treatment_method : str, optional
        Method for post-processing (default: 'Grandjean').
    band : str, optional
        Temporal filter band (default: '0.01 0.1').
    blur : float, optional
        Gaussian smoothing kernel in mm.
    do_not_correct_signal : bool, optional
        Skip temporal or spatial signal corrections.
    motion_method : str, optional
        Program that estimates the head motion in step 1: 'ants' (default, ANTs rigid), '3dvolreg' (AFNI) or
        'mcflirt' (FSL). The transforms of 3dvolreg and MCFLIRT are converted into the ANTs matrices used
        by step 5 (one interpolation per volume with the registration and distortion correction), after a
        check that the converted matrices reproduce the registered volumes written by the program (the run
        stops otherwise). Needs the python package nitransforms. The brain mask (motion_mask) is used with 'ants' only.
    motion_mask : str, optional
        Brain mask for the motion estimation of step 1. 'None' (default): the whole image is used (ANTs
        rigid, 'BOLDRigid'). 'previous': the functional mask written by step 3 of an earlier execution
        (run the pipeline once, then again from step 1). Or the path of a mask in the frame of the
        reoriented functional series. Needed when tissue that does not move with the brain (skull, muscle,
        skin) dominates the image: the log of step 1 compares the two estimations.
    motion_mask_erode : int, optional
        Voxels removed from the border of the motion mask (default 3): it must stay inside the brain
        whatever the head does, so use at least the largest head displacement in voxels.
    extract_WM, extract_Vc, extract_GS : bool, optional
        Nuisance signals of the white matter (Wmask), the ventricles (Vmask) and the brain (global
        signal). They select the regressors of the 'AFNI' method only. 'Grandjean' always uses the mean
        signal of the white-matter and ventricle masks (when the masks exist) and 'fMRIPrep_style'
        always uses aCompCor of both; extract_GS adds the global signal to 'AFNI' and 'fMRIPrep_style'.
    acompcor_mode : str, optional
        'fixed' (default): acompcor_n_fixed components per tissue (Behzadi et al. 2007); 'parallel':
        the components that stand above chance (parallel analysis against a phase-randomised null),
        at most acompcor_max; 'adaptive': as many components as needed to explain acompcor_variance
        of the variance, at most acompcor_max (the cap is usually reached: it then acts as a fixed
        number). Only for post_treatment_method='fMRIPrep_style'. The log reports, for every tissue,
        the number found by parallel analysis.
    acompcor_n_fixed : int, optional
        Number of aCompCor components per tissue in 'fixed' mode (default 5).
    acompcor_variance : float, optional
        Fraction of the variance to explain in 'adaptive' mode (default 0.50).
    acompcor_max : int, optional
        Cap on the number of aCompCor components per tissue in 'adaptive' and 'parallel' modes (default 15).
    use_erode_WM_func_masks, use_erode_V_func_masks : bool, optional
        Apply erosion to WM and ventricle masks.
    normalize : str, optional
        Normalization method or 'Skip'.
    selected_atlases_matrix, wanted_level_matrix : str, optional
        Atlas selection for matrix-based analyses ('all', 'default', etc.).
    selected_atlases_SBA, panda_files_SBA : str, optional
        Atlas and file settings for seed-based analyses.
    SBAspace : list, optional
        Space(s) for SBA outputs ('func', 'anat', 'atlas').
    erod_seed : bool, optional
        Erode seed masks for SBA.
    smoothSBA : bool, optional
        Smooth SBA results.
    intra_thresh, delta_thresh : float, optional
        Connectivity-specificity thresholds of step 13 (Grandjean et al. 2020): a run is 'Specific'
        when homotopic connectivity exceeds intra_thresh and exceeds the cross-network connectivity
        by more than delta_thresh; 'Unspecific' when it exceeds intra_thresh but not by delta_thresh;
        'Spurious' when only the cross-network connectivity exceeds intra_thresh; 'No' otherwise.
    oversample_map : bool, optional
        Oversample maps to higher resolution.
    use_cortical_mask_func : bool, optional
        Restrict analyses to cortical regions.
    n_cut : int, optional
        Number of voxels to exclude at edges.
    threshold_val : float, optional
        Intensity threshold for masking.
    **kwargs : dict
        Additional parameters to pass to custom functions or overrides.

    Workflow
    --------
    1. Load BIDS dataset and extract fMRI sessions.
    2. Filter subjects and sessions according to include/exclude lists.
    3. Apply slice timing correction, motion correction, and temporal filtering.
    4. Register functional images to anat and template space.
    5. Generate functional masks (WM, CSF, ventricles, brain).
    6. Optional ICA-based denoising and post-processing.
    7. Surface-based analyses if enabled (SBA).
    8. Spatial smoothing and normalization.
    9. Apply atlas transformations and compute derived metrics.
    10. Optional QC checks and overwriting of existing outputs.

    Notes
    -----
    - Supports both task and resting-state fMRI.
    - Modular: steps can be skipped with Skip_step.
    - Compatible with multi-species studies and longitudinal designs.
    - Integrates fully with anat preprocessing to ensure accurate coregistration.
    """
    if kwargs:
        run_cmd.printcolor('WARNING: unknown argument(s) ignored by preprocess_fmri: ' +
                           ', '.join(sorted(kwargs)) + ' — check the spelling in your launcher.', 'WARNING')

    (FS_refs, template_dir, reference, balsa_folder, BALSAname, balsa_brainT1, BASE_atlas_folder, BASE_template, BASE_SS,
    BASE_mask, BASE_Gmask, BASE_WBGmask, BASE_Wmask, BASE_Vmask, CSF, GM, WM, Aseg_ref, list_atlas, path_label_code, all_ID,
    all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max,
    fs_tools, reftemplate_path, MNIBcorrect_indiv, masking_img) = set_launcher.get(MAIN_PATH, bids_dir, allinfo_study_c,
                                                                                   species, list_to_keep,
                                                                                   list_to_remove, reference, type_norm,
                                                                                   '', '', atlas_followers, doBALSAspace)

    sing_afni, sing_fsl, sing_fs, sing_itk, sing_wb, _, sing_synstrip, Unetpath = Load_EDNiX_requirement.load_requirement(
        MAIN_PATH, reftemplate_path, bids_dir, 'yes')

    #######for matrix analysis (step 10)
    #### name of the atlases  you want to use for the matrix analysis
    Lut_dir, selected_atlases_matrix, wanted_level, segmentation_name_list, selected_atlases, panda_files = (
        atlas4func.setup(MAIN_PATH, species, reference, selected_atlases_matrix, wanted_level_matrix,
                         selected_atlases_SBA, panda_files_SBA, template_dir))

    if species in ['Human', 'Chimpanzee']:
        config_f = opj(MAIN_PATH, 'Singularity_library', 'config', 'b02b0Human.cnf')
    else:
        config_f = opj(MAIN_PATH, 'Singularity_library', 'config', 'b02b0Macaque.cnf')

    # Set the environment variable for the current process
    os.environ["AFNI_NIFTI_TYPE_WARN"] = "NO"
    overwrite = ''
    if overwrite_option == True:
        overwrite = ' -overwrite'
    # Usage
    type_norm = check_nii.normalize_anat_type(type_norm)
    TfMRI = check_nii.normalize_anat_type(TfMRI)

    for ID, Session, data_path, max_ses in zip(all_ID, all_Session, all_data_path, all_Session_max):
        nl = 'INFO: Work on ' + str(ID) + ' session ' + str(Session)
        run_cmd.printcolor(nl, 'HEADER')

        # Resting data organization folders ###########################################################################
        (dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_raw_matrices, path_func, dir_fmap, dir_prepro_fmap,
            dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
            dir_prepro_orig_process, dir_prepro_orig_postprocessed,
            dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
            dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
            dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
            dir_prepro_template_process, dir_prepro_template_postprocessed, dir_prepro_acpc_matrices, dir_prepro_orig_matrices) \
            = getpath.func(data_path, reference, resting_or_task)

        for path in [dir_prepro_raw, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_raw_matrices, path_func,
            dir_prepro_orig, dir_prepro_orig_labels, dir_prepro_orig_masks,
            dir_prepro_orig_process, dir_prepro_orig_postprocessed,
            dir_prepro_acpc, dir_prepro_acpc_labels, dir_prepro_acpc_masks,
            dir_prepro_acpc_process, dir_prepro_acpc_postprocessed,
            dir_prepro_template, dir_prepro_template_labels, dir_prepro_template_masks,
            dir_prepro_template_process, dir_prepro_template_postprocessed, dir_prepro_acpc_matrices, dir_prepro_orig_matrices]:

            if not os.path.exists(path):
                os.makedirs(path)

        diary_file = diaryfile.create(opj(path_func, str(ID) + ' session ' + str(Session)), nl)
        launcher_parameters = diaryfile.create(opj(path_func, str(ID) + 'launcher_parameters session ' + str(Session)), nl)

        frame = inspect.currentframe()
        args_dict = inspect.getargvalues(frame).locals

        with open(launcher_parameters, "w") as f:
            f.write("Function: preprocess_data\n")
            f.write("Effective parameters at runtime:\n\n")
            for k, v in args_dict.items():
                if k == "kwargs":
                    for kk, vv in v.items():
                        f.write(f"- {kk} = {vv}\n")
                else:
                    f.write(f"- {k} = {v}\n")
        # ----------------------------------------

        diary_WARNING = opj(path_func, 'MAJOR_WARNING.txt')
        dir_RS_ICA_native = opj(dir_prepro_orig_process, 'Melodic')
        # link with the individual anat template ################################################################

        anat_session_used = Session
        anat_data_path_used = data_path

        if IhaveanANAT == False:
            if folderforTemplate_Anat == '':
                nl = 'ERROR: You need to provide either an anat image or a folder with template images to perform the registration'
                raise ValueError(run_cmd.error(nl, diary_file))
            else:
                (anat_subject, brainmask, G_mask, WBG_mask, V_mask, W_mask, dir_transfo, FS_dir,
                 dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.create(folderforTemplate_Anat,
                                                                                     use_erode_V_func_masks,
                                                                                     use_erode_WM_func_masks, species,
                                                                                     TfMRI, diary_file)
        else:
            (anat_subject, brainmask, G_mask, WBG_mask, V_mask, W_mask, dir_transfo, FS_dir,
             dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.retrieve(ID, data_path,
                                                                                   Session, anat_func_same_space,
                                                                                   use_erode_V_func_masks,
                                                                                   use_erode_WM_func_masks,
                                                                                   TfMRI, diary_file)

            if anat_subject is None or not ope(anat_subject):
                nl = (f'INFO: No anat found for {ID} session {Session}, '
                      f'searching across all sessions for this subject...')
                run_cmd.msg(nl, diary_file, 'WARNING')

                anat_found = False
                for ses in sorted(set(all_Session)):
                    candidate_path = opj(bids_dir, f'sub-{ID}', f'ses-{ses}')
                    if not ope(candidate_path):
                        continue
                    try:
                        (anat_subject, brainmask, G_mask, WBG_mask, V_mask, W_mask, dir_transfo, FS_dir,
                         dir_prepro, volumes_dir, labels_dir, masks_dir) = chooseanat.retrieve(ID, candidate_path,
                                                                                               ses,
                                                                                               anat_func_same_space,
                                                                                               use_erode_V_func_masks,
                                                                                               use_erode_WM_func_masks,
                                                                                               TfMRI, diary_file)
                        if anat_subject is not None and ope(anat_subject):
                            nl = f'INFO: Found anat for {ID} in session {ses}, using it as fallback.'
                            run_cmd.msg(nl, diary_file, 'OKGREEN')
                            anat_session_used = ses  # ← on mémorise la session anat réelle
                            anat_data_path_used = candidate_path
                            anat_found = True
                            break
                    except Exception:
                        continue

                if not anat_found:
                    nl = f'ERROR: No anat found for {ID} in any session. Cannot proceed.'
                    raise FileNotFoundError(run_cmd.error(nl, diary_file))

        study_template_atlas_folder = opj(bids_dir, 'sty_template')
        targetsuffix = 'space-acpc_desc-SS'
        info = backtonative.get(ID, anat_data_path_used, bids_dir, anat_session_used, max_ses,
                                # ← utilise le path/session anat réels
                                targetsuffix, type_norm, BASE_SS, BASE_mask, BASE_atlas_folder,
                                study_template_atlas_folder, creat_study_template,
                                coregistration_longitudinal, reference, 'Final', species, diary_file)
        # Check the func runs
        list_RS = sorted(glob.glob(opj(path_func, endfmri)))
        RS = [opb(i) for i in list_RS]

        if len(list_RS) == 0:
            nl = ('WARNING : No func image found, we are look for an image define such as opj(dir_fMRI_Refth_RS, endfmri) and here it is ' +
                        str(opj(path_func, endfmri)) + ' I would check how you define "endfmri", let see if we can found this pattern in other subject to continu')
            run_cmd.msg(nl, diary_file, 'WARNING')

        else:
            nl = ("INFO: now let's check that this is a real a 4D fMRI image with enough time point as in define with the variable ntimepoint_treshold=" +
                        str(ntimepoint_treshold))
            run_cmd.msg(nl, diary_file, 'OKGREEN')

            list_RS_list = list_RS.copy()
            list_pop_index = []

            for imageF in list_RS_list:
                # Load the fMRI NIfTI image
                fmri_image = nib.load(imageF)
                # Get the shape of the image (x, y, z, t)
                image_shape = fmri_image.shape
                # Check the number of time points (4th dimension)
                if len(image_shape) == 4:
                    ntimepoint = image_shape[3]  # The 4th dimension represents time
                    nl = "INFO: " + f"Number of time points: {ntimepoint}"
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    if int(ntimepoint) < ntimepoint_treshold:
                        index_of_imageF = list_RS.index(imageF)
                        nl = "INFO: We will not analyze " + str(imageF) + " because there are not enough time points"
                        run_cmd.msg(nl, diary_file, 'WARNING')

                        diary_WARNING_file = diaryfile.create(diary_WARNING, nl)

                        list_RS.pop(index_of_imageF)
                        list_pop_index.append(index_of_imageF)
                else:
                    nl = "INFO: " + str(imageF) + " is not a 4D fMRI image"
                    run_cmd.msg(nl, diary_file, 'WARNING')
                    index_of_imageF = list_RS.index(imageF)
                    list_RS.pop(index_of_imageF)
                    list_pop_index.append(index_of_imageF)
                    diary_WARNING_file = diaryfile.create(diary_WARNING, nl)

            nb_run = len(list_RS)
            RS = [opb(i) for i in list_RS]

            # Diagnostic print to verify run mapping
            run_cmd.msg('--- Run mapping after filtering ---', diary_file, 'OKGREEN')
            for idx, (in_file, out_base) in enumerate(zip(list_RS, RS)):
                run_cmd.msg(f'  Index {idx:2d} : Input = {opb(in_file):40s}  →  Output base = {out_base}',
                            diary_file, 'OKGREEN')
            run_cmd.msg('--- End of mapping ---', diary_file, 'OKGREEN')

            # Setup for distortion correction using Fieldmaps
            if nb_run>0:
                #### find and correct the confound files
                for imageF in list_RS:
                    # Load the fMRI NIfTI image
                    fmri_image = nib.load(imageF)
                    # Get the shape of the image (x, y, z, t)
                    image_shape = fmri_image.shape
                    # Check the number of time points (4th dimension)
                    ntimepoint = image_shape[3]  # The 4th dimension represents time

                    if ope(opj(path_func, extract_filename(imageF) + '_confounds.tsv')):
                        confounds_df = pd.read_csv(opj(path_func, extract_filename(imageF) + '_confounds.tsv'), sep='\t')

                        # Check if the number of rows in the confounds file matches the number of time points
                        if len(confounds_df) != ntimepoint:
                            nl = ('Mismatch in the number of time points: ' + str(len(confounds_df)) +
                                  ' in confounds file vs ' + str(ntimepoint) + ' in fMRI image.')

                            diary_WARNING_file = diaryfile.create(diary_WARNING, nl)
                            run_cmd.msg(nl, diary_WARNING_file, 'WARNING')

                        # Remove the first T1_eq rows from the confounds DataFrame
                        if T1_eq > 0:
                            # Preserve the first column (assuming it's a label or name)
                            first_column = confounds_df.iloc[:, 0]
                            remaining_columns = confounds_df.iloc[:, 1:]

                            # Remove the first T1_eq rows from the remaining columns
                            remaining_columns = remaining_columns.iloc[T1_eq:]

                            # Concatenate the first column with the remaining columns
                            confounds_df = pd.concat([first_column, remaining_columns], axis=1)

                            # Save the modified confounds DataFrame back to the TSV file
                            confounds_df.to_csv(opj(dir_prepro_orig, extract_filename(imageF) + '_confounds_correct.tsv'), sep='\t', index=False)
                            print(f"Removed the first {T1_eq} TRs from the confounds file.")
                        else:
                            print("T1_eq is not positive; no rows removed.")
                    else:
                        print(f"Confounds file not found: {opj(opd(imageF), extract_filename(imageF) + '_confounds.tsv')}")

                # Find the fmap images
                list_map = sorted(glob.glob(opj(dir_fmap, endmap)))
                nl = "looking for fmap image with the command glob.glob(" + str(opj(dir_fmap, endmap))
                run_cmd.msg(nl, diary_file, 'OKGREEN')
                if len(list_map) == 0:
                    list_map = sorted(glob.glob(opj(path_func, endmap)))
                    nl = "No fmap found in fmap folder, looking now in func folder with the command glob.glob(" + str(
                        opj(path_func, endmap))
                    run_cmd.msg(nl, diary_file, 'OKGREEN')
                    if len(list_map) > 0:
                        nl = "WARNING: We found some fieldmap images in the func folder, be sure that this is what you want"
                        run_cmd.msg(nl, diary_file, 'WARNING')
                    else:
                        nl = "No fmap found in func folder either"
                        run_cmd.msg(nl, diary_file, 'WARNING')
                else:
                    nl = "We found " + str(list_map)
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    for path in [dir_fmap, dir_prepro_fmap]:
                        if not os.path.exists(path):
                            os.makedirs(path)

                if len(list_RS) > 0:
                    RS_map = [opb(i) for i in list_map]

                    ######### choose TOPUP strategy ##########
                    if len(list_map) == 0:
                        recordings = 'very_old'
                        nl = 'WARNING: There is no image available for building a fieldmaps'
                        run_cmd.msg(nl, diary_file, 'WARNING')

                    elif len(list_map) == 1:
                        comp = check_nii.comphd(list_map[0], opj(path_func, RS[int(REF_int) - 1]))  # fieldmap may live in fmap/, not func/

                        if comp == True:
                            recordings = 'old'  # there is only one AP recording to correct for field distortion
                        else:
                            recordings = 'very_old'  # the only one AP recording to correct for field distortion is useless !!!
                            nl = 'WARNING : Before moving on, check the quality of the AP image, you may decide to NOT use it for correction'
                            run_cmd.msg(nl, diary_file, 'WARNING')

                    elif len(list_map) == 2:
                        recordings = '2_mapdir'  # there one AP per PA recordings in total
                    elif len(list_map) > 2:
                        recordings = 'new'  # there is one AP per PA recordings to correct for field distortion
                        if len(list_pop_index) > 0:
                            list_pop_index.sort(reverse=True)
                            for index in list_pop_index:
                                list_map.pop(index)
                        if not len(list_map) == len(list_RS):
                            nl = 'ERROR: Check the runs. There is probably one or two broken files that has been repeated and that you should remove !'
                            raise NameError(run_cmd.error(nl, diary_file))

                    nl = 'INFO: recordings type detected: ' + str(recordings)
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    list_json = sorted(glob.glob(opj(path_func, endjson)))

                    # get useful information about the func images #############################################################
                    ### check if we found some .json file
                    if not list_json:
                        nl = 'WARNING: no .json found!!, you will need to at least provide the TR. Beware that no TOPUP correction can be applied if you do not provide the DwellT as well.'
                        run_cmd.msg(nl, diary_file, 'WARNING')

                    else:
                        f = open(list_json[0])
                        info_RS = json.load(f)

                    ## find metrics in header (for fun)
                    nl = 'Get infos about the func image'
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    ## TR (always a float from here on: it was left as a string when given in the launcher)
                    try:
                        TR_val, nl, level = resolve_TR(TR, info_RS if list_json else None)
                    except ValueError as exc:
                        raise Exception(run_cmd.error('ERROR: ' + str(exc), diary_file))
                    run_cmd.msg(nl, diary_file, level)

                    if Slice_timing_info == 'Auto':
                        try:
                            slice_timing = info_RS["SliceTiming"]
                            nl = "INFO: SliceTiming = " + str(slice_timing)
                            run_cmd.msg(nl, diary_file, 'OKGREEN')
                            kind, nl = check_slice_timing(slice_timing, TR_val)
                            if kind == 'error':
                                raise SliceTimingError(run_cmd.error('ERROR: ' + nl, diary_file))
                            if kind == 'warning':
                                run_cmd.msg(nl, diary_file, 'WARNING')

                            STC = map(str, slice_timing)
                            stc = ' '.join(STC)
                            with open(opj(dir_prepro_raw_process, 'stc.txt'), 'w') as f:
                                f.write(stc)
                            f.close()
                        except SliceTimingError:
                            raise
                        except:
                            cmd = 'export SINGULARITYENV_AFNI_NIFTI_TYPE_WARN="NO";' + sing_afni + '3dinfo -slice_timing ' + \
                                  list_RS[0]
                            nl = spgo(cmd).split('\n')
                            STC = nl[-1].split('|')
                            STC = list(map(float, STC))
                            if np.sum(STC) > 0:
                                nl = "INFO: SliceTiming = " + str(STC)
                                run_cmd.msg(nl, diary_file, 'OKGREEN')

                            else:
                                nl = "WARNING: Slice Timing not found, this will be particularly DANGEROUS, you SHOULD PROVIDE MANUALLY ONE!"
                                run_cmd.msg(nl, diary_file, 'WARNING')

                                diary_WARNING_file = diaryfile.create(diary_WARNING, nl)

                    elif isinstance(Slice_timing_info, list) == True:
                        nl = "INFO: SliceTiming = " + str(Slice_timing_info)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                        kind, nl = check_slice_timing(Slice_timing_info, TR_val)
                        if kind == 'error':
                            raise Exception(run_cmd.error('ERROR: ' + nl, diary_file))
                        if kind == 'warning':
                            run_cmd.msg(nl, diary_file, 'WARNING')

                        STC = map(str, Slice_timing_info)
                        stc = ' '.join(STC)
                        with open(opj(dir_prepro_orig, 'stc.txt'), 'w') as f:
                            f.write(stc)
                        f.close()

                    try:
                        TE = info_RS["EchoTime"]
                        nl = "INFO: EchoTime = " + str(TE)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                    except:
                        nl = "INFO: EchoTime not found in header"
                        run_cmd.msg(nl, diary_file, 'WARNING')

                    try:
                        EES = info_RS["EffectiveEchoSpacing"]
                        nl = "INFO: Effective Echo Spacing = " + str(EES)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                    except:
                        nl = "INFO: Effective Echo Spacing not found in header"
                        run_cmd.msg(nl, diary_file, 'WARNING')

                    TRT_val, nl, level = resolve_TRT(TRT, info_RS if list_json else None)
                    run_cmd.msg(nl, diary_file, level)

                    # Find correction_direction
                    if correction_direction == 'Auto':
                        nl = 'INFO: input correction_direction was empty, let s try to find what is with the header'
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                        try:
                            PE_d2 = info_RS['PhaseEncodingDirection']
                        except:
                            nl = "INFO: Phase Encoding Direction not found in header"
                            run_cmd.msg(nl, diary_file, 'WARNING')

                            recordings = 'very_old'
                            nl = 'WARNING : No distortion correction will be applied with fugue'
                            run_cmd.msg(nl, diary_file, 'WARNING')

                            correction_direction_val = 'None'
                            dmap = ''
                            dbold = ''
                            PE_d2 = 'None'
                    else:
                        nl = 'INFO: input correction_direction is the launcher was determined as' + str(correction_direction)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                        PE_d2 = 'None'

                    if TRT_val != 'None':
                        if PE_d2 == 'j' or correction_direction == 'y-':
                            dmap = '0 1 0 ' + str(TRT_val)
                            dbold = '0 -1 0 ' + str(TRT_val)
                            correction_direction_val = 'y-'
                        elif PE_d2 == 'j-' or correction_direction == 'y':
                            dmap = '0 -1 0 ' + str(TRT_val)
                            dbold = '0 1 0 ' + str(TRT_val)
                            correction_direction_val = 'y'
                        elif PE_d2 == 'i' or correction_direction == 'x-':
                            dmap = '1 0 0 ' + str(TRT_val)
                            dbold = '-1 0 0 ' + str(TRT_val)
                            correction_direction_val = 'x-'
                        elif PE_d2 == 'i-' or correction_direction == 'x':
                            dmap = '-1 0 0 ' + str(TRT_val)
                            dbold = '1 0 0 ' + str(TRT_val)
                            correction_direction_val = 'x'
                        else:
                            recordings = 'very_old'
                            nl = 'WARNING : No distortion correction will be applied with fugue'
                            run_cmd.msg(nl, diary_file, 'WARNING')
                            dmap = ''
                            dbold = ''
                            correction_direction_val = 'None'

                    else:
                        nl = 'WARNING: TRT not found'
                        run_cmd.msg(nl, diary_file, 'WARNING')

                        nl = 'WARNING : No distortion correction will be applied with fugue'
                        run_cmd.msg(nl, diary_file, 'WARNING')

                        recordings = 'very_old'
                        dmap = ''
                        dbold = ''
                        correction_direction_val = ''

                    # Find SliceEncodingDirection
                    if SED == 'Auto':
                        try:
                            SED_val = info_RS["SliceEncodingDirection"]
                        except KeyError:
                            try:
                                if info_RS["ImageOrientationPatientDICOM"][0] == 1:
                                    SED_val = "i"
                                elif info_RS["ImageOrientationPatientDICOM"][0] == -1:
                                    SED_val = "i-"
                                elif info_RS["ImageOrientationPatientDICOM"][1] == 1:
                                    SED_val = "j"
                                elif info_RS["ImageOrientationPatientDICOM"][1] == -1:
                                    SED_val = "j-"
                                elif info_RS["ImageOrientationPatientDICOM"][2] == 1:
                                    SED_val = "k"
                                elif info_RS["ImageOrientationPatientDICOM"][2] == -1:
                                    SED_val = "k-"
                                else:
                                    SED_val = 'None'
                                    nl = 'WARNING !!!! Can not find the Slice Encoding Direction. No restriction of deformation can be applied (not a big deal)'
                                    run_cmd.msg(nl, diary_file, 'WARNING')

                                    nl = 'SED = ' + str(SED_val)
                                    run_cmd.msg(nl, diary_file, 'WARNING')

                            except KeyError:
                                SED_val = 'None'
                                nl = 'WARNING !!!! Can not find the Slice Encoding Direction. No restriction of deformation can be applied (not a big deal)'
                                run_cmd.msg(nl, diary_file, 'WARNING')

                                nl = 'SED = ' + str(SED_val)
                                run_cmd.msg(nl, diary_file, 'OKGREEN')
                    else:
                        SED_val = SED
                        nl = 'SED = ' + str(SED_val)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    ### Effective echo spacing for FUGUE (--dwell)
                    DwellT_val, nl, level = resolve_dwell(DwellT, info_RS if list_json else None)
                    run_cmd.msg(nl, diary_file, level)

                    DIR = os.getcwd()
                    nl = 'Working path : ' + DIR
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                    # Run the preprocessing ####################################################################################

                    if 1 in Skip_step:
                        nl = 'skip step ' + str(1)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        # brain mask for the motion estimation: 'previous' = the functional mask written by step 3
                        # of an earlier execution (it does not exist yet at the first one), or a path
                        motion_mask_path = None
                        if str(motion_mask) == 'previous':
                            _prev = opj(dir_prepro_raw_masks, ID + '_fMRI_mask.nii.gz')
                            if os.path.exists(_prev):
                                motion_mask_path = _prev
                            else:
                                run_cmd.msg('INFO: motion_mask=previous but ' + _prev + ' does not exist yet (first execution): '
                                            'the motion is estimated without a mask. Run the pipeline again from step 1 once '
                                            'the functional mask of step 3 exists.', diary_file, 'WARNING')
                        elif str(motion_mask) not in ('None', ''):
                            if not os.path.exists(str(motion_mask)):
                                raise Exception(run_cmd.error('ERROR: motion_mask ' + str(motion_mask) + ' does not exist.', diary_file))
                            motion_mask_path = str(motion_mask)
                        _1_fMRI_preTTT_in_fMRIspace.preprocess_data(dir_prepro_raw_process, RS, list_RS, nb_run, T1_eq, TR_val, Slice_timing_info, dir_prepro_raw_matrices, n_for_ANTS,
                                                                    overwrite, sing_afni, diary_file, animalPosition, humanPosition, orientation, doWARPonfunc, diary_WARNING,
                                                                    motion_mask_path=motion_mask_path, motion_mask_erode=motion_mask_erode,
                                                                    motion_method=motion_method, sing_fsl=sing_fsl)

                    if 2 in Skip_step:
                        nl = 'skip step ' + str(2)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _2_coregistration_to_norm.coregist_to_norm(correction_direction, list_RS, dir_prepro_fmap, dir_prepro_raw_process, dir_prepro_orig_process, RS, RS_map, nb_run, recordings, bids_dir,
                                                                   REF_int, list_map, animalPosition, humanPosition, doWARPonfunc, dir_prepro_raw_matrices, orientation, DwellT, n_for_ANTS,
                                                                   overwrite, sing_afni, sing_fsl, dmap, dbold, config_f, diary_file)

                    if 3 in Skip_step:
                        nl = 'skip step ' + str(3)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _3_mask_fMRI.Refimg_to_meanfMRI(MAIN_PATH, anat_func_same_space, BASE_SS, TfMRI, dir_prepro_raw_process, dir_prepro_raw_masks, dir_prepro_acpc_masks, dir_prepro_acpc_process,
                                                        dir_prepro_template_process, RS, nb_run, REF_int, ID, dir_transfo, brainmask, V_mask, W_mask, G_mask, WBG_mask, dilate_mask, n_for_ANTS, bids_dir,
                                                        costAllin, anat_subject, Method_mask_func, overwrite, type_of_transform, aff_metric_ants, extra_erode,
                                                        sing_afni, sing_fs, sing_fsl, sing_itk, diary_file)

                    if 'itk_1' in Skip_step:
                        nl = 'skip step ' + str('itk_1')
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _4_check_mask._itk_check_masks(dir_prepro_raw_masks, dir_prepro_raw_process, ID, sing_itk, diary_file, sing_afni, overwrite)

                    if 5 in Skip_step:
                        nl = 'skip step ' + str(5)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        # phase-encoding direction (JSON, or correction_direction given in the launcher):
                        # the non-linear func->anat registration is restricted to this axis
                        pe_for_restriction = PE_d2 if str(PE_d2) not in ('None', '') else correction_direction
                        _5_anat_to_fMRI.Refimg_to_meanfMRI(pe_for_restriction, anat_func_same_space, TfMRI, dir_prepro_raw_process, RS, nb_run, ID, bids_dir, dir_prepro_raw_masks, REF_int, dir_prepro_raw_matrices, recordings,
                                                           n_for_ANTS, aff_metric_ants, aff_metric_ants_Transl, list_atlas, labels_dir, anat_subject, dir_transfo, IhaveanANAT, do_anat_to_func, TR_val,
                                                           type_of_transform, registration_type, dir_prepro_acpc_masks, dir_prepro_acpc_process, dir_prepro_orig_masks, dir_prepro_acpc_labels,
                                                           dir_prepro_orig_labels, BASE_atlas_folder, dir_prepro_orig_process, doWARPonfunc, BASE_atlas_folder, species, opd(BASE_mask),
                                                           overwrite, sing_afni, diary_file)

                    if 'itk_2' in Skip_step:
                        nl = 'skip step ' + str('itk_2')
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                    else:
                        itk_2._itk_check_func_in_template(TfMRI, dir_prepro_raw_process, dir_prepro_raw_masks, sing_itk, diary_file)

                    if 6 in Skip_step or ICA_cleaning == 'Skip':
                        nl = 'skip step ' + str(6)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')


                    else:
                        _6_Melodic.Melodic_correct(dir_RS_ICA_native, dir_RS_ICA_native, dir_prepro_orig_process, dir_prepro_orig_masks, ID,
                                                   nb_ICA_run, nb_run, RS, ICA_cleaning, sing_fsl, sing_itk, TR_val, TfMRI, diary_file)


                    if 7 in Skip_step:
                        nl = 'skip step ' + str(7)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _7_post_TTT.signal_regression(dir_prepro_orig_process, dir_RS_ICA_native, dir_prepro_orig_masks, dir_prepro_raw_process,
                                                      nb_run, RS, blur, TR_val, ICA_cleaning, extract_WM, normalize, ID, post_treatment_method,
                                                      do_not_correct_signal, band, extract_Vc, extract_GS, dir_prepro_orig_postprocessed, dir_prepro_raw_matrices, overwrite, sing_afni, sing_fsl, diary_file,
                                                      acompcor_mode=acompcor_mode, acompcor_n_fixed=acompcor_n_fixed,
                                                      acompcor_variance=acompcor_variance, acompcor_max=acompcor_max)

                    if 8 in Skip_step:
                        nl = 'skip step ' + str(8)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _8_fMRI_to_anat.to_anat_space(dir_prepro_acpc_process, dir_prepro_orig_process, bids_dir, ID, TfMRI,
                                                      nb_run, RS, n_for_ANTS, do_anat_to_func, anat_func_same_space, dir_prepro_acpc_postprocessed, dir_prepro_orig_postprocessed, diary_file)

                    if 9 in Skip_step:
                        nl = 'skip step ' + str(9)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _9_coregistration_to_template_space.to_common_template_space(dir_prepro_template_process, bids_dir, ID, dir_prepro_template_labels, n_for_ANTS,
                                                                                     dir_prepro_orig_postprocessed, dir_prepro_acpc_postprocessed, dir_prepro_template_postprocessed,
                                                                                     nb_run, RS, do_anat_to_func, list_atlas, info, dir_prepro_orig_process, species,
                                                                                     BASE_atlas_folder, opd(BASE_mask), anat_func_same_space, dir_prepro_acpc_process,
                                                                                     dir_prepro_template_masks, IhaveanANAT, use_erode_WM_func_masks, overwrite, sing_afni, diary_file)

                    if 10 in Skip_step:
                        nl = 'skip step ' + str(10)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _10_Correl_matrix.correl_matrix(dir_prepro_orig_postprocessed, RS, nb_run, selected_atlases_matrix, segmentation_name_list,
                                                        ID, Session, TR_val, dir_prepro_orig_labels, dir_prepro_orig, matrix_fit,
                                                        sing_afni, diary_file)

                    if 11 in Skip_step:
                        nl = 'skip step ' + str(11)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _11_Seed_base_many_regionsatlas.SBA(SBAspace, BASE_SS, erod_seed, dir_prepro_orig_labels, dir_prepro_orig, dir_prepro_orig_process,
                                                            dir_prepro_acpc_labels, dir_prepro_acpc, dir_prepro_acpc_postprocessed, anat_subject, dir_prepro_acpc_process,
                                                            RS, nb_run, selected_atlases, panda_files, oversample_map, use_cortical_mask_func, dir_prepro_acpc_masks, TfMRI, ID,
                                                            dir_prepro_template_postprocessed, dir_prepro_template_labels, dir_prepro_template_masks, dir_prepro_orig_postprocessed,
                                                            n_cut, threshold_val, sing_afni, diary_file, smoothSBA, TR_val, dir_prepro_template, dir_prepro_template_process)

                    if 12 in Skip_step:
                        nl = 'skip step ' + str(12)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _12_fMRI_QC.fMRI_QC(correction_direction, path_func, ID, dir_prepro_template_process, dir_prepro_template_labels, dir_prepro_orig_masks, dir_prepro_orig_process, dir_prepro_orig_postprocessed, dir_prepro_raw_matrices,
                                            dir_prepro_template_postprocessed, dir_prepro_raw_process, dir_prepro_orig_labels, RS, nb_run, sing_afni, diary_file)

                    if 13 in Skip_step:
                        nl = 'skip step ' + str(13)
                        run_cmd.msg(nl, diary_file, 'OKGREEN')

                    else:
                        _14_fMRI_QC_matrix.fMRI_QC_matrix(path_func, dir_prepro_orig, intra_thresh, delta_thresh, RS, nb_run, diary_file)

                    if 'itk_3' in Skip_step:
                        nl = 'skip step ' + str('itk_3')
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                    else:
                        _200_Data_QC._itk_check_func_in_template(dir_prepro_template_postprocessed, dir_prepro_template_masks,
                                                                 dir_prepro_template_process, sing_itk, diary_file, sing_afni, overwrite)

                    if 'Clean' in Skip_step:
                        nl = 'skip step ' + str("Clean")
                        run_cmd.msg(nl, diary_file, 'OKGREEN')
                    else:
                        _100_Data_Clean.clean(dir_prepro_raw_process, dir_prepro_fmap, dir_prepro_acpc_process, dir_prepro_orig_process,
                                              dir_prepro_template_process, diary_file)
            else:
                nl = 'WARNING: No fMRI images were above the number of time point (' + str(ntimepoint) + '), we skip the subject'
                run_cmd.msg(nl, diary_file, 'WARNING')