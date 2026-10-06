import os
import json
import traceback
import numpy as np
import ants
import pandas as pd
import nibabel as nib
from scipy import signal as scipy_signal
from sklearn.decomposition import PCA

opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

from Tools import run_cmd
from modalities.fMRI.extract_filename import extract_filename


# =============================================================================
# Module-level helpers
# =============================================================================

def _write_json(nii_path, sources, description, command=''):
    """Write a minimal BIDS-style JSON sidecar next to a NIfTI file."""
    meta = {"Sources": sources, "Description": description}
    if command:
        meta["Command"] = command
    with open(nii_path.replace('.nii.gz', '.json'), 'w') as f:
        json.dump(meta, f, indent=3)


def _write_1d(lines_bytes, out_path):
    """
    Decode AFNI stdout bytes and write as a clean .1D file.
    Handles the repeated pattern of decoding + whitespace normalisation.
    """
    lines = lines_bytes.decode('utf-8').split('\n')
    with open(out_path, 'w') as f:
        for line in lines:
            f.write(' '.join(line.split()) + '\n')


def _load_4d(func):
    """4D data as float32 (a path or an array); float64 doubled the memory for nothing."""
    if isinstance(func, np.ndarray):
        return func
    return np.asarray(nib.load(func).dataobj, dtype=np.float32)


def _dct_basis(n_trs, tr_sec, hp_cutoff):
    """
    Cosine basis of the frequencies below hp_cutoff (SPM / fMRIPrep convention):
    floor(2 * n * TR * hp) regressors, without the constant. Empty if hp_cutoff <= 0.
    """
    if hp_cutoff <= 0:
        return np.zeros((n_trs, 0))
    n_cos = int(np.floor(2 * n_trs * tr_sec * hp_cutoff))
    t = np.arange(n_trs)
    return np.column_stack([np.cos(np.pi * k * (2 * t + 1) / (2 * n_trs))
                            for k in range(1, n_cos + 1)]) if n_cos else np.zeros((n_trs, 0))


def _tissue_timeseries(data, mask_path, tr_sec, hp_cutoff):
    """
    Voxel time series of a tissue mask, (T, n_voxels), prepared for CompCor as in
    fMRIPrep (pre_filter='cosine'): high-pass with the cosine basis of the analysis
    band (linear detrend if there is no high-pass), then standardised per voxel.
    """
    mask = np.asarray(nib.load(mask_path).dataobj) > 0
    if mask.shape != data.shape[:3]:
        raise ValueError(f'aCompCor: mask {mask_path} {mask.shape} does not match the data {data.shape[:3]}')
    ts = data[mask].T.astype(np.float64)
    if ts.shape[1] == 0:
        raise ValueError(f'aCompCor: mask {mask_path} contains no voxels.')
    n_trs = ts.shape[0]
    basis = np.column_stack([np.ones(n_trs), _dct_basis(n_trs, tr_sec, hp_cutoff)])
    if basis.shape[1] > 1:
        ts = ts - basis @ np.linalg.lstsq(basis, ts, rcond=None)[0]
    else:
        ts = scipy_signal.detrend(ts, axis=0, type='linear')
    std = ts.std(axis=0)
    std[std == 0] = 1.0
    return ts / std


def _acompcor(func, mask_path, n_components=3, tr_sec=None, hp_cutoff=0.0):
    """
    Anatomical CompCor with a fixed number of components (Behzadi et al. 2007;
    fMRIPrep: Esteban et al. 2019). `func` is a path or a (X, Y, Z, T) array.
    Returns (components (T, n), fraction of variance of each component).
    """
    ts = _tissue_timeseries(_load_4d(func), mask_path, tr_sec or 1.0, hp_cutoff if tr_sec else 0.0)
    n_comp = min(n_components, ts.shape[1], ts.shape[0] - 1)
    U, S, _ = np.linalg.svd(ts, full_matrices=False)
    return U[:, :n_comp], (S[:n_comp] ** 2) / np.sum(S ** 2)


def _acompcor_adaptive(func, mask_path, variance_threshold=0.50, max_components=5,
                       min_components=1, tr_sec=None, hp_cutoff=0.0):
    """
    Adaptive aCompCor (aCompCor50-like): as many components as needed to explain
    `variance_threshold` of the variance in the mask, capped at `max_components`.
    One SVD (the previous version fitted the PCA twice).
    Returns (components (T, n), fraction of variance of each, n).
    """
    ts = _tissue_timeseries(_load_4d(func), mask_path, tr_sec or 1.0, hp_cutoff if tr_sec else 0.0)
    U, S, _ = np.linalg.svd(ts, full_matrices=False)
    var = (S ** 2) / np.sum(S ** 2)
    n_max = min(ts.shape[0] - 1, ts.shape[1])
    n_needed = int(np.searchsorted(np.cumsum(var), variance_threshold) + 1)
    n_comp = max(min_components, min(n_needed, max_components, n_max))
    return U[:, :n_comp], var[:n_comp], n_comp


def _phase_randomized(ts, rng):
    """
    Surrogate of a (T, V) matrix with the same power spectrum in every column (so the same
    variance and temporal autocorrelation) and independent random phases across columns:
    only the structure shared between voxels is destroyed.
    """
    n = ts.shape[0]
    spec = np.fft.rfft(ts, axis=0)
    phase = rng.uniform(0.0, 2.0 * np.pi, spec.shape)
    phase[0] = 0.0                                  # constant term
    if n % 2 == 0:
        phase[-1] = 0.0                             # Nyquist stays real
    return np.fft.irfft(spec * np.exp(1j * phase), n=n, axis=0)


def _n_components_parallel(ts, n_surrogates=20, quantile=0.99, max_voxels=5000, seed=0):
    """
    Number of principal components of a (T, V) matrix of tissue time series that stand above
    chance (parallel analysis, Horn 1965): component k is kept, together with all the previous
    ones, if its eigenvalue exceeds the `quantile` of the k-th eigenvalue of phase-randomised
    surrogates. Unlike a variance threshold, it does not depend on how the noise variance is
    spread over components: with spatially independent noise a 50 % variance criterion needs
    many noise components, while this one returns the number of shared signals.
    Large masks are subsampled to `max_voxels` columns (the leading components are unchanged).
    Returns (n, observed eigenvalues, null threshold).
    """
    rng = np.random.default_rng(seed)
    if ts.shape[1] > max_voxels:
        ts = ts[:, rng.choice(ts.shape[1], max_voxels, replace=False)]

    def eigenvalues(x):
        x = x - x.mean(axis=0)
        return np.sort(np.linalg.eigvalsh(x @ x.T))[::-1]          # squared singular values, via the T x T matrix

    observed = eigenvalues(ts)
    null = np.array([eigenvalues(_phase_randomized(ts, rng)) for _ in range(n_surrogates)])
    threshold = np.quantile(null, quantile, axis=0)
    n_max = min(ts.shape[0] - 1, ts.shape[1])
    n = 0
    while n < n_max and observed[n] > threshold[n]:
        n += 1
    return n, observed, threshold


def _acompcor_parallel(func, mask_path, max_components=15, tr_sec=None, hp_cutoff=0.0,
                       n_surrogates=20, quantile=0.99, seed=0):
    """
    aCompCor with the number of components chosen against a null (see _n_components_parallel).
    max_components is a safety cap only; `capped` tells whether it was reached.
    Returns (components (T, n), fraction of variance of each, n, capped). n may be 0.
    """
    ts = _tissue_timeseries(_load_4d(func), mask_path, tr_sec or 1.0, hp_cutoff if tr_sec else 0.0)
    n, _, _ = _n_components_parallel(ts, n_surrogates, quantile, seed=seed)
    capped = n > max_components
    n = min(n, max_components)
    U, S, _ = np.linalg.svd(ts, full_matrices=False)
    var = (S ** 2) / np.sum(S ** 2)
    return U[:, :n], var[:n], n, capped


def _parse_band(band, tr_sec):
    """
    band: 'hp' (high-pass) or 'hp lp' (band-pass), in Hz; hp may be 0 (low-pass only).
    Returns hp, lp (None for high-pass only), nyquist, 1dBport argument, 3dTproject argument.
    """
    try:
        vals = [float(v) for v in str(band).split()]
    except ValueError:
        raise ValueError(f'band={band!r}: expected one or two frequencies in Hz, e.g. "0.01 0.1".')
    if len(vals) not in (1, 2) or vals[0] < 0:
        raise ValueError(f'band={band!r}: expected one or two non-negative frequencies in Hz.')
    hp = vals[0]
    lp = vals[1] if len(vals) == 2 else None
    nyquist = 1.0 / (2.0 * tr_sec)
    if lp is not None and lp <= hp:
        raise ValueError(f'band={band!r}: the low-pass cutoff must be above the high-pass cutoff.')
    if lp is None and hp <= 0:
        raise ValueError(f'band={band!r}: a single value is a high-pass cutoff and must be > 0.')
    if lp is not None and lp >= nyquist:
        lp = None if hp > 0 else lp          # cutoff at or above Nyquist: high-pass only
    if lp is not None:
        return hp, lp, nyquist, f' -band {hp} {lp} -invert -nozero', f' -passband {hp} {lp}'
    return hp, None, nyquist, f' -band {hp} {nyquist} -invert -nozero', f' -stopband 0 {hp}'


def _n_filter_regressors(n_trs, tr_sec, hp, lp):
    """Number of sine / cosine regressors the temporal filter removes (1dBport, 3dTproject)."""
    f = np.arange(1, n_trs // 2 + 1) / (n_trs * tr_sec)
    removed = (f < hp) | ((f > lp) if lp is not None else np.zeros_like(f, dtype=bool))
    n = 2 * int(removed.sum())
    if n_trs % 2 == 0 and removed[-1]:
        n -= 1                                # the Nyquist frequency has no sine
    return n


MIN_TISSUE_VOXELS = 10      # a mean or a PCA over fewer voxels is not a tissue signal (small brains, eroded masks)


def _mask_voxels(path):
    """Number of voxels of a binary mask."""
    return int((np.asarray(nib.load(path).dataobj) > 0).sum())


def _tissue_too_small(label, mask_path, root_RS, diary_file, report):
    """
    True (with a warning) when a tissue mask has fewer than MIN_TISSUE_VOXELS voxels. The voxel
    count and the decision are recorded in `report`, written next to the residual, because the
    regressors then differ between subjects and species and must be reported.
    """
    n_vox = _mask_voxels(mask_path)
    too_small = n_vox < MIN_TISSUE_VOXELS
    report[label] = dict(mask=opb(mask_path), voxels=n_vox, used=not too_small)
    if too_small:
        run_cmd.msg(f'WARNING: {label} mask has {n_vox} voxel(s) (< {MIN_TISSUE_VOXELS}): the {label} regressor is NOT '
                    f'used for {root_RS}. Small brains and eroded masks can leave no tissue at the functional '
                    f'resolution: report which regressors each subject has.', diary_file, 'WARNING')
    else:
        run_cmd.msg(f'INFO: {label} mask: {n_vox} voxels', diary_file, 'OKCYAN')
    return too_small


def _header_tr(path):
    """TR (s) stored in the header of a 4D NIfTI, as AFNI reads it (unknown time unit = seconds)."""
    hdr = nib.load(path).header
    zooms = hdr.get_zooms()
    if len(zooms) < 4:
        raise ValueError(f'{opb(path)} is not a 4D image')
    unit = hdr.get_xyzt_units()[1]
    return float(zooms[3]) * {'sec': 1.0, 'msec': 1e-3, 'usec': 1e-6}.get(unit, 1.0)


def _check_header_tr(path, tr_sec, tol=1e-3):
    """
    3dTproject takes the TR from the header of its input for -passband / -stopband, whereas
    1dBport receives it on the command line. A header TR different from the TR of the analysis
    makes the filter remove the wrong frequencies, silently: refuse it.
    """
    tr_hdr = _header_tr(path)
    if abs(tr_hdr - tr_sec) > tol * max(tr_sec, 1e-9):
        raise ValueError(f'{opb(path)} has TR = {tr_hdr:g} s in its header but the TR of the analysis is {tr_sec:g} s: '
                         f'3dTproject reads the TR in the header and would filter the wrong frequencies. First decide which '
                         f'value is the TR of the acquisition protocol. If it is the TR of the analysis and this file was written '
                         f'by an earlier run of step 5 with another TR, run step 5 again (remove 5 from Skip_step), which writes '
                         f'the TR in the header; otherwise correct the TR given to the pipeline (TR=... in the launcher), '
                         f'or the header (3drefit -TR {tr_sec:g} {path}).')


def _check_dof(n_kept, n_regressors, n_filter, label, diary_file):
    """
    Degrees of freedom left after regression. Raises if none are left (the residual
    would be meaningless); warns below 20 % of the uncensored volumes.
    """
    dof = n_kept - n_regressors - n_filter
    msg = (f'{label}: {n_kept} uncensored volumes - {n_regressors} regressors - '
           f'{n_filter} filter regressors = {dof} degrees of freedom left')
    if dof <= 0:
        raise ValueError(msg + '. Reduce the regressors (fewer aCompCor components, no global signal, '
                               'high-pass only) or use longer runs.')
    run_cmd.msg(('WARNING: ' if dof < 0.2 * n_kept else 'INFO: ') + msg,
                diary_file, 'WARNING' if dof < 0.2 * n_kept else 'OKGREEN')
    return dof


def _n_kept(censor_1d, n_trs):
    if censor_1d and ope(censor_1d):
        return int(np.loadtxt(censor_1d).astype(bool).sum())
    return n_trs


def _read_confounds(tsv_path, n_trs, diary_file=None):
    """
    Custom confounds table (one column per regressor, one row per volume) as a numeric DataFrame.
    Checks the number of rows. fMRIPrep tables have NaN in the first row of the derivative columns,
    and 3dTproject / 3dDeconvolve refuse NaN: they are replaced by 0 (the value of a derivative
    before the first volume), with a warning naming the columns. Text columns are dropped.
    """
    df = pd.read_csv(tsv_path, sep='\t')
    if len(df) != n_trs:
        raise ValueError(f'{opb(tsv_path)} has {len(df)} rows, expected {n_trs} (one per volume)')
    num = df.apply(pd.to_numeric, errors='coerce')
    text_cols = [c for c in df.columns if num[c].isna().all() and not df[c].isna().all()]
    if text_cols and diary_file is not None:
        run_cmd.msg(f'WARNING: non-numeric columns of {opb(tsv_path)} ignored: {text_cols}', diary_file, 'WARNING')
    num = num.drop(columns=text_cols)
    nan_cols = [c for c in num.columns if num[c].isna().any()]
    if nan_cols:
        if diary_file is not None:
            run_cmd.msg(f'WARNING: NaN replaced by 0 in {opb(tsv_path)}: {nan_cols}', diary_file, 'WARNING')
        num = num.fillna(0.0)
    return num


def _save_regressors_tsv(regressors_dict, out_path):
    """
    Save a dictionary of {name: array} regressors as a TSV file.
    Arrays must all have the same length.
    """
    df = pd.DataFrame(regressors_dict)
    df.to_csv(out_path, sep='\t', index=False, float_format='%.8f')


# =============================================================================
# Main signal regression function
# =============================================================================

def signal_regression(dir_prepro_orig_process, dir_RS_ICA_native,
                      dir_prepro_orig_masks, dir_prepro_raw_process,
                      nb_run, RS, blur, TR, ICA_cleaning,
                      extract_WM, normalize, ID,
                      post_treatment_method, do_not_correct_signal, band,
                      extract_Vc, extract_GS,
                      dir_prepro_orig_postprocessed, dir_prepro_raw_matrices,
                      overwrite, sing_afni, sing_fsl, diary_file,
                      acompcor_mode='fixed',
                      acompcor_n_fixed=5,
                      acompcor_variance=0.50,
                      acompcor_max=15):

    run_cmd.msg(
        '##  Working on step 7 (function: _7_post_TTT).  ##',
        diary_file, 'HEADER')

    os.makedirs(dir_RS_ICA_native, exist_ok=True)
    os.makedirs(dir_prepro_orig_postprocessed, exist_ok=True)

    # =========================================================================
    # Parse band argument once — supports 1 value (highpass) or 2 (bandpass).
    # The mode is determined automatically here and reused by every method.
    #
    #   1) highpass only (1 value, e.g. band='0.01')
    #        1dBport  : -band hp nyquist -invert -nozero
    #        3dTproject: -stopband 0 hp   (AFNI has no -highpass)
    #
    #   2) bandpass (2 values, e.g. band='0.01 0.1')
    #        1dBport  : -band hp lp -invert -nozero
    #        3dTproject: -passband hp lp
    #
    # Note: AFNI's -passband requires exactly TWO values. Passing a single
    # value is interpreted as an illegal -stopband range.
    # =========================================================================
    tr_sec = float(TR)
    if acompcor_mode not in ('fixed', 'adaptive', 'parallel'):
        raise ValueError(f"acompcor_mode={acompcor_mode!r}: expected 'fixed', 'adaptive' or 'parallel'.")
    hp_cutoff, lp_cutoff, nyquist, bport_band_arg, tproject_filter_arg = _parse_band(band, tr_sec)
    if lp_cutoff is not None:
        run_cmd.msg(f'INFO: band-pass {hp_cutoff}-{lp_cutoff} Hz (Nyquist {nyquist:.4f} Hz).',
                    diary_file, 'OKGREEN')
    else:
        run_cmd.msg(f'INFO: high-pass {hp_cutoff} Hz (Nyquist {nyquist:.4f} Hz; a low-pass cutoff '
                    f'at or above Nyquist is ignored).', diary_file, 'OKGREEN')

    # Masks of the functional AC-PC space (written at step 5)
    brain_mask = opj(dir_prepro_orig_masks, 'mask_ref.nii.gz')
    if not ope(brain_mask):
        brain_mask = opj(dir_prepro_orig_masks, ID + '_space-acpc-func_desc-fMRI_mask_dilated.nii.gz')
    tissue_masks = {
        # exterior CSF needs its own mask; the ventricle mask is extract_Vc
        '_Wc':   (extract_WM,           opj(dir_prepro_orig_masks, 'Wmask.nii.gz')),
        '-GS':   (extract_GS,           brain_mask),
        '_Vc':   (extract_Vc,           opj(dir_prepro_orig_masks, 'Vmask.nii.gz')),
    }
    failed_runs = []

    # =========================================================================
    # LOOP 1 — QC statistics + nuisance signal extraction + bandpass regressors
    # =========================================================================

    for run_idx in range(int(nb_run)):

        root_RS = extract_filename(RS[run_idx])
        run_cmd.msg(f'INFO: QC + regressor extraction — run {run_idx + 1}/{nb_run}',
                    diary_file, 'OKGREEN')

        # Input selection
        if ICA_cleaning == 'Skip':
            func_input = opj(dir_prepro_orig_process,
                             root_RS + '_space-acpc-func_desc-fMRI_run_inRef_SS.nii.gz')
        else:
            func_input = opj(dir_RS_ICA_native,
                             root_RS + '_norm_final_clean.nii.gz')

        # ---- QC statistics --------------------------------------------------
        # cvarinv = fabs(mean)/stdev with detrend  (pseudo-tSNR)
        # cvar    = stdev/fabs(mean) with detrend  (CoV)
        # tsnr    = fabs(mean)/stdev NOT detrended (raw tSNR)
        # stdev   = standard deviation with detrend
        for option, suffix, description in [
            (' -cvarinv', '_tsnr1', 'fabs(mean)/stdev with detrend (3dTstat, AFNI)'),
            (' -cvar',    '_cvar',  'stdev/fabs(mean) with detrend (3dTstat, AFNI)'),
            (' -tsnr',    '_tsnr2', 'fabs(mean)/stdev NOT detrended — raw tSNR (3dTstat, AFNI)'),
            (' -stdev',   '_stdev', 'Standard deviation with detrend (3dTstat, AFNI)'),
        ]:
            out = opj(dir_prepro_orig_postprocessed,
                      root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.nii.gz')
            command = (sing_afni + '3dTstat' + overwrite + option +
                       ' -prefix ' + out + ' ' + func_input)
            run_cmd.run(command, diary_file)
            _write_json(out, func_input, description, command)

        # ---- Nuisance signal extraction via 3dmaskSVD -----------------------
        # First singular vector of the voxels of each mask (more robust than
        # the mean to partial volume). A requested regressor whose mask does not
        # exist is reported, not replaced by another mask.
        for suffix, (flag, mask_path) in tissue_masks.items():
            if not flag:
                continue
            out_1d = opj(dir_prepro_orig_postprocessed,
                         root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.1D')
            if ope(out_1d):
                os.remove(out_1d)            # never reuse a regressor of a previous execution
            if not ope(mask_path):
                run_cmd.msg(f'WARNING: {opb(mask_path)} not found: regressor {suffix} not extracted.',
                            diary_file, 'WARNING')
                continue
            if suffix in ('_Wc', '_Vc') and _tissue_too_small({'_Wc': 'WM', '_Vc': 'CSF'}[suffix], mask_path,
                                                              root_RS, diary_file, {}):
                continue
            command = (sing_afni + '3dmaskSVD' + overwrite +
                       ' -polort 2 -vnorm -mask ' + mask_path + ' ' + func_input)
            val, _ = run_cmd.get(command, diary_file)
            _write_1d(val, out_1d)

        # ---- Bandpass regressors (sinusoidal basis, for GLM approach) -------
        # Using regressors rather than direct filtering preserves DOF around
        # censored TRs (Power et al. 2014). This is the AFNI / fMRIPrep-aligned
        # approach; the Grandjean and fMRIPrep_style methods use 3dTproject
        # filter options directly.
        hd    = ants.image_header_info(func_input)
        n_trs = int(hd['dimensions'][3])
        command = (sing_afni + '1dBport' + overwrite +
                   ' -nodata ' + str(n_trs) +
                   ' ' + str(TR) +
                   bport_band_arg)
        bp_bytes, _ = run_cmd.get(command, diary_file)
        bp_out = opj(dir_prepro_orig_postprocessed,
                     root_RS + '_bandpass_rall.1D')
        _write_1d(bp_bytes, bp_out)

    # =========================================================================
    # LOOP 2 — Signal regression (or bypass if do_not_correct_signal)
    # =========================================================================

    if do_not_correct_signal:
        # Bypass: copy input to residual slot so downstream steps have a
        # consistent file to read regardless of correction setting.
        for run_idx in range(int(nb_run)):
            root_RS = extract_filename(RS[run_idx])
            func_input = (
                opj(dir_prepro_orig_process,
                    root_RS + '_space-acpc-func_desc-fMRI_run_inRef.nii.gz')
                if ICA_cleaning == 'Skip' else
                opj(dir_RS_ICA_native, root_RS + '_norm_final_clean.nii.gz'))
            residual = opj(dir_prepro_orig_postprocessed,
                           root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
            command = (sing_afni + '3dcalc' + overwrite +
                       ' -a ' + func_input +
                       ' -prefix ' + residual + ' -expr "a"')
            run_cmd.do(command, diary_file)
            _write_json(residual, func_input,
                        'Bypass copy — do_not_correct_signal=True.', command)
        return

    for run_idx in range(int(nb_run)):

        root_RS = extract_filename(RS[run_idx])
        run_cmd.msg(
            f'INFO: signal regression ({post_treatment_method}) '
            f'— run {run_idx + 1}/{nb_run}',
            diary_file, 'OKGREEN')

        if ICA_cleaning == 'Skip':
            func_input = opj(dir_prepro_orig_process,
                             root_RS + '_space-acpc-func_desc-fMRI_run_inRef.nii.gz')
        else:
            func_input = opj(dir_RS_ICA_native,
                             root_RS + '_norm_final_clean.nii.gz')

        # Paths shared across all methods
        Mean               = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_Mean.nii.gz')
        Sdev               = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_run_inRef_Sdev.nii.gz')
        maskDilat          = opj(dir_prepro_orig_masks, ID + '_space-acpc-func_desc-fMRI_mask_dilated.nii.gz')
        bandpass           = opj(dir_prepro_orig_postprocessed, root_RS + '_bandpass_rall.1D')
        censore1D          = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-censor.1D')
        demean             = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-demean.1D')
        deriv              = opj(dir_prepro_raw_process, root_RS + '_space-func_desc-deriv.1D')
        motion_correction  = opj(dir_prepro_raw_matrices, root_RS + '_space-func_desc-motion_correction.1D')
        residual           = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')
        xmat               = opj(dir_prepro_orig_postprocessed, root_RS + '_X.xmat.1D')
        xmat_nocensor      = opj(dir_prepro_orig_postprocessed, root_RS + '_X.nocensor.xmat.1D')

        # Preserve working directory — restored unconditionally in finally
        original_dir = os.getcwd()

        try:
            os.chdir(dir_prepro_orig_process)

            # -----------------------------------------------------------------
            # METHOD: AFNI
            # Strategy: estimate GLM design matrix with 3dDeconvolve -x1D_stop
            # (no actual fitting), then project out the full regressor matrix
            # with 3dTproject. Two-step approach is more numerically stable
            # and handles censored TRs correctly via -cenmode ZERO.
            # Matches the recommended AFNI rs-fMRI denoising pipeline.
            # -----------------------------------------------------------------
            if post_treatment_method in ('Grandjean', 'fMRIPrep_style'):
                _check_header_tr(func_input, tr_sec)

            if post_treatment_method == 'AFNI':

                # Build 3dDeconvolve command (design matrix estimation only)
                command = (
                    sing_afni + '3dDeconvolve'
                    ' -input '               + func_input +
                    ' -mask '                + maskDilat +
                    ' -ortvec '              + bandpass + ' bandpass_rall' +
                    ' -ortvec '              + demean   + ' mot_demean' +
                    ' -ortvec '              + deriv    + ' mot_deriv' +
                    ' -censor '              + censore1D +
                    ' -polort A -float'
                    ' -num_stimts 0'         + overwrite +
                    ' -fout -tout'
                    ' -x1D '                 + xmat +
                    ' -xjpeg '               + opj(dir_prepro_orig_postprocessed, root_RS + '_X.jpg') +
                    ' -x1D_uncensored '      + xmat_nocensor +
                    ' -fitts '               + opj(dir_prepro_orig_postprocessed, root_RS + '_Xfitts') +
                    ' -errts '               + opj(dir_prepro_orig_postprocessed, root_RS + '_errts') +
                    ' -x1D_stop'
                    ' -bucket '              + opj(dir_prepro_orig_postprocessed, root_RS + '_stats'))

                # Tissue regressors actually extracted in loop 1
                n_tissue = 0
                for suffix, (flag, _mask) in tissue_masks.items():
                    reg_path = opj(dir_prepro_orig_postprocessed,
                                   root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.1D')
                    if flag and ope(reg_path) and os.path.getsize(reg_path) > 0:
                        command += f' -ortvec {reg_path} residual_norm{suffix}'
                        n_tissue += 1

                # Add custom confounds TSV columns if present
                confounds_tsv = opj(dir_prepro_orig_postprocessed,
                                    root_RS + '_confounds_correct.tsv')
                if ope(confounds_tsv):
                    confounds_df = _read_confounds(confounds_tsv, int(ants.image_header_info(func_input)['dimensions'][3]),
                                                   diary_file)
                    run_cmd.msg(
                        f'INFO: loading confounds TSV — columns: '
                        f'{list(confounds_df.columns)}',
                        diary_file, 'OKGREEN')
                    for col in confounds_df.columns:
                        col_file = opj(dir_prepro_orig_postprocessed,
                                       root_RS + f'_{col}.1D')
                        confounds_df[[col]].to_csv(col_file, sep=' ',
                                                   index=False, header=False)
                        command += f' -ortvec {col_file} residual_norm_{col}'
                        n_tissue += 1

                # degrees of freedom: filter (1dBport) + 12 motion + tissue + polynomials
                n_trs = int(ants.image_header_info(func_input)['dimensions'][3])
                n_poly = 1 + int(1 + n_trs * tr_sec / 150.0)       # 3dDeconvolve -polort A
                _check_dof(_n_kept(censore1D, n_trs), 12 + n_tissue + n_poly,
                           _n_filter_regressors(n_trs, tr_sec, hp_cutoff, lp_cutoff),
                           f'{root_RS} (AFNI)', diary_file)

                # logged and checked like every other command (subprocess bypassed the log)
                for f_ in (xmat, xmat_nocensor):
                    if ope(f_):
                        os.remove(f_)
                run_cmd.run(command, diary_file)
                if not ope(xmat_nocensor):
                    raise RuntimeError(f'3dDeconvolve did not write the design matrix {xmat_nocensor}')

                # Project out design matrix with optional spatial smoothing
                command = (
                    sing_afni + '3dTproject -polort 0' + overwrite +
                    ' -input '   + func_input +
                    ' -censor '  + censore1D +
                    ' -cenmode ZERO'
                    ' -ort '     + xmat_nocensor +
                    ' -prefix '  + residual)
                if blur > 0:
                    command += ' -blur ' + str(blur)
                run_cmd.run(command, diary_file)
                _write_json(residual, func_input,
                            'Signal regression (AFNI: 3dDeconvolve + 3dTproject).',
                            command)

                # Normalisation
                _apply_normalisation(residual, Mean, Sdev, normalize,
                                     sing_afni, diary_file,
                                     func_input=func_input, censor_1d=censore1D)

            # -----------------------------------------------------------------
            # METHOD: Grandjean
            # Strategy: simultaneous bandpass + regressor removal in 3dTproject,
            # then add back the temporal mean to restore baseline intensity.
            # Matches the Grandjean lab rodent rs-fMRI protocol (single-step,
            # mean-restoring). Standard for animal neuroimaging pipelines.
            # -----------------------------------------------------------------
            elif post_treatment_method == 'Grandjean':

                regressors_file = opj(dir_prepro_orig_postprocessed,
                                      root_RS + '_regressors.1D')
                wm_signal  = opj(dir_prepro_orig_postprocessed, root_RS + '_wm.1D')
                csf_signal = opj(dir_prepro_orig_postprocessed, root_RS + '_csf.1D')
                wm_mask    = opj(dir_prepro_orig_masks, 'Wmask.nii.gz')
                csf_mask   = opj(dir_prepro_orig_masks, 'Vmask.nii.gz')

                # Extract available tissue signals with fslmeants (simple mean)
                tissue_signals = [motion_correction]  # always include motion
                tissue_report = {}
                for mask_path, sig_path, label in [
                    (wm_mask,  wm_signal,  'WM'),
                    (csf_mask, csf_signal, 'CSF'),
                ]:
                    if ope(mask_path):
                        if _tissue_too_small(label, mask_path, root_RS, diary_file, tissue_report):
                            continue
                        run_cmd.msg(f'INFO: extracting {label} signal',
                                    diary_file, 'OKCYAN')
                        cmd = (sing_fsl + 'fslmeants'
                               ' -i ' + func_input +
                               ' -o ' + sig_path +
                               ' -m ' + mask_path)
                        run_cmd.run(cmd, diary_file)
                        tissue_signals.append(sig_path)
                    else:
                        tissue_report[label] = dict(mask=opb(mask_path), voxels=0, used=False)
                        run_cmd.msg(f'WARNING: {label} mask not found — skipped.',
                                    diary_file, 'WARNING')
                with open(opj(dir_prepro_orig_postprocessed, root_RS + '_tissue_regressors.json'), 'w') as fh:
                    json.dump(dict(method='Grandjean', tissues=tissue_report), fh, indent=2)

                run_cmd.msg(
                    f'INFO: concatenating regressors: '
                    f'{[opb(p) for p in tissue_signals]}',
                    diary_file, 'OKCYAN')
                if ope(regressors_file):
                    os.remove(regressors_file)
                run_cmd.do(sing_afni + '1dcat ' + ' '.join(tissue_signals) + ' > ' + regressors_file,
                           diary_file)
                if not ope(regressors_file) or os.path.getsize(regressors_file) == 0:
                    raise RuntimeError(f'1dcat did not write {regressors_file}')
                n_trs = int(ants.image_header_info(func_input)['dimensions'][3])
                n_reg = np.atleast_2d(np.loadtxt(regressors_file, comments='#')).shape
                n_reg = n_reg[1] if n_reg[0] == n_trs else n_reg[0]
                _check_dof(_n_kept(censore1D, n_trs), 3 + n_reg,          # 3dTproject default polort 2
                           _n_filter_regressors(n_trs, tr_sec, hp_cutoff, lp_cutoff),
                           f'{root_RS} (Grandjean)', diary_file)

                # Simultaneous bandpass + regression + optional smoothing
                fwhm_str    = str(blur) if blur > 0 else '0'
                prefiltered = opj(dir_prepro_orig_postprocessed,
                                  root_RS + '_prefiltered_func_data_tempfilt.nii.gz')
                cmd = (sing_afni + '3dTproject' +
                       ' -blur '    + fwhm_str +
                       tproject_filter_arg +
                       ' -censor '  + censore1D +
                       ' -cenmode ZERO' +
                       ' -ort '     + regressors_file +
                       ' -prefix '  + prefiltered +
                       ' -input '   + func_input +
                       ' '          + overwrite)
                run_cmd.run(cmd, diary_file)

                if normalize in ('zscore', 'psc'):
                    # normalisation was silently ignored for this method: it is applied to the
                    # residual itself (restoring the mean first would make PSC wrong)
                    cmd = sing_afni + '3dcalc' + overwrite + ' -a ' + prefiltered + ' -prefix ' + residual + ' -expr "a"'
                    run_cmd.run(cmd, diary_file)
                    _write_json(residual, func_input, 'Signal regression (Grandjean: 3dTproject).', cmd)
                    _apply_normalisation(residual, Mean, Sdev, normalize, sing_afni, diary_file,
                                         func_input=func_input, censor_1d=censore1D)
                else:
                    # Restore mean intensity (Grandjean convention)
                    cmd = sing_fsl + 'fslmaths ' + func_input + ' -Tmean ' + Mean
                    run_cmd.run(cmd, diary_file)
                    cmd = (sing_fsl + 'fslmaths ' + prefiltered +
                           ' -add ' + Mean + ' ' + residual)
                    run_cmd.run(cmd, diary_file)
                    _write_json(residual, func_input,
                                'Signal regression (Grandjean: 3dTproject + mean restoration).',
                                cmd)

            # -----------------------------------------------------------------
            # METHOD: fMRIPrep_style
            # Strategy: anatomical CompCor (aCompCor) with fixed or adaptive
            # number of PCs per tissue mask, 6 motion parameters + derivatives
            # + squares (24-parameter Friston model), DCT drift regressors,
            # all projected out simultaneously via 3dTproject.
            #
            # Follows Esteban et al. 2019 (fMRIPrep), Behzadi et al. 2007
            # (CompCor), Friston et al. 1996 (24-parameter motion model).
            #
            # aCompCor mode:
            #   'fixed'    : acompcor_n_fixed PCs per tissue (default 5)
            #   'adaptive' : PCs needed to reach acompcor_variance of the
            #                variance, capped at acompcor_max (aCompCor50-like)
            #
            # Temporal filter: applied via tproject_filter_arg (auto-selected
            # as -passband or -stopband 0 hp).
            # -----------------------------------------------------------------
            elif post_treatment_method == 'fMRIPrep_style':

                regressors = {}   # {name: np.array(T,)} — built up below
                data = _load_4d(func_input)          # loaded once, float32
                n_trs = data.shape[3]

                # ---- Motion: 24-parameter Friston model ---------------------
                # params(6) + derivatives(6) + squares(6) + deriv_squares(6)
                motion_params = np.atleast_2d(np.loadtxt(motion_correction))
                if motion_params.shape != (n_trs, 6):
                    raise ValueError(f'{opb(motion_correction)} has shape {motion_params.shape}, '
                                     f'expected ({n_trs}, 6) for {opb(func_input)}')
                param_names   = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']

                for k, name in enumerate(param_names):
                    regressors[f'motion_{name}'] = motion_params[:, k]

                motion_deriv = np.vstack([np.zeros((1, 6)),
                                          np.diff(motion_params, axis=0)])
                for k, name in enumerate(param_names):
                    regressors[f'motion_{name}_deriv'] = motion_deriv[:, k]

                for k, name in enumerate(param_names):
                    regressors[f'motion_{name}_sq']       = motion_params[:, k] ** 2
                    regressors[f'motion_{name}_deriv_sq'] = motion_deriv[:, k] ** 2

                # ---- aCompCor: WM and CSF separately -----------------------
                # Voxel time series high-passed with the cosine basis of the band
                # before the PCA, as in fMRIPrep (pre_filter='cosine').
                wm_mask  = opj(dir_prepro_orig_masks, 'Wmask.nii.gz')
                csf_mask = opj(dir_prepro_orig_masks, 'Vmask.nii.gz')

                tissue_report = {}
                for tissue_label, mask_path in [('WM', wm_mask),
                                                ('CSF', csf_mask)]:
                    if not ope(mask_path):
                        tissue_report[tissue_label] = dict(mask=opb(mask_path), voxels=0, used=False)
                        run_cmd.msg(
                            f'WARNING: {tissue_label} mask not found — '
                            f'aCompCor skipped for {tissue_label}.',
                            diary_file, 'WARNING')
                        continue
                    if _tissue_too_small(tissue_label, mask_path, root_RS, diary_file, tissue_report):
                        continue
                    # an error here stops the run (it was downgraded to a warning, so a run could
                    # silently lose its aCompCor regressors)
                    if acompcor_mode == 'parallel':
                        comps, var_exp, n_used, capped = _acompcor_parallel(
                            data, mask_path, max_components=acompcor_max,
                            tr_sec=tr_sec, hp_cutoff=hp_cutoff)
                        run_cmd.msg(
                            f'INFO: aCompCor {tissue_label} (parallel analysis) — {n_used} component(s) above '
                            f'chance' + (f', capped at {acompcor_max}' if capped else ''),
                            diary_file, 'WARNING' if capped else 'OKGREEN')
                    else:
                        if acompcor_mode == 'adaptive':
                            comps, var_exp, n_used = _acompcor_adaptive(
                                data, mask_path,
                                variance_threshold=acompcor_variance,
                                max_components=acompcor_max,
                                tr_sec=tr_sec, hp_cutoff=hp_cutoff)
                            reached = float(np.cumsum(var_exp)[-1])
                            run_cmd.msg(
                                f'INFO: aCompCor {tissue_label} (variance criterion) — {n_used} components, '
                                f'cumulative variance {reached:.3f} (target {acompcor_variance}, cap {acompcor_max})',
                                diary_file, 'OKGREEN')
                            if reached < acompcor_variance:
                                run_cmd.msg(
                                    f'WARNING: aCompCor {tissue_label}: the {acompcor_variance:.0%} variance target is '
                                    f'not reached with {acompcor_max} components, so the cap decides the number of '
                                    f'regressors (this mode behaves as a fixed {acompcor_max}). Spatially independent '
                                    f'noise spreads the variance over many components.',
                                    diary_file, 'WARNING')
                        else:
                            comps, var_exp = _acompcor(
                                data, mask_path,
                                n_components=acompcor_n_fixed,
                                tr_sec=tr_sec, hp_cutoff=hp_cutoff)
                            n_used = comps.shape[1]
                            run_cmd.msg(
                                f'INFO: aCompCor {tissue_label} (fixed) — {n_used} components, '
                                f'variance explained: {[f"{v:.3f}" for v in var_exp]}',
                                diary_file, 'OKGREEN')
                        # diagnostic: how many components stand above chance in this mask
                        n_chance = _n_components_parallel(
                            _tissue_timeseries(data, mask_path, tr_sec, hp_cutoff))[0]
                        run_cmd.msg(
                            f'INFO: aCompCor {tissue_label}: parallel analysis finds {n_chance} component(s) above '
                            f'chance ({n_used} used)', diary_file, 'OKCYAN')

                    for k in range(n_used):
                        regressors[f'aCompCor_{tissue_label}_PC{k+1:02d}'] = comps[:, k]

                with open(opj(dir_prepro_orig_postprocessed, root_RS + '_tissue_regressors.json'), 'w') as fh:
                    json.dump(dict(method='fMRIPrep_style', tissues=tissue_report), fh, indent=2)

                # ---- Optional global signal ---------------------------------
                # fMRIPrep definition: mean signal within the brain mask.
                if extract_GS:
                    if not ope(brain_mask):
                        raise FileNotFoundError(f'global signal requested but no brain mask: {brain_mask}')
                    bm = np.asarray(nib.load(brain_mask).dataobj) > 0
                    regressors['global_signal'] = data[bm].mean(axis=0)
                    run_cmd.msg(f'INFO: global signal regressor added (mean in {opb(brain_mask)}).',
                                diary_file, 'OKGREEN')

                # ---- No cosine drift regressors -----------------------------
                # The temporal filter of 3dTproject (tproject_filter_arg) already
                # removes every frequency below hp_cutoff: the DCT set spanned the
                # same frequencies and only used degrees of freedom.

                # ---- Optional custom confounds TSV --------------------------
                confounds_tsv = opj(dir_prepro_orig_postprocessed,
                                    root_RS + '_confounds_correct.tsv')
                if ope(confounds_tsv):
                    extra_df = _read_confounds(confounds_tsv, n_trs, diary_file)
                    run_cmd.msg(
                        f'INFO: loading confounds TSV — columns: '
                        f'{list(extra_df.columns)}',
                        diary_file, 'OKGREEN')
                    for col in extra_df.columns:
                        regressors[f'custom_{col}'] = extra_df[col].values

                _check_dof(_n_kept(censore1D, n_trs), 1 + len(regressors),
                           _n_filter_regressors(n_trs, tr_sec, hp_cutoff, lp_cutoff),
                           f'{root_RS} (fMRIPrep_style)', diary_file)
                del data

                # ---- Save full regressor matrix as TSV (for QC / audit) ----
                regressor_tsv = opj(dir_prepro_orig_postprocessed,
                                    root_RS + '_confounds_fMRIPrep_style.tsv')
                _save_regressors_tsv(regressors, regressor_tsv)
                run_cmd.msg(
                    f'INFO: {len(regressors)} regressors saved → '
                    f'{opb(regressor_tsv)}',
                    diary_file, 'OKGREEN')

                # ---- Write each column as .1D for 3dTproject ----------------
                reg_1d_paths = []
                for reg_name, reg_vals in regressors.items():
                    reg_path = opj(dir_prepro_orig_postprocessed,
                                   root_RS + f'_reg_{reg_name}.1D')
                    np.savetxt(reg_path, reg_vals, fmt='%.8f')
                    reg_1d_paths.append(reg_path)

                # ---- Single 3dTproject step ---------------------------------
                # All confounds projected simultaneously; optional smoothing;
                # temporal filtering via tproject_filter_arg (auto-selected).
                # Censored volumes are zeroed (-cenmode ZERO) so they do not
                # contribute to the regression but are retained as placeholders
                # for downstream TR-counting consistency.
                ort_args = ' '.join(f'-ort {p}' for p in reg_1d_paths)
                fwhm_str = str(blur) if blur > 0 else '0'
                command = (
                    sing_afni + '3dTproject -polort 0' + overwrite +
                    ' -input '   + func_input +
                    ' -censor '  + censore1D +
                    ' -cenmode ZERO' +
                    tproject_filter_arg +
                    ' -blur '    + fwhm_str +
                    ' ' + ort_args +
                    ' -prefix '  + residual)
                run_cmd.msg(f'INFO: 3dTproject command:\n{command}',
                            diary_file, 'OKGREEN')
                run_cmd.run(command, diary_file)
                _write_json(residual, func_input,
                            'Signal regression (fMRIPrep_style: aCompCor + '
                            '24-param motion + temporal filter, 3dTproject).',
                            command)

                # Normalisation
                _apply_normalisation(residual, Mean, Sdev, normalize,
                                     sing_afni, diary_file,
                                     func_input=func_input, censor_1d=censore1D)

            else:
                raise ValueError(
                    f'Unknown post_treatment_method: "{post_treatment_method}". '
                    f'Expected one of: AFNI, Grandjean, fMRIPrep_style.')

        except Exception as exc:
            # The error is logged with its traceback and the run is reported at the end.
            # (A copy of the input named *_3dDeconvolve_failed was written instead: no step
            # uses it, so the failure only showed up later as a missing residual.)
            run_cmd.msg(
                f'ERROR: signal regression failed for run {RS[run_idx]}: {exc}\n'
                + traceback.format_exc(), diary_file, 'FAIL')
            failed_runs.append((RS[run_idx], str(exc)))
            if ope(residual):
                os.remove(residual)          # no stale residual of a previous execution

        finally:
            # Always restore working directory, success or failure
            os.chdir(original_dir)

        # Warn if residual was not produced
        if not ope(residual):
            run_cmd.msg(
                f'WARNING: residual not found after processing run '
                f'{RS[run_idx]} — check logs.',
                diary_file, 'WARNING')

    if failed_runs:
        summary = '; '.join(f'{r}: {e}' for r, e in failed_runs)
        if os.environ.get('EDNIX_STRICT', '0') == '1':
            raise RuntimeError(f'step 7 failed for {len(failed_runs)} run(s): {summary}')
        run_cmd.msg(f'ERROR: step 7 failed for {len(failed_runs)} run(s), which have no residual '
                    f'and will fail in the next steps: {summary}', diary_file, 'FAIL')


# =============================================================================
# Normalisation helper (shared by AFNI and fMRIPrep_style methods)
# =============================================================================

def _apply_normalisation(residual, Mean, Sdev, normalize, sing_afni, diary_file,
                         func_input=None, censor_1d=None):
    """
    Post-regression normalisation, in place on `residual`.

    zscore : (x - mean_kept) / sd_kept, statistics over uncensored volumes only
             (ddof=1); censored volumes are kept at 0.
    psc    : 100 * residual / baseline, baseline = temporal mean of the
             PRE-regression image (`func_input`) over uncensored volumes.
             The residual has zero mean after projection (polort / passband),
             so it cannot serve as its own baseline.
    other  : no normalisation.
    """
    if normalize not in ('zscore', 'psc'):
        run_cmd.msg('INFO: no normalisation applied.', diary_file, 'OKCYAN')
        return

    img = nib.load(residual)
    data = np.asarray(img.dataobj, dtype=np.float32)
    n_t = data.shape[-1]

    if censor_1d is not None and ope(censor_1d):
        keep = np.loadtxt(censor_1d).astype(bool).reshape(-1)
        if keep.size != n_t:
            raise ValueError(f'Censor file has {keep.size} rows, residual has {n_t} volumes: {censor_1d}')
    else:
        keep = np.ones(n_t, dtype=bool)
    if keep.sum() < 3:
        raise ValueError(f'Only {keep.sum()} uncensored volumes in {residual}')

    out = np.zeros_like(data)

    if normalize == 'zscore':
        mu = data[..., keep].mean(axis=-1)
        sd = data[..., keep].std(axis=-1, ddof=1)
        valid = sd > 0
        out[valid] = (data[valid] - mu[valid, None]) / sd[valid, None]
        out[..., ~keep] = 0.0
        nib.save(nib.Nifti1Image(mu, img.affine, img.header), Mean)
        nib.save(nib.Nifti1Image(sd, img.affine, img.header), Sdev)
        desc = 'Z-score over uncensored volumes: (x - mean) / sd (ddof=1); censored volumes = 0.'

    else:  # psc
        if func_input is None or not ope(func_input):
            raise ValueError('PSC requires the pre-regression image (func_input) as baseline.')
        base = np.asarray(nib.load(func_input).dataobj, dtype=np.float32)
        if base.shape != data.shape:
            raise ValueError(f'Shape mismatch residual {data.shape} vs func_input {base.shape}')
        mu = base[..., keep].mean(axis=-1)
        pos = mu[mu > 0]
        floor = 0.01 * np.percentile(pos, 99) if pos.size else np.inf
        valid = mu > floor
        out[valid] = 100.0 * data[valid] / mu[valid, None]
        out[..., ~keep] = 0.0
        nib.save(nib.Nifti1Image(mu, img.affine, img.header), Mean)
        desc = 'PSC: 100 * residual / temporal mean of pre-regression data (uncensored volumes).'

    hdr = img.header.copy()
    hdr.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(out, img.affine, hdr), residual)
    _write_json(residual, [residual] + ([func_input] if func_input else []), desc,
                f'_apply_normalisation(normalize={normalize})')
    run_cmd.msg(f'INFO: {normalize} normalisation applied ({int(keep.sum())}/{n_t} volumes used).',
                diary_file, 'OKGREEN')
