import os
import subprocess
import json
import numpy as np
import ants
import pandas as pd
import nibabel as nib
from scipy import signal as scipy_signal

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


def _acompcor(func_img_path, mask_path, n_components=5):
    """
    Anatomical CompCor: PCA on voxels within a tissue mask.

    Follows the fMRIPrep implementation (Behzadi et al. 2007, Muschelli et al.
    2014): detrend each voxel time series, run PCA on the masked voxels, return
    the first `n_components` principal components as nuisance regressors.

    Parameters
    ----------
    func_img_path : str
        Path to 4D functional NIfTI.
    mask_path : str
        Path to binary tissue mask (WM or CSF).
    n_components : int
        Number of principal components to extract (default 5, as in fMRIPrep).

    Returns
    -------
    components : np.ndarray, shape (n_timepoints, n_components)
        aCompCor regressors, or fewer columns if the mask has < n_components
        voxels.
    variance_explained : np.ndarray
        Fraction of variance explained by each component.
    """
    img   = nib.load(func_img_path)
    data  = img.get_fdata()                      # (X, Y, Z, T)
    mask  = nib.load(mask_path).get_fdata() > 0  # (X, Y, Z)

    # Extract masked voxels → shape (T, n_voxels)
    masked = data[mask].T
    if masked.shape[1] == 0:
        raise ValueError(f"aCompCor: mask {mask_path} contains no voxels.")

    n_comp = min(n_components, masked.shape[1], masked.shape[0] - 1)

    # Linear detrend each voxel time series
    masked = scipy_signal.detrend(masked, axis=0, type='linear')

    # Standardise (zero mean, unit variance) per voxel
    masked -= masked.mean(axis=0)
    std = masked.std(axis=0)
    std[std == 0] = 1.0
    masked /= std

    # SVD (economy)
    U, S, Vt = np.linalg.svd(masked, full_matrices=False)
    components = U[:, :n_comp]
    variance_explained = (S[:n_comp] ** 2) / np.sum(S ** 2)

    return components, variance_explained


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
                      extract_exterior_CSF, extract_WM, normalize, ID,
                      post_treatment_method, do_not_correct_signal, band,
                      extract_Vc, extract_GS,
                      dir_prepro_orig_postprocessed, dir_prepro_raw_matrices,
                      overwrite, sing_afni, sing_fsl, diary_file):

    run_cmd.msg(
        '##  Working on step 7 (function: _7_post_TTT).  ##',
        diary_file, 'HEADER')

    os.makedirs(dir_RS_ICA_native, exist_ok=True)
    os.makedirs(dir_prepro_orig_postprocessed, exist_ok=True)

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
        # First singular vector of voxels within each tissue mask.
        # More robust than the mean against partial-volume effects.
        for extract_flag, img_name, suffix in [
            (extract_exterior_CSF, 'exterior_ligne.nii.gz', '_NonB'),
            (extract_WM,           'Wmask.nii.gz',          '_Wc'),
            (extract_GS,           'maskDilat.nii.gz',      '-GS'),
        ]:
            if not extract_flag:
                continue
            out_1d = opj(dir_prepro_orig_postprocessed,
                         root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.1D')
            command = (sing_afni + '3dmaskSVD' + overwrite +
                       ' -polort 2 -vnorm -mask ' +
                       opj(dir_prepro_orig_masks, img_name) +
                       ' ' + func_input)
            val, _ = run_cmd.get(command, diary_file)
            _write_1d(val, out_1d)

        # ---- Ventricle signal (size-gated) ----------------------------------
        if extract_Vc:
            vmask_path = opj(dir_prepro_orig_masks, 'Vmask.nii.gz')
            msk        = ants.image_read(vmask_path)
            if int(msk.sum()) > 10:
                out_1d = opj(dir_prepro_orig_postprocessed,
                             root_RS + '_space-acpc-func_desc-fMRI_run_inRef_Vc.1D')
                command = (sing_afni + '3dmaskSVD' + overwrite +
                           ' -polort 2 -vnorm -mask ' + vmask_path +
                           ' ' + func_input)
                val, _ = run_cmd.get(command, diary_file)
                _write_1d(val, out_1d)
            else:
                run_cmd.msg(
                    'WARNING: Vmask has ≤10 voxels — ventricle regressor skipped.',
                    diary_file, 'WARNING')

        # ---- Bandpass regressors (sinusoidal basis, for GLM approach) -------
        # Using regressors rather than direct filtering preserves DOF around
        # censored TRs (Power et al. 2014). This is the AFNI / fMRIPrep-aligned
        # approach; the Grandjean method uses 3dTproject -passband directly.
        hd = ants.image_header_info(func_input)
        n_trs = int(hd['dimensions'][3])
        command = (sing_afni + '1dBport' + overwrite +
                   ' -nodata ' + str(n_trs) +
                   ' ' + str(TR) +
                   ' -band ' + band +
                   ' -invert -nozero')
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

                # Add optional tissue regressors
                for extract_flag, suffix in [
                    (extract_exterior_CSF, '_NonB'),
                    (extract_WM,           '_Wc'),
                    (extract_GS,           '-GS'),
                    (extract_Vc,           '_Vc'),
                ]:
                    if extract_flag:
                        reg_path = opj(dir_prepro_orig_postprocessed,
                                       root_RS + '_space-acpc-func_desc-fMRI_run_inRef' + suffix + '.1D')
                        command += f' -ortvec {reg_path} residual_norm{suffix}'

                # Add custom confounds TSV columns if present
                confounds_tsv = opj(dir_prepro_orig_postprocessed,
                                    root_RS + '_confounds_correct.tsv')
                if ope(confounds_tsv):
                    confounds_df = pd.read_csv(confounds_tsv, sep='\t')
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

                run_cmd.msg(f'INFO: 3dDeconvolve command:\n{command}',
                            diary_file, 'OKGREEN')
                subprocess.run(command, shell=True, check=True)

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
                                     sing_afni, diary_file)

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
                for mask_path, sig_path, label in [
                    (wm_mask,  wm_signal,  'WM'),
                    (csf_mask, csf_signal, 'CSF'),
                ]:
                    if ope(mask_path):
                        run_cmd.msg(f'INFO: extracting {label} signal',
                                    diary_file, 'OKCYAN')
                        cmd = (sing_fsl + 'fslmeants'
                               ' -i ' + func_input +
                               ' -o ' + sig_path +
                               ' -m ' + mask_path)
                        run_cmd.run(cmd, diary_file)
                        tissue_signals.append(sig_path)
                    else:
                        run_cmd.msg(f'WARNING: {label} mask not found — skipped.',
                                    diary_file, 'WARNING')

                run_cmd.msg(
                    f'INFO: concatenating regressors: '
                    f'{[opb(p) for p in tissue_signals]}',
                    diary_file, 'OKCYAN')
                subprocess.run(
                    sing_afni + '1dcat ' +
                    ' '.join(tissue_signals) + ' > ' + regressors_file,
                    shell=True, check=True)

                # Simultaneous bandpass + regression + optional smoothing
                fwhm_str    = str(blur) if blur > 0 else '0'
                prefiltered = opj(dir_prepro_orig_postprocessed,
                                  root_RS + '_prefiltered_func_data_tempfilt.nii.gz')
                cmd = (sing_afni + '3dTproject' +
                       ' -blur '    + fwhm_str +
                       ' -passband ' + band +
                       ' -ort '     + regressors_file +
                       ' -prefix '  + prefiltered +
                       ' -input '   + func_input +
                       ' '          + overwrite)
                run_cmd.run(cmd, diary_file)

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
            # Strategy: anatomical CompCor (aCompCor) with 5 PCs from WM and
            # CSF masks separately, 6 motion parameters + derivatives + squares
            # (24-parameter Friston model), bandpass via sinusoidal regressors,
            # all projected out simultaneously via 3dTproject.
            #
            # Follows Esteban et al. 2019 (fMRIPrep), Behzadi et al. 2007
            # (CompCor), Friston et al. 1996 (24-parameter motion model).
            #
            # Key differences from AFNI method:
            #   • aCompCor (5 PCs per tissue) instead of single SVD component
            #   • 24-parameter motion model (params + deriv + squares)
            #   • All regressors assembled in Python then passed to 3dTproject
            #     → single projection step, avoids repeated interpolation
            #   • No -polort A (polynomial detrending folded into regressors)
            # -----------------------------------------------------------------
            elif post_treatment_method == 'fMRIPrep_style':

                regressors = {}   # {name: np.array(T,)} — built up below

                # ---- Motion: 24-parameter Friston model ---------------------
                # params(6) + derivatives(6) + squares(6) + deriv_squares(6)
                motion_params = np.loadtxt(motion_correction)  # (T, 6)
                param_names   = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz']

                for k, name in enumerate(param_names):
                    regressors[f'motion_{name}']       = motion_params[:, k]

                motion_deriv = np.vstack([np.zeros((1, 6)),
                                          np.diff(motion_params, axis=0)])
                for k, name in enumerate(param_names):
                    regressors[f'motion_{name}_deriv'] = motion_deriv[:, k]

                for k, name in enumerate(param_names):
                    regressors[f'motion_{name}_sq']       = motion_params[:, k] ** 2
                    regressors[f'motion_{name}_deriv_sq'] = motion_deriv[:, k] ** 2

                # ---- aCompCor: WM and CSF separately, 5 PCs each ------------
                wm_mask  = opj(dir_prepro_orig_masks, 'Wmask.nii.gz')
                csf_mask = opj(dir_prepro_orig_masks, 'Vmask.nii.gz')

                for tissue_label, mask_path in [('WM', wm_mask),
                                                ('CSF', csf_mask)]:
                    if not ope(mask_path):
                        run_cmd.msg(
                            f'WARNING: {tissue_label} mask not found — '
                            f'aCompCor skipped for {tissue_label}.',
                            diary_file, 'WARNING')
                        continue

                    try:
                        comps, var_exp = _acompcor(func_input, mask_path,
                                                   n_components=5)
                        run_cmd.msg(
                            f'INFO: aCompCor {tissue_label} — '
                            f'{comps.shape[1]} components, '
                            f'variance explained: '
                            f'{[f"{v:.3f}" for v in var_exp]}',
                            diary_file, 'OKGREEN')
                        for k in range(comps.shape[1]):
                            regressors[f'aCompCor_{tissue_label}_PC{k+1:02d}'] = comps[:, k]
                    except Exception as exc:
                        run_cmd.msg(
                            f'WARNING: aCompCor failed for {tissue_label}: {exc}',
                            diary_file, 'WARNING')

                # ---- Optional global signal ---------------------------------
                if extract_GS:
                    gs_path = opj(dir_prepro_orig_masks, 'maskDilat.nii.gz')
                    if ope(gs_path):
                        command = (sing_afni + '3dmaskSVD' + overwrite +
                                   ' -polort 2 -vnorm -mask ' + gs_path +
                                   ' ' + func_input)
                        val, _ = run_cmd.get(command, diary_file)
                        gs_lines = val.decode('utf-8').strip().split('\n')
                        gs_vals  = np.array([float(l.strip())
                                             for l in gs_lines if l.strip()])
                        regressors['global_signal'] = gs_vals
                        run_cmd.msg('INFO: global signal regressor added.',
                                    diary_file, 'OKGREEN')

                # ---- Cosine drift regressors (DCT set, fMRIPrep style) ------
                # Discrete cosine transform basis replaces -polort; captures
                # low-frequency drift without assuming polynomial form.
                # High-pass cutoff = band.split()[0] Hz (lower edge of passband).
                n_trs      = motion_params.shape[0]
                try:
                    hp_cutoff  = float(band.split()[0])
                    period_max = 1.0 / hp_cutoff       # seconds
                    n_cosines  = int(np.floor(2 * n_trs * TR / period_max))
                    if n_cosines > 0:
                        t = np.arange(n_trs)
                        for k in range(1, n_cosines + 1):
                            cosine = np.cos(np.pi * k * (2 * t + 1) / (2 * n_trs))
                            regressors[f'cosine_{k:02d}'] = cosine
                        run_cmd.msg(
                            f'INFO: {n_cosines} DCT drift regressors added '
                            f'(HP cutoff {hp_cutoff} Hz).',
                            diary_file, 'OKGREEN')
                except (ValueError, IndexError):
                    run_cmd.msg(
                        f'WARNING: could not parse HP cutoff from band="{band}" '
                        f'— DCT regressors skipped.',
                        diary_file, 'WARNING')

                # ---- Optional custom confounds TSV --------------------------
                confounds_tsv = opj(dir_prepro_orig_postprocessed,
                                    root_RS + '_confounds_correct.tsv')
                if ope(confounds_tsv):
                    extra_df = pd.read_csv(confounds_tsv, sep='\t')
                    run_cmd.msg(
                        f'INFO: loading confounds TSV — columns: '
                        f'{list(extra_df.columns)}',
                        diary_file, 'OKGREEN')
                    for col in extra_df.columns:
                        regressors[f'custom_{col}'] = extra_df[col].values

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
                # bandpass via -passband (equivalent to sinusoidal regressors
                # when no censoring, but simpler to implement here since all
                # regressors are explicit).
                # NOTE: censored volumes are zeroed (-cenmode ZERO) so they
                # do not contribute to the regression but are retained as
                # placeholders for downstream TR-counting consistency.
                ort_args = ' '.join(f'-ort {p}' for p in reg_1d_paths)
                fwhm_str = str(blur) if blur > 0 else '0'
                command = (
                    sing_afni + '3dTproject -polort 0' + overwrite +
                    ' -input '   + func_input +
                    ' -censor '  + censore1D +
                    ' -cenmode ZERO'
                    ' -passband ' + band +
                    ' -blur '    + fwhm_str +
                    ' ' + ort_args +
                    ' -prefix '  + residual)
                run_cmd.msg(f'INFO: 3dTproject command:\n{command}',
                            diary_file, 'OKGREEN')
                run_cmd.run(command, diary_file)
                _write_json(residual, func_input,
                            'Signal regression (fMRIPrep_style: aCompCor + '
                            '24-param motion + DCT drift + 3dTproject).',
                            command)

                # Normalisation
                _apply_normalisation(residual, Mean, Sdev, normalize,
                                     sing_afni, diary_file)

            else:
                raise ValueError(
                    f'Unknown post_treatment_method: "{post_treatment_method}". '
                    f'Expected one of: AFNI, Grandjean, fMRIPrep_style.')

        except Exception as exc:
            run_cmd.msg(
                f'ERROR: signal regression failed for run {RS[run_idx]}: {exc}',
                diary_file, 'WARNING')

            # Fall back: copy input to residual so pipeline can continue
            deconvolve_failed = opj(
                dir_prepro_orig_process,
                root_RS + '_space-acpc-func_desc-3dDeconvolve_failed.nii.gz')
            command = (sing_afni + '3dcalc' + overwrite +
                       ' -a ' + func_input +
                       ' -prefix ' + deconvolve_failed + ' -expr "a"')
            run_cmd.do(command, diary_file)
            _write_json(deconvolve_failed, func_input,
                        f'Fallback copy — signal regression failed: {exc}',
                        command)

        finally:
            # Always restore working directory, success or failure
            os.chdir(original_dir)

        # Warn if residual was not produced
        if not ope(residual):
            run_cmd.msg(
                f'WARNING: residual not found after processing run '
                f'{RS[run_idx]} — check logs.',
                diary_file, 'WARNING')


# =============================================================================
# Normalisation helper (shared by AFNI and fMRIPrep_style methods)
# =============================================================================

def _apply_normalisation(residual, Mean, Sdev, normalize, sing_afni, diary_file):
    """
    Apply post-regression normalisation in-place on the residual image.

    Parameters
    ----------
    residual : str
        Path to residual image (modified in-place).
    Mean : str
        Output path for temporal mean image.
    Sdev : str
        Output path for temporal stdev image (zscore only).
    normalize : str
        'zscore', 'psc', or anything else (no normalisation).
    """
    if normalize == 'zscore':
        # (x - mean) / stdev  — zero mean, unit variance per voxel
        command = (
            sing_afni + '3dTstat -mean -overwrite -prefix ' + Mean +
            ' ' + residual +
            ' && ' + sing_afni + '3dTstat -stdev -overwrite -prefix ' + Sdev +
            ' ' + residual +
            ' && ' + sing_afni + '3dcalc -overwrite'
            ' -a ' + residual +
            ' -b ' + Mean +
            ' -c ' + Sdev +
            ' -expr "(a-b)/c" -prefix ' + residual)
        run_cmd.do(command, diary_file)
        _write_json(residual, residual,
                    'Z-score normalisation: (x - mean) / stdev.', command)

    elif normalize == 'psc':
        # percent signal change: ((x - mean) / mean) * 100
        command = (
            sing_afni + '3dTstat -mean -overwrite -prefix ' + Mean +
            ' ' + residual +
            ' && ' + sing_afni + '3dcalc -overwrite'
            ' -a ' + residual +
            ' -b ' + Mean +
            ' -expr "((a-b)/b)*100" -prefix ' + residual)
        run_cmd.do(command, diary_file)
        _write_json(residual, residual,
                    'PSC normalisation: ((x - mean) / mean) × 100.', command)

    else:
        run_cmd.msg('INFO: no normalisation applied.', diary_file, 'OKCYAN')
