import os
import nilearn
from nilearn import image
from nilearn.input_data import NiftiLabelsMasker
import shutil
from shutil import copyfile
import math
import scipy
from scipy.stats import pearsonr, entropy
from nitime.lazy import scipy_linalg as linalg
import nitime.utils as utils
import nibabel as nib
import pandas as pd
import json
import time
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.metrics import mutual_info_score
from sklearn.preprocessing import StandardScaler
import ants
import numpy as np
from atlases import atlas4func

opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists
ops = os.path.splitext
opa = os.path.abspath

from Tools import run_cmd
from fonctions.extract_filename import extract_filename

def create_surrogate_from_label_atlas(anat_img,diary_file):
    """Create a surrogate T1/T2-like image from a label atlas using inversion trick with safer handling."""
    atlas_data = anat_img.numpy()

    valid_mask = atlas_data > 0
    if not np.any(valid_mask):
        nl = "Input label atlas contains no non-zero values."
        run_cmd.msg(nl, diary_file, 'WARNING')
        raise ValueError(nl)

    # Get only valid data
    atlas_data_valid = atlas_data[valid_mask]
    unique_values = np.unique(atlas_data_valid)

    if len(unique_values) < 2:
        # Instead of failing, fallback to flat image with small noise
        nl="WARNING: Only one label found, generating flat surrogate image with noise."
        run_cmd.msg(nl,diary_file,'WARNING')
        surrogate_data = np.zeros_like(atlas_data, dtype=np.float32)
        surrogate_data[valid_mask] = 50 + np.random.normal(0, 1, size=np.sum(valid_mask))
    else:
        # Safe linear inversion
        max_val = atlas_data_valid.max()
        min_val = atlas_data_valid.min()
        R = 100 / (max_val - min_val)
        surrogate_data = (max_val - atlas_data) * R
        surrogate_data = np.where(valid_mask, surrogate_data, 0.0)

    # Convert back to ANTs
    surrogate_img = ants.from_numpy(surrogate_data.astype(np.float32),
                                    origin=anat_img.origin,
                                    spacing=anat_img.spacing,
                                    direction=anat_img.direction)
    return surrogate_img


def detect_image_type(anat_img, atlas_img,diary_file):
    """Detect whether anatomical image is T1 or T2 weighted based on tissue contrast."""
    try:
        # Get numpy arrays
        anat_data = anat_img.numpy()
        atlas_data = atlas_img.numpy()
        nonzero = anat_data > 0  # Create mask of non-zero values

        # Define regions of interest (using simplified atlas values)
        wm_mask = atlas_data == 6  # White matter
        gm_mask = (atlas_data == 8) #| (atlas_data == 9) | (atlas_data == 10)  | (atlas_data == 11) | (atlas_data == 12)# Gray regions
        csf_mask = (atlas_data == 1) | (atlas_data == 2) | (atlas_data == 3)

        # Calculate mean intensities
        wm_mean = np.mean(anat_data[wm_mask & nonzero]) if np.any(wm_mask & nonzero) else np.nan
        gm_mean = np.mean(anat_data[gm_mask & nonzero]) if np.any(gm_mask & nonzero) else np.nan
        csf_mean = np.mean(anat_data[csf_mask & nonzero]) if np.any(csf_mask & nonzero) else np.nan

        # Calculate contrast ratios
        if not np.isnan(wm_mean) and not np.isnan(gm_mean) and not np.isnan(csf_mean):
            wm_gm_ratio = wm_mean / gm_mean
            wm_csf_ratio = wm_mean / csf_mean
            gm_csf_ratio = gm_mean / csf_mean

            # Typical T1 contrasts: WM > GM > CSF
            t1_likelihood = 0
            if wm_mean > gm_mean: t1_likelihood += 1

            # Typical T2 contrasts: CSF > GM > WM
            t2_likelihood = 0
            if gm_mean > wm_mean: t2_likelihood += 1
            nl = "t1_likelihood " + str(t1_likelihood)
            run_cmd.msg(nl, diary_file, 'ENDC')
            nl = "t2_likelihood " + str(t2_likelihood)
            run_cmd.msg(nl, diary_file, 'ENDC')
            # Determine image type
            if t1_likelihood > t2_likelihood:
                image_type = "T1"
            elif t2_likelihood > t1_likelihood:
                image_type = "T2"
            else:
                image_type = "Unknown"

            return {
                "image_type": image_type,
                "wm_mean": wm_mean,
                "gm_mean": gm_mean,
                "csf_mean": csf_mean,
                "wm_gm_ratio": wm_gm_ratio,
                "wm_csf_ratio": wm_csf_ratio,
                "gm_csf_ratio": gm_csf_ratio,
                "t1_likelihood": t1_likelihood,
                "t2_likelihood": t2_likelihood
            }
        else:
            return {
                "image_type": "Unknown",
                "wm_mean": np.nan,
                "gm_mean": np.nan,
                "csf_mean": np.nan,
                "wm_gm_ratio": np.nan,
                "wm_csf_ratio": np.nan,
                "gm_csf_ratio": np.nan,
                "t1_likelihood": np.nan,
                "t2_likelihood": np.nan
            }
    except Exception as e:
        nl = 'ERROR in detect_image_type: ' + str(e)
        run_cmd.msg(nl,diary_file,'FAIL')
        return {
            "image_type": "Unknown",
            "wm_mean": np.nan,
            "gm_mean": np.nan,
            "csf_mean": np.nan,
            "wm_gm_ratio": np.nan,
            "wm_csf_ratio": np.nan,
            "gm_csf_ratio": np.nan,
            "t1_likelihood": np.nan,
            "t2_likelihood": np.nan
        }

def compute_coverage_stats(anat_img, fmri_img):
    # Load images
    anat_img = anat_img.numpy()
    fmri_img = fmri_img.numpy()

    # Generate masks
    anat_mask = anat_img != 0
    fmri_mask = fmri_img != 0

    # Find uncovered voxels
    anat_uncovered = anat_mask & ~fmri_mask  # Anat non-zero not covered by fMRI
    fmri_outside_anat = fmri_mask & ~anat_mask  # fMRI non-zero outside anat

    # Compute stats
    anat_total = np.count_nonzero(anat_mask)
    func_total = np.count_nonzero(fmri_mask)
    anat_uncovered_count = np.count_nonzero(anat_uncovered)
    fmri_outside_count = np.count_nonzero(fmri_outside_anat)

    uncovered_pct_anat = (anat_uncovered_count / anat_total) * 100 if anat_total > 0 else 0
    uncovered_pct_func = (fmri_outside_count / func_total) * 100 if func_total > 0 else 0
    return {
        "pct_anat_uncovered": uncovered_pct_anat,
        "pct_fmri_uncovered": uncovered_pct_func}

def create_qc_figure(output_results, root_RS, snr_results, fd_jenkinson_array,
                     calc_dvars_array, motion_params, image_metrics, qc_values,
                     anat_img_path, fmri_img_path):
    """
    Create a comprehensive QC figure with multiple subplots including:
    - SNR plots
    - Motion parameters
    - DVARS/FD
    - Coverage overlay
    - Comprehensive metrics table
    """
    try:
        # Create figure with constrained layout
        fig = plt.figure(figsize=(20, 24), constrained_layout=True)
        plt.style.use('seaborn' if 'seaborn' in plt.style.available else 'ggplot')

        # Set font sizes
        plt.rcParams.update({
            'font.size': 10,
            'axes.titlesize': 12,
            'axes.labelsize': 10,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 9
        })

        # Create grid layout (4 rows, 3 columns)
        gs = GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 0.7])

        # Subplot A: Regional SNR Over Time
        ax1 = fig.add_subplot(gs[0, 0])
        if 'Gray_Matter' in snr_results:
            for region in snr_results['Gray_Matter'].keys():
                left_data = snr_results['Gray_Matter'][region].get('Left', [])
                right_data = snr_results['Gray_Matter'][region].get('Right', [])
                if len(left_data) > 0 and not all(np.isnan(x) for x in left_data):
                    ax1.plot(left_data, label=f"{region} Left", linestyle='--', alpha=0.7)
                if len(right_data) > 0 and not all(np.isnan(x) for x in right_data):
                    ax1.plot(right_data, label=f"{region} Right", alpha=0.7)
        ax1.set_xlabel('Time (volumes)')
        ax1.set_ylabel('SNR')
        ax1.set_title('A: Regional SNR Over Time')
        ax1.grid(True, linestyle=':', alpha=0.5)
        ax1.legend()

        # Subplot B: Gray vs White Matter SNR
        ax2 = fig.add_subplot(gs[0, 1])
        gray_data = snr_results.get('Avg_Gray_SNR', [])
        white_data = snr_results.get('Avg_White_SNR', [])
        if len(gray_data) > 0 and not all(np.isnan(x) for x in gray_data):
            ax2.plot(gray_data, label='Gray Matter', color='blue')
        if len(white_data) > 0 and not all(np.isnan(x) for x in white_data):
            ax2.plot(white_data, label='White Matter', color='green')
        ax2.set_xlabel('Time (volumes)')
        ax2.set_ylabel('Average SNR')
        ax2.set_title('B: Gray vs White Matter SNR')
        ax2.legend()
        ax2.grid(True, linestyle=':', alpha=0.5)

        # Subplot C: Framewise Displacement
        ax3 = fig.add_subplot(gs[0, 2])
        if len(fd_jenkinson_array) > 0:
            timepoints = np.arange(len(fd_jenkinson_array))
            ax3.plot(timepoints, fd_jenkinson_array, label='FD Jenkinson', color='b')
            ax3.axhline(y=0.5, color='r', linestyle='--', label='Threshold (0.5mm)')
            ax3.set_xlabel('Time (volumes)')
            ax3.set_ylabel('Framewise Displacement (mm)')
            ax3.set_title('C: Framewise Displacement')
            ax3.legend()
            ax3.grid(True, linestyle=':', alpha=0.5)

        # Subplot D: DVARS
        ax4 = fig.add_subplot(gs[1, 0])
        if len(calc_dvars_array) > 0:
            timepoints = np.arange(len(calc_dvars_array))
            ax4.plot(timepoints, calc_dvars_array, label='DVARS', color='b')
            ax4.set_xlabel('Time (volumes)')
            ax4.set_ylabel('DVARS')
            ax4.set_title('D: DVARS Over Time')
            ax4.grid(True, linestyle=':', alpha=0.5)

        # Subplot E: Motion Parameters
        ax5 = fig.add_subplot(gs[1, 1:3])
        if motion_params is not None and motion_params.size > 0 and motion_params.ndim == 2:
            labels = ['X', 'Y', 'Z', 'Pitch', 'Roll', 'Yaw']
            colors = plt.cm.tab10(np.linspace(0, 1, 6))
            for i in range(min(6, motion_params.shape[1])):
                ax5.plot(motion_params[:, i], label=labels[i], color=colors[i])
            ax5.set_xlabel('Time (volumes)')
            ax5.set_ylabel('Motion (mm/degrees)')
            ax5.set_title('E: Motion Parameters')
            ax5.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax5.grid(True, linestyle=':', alpha=0.5)

        # Subplot F: Image Intensity Distribution
        ax6 = fig.add_subplot(gs[2, 0])
        vals1 = image_metrics.get('vals1', [])
        vals2 = image_metrics.get('vals2', [])
        if len(vals1) > 0 and len(vals2) > 0:
            ax6.hist(vals1, bins=50, alpha=0.5,
                     label='Template', density=True, color='blue')
            ax6.hist(vals2, bins=50, alpha=0.5,
                     label='fMRI', density=True, color='orange')
            ax6.set_xlabel('Intensity')
            ax6.set_ylabel('Density')
            ax6.set_title('F: Intensity Distribution')
            ax6.legend()
            ax6.grid(True, linestyle=':', alpha=0.5)

        # Subplot G: Image Correlation
        ax7 = fig.add_subplot(gs[2, 1])
        if len(vals1) > 0 and len(vals2) > 0:
            ax7.scatter(vals2[::10], vals1[::10],
                        alpha=0.3, s=5, color='green')
            cc = image_metrics.get('cc', np.nan)
            ax7.set_xlabel('fMRI Intensity')
            ax7.set_ylabel('Template Intensity')
            ax7.set_title(f'G: Intensity Correlation (r={cc:.2f})')
            ax7.grid(True, linestyle=':', alpha=0.5)

        # Subplot H: Coverage Overlay (using nilearn)
        ax8 = fig.add_subplot(gs[2, 2])
        if anat_img_path is not None and fmri_img_path is not None:
            try:
                from nilearn import plotting
                display = plotting.plot_anat(anat_img_path, axes=ax8,
                                             title='H: Anat-fMRI Coverage',
                                             draw_cross=False)
                display.add_contours(fmri_img_path, colors='r', linewidths=0.5)

                # Add coverage text
                coverage_text = (
                    f"fMRI uncovered: {qc_values.get('pct_fmri_uncovered', np.nan):.1f}%\n"
                    f"Anat uncovered: {qc_values.get('pct_anat_uncovered', np.nan):.1f}%"
                )
                ax8.text(0.05, 0.05, coverage_text,
                         transform=ax8.transAxes, color='white',
                         bbox=dict(facecolor='black', alpha=0.5))
            except Exception as e:
                print(f"ERROR in coverage visualization: {str(e)}")
                ax8.axis('off')
                ax8.set_title('H: Coverage (Failed to render)')
        else:
            ax8.axis('off')
            ax8.set_title('H: Coverage (No data)')

        # Subplot I: Comprehensive Metrics Table
        ax9 = fig.add_subplot(gs[3, :])
        ax9.axis('off')

        # Safe formatting function
        def safe_format(value, format_spec=".2f", default="N/A", unit=""):
            if value is None:
                return f"{default}{unit}"
            try:
                if isinstance(value, (np.ndarray, list)):
                    if len(value) == 0:
                        return f"{default}{unit}"
                    value = np.nanmean(value)
                float_val = float(value)
                return f"{float_val:{format_spec}}{unit}" if not np.isnan(float_val) else f"{default}{unit}"
            except (TypeError, ValueError):
                return f"{default}{unit}"

        # Organize metrics into categories
        metrics_data = [
            ["SNR Metrics",
             f"Gray Matter SNR: {safe_format(qc_values.get('avg_snr_gray'), '.1f')}",
             f"White Matter SNR: {safe_format(qc_values.get('avg_snr_white'), '.1f')}",
             f"Stdev: {safe_format(qc_values.get('stdev'))}",
             f"TSNR: {safe_format(qc_values.get('TSNR'))}",
             f"TSNRcvarinv: {safe_format(qc_values.get('TSNRcvarinv'))}"],

            ["Motion Metrics",
             f"Mean FD: {safe_format(qc_values.get('mean_fd'), '.3f', unit=' mm')}",
             f"Mean DVARS: {safe_format(qc_values.get('mean_dvars'), '.1f')}",
             f"Motion Velocity: {safe_format(qc_values.get('avg_velocity'))}",
             f"Censor Fraction: {safe_format(qc_values.get('censor_fraction'))}",
             f"Outliers: {safe_format(qc_values.get('avg_outcount'))}"],

            ["Image Quality",
             f"Ghost Ratio: {safe_format(qc_values.get('ghost_ratio'))}",
             f"FWHM: {safe_format(qc_values.get('fwhm'), '.2f', unit=' mm')}",
             f"Global Correlation: {safe_format(qc_values.get('gcor'))}",
             f"Mutual Info: {safe_format(qc_values.get('mi'))}",
             f"NMI: {safe_format(qc_values.get('nmi'))}"],

            ["Coverage & Image Type",
             f"Comp. modal.: {qc_values.get('Comp modal', 'N/A')}",
             f"fMRI not covered: {safe_format(qc_values.get('pct_fmri_uncovered'), '.1f', unit='%')}",
             f"anat not covered: {safe_format(qc_values.get('pct_anat_uncovered'), '.1f', unit='%')}",
             f"WM/GM Contrast: {safe_format(qc_values.get('wm_gm_contrast'))}",
             f"WM/CSF Contrast: {safe_format(qc_values.get('wm_csf_contrast'))}"]
        ]

        # Create table
        table = ax9.table(
            cellText=metrics_data,
            loc='center',
            cellLoc='left',
            colWidths=[0.2] * 6
        )

        # Style table
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.8)

        # Add title to the table
        ax9.set_title('I: Summary of QC Metrics', y=1.1)

        # Save figure
        plt.savefig(
            opj(output_results, root_RS + '_QC_summary.png'),
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()

    except Exception as e:
        print(f"ERROR in create_qc_figure: {str(e)}")
        raise
def fd_jenkinson(in_file, rmax=80., out_file=None, out_array=True):
    """Calculate Jenkinson's Framewise Displacement."""
    if out_file is None:
        fname, ext = ops(opb(in_file))
        out_file = opa('%s_fdfile%s' % (fname, ext))

    if 'rel.rms' in in_file:
        copyfile(in_file, out_file)
        return out_file

    try:
        pm_ = np.genfromtxt(in_file)
    except:
        raise Exception('ERROR reading motion file')

    original_shape = pm_.shape
    pm = np.zeros((pm_.shape[0], pm_.shape[1] + 4))
    pm[:, :original_shape[1]] = pm_
    pm[:, original_shape[1]:] = [0.0, 0.0, 0.0, 1.0]

    T_rb_prev = np.matrix(np.eye(4))
    flag = 0
    X = [0]

    for i in range(pm.shape[0]):
        T_rb = np.matrix(pm[i].reshape(4, 4))

        if flag == 0:
            flag = 1
        else:
            M = np.dot(T_rb, T_rb_prev.I) - np.eye(4)
            A = M[0:3, 0:3]
            b = M[0:3, 3]
            FD_J = math.sqrt((rmax * rmax / 5) * np.trace(np.dot(A.T, A)) + np.dot(b.T, b))
            X.append(FD_J)

        T_rb_prev = T_rb

    try:
        X = np.array(X).reshape(-1)
        np.savetxt(out_file, X)
    except:
        raise Exception('ERROR saving FD file')

    return np.array(X).reshape(-1) if out_array else out_file


def global_correlation(func_reorient, func_mask):
    """Calculate global correlation (GCOR)."""
    func = load(func_reorient, func_mask)
    list_of_ts = func.transpose()
    demeaned_normed = [scipy.stats.mstats.zscore(ts) for ts in list_of_ts]
    demeaned_normed = np.asarray(demeaned_normed)
    volume_list = demeaned_normed.transpose()
    avg_ts = np.asarray([voxel.mean() for voxel in volume_list])
    return (avg_ts.transpose().dot(avg_ts)) / len(avg_ts)


def load(func_file, mask_file):
    """Load functional data with mask."""
    func_img = nib.load(func_file)
    mask_img = nib.load(mask_file)
    mask = mask_img.get_fdata()
    mask_var_filtered = mask
    func = func_img.get_fdata().astype(float)
    variance_map = func.var(axis=-1)
    mask_var_filtered[variance_map < 1e-6] = 0
    return func[mask_var_filtered.nonzero()].T


def calc_dvars(func_file, mask_file, output_all=False):
    """Calculate DVARS metric."""
    func = load(func_file, mask_file)
    func_sd = robust_stdev(func)
    func_ar1 = ar1(func)
    func_sd_pd = np.sqrt(2 * (1 - func_ar1)) * func_sd
    diff_sd_mean = func_sd_pd.mean()
    func_deriv = np.diff(func, axis=0)
    dvars_stdz = func_deriv.std(1, ddof=1) / diff_sd_mean
    return dvars_stdz.reshape(len(dvars_stdz), 1)


def robust_stdev(func):
    """Compute robust standard deviation."""
    lower_qs = np.percentile(func, 25, axis=0)
    upper_qs = np.percentile(func, 75, axis=0)
    return (upper_qs - lower_qs) / 1.349


def ar_nitime(x, order=1, center=False):
    """AR model estimation."""
    if center:
        x = x.copy()
        x = x - x.mean()
    r_m = utils.autocorr(x)[:order + 1]
    Tm = linalg.toeplitz(r_m[:order])
    y = r_m[1:]
    return linalg.solve(Tm, y)[0]


def ar1(func, method=ar_nitime):
    """Apply AR model across timeseries."""
    func_centered = func - func.mean(0)
    return np.apply_along_axis(method, 0, func_centered)


def ghost_direction(epi_data_path, mask_data_path, direction="y", ref_file=None, out_file=None):
    """Calculate Ghost to Signal Ratio."""
    epi_img = nib.load(epi_data_path)
    epi_data = epi_img.get_fdata()
    mask_img = nib.load(mask_data_path)
    mask_data = mask_img.get_fdata()
    n2_mask_data = np.zeros_like(mask_data)

    if direction == "x":
        n2 = int(np.floor(mask_data.shape[0] / 2))
        n2_mask_data[:n2, :, :] = mask_data[n2:(n2 * 2), :, :]
        n2_mask_data[n2:(n2 * 2), :, :] = mask_data[:n2, :, :]
    elif direction == "y":
        n2 = int(np.floor(mask_data.shape[1] / 2))
        n2_mask_data[:, :n2, :] = mask_data[:, n2:(n2 * 2), :]
        n2_mask_data[:, n2:(n2 * 2), :] = mask_data[:, :n2, :]
    elif direction == "z":
        n2 = int(np.floor(mask_data.shape[2] / 2))
        n2_mask_data[:, :, :n2] = mask_data[:, :, n2:(n2 * 2)]
        n2_mask_data[:, :, n2:(n2 * 2)] = mask_data[:, :, :n2]
    else:
        raise Exception("Unknown direction %s" % direction)

    n2_mask_data = n2_mask_data * (1 - mask_data)
    n2_mask_data = n2_mask_data + 2 * (1 - n2_mask_data - mask_data)

    if ref_file is not None and out_file is not None:
        ref = nib.load(ref_file)
        out = nib.Nifti1Image(n2_mask_data, ref.affine, ref.header)
        out.to_filename(out_file)

    return ((epi_data[n2_mask_data == 1].mean() - epi_data[n2_mask_data == 2].mean()) /
            epi_data[n2_mask_data == 0].mean())


def safe_format(value, format_spec=".2f", default="N/A", unit=""):
    """Safely format values that might be numpy types, missing, or invalid"""
    if value is None:
        return f"{default}{unit}"
    try:
        # Handle numpy arrays/lists by taking mean
        if isinstance(value, (np.ndarray, list)):
            if len(value) == 0:
                return f"{default}{unit}"
            value = np.nanmean(value)
        # Convert numpy types to native float
        float_val = float(value)
        return f"{float_val:{format_spec}}{unit}" if not np.isnan(float_val) else f"{default}{unit}"
    except (TypeError, ValueError):
        return f"{default}{unit}"

def calculate_jaccard(img1, img2):
    """Calculate Jaccard index."""
    intersection = np.sum(img1 & img2)
    union = np.sum(img1 | img2)
    return intersection / union if union > 0 else 0.0

def save_qc_values(output_results, root_RS, qc_values):
    """Save QC values to JSON file with proper serialization."""

    def numpy_to_python(obj):
        if isinstance(obj, (np.floating, np.integer)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (list, tuple)):
            return [numpy_to_python(x) for x in obj]
        elif isinstance(obj, dict):
            return {k: numpy_to_python(v) for k, v in obj.items()}
        return obj

    try:
        # Convert all numpy types to native Python types
        py_qc_values = numpy_to_python(qc_values)

        # Save to JSON file
        with open(opj(output_results, root_RS + '_QC_values.json'), 'w') as f:
            json.dump(py_qc_values, f, indent=2)

    except Exception as e:
        print(f"ERROR saving QC values: {str(e)}")
        # Try saving with simplified values
        try:
            simple_values = {}
            for k, v in qc_values.items():
                if isinstance(v, (np.floating, np.integer)):
                    simple_values[k] = float(v)
                elif isinstance(v, (list, tuple, np.ndarray)):
                    simple_values[k] = [float(x) if isinstance(x, (np.floating, np.integer)) else str(x)
                                        for x in v]
                else:
                    simple_values[k] = str(v)

            with open(opj(output_results, root_RS + '_QC_values_simple.json'), 'w') as f:
                json.dump(simple_values, f, indent=2)
        except Exception as e2:
            print(f"ERROR saving simplified QC values: {str(e2)}")


def fMRI_QC(correction_direction, dir_prepro_orig, ID, dir_prepro_template_process, dir_prepro_template_labels,
            dir_prepro_template_postprocessed, species, dir_prepro_orig_labels, template_dir, MAIN_PATH, RS, nb_run, sing_afni, diary_file):

    """
    Main function for fMRI quality control analysis with enhanced coverage analysis and image type detection.
    """

    nl = '##  Working on step ' + str(12) + '(function: _12_fMRI_QC).  ##'
    run_cmd.msg(nl, diary_file, 'HEADER')

    output_results = opj(dir_prepro_orig, 'fMRI_QC')

    if not ope(output_results):
        os.mkdir(output_results)

    # Try to remove directory if it exists
    if ope(output_results):
        shutil.rmtree(output_results, ignore_errors=True)
        # Wait a bit to avoid race condition (especially on network filesystems)
        time.sleep(0.1)

    # Ensure it's gone before recreating
    if not ope(output_results):
        os.makedirs(output_results)

    for i in range(0, int(nb_run)):
        qc_values = {}  # Dictionary to store all QC metrics
        root_RS = extract_filename(RS[i])

        ######### QC in func space #######
        func_filename = opj(dir_prepro_orig, root_RS + '_xdtrf_2ref.nii.gz')
        func_filename = opj(dir_prepro_orig_postprocessed, root_RS + '_space-acpc-func_desc-fMRI_residual.nii.gz')

        selected_atlases_matrix = ['EDNIxCSCLR', 1]
        atlas = [selected_atlases_matrix[0][0],selected_atlases_matrix[1][0]]
        atlas_filename = opj(dir_prepro_orig_labels, ID + '_seg-' + atlas[0] + '_dseg.nii.gz')

        if not ope(func_filename):
            nl = 'WARNING: ' + str(func_filename) + ' not found!!'
            run_cmd.msg(nl, diary_file, 'WARNING')
            continue

        if not ope(atlas_filename):
            nl = 'WARNING: no atlas lvl 1 LR found, this is a requirement for some of QC analysis'
            run_cmd.msg(nl, diary_file, 'WARNING')
            continue

        # Load and check image headers
        try:
            img1 = nib.load(atlas_filename)
            img2 = nib.load(func_filename)
            header1 = img1.header
            header2 = img2.header

            # Compare specific fields
            fields_to_compare = ['dim', 'pixdim']
            tolerance = 0.0000006
            differences = {}
            for field in fields_to_compare:
                value1 = header1[field][1:4]
                value2 = header2[field][1:4]
                if not np.allclose(value1, value2, atol=tolerance):
                    differences[field] = (value1, value2)

            if differences:
                for field, values in differences.items():
                    nl = f"{field}: Image 1: {values[0]}, Image 2: {values[1]}"
                    run_cmd.msg(nl, diary_file, 'OKGREEN')

                # Resample atlas to match fMRI
                dummy = nilearn.image.resample_to_img(atlas_filename, func_filename, interpolation='nearest')
                dummy.to_filename(atlas_filename)
                extracted_data = nib.load(atlas_filename).get_fdata()
                extracted_data = np.rint(extracted_data).astype(np.int32)
                labeled_img2 = nilearn.image.new_img_like(func_filename, extracted_data, copy_header=True)
                labeled_img2.to_filename(atlas_filename)
            else:
                nl = "No differences found in the specified fields."
                run_cmd.msg(nl, diary_file, 'OKGREEN')
        except Exception as e:
            nl = f"ERROR loading/checking images: {str(e)}"
            run_cmd.msg(nl, diary_file, 'FAIL')
            continue

        # Atlas labels
        labels = {
            2: '3rd ventricles', 3: '4th ventricles', 7: 'Cerebellum',
            5: 'Cerebellum White', 6: 'Cortical White matter', 1: 'CSF',
            4: 'Lateral ventricles', 8: 'Isocortex', 9: 'Allocortex',
            10: 'Periallocortex', 11: 'Subcortical areas',
            12: 'Diencephalon', 13: 'Brain stem'
        }

        # Load data
        try:
            atlas_img = nib.load(atlas_filename)
            atlas_data = atlas_img.get_fdata()
            atlas_data = np.rint(atlas_data).astype(np.int32)
            fmri_img = nib.load(func_filename)
            fmri_data = fmri_img.get_fdata()
        except Exception as e:
            nl = f"ERROR loading atlas/fMRI data: {str(e)}"
            run_cmd.msg(nl, diary_file, 'FAIL')
            continue

        # Compute SNR for each region
        def compute_snr(atlas_data, fmri_data, labels, hemisphere_offset=1000):
            time_points = fmri_data.shape[-1]
            snr_values = {}

            for label, region_name in labels.items():
                snr_values[region_name] = {'Left': [], 'Right': []}

                # Left hemisphere
                left_mask = atlas_data == label
                left_signal = fmri_data[left_mask]
                if left_signal.ndim == 1:
                    left_signal = left_signal[:, np.newaxis]

                # Right hemisphere
                right_mask = atlas_data == label + hemisphere_offset
                right_signal = fmri_data[right_mask]
                if right_signal.ndim == 1:
                    right_signal = right_signal[:, np.newaxis]

                for t in range(time_points):
                    # Skip if no voxels in region
                    if left_signal.size > 0:
                        left_signal_t = np.mean(left_signal[:, t]) if left_signal[:, t].size > 0 else np.nan
                    else:
                        left_signal_t = np.nan

                    if right_signal.size > 0:
                        right_signal_t = np.mean(right_signal[:, t]) if right_signal[:, t].size > 0 else np.nan
                    else:
                        right_signal_t = np.nan

                    # Compute noise from bottom 10% of the histogram
                    brain_values = fmri_data[:, :, :, t].flatten()
                    non_zero_values = brain_values[brain_values > 0]
                    if len(non_zero_values) > 0:
                        noise_threshold = np.percentile(non_zero_values, 10)
                        noise_values = non_zero_values[non_zero_values <= noise_threshold]
                        noise = np.std(noise_values) if len(noise_values) > 0 else np.nan
                    else:
                        noise = np.nan

                    # Calculate SNR if valid values
                    left_snr = left_signal_t / noise if not np.isnan(left_signal_t) and noise > 0 else np.nan
                    right_snr = right_signal_t / noise if not np.isnan(right_signal_t) and noise > 0 else np.nan

                    snr_values[region_name]['Left'].append(left_snr)
                    snr_values[region_name]['Right'].append(right_snr)

            return snr_values, noise

        try:
            snr_results, noise = compute_snr(atlas_data, fmri_data, labels)
        except Exception as e:
            nl = f"ERROR computing SNR: {str(e)}"
            run_cmd.msg(nl, diary_file, 'FAIL')
            snr_results = {}
            noise = np.nan

        # Store SNR results for plotting
        plot_snr_results = {
            'Gray_Matter': {k: v for k, v in snr_results.items()
                            if k in ['Isocortex', 'Allocortex', 'Periallocortex',
                                     'Subcortical areas', 'Diencephalon']},
            'White_Matter': {k: v for k, v in snr_results.items()
                             if k in ['Cortical White matter']}
        }

        # Compute average SNR for gray and white matter
        def compute_average_snr(snr_results, regions):
            avg_snr = []
            if not snr_results:
                return avg_snr

            time_points = len(next(iter(snr_results.values()))['Left'])

            for t in range(time_points):
                region_snr = []
                for region in regions:
                    if region in snr_results:
                        left = snr_results[region]['Left'][t] if not np.isnan(snr_results[region]['Left'][t]) else 0
                        right = snr_results[region]['Right'][t] if not np.isnan(snr_results[region]['Right'][t]) else 0
                        if not np.isnan(left) and not np.isnan(right):
                            region_snr.append(np.mean([left, right]))
                if region_snr:
                    avg_snr.append(np.mean(region_snr))
            return avg_snr

        gray_regions = ['Isocortex', 'Allocortex', 'Periallocortex',
                        'Subcortical areas', 'Diencephalon']
        white_regions = ['Cortical White matter']

        avg_gray_snr = compute_average_snr(snr_results, gray_regions)
        avg_white_snr = compute_average_snr(snr_results, white_regions)

        plot_snr_results['Avg_Gray_SNR'] = avg_gray_snr if avg_gray_snr else [np.nan]
        plot_snr_results['Avg_White_SNR'] = avg_white_snr if avg_white_snr else [np.nan]

        # Calculate average SNR values for table
        avg_snr_gray = np.nanmean(avg_gray_snr) if avg_gray_snr else np.nan
        avg_snr_white = np.nanmean(avg_white_snr) if avg_white_snr else np.nan
        qc_values['avg_snr_gray'] = avg_snr_gray
        qc_values['avg_snr_white'] = avg_snr_white

        # Calculate CNR and cortical contrast if we have valid values
        if not np.isnan(avg_snr_gray) and not np.isnan(avg_snr_white) and not np.isnan(noise):
            cnr_val = np.abs(avg_snr_gray - avg_snr_white) / noise
            cortical_contrast = ((avg_snr_white - avg_snr_gray) /
                                 ((avg_snr_white + avg_snr_gray) / 2))
        else:
            cnr_val = np.nan
            cortical_contrast = np.nan

        qc_values['cnr'] = cnr_val
        qc_values['cortical_contrast'] = cortical_contrast

        for imageQC, QCexplain in zip(
                [opj(dir_prepro_orig, root_RS + '_xdtrfwS_stdev.nii.gz'),
                 opj(dir_prepro_orig, root_RS + '_xdtrfwS_tsnr1.nii.gz'),
                 opj(dir_prepro_orig, root_RS + '_xdtrfwS_tsnr2.nii.gz')],
                ['stdev', 'TSNRcvarinv', 'TSNR']):

            try:
                nifti_img = image.load_img(imageQC)
                masker = NiftiLabelsMasker(
                    labels_img=atlas_filename, detrend=False, smoothing_fwhm=None,
                    standardize=False, low_pass=None, high_pass=None, t_r=None,
                    memory=None, verbose=5)

                time_series = masker.fit_transform(nifti_img)
                mean_signal_per_region = np.mean(time_series, axis=1)

                if QCexplain == 'stdev':
                    stdev = mean_signal_per_region
                    qc_values['stdev'] = stdev
                elif QCexplain == 'TSNRcvarinv':
                    TSNRcvarinv = mean_signal_per_region
                    qc_values['TSNRcvarinv'] = TSNRcvarinv
                elif QCexplain == 'TSNR':
                    TSNR = mean_signal_per_region
                    qc_values['TSNR'] = TSNR

            except Exception as e:
                nl = f"ERROR processing {QCexplain}: {str(e)}"
                run_cmd.msg(nl, diary_file, 'WARNING')

        # Motion analysis
        try:
            fd_jenkinson_array = fd_jenkinson(
                opj(dir_prepro_orig, root_RS + '.aff12.1D'), rmax=80.,
                out_file=opj(dir_prepro_orig, root_RS + '.aff12_fdfile.1D'),
                out_array=True)

            qc_values['mean_fd'] = np.nanmean(fd_jenkinson_array) if len(fd_jenkinson_array) > 0 else np.nan
        except Exception as e:
            nl = f"ERROR computing FD: {str(e)}"
            run_cmd.msg(nl, diary_file, 'WARNING')
            fd_jenkinson_array = []
            qc_values['mean_fd'] = np.nan

        # DVARS calculation
        try:
            calc_dvars_array = calc_dvars(
                opj(dir_prepro_orig, root_RS + '_xdtr_deob.nii.gz'),
                opj(dir_prepro_orig, root_RS + '_mask_final_in_fMRI_orig.nii.gz'))

            qc_values['mean_dvars'] = np.nanmean(calc_dvars_array) if len(calc_dvars_array) > 0 else np.nan
        except Exception as e:
            nl = f"ERROR computing DVARS: {str(e)}"
            run_cmd.msg(nl, diary_file, 'WARNING')
            calc_dvars_array = []
            qc_values['mean_dvars'] = np.nan

        # Motion parameters
        try:
            motion_params = np.genfromtxt(opj(dir_prepro_orig, root_RS + '_dfile.1D'))
            if motion_params.size == 0:
                raise ValueError("Empty motion parameters file")
        except Exception as e:
            nl = f"ERROR loading motion parameters: {str(e)}"
            run_cmd.msg(nl, diary_file, 'WARNING')
            motion_params = np.array([])

        # Global correlation
        try:
            gcor_val = global_correlation(
                opj(dir_prepro_orig, root_RS + '_xdtr_deob.nii.gz'),
                opj(dir_prepro_orig, root_RS + '_mask_final_in_fMRI_orig.nii.gz'))

            qc_values['gcor'] = gcor_val
        except Exception as e:
            nl = f"ERROR computing global correlation: {str(e)}"
            run_cmd.msg(nl, diary_file, 'WARNING')
            qc_values['gcor'] = np.nan

        # Ghost ratio
        if correction_direction in ['y', 'y-', 'z', 'z-', 'x', 'x-']:
            try:
                direct_aqc = 'y' if correction_direction in ['y', 'y-'] else (
                    'z' if correction_direction in ['z', 'z-'] else 'x')

                ghost_ratio = ghost_direction(
                    opj(dir_prepro_orig, root_RS + '_xdtr_deob.nii.gz'),
                    opj(dir_prepro_orig, root_RS + '_mask_final_in_fMRI_orig.nii.gz'),
                    direction=direct_aqc,
                    ref_file=opj(dir_prepro_orig, root_RS + '_xdtr_deob.nii.gz'),
                    out_file=opj(output_results, root_RS + '_ghost_mask.nii.gz'))

                qc_values['ghost_ratio'] = float(ghost_ratio)
            except Exception as e:
                nl = f"ERROR computing ghost ratio: {str(e)}"
                run_cmd.msg(nl, diary_file, 'WARNING')
                qc_values['ghost_ratio'] = np.nan
        else:
            qc_values['ghost_ratio'] = np.nan

        # FWHM estimation
        try:
            anat_file = opj(dir_prepro_orig, root_RS + '_residual.nii.gz')
            mask_file = opj(dir_prepro_orig, 'maskDilat.nii.gz')

            if ope(anat_file) and ope(mask_file):
                original_dir = os.getcwd()
                os.chdir(dir_prepro_orig)

                command = (sing_afni + '3dFWHMx -overwrite -combined -mask ' + mask_file +
                           ' -input ' + anat_file + ' -acf ' + anat_file[:-7] + '.acf.txt > ' +
                           anat_file[:-7] + '_FWHMx.txt')
                run_cmd.do(command, diary_file)
                os.chdir(original_dir)

                if ope(anat_file[:-7] + '.acf.txt'):
                    df = pd.read_csv(anat_file[:-7] + '.acf.txt', sep='\s+',
                                     names=["a", "b", "c", "combined_estimation"])
                    fwhm_val = np.float64(df["combined_estimation"][1])

                    if not np.isnan(fwhm_val):
                        qc_values['fwhm'] = fwhm_val
                    else:
                        qc_values['fwhm'] = np.nan
                else:
                    qc_values['fwhm'] = np.nan
            else:
                qc_values['fwhm'] = np.nan
        except Exception as e:
            nl = f"ERROR computing FWHM: {str(e)}"
            run_cmd.msg(nl, diary_file,'WARNING')
            qc_values['fwhm'] = np.nan

        # Image similarity metrics
        image_metrics = {
            'vals1': [],
            'vals2': [],
            'mi': np.nan,
            'cc': np.nan,
            'jaccard': np.nan}

        ############################################ QC in template space ################
        # Load images for coverage analysis and surrogate creation
        anat_img_path = opj(dir_prepro_template_process, 'BASE_SS_fMRI.nii.gz')
        fmri_img_path = opj(dir_prepro_template_postprocessed, 'all_runs_space-template-func_desc-fMRI_Mean_Image_SS.nii.gz')
        atlasfile = ID + '_seg-EDNIxCSC_dseg.nii.gz'
        atlas_filename = opj(dir_prepro_template_labels, atlasfile)
        atlas_img_4d = ants.image_read(atlas_filename)
        atlas_filename_template = ants.slice_image(atlas_img_4d, axis=3, idx=3)  # idx=3 → 4ème volume

        if ope(anat_img_path) and ope(fmri_img_path):
            try:
                # Load images with ANTs
                anat_img = ants.image_read(anat_img_path)
                fmri_img = ants.image_read(fmri_img_path)
                atlas_img = ants.image_read(atlas_filename_template)

                # Compute coverage statistics
                coverage_stats = compute_coverage_stats(anat_img, fmri_img)
                qc_values.update({
                    'pct_anat_uncovered': float(coverage_stats['pct_anat_uncovered']),
                    'pct_fmri_uncovered': float(coverage_stats['pct_fmri_uncovered'])
                })

                # Detect image type
                anat_type_stats = detect_image_type(anat_img, atlas_img,diary_file)
                image_type_stats = detect_image_type(fmri_img, atlas_img,diary_file)

                if image_type_stats['image_type'] == anat_type_stats['image_type']:
                    modality = 'Same (' + str(image_type_stats['image_type']) + '/' + str(anat_type_stats['image_type']) + ')'
                    run_cmd.msg('anat and func same modality', diary_file, 'OKGREEN')
                else:
                    modality = 'Different ! (' + str(image_type_stats['image_type']) + '/' + str(
                        anat_type_stats['image_type']) + ')'
                    run_cmd.msg('create a surogate anat', diary_file, 'WARNING')
                    # Create surrogate anatomical image matching fMRI contrast
                    surrogate_img = create_surrogate_from_label_atlas(anat_img,diary_file)

                    if surrogate_img is not None:
                        # Save surrogate image
                        surrogate_path = opj(output_results, root_RS + '_surrogate_anat.nii.gz')
                        ants.image_write(surrogate_img, surrogate_path)

                        # Update anatomical image path to use surrogate for correlation analysis
                        anat_img_path = surrogate_path
                        anat_img = surrogate_img

                        nl = "Created surrogate anatomical image matching fMRI contrast"
                        run_cmd.msg(nl, diary_file, 'OKBLUE')

                qc_values.update({
                    'Comp modal': modality,
                    'anat_image_type': image_type_stats['image_type'],
                    'wm_mean_intensity': image_type_stats['wm_mean'],
                    'gm_mean_intensity': image_type_stats['gm_mean'],
                    'csf_mean_intensity': image_type_stats['csf_mean'],
                    'wm_gm_contrast': image_type_stats['wm_gm_ratio'],
                    'wm_csf_contrast': image_type_stats['wm_csf_ratio']
                })

                # Print detection results
                nl = f"Anatomical image type detected as: {image_type_stats['image_type']}"
                run_cmd.msg(nl, diary_file, 'OKBLUE')

                nl = f"Coverage: {coverage_stats['pct_anat_uncovered']:.1f}% of anatomical not covered by fMRI"
                run_cmd.msg(nl, diary_file, 'OKBLUE')

            except Exception as e:
                nl = f"ERROR in coverage analysis: {str(e)}"
                run_cmd.msg(nl, diary_file, 'WARNING')
        else:
            nl = "WARNING: Missing anatomical or fMRI image for coverage analysis"
            run_cmd.msg(nl, diary_file, 'WARNING')

        if ope(anat_img_path) and ope(fmri_img_path):
            try:
                # Load images for metrics
                img1_data = ants.image_read(anat_img_path).numpy()
                img2_data = ants.image_read(fmri_img_path).numpy()
                mask = (img1_data != 0) & (img2_data != 0) & np.isfinite(img1_data) & np.isfinite(img2_data)
                vals1 = img1_data[mask]
                vals2 = img2_data[mask]

                vals1 = StandardScaler().fit_transform(vals1.reshape(-1, 1)).flatten()
                vals2 = StandardScaler().fit_transform(vals2.reshape(-1, 1)).flatten()
                if len(vals1) > 0 and len(vals2) > 0:
                    image_metrics['vals1'] = vals1
                    image_metrics['vals2'] = vals2

                    # Calculate metrics
                    hist_2d = np.histogram2d(vals1, vals2, bins=32)[0]
                    mi = mutual_info_score(None, None, contingency=hist_2d)
                    mi /= np.mean([entropy(np.sum(hist_2d, axis=1)),
                                   entropy(np.sum(hist_2d, axis=0))])
                    cc = pearsonr(vals1, vals2)[0]
                    jaccard = calculate_jaccard(img1_data != 0, img2_data != 0)

                    image_metrics['mi'] = mi
                    image_metrics['cc'] = cc
                    image_metrics['jaccard'] = jaccard

                    qc_values['mi'] = mi
                    qc_values['cc'] = cc
                    qc_values['jaccard'] = jaccard

                    # NMI index
                    fixed = ants.image_read(anat_img_path)
                    moving = ants.image_read(fmri_img_path)
                    nmi = ants.image_mutual_information(fixed, moving)
                    qc_values['nmi'] = nmi
                else:
                    qc_values['mi'] = np.nan
                    qc_values['cc'] = np.nan
                    qc_values['jaccard'] = np.nan
                    qc_values['nmi'] = np.nan
            except Exception as e:
                nl = f"ERROR computing image metrics: {str(e)}"
                run_cmd.msg(nl, diary_file, 'WARNING')
                qc_values['mi'] = np.nan
                qc_values['cc'] = np.nan
                qc_values['jaccard'] = np.nan
                qc_values['nmi'] = np.nan
        else:
            qc_values['mi'] = np.nan
            qc_values['cc'] = np.nan
            qc_values['jaccard'] = np.nan
            qc_values['nmi'] = np.nan

        # Motion metrics
        try:
            motion_enorm = np.loadtxt(opj(dir_prepro_orig, root_RS + 'motion_enorm.1D'))
            derivatives = np.loadtxt(opj(dir_prepro_orig, root_RS + '_xdtr_deriv.1D'))
            censor_1d = np.loadtxt(opj(dir_prepro_orig, root_RS + '_xdtr_censor.1D'))
            outcount = np.loadtxt(opj(dir_prepro_orig, root_RS + '_xdt_outcount.r.1D'),
                                  skiprows=2)

            qc_values['avg_enorm'] = np.nanmean(motion_enorm) if motion_enorm.size > 0 else np.nan
            qc_values['avg_velocity'] = np.nanmean(np.abs(derivatives), axis=0) if derivatives.size > 0 else np.nan
            qc_values['censor_fraction'] = 1 - np.nanmean(censor_1d) if censor_1d.size > 0 else np.nan
            qc_values['avg_outcount'] = np.nanmean(outcount) if outcount.size > 0 else np.nan
        except Exception as e:
            nl = f"ERROR loading motion metrics: {str(e)}"
            run_cmd.msg(nl, diary_file, 'WARNING')
            qc_values['avg_enorm'] = np.nan
            qc_values['avg_velocity'] = np.nan
            qc_values['censor_fraction'] = np.nan
            qc_values['avg_outcount'] = np.nan



        ######################################
        # Create the comprehensive QC figure #
        ######################################
        try:
            create_qc_figure(output_results, root_RS, plot_snr_results,
                             fd_jenkinson_array, calc_dvars_array, motion_params,
                             image_metrics, qc_values,
                             anat_img_path, fmri_img_path)
        except Exception as e:
            nl = f"ERROR creating QC figure: {str(e)}"
            run_cmd.msg(nl, diary_file, 'FAIL')

        # Save QC values to text file
        try:
            save_qc_values(output_results, root_RS, qc_values)
        except Exception as e:
            nl = f"ERROR saving QC values: {str(e)}"
            run_cmd.msg(nl, diary_file, 'FAIL')
