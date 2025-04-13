import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.gridspec import GridSpec
import seaborn as sns
import json
import datetime
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
import shutil
from fonctions.extract_filename import extract_filename


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# Path utilities
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists


def validate_matrix(matrix, roi_names):
    """Validate matrix with detailed asymmetry reporting"""
    if matrix.shape[0] != matrix.shape[1]:
        if isinstance(roi_names, (list, np.ndarray)):
            if len(roi_names) == matrix.shape[0]:
                matrix = matrix[:, :len(roi_names)]
            elif len(roi_names) == matrix.shape[1]:
                matrix = matrix[:len(roi_names), :]
            else:
                raise ValueError(
                    f"ROI names count ({len(roi_names)}) doesn't match "
                    f"matrix dimensions ({matrix.shape[0]} rows, {matrix.shape[1]} columns)"
                )

    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Failed to make matrix square. Final shape: {matrix.shape}")

    if isinstance(roi_names, (list, np.ndarray)) and len(roi_names) != matrix.shape[0]:
        raise ValueError(
            f"ROI names count ({len(roi_names)}) doesn't match "
            f"matrix dimension ({matrix.shape[0]}) after adjustment"
        )

    return matrix, roi_names


def analyze_specificity(correlation_matrix, roi_names, specific_thresh, unspecific_thresh):
    try:
        corr_matrix, roi_names = validate_matrix(correlation_matrix, roi_names)
    except Exception as e:
        return {'error': str(e)}

    results = {
        'target_specificity': None,
        'roi_specificity': [],
        'matrix_shape': corr_matrix.shape,
        'n_rois': len(roi_names)
    }

    # Find target ROIs
    l_sensory_idx, r_sensory_idx, l_mpfc_idx, r_mpfc_idx = None, None, None, None
    sensory_patterns = ['Somatosensory', 'Motor', 'Visual', 'Auditory']
    mpfc_patterns = ['mPFC', 'Prefrontal', 'dlPFC', 'vlPFC']

    # Find sensory ROIs
    for pattern in sensory_patterns:
        l_matches = [i for i, name in enumerate(roi_names) if name.startswith('L_') and pattern in name]
        r_matches = [i for i, name in enumerate(roi_names) if name.startswith('R_') and pattern in name]
        if l_matches and r_matches:
            l_sensory_idx, r_sensory_idx = l_matches[0], r_matches[0]
            break

    # Find mPFC ROIs
    for pattern in mpfc_patterns:
        l_matches = [i for i, name in enumerate(roi_names) if name.startswith('L_') and pattern in name]
        r_matches = [i for i, name in enumerate(roi_names) if name.startswith('R_') and pattern in name]
        if l_matches:
            l_mpfc_idx = l_matches[0]
        if r_matches:
            r_mpfc_idx = r_matches[0]
        if l_mpfc_idx is not None or r_mpfc_idx is not None:
            break

    # Target specificity analysis
    if l_sensory_idx is not None and r_sensory_idx is not None:
        specific_corr = corr_matrix[l_sensory_idx, r_sensory_idx]
        non_specific_corrs = []

        connections = [
            (l_mpfc_idx, l_sensory_idx),
            (r_mpfc_idx, l_sensory_idx),
            (l_mpfc_idx, r_sensory_idx),
            (r_mpfc_idx, r_sensory_idx)]

        for idx1, idx2 in connections:
            if idx1 is not None and idx2 is not None:
                non_specific_corrs.append(corr_matrix[idx1, idx2])

        non_specific_corr = np.mean(non_specific_corrs) if non_specific_corrs else np.nan

        if (specific_corr >= specific_thresh) and (non_specific_corr < unspecific_thresh):
            category = 'Specific'
        elif (specific_corr >= specific_thresh) and (non_specific_corr >= unspecific_thresh):
            category = 'Unspecific'
        elif (abs(specific_corr) < specific_thresh) and (abs(non_specific_corr) < unspecific_thresh):
            category = 'No'
        else:
            category = 'Spurious'

        results['target_specificity'] = {
            'Specific_ROI_pair': f"{'R-' + str(roi_names[l_sensory_idx])}",
            'Specific_Correlation': round(specific_corr, 3),
            'NonSpecific_ROI_pair': f"{'R-' + str(roi_names[l_mpfc_idx])}",
            'NonSpecific_Correlation': round(non_specific_corr, 3),
            'Category': category
        }

    return results


def calculate_network_metrics(correlation_matrix):
    """Calculate various network metrics with safety checks"""
    metrics = {}
    n_nodes = correlation_matrix.shape[0]

    # Eigenvalue analysis
    try:
        eigenvalues = np.sort(np.abs(np.linalg.eig(correlation_matrix)[0]))[::-1]
        metrics['Top_Eigenvalue'] = round(eigenvalues[0], 3)
        metrics['Eigenvalue_Ratio'] = round(eigenvalues[0] / eigenvalues[1], 3) if len(eigenvalues) > 1 else np.nan
    except Exception as e:
        metrics.update({'Top_Eigenvalue': np.nan, 'Eigenvalue_Ratio': np.nan})

    # Clustering metrics
    if n_nodes > 1:
        n_clusters = min(7, n_nodes - 1)
        try:
            clustering = SpectralClustering(
                n_clusters=n_clusters,
                affinity='precomputed',
                assign_labels='kmeans',
                random_state=42
            )
            labels = clustering.fit_predict((correlation_matrix + 1) / 2)

            metrics.update({
                'Silhouette_Score': round(silhouette_score(correlation_matrix, labels), 3),
                'Davies_Bouldin': round(davies_bouldin_score(correlation_matrix, labels), 3),
                'Calinski_Harabasz': round(calinski_harabasz_score(correlation_matrix, labels), 3)
            })
        except Exception:
            metrics.update({
                'Silhouette_Score': np.nan,
                'Davies_Bouldin': np.nan,
                'Calinski_Harabasz': np.nan
            })

    # Correlation statistics
    mask = np.tri(n_nodes, k=-1).astype(bool)
    if mask.any():
        metrics.update({
            'Mean_Correlation': round(np.mean(correlation_matrix[mask]), 3),
            'Median_Correlation': round(np.median(correlation_matrix[mask]), 3),
            'Std_Correlation': round(np.std(correlation_matrix[mask]), 3)
        })
    else:
        metrics.update({
            'Mean_Correlation': np.nan,
            'Median_Correlation': np.nan,
            'Std_Correlation': np.nan
        })

    return metrics


def analyze_hemisphere_comparisons(correlation_matrix, roi_names):
    """Compare intra- vs inter-hemisphere connectivity with detailed stats"""
    left_rois = [i for i, name in enumerate(roi_names) if name.startswith('L_')]
    right_rois = [i for i, name in enumerate(roi_names) if name.startswith('R_')]

    if not left_rois or not right_rois:
        return None

    intra_left = correlation_matrix[np.ix_(left_rois, left_rois)]
    intra_right = correlation_matrix[np.ix_(right_rois, right_rois)]
    inter_lr_vals = []
    for name in roi_names:
        if name.startswith('L_'):
            match_name = 'R_' + name[2:]  # Replace 'L_' with 'R_'
            if match_name in roi_names:
                i = roi_names.index(name)
                j = roi_names.index(match_name)
                inter_lr_vals.append(correlation_matrix[i, j])

    # Get values (excluding diagonal for intra)
    intra_left_vals = intra_left[np.triu_indices_from(intra_left, k=1)]
    intra_right_vals = intra_right[np.triu_indices_from(intra_right, k=1)]
    # Statistical tests
    t_left_right, p_left_right = stats.ttest_rel(intra_left_vals, intra_right_vals)
    t_intra_inter, p_intra_inter = stats.ttest_ind(
        np.concatenate([intra_left_vals, intra_right_vals]),
        inter_lr_vals)

    return {
        'intra_left': intra_left_vals,
        'intra_right': intra_right_vals,
        'inter': inter_lr_vals,
        'p_intra_vs_inter': p_intra_inter,
        'p_left_vs_right': p_left_right,
        'intra_left_mean': np.mean(intra_left_vals),
        'intra_right_mean': np.mean(intra_right_vals),
        'inter_mean': np.mean(inter_lr_vals),
        't_intra_vs_inter': t_intra_inter,
        't_left_vs_right': t_left_right
    }


def generate_qc_plots(corr_matrix, roi_names, output_dir, prefix,
                      specific_thresh, unspecific_thresh):
    """Generate QC plots with matrix stats in table"""
    fig = plt.figure(figsize=(18, 10))
    gs = GridSpec(2, 3, figure=fig,
                  width_ratios=[1, 1, 1.5],
                  height_ratios=[1, 1.2],
                  wspace=0.3, hspace=0.4)

    # --- Top Row: Three Panels ---
    # Panel A: Correlation Matrix
    ax1 = fig.add_subplot(gs[0, 0])
    img = ax1.imshow(corr_matrix, cmap='cold_hot', vmin=-0.8, vmax=0.8)
    plt.colorbar(img, ax=ax1, fraction=0.04, pad=0.01)
    ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # Panel B: Correlation Distribution
    ax2 = fig.add_subplot(gs[0, 1])
    mask = np.tri(corr_matrix.shape[0], k=-1).astype(bool)
    corr_values = corr_matrix[mask]
    ax2.hist(corr_values, bins=50, color='#4e79a7', edgecolor='white', alpha=0.8)
    ax2.set_xlabel('Correlation Value')
    ax2.set_ylabel('Count')
    ax2.text(-0.1, 1.1, 'B', transform=ax2.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # Panel C: Hemisphere Connectivity
    ax3 = fig.add_subplot(gs[0, 2])
    hemi_results = analyze_hemisphere_comparisons(corr_matrix, roi_names) or {}

    # Initialize hemisphere data with default values
    left_data = hemi_results.get('intra_left', [np.nan])
    right_data = hemi_results.get('intra_right', [np.nan])
    inter_data = hemi_results.get('inter', [np.nan])

    # Create dataframe for plotting
    plot_data = pd.DataFrame({
        'Connectivity': np.concatenate([left_data, right_data, inter_data]),
        'Type': (['Intra-Left'] * len(left_data) +
                 ['Intra-Right'] * len(right_data) +
                 ['Inter-Hemi'] * len(inter_data))
    })

    # Violin plot with connections
    palette = ['#4e79a7', '#f28e2b', '#59a14f']
    sns.violinplot(x='Type', y='Connectivity', data=plot_data,
                   palette=palette, ax=ax3, inner=None, saturation=0.8)

    # Connect left-right pairs if we have both
    if len(left_data) == len(right_data):
        for l_val, r_val in zip(left_data, right_data):
            jitter = np.random.uniform(-0.1, 0.1)
            ax3.plot([0 + jitter, 1 + jitter], [l_val, r_val],
                     color='#7f7f7f', alpha=0.3, linewidth=1, linestyle='--')

    # Add points
    sns.stripplot(x='Type', y='Connectivity', data=plot_data,
                  palette=palette, alpha=0.7, size=5, jitter=0.1, ax=ax3)

    # Significance markers
    y_max = max(np.nanmax(left_data), np.nanmax(right_data), np.nanmax(inter_data))
    if 'p_left_vs_right' in hemi_results:
        ax3.text(0.5, y_max * 1.05, f"p = {hemi_results['p_left_vs_right']:.3f}",
                 ha='center', fontsize=10)
    if 'p_intra_vs_inter' in hemi_results:
        ax3.text(1.5, y_max * 1.1, f"p = {hemi_results['p_intra_vs_inter']:.3f}",
                 ha='center', fontsize=10)

    ax3.set_ylim(min(np.nanmin(left_data), np.nanmin(right_data), np.nanmin(inter_data)) * 1.1,
                 y_max * 1.15)
    ax3.text(-0.1, 1.1, 'C', transform=ax3.transAxes,
             fontsize=16, fontweight='bold', va='top')

    # --- Bottom Row: Consolidated Table ---
    ax4 = fig.add_subplot(gs[1, :])
    ax4.axis('off')

    # Get all results
    metrics = calculate_network_metrics(corr_matrix) or {}
    spec_results = analyze_specificity(corr_matrix, roi_names, specific_thresh, unspecific_thresh) or {}
    target_data = spec_results.get('target_specificity', {})

    # Calculate matrix statistics (using only lower triangle)
    matrix_stats = {
        'Min': np.nanmin(corr_values),
        'Max': np.nanmax(corr_values),
        'Mean': np.nanmean(corr_values),
        'Std': np.nanstd(corr_values)
    }

    # Additional metrics
    pos_frac = np.mean(corr_values > 0)
    neg_frac = np.mean(corr_values < 0)
    n_rois = corr_matrix.shape[0]

    # Color coding for specificity categories
    category_colors = {
        'Specific': '#2ca02c',
        'Unspecific': '#ff7f0e',
        'No': '#1f77b4',
        'Spurious': '#d62728'
    }

    # Prepare table content
    table_content = [
        ["Network Metrics", "", "Hemisphere Results", "", "Matrix Stats", ""],
        ["Top Eigenvalue", f"{metrics.get('Top_Eigenvalue', np.nan):.2f}",
         "Left vs Right p-value", f"{hemi_results.get('p_left_vs_right', np.nan):.3f}",
         "Min", f"{matrix_stats['Min']:.3f}"],
        ["Eigenvalue Ratio", f"{metrics.get('Eigenvalue_Ratio', np.nan):.3f}",
         "Intra vs Inter p-value", f"{hemi_results.get('p_intra_vs_inter', np.nan):.3f}",
         "Max", f"{matrix_stats['Max']:.3f}"],
        ["Silhouette Score", f"{metrics.get('Silhouette_Score', np.nan):.3f}",
         "Mean Intra-Left", f"{np.nanmean(left_data):.3f}",
         "Mean", f"{matrix_stats['Mean']:.3f}"],
        ["Davies-Bouldin", f"{metrics.get('Davies_Bouldin', np.nan):.3f}",
         "Mean Intra-Right", f"{np.nanmean(right_data):.3f}",
         "Std", f"{matrix_stats['Std']:.3f}"],
        ["Calinski-Harabasz", f"{metrics.get('Calinski_Harabasz', np.nan):.0f}",
         "Mean Inter-Hemi", f"{np.nanmean(inter_data):.3f}",
         "Median", f"{np.nanmedian(corr_values):.3f}"],
        ["Positive Correlations", f"{pos_frac:.1%}",
         "Negative Correlations", f"{neg_frac:.1%}",
         "IQR", f"{np.nanpercentile(corr_values, 75) - np.nanpercentile(corr_values, 25):.3f}"],
        ["Specificity Results", "", "", "", "", ""],
        ["Specific Pair", target_data.get('Specific_ROI_pair', 'N/A'),
         "Correlation", f"{target_data.get('Specific_Correlation', np.nan):.3f}",
         "Effect Size",
         f"{target_data.get('Specific_Correlation', 0) - target_data.get('NonSpecific_Correlation', 0):.3f}"],
        ["Non-Specific Pair", target_data.get('NonSpecific_ROI_pair', 'N/A'),
         "Correlation", f"{target_data.get('NonSpecific_Correlation', np.nan):.3f}",
         "Classification", target_data.get('Category', 'N/A')],
    ]

    # Create table
    table = ax4.table(
        cellText=table_content,
        loc='center',
        cellLoc='center',
        colWidths=[0.2, 0.15, 0.2, 0.15, 0.15, 0.15]
    )

    # Style table
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    # Apply colors and formatting
    for (row, col), cell in table.get_celld().items():
        # Section headers
        if row in [0, 7]:
            cell.set_facecolor('#404040')
            cell.set_text_props(color='white', weight='bold')
        # Color p-values
        elif row == 1 and col == 3:
            p_val = hemi_results.get('p_left_vs_right', 1)
            cell.set_facecolor('#d62728' if p_val < 0.05 else '#2ca02c')
            cell.set_text_props(color='white')
        elif row == 2 and col == 3:
            p_val = hemi_results.get('p_intra_vs_inter', 1)
            cell.set_facecolor('#d62728' if p_val > 0.05 else '#2ca02c')
            cell.set_text_props(color='white')
        # Classification cell
        elif row == 9 and col == 5:
            classification = target_data.get('Category', 'N/A')
            if classification in category_colors:
                cell.set_facecolor(category_colors[classification])
                cell.set_text_props(color='white' if classification != 'No' else 'black')
        # Metric names
        elif row > 0 and col in [0, 2, 4]:
            cell.set_facecolor('#f0f0f0')

    ax4.text(-0.05, 1.05, 'D', transform=ax4.transAxes,
             fontsize=16, fontweight='bold', va='top')

    plt.savefig(os.path.join(output_dir, f'{prefix}_qc_report.png'),
                dpi=300, bbox_inches='tight')
    plt.close()

    return metrics, hemi_results


def filter_bilateral_rois(matrix, roi_names):
    """Keeps only bilateral ROIs (L_ and R_ pairs) and reports removed ones."""
    bilateral_suffixes = set()
    left_rois = {name[2:]: i for i, name in enumerate(roi_names) if name.startswith("L_")}
    right_rois = {name[2:]: i for i, name in enumerate(roi_names) if name.startswith("R_")}

    # Identify suffixes that have both L_ and R_ counterparts
    for suffix in left_rois:
        if suffix in right_rois:
            bilateral_suffixes.add(suffix)

    # Keep only indices that are bilateral
    keep_indices = []
    for i, name in enumerate(roi_names):
        if name.startswith(("L_", "R_")):
            suffix = name[2:]
            if suffix in bilateral_suffixes:
                keep_indices.append(i)

    # Apply filtering
    filtered_matrix = matrix[np.ix_(keep_indices, keep_indices)]
    filtered_names = [roi_names[i] for i in keep_indices]

    return filtered_matrix, filtered_names


def load_and_validate_matrix(matrix_path):
    """Load and validate matrix from file with robust error handling and bilateral checking"""
    try:
        if not os.path.exists(matrix_path):
            raise FileNotFoundError(f"Matrix file not found: {matrix_path}")

        # Load with pandas
        df = pd.read_csv(matrix_path, index_col=0)

        # Clean headers and index
        df.columns = df.columns.astype(str).str.strip()
        df.index = df.index.astype(str).str.strip()

        # Identify bilateral ROIs (those with L_ and R_ versions)
        all_rois = set(df.columns)
        bilateral_pairs = []

        for roi in all_rois:
            if roi.startswith('L_'):
                counterpart = 'R_' + roi[2:]
                if counterpart in all_rois:
                    bilateral_pairs.extend([roi, counterpart])
            elif roi.startswith('R_'):
                counterpart = 'L_' + roi[2:]
                if counterpart in all_rois:
                    bilateral_pairs.extend([roi, counterpart])

        # Remove duplicates while preserving order
        bilateral_rois = []
        seen = set()
        for roi in bilateral_pairs:
            if roi not in seen:
                seen.add(roi)
                bilateral_rois.append(roi)

        # If we found bilateral pairs, filter the matrix
        if bilateral_rois:
            df = df.loc[bilateral_rois, bilateral_rois]

        # Final validation
        if not list(df.columns) == list(df.index):
            raise ValueError("Header and index mismatch after bilateral filtering")

        matrix = df.values.astype(float)
        roi_names = list(df.columns)

        return matrix, roi_names

    except Exception as e:
        raise ValueError(f"Failed to load and validate matrix: {str(e)}")


def fMRI_QC_matrix(dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3,
                   specific_roi_tresh, unspecific_ROI_thresh, RS, nb_run, diary_file):
    """Optimized QC analysis with updated output"""
    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '## Working on fMRI QC matrix analysis (optimized) ##'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    for direction in [dir_fMRI_Refth_RS_prepro1, dir_fMRI_Refth_RS_prepro2, dir_fMRI_Refth_RS_prepro3]:
        out_results = opj(direction, '10_Results')
        os.makedirs(out_results, exist_ok=True)

        out_results_V = opj(out_results, 'fMRI_QC_matrix')
        if ope(out_results_V):
            shutil.rmtree(out_results_V)
        os.makedirs(out_results_V)

        for i in range(int(nb_run)):
            root_RS = extract_filename(RS[i])
            matrix_file = opj(direction, '10_Results', 'correl_matrix',
                              f'atlaslvl3_LR_run_{i}_matrix.csv')

            if not ope(matrix_file):
                nl = f'WARNING: Missing matrix file {matrix_file}'
                print(bcolors.WARNING + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                continue

            try:
                # Load and validate matrix with new robust function
                full_corr, roi_names = load_and_validate_matrix(matrix_file)

                # Filter to keep only bilateral ROIs
                full_corr, roi_names = filter_bilateral_rois(full_corr, roi_names)

                # Generate QC plots and get metrics
                metrics, hemi_results = generate_qc_plots(
                    full_corr, roi_names, out_results_V, root_RS,
                    specific_roi_tresh, unspecific_ROI_thresh)

                # Save results
                results = {
                    'network_metrics': metrics,
                    'hemisphere_results': hemi_results,
                    'specificity_results': analyze_specificity(
                        full_corr, roi_names, specific_roi_tresh, unspecific_ROI_thresh),
                    'timestamp': str(datetime.datetime.now())
                }

                # Helper to convert ndarrays to lists recursively
                def convert_ndarrays(obj):
                    if isinstance(obj, dict):
                        return {k: convert_ndarrays(v) for k, v in obj.items()}
                    elif isinstance(obj, list):
                        return [convert_ndarrays(v) for v in obj]
                    elif isinstance(obj, np.ndarray):
                        return obj.tolist()
                    else:
                        return obj

                # Save to JSON
                with open(opj(out_results_V, f'{root_RS}_full_results.json'), 'w') as f:
                    json.dump(convert_ndarrays(results), f, indent=2)

                # Save CSV versions
                pd.DataFrame(full_corr, index=roi_names, columns=roi_names).to_csv(
                    opj(out_results_V, f'{root_RS}_corr_matrix.csv'))

                diary.write(f'\n{root_RS} analysis complete\n')

            except Exception as e:
                nl = f'Error processing {matrix_file}: {str(e)}'
                print(bcolors.FAIL + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
                continue

    diary.write('\nProcessing complete.\n')
    diary.close()
    print(bcolors.OKGREEN + 'QC analysis completed successfully' + bcolors.ENDC)