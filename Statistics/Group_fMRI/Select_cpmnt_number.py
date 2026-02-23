import numpy as np
import pandas as pd
import os
from nilearn.decomposition import DictLearning
from scipy.optimize import linear_sum_assignment
from scipy.ndimage import gaussian_filter1d
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from Tools import Load_subject_with_BIDS, load_bids, Load_BIDS_data_for_analysis

def spatial_corr(a, b):
    """Spatial correlation between two components"""
    a = a - a.mean()
    b = b - b.mean()
    return np.dot(a, b) / (norm(a) * norm(b) + 1e-8)


def match_components(components_A, components_B):
    """Match components using Hungarian algorithm."""
    nA, nB = components_A.shape[0], components_B.shape[0]
    corr_mat = np.zeros((nA, nB))

    for i in range(nA):
        for j in range(nB):
            corr_mat[i, j] = spatial_corr(components_A[i], components_B[j])

    row_ind, col_ind = linear_sum_assignment(-corr_mat)
    matched_corrs = corr_mat[row_ind, col_ind]

    return {
        'mean_correlation': matched_corrs.mean(),
        'std_correlation': matched_corrs.std(),
        'mean_instability': 1 - matched_corrs.mean(),
        'all_correlations': matched_corrs
    }


def normalized_stability_analysis(
        images_dir_list,  # List of image lists for different specificity groups
        mask_img,
        component_list,
        alpha_dic,
        smoothing,
        TR,
        n_repeats=10,
        test_size=0.5,
        random_state=0,
        min_subjects_per_group=None,  # Minimum subjects to use for all groups
        use_instability=True
):
    """
    Compare stability across different specificity groups with normalization.

    Parameters:
    -----------
    images_dir_list : list of lists
        Each element is a list of image paths for a specificity group
    min_subjects_per_group : int or None
        If provided, subsample each group to this number of subjects
        If None, use the minimum number across all groups
    """

    rng = np.random.RandomState(random_state)

    # Determine minimum number of subjects across groups
    group_sizes = [len(imgs) for imgs in images_dir_list]
    if min_subjects_per_group is None:
        min_subjects_per_group = min(group_sizes)

    print(f"Original group sizes: {group_sizes}")
    print(f"Using {min_subjects_per_group} subjects per group for fair comparison")

    # Prepare results storage
    all_results = []

    for group_idx, images_dir in enumerate(images_dir_list):
        group_results = []

        # Subsample to have equal number of subjects per group
        if len(images_dir) > min_subjects_per_group:
            # Randomly select subjects for this group
            indices = np.arange(len(images_dir))
            rng.shuffle(indices)
            selected_indices = indices[:min_subjects_per_group]
            selected_images = [images_dir[i] for i in selected_indices]
        else:
            selected_images = images_dir

        for n_comp in component_list:
            print(f"Group {group_idx}, n_components = {n_comp}")
            correlations = []
            instabilities = []

            for r in range(n_repeats):
                # Split-half for this repeat
                idx = np.arange(len(selected_images))
                rng.shuffle(idx)
                split = int(len(idx) * test_size)

                imgs_A = [selected_images[i] for i in idx[:split]]
                imgs_B = [selected_images[i] for i in idx[split:]]

                # Fit Dictionary Learning
                dl_A = DictLearning(
                    n_components=n_comp,
                    alpha=alpha_dic,
                    mask=mask_img,
                    standardize="zscore_sample",
                    smoothing_fwhm=smoothing,
                    detrend=False,
                    t_r=TR,
                    n_epochs=1,
                    random_state=r,
                    n_jobs=-1
                ).fit(imgs_A)

                dl_B = DictLearning(
                    n_components=n_comp,
                    alpha=alpha_dic,
                    mask=mask_img,
                    standardize="zscore_sample",
                    smoothing_fwhm=smoothing,
                    detrend=False,
                    t_r=TR,
                    n_epochs=1,
                    random_state=r + 100,
                    n_jobs=-1
                ).fit(imgs_B)

                result = match_components(dl_A.components_, dl_B.components_)
                correlations.append(result["mean_correlation"])
                instabilities.append(result["mean_instability"])

            group_results.append({
                "n_components": n_comp,
                "stability_mean": np.mean(correlations),
                "stability_std": np.std(correlations),
                "instability_mean": np.mean(instabilities),
                "instability_std": np.std(instabilities),
                "group_size": len(selected_images)
            })

        all_results.append(pd.DataFrame(group_results))

    return all_results


def compare_specificity_groups(
        bids_dir,
        specificity_groups=['Specific', 'Unspecific', 'Spurious', 'No', 'all'],
        component_list=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
        alpha_dic=10,
        smoothing=3,
        TR=2.0,
        mask_func=None,
        n_repeats=10,
        output_dir='./results'
):
    """
    Main function to compare stability across specificity groups.
    """

    os.makedirs(output_dir, exist_ok=True)

    # Load data for each specificity group
    images_dir_list = []
    group_labels = []
    group_sizes = []

    for spec in specificity_groups:
        print(f"\nLoading data for specificity: {spec}")

        # Use your existing function to load data
        all_ID, all_Session, all_data_path, all_ID_max, all_Session_max, all_data_path_max, mean_imgs, images_dir = Load_BIDS_data_for_analysis.reverse_load_data_bids(
            bids_dir, spec, file_pattern="_space-template_desc-fMRI_residual.nii.gz"
        )

        if len(images_dir) > 0:
            images_dir_list.append(images_dir)
            group_labels.append(spec)
            group_sizes.append(len(images_dir))
            print(f"  Found {len(images_dir)} subjects")
        else:
            print(f"  WARNING: No subjects found for {spec}")

    print(f"\nGroup sizes: {dict(zip(group_labels, group_sizes))}")

    # Run normalized stability analysis
    all_results = normalized_stability_analysis(
        images_dir_list=images_dir_list,
        mask_img=mask_func,
        component_list=component_list,
        alpha_dic=alpha_dic,
        smoothing=smoothing,
        TR=TR,
        n_repeats=n_repeats,
        test_size=0.5,
        min_subjects_per_group=min(group_sizes)  # Normalize to smallest group
    )

    # Combine results into a single DataFrame for comparison
    comparison_data = []
    for label, df in zip(group_labels, all_results):
        df['specificity'] = label
        comparison_data.append(df)

    comparison_df = pd.concat(comparison_data, ignore_index=True)

    # Save results
    csv_path = os.path.join(output_dir, 'specificity_comparison.csv')
    comparison_df.to_csv(csv_path, index=False)
    print(f"\nResults saved to: {csv_path}")

    # Create comparison plots
    plot_comparison_results(comparison_df, output_dir, component_list)

    return comparison_df


def plot_comparison_results(comparison_df, output_dir, component_list):
    """Plot comparison of stability across specificity groups."""

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Instability vs Components for each group
    ax1 = axes[0]
    for spec in comparison_df['specificity'].unique():
        group_data = comparison_df[comparison_df['specificity'] == spec]
        ax1.errorbar(
            group_data['n_components'],
            group_data['instability_mean'],
            yerr=group_data['instability_std'],
            label=spec,
            marker='o',
            capsize=4
        )

    ax1.set_xlabel('Number of Components')
    ax1.set_ylabel('Instability (1 - correlation)')
    ax1.set_title('Instability Comparison Across Specificity Groups')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Stability vs Components for each group
    ax2 = axes[1]
    for spec in comparison_df['specificity'].unique():
        group_data = comparison_df[comparison_df['specificity'] == spec]
        ax2.errorbar(
            group_data['n_components'],
            group_data['stability_mean'],
            yerr=group_data['stability_std'],
            label=spec,
            marker='s',
            capsize=4
        )

    ax2.set_xlabel('Number of Components')
    ax2.set_ylabel('Stability (correlation)')
    ax2.set_title('Stability Comparison Across Specificity Groups')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Plot 3: Best stability per group
    ax3 = axes[2]
    best_results = []
    for spec in comparison_df['specificity'].unique():
        group_data = comparison_df[comparison_df['specificity'] == spec]
        best_idx = group_data['stability_mean'].idxmax()
        best_results.append({
            'specificity': spec,
            'n_components': group_data.loc[best_idx, 'n_components'],
            'stability': group_data.loc[best_idx, 'stability_mean'],
            'instability': group_data.loc[best_idx, 'instability_mean']
        })

    best_df = pd.DataFrame(best_results)
    x_pos = np.arange(len(best_df))
    bars = ax3.bar(x_pos, best_df['stability'], alpha=0.7)

    # Color bars by instability (red for high instability, green for low)
    colors = plt.cm.RdYlGn_r(best_df['instability'])  # Reverse colormap
    for bar, color in zip(bars, colors):
        bar.set_color(color)

    ax3.set_xlabel('Specificity Group')
    ax3.set_ylabel('Best Stability Achieved')
    ax3.set_title('Best Stability per Group (colored by instability)')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(best_df['specificity'], rotation=45)

    # Add text labels
    for i, row in best_df.iterrows():
        ax3.text(i, row['stability'] + 0.02,
                 f"n={int(row['n_components'])}",
                 ha='center', fontsize=9)

    # Add colorbar for instability
    sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlGn_r,
                               norm=plt.Normalize(vmin=0, vmax=0.5))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax3, orientation='vertical', pad=0.1)
    cbar.set_label('Instability (lower = better)')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'specificity_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Create a summary table
    print("\n" + "=" * 80)
    print("SPECIFICITY COMPARISON SUMMARY")
    print("=" * 80)

    for spec in comparison_df['specificity'].unique():
        group_data = comparison_df[comparison_df['specificity'] == spec]
        print(f"\n{spec}:")
        print("-" * 40)
        print(f"  Best stability: {group_data['stability_mean'].max():.3f} "
              f"at n={group_data.loc[group_data['stability_mean'].idxmax(), 'n_components']}")
        print(f"  Best instability: {group_data['instability_mean'].min():.3f} "
              f"at n={group_data.loc[group_data['instability_mean'].idxmin(), 'n_components']}")

        # Find elbow point
        n_comps = np.array(group_data['n_components'])
        values = np.array(group_data['instability_mean'])
        if len(values) > 3:
            values_smooth = gaussian_filter1d(values, sigma=1.0)
            second_deriv = np.gradient(np.gradient(values_smooth))
            # Avoid edges
            search_region = second_deriv[2:-2]
            elbow_idx_in_region = np.argmax(search_region) + 2
            elbow_n = n_comps[elbow_idx_in_region]
            print(f"  Elbow (Yeo method): n={elbow_n}")