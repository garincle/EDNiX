import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from nilearn.decomposition import DictLearning
from scipy.optimize import linear_sum_assignment
from scipy.ndimage import gaussian_filter1d
from numpy.linalg import norm
import pandas as pd
import seaborn as sns


def spatial_corr(a, b):
    """Spatial correlation between two components"""
    a = a - a.mean()
    b = b - b.mean()
    return np.dot(a, b) / (norm(a) * norm(b) + 1e-8)


def match_components(components_A, components_B):
    """
    Match components between two splits using Hungarian algorithm.

    Returns mean correlation of best matches.
    """
    nA, nB = components_A.shape[0], components_B.shape[0]
    corr_mat = np.zeros((nA, nB))

    for i in range(nA):
        for j in range(nB):
            corr_mat[i, j] = spatial_corr(
                components_A[i], components_B[j]
            )

    row_ind, col_ind = linear_sum_assignment(-corr_mat)
    matched_corrs = corr_mat[row_ind, col_ind]

    return {
        'mean': matched_corrs.mean(),
        'std': matched_corrs.std(),
        'min': matched_corrs.min(),
        'max': matched_corrs.max(),
        'matched_corrs': matched_corrs
    }


def stability_analysis(
        images_dir,
        mask_img,
        component_list,
        alpha_dic,
        smoothing,
        TR,
        n_repeats=10,
        test_size=0.5,
        random_state=0,
        verbose=True
):
    """
    Split-half stability analysis for Dictionary Learning.

    Parameters:
    -----------
    images_dir : list
        List of functional image paths
    mask_img : str or Nifti1Image
        Brain mask
    component_list : list
        List of component numbers to test
    alpha_dic : float
        Dictionary Learning alpha parameter
    smoothing : float
        Spatial smoothing FWHM in mm
    TR : float
        Repetition time
    n_repeats : int
        Number of split-half repeats (default: 10)
    test_size : float
        Proportion for first split (default: 0.5)
    random_state : int
        Random seed
    verbose : bool
        Print progress

    Returns:
    --------
    dict : Stability scores for each component number
    """

    rng = np.random.RandomState(random_state)
    stability_scores = {}

    for n_comp in component_list:
        if verbose:
            print(f"\n[STABILITY] n_components = {n_comp}")
            print(f"  Running {n_repeats} split-half iterations...")

        scores = []
        component_level_scores = []  # Track individual component matches

        for r in range(n_repeats):
            if verbose and (r + 1) % 3 == 0:
                print(f"    Iteration {r + 1}/{n_repeats}...", end=" ")

            # --- Split subjects ---
            idx = np.arange(len(images_dir))
            rng.shuffle(idx)
            split = int(len(idx) * test_size)

            imgs_A = [images_dir[i] for i in idx[:split]]
            imgs_B = [images_dir[i] for i in idx[split:]]

            # --- DL on split A ---
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
                n_jobs=-1  # Use all CPUs
            )
            dl_A.fit(imgs_A)

            # --- DL on split B ---
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
            )
            dl_B.fit(imgs_B)

            # --- Stability score ---
            match_result = match_components(
                dl_A.components_,
                dl_B.components_
            )

            scores.append(match_result['mean'])
            component_level_scores.append(match_result['matched_corrs'])

            if verbose and (r + 1) % 3 == 0:
                print(f"stability = {match_result['mean']:.3f}")

        # Aggregate results
        stability_scores[n_comp] = {
            "mean": np.mean(scores),
            "std": np.std(scores),
            "median": np.median(scores),
            "min": np.min(scores),
            "max": np.max(scores),
            "all_scores": scores,
            "component_scores": component_level_scores  # For detailed analysis
        }

        if verbose:
            print(
                f"\n  ✓ Overall Stability = "
                f"{stability_scores[n_comp]['mean']:.3f} ± "
                f"{stability_scores[n_comp]['std']:.3f}"
            )
            print(f"    Range: [{stability_scores[n_comp]['min']:.3f}, "
                  f"{stability_scores[n_comp]['max']:.3f}]")

    return stability_scores


def plot_stability_comprehensive(stability_scores, output_dir,
                                 smooth_sigma=1.0, alpha=None, smoothing=None):
    """
    Comprehensive stability visualization with elbow detection.

    Creates a 3-panel figure showing:
    - Stability curve with elbow
    - Distribution of stability scores
    - Component-level stability breakdown
    """

    n_comps = np.array(sorted(stability_scores.keys()))
    means = np.array([stability_scores[k]["mean"] for k in n_comps])
    stds = np.array([stability_scores[k]["std"] for k in n_comps])
    mins = np.array([stability_scores[k]["min"] for k in n_comps])
    maxs = np.array([stability_scores[k]["max"] for k in n_comps])

    # --- Smooth curve for derivatives ---
    means_smooth = gaussian_filter1d(means, sigma=smooth_sigma)

    # --- Curvature-based elbow ---
    first_deriv = np.gradient(means_smooth)
    second_deriv = np.gradient(first_deriv)

    # Avoid edges
    if len(second_deriv) > 4:
        elbow_idx = np.argmax(-second_deriv[2:-2]) + 2
    else:
        elbow_idx = np.argmax(means)

    elbow_n = n_comps[elbow_idx]

    # --- Find stability plateau (where improvement drops below threshold) ---
    improvements = np.diff(means_smooth)
    plateau_threshold = 0.01  # Less than 1% improvement
    plateau_idx = np.where(improvements < plateau_threshold)[0]
    if len(plateau_idx) > 0:
        plateau_n = n_comps[plateau_idx[0] + 1]
    else:
        plateau_n = n_comps[-1]

    # --- Create figure ---
    fig = plt.figure(figsize=(15, 5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.3)

    # Panel 1: Main stability curve with elbow
    ax1 = fig.add_subplot(gs[0, 0])

    # Plot with error bars
    ax1.fill_between(n_comps, mins, maxs, alpha=0.2, color='#0173B2',
                     label='Min-Max range')
    ax1.errorbar(n_comps, means, yerr=stds, fmt='o-', capsize=4,
                 color='#0173B2', linewidth=2, markersize=6,
                 label="Mean ± SD", zorder=3)

    # Smoothed curve
    ax1.plot(n_comps, means_smooth, linestyle='--', linewidth=2,
             color='#DE8F05', alpha=0.7, label="Smoothed")

    # Thresholds
    ax1.axhline(0.7, color='green', linestyle='--', linewidth=1.5,
                alpha=0.6, label='Excellent (>0.7)')
    ax1.axhline(0.6, color='orange', linestyle='--', linewidth=1.5,
                alpha=0.6, label='Good (>0.6)')

    # Elbow
    ax1.axvline(elbow_n, color='red', linestyle=':', linewidth=2.5,
                label=f'Elbow (n={elbow_n})', zorder=4)
    ax1.plot(elbow_n, means[elbow_idx], 'rD', markersize=12,
             markeredgewidth=2, markeredgecolor='white', zorder=5)

    # Plateau
    if plateau_n != n_comps[-1]:
        ax1.axvline(plateau_n, color='purple', linestyle=':', linewidth=2,
                    alpha=0.7, label=f'Plateau (n={plateau_n})')

    ax1.set_xlabel('Number of Components', fontweight='bold', fontsize=11)
    ax1.set_ylabel('Stability (Spatial Correlation)', fontweight='bold', fontsize=11)
    ax1.set_title('A. Stability Curve with Elbow Detection', fontweight='bold', fontsize=12)
    ax1.set_ylim([0, 1.05])
    ax1.legend(fontsize=8, frameon=True, edgecolor='black', loc='lower right')
    ax1.grid(alpha=0.3)

    # Panel 2: Distribution of scores
    ax2 = fig.add_subplot(gs[0, 1])

    positions = []
    data_for_box = []
    colors_box = []

    for i, n in enumerate(n_comps):
        positions.append(n)
        data_for_box.append(stability_scores[n]['all_scores'])

        # Color code by stability
        mean_val = stability_scores[n]['mean']
        if mean_val >= 0.7:
            colors_box.append('#029E73')  # Green
        elif mean_val >= 0.6:
            colors_box.append('#DE8F05')  # Orange
        else:
            colors_box.append('#D55E00')  # Red

    bp = ax2.boxplot(data_for_box, positions=positions, widths=0.6,
                     patch_artist=True, showfliers=True)

    for patch, color in zip(bp['boxes'], colors_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # Highlight elbow
    elbow_pos = list(n_comps).index(elbow_n)
    bp['boxes'][elbow_pos].set_edgecolor('red')
    bp['boxes'][elbow_pos].set_linewidth(3)

    ax2.axhline(0.7, color='gray', linestyle='--', alpha=0.5)
    ax2.axhline(0.6, color='gray', linestyle='--', alpha=0.5)

    ax2.set_xlabel('Number of Components', fontweight='bold', fontsize=11)
    ax2.set_ylabel('Stability Distribution', fontweight='bold', fontsize=11)
    ax2.set_title('B. Stability Across Repeats', fontweight='bold', fontsize=12)
    ax2.set_ylim([0, 1.05])
    ax2.grid(alpha=0.3, axis='y')

    # Panel 3: Component-level analysis (show for elbow n)
    ax3 = fig.add_subplot(gs[0, 2])

    # Get component-level scores for elbow_n
    comp_scores = stability_scores[elbow_n]['component_scores']
    comp_means = np.mean(comp_scores, axis=0)
    comp_stds = np.std(comp_scores, axis=0)

    component_ids = np.arange(1, len(comp_means) + 1)

    ax3.bar(component_ids, comp_means, yerr=comp_stds,
            capsize=4, alpha=0.7, color='#0173B2', edgecolor='black')
    ax3.axhline(0.7, color='green', linestyle='--', alpha=0.6)
    ax3.axhline(0.6, color='orange', linestyle='--', alpha=0.6)

    # Highlight unreliable components
    unreliable = comp_means < 0.5
    if unreliable.any():
        ax3.bar(component_ids[unreliable], comp_means[unreliable],
                color='red', alpha=0.8, edgecolor='black')

    ax3.set_xlabel('Component ID', fontweight='bold', fontsize=11)
    ax3.set_ylabel('Mean Stability', fontweight='bold', fontsize=11)
    ax3.set_title(f'C. Component Reliability (n={elbow_n})', fontweight='bold', fontsize=12)
    ax3.set_ylim([0, 1.05])
    ax3.grid(alpha=0.3, axis='y')

    # Add info box
    n_reliable = np.sum(comp_means >= 0.6)
    info_text = f"{n_reliable}/{len(comp_means)} components\nstability ≥ 0.6"
    ax3.text(0.98, 0.98, info_text, transform=ax3.transAxes,
             ha='right', va='top', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    # Overall title
    if alpha is not None and smoothing is not None:
        plt.suptitle(f'Dictionary Learning Stability Analysis (α={alpha}, smoothing={smoothing}mm)',
                     fontsize=13, fontweight='bold', y=0.98)
    else:
        plt.suptitle('Dictionary Learning Stability Analysis',
                     fontsize=13, fontweight='bold', y=0.98)

    # Save
    out_fig = os.path.join(output_dir, "DL_stability_comprehensive.png")
    plt.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n[INFO] Figure saved: {out_fig}")
    print(f"[INFO] Elbow-based suggestion: n = {elbow_n} (stability = {means[elbow_idx]:.3f})")
    print(f"[INFO] Plateau point: n = {plateau_n}")

    # Print quality assessment
    print(f"\n[QUALITY ASSESSMENT]")
    for n in n_comps:
        mean_stab = stability_scores[n]['mean']
        if mean_stab >= 0.7:
            quality = "EXCELLENT ✓✓"
        elif mean_stab >= 0.6:
            quality = "GOOD ✓"
        elif mean_stab >= 0.5:
            quality = "ACCEPTABLE"
        else:
            quality = "POOR ✗"

        print(f"  n={n:2d}: {mean_stab:.3f} ± {stability_scores[n]['std']:.3f} - {quality}")

    return elbow_n, plateau_n


def optimize_parameters_with_stability(
        bids_dir, images_dir, mask_func,
        component_range=[7, 12, 16, 20],
        smoothing_range=[0, 2, 2.5, 3, 3.5, 4],
        alpha_range=[9, 9.5, 10, 10.5, 11],
        TR=1.0,
        n_repeats=10,
        test_size=0.5,
        random_state=0,
        label_specific='notprovided'):
    """
    Comprehensive parameter optimization using stability analysis.

    Tests all combinations of smoothing × alpha using split-half stability.

    Returns:
    --------
    DataFrame : All results
    dict : Best parameters
    """

    output_dir = os.path.join(bids_dir, 'DicL', 'QC_specificity_' + label_specific)
    os.makedirs(output_dir, exist_ok=True)

    all_results = []
    total_runs = len(smoothing_range) * len(alpha_range)
    run_count = 0

    print("\n" + "=" * 70)
    print("PARAMETER OPTIMIZATION WITH STABILITY ANALYSIS")
    print(f"Testing {total_runs} parameter combinations")
    print(f"Components: {component_range}")
    print(f"Smoothing: {smoothing_range}")
    print(f"Alpha: {alpha_range}")
    print(f"Repeats per config: {n_repeats}")
    print("=" * 70)

    for smooth in smoothing_range:
        for alpha in alpha_range:
            run_count += 1
            print(f"\n{'=' * 70}")
            print(f"[{run_count}/{total_runs}] Testing smoothing={smooth}mm, alpha={alpha}")
            print(f"{'=' * 70}")

            try:
                # Run stability analysis
                stability_scores = stability_analysis(
                    images_dir=images_dir,
                    mask_img=mask_func,
                    component_list=component_range,
                    alpha_dic=alpha,
                    smoothing=smooth,
                    TR=TR,
                    n_repeats=n_repeats,
                    test_size=test_size,
                    random_state=random_state,
                    verbose=True
                )

                # Create plots
                elbow_n, plateau_n = plot_stability_comprehensive(
                    stability_scores,
                    output_dir,
                    smooth_sigma=1.0,
                    alpha=alpha,
                    smoothing=smooth
                )

                # Rename plot with parameters
                old_name = os.path.join(output_dir, "DL_stability_comprehensive.png")
                new_name = os.path.join(output_dir,
                                        f"stability_alpha{alpha}_smooth{smooth}.png")
                if os.path.exists(old_name):
                    os.rename(old_name, new_name)

                # Store results for each component
                for n in component_range:
                    all_results.append({
                        'smoothing': smooth,
                        'alpha': alpha,
                        'n_components': n,
                        'stability_mean': stability_scores[n]['mean'],
                        'stability_std': stability_scores[n]['std'],
                        'stability_min': stability_scores[n]['min'],
                        'stability_max': stability_scores[n]['max'],
                        'is_elbow': (n == elbow_n),
                        'is_plateau': (n == plateau_n)
                    })

                print(f"\n✓ Completed: α={alpha}, s={smooth}mm")
                print(f"  Elbow: n={elbow_n}")
                print(f"  Best stability: {max(stability_scores[n]['mean'] for n in component_range):.3f}")

            except Exception as e:
                print(f"\n✗ Error: {e}")
                import traceback
                traceback.print_exc()
                continue

    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(os.path.join(output_dir, 'stability_optimization_results.csv'),
                      index=False)

    # Find best parameters
    best_params = find_best_parameters_from_stability(results_df)

    # Create summary plot
    create_stability_optimization_plot(results_df, best_params, output_dir)

    print("\n" + "=" * 70)
    print("OPTIMIZATION COMPLETE!")
    print("=" * 70)
    print(f"\nBest parameters:")
    print(f"  Smoothing: {best_params['smoothing']} mm")
    print(f"  Alpha: {best_params['alpha']}")
    print(f"  Best n: {best_params['n_components']}")
    print(f"  Stability: {best_params['stability']:.3f} ± {best_params['stability_std']:.3f}")
    print("=" * 70)

    return results_df, best_params


def find_best_parameters_from_stability(results_df):
    """
    Find best parameter combination based on stability.

    Priority:
    1. Highest mean stability across all n
    2. Prefer elbow points
    3. Prefer lower alpha (less regularization) if tie
    """

    # Aggregate by (smoothing, alpha)
    grouped = results_df.groupby(['smoothing', 'alpha']).agg({
        'stability_mean': 'mean',
        'stability_std': 'mean',
        'stability_min': 'min'
    }).reset_index()

    grouped.columns = ['smoothing', 'alpha', 'mean_stability', 'mean_std', 'min_stability']

    # Find best mean stability
    best_idx = grouped['mean_stability'].idxmax()
    best_row = grouped.loc[best_idx]

    # Get the specific n for this combo
    subset = results_df[
        (results_df['smoothing'] == best_row['smoothing']) &
        (results_df['alpha'] == best_row['alpha'])
        ]

    # Prefer elbow if available
    elbows = subset[subset['is_elbow']]
    if len(elbows) > 0:
        best_n_row = elbows.iloc[0]
    else:
        # Otherwise take highest stability
        best_n_row = subset.loc[subset['stability_mean'].idxmax()]

    return {
        'smoothing': best_row['smoothing'],
        'alpha': best_row['alpha'],
        'n_components': int(best_n_row['n_components']),
        'stability': best_n_row['stability_mean'],
        'stability_std': best_n_row['stability_std'],
        'mean_across_n': best_row['mean_stability']
    }


def create_stability_optimization_plot(results_df, best_params, output_dir):
    """
    Create heatmap visualization of parameter optimization results.
    """

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # For each n, create a heatmap
    n_values = sorted(results_df['n_components'].unique())

    for idx, n in enumerate(n_values[:4]):  # Show up to 4
        if idx >= 4:
            break

        ax = axes[idx // 2, idx % 2]

        # Pivot for heatmap
        subset = results_df[results_df['n_components'] == n]
        pivot = subset.pivot(index='alpha', columns='smoothing', values='stability_mean')

        # Plot
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn',
                    vmin=0, vmax=1, cbar_kws={'label': 'Stability'},
                    ax=ax)

        # Mark best if this is the selected n
        if n == best_params['n_components']:
            best_row_idx = list(pivot.index).index(best_params['alpha'])
            best_col_idx = list(pivot.columns).index(best_params['smoothing'])
            ax.add_patch(plt.Rectangle((best_col_idx, best_row_idx), 1, 1,
                                       fill=False, edgecolor='blue', lw=4))
            ax.set_title(f'n={n} Components ★ SELECTED', fontweight='bold', fontsize=12)
        else:
            ax.set_title(f'n={n} Components', fontweight='bold', fontsize=12)

        ax.set_xlabel('Smoothing (mm)', fontweight='bold')
        ax.set_ylabel('Alpha', fontweight='bold')

    plt.suptitle('Stability Across All Parameter Combinations',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, 'stability_optimization_heatmaps.png'),
                dpi=300, bbox_inches='tight')
    plt.close()