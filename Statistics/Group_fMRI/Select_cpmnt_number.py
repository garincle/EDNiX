from nilearn.decomposition import DictLearning
import os
import numpy as np

opn = os.path.normpath
opj = os.path.join

def select_optimal_components_publication(bids_dir, images_dir, component_range=range(2, 41),
                                          alpha_dic=10, smoothing=6, TR=2, n_top_candidates=3):
    """
    Automatically select optimal number of components with publication-quality visualization.

    This function evaluates multiple complementary metrics to identify not just one, but
    several candidate solutions that represent different trade-offs in the atlas construction.

    Parameters:
    -----------
    n_top_candidates : int
        Number of top candidate solutions to identify (default: 3)

    Returns:
    --------
    dict: Results containing optimal n_components candidates and all metrics
    """
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import MinMaxScaler
    from scipy.signal import find_peaks
    from scipy.ndimage import gaussian_filter1d
    import numpy as np
    import pandas as pd

    output_results = opj(bids_dir, 'Group_Stats')
    mask_path = opj(output_results, 'mask_mean_func_overlapp.nii.gz')

    results = {
        'n_components': [],
        'explained_variance': [],
        'explained_variance_ratio': [],
        'reconstruction_error': [],
        'sparsity': [],
        'stability': [],
        'spatial_smoothness': [],
        'silhouette': [],
        'component_overlap': []
    }

    for n_comp in component_range:
        print(f"\n--- Testing {n_comp} components ---")

        # Fit dictionary learning
        dict_learning = DictLearning(
            mask=mask_path,
            n_components=n_comp,
            alpha=alpha_dic,
            n_epochs=1,
            standardize="zscore_sample",
            random_state=0,
            n_jobs=-1,
            smoothing_fwhm=smoothing,
            t_r=TR
        )
        dict_learning.fit(images_dir)

        # 1. EXPLAINED VARIANCE
        # Cumulative variance explained by components
        per_component_var = dict_learning.score(images_dir, per_component=True)
        total_var = np.sum(per_component_var)

        # Calculate incremental variance (diminishing returns)
        if len(results['explained_variance']) > 0:
            prev_var = results['explained_variance'][-1]
            var_ratio = (total_var - prev_var) / prev_var if prev_var > 0 else 0
        else:
            var_ratio = 0

        # 2. RECONSTRUCTION ERROR
        # How well can we reconstruct original data?
        transformed = dict_learning.transform(images_dir)
        mse_list = []
        for i, img_path in enumerate(images_dir):
            original = dict_learning.masker_.transform(img_path)
            reconstructed = np.dot(transformed[i:i + 1], dict_learning.components_)
            mse = np.mean((original - reconstructed) ** 2)
            mse_list.append(mse)
        avg_mse = np.mean(mse_list)

        # 3. SPARSITY (component distinctiveness)
        # Components should be sparse (activate in limited regions)
        # Higher sparsity = more localized, interpretable regions
        component_sparsity = np.mean(np.sum(np.abs(dict_learning.components_) > 0.1, axis=1))

        # 4. SPATIAL SMOOTHNESS
        # Biological plausibility - avoid salt-and-pepper noise
        components_img = dict_learning.masker_.inverse_transform(dict_learning.components_)
        components_data = components_img.get_fdata()

        smoothness_scores = []
        for comp_idx in range(n_comp):
            comp_data = components_data[:, :, :, comp_idx]
            # Lower gradient = smoother = more biologically plausible
            grad_mag = np.sqrt(
                np.sum(np.gradient(comp_data)[i] ** 2 for i in range(3))
            )
            smoothness_scores.append(np.mean(grad_mag))
        avg_smoothness = np.mean(smoothness_scores)

        # 5. STABILITY (split-half reliability)
        # Do we get similar components with different data subsets?
        # Critical for reproducibility
        n_imgs = len(images_dir)
        half = n_imgs // 2

        dict_learning_1 = DictLearning(
            mask=mask_path, n_components=n_comp, alpha=alpha_dic,
            n_epochs=1, standardize="zscore_sample", random_state=0,
            smoothing_fwhm=smoothing, t_r=TR, n_jobs=-1
        )
        dict_learning_1.fit(images_dir[:half])

        dict_learning_2 = DictLearning(
            mask=mask_path, n_components=n_comp, alpha=alpha_dic,
            n_epochs=1, standardize="zscore_sample", random_state=1,
            smoothing_fwhm=smoothing, t_r=TR, n_jobs=-1
        )
        dict_learning_2.fit(images_dir[half:])

        # Hungarian algorithm would be better, but max correlation is simpler
        corr_matrix = np.corrcoef(dict_learning_1.components_, dict_learning_2.components_)[:n_comp, n_comp:]
        # For each component in set 1, find best match in set 2
        stability = np.mean(np.max(np.abs(corr_matrix), axis=1))

        # 6. SILHOUETTE SCORE (clustering quality)
        # How well-separated are the functional networks?
        try:
            labels = np.argmax(np.abs(transformed), axis=1)
            unique_labels = np.unique(labels)
            if len(unique_labels) > 1 and len(unique_labels) < len(labels) - 1:
                silhouette = silhouette_score(transformed, labels, metric='euclidean',
                                              sample_size=min(1000, len(labels)))
            else:
                silhouette = 0
        except:
            silhouette = 0

        # 7. COMPONENT OVERLAP
        # Measure spatial overlap between components (lower is better)
        overlap_scores = []
        for i in range(n_comp):
            for j in range(i + 1, n_comp):
                comp_i = np.abs(components_data[:, :, :, i]) > 0.1
                comp_j = np.abs(components_data[:, :, :, j]) > 0.1
                intersection = np.sum(comp_i & comp_j)
                union = np.sum(comp_i | comp_j)
                if union > 0:
                    overlap_scores.append(intersection / union)
        avg_overlap = np.mean(overlap_scores) if overlap_scores else 0

        # Store results
        results['n_components'].append(n_comp)
        results['explained_variance'].append(total_var)
        results['explained_variance_ratio'].append(var_ratio)
        results['reconstruction_error'].append(avg_mse)
        results['sparsity'].append(component_sparsity)
        results['stability'].append(stability)
        results['spatial_smoothness'].append(avg_smoothness)
        results['silhouette'].append(silhouette)
        results['component_overlap'].append(avg_overlap)

        print(f"Variance: {total_var:.4f}, Stability: {stability:.4f}, Silhouette: {silhouette:.4f}")

    # Convert to DataFrame
    metrics_df = pd.DataFrame(results)

    # === INTELLIGENT COMPONENT SELECTION ===
    # We identify multiple candidates based on different criteria

    # 1. Normalize metrics to [0, 1]
    scaler = MinMaxScaler()
    normalized = pd.DataFrame()
    normalized['n_components'] = metrics_df['n_components']

    # Invert metrics where lower is better
    metrics_df['reconstruction_error_inv'] = 1 / (1 + metrics_df['reconstruction_error'])
    metrics_df['smoothness_score'] = 1 / (1 + metrics_df['spatial_smoothness'])
    metrics_df['overlap_inv'] = 1 - metrics_df['component_overlap']

    # Normalize all relevant metrics
    for col in ['explained_variance', 'reconstruction_error_inv', 'sparsity',
                'stability', 'smoothness_score', 'silhouette', 'overlap_inv']:
        if col in metrics_df.columns:
            normalized[col] = scaler.fit_transform(metrics_df[[col]])

    # 2. Calculate multiple composite scores with different philosophies

    # BALANCED SCORE: Equal weight to all criteria
    normalized['balanced_score'] = normalized[[
        'explained_variance', 'reconstruction_error_inv', 'stability',
        'sparsity', 'smoothness_score', 'silhouette', 'overlap_inv'
    ]].mean(axis=1)

    # INTERPRETABILITY SCORE: Prioritize stable, distinct, smooth components
    weights_interp = {
        'stability': 0.30,
        'sparsity': 0.20,
        'smoothness_score': 0.20,
        'overlap_inv': 0.15,
        'silhouette': 0.10,
        'explained_variance': 0.05
    }
    normalized['interpretability_score'] = sum(
        normalized[col] * weight for col, weight in weights_interp.items()
    )

    # PERFORMANCE SCORE: Prioritize variance explained and reconstruction
    weights_perf = {
        'explained_variance': 0.40,
        'reconstruction_error_inv': 0.30,
        'stability': 0.15,
        'silhouette': 0.10,
        'sparsity': 0.05
    }
    normalized['performance_score'] = sum(
        normalized[col] * weight for col, weight in weights_perf.items()
    )

    # 3. ELBOW DETECTION on explained variance
    variances = np.array(results['explained_variance'])
    n_comps = np.array(results['n_components'])

    # Smooth the curve
    variances_smooth = gaussian_filter1d(variances, sigma=1.5)

    # Calculate curvature (second derivative)
    first_deriv = np.gradient(variances_smooth)
    second_deriv = np.gradient(first_deriv)

    # Find elbow as point of maximum curvature
    elbow_idx = np.argmax(-second_deriv[2:-2]) + 2  # Avoid edges
    elbow_n = n_comps[elbow_idx]

    # 4. KNEE DETECTION (Kneedle algorithm approximation)
    # Normalize curve to [0,1] x [0,1]
    x_norm = (n_comps - n_comps.min()) / (n_comps.max() - n_comps.min())
    y_norm = (variances - variances.min()) / (variances.max() - variances.min())

    # Distance from diagonal (y = x line)
    distances = y_norm - x_norm
    knee_idx = np.argmax(distances)
    knee_n = n_comps[knee_idx]

    # 5. Find top candidates from each score
    candidates = {}

    candidates['balanced'] = int(normalized.loc[normalized['balanced_score'].idxmax(), 'n_components'])
    candidates['interpretability'] = int(normalized.loc[normalized['interpretability_score'].idxmax(), 'n_components'])
    candidates['performance'] = int(normalized.loc[normalized['performance_score'].idxmax(), 'n_components'])
    candidates['elbow'] = int(elbow_n)
    candidates['knee'] = int(knee_n)

    # Find stability plateau (where stability stops improving significantly)
    stability_arr = np.array(results['stability'])
    stability_diff = np.diff(stability_arr)
    # Find where improvement drops below threshold
    plateau_idx = np.where(stability_diff < 0.01)[0]
    if len(plateau_idx) > 0:
        candidates['stability_plateau'] = int(n_comps[plateau_idx[0]])
    else:
        candidates['stability_plateau'] = int(n_comps[-1])

    # Get unique candidates sorted
    unique_candidates = sorted(list(set(candidates.values())))

    # Rank all candidates by composite balanced score
    candidate_scores = []
    for n in unique_candidates:
        idx = normalized[normalized['n_components'] == n].index[0]
        score = normalized.loc[idx, 'balanced_score']
        candidate_scores.append((n, score))

    candidate_scores.sort(key=lambda x: x[1], reverse=True)
    top_candidates = [x[0] for x in candidate_scores[:n_top_candidates]]

    print(f"\n{'=' * 60}")
    print(f"CANDIDATE SELECTION RESULTS:")
    print(f"{'=' * 60}")
    print(f"\nTop {n_top_candidates} candidates (by balanced score):")
    for i, n in enumerate(top_candidates, 1):
        idx = normalized[normalized['n_components'] == n].index[0]
        score = normalized.loc[idx, 'balanced_score']
        print(f"  {i}. n={n:2d} components (balanced score: {score:.3f})")

    print(f"\nDetection methods:")
    for method, n in candidates.items():
        print(f"  {method:20s}: n={n:2d}")
    print(f"{'=' * 60}\n")

    # === CREATE PUBLICATION-QUALITY FIGURE ===
    create_elife_style_figure(results, normalized, candidates, top_candidates,
                              output_results, elbow_idx, knee_idx)

    # Save detailed results
    metrics_df.to_csv(opj(output_results, 'component_selection_metrics.csv'), index=False)
    normalized.to_csv(opj(output_results, 'component_selection_normalized.csv'), index=False)

    # Save candidate summary
    summary = pd.DataFrame({
        'rank': range(1, len(top_candidates) + 1),
        'n_components': top_candidates,
        'balanced_score': [normalized.loc[normalized['n_components'] == n, 'balanced_score'].values[0]
                           for n in top_candidates],
        'stability': [results['stability'][results['n_components'].index(n)]
                      for n in top_candidates],
        'explained_variance': [results['explained_variance'][results['n_components'].index(n)]
                               for n in top_candidates]
    })
    summary.to_csv(opj(output_results, 'top_candidates_summary.csv'), index=False)

    return {
        'top_candidates': top_candidates,
        'all_candidates': candidates,
        'all_metrics': metrics_df,
        'normalized_metrics': normalized,
        'recommended': top_candidates[0]
    }


def create_elife_style_figure(results, normalized, candidates, top_candidates,
                              output_dir, elbow_idx, knee_idx):
    """
    Create publication-quality eLife-style figure for component selection.

    eLife style characteristics:
    - Clean, minimal design
    - Sans-serif fonts (Arial/Helvetica)
    - High contrast
    - Clear panel labels (A, B, C, etc.)
    - Colorblind-friendly palette
    """
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    import matplotlib.patches as mpatches

    # Set eLife style parameters
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    rcParams['font.size'] = 8
    rcParams['axes.labelsize'] = 9
    rcParams['axes.titlesize'] = 10
    rcParams['xtick.labelsize'] = 8
    rcParams['ytick.labelsize'] = 8
    rcParams['legend.fontsize'] = 7
    rcParams['axes.linewidth'] = 0.8
    rcParams['xtick.major.width'] = 0.8
    rcParams['ytick.major.width'] = 0.8
    rcParams['xtick.major.size'] = 3
    rcParams['ytick.major.size'] = 3

    # Colorblind-friendly palette (Wong 2011)
    colors = {
        'primary': '#0173B2',  # Blue
        'secondary': '#DE8F05',  # Orange
        'tertiary': '#029E73',  # Green
        'quaternary': '#CC78BC',  # Purple
        'candidate1': '#CA9161',  # Tan
        'candidate2': '#949494',  # Gray
        'candidate3': '#ECE133',  # Yellow
        'elbow': '#D55E00',  # Vermillion
        'error': '#999999'  # Light gray
    }

    # Create figure with specific dimensions (eLife typically uses 85mm or 180mm width)
    fig = plt.figure(figsize=(7.08, 8))  # 180mm width, ~200mm height

    # Create grid for subplots
    from matplotlib.gridspec import GridSpec
    gs = GridSpec(4, 2, figure=fig, hspace=0.4, wspace=0.35,
                  left=0.08, right=0.96, top=0.95, bottom=0.06)

    n_comps = np.array(results['n_components'])

    # === PANEL A: Explained Variance with Elbow Detection ===
    ax_a = fig.add_subplot(gs[0, :])

    variances = np.array(results['explained_variance'])
    ax_a.plot(n_comps, variances, 'o-', color=colors['primary'],
              linewidth=1.5, markersize=4, label='Explained variance', zorder=3)

    # Show elbow point
    ax_a.plot(n_comps[elbow_idx], variances[elbow_idx], 'D',
              color=colors['elbow'], markersize=8, label=f'Elbow (n={n_comps[elbow_idx]})',
              zorder=4, markeredgewidth=1, markeredgecolor='white')

    # Highlight top candidates
    for i, n in enumerate(top_candidates[:3]):
        idx = list(n_comps).index(n)
        marker_colors = [colors['candidate1'], colors['candidate2'], colors['candidate3']]
        ax_a.plot(n, variances[idx], 's', color=marker_colors[i],
                  markersize=7, zorder=5, markeredgewidth=1, markeredgecolor='white',
                  label=f'Candidate {i + 1} (n={n})')

    ax_a.set_xlabel('Number of components', fontweight='bold')
    ax_a.set_ylabel('Explained variance', fontweight='bold')
    ax_a.set_xlim([n_comps.min() - 1, n_comps.max() + 1])
    ax_a.grid(True, alpha=0.3, linewidth=0.5)
    ax_a.legend(loc='lower right', frameon=True, fancybox=False,
                edgecolor='black', framealpha=1)

    # Panel label
    ax_a.text(-0.08, 1.05, 'A', transform=ax_a.transAxes,
              fontsize=12, fontweight='bold', va='top')

    # === PANEL B: Reconstruction Error ===
    ax_b = fig.add_subplot(gs[1, 0])

    errors = np.array(results['reconstruction_error'])
    ax_b.plot(n_comps, errors, 'o-', color=colors['secondary'],
              linewidth=1.5, markersize=4, zorder=3)

    # Highlight top candidates
    for i, n in enumerate(top_candidates[:3]):
        idx = list(n_comps).index(n)
        marker_colors = [colors['candidate1'], colors['candidate2'], colors['candidate3']]
        ax_b.plot(n, errors[idx], 's', color=marker_colors[i],
                  markersize=6, zorder=5, markeredgewidth=1, markeredgecolor='white')

    ax_b.set_xlabel('Number of components', fontweight='bold')
    ax_b.set_ylabel('Reconstruction error (MSE)', fontweight='bold')
    ax_b.set_xlim([n_comps.min() - 1, n_comps.max() + 1])
    ax_b.grid(True, alpha=0.3, linewidth=0.5)

    ax_b.text(-0.15, 1.05, 'B', transform=ax_b.transAxes,
              fontsize=12, fontweight='bold', va='top')

    # === PANEL C: Stability (Split-Half Reliability) ===
    ax_c = fig.add_subplot(gs[1, 1])

    stability = np.array(results['stability'])
    ax_c.plot(n_comps, stability, 'o-', color=colors['tertiary'],
              linewidth=1.5, markersize=4, zorder=3)

    # Add stability threshold line
    ax_c.axhline(y=0.7, color=colors['error'], linestyle='--',
                 linewidth=1, alpha=0.7, label='Good stability (>0.7)')

    # Highlight top candidates
    for i, n in enumerate(top_candidates[:3]):
        idx = list(n_comps).index(n)
        marker_colors = [colors['candidate1'], colors['candidate2'], colors['candidate3']]
        ax_c.plot(n, stability[idx], 's', color=marker_colors[i],
                  markersize=6, zorder=5, markeredgewidth=1, markeredgecolor='white')

    ax_c.set_xlabel('Number of components', fontweight='bold')
    ax_c.set_ylabel('Stability (split-half r)', fontweight='bold')
    ax_c.set_xlim([n_comps.min() - 1, n_comps.max() + 1])
    ax_c.set_ylim([0, 1.05])
    ax_c.grid(True, alpha=0.3, linewidth=0.5)
    ax_c.legend(loc='lower right', frameon=True, fancybox=False,
                edgecolor='black', framealpha=1)

    ax_c.text(-0.15, 1.05, 'C', transform=ax_c.transAxes,
              fontsize=12, fontweight='bold', va='top')

    # === PANEL D: Silhouette Score ===
    ax_d = fig.add_subplot(gs[2, 0])

    silhouette = np.array(results['silhouette'])
    ax_d.plot(n_comps, silhouette, 'o-', color=colors['quaternary'],
              linewidth=1.5, markersize=4, zorder=3)

    # Highlight top candidates
    for i, n in enumerate(top_candidates[:3]):
        idx = list(n_comps).index(n)
        marker_colors = [colors['candidate1'], colors['candidate2'], colors['candidate3']]
        ax_d.plot(n, silhouette[idx], 's', color=marker_colors[i],
                  markersize=6, zorder=5, markeredgewidth=1, markeredgecolor='white')

    ax_d.set_xlabel('Number of components', fontweight='bold')
    ax_d.set_ylabel('Silhouette score', fontweight='bold')
    ax_d.set_xlim([n_comps.min() - 1, n_comps.max() + 1])
    ax_d.grid(True, alpha=0.3, linewidth=0.5)

    ax_d.text(-0.15, 1.05, 'D', transform=ax_d.transAxes,
              fontsize=12, fontweight='bold', va='top')

    # === PANEL E: Sparsity and Overlap ===
    ax_e = fig.add_subplot(gs[2, 1])

    sparsity = np.array(results['sparsity'])
    overlap = np.array(results['component_overlap'])

    ax_e_twin = ax_e.twinx()

    line1 = ax_e.plot(n_comps, sparsity, 'o-', color=colors['primary'],
                      linewidth=1.5, markersize=4, label='Sparsity', zorder=3)
    line2 = ax_e_twin.plot(n_comps, overlap, 's-', color=colors['secondary'],
                           linewidth=1.5, markersize=4, label='Overlap', zorder=3)

    ax_e.set_xlabel('Number of components', fontweight='bold')
    ax_e.set_ylabel('Sparsity (active voxels)', fontweight='bold', color=colors['primary'])
    ax_e_twin.set_ylabel('Component overlap (IoU)', fontweight='bold', color=colors['secondary'])
    ax_e.set_xlim([n_comps.min() - 1, n_comps.max() + 1])

    ax_e.tick_params(axis='y', labelcolor=colors['primary'])
    ax_e_twin.tick_params(axis='y', labelcolor=colors['secondary'])

    ax_e.grid(True, alpha=0.3, linewidth=0.5)

    # Combined legend
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax_e.legend(lines, labels, loc='upper right', frameon=True,
                fancybox=False, edgecolor='black', framealpha=1)

    ax_e.text(-0.15, 1.05, 'E', transform=ax_e.transAxes,
              fontsize=12, fontweight='bold', va='top')

    # === PANEL F: Composite Scores Comparison ===
    ax_f = fig.add_subplot(gs[3, :])

    # Plot three different composite scores
    ax_f.plot(normalized['n_components'], normalized['balanced_score'],
              'o-', color=colors['primary'], linewidth=1.5, markersize=4,
              label='Balanced', zorder=3)
    ax_f.plot(normalized['n_components'], normalized['interpretability_score'],
              's-', color=colors['tertiary'], linewidth=1.5, markersize=4,
              label='Interpretability', zorder=3)
    ax_f.plot(normalized['n_components'], normalized['performance_score'],
              '^-', color=colors['secondary'], linewidth=1.5, markersize=4,
              label='Performance', zorder=3)

    # Highlight top 3 candidates with vertical spans
    colors_span = [colors['candidate1'], colors['candidate2'], colors['candidate3']]
    for i, n in enumerate(top_candidates[:3]):
        ax_f.axvline(x=n, color=colors_span[i], linestyle='--',
                     linewidth=1.5, alpha=0.6, zorder=2)
        ax_f.text(n, 0.95, f'{n}', ha='center', va='top',
                  fontsize=7, color=colors_span[i], fontweight='bold',
                  transform=ax_f.get_xaxis_transform())

    ax_f.set_xlabel('Number of components', fontweight='bold')
    ax_f.set_ylabel('Normalized composite score', fontweight='bold')
    ax_f.set_xlim([n_comps.min() - 1, n_comps.max() + 1])
    ax_f.set_ylim([0, 1.05])
    ax_f.grid(True, alpha=0.3, linewidth=0.5)
    ax_f.legend(loc='lower right', frameon=True, fancybox=False,
                edgecolor='black', framealpha=1, ncol=3)

    ax_f.text(-0.04, 1.05, 'F', transform=ax_f.transAxes,
              fontsize=12, fontweight='bold', va='top')

    # Save figure
    plt.savefig(opj(output_dir, 'component_selection_elife_style.pdf'),
                dpi=300, bbox_inches='tight')
    plt.savefig(opj(output_dir, 'component_selection_elife_style.png'),
                dpi=300, bbox_inches='tight')

    print(f"Publication figure saved to: {output_dir}")

    plt.close()