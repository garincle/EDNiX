import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import string

# For subplot labels (A, B, C...)
# Path setup
opj = os.path.join
opd = os.path.dirname

# Load data
source_csv = "/home/cgarin/Documents/EDNiX_study/MRI/QC_func_inter-species/QC_results_summary.csv"
df = pd.read_csv(source_csv)

# Clean species names
df['species'] = df['species'].str.replace('BIDS_Cdt_Garin', '').str.replace('BIDS_BenHamed', '').str.strip()

# Select 6 key metrics to visualize
metrics = ['Median', 'Variance', 'Density (Positive Correlations)',
           'Density (Negative Correlations)', 'W-statistic', 'P-value']

# Melt dataframe
melted = df.melt(id_vars=['species'], value_vars=metrics,
                 var_name='Metric', value_name='Value')

# Set style
sns.set_style("whitegrid")
palette = sns.color_palette("husl", n_colors=len(df['species'].unique()))

# Create figure with 2 rows and 3 columns
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(24, 16))
plt.subplots_adjust(hspace=0.2, wspace=0.3, top=0.9)  # Adjusted hspace to bring rows closer

# Flatten axes for easy iteration
axes = axes.flatten()

# Create a list to store legend handles and labels
legend_handles = []
legend_labels = df['species'].unique()

# Plot each metric
for i, metric in enumerate(metrics):
    ax = axes[i]
    metric_data = melted[melted['Metric'] == metric]

    # Violin plot
    vplot = sns.violinplot(x='species', y='Value', data=metric_data,
                           inner=None, cut=0, bw=0.2, ax=ax, palette=palette)

    # Boxplot
    sns.boxplot(x='species', y='Value', data=metric_data, ax=ax,
                showcaps=False, boxprops=dict(alpha=0.5), width=0.2, color='black')

    # Strip plot instead of swarm plot
    sns.stripplot(x='species', y='Value', data=metric_data,
                  size=4, palette=palette, alpha=0.6, ax=ax, jitter=True)

    # Add subplot label (A, B, C...)
    ax.text(-0.1, 1.05, string.ascii_uppercase[i], transform=ax.transAxes,
            size=20, weight='bold')

    # Titles and labels
    ax.set_title('')  # Remove title since we have letters
    ax.set_ylabel(metric, fontsize=14)  # Changed from 'Value' to metric name
    ax.tick_params(axis='y', labelsize=12)
    ax.grid(True, alpha=0.3)

    # Remove species names from x-axis
    ax.set_xticklabels([])
    ax.set_xlabel('')  # Remove 'species' label from x-axis

    # Remove individual legends from each subplot if they exist
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

    # Store handles for the legend (use violin plot patches)
    if i == 0:  # Only need to do this once
        for j, patch in enumerate(vplot.collections):
            if j < len(legend_labels):  # Only take the first n patches
                legend_handles.append(patch)

# Create a single legend for the entire figure
fig.legend(legend_handles, legend_labels,
           title='Species',
           loc='upper center',
           bbox_to_anchor=(0.5, 0.98),  # Adjusted vertical position (0.98 instead of 1.05)
           ncol=len(legend_labels),
           fontsize=12,
           frameon=True)

# Save the figure
output_path = opj(opd(source_csv), "QC_comparison_6metrics_2rows.png")
plt.savefig(output_path, bbox_inches='tight', dpi=300)
plt.close()

print(f"Plot saved to: {output_path}")
