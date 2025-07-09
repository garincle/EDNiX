import re
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
from math import pi
from collections import defaultdict
from scipy import stats

# Set matplotlib to non-interactive backend
plt.switch_backend('Agg')

# Define base path and output path
base_path = '/srv/projects/easymribrain/data/MRI/'
output_path = '/srv/projects/easymribrain/data/MRI/QC_func_inter-species'

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)
print(f"Output directory created at: {output_path}")


def process_numeric_value(value):
    """Convert values to numeric, handling lists and other formats"""
    if isinstance(value, list):
        return float(value[0]) if len(value) == 1 else np.mean(value)
    elif isinstance(value, (int, float)):
        return float(value)
    elif isinstance(value, str) and value.lower() in ['nan', 'null', 'none']:
        return np.nan
    else:
        try:
            return float(value)
        except (ValueError, TypeError):
            return np.nan


def get_base_name(file_path):
    """Extract the base name for file pairing"""
    return os.path.dirname(os.path.dirname(file_path))


def extract_run_number(filepath):
    """
    Extracts run identifier from filename, normalizing numbers and handling complex run IDs.
    Returns 'run-1' if no run found.
    """
    filename = os.path.basename(filepath)
    # Find run entity: 'run-' followed by anything non-underscore
    match = re.search(r'run-([^_]+)', filename)
    if not match:
        return 'run-1'  # Default to run-1 if no run found

    run_val = match.group(1)
    # Normalize numbers like '01' -> '1', and embedded ones like 'RL03' -> 'RL3'
    run_val = re.sub(r'\b0+(\d)\b', r'\1', run_val)  # e.g. "01" -> "1"
    run_val = re.sub(r'(\D)0+(\d)', r'\1\2', run_val)  # e.g. "RL03" -> "RL3"

    return f"run-{run_val}"


def extract_pairing_key(filepath):
    filename = os.path.basename(filepath)

    # Match BIDS-style entities: sub-XX, ses-YY, task-*, run-*, acq-*, bold, etc.
    entities = re.findall(r'(sub-[^_]+|ses-[^_]+|task-[^_]+|acq-[^_]+|run-[^_]+|bold)', filename)

    normalized = []
    run_found = False

    for e in entities:
        parts = e.split('-', maxsplit=1)  # Only split at the first dash
        if len(parts) != 2:
            continue  # skip malformed entries
        key, val = parts

        # Normalize numerical parts like 01 -> 1, and embedded ones like RL03 -> RL3
        val = re.sub(r'\b0+(\d)\b', r'\1', val)  # e.g. "01" -> "1"
        val = re.sub(r'(\D)0+(\d)', r'\1\2', val)  # e.g. "RL03" -> "RL3"

        if key == 'run':
            run_found = True

        normalized.append(f"{key}-{val}")

    # If no run entity found, add default run-1
    if not run_found:
        normalized.append('run-1')

    # Sort entities to allow flexible ordering
    pairing_key = '_'.join(sorted(normalized))
    return pairing_key if normalized else None


def files_match(key1, key2):
    return key1 == key2 or key1 in key2 or key2 in key1


# Get all functional QC JSON files
print("Scanning for QC files...")
QC_files = glob.glob(f"{base_path}/**/**/**/**/func/01_prepro/01_funcspace/10_Results/fMRI_QC_SNR/*_QC_values.json",
                     recursive=False)
full_results_files = glob.glob(
    f"{base_path}/**/**/**/**/func/01_prepro/01_funcspace/10_Results/fMRI_QC_matrix/*_full_results.json",
    recursive=False)
print(f"Found {len(QC_files)} QC values files and {len(full_results_files)} full results files")

# Create a mapping between QC values and full results files
print("Pairing QC files with full results...")
from collections import defaultdict

file_pairs = defaultdict(dict)
unmatched_qc = []
unmatched_fr = []

for qc_file in QC_files:
    base_name = get_base_name(qc_file)
    pairing_key = extract_pairing_key(qc_file)
    if base_name and pairing_key:
        file_pairs[(pairing_key, base_name)]['qc_values'] = qc_file
        file_pairs[(pairing_key, base_name)]['base_name'] = base_name
    else:
        unmatched_qc.append(qc_file)
        print(f"[WARNING] QC file could not be used for pairing: {qc_file}")

for fr_file in full_results_files:
    base_name = get_base_name(fr_file)
    pairing_key = extract_pairing_key(fr_file)
    if base_name and pairing_key:
        file_pairs[(pairing_key, base_name)]['full_results'] = fr_file
        if 'base_name' not in file_pairs[(pairing_key, base_name)]:
            file_pairs[(pairing_key, base_name)]['base_name'] = base_name
    else:
        unmatched_fr.append(fr_file)
        print(f"[WARNING] Full results file could not be used for pairing: {fr_file}")

# Filter only complete pairs
complete_pairs = {k: v for k, v in file_pairs.items() if 'qc_values' in v and 'full_results' in v}
incomplete_pairs = {k: v for k, v in file_pairs.items() if not ('qc_values' in v and 'full_results' in v)}

print(f"\n[SUMMARY]")
print(f"✔️  Complete file pairs found: {len(complete_pairs)}")
print(f"❌  Incomplete pairs (missing QC or full result): {len(incomplete_pairs)}")
print(f"❌  Unmatched QC files (missing key or base): {len(unmatched_qc)}")
print(f"❌  Unmatched full results files: {len(unmatched_fr)}")

if incomplete_pairs:
    print("\nExample incomplete pairs:")
    for k, v in list(incomplete_pairs.items())[:5]:  # show only first 5
        print(f"  Pairing key: {k}")
        print(f"    QC file: {v.get('qc_values', '❌ Missing')}")
        print(f"    Full results: {v.get('full_results', '❌ Missing')}")

print("\nPairing process completed.")

if incomplete_pairs:
    print("\nExample incomplete pairs:")
    for k, v in list(incomplete_pairs.items())[:5]:  # show only first 5
        print(f"  Pairing key: {k}")
        print(f"    QC file: {v.get('qc_values', '❌ Missing')}")
        print(f"    Full results: {v.get('full_results', '❌ Missing')}")

print("\nPairing process completed.")

data = []
# Process each pair of files
print("\nProcessing files...")
for i, (base_name, files) in enumerate(complete_pairs.items()):
    if i % 10 == 0 or i == len(complete_pairs) - 1:
        print(f"Processing file pair {i + 1}/{len(complete_pairs)}: {base_name}")

    try:

        path_parts = files['qc_values'].split('/')
        species = path_parts[6]
        bids_dataset = path_parts[7]
        subject = path_parts[8]
        session = path_parts[9] if len(path_parts) > 9 else None
        task = path_parts[10] if len(path_parts) > 10 else None

        # Extract run number from qc_values filename
        run_number = extract_run_number(files['qc_values'])

        metrics = {
            'species': species,
            'bids_dataset': bids_dataset,
            'subject': subject,
            'session': session,
            'task': task,
            'run': run_number
        }

        # Load QC values JSON
        with open(files['qc_values'], 'r') as f:
            qc_values = json.load(f)
            for k, v in qc_values.items():
                metrics[k] = process_numeric_value(v)

        # Load full results JSON
        with open(files['full_results'], 'r') as f:
            full_results = json.load(f)

            for k, v in full_results['network_metrics'].items():
                metrics[f'network_{k}'] = process_numeric_value(v)

            for k, v in full_results['hemisphere_results'].items():
                if isinstance(v, list):
                    metrics[f'hemisphere_{k}_mean'] = np.mean([process_numeric_value(x) for x in v]) if v else np.nan
                else:
                    metrics[f'hemisphere_{k}'] = process_numeric_value(v)

            if 'specificity_results' in full_results:
                for k, v in full_results['specificity_results'].items():
                    if k == 'target_specificity':
                        for tk, tv in v.items():
                            metrics[f'specificity_{tk}'] = process_numeric_value(tv)
                    elif isinstance(v, (list, dict)):
                        continue
                    else:
                        metrics[f'specificity_{k}'] = process_numeric_value(v)

        data.append(metrics)
    except Exception as e:
        print(f"Error processing {base_name}: {str(e)}")
        continue

# Convert to DataFrame
print("\nCreating DataFrame...")
if not data:
    print("ERROR: No data was processed. Check file contents and structure.")
    exit(1)

df = pd.DataFrame(data)
df = df[df['bids_dataset'] != 'BIDS_Ranft']
df = df[df['bids_dataset'] != 'ARITEP-PNH']
df = df[df['bids_dataset'] != 'BIDS_merg']
df = df[df['bids_dataset'] != 'BIDS_GdGSR']
df = df[df['bids_dataset'] != 'ds004856']

# Define metrics to plot
exclude_cols = ['species', 'bids_dataset', 'subject', 'session', 'task',
                'Comp modal', 'anat_image_type', 'specificity_Specific_ROI_pair',
                'specificity_NonSpecific_ROI_pair', 'specificity_Category',
                'timestamp']
plot_metrics = [col for col in df.columns if col not in exclude_cols and pd.api.types.is_numeric_dtype(df[col])]

if not plot_metrics:
    print("ERROR: No numeric metrics found to plot. Check data processing.")
    exit(1)

print(f"\nFound {len(plot_metrics)} numeric metrics to plot:")


# =============================================
# 1. Add Specificity Classification and Delta
# =============================================
def classify_specificity(row):
    try:
        specific = row.get('specificity_Specific_Correlation', np.nan)
        nonspecific = row.get('specificity_NonSpecific_Correlation', np.nan)
        specific_thresh = 0.2
        delta_thresh = 0.2

        if pd.isna(specific) or pd.isna(nonspecific):
            return "Missing"

        delta = specific - nonspecific

        if delta > delta_thresh and specific > specific_thresh:
            return "Specific"
        elif specific < specific_thresh:
            return "No"
        elif delta <= delta_thresh:
            return "Non-Specific"
        else:
            return "Spurious"
    except Exception as e:
        return "Error"


def calculate_delta(row):
    try:
        specific = row.get('specificity_Specific_Correlation', np.nan)
        nonspecific = row.get('specificity_NonSpecific_Correlation', np.nan)
        if pd.isna(specific) or pd.isna(nonspecific):
            return np.nan
        return specific - nonspecific
    except Exception as e:
        return np.nan


df['Specificity_Class'] = df.apply(classify_specificity, axis=1)
df['delta'] = df.apply(calculate_delta, axis=1)

# Add delta to plot_metrics if it's not already there
if 'delta' not in plot_metrics and 'delta' in df.columns:
    plot_metrics.append('delta')

# Save to CSV
output_csv = os.path.join(output_path, "QC_results_summary.csv")
# Step 1: Replace empty strings with NaN (to treat them the same)
df = df.replace("", np.nan)
# Step 2: Drop columns with ANY NaN (now includes original NaN + former empty strings)
df = df.dropna(axis=1)

df.to_csv(output_csv, index=False)
print(f"\nData saved to {output_csv}")
print(f"Data shape: {df.shape}")
print("Sample data:\n", df.head())

# Create output subdirectories
os.makedirs(os.path.join(output_path, "violin_plots"), exist_ok=True)
os.makedirs(os.path.join(output_path, "radar_plots"), exist_ok=True)
os.makedirs(os.path.join(output_path, "correlation_plots"), exist_ok=True)
os.makedirs(os.path.join(output_path, "species_heatmaps"), exist_ok=True)
os.makedirs(os.path.join(output_path, "swarm_plots"), exist_ok=True)
os.makedirs(os.path.join(output_path, "specificty_plot"), exist_ok=True)

# =============================================
# 5. Create Specificity Distribution Plots
# =============================================
print("\nCreating specificity distribution plots...")


def plot_specificity(data, group_name, plot_type="percentage"):
    counts = data['Specificity_Class'].value_counts(normalize=(plot_type == "percentage"))
    if plot_type == "percentage":
        counts *= 100
        ylabel = "Percentage"
    else:
        ylabel = "Count"

    colors = {
        "Specific": "green",
        "Non-Specific": "blue",
        "No": "red",
        "Spurious": "orange",
        "Missing": "gray",
        "Error": "black"
    }

    plt.figure(figsize=(10, 6))
    counts.plot(kind='bar', color=[colors.get(x, 'gray') for x in counts.index])
    plt.title(f"Specificity Distribution - {group_name} ({plot_type})")
    plt.ylabel(ylabel)
    plt.xlabel("Specificity Class")
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    filename = f"specificity_{group_name.lower().replace(' ', '_')}_{plot_type}.png"
    plt.savefig(os.path.join(output_path, "specificty_plot", filename), dpi=300, bbox_inches='tight')
    plt.close()


# Plot for each species and dataset
for species in df['species'].unique():
    species_data = df[df['species'] == species]
    plot_specificity(species_data, f"Species: {species}", "percentage")
    plot_specificity(species_data, f"Species: {species}", "count")

for dataset in df['bids_dataset'].unique():
    dataset_data = df[df['bids_dataset'] == dataset]
    plot_specificity(dataset_data, f"Dataset: {dataset}", "percentage")
    plot_specificity(dataset_data, f"Dataset: {dataset}", "count")

# =============================================
# 1. Violin + Swarm Plots for All Metrics
# =============================================
print("\nCreating violin + swarm plots...")
for metric in plot_metrics:
    try:
        print(f"  Creating plot for {metric}...")
        plt.figure(figsize=(12, 8))

        # Violin plot with quartile lines
        ax = sns.violinplot(data=df, x='species', y=metric,
                            inner='quartile', palette='Set2')

        # Add swarm plot for individual data points
        sns.swarmplot(data=df, x='species', y=metric,
                      color='black', alpha=0.5, size=3, ax=ax)

        plt.title(f"Distribution of {metric} by Species", pad=20)
        plt.xlabel('Species')
        plt.ylabel(metric)
        plt.xticks(rotation=45)
        plt.tight_layout()

        # Save plot
        plot_path = os.path.join(output_path, "violin_plots", f"{metric}_violin_swarm.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()

    except Exception as e:
        print(f"  Could not create violin plot for {metric}: {str(e)}")

print(f"Violin plots saved to {'violin_plots'}")

# =============================================
# 2. Radar (Spider) Plots for All Metrics
# =============================================
print("\nCreating radar plots for all metrics...")

# Group metrics by category for better radar plot organization
metric_categories = {
    'Signal Quality': ['avg_snr_gray', 'avg_snr_white', 'cnr', 'cortical_contrast',
                       'wm_gm_contrast', 'TSNR', 'TSNRcvarinv'],
    'Motion': ['mean_fd', 'mean_dvars', 'gcor', 'ghost_ratio', 'avg_enorm'],
    'Spatial': ['jaccard', 'nmi', 'mi', 'cc']}

# Create radar plots for each category
for category_name, metrics in metric_categories.items():
    if not metrics:
        continue

    print(f"  Creating radar plot for {category_name} ({len(metrics)} metrics)...")

    # Filter metrics that exist in the data
    existing_metrics = [m for m in metrics if m in df.columns]
    if len(existing_metrics) < 3:  # Need at least 3 metrics for radar plot
        print(f"    Skipping {category_name} - not enough metrics (found {len(existing_metrics)})")
        continue

    # Prepare data - mean values per species
    radar_df = df.groupby('species')[existing_metrics].mean().reset_index()

    # Normalize each metric to 0-1 scale for better visualization
    for metric in existing_metrics:
        min_val = radar_df[metric].min()
        max_val = radar_df[metric].max()
        if max_val != min_val:  # Avoid division by zero
            radar_df[metric] = (radar_df[metric] - min_val) / (max_val - min_val)
        else:
            radar_df[metric] = 0.5  # If all values are identical

    # Number of variables we're plotting
    num_vars = len(existing_metrics)

    # Compute angle of each axis
    angles = [n / float(num_vars) * 2 * pi for n in range(num_vars)]
    angles += angles[:1]  # Close the loop

    # Initialize radar plot
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'polar': True})

    # Draw one axis per variable and add labels
    plt.xticks(angles[:-1], existing_metrics, color='grey', size=10)

    # Draw ylabels
    ax.set_rlabel_position(30)
    plt.yticks([0.25, 0.5, 0.75], ["0.25", "0.5", "0.75"], color="grey", size=8)
    plt.ylim(0, 1)

    # Plot data for each species
    for idx, row in radar_df.iterrows():
        values = row[existing_metrics].tolist()
        values += values[:1]  # Close the loop
        ax.plot(angles, values, linewidth=2, linestyle='solid',
                label=row['species'])
        ax.fill(angles, values, alpha=0.1)

    # Add legend and title
    plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    plt.title(f"{category_name} Metrics Comparison", size=15, y=1.1)

    # Save plot
    plot_path = os.path.join(output_path, "radar_plots", f"radar_{category_name.lower().replace(' ', '_')}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
print(f"Radar plots saved to {os.path.join(output_path, 'radar_plots')}")
# =============================================
# 6. Enhanced Swarm Plots with Specificity
# =============================================
print("\nCreating enhanced swarm plots...")


def create_swarm_plot(metric, data):
    plt.figure(figsize=(12, 8))
    palette = {
        "Specific": "green",
        "Non-Specific": "blue",
        "No": "red",
        "Spurious": "orange",
        "Missing": "gray",
        "Error": "black"
    }

    ax = sns.swarmplot(
        data=data,
        x='species',
        y=metric,
        hue='Specificity_Class',
        palette=palette,
        size=5
    )

    plt.title(f"{metric} by Species with Specificity", pad=20)
    plt.xlabel('Species')
    plt.ylabel(metric)
    plt.xticks(rotation=45)
    plt.legend(title='Specificity Class', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    plot_path = os.path.join(output_path, "swarm_plots", f"swarm_specificity_{metric}.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()


# Create swarm plots for key metrics
key_metrics = ['avg_snr_gray', 'mean_fd', 'network_Mean_Correlation', 'hemisphere_inter_mean']
for metric in key_metrics:
    if metric in df.columns:
        create_swarm_plot(metric, df)

print("\nAll visualizations completed successfully!")
print(f"Final outputs saved in: {os.path.join(output_path, 'swarm_plots')}")


# =============================================
# 2. Create Comprehensive Heatmap (All Metrics)
# =============================================
def create_heatmap(data, title, filename, columns_order=None, pvalue_threshold=None):
    """Helper function to create heatmaps"""
    # Use specified column order or default to data columns
    if columns_order is not None:
        # Only keep columns that actually exist in the data
        columns_order = [col for col in columns_order if col in data.columns]
        data = data[columns_order]

    # Calculate correlation matrix with NaN handling
    corr_matrix = data.corr()

    # Replace infinite values with NaN and drop columns/rows with all NaN
    corr_matrix = corr_matrix.replace([np.inf, -np.inf], np.nan)
    corr_matrix = corr_matrix.dropna(axis=0, how='all').dropna(axis=1, how='all')

    if pvalue_threshold is not None:
        # Calculate p-values for correlations
        n = data.shape[0]
        p_matrix = np.zeros_like(corr_matrix)
        for i in range(len(corr_matrix.columns)):
            for j in range(len(corr_matrix.columns)):
                if i != j:
                    valid = data.iloc[:, [i, j]].dropna()
                    if len(valid) >= 3:
                        _, p_matrix[i, j] = stats.pearsonr(valid.iloc[:, 0], valid.iloc[:, 1])
                    else:
                        p_matrix[i, j] = 1  # Treat as non-significant if not enough data

        # Mask non-significant correlations
        mask = p_matrix > pvalue_threshold
        corr_matrix[mask] = 0

    if corr_matrix.empty:
        print(f"Skipping {title} - no valid correlations after NaN removal")
        return None

    try:
        # Create regular heatmap with fixed order
        plt.figure(figsize=(20, 18))
        sns.heatmap(corr_matrix,
                    cmap='coolwarm',
                    annot=True,
                    fmt=".1f",
                    center=0,
                    vmin=-1,
                    vmax=1,
                    mask=np.isnan(corr_matrix) if pvalue_threshold is None else (corr_matrix == 0))
        plt.title(title, pad=20)
        plot_path = os.path.join(output_path, filename)
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        return list(corr_matrix.columns)
    except Exception as e:
        print(f"Could not create heatmap for {title}: {str(e)}")
        plt.close()
        return None


print("\nCreating comprehensive heatmap with all metrics...")
valid_metrics = [m for m in plot_metrics if m in df.columns and df[m].notna().sum() >= 3 and df[m].std() > 0]

# First create regular heatmap and get column order
col_order = create_heatmap(df[valid_metrics],
                           "All Metrics Correlation",
                           "correlation_plots/all_metrics_heatmap.png")

# Then create statistically significant heatmap using same order
if col_order:
    _ = create_heatmap(df[valid_metrics],
                       "All Metrics Correlation (p < 0.05)",
                       "correlation_plots/all_metrics_heatmap_significant.png",
                       columns_order=col_order,
                       pvalue_threshold=0.05)

# Fallback to alphabetical order if clustering failed
if not col_order:
    col_order = sorted(valid_metrics)
    print("Using alphabetical column order for consistency")

# =============================================
# 3. Create Heatmaps per Species with Consistent Order
# =============================================
print("\nCreating heatmaps per species with consistent column order...")
for species in df['species'].unique():
    species_data = df[df['species'] == species]
    # Only use columns that exist in both the data and our column order
    valid_metrics_species = [m for m in col_order if m in species_data.columns]

    if len(valid_metrics_species) >= 2:
        print(f"  Creating heatmap for {species}...")
        # Use the same column order as group heatmap
        create_heatmap(species_data[valid_metrics_species],
                       f"{species} Metrics Correlation",
                       f"species_heatmaps/{species.lower()}_heatmap.png",
                       columns_order=col_order)

print("\nCreating combined specificity distribution plots...")


def create_combined_specificity_plot(df, group_column, title_suffix):
    # Calculate percentages for each group
    grouped = df.groupby([group_column, 'Specificity_Class']).size().unstack(fill_value=0)
    percentages = grouped.div(grouped.sum(axis=1), axis=0) * 100

    # Define colors and order
    specificity_order = ["Specific", "Non-Specific", "No", "Spurious", "Missing", "Error"]
    colors = {
        "Specific": "#4CAF50",  # Green
        "Non-Specific": "#2196F3",  # Blue
        "No": "#F44336",  # Red
        "Spurious": "#FF9800",  # Orange
        "Missing": "#9E9E9E",  # Gray
        "Error": "#000000"  # Black
    }

    # Reorder columns
    percentages = percentages.reindex(columns=specificity_order)

    # Plot setup
    fig, ax = plt.subplots(figsize=(12, 6))
    bar_width = 0.8 / len(specificity_order)
    x = np.arange(len(percentages.index))

    # Plot each specificity class
    for i, spec_class in enumerate(specificity_order):
        if spec_class in percentages.columns:
            ax.bar(x + i * bar_width, percentages[spec_class],
                   width=bar_width,
                   label=spec_class,
                   color=colors[spec_class])

    # Formatting
    ax.set_title(f"Specificity Distribution by {title_suffix} (Percentage)", pad=20)
    ax.set_ylabel("Percentage (%)")
    ax.set_xlabel(title_suffix)
    ax.set_xticks(x + bar_width * (len(specificity_order) - 1) / 2)
    ax.set_xticklabels(percentages.index, rotation=45 if len(percentages.index) > 5 else 0)
    ax.legend(title="Specificity Class", bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()

    # Save plot
    filename = f"combined_specificity_{group_column.lower()}.png"
    plt.savefig(os.path.join(output_path, "specificty_plot", filename), dpi=300, bbox_inches='tight')
    plt.close()


# Create output directory
os.makedirs(os.path.join(output_path, "specificty_plot"), exist_ok=True)

# Generate combined plots
create_combined_specificity_plot(df, 'species', "Species")
create_combined_specificity_plot(df, 'bids_dataset', "BIDS Dataset")

print("Combined specificity plots created successfully!")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix

# Column names
target_col = "Specificity_Class"
predictors = [
    "avg_snr_gray",
    "avg_snr_white",
    "cnr",
    "cortical_contrast",
    "stdev",
    "TSNR",
    "mean_fd",
    "mean_dvars",
    "gcor",
    "fwhm",
    "mi",
    "cc",
    "nmi",
    "avg_enorm",
    "avg_velocity",
    "censor_fraction",
    "avg_outcount",
    "network_Davies_Bouldin",
    "network_Calinski_Harabasz",
    "network_Mean_Correlation",
    "network_Median_Correlation",
    "network_Std_Correlation",
    "hemisphere_intra_left_mean",
    "hemisphere_inter_mean",
    "hemisphere_t_intra_vs_inter",
    "hemisphere_t_left_vs_right"]

# Drop rows with missing values
df_clean = df.dropna(subset=predictors + [target_col])

# Features and target
X = df_clean[predictors]
y = df_clean[target_col]

# Optional: encode target labels (RandomForestClassifier can handle strings, but encoding is good for other classifiers too)
le = LabelEncoder()
y_encoded = le.fit_transform(y)

# Train classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
clf.fit(X, y_encoded)

# Feature importances
importances = pd.Series(clf.feature_importances_, index=predictors).sort_values(ascending=False)

# Plot feature importances
plt.figure(figsize=(10, 5))
sns.barplot(x=importances.values, y=importances.index)
plt.title("Top predictors of specificity label (Random Forest)")
plt.xlabel("Importance")
plt.tight_layout()
plt.savefig(os.path.join(output_path, "predictors_specificity.png"), dpi=300, bbox_inches='tight')
plt.close()

# Optional: Show classification report
X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.3, random_state=42, stratify=True)
clf_eval = RandomForestClassifier(n_estimators=100, random_state=42)
clf_eval.fit(X_train, y_train)
y_pred = clf_eval.predict(X_test)

print("\nClassification Report:")
print(classification_report(le.inverse_transform(y_test), le.inverse_transform(y_pred)))

# Optional: Confusion matrix
cm = pd.DataFrame(confusion_matrix(le.inverse_transform(y_test), le.inverse_transform(y_pred), normalize=True),
                  index=le.classes_, columns=le.classes_)
plt.figure(figsize=(6, 5))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
plt.title("Confusion Matrix")
plt.ylabel("True label")
plt.xlabel("Predicted label")
plt.tight_layout()
plt.savefig(os.path.join(output_path, "Classification_specificity.png"), dpi=300, bbox_inches='tight')
plt.close()

