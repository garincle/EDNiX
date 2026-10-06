import os
import json
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
import math

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Define path operations
opj = os.path.join
opb = os.path.basename
opd = os.path.dirname
ope = os.path.exists

# Define plotting style
plt.style.use('seaborn-v0_8')
sns.set_context("notebook", font_scale=1.1)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

# Define metrics of interest
SELECTED_METRICS = [
    't_intra_vs_inter',
    't_left_vs_right',
    'Mean_Correlation',
    'Std_Correlation',
    'Silhouette_Score',
    'Davies_Bouldin']

def load_qc_data(bids_dir):
    """Load all QC JSON files from BIDS directory structure."""
    qc_files = glob.glob(f"{bids_dir}/**/**/func/01_prepro/01_funcspace/10_Results/fMRI_QC_matrix/**full_results.json")

    if not qc_files:
        raise ValueError(f"No QC JSON files found in {bids_dir}")

    all_data = []
    seen_subjects = set()  # To track seen subjects and avoid duplicates

    for qc_file in qc_files:
        try:
            path_parts = Path(qc_file).parts
            subject = next((p for p in path_parts if p.startswith("sub-")), None)
            session = next((p for p in path_parts if p.startswith("ses-")), None)

            # Create a unique identifier for each subject-session combination
            subject_session_id = f"{subject}_{session}" if session else subject

            # Skip if we've already seen this subject-session combination
            if subject_session_id in seen_subjects:
                print(f"Skipping duplicate subject-session: {subject_session_id}")
                continue

            seen_subjects.add(subject_session_id)

            with open(qc_file, 'r') as f:
                qc_data = json.load(f)

            # Create flat dictionary for this subject
            subject_data = {
                'subject': subject,
                'session': session if session else 'no_session'}

            # Extract network metrics
            if 'network_metrics' in qc_data:
                for metric, value in qc_data['network_metrics'].items():
                    subject_data[metric] = value

            # Extract hemisphere results
            if 'hemisphere_results' in qc_data:
                for metric, value in qc_data['hemisphere_results'].items():
                    if metric in ['p_intra_vs_inter', 'p_left_vs_right']:
                        continue  # Skip p-values
                    if isinstance(value, list):
                        subject_data[metric] = np.nanmean(value) if value else np.nan
                    else:
                        subject_data[metric] = value

            # Extract specificity category
            if 'specificity_results' in qc_data and 'target_specificity' in qc_data['specificity_results']:
                subject_data['Correlation_Category'] = qc_data['specificity_results']['target_specificity'].get('Category', 'Unknown')
                subject_data['Specific_Correlation'] = qc_data['specificity_results']['target_specificity'].get('Specific_Correlation', 'Unknown')
                subject_data['NonSpecific_Correlation'] = qc_data['specificity_results']['target_specificity'].get('NonSpecific_Correlation', 'Unknown')

            all_data.append(subject_data)

        except Exception as e:
            print(f"Error processing {qc_file}: {str(e)}")
            continue

    if not all_data:
        raise ValueError("No valid QC data found in any files")

    df = pd.DataFrame(all_data)

    # Ensure all selected metrics exist
    for metric in SELECTED_METRICS:
        if metric not in df.columns:
            df[metric] = np.nan

    print("\n=== Data Loading Debug Info ===")
    print(f"Found {len(df)} records")
    print("Columns found:", df.columns.tolist())
    print("\nSample data:")
    print(df[['subject', 'session'] + SELECTED_METRICS + ['Correlation_Category']].head())

    return df

def identify_outliers(series, method='iqr', threshold=1.5):
    """Identify outliers in a pandas Series."""
    if method == 'iqr':
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - threshold * iqr
        upper_bound = q3 + threshold * iqr
        return (series < lower_bound) | (series > upper_bound)
    elif method == 'zscore':
        z = np.abs(stats.zscore(series, nan_policy='omit'))
        return z > threshold
    else:
        raise ValueError(f"Unknown outlier detection method: {method}")

def create_comprehensive_report(df, output_dir):
    """Create a single-page report with all requested plots."""
    # Create figure
    fig = plt.figure(figsize=(20, 24))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.5, 1, 1])

    # --- Plot A: Distribution plots ---
    ax_a = fig.add_subplot(gs[0])
    n_metrics = len(SELECTED_METRICS)
    n_cols = 3
    n_rows = math.ceil(n_metrics / n_cols)
    sub_gs = gs[0].subgridspec(n_rows, n_cols)

    outlier_reports = {}

    for i, metric in enumerate(SELECTED_METRICS):
        ax = fig.add_subplot(sub_gs[i//n_cols, i%n_cols])
        data = df[metric].dropna()

        if len(data) > 0:
            sns.boxplot(x=data, width=0.3, fliersize=0, ax=ax)
            sns.swarmplot(x=data, color=".25", size=4, ax=ax)

            # Identify and label outliers
            outliers = identify_outliers(data)
            outlier_data = df.loc[outliers, ['subject', 'session', metric]].dropna()
            outlier_reports[metric] = outlier_data

            for _, row in outlier_data.iterrows():
                label = f"{row['subject']}_{row['session']}" if row['session'] != 'no_session' else row['subject']
                ax.annotate(label, xy=(row[metric], 0), xytext=(0, 5),
                           textcoords='offset points', ha='center', va='bottom',
                           fontsize=8, bbox=dict(boxstyle='round,pad=0.2',
                                               facecolor='white', alpha=0.8,
                                               edgecolor='gray'))

            mean_val = data.mean()
            ax.axvline(mean_val, color='red', linestyle='--', linewidth=1,
                       label=f'Mean: {mean_val:.2f}')

        ax.set_title(metric, fontsize=10)
        ax.set_xlabel('')
        ax.legend(loc='upper right', fontsize=8)

    # Remove empty subplots
    for i in range(len(SELECTED_METRICS), n_rows*n_cols):
        fig.delaxes(fig.axes[-1])

    # --- Plot B: Outlier distribution ---
    ax_b = fig.add_subplot(gs[1])

    # Calculate outlier counts per subject
    outlier_counts = []
    for metric in SELECTED_METRICS:
        data = df[metric].dropna()
        if len(data) > 0:
            outliers = identify_outliers(data)
            outlier_data = df.loc[outliers, ['subject', 'session']]
            outlier_counts.extend(outlier_data.to_dict('records'))

    if outlier_counts:
        outlier_df = pd.DataFrame(outlier_counts)

        subject_outliers = outlier_df.groupby(['subject', 'session']).size().reset_index(name='outlier_count')

        # Include subjects with 0 outliers
        all_subjects = df[['subject', 'session']].drop_duplicates()
        subject_outliers = pd.merge(all_subjects, subject_outliers,
                                   on=['subject', 'session'], how='left').fillna(0)

        # Create labels
        labels = []
        for _, row in subject_outliers.iterrows():
            label = row['subject']
            if row['session'] != 'no_session':
                label += f"\n{row['session']}"
            labels.append(label)

        summary_path = opj(output_dir, "outlier_summary_report.csv")
        subject_outliers.to_csv(summary_path, index=False)

        # Plot
        sns.barplot(x='subject', y='outlier_count', hue='session',
                    data=subject_outliers, palette='viridis', dodge=False, ax=ax_b)

        ax_b.set_title('B: Outlier Distribution Across Subjects', fontsize=12)
        ax_b.set_xlabel('Subject')
        ax_b.set_ylabel('Number of Outliers')
        plt.setp(ax_b.get_xticklabels(), rotation=45, ha='right')

        if len(subject_outliers['session'].unique()) > 1:
            ax_b.legend(title='Session', bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            ax_b.get_legend().remove()
    else:
        ax_b.text(0.5, 0.5, 'No outliers detected', ha='center', va='center')
        ax_b.set_title('B: Outlier Distribution Across Subjects', fontsize=12)

    # --- Plot C: Correlation Category Distribution ---
    ax_c = fig.add_subplot(gs[2])

    if 'Correlation_Category' in df.columns:
        category_counts = df['Correlation_Category'].value_counts()
        if not category_counts.empty:
            # Define a color palette with different colors for each category
            palette = sns.color_palette("husl", len(category_counts))

            sns.barplot(x=category_counts.index, y=category_counts.values, ax=ax_c, palette=palette)
            ax_c.set_title('C: Distribution of Correlation Categories', fontsize=12)
            ax_c.set_xlabel('Category')
            ax_c.set_ylabel('Count')

            # Add count labels
            for p in ax_c.patches:
                ax_c.annotate(f"{int(p.get_height())}",
                             (p.get_x() + p.get_width() / 2., p.get_height()),
                             ha='center', va='center', xytext=(0, 5),
                             textcoords='offset points')
        else:
            ax_c.text(0.5, 0.5, 'No category data available', ha='center', va='center')
            ax_c.set_title('C: Distribution of Correlation Categories', fontsize=12)
    else:
        ax_c.text(0.5, 0.5, 'No category data available', ha='center', va='center')
        ax_c.set_title('C: Distribution of Correlation Categories', fontsize=12)

    # Adjust layout and save
    plt.tight_layout()
    report_path = opj(output_dir, "comprehensive_qc_report.png")
    plt.savefig(report_path, dpi=300, bbox_inches='tight')
    plt.close()

    # Create and save category summary
    if 'Correlation_Category' in df.columns:
        category_summary = df[['subject', 'session', 'Correlation_Category']].copy()
        category_summary.to_csv(opj(output_dir, "subject_categories.csv"), index=False)
    else:
        pd.DataFrame(columns=['subject', 'session', 'Correlation_Category']).to_csv(
            opj(output_dir, "subject_categories.csv"), index=False)

    # Save summary statistics
    summary_stats = pd.DataFrame(index=SELECTED_METRICS,
                                columns=['mean', 'std', 'min', '25%', '50%', '75%', 'max'])
    for metric in SELECTED_METRICS:
        if metric in df.columns:
            summary_stats.loc[metric] = df[metric].describe()[['mean', 'std', 'min', '25%', '50%', '75%', 'max']]
    summary_stats.to_csv(opj(output_dir, "qc_metrics_summary.csv"))

    print("\n=== Report Generation Complete ===")
    print(f"Saved comprehensive report to: {report_path}")

