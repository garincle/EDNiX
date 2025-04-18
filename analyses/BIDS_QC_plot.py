import glob
import ast
import os
import pandas as pd
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns


#Path to the excels files and data structure
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists

# Define paths
MAIN_PATH = r'/mnt/c/Users/cgarin/Documents/EDNiX'
sys.path.append('/mnt/c/Users/cgarin/PycharmProjects/EDNiX')
import Tools.Load_subject_with_BIDS
bids_dir = Tools.Load_subject_with_BIDS.linux_path(opj(r"C:\Users\cgarin\Desktop\BIDS_k9"))
output_path = Tools.Load_subject_with_BIDS.linux_path(opj(r"C:\Users\cgarin\Desktop\BIDS_k9\QC_report"))


# Create output directories
os.makedirs(output_path, exist_ok=True)

# Function to load and merge QC data
def load_qc_data(bids_dir, output_path):
    # Get all functional QC files
    QC_files_matrix = glob.glob(f"{bids_dir}/**/**/func/01_prepro/03_atlas_space/10_Results/fMRI_QC_matrix/**QC_result.txt", recursive=False)
    QC_files_SNR = glob.glob(f"{bids_dir}/**/**/func/01_prepro/01_funcspace/10_Results/fMRI_QC_SNR/**QC_result.txt", recursive=False)

    data = []

    # Process matrix files
    for file in QC_files_matrix:
        path_parts = file.split('/')
        species = path_parts[6]
        bids_dataset = path_parts[7]
        subject = path_parts[8]
        session = path_parts[9] if len(path_parts) > 9 else None

        # Extract subject and session from filename (e.g., sub-01_ses-1_task-rest_boldQC_result.txt)
        filename = Path(file).name
        parts = filename.split('_')
        sub = parts[0]
        ses = parts[1] if len(parts) > 1 else None

        metrics = {
            'species': species,
            'bids_dataset': bids_dataset,
            'subject': subject,
            'session': session,
            'filename': filename,
            'sub': sub,
            'ses': ses
        }

        with open(file, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()

                    if value.startswith('[') and value.endswith(']'):
                        try:
                            value = ast.literal_eval(value)[0]
                        except (ValueError, SyntaxError):
                            value = None

                    if value == 'nan' or value == '':
                        metrics[key] = None
                    else:
                        try:
                            metrics[key] = float(value)
                        except ValueError:
                            metrics[key] = value

                if "specific_roi" in line:
                    metrics["specific_roi"] = lines[i + 1].strip()
                elif "unspecific_ROI" in line:
                    metrics["unspecific_ROI"] = lines[i + 1].strip()
                elif "Result" in line:
                    metrics["Result"] = lines[i + 1].strip()

        data.append(metrics)

    # Process SNR files
    for file in QC_files_SNR:
        path_parts = file.split('/')
        species = path_parts[6]
        bids_dataset = path_parts[7]
        subject = path_parts[8]
        session = path_parts[9] if len(path_parts) > 9 else None

        filename = Path(file).name
        parts = filename.split('_')
        sub = parts[0]
        ses = parts[1] if len(parts) > 1 else None

        # Find matching entry in data or create new one
        existing_entry = next((item for item in data if item['filename'] == filename), None)
        if existing_entry:
            metrics = existing_entry
        else:
            metrics = {
                'species': species,
                'bids_dataset': bids_dataset,
                'subject': subject,
                'session': session,
                'filename': filename,
                'sub': sub,
                'ses': ses
            }
            data.append(metrics)

        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()

                    if value.startswith('[') and value.endswith(']'):
                        try:
                            value = ast.literal_eval(value)[0]
                        except (ValueError, SyntaxError):
                            value = None

                    if value == 'nan' or value == '':
                        metrics[key] = None
                    else:
                        try:
                            metrics[key] = float(value)
                        except ValueError:
                            metrics[key] = value

    df = pd.DataFrame(data)

    # Save raw combined data
    raw_data_path = os.path.join(output_path, "QC_raw_combined_data.csv")
    df.to_csv(raw_data_path, index=False)

    return df


# Function to identify outliers using IQR method
def identify_outliers(series):
    Q1 = series.quantile(0.25)
    Q3 = series.quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return (series < lower_bound) | (series > upper_bound)


# Function to create box plots for each variable
def create_boxplots(df, output_path):
    plot_dir = Path(output_path) / "boxplots"
    plot_dir.mkdir(exist_ok=True)

    outlier_report = {}

    # Get numeric columns only (excluding metadata columns)
    exclude_cols = ['species', 'bids_dataset', 'subject', 'session', 'filename', 'sub', 'ses', 'Result']
    numeric_cols = [col for col in df.select_dtypes(include=['float64', 'int64']).columns
                    if col not in exclude_cols]

    # Set style for better looking plots
    sns.set(style="whitegrid", palette="pastel")

    for col in numeric_cols:
        plt.figure(figsize=(14, 8))
        ax = sns.boxplot(x=df[col], width=0.3, fliersize=0)  # fliersize=0 hides default outlier markers

        # Add swarmplot to show all data points
        sns.swarmplot(x=df[col], color=".25", size=6, ax=ax)

        # Identify outliers
        outliers = identify_outliers(df[col])
        outlier_data = df.loc[outliers, ['sub', 'ses', col]]

        # Add to outlier report
        outlier_report[col] = outlier_data.values.tolist()

        # Annotate only the outliers with their labels
        for i, (sub, ses, value) in enumerate(zip(df['sub'], df['ses'], df[col])):
            if outliers[i]:
                label = f"{sub}_{ses}" if pd.notna(ses) else sub
                # Position the text slightly above the point
                ax.annotate(label,
                            xy=(0, value),
                            xytext=(0, 10),
                            textcoords='offset points',
                            ha='center',
                            va='bottom',
                            fontsize=9,
                            bbox=dict(boxstyle='round,pad=0.3',
                                      facecolor='white',
                                      alpha=0.8,
                                      edgecolor='gray'))

        # Add mean line
        mean_val = df[col].mean()
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=1, label=f'Mean: {mean_val:.2f}')

        # Customize plot appearance
        plt.title(f'Distribution of {col}\n(Outliers labeled)', fontsize=14, pad=20)
        plt.xlabel(col, fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)

        # Add legend for mean line
        plt.legend(loc='upper right')

        # Adjust layout to prevent label clipping
        plt.tight_layout()

        plot_path = plot_dir / f"{col}_boxplot.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

    return outlier_report

# Function to generate outlier report
def generate_outlier_report(outlier_report, output_path):
    # Create a list of all outlier occurrences
    outlier_list = []
    for metric, occurrences in outlier_report.items():
        for occ in occurrences:
            sub = occ[0]
            ses = occ[1] if len(occ) > 1 else None
            value = occ[2]
            outlier_list.append({
                'subject': sub,
                'session': ses,
                'metric': metric,
                'value': value
            })

    # Convert to DataFrame
    outlier_df = pd.DataFrame(outlier_list)

    # Count outliers per subject/session
    if 'session' in outlier_df.columns:
        count_df = outlier_df.groupby(['subject', 'session']).size().reset_index(name='outlier_count')
        count_df = count_df.sort_values('outlier_count', ascending=False)
    else:
        count_df = outlier_df.groupby('subject').size().reset_index(name='outlier_count')
        count_df = count_df.sort_values('outlier_count', ascending=False)

    # Save reports
    report_dir = Path(output_path) / "reports"
    report_dir.mkdir(exist_ok=True)

    detailed_path = report_dir / "outlier_detailed_report.csv"
    summary_path = report_dir / "outlier_summary_report.csv"

    outlier_df.to_csv(detailed_path, index=False)
    count_df.to_csv(summary_path, index=False)

    # Print summary
    print("\n=== Outlier Report Summary ===")
    print(f"Total outlier occurrences: {len(outlier_df)}")
    print(f"Total unique subject/session combinations with outliers: {len(count_df)}")
    print("\nTop 5 subject/sessions with most outliers:")
    print(count_df.head().to_string(index=False))

    return outlier_df, count_df


# Main function
def main():
    # Load and combine QC data
    print("Loading QC data...")
    df = load_qc_data(bids_dir, output_path)

    # Create boxplots and identify outliers
    print("Creating boxplots and identifying outliers...")
    outlier_report = create_boxplots(df, output_path)

    # Generate and save reports
    print("Generating outlier reports...")
    detailed_df, summary_df = generate_outlier_report(outlier_report, output_path)

    print(f"\nProcessing complete. Results saved to: {output_path}")


if __name__ == "__main__":
    main()