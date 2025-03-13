import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
df = pd.read_csv("/scratch/cgarin/QC_func_matrix/QC_results_summary.csv")
output_path =  '/scratch/cgarin/QC_func_matrix/'
# Filter out rows where 'Result' is None
filtered_df = df.dropna(subset=['Result'])

# Group by 'species' and 'Result' to count occurrences
result_counts = filtered_df.groupby(['species', 'Result']).size().unstack(fill_value=0)

# Ensure the data is in integer format
result_counts = result_counts.astype(int)

# Plotting with matplotlib
fig, ax = plt.subplots(figsize=(12, 8))

# Stacked bar plot
result_counts.plot(kind='bar', stacked=True, ax=ax, color=sns.color_palette("Set2", len(result_counts.columns)))

# Customize plot
ax.set_title("Frequency of 'Result' Categories Across Species")
ax.set_xlabel("Species")
ax.set_ylabel("Frequency")
ax.legend(title="Result")
plt.tight_layout()

# Save the Result frequency plot
result_plot_path = os.path.join(output_path, "Result_Frequency_QC_plot.png")
plt.savefig(result_plot_path)
plt.close()

print(f"Result frequency plot saved in {output_path}")




# Convert the collected data into a DataFrame
df = pd.read_csv("/scratch/cgarin/QC_func_matrix/QC_results_summary.csv")
# Filter out 'BIDS_BENHamed' BIDS dataset
df = df[df['bids_dataset'] != 'BIDS_BENHamed']

# Define the metrics of interest (for now, use example metrics from the func QC file)
#metrics_of_interest = ['species', 'bids_dataset', 'average_snr_Gray_Matter', 'avg_motion_velocity', 'W-statistic', 'Clustering Stability']
#df = df[metrics_of_interest].dropna()  # Drop rows with NaN values

# Plotting example metrics for functional QC
sns.set(style='whitegrid')
plt.figure(figsize=(12, 8))

# Plot each metric (this is an example for one metric, you can repeat this for others)
for metric in df[2:]:
    plt.figure(figsize=(12, 8))
    sns.boxplot(data=df, x='species', y=metric, hue='bids_dataset', palette='Set2')
    plt.title(f"Boxplot of {metric} Across Species and BIDS Datasets")
    plt.xlabel('Species')
    plt.ylabel(f'{metric}')
    plt.legend(title='BIDS Dataset')
    plt.savefig(f"{output_path}{metric}_boxplot.png")
    plt.close()




# Example: Show a plot for average SNR in gray matter
sns.boxplot(data=df, x='species', y='average_snr_Gray_Matter', hue='bids_dataset', palette='Set2')
plt.title("Boxplot of Average SNR in Gray Matter Across Species and BIDS Datasets")
plt.xlabel('Species')
plt.ylabel('Average SNR (Gray Matter)')
plt.legend(title='BIDS Dataset')
plt.show()