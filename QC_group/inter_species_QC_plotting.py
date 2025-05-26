import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import ast

# Define base path and output path
base_path = '/home/cgarin/Documents/EDNiX_study/MRI'
output_path = '/home/cgarin/Documents/EDNiX_study/MRI/QC_func_inter-species'

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# Get all functional QC files
QC_files = glob.glob(f"{base_path}/**/**/**/**/func/01_prepro/03_atlas_space/10_Results/fMRI_QC_matrix/**QC_result.txt", recursive=False)

data = []
# Traverse through each QC file
for file in QC_files:
    # Extract species, BIDS dataset, and subject based on the path structure
    path_parts = file.split('/')
    species = path_parts[6]
    bids_dataset = path_parts[7]
    subject = path_parts[8]

    # Initialize a dictionary to store metrics for this file
    metrics = {'species': species, 'bids_dataset': bids_dataset, 'subject': subject}

    # Read each file and extract metrics
    with open(file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # Standard key-value pairs
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()

                # Handle lists in the form '[value]' by extracting the first item
                if value.startswith('[') and value.endswith(']'):
                    try:
                        value = ast.literal_eval(value)[0]
                    except (ValueError, SyntaxError):
                        value = None

                # Convert 'nan' values and try to cast to float
                if value == 'nan' or value == '':
                    metrics[key] = None
                else:
                    try:
                        metrics[key] = float(value)
                    except ValueError:
                        metrics[key] = value

            # Handle specific keys that are on separate lines
            if "specific_roi" in line:
                metrics["specific_roi"] = lines[i + 1].strip()
            elif "unspecific_ROI" in line:
                metrics["unspecific_ROI"] = lines[i + 1].strip()
            elif "Result" in line:
                metrics["Result"] = lines[i + 1].strip()

    # Append the collected metrics for this file
    data.append(metrics)

# Convert the data into a DataFrame
df = pd.DataFrame(data)

# Display the DataFrame to confirm extraction
print(df.head())

# Save DataFrame to a CSV for easy reference if needed
output_csv = os.path.join(output_path, "QC_results_summary.csv")
df.to_csv(output_csv, index=False)
print(f"Data saved to {output_csv}")



# Define base path and output path
base_path = '/home/cgarin/Documents/EDNiX_study/MRI'
output_path = '/home/cgarin/Documents/EDNiX_study/MRI/QC_func_inter-species'

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# Get all functional QC files
QC_files = glob.glob(f"{base_path}/**/**/**/**/func/01_prepro/03_atlas_space/10_Results/fMRI_QC_matrix/**QC_result.txt", recursive=False)

data = []
# Traverse through each QC file
for file in QC_files:
    # Extract species, BIDS dataset, and subject based on the path structure
    path_parts = file.split('/')
    species = path_parts[6]
    bids_dataset = path_parts[7]
    subject = path_parts[8]

    # Initialize a dictionary to store metrics for this file
    metrics = {'species': species, 'bids_dataset': bids_dataset, 'subject': subject}

    # Read each file and extract metrics
    with open(file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            # Standard key-value pairs
            if ':' in line:
                key, value = line.split(':', 1)
                key = key.strip()
                value = value.strip()

                # Handle lists in the form '[value]' by extracting the first item
                if value.startswith('[') and value.endswith(']'):
                    try:
                        value = ast.literal_eval(value)[0]
                    except (ValueError, SyntaxError):
                        value = None

                # Convert 'nan' values and try to cast to float
                if value == 'nan' or value == '':
                    metrics[key] = None
                else:
                    try:
                        metrics[key] = float(value)
                    except ValueError:
                        metrics[key] = value

            # Handle specific keys that are on separate lines
            if "specific_roi" in line:
                metrics["specific_roi"] = lines[i + 1].strip()
            elif "unspecific_ROI" in line:
                metrics["unspecific_ROI"] = lines[i + 1].strip()
            elif "Result" in line:
                metrics["Result"] = lines[i + 1].strip()

    # Append the collected metrics for this file
    data.append(metrics)

# Convert the data into a DataFrame
df = pd.DataFrame(data)

# Display the DataFrame to confirm extraction
print(df.head())

# Save DataFrame to a CSV for easy reference if needed
output_csv = os.path.join(output_path, "QC_results_summary.csv")
df.to_csv(output_csv, index=False)
print(f"Data saved to {output_csv}")
