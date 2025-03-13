import matplotlib.pyplot as plt

# Filter data for Rat
rat_data = df[df["species"] == "Rat"]

# Count occurrences of each category in the "Result" column
result_counts = rat_data["Result"].value_counts(normalize=True) * 100  # Convert to percentage

# Plot bar chart
plt.figure(figsize=(8, 5))
result_counts.plot(kind="bar", color=["blue", "green", "red", "orange"])
plt.title("Percentage Distribution of Result Categories for Rat")
plt.ylabel("Percentage")
plt.xlabel("Result Category")
plt.xticks(rotation=45)
plt.grid(axis="y", linestyle="--", alpha=0.7)

# Show plot
plt.show()


import pandas as pd

# Load the dataset
file_path = "/scratch/cgarin/QC_func_matrix/QC_results_summary.csv"  # Update with the correct file path
df = pd.read_csv(file_path)


def classify_fc(row):
    specific = row["specific_roi"]  # rS1bf_left_to_right
    unspecific = row["unspecific_ROI"]  # rS1bf_left_to_ACA

    if specific > 0.1 and unspecific < 0.1:
        return "Specific"
    elif specific < 0.1 and unspecific > 0.1:
        return "Unspecific"
    elif abs(specific) < 0.1 and abs(unspecific) < 0.1:
        return "No"
    else:
        return "Spurious"

# Apply function to create a new column
df["FC_Category"] = df.apply(classify_fc, axis=1)

import matplotlib.pyplot as plt

# Filter only Rat data
rat_data = df[df["species"] == "Rat"]

# Count occurrences of each category
category_counts = rat_data["FC_Category"].value_counts()
print(category_counts)
# Create a single stacked bar with raw counts
plt.figure(figsize=(4, 6))
plt.bar(["Rat"], [category_counts.get("Non-Specific", 0)], color="blue", label="Non-Specific")
plt.bar(["Rat"], [category_counts.get("Specific", 0)], bottom=[category_counts.get("Non-Specific", 0)], color="green", label="Specific")
plt.bar(["Rat"], [category_counts.get("No", 0)], bottom=[category_counts.get("Non-Specific", 0) + category_counts.get("Specific", 0)], color="red", label="No")
plt.bar(["Rat"], [category_counts.get("Spurious", 0)], bottom=[category_counts.get("Non-Specific", 0) + category_counts.get("Specific", 0) + category_counts.get("No", 0)], color="orange", label="Spurious")

# Labels and legend
plt.title("FC Category Distribution for Rat (Raw Counts)")
plt.ylabel("Number of Subjects")
plt.legend()

# Show plot
plt.show()
