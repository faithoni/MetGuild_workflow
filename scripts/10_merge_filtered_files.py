import pandas as pd
import os

# Define the paths to the three CSV files
folder_paths = [
    "06_correlation/01_featurevsfeature_correlation/02_hierarchical/14_clusters/05_DSK_DB_14_HC_filtered_corr_pvals.csv",
    "06_correlation/02_featurevscluster_correlation/02_hierarchical/14_clusters/05_DSK_DB_14_HC_filtered_corr_pvalues_matrix.csv",
    "06_correlation/03_clustervscluster_correlation/02_hierarchical/14_clusters/05_DSK_DB_14_HC_filtered_results.csv"
]

# Read all CSV files into a list of DataFrames
dataframes = [pd.read_csv(path) for path in folder_paths]

# Concatenate the DataFrames and keep only one header
combined_df = pd.concat(dataframes, ignore_index=True)

# Function to remove duplicate rows based on matching columns
def remove_duplicate_rows(df):
    # Create a set to track unique pairs (column1, column2)
    unique_pairs = set()
    rows_to_keep = []

    for _, row in df.iterrows():
        # Create pairs from the first two columns
        pair_1 = (row.iloc[0], row.iloc[1])
        pair_2 = (row.iloc[1], row.iloc[0])

        # If neither pair is in the set, keep the row and add pair_1 to the set
        if pair_1 not in unique_pairs and pair_2 not in unique_pairs:
            rows_to_keep.append(row)
            unique_pairs.add(pair_1)

    # Convert the list of rows to a DataFrame
    return pd.DataFrame(rows_to_keep, columns=df.columns)

# Apply the function to remove duplicate rows
cleaned_df = remove_duplicate_rows(combined_df)

# Extract the prefix from the input file paths
prefix = "05_DSK_DB_14_HC"

# Define output directory and file name
output_dir = "07_network_construction/edge_lists"
output_file_name = f"{prefix}_edge_list.csv"
output_file_path = os.path.join(output_dir, output_file_name)

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Save the cleaned DataFrame to a new CSV file
cleaned_df.to_csv(output_file_path, index=False)

print("Merging complete, duplicates removed, and file saved as:", output_file_path)