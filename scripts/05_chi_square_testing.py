import os
import pandas as pd
from scipy.stats import chi2_contingency

# File paths - edit as necessary
input_csv = "03_04_clustering/02_hierarchical/14_clusters/04_01_DSK_DB_cleaned_14_HC.csv"
output_chi2_csv = "05_stat_analysis/02_hierarchical/14_clusters/05_DSK_DB_14_HC_chi_squared_results.csv"
output_significant_csv = "05_stat_analysis/02_hierarchical/14_clusters/05_DSK_DB_14_HC_chi_stat_sig.csv"
output_non_significant_csv = "05_stat_analysis/02_hierarchical/14_clusters/05_DSK_DB_14_HC_chi_non_sig.csv"

# Ensure the output directories exist
def ensure_dir(file_path):
    """Ensures that the directory for the file_path exists, creating it if necessary."""
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

# Ensure all output directories exist
ensure_dir(output_chi2_csv)
ensure_dir(output_significant_csv)
ensure_dir(output_non_significant_csv)

# Read the CSV file
data = pd.read_csv(input_csv)

# Check if the column "Cluster" exists and convert it to a category
if 'Cluster' in data.columns:
    data['Cluster'] = data['Cluster'].astype('category')
else:
    raise ValueError("Column 'Cluster' not found in the data.")

# Check for missing values in the Cluster column
if data['Cluster'].isna().any():
    raise ValueError("Missing values found in 'Cluster' column.")

# Function to perform chi-squared test for each gene/pathway/module
def perform_chi_squared_test(df, gene_column):
    # Check for missing values in the gene column
    if df[gene_column].isna().any():
        raise ValueError(f"Missing values found in gene column: {gene_column}")
    
    # Ensure gene column contains only 0s and 1s
    if not df[gene_column].isin([0, 1]).all():
        raise ValueError(f"Non-binary values found in gene column: {gene_column}")
    
    # Create a contingency table
    contingency_table = pd.crosstab(df['Cluster'], df[gene_column])
    
    # Perform the chi-squared test
    chi2_stat, p_value, _, _ = chi2_contingency(contingency_table)
    
    return p_value

# Get the list of genes (assuming they start from the 3rd column onwards)
gene_columns = data.columns[2:]

# Apply the chi-squared test for each gene and store the results
results = []
for gene in gene_columns:
    try:
        p_value = perform_chi_squared_test(data, gene)
        p_value_rounded = round(p_value, 2)
        gene_name = f"*{gene}" if p_value < 0.05 else gene
        results.append({'gene': gene_name, 'p_value': p_value_rounded})
    except ValueError as e:
        print(e)

# Convert the results to a DataFrame
results_df = pd.DataFrame(results)

# Save the results to a CSV file
results_df.to_csv(output_chi2_csv, index=False)

# Identify significant genes (p-value < 0.05)
significant_genes = results_df['p_value'] < 0.05
num_significant = significant_genes.sum()
num_non_significant = len(significant_genes) - num_significant

print(f"Number of significant genes: {num_significant}")
print(f"Number of non-significant genes: {num_non_significant}")

# Filter the original data to exclude non-significant genes
non_significant_gene_columns = results_df[~significant_genes]['gene'].str.lstrip('*')
filtered_data = data.drop(columns=non_significant_gene_columns)

# Save the data of statistically significant genes to a new file
filtered_data.to_csv(output_significant_csv, index=False)

# Save non-significant genes to a file
pd.DataFrame({'gene': non_significant_gene_columns.tolist()}).to_csv(output_non_significant_csv, index=False)
