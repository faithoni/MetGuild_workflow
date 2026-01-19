import pandas as pd
import numpy as np
from scipy import stats
import os

# Set directory for data
os.makedirs("06_correlation", exist_ok=True)

# Function to read data
def read_data(file_path):
    if os.path.exists(file_path):
        return pd.read_csv(file_path, index_col=0)
    else:
        raise FileNotFoundError(f"File does not exist: {file_path}")

# Function to calculate correlation p-values
def cor_prob(X):
    corr = X.corr(method='spearman')
    p_values = pd.DataFrame(np.zeros(corr.shape), columns=corr.columns, index=corr.index)

    for i in range(len(corr.columns)):
        for j in range(i):
            # Calculate p-value using Spearman correlation
            _, p_value = stats.spearmanr(X[corr.columns[i]], X[corr.columns[j]])
            p_values.iloc[j, i] = p_value
            p_values.iloc[i, j] = p_value  # Symmetric

    np.fill_diagonal(p_values.values, np.nan)  # Fill diagonal with NaN
    return p_values

# Load data
data_05 = read_data("05_stat_analysis/05_DSK_DB_7_HC_chi_stat_sig.csv")

# Calculate the correlation matrix
correlation_matrix_05 = data_05.corr(method='spearman')

# Write the correlation matrix to a file
correlation_matrix_05.to_csv("06_correlation/02_correlationmatrix_05.txt", sep="\t")

# Calculate the P-value matrix
Pval_matrix_05 = cor_prob(data_05)

# Write the P-value matrix to a file
Pval_matrix_05.to_csv("06_correlation/03_Pval_matrix_05.txt", sep="\t")

# Flatten the P-value matrix
flat_matrix_05 = Pval_matrix_05.stack().reset_index()
flat_matrix_05.columns = ['i', 'j', 'p']

# Remove NaNs
flat_matrix_05.dropna(inplace=True)

# Add correlation values to the flattened matrix
# Ensure the correlation matrix is stacked similarly
flat_matrix_05['cor'] = correlation_matrix_05.stack().values

# Check lengths before filtering
print(f"Length of flat_matrix_05: {len(flat_matrix_05)}")
print(f"Length of correlation values: {len(correlation_matrix_05.stack().values)}")

# Filter relevant correlations (absolute values above 0.8 and p-values lower than 0.05)
select_cor_05 = flat_matrix_05[abs(flat_matrix_05['cor']) > 0.8]
select_cor_pval_05 = select_cor_05[select_cor_05['p'] < 0.05]

# Sort by correlation
select_cor_pval_ordered_05 = select_cor_pval_05.sort_values(by='cor', ascending=False)

# Write the flattened matrix to a file
flat_matrix_05.to_csv("06_correlation/04_Correlations_05.txt", sep="\t", index=False)

# Write selected correlations to a file
select_cor_pval_ordered_05.to_csv("06_correlation/05_Selected_Correlations_05_AllSamples.txt", sep="\t", index=False)

# Output the head of the selected correlations
print(select_cor_pval_ordered_05.head())
