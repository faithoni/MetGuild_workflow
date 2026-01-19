# Load necessary libraries
if (!require(Hmisc)) install.packages('Hmisc', dependencies=TRUE)
if (!require(reshape2)) install.packages('reshape2', dependencies=TRUE)
library(Hmisc)
library(reshape2)

# Read the dataset from CSV file
data <- read.csv("05_stat_analysis/02_hierarchical/14_clusters/05_DSK_DB_14_HC_chi_stat_sig.csv", header=TRUE)

# Remove the first column which contains species/genomes names
data <- data[,-1]

# Exclude the second column from the analysis
data <- data[,-1]  # Exclude second column

# Calculate the Spearman correlation matrix
spearman_matrix <- rcorr(as.matrix(data), type="spearman")

# Extract the correlation coefficients and p-values
spearman_corr <- spearman_matrix$r  # Correlation coefficients
spearman_pvals <- spearman_matrix$P # p-values

# Flatten the correlation matrix and p-value matrix
spearman_corr_flat <- melt(spearman_corr)
spearman_pvals_flat <- melt(spearman_pvals)

# Combine the correlation coefficients and p-values into one data frame
flat_matrix <- data.frame(
  Var1 = spearman_corr_flat$Var1,
  Var2 = spearman_corr_flat$Var2,
  Correlation = spearman_corr_flat$value,
  P_value = spearman_pvals_flat$value
)

# Remove rows where Var1 equals Var2 (diagonal elements)
flat_matrix <- flat_matrix[flat_matrix$Var1 != flat_matrix$Var2,]

flat_matrix <- flat_matrix

# Define output paths
output_corr_matrix <- "/Users/faith/Documents/Oni_PhD/01_data_analyses_per_chapter/01_chapter_1_redone/01_DSK_DB/06_correlation/01_featurevsfeature_correlation/02_hierarchical/14_clusters/05_DSK_DB_14_HC_corr_matrix.csv"
output_flattened_matrix <- "/Users/faith/Documents/Oni_PhD/01_data_analyses_per_chapter/01_chapter_1_redone/01_DSK_DB/06_correlation/01_featurevsfeature_correlation/02_hierarchical/14_clusters/05_DSK_DB_14_HC_flattened_corr_pvals.csv"
output_filtered_matrix <- "/Users/faith/Documents/Oni_PhD/01_data_analyses_per_chapter/01_chapter_1_redone/01_DSK_DB/06_correlation/01_featurevsfeature_correlation/02_hierarchical/14_clusters/05_DSK_DB_14_HC_filtered_corr_pvals.csv"

# Create directories if they don't exist
output_dir <- dirname(output_corr_matrix)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the correlation matrix
write.csv(spearman_corr, output_corr_matrix, row.names = FALSE)

# Save the flattened matrix
write.csv(flat_matrix, output_flattened_matrix, row.names = FALSE)

# Filter the observations: Correlation >= 0.5 and P-value < 0.05
filtered_matrix <- flat_matrix[flat_matrix$Correlation >= 0.5 | flat_matrix$Correlation <= -0.5 & flat_matrix$P_value < 0.05,]

# Save the filtered matrix to a new CSV file, excluding row numbers
write.csv(filtered_matrix, output_filtered_matrix, row.names = FALSE)