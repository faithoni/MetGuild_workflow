# Load necessary libraries
library(dplyr)
library(tidyr)

# Load the dataset
df <- read.csv("/Users/faith/Documents/Oni_PhD/01_data_analyses_per_chapter/01_chapter_1_redone/01_DSK_DB/05_stat_analysis/01_kmeans/3_clusters/05_DSK_DB_3_KM_chi_stat_sig.csv", header = TRUE)

# Create a binary presence-absence matrix for proteins by Cluster
binary_matrix <- df %>%
  pivot_longer(cols = -c(Genomes, Cluster), names_to = "Protein", values_to = "Presence") %>%
  group_by(Cluster, Protein) %>%
  summarise(Presence = max(Presence, na.rm = TRUE)) %>%
  pivot_wider(names_from = Cluster, values_from = Presence, values_fill = 0)  # Fill NA with 0

# Calculate the correlation matrix between Cluster groups
correlation_matrix <- cor(binary_matrix[-1])  # Exclude the first column (Protein) for correlation

# Calculate p-value matrix
p_value_matrix <- matrix(NA, nrow = ncol(correlation_matrix), ncol = ncol(correlation_matrix))
rownames(p_value_matrix) <- colnames(correlation_matrix)
colnames(p_value_matrix) <- colnames(correlation_matrix)

for (i in 1:(ncol(binary_matrix) - 1)) {
  for (j in 1:(ncol(binary_matrix) - 1)) {
    if (i != j) {
      test <- cor.test(binary_matrix[[i + 1]], binary_matrix[[j + 1]], method = "spearman")
      p_value_matrix[i, j] <- test$p.value
    } else {
      p_value_matrix[i, j] <- NA  # No p-value for correlation with itself
    }
  }
}

# Print the correlation matrix and p-value matrix
print(correlation_matrix)
print(p_value_matrix)

# Define output directories
output_dir <- "/Users/faith/Documents/Oni_PhD/01_data_analyses_per_chapter/01_chapter_1_redone/01_DSK_DB/06_correlation/03_clustervscluster_correlation/01_kmeans/3_clusters"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Write the correlation matrix to a CSV file
write.csv(correlation_matrix, file.path(output_dir, "05_DSK_DB_3_KM_corr_matrix.csv"), row.names = TRUE)

# Write the p-value matrix to a CSV file
write.csv(p_value_matrix, file.path(output_dir, "05_DSK_DB_3_KM_p_value_matrix.csv"), row.names = TRUE)

# Flatten the correlation matrix and p-value matrix to long format
correlation_long <- as.data.frame(as.table(correlation_matrix))
colnames(correlation_long) <- c("Var1", "Var2", "Correlation")

p_value_long <- as.data.frame(as.table(p_value_matrix))
colnames(p_value_long) <- c("Var1", "Var2", "P_value")

# Merge the two long format data frames on Var1 and Var2
merged_long <- correlation_long %>%
  left_join(p_value_long, by = c("Var1", "Var2"))


# Add 'Cluster_' prefix to Var1 and Var2 without a space
merged_long <- merged_long %>%
  mutate(Var1 = paste0("Cluster_", Var1),
         Var2 = paste0("Cluster_", Var2))


# Write the flattened data to a new CSV file
write.csv(merged_long, file.path(output_dir, "05_DSK_DB_3_KM_flattened_results.csv"), row.names = FALSE)

# Filter for correlations greater than 0.5 or less than -0.5 and p-values less than or equal to 0.05
filtered_results <- merged_long %>%
  filter((Correlation >= 0.5 | Correlation <= -0.5) & P_value < 0.05)

# Write the filtered results to a new CSV file
write.csv(filtered_results, file.path(output_dir, "05_DSK_DB_3_KM_filtered_results.csv"), row.names = FALSE)

# Print counts of correlations
total_vars <- nrow(merged_long)
significant_vars <- nrow(filtered_results)
non_significant_vars <- total_vars - significant_vars

cat("Total number of variables:", total_vars, "\n")
cat("Number of variables that satisfy the conditions (Correlation >= 0.5 or Correlation <= -0.5 & P_value < 0.05):", significant_vars, "\n")
cat("Number of variables that do not satisfy the conditions:", non_significant_vars, "\n")
