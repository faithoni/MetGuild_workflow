# Load necessary libraries
library(tidyverse)
if (!require(Hmisc)) install.packages('Hmisc', dependencies=TRUE)
if (!require(reshape2)) install.packages('reshape2', dependencies=TRUE)
library(Hmisc)
library(reshape2)

# Read the data
data <- read.csv("/Users/faith/Documents/Oni_PhD/01_data_analyses_per_chapter/01_chapter_1_redone/01_DSK_DB/05_stat_analysis/01_kmeans/14_clusters/05_DSK_DB_14_KM_chi_stat_sig.csv", header=TRUE)

# Extract the clusters and pathways data
clusters <- data$Cluster
pathways_data <- data %>% select(-Genomes, -Cluster)

# Initialize data frames to store the correlation coefficients and p-values
correlation_results <- data.frame(Pathway = colnames(pathways_data))
pvalue_results <- data.frame(Pathway = colnames(pathways_data))

# Calculate the correlation and p-values for each cluster
for (cluster_num in unique(clusters)) {
  # Create a binary indicator for the current cluster
  cluster_indicator <- ifelse(clusters == cluster_num, 1, 0)
  
  # Initialize vectors to store correlations and p-values
  correlations <- numeric(ncol(pathways_data))
  pvalues <- numeric(ncol(pathways_data))
  
  # Calculate the Spearman correlation and p-values between the cluster indicator and each pathway
  for (i in 1:ncol(pathways_data)) {
    result <- rcorr(cluster_indicator, pathways_data[[i]], type = "spearman")
    correlations[i] <- result$r[1, 2]  # Correlation value
    pvalues[i] <- result$P[1, 2]       # p-value
  }
  
  # Add the correlations and p-values to the respective data frames
  correlation_results <- cbind(correlation_results, correlations)
  pvalue_results <- cbind(pvalue_results, pvalues)
}

# Rename the columns to indicate the cluster numbers
colnames(correlation_results)[-1] <- paste0("Cluster_", unique(clusters))
colnames(pvalue_results)[-1] <- paste0("Cluster_", unique(clusters))

# Flatten (melt) both the correlation and p-value matrices
correlation_long <- melt(correlation_results, id.vars = "Pathway", variable.name = "Cluster", value.name = "Correlation")
pvalue_long <- melt(pvalue_results, id.vars = "Pathway", variable.name = "Cluster", value.name = "P_value")

# Merge the two long-form data frames on Pathway and Cluster
merged_results <- merge(correlation_long, pvalue_long, by = c("Pathway", "Cluster"))

# Rename columns to Var1, Var2, Correlation, and P_value
colnames(merged_results) <- c("Pathway", "Cluster", "Correlation", "P_value")

# Reorder columns to match the desired output
merged_results <- merged_results %>%
  select(Var1 = Pathway, Var2 = Cluster, Correlation, P_value)

# Display the merged results
#print(merged_results)

# Define output paths
output_dir <- "/Users/faith/Documents/Oni_PhD/01_data_analyses_per_chapter/01_chapter_1_redone/01_DSK_DB/06_correlation/02_featurevscluster_correlation/01_kmeans/14_clusters"
output_file_full <- file.path(output_dir, "05_DSK_DB_14_KM_corr_pvalues_matrix.csv")
output_file_filtered <- file.path(output_dir, "05_DSK_DB_14_KM_filtered_corr_pvalues_matrix.csv")
output_file_corr <- file.path(output_dir, "05_DSK_DB_14_KM_corr_matrix.csv")
output_file_pvals <- file.path(output_dir, "05_DSK_DB_14_KM_pvalues_matrix.csv")

# Create directories if they don't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save the full merged results to a CSV file
write.csv(merged_results, file = output_file_full, row.names = FALSE)

# Filter the results based on correlation and p-value criteria
filtered_results <- merged_results %>%
  filter(Correlation >= 0.5 | Correlation <= -0.5, P_value < 0.05)

# Display the filtered results
#print(filtered_results)

# Save the filtered results to a new CSV file
write.csv(filtered_results, file = output_file_filtered, row.names = FALSE)

# Save the separate correlation and p-value matrices for reference
write.csv(correlation_results, file = output_file_corr, row.names = FALSE)
write.csv(pvalue_results, file = output_file_pvals, row.names = FALSE)
