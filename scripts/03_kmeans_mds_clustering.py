import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import plotly.express as px
import os

# Step 1: Import the data from CSV
file_path = "01_DSK_DB_cleaned.csv"  # Update this path if needed
data = pd.read_csv(file_path, encoding='latin1', index_col=0)

# Extract genome IDs and gene presence/absence data
genome_ids = data.index
gene_data = data.values

# Step 2: Calculate distance matrix using Jaccard distance for binary values
distance_matrix = pairwise_distances(gene_data, metric='jaccard')

# Step 3: Apply Classical MDS to reduce to 3 dimensions
mds = MDS(n_components=3, dissimilarity="precomputed", random_state=42)
mds_transformed = mds.fit_transform(distance_matrix)

# Step 4: Manual input for the number of clusters
try:
    num_clusters = int(input("Enter the number of clusters: "))
except ValueError:
    print("Invalid input. Using default value of 5.")
    num_clusters = 5

# Perform KMeans clustering with the specified number of clusters
kmeans = KMeans(n_clusters=num_clusters, n_init=25, random_state=42)
kmeans.fit(mds_transformed)
clusters = kmeans.labels_

# Adjust clusters to start from 1 instead of 0
clusters = clusters + 1

# Prepare DataFrame for visualization
mds_df = pd.DataFrame(mds_transformed, columns=['MDS1', 'MDS2', 'MDS3'])
mds_df['Cluster'] = clusters
mds_df['Genome_ID'] = genome_ids

# Create a DataFrame to combine gene data with cluster assignments
combined_data = pd.DataFrame(data.values, index=genome_ids, columns=data.columns)
combined_data['Cluster'] = clusters
combined_data = combined_data.reset_index()  # Reset index to include genome IDs as a column

# Rename the index column to 'Genomes'
combined_data.rename(columns={'index': 'Genomes'}, inplace=True)

# Reorder columns: 'Cluster' first, then the original gene dataset columns
column_order = ['Genomes', 'Cluster'] + [col for col in combined_data.columns if col not in ['Genomes', 'Cluster']]
combined_data = combined_data[column_order]

# Step 5: Set the output directory dynamically based on the number of clusters
output_base_dir = "03_04_clustering/01_k_means"
output_directory = os.path.join(output_base_dir, f"{num_clusters}_clusters")

# Create directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Save the combined dataset to a single CSV file in the dynamically created directory
base_name, _ = os.path.splitext(file_path)
combined_output_file = os.path.join(output_directory, f"03_{base_name}_{num_clusters}_KM.csv")
combined_data.to_csv(combined_output_file, index=False)
print(f"Combined gene and cluster assignments saved to: {combined_output_file}")

# Generate a consistent color map
viridis = cm.get_cmap('viridis', num_clusters)
color_map_hex = [mcolors.to_hex(viridis(i / (num_clusters - 1))) for i in range(num_clusters)]

# Step 6: Visualize the clusters in 2D and save the plot
plt.figure(figsize=(10, 8))
scatter = plt.scatter(mds_df['MDS1'], mds_df['MDS2'], c=mds_df['Cluster'], cmap='viridis', s=50)
plt.title(f'KMeans Clustering with {num_clusters} Clusters')
plt.xlabel('MDS Dimension 1')
plt.ylabel('MDS Dimension 2')
plt.colorbar(scatter, label='Cluster')

# Save the 2D plot before showing it
plot_output_file = os.path.join(output_directory, f"03_{base_name}_{num_clusters}_KM_2D.png")
plt.savefig(plot_output_file)
print(f"2D plot saved to: {plot_output_file}")

# Show the plot after saving it
plt.show()

# Step 7: Visualize and save the 3D plot using Plotly
fig = px.scatter_3d(mds_df, x='MDS1', y='MDS2', z='MDS3', color='Cluster', 
                    color_discrete_sequence=color_map_hex,
                    title=f'KMeans Clustering with {num_clusters} Clusters',
                    labels={'Cluster': 'Cluster'})

# Save the interactive 3D plot as an HTML file in the dynamically created directory
html_output_file = os.path.join(output_directory, f"03_{base_name}_{num_clusters}_KM.html")
fig.write_html(html_output_file)
print(f"Interactive 3D plot saved to: {html_output_file}")
