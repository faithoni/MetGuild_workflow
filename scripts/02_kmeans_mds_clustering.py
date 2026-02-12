import os
import json
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.manifold import MDS
from sklearn.cluster import KMeans

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import plotly.express as px

# -----------------------------
# Auto-load manifest from Script 1
# -----------------------------
manifest_path = os.path.join("manifests", "optimal_k_values.json")
if not os.path.exists(manifest_path):
    raise FileNotFoundError(
        f"Manifest not found at {manifest_path}. "
        "Run Script 1 first (it writes this file automatically)."
    )

with open(manifest_path, "r", encoding="utf-8") as f:
    m = json.load(f)

k_values = m["unique_k_values"]
file_path = m["cleaned_output_file"]

print(f"Loaded manifest: {manifest_path}")
print(f"Cleaned CSV: {file_path}")
print(f"Will run clustering for k values: {k_values}")

if not os.path.exists(file_path):
    raise FileNotFoundError(f"Cleaned CSV not found: {file_path}")

# -----------------------------
# Load data
# -----------------------------
data = pd.read_csv(file_path, encoding="latin1", index_col=0)

genome_ids = data.index
gene_data = data.values

# -----------------------------
# Distance + MDS (compute once, reuse for each k)
# -----------------------------
distance_matrix = pairwise_distances(gene_data, metric="jaccard")

mds = MDS(n_components=3, dissimilarity="precomputed", random_state=42)
mds_transformed = mds.fit_transform(distance_matrix)

mds_df_base = pd.DataFrame(mds_transformed, columns=["MDS1", "MDS2", "MDS3"])
mds_df_base["Genome_ID"] = genome_ids

base_name = os.path.splitext(os.path.basename(file_path))[0]

output_base_dir = "03_04_clustering/01_kmeans"
os.makedirs(output_base_dir, exist_ok=True)

# -----------------------------
# Loop over K values
# -----------------------------
for num_clusters in k_values:
    print(f"\n=== Running KMeans for k={num_clusters} ===")

    kmeans = KMeans(n_clusters=num_clusters, n_init=25, random_state=42)
    clusters = kmeans.fit_predict(mds_transformed) + 1  # start at 1

    # MDS df for plotting
    mds_df = mds_df_base.copy()
    mds_df["Cluster"] = clusters

    # Combine gene data with cluster assignments
    combined_data = data.copy()
    combined_data["Cluster"] = clusters
    combined_data = combined_data.reset_index()
    combined_data.rename(columns={combined_data.columns[0]: "Genomes"}, inplace=True)

    # Reorder columns
    gene_cols = [c for c in combined_data.columns if c not in ["Genomes", "Cluster"]]
    combined_data = combined_data[["Genomes", "Cluster"] + gene_cols]

    # Output directory for this k
    output_directory = os.path.join(output_base_dir, f"{num_clusters}_clusters")
    os.makedirs(output_directory, exist_ok=True)

    # Save combined CSV
    combined_output_file = os.path.join(output_directory, f"03_{base_name}_{num_clusters}_KM.csv")
    combined_data.to_csv(combined_output_file, index=False)
    print(f"Combined gene + cluster assignments saved to: {combined_output_file}")

    # Colors
    viridis = cm.get_cmap("viridis", num_clusters)
    if num_clusters == 1:
        color_map_hex = [mcolors.to_hex(viridis(0))]
    else:
        color_map_hex = [mcolors.to_hex(viridis(i / (num_clusters - 1))) for i in range(num_clusters)]

    # 2D plot (save, don’t block loop)
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(mds_df["MDS1"], mds_df["MDS2"], c=mds_df["Cluster"], cmap="viridis", s=50)
    plt.title(f"KMeans Clustering with {num_clusters} Clusters")
    plt.xlabel("MDS Dimension 1")
    plt.ylabel("MDS Dimension 2")
    plt.colorbar(scatter, label="Cluster")

    plot_output_file = os.path.join(output_directory, f"03_{base_name}_{num_clusters}_KM_2D.png")
    plt.savefig(plot_output_file)
    plt.close()
    print(f"2D plot saved to: {plot_output_file}")

    # 3D plot (Plotly)
    fig = px.scatter_3d(
        mds_df, x="MDS1", y="MDS2", z="MDS3",
        color="Cluster",
        color_discrete_sequence=color_map_hex,
        title=f"KMeans Clustering with {num_clusters} Clusters",
        labels={"Cluster": "Cluster"}
    )

    html_output_file = os.path.join(output_directory, f"03_{base_name}_{num_clusters}_KM.html")
    fig.write_html(html_output_file)
    print(f"Interactive 3D plot saved to: {html_output_file}")

print("\nAll clustering runs completed.")
