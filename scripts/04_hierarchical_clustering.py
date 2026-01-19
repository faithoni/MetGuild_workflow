import pandas as pd
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import fcluster
import matplotlib.pyplot as plt
from sklearn.metrics import jaccard_score
import seaborn as sns
import os

# Load the data
file_path = '01_DSK_DB_cleaned.csv'
data = pd.read_csv(file_path, header=0, index_col=0)

# Check if there are enough rows for clustering
if len(data) > 1:

    # Prompt user to enter the desired number of clusters
    try:
        k_optimal = int(input("Enter the desired number of clusters: "))
    except ValueError:
        print("Please enter a valid integer for the number of clusters.")
        exit()

    # Compute the Jaccard distance matrix (assuming binary data)
    def jaccard_distance(x, y):
        return 1 - jaccard_score(x, y, average='binary')

    # Calculate the distance matrix using pdist with the custom Jaccard distance
    dist_matrix = pdist(data.values, metric=jaccard_distance)
    
    # Perform hierarchical clustering using the average linkage method
    hclust_result = sch.linkage(dist_matrix, method='average')
    
    # Assign clusters based on the linkage matrix and the user-defined number of clusters
    clusters = fcluster(hclust_result, k_optimal, criterion='maxclust')
    
    # Generate a color palette with enough colors
    colors = sns.color_palette("Set3", n_colors=min(12, k_optimal)) + sns.color_palette("hsv", n_colors=max(k_optimal - 12, 0))

    # Plot the dendrogram and assign colors to clusters
    plt.figure(figsize=(10, 7))
    dend = sch.dendrogram(hclust_result, labels=data.index, color_threshold=0, above_threshold_color='black')

    # Add the cluster assignments as the first column in the dataset
    data.insert(0, 'Cluster', clusters)

    # Set the output directory based on the number of clusters (k_optimal)
    output_base_dir = "03_04_clustering/02_hierarchical"
    output_directory = os.path.join(output_base_dir, f"{k_optimal}_clusters")

    # Create directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    # Save the dataset with cluster assignments to a CSV file
    output_file = os.path.join(output_directory, f"04_{os.path.basename(file_path).split('.')[0]}_{k_optimal}_HC.csv")
    data.to_csv(output_file)
    print(f"Clustered data saved to: {output_file}")

    # Color labels according to cluster assignment
    label_colors = {str(data.index[i]): colors[clusters[i]-1] for i in range(len(clusters))}
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(label_colors[lbl.get_text()])    

    plt.title('Dendrogram of Hierarchical Clustering')
    plt.xlabel('Samples')

    # Save the plot image before displaying it
    plot_output_file = os.path.join(output_directory, f"04_{os.path.basename(file_path).split('.')[0]}_{k_optimal}_HC1.png")
    plt.savefig(plot_output_file)
    print(f"Plot saved to: {plot_output_file}")

    # Now show the plot
    plt.show()

else:
    print("Not enough data for clustering.")
