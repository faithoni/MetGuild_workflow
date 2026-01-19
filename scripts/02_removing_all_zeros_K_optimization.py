import pandas as pd
import os
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)

# Load the dataset
input_file = 'DSK_DB.csv'  # Replace with your actual CSV file path

try:
    # Read the dataset without setting the index for debugging
    df = pd.read_csv(input_file, encoding='latin1')
except Exception as e:
    logging.error(f"Error reading {input_file}: {e}")
    raise

# Print the columns to check for 'Genomes'
logging.info(f"Columns in the dataset: {df.columns.tolist()}")

# Check if 'Genomes' column exists
if 'Genomes' not in df.columns:
    logging.error("'Genomes' column not found in the dataset.")
else:
    # Store the "Genomes" column and drop it from the original dataframe for later use
    genomes_column = df['Genomes']
    df = df.drop(columns=['Genomes'])

    # Drop columns where all values are zero
    df_cleaned = df.loc[:, (df != 0).any(axis=0)]

    # Optionally, drop columns where all values are NaN
    df_cleaned = df_cleaned.dropna(axis=1, how='all')

    # Print shapes of original and cleaned DataFrames
    logging.info(f"Original shape: {df.shape}, Cleaned shape: {df_cleaned.shape}")

    # Add the "Genomes" column back to the cleaned DataFrame
    df_cleaned.insert(0, 'Genomes', genomes_column)

    # Generate the output file name by adding '_cleaned' suffix before the file extension
    base, ext = os.path.splitext(input_file)
    output_file = f"01_{base}_cleaned{ext}"

    # Save the cleaned dataset to a new CSV file
    try:
        df_cleaned.to_csv(output_file, index=False)
        logging.info(f"Cleaned dataset saved to {output_file}")
    except Exception as e:
        logging.error(f"Error saving {output_file}: {e}")


#script for k_optimization begins!!!!!!!!!!!!


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.preprocessing import StandardScaler

# Load data from CSV file
file_path = '01_DSK_DB_cleaned.csv'  # Update this to your CSV file path
data = pd.read_csv(file_path, index_col=0)  # Skip the first column for row indices

# Extract presence/absence matrix
X = data.values

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

def evaluate_clustering(X, max_k):
    """
    Evaluate clustering performance using different methods.
    
    Parameters:
    X (array-like): Input data for clustering.
    max_k (int): Maximum number of clusters to test.
    
    Returns:
    tuple: Lists containing scores for each method.
    """
    elbow_scores = []
    silhouette_scores = []
    calinski_harabasz_scores = []
    davies_bouldin_scores = []
    
    for k in range(2, max_k):
        kmeans = KMeans(n_clusters=k, random_state=0)
        labels = kmeans.fit_predict(X)
        
        # Elbow method: within-cluster sum of squares
        elbow_scores.append(kmeans.inertia_)
        
        # Silhouette method
        if k > 1:
            silhouette_avg = silhouette_score(X, labels)
            silhouette_scores.append(silhouette_avg)
        else:
            silhouette_scores.append(None)
        
        # Calinski-Harabasz score
        calinski_harabasz_avg = calinski_harabasz_score(X, labels)
        calinski_harabasz_scores.append(calinski_harabasz_avg)
        
        # Davies-Bouldin index
        davies_bouldin_avg = davies_bouldin_score(X, labels)
        davies_bouldin_scores.append(davies_bouldin_avg)
    
    return (elbow_scores, silhouette_scores, calinski_harabasz_scores, davies_bouldin_scores)

def find_optimal_k(scores, method):
    """
    Find the optimal number of clusters based on the given method.
    
    Parameters:
    scores (list): List of scores for different cluster numbers.
    method (str): Method name for printing.
    
    Returns:
    int: Optimal number of clusters.
    """
    if method == 'elbow':
        # Use the "elbow" heuristic
        diffs = np.diff(scores)
        second_diffs = np.diff(diffs)
        optimal_k = np.argmin(second_diffs) + 3  # Add 2 to account for zero-based index and initial k=2
    elif method == 'silhouette':
        # Silhouette score is maximized for better clustering
        if all(score is None for score in scores):
            return None
        optimal_k = np.argmax(scores) + 2  # +2 because range starts from 2
    elif method == 'calinski_harabasz':
        # Higher values are better
        optimal_k = np.argmax(scores) + 2  # +2 because range starts from 2
    elif method == 'davies_bouldin':
        # Lower values are better
        optimal_k = np.argmin(scores) + 2  # +2 because range starts from 2
    else:
        raise ValueError("Unknown method")
    
    return optimal_k

# Parameters
max_k = 15  # You can adjust this to explore more cluster numbers if needed

# Evaluate clustering
elbow_scores, silhouette_scores, calinski_harabasz_scores, davies_bouldin_scores = evaluate_clustering(X_scaled, max_k)

# Find optimal clusters for each method
optimal_elbow_k = find_optimal_k(elbow_scores, 'elbow')
optimal_silhouette_k = find_optimal_k(silhouette_scores, 'silhouette')
optimal_calinski_k = find_optimal_k(calinski_harabasz_scores, 'calinski_harabasz')
optimal_davies_k = find_optimal_k(davies_bouldin_scores, 'davies_bouldin')

# Print optimal k for each method
print(f"The optimal number of clusters according to the Elbow method is: {optimal_elbow_k}")
print(f"The optimal number of clusters according to the Silhouette method is: {optimal_silhouette_k}")
print(f"The optimal number of clusters according to the Calinski-Harabasz index is: {optimal_calinski_k}")
print(f"The optimal number of clusters according to the Davies-Bouldin index is: {optimal_davies_k}")

# Plotting for visual confirmation (optional)
plt.figure(figsize=(14, 10))

# Plot Elbow Method
plt.subplot(2, 2, 1)
plt.plot(range(2, max_k), elbow_scores, marker='o', label='Inertia')
plt.scatter(optimal_elbow_k, elbow_scores[optimal_elbow_k - 2], color='red', s=100, edgecolor='black', zorder=5, label='Optimal k')
plt.title('Elbow Method')
plt.xlabel('Number of clusters')
plt.ylabel('Within-cluster Sum of Squares (Inertia)')
plt.legend()

# Plot Silhouette Method
plt.subplot(2, 2, 2)
plt.plot(range(2, max_k), silhouette_scores, marker='o', label='Silhouette Score')
plt.scatter(optimal_silhouette_k, silhouette_scores[optimal_silhouette_k - 2], color='red', s=100, edgecolor='black', zorder=5, label='Optimal k')
plt.title('Silhouette Method')
plt.xlabel('Number of clusters')
plt.ylabel('Silhouette Score')
plt.legend()

# Plot Calinski-Harabasz Method
plt.subplot(2, 2, 3)
plt.plot(range(2, max_k), calinski_harabasz_scores, marker='o', label='Calinski-Harabasz Score')
plt.scatter(optimal_calinski_k, calinski_harabasz_scores[optimal_calinski_k - 2], color='red', s=100, edgecolor='black', zorder=5, label='Optimal k')
plt.title('Calinski-Harabasz Method')
plt.xlabel('Number of clusters')
plt.ylabel('Calinski-Harabasz Score')
plt.legend()

# Plot Davies-Bouldin Method
plt.subplot(2, 2, 4)
plt.plot(range(2, max_k), davies_bouldin_scores, marker='o', label='Davies-Bouldin Score')
plt.scatter(optimal_davies_k, davies_bouldin_scores[optimal_davies_k - 2], color='red', s=100, edgecolor='black', zorder=5, label='Optimal k')
plt.title('Davies-Bouldin Method')
plt.xlabel('Number of clusters')
plt.ylabel('Davies-Bouldin Score')
plt.legend()

#Generate the output file name for the plot by adding '_k_opt' suffix
base_name, _ = os.path.splitext(file_path)
plot_output_file = f"02_{base_name}_k_opt.png"


# Ensure the output directory exists
os.makedirs('02_k_optimization', exist_ok=True)

# Save the plot
plt.savefig(os.path.join('02_k_optimization', plot_output_file))



plt.tight_layout()
plt.show()

###ends#####
