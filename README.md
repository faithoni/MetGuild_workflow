# Clustering → Statistics → Correlation → Networks (Cytoscape) Workflow

This repository contains a stepwise workflow to go from a genome feature presence/absence table to:
1) optimized clustering choices,  
2) clustering outputs (KMeans + hierarchical),  
3) statistical testing (Chi-square + ANOVA),  
4) correlation analyses (feature–feature, feature–cluster, cluster–cluster),  
5) network construction (edge lists + node attributes), and  
6) Cytoscape network export (PNG).

The workflow is designed so that **each script saves outputs in a predictable folder structure**, and downstream scripts automatically discover required inputs via a **manifest**.

---

## Folder Structure (high level)

- `manifests/`
  - `optimal_k_values.json` (created in Step 1; used by downstream scripts)
- `02_k_optimization/`
  - k-optimization plot(s) from Step 1
- `03_04_clustering/`
  - `01_kmeans/<k>_clusters/` (KMeans cluster outputs)
  - `02_hierarchical/<k>_clusters/` (Hierarchical cluster outputs)
- `05_stat_analysis/`
  - `01_kmeans/<k>_clusters/` (Chi-square + ANOVA outputs on KMeans results)
  - `02_hierarchical/<k>_clusters/` (Chi-square + ANOVA outputs on hierarchical results)
- `06_correlation/`
  - `01_featurevsfeature_correlation/`
  - `02_featurevscluster_correlation/`
  - `03_clustervscluster_correlation/`
- `07_network_construction/`
  - `edge_lists/` (network edge lists per dataset)
  - `nodes_attributes/` (node attribute tables per dataset)
  - `cytoscape_pngs/` (exported network PNGs)

---

## Input Data

Expected input is a CSV containing:
- a column named `Genomes`
- remaining columns are numeric/binary presence-absence features

Example: `01_DSK_DB.csv`

---

# Step-by-step Workflow

## Step 1 — Cleaning + K optimization (KMeans metrics)
**Script:** `01_<your_script_name>.py` (Cleaning + K optimization combined)

### What it does
1. Prompts for an input CSV filename.
2. Cleans the dataset:
   - keeps `Genomes`
   - drops all-zero columns
   - drops all-NaN columns
3. Saves a cleaned CSV:  
   `01_DSK_DB_cleaned.csv` (name depends on input)
4. Performs KMeans clustering evaluation (k = 2..14 by default) using:
   - Elbow (inertia)
   - Silhouette
   - Calinski–Harabasz
   - Davies–Bouldin
5. Saves the optimization plot under:
   - `02_k_optimization/`

### Key output for downstream steps
A manifest file is written:
- `manifests/optimal_k_values.json`

This manifest contains:
- the cleaned input file path
- the unique optimal k values to use downstream (duplicates removed)

> Example: if metrics yield `4, 2, 2, 11`, the unique set is `[2, 4, 11]`.

---

## Step 2 — KMeans clustering (for each k from the manifest)
**Script:** `03_kmeans_clustering.py` (your updated Script 2)

### What it does
For each `k` listed in `manifests/optimal_k_values.json`:
1. Loads the cleaned CSV saved in Step 1.
2. Computes a distance matrix (Jaccard recommended for binary features).
3. Runs MDS to 3 dimensions.
4. Runs KMeans clustering with `k` clusters.
5. Saves outputs to:
   - `03_04_clustering/01_kmeans/<k>_clusters/`
   - including:
     - clustered CSV (features + cluster assignments)
     - 2D plot PNG
     - 3D plot HTML

---

## Step 3 — Hierarchical clustering (for each k from the manifest)
**Script:** `04_hierarchical_clustering.py` (your updated Script 3)

### What it does
For each `k` in the manifest:
1. Loads cleaned CSV.
2. Computes a distance matrix (Euclidean by default; can be changed).
3. Performs hierarchical clustering (average linkage).
4. Cuts tree into `k` clusters.
5. Saves outputs to:
   - `03_04_clustering/02_hierarchical/<k>_clusters/`
   - including:
     - clustered CSV (features + cluster assignments)
     - dendrogram PNG

---

## Step 4 — Chi-square tests (run on each clustering output)
**Script:** `05_chi_square.py` (your updated chi-square script)

### What it does
For each clustering output CSV (KMeans + hierarchical, each `k`):
1. Performs Chi-square test of independence for each feature vs Cluster.
2. Saves:
   - full chi-square results
   - significant-only dataset (p < 0.05)
   - non-significant gene list

Outputs are written under:
- `05_stat_analysis/01_kmeans/<k>_clusters/`
- `05_stat_analysis/02_hierarchical/<k>_clusters/`

---

## Step 5 — ANOVA (run on each clustering output)
**Script:** `05_anova.py` (your updated ANOVA script)

### What it does
For each clustering output CSV (KMeans + hierarchical, each `k`):
1. Runs ANOVA per feature across cluster groups.
2. Saves:
   - full ANOVA results
   - filtered significant results (p < 0.05)
   - dataset with only statistically significant features

Outputs are written under:
- `05_stat_analysis/01_kmeans/<k>_clusters/`
- `05_stat_analysis/02_hierarchical/<k>_clusters/`

---

## Step 6 — Correlations
These are run on statistically significant datasets (from Step 4/5).

### Step 6A — Feature vs Feature correlation
**Script:** `06_feature_vs_feature_corr.py`

**What it does**
- Loads significant dataset (drops `Genomes`, `Cluster`)
- Computes Spearman correlation matrix + p-values
- Saves:
  - correlation matrix
  - flattened unique pairs
  - filtered pairs (|corr| >= 0.5 and p < 0.05)

Outputs:
- `06_correlation/01_featurevsfeature_correlation/<method>/<k>_clusters/`

Also writes summary files (overwritten each run):
- `06_correlation/01_featurevsfeature_correlation/01_fvf_correlation_summary.csv`
- `06_correlation/01_featurevsfeature_correlation/01_fvf_correlation_summary.json`

---

### Step 6B — Feature vs Cluster correlation
**Script:** `06_feature_vs_cluster_corr.py`

**What it does**
- Computes Spearman correlation between each feature and cluster membership indicator (one-vs-rest for each cluster)
- Uses cluster labels formatted as:
  - `C1`, `C2`, …
- Saves:
  - correlation matrix
  - p-value matrix
  - long-format merged table
  - filtered results (corr >= 0.5 and p < 0.05)

Outputs:
- `06_correlation/02_featurevscluster_correlation/<method>/<k>_clusters/`

---

### Step 6C — Cluster vs Cluster correlation
**Script:** `07_cluster_vs_cluster_corr.py`

**What it does**
- Builds a **cluster × feature** presence-absence summary matrix (max presence per cluster)
- Computes Spearman correlations between cluster columns
- Saves:
  - cluster correlation matrix
  - p-value matrix
  - flattened unique cluster pairs
  - filtered results (corr >= 0.5 and p < 0.05)

Outputs:
- `06_correlation/03_clustervscluster_correlation/<method>/<k>_clusters/`

---

## Step 7 — Network construction (edge lists + node attributes)
**Script:** `08_merge_filtered_files.py`

### What it does
For each dataset (each `k`, each method KM/HC):
1. Reads the filtered outputs from:
   - feature–feature
   - feature–cluster
   - cluster–cluster
2. Normalizes them into a common edge list schema:
   - `node1`, `node2`, `weight`, `p_value`, `source`
3. Removes duplicate edges (undirected dedupe, so A–B == B–A).
4. Writes edge lists to:
   - `07_network_construction/edge_lists/<prefix>_edge_list.csv`
5. Writes node attributes to:
   - `07_network_construction/nodes_attributes/<prefix>_node_attributes.csv`

### Node attributes format
Two columns:
- `Node`
- `Group`

Rules:
- Any cluster node is labeled `Group = cluster`
- everything else is `Group = feature`

Cluster nodes are uniquely named per dataset by appending the dataset label:
- `C1_KM2`, `C2_KM2`, …  
- `C1_HC4`, `C2_HC4`, …

This prevents cluster IDs from different networks colliding.

---

## Step 8 — Cytoscape import + PNG export (this scripts doesn't work! switch to the cytoscape app)
**Script:** `09_network_construction_cytoscape.py` 

### What it does
For each edge list in:
- `07_network_construction/edge_lists/`

It:
1. Imports the edge list as a Cytoscape network (`node1` → source, `node2` → target).
2. Imports the matching node attribute file into the node table (maps `Node` → Cytoscape node `name`).
3. Applies a style:
   - `Group=feature` vs `Group=cluster` mapped to different colors/shapes
4. Runs a layout (force-directed with fallback).
5. Exports a PNG to:
   - `07_network_construction/cytoscape_pngs/<prefix>.png`

### Requirements
- Cytoscape must be running
- CyREST must be enabled (default)
- If your CyREST port differs, update:
  - `CYREST = "http://127.0.0.1:1234/v1"`

---

# Running the Pipeline

Typical run order:

1. Step 1: Cleaning + K optimization  
2. Step 2: KMeans clustering (all k values)  
3. Step 3: Hierarchical clustering (all k values)  
4. Step 4: Chi-square (all outputs)  
5. Step 5: ANOVA (all outputs)  
6. Step 6: Correlations (all significant datasets)  
7. Step 7: Merge filtered correlations → edge list + node attributes  
8. Step 8: Cytoscape → PNG export (batch)

---

# Notes & Tips

- **Windows file locks:** If you open CSV outputs in Excel while rerunning scripts, Windows may lock the file. Close Excel before rerunning steps that overwrite outputs.
- **Overwriting behavior:** Most outputs are overwritten when filenames match; some scripts include fallbacks if locked.
- **Consistency of node names:** Cluster nodes are renamed with `_KMk` / `_HCk` suffix to prevent collisions.

---

# Outputs You’ll Typically Use Downstream

For each dataset (e.g., `01_DSK_DB_cleaned_4_KM`):

- Edge list:  
  `07_network_construction/edge_lists/01_DSK_DB_cleaned_4_KM_edge_list.csv`

- Node attributes:  
  `07_network_construction/nodes_attributes/01_DSK_DB_cleaned_4_KM_node_attributes.csv`

- Cytoscape PNG:  
  `07_network_construction/cytoscape_pngs/01_DSK_DB_cleaned_4_KM.png`
