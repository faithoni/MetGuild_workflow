import os
try:
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
except NameError:
    pass

import json
import pandas as pd
from scipy.stats import chi2_contingency

# -----------------------------
# Helpers
# -----------------------------
def ensure_dir_for_file(file_path: str) -> None:
    """Ensure directory for a file path exists."""
    directory = os.path.dirname(file_path)
    if directory and not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)

def load_manifest() -> dict:
    manifest_path = os.path.join("manifests", "optimal_k_values.json")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"Manifest not found at {manifest_path}. "
            "Run Script 1 first (it creates this file)."
        )
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)

def normalize_clustered_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize different clustering outputs into a consistent shape:
    - Ensure there's a 'Genomes' column (if present as index or unnamed col)
    - Ensure there's a 'Cluster' column
    - Keep only ['Genomes','Cluster', gene columns...]
    """
    # If 'Genomes' is not a column, try to recover it from index or unnamed first column
    if "Genomes" not in df.columns:
        # common case: saved with index -> first col is 'Unnamed: 0'
        if "Unnamed: 0" in df.columns:
            df = df.rename(columns={"Unnamed: 0": "Genomes"})
        else:
            # if index looks meaningful, bring it out
            if df.index.name == "Genomes" or df.index.dtype == object:
                df = df.reset_index().rename(columns={"index": "Genomes"})

    if "Cluster" not in df.columns:
        raise ValueError("Column 'Cluster' not found in the clustered data.")

    # Make sure Cluster is category (as your original script did)
    df["Cluster"] = df["Cluster"].astype("category")

    # Reorder columns if possible
    base_cols = []
    if "Genomes" in df.columns:
        base_cols.append("Genomes")
    base_cols.append("Cluster")

    gene_cols = [c for c in df.columns if c not in base_cols]
    return df[base_cols + gene_cols]

def perform_chi_squared_test(df: pd.DataFrame, gene_column: str) -> float:
    # Missing check
    if df[gene_column].isna().any():
        raise ValueError(f"Missing values found in gene column: {gene_column}")

    # Ensure binary (0/1)
    if not df[gene_column].isin([0, 1]).all():
        raise ValueError(f"Non-binary values found in gene column: {gene_column}")

    contingency_table = pd.crosstab(df["Cluster"], df[gene_column])

    # If contingency table is degenerate (e.g., only one column), chi2 can fail or be meaningless.
    # chi2_contingency handles many cases, but we can skip fully-constant genes safely.
    if contingency_table.shape[1] < 2:
        # gene always 0 or always 1 across all clusters
        return 1.0

    chi2_stat, p_value, _, _ = chi2_contingency(contingency_table)
    return float(p_value)

def run_chi_square_for_file(input_csv: str, output_dir: str, prefix: str) -> None:
    """
    Run chi-square workflow for one clustered CSV.
    Writes:
      - <prefix>_chi_squared_results.csv
      - <prefix>_chi_stat_sig.csv
      - <prefix>_chi_non_sig.csv
    """
    if not os.path.exists(input_csv):
        print(f"SKIP (missing): {input_csv}")
        return

    df = pd.read_csv(input_csv, encoding="latin1")
    df = normalize_clustered_df(df)

    # gene columns = everything except Genomes/Cluster
    gene_columns = [c for c in df.columns if c not in ["Genomes", "Cluster"]]

    # If your gene columns should start at "3rd column onwards" like before,
    # the above is equivalent but more robust across formats.

    results = []
    for gene in gene_columns:
        try:
            p_value = perform_chi_squared_test(df, gene)
            p_value_rounded = round(p_value, 2)
            gene_name = f"*{gene}" if p_value < 0.05 else gene
            results.append({"gene": gene_name, "p_value": p_value_rounded})
        except ValueError as e:
            print(f"[WARN] {e} (file={input_csv})")

    results_df = pd.DataFrame(results)

    # Output paths
    output_chi2_csv = os.path.join(output_dir, f"{prefix}_chi_squared_results.csv")
    output_significant_csv = os.path.join(output_dir, f"{prefix}_chi_stat_sig.csv")
    output_non_significant_csv = os.path.join(output_dir, f"{prefix}_chi_non_sig.csv")

    ensure_dir_for_file(output_chi2_csv)
    ensure_dir_for_file(output_significant_csv)
    ensure_dir_for_file(output_non_significant_csv)

    # Save chi2 p-values
    results_df.to_csv(output_chi2_csv, index=False)

    # Significant mask
    if not results_df.empty:
        significant_genes_mask = results_df["p_value"] < 0.05
        num_significant = int(significant_genes_mask.sum())
        num_non_significant = int(len(significant_genes_mask) - num_significant)
    else:
        significant_genes_mask = pd.Series([], dtype=bool)
        num_significant = 0
        num_non_significant = 0

    print(f"  Significant genes: {num_significant}")
    print(f"  Non-significant genes: {num_non_significant}")

    # Filter original data to keep only significant genes (+ Genomes/Cluster)
    # Note: results_df["gene"] may have '*' prefix, strip it
    non_significant_gene_columns = results_df.loc[~significant_genes_mask, "gene"].astype(str).str.lstrip("*").tolist()

    filtered_df = df.drop(columns=[c for c in non_significant_gene_columns if c in df.columns], errors="ignore")
    filtered_df.to_csv(output_significant_csv, index=False)

    pd.DataFrame({"gene": non_significant_gene_columns}).to_csv(output_non_significant_csv, index=False)

    print(f"  Saved: {output_chi2_csv}")
    print(f"  Saved: {output_significant_csv}")
    print(f"  Saved: {output_non_significant_csv}")

# -----------------------------
# Main
# -----------------------------
m = load_manifest()
k_values = m["unique_k_values"]

cleaned_csv = m["cleaned_output_file"]
base_stub = os.path.splitext(os.path.basename(cleaned_csv))[0]

print(f"Loaded k values from manifest: {k_values}")
print(f"Base stub from cleaned CSV: {base_stub}")

# Patterns from your scripts:
# KMeans: 03_04_clustering/01_kmeans/{k}_clusters/03_{base_stub}_{k}_KM.csv
# Hier:   03_04_clustering/02_hierarchical/{k}_clusters/04_{base_stub}_{k}_HC.csv

for k in k_values:
    print(f"\n=== Chi-square for k={k} ===")

    # ----- KMEANS -----
    kmeans_input = os.path.join(
        "03_04_clustering", "01_kmeans", f"{k}_clusters", f"03_{base_stub}_{k}_KM.csv"
    )
    kmeans_outdir = os.path.join(
        "05_stat_analysis", "01_kmeans", f"{k}_clusters"
    )
    kmeans_prefix = f"05_{base_stub}_{k}_KM"
    print(f"[KMeans] Input: {kmeans_input}")
    run_chi_square_for_file(kmeans_input, kmeans_outdir, kmeans_prefix)

    # ----- HIERARCHICAL -----
    hier_input = os.path.join(
        "03_04_clustering", "02_hierarchical", f"{k}_clusters", f"04_{base_stub}_{k}_HC.csv"
    )
    hier_outdir = os.path.join(
        "05_stat_analysis", "02_hierarchical", f"{k}_clusters"
    )
    hier_prefix = f"05_{base_stub}_{k}_HC"
    print(f"[Hierarchical] Input: {hier_input}")
    run_chi_square_for_file(hier_input, hier_outdir, hier_prefix)

print("\nAll chi-square runs completed.")
