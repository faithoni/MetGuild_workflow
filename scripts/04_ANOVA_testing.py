import os
import json
import pandas as pd
import numpy as np
from scipy.stats import f_oneway

# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

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
    - Ensure there's a 'Genomes' column
    - Ensure there's a 'Cluster' column
    - Return columns ordered: Genomes, Cluster, genes...
    """
    # Recover Genomes if it's saved as unnamed index column
    if "Genomes" not in df.columns:
        if "Unnamed: 0" in df.columns:
            df = df.rename(columns={"Unnamed: 0": "Genomes"})
        else:
            # if index is genomes, bring it out
            if df.index.name == "Genomes" or df.index.dtype == object:
                df = df.reset_index().rename(columns={"index": "Genomes"})

    if "Cluster" not in df.columns:
        raise ValueError("Column 'Cluster' not found in the clustered data.")

    # Reorder
    gene_cols = [c for c in df.columns if c not in ["Genomes", "Cluster"]]
    cols = (["Genomes"] if "Genomes" in df.columns else []) + ["Cluster"] + gene_cols
    return df[cols]

def run_anova_for_file(input_csv: str, output_dir: str, prefix: str) -> None:
    """
    Runs ANOVA per gene and writes:
      - <prefix>_anova_results.csv
      - <prefix>_anova_filtered_results.csv
      - <prefix>_anova_stat_significant.csv
      - <prefix>_anova_non_significant_genes.csv
    """
    if not os.path.exists(input_csv):
        print(f"SKIP (missing): {input_csv}")
        return

    df = pd.read_csv(input_csv, encoding="latin1")
    df = normalize_clustered_df(df)

    independent_variable = "Cluster"
    if independent_variable not in df.columns:
        raise ValueError(f"'{independent_variable}' column not found in {input_csv}")

    # dependent variables = everything except Genomes/Cluster
    dependent_variables = [c for c in df.columns if c not in ["Genomes", "Cluster"]]

    anova_results = []

    for column in dependent_variables:
        # Replace inf, drop NA
        df_filtered = df[[independent_variable, column]].replace([np.inf, -np.inf], np.nan).dropna()

        # Need at least two groups
        groups = [group[column].values for _, group in df_filtered.groupby(independent_variable)]
        if len(groups) > 1:
            try:
                f_value, p_value = f_oneway(*groups)
                anova_results.append({
                    "Dependent_Variable": column,
                    "F_Value": round(float(f_value), 2),
                    "Pr_Value": round(float(p_value), 2)
                })
            except Exception as e:
                print(f"[WARN] ANOVA failed for {column} in {input_csv}: {e}")

    anova_df = pd.DataFrame(anova_results)

    ensure_dir(output_dir)

    # Save full ANOVA results
    full_path = os.path.join(output_dir, f"{prefix}_anova_results.csv")
    anova_df.to_csv(full_path, index=False)

    # Filter significant
    filtered_anova_df = anova_df[anova_df["Pr_Value"] < 0.05].copy()
    filtered_path = os.path.join(output_dir, f"{prefix}_anova_filtered_results.csv")
    filtered_anova_df.to_csv(filtered_path, index=False)

    num_significant = int(len(filtered_anova_df))
    num_non_significant = int(len(anova_df) - num_significant)

    print(f"  Significant (p<0.05): {num_significant}")
    print(f"  Non-significant:      {num_non_significant}")

    # Save dataset with only significant genes
    significant_genes = filtered_anova_df["Dependent_Variable"].tolist()
    base_cols = ["Genomes", "Cluster"] if "Genomes" in df.columns else ["Cluster"]
    df_significant = df[base_cols + significant_genes] if significant_genes else df[base_cols].copy()

    sig_dataset_path = os.path.join(output_dir, f"{prefix}_anova_stat_significant.csv")
    df_significant.to_csv(sig_dataset_path, index=False)

    # Save list of non-significant genes
    non_significant_genes = sorted(set(anova_df["Dependent_Variable"]) - set(significant_genes))
    non_sig_path = os.path.join(output_dir, f"{prefix}_anova_non_significant_genes.csv")
    pd.DataFrame({"gene": non_significant_genes}).to_csv(non_sig_path, index=False)

    print(f"  Saved: {full_path}")
    print(f"  Saved: {filtered_path}")
    print(f"  Saved: {sig_dataset_path}")
    print(f"  Saved: {non_sig_path}")

# -----------------------------
# Main
# -----------------------------
m = load_manifest()
k_values = m["unique_k_values"]

cleaned_csv = m["cleaned_output_file"]
base_stub = os.path.splitext(os.path.basename(cleaned_csv))[0]

print(f"Loaded k values from manifest: {k_values}")
print(f"Base stub from cleaned CSV: {base_stub}")

for k in k_values:
    print(f"\n=== ANOVA for k={k} ===")

    # ----- KMEANS -----
    kmeans_input = os.path.join(
        "03_04_clustering", "01_kmeans", f"{k}_clusters", f"03_{base_stub}_{k}_KM.csv"
    )
    kmeans_outdir = os.path.join(
        "05_stat_analysis", "01_kmeans", f"{k}_clusters"
    )
    kmeans_prefix = f"05_{base_stub}_{k}_KM"
    print(f"[KMeans] Input: {kmeans_input}")
    run_anova_for_file(kmeans_input, kmeans_outdir, kmeans_prefix)

    # ----- HIERARCHICAL -----
    hier_input = os.path.join(
        "03_04_clustering", "02_hierarchical", f"{k}_clusters", f"04_{base_stub}_{k}_HC.csv"
    )
    hier_outdir = os.path.join(
        "05_stat_analysis", "02_hierarchical", f"{k}_clusters"
    )
    hier_prefix = f"05_{base_stub}_{k}_HC"
    print(f"[Hierarchical] Input: {hier_input}")
    run_anova_for_file(hier_input, hier_outdir, hier_prefix)

print("\nAll ANOVA runs completed.")
