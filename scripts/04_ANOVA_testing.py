import os
try:
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
except NameError:
    pass

import json
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import f_oneway

# -----------------------------
# Config
# -----------------------------
P_THRESHOLD = 0.05

MANIFEST_PATH = os.path.join("manifests", "optimal_k_values.json")

# Clustering outputs
KMEANS_CLUSTER_DIR = os.path.join("03_04_clustering", "01_kmeans")
HIER_CLUSTER_DIR = os.path.join("03_04_clustering", "02_hierarchical")

# Statistical outputs
KMEANS_STAT_DIR = os.path.join("05_stat_analysis", "01_kmeans")
HIER_STAT_DIR = os.path.join("05_stat_analysis", "02_hierarchical")


# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_manifest() -> dict:
    if not os.path.exists(MANIFEST_PATH):
        raise FileNotFoundError(f"Manifest not found at {MANIFEST_PATH}. Run Script 1 first.")
    with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


def find_clustered_csv(base_stub: str, k: int, method: str) -> str:
    """
    method: 'KM' or 'HC'
    Expected filenames from your clustering scripts:
      - KMeans: 03_<base_stub>_<k>_KM.csv in 03_04_clustering/01_kmeans/<k>_clusters/
      - Hier:  04_<base_stub>_<k>_HC.csv in 03_04_clustering/02_hierarchical/<k>_clusters/
    """
    if method == "KM":
        folder = os.path.join(KMEANS_CLUSTER_DIR, f"{k}_clusters")
        fname = f"03_{base_stub}_{k}_KM.csv"
    else:
        folder = os.path.join(HIER_CLUSTER_DIR, f"{k}_clusters")
        fname = f"04_{base_stub}_{k}_HC.csv"

    path = os.path.join(folder, fname)
    return path


def run_anova(df: pd.DataFrame, cluster_col: str = "Cluster") -> pd.DataFrame:
    """
    Runs ANOVA for each feature column (assumes 0/1 or numeric) grouped by Cluster.
    Returns a dataframe with columns:
      Dependent_Variable, F_Value, Pr_Value
    """
    if cluster_col not in df.columns:
        raise ValueError("Column 'Cluster' not found in input dataframe.")

    # Keep Genomes if present; but analysis uses only Cluster + features
    drop_cols = [c for c in ["Genomes"] if c in df.columns]
    work = df.drop(columns=drop_cols, errors="ignore").copy()

    # Ensure Cluster is treated as categorical-ish
    work[cluster_col] = work[cluster_col].astype("category")

    feature_cols = [c for c in work.columns if c != cluster_col]

    results = []
    for col in feature_cols:
        # numeric clean
        tmp = work[[cluster_col, col]].replace([np.inf, -np.inf], np.nan).dropna()
        tmp[col] = pd.to_numeric(tmp[col], errors="coerce")
        tmp = tmp.dropna()

        groups = [g[col].values for _, g in tmp.groupby(cluster_col)]
        if len(groups) > 1:
            f_value, p_value = f_oneway(*groups)
            results.append({
                "Dependent_Variable": col,
                "F_Value": round(float(f_value), 4),
                "Pr_Value": round(float(p_value), 6),
            })

    return pd.DataFrame(results)


def save_outputs(base_stub: str, k: int, method: str, df_original: pd.DataFrame, anova_df: pd.DataFrame) -> dict:
    """
    Writes outputs with the EXACT names your downstream scripts expect.

    Outputs written to:
      05_stat_analysis/<01_kmeans or 02_hierarchical>/<k>_clusters/

    Files:
      05_<base_stub>_<k>_<KM|HC>_anova_results.csv
      05_<base_stub>_<k>_<KM|HC>_anova_filtered_results.csv
      05_<base_stub>_<k>_<KM|HC>_anova_stat_significant.csv   <-- THIS is the key file used downstream
      05_<base_stub>_<k>_<KM|HC>_anova_non_sig_features.csv
    """
    out_root = KMEANS_STAT_DIR if method == "KM" else HIER_STAT_DIR
    out_dir = os.path.join(out_root, f"{k}_clusters")
    ensure_dir(out_dir)

    # Base prefix for filenames
    prefix = f"05_{base_stub}_{k}_{method}"

    out_results = os.path.join(out_dir, f"{prefix}_anova_results.csv")
    out_filtered = os.path.join(out_dir, f"{prefix}_anova_filtered_results.csv")
    out_sig_dataset = os.path.join(out_dir, f"{prefix}_anova_stat_significant.csv")
    out_non_sig = os.path.join(out_dir, f"{prefix}_anova_non_sig_features.csv")

    # Filter significant ANOVA rows
    filtered = anova_df[anova_df["Pr_Value"] < P_THRESHOLD].copy()

    # Save results tables (overwrite)
    anova_df.to_csv(out_results, index=False)
    filtered.to_csv(out_filtered, index=False)

    # Build significant-only dataset (keep Genomes+Cluster if present)
    sig_features = filtered["Dependent_Variable"].tolist()

    keep_cols = []
    if "Genomes" in df_original.columns:
        keep_cols.append("Genomes")
    if "Cluster" in df_original.columns:
        keep_cols.append("Cluster")

    df_sig = df_original[keep_cols + sig_features].copy() if sig_features else df_original[keep_cols].copy()
    df_sig.to_csv(out_sig_dataset, index=False)

    # Save non-significant feature list
    non_sig = sorted(set(anova_df["Dependent_Variable"]) - set(sig_features))
    pd.DataFrame({"feature": non_sig}).to_csv(out_non_sig, index=False)

    return {
        "out_dir": out_dir,
        "anova_results": out_results,
        "anova_filtered": out_filtered,
        "anova_stat_significant": out_sig_dataset,
        "anova_non_sig_features": out_non_sig,
        "n_features_tested": int(len(anova_df)),
        "n_significant": int(len(filtered)),
        "n_non_significant": int(len(anova_df) - len(filtered)),
    }


# -----------------------------
# Main
# -----------------------------
def main():
    m = load_manifest()
    k_values = m["unique_k_values"]
    cleaned_csv = m["cleaned_output_file"]
    base_stub = os.path.splitext(os.path.basename(cleaned_csv))[0]

    print(f"Loaded k values: {k_values}")
    print(f"Base stub: {base_stub}")
    print(f"p-threshold: {P_THRESHOLD}")

    for k in k_values:
        for method in ["KM", "HC"]:
            clustered_csv = find_clustered_csv(base_stub, k, method)
            tag = "KMeans" if method == "KM" else "Hierarchical"

            if not os.path.exists(clustered_csv):
                print(f"\n[{tag} k={k}] SKIP (missing clustering output): {clustered_csv}")
                continue

            print(f"\n[{tag} k={k}] Input: {clustered_csv}")

            df = pd.read_csv(clustered_csv, encoding="latin1")

            # Run ANOVA
            anova_df = run_anova(df, cluster_col="Cluster")
            if anova_df.empty:
                print(f"  WARNING: ANOVA produced no results (possibly single cluster or no usable features).")
                # Still write empty outputs to keep downstream predictable
                anova_df = pd.DataFrame(columns=["Dependent_Variable", "F_Value", "Pr_Value"])

            out = save_outputs(base_stub, k, method, df, anova_df)

            print(f"  Tested features:      {out['n_features_tested']}")
            print(f"  Significant (p<0.05): {out['n_significant']}")
            print(f"  Non-significant:      {out['n_non_significant']}")
            print(f"  Wrote: {out['anova_stat_significant']}")

    print("\nDone.")


if __name__ == "__main__":
    main()
