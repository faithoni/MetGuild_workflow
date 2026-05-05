import os
try:
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
except NameError:
    pass

import json
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import spearmanr

# -----------------------------
# Config
# -----------------------------
# Filter: keep correlations >= 0.5 and p < 0.05 
CORR_THRESHOLD = 0.5
P_THRESHOLD = 0.05

# Output root requested
OUTPUT_ROOT = os.path.join("06_correlation", "02_featurevscluster_correlation")

# Summary files (overwrite if same names)
SUMMARY_CSV = os.path.join(OUTPUT_ROOT, "01_fvc_correlation_summary.csv")
SUMMARY_JSON = os.path.join(OUTPUT_ROOT, "01_fvc_correlation_summary.json")

# Inputs we will use from 05_stat_analysis (chi + anova)
SIGNIFICANT_INPUT_SUFFIXES = [
    "_chi_stat_sig.csv",              # Script 4
    "_anova_stat_significant.csv",     # Script 5
]

# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def load_manifest() -> dict:
    manifest_path = os.path.join("manifests", "optimal_k_values.json")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"Manifest not found at {manifest_path}. Run Script 1 first."
        )
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)

def load_cluster_feature_df(input_file: str) -> pd.DataFrame:
    """
    Load significant dataset CSV and ensure it has Cluster + features.
    Drops Genomes if present.
    Forces numeric for features.
    """
    df = pd.read_csv(input_file, encoding="latin1")

    if "Cluster" not in df.columns:
        raise ValueError(f"Missing 'Cluster' column in {input_file}")

    if "Genomes" in df.columns:
        df = df.drop(columns=["Genomes"])

    clusters = df["Cluster"]

    feature_cols = [c for c in df.columns if c != "Cluster"]
    features = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    features = features.dropna(axis=1, how="all")

    out = pd.concat([clusters, features], axis=1)
    out = out.dropna(subset=["Cluster"])
    return out

def spearman_corr_pval_binary(x: np.ndarray, y: np.ndarray):
    """
    Spearman correlation + p-value between binary indicator x (0/1) and feature y.
    Pairwise drops NaNs.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    x2 = x[mask]
    y2 = y[mask]

    if x2.size < 3:
        return np.nan, np.nan
    if np.nanstd(y2) == 0:
        return np.nan, np.nan

    corr, pval = spearmanr(x2, y2, nan_policy="omit")
    return float(corr), float(pval)

def _format_cluster_col(level) -> str:
    """Format cluster level as C1, C5, ..."""
    # If level is numeric-ish, format without .0
    try:
        lv = float(level)
        if lv.is_integer():
            return f"C{int(lv)}"
        return f"C{lv}"
    except Exception:
        # string cluster label
        return f"C{str(level)}"

def run_feature_vs_cluster(input_file: str, out_dir: str, out_prefix: str) -> dict:
    """
    Writes:
      - <prefix>_corr_pvalues_long.csv
      - <prefix>_filtered_corr_pvalues_long.csv
      - <prefix>_corr_matrix.csv
      - <prefix>_pvalues_matrix.csv

    Uses:
      - feature name column = 'feature'
      - cluster columns = C1, C2, ... (instead of Cluster_1, ...)
      - long-format columns = feature, cluster, Correlation, P_value

    Returns summary dict.
    """
    df = load_cluster_feature_df(input_file)

    clusters = df["Cluster"].to_numpy()
    feature_cols = [c for c in df.columns if c != "Cluster"]
    features_df = df[feature_cols]

    # Determine cluster levels (sorted if possible)
    unique_levels = pd.unique(clusters)
    try:
        cluster_levels = sorted(unique_levels.astype(float))
        cluster_levels = [int(x) if float(x).is_integer() else float(x) for x in cluster_levels]
        cluster_cmp_mode = "numeric"
    except Exception:
        cluster_levels = sorted([str(x) for x in unique_levels])
        cluster_cmp_mode = "string"

    feature_names = features_df.columns.to_list()

    if len(cluster_levels) == 0 or len(feature_names) == 0:
        return {
            "status": "skipped",
            "reason": "no clusters or no features",
            "clusters": int(len(cluster_levels)),
            "features": int(len(feature_names)),
            "rows_in_long_table": 0,
            "filtered_rows": 0,
        }

    cluster_cols = [_format_cluster_col(c) for c in cluster_levels]

    # Matrices: rows=feature, cols=C1/C2/...
    corr_mat = pd.DataFrame(data=np.nan, index=feature_names, columns=cluster_cols)
    pval_mat = pd.DataFrame(data=np.nan, index=feature_names, columns=cluster_cols)

    for level, colname in zip(cluster_levels, cluster_cols):
        if cluster_cmp_mode == "numeric":
            indicator = (clusters == level).astype(float)
        else:
            indicator = (clusters.astype(str) == str(level)).astype(float)

        for feat in feature_names:
            y = features_df[feat].to_numpy(dtype=float, copy=False)
            corr, pval = spearman_corr_pval_binary(indicator, y)
            corr_mat.loc[feat, colname] = corr
            pval_mat.loc[feat, colname] = pval

    # Long format with cluster column like C1, C2, ...
    corr_long = corr_mat.reset_index(names="feature").melt(
        id_vars="feature", var_name="cluster", value_name="Correlation"
    )
    p_long = pval_mat.reset_index(names="feature").melt(
        id_vars="feature", var_name="cluster", value_name="P_value"
    )
    merged = corr_long.merge(p_long, on=["feature", "cluster"], how="inner")

    filtered = merged[(merged["Correlation"] >= CORR_THRESHOLD) & (merged["P_value"] < P_THRESHOLD)].copy()

    ensure_dir(out_dir)

    out_full_long = os.path.join(out_dir, f"{out_prefix}_corr_pvalues_long.csv")
    out_filtered_long = os.path.join(out_dir, f"{out_prefix}_filtered_corr_pvalues_long.csv")
    out_corr_matrix = os.path.join(out_dir, f"{out_prefix}_corr_matrix.csv")
    out_pvals_matrix = os.path.join(out_dir, f"{out_prefix}_pvalues_matrix.csv")

    # Overwrite if exists (default behavior)
    merged.to_csv(out_full_long, index=False)
    filtered.to_csv(out_filtered_long, index=False)

    # Save matrices with first column named 'feature'
    corr_out = corr_mat.reset_index().rename(columns={"index": "feature"})
    pval_out = pval_mat.reset_index().rename(columns={"index": "feature"})
    corr_out.to_csv(out_corr_matrix, index=False)
    pval_out.to_csv(out_pvals_matrix, index=False)

    return {
        "status": "ok",
        "clusters": int(len(cluster_levels)),
        "features": int(len(feature_names)),
        "rows_in_long_table": int(len(merged)),
        "filtered_rows": int(len(filtered)),
        "out_full_long": out_full_long,
        "out_filtered_long": out_filtered_long,
        "out_corr_matrix": out_corr_matrix,
        "out_pvals_matrix": out_pvals_matrix,
    }

# -----------------------------
# Main
# -----------------------------
ensure_dir(OUTPUT_ROOT)

m = load_manifest()
k_values = m["unique_k_values"]
cleaned_csv = m["cleaned_output_file"]
base_stub = os.path.splitext(os.path.basename(cleaned_csv))[0]

print(f"Loaded k values from manifest: {k_values}")
print(f"Base stub: {base_stub}")
print(f"Thresholds: corr >= {CORR_THRESHOLD}, p < {P_THRESHOLD}")
print(f"Outputs will be under: {OUTPUT_ROOT}")

summary_rows = []
summary_json_entries = []

for k in k_values:
    print(f"\n=== Feature vs Cluster correlation for k={k} ===")

    # -----------------------------
    # KMEANS
    # -----------------------------
    kmeans_stat_dir = os.path.join("05_stat_analysis", "01_kmeans", f"{k}_clusters")
    kmeans_out_dir = os.path.join(OUTPUT_ROOT, "01_kmeans", f"{k}_clusters")
    ensure_dir(kmeans_out_dir)

    for suffix in SIGNIFICANT_INPUT_SUFFIXES:
        in_file = os.path.join(kmeans_stat_dir, f"05_{base_stub}_{k}_KM{suffix}")
        if not os.path.exists(in_file):
            continue

        # No "_chi" label; keep "_anova" label
        label = "" if suffix == "_chi_stat_sig.csv" else "_anova"
        out_prefix = f"{base_stub}_{k}_KM{label}"

        print(f"  [KMeans] Input: {in_file}")
        result = run_feature_vs_cluster(in_file, kmeans_out_dir, out_prefix)

        row = {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "method": "kmeans",
            "k": int(k),
            "source": "chi" if label == "" else "anova",
            "input_file": in_file,
            "output_dir": kmeans_out_dir,
            "output_prefix": out_prefix,
            "status": result.get("status"),
            "clusters": result.get("clusters"),
            "features": result.get("features"),
            "rows_in_long_table": result.get("rows_in_long_table"),
            "filtered_rows": result.get("filtered_rows"),
        }
        summary_rows.append(row)

        summary_json_entries.append({
            **row,
            "thresholds": {"corr": CORR_THRESHOLD, "p_value": P_THRESHOLD},
            "outputs": {
                "full_long": result.get("out_full_long"),
                "filtered_long": result.get("out_filtered_long"),
                "corr_matrix": result.get("out_corr_matrix"),
                "pvalues_matrix": result.get("out_pvals_matrix"),
            },
            "skip_reason": result.get("reason"),
        })

        if result.get("status") == "ok":
            print(f"    Clusters:           {result['clusters']}")
            print(f"    Features:           {result['features']}")
            print(f"    Rows in long table: {result['rows_in_long_table']}")
            print(f"    Filtered rows:      {result['filtered_rows']}")
        else:
            print(f"    Skipped: {result.get('reason', 'unknown reason')}")

    # -----------------------------
    # HIERARCHICAL
    # -----------------------------
    hier_stat_dir = os.path.join("05_stat_analysis", "02_hierarchical", f"{k}_clusters")
    hier_out_dir = os.path.join(OUTPUT_ROOT, "02_hierarchical", f"{k}_clusters")
    ensure_dir(hier_out_dir)

    for suffix in SIGNIFICANT_INPUT_SUFFIXES:
        in_file = os.path.join(hier_stat_dir, f"05_{base_stub}_{k}_HC{suffix}")
        if not os.path.exists(in_file):
            continue

        label = "" if suffix == "_chi_stat_sig.csv" else "_anova"
        out_prefix = f"{base_stub}_{k}_HC{label}"

        print(f"  [Hierarchical] Input: {in_file}")
        result = run_feature_vs_cluster(in_file, hier_out_dir, out_prefix)

        row = {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "method": "hierarchical",
            "k": int(k),
            "source": "chi" if label == "" else "anova",
            "input_file": in_file,
            "output_dir": hier_out_dir,
            "output_prefix": out_prefix,
            "status": result.get("status"),
            "clusters": result.get("clusters"),
            "features": result.get("features"),
            "rows_in_long_table": result.get("rows_in_long_table"),
            "filtered_rows": result.get("filtered_rows"),
        }
        summary_rows.append(row)

        summary_json_entries.append({
            **row,
            "thresholds": {"corr": CORR_THRESHOLD, "p_value": P_THRESHOLD},
            "outputs": {
                "full_long": result.get("out_full_long"),
                "filtered_long": result.get("out_filtered_long"),
                "corr_matrix": result.get("out_corr_matrix"),
                "pvalues_matrix": result.get("out_pvals_matrix"),
            },
            "skip_reason": result.get("reason"),
        })

        if result.get("status") == "ok":
            print(f"    Clusters:           {result['clusters']}")
            print(f"    Features:           {result['features']}")
            print(f"    Rows in long table: {result['rows_in_long_table']}")
            print(f"    Filtered rows:      {result['filtered_rows']}")
        else:
            print(f"    Skipped: {result.get('reason', 'unknown reason')}")

# -----------------------------
# Write summary files (overwrite)
# -----------------------------
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(SUMMARY_CSV, index=False)

summary_payload = {
    "created_at": datetime.now().isoformat(timespec="seconds"),
    "output_root": OUTPUT_ROOT,
    "thresholds": {"corr": CORR_THRESHOLD, "p_value": P_THRESHOLD},
    "runs": summary_json_entries,
}

with open(SUMMARY_JSON, "w", encoding="utf-8") as f:
    json.dump(summary_payload, f, indent=2)

print("\nAll feature-vs-cluster correlation runs completed.")
print(f"Summary CSV saved to:  {SUMMARY_CSV}")
print(f"Summary JSON saved to: {SUMMARY_JSON}")
