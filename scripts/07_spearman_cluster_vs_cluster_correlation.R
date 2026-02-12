import os
import json
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.stats import spearmanr

# -----------------------------
# Config
# -----------------------------
CORR_THRESHOLD = 0.5
P_THRESHOLD = 0.05

OUTPUT_ROOT = os.path.join("06_correlation", "03_clustervscluster_correlation")
SUMMARY_CSV = os.path.join(OUTPUT_ROOT, "01_cvc_correlation_summary.csv")
SUMMARY_JSON = os.path.join(OUTPUT_ROOT, "01_cvc_correlation_summary.json")

# Use significant datasets from 05_stat_analysis
SIGNIFICANT_INPUT_SUFFIXES = [
    "_chi_stat_sig.csv",
    "_anova_stat_significant.csv",
]

# -----------------------------
# Helpers
# -----------------------------
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def load_manifest() -> dict:
    manifest_path = os.path.join("manifests", "optimal_k_values.json")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Manifest not found at {manifest_path}. Run Script 1 first.")
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)

def load_sig_df(input_file: str) -> pd.DataFrame:
    df = pd.read_csv(input_file, encoding="latin1")
    if "Cluster" not in df.columns:
        raise ValueError(f"Missing 'Cluster' column in {input_file}")
    return df

def build_cluster_feature_binary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Build a binary presence/absence matrix per (feature x cluster):
      long melt on Cluster -> max presence per (Cluster, feature) -> pivot wide
    Output: rows=feature, cols=C1,C2,... values 0/1
    """
    drop_cols = [c for c in ["Genomes", "Cluster"] if c in df.columns]
    feature_cols = [c for c in df.columns if c not in drop_cols]

    # numeric -> fill NA -> presence/absence
    features = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    features = (features > 0).astype(int)

    tmp = pd.concat([df["Cluster"], features], axis=1)

    long_df = tmp.melt(id_vars=["Cluster"], var_name="feature", value_name="Presence")
    grouped = long_df.groupby(["Cluster", "feature"], as_index=False)["Presence"].max()
    wide = grouped.pivot(index="feature", columns="Cluster", values="Presence").fillna(0).astype(int)

    # sort cluster columns if possible
    cols = list(wide.columns)
    try:
        cols_sorted = sorted(cols, key=lambda x: float(x))
    except Exception:
        cols_sorted = sorted(cols, key=lambda x: str(x))
    wide = wide[cols_sorted]

    # rename to C1, C2, ...
    def fmt_c(c):
        try:
            v = float(c)
            if v.is_integer():
                return f"C{int(v)}"
            return f"C{v}"
        except Exception:
            return f"C{str(c)}"

    wide.columns = [fmt_c(c) for c in wide.columns]
    return wide

def cluster_cluster_corr_pvals(wide: pd.DataFrame):
    """
    Compute Spearman correlation matrix + p-value matrix for cluster columns.
    Robust to:
      - constant columns (dropped)
      - too few feature rows
      - scipy returning scalar outputs
      - read-only backing arrays when setting diagonal
    Returns: (corr_df, pval_df, used_cluster_cols)
    """
    if wide.shape[1] < 2:
        raise ValueError("Need at least 2 clusters to compute cluster-vs-cluster correlation.")

    # Drop constant columns (no variance => Spearman undefined)
    nunique = wide.nunique(dropna=True)
    keep_cols = nunique[nunique > 1].index.tolist()
    wide2 = wide[keep_cols]

    if wide2.shape[1] < 2:
        raise ValueError(
            "After dropping constant cluster columns, <2 clusters remain. "
            "This usually happens when one or more clusters are all-0 or all-1 across features."
        )

    # Need at least 2 features (rows)
    if wide2.shape[0] < 2:
        raise ValueError(
            "Not enough features (rows) to compute cluster-vs-cluster correlation. "
            "Need at least 2 features after filtering."
        )

    corr, pvals = spearmanr(wide2.to_numpy(), axis=0, nan_policy="omit")
    corr = np.asarray(corr)
    pvals = np.asarray(pvals)

    if corr.ndim != 2 or pvals.ndim != 2:
        raise ValueError("Spearman returned non-matrix outputs (likely due to insufficient variability).")

    cols = wide2.columns.to_list()
    corr_df = pd.DataFrame(corr, index=cols, columns=cols)
    pval_df = pd.DataFrame(pvals, index=cols, columns=cols)

    # Make diagonal p-values NaN (avoid 'underlying array is read-only')
    p_arr = pval_df.to_numpy(copy=True)  # writable
    np.fill_diagonal(p_arr, np.nan)
    pval_df.iloc[:, :] = p_arr

    return corr_df, pval_df, cols

def flatten_unique_pairs(corr_df: pd.DataFrame, pval_df: pd.DataFrame) -> pd.DataFrame:
    """
    Unique pairs only; Var1/Var2 formatted as 'cluster C1', 'cluster C2', ...
    """
    cols = corr_df.columns.to_list()
    rows = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            v1, v2 = cols[i], cols[j]
            rows.append({
                "Var1": f"cluster {v1}",
                "Var2": f"cluster {v2}",
                "Correlation": float(corr_df.loc[v1, v2]),
                "P_value": float(pval_df.loc[v1, v2]),
            })
    return pd.DataFrame(rows)

def run_cluster_vs_cluster(input_file: str, out_dir: str, out_prefix: str) -> dict:
    """
    Outputs (overwrite by default):
      - <prefix>_corr_matrix.csv
      - <prefix>_p_value_matrix.csv
      - <prefix>_flattened_results.csv
      - <prefix>_filtered_results.csv
    """
    df = load_sig_df(input_file)
    wide = build_cluster_feature_binary(df)

    ensure_dir(out_dir)

    try:
        corr_df, pval_df, used_cols = cluster_cluster_corr_pvals(wide)
    except ValueError as e:
        return {
            "status": "skipped",
            "reason": str(e),
            "clusters_total": int(wide.shape[1]),
            "clusters_used": 0,
            "features": int(wide.shape[0]),
            "pairs": 0,
            "filtered_pairs": 0,
        }

    flat_df = flatten_unique_pairs(corr_df, pval_df)

    filtered_df = flat_df[
        (flat_df["Correlation"] >= CORR_THRESHOLD) &
        (flat_df["P_value"] < P_THRESHOLD)
    ].copy()

    out_corr = os.path.join(out_dir, f"{out_prefix}_corr_matrix.csv")
    out_pvals = os.path.join(out_dir, f"{out_prefix}_p_value_matrix.csv")
    out_flat = os.path.join(out_dir, f"{out_prefix}_flattened_results.csv")
    out_filt = os.path.join(out_dir, f"{out_prefix}_filtered_results.csv")

    corr_df.to_csv(out_corr, index=True)
    pval_df.to_csv(out_pvals, index=True)
    flat_df.to_csv(out_flat, index=False)
    filtered_df.to_csv(out_filt, index=False)

    return {
        "status": "ok",
        "clusters_total": int(wide.shape[1]),
        "clusters_used": int(len(used_cols)),
        "clusters_used_names": used_cols,
        "features": int(wide.shape[0]),
        "pairs": int(len(flat_df)),
        "filtered_pairs": int(len(filtered_df)),
        "out_corr_matrix": out_corr,
        "out_p_value_matrix": out_pvals,
        "out_flattened": out_flat,
        "out_filtered": out_filt,
    }

# -----------------------------
# Main
# -----------------------------
def main():
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
        print(f"\n=== Cluster vs Cluster correlation for k={k} ===")

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

            label = "" if suffix == "_chi_stat_sig.csv" else "_anova"
            out_prefix = f"{base_stub}_{k}_KM{label}"

            print(f"  [KMeans] Input: {in_file}")
            result = run_cluster_vs_cluster(in_file, kmeans_out_dir, out_prefix)

            row = {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "method": "kmeans",
                "k": int(k),
                "source": "chi" if label == "" else "anova",
                "input_file": in_file,
                "output_dir": kmeans_out_dir,
                "output_prefix": out_prefix,
                "status": result.get("status"),
                "clusters_total": result.get("clusters_total"),
                "clusters_used": result.get("clusters_used"),
                "features": result.get("features"),
                "pairs": result.get("pairs"),
                "filtered_pairs": result.get("filtered_pairs"),
            }
            summary_rows.append(row)
            summary_json_entries.append({
                **row,
                "thresholds": {"corr": CORR_THRESHOLD, "p_value": P_THRESHOLD},
                "outputs": {
                    "corr_matrix": result.get("out_corr_matrix"),
                    "p_value_matrix": result.get("out_p_value_matrix"),
                    "flattened": result.get("out_flattened"),
                    "filtered": result.get("out_filtered"),
                },
                "clusters_used_names": result.get("clusters_used_names"),
                "skip_reason": result.get("reason"),
            })

            if result.get("status") == "ok":
                print(f"    Clusters total:  {result['clusters_total']}")
                print(f"    Clusters used:   {result['clusters_used']} ({result.get('clusters_used_names')})")
                print(f"    Features:        {result['features']}")
                print(f"    Unique pairs:    {result['pairs']}")
                print(f"    Filtered pairs:  {result['filtered_pairs']}")
            else:
                print(f"    Skipped: {result.get('reason')}")

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
            result = run_cluster_vs_cluster(in_file, hier_out_dir, out_prefix)

            row = {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "method": "hierarchical",
                "k": int(k),
                "source": "chi" if label == "" else "anova",
                "input_file": in_file,
                "output_dir": hier_out_dir,
                "output_prefix": out_prefix,
                "status": result.get("status"),
                "clusters_total": result.get("clusters_total"),
                "clusters_used": result.get("clusters_used"),
                "features": result.get("features"),
                "pairs": result.get("pairs"),
                "filtered_pairs": result.get("filtered_pairs"),
            }
            summary_rows.append(row)
            summary_json_entries.append({
                **row,
                "thresholds": {"corr": CORR_THRESHOLD, "p_value": P_THRESHOLD},
                "outputs": {
                    "corr_matrix": result.get("out_corr_matrix"),
                    "p_value_matrix": result.get("out_p_value_matrix"),
                    "flattened": result.get("out_flattened"),
                    "filtered": result.get("out_filtered"),
                },
                "clusters_used_names": result.get("clusters_used_names"),
                "skip_reason": result.get("reason"),
            })

            if result.get("status") == "ok":
                print(f"    Clusters total:  {result['clusters_total']}")
                print(f"    Clusters used:   {result['clusters_used']} ({result.get('clusters_used_names')})")
                print(f"    Features:        {result['features']}")
                print(f"    Unique pairs:    {result['pairs']}")
                print(f"    Filtered pairs:  {result['filtered_pairs']}")
            else:
                print(f"    Skipped: {result.get('reason')}")

    # Write summary (overwrite)
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

    print("\nAll cluster-vs-cluster correlation runs completed.")
    print(f"Summary CSV saved to:  {SUMMARY_CSV}")
    print(f"Summary JSON saved to: {SUMMARY_JSON}")

if __name__ == "__main__":
    main()
