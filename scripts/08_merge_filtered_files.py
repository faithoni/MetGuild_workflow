import os
try:
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
except NameError:
    pass

import json
import re
import pandas as pd
from datetime import datetime

# -----------------------------
# Config
# -----------------------------
EDGE_OUTPUT_DIR = os.path.join("07_network_construction", "edge_lists")
NODE_OUTPUT_DIR = os.path.join("07_network_construction", "nodes_attributes")
os.makedirs(EDGE_OUTPUT_DIR, exist_ok=True)
os.makedirs(NODE_OUTPUT_DIR, exist_ok=True)

SUMMARY_DIR = os.path.join("07_network_construction")
os.makedirs(SUMMARY_DIR, exist_ok=True)

MERGE_SUMMARY_CSV = os.path.join(SUMMARY_DIR, "01_edge_list_merge_summary.csv")
MERGE_SUMMARY_JSON = os.path.join(SUMMARY_DIR, "01_edge_list_merge_summary.json")

NODE_SUMMARY_CSV = os.path.join(SUMMARY_DIR, "02_node_attributes_summary.csv")
NODE_SUMMARY_JSON = os.path.join(SUMMARY_DIR, "02_node_attributes_summary.json")

FVSF_ROOT = os.path.join("06_correlation", "01_featurevsfeature_correlation")
FVSCL_ROOT = os.path.join("06_correlation", "02_featurevscluster_correlation")
CLVSCL_ROOT = os.path.join("06_correlation", "03_clustervscluster_correlation")

FVSF_FILES = [
    "{base}_{k}_{method}_filtered_corr_pvals.csv",
    "{base}_{k}_{method}_anova_filtered_corr_pvals.csv",
]
FVSCL_FILES = [
    "{base}_{k}_{method}_filtered_corr_pvalues_long.csv",
    "{base}_{k}_{method}_anova_filtered_corr_pvalues_long.csv",
]
CLVSCL_FILES = [
    "{base}_{k}_{method}_filtered_results.csv",
    "{base}_{k}_{method}_anova_filtered_results.csv",
]

CLUSTER_CORE_RE = re.compile(r"^C\d+(?:\.\d+)?$")

# -----------------------------
# Helpers
# -----------------------------
def load_manifest() -> dict:
    manifest_path = os.path.join("manifests", "optimal_k_values.json")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Manifest not found at {manifest_path}. Run Script 1 first.")
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def safe_read_csv(path: str):
    if not os.path.exists(path):
        return None, f"missing: {path}"
    try:
        return pd.read_csv(path), None
    except Exception as e:
        return None, f"read_error: {path} ({e})"

def build_paths(base: str, k: int, method: str):
    method_folder = "01_kmeans" if method == "KM" else "02_hierarchical"
    cluster_folder = f"{k}_clusters"

    fvsf_paths = [
        os.path.join(FVSF_ROOT, method_folder, cluster_folder, pat.format(base=base, k=k, method=method))
        for pat in FVSF_FILES
    ]
    fvscl_paths = [
        os.path.join(FVSCL_ROOT, method_folder, cluster_folder, pat.format(base=base, k=k, method=method))
        for pat in FVSCL_FILES
    ]
    clvscl_paths = [
        os.path.join(CLVSCL_ROOT, method_folder, cluster_folder, pat.format(base=base, k=k, method=method))
        for pat in CLVSCL_FILES
    ]
    return fvsf_paths, fvscl_paths, clvscl_paths

def normalize_edge_df(df: pd.DataFrame, source: str) -> pd.DataFrame:
    cols = {c.lower(): c for c in df.columns}

    if "feature" in cols and "cluster" in cols:
        out = pd.DataFrame({
            "node1": df[cols["feature"]].astype(str),
            "node2": df[cols["cluster"]].astype(str),
        })
        out["weight"] = df[cols.get("correlation", "Correlation")]
        out["p_value"] = df[cols.get("p_value", "P_value")]
        out["source"] = source
        return out

    if "var1" in cols and "var2" in cols:
        out = pd.DataFrame({
            "node1": df[cols["var1"]].astype(str),
            "node2": df[cols["var2"]].astype(str),
        })
        out["weight"] = df[cols.get("correlation", "Correlation")]
        out["p_value"] = df[cols.get("p_value", "P_value")]
        out["source"] = source
        return out

    out = pd.DataFrame({
        "node1": df.iloc[:, 0].astype(str),
        "node2": df.iloc[:, 1].astype(str),
    })
    out["weight"] = df[cols["correlation"]] if "correlation" in cols else pd.NA
    out["p_value"] = df[cols["p_value"]] if "p_value" in cols else pd.NA
    out["source"] = source
    return out

def filter_out_perfect_corr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove perfect correlations from the edge list.
    - Drops rows where abs(weight) == 1 (covers 1 and -1).
    - Also drops rows with missing/non-numeric weights.
    """
    out = df.copy()
    out["weight"] = pd.to_numeric(out["weight"], errors="coerce")
    out = out.dropna(subset=["weight"])
    out = out[out["weight"].abs() < 1.0]
    return out

def _strip_cluster_prefix(s: str) -> str:
    s = str(s).strip()
    if s.lower().startswith("cluster "):
        return s.split(" ", 1)[1].strip()
    return s

def append_cluster_suffix(node: str, method: str, k: int) -> str:
    core = _strip_cluster_prefix(node)
    if CLUSTER_CORE_RE.match(core):
        return f"{core}_{method}{k}"
    return str(node)

def apply_cluster_suffix_to_edges(df: pd.DataFrame, method: str, k: int) -> pd.DataFrame:
    out = df.copy()
    out["node1"] = out["node1"].map(lambda x: append_cluster_suffix(x, method, k))
    out["node2"] = out["node2"].map(lambda x: append_cluster_suffix(x, method, k))
    return out

def dedupe_undirected_edges(df: pd.DataFrame) -> pd.DataFrame:
    a = df["node1"].astype(str)
    b = df["node2"].astype(str)
    key_u = a.where(a <= b, b)
    key_v = b.where(a <= b, a)

    df2 = df.copy()
    df2["_u"] = key_u
    df2["_v"] = key_v
    df2 = df2.drop_duplicates(subset=["_u", "_v"], keep="first").drop(columns=["_u", "_v"])
    return df2

def build_node_attributes(edge_df: pd.DataFrame) -> pd.DataFrame:
    nodes = pd.unique(pd.concat([edge_df["node1"], edge_df["node2"]], ignore_index=True).astype(str))
    nodes = [n for n in nodes if n and n.lower() != "nan"]

    def group_of(n: str) -> str:
        return "cluster" if re.match(r"^C\d+(?:\.\d+)?_(?:KM|HC)\d+$", n) else "feature"

    attr = pd.DataFrame({"Node": nodes, "Group": [group_of(n) for n in nodes]})
    return attr.sort_values(by=["Group", "Node"]).reset_index(drop=True)

def safe_to_csv(df: pd.DataFrame, target_path: str, index: bool = False) -> str:
    try:
        df.to_csv(target_path, index=index)
        return target_path
    except PermissionError:
        base, ext = os.path.splitext(target_path)
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        fallback = f"{base}__{ts}{ext}"
        df.to_csv(fallback, index=index)
        return fallback

# -----------------------------
# Main
# -----------------------------
def main():
    ensure_dir(SUMMARY_DIR)

    m = load_manifest()
    k_values = m["unique_k_values"]
    cleaned_csv = m["cleaned_output_file"]
    base_stub = os.path.splitext(os.path.basename(cleaned_csv))[0]

    print(f"Loaded k values from manifest: {k_values}")
    print(f"Base stub: {base_stub}")
    print(f"Edge lists: {EDGE_OUTPUT_DIR}")
    print(f"Node attributes: {NODE_OUTPUT_DIR}")

    merge_summary_rows = []
    merge_summary_json_runs = []

    node_summary_rows = []
    node_summary_json_runs = []

    for k in k_values:
        for method in ["KM", "HC"]:
            prefix = f"{base_stub}_{k}_{method}"
            intended_edge_path = os.path.join(EDGE_OUTPUT_DIR, f"{prefix}_edge_list.csv")
            intended_node_path = os.path.join(NODE_OUTPUT_DIR, f"{prefix}_node_attributes.csv")

            fvsf_paths, fvscl_paths, clvscl_paths = build_paths(base_stub, k, method)

            dfs = []
            issues = []
            inputs_used = {"feature_vs_feature": [], "feature_vs_cluster": [], "cluster_vs_cluster": []}

            for p in fvsf_paths:
                df, err = safe_read_csv(p)
                if err:
                    issues.append(err)
                else:
                    src = "feature_vs_feature_anova" if "_anova_" in os.path.basename(p) else "feature_vs_feature"
                    dfs.append(normalize_edge_df(df, source=src))
                    inputs_used["feature_vs_feature"].append(p)

            for p in fvscl_paths:
                df, err = safe_read_csv(p)
                if err:
                    issues.append(err)
                else:
                    src = "feature_vs_cluster_anova" if "_anova_" in os.path.basename(p) else "feature_vs_cluster"
                    dfs.append(normalize_edge_df(df, source=src))
                    inputs_used["feature_vs_cluster"].append(p)

            for p in clvscl_paths:
                df, err = safe_read_csv(p)
                if err:
                    issues.append(err)
                else:
                    src = "cluster_vs_cluster_anova" if "_anova_" in os.path.basename(p) else "cluster_vs_cluster"
                    dfs.append(normalize_edge_df(df, source=src))
                    inputs_used["cluster_vs_cluster"].append(p)

            ts_now = datetime.now().isoformat(timespec="seconds")

            if not dfs:
                print(f"[{prefix}] SKIP: no input files found/readable.")
                merge_row = {
                    "timestamp": ts_now,
                    "k": int(k),
                    "method": "kmeans" if method == "KM" else "hierarchical",
                    "prefix": prefix,
                    "status": "skipped",
                    "reason": "no inputs",
                    "issues": "; ".join(issues) if issues else "",
                    "output_edge_list": intended_edge_path,
                    "output_node_attributes": intended_node_path,
                    "rows_combined": 0,
                    "rows_after_filter": 0,
                    "rows_deduped": 0,
                    "rows_removed_weight_eq_1": 0,
                    "nodes": 0,
                    "clusters": 0,
                    "features": 0,
                    "inputs_used": inputs_used,
                }
                merge_summary_rows.append(merge_row)
                merge_summary_json_runs.append(merge_row)

                node_row = {
                    "timestamp": ts_now,
                    "k": int(k),
                    "method": "kmeans" if method == "KM" else "hierarchical",
                    "prefix": prefix,
                    "status": "skipped",
                    "reason": "no inputs",
                    "node_attributes_file": intended_node_path,
                    "nodes_total": 0,
                    "nodes_feature": 0,
                    "nodes_cluster": 0,
                }
                node_summary_rows.append(node_row)
                node_summary_json_runs.append(node_row)
                continue

            combined = pd.concat(dfs, ignore_index=True).dropna(subset=["node1", "node2"])
            n_before = len(combined)

            # Remove perfect correlations (abs(weight) == 1)
            filtered = filter_out_perfect_corr(combined)
            n_after = len(filtered)
            n_removed = n_before - n_after

            filtered = apply_cluster_suffix_to_edges(filtered, method=method, k=int(k))
            deduped = dedupe_undirected_edges(filtered)

            edge_path_used = safe_to_csv(deduped, intended_edge_path, index=False)

            node_attr = build_node_attributes(deduped)
            node_path_used = safe_to_csv(node_attr, intended_node_path, index=False)

            n_nodes = int(len(node_attr))
            n_clusters = int((node_attr["Group"] == "cluster").sum())
            n_features = int((node_attr["Group"] == "feature").sum())

            print(
                f"[{prefix}] removed(|w|==1)={n_removed} | "
                f"edge_list: {edge_path_used} (combined={n_before}, after_filter={n_after}, deduped={len(deduped)}) | "
                f"node_attributes: {node_path_used} (nodes={n_nodes}, clusters={n_clusters}, features={n_features})"
            )

            merge_row = {
                "timestamp": ts_now,
                "k": int(k),
                "method": "kmeans" if method == "KM" else "hierarchical",
                "prefix": prefix,
                "status": "ok",
                "issues": "; ".join(issues) if issues else "",
                "output_edge_list": edge_path_used,
                "output_node_attributes": node_path_used,
                "rows_combined": int(n_before),
                "rows_after_filter": int(n_after),
                "rows_deduped": int(len(deduped)),
                "rows_removed_weight_eq_1": int(n_removed),
                "nodes": n_nodes,
                "clusters": n_clusters,
                "features": n_features,
                "inputs_used": inputs_used,
            }
            merge_summary_rows.append(merge_row)
            merge_summary_json_runs.append(merge_row)

            node_row = {
                "timestamp": ts_now,
                "k": int(k),
                "method": "kmeans" if method == "KM" else "hierarchical",
                "prefix": prefix,
                "status": "ok",
                "node_attributes_file": node_path_used,
                "nodes_total": n_nodes,
                "nodes_feature": n_features,
                "nodes_cluster": n_clusters,
            }
            node_summary_rows.append(node_row)
            node_summary_json_runs.append(node_row)

    pd.DataFrame(merge_summary_rows).to_csv(MERGE_SUMMARY_CSV, index=False)
    with open(MERGE_SUMMARY_JSON, "w", encoding="utf-8") as f:
        json.dump({
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "edge_output_dir": EDGE_OUTPUT_DIR,
            "node_output_dir": NODE_OUTPUT_DIR,
            "runs": merge_summary_json_runs,
        }, f, indent=2)

    pd.DataFrame(node_summary_rows).to_csv(NODE_SUMMARY_CSV, index=False)
    with open(NODE_SUMMARY_JSON, "w", encoding="utf-8") as f:
        json.dump({
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "node_output_dir": NODE_OUTPUT_DIR,
            "runs": node_summary_json_runs,
        }, f, indent=2)

    print("\nDone.")
    print(f"Merge summary CSV:  {MERGE_SUMMARY_CSV}")
    print(f"Merge summary JSON: {MERGE_SUMMARY_JSON}")
    print(f"Node summary CSV:   {NODE_SUMMARY_CSV}")
    print(f"Node summary JSON:  {NODE_SUMMARY_JSON}")

if __name__ == "__main__":
    main()
