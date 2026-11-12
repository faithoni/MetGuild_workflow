#!/usr/bin/env python3
"""MetGuild Step 8: merge filtered correlation edges into network files.

If Step 3b taxonomy annotation is available, cluster nodes receive one additional
MAG_Phyla attribute containing every MAG-phylum pair represented in that cluster.

This stage combines the filtered outputs from:
  Step 5 - feature-feature correlations
  Step 6 - feature-cluster correlations
  Step 7 - cluster-cluster correlations

Naming convention
-----------------
Correlation-stage outputs use:
  * features: exact original feature names
  * clusters: C1, C2, ..., Ck

Only at final network construction are cluster nodes made solution-specific:
  * C1_KM3, C2_KM3, ... for KMeans k=3
  * C1_HC3, C2_HC3, ... for hierarchical k=3

Feature names are never modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PERFECT_ATOL = 1e-12


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"Required manifest not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def expected_cluster_names(k: int) -> set[str]:
    return {f"C{i}" for i in range(1, k + 1)}


def validate_cluster_name(value: object, k: int, context: str) -> str:
    name = str(value).strip()
    if name not in expected_cluster_names(k):
        raise ValueError(
            f"{context}: invalid cluster name {name!r}. "
            f"Expected one of {sorted(expected_cluster_names(k))}."
        )
    return name


def cluster_network_name(cluster: str, method_code: str, k: int) -> str:
    validate_cluster_name(cluster, k, "cluster naming")
    return f"{cluster}_{method_code}{k}"


def read_filtered_table(path: Path, required: set[str]) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Expected filtered correlation file not found: {path}")
    df = pd.read_csv(path)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")
    return df


def numeric_series(df: pd.DataFrame, column: str, path: Path) -> pd.Series:
    s = pd.to_numeric(df[column], errors="coerce")
    if s.isna().any() or not np.isfinite(s.to_numpy(dtype=float)).all():
        raise ValueError(f"{path}: column {column!r} contains missing/non-numeric/non-finite values.")
    return s.astype(float)


def normalize_fvf(path: Path) -> pd.DataFrame:
    """Normalize Step 5 feature-feature filtered edges."""
    df = read_filtered_table(
        path,
        {"Var1", "Var2", "Correlation", "P_value", "Adjusted_P_Value"},
    )
    out = pd.DataFrame({
        "node1": df["Var1"].astype(str),
        "node2": df["Var2"].astype(str),
        "node1_type": "feature",
        "node2_type": "feature",
        "correlation": numeric_series(df, "Correlation", path),
        "p_value": numeric_series(df, "P_value", path),
        "adjusted_p_value": numeric_series(df, "Adjusted_P_Value", path),
        "edge_type": "feature_feature",
    })
    return out


def normalize_fvc(path: Path, k: int) -> pd.DataFrame:
    """Normalize Step 6 feature-cluster filtered edges."""
    df = read_filtered_table(
        path,
        {"feature", "cluster", "Correlation", "P_value", "Adjusted_P_Value"},
    )
    clusters = [validate_cluster_name(v, k, str(path)) for v in df["cluster"]]
    out = pd.DataFrame({
        "node1": df["feature"].astype(str),
        "node2": clusters,
        "node1_type": "feature",
        "node2_type": "cluster",
        "correlation": numeric_series(df, "Correlation", path),
        "p_value": numeric_series(df, "P_value", path),
        "adjusted_p_value": numeric_series(df, "Adjusted_P_Value", path),
        "edge_type": "feature_cluster",
    })
    return out


def normalize_cvc(path: Path, k: int) -> pd.DataFrame:
    """Normalize Step 7 cluster-cluster filtered edges."""
    df = read_filtered_table(
        path,
        {"Var1", "Var2", "Correlation", "P_value", "Adjusted_P_Value"},
    )
    node1 = [validate_cluster_name(v, k, str(path)) for v in df["Var1"]]
    node2 = [validate_cluster_name(v, k, str(path)) for v in df["Var2"]]
    out = pd.DataFrame({
        "node1": node1,
        "node2": node2,
        "node1_type": "cluster",
        "node2_type": "cluster",
        "correlation": numeric_series(df, "Correlation", path),
        "p_value": numeric_series(df, "P_value", path),
        "adjusted_p_value": numeric_series(df, "Adjusted_P_Value", path),
        "edge_type": "cluster_cluster",
    })
    return out


def suffix_cluster_nodes(df: pd.DataFrame, method_code: str, k: int) -> pd.DataFrame:
    out = df.copy()
    mask1 = out["node1_type"].eq("cluster")
    mask2 = out["node2_type"].eq("cluster")
    out.loc[mask1, "node1"] = out.loc[mask1, "node1"].map(
        lambda x: cluster_network_name(str(x), method_code, k)
    )
    out.loc[mask2, "node2"] = out.loc[mask2, "node2"].map(
        lambda x: cluster_network_name(str(x), method_code, k)
    )
    return out


def remove_invalid_edges(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Remove missing, self-loop, and perfect-correlation edges as a safety check."""
    out = df.copy()
    before = len(out)

    missing_node = (
        out["node1"].isna() | out["node2"].isna()
        | out["node1"].astype(str).str.strip().eq("")
        | out["node2"].astype(str).str.strip().eq("")
    )
    n_missing = int(missing_node.sum())
    out = out.loc[~missing_node].copy()

    self_loop = out["node1"].astype(str).eq(out["node2"].astype(str))
    n_self = int(self_loop.sum())
    out = out.loc[~self_loop].copy()

    perfect = np.isclose(out["correlation"].abs(), 1.0, atol=PERFECT_ATOL, rtol=0.0)
    n_perfect = int(perfect.sum())
    out = out.loc[~perfect].copy()

    return out, {
        "rows_input": int(before),
        "removed_missing_node": n_missing,
        "removed_self_loops": n_self,
        "removed_perfect_correlations": n_perfect,
    }


def canonical_pair_key(df: pd.DataFrame) -> pd.Series:
    a = df["node1"].astype(str)
    b = df["node2"].astype(str)
    left = a.where(a <= b, b)
    right = b.where(a <= b, a)
    return left + "\x1f" + right


def assert_no_duplicate_undirected_edges(df: pd.DataFrame) -> None:
    if df.empty:
        return
    key = canonical_pair_key(df)
    duplicated = key.duplicated(keep=False)
    if duplicated.any():
        examples = df.loc[duplicated, ["node1", "node2", "edge_type"]].head(10).to_dict("records")
        raise ValueError(
            "Duplicate undirected endpoint pairs were produced while merging correlation sources. "
            f"This should not happen with type-aware naming. Examples: {examples}"
        )


def taxonomy_run_map(taxonomy_manifest: dict) -> dict[tuple[str, int], dict]:
    out: dict[tuple[str, int], dict] = {}
    for run in taxonomy_manifest.get("runs", []):
        code = str(run.get("method_code"))
        k = int(run.get("requested_k"))
        out[(code, k)] = run
    return out


def load_cluster_taxonomy(path: Path, k: int) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Expected taxonomy annotation file not found: {path}")
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"Cluster", "MAG_Phyla"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing taxonomy columns: {sorted(missing)}")
    if frame["Cluster"].duplicated().any():
        dup = frame.loc[frame["Cluster"].duplicated(keep=False), "Cluster"].tolist()
        raise ValueError(f"{path}: duplicate cluster rows: {dup}")

    lookup: dict[str, str] = {}
    for cluster, mag_phyla in frame[["Cluster", "MAG_Phyla"]].itertuples(index=False, name=None):
        cluster_name = validate_cluster_name(cluster, k, str(path))
        value = str(mag_phyla).strip()
        if not value:
            raise ValueError(f"{path}: MAG_Phyla is empty for {cluster_name}.")
        lookup[cluster_name] = value

    expected = expected_cluster_names(k)
    if set(lookup) != expected:
        raise ValueError(
            f"{path}: taxonomy annotation must contain exactly {sorted(expected)}; "
            f"observed {sorted(lookup)}"
        )
    return lookup


def build_node_attributes(
    edges: pd.DataFrame,
    method_code: str,
    k: int,
    taxonomy_lookup: dict[str, str] | None = None,
) -> pd.DataFrame:
    records: dict[str, dict] = {}
    for endpoint, type_col in [("node1", "node1_type"), ("node2", "node2_type")]:
        for name, node_type in zip(edges[endpoint].astype(str), edges[type_col].astype(str)):
            if name in records and records[name]["Group"] != node_type:
                raise ValueError(
                    f"Node-name collision: {name!r} is being used as both a feature and a cluster. "
                    "Feature names are intentionally preserved exactly; rename the conflicting feature upstream."
                )
            if node_type == "cluster":
                suffix = f"_{method_code}{k}"
                if not name.endswith(suffix):
                    raise ValueError(f"Unexpected final cluster node name: {name!r}")
                original = name[:-len(suffix)]
                mag_phyla = ""
                if taxonomy_lookup is not None:
                    if original not in taxonomy_lookup:
                        raise ValueError(
                            f"Taxonomy annotation is missing cluster {original} for {method_code}, k={k}."
                        )
                    mag_phyla = taxonomy_lookup[original]
                records[name] = {
                    "Node": name,
                    "Group": "cluster",
                    "Original_Name": original,
                    "Method": method_code,
                    "K": int(k),
                    "MAG_Phyla": mag_phyla,
                }
            else:
                records[name] = {
                    "Node": name,
                    "Group": "feature",
                    "Original_Name": name,
                    "Method": "",
                    "K": pd.NA,
                    "MAG_Phyla": "",
                }
    columns = ["Node", "Group", "Original_Name", "Method", "K", "MAG_Phyla"]
    if not records:
        return pd.DataFrame(columns=columns)
    return (
        pd.DataFrame(records.values(), columns=columns)
        .sort_values(["Group", "Node"], kind="stable")
        .reset_index(drop=True)
    )


def stage_run_map(stage_manifest: dict) -> dict[tuple[str, int], dict]:
    out: dict[tuple[str, int], dict] = {}
    for r in stage_manifest.get("runs", []):
        code = str(r.get("method_code"))
        k = int(r.get("requested_k"))
        out[(code, k)] = r
    return out


def source_file_for(
    *,
    source: str,
    run_dir: Path,
    base_stub: str,
    method_dir: str,
    method_code: str,
    k: int,
) -> Path:
    if source == "feature_feature":
        return (
            run_dir / "06_correlation" / "01_featurevsfeature_correlation"
            / method_dir / f"{k}_clusters"
            / f"{base_stub}_{k}_{method_code}_filtered_corr_pvals.csv"
        )
    if source == "feature_cluster":
        return (
            run_dir / "06_correlation" / "02_featurevscluster_correlation"
            / method_dir / f"{k}_clusters"
            / f"{base_stub}_{k}_{method_code}_filtered_corr_pvalues_long.csv"
        )
    if source == "cluster_cluster":
        return (
            run_dir / "06_correlation" / "03_clustervscluster_correlation"
            / method_dir / f"{k}_clusters"
            / f"{base_stub}_{k}_{method_code}_filtered_results.csv"
        )
    raise ValueError(f"Unknown source: {source}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="MetGuild Step 8: merge filtered correlation outputs into network edge/node files."
    )
    p.add_argument(
        "--run-dir",
        required=True,
        help="MetGuild run directory containing completed Steps 1 and 5-7.",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing 07_network_construction directory.",
    )
    p.add_argument(
        "--require-taxonomy",
        action="store_true",
        help=(
            "Require Step 3b taxonomy annotation and attach MAG_Phyla to cluster nodes. "
            "The runner enables this automatically when --taxonomy is supplied."
        ),
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    step1_path = run_dir / "manifests" / "optimal_k_values.json"
    fvf_manifest_path = run_dir / "06_correlation" / "01_featurevsfeature_correlation" / "01_fvf_correlation_summary.json"
    fvc_manifest_path = run_dir / "06_correlation" / "02_featurevscluster_correlation" / "01_fvc_correlation_summary.json"
    cvc_manifest_path = run_dir / "06_correlation" / "03_clustervscluster_correlation" / "01_cvc_correlation_summary.json"
    taxonomy_manifest_path = run_dir / "03_04_clustering" / "taxonomy" / "taxonomy_manifest.json"

    step1 = read_json(step1_path)
    fvf_manifest = read_json(fvf_manifest_path)
    fvc_manifest = read_json(fvc_manifest_path)
    cvc_manifest = read_json(cvc_manifest_path)
    taxonomy_manifest = read_json(taxonomy_manifest_path) if taxonomy_manifest_path.is_file() else None
    if args.require_taxonomy and taxonomy_manifest is None:
        raise FileNotFoundError(
            f"--require-taxonomy was requested but Step 3b manifest was not found: {taxonomy_manifest_path}"
        )
    if taxonomy_manifest is not None and taxonomy_manifest.get("status") != "success":
        raise ValueError("Step 3b taxonomy manifest does not report status='success'.")

    if step1.get("status") != "success":
        raise ValueError("Step 1 manifest does not report status='success'.")
    for name, manifest in [
        ("Step 5 feature-feature", fvf_manifest),
        ("Step 6 feature-cluster", fvc_manifest),
        ("Step 7 cluster-cluster", cvc_manifest),
    ]:
        if manifest.get("status") != "success":
            raise ValueError(f"{name} stage manifest does not report status='success'.")

    k_values = [int(k) for k in step1.get("unique_k_values", [])]
    if not k_values:
        raise ValueError("Step 1 manifest contains no unique_k_values.")

    base_stub = Path(step1.get("cleaned_output_file", "cleaned.tsv")).stem
    matrix_type = step1.get("preprocessing", {}).get("matrix_type")

    output_root = run_dir / "07_network_construction"
    if output_root.exists() and any(output_root.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Step 8 output already exists: {output_root}. Use --overwrite to replace it."
            )
        shutil.rmtree(output_root)
    edge_dir = output_root / "edge_lists"
    node_dir = output_root / "nodes_attributes"
    edge_dir.mkdir(parents=True, exist_ok=True)
    node_dir.mkdir(parents=True, exist_ok=True)

    stage_maps = {
        "feature_feature": stage_run_map(fvf_manifest),
        "feature_cluster": stage_run_map(fvc_manifest),
        "cluster_cluster": stage_run_map(cvc_manifest),
    }
    taxonomy_map = taxonomy_run_map(taxonomy_manifest) if taxonomy_manifest is not None else {}

    method_info = [
        ("KM", "01_kmeans", "kmeans"),
        ("HC", "02_hierarchical", "hierarchical"),
    ]

    summaries: list[dict] = []
    print(f"Run directory: {run_dir}")
    print(f"Matrix type: {matrix_type}")
    print("Correlation naming convention: features unchanged; clusters C1..Ck")
    print("Final network cluster naming: Cn_KMk or Cn_HCk")
    print(f"Candidate k values: {k_values}")
    print(f"Taxonomy annotation: {'enabled' if taxonomy_manifest is not None else 'not supplied'}")

    for k in k_values:
        for method_code, method_dir, method_name in method_info:
            print(f"\n=== Network merge: {method_name}, k={k} ===")

            taxonomy_lookup = None
            taxonomy_input = None
            if taxonomy_manifest is not None:
                tax_entry = taxonomy_map.get((method_code, k))
                if tax_entry is None or tax_entry.get("status") != "success":
                    raise ValueError(
                        f"Taxonomy manifest has no successful entry for {method_code}, k={k}."
                    )
                taxonomy_input = Path(tax_entry["output_file"])
                if not taxonomy_input.is_absolute():
                    taxonomy_input = (run_dir / taxonomy_input).resolve()
                taxonomy_lookup = load_cluster_taxonomy(taxonomy_input, k)

            frames: list[pd.DataFrame] = []
            inputs: dict[str, object] = {}
            source_counts: dict[str, int] = {}

            for source in ["feature_feature", "feature_cluster", "cluster_cluster"]:
                entry = stage_maps[source].get((method_code, k))
                if entry is None:
                    raise ValueError(f"{source}: no run-manifest entry for {method_code}, k={k}.")

                status = entry.get("status")
                if status == "skipped":
                    inputs[source] = {
                        "status": "skipped",
                        "reason": entry.get("reason"),
                    }
                    source_counts[source] = 0
                    print(f"  {source}: skipped upstream ({entry.get('reason')})")
                    continue
                if status != "success":
                    raise ValueError(f"{source}: unexpected upstream status {status!r} for {method_code}, k={k}.")

                path = source_file_for(
                    source=source,
                    run_dir=run_dir,
                    base_stub=base_stub,
                    method_dir=method_dir,
                    method_code=method_code,
                    k=k,
                )
                if source == "feature_feature":
                    frame = normalize_fvf(path)
                elif source == "feature_cluster":
                    frame = normalize_fvc(path, k)
                else:
                    frame = normalize_cvc(path, k)

                frames.append(frame)
                source_counts[source] = int(len(frame))
                inputs[source] = {
                    "status": "success",
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "edges": int(len(frame)),
                }
                print(f"  {source}: {len(frame)} filtered edge(s)")

            if frames:
                combined = pd.concat(frames, ignore_index=True)
            else:
                combined = pd.DataFrame(columns=[
                    "node1", "node2", "node1_type", "node2_type",
                    "correlation", "p_value", "adjusted_p_value", "edge_type",
                ])

            combined, removal = remove_invalid_edges(combined)
            combined = suffix_cluster_nodes(combined, method_code, k)
            assert_no_duplicate_undirected_edges(combined)

            # Preserve signed correlation and also provide a non-negative magnitude for MCL.
            combined["weight"] = combined["correlation"]
            combined["abs_weight"] = combined["correlation"].abs()
            combined["mcl_weight"] = combined["abs_weight"]

            final_edges = combined[[
                "node1", "node2", "weight", "abs_weight", "mcl_weight",
                "p_value", "adjusted_p_value", "edge_type",
            ]].copy()
            final_edges = final_edges.sort_values(
                ["edge_type", "node1", "node2"], kind="stable"
            ).reset_index(drop=True)

            prefix = f"{base_stub}_{k}_{method_code}"
            edge_path = edge_dir / f"{prefix}_edge_list.csv"
            node_path = node_dir / f"{prefix}_node_attributes.csv"
            final_edges.to_csv(edge_path, index=False)

            node_attr = build_node_attributes(combined, method_code, k, taxonomy_lookup=taxonomy_lookup)
            node_attr.to_csv(node_path, index=False)

            n_features = int((node_attr["Group"] == "feature").sum()) if len(node_attr) else 0
            n_clusters = int((node_attr["Group"] == "cluster").sum()) if len(node_attr) else 0

            result = {
                "timestamp_utc": utc_now(),
                "method_code": method_code,
                "method": method_name,
                "requested_k": int(k),
                "status": "success",
                "inputs": inputs,
                "source_edge_counts": source_counts,
                **removal,
                "final_edges": int(len(final_edges)),
                "nodes_total": int(len(node_attr)),
                "feature_nodes": n_features,
                "cluster_nodes": n_clusters,
                "cluster_name_format_upstream": "C1..Ck",
                "cluster_name_format_network": f"Cn_{method_code}{k}",
                "feature_names_modified": False,
                "taxonomy_annotation_used": taxonomy_lookup is not None,
                "taxonomy_input": str(taxonomy_input) if taxonomy_input is not None else None,
                "outputs": {
                    edge_path.name: {"path": str(edge_path), "sha256": sha256_file(edge_path)},
                    node_path.name: {"path": str(node_path), "sha256": sha256_file(node_path)},
                },
            }
            summaries.append(result)
            print(
                f"  Final network: {len(final_edges)} edges; "
                f"{len(node_attr)} nodes ({n_features} features, {n_clusters} clusters)"
            )

    summary_rows = []
    for r in summaries:
        summary_rows.append({
            "method": r["method_code"],
            "k": r["requested_k"],
            "final_edges": r["final_edges"],
            "nodes_total": r["nodes_total"],
            "feature_nodes": r["feature_nodes"],
            "cluster_nodes": r["cluster_nodes"],
            "feature_feature_edges": r["source_edge_counts"].get("feature_feature", 0),
            "feature_cluster_edges": r["source_edge_counts"].get("feature_cluster", 0),
            "cluster_cluster_edges": r["source_edge_counts"].get("cluster_cluster", 0),
            "removed_perfect_correlations": r["removed_perfect_correlations"],
            "removed_self_loops": r["removed_self_loops"],
        })

    summary_csv = output_root / "01_edge_list_merge_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

    manifest_path = output_root / "network_construction_manifest.json"
    manifest = {
        "manifest_schema_version": 1,
        "stage": "network_construction",
        "status": "success",
        "timestamp_utc": utc_now(),
        "matrix_type": matrix_type,
        "candidate_k_values": k_values,
        "naming_convention": {
            "feature_nodes": "exact original feature names; never modified",
            "correlation_cluster_nodes": "C1, C2, ..., Ck across Steps 6 and 7",
            "final_network_cluster_nodes": "Cn_KMk for KMeans or Cn_HCk for hierarchical clustering",
        },
        "perfect_correlations_excluded": True,
        "taxonomy_annotation": {
            "used": taxonomy_manifest is not None,
            "node_attribute_field": "MAG_Phyla",
            "description": "semicolon-separated MAG-Phylum pairs for each cluster node",
        },
        "edge_weights": {
            "weight": "signed Spearman correlation",
            "abs_weight": "absolute Spearman correlation",
            "mcl_weight": "non-negative absolute correlation magnitude for downstream MCL compatibility",
        },
        "runs": summaries,
        "provenance": {
            "step1_manifest": str(step1_path),
            "step1_manifest_sha256": sha256_file(step1_path),
            "step5_manifest": str(fvf_manifest_path),
            "step5_manifest_sha256": sha256_file(fvf_manifest_path),
            "step6_manifest": str(fvc_manifest_path),
            "step6_manifest_sha256": sha256_file(fvc_manifest_path),
            "step7_manifest": str(cvc_manifest_path),
            "step7_manifest_sha256": sha256_file(cvc_manifest_path),
            "taxonomy_manifest": str(taxonomy_manifest_path) if taxonomy_manifest is not None else None,
            "taxonomy_manifest_sha256": sha256_file(taxonomy_manifest_path) if taxonomy_manifest is not None else None,
            "software_versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
            },
            "script_sha256": sha256_file(Path(__file__)),
        },
        "outputs": {
            summary_csv.name: {"path": str(summary_csv), "sha256": sha256_file(summary_csv)},
        },
    }
    write_json(manifest_path, manifest)

    print("\nStep 8 network construction completed.")
    print(f"Edge lists:      {edge_dir}")
    print(f"Node attributes: {node_dir}")
    print(f"Summary:         {summary_csv}")
    print(f"Manifest:        {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
