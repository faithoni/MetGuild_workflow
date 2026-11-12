#!/usr/bin/env python3
"""MetGuild Step 6: feature-cluster correlation screening.

For each clustering method and candidate k, this stage correlates every feature
retained by Step 4 with one-vs-rest cluster membership indicators using
Spearman rank correlation. Pairwise p-values are adjusted with the
Benjamini-Hochberg procedure within each clustering solution.

By default, only positive associations are retained for network construction
(rho >= threshold and BH-adjusted p < alpha). Negative associations remain in
the full output and can be retained instead with --direction absolute.
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
import scipy
from scipy.stats import spearmanr

DEFAULT_CORR_THRESHOLD = 0.50
DEFAULT_ALPHA = 0.05
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


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, preserving input order."""
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out

    vals = p[finite]
    if np.any((vals < 0) | (vals > 1)):
        raise ValueError("P-values must lie in [0, 1].")

    m = len(vals)
    order = np.argsort(vals, kind="mergesort")
    ranked = vals[order]
    adjusted_ranked = ranked * m / np.arange(1, m + 1, dtype=float)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0.0, 1.0)

    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    out[finite] = adjusted
    return out


def load_retained_dataset(
    path: Path,
    *,
    expected_genomes: list[str],
    expected_features: list[str],
    matrix_type: str,
    expected_k: int,
) -> tuple[pd.DataFrame, list[str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Step 4 significant-feature dataset not found: {path}")

    df = pd.read_csv(path, dtype={"Genomes": "string"}, keep_default_na=False)
    required = {"Genomes", "Cluster"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")

    ids = df["Genomes"].astype(str).tolist()
    if len(ids) != len(set(ids)):
        raise ValueError(f"{path}: duplicate genome identifiers found.")
    if ids != expected_genomes:
        if set(ids) != set(expected_genomes):
            missing_ids = sorted(set(expected_genomes) - set(ids))
            extra_ids = sorted(set(ids) - set(expected_genomes))
            raise ValueError(
                f"{path}: genome identifiers do not match Step 1. "
                f"Missing={missing_ids[:10]}, extra={extra_ids[:10]}"
            )
        order = pd.Categorical(df["Genomes"].astype(str), categories=expected_genomes, ordered=True)
        df = (
            df.assign(_order=order)
            .sort_values("_order")
            .drop(columns="_order")
            .reset_index(drop=True)
        )

    clusters = pd.to_numeric(df["Cluster"], errors="coerce")
    if clusters.isna().any() or not np.allclose(clusters, np.round(clusters)):
        raise ValueError(f"{path}: cluster labels must be finite integers.")
    clusters = clusters.astype(int)
    observed = sorted(clusters.unique().tolist())
    expected_levels = list(range(1, expected_k + 1))
    if observed != expected_levels:
        raise ValueError(
            f"{path}: expected cluster labels {expected_levels}, observed {observed}."
        )
    df["Cluster"] = clusters

    feature_cols = [c for c in df.columns if c not in {"Genomes", "Cluster"}]
    if len(feature_cols) != len(set(feature_cols)):
        raise ValueError(f"{path}: duplicate feature names found.")
    if feature_cols != expected_features:
        if set(feature_cols) != set(expected_features):
            raise ValueError(
                f"{path}: feature set differs from Step 4 manifest. "
                f"CSV has {len(feature_cols)}, manifest has {len(expected_features)}."
            )
        df = df[["Genomes", "Cluster", *expected_features]]
        feature_cols = expected_features.copy()

    for col in feature_cols:
        converted = pd.to_numeric(df[col], errors="coerce")
        if converted.isna().any():
            raise ValueError(f"{path}: feature {col!r} contains missing/non-numeric values.")
        arr = converted.to_numpy(dtype=float)
        if not np.isfinite(arr).all():
            raise ValueError(f"{path}: feature {col!r} contains non-finite values.")
        if matrix_type == "binary" and not np.isin(arr, [0.0, 1.0]).all():
            raise ValueError(f"{path}: binary feature {col!r} contains values other than 0/1.")
        if matrix_type == "proportion" and (np.any(arr < 0.0) or np.any(arr > 1.0)):
            raise ValueError(f"{path}: proportional feature {col!r} contains values outside [0,1].")
        df[col] = converted.astype(float)

    return df, feature_cols


def cluster_label(level: int) -> str:
    return f"C{int(level)}"


def run_one(
    *,
    input_file: Path,
    out_dir: Path,
    out_prefix: str,
    expected_genomes: list[str],
    expected_features: list[str],
    matrix_type: str,
    method_code: str,
    k: int,
    corr_threshold: float,
    alpha: float,
    direction: str,
) -> dict:
    df, feature_cols = load_retained_dataset(
        input_file,
        expected_genomes=expected_genomes,
        expected_features=expected_features,
        matrix_type=matrix_type,
        expected_k=k,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    run_json = out_dir / f"{out_prefix}_run.json"

    if not feature_cols:
        payload = {
            "manifest_schema_version": 1,
            "stage": "feature_cluster_correlation",
            "status": "skipped",
            "timestamp_utc": utc_now(),
            "method_code": method_code,
            "requested_k": int(k),
            "matrix_type": matrix_type,
            "input_file": str(input_file),
            "input_sha256": sha256_file(input_file),
            "features_from_step4": 0,
            "reason": "Step 4 retained no features for this clustering solution.",
        }
        write_json(run_json, payload)
        return payload

    clusters = df["Cluster"].to_numpy(dtype=int)
    levels = list(range(1, k + 1))
    cluster_cols = [cluster_label(level) for level in levels]
    n = len(df)

    corr_mat = pd.DataFrame(np.nan, index=feature_cols, columns=cluster_cols, dtype=float)
    pval_mat = pd.DataFrame(np.nan, index=feature_cols, columns=cluster_cols, dtype=float)

    diagnostics_rows = []
    long_rows = []

    for level, cname in zip(levels, cluster_cols):
        indicator = (clusters == level).astype(float)
        cluster_size = int(indicator.sum())
        other_size = int(n - cluster_size)
        small_group = bool(cluster_size < 2 or other_size < 2)

        diagnostics_rows.append({
            "cluster": cname,
            "Cluster_Size": cluster_size,
            "Other_Size": other_size,
            "Small_Group_Flag": small_group,
        })

        for feature in feature_cols:
            y = df[feature].to_numpy(dtype=float)
            if np.nanstd(y) == 0.0:
                rho = np.nan
                p_value = np.nan
                status = "constant_feature"
            else:
                result = spearmanr(indicator, y, nan_policy="raise")
                rho = float(result.statistic)
                p_value = float(result.pvalue)
                status = "small_cluster_or_complement" if small_group else "ok"

            corr_mat.loc[feature, cname] = rho
            pval_mat.loc[feature, cname] = p_value
            long_rows.append({
                "feature": feature,
                "cluster": cname,
                "Correlation": rho,
                "P_value": p_value,
                "Status": status,
            })

    long_df = pd.DataFrame(long_rows)
    long_df["Adjusted_P_Value"] = bh_adjust(long_df["P_value"].to_numpy(dtype=float))
    long_df["Significant_FDR"] = long_df["Adjusted_P_Value"] < alpha
    long_df["Perfect_Correlation"] = np.isclose(
        long_df["Correlation"].abs(), 1.0, atol=PERFECT_ATOL, rtol=0.0
    )

    if direction == "positive":
        pass_corr = long_df["Correlation"] >= corr_threshold
    else:
        pass_corr = long_df["Correlation"].abs() >= corr_threshold

    long_df["Pass_Correlation_Threshold"] = pass_corr
    retained = (
        long_df["Pass_Correlation_Threshold"]
        & long_df["Significant_FDR"]
        & (long_df["Status"] != "constant_feature")
        & ~long_df["Perfect_Correlation"]
    )
    long_df["Retained_Edge"] = retained

    filtered = long_df.loc[
        long_df["Retained_Edge"],
        ["feature", "cluster", "Correlation", "P_value", "Adjusted_P_Value"],
    ].copy()

    # Keep the primary long table compact; technical flags live in diagnostics/run JSON.
    compact_long = long_df[[
        "feature", "cluster", "Correlation", "P_value",
        "Adjusted_P_Value", "Significant_FDR", "Retained_Edge",
    ]].copy()

    adj_mat = pd.DataFrame(np.nan, index=feature_cols, columns=cluster_cols, dtype=float)
    for row in long_df.itertuples(index=False):
        adj_mat.loc[row.feature, row.cluster] = row.Adjusted_P_Value

    diagnostics_df = pd.DataFrame(diagnostics_rows)

    out_full = out_dir / f"{out_prefix}_corr_pvalues_long.csv"
    out_filtered = out_dir / f"{out_prefix}_filtered_corr_pvalues_long.csv"
    out_corr = out_dir / f"{out_prefix}_corr_matrix.csv"
    out_p = out_dir / f"{out_prefix}_pvalue_matrix.csv"
    out_adj = out_dir / f"{out_prefix}_adjusted_pvalue_matrix.csv"
    out_diag = out_dir / f"{out_prefix}_cluster_diagnostics.csv"

    compact_long.to_csv(out_full, index=False)
    filtered.to_csv(out_filtered, index=False)
    corr_mat.reset_index(names="feature").to_csv(out_corr, index=False)
    pval_mat.reset_index(names="feature").to_csv(out_p, index=False)
    adj_mat.reset_index(names="feature").to_csv(out_adj, index=False)
    diagnostics_df.to_csv(out_diag, index=False)

    warning = None
    if n <= 500:
        warning = (
            "SciPy documents the asymptotic Spearman p-value approximation as most accurate for very large samples "
            "(>500 observations). Treat p-values/FDR as exploratory at this sample size; correlation coefficients "
            "are the primary network association weights."
        )

    payload = {
        "manifest_schema_version": 1,
        "stage": "feature_cluster_correlation",
        "status": "success",
        "timestamp_utc": utc_now(),
        "method_code": method_code,
        "requested_k": int(k),
        "matrix_type": matrix_type,
        "association_measure": "Spearman rank correlation between feature values and one-vs-rest cluster indicators",
        "cluster_encoding": "one_vs_rest_binary_indicator",
        "pvalue_method": "scipy_spearmanr_asymptotic_two_sided",
        "pvalue_warning": warning,
        "multiple_testing_correction": "Benjamini-Hochberg across all feature x cluster tests within this clustering solution",
        "alpha": float(alpha),
        "correlation_threshold": float(corr_threshold),
        "direction_filter": direction,
        "exclude_perfect_correlations": True,
        "perfect_correlation_tolerance": PERFECT_ATOL,
        "input_file": str(input_file),
        "input_sha256": sha256_file(input_file),
        "n_genomes": int(n),
        "features_from_step4": int(len(feature_cols)),
        "clusters": int(k),
        "feature_cluster_tests": int(len(long_df)),
        "fdr_significant_tests": int(long_df["Significant_FDR"].sum()),
        "tests_passing_correlation_threshold": int(long_df["Pass_Correlation_Threshold"].sum()),
        "perfect_correlation_tests": int(long_df["Perfect_Correlation"].sum()),
        "retained_edges": int(len(filtered)),
        "cluster_sizes": {
            row["cluster"]: int(row["Cluster_Size"])
            for row in diagnostics_rows
        },
        "interpretation": (
            "Exploratory feature-to-cluster association screening for network construction. Positive correlations "
            "identify features enriched in a cluster relative to the remaining genomes. Negative correlations are "
            "retained in the full table but are excluded from the default positive-edge network filter. Perfect "
            "correlations (|rho| = 1) are also always excluded from retained network edges. Because "
            "features and clusters derive from the same original matrix and Step 4 already selected cluster-discriminating "
            "features, these p-values are not independent confirmatory inference."
        ),
        "outputs": {
            out_full.name: {"path": str(out_full), "sha256": sha256_file(out_full)},
            out_filtered.name: {"path": str(out_filtered), "sha256": sha256_file(out_filtered)},
            out_corr.name: {"path": str(out_corr), "sha256": sha256_file(out_corr)},
            out_p.name: {"path": str(out_p), "sha256": sha256_file(out_p)},
            out_adj.name: {"path": str(out_adj), "sha256": sha256_file(out_adj)},
            out_diag.name: {"path": str(out_diag), "sha256": sha256_file(out_diag)},
        },
    }
    write_json(run_json, payload)
    return payload


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="MetGuild Step 6: feature-cluster Spearman correlations."
    )
    p.add_argument(
        "--run-dir",
        required=True,
        help="Completed MetGuild run directory containing Steps 1-4.",
    )
    p.add_argument(
        "--corr-threshold",
        type=float,
        default=DEFAULT_CORR_THRESHOLD,
        help=f"Minimum Spearman correlation magnitude/positive value for retained edges (default: {DEFAULT_CORR_THRESHOLD}).",
    )
    p.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help=f"BH-adjusted p-value cutoff (default: {DEFAULT_ALPHA}).",
    )
    p.add_argument(
        "--direction",
        choices=["positive", "absolute"],
        default="positive",
        help=(
            "Correlation filter for network edges: 'positive' keeps rho >= threshold; "
            "'absolute' keeps |rho| >= threshold. Default: positive."
        ),
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing Step 6 feature-cluster output directory.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not (0.0 <= args.corr_threshold <= 1.0):
        raise ValueError("--corr-threshold must be between 0 and 1.")
    if not (0.0 < args.alpha < 1.0):
        raise ValueError("--alpha must lie strictly between 0 and 1.")

    run_dir = Path(args.run_dir).expanduser().resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    step1_path = run_dir / "manifests" / "optimal_k_values.json"
    step4_path = run_dir / "05_stat_analysis" / "statistical_testing_manifest.json"
    step1 = read_json(step1_path)
    step4 = read_json(step4_path)

    if step1.get("status") != "success":
        raise ValueError("Step 1 manifest does not report status='success'.")
    if step4.get("status") != "success":
        raise ValueError("Step 4 manifest does not report status='success'.")

    matrix_type = step1.get("preprocessing", {}).get("matrix_type")
    if matrix_type not in {"binary", "proportion"}:
        raise ValueError(f"Unsupported Step 1 matrix type: {matrix_type!r}")
    if step4.get("matrix_type") != matrix_type:
        raise ValueError("Step 1 and Step 4 matrix types disagree.")

    selected_test = step4.get("selected_test")
    expected_test = "chi-square" if matrix_type == "binary" else "anova"
    if selected_test != expected_test:
        raise ValueError(
            f"Step 4 selected_test={selected_test!r}; expected {expected_test!r} for {matrix_type} data."
        )

    k_values = [int(k) for k in step1.get("unique_k_values", [])]
    if not k_values:
        raise ValueError("Step 1 manifest contains no unique_k_values.")
    step4_k = [int(k) for k in step4.get("tested_k_values", [])]
    if sorted(step4_k) != sorted(k_values):
        raise ValueError(
            f"Step 4 tested_k_values {step4_k} do not match Step 1 unique_k_values {k_values}."
        )

    expected_genomes = [str(x) for x in step1.get("genome_order", [])]
    if not expected_genomes:
        raise ValueError("Step 1 manifest contains no genome_order.")

    base_stub = Path(step1.get("cleaned_output_file", "cleaned.tsv")).stem
    output_root = run_dir / "06_correlation" / "02_featurevscluster_correlation"
    if output_root.exists() and any(output_root.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Step 6 output already exists: {output_root}. Use --overwrite to replace it."
            )
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    step4_runs = {
        (str(r.get("method_code")), int(r.get("requested_k"))): r
        for r in step4.get("runs", [])
        if r.get("status") == "success"
    }

    suffix = "chi_stat_sig.csv" if selected_test == "chi-square" else "anova_stat_significant.csv"
    method_info = [
        ("KM", "01_kmeans", "kmeans"),
        ("HC", "02_hierarchical", "hierarchical"),
    ]

    summaries = []
    print(f"Run directory: {run_dir}")
    print(f"Matrix type: {matrix_type}")
    print("Association: Spearman rank correlation against one-vs-rest cluster indicators")
    if args.direction == "positive":
        print(f"Thresholds: rho >= {args.corr_threshold}; BH-adjusted p < {args.alpha}")
    else:
        print(f"Thresholds: |rho| >= {args.corr_threshold}; BH-adjusted p < {args.alpha}")
    print(f"Candidate k values: {k_values}")

    for k in k_values:
        print(f"\n=== Feature-cluster correlations for k={k} ===")
        for code, method_dir, method_name in method_info:
            entry = step4_runs.get((code, k))
            if entry is None:
                raise ValueError(f"Step 4 manifest has no successful run for {code}, k={k}.")

            stat_dir = run_dir / "05_stat_analysis" / method_dir / f"{k}_clusters"
            input_file = stat_dir / f"05_{base_stub}_{k}_{code}_{suffix}"
            out_dir = output_root / method_dir / f"{k}_clusters"
            out_prefix = f"{base_stub}_{k}_{code}"

            print(f"[{method_name} k={k}] Input: {input_file}")
            result = run_one(
                input_file=input_file,
                out_dir=out_dir,
                out_prefix=out_prefix,
                expected_genomes=expected_genomes,
                expected_features=[str(x) for x in entry.get("significant_features", [])],
                matrix_type=matrix_type,
                method_code=code,
                k=k,
                corr_threshold=args.corr_threshold,
                alpha=args.alpha,
                direction=args.direction,
            )
            summaries.append(result)

            if result["status"] == "success":
                print(f"  Features: {result['features_from_step4']}")
                print(f"  Clusters: {result['clusters']}")
                print(f"  Feature-cluster tests: {result['feature_cluster_tests']}")
                print(f"  Retained edges: {result['retained_edges']}")
                if result.get("pvalue_warning"):
                    print("  WARNING: asymptotic Spearman p-values are exploratory at this sample size.")
            else:
                print(f"  SKIPPED: {result['reason']}")

    summary_rows = []
    for r in summaries:
        summary_rows.append({
            "timestamp_utc": r.get("timestamp_utc"),
            "method_code": r.get("method_code"),
            "k": r.get("requested_k"),
            "status": r.get("status"),
            "features_from_step4": r.get("features_from_step4", 0),
            "clusters": r.get("clusters", 0),
            "feature_cluster_tests": r.get("feature_cluster_tests", 0),
            "fdr_significant_tests": r.get("fdr_significant_tests", 0),
            "tests_passing_correlation_threshold": r.get("tests_passing_correlation_threshold", 0),
            "perfect_correlation_tests": r.get("perfect_correlation_tests", 0),
            "retained_edges": r.get("retained_edges", 0),
            "reason": r.get("reason"),
        })

    summary_csv = output_root / "01_fvc_correlation_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

    stage_manifest = {
        "manifest_schema_version": 1,
        "stage": "feature_cluster_correlation",
        "status": "success",
        "timestamp_utc": utc_now(),
        "matrix_type": matrix_type,
        "step4_selected_test": selected_test,
        "candidate_k_values": k_values,
        "association_measure": "Spearman rank correlation against one-vs-rest cluster indicators",
        "cluster_encoding": "one_vs_rest_binary_indicator",
        "pvalue_method": "scipy_spearmanr_asymptotic_two_sided",
        "multiple_testing_correction": "Benjamini-Hochberg across feature x cluster tests within each clustering method x k solution",
        "alpha": float(args.alpha),
        "correlation_threshold": float(args.corr_threshold),
        "direction_filter": args.direction,
        "exclude_perfect_correlations": True,
        "provenance": {
            "step1_manifest": str(step1_path),
            "step1_manifest_sha256": sha256_file(step1_path),
            "step4_manifest": str(step4_path),
            "step4_manifest_sha256": sha256_file(step4_path),
            "software_versions": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
            },
        },
        "interpretation": (
            "Exploratory feature-cluster association screening for network construction. The default positive filter "
            "retains features enriched in a cluster; negative correlations remain available in the full table. "
            "FDR correction addresses multiplicity within each clustering solution but not selection induced by prior "
            "clustering and Step 4 feature screening."
        ),
        "runs": summaries,
    }
    summary_json = output_root / "01_fvc_correlation_summary.json"
    write_json(summary_json, stage_manifest)

    print("\nStep 6 feature-cluster correlation completed.")
    print(f"Summary CSV:  {summary_csv}")
    print(f"Summary JSON: {summary_json}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
