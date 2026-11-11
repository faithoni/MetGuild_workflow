#!/usr/bin/env python3
"""MetGuild Step 7: cluster-cluster correlation screening.

For each clustering method and candidate k, this stage builds one metabolic
profile per cluster from the Step 4 retained features and compares cluster
profiles using Spearman rank correlation.

Pairwise p-values are adjusted with Benjamini-Hochberg FDR within each
clustering solution. By default, positive cluster-cluster associations are
retained when rho >= threshold and BH-adjusted p < alpha. Perfect correlations
(|rho| = 1) are always excluded from retained network edges, while remaining
visible in the full results for traceability.
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
MIN_PROFILE_FEATURES = 3


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
        value = json.load(fh)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


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

    df = pd.read_csv(path, dtype={"Genomes": "string"})
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


def cluster_name(level: int) -> str:
    return f"C{int(level)}"


def build_cluster_profiles(
    df: pd.DataFrame,
    feature_cols: list[str],
    *,
    matrix_type: str,
    k: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Return feature x cluster profile matrix and diagnostics.

    A cluster profile is the cluster-wise mean for every retained feature.
    For binary data this is feature prevalence; for proportional data it is
    mean feature proportion. Features with identical profile values across all
    clusters are removed from correlation calculations but retained in the
    exported full profile table.
    """
    levels = list(range(1, k + 1))
    cluster_sizes = df.groupby("Cluster", sort=True).size().reindex(levels)
    if cluster_sizes.isna().any() or (cluster_sizes <= 0).any():
        raise ValueError("One or more requested clusters contain no genomes.")

    if not feature_cols:
        empty = pd.DataFrame(index=pd.Index([], name="feature"), columns=[cluster_name(x) for x in levels])
        diagnostics = pd.DataFrame({
            "cluster": [cluster_name(x) for x in levels],
            "Cluster_Size": [int(cluster_sizes.loc[x]) for x in levels],
            "Profile_Aggregation": ["prevalence" if matrix_type == "binary" else "mean_proportion"] * k,
        })
        return empty, diagnostics, []

    means = df.groupby("Cluster", sort=True)[feature_cols].mean().reindex(levels)
    profiles = means.T
    profiles.columns = [cluster_name(x) for x in levels]
    profiles.index.name = "feature"

    if matrix_type == "binary":
        aggregation = "feature_prevalence_within_cluster"
    else:
        aggregation = "mean_proportional_feature_value_within_cluster"

    diagnostics_rows = []
    for level, cname in zip(levels, profiles.columns):
        vals = profiles[cname].to_numpy(dtype=float)
        diagnostics_rows.append({
            "cluster": cname,
            "Cluster_Size": int(cluster_sizes.loc[level]),
            "Profile_Aggregation": aggregation,
            "Profile_Min": float(np.min(vals)) if len(vals) else np.nan,
            "Profile_Max": float(np.max(vals)) if len(vals) else np.nan,
            "Profile_Mean": float(np.mean(vals)) if len(vals) else np.nan,
            "Nonzero_Profile_Features": int(np.count_nonzero(vals)) if len(vals) else 0,
            "Constant_Cluster_Profile": bool(len(vals) > 0 and np.allclose(vals, vals[0])),
        })

    profile_variation = profiles.max(axis=1) - profiles.min(axis=1)
    invariant_features = profiles.index[np.isclose(profile_variation.to_numpy(dtype=float), 0.0)].tolist()
    diagnostics = pd.DataFrame(diagnostics_rows)
    return profiles, diagnostics, invariant_features


def spearman_cluster_matrices(values: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute Spearman correlation/p-value matrices for cluster profile columns."""
    cols = values.columns.tolist()
    m = len(cols)
    if m < 2:
        raise ValueError("Need at least two non-constant cluster profiles.")

    x = values.to_numpy(dtype=float)
    if m == 2:
        result = spearmanr(x[:, 0], x[:, 1], nan_policy="raise")
        rho = float(result.statistic)
        p = float(result.pvalue)
        corr = np.array([[1.0, rho], [rho, 1.0]], dtype=float)
        pvals = np.array([[0.0, p], [p, 0.0]], dtype=float)
    else:
        result = spearmanr(x, axis=0, nan_policy="raise")
        corr = np.asarray(result.statistic, dtype=float)
        pvals = np.asarray(result.pvalue, dtype=float)
        if corr.shape != (m, m) or pvals.shape != (m, m):
            raise RuntimeError(
                f"Unexpected scipy.stats.spearmanr output shapes: corr={corr.shape}, "
                f"p={pvals.shape}, clusters={m}"
            )

    corr_df = pd.DataFrame(corr, index=cols, columns=cols)
    pval_df = pd.DataFrame(pvals, index=cols, columns=cols)
    np.fill_diagonal(pval_df.values, np.nan)
    return corr_df, pval_df


def flatten_pairs(corr_df: pd.DataFrame, pval_df: pd.DataFrame) -> pd.DataFrame:
    cols = np.asarray(corr_df.columns, dtype=object)
    ii, jj = np.triu_indices(len(cols), k=1)
    corr = corr_df.to_numpy(dtype=float)[ii, jj]
    pvals = pval_df.to_numpy(dtype=float)[ii, jj]
    out = pd.DataFrame({
        "Var1": cols[ii],
        "Var2": cols[jj],
        "Correlation": corr,
        "P_value": pvals,
    })
    out["Adjusted_P_Value"] = bh_adjust(out["P_value"].to_numpy(dtype=float))
    out["Significant_FDR"] = out["Adjusted_P_Value"] < DEFAULT_ALPHA  # overwritten by caller
    out["Perfect_Correlation"] = np.isclose(
        out["Correlation"].abs(), 1.0, atol=PERFECT_ATOL, rtol=0.0
    )
    return out


def adjusted_p_matrix(flat: pd.DataFrame, clusters: list[str]) -> pd.DataFrame:
    mat = pd.DataFrame(np.nan, index=clusters, columns=clusters, dtype=float)
    np.fill_diagonal(mat.values, 0.0)
    for row in flat.itertuples(index=False):
        c1 = str(row.Var1)
        c2 = str(row.Var2)
        mat.loc[c1, c2] = row.Adjusted_P_Value
        mat.loc[c2, c1] = row.Adjusted_P_Value
    return mat


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

    profiles, diagnostics_df, invariant_features = build_cluster_profiles(
        df,
        feature_cols,
        matrix_type=matrix_type,
        k=k,
    )

    out_profiles = out_dir / f"{out_prefix}_cluster_feature_profiles.csv"
    profiles.reset_index().to_csv(out_profiles, index=False)

    out_diag = out_dir / f"{out_prefix}_cluster_profile_diagnostics.csv"
    diagnostics_df.to_csv(out_diag, index=False)

    usable_features = [f for f in profiles.index.tolist() if f not in set(invariant_features)]
    usable = profiles.loc[usable_features].copy()

    constant_cluster_profiles = [
        c for c in usable.columns
        if usable.shape[0] == 0 or usable[c].nunique(dropna=False) < 2
    ]
    usable_cluster_cols = [c for c in usable.columns if c not in constant_cluster_profiles]
    usable = usable[usable_cluster_cols]

    base_payload = {
        "manifest_schema_version": 1,
        "stage": "cluster_cluster_correlation",
        "timestamp_utc": utc_now(),
        "method_code": method_code,
        "requested_k": int(k),
        "matrix_type": matrix_type,
        "input_file": str(input_file),
        "input_sha256": sha256_file(input_file),
        "profile_aggregation": (
            "feature prevalence within cluster (cluster-wise mean of 0/1 values)"
            if matrix_type == "binary"
            else "mean proportional feature value within cluster"
        ),
        "features_from_step4": int(len(feature_cols)),
        "invariant_profile_features_removed": invariant_features,
        "features_used_for_cluster_correlations": int(len(usable_features)),
        "constant_cluster_profiles_removed": constant_cluster_profiles,
        "clusters_total": int(k),
        "clusters_used": int(len(usable_cluster_cols)),
        "exclude_perfect_correlations": True,
        "perfect_correlation_tolerance": PERFECT_ATOL,
        "outputs": {
            out_profiles.name: {"path": str(out_profiles), "sha256": sha256_file(out_profiles)},
            out_diag.name: {"path": str(out_diag), "sha256": sha256_file(out_diag)},
        },
    }

    if len(feature_cols) == 0:
        payload = {
            **base_payload,
            "status": "skipped",
            "reason": "Step 4 retained no features for this clustering solution.",
        }
        write_json(run_json, payload)
        return payload

    if len(usable_features) < MIN_PROFILE_FEATURES:
        payload = {
            **base_payload,
            "status": "skipped",
            "reason": (
                f"Need at least {MIN_PROFILE_FEATURES} varying retained features to compare cluster profiles; "
                f"only {len(usable_features)} remain after removing invariant profile features."
            ),
        }
        write_json(run_json, payload)
        return payload

    if len(usable_cluster_cols) < 2:
        payload = {
            **base_payload,
            "status": "skipped",
            "reason": "Fewer than two non-constant cluster profiles remain for correlation.",
        }
        write_json(run_json, payload)
        return payload

    corr_df, pval_df = spearman_cluster_matrices(usable)
    flat = flatten_pairs(corr_df, pval_df)
    flat["Significant_FDR"] = flat["Adjusted_P_Value"] < alpha
    flat["Pass_Correlation_Threshold"] = flat["Correlation"] >= corr_threshold
    flat["Retained_Edge"] = (
        flat["Pass_Correlation_Threshold"]
        & flat["Significant_FDR"]
        & ~flat["Perfect_Correlation"]
    )

    filtered = flat.loc[
        flat["Retained_Edge"],
        ["Var1", "Var2", "Correlation", "P_value", "Adjusted_P_Value"],
    ].copy()

    compact_full = flat[[
        "Var1", "Var2", "Correlation", "P_value",
        "Adjusted_P_Value", "Significant_FDR", "Retained_Edge",
    ]].copy()

    adj_df = adjusted_p_matrix(flat, corr_df.columns.tolist())

    out_corr = out_dir / f"{out_prefix}_corr_matrix.csv"
    out_p = out_dir / f"{out_prefix}_p_value_matrix.csv"
    out_adj = out_dir / f"{out_prefix}_adjusted_p_value_matrix.csv"
    out_flat = out_dir / f"{out_prefix}_flattened_results.csv"
    out_filt = out_dir / f"{out_prefix}_filtered_results.csv"

    corr_df.to_csv(out_corr, index=True)
    pval_df.to_csv(out_p, index=True)
    adj_df.to_csv(out_adj, index=True)
    compact_full.to_csv(out_flat, index=False)
    filtered.to_csv(out_filt, index=False)

    n_profile_features = len(usable_features)
    pvalue_warning = None
    if n_profile_features <= 500:
        pvalue_warning = (
            "SciPy documents the asymptotic Spearman p-value approximation as most accurate for very large samples "
            "(>500 observations). Here the observations are retained profile features, so p-values/FDR should be "
            "treated as exploratory; correlation coefficients are the primary network weights."
        )

    output_records = dict(base_payload["outputs"])
    for p in [out_corr, out_p, out_adj, out_flat, out_filt]:
        output_records[p.name] = {"path": str(p), "sha256": sha256_file(p)}

    payload = {
        **base_payload,
        "status": "success",
        "association_measure": "Spearman rank correlation between cluster-level retained-feature profiles",
        "pvalue_method": "scipy_spearmanr_asymptotic_two_sided",
        "pvalue_warning": pvalue_warning,
        "multiple_testing_correction": "Benjamini-Hochberg across unique cluster pairs within this clustering solution",
        "alpha": float(alpha),
        "correlation_threshold": float(corr_threshold),
        "direction_filter": "positive",
        "cluster_pairs_tested": int(len(flat)),
        "fdr_significant_pairs": int(flat["Significant_FDR"].sum()),
        "pairs_passing_correlation_threshold": int(flat["Pass_Correlation_Threshold"].sum()),
        "perfect_correlation_pairs": int(flat["Perfect_Correlation"].sum()),
        "retained_edges": int(len(filtered)),
        "interpretation": (
            "Exploratory similarity screening between cluster-level metabolic profiles. Positive correlations indicate "
            "clusters whose retained functional features show similar relative patterns. Perfect correlations are always "
            "excluded from retained network edges. The retained features were selected using the same cluster solution in "
            "Step 4, and functional features are not independent observations, so p-values/FDR are not independent "
            "confirmatory inference."
        ),
        "outputs": output_records,
    }
    write_json(run_json, payload)
    return payload


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="MetGuild Step 7: cluster-cluster Spearman correlations."
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
        help=f"Minimum positive Spearman rho for retained cluster-cluster edges (default: {DEFAULT_CORR_THRESHOLD}).",
    )
    p.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help=f"BH-adjusted p-value cutoff (default: {DEFAULT_ALPHA}).",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing Step 7 cluster-cluster output directory.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not (0.0 <= args.corr_threshold < 1.0):
        raise ValueError("--corr-threshold must be >=0 and <1 because perfect correlations are excluded.")
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
    output_root = run_dir / "06_correlation" / "03_clustervscluster_correlation"
    if output_root.exists() and any(output_root.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Step 7 output already exists: {output_root}. Use --overwrite to replace it."
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
    if matrix_type == "binary":
        print("Cluster profiles: feature prevalence within each cluster")
    else:
        print("Cluster profiles: mean proportional feature value within each cluster")
    print("Association: Spearman rank correlation between cluster profiles")
    print(f"Thresholds: rho >= {args.corr_threshold}; BH-adjusted p < {args.alpha}; |rho|=1 excluded")
    print(f"Candidate k values: {k_values}")

    for k in k_values:
        print(f"\n=== Cluster-cluster correlations for k={k} ===")
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
            )
            summaries.append(result)

            if result["status"] == "success":
                print(f"  Step 4 features:       {result['features_from_step4']}")
                print(f"  Profile features used: {result['features_used_for_cluster_correlations']}")
                print(f"  Clusters used:         {result['clusters_used']}")
                print(f"  Cluster pairs tested:  {result['cluster_pairs_tested']}")
                print(f"  Perfect pairs removed: {result['perfect_correlation_pairs']}")
                print(f"  Retained edges:        {result['retained_edges']}")
            else:
                print(f"  SKIPPED: {result['reason']}")

    summary_rows = []
    for r in summaries:
        summary_rows.append({
            "method": r.get("method_code"),
            "k": r.get("requested_k"),
            "matrix_type": r.get("matrix_type"),
            "status": r.get("status"),
            "features_from_step4": r.get("features_from_step4"),
            "features_used": r.get("features_used_for_cluster_correlations"),
            "clusters_used": r.get("clusters_used"),
            "cluster_pairs_tested": r.get("cluster_pairs_tested", 0),
            "retained_edges": r.get("retained_edges", 0),
            "reason": r.get("reason"),
        })

    summary_csv = output_root / "01_cvc_correlation_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

    summary_json = output_root / "01_cvc_correlation_summary.json"
    stage_payload = {
        "manifest_schema_version": 1,
        "stage": "cluster_cluster_correlation",
        "status": "success",
        "timestamp_utc": utc_now(),
        "matrix_type": matrix_type,
        "profile_aggregation": (
            "feature prevalence within cluster" if matrix_type == "binary"
            else "mean proportional feature value within cluster"
        ),
        "association_measure": "Spearman rank correlation",
        "correlation_threshold": float(args.corr_threshold),
        "alpha": float(args.alpha),
        "multiple_testing_correction": "Benjamini-Hochberg within each method x k solution",
        "exclude_perfect_correlations": True,
        "candidate_k_values": k_values,
        "runs": summaries,
        "software_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "provenance": {
            "step1_manifest": str(step1_path),
            "step1_manifest_sha256": sha256_file(step1_path),
            "step4_manifest": str(step4_path),
            "step4_manifest_sha256": sha256_file(step4_path),
            "script_sha256": sha256_file(Path(__file__)),
        },
        "outputs": {
            summary_csv.name: {"path": str(summary_csv), "sha256": sha256_file(summary_csv)},
        },
    }
    write_json(summary_json, stage_payload)

    print("\nStep 7 cluster-cluster correlation completed.")
    print(f"Summary CSV:  {summary_csv}")
    print(f"Summary JSON: {summary_json}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
