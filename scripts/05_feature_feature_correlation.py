#!/usr/bin/env python3
"""MetGuild Step 5: feature-feature correlation screening.

Takes the statistically retained feature datasets from Step 4 and computes
Spearman feature-feature correlations separately for each clustering method and
candidate k. Pairwise p-values are adjusted with Benjamini-Hochberg FDR within
each clustering solution. Perfect correlations (|rho| = 1) are always excluded
from the retained network-edge table, while remaining visible in the full results.
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

DEFAULT_ABS_CORR_THRESHOLD = 0.50
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


def load_retained_dataset(path: Path, expected_genomes: list[str]) -> tuple[pd.DataFrame, list[str]]:
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
        df = df.assign(_order=order).sort_values("_order").drop(columns="_order").reset_index(drop=True)

    feature_cols = [c for c in df.columns if c not in {"Genomes", "Cluster"}]
    if len(feature_cols) != len(set(feature_cols)):
        raise ValueError(f"{path}: duplicate feature names found.")

    for col in feature_cols:
        converted = pd.to_numeric(df[col], errors="coerce")
        if converted.isna().any():
            bad_n = int(converted.isna().sum())
            raise ValueError(f"{path}: feature {col!r} contains {bad_n} missing/non-numeric value(s).")
        arr = converted.to_numpy(dtype=float)
        if not np.isfinite(arr).all():
            raise ValueError(f"{path}: feature {col!r} contains non-finite values.")
        df[col] = converted.astype(float)

    return df, feature_cols


def spearman_matrices(values: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = values.columns.tolist()
    m = len(cols)
    if m < 2:
        raise ValueError("Need at least two non-constant features for feature-feature correlation.")

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
                f"Unexpected scipy.stats.spearmanr output shapes: corr={corr.shape}, p={pvals.shape}, features={m}"
            )

    return (
        pd.DataFrame(corr, index=cols, columns=cols),
        pd.DataFrame(pvals, index=cols, columns=cols),
    )


def flatten_pairs(corr_df: pd.DataFrame, pval_df: pd.DataFrame, n_genomes: int) -> pd.DataFrame:
    cols = np.asarray(corr_df.columns, dtype=object)
    ii, jj = np.triu_indices(len(cols), k=1)
    corr = corr_df.to_numpy(dtype=float)[ii, jj]
    pvals = pval_df.to_numpy(dtype=float)[ii, jj]
    out = pd.DataFrame(
        {
            "Var1": cols[ii],
            "Var2": cols[jj],
            "N": int(n_genomes),
            "Correlation": corr,
            "Abs_Correlation": np.abs(corr),
            "P_value": pvals,
        }
    )
    out["Adjusted_P_Value"] = bh_adjust(out["P_value"].to_numpy(dtype=float))
    out["Direction"] = np.where(out["Correlation"] > 0, "positive", np.where(out["Correlation"] < 0, "negative", "zero"))
    out["Perfect_Correlation"] = np.isclose(out["Abs_Correlation"], 1.0, atol=PERFECT_ATOL, rtol=0.0)
    return out


def adjusted_p_matrix(flat: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    mat = pd.DataFrame(np.nan, index=features, columns=features, dtype=float)
    np.fill_diagonal(mat.values, 0.0)
    for row in flat.itertuples(index=False):
        mat.loc[row.Var1, row.Var2] = row.Adjusted_P_Value
        mat.loc[row.Var2, row.Var1] = row.Adjusted_P_Value
    return mat


def run_one(
    *,
    input_file: Path,
    out_dir: Path,
    out_prefix: str,
    expected_genomes: list[str],
    expected_features: list[str],
    abs_threshold: float,
    alpha: float,
    method_code: str,
    k: int,
    matrix_type: str,
) -> dict:
    df, feature_cols = load_retained_dataset(input_file, expected_genomes)

    if feature_cols != expected_features:
        if set(feature_cols) != set(expected_features):
            raise ValueError(
                f"{input_file}: feature set differs from Step 4 manifest. "
                f"CSV has {len(feature_cols)}, manifest has {len(expected_features)}."
            )
        df = df[["Genomes", "Cluster"] + expected_features]
        feature_cols = expected_features.copy()

    constant_features = [c for c in feature_cols if df[c].nunique(dropna=False) < 2]
    usable_features = [c for c in feature_cols if c not in constant_features]

    out_dir.mkdir(parents=True, exist_ok=True)
    run_json = out_dir / f"{out_prefix}_run.json"

    if len(usable_features) < 2:
        payload = {
            "manifest_schema_version": 1,
            "stage": "feature_feature_correlation",
            "status": "skipped",
            "timestamp_utc": utc_now(),
            "method_code": method_code,
            "requested_k": int(k),
            "matrix_type": matrix_type,
            "input_file": str(input_file),
            "input_sha256": sha256_file(input_file),
            "features_from_step4": len(feature_cols),
            "constant_features_removed": constant_features,
            "usable_features": len(usable_features),
            "reason": "Fewer than two non-constant statistically retained features.",
        }
        write_json(run_json, payload)
        return payload

    values = df[usable_features]
    corr_df, pval_df = spearman_matrices(values)
    flat = flatten_pairs(corr_df, pval_df, len(df))

    flat["Pass_Correlation_Threshold"] = flat["Abs_Correlation"] >= abs_threshold
    flat["Significant_FDR"] = flat["Adjusted_P_Value"] < alpha
    retained = (
        flat["Pass_Correlation_Threshold"]
        & flat["Significant_FDR"]
        & ~flat["Perfect_Correlation"]
    )
    flat["Retained_Edge"] = retained
    filtered = flat[flat["Retained_Edge"]].copy()

    adj_df = adjusted_p_matrix(flat, usable_features)

    out_corr = out_dir / f"{out_prefix}_corr_matrix.csv"
    out_p = out_dir / f"{out_prefix}_pvalue_matrix.csv"
    out_adj = out_dir / f"{out_prefix}_adjusted_pvalue_matrix.csv"
    out_flat = out_dir / f"{out_prefix}_flattened_corr_pvals.csv"
    out_filtered = out_dir / f"{out_prefix}_filtered_corr_pvals.csv"

    corr_df.to_csv(out_corr, index=True)
    pval_df.to_csv(out_p, index=True)
    adj_df.to_csv(out_adj, index=True)
    flat.to_csv(out_flat, index=False)
    filtered.to_csv(out_filtered, index=False)

    warning = None
    if len(df) <= 500:
        warning = (
            "SciPy's asymptotic Spearman p-value approximation is documented as most accurate for very large samples "
            "(>500 observations). Treat p-values/FDR as exploratory at this sample size; the correlation coefficients "
            "remain the primary network weights."
        )

    payload = {
        "manifest_schema_version": 1,
        "stage": "feature_feature_correlation",
        "status": "success",
        "timestamp_utc": utc_now(),
        "method_code": method_code,
        "requested_k": int(k),
        "matrix_type": matrix_type,
        "association_measure": "Spearman rank correlation",
        "pvalue_method": "scipy_spearmanr_asymptotic_two_sided",
        "pvalue_warning": warning,
        "multiple_testing_correction": "Benjamini-Hochberg",
        "alpha": float(alpha),
        "absolute_correlation_threshold": float(abs_threshold),
        "exclude_perfect_correlations": True,
        "perfect_correlation_tolerance": PERFECT_ATOL,
        "input_file": str(input_file),
        "input_sha256": sha256_file(input_file),
        "n_genomes": int(len(df)),
        "features_from_step4": int(len(feature_cols)),
        "constant_features_removed": constant_features,
        "features_used": int(len(usable_features)),
        "all_unique_pairs": int(len(flat)),
        "fdr_significant_pairs": int(flat["Significant_FDR"].sum()),
        "pairs_passing_abs_correlation": int(flat["Pass_Correlation_Threshold"].sum()),
        "perfect_correlation_pairs": int(flat["Perfect_Correlation"].sum()),
        "filtered_pairs": int(len(filtered)),
        "interpretation": (
            "Exploratory association/network screening among Step 4 retained features. FDR correction addresses "
            "multiplicity of pairwise tests but does not make the analysis independent of the preceding data-driven "
            "feature selection and clustering."
        ),
        "outputs": {
            out_corr.name: {"path": str(out_corr), "sha256": sha256_file(out_corr)},
            out_p.name: {"path": str(out_p), "sha256": sha256_file(out_p)},
            out_adj.name: {"path": str(out_adj), "sha256": sha256_file(out_adj)},
            out_flat.name: {"path": str(out_flat), "sha256": sha256_file(out_flat)},
            out_filtered.name: {"path": str(out_filtered), "sha256": sha256_file(out_filtered)},
        },
    }
    write_json(run_json, payload)
    return payload


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="MetGuild Step 5: feature-feature Spearman correlations.")
    p.add_argument("--run-dir", required=True, help="Completed MetGuild run directory containing Steps 1-4.")
    p.add_argument("--abs-corr-threshold", type=float, default=DEFAULT_ABS_CORR_THRESHOLD,
                   help=f"Minimum absolute Spearman correlation for retained edges (default: {DEFAULT_ABS_CORR_THRESHOLD}).")
    p.add_argument("--alpha", type=float, default=DEFAULT_ALPHA,
                   help=f"BH-adjusted p-value cutoff (default: {DEFAULT_ALPHA}).")
    p.add_argument("--overwrite", action="store_true", help="Replace an existing Step 5 feature-feature output directory.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not (0.0 <= args.abs_corr_threshold <= 1.0):
        raise ValueError("--abs-corr-threshold must be between 0 and 1.")
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
        raise ValueError(f"Step 4 selected_test={selected_test!r}; expected {expected_test!r} for {matrix_type} data.")

    k_values = [int(k) for k in step1.get("unique_k_values", [])]
    if not k_values:
        raise ValueError("Step 1 manifest contains no unique_k_values.")
    step4_k = [int(k) for k in step4.get("tested_k_values", [])]
    if sorted(step4_k) != sorted(k_values):
        raise ValueError(f"Step 4 tested_k_values {step4_k} do not match Step 1 unique_k_values {k_values}.")

    expected_genomes = [str(x) for x in step1.get("genome_order", [])]
    if not expected_genomes:
        raise ValueError("Step 1 manifest contains no genome_order.")

    base_stub = Path(step1.get("cleaned_output_file", "cleaned.tsv")).stem
    output_root = run_dir / "06_correlation" / "01_featurevsfeature_correlation"
    if output_root.exists() and any(output_root.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Step 5 output already exists: {output_root}. Use --overwrite to replace it."
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
    print("Association: Spearman rank correlation")
    print(f"Thresholds: |rho| >= {args.abs_corr_threshold}; BH-adjusted p < {args.alpha}")
    print(f"Candidate k values: {k_values}")

    for k in k_values:
        print(f"\n=== Feature-feature correlations for k={k} ===")
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
                abs_threshold=args.abs_corr_threshold,
                alpha=args.alpha,
                method_code=code,
                k=k,
                matrix_type=matrix_type,
            )
            summaries.append(result)
            if result["status"] == "success":
                print(f"  Features used: {result['features_used']}")
                print(f"  Unique pairs:  {result['all_unique_pairs']}")
                print(f"  Retained edges: {result['filtered_pairs']}")
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
            "features_from_step4": r.get("features_from_step4"),
            "features_used": r.get("features_used", r.get("usable_features")),
            "all_unique_pairs": r.get("all_unique_pairs", 0),
            "fdr_significant_pairs": r.get("fdr_significant_pairs", 0),
            "pairs_passing_abs_correlation": r.get("pairs_passing_abs_correlation", 0),
            "perfect_correlation_pairs": r.get("perfect_correlation_pairs", 0),
            "filtered_pairs": r.get("filtered_pairs", 0),
            "reason": r.get("reason"),
        })

    summary_csv = output_root / "01_fvf_correlation_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)

    stage_manifest = {
        "manifest_schema_version": 1,
        "stage": "feature_feature_correlation",
        "status": "success",
        "timestamp_utc": utc_now(),
        "matrix_type": matrix_type,
        "step4_selected_test": selected_test,
        "candidate_k_values": k_values,
        "association_measure": "Spearman rank correlation",
        "pvalue_method": "scipy_spearmanr_asymptotic_two_sided",
        "multiple_testing_correction": "Benjamini-Hochberg within each clustering method x k solution",
        "alpha": float(args.alpha),
        "absolute_correlation_threshold": float(args.abs_corr_threshold),
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
            "Exploratory feature-feature association screening among features retained in Step 4. "
            "BH FDR correction addresses multiplicity within each solution but not selection induced by prior clustering/statistical screening."
        ),
        "runs": summaries,
    }
    summary_json = output_root / "01_fvf_correlation_summary.json"
    write_json(summary_json, stage_manifest)

    print("\nStep 5 feature-feature correlation completed.")
    print(f"Summary CSV:  {summary_csv}")
    print(f"Summary JSON: {summary_json}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
