#!/usr/bin/env python3
"""MetGuild Step 4: automatic statistical feature screening.

The statistical method is selected from the matrix type recorded by the
successful Step 1 manifest:

- binary      -> Pearson chi-square test of Cluster x feature presence/absence
- proportion  -> one-way ANOVA of feature values across clusters

Each selected candidate k is screened separately for the final KMeans-MDS and
average-linkage hierarchical partitions. Benjamini-Hochberg FDR correction is
applied within each clustering solution by default. Effect sizes and basic
assumption diagnostics are exported alongside the p-values.

Important interpretation note: the same feature matrix was used to derive the
clusters and to screen cluster-discriminating features. These tests are
therefore exploratory/post-clustering feature-screening statistics, not
independent confirmatory validation of the clusters.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, f_oneway, levene

from metguild_common import (
    atomic_write_bundle,
    load_analysis,
    sha256_file,
    software_versions,
)

LOG = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Screen cluster-discriminating features using the statistical test "
            "appropriate for the Step 1 matrix type."
        ),
        epilog=(
            "Default method selection: binary -> Pearson chi-square; "
            "proportion -> one-way ANOVA. Existing output files are not overwritten."
        ),
    )
    parser.add_argument(
        "--run-dir",
        help=(
            "Dataset-specific MetGuild run directory. Derives the Step 1 manifest, "
            "clustering manifests, clustering tables, and Step 4 output directory."
        ),
    )
    parser.add_argument(
        "--manifest",
        help=(
            "Step 1 manifest when --run-dir is not used. "
            "Default: manifests/optimal_k_values.json"
        ),
    )
    parser.add_argument(
        "--output-dir",
        help=(
            "Step 4 output directory when --run-dir is not used. "
            "Default: 05_stat_analysis"
        ),
    )
    parser.add_argument(
        "--test",
        choices=("auto", "chi-square", "anova"),
        default="auto",
        help=(
            "Statistical test. auto selects from the Step 1 matrix type and is "
            "recommended. Explicit incompatible selections are rejected."
        ),
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance threshold after the selected multiple-testing correction (default: 0.05).",
    )
    parser.add_argument(
        "--correction",
        choices=("bh", "none"),
        default="bh",
        help="Multiple-testing correction within each clustering solution (default: bh).",
    )
    parser.add_argument(
        "--require-expected-ge5",
        action="store_true",
        help=(
            "Binary only: require every expected contingency-table count to be >=5 "
            "in addition to the significance threshold. By default sparse tables are "
            "flagged but not automatically excluded."
        ),
    )
    return parser


def resolve_cli_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    if args.run_dir:
        if args.manifest or args.output_dir:
            raise ValueError("Use --run-dir by itself, or --manifest/--output-dir without --run-dir.")
        run_dir = Path(args.run_dir).expanduser().resolve()
        manifest = run_dir / "manifests" / "optimal_k_values.json"
        output_dir = run_dir / "05_stat_analysis"
    else:
        manifest = Path(args.manifest or "manifests/optimal_k_values.json").expanduser().resolve()
        run_dir = manifest.parent.parent.resolve()
        output_dir = Path(args.output_dir or "05_stat_analysis").expanduser().resolve()
    if not run_dir.is_dir():
        raise ValueError(f"Run directory does not exist: {run_dir}")
    return run_dir, manifest, output_dir


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise ValueError(f"Required file is missing: {path}")
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def validate_clustering_stage(
    stage_manifest_path: Path,
    expected_stage: str,
    expected_ks: list[int],
    step1_manifest_sha256: str,
) -> dict:
    stage = read_json(stage_manifest_path)
    if stage.get("status") != "success" or stage.get("stage") != expected_stage:
        raise ValueError(f"Clustering manifest is not a successful {expected_stage} result: {stage_manifest_path}")
    if stage.get("tested_k_values") != expected_ks:
        raise ValueError(
            f"{expected_stage} candidate k values differ from Step 1: "
            f"{stage.get('tested_k_values')} != {expected_ks}"
        )
    provenance = stage.get("provenance", {})
    if provenance.get("step1_manifest_sha256") != step1_manifest_sha256:
        raise ValueError(
            f"{expected_stage} was not produced from the current Step 1 manifest. "
            "Rerun clustering before Step 4."
        )
    runs = stage.get("runs")
    if not isinstance(runs, list) or len(runs) != len(expected_ks):
        raise ValueError(f"{expected_stage} manifest has an incomplete run list.")
    observed = []
    for item in runs:
        if item.get("status") != "success":
            raise ValueError(f"{expected_stage} contains a failed run.")
        requested = item.get("requested_k")
        realized = item.get("realized_clusters")
        if type(requested) is not int or requested != realized:
            raise ValueError(
                f"{expected_stage} requested/realized cluster counts disagree: "
                f"requested={requested}, realized={realized}"
            )
        observed.append(requested)
    if observed != expected_ks:
        raise ValueError(f"{expected_stage} run order/counts disagree with Step 1.")
    return stage


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, preserving input order."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("BH correction requires finite p-values within [0,1].")
    m = len(p)
    if m == 0:
        return p.copy()
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted_ranked = ranked * m / np.arange(1, m + 1, dtype=float)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0.0, 1.0)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    return adjusted


def apply_multiplicity(results: pd.DataFrame, raw_column: str, alpha: float, correction: str) -> pd.DataFrame:
    out = results.copy()
    raw = out[raw_column].to_numpy(dtype=float)
    if correction == "bh":
        adjusted = benjamini_hochberg(raw)
        method = "Benjamini-Hochberg"
    elif correction == "none":
        adjusted = raw.copy()
        method = "none"
    else:  # argparse protects this path
        raise ValueError(f"Unknown correction: {correction}")
    out["Adjusted_P_Value"] = adjusted
    return out


def load_clustered_table(path: Path, expected_features: pd.DataFrame, expected_genomes: list[str], k: int) -> pd.DataFrame:
    if not path.is_file():
        raise ValueError(f"Required clustering output is missing: {path}")
    frame = pd.read_csv(
        path,
        encoding="utf-8",
        dtype={"Genomes": str},
        keep_default_na=False,
    )
    expected_columns = ["Genomes", "Cluster", *expected_features.columns.tolist()]
    if frame.columns.tolist() != expected_columns:
        raise ValueError(
            f"Unexpected columns in {path}. Clustering exports must contain exactly "
            "Genomes, Cluster, and the Step 1 cleaned features in the original order."
        )
    genomes = frame["Genomes"].astype(str).tolist()
    if genomes != expected_genomes or len(set(genomes)) != len(genomes):
        raise ValueError(f"Genome identifiers/order in {path} differ from Step 1.")
    clusters = pd.to_numeric(frame["Cluster"], errors="coerce")
    if clusters.isna().any() or not np.allclose(clusters, np.round(clusters)):
        raise ValueError(f"Cluster labels in {path} are not integers.")
    clusters = clusters.astype(int)
    if sorted(clusters.unique().tolist()) != list(range(1, k + 1)):
        raise ValueError(f"Expected cluster labels 1..{k} in {path}; observed {sorted(clusters.unique().tolist())}.")
    frame["Cluster"] = clusters
    values = frame[expected_features.columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    reference = expected_features.to_numpy(dtype=float)
    if values.shape != reference.shape or not np.isfinite(values).all():
        raise ValueError(f"Feature values in {path} are missing, nonnumeric, nonfinite, or have incompatible dimensions.")
    if not np.allclose(values, reference, rtol=1e-12, atol=1e-15):
        raise ValueError(f"Feature values in {path} differ from the Step 1 cleaned matrix.")
    frame.loc[:, expected_features.columns] = values
    return frame


def chi_square_screen(frame: pd.DataFrame, feature_names: list[str], alpha: float, correction: str, require_expected_ge5: bool) -> pd.DataFrame:
    rows = []
    n = len(frame)
    for feature in feature_names:
        values = frame[feature].to_numpy(dtype=float)
        if not np.isin(values, [0.0, 1.0]).all():
            raise ValueError(f"Binary Step 4 received non-binary values in feature {feature!r}.")
        contingency = pd.crosstab(frame["Cluster"], frame[feature], dropna=False)
        constant = contingency.shape[1] < 2
        if constant:
            chi2_value = 0.0
            p_value = 1.0
            dof = 0
            expected = contingency.to_numpy(dtype=float)
            status = "constant_feature"
        else:
            chi2_value, p_value, dof, expected = chi2_contingency(contingency, correction=False)
            chi2_value = float(chi2_value)
            p_value = float(p_value)
            dof = int(dof)
            expected = np.asarray(expected, dtype=float)
            status = "ok"
        min_expected = float(np.min(expected)) if expected.size else math.nan
        cells_lt5 = int(np.sum(expected < 5.0)) if expected.size else 0
        fraction_lt5 = float(cells_lt5 / expected.size) if expected.size else math.nan
        all_ge5 = bool(expected.size and np.all(expected >= 5.0))
        r, c = contingency.shape
        denominator = n * min(max(r - 1, 0), max(c - 1, 0))
        cramers_v = float(math.sqrt(chi2_value / denominator)) if denominator > 0 else 0.0
        if not constant and not all_ge5:
            status = "sparse_expected_counts"
        rows.append({
            "gene": feature,
            "Chi2_Value": chi2_value,
            "Degrees_of_Freedom": dof,
            "p_value": p_value,
            "Cramers_V": cramers_v,
            "Min_Expected_Count": min_expected,
            "Cells_Expected_LT5": cells_lt5,
            "Fraction_Cells_Expected_LT5": fraction_lt5,
            "Expected_Counts_All_GE5": all_ge5,
            "Status": status,
        })
    result = pd.DataFrame(rows)
    result = apply_multiplicity(result, "p_value", alpha, correction)
    significant = result["Adjusted_P_Value"] < alpha
    if require_expected_ge5:
        significant &= result["Expected_Counts_All_GE5"]
    significant &= result["Status"] != "constant_feature"
    result["Significant"] = significant.astype(bool)
    return result


def eta_squared(groups: list[np.ndarray]) -> float:
    all_values = np.concatenate(groups)
    grand_mean = float(np.mean(all_values))
    ss_total = float(np.sum((all_values - grand_mean) ** 2))
    if ss_total <= 0:
        return 0.0
    ss_between = float(sum(len(group) * (float(np.mean(group)) - grand_mean) ** 2 for group in groups))
    return max(0.0, min(1.0, ss_between / ss_total))


def safe_levene(groups: list[np.ndarray]) -> float:
    # Levene requires at least two nonempty groups and can be undefined when
    # all groups have zero within-group variance. Report NaN as JSON/CSV blank.
    try:
        with np.errstate(all="ignore"):
            statistic, p_value = levene(*groups, center="median")
        return float(p_value) if np.isfinite(p_value) else math.nan
    except ValueError:
        return math.nan


def anova_screen(frame: pd.DataFrame, feature_names: list[str], alpha: float, correction: str) -> pd.DataFrame:
    rows = []
    cluster_ids = sorted(frame["Cluster"].unique().tolist())
    for feature in feature_names:
        groups = [frame.loc[frame["Cluster"] == cluster, feature].to_numpy(dtype=float) for cluster in cluster_ids]
        sizes = [int(len(group)) for group in groups]
        means = [float(np.mean(group)) for group in groups]
        values = np.concatenate(groups)
        constant = bool(np.allclose(values, values[0], rtol=0.0, atol=0.0))
        if constant:
            f_value = 0.0
            p_value = 1.0
            effect = 0.0
            status = "constant_feature"
        else:
            with np.errstate(all="ignore"):
                f_value, p_value = f_oneway(*groups)
            f_value = float(f_value)
            p_value = float(p_value)
            if not np.isfinite(f_value) or not np.isfinite(p_value):
                # This can occur with degenerate within-cluster variance. Keep a
                # conservative p-value instead of silently dropping the feature.
                f_value = 0.0 if not np.isfinite(f_value) else f_value
                p_value = 1.0
                status = "degenerate_anova"
            else:
                status = "ok"
            effect = eta_squared(groups)
        rows.append({
            "Dependent_Variable": feature,
            "F_Value": f_value,
            "Pr_Value": p_value,
            "Eta_Squared": effect,
            "Group_Count": len(groups),
            "Min_Group_Size": min(sizes),
            "Max_Group_Size": max(sizes),
            "Group_Sizes": json.dumps({str(cluster): size for cluster, size in zip(cluster_ids, sizes)}, sort_keys=True),
            "Group_Means": json.dumps({str(cluster): mean for cluster, mean in zip(cluster_ids, means)}, sort_keys=True),
            "Levene_P_Value": safe_levene(groups),
            "Status": status,
        })
    result = pd.DataFrame(rows)
    result = apply_multiplicity(result, "Pr_Value", alpha, correction)
    result["Significant"] = ((result["Adjusted_P_Value"] < alpha) & (result["Status"] != "constant_feature")).astype(bool)
    return result


def dataframe_csv_writer(frame: pd.DataFrame):
    def writer(path: Path) -> None:
        frame.to_csv(path, index=False, encoding="utf-8", float_format="%.17g", na_rep="")
    return writer


def json_writer(payload_factory):
    def writer(path: Path) -> None:
        payload = payload_factory() if callable(payload_factory) else payload_factory
        Path(path).write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    return writer


def run(args: argparse.Namespace) -> Path:
    if not 0.0 < args.alpha < 1.0:
        raise ValueError("--alpha must be between 0 and 1.")
    run_dir, manifest_path, output_dir = resolve_cli_paths(args)
    analysis = load_analysis(manifest_path)
    manifest = analysis["manifest"]
    step1_hash = sha256_file(Path(analysis["manifest_path"]))
    matrix_type = manifest.get("preprocessing", {}).get("matrix_type")
    if matrix_type not in {"binary", "proportion"}:
        raise ValueError(
            "Step 4 currently supports the two MetGuild input types only: binary and proportion. "
            f"The Step 1 manifest records matrix_type={matrix_type!r}."
        )
    expected_test = "chi-square" if matrix_type == "binary" else "anova"
    selected_test = expected_test if args.test == "auto" else args.test
    if selected_test != expected_test:
        raise ValueError(
            f"Incompatible test selection: matrix_type={matrix_type!r} requires {expected_test!r}, "
            f"but {selected_test!r} was requested."
        )
    if args.require_expected_ge5 and selected_test != "chi-square":
        raise ValueError("--require-expected-ge5 applies only to binary chi-square screening.")

    k_values = list(analysis["k_values"])
    genomes = list(analysis["genomes"])
    features = analysis["features"].copy()
    feature_names = features.columns.tolist()
    base_name = Path(manifest["cleaned_output_file"]).stem

    km_stage_path = run_dir / "03_04_clustering" / "01_kmeans" / "kmeans_manifest.json"
    hc_stage_path = run_dir / "03_04_clustering" / "02_hierarchical" / "hierarchical_manifest.json"
    validate_clustering_stage(km_stage_path, "kmeans_mds", k_values, step1_hash)
    validate_clustering_stage(hc_stage_path, "hierarchical", k_values, step1_hash)

    stage_manifest_path = output_dir / "statistical_testing_manifest.json"
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    writers: dict[Path, object] = {}
    artifact_records: dict[str, dict] = {}
    run_records: list[dict] = []

    def register(path: Path, callback) -> None:
        def tracked(temp_path: Path) -> None:
            callback(temp_path)
            artifact_records[str(path)] = {
                "path": str(path),
                "sha256": sha256_file(temp_path),
                "size_bytes": Path(temp_path).stat().st_size,
            }
        writers[path] = tracked

    methods = [
        ("KM", "01_kmeans", "03", km_stage_path),
        ("HC", "02_hierarchical", "04", hc_stage_path),
    ]

    for k in k_values:
        for method, folder, input_prefix, stage_path in methods:
            clustered_path = (
                run_dir
                / "03_04_clustering"
                / folder
                / f"{k}_clusters"
                / f"{input_prefix}_{base_name}_{k}_{method}.csv"
            )
            frame = load_clustered_table(clustered_path, features, genomes, k)
            output_subdir = output_dir / folder / f"{k}_clusters"
            prefix = f"05_{base_name}_{k}_{method}"

            if selected_test == "chi-square":
                results = chi_square_screen(
                    frame,
                    feature_names,
                    args.alpha,
                    args.correction,
                    args.require_expected_ge5,
                )
                significant_features = results.loc[results["Significant"], "gene"].tolist()
                results_path = output_subdir / f"{prefix}_chi_squared_results.csv"
                filtered_path = output_subdir / f"{prefix}_chi_filtered_results.csv"
                diagnostics_path = output_subdir / f"{prefix}_chi_diagnostics.csv"
                sig_data_path = output_subdir / f"{prefix}_chi_stat_sig.csv"
                non_sig_path = output_subdir / f"{prefix}_chi_non_sig.csv"
                run_json_path = output_subdir / f"{prefix}_chi_run.json"

                # Keep the main result table concise. Assumption checks and other
                # technical details are written to a separate diagnostics table.
                compact_results = results[[
                    "gene", "Chi2_Value", "p_value", "Adjusted_P_Value",
                    "Cramers_V", "Significant",
                ]].copy()
                filtered = compact_results.loc[compact_results["Significant"]].copy()
                diagnostic_results = results[[
                    "gene", "Degrees_of_Freedom", "Min_Expected_Count",
                    "Cells_Expected_LT5", "Fraction_Cells_Expected_LT5",
                    "Expected_Counts_All_GE5", "Status",
                ]].copy()

                non_significant = [name for name in feature_names if name not in set(significant_features)]
                non_sig_frame = pd.DataFrame({"gene": non_significant})
                diagnostics = {
                    "features_with_sparse_expected_counts": int((results["Expected_Counts_All_GE5"] == False).sum()),
                    "require_expected_ge5_for_significance": bool(args.require_expected_ge5),
                    "diagnostics_table": str(diagnostics_path),
                }
                effect_size = "Cramer's V"
            else:
                results = anova_screen(frame, feature_names, args.alpha, args.correction)
                significant_features = results.loc[results["Significant"], "Dependent_Variable"].tolist()
                results_path = output_subdir / f"{prefix}_anova_results.csv"
                filtered_path = output_subdir / f"{prefix}_anova_filtered_results.csv"
                diagnostics_path = output_subdir / f"{prefix}_anova_diagnostics.csv"
                sig_data_path = output_subdir / f"{prefix}_anova_stat_significant.csv"
                non_sig_path = output_subdir / f"{prefix}_anova_non_sig_features.csv"
                run_json_path = output_subdir / f"{prefix}_anova_run.json"

                compact_results = results[[
                    "Dependent_Variable", "F_Value", "Pr_Value",
                    "Adjusted_P_Value", "Eta_Squared", "Significant",
                ]].copy()
                filtered = compact_results.loc[compact_results["Significant"]].copy()
                diagnostic_results = results[[
                    "Dependent_Variable", "Group_Count", "Min_Group_Size",
                    "Max_Group_Size", "Group_Sizes", "Group_Means",
                    "Levene_P_Value", "Status",
                ]].copy()

                non_significant = [name for name in feature_names if name not in set(significant_features)]
                non_sig_frame = pd.DataFrame({"feature": non_significant})
                finite_levene = pd.to_numeric(results["Levene_P_Value"], errors="coerce")
                diagnostics = {
                    "features_with_levene_p_below_0_05": int((finite_levene < 0.05).sum()),
                    "features_with_undefined_levene_p": int(finite_levene.isna().sum()),
                    "levene_note": "Median-centered Levene p-values are diagnostic only; standard one-way ANOVA is still the configured test.",
                    "diagnostics_table": str(diagnostics_path),
                }
                effect_size = "eta squared"

            sig_data = frame[["Genomes", "Cluster", *significant_features]].copy()
            register(results_path, dataframe_csv_writer(compact_results))
            register(filtered_path, dataframe_csv_writer(filtered))
            register(diagnostics_path, dataframe_csv_writer(diagnostic_results))
            register(sig_data_path, dataframe_csv_writer(sig_data))
            register(non_sig_path, dataframe_csv_writer(non_sig_frame))

            record = {
                "manifest_schema_version": 1,
                "stage": "statistical_feature_screening",
                "status": "success",
                "timestamp_utc": timestamp,
                "matrix_type": matrix_type,
                "statistical_method": "pearson_chi_square" if selected_test == "chi-square" else "one_way_anova",
                "method_code": method,
                "requested_k": k,
                "n_genomes": len(genomes),
                "n_features_tested": len(feature_names),
                "n_significant": len(significant_features),
                "n_non_significant": len(feature_names) - len(significant_features),
                "alpha": args.alpha,
                "multiple_testing_correction": "Benjamini-Hochberg" if args.correction == "bh" else "none",
                "effect_size": effect_size,
                "diagnostics": diagnostics,
                "input_clustered_table": str(clustered_path),
                "input_clustered_table_sha256": sha256_file(clustered_path),
                "clustering_stage_manifest": str(stage_path),
                "clustering_stage_manifest_sha256": sha256_file(stage_path),
                "significant_features": significant_features,
                "interpretation": (
                    "Exploratory/post-clustering feature screening. Because the same feature matrix was used to derive "
                    "the clusters, these p-values are not independent confirmatory validation of the clusters."
                ),
                "run_manifest": str(run_json_path),
            }

            def make_run_payload(info=record, files=(results_path, filtered_path, diagnostics_path, sig_data_path, non_sig_path)):
                payload = dict(info)
                payload["outputs"] = {p.name: artifact_records[str(p)] for p in files}
                return payload

            register(run_json_path, json_writer(make_run_payload))
            run_records.append(record)
            LOG.info(
                "%s k=%d: tested %d features; retained %d after %s correction.",
                method,
                k,
                len(feature_names),
                len(significant_features),
                "BH" if args.correction == "bh" else "no",
            )

    provenance = {
        "step1_manifest": str(analysis["manifest_path"]),
        "step1_manifest_sha256": step1_hash,
        "kmeans_manifest": str(km_stage_path),
        "kmeans_manifest_sha256": sha256_file(km_stage_path),
        "hierarchical_manifest": str(hc_stage_path),
        "hierarchical_manifest_sha256": sha256_file(hc_stage_path),
        "cleaned_matrix": manifest["cleaned_output_file"],
        "matrix_type": matrix_type,
        "candidate_k_values": k_values,
        "software_versions": software_versions(),
        "script_sha256": sha256_file(Path(__file__)),
    }

    def stage_payload() -> dict:
        return {
            "manifest_schema_version": 1,
            "stage": "statistical_feature_screening",
            "status": "success",
            "timestamp_utc": timestamp,
            "matrix_type": matrix_type,
            "selected_test": selected_test,
            "selection_policy": "automatic_from_step1_matrix_type" if args.test == "auto" else "explicit_compatible_selection",
            "alpha": args.alpha,
            "multiple_testing_correction": "Benjamini-Hochberg" if args.correction == "bh" else "none",
            "require_expected_ge5": bool(args.require_expected_ge5),
            "tested_k_values": k_values,
            "provenance": provenance,
            "runs": run_records,
            "outputs": list(artifact_records.values()),
            "interpretation_warning": (
                "The clusters and screened features derive from the same matrix. Treat Step 4 as cluster-discriminating "
                "feature screening/prioritization, not independent statistical validation of the cluster solution."
            ),
        }

    writers[stage_manifest_path] = json_writer(stage_payload)
    atomic_write_bundle(writers, commit_last=stage_manifest_path)
    LOG.info("Step 4 completed with %s for matrix type %s.", selected_test, matrix_type)
    print(f"MATRIX_TYPE={matrix_type}")
    print(f"STATISTICAL_TEST={selected_test}")
    print(f"STATISTICAL_MANIFEST={stage_manifest_path}")
    return stage_manifest_path


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = build_parser().parse_args(argv)
    try:
        run(args)
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        ImportError,
        FloatingPointError,
        csv.Error,
        json.JSONDecodeError,
    ) as error:
        LOG.error("%s", error)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
