#!/usr/bin/env python3
"""Clean CSV/TSV matrices and select candidate k on the final MDS geometry.

Accepted matrices contain binary annotations, pathway completeness fractions,
percentages, or other nonnegative quantitative measurements. Binary defaults
to Jaccard distance; quantitative data use Ruzicka distance without binarizing
positive values. Classical MDS coordinates are shared with Step 2. Candidate
cluster counts are selected from four KMeans diagnostics plus an independent
average-linkage hierarchical silhouette diagnostic. Their unique union is
forwarded to both clustering branches by default so the methods can be compared
at equivalent k. These unknown-k selectors are MetGuild additions to the
referenced paper, which studied binary data with known k:
https://arxiv.org/html/2108.09243v2.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score

from metguild_common import (
    DISTANCES, MATRIX_TYPES, MISSING_TOKENS, PAPER_URL, SCHEMA_VERSION,
    atomic_write_bundle, clean_matrix, cut_hierarchy, fit_kmeans,
    output_separator, prepare_geometry, read_matrix, resolve_format,
    sha256_file, software_versions, validate_output_paths,
)

LOG = logging.getLogger(__name__)
SCORE_NAMES = ("inertia", "silhouette", "calinski_harabasz", "davies_bouldin")
ELBOW_RELATIVE_TOLERANCE = 1e-10


def build_parser():
    parser = argparse.ArgumentParser(
        description="Clean a CSV/TSV annotation or pathway-completeness matrix and evaluate candidate k on classical MDS coordinates.",
        epilog="Fresh outputs are required. No automatic binarization or feature/MDS-axis standardization. The quantitative branch is a MetGuild extension to the binary-data study.")
    parser.add_argument("--input", required=True, help="CSV or TSV with a unique nonempty Genomes identifier and numeric feature columns.")
    parser.add_argument("--input-format", choices=("auto", "csv", "tsv"), default="auto", help="auto uses .csv/.tsv suffix, otherwise detects comma/tab; explicit values override the suffix.")
    parser.add_argument("--run-dir", help="Fresh dataset-specific run directory. If omitted together with the legacy output arguments, the next <dataset>_run_XX directory is selected automatically.")
    parser.add_argument("--output-root", default=".", help="Parent directory for automatic run directories (default: current directory).")
    parser.add_argument("--cleaned-output", help="Legacy/custom cleaned output path. Supply this together with --plot-dir and --manifest, or omit all three for automatic run-directory mode.")
    parser.add_argument("--plot-dir", help="Legacy/custom directory for candidate-score PNGs. Supply together with --cleaned-output and --manifest.")
    parser.add_argument("--manifest", help="Legacy/custom output JSON manifest. Supply together with --cleaned-output and --plot-dir.")
    parser.add_argument("--scores-output", help="Optional KMeans score TSV; default <manifest stem>_scores.tsv.")
    parser.add_argument("--max-k", type=int, default=15, help="Exclusive upper bound: k=2..max_k-1, limited by available profiles (default15).")
    parser.add_argument("--missing-values", choices=("error", "zero"), default="error", help="Reject partial missing values by default; zero explicitly treats missing measurements as zero.")
    parser.add_argument("--matrix-type", choices=MATRIX_TYPES, default="auto", help="binary:0/1; proportion:0..1; percentage:0..100; quantitative:nonnegative. Auto detects exact binary vs other numeric values, not biological units.")
    parser.add_argument("--distance", choices=DISTANCES, default="auto", help="auto: Jaccard for binary, Ruzicka/generalized Jaccard for quantitative values. Ordinary Jaccard rejects nonbinary data.")
    parser.add_argument("--mds-components", type=int, default=3, help="Requested classical MDS axes (default3); use available positive axes if fewer exist.")
    parser.add_argument("--n-init", type=int, default=100, help="KMeans initializations per candidate (default100).")
    parser.add_argument("--kmeans-max-iter", type=int, default=100, help="Maximum KMeans iterations per initialization (default100).")
    parser.add_argument("--random-state", type=int, default=0, help="KMeans seed, 0..4294967295 (default0).")
    parser.add_argument("--candidate-source", choices=("union", "kmeans"), default="union", help="union (default): unique union of the four KMeans recommendations and the hierarchy's best original-distance silhouette candidate. kmeans retains the legacy four-KMeans-only forwarding policy.")
    return parser


def next_run_directory(source: Path, output_root: Path) -> Path:
    """Return the next unused <dataset>_run_XX path without creating it."""
    stem = source.stem.strip()
    if not stem:
        raise ValueError("The input filename must have a nonempty dataset prefix.")
    root = output_root.expanduser().resolve()
    for number in range(1, 10000):
        candidate = root / f"{stem}_run_{number:02d}"
        if not candidate.exists() and not candidate.is_symlink():
            return candidate
    raise ValueError(f"Could not allocate a fresh run directory beneath {root}.")


def resolve_paths(args):
    absolute = lambda value: Path(value).expanduser().resolve()
    source = absolute(args.input)
    if not source.is_file():
        raise ValueError(f"Input is not a regular file: {source}")

    legacy_values = (args.cleaned_output, args.plot_dir, args.manifest)
    supplied_legacy = [value is not None for value in legacy_values]
    if any(supplied_legacy) and not all(supplied_legacy):
        raise ValueError("Supply --cleaned-output, --plot-dir, and --manifest together, or omit all three for automatic run-directory mode.")
    if args.run_dir and any(supplied_legacy):
        raise ValueError("Use either --run-dir or the legacy/custom output arguments, not both.")

    if all(supplied_legacy):
        cleaned = absolute(args.cleaned_output)
        manifest = absolute(args.manifest)
        plot_dir = absolute(args.plot_dir)
        run_dir = None
    else:
        if args.run_dir:
            run_dir = absolute(args.run_dir)
            if run_dir.exists() or run_dir.is_symlink():
                raise ValueError(f"Run directory already exists: {run_dir}. Choose a fresh run directory.")
        else:
            run_dir = next_run_directory(source, Path(args.output_root))
        cleaned = run_dir / f"{source.stem}_cleaned.tsv"
        manifest = run_dir / "manifests" / "optimal_k_values.json"
        plot_dir = run_dir / "02_k_optimization"

    output_separator(cleaned)
    paths = {
        "input": source, "cleaned_matrix": cleaned,
        "kmeans_scores": absolute(args.scores_output) if args.scores_output else manifest.with_name(manifest.stem + "_scores.tsv"),
        "hierarchical_scores": manifest.with_name(manifest.stem + "_hierarchical_scores.tsv"),
        "mds_coordinates": manifest.with_name(manifest.stem + "_mds.tsv"),
        "analysis_cache": manifest.with_name(manifest.stem + "_analysis.npz"),
        "kmeans_plot": plot_dir / f"02_{cleaned.stem}_k_opt.png",
        "hierarchical_plot": plot_dir / f"02_{cleaned.stem}_hc_k_opt.png",
        "manifest": manifest,
    }
    if run_dir is not None:
        paths["run_dir"] = run_dir
    validate_output_paths([p for name, p in paths.items() if name not in {"input", "run_dir"}], protected=[source])
    return paths


def evaluate_clustering(embedding, distances, hierarchy, ks, settings):
    records, hc_records = [], []
    for k in ks:
        LOG.info("Evaluating candidate k=%d on the saved-fit representation.", k)
        record = {"requested_k": k, "realized_clusters": None,
                  **{name: None for name in SCORE_NAMES}, "status": "failed", "message": "", "warnings": []}
        hc = {"requested_k": k, "realized_clusters": None, "silhouette": None,
              "status": "failed", "message": "", "metric": "precomputed_original_distances"}
        try:
            labels, fit = fit_kmeans(embedding, k, settings)
            record.update(fit)
            scores = {
                "silhouette": float(silhouette_score(embedding, labels, metric="euclidean")),
                "calinski_harabasz": float(calinski_harabasz_score(embedding, labels)),
                "davies_bouldin": float(davies_bouldin_score(embedding, labels)),
            }
            if not all(np.isfinite(value) for value in scores.values()):
                raise ValueError("A KMeans score is nonfinite.")
            record.update(scores, status="ok")
            # Additional original-distance silhouette is diagnostic; the four
            # KM selectors above all evaluate the actual Euclidean fit space.
            record["silhouette_original_distance"] = float(silhouette_score(distances, labels, metric="precomputed"))
        except (ValueError, FloatingPointError, OverflowError) as exc:
            record.update(status="failed", message=str(exc))
            LOG.warning("KMeans k=%d excluded: %s", k, exc)
        try:
            labels, cut = cut_hierarchy(hierarchy, k)
            hc.update(cut)
            hc["silhouette"] = float(silhouette_score(distances, labels, metric="precomputed"))
            if not np.isfinite(hc["silhouette"]):
                raise ValueError("Hierarchical silhouette is nonfinite.")
            hc["status"] = "ok"
        except (ValueError, FloatingPointError, OverflowError) as exc:
            hc.update(status="failed", message=str(exc))
            LOG.warning("Hierarchy k=%d excluded: %s", k, exc)
        records.append(record)
        hc_records.append(hc)
    return records, hc_records


def find_elbow(k_values, inertia_values):
    """Consecutive-k positive second-difference heuristic; no forced elbow."""
    if len(k_values) != len(inertia_values):
        raise ValueError("Elbow k and inertia arrays must have equal lengths.")
    if len(k_values) < 3:
        return None, "Fewer than three valid k values; no second difference is available."
    ks, scores = np.asarray(k_values, dtype=int), np.asarray(inertia_values, dtype=float)
    if not np.isfinite(scores).all():
        return None, "The inertia curve contains nonfinite values."
    if not np.all(np.diff(ks) == 1):
        return None, "Valid k values are not consecutive."
    tolerance = ELBOW_RELATIVE_TOLERANCE * float(np.max(np.abs(scores)))
    differences = np.diff(scores)
    if np.any(differences > tolerance):
        return None, "Inertia increases beyond tolerance; inspect fitting stability."
    if np.max(np.abs(differences)) <= tolerance:
        return None, "The inertia curve is flat within numerical tolerance."
    curvature = np.diff(differences)
    best = int(np.argmax(curvature))
    if curvature[best] <= tolerance:
        return None, "No positive change in slope exceeds numerical tolerance."
    return int(ks[best + 1]), "Largest positive second difference on the valid consecutive curve."


def select_candidates(records, hc_records, source):
    hc_by_k = {r["requested_k"]: r for r in hc_records if r["status"] == "ok"}
    valid = sorted((r for r in records if r["status"] == "ok" and r["requested_k"] in hc_by_k), key=lambda r: r["requested_k"])
    if not valid:
        raise ValueError("No valid candidate can be run by both clustering procedures.")
    elbow, reason = find_elbow([r["requested_k"] for r in valid], [r["inertia"] for r in valid])
    optimal = {
        "elbow": elbow,
        "silhouette": max(valid, key=lambda r: r["silhouette"])["requested_k"],
        "calinski_harabasz": max(valid, key=lambda r: r["calinski_harabasz"])["requested_k"],
        "davies_bouldin": min(valid, key=lambda r: r["davies_bouldin"])["requested_k"],
    }
    km_candidates = sorted({k for k in optimal.values() if k is not None})
    hc_valid = [hc_by_k[r["requested_k"]] for r in valid]
    hc_candidate = max(hc_valid, key=lambda r: r["silhouette"])["requested_k"]
    forwarded = sorted(set(km_candidates + ([hc_candidate] if source == "union" else [])))
    details = {
        "candidate_source": source,
        "interpretation": "Candidates for inspection, not a proven biological number of guilds.",
        "kmeans_score_space": "Unstandardized classical MDS coordinates, the same embedding used for final KMeans.",
        "hierarchical_score_space": "Original shared Jaccard/Ruzicka/Euclidean distances via precomputed silhouette.",
        "hierarchical_selection_domain": "Tested counts with valid KMeans and hierarchical fits, so all forwarded counts are feasible for both branches.",
        "tie_policy": "Smallest k among exactly equal scores.",
        "elbow_method": "maximum_positive_second_difference",
        "elbow_status": "selected" if elbow is not None else "unavailable",
        "elbow_reason": reason, "elbow_relative_numerical_tolerance": ELBOW_RELATIVE_TOLERANCE,
        "paper_scope": "The reference used fixed known k; these candidate selectors are MetGuild additions.",
        "forwarding_policy": "Both final algorithms run unique_k_values for comparable counts and legacy downstream compatibility.",
    }
    return optimal, {"silhouette": hc_candidate}, {"kmeans_mds": km_candidates, "hierarchical": [hc_candidate]}, forwarded, details


def make_plot(path, records, optimal, hc_records, hc_optimal, forwarded, candidate_source, distance_metric):
    """Create the combined six-panel candidate-k diagnostic figure.

    Panels 1-4 show the four KMeans/MDS criteria, panel 5 shows the
    average-linkage hierarchical silhouette evaluated on the original
    dissimilarity matrix, and panel 6 summarizes all selected counts and the
    unique set forwarded downstream.
    """
    hc_ok = {r["requested_k"]: r for r in hc_records if r["status"] == "ok"}
    valid = sorted(
        (r for r in records if r["status"] == "ok" and r["requested_k"] in hc_ok),
        key=lambda r: r["requested_k"],
    )
    if not valid:
        raise ValueError("No common valid candidate k values are available for plotting.")

    ks = [r["requested_k"] for r in valid]
    km_panels = (
        ("inertia", "elbow", "KMeans â€” inertia / elbow", "Within-cluster squared distance (MDS)"),
        ("silhouette", "silhouette", "KMeans â€” silhouette", "Silhouette (MDS Euclidean)"),
        ("calinski_harabasz", "calinski_harabasz", "KMeans â€” Calinskiâ€“Harabasz", "Calinskiâ€“Harabasz"),
        ("davies_bouldin", "davies_bouldin", "KMeans â€” Daviesâ€“Bouldin", "Daviesâ€“Bouldin"),
    )

    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    try:
        for ax, (score, method, title, ylabel) in zip(axes.flat[:4], km_panels):
            values = [r[score] for r in valid]
            ax.plot(ks, values, marker="o", label="Evaluated k")
            selected = optimal[method]
            if selected is None:
                ax.text(0.03, 0.95, "No elbow reported", transform=ax.transAxes, va="top")
            else:
                value = next(r[score] for r in valid if r["requested_k"] == selected)
                ax.scatter(
                    selected, value, color="red", s=100, edgecolor="black",
                    zorder=5, label=f"Selected k={selected}",
                )
            ax.set(title=title, xlabel="Requested clusters (k)", ylabel=ylabel)
            ax.set_xticks(ks)
            ax.legend()
            ax.grid(alpha=0.2)

        # Panel 5: hierarchical silhouette on the original dissimilarity matrix.
        ax = axes.flat[4]
        hc_values = [hc_ok[k]["silhouette"] for k in ks]
        ax.plot(ks, hc_values, marker="o", label="Evaluated k")
        hc_selected = hc_optimal["silhouette"]
        hc_value = hc_ok[hc_selected]["silhouette"]
        ax.scatter(
            hc_selected, hc_value, color="red", s=100, edgecolor="black",
            zorder=5, label=f"Selected k={hc_selected}",
        )
        metric_label = str(distance_metric).replace("_", " ").title()
        ax.set(
            title="Hierarchical â€” silhouette",
            xlabel="Requested clusters (k)",
            ylabel=f"Silhouette ({metric_label} distance)",
        )
        ax.set_xticks(ks)
        ax.legend()
        ax.grid(alpha=0.2)

        # Panel 6: selections contributing to the final candidate set.
        ax = axes.flat[5]
        ax.axis("off")
        candidate_text = ", ".join(str(k) for k in forwarded) if forwarded else "None"
        elbow_text = str(optimal["elbow"]) if optimal["elbow"] is not None else "unavailable"
        summary_lines = [
            "Selected candidate k",
            "",
            "KMeans on classical MDS",
            f"  Elbow:                 {elbow_text}",
            f"  Silhouette:            {optimal['silhouette']}",
            f"  Calinskiâ€“Harabasz:     {optimal['calinski_harabasz']}",
            f"  Daviesâ€“Bouldin:        {optimal['davies_bouldin']}",
            "",
            "Average-linkage hierarchical",
            f"  Silhouette:            {hc_selected}",
            "",
            f"Forwarding policy: {candidate_source}",
            f"Final unique candidate set: {candidate_text}",
        ]
        ax.text(
            0.03, 0.97, "\n".join(summary_lines), transform=ax.transAxes,
            va="top", ha="left", fontsize=11, family="monospace",
            bbox={"boxstyle": "round,pad=0.6", "facecolor": "white", "edgecolor": "0.75"},
        )

        fig.suptitle(
            "MetGuild candidate-k selection: KMeans/MDS and hierarchical diagnostics",
            fontsize=16,
        )
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(path, dpi=200, format="png")
    finally:
        plt.close(fig)


def make_hc_plot(path, records, optimal):
    valid = [r for r in records if r["status"] == "ok"]
    fig, ax = plt.subplots(figsize=(9, 6))
    try:
        ax.plot([r["requested_k"] for r in valid], [r["silhouette"] for r in valid], marker="o")
        row = next(r for r in valid if r["requested_k"] == optimal["silhouette"])
        ax.scatter(row["requested_k"], row["silhouette"], color="red", edgecolor="black", s=100,
                   label=f"Selected k={row['requested_k']}", zorder=5)
        ax.set(title="Average-linkage cuts: original-distance silhouette",
               xlabel="Requested clusters (exact-count cuts)", ylabel="Silhouette")
        ax.set_xticks([r["requested_k"] for r in valid])
        ax.legend()
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(path, dpi=200, format="png")
    finally:
        plt.close(fig)


def run(args):
    for name in ("n_init", "kmeans_max_iter", "mds_components"):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive.")
    if args.max_k < 3 or not 0 <= args.random_state <= 2**32 - 1:
        raise ValueError("--max-k must be at least3 and --random-state must be within0..4294967295.")
    paths = resolve_paths(args)
    input_hash = sha256_file(paths["input"])
    fmt = resolve_format(paths["input"], args.input_format)
    genomes, features, lines = read_matrix(paths["input"], fmt)
    if input_hash != sha256_file(paths["input"]):
        raise ValueError("Input changed while being read; retry with a stable matrix.")
    cleaned, cleaning = clean_matrix(genomes, features, lines, args.missing_values)
    D, embedding, hierarchy, preprocessing, mds = prepare_geometry(cleaned, args.matrix_type, args.distance, args.mds_components)
    preprocessing.update({"missing_values_policy": args.missing_values,
                          "recognized_feature_missing_tokens_case_insensitive": sorted(MISSING_TOKENS),
                          "cleaned_output_values": "Original numeric measurements after only the recorded cleaning/missing-value policy.",
                          "genome_identifiers": "Exact nonempty unique strings; missing tokens apply only to feature cells."})
    maximum = min(args.max_k - 1, len(genomes) - 1, preprocessing["distinct_input_profiles"], preprocessing["distinct_profiles_for_clustering"])
    ks = list(range(2, maximum + 1))
    if not ks:
        raise ValueError("No feasible candidate k remains after the chosen representation.")
    settings = {"init": "k-means++", "n_init": args.n_init, "random_state": args.random_state,
                "max_iter": args.kmeans_max_iter, "tol": 1e-4, "algorithm": "lloyd"}
    LOG.info("Input %s; matrix type %s; distance %s; classical MDS axes %d.", fmt.upper(), preprocessing["matrix_type"], preprocessing["distance_metric"], embedding.shape[1])
    LOG.info("Features retained %d/%d; candidate counts %s.", cleaned.shape[1], features.shape[1], ks)
    if preprocessing["all_zero_row_indices"]:
        LOG.warning("%d genomes have all-zero profiles; they are retained and identified in the manifest.", len(preprocessing["all_zero_row_indices"]))
    records, hc_records = evaluate_clustering(embedding, D, hierarchy, ks, settings)
    optimal, hc_optimal, by_algorithm, forwarded, selection = select_candidates(records, hc_records, args.candidate_source)
    manifest = {
        "manifest_schema_version": SCHEMA_VERSION, "status": "success",
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input_file": str(paths["input"]), "input_format": fmt, "input_sha256": input_hash,
        "run_directory": str(paths.get("run_dir", paths["manifest"].parent.parent)),
        "cleaned_output_file": str(paths["cleaned_matrix"]),
        "cleaned_output_format": "csv" if output_separator(paths["cleaned_matrix"]) == "," else "tsv",
        "kopt_plot": str(paths["kmeans_plot"]), "kopt_scores": str(paths["kmeans_scores"]),
        "hierarchical_kopt_plot": str(paths["hierarchical_plot"]),
        "mds_coordinates_file": str(paths["mds_coordinates"]), "analysis_cache_file": str(paths["analysis_cache"]),
        "max_k": maximum + 1, "requested_max_k": args.max_k, "max_k_is_exclusive": True,
        "tested_k_values": ks, "valid_k_values": [r["requested_k"] for r in records if r["status"] == "ok"],
        "optimal_k_by_method": optimal, "hierarchical_optimal_k_by_method": hc_optimal,
        "candidate_k_by_algorithm": by_algorithm, "unique_k_values": forwarded,
        "n_genomes": len(genomes), "genome_order": genomes.tolist(),
        "cleaning": cleaning, "preprocessing": preprocessing,
        "distance": {"metric": preprocessing["distance_metric"], "formula": preprocessing["formula"],
                     "quantitative_extension": preprocessing["quantitative_extension"]},
        "mds": mds, "kmeans": settings,
        "hierarchical": {"linkage": "average", "distance_metric": preprocessing["distance_metric"],
                         "cut_method": "exact_count_cut_tree", "input": "Original distance matrix; not MDS coordinates."},
        "selection": selection, "score_records": records, "hierarchical_score_records": hc_records,
        "method_reference": {"url": PAPER_URL, "binary_default": "Jaccard -> three-dimensional classical MDS -> KMeans",
                             "implementation_note": "Python uses KMeans++/Lloyd; the paper used R kmeans. This follows the method, not bit-for-bit R output.",
                             "scope": "Quantitative distances and unknown-k selection are MetGuild extensions."},
        "software_versions": software_versions(), "step1_script_sha256": sha256_file(Path(__file__)),
        "shared_module_sha256": sha256_file(Path(__file__).with_name("metguild_common.py")),
    }
    artifacts, writers = {}, {}

    def register(name, callback):
        def tracked(temp):
            callback(temp)
            artifacts[name] = {"path": str(paths[name]), "sha256": sha256_file(temp), "size_bytes": Path(temp).stat().st_size}
        writers[paths[name]] = tracked

    exported = cleaned.copy()
    exported.insert(0, "Genomes", genomes)
    register("cleaned_matrix", lambda temp: exported.to_csv(temp, sep=output_separator(paths["cleaned_matrix"]), index=False, encoding="utf-8", float_format="%.17g"))
    for name, rows in (("kmeans_scores", records), ("hierarchical_scores", hc_records)):
        table = pd.DataFrame([{key: json.dumps(value) if isinstance(value, (list, dict)) else value for key, value in row.items()} for row in rows])
        register(name, lambda temp, frame=table: frame.to_csv(temp, sep="\t", index=False, encoding="utf-8", float_format="%.17g"))
    coordinates = pd.DataFrame(embedding, columns=[f"MDS{i+1}" for i in range(embedding.shape[1])])
    coordinates.insert(0, "Genomes", genomes)
    register("mds_coordinates", lambda temp: coordinates.to_csv(temp, sep="\t", index=False, encoding="utf-8", float_format="%.17g"))

    def write_cache(temp):
        with Path(temp).open("wb") as handle:
            np.savez_compressed(handle, distances=D, embedding=embedding, linkage=hierarchy,
                                genome_ids=np.asarray(genomes.tolist(), dtype=str),
                                feature_names=np.asarray(cleaned.columns.tolist(), dtype=str))
    register("analysis_cache", write_cache)
    register("kmeans_plot", lambda temp: make_plot(
        temp, records, optimal, hc_records, hc_optimal, forwarded,
        args.candidate_source, preprocessing["distance_metric"],
    ))
    register("hierarchical_plot", lambda temp: make_hc_plot(temp, hc_records, hc_optimal))

    def write_manifest(temp):
        manifest["artifacts"] = dict(artifacts)
        Path(temp).write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    writers[paths["manifest"]] = write_manifest
    atomic_write_bundle(writers, commit_last=paths["manifest"])
    for method, value in optimal.items():
        print(f"KMeans candidate by {method}: {value if value is not None else 'unavailable'}")
    print(f"Hierarchical candidate by original-distance silhouette: {hc_optimal['silhouette']}")
    print(f"Candidate source: {args.candidate_source}; values forwarded to both final methods: {forwarded}")
    print(f"RUN_DIR={paths.get('run_dir', paths['manifest'].parent.parent)}")
    print(f"MANIFEST={paths['manifest']}")
    return manifest


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = build_parser().parse_args(argv)
    try:
        run(args)
    except (ValueError, OSError, UnicodeError, csv.Error, FloatingPointError) as exc:
        LOG.error("%s", exc)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
