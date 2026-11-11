#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
import pandas as pd

from metguild_common import (
    atomic_write_bundle,
    fit_kmeans,
    load_analysis,
    sha256_file,
    software_versions,
)


LOG = logging.getLogger(__name__)


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Produce KMeans assignments and MDS plots from the validated Step 1 "
            "manifest and its saved embedding. Supports cleaned CSV/TSV matrices "
            "with binary or quantitative pathway/completeness values."
        ),
        epilog=(
            "Run the updated Step 1 first. Use fresh output paths: existing files "
            "are never overwritten. Relative paths are resolved from the current directory."
        ),
    )
    parser.add_argument("--run-dir", help="Dataset-specific Step 1 run directory. Derives the manifest and KMeans output directory automatically.")
    parser.add_argument("--manifest", help="Step 1 JSON manifest. Default: manifests/optimal_k_values.json when --run-dir is not used.")
    parser.add_argument("--output-dir", help="KMeans output directory. Default: 03_04_clustering/01_kmeans when --run-dir is not used.")
    return parser


def cluster_colors(k):
    name = "tab10" if k <= 10 else ("tab20" if k <= 20 else "hsv")
    cmap = matplotlib.colormaps[name]
    if k > 10:
        cmap = cmap.resampled(k)
    return {str(cluster): to_hex(cmap(index)) for index, cluster in enumerate(range(1, k + 1))}


def check_labels(labels, n, requested_k):
    labels = np.asarray(labels)
    if labels.shape != (n,):
        raise ValueError(f"Assignment shape {labels.shape} does not match {n} genomes.")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("KMeans assignments must be integer cluster labels.")
    expected = np.arange(1, requested_k + 1)
    if not np.array_equal(np.unique(labels), expected):
        raise ValueError(
            f"Requested k={requested_k}, but obtained labels {np.unique(labels).tolist()}; "
            "no successful clustering output will be written."
        )
    return labels


def verify_selection_fit(manifest, k, fit_info):
    """The same saved coordinates and settings should reproduce the trial fit."""
    matches = [record for record in manifest.get("score_records", [])
               if record.get("requested_k") == k and record.get("labels_sha256")]
    for record in matches:
        if record["labels_sha256"] != fit_info.get("labels_sha256"):
            raise ValueError(
                f"KMeans labels for k={k} differ from the Step 1 trial fit. "
                "Use the same software environment or rerun Step 1 with the current environment."
            )


def save_2d(path, coordinates, labels, colors, k):
    fitted_dimensions = coordinates.shape[1]
    display = np.zeros((len(coordinates), 2), dtype=float)
    display[:, :min(2, fitted_dimensions)] = coordinates[:, :2]
    fig, ax = plt.subplots(figsize=(11, 8))
    try:
        for cluster in range(1, k + 1):
            selected = labels == cluster
            ax.scatter(display[selected, 0], display[selected, 1],
                       s=50, color=colors[str(cluster)], label=f"Cluster {cluster}",
                       edgecolors="black", linewidths=0.25, alpha=0.9)
        ax.set_title(
            f"KMeans: {k} clusters\n"
            f"First two axes of the {fitted_dimensions}-dimensional classical MDS embedding"
        )
        ax.set_xlabel("MDS 1")
        ax.set_ylabel("MDS 2" if fitted_dimensions >= 2 else "MDS 2 (zero-padded for display)")
        ax.legend(title="Cluster membership", loc="center left", bbox_to_anchor=(1.02, 0.5),
                  fontsize="small", ncol=1 if k <= 15 else 2)
        ax.grid(alpha=0.15)
        fig.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches="tight", format="png")
    finally:
        plt.close(fig)


def save_3d(path, coordinates, genomes, labels, colors, k):
    import plotly.express as px

    fitted_dimensions = coordinates.shape[1]
    display = np.zeros((len(coordinates), 3), dtype=float)
    display[:, :min(3, fitted_dimensions)] = coordinates[:, :3]
    frame = pd.DataFrame(display, columns=["MDS1", "MDS2", "MDS3"])
    frame["Genome_ID"] = list(genomes)
    frame["Cluster"] = labels.astype(str)
    projection_note = f"first 3 axes of {fitted_dimensions} fitted dimensions"
    if fitted_dimensions < 3:
        projection_note = f"{fitted_dimensions} fitted dimensions; remaining display axes are zero"
    fig = px.scatter_3d(
        frame, x="MDS1", y="MDS2", z="MDS3", color="Cluster",
        color_discrete_map=colors,
        category_orders={"Cluster": [str(value) for value in range(1, k + 1)]},
        hover_name="Genome_ID",
        title=f"KMeans: {k} clusters â€” classical MDS ({projection_note})",
        labels={"Cluster": "Cluster membership", "MDS1": "MDS 1", "MDS2": "MDS 2", "MDS3": "MDS 3"},
    )
    fig.update_traces(marker={"size": 5, "opacity": 0.9})
    fig.update_layout(legend_title_text="Cluster membership")
    # A self-contained file works on HPC output downloads without internet access.
    fig.write_html(str(path), include_plotlyjs=True, full_html=True, auto_open=False)


def main(argv=None):
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        if args.run_dir:
            if args.manifest or args.output_dir:
                raise ValueError("Use --run-dir by itself, or use --manifest/--output-dir without --run-dir.")
            run_dir = Path(args.run_dir).expanduser().resolve()
            args.manifest = str(run_dir / "manifests" / "optimal_k_values.json")
            args.output_dir = str(run_dir / "03_04_clustering" / "01_kmeans")
        else:
            args.manifest = args.manifest or "manifests/optimal_k_values.json"
            args.output_dir = args.output_dir or "03_04_clustering/01_kmeans"
        # Fail before computing or staging outputs if a required plotting dependency is absent.
        import plotly.express  # noqa: F401

        analysis = load_analysis(args.manifest)
        manifest = analysis["manifest"]
        manifest_path = Path(analysis["manifest_path"])
        genomes = list(analysis["genomes"])
        features = analysis["features"]
        coordinates = np.asarray(analysis["embedding"], dtype=float)
        k_values = list(analysis["k_values"])
        if features.index.tolist() != genomes or coordinates.shape[0] != len(genomes):
            raise ValueError("Feature, identifier, and saved embedding row order do not agree.")
        if coordinates.ndim != 2 or coordinates.shape[1] < 1 or not np.isfinite(coordinates).all():
            raise ValueError("The saved MDS embedding is empty or nonfinite.")
        if "Cluster" in features.columns:
            raise ValueError("The feature name 'Cluster' is reserved for assignment output.")

        output_dir = Path(args.output_dir).expanduser().resolve()
        base_name = Path(manifest["cleaned_output_file"]).stem
        stage_manifest_path = output_dir / "kmeans_manifest.json"
        embedding_path = output_dir / "embedding_coordinates.tsv"
        settings = manifest["kmeans"]
        run_time = datetime.now(timezone.utc).isoformat(timespec="seconds")
        provenance = {
            "step1_manifest": str(manifest_path),
            "step1_manifest_sha256": sha256_file(manifest_path),
            "step1_cleaned_output_file": manifest["cleaned_output_file"],
            "checked_input_artifacts": {
                name: manifest["artifacts"][name] for name in analysis["artifact_paths"]
            },
            "candidate_k_by_algorithm": manifest.get("candidate_k_by_algorithm"),
            "preprocessing": manifest.get("preprocessing"),
            "distance": manifest.get("distance"),
            "mds": manifest.get("mds"),
            "kmeans": settings,
            "software_versions": software_versions(),
            "script_sha256": sha256_file(Path(__file__)),
            "shared_module_sha256": sha256_file(Path(__file__).with_name("metguild_common.py")),
            "method_reference": manifest.get("method_reference"),
        }
        writers = {}
        artifacts = {}

        def register(path, writer):
            def tracked(temp_path):
                writer(temp_path)
                artifacts[str(path)] = {
                    "path": str(path), "sha256": sha256_file(temp_path),
                    "size_bytes": Path(temp_path).stat().st_size,
                }
            writers[path] = tracked

        embedding_frame = pd.DataFrame(
            coordinates, columns=[f"MDS{i + 1}" for i in range(coordinates.shape[1])]
        )
        embedding_frame.insert(0, "Genomes", genomes)
        register(embedding_path, lambda temp: embedding_frame.to_csv(
            temp, sep="\t", index=False, encoding="utf-8", float_format="%.17g"
        ))

        runs = []
        for k in k_values:
            LOG.info("Fitting final KMeans for k=%d on the saved MDS embedding.", k)
            labels, fit_info = fit_kmeans(coordinates, k, settings)
            labels = check_labels(labels, len(genomes), k)
            verify_selection_fit(manifest, k, fit_info)
            exported = features.copy()
            exported.index.name = "Genomes"
            exported.insert(0, "Cluster", labels)
            exported = exported.reset_index()
            colors = cluster_colors(k)
            run_dir = output_dir / f"{k}_clusters"
            prefix = f"03_{base_name}_{k}_KM"
            csv_path = run_dir / f"{prefix}.csv"
            png_path = run_dir / f"{prefix}_2D.png"
            html_path = run_dir / f"{prefix}.html"
            json_path = run_dir / f"{prefix}_run.json"
            register(csv_path, lambda temp, frame=exported: frame.to_csv(
                temp, index=False, encoding="utf-8", float_format="%.17g"
            ))
            register(png_path, lambda temp, lab=labels.copy(), palette=colors, count=k:
                     save_2d(temp, coordinates, lab, palette, count))
            register(html_path, lambda temp, lab=labels.copy(), palette=colors, count=k:
                     save_3d(temp, coordinates, genomes, lab, palette, count))
            run = {
                "manifest_schema_version": 3,
                "stage": "kmeans_mds",
                "status": "success",
                "timestamp_utc": run_time,
                "requested_k": k,
                "realized_clusters": int(np.unique(labels).size),
                "n_genomes": len(genomes),
                "n_features": int(features.shape[1]),
                "cluster_sizes": {str(value): int(np.sum(labels == value)) for value in range(1, k + 1)},
                "fit": fit_info,
                "fitted_mds_dimensions": int(coordinates.shape[1]),
                "feature_values_in_export": "Original cleaned numeric values; no binarization or MDS replacement.",
                "cluster_colors": colors,
                "provenance": provenance,
                "run_manifest": str(json_path),
            }

            def write_run(temp, info=run, files=(csv_path, png_path, html_path)):
                payload = dict(info)
                payload["outputs"] = {"embedding": artifacts[str(embedding_path)],
                                      **{p.name: artifacts[str(p)] for p in files}}
                Path(temp).write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            register(json_path, write_run)
            runs.append(run)

        def write_stage_manifest(temp):
            payload = {
                "manifest_schema_version": 3,
                "stage": "kmeans_mds",
                "status": "success",
                "timestamp_utc": run_time,
                "tested_k_values": k_values,
                "provenance": provenance,
                "runs": runs,
                "outputs": list(artifacts.values()),
            }
            Path(temp).write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        writers[stage_manifest_path] = write_stage_manifest
        atomic_write_bundle(writers, commit_last=stage_manifest_path)
        LOG.info("Completed %d KMeans runs. Manifest: %s", len(runs), stage_manifest_path)
        print(f"KMEANS_MANIFEST={stage_manifest_path}")
        return 0
    except (OSError, ValueError, KeyError, TypeError, ImportError, FloatingPointError, csv.Error) as error:
        LOG.error("%s", error)
        return 2


if __name__ == "__main__":
    sys.exit(main())
