#!/usr/bin/env python3
"""Cut the validated Step 1 hierarchy and export original feature measurements.

"""
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
from matplotlib.patches import Patch
import numpy as np
from scipy.cluster.hierarchy import dendrogram

from metguild_common import (
    atomic_write_bundle,
    cut_hierarchy,
    load_analysis,
    sha256_file,
    software_versions,
)


LOG = logging.getLogger(__name__)


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Produce hierarchical cluster assignments and dendrograms from the "
            "validated Step 1 manifest. Original binary or quantitative measurements "
            "are retained in the output CSV files."
        ),
        epilog=(
            "Run the updated Step 1 first. Exact-k cuts follow the stored merge order "
            "when heights tie; these ties are reported. Existing output files are never overwritten."
        ),
    )
    parser.add_argument("--run-dir", help="Dataset-specific Step 1 run directory. Derives the manifest and hierarchical output directory automatically.")
    parser.add_argument("--manifest", help="Step 1 JSON manifest. Default: manifests/optimal_k_values.json when --run-dir is not used.")
    parser.add_argument("--output-dir", help="Hierarchical output directory. Default: 03_04_clustering/02_hierarchical when --run-dir is not used.")
    return parser


def cluster_colors(k):
    name = "tab10" if k <= 10 else ("tab20" if k <= 20 else "hsv")
    cmap = matplotlib.colormaps[name]
    if k > 10:
        cmap = cmap.resampled(k)
    return {str(value): to_hex(cmap(index)) for index, value in enumerate(range(1, k + 1))}


def check_labels(labels, n, k):
    labels = np.asarray(labels)
    if labels.shape != (n,) or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("The hierarchy did not return one integer cluster assignment per genome.")
    if not np.array_equal(np.unique(labels), np.arange(1, k + 1)):
        raise ValueError(
            f"Requested k={k}, but obtained labels {np.unique(labels).tolist()}; "
            "no successful clustering output will be written."
        )
    return labels


def verify_selection_cut(manifest, k, cut_info):
    """Require the saved hierarchy to reproduce its Step 1 partition."""
    matches = [record for record in manifest.get("hierarchical_score_records", [])
               if record.get("requested_k") == k and record.get("labels_sha256")]
    for record in matches:
        if record["labels_sha256"] != cut_info.get("labels_sha256"):
            raise ValueError(
                f"Hierarchical labels for k={k} differ from the Step 1 trial cut. "
                "Use the same software environment or rerun Step 1 with the current environment."
            )


def save_dendrogram(path, linkage_matrix, genomes, labels, colors, k, tied):
    # All leaves remain visible; width grows with the number of genomes.
    width = min(40.0, max(11.0, len(genomes) * 0.16))
    fig, ax = plt.subplots(figsize=(width, 8))
    try:
        dendrogram(
            linkage_matrix,
            labels=genomes,
            color_threshold=0,
            above_threshold_color="black",
            leaf_rotation=90,
            leaf_font_size=7 if len(genomes) > 80 else 9,
            ax=ax,
        )
        mapping = {genome: colors[str(int(label))] for genome, label in zip(genomes, labels)}
        for tick in ax.get_xmajorticklabels():
            tick.set_color(mapping[tick.get_text()])
        title = f"Average-linkage hierarchical clustering: {k} clusters"
        if tied:
            title += "\nExact-count cut resolves tied merge heights"
        ax.set_title(title)
        ax.set_xlabel("Genomes")
        ax.set_ylabel("Linkage distance")
        handles = [Patch(facecolor=colors[str(value)], label=f"Cluster {value}")
                   for value in range(1, k + 1)]
        ax.legend(handles=handles, title="Leaf-label membership", loc="upper left",
                  bbox_to_anchor=(1.01, 1), fontsize="small", ncol=1 if k <= 15 else 2)
        fig.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches="tight", format="png")
    finally:
        plt.close(fig)


def main(argv=None):
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        if args.run_dir:
            if args.manifest or args.output_dir:
                raise ValueError("Use --run-dir by itself, or use --manifest/--output-dir without --run-dir.")
            run_dir = Path(args.run_dir).expanduser().resolve()
            args.manifest = str(run_dir / "manifests" / "optimal_k_values.json")
            args.output_dir = str(run_dir / "03_04_clustering" / "02_hierarchical")
        else:
            args.manifest = args.manifest or "manifests/optimal_k_values.json"
            args.output_dir = args.output_dir or "03_04_clustering/02_hierarchical"
        analysis = load_analysis(args.manifest)
        manifest = analysis["manifest"]
        manifest_path = Path(analysis["manifest_path"])
        genomes = list(analysis["genomes"])
        features = analysis["features"]
        linkage_matrix = np.asarray(analysis["linkage"], dtype=float)
        k_values = list(analysis["k_values"])
        if features.index.tolist() != genomes or linkage_matrix.shape != (len(genomes) - 1, 4):
            raise ValueError("Feature, identifier, and saved linkage dimensions do not agree.")
        if "Cluster" in features.columns:
            raise ValueError("The feature name 'Cluster' is reserved for assignment output.")

        output_dir = Path(args.output_dir).expanduser().resolve()
        base_name = Path(manifest["cleaned_output_file"]).stem
        stage_manifest_path = output_dir / "hierarchical_manifest.json"
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
            "hierarchical": manifest.get("hierarchical"),
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

        runs = []
        for k in k_values:
            LOG.info("Cutting the saved average-linkage hierarchy at k=%d.", k)
            labels, cut_info = cut_hierarchy(linkage_matrix, k)
            labels = check_labels(labels, len(genomes), k)
            verify_selection_cut(manifest, k, cut_info)
            tied = bool(cut_info.get("tied_cut_boundary", False))
            if tied:
                LOG.warning("k=%d resolves tied merge heights by stored merge order.", k)
            exported = features.copy()
            exported.index.name = "Genomes"
            exported.insert(0, "Cluster", labels)
            exported = exported.reset_index()
            colors = cluster_colors(k)
            run_dir = output_dir / f"{k}_clusters"
            prefix = f"04_{base_name}_{k}_HC"
            csv_path = run_dir / f"{prefix}.csv"
            png_path = run_dir / f"{prefix}.png"
            json_path = run_dir / f"{prefix}_run.json"
            register(csv_path, lambda temp, frame=exported: frame.to_csv(
                temp, index=False, encoding="utf-8", float_format="%.17g"
            ))
            register(png_path, lambda temp, lab=labels.copy(), palette=colors, count=k, ties=tied:
                     save_dendrogram(temp, linkage_matrix, genomes, lab, palette, count, ties))
            run = {
                "manifest_schema_version": 3,
                "stage": "hierarchical",
                "status": "success",
                "timestamp_utc": run_time,
                "requested_k": k,
                "realized_clusters": int(np.unique(labels).size),
                "n_genomes": len(genomes),
                "n_features": int(features.shape[1]),
                "cluster_sizes": {str(value): int(np.sum(labels == value)) for value in range(1, k + 1)},
                "cut": cut_info,
                "feature_values_in_export": "Original cleaned numeric values; no binarization.",
                "cluster_colors": colors,
                "provenance": provenance,
                "run_manifest": str(json_path),
            }

            def write_run(temp, info=run, files=(csv_path, png_path)):
                payload = dict(info)
                payload["outputs"] = {p.name: artifacts[str(p)] for p in files}
                Path(temp).write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            register(json_path, write_run)
            runs.append(run)

        def write_stage_manifest(temp):
            payload = {
                "manifest_schema_version": 3,
                "stage": "hierarchical",
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
        LOG.info("Completed %d hierarchical runs. Manifest: %s", len(runs), stage_manifest_path)
        print(f"HIERARCHICAL_MANIFEST={stage_manifest_path}")
        return 0
    except (OSError, ValueError, KeyError, TypeError, ImportError, FloatingPointError, csv.Error) as error:
        LOG.error("%s", error)
        return 2


if __name__ == "__main__":
    sys.exit(main())
