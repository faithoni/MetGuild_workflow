#!/usr/bin/env python3
"""Shared IO, distances and classical-MDS clustering for MetGuild.

Binary default: Jaccard -> classical MDS -> KMeans, following the approach
of d'Angella & Hennig (2021), https://arxiv.org/html/2108.09243v2.
Quantitative default: Ruzicka distance -> classical MDS -> KMeans. This is a
MetGuild extension for nonnegative pathway completeness/proportions, not a
quantitative-data result established by that binary-data study. No implicit
binarization, row normalization or feature/coordinate standardization is used.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import logging
import os
import platform
import tempfile
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, is_valid_linkage, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans

LOG = logging.getLogger(__name__)
SCHEMA_VERSION = 3
MISSING_TOKENS = frozenset({"", "na", "nan", "n/a", "null", "none"})
MATRIX_TYPES = ("auto", "binary", "proportion", "percentage", "quantitative")
DISTANCES = ("auto", "jaccard", "ruzicka", "euclidean")
PAPER_URL = "https://arxiv.org/html/2108.09243v2"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def software_versions():
    result = {"python": platform.python_version()}
    for package in ("numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "plotly"):
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = None
    return result


def resolve_format(path, input_format="auto"):
    if input_format not in {"auto", "csv", "tsv"}:
        raise ValueError("Input format must be auto, csv, or tsv.")
    if input_format != "auto":
        return input_format
    suffix = Path(path).suffix.lower()
    if suffix in {".csv", ".tsv"}:
        return suffix[1:]
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error as exc:
        raise ValueError("Could not identify comma/tab delimiters; set --input-format csv or tsv.") from exc
    return "csv" if delimiter == "," else "tsv"


def output_separator(path):
    suffix = Path(path).suffix.lower()
    if suffix not in {".csv", ".tsv"}:
        raise ValueError(f"Matrix output requires a .csv or .tsv extension: {path}")
    return "," if suffix == ".csv" else "\t"


def read_matrix(path, input_format="auto"):
    """Parse numeric features without changing IDs or duplicate headers.

    Accepted values include binary annotations, nonnegative quantitative
    measurements, pathway completeness fractions and completeness percentages.
    Missing feature tokens are represented as NaN for the explicit cleaning
    policy. Identifier strings, including 'NA' and leading zeros, are preserved.
    """
    fmt = resolve_format(path, input_format)
    genome_ids, rows, lines, seen = [], [], [], set()
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="," if fmt == "csv" else "\t", strict=True)
        headers = next(reader, None)
        if not headers:
            raise ValueError("Input matrix is empty or has no header.")
        if any(not name.strip() or name != name.strip() for name in headers):
            raise ValueError("Column headers must be nonempty and have no surrounding whitespace.")
        duplicates = [name for name, count in Counter(headers).items() if count > 1]
        if duplicates:
            raise ValueError(f"Duplicate column headers: {duplicates}")
        if "Genomes" not in headers:
            raise ValueError(f"A column named 'Genomes' is required; reading the input as {fmt.upper()}.")
        if "Cluster" in headers:
            raise ValueError("'Cluster' is reserved for output assignments; rename that input feature.")
        genome_index = headers.index("Genomes")
        indices = [i for i in range(len(headers)) if i != genome_index]
        if not indices:
            raise ValueError("No feature columns were supplied.")
        for row in reader:
            number = reader.line_num
            if len(row) != len(headers):
                raise ValueError(f"Matrix line {number} has {len(row)} fields; expected {len(headers)}.")
            genome = row[genome_index]
            if not genome.strip():
                raise ValueError(f"Blank genome identifier on matrix line {number}.")
            if genome in seen:
                raise ValueError(f"Duplicate genome identifier {genome!r} on matrix line {number}.")
            seen.add(genome)
            values = []
            for i in indices:
                token = row[i].strip()
                if token.casefold() in MISSING_TOKENS:
                    values.append(np.nan)
                    continue
                try:
                    if "_" in token:
                        raise ValueError
                    value = float(token)
                except ValueError as exc:
                    raise ValueError(f"Invalid number {row[i]!r}: genome {genome!r}, feature {headers[i]!r}, line {number}.") from exc
                if not np.isfinite(value):
                    raise ValueError(f"Nonfinite number {row[i]!r}: genome {genome!r}, feature {headers[i]!r}, line {number}.")
                values.append(value)
            genome_ids.append(genome)
            rows.append(values)
            lines.append(number)
    if not rows:
        raise ValueError("Input matrix contains no genome rows.")
    return (pd.Series(genome_ids, name="Genomes", dtype=object),
            pd.DataFrame(rows, columns=[headers[i] for i in indices], dtype=float), lines)


def clean_matrix(genomes, features, line_numbers, missing_values="error"):
    if missing_values not in {"error", "zero"}:
        raise ValueError("Missing-value policy must be error or zero.")
    all_missing = features.columns[features.isna().all()].tolist()
    all_zero = features.columns[features.eq(0).all()].tolist()
    cleaned = features.drop(columns=all_missing + all_zero).copy()
    missing = int(cleaned.isna().sum().sum())
    removed_after_imputation = []
    if missing:
        if missing_values == "error":
            examples = [f"{genomes.iloc[int(r)]!r}/{cleaned.columns[int(c)]!r} (line {line_numbers[int(r)]})"
                        for r, c in np.argwhere(cleaned.isna().to_numpy())[:5]]
            raise ValueError(f"Found {missing} partial missing feature values: {'; '.join(examples)}. "
                             "Correct the matrix or explicitly select --missing-values zero.")
        cleaned = cleaned.fillna(0.0)
        removed_after_imputation = cleaned.columns[cleaned.eq(0).all()].tolist()
        cleaned = cleaned.drop(columns=removed_after_imputation)
        LOG.warning("Explicitly replaced %d missing measurements with zero.", missing)
    if cleaned.shape[1] == 0:
        raise ValueError("No feature columns remain after cleaning.")
    if not np.isfinite(cleaned.to_numpy()).all():
        raise ValueError("Cleaned features contain nonfinite values.")
    return cleaned, {
        "original_n_features": int(features.shape[1]),
        "cleaned_n_features": int(cleaned.shape[1]),
        "missing_cells_in_input": int(features.isna().sum().sum()),
        "imputed_cells": missing if missing_values == "zero" else 0,
        "missing_values_policy": missing_values,
        "removed_all_missing_columns": all_missing,
        "removed_all_zero_columns": all_zero,
        "removed_all_zero_after_imputation_columns": removed_after_imputation,
        "retained_features": cleaned.columns.tolist(),
        "retained_constant_features": cleaned.columns[cleaned.nunique() == 1].tolist(),
    }


def validate_values(values, matrix_type="auto"):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or min(values.shape) == 0 or not np.isfinite(values).all():
        raise ValueError("Expected a nonempty, finite numeric feature matrix.")
    if matrix_type not in MATRIX_TYPES:
        raise ValueError(f"Unknown matrix type: {matrix_type}")
    if np.any(values < 0):
        raise ValueError("MetGuild annotation/completeness values must be nonnegative; negative values are not interpreted as absence.")
    is_binary = bool(np.all((values == 0) | (values == 1)))
    if matrix_type == "binary" and not is_binary:
        raise ValueError("--matrix-type binary requires exact 0/1 values; use proportion, percentage, or quantitative for other measurements.")
    if matrix_type == "proportion" and np.any(values > 1):
        raise ValueError("Proportion/pathway-completeness fractions must be within [0, 1].")
    if matrix_type == "percentage" and np.any(values > 100):
        raise ValueError("Percentage pathway completeness must be within [0, 100].")
    effective = ("binary" if is_binary else "quantitative") if matrix_type == "auto" else matrix_type
    return {"matrix_type_requested": matrix_type, "matrix_type": effective,
            "observed_values_are_binary": is_binary,
            "value_min": float(values.min()), "value_max": float(values.max()),
            "type_detection_note": "Auto distinguishes exact 0/1 from other nonnegative values; it does not infer biological units.",
            "all_zero_row_indices": np.flatnonzero(np.all(values == 0, axis=1)).tolist()}


def distance_matrix(values, metric="auto", matrix_type="auto"):
    X = np.asarray(values, dtype=float)
    info = validate_values(X, matrix_type)
    if metric not in DISTANCES:
        raise ValueError(f"Unknown distance: {metric}")
    if metric == "auto":
        metric = "jaccard" if info["matrix_type"] == "binary" else "ruzicka"
    if metric == "jaccard":
        if not info["observed_values_are_binary"]:
            raise ValueError("Jaccard requires exact binary values here; select ruzicka to preserve quantitative magnitudes.")
        D = squareform(pdist(X.astype(bool), metric="jaccard"))
    elif metric == "ruzicka":
        # This common positive scaling is algebraically neutral and avoids
        # overflow when summing nonnegative measurements with large units.
        scale = float(X.max())
        safe = X / scale if scale else X.copy()
        if np.any((X > 0) & (safe < np.finfo(float).tiny)):
            raise ValueError("Numeric dynamic range is too wide for stable quantitative Ruzicka distances; positive measurements would underflow.")
        difference = squareform(pdist(safe, metric="cityblock"))
        sums = safe.sum(axis=1)
        denominator = sums[:, None] + sums[None, :] + difference
        D = np.divide(2 * difference, denominator, out=np.zeros_like(difference), where=denominator > 0)
        D = np.clip(D, 0.0, 1.0)
    else:
        D = squareform(pdist(X, metric="euclidean"))
    np.fill_diagonal(D, 0.0)
    if not np.isfinite(D).all() or np.any(D < 0):
        raise ValueError("Distances are nonfinite or negative; check numeric magnitudes.")
    info.update({"distance_metric": metric, "scaling": "none", "binarization": "none",
                 "quantitative_extension": metric == "ruzicka",
                 "zero_zero_distance": 0.0,
                 "formula": {"jaccard": "1 - joint_presences / union_presences",
                             "ruzicka": "1 - sum(min(x_j,y_j)) / sum(max(x_j,y_j))",
                             "euclidean": "sqrt(sum((x_j-y_j)^2))"}[metric]})
    return D, info


def distinct_embedding_profiles(embedding):
    X = np.asarray(embedding, dtype=float)
    scale = float(np.max(np.abs(X)))
    return 1 if not scale else int(np.unique(np.round(X / scale, decimals=12), axis=0).shape[0])


def classical_mds(distances, n_components=3):
    """Classical scaling / principal coordinates, without additive correction.

    Center squared distances, eigendecompose and retain positive axes. MDS
    coordinates are not standardized. Zero padding belongs only in displays.
    """
    D = np.asarray(distances, dtype=float)
    if type(n_components) is not int or n_components < 1:
        raise ValueError("MDS components must be a positive integer.")
    if D.ndim != 2 or D.shape[0] != D.shape[1] or D.shape[0] < 2:
        raise ValueError("MDS requires a square distance matrix for at least two genomes.")
    if not np.isfinite(D).all() or np.any(D < 0) or not np.allclose(D, D.T, rtol=1e-12, atol=0):
        raise ValueError("MDS distances must be finite, nonnegative and symmetric.")
    if not np.all(np.diag(D) == 0):
        raise ValueError("MDS distances must have a zero diagonal.")
    scale = float(D.max())
    if scale == 0:
        raise ValueError("All pairwise distances are zero; no distinct profiles can be clustered.")
    squared = (D / scale) ** 2
    B = -0.5 * (squared - squared.mean(axis=0)[None, :]
                - squared.mean(axis=1)[:, None] + squared.mean())
    eigenvalues, eigenvectors = np.linalg.eigh((B + B.T) / 2)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    tolerance = float(np.max(np.abs(eigenvalues)) * np.finfo(float).eps * len(D) * 100)
    positive = eigenvalues > tolerance
    count = min(n_components, len(D) - 1, int(positive.sum()))
    if count == 0:
        raise ValueError("Classical MDS found no positive axes for clustering.")
    embedding = eigenvectors[:, :count] * np.sqrt(eigenvalues[:count]) * scale
    for j in range(count):
        pivot = int(np.argmax(np.abs(embedding[:, j])))
        if embedding[pivot, j] < 0:
            embedding[:, j] *= -1
    if not np.isfinite(embedding).all():
        raise ValueError("MDS coordinates are nonfinite; check numeric magnitudes.")
    # Report fit to the original dissimilarities, with a stated denominator;
    # this diagnostic is not SMACOF's optimized stress or a goodness cutoff.
    projected = squareform(pdist(embedding / scale))
    upper = np.triu_indices(len(D), 1)
    residual = (D / scale)[upper] - projected[upper]
    denominator = float(np.sum((D / scale)[upper] ** 2))
    normalized_error = float(np.sqrt(np.sum(residual ** 2) / denominator))
    positive_sum = float(eigenvalues[positive].sum())
    abs_sum = float(np.abs(eigenvalues).sum())
    info = {"method": "classical", "requested_components": n_components,
            "effective_components": count, "additive_correction": False,
            "coordinate_standardization": "none", "distance_scale_for_eigendecomposition": scale,
            "eigenvalues_of_scaled_distances": eigenvalues.tolist(),
            "eigenvalue_relative_tolerance": tolerance,
            "retained_positive_eigenvalue_fraction": float(eigenvalues[:count].sum() / positive_sum),
            "negative_eigenvalue_absolute_fraction": float(np.abs(eigenvalues[eigenvalues < -tolerance]).sum() / abs_sum),
            "normalized_distance_rmse": normalized_error,
            "normalized_distance_rmse_definition": "sqrt(sum((original_distance-embedded_distance)^2)/sum(original_distance^2)), i<j",
            "truncated_eigenvalue_tie": bool(count < int(positive.sum()) and np.isclose(eigenvalues[count-1], eigenvalues[count], rtol=1e-10, atol=tolerance)),
            "interpretation": "Low-dimensional distance approximation; diagnostic fractions are not evidence of biological guild validity."}
    if count < n_components:
        LOG.warning("Requested %d MDS axes; %d informative positive axes are available.", n_components, count)
    if info["truncated_eigenvalue_tie"]:
        LOG.warning("MDS truncation intersects tied eigenvalues; the selected subspace may depend on numerical eigensolver choices.")
    return embedding, info


def prepare_geometry(features, matrix_type="auto", metric="auto", n_components=3):
    X = features.to_numpy(dtype=float)
    if X.shape[0] < 3:
        raise ValueError("At least three genomes are required for candidate-k evaluation.")
    D, preparation = distance_matrix(X, metric, matrix_type)
    embedding, mds_info = classical_mds(D, n_components)
    # Identical measurements must not acquire numerical differences in MDS.
    _, inverse = np.unique(X, axis=0, return_inverse=True)
    for group in np.unique(inverse):
        rows = inverse == group
        if rows.sum() > 1:
            embedding[rows] = embedding[rows].mean(axis=0)
    preparation.update({"feature_order": features.columns.tolist(),
                        "distinct_input_profiles": int(np.unique(inverse).size),
                        "distinct_profiles_for_clustering": distinct_embedding_profiles(embedding),
                        "embedding_duplicate_policy": "Identical input profiles share their mean MDS coordinates."})
    if preparation["distinct_profiles_for_clustering"] < 2:
        raise ValueError("The MDS representation has fewer than two distinct profiles; increase --mds-components or inspect the matrix.")
    Z = linkage(squareform(D, checks=True), method="average", optimal_ordering=False)
    return D, embedding, Z, preparation, mds_info


def canonical_labels(labels):
    mapping = {}
    return np.array([mapping.setdefault(int(label), len(mapping) + 1) for label in labels], dtype=np.int64)


def label_metadata(labels, k):
    ids, counts = np.unique(labels, return_counts=True)
    actual = int(len(ids))
    if actual != k:
        raise ValueError(f"Requested k={k}, but the fit realized {actual} clusters.")
    return {"requested_k": int(k), "realized_clusters": actual,
            "cluster_sizes": {str(int(i)): int(n) for i, n in zip(ids, counts)},
            "labels_sha256": hashlib.sha256(np.asarray(labels, dtype="<i8").tobytes()).hexdigest()}


def fit_kmeans(embedding, k, settings):
    X = np.asarray(embedding, dtype=float)
    if type(k) is not int or not 2 <= k < len(X):
        raise ValueError("Candidate k must be an integer between 2 and n_genomes-1.")
    if not np.isfinite(X).all() or distinct_embedding_profiles(X) < k:
        raise ValueError(f"The embedding does not contain {k} distinct finite profiles.")
    allowed = {"init", "n_init", "random_state", "max_iter", "tol", "algorithm"}
    if set(settings) != allowed:
        raise ValueError(f"KMeans settings must contain exactly {sorted(allowed)}.")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = KMeans(n_clusters=k, **settings)
        labels = canonical_labels(model.fit_predict(X))
    info = label_metadata(labels, k)
    inertia = float(model.inertia_)
    if not np.isfinite(inertia):
        raise ValueError("KMeans returned nonfinite inertia.")
    info.update({"inertia": inertia, "n_iter": int(model.n_iter_),
                 "iteration_limit_reached": bool(model.n_iter_ >= settings["max_iter"]),
                 "warnings": [str(w.message) for w in caught]})
    return labels, info


def cut_hierarchy(Z, k):
    Z = np.asarray(Z, dtype=float)
    is_valid_linkage(Z, throw=True)
    n = len(Z) + 1
    if type(k) is not int or not 2 <= k < n:
        raise ValueError("Candidate k must be an integer between 2 and n_genomes-1.")
    labels = canonical_labels(cut_tree(Z, n_clusters=[k]).ravel())
    info = label_metadata(labels, k)
    merges = n - k
    lower, upper = float(Z[merges - 1, 2]), float(Z[merges, 2])
    tied = bool(np.isclose(lower, upper, rtol=1e-12, atol=1e-15 * float(np.max(Z[:, 2]))))
    info.update({"cut_method": "exact_count_cut_tree", "tied_cut_boundary": tied,
                 "last_included_merge_distance": lower, "next_merge_distance": upper,
                 "tie_note": "An exact-count cut through tied heights uses stored merge order; it does not establish distinct natural groups." if tied else None})
    return labels, info


def validate_output_paths(paths, protected=()):
    resolved = [Path(p).expanduser().resolve() for p in paths]
    protected = [Path(p).expanduser().resolve() for p in protected]
    for i, path in enumerate(resolved):
        if path.exists() or path.is_symlink():
            raise ValueError(f"Output already exists: {path}. Choose a fresh output location.")
        if path in protected:
            raise ValueError(f"Output would replace an input: {path}")
        for other in resolved[i + 1:] + protected:
            if path == other or path in other.parents or other in path.parents:
                raise ValueError(f"Conflicting paths: {path} and {other}")
        for parent in path.parents:
            if parent.exists() and not parent.is_dir():
                raise ValueError(f"Output parent is not a directory: {parent}")
    return resolved


def atomic_write_bundle(writers, commit_last=None):
    """Run callbacks in insertion order; publish only after all staging succeeds.

    Hard-link publication fails if a target appears concurrently. On an error
    remove only files published by this call. The success manifest is last.
    """
    targets = validate_output_paths(writers)
    callbacks = dict(zip(targets, writers.values()))
    last = Path(commit_last).expanduser().resolve() if commit_last is not None else None
    if last is not None and last not in callbacks:
        raise ValueError("commit_last must be a target in the output bundle.")
    staged, published = {}, []
    try:
        for target, writer in callbacks.items():
            target.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(dir=target.parent, prefix=f".{target.stem}.", suffix=target.suffix)
            os.close(fd)
            temporary = Path(name)
            staged[target] = temporary
            writer(temporary)
        order = [p for p in targets if p != last] + ([last] if last is not None else [])
        for target in order:
            os.link(staged[target], target)
            published.append(target)
    except BaseException:
        for target in published:
            target.unlink(missing_ok=True)
        raise
    finally:
        for temporary in staged.values():
            temporary.unlink(missing_ok=True)


def load_analysis(manifest_path):
    """Load the exact checked Step 1 distances/coordinates, never refit MDS."""
    path = Path(manifest_path).expanduser().resolve()
    with path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("manifest_schema_version") != SCHEMA_VERSION or manifest.get("status") != "success":
        raise ValueError("An updated successful schema-3 Step 1 manifest is required. Rerun the updated Step 1; older standardized-feature results are incompatible.")
    artifacts = manifest.get("artifacts", {})
    required = ("cleaned_matrix", "analysis_cache", "mds_coordinates")
    paths = {}
    for name in required:
        item = artifacts.get(name)
        if not isinstance(item, dict) or not isinstance(item.get("path"), str) or not isinstance(item.get("sha256"), str):
            raise ValueError(f"Manifest is missing the checked {name} artifact.")
        artifact = Path(item["path"]).expanduser()
        if not artifact.is_absolute():
            artifact = path.parent / artifact
        artifact = artifact.resolve()
        if not artifact.is_file() or sha256_file(artifact) != item["sha256"]:
            raise ValueError(f"Missing or modified {name} artifact: {artifact}. Rerun Step 1 with stable files.")
        paths[name] = artifact
    cleaned_legacy = Path(manifest.get("cleaned_output_file", "")).expanduser()
    if not cleaned_legacy.is_absolute():
        cleaned_legacy = path.parent / cleaned_legacy
    if cleaned_legacy.resolve() != paths["cleaned_matrix"]:
        raise ValueError("Manifest cleaned_output_file disagrees with its checked matrix artifact.")
    genomes, features, _ = read_matrix(paths["cleaned_matrix"], manifest["cleaned_output_format"])
    if not np.isfinite(features.to_numpy()).all():
        raise ValueError("The cleaned artifact contains missing/nonfinite measurements.")
    validate_values(features.to_numpy(), manifest["preprocessing"]["matrix_type"])
    n = len(genomes)
    if genomes.tolist() != manifest.get("genome_order") or features.columns.tolist() != manifest["preprocessing"].get("feature_order"):
        raise ValueError("Genome/feature order differs from the Step 1 manifest.")
    with np.load(paths["analysis_cache"], allow_pickle=False) as cache:
        if not {"distances", "embedding", "linkage", "genome_ids", "feature_names"}.issubset(cache.files):
            raise ValueError("Analysis cache is incomplete.")
        D, E, Z = (np.asarray(cache[name], dtype=float).copy() for name in ("distances", "embedding", "linkage"))
        if cache["genome_ids"].tolist() != genomes.tolist() or cache["feature_names"].tolist() != features.columns.tolist():
            raise ValueError("Analysis cache identifiers do not match the cleaned matrix.")
    if D.shape != (n, n) or E.ndim != 2 or E.shape[0] != n or E.shape[1] < 1 or Z.shape != (n-1, 4):
        raise ValueError("Cached geometry has incompatible dimensions.")
    if not all(np.isfinite(a).all() for a in (D, E, Z)) or np.any(D < 0) or not np.allclose(D, D.T) or not np.all(np.diag(D) == 0):
        raise ValueError("Cached geometry contains invalid distances/coordinates.")
    is_valid_linkage(Z, throw=True)
    ks = manifest.get("unique_k_values")
    if not isinstance(ks, list) or not ks or any(type(k) is not int or not 2 <= k < n for k in ks) or ks != sorted(set(ks)):
        raise ValueError("Manifest contains invalid candidate counts.")
    valid = {r["requested_k"] for r in manifest.get("score_records", []) if r.get("status") == "ok"}
    if not set(ks).issubset(valid):
        raise ValueError("A selected candidate was not a valid Step 1 KMeans fit.")
    for name, artifact in paths.items():
        if sha256_file(artifact) != artifacts[name]["sha256"]:
            raise ValueError(f"The {name} artifact changed while it was being loaded.")
    features.index = pd.Index(genomes.tolist(), name="Genomes")
    return {"manifest": manifest, "manifest_path": path, "genomes": genomes.tolist(),
            "features": features, "distances": D, "embedding": E, "linkage": Z,
            "k_values": ks, "artifact_paths": paths}
