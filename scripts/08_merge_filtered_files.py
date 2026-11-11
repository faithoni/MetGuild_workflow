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
        raise ValueError(
            f"{path}: column {column!r} contains missing/non-numeric/non-finite values."
        )
    return s.astype(float)
