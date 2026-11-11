#!/usr/bin/env python3
"""MetGuild Step 3b: annotate each clustering solution with MAG-phylum membership.

For every candidate k and for both KMeans and hierarchical clustering, this script
joins clustered genome IDs to a taxonomy table containing at least:
Genomes, Phylum
It writes one compact two-column table per clustering solution:
    Cluster, MAG_Phyla
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


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


def detect_delimiter(path: Path, requested: str) -> tuple[str, str]:
    """Return (pandas separator, resolved format).

    Auto-detection deliberately inspects file content before considering the
    filename extension. This allows a file named ``taxonomy.csv`` to contain
    tab- or whitespace-separated columns, as in many taxonomy exports.
    """
    if requested == "csv":
        return ",", "csv"
    if requested == "tsv":
        return "\t", "tsv"
    if requested == "whitespace":
        return r"\s+", "whitespace"

    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        sample = fh.read(65536)

    first_nonempty = next((line for line in sample.splitlines() if line.strip()), "")
    if not first_nonempty:
        raise ValueError(f"Taxonomy file is empty: {path}")

    # Content takes precedence over the extension.
    if "," in first_nonempty:
        return ",", "csv"
    if "\t" in first_nonempty:
        return "\t", "tsv"

    # No comma or tab in the header: treat runs of whitespace as separators.
    if len(first_nonempty.split()) >= 2:
        return r"\s+", "whitespace"

    raise ValueError(
        f"Could not determine the taxonomy delimiter from: {path}. "
        "Use --taxonomy-format csv, tsv, or whitespace."
    )


def load_taxonomy(path: Path, requested_format: str) -> tuple[pd.DataFrame, str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Taxonomy file not found: {path}")

    separator, resolved_format = detect_delimiter(path, requested_format)
    read_kwargs = {
        "sep": separator,
        "encoding": "utf-8-sig",
        "dtype": str,
        "keep_default_na": False,
    }
    if resolved_format == "whitespace":
        read_kwargs["engine"] = "python"

    frame = pd.read_csv(path, **read_kwargs)
    frame.columns = [str(c).strip() for c in frame.columns]

    if "Genomes" not in frame.columns:
        raise ValueError(
            f"Taxonomy file must contain a 'Genomes' column. Observed columns: {frame.columns.tolist()}"
        )

    phylum_candidates = [c for c in ("Phylum", "Phyla") if c in frame.columns]
    if not phylum_candidates:
        raise ValueError(
            "Taxonomy file must contain either a 'Phylum' or 'Phyla' column. "
            f"Observed columns: {frame.columns.tolist()}"
        )
    if len(phylum_candidates) > 1:
        raise ValueError(
            "Taxonomy file contains both 'Phylum' and 'Phyla'. Keep only one taxonomy column "
            "to avoid ambiguous assignments."
        )

    source_phylum_col = phylum_candidates[0]
    tax = frame[["Genomes", source_phylum_col]].copy()
    tax = tax.rename(columns={source_phylum_col: "Phylum"})
    tax["Genomes"] = tax["Genomes"].astype(str).str.strip()
    tax["Phylum"] = tax["Phylum"].astype(str).str.strip()

    if tax["Genomes"].eq("").any():
        raise ValueError("Taxonomy file contains empty Genomes values.")
    if tax["Phylum"].eq("").any():
        bad = tax.loc[tax["Phylum"].eq(""), "Genomes"].tolist()
        raise ValueError(
            "Taxonomy file contains empty phylum assignments for: " + ", ".join(bad)
        )
    if tax["Genomes"].duplicated().any():
        dup = (
            tax.loc[tax["Genomes"].duplicated(keep=False), "Genomes"]
            .drop_duplicates()
            .tolist()
        )
        raise ValueError(
            "Taxonomy file must contain exactly one row per genome. Duplicate genome IDs: "
            + ", ".join(dup)
        )

    return tax, resolved_format, source_phylum_col


def normalize_cluster(value: object, k: int, context: str) -> str:
    text = str(value).strip()
    if text.startswith("C") and text[1:].isdigit():
        number = int(text[1:])
    else:
        try:
            numeric = float(text)
        except ValueError as exc:
            raise ValueError(f"{context}: invalid cluster label {value!r}") from exc
        if not numeric.is_integer():
            raise ValueError(f"{context}: cluster label must be an integer or C<number>, got {value!r}")
        number = int(numeric)

    if number < 1 or number > k:
        raise ValueError(f"{context}: cluster label {value!r} is outside 1..{k}")
    return f"C{number}"


def load_cluster_assignments(path: Path, k: int) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Expected clustering output not found: {path}")

    frame = pd.read_csv(path, encoding="utf-8-sig", dtype={"Genomes": str})
    missing = {"Genomes", "Cluster"} - set(frame.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")

    out = frame[["Genomes", "Cluster"]].copy()
    out["Genomes"] = out["Genomes"].astype(str).str.strip()
    if out["Genomes"].eq("").any():
        raise ValueError(f"{path}: contains empty genome IDs.")
    if out["Genomes"].duplicated().any():
        dup = out.loc[out["Genomes"].duplicated(keep=False), "Genomes"].drop_duplicates().tolist()[:10]
        raise ValueError(f"{path}: duplicate genome IDs: {dup}")

    out["Cluster"] = [normalize_cluster(v, k, str(path)) for v in out["Cluster"]]
    expected = {f"C{i}" for i in range(1, k + 1)}
    observed = set(out["Cluster"])
    if observed != expected:
        raise ValueError(
            f"{path}: expected exactly clusters {sorted(expected)}, observed {sorted(observed)}"
        )
    return out


def build_cluster_taxonomy(assignments: pd.DataFrame, taxonomy: pd.DataFrame, k: int) -> pd.DataFrame:
    missing_genomes = sorted(set(assignments["Genomes"]) - set(taxonomy["Genomes"]))
    if missing_genomes:
        raise ValueError(
            f"Taxonomy is incomplete: {len(missing_genomes)} clustered genome(s) have no phylum assignment: "
            + ", ".join(missing_genomes)
        )

    joined = assignments.merge(taxonomy, on="Genomes", how="left", validate="one_to_one")
    joined = joined.sort_values(["Cluster", "Genomes"], key=lambda s: s.map(
        lambda x: int(str(x)[1:]) if str(x).startswith("C") and str(x)[1:].isdigit() else str(x)
    ) if s.name == "Cluster" else s, kind="stable")

    rows = []
    for number in range(1, k + 1):
        cluster = f"C{number}"
        part = joined.loc[joined["Cluster"].eq(cluster), ["Genomes", "Phylum"]].copy()
        if part.empty:
            raise ValueError(f"Cluster {cluster} contains no genomes.")
        pairs = [f"{genome}-{phylum}" for genome, phylum in part.itertuples(index=False, name=None)]
        rows.append({"Cluster": cluster, "MAG_Phyla": "; ".join(pairs)})

    return pd.DataFrame(rows, columns=["Cluster", "MAG_Phyla"])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="MetGuild Step 3b: annotate clustered MAGs with phylum assignments."
    )
    p.add_argument("--run-dir", required=True, help="MetGuild run directory containing completed Steps 1-3.")
    p.add_argument(
        "--taxonomy",
        required=True,
        help=(
            "Taxonomy file containing Genomes plus either Phylum or Phyla. "
            "Comma-, tab-, and whitespace-delimited files are supported."
        ),
    )
    p.add_argument(
        "--taxonomy-format",
        choices=["auto", "csv", "tsv", "whitespace"],
        default="auto",
        help="Taxonomy file format (default: auto; detects comma, tab, or whitespace).",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing taxonomy annotation directory.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    taxonomy_path = Path(args.taxonomy).expanduser().resolve()
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    step1_path = run_dir / "manifests" / "optimal_k_values.json"
    step1 = read_json(step1_path)
    if step1.get("status") != "success":
        raise ValueError("Step 1 manifest does not report status='success'.")

    k_values = [int(k) for k in step1.get("unique_k_values", [])]
    if not k_values:
        raise ValueError("Step 1 manifest contains no unique_k_values.")
    base_name = Path(step1.get("cleaned_output_file", "cleaned.tsv")).stem

    # Require successful clustering stages before annotating them.
    km_manifest_path = run_dir / "03_04_clustering" / "01_kmeans" / "kmeans_manifest.json"
    hc_manifest_path = run_dir / "03_04_clustering" / "02_hierarchical" / "hierarchical_manifest.json"
    for label, path in [("KMeans", km_manifest_path), ("hierarchical", hc_manifest_path)]:
        manifest = read_json(path)
        if manifest.get("status") != "success":
            raise ValueError(f"{label} clustering manifest does not report status='success'.")

    taxonomy, resolved_format, source_phylum_col = load_taxonomy(
        taxonomy_path, args.taxonomy_format
    )

    output_root = run_dir / "03_04_clustering" / "taxonomy"
    if output_root.exists() and any(output_root.iterdir()):
        if not args.overwrite:
            raise FileExistsError(
                f"Taxonomy annotation output already exists: {output_root}. Use --overwrite to replace it."
            )
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    methods = [
        ("KM", "01_kmeans", "03", "kmeans"),
        ("HC", "02_hierarchical", "04", "hierarchical"),
    ]
    runs = []

    print(f"Run directory: {run_dir}")
    print(f"Taxonomy file: {taxonomy_path}")
    print(f"Taxonomy format: {resolved_format}")
    print(f"Taxonomy phylum column: {source_phylum_col} -> Phylum")
    print(f"Candidate k values: {k_values}")

    for k in k_values:
        for method_code, method_dir, input_prefix, method_name in methods:
            clustered_path = (
                run_dir / "03_04_clustering" / method_dir / f"{k}_clusters"
                / f"{input_prefix}_{base_name}_{k}_{method_code}.csv"
            )
            assignments = load_cluster_assignments(clustered_path, k)
            cluster_taxonomy = build_cluster_taxonomy(assignments, taxonomy, k)

            out_dir = output_root / method_dir / f"{k}_clusters"
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / f"{base_name}_{k}_{method_code}_cluster_taxonomy.csv"
            cluster_taxonomy.to_csv(out_path, index=False)

            runs.append({
                "method_code": method_code,
                "method": method_name,
                "requested_k": int(k),
                "status": "success",
                "clustered_input": str(clustered_path),
                "clustered_input_sha256": sha256_file(clustered_path),
                "output_file": str(out_path),
                "output_sha256": sha256_file(out_path),
                "clusters": int(len(cluster_taxonomy)),
                "genomes": int(len(assignments)),
            })
            print(f"  {method_code} k={k}: {out_path}")

    manifest_path = output_root / "taxonomy_manifest.json"
    manifest = {
        "manifest_schema_version": 1,
        "stage": "taxonomy_annotation",
        "status": "success",
        "timestamp_utc": utc_now(),
        "taxonomy_source": str(taxonomy_path),
        "taxonomy_source_sha256": sha256_file(taxonomy_path),
        "taxonomy_format": resolved_format,
        "taxonomy_phylum_source_column": source_phylum_col,
        "accepted_taxonomy_columns": {
            "genome": "Genomes",
            "phylum": ["Phylum", "Phyla"],
        },
        "output_columns": ["Cluster", "MAG_Phyla"],
        "candidate_k_values": k_values,
        "runs": runs,
        "provenance": {
            "step1_manifest": str(step1_path),
            "step1_manifest_sha256": sha256_file(step1_path),
            "kmeans_manifest": str(km_manifest_path),
            "kmeans_manifest_sha256": sha256_file(km_manifest_path),
            "hierarchical_manifest": str(hc_manifest_path),
            "hierarchical_manifest_sha256": sha256_file(hc_manifest_path),
            "script_sha256": sha256_file(Path(__file__)),
            "software_versions": {
                "python": platform.python_version(),
                "pandas": pd.__version__,
            },
        },
    }
    write_json(manifest_path, manifest)

    print("\nStep 3b taxonomy annotation completed.")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
