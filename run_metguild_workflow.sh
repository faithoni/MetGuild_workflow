#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------
# run_niches_workflow.sh
#
# Runs the full NICHEs workflow, saves all outputs inside
# a folder named after the original input file prefix, and creates
# a populated SQLite database linking inputs, scripts, workflow steps,
# intermediate outputs, final network outputs, and the actual contents
# of every generated CSV, TSV, and JSON scientific result.
#
# Example:
#   Input file:
#     data/your_dataset.tsv
#
#   Command:
#     bash run_niches_workflow.sh data/your_dataset.tsv
#
#   Output folder:
#     your_dataset/
#
# Usage:
#   bash run_niches_workflow.sh path/to/your_input_file.tsv
#   bash run_niches_workflow.sh path/to/your_input_file.csv
#
# Optional:
#   bash run_niches_workflow.sh path/to/your_input_file.tsv 20
#
# Step 4 statistical-test selection:
#   - Inspect feature values first: only 0/1 -> recommend chi-square;
#     other numeric values -> recommend ANOVA (a data-type heuristic only).
#   - The user is prompted to choose chi-square or ANOVA.
#   - Enter 1 for chi-square or 2 for ANOVA.
#
# After clustering, supply a CSV/TSV taxonomy file with Genomes and Phyla
# columns. Relative paths are resolved from the workflow root. Suffixes are
# extracted automatically: *_2_KM.csv -> KM2; *_4_HC.csv -> HC4.
# Output: <input_prefix>/03_04_clustering/cluster_group_taxonomy.csv
#
# Provenance and scientific-results outputs:
#   <input_prefix>/08_workflow_provenance/
#       <input_prefix>_workflow_provenance.db
#       <input_prefix>_workflow_provenance.sql
#
# The .db file contains populated result tables. The .sql file is a full
# text dump that can recreate the same populated database.
#
# This file is self-contained: no separate 09_generate_workflow_sql.py
# is required.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Locate workflow root
# ------------------------------------------------------------
WORKFLOW_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------------------------------------------
# Check Python and required packages
# ------------------------------------------------------------
if ! command -v python >/dev/null 2>&1; then
    echo "Error: Python was not found."
    echo ""
    echo "Please activate the conda environment first:"
    echo "  conda activate niches_workflow"
    exit 1
fi

python - <<'PY'
missing = []

packages = {
    "pandas": "pandas",
    "numpy": "numpy",
    "matplotlib": "matplotlib",
    "sklearn": "scikit-learn",
    "scipy": "scipy",
    "seaborn": "seaborn",
    "plotly": "plotly",
}

for module_name, package_name in packages.items():
    try:
        __import__(module_name)
    except ImportError:
        missing.append(package_name)

if missing:
    print("Error: Missing required Python packages:")
    for package in missing:
        print(f"  - {package}")
    print("")
    print("Install them with:")
    print("  conda env create -f environment.yml")
    print("  conda activate niches_workflow")
    raise SystemExit(1)

print("Python package check passed.")
PY

# ------------------------------------------------------------
# Read command-line arguments
# ------------------------------------------------------------
INPUT_FILE=${1:-}
MAX_K=${2:-15}

if [ -z "$INPUT_FILE" ]; then
    echo "Error: Please provide the path to your input dataset."
    echo ""
    echo "Usage:"
    echo "  bash run_niches_workflow.sh path/to/your_input_file.tsv"
    echo "  bash run_niches_workflow.sh path/to/your_input_file.csv"
    echo ""
    echo "Examples:"
    echo "  bash run_niches_workflow.sh data/your_dataset.tsv"
    echo "  bash run_niches_workflow.sh data/your_dataset.csv"
    echo "  bash run_niches_workflow.sh data/your_dataset.tsv 20"
    echo ""
    echo "Notes:"
    echo "  - Replace 'path/to/your_input_file.tsv' with the actual path to your dataset."
    echo "  - The input file can have any name."
    echo "  - Outputs are saved in a folder named after the input file prefix."
    echo "  - Example: data/my_dataset.tsv creates output folder my_dataset/"
    echo "  - The input values are checked to recommend a statistical test before you choose."
    echo "  - Enter 1 for chi-square or 2 for ANOVA."
    exit 1
fi

# Convert input path to absolute path
if [[ "$INPUT_FILE" = /* ]]; then
    ABS_INPUT_FILE="$INPUT_FILE"
else
    ABS_INPUT_FILE="$WORKFLOW_ROOT/$INPUT_FILE"
fi

if [ ! -f "$ABS_INPUT_FILE" ]; then
    echo "Error: Input file does not exist:"
    echo "  $ABS_INPUT_FILE"
    echo ""
    echo "Please check the dataset path and run again, for example:"
    echo "  bash run_niches_workflow.sh data/your_dataset.tsv"
    exit 1
fi

# ------------------------------------------------------------
# Inspect feature values, then ask which statistical test Step 4 should use
# ------------------------------------------------------------
echo "============================================================"
echo "Select the statistical test for Step 4"
echo "============================================================"
python - "$ABS_INPUT_FILE" <<'NICHE_RECOMMEND_TEST_PY'
import csv
import sys
from decimal import Decimal, InvalidOperation


def recommend(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel_tab if str(path).lower().endswith(".tsv") else csv.excel
        reader = csv.reader(handle, dialect=dialect)
        headers = [name.strip() for name in next(reader, [])]
        if not headers:
            print("No recommendation: the input file is empty.")
            return
        id_names = {"genomes", "genome", "genome_id", "mag_id"}
        id_indices = {i for i, name in enumerate(headers) if name.lower() in id_names}
        if not id_indices:
            id_indices = {0}
            print(f"No named genome-ID column found; treating the first column ({headers[0]!r}) as genome IDs.")
        indices = [i for i in range(len(headers)) if i not in id_indices]
        print("Ignoring genome-ID column(s): " + ", ".join(headers[i] or "(unnamed)" for i in sorted(id_indices)))
        if not indices:
            print("No recommendation: no feature columns found.")
            return
        observed = set()
        nonbinary_columns = set()
        finite_count = missing_count = invalid_count = row_count = 0
        noninteger = False
        invalid_examples = []
        missing_tokens = {"", "na", "nan", "n/a", "null", "none"}
        for line_number, row in enumerate(reader, start=2):
            if not any(cell.strip() for cell in row):
                continue
            row_count += 1
            if len(row) != len(headers):
                print(f"No recommendation: row {line_number} has {len(row)} fields; expected {len(headers)}.")
                return
            for i in indices:
                raw = row[i].strip()
                if raw.lower() in missing_tokens:
                    missing_count += 1
                    continue
                try:
                    value = Decimal(raw)
                    if not value.is_finite():
                        raise InvalidOperation
                except InvalidOperation:
                    invalid_count += 1
                    if len(invalid_examples) < 3:
                        invalid_examples.append(f"row {line_number}, {headers[i]!r}: {raw!r}")
                    continue
                finite_count += 1
                observed.add(i)
                if value not in (0, 1):
                    nonbinary_columns.add(i)
                if value != value.to_integral_value():
                    noninteger = True
        print(f"Checked {row_count} rows, {len(observed)} populated numeric feature columns, {finite_count} numeric values.")
        if missing_count:
            print(f"Note: ignored {missing_count} blank/missing cells for this recommendation only; the input was not modified.")
        if invalid_count:
            print(f"No recommendation: found {invalid_count} nonnumeric or nonfinite feature values.")
            print("Examples: " + "; ".join(invalid_examples))
            print("Check for metadata columns or invalid values before continuing.")
        elif not finite_count:
            print("No recommendation: no usable numeric feature values found.")
        elif not nonbinary_columns:
            print("Detected binary data: all observed feature values are 0 or 1.")
            print("Recommended choice: 1) Chi-square.")
        else:
            if noninteger:
                print("Detected nonbinary numeric data, including fractional/continuous-valued observations.")
            else:
                print("Detected nonbinary integer values (possibly counts, not necessarily continuous data).")
            print("Recommended choice: 2) ANOVA, based on the nonbinary values.")
            if observed - nonbinary_columns:
                print("Warning: binary-only and nonbinary feature columns are mixed; one test may not suit every feature.")
        print("This is a value-type recommendation, not a check of statistical-test assumptions.")
        print("You may choose either test below; no test is selected automatically.")


if __name__ == "__main__":
    try:
        recommend(sys.argv[1])
    except (OSError, UnicodeError, csv.Error) as error:
        print(f"No recommendation: could not inspect the input file: {error}")
NICHE_RECOMMEND_TEST_PY
echo ""
echo "1) Chi-square"
echo "2) ANOVA"

while true; do
    if ! read -r -p "Enter 1 or 2: " STATISTICAL_TEST_CHOICE; then
        echo ""
        echo "Error: No statistical-test selection was received."
        exit 1
    fi

    case "$STATISTICAL_TEST_CHOICE" in
        1)
            STATISTICAL_TEST="chi-square"
            break
            ;;
        2)
            STATISTICAL_TEST="ANOVA"
            break
            ;;
        *)
            echo "Invalid selection. Please enter 1 for chi-square or 2 for ANOVA."
            ;;
    esac
done

echo "Selected statistical test: $STATISTICAL_TEST"
echo ""

INPUT_BASENAME=$(basename "$ABS_INPUT_FILE")
INPUT_PREFIX="${INPUT_BASENAME%.*}"

OUTPUT_DIR="$WORKFLOW_ROOT/$INPUT_PREFIX"
RUN_SCRIPTS_DIR="$OUTPUT_DIR/scripts"
RUN_INPUT_FILE="$OUTPUT_DIR/$INPUT_BASENAME"

# ------------------------------------------------------------
# Prepare dataset-specific output folder
# ------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"
mkdir -p "$RUN_SCRIPTS_DIR"

cp "$ABS_INPUT_FILE" "$RUN_INPUT_FILE"
cp "$WORKFLOW_ROOT"/scripts/*.py "$RUN_SCRIPTS_DIR"/
cp "${BASH_SOURCE[0]}" "$RUN_SCRIPTS_DIR/run_niches_workflow.sh"

cd "$OUTPUT_DIR"

RUN_STARTED_AT=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

echo "============================================================"
echo "Starting NICHEs workflow"
echo "============================================================"
echo "Original input file: $ABS_INPUT_FILE"
echo "Workflow output folder: $OUTPUT_DIR"
echo "Working input file: $INPUT_BASENAME"
echo "Input file prefix: $INPUT_PREFIX"
echo "max_k: $MAX_K"
echo "Statistical test: $STATISTICAL_TEST"
echo ""
echo "All workflow outputs for this run will be saved inside:"
echo "  $OUTPUT_DIR"
echo ""

# ------------------------------------------------------------
# Step 1: Remove zero/NaN columns and optimize k
# ------------------------------------------------------------
echo "============================================================"
echo "Step 1/9: Removing zero columns and running k-optimization"
echo "============================================================"

printf "%s\n%s\n" "$INPUT_BASENAME" "$MAX_K" | python scripts/01_removing_zero_columns_K_Optimization.py

echo ""
echo "Step 1 complete."
echo "Expected outputs:"
echo "  - Cleaned file in: $OUTPUT_DIR"
echo "  - K-optimization plot folder: 02_k_optimization/"
echo "  - Manifest: manifests/optimal_k_values.json"
echo ""

if [ ! -f "manifests/optimal_k_values.json" ]; then
    echo "Error: Script 1 did not create manifests/optimal_k_values.json"
    exit 1
fi

# ------------------------------------------------------------
# Step 2: KMeans clustering and MDS visualization
# ------------------------------------------------------------
echo "============================================================"
echo "Step 2/9: Running KMeans/MDS clustering"
echo "============================================================"

python scripts/02_kmeans_mds_clustering.py

echo ""
echo "Step 2 complete."
echo "Expected outputs:"
echo "  - Clustered files in: 03_04_clustering/01_kmeans/"
echo "  - 2D cluster plots in: 03_04_clustering/01_kmeans/"
echo "  - Interactive 3D HTML plots in: 03_04_clustering/01_kmeans/"
echo ""

# ------------------------------------------------------------
# Step 3: Hierarchical clustering
# ------------------------------------------------------------
echo "============================================================"
echo "Step 3/9: Running hierarchical clustering"
echo "============================================================"

python scripts/03_hierarchical_clustering.py

echo ""
echo "Step 3 complete."
echo "Expected outputs:"
echo "  - Hierarchical clustering files in: 03_04_clustering/02_hierarchical/"
echo "  - Dendrogram plots in: 03_04_clustering/02_hierarchical/"
echo ""

# ------------------------------------------------------------
# Step 3b: Join all clustering assignments to genome taxonomy
# ------------------------------------------------------------
echo "Creating cluster-group taxonomy summary after clustering."
# Preserve interactive input on fd 3: stdin supplies the embedded Python code.
python - "$WORKFLOW_ROOT" "$OUTPUT_DIR" 3<&0 <<'NICHE_TAXONOMY_PY'
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def ask(message):
    print(message, end="", flush=True)
    line = user_input.readline()
    if not line:
        raise ValueError("No input received. Run interactively and provide the requested selection.")
    return line.strip()


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(65536)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel_tab if path.suffix.lower() == ".tsv" else csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        headers = [h.strip() for h in (reader.fieldnames or [])]
        if len([h for h in headers if h]) != len(set(h for h in headers if h)):
            raise ValueError(f"Duplicate column names in {path}")
        rows = [{str(k).strip(): (v or "").strip() for k, v in row.items()
                 if k is not None and k.strip()} for row in reader]
        return set(headers), rows


def natural(value):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", value)]


def export(root, run_dir):
    while True:
        name = ask("Taxonomy filename/path (CSV or TSV; columns Genomes, Phyla): ").strip("\"'")
        path = Path(name).expanduser()
        if not path.is_absolute():
            path = root / path
        if name and path.is_file():
            break
        print(f"File not found: {path}. Please try again.")
    headers, rows = read_table(path)
    if not {"Genomes", "Phyla"} <= headers:
        raise ValueError("Taxonomy file must have columns named Genomes and Phyla.")
    taxonomy = {}
    for row in rows:
        genome, phylum = row["Genomes"], row["Phyla"]
        if not genome:
            continue
        phylum = phylum or "Unknown_taxonomy"
        if genome in taxonomy and taxonomy[genome] != phylum:
            raise ValueError(f"Conflicting taxonomy entries for {genome}.")
        taxonomy[genome] = phylum
    if not taxonomy:
        raise ValueError("Taxonomy file contains no genome entries.")

    base = run_dir / "03_04_clustering"
    groups, group_sources = {}, {}
    missing, processed = set(), 0
    for folder, method in [("01_kmeans", "KM"), ("02_hierarchical", "HC")]:
        candidates = sorted((p for p in (base / folder).rglob("*")
                             if p.suffix.lower() in {".csv", ".tsv"}),
                            key=lambda p: natural(str(p)))
        method_files = 0
        for file in candidates:
            headers, rows = read_table(file)
            if not {"Genomes", "Cluster"} <= headers:
                continue
            rows = [r for r in rows if r["Genomes"] or r["Cluster"]]
            if not rows:
                continue
            if any(not r["Genomes"] or not r["Cluster"] for r in rows):
                raise ValueError(f"Missing genome or cluster label in {file}.")
            labels = {r["Cluster"] for r in rows}
            full_pattern = rf"C\d+_{method}\d+"
            full_labels = all(re.fullmatch(full_pattern, label) for label in labels)
            filename_match = re.search(r"_(\d+)_(KM|HC)$", file.stem, re.IGNORECASE)
            if not filename_match:
                raise ValueError(f"Cannot extract clustering suffix from {file.name}; "
                                 "expected a filename ending in _<number>_KM.csv or _<number>_HC.csv "
                                 "(or .tsv).")
            count, filename_method = filename_match.groups()
            if filename_method.upper() != method:
                raise ValueError(f"Clustering method in {file.name} does not match its {folder} folder.")
            if int(count) < 1:
                raise ValueError(f"Invalid cluster count in {file.name}.")
            suffix = f"{method}{int(count)}"
            if full_labels and any(not label.endswith("_" + suffix) for label in labels):
                raise ValueError(f"Cluster labels in {file.name} disagree with filename suffix {suffix}.")
            if not full_labels:
                if not all(re.fullmatch(r"(?:C)?\d+(?:\.0+)?", label) for label in labels):
                    raise ValueError(f"Unsupported or mixed cluster labels in {file}: {sorted(labels)}")
            print(f"Clustering file: {file.relative_to(run_dir)} -> {suffix}")
            members = defaultdict(set)
            assignments = {}
            for row in rows:
                genome, label = row["Genomes"], row["Cluster"]
                group = label if full_labels else f"C{int(float(label.removeprefix('C')))}_{suffix}"
                if genome in assignments and assignments[genome] != group:
                    raise ValueError(f"Conflicting cluster assignments for {genome} in {file}.")
                assignments[genome] = group
                if genome not in taxonomy or taxonomy[genome] == "Unknown_taxonomy":
                    missing.add(genome)
                members[group].add(f"{genome}-{taxonomy.get(genome, 'Unknown_taxonomy')}")
            for group, values in members.items():
                if group in groups:
                    raise ValueError(f"Duplicate group {group} in {file} and {group_sources[group]}; "
                                     "use distinct clustering suffixes or remove duplicate output files.")
                groups[group] = values
                group_sources[group] = file
            processed += 1
            method_files += 1
        if not method_files:
            raise ValueError(f"No nonempty Genomes/Cluster CSV or TSV files found in {base / folder}.")

    output = base / "cluster_group_taxonomy.csv"
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Cluster_group", "Genomes_taxonomy"])
        for group in sorted(groups, key=natural):
            writer.writerow([group, "; ".join(sorted(groups[group], key=natural))])
    temporary.replace(output)
    print(f"Saved {len(groups)} groups from {processed} clustering files to {output}")
    if missing:
        print("WARNING: Unknown_taxonomy used for: " + "; ".join(sorted(missing)))


if __name__ == "__main__":
    try:
        with os.fdopen(3, "r") as user_input:
            export(Path(sys.argv[1]), Path(sys.argv[2]))
    except (ValueError, OSError, csv.Error) as error:
        print(f"Cluster taxonomy export failed: {error}", file=sys.stderr)
        sys.exit(1)
NICHE_TAXONOMY_PY

# ------------------------------------------------------------
# Step 4: Statistical testing
# ------------------------------------------------------------
echo "============================================================"
echo "Step 4/9: Running statistical testing"
echo "============================================================"

if [ "$STATISTICAL_TEST" = "chi-square" ]; then
    echo "Using chi-square testing..."
    python scripts/04_chi-square_testing.py
else
    echo "Using ANOVA testing..."
    python scripts/04_ANOVA_testing.py
fi

echo ""
echo "Step 4 complete."
echo "Expected outputs:"
echo "  - Statistical analysis files in: 05_stat_analysis/01_kmeans/"
echo "  - Statistical analysis files in: 05_stat_analysis/02_hierarchical/"
echo ""

# ------------------------------------------------------------
# Step 5: Feature-vs-feature Spearman correlation
# ------------------------------------------------------------
echo "============================================================"
echo "Step 5/9: Running feature-vs-feature Spearman correlation"
echo "============================================================"

python scripts/05_spearman_feature_vs_feature_correlation.py

echo ""
echo "Step 5 complete."
echo "Expected outputs:"
echo "  - Feature-vs-feature correlation files in:"
echo "    06_correlation/01_featurevsfeature_correlation/"
echo "  - Summary CSV:"
echo "    06_correlation/01_featurevsfeature_correlation/01_fvf_correlation_summary.csv"
echo "  - Summary JSON:"
echo "    06_correlation/01_featurevsfeature_correlation/01_fvf_correlation_summary.json"
echo ""

# ------------------------------------------------------------
# Step 6: Feature-vs-cluster Spearman correlation
# ------------------------------------------------------------
echo "============================================================"
echo "Step 6/9: Running feature-vs-cluster Spearman correlation"
echo "============================================================"

python scripts/06_spearman_feature_vs_cluster_correlation.py

echo ""
echo "Step 6 complete."
echo "Expected outputs:"
echo "  - Feature-vs-cluster correlation files in:"
echo "    06_correlation/02_featurevscluster_correlation/"
echo "  - Summary CSV:"
echo "    06_correlation/02_featurevscluster_correlation/01_fvc_correlation_summary.csv"
echo "  - Summary JSON:"
echo "    06_correlation/02_featurevscluster_correlation/01_fvc_correlation_summary.json"
echo ""

# ------------------------------------------------------------
# Step 7: Cluster-vs-cluster Spearman correlation
# ------------------------------------------------------------
echo "============================================================"
echo "Step 7/9: Running cluster-vs-cluster Spearman correlation"
echo "============================================================"

python scripts/07_spearman_cluster_vs_cluster_correlation.py

echo ""
echo "Step 7 complete."
echo "Expected outputs:"
echo "  - Cluster-vs-cluster correlation files in:"
echo "    06_correlation/03_clustervscluster_correlation/"
echo "  - Summary CSV:"
echo "    06_correlation/03_clustervscluster_correlation/01_cvc_correlation_summary.csv"
echo "  - Summary JSON:"
echo "    06_correlation/03_clustervscluster_correlation/01_cvc_correlation_summary.json"
echo ""

# ------------------------------------------------------------
# Step 8: Merge filtered files into network edge/node files
# ------------------------------------------------------------
echo "============================================================"
echo "Step 8/9: Generating network edge lists and node attributes"
echo "============================================================"

python scripts/08_merge_filtered_files.py

echo ""
echo "Step 8 complete."
echo "Expected outputs:"
echo "  - Edge lists in:"
echo "    07_network_construction/edge_lists/"
echo "  - Node attributes in:"
echo "    07_network_construction/nodes_attributes/"
echo "  - Merge summary:"
echo "    07_network_construction/01_edge_list_merge_summary.csv"
echo "    07_network_construction/01_edge_list_merge_summary.json"
echo "  - Node summary:"
echo "    07_network_construction/02_node_attributes_summary.csv"
echo "    07_network_construction/02_node_attributes_summary.json"
echo ""
# ------------------------------------------------------------
# Step 9: Generate SQL workflow provenance
# ------------------------------------------------------------
echo "============================================================"
echo "Step 9/9: Building populated SQLite workflow database"
echo "============================================================"

RUN_COMPLETED_AT=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
PROVENANCE_DIR="08_workflow_provenance"
PROVENANCE_SQL="$PROVENANCE_DIR/${INPUT_PREFIX}_workflow_provenance.sql"
PROVENANCE_DB="$PROVENANCE_DIR/${INPUT_PREFIX}_workflow_provenance.db"

# Remove previous exports so a rerun cannot accidentally import stale tables.
rm -f "$PROVENANCE_SQL" "$PROVENANCE_DB"

# The SQL generator is embedded below, so this runner is self-contained.
python - \
    --run-dir "$OUTPUT_DIR" \
    --original-input "$ABS_INPUT_FILE" \
    --working-input "$INPUT_BASENAME" \
    --max-k "$MAX_K" \
    --statistical-test "$STATISTICAL_TEST" \
    --started-at "$RUN_STARTED_AT" \
    --completed-at "$RUN_COMPLETED_AT" \
    --output "$PROVENANCE_SQL" <<'NICHE_SQL_PY'
"""Generate a SQLite-compatible SQL provenance file for one NICHEs workflow run."""

from __future__ import annotations

import argparse
import hashlib
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional


@dataclass(frozen=True)
class Step:
    step_id: int
    number: int
    name: str
    script_name: Optional[str]
    description: str


@dataclass
class Artifact:
    artifact_id: int
    relative_path: str
    artifact_name: str
    artifact_role: str
    artifact_type: str
    extension: Optional[str]
    size_bytes: Optional[int]
    sha256: Optional[str]
    modified_at: Optional[str]
    producing_step_id: Optional[int]
    is_external: int = 0


STEPS = [
    Step(1, 0, "Run setup and input staging", "run_niches_workflow.sh", "Validate the environment, create the run directory, copy the input dataset, and stage workflow scripts."),
    Step(2, 1, "Zero-column removal and k optimization", "01_removing_zero_columns_K_Optimization.py", "Clean the feature matrix and determine candidate or optimal values of k."),
    Step(3, 2, "KMeans and MDS clustering", "02_kmeans_mds_clustering.py", "Generate KMeans cluster assignments and two- and three-dimensional MDS visualizations."),
    Step(4, 3, "Hierarchical clustering", "03_hierarchical_clustering.py", "Generate hierarchical cluster assignments and dendrograms."),
    Step(5, 4, "Statistical testing", None, "Run the chi-square or ANOVA test selected interactively by the user."),
    Step(6, 5, "Feature-versus-feature correlation", "05_spearman_feature_vs_feature_correlation.py", "Calculate and filter Spearman correlations between features."),
    Step(7, 6, "Feature-versus-cluster correlation", "06_spearman_feature_vs_cluster_correlation.py", "Calculate and filter Spearman correlations between features and clusters."),
    Step(8, 7, "Cluster-versus-cluster correlation", "07_spearman_cluster_vs_cluster_correlation.py", "Calculate and filter Spearman correlations between clusters."),
    Step(9, 8, "Network file construction", "08_merge_filtered_files.py", "Merge filtered correlations into network edge lists and node-attribute tables."),
    Step(10, 9, "Workflow provenance export", "run_niches_workflow.sh", "Use the SQL generator embedded in the main runner to scan the completed run and export workflow structure and file lineage."),
]

# (upstream step number, downstream step number, description)
STEP_DEPENDENCIES = [
    (0, 1, "The staged input matrix is cleaned and evaluated for k."),
    (1, 2, "The cleaned matrix and k manifest drive KMeans clustering."),
    (1, 3, "The cleaned matrix and k manifest drive hierarchical clustering."),
    (1, 4, "The cleaned feature matrix is used in statistical testing."),
    (2, 4, "KMeans cluster assignments are tested for associated features."),
    (3, 4, "Hierarchical cluster assignments are tested for associated features."),
    (1, 5, "The cleaned feature matrix supplies feature values."),
    (4, 5, "Statistically retained features are used for feature-feature correlations."),
    (1, 6, "The cleaned feature matrix supplies feature values."),
    (2, 6, "KMeans assignments supply cluster variables."),
    (3, 6, "Hierarchical assignments supply cluster variables."),
    (4, 6, "Statistically retained features are used for feature-cluster correlations."),
    (2, 7, "KMeans assignments supply cluster variables."),
    (3, 7, "Hierarchical assignments supply cluster variables."),
    (5, 8, "Filtered feature-feature correlations contribute network edges."),
    (6, 8, "Filtered feature-cluster correlations contribute network edges."),
    (7, 8, "Filtered cluster-cluster correlations contribute network edges."),
    (0, 9, "Run metadata and staged files are included in the provenance export."),
    (1, 9, "Preprocessing outputs are included in the provenance export."),
    (2, 9, "KMeans outputs are included in the provenance export."),
    (3, 9, "Hierarchical outputs are included in the provenance export."),
    (4, 9, "Statistical-testing outputs are included in the provenance export."),
    (5, 9, "Feature-feature outputs are included in the provenance export."),
    (6, 9, "Feature-cluster outputs are included in the provenance export."),
    (7, 9, "Cluster-cluster outputs are included in the provenance export."),
    (8, 9, "Network outputs are included in the provenance export."),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="Dataset-specific workflow output directory.")
    parser.add_argument("--original-input", required=True, help="Absolute path of the source input dataset.")
    parser.add_argument("--working-input", required=True, help="Filename of the copied input within the run directory.")
    parser.add_argument("--max-k", required=True, type=int, help="Maximum k supplied to the workflow.")
    parser.add_argument("--statistical-test", required=True, choices=["ANOVA", "chi-square"], help="Statistical test selected by the runner.")
    parser.add_argument("--started-at", required=True, help="UTC workflow start time in ISO 8601 format.")
    parser.add_argument("--completed-at", required=True, help="UTC workflow completion time in ISO 8601 format.")
    parser.add_argument("--output", required=True, help="SQL output path, absolute or relative to --run-dir.")
    return parser.parse_args()


def iso_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sql_value(value: object) -> str:
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    return "'" + str(value).replace("'", "''") + "'"


def classify_type(path: Path) -> str:
    suffix = path.suffix.lower()
    parts = set(path.parts)
    if suffix in {".py", ".sh"}:
        return "script"
    if suffix == ".sql":
        return "database_sql"
    if suffix in {".csv", ".tsv"}:
        return "table"
    if suffix == ".json":
        return "manifest" if "manifests" in parts else "summary_json"
    if suffix in {".png", ".svg", ".pdf", ".jpg", ".jpeg"}:
        return "plot_or_figure"
    if suffix == ".html":
        return "interactive_plot"
    if suffix in {".log", ".out", ".err"}:
        return "log"
    if suffix in {".txt", ".md"}:
        return "text"
    return "other"


def producing_step_number(relative_path: str, working_input: str, output_relative: str) -> Optional[int]:
    path = Path(relative_path)
    text = path.as_posix()

    if text == output_relative:
        return 9
    if text == working_input or text.startswith("scripts/"):
        return 0
    if text.startswith("02_k_optimization/") or text.startswith("manifests/"):
        return 1
    if text.startswith("03_04_clustering/01_kmeans/"):
        return 2
    if text.startswith("03_04_clustering/02_hierarchical/"):
        return 3
    if text == "03_04_clustering/cluster_group_taxonomy.csv":
        return 3
    if text.startswith("05_stat_analysis/"):
        return 4
    if text.startswith("06_correlation/01_featurevsfeature_correlation/"):
        return 5
    if text.startswith("06_correlation/02_featurevscluster_correlation/"):
        return 6
    if text.startswith("06_correlation/03_clustervscluster_correlation/"):
        return 7
    if text.startswith("07_network_construction/"):
        return 8
    if text.startswith("08_workflow_provenance/"):
        return 9

    # Data files produced in the run root after the staged input are treated as
    # Step 1 cleaning/preprocessing outputs.
    if len(path.parts) == 1 and path.suffix.lower() in {".csv", ".tsv", ".txt", ".json"}:
        return 1
    return None


def role_for(relative_path: str, working_input: str, output_relative: str) -> str:
    if relative_path == working_input:
        return "staged_input"
    if relative_path == output_relative:
        return "provenance_output"
    if relative_path.startswith("scripts/"):
        return "script"
    return "workflow_output"


def is_tabular_or_manifest(artifact: Artifact) -> bool:
    return artifact.artifact_type in {"table", "manifest", "summary_json", "text"}


def script_for_step(step_number: int, statistical_test: str) -> Optional[str]:
    if step_number == 4:
        return "04_chi-square_testing.py" if statistical_test == "chi-square" else "04_ANOVA_testing.py"
    for step in STEPS:
        if step.number == step_number:
            return step.script_name
    return None


def scan_artifacts(run_dir: Path, original_input: Path, working_input: str, output_path: Path) -> list[Artifact]:
    output_relative = output_path.relative_to(run_dir).as_posix()
    records: list[Artifact] = []
    next_id = 1

    if original_input.is_file():
        records.append(
            Artifact(
                artifact_id=next_id,
                relative_path=str(original_input),
                artifact_name=original_input.name,
                artifact_role="external_input",
                artifact_type="table" if original_input.suffix.lower() in {".csv", ".tsv"} else classify_type(original_input),
                extension=original_input.suffix.lower() or None,
                size_bytes=original_input.stat().st_size,
                sha256=sha256sum(original_input),
                modified_at=iso_mtime(original_input),
                producing_step_id=None,
                is_external=1,
            )
        )
        next_id += 1

    for path in sorted((p for p in run_dir.rglob("*") if p.is_file()), key=lambda p: p.relative_to(run_dir).as_posix()):
        relative = path.relative_to(run_dir).as_posix()
        if relative == output_relative:
            continue
        step_number = producing_step_number(relative, working_input, output_relative)
        step_id = next((step.step_id for step in STEPS if step.number == step_number), None)
        records.append(
            Artifact(
                artifact_id=next_id,
                relative_path=relative,
                artifact_name=path.name,
                artifact_role=role_for(relative, working_input, output_relative),
                artifact_type=classify_type(path.relative_to(run_dir)),
                extension=path.suffix.lower() or None,
                size_bytes=path.stat().st_size,
                sha256=sha256sum(path),
                modified_at=iso_mtime(path),
                producing_step_id=step_id,
                is_external=0,
            )
        )
        next_id += 1

    # The SQL file describes itself, but its checksum and size cannot be stored
    # inside itself without changing the file. They are intentionally NULL.
    records.append(
        Artifact(
            artifact_id=next_id,
            relative_path=output_relative,
            artifact_name=output_path.name,
            artifact_role="provenance_output",
            artifact_type="database_sql",
            extension=".sql",
            size_bytes=None,
            sha256=None,
            modified_at=None,
            producing_step_id=next(step.step_id for step in STEPS if step.number == 9),
            is_external=0,
        )
    )
    return records


def artifacts_by_step(artifacts: Iterable[Artifact], step_number: int) -> list[Artifact]:
    step_id = next(step.step_id for step in STEPS if step.number == step_number)
    return [artifact for artifact in artifacts if artifact.producing_step_id == step_id]


def choose_step_inputs(artifacts: list[Artifact], step_number: int, working_input: str) -> list[Artifact]:
    external = [a for a in artifacts if a.artifact_role == "external_input"]
    staged = [a for a in artifacts if a.relative_path == working_input]
    root_cleaned = [
        a for a in artifacts
        if a.producing_step_id == next(step.step_id for step in STEPS if step.number == 1)
        and Path(a.relative_path).parent == Path(".")
        and is_tabular_or_manifest(a)
    ]
    step1_key = [a for a in artifacts_by_step(artifacts, 1) if is_tabular_or_manifest(a)]
    step2_tables = [a for a in artifacts_by_step(artifacts, 2) if is_tabular_or_manifest(a)]
    step3_tables = [a for a in artifacts_by_step(artifacts, 3) if is_tabular_or_manifest(a)]
    step4_tables = [a for a in artifacts_by_step(artifacts, 4) if is_tabular_or_manifest(a)]
    step5_data = [a for a in artifacts_by_step(artifacts, 5) if is_tabular_or_manifest(a)]
    step6_data = [a for a in artifacts_by_step(artifacts, 6) if is_tabular_or_manifest(a)]
    step7_data = [a for a in artifacts_by_step(artifacts, 7) if is_tabular_or_manifest(a)]

    selections = {
        0: external,
        1: staged,
        2: step1_key,
        3: step1_key,
        4: root_cleaned + step2_tables + step3_tables,
        5: root_cleaned + step4_tables,
        6: root_cleaned + step2_tables + step3_tables + step4_tables,
        7: step2_tables + step3_tables,
        8: step5_data + step6_data + step7_data,
        9: [a for a in artifacts if a.artifact_role != "provenance_output"],
    }

    # Preserve order while removing duplicates.
    seen: set[int] = set()
    unique: list[Artifact] = []
    for artifact in selections.get(step_number, []):
        if artifact.artifact_id not in seen:
            seen.add(artifact.artifact_id)
            unique.append(artifact)
    return unique


def generate_sql(args: argparse.Namespace) -> str:
    run_dir = Path(args.run_dir).resolve()
    original_input = Path(args.original_input).resolve()
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = run_dir / output_path
    output_path = output_path.resolve()

    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory does not exist: {run_dir}")
    try:
        output_relative = output_path.relative_to(run_dir).as_posix()
    except ValueError as exc:
        raise ValueError("--output must be inside --run-dir so it can be represented by a relative path.") from exc

    artifacts = scan_artifacts(run_dir, original_input, args.working_input, output_path)
    step_by_number = {step.number: step for step in STEPS}

    lines: list[str] = []
    lines.extend([
        "-- NICHEs workflow provenance export",
        f"-- Dataset: {Path(args.working_input).stem}",
        f"-- Generated: {args.completed_at}",
        "-- SQLite-compatible SQL. Import with:",
        f"--   sqlite3 {Path(args.working_input).stem}_workflow_provenance.db < {Path(output_relative).name}",
        "",
        "PRAGMA foreign_keys = ON;",
        "BEGIN TRANSACTION;",
        "",
        "DROP VIEW IF EXISTS artifact_connections;",
        "DROP VIEW IF EXISTS step_file_inventory;",
        "DROP TABLE IF EXISTS artifact_dependencies;",
        "DROP TABLE IF EXISTS step_artifacts;",
        "DROP TABLE IF EXISTS artifacts;",
        "DROP TABLE IF EXISTS step_dependencies;",
        "DROP TABLE IF EXISTS workflow_steps;",
        "DROP TABLE IF EXISTS workflow_runs;",
        "",
        "CREATE TABLE workflow_runs (",
        "    run_id INTEGER PRIMARY KEY,",
        "    run_name TEXT NOT NULL,",
        "    dataset_name TEXT NOT NULL,",
        "    original_input_path TEXT NOT NULL,",
        "    run_directory TEXT NOT NULL,",
        "    max_k INTEGER NOT NULL,",
        "    statistical_test TEXT NOT NULL,",
        "    started_at TEXT NOT NULL,",
        "    completed_at TEXT NOT NULL,",
        "    status TEXT NOT NULL CHECK (status IN ('completed', 'failed', 'running'))",
        ");",
        "",
        "CREATE TABLE workflow_steps (",
        "    step_id INTEGER PRIMARY KEY,",
        "    step_number INTEGER NOT NULL UNIQUE,",
        "    step_name TEXT NOT NULL,",
        "    script_name TEXT,",
        "    description TEXT NOT NULL",
        ");",
        "",
        "CREATE TABLE step_dependencies (",
        "    upstream_step_id INTEGER NOT NULL,",
        "    downstream_step_id INTEGER NOT NULL,",
        "    dependency_description TEXT NOT NULL,",
        "    PRIMARY KEY (upstream_step_id, downstream_step_id),",
        "    FOREIGN KEY (upstream_step_id) REFERENCES workflow_steps(step_id),",
        "    FOREIGN KEY (downstream_step_id) REFERENCES workflow_steps(step_id)",
        ");",
        "",
        "CREATE TABLE artifacts (",
        "    artifact_id INTEGER PRIMARY KEY,",
        "    run_id INTEGER NOT NULL,",
        "    relative_path TEXT NOT NULL,",
        "    artifact_name TEXT NOT NULL,",
        "    artifact_role TEXT NOT NULL,",
        "    artifact_type TEXT NOT NULL,",
        "    extension TEXT,",
        "    size_bytes INTEGER,",
        "    sha256 TEXT,",
        "    modified_at TEXT,",
        "    producing_step_id INTEGER,",
        "    is_external INTEGER NOT NULL DEFAULT 0 CHECK (is_external IN (0, 1)),",
        "    UNIQUE (run_id, relative_path),",
        "    FOREIGN KEY (run_id) REFERENCES workflow_runs(run_id),",
        "    FOREIGN KEY (producing_step_id) REFERENCES workflow_steps(step_id)",
        ");",
        "",
        "CREATE TABLE step_artifacts (",
        "    step_id INTEGER NOT NULL,",
        "    artifact_id INTEGER NOT NULL,",
        "    io_role TEXT NOT NULL CHECK (io_role IN ('input', 'output', 'script')),",
        "    PRIMARY KEY (step_id, artifact_id, io_role),",
        "    FOREIGN KEY (step_id) REFERENCES workflow_steps(step_id),",
        "    FOREIGN KEY (artifact_id) REFERENCES artifacts(artifact_id)",
        ");",
        "",
        "CREATE TABLE artifact_dependencies (",
        "    parent_artifact_id INTEGER NOT NULL,",
        "    child_artifact_id INTEGER NOT NULL,",
        "    relation TEXT NOT NULL,",
        "    PRIMARY KEY (parent_artifact_id, child_artifact_id, relation),",
        "    FOREIGN KEY (parent_artifact_id) REFERENCES artifacts(artifact_id),",
        "    FOREIGN KEY (child_artifact_id) REFERENCES artifacts(artifact_id)",
        ");",
        "",
    ])

    run_name = f"{Path(args.working_input).stem}_{args.started_at.replace(':', '').replace('-', '')}"
    lines.append(
        "INSERT INTO workflow_runs VALUES (" + ", ".join([
            "1",
            sql_value(run_name),
            sql_value(Path(args.working_input).stem),
            sql_value(str(original_input)),
            sql_value(str(run_dir)),
            sql_value(args.max_k),
            sql_value(args.statistical_test),
            sql_value(args.started_at),
            sql_value(args.completed_at),
            sql_value("completed"),
        ]) + ");"
    )
    lines.append("")

    for step in STEPS:
        script_name = script_for_step(step.number, args.statistical_test)
        lines.append(
            "INSERT INTO workflow_steps VALUES (" + ", ".join([
                sql_value(step.step_id),
                sql_value(step.number),
                sql_value(step.name),
                sql_value(script_name),
                sql_value(step.description),
            ]) + ");"
        )
    lines.append("")

    for upstream_number, downstream_number, description in STEP_DEPENDENCIES:
        lines.append(
            "INSERT INTO step_dependencies VALUES (" + ", ".join([
                sql_value(step_by_number[upstream_number].step_id),
                sql_value(step_by_number[downstream_number].step_id),
                sql_value(description),
            ]) + ");"
        )
    lines.append("")

    for artifact in artifacts:
        lines.append(
            "INSERT INTO artifacts VALUES (" + ", ".join([
                sql_value(artifact.artifact_id),
                "1",
                sql_value(artifact.relative_path),
                sql_value(artifact.artifact_name),
                sql_value(artifact.artifact_role),
                sql_value(artifact.artifact_type),
                sql_value(artifact.extension),
                sql_value(artifact.size_bytes),
                sql_value(artifact.sha256),
                sql_value(artifact.modified_at),
                sql_value(artifact.producing_step_id),
                sql_value(artifact.is_external),
            ]) + ");"
        )
    lines.append("")

    artifact_by_path = {artifact.relative_path: artifact for artifact in artifacts}

    # Link scripts to their execution steps.
    for step in STEPS:
        script_name = script_for_step(step.number, args.statistical_test)
        if not script_name:
            continue
        script_path = "scripts/" + script_name
        script_artifact = artifact_by_path.get(script_path)
        if script_artifact:
            lines.append(
                f"INSERT OR IGNORE INTO step_artifacts VALUES ({step.step_id}, {script_artifact.artifact_id}, 'script');"
            )

    # Link data inputs and outputs to each step.
    for step in STEPS:
        for artifact in choose_step_inputs(artifacts, step.number, args.working_input):
            lines.append(
                f"INSERT OR IGNORE INTO step_artifacts VALUES ({step.step_id}, {artifact.artifact_id}, 'input');"
            )
        for artifact in artifacts_by_step(artifacts, step.number):
            if artifact.artifact_role == "script":
                continue
            lines.append(
                f"INSERT OR IGNORE INTO step_artifacts VALUES ({step.step_id}, {artifact.artifact_id}, 'output');"
            )
    lines.append("")

    # Explicit copy relationship between the source dataset and staged dataset.
    external_input = next((a for a in artifacts if a.artifact_role == "external_input"), None)
    staged_input = next((a for a in artifacts if a.relative_path == args.working_input), None)
    if external_input and staged_input:
        lines.append(
            f"INSERT OR IGNORE INTO artifact_dependencies VALUES ({external_input.artifact_id}, {staged_input.artifact_id}, 'copied_as');"
        )

    # Direct file lineage for every output. Scripts are included as generators;
    # selected input tables/manifests are included as data dependencies.
    for step in STEPS:
        inputs = choose_step_inputs(artifacts, step.number, args.working_input)
        script_name = script_for_step(step.number, args.statistical_test)
        script_artifact = artifact_by_path.get("scripts/" + script_name) if script_name else None
        outputs = [a for a in artifacts_by_step(artifacts, step.number) if a.artifact_role != "script"]
        for child in outputs:
            for parent in inputs:
                if parent.artifact_id != child.artifact_id:
                    lines.append(
                        f"INSERT OR IGNORE INTO artifact_dependencies VALUES ({parent.artifact_id}, {child.artifact_id}, 'used_to_generate');"
                    )
            if script_artifact and script_artifact.artifact_id != child.artifact_id:
                lines.append(
                    f"INSERT OR IGNORE INTO artifact_dependencies VALUES ({script_artifact.artifact_id}, {child.artifact_id}, 'generated_by_script');"
                )
    lines.extend([
        "",
        "CREATE VIEW artifact_connections AS",
        "SELECT",
        "    p.artifact_id AS parent_artifact_id,",
        "    p.relative_path AS parent_path,",
        "    ps.step_number AS parent_step_number,",
        "    ps.step_name AS parent_step_name,",
        "    d.relation,",
        "    c.artifact_id AS child_artifact_id,",
        "    c.relative_path AS child_path,",
        "    cs.step_number AS child_step_number,",
        "    cs.step_name AS child_step_name",
        "FROM artifact_dependencies d",
        "JOIN artifacts p ON p.artifact_id = d.parent_artifact_id",
        "JOIN artifacts c ON c.artifact_id = d.child_artifact_id",
        "LEFT JOIN workflow_steps ps ON ps.step_id = p.producing_step_id",
        "LEFT JOIN workflow_steps cs ON cs.step_id = c.producing_step_id;",
        "",
        "CREATE VIEW step_file_inventory AS",
        "SELECT",
        "    s.step_number,",
        "    s.step_name,",
        "    sa.io_role,",
        "    a.relative_path,",
        "    a.artifact_type,",
        "    a.size_bytes,",
        "    a.sha256",
        "FROM step_artifacts sa",
        "JOIN workflow_steps s ON s.step_id = sa.step_id",
        "JOIN artifacts a ON a.artifact_id = sa.artifact_id;",
        "",
        "CREATE INDEX idx_artifacts_step ON artifacts(producing_step_id);",
        "CREATE INDEX idx_artifacts_path ON artifacts(relative_path);",
        "CREATE INDEX idx_dependencies_parent ON artifact_dependencies(parent_artifact_id);",
        "CREATE INDEX idx_dependencies_child ON artifact_dependencies(child_artifact_id);",
        "",
        "COMMIT;",
        "",
        "-- Example query 1: list every workflow output with its producing step.",
        "-- SELECT step_number, step_name, relative_path, artifact_type",
        "-- FROM step_file_inventory",
        "-- WHERE io_role = 'output'",
        "-- ORDER BY step_number, relative_path;",
        "",
        "-- Example query 2: direct inputs for the final network edge lists.",
        "-- SELECT parent_path, relation, child_path",
        "-- FROM artifact_connections",
        "-- WHERE child_path LIKE '07_network_construction/edge_lists/%'",
        "-- ORDER BY child_path, parent_path;",
        "",
        "-- Example query 3: recursively trace all upstream files for one output.",
        "-- WITH RECURSIVE lineage(parent_id, child_id, relation, depth) AS (",
        "--     SELECT parent_artifact_id, child_artifact_id, relation, 1",
        "--     FROM artifact_dependencies",
        "--     WHERE child_artifact_id = (",
        "--         SELECT artifact_id FROM artifacts",
        "--         WHERE relative_path = '07_network_construction/01_edge_list_merge_summary.csv'",
        "--     )",
        "--     UNION ALL",
        "--     SELECT d.parent_artifact_id, d.child_artifact_id, d.relation, l.depth + 1",
        "--     FROM artifact_dependencies d",
        "--     JOIN lineage l ON d.child_artifact_id = l.parent_id",
        "-- )",
        "-- SELECT DISTINCT a.relative_path, lineage.relation, lineage.depth",
        "-- FROM lineage",
        "-- JOIN artifacts a ON a.artifact_id = lineage.parent_id",
        "-- ORDER BY lineage.depth, a.relative_path;",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = run_dir / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sql_text = generate_sql(args)
    output_path.write_text(sql_text, encoding="utf-8")
    print(f"Workflow provenance SQL written to: {output_path}")


if __name__ == "__main__":
    main()
NICHE_SQL_PY

if [ ! -f "$PROVENANCE_SQL" ]; then
    echo "Error: SQL provenance file was not created:"
    echo "  $OUTPUT_DIR/$PROVENANCE_SQL"
    exit 1
fi

# Materialize the provenance SQL as a SQLite database without requiring the
# sqlite3 command-line application.
python - "$PROVENANCE_SQL" "$PROVENANCE_DB" <<'NICHE_CREATE_DB_PY'
from pathlib import Path
import sqlite3
import sys

sql_path = Path(sys.argv[1])
db_path = Path(sys.argv[2])
db_path.parent.mkdir(parents=True, exist_ok=True)
if db_path.exists():
    db_path.unlink()

connection = sqlite3.connect(db_path)
try:
    connection.execute("PRAGMA foreign_keys = ON")
    connection.executescript(sql_path.read_text(encoding="utf-8"))
    connection.commit()
finally:
    connection.close()

print(f"SQLite provenance database created: {db_path}")
NICHE_CREATE_DB_PY

# Import the actual contents of every CSV, TSV, and JSON artifact into the
# database. A temporary helper is generated from this self-contained runner
# and removed after use.
RESULT_IMPORTER="$PROVENANCE_DIR/.niches_import_results.py"
cat > "$RESULT_IMPORTER" <<'NICHE_RESULT_IMPORTER_PY'
"""Import NICHEs CSV/TSV/JSON artifacts into a provenance SQLite database.

The script augments an existing NICHEs provenance database with:
- one populated SQLite table per CSV/TSV/JSON artifact;
- metadata describing each imported table and its columns;
- table-level lineage derived from artifact_dependencies;
- foreign-key links from every imported result table to dataset_tables.

It can also write a full SQL dump that recreates both provenance metadata and
all imported scientific result data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence


RESERVED_METADATA_COLUMNS = {
    "__niches_row_id",
    "__dataset_table_id",
    "__source_row_number",
}


@dataclass(frozen=True)
class ArtifactRecord:
    artifact_id: int
    relative_path: str
    artifact_name: str
    artifact_role: str
    artifact_type: str
    extension: str | None
    producing_step_id: int | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="Dataset-specific NICHEs output directory.")
    parser.add_argument("--db", required=True, help="Existing provenance SQLite database to augment.")
    parser.add_argument(
        "--dump-sql",
        help="Optional full SQL dump path. The dump includes all imported result rows.",
    )
    parser.add_argument(
        "--replace-imports",
        action="store_true",
        help="Drop previously imported scientific-result tables before importing again.",
    )
    parser.add_argument(
        "--include-json",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Import JSON manifests and summaries as flattened tables (default: yes).",
    )
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def slugify(text: str, max_length: int = 48) -> str:
    slug = re.sub(r"[^0-9A-Za-z_]+", "_", text).strip("_").lower()
    slug = re.sub(r"_+", "_", slug)
    if not slug:
        slug = "result"
    if slug[0].isdigit():
        slug = "t_" + slug
    return slug[:max_length].rstrip("_") or "result"


def make_table_name(artifact: ArtifactRecord, existing: set[str]) -> str:
    step = artifact.producing_step_id if artifact.producing_step_id is not None else 0
    base = slugify(Path(artifact.artifact_name).stem)
    digest = hashlib.sha1(artifact.relative_path.encode("utf-8")).hexdigest()[:8]
    candidate = f"result_s{step:02d}_{base}_{digest}"
    suffix = 2
    while candidate in existing:
        candidate = f"result_s{step:02d}_{base}_{digest}_{suffix}"
        suffix += 1
    existing.add(candidate)
    return candidate


def make_unique_column_names(source_names: Sequence[str]) -> list[str]:
    used = set(RESERVED_METADATA_COLUMNS)
    result: list[str] = []
    for index, source_name in enumerate(source_names, start=1):
        base = slugify(str(source_name or f"column_{index}"), max_length=54)
        if base in RESERVED_METADATA_COLUMNS:
            base = "data_" + base.lstrip("_")
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used.add(candidate)
        result.append(candidate)
    return result


def resolve_artifact_path(run_dir: Path, relative_path: str) -> Path:
    path = Path(relative_path)
    return path if path.is_absolute() else run_dir / path


def detect_delimiter(path: Path, extension: str | None) -> str:
    if (extension or path.suffix).lower() == ".tsv":
        return "\t"
    try:
        sample = path.read_text(encoding="utf-8-sig", errors="replace")[:65536]
        if not sample:
            return ","
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        return dialect.delimiter
    except (csv.Error, OSError):
        return ","


def value_kind(value: str) -> str | None:
    text = value.strip()
    if text == "":
        return None
    if re.fullmatch(r"[+-]?\d+", text):
        return "INTEGER"
    try:
        float(text)
        return "REAL"
    except ValueError:
        return "TEXT"


def merge_sqlite_type(current: str | None, observed: str | None) -> str | None:
    if observed is None:
        return current
    if current is None:
        return observed
    if current == observed:
        return current
    if {current, observed} <= {"INTEGER", "REAL"}:
        return "REAL"
    return "TEXT"


def infer_types(rows: Sequence[Sequence[str]], column_count: int) -> list[str]:
    inferred: list[str | None] = [None] * column_count
    for row in rows:
        padded = list(row[:column_count]) + [""] * max(0, column_count - len(row))
        for index, value in enumerate(padded):
            inferred[index] = merge_sqlite_type(inferred[index], value_kind(value))
    return [item or "TEXT" for item in inferred]


def convert_value(value: str, sqlite_type: str) -> Any:
    text = value.strip()
    if text == "":
        return None
    if sqlite_type == "INTEGER":
        try:
            return int(text)
        except ValueError:
            return value
    if sqlite_type == "REAL":
        try:
            return float(text)
        except ValueError:
            return value
    return value


def batched(iterator: Iterable[tuple[Any, ...]], size: int = 1000) -> Iterator[list[tuple[Any, ...]]]:
    batch: list[tuple[Any, ...]] = []
    for item in iterator:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


def ensure_catalog_schema(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS dataset_tables (
            dataset_table_id INTEGER PRIMARY KEY,
            artifact_id INTEGER NOT NULL UNIQUE,
            sqlite_table_name TEXT NOT NULL UNIQUE,
            display_name TEXT NOT NULL,
            source_path TEXT NOT NULL,
            source_format TEXT NOT NULL,
            producing_step_id INTEGER,
            row_count INTEGER NOT NULL DEFAULT 0,
            column_count INTEGER NOT NULL DEFAULT 0,
            import_status TEXT NOT NULL,
            import_error TEXT,
            imported_at TEXT NOT NULL,
            FOREIGN KEY (artifact_id) REFERENCES artifacts(artifact_id),
            FOREIGN KEY (producing_step_id) REFERENCES workflow_steps(step_id)
        );

        CREATE TABLE IF NOT EXISTS dataset_columns (
            dataset_table_id INTEGER NOT NULL,
            ordinal_position INTEGER NOT NULL,
            source_column_name TEXT NOT NULL,
            sqlite_column_name TEXT NOT NULL,
            sqlite_type TEXT NOT NULL,
            PRIMARY KEY (dataset_table_id, ordinal_position),
            UNIQUE (dataset_table_id, sqlite_column_name),
            FOREIGN KEY (dataset_table_id) REFERENCES dataset_tables(dataset_table_id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS table_dependencies (
            parent_dataset_table_id INTEGER NOT NULL,
            child_dataset_table_id INTEGER NOT NULL,
            relation TEXT NOT NULL,
            PRIMARY KEY (parent_dataset_table_id, child_dataset_table_id, relation),
            FOREIGN KEY (parent_dataset_table_id) REFERENCES dataset_tables(dataset_table_id),
            FOREIGN KEY (child_dataset_table_id) REFERENCES dataset_tables(dataset_table_id)
        );

        CREATE INDEX IF NOT EXISTS idx_dataset_tables_artifact
            ON dataset_tables(artifact_id);
        CREATE INDEX IF NOT EXISTS idx_dataset_tables_step
            ON dataset_tables(producing_step_id);
        CREATE INDEX IF NOT EXISTS idx_table_dependencies_parent
            ON table_dependencies(parent_dataset_table_id);
        CREATE INDEX IF NOT EXISTS idx_table_dependencies_child
            ON table_dependencies(child_dataset_table_id);
        """
    )

    connection.execute("DROP VIEW IF EXISTS scientific_table_inventory")
    connection.execute("DROP VIEW IF EXISTS scientific_table_connections")
    connection.executescript(
        """
        CREATE VIEW scientific_table_inventory AS
        SELECT
            dt.dataset_table_id,
            dt.sqlite_table_name,
            dt.display_name,
            dt.source_path,
            dt.source_format,
            ws.step_number,
            ws.step_name,
            dt.row_count,
            dt.column_count,
            dt.import_status,
            dt.import_error,
            dt.imported_at
        FROM dataset_tables dt
        LEFT JOIN workflow_steps ws
            ON ws.step_id = dt.producing_step_id;

        CREATE VIEW scientific_table_connections AS
        SELECT
            p.dataset_table_id AS parent_dataset_table_id,
            p.sqlite_table_name AS parent_sqlite_table,
            p.source_path AS parent_source_path,
            td.relation,
            c.dataset_table_id AS child_dataset_table_id,
            c.sqlite_table_name AS child_sqlite_table,
            c.source_path AS child_source_path
        FROM table_dependencies td
        JOIN dataset_tables p
            ON p.dataset_table_id = td.parent_dataset_table_id
        JOIN dataset_tables c
            ON c.dataset_table_id = td.child_dataset_table_id;
        """
    )


def existing_imported_table_names(connection: sqlite3.Connection) -> set[str]:
    return {
        row[0]
        for row in connection.execute(
            "SELECT sqlite_table_name FROM dataset_tables WHERE import_status = 'imported'"
        )
    }


def drop_previous_imports(connection: sqlite3.Connection) -> None:
    names = [row[0] for row in connection.execute("SELECT sqlite_table_name FROM dataset_tables")]
    connection.execute("PRAGMA foreign_keys = OFF")
    try:
        for name in names:
            connection.execute(f"DROP TABLE IF EXISTS {quote_identifier(name)}")
        connection.execute("DELETE FROM table_dependencies")
        connection.execute("DELETE FROM dataset_columns")
        connection.execute("DELETE FROM dataset_tables")
        connection.commit()
    finally:
        connection.execute("PRAGMA foreign_keys = ON")


def list_importable_artifacts(connection: sqlite3.Connection, include_json: bool) -> list[ArtifactRecord]:
    extensions = [".csv", ".tsv"] + ([".json"] if include_json else [])
    placeholders = ",".join("?" for _ in extensions)
    query = f"""
        SELECT
            artifact_id,
            relative_path,
            artifact_name,
            artifact_role,
            artifact_type,
            extension,
            producing_step_id
        FROM artifacts
        WHERE lower(COALESCE(extension, '')) IN ({placeholders})
          AND artifact_role != 'script'
        ORDER BY artifact_id
    """
    return [ArtifactRecord(*row) for row in connection.execute(query, extensions)]


def register_dataset_table(
    connection: sqlite3.Connection,
    artifact: ArtifactRecord,
    sqlite_table_name: str,
    source_format: str,
    status: str = "pending",
    error: str | None = None,
) -> int:
    cursor = connection.execute(
        """
        INSERT INTO dataset_tables (
            artifact_id,
            sqlite_table_name,
            display_name,
            source_path,
            source_format,
            producing_step_id,
            row_count,
            column_count,
            import_status,
            import_error,
            imported_at
        ) VALUES (?, ?, ?, ?, ?, ?, 0, 0, ?, ?, ?)
        ON CONFLICT(artifact_id) DO UPDATE SET
            sqlite_table_name = excluded.sqlite_table_name,
            display_name = excluded.display_name,
            source_path = excluded.source_path,
            source_format = excluded.source_format,
            producing_step_id = excluded.producing_step_id,
            import_status = excluded.import_status,
            import_error = excluded.import_error,
            imported_at = excluded.imported_at
        RETURNING dataset_table_id
        """,
        (
            artifact.artifact_id,
            sqlite_table_name,
            artifact.artifact_name,
            artifact.relative_path,
            source_format,
            artifact.producing_step_id,
            status,
            error,
            utc_now(),
        ),
    )
    row = cursor.fetchone()
    if row is None:
        raise RuntimeError(f"Could not register dataset table for artifact {artifact.artifact_id}")
    return int(row[0])


def create_result_table(
    connection: sqlite3.Connection,
    table_name: str,
    dataset_table_id: int,
    sqlite_columns: Sequence[str],
    sqlite_types: Sequence[str],
) -> None:
    connection.execute(f"DROP TABLE IF EXISTS {quote_identifier(table_name)}")
    definitions = [
        '"__niches_row_id" INTEGER PRIMARY KEY AUTOINCREMENT',
        (
            '"__dataset_table_id" INTEGER NOT NULL '
            f'DEFAULT {dataset_table_id} CHECK ("__dataset_table_id" = {dataset_table_id}) '
            'REFERENCES dataset_tables(dataset_table_id)'
        ),
        '"__source_row_number" INTEGER NOT NULL',
    ]
    definitions.extend(
        f"{quote_identifier(name)} {sqlite_type}"
        for name, sqlite_type in zip(sqlite_columns, sqlite_types)
    )
    sql = f"CREATE TABLE {quote_identifier(table_name)} (\n    " + ",\n    ".join(definitions) + "\n)"
    connection.execute(sql)
    connection.execute(
        f"CREATE INDEX {quote_identifier('idx_' + table_name + '_dataset')} "
        f"ON {quote_identifier(table_name)}(\"__dataset_table_id\")"
    )


def import_delimited_file(
    connection: sqlite3.Connection,
    run_dir: Path,
    artifact: ArtifactRecord,
    table_name: str,
    dataset_table_id: int,
) -> tuple[int, int]:
    path = resolve_artifact_path(run_dir, artifact.relative_path)
    delimiter = detect_delimiter(path, artifact.extension)

    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        try:
            header = next(reader)
        except StopIteration:
            create_result_table(connection, table_name, dataset_table_id, [], [])
            return 0, 0

        source_columns = [str(item) for item in header]
        sqlite_columns = make_unique_column_names(source_columns)

        sample_rows: list[list[str]] = []
        for _ in range(5000):
            try:
                sample_rows.append(next(reader))
            except StopIteration:
                break
        sqlite_types = infer_types(sample_rows, len(source_columns))

        create_result_table(
            connection,
            table_name,
            dataset_table_id,
            sqlite_columns,
            sqlite_types,
        )

        connection.execute("DELETE FROM dataset_columns WHERE dataset_table_id = ?", (dataset_table_id,))
        connection.executemany(
            """
            INSERT INTO dataset_columns (
                dataset_table_id,
                ordinal_position,
                source_column_name,
                sqlite_column_name,
                sqlite_type
            ) VALUES (?, ?, ?, ?, ?)
            """,
            [
                (dataset_table_id, index, source, sqlite_name, sqlite_type)
                for index, (source, sqlite_name, sqlite_type) in enumerate(
                    zip(source_columns, sqlite_columns, sqlite_types), start=1
                )
            ],
        )

        insert_columns = ["__dataset_table_id", "__source_row_number", *sqlite_columns]
        placeholders = ",".join("?" for _ in insert_columns)
        insert_sql = (
            f"INSERT INTO {quote_identifier(table_name)} "
            f"({', '.join(quote_identifier(item) for item in insert_columns)}) "
            f"VALUES ({placeholders})"
        )

        def normalized_rows() -> Iterator[tuple[Any, ...]]:
            row_number = 0
            for row in [*sample_rows, *reader]:
                row_number += 1
                values = list(row[: len(source_columns)])
                if len(values) < len(source_columns):
                    values.extend([""] * (len(source_columns) - len(values)))
                converted = [
                    convert_value(value, sqlite_type)
                    for value, sqlite_type in zip(values, sqlite_types)
                ]
                yield (dataset_table_id, row_number, *converted)

        row_count = 0
        for batch in batched(normalized_rows(), size=1000):
            connection.executemany(insert_sql, batch)
            row_count += len(batch)

    return row_count, len(source_columns)


def flatten_json(value: Any, path: str = "$") -> Iterator[tuple[str, str | None, str]]:
    if isinstance(value, dict):
        if not value:
            yield path, "{}", "object"
        for key, child in value.items():
            escaped = str(key).replace("~", "~0").replace("/", "~1")
            yield from flatten_json(child, f"{path}/{escaped}")
    elif isinstance(value, list):
        if not value:
            yield path, "[]", "array"
        for index, child in enumerate(value):
            yield from flatten_json(child, f"{path}/{index}")
    elif value is None:
        yield path, None, "null"
    elif isinstance(value, bool):
        yield path, "true" if value else "false", "boolean"
    elif isinstance(value, int):
        yield path, str(value), "integer"
    elif isinstance(value, float):
        yield path, repr(value), "real"
    else:
        yield path, str(value), "text"


def import_json_file(
    connection: sqlite3.Connection,
    run_dir: Path,
    artifact: ArtifactRecord,
    table_name: str,
    dataset_table_id: int,
) -> tuple[int, int]:
    path = resolve_artifact_path(run_dir, artifact.relative_path)
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        document = json.load(handle)

    source_columns = ["json_path", "value", "value_type"]
    sqlite_columns = source_columns[:]
    sqlite_types = ["TEXT", "TEXT", "TEXT"]
    create_result_table(connection, table_name, dataset_table_id, sqlite_columns, sqlite_types)

    connection.execute("DELETE FROM dataset_columns WHERE dataset_table_id = ?", (dataset_table_id,))
    connection.executemany(
        """
        INSERT INTO dataset_columns (
            dataset_table_id,
            ordinal_position,
            source_column_name,
            sqlite_column_name,
            sqlite_type
        ) VALUES (?, ?, ?, ?, ?)
        """,
        [
            (dataset_table_id, index, source, sqlite_name, sqlite_type)
            for index, (source, sqlite_name, sqlite_type) in enumerate(
                zip(source_columns, sqlite_columns, sqlite_types), start=1
            )
        ],
    )

    insert_sql = f"""
        INSERT INTO {quote_identifier(table_name)} (
            "__dataset_table_id",
            "__source_row_number",
            "json_path",
            "value",
            "value_type"
        ) VALUES (?, ?, ?, ?, ?)
    """
    rows = [
        (dataset_table_id, row_number, json_path, value, value_type)
        for row_number, (json_path, value, value_type) in enumerate(flatten_json(document), start=1)
    ]
    connection.executemany(insert_sql, rows)
    return len(rows), 3


def update_import_result(
    connection: sqlite3.Connection,
    dataset_table_id: int,
    row_count: int,
    column_count: int,
    status: str,
    error: str | None = None,
) -> None:
    connection.execute(
        """
        UPDATE dataset_tables
        SET row_count = ?,
            column_count = ?,
            import_status = ?,
            import_error = ?,
            imported_at = ?
        WHERE dataset_table_id = ?
        """,
        (row_count, column_count, status, error, utc_now(), dataset_table_id),
    )


def rebuild_table_dependencies(connection: sqlite3.Connection) -> int:
    connection.execute("DELETE FROM table_dependencies")
    connection.execute(
        """
        INSERT OR IGNORE INTO table_dependencies (
            parent_dataset_table_id,
            child_dataset_table_id,
            relation
        )
        SELECT
            parent.dataset_table_id,
            child.dataset_table_id,
            dependencies.relation
        FROM artifact_dependencies AS dependencies
        JOIN dataset_tables AS parent
            ON parent.artifact_id = dependencies.parent_artifact_id
        JOIN dataset_tables AS child
            ON child.artifact_id = dependencies.child_artifact_id
        WHERE parent.import_status = 'imported'
          AND child.import_status = 'imported'
          AND parent.dataset_table_id != child.dataset_table_id
        """
    )
    return int(connection.execute("SELECT COUNT(*) FROM table_dependencies").fetchone()[0])


def write_full_sql_dump(connection: sqlite3.Connection, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("-- Complete NICHEs SQLite export: provenance plus scientific result data\n")
        handle.write(f"-- Generated: {utc_now()}\n\n")
        for line in connection.iterdump():
            handle.write(line)
            handle.write("\n")
    temporary.replace(destination)


def validate_provenance_schema(connection: sqlite3.Connection) -> None:
    required = {"artifacts", "artifact_dependencies", "workflow_steps"}
    existing = {
        row[0]
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }
    missing = required - existing
    if missing:
        raise RuntimeError(
            "The database is not a compatible NICHEs provenance database. "
            f"Missing tables: {', '.join(sorted(missing))}"
        )


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    db_path = Path(args.db)
    if not db_path.is_absolute():
        db_path = run_dir / db_path
    dump_path = None
    if args.dump_sql:
        dump_path = Path(args.dump_sql)
        if not dump_path.is_absolute():
            dump_path = run_dir / dump_path

    if not run_dir.is_dir():
        raise SystemExit(f"Run directory does not exist: {run_dir}")
    if not db_path.is_file():
        raise SystemExit(f"Provenance database does not exist: {db_path}")

    connection = sqlite3.connect(db_path)
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA journal_mode = DELETE")
    connection.execute("PRAGMA synchronous = NORMAL")

    imported = 0
    missing = 0
    failed = 0
    total_rows = 0

    try:
        validate_provenance_schema(connection)
        ensure_catalog_schema(connection)
        if args.replace_imports:
            drop_previous_imports(connection)
            ensure_catalog_schema(connection)

        existing_names = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type IN ('table', 'view')"
            )
        }
        existing_names.update(existing_imported_table_names(connection))

        artifacts = list_importable_artifacts(connection, args.include_json)
        print(f"Found {len(artifacts)} importable CSV/TSV/JSON artifacts.")

        for artifact in artifacts:
            path = resolve_artifact_path(run_dir, artifact.relative_path)
            previous = connection.execute(
                "SELECT dataset_table_id, sqlite_table_name FROM dataset_tables WHERE artifact_id = ?",
                (artifact.artifact_id,),
            ).fetchone()
            if previous and not args.replace_imports:
                table_name = str(previous[1])
                dataset_table_id = int(previous[0])
                table_exists = connection.execute(
                    "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
                    (table_name,),
                ).fetchone()
                if table_exists:
                    print(f"SKIP    {artifact.relative_path} -> {table_name}")
                    continue
            else:
                table_name = make_table_name(artifact, existing_names)
                dataset_table_id = register_dataset_table(
                    connection,
                    artifact,
                    table_name,
                    (artifact.extension or path.suffix).lower().lstrip("."),
                )

            if not path.is_file():
                update_import_result(
                    connection,
                    dataset_table_id,
                    0,
                    0,
                    "missing_source",
                    f"Source file not found: {path}",
                )
                missing += 1
                print(f"MISSING {artifact.relative_path}")
                continue

            try:
                extension = (artifact.extension or path.suffix).lower()
                if extension in {".csv", ".tsv"}:
                    row_count, column_count = import_delimited_file(
                        connection,
                        run_dir,
                        artifact,
                        table_name,
                        dataset_table_id,
                    )
                elif extension == ".json":
                    row_count, column_count = import_json_file(
                        connection,
                        run_dir,
                        artifact,
                        table_name,
                        dataset_table_id,
                    )
                else:
                    continue
                update_import_result(
                    connection,
                    dataset_table_id,
                    row_count,
                    column_count,
                    "imported",
                )
                connection.commit()
                imported += 1
                total_rows += row_count
                print(
                    f"IMPORTED {artifact.relative_path} -> {table_name} "
                    f"({row_count} rows, {column_count} columns)"
                )
            except Exception as exc:  # Continue so one malformed file does not erase the run database.
                connection.rollback()
                ensure_catalog_schema(connection)
                # Re-register because a rollback may have removed a pending catalog row.
                dataset_table_id = register_dataset_table(
                    connection,
                    artifact,
                    table_name,
                    (artifact.extension or path.suffix).lower().lstrip("."),
                    status="failed",
                    error=str(exc),
                )
                update_import_result(connection, dataset_table_id, 0, 0, "failed", str(exc))
                connection.commit()
                failed += 1
                print(f"FAILED  {artifact.relative_path}: {exc}", file=sys.stderr)

        dependency_count = rebuild_table_dependencies(connection)
        connection.commit()

        # Record a compact import summary in the database itself.
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS result_import_runs (
                import_run_id INTEGER PRIMARY KEY,
                imported_at TEXT NOT NULL,
                run_directory TEXT NOT NULL,
                imported_tables INTEGER NOT NULL,
                missing_sources INTEGER NOT NULL,
                failed_tables INTEGER NOT NULL,
                imported_rows INTEGER NOT NULL,
                table_dependencies INTEGER NOT NULL
            )
            """
        )
        connection.execute(
            """
            INSERT INTO result_import_runs (
                imported_at,
                run_directory,
                imported_tables,
                missing_sources,
                failed_tables,
                imported_rows,
                table_dependencies
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (utc_now(), str(run_dir), imported, missing, failed, total_rows, dependency_count),
        )
        connection.commit()

        if dump_path is not None:
            write_full_sql_dump(connection, dump_path)
            print(f"Full SQL dump written to: {dump_path}")

        print("\nImport summary")
        print(f"  Imported tables:   {imported}")
        print(f"  Imported rows:     {total_rows}")
        print(f"  Missing sources:   {missing}")
        print(f"  Failed imports:    {failed}")
        print(f"  Table connections: {dependency_count}")
        print(f"  SQLite database:   {db_path}")
        return 0 if failed == 0 else 2
    finally:
        connection.close()


if __name__ == "__main__":
    raise SystemExit(main())

NICHE_RESULT_IMPORTER_PY

python "$RESULT_IMPORTER" \
    --run-dir "$OUTPUT_DIR" \
    --db "$PROVENANCE_DB" \
    --dump-sql "$PROVENANCE_SQL" \
    --replace-imports

rm -f "$RESULT_IMPORTER"

if [ ! -f "$PROVENANCE_DB" ]; then
    echo "Error: Populated SQLite database was not created:"
    echo "  $OUTPUT_DIR/$PROVENANCE_DB"
    exit 1
fi

# Validate that real scientific-result tables were imported.
python - "$PROVENANCE_DB" <<'NICHE_VALIDATE_DB_PY'
import sqlite3
import sys

connection = sqlite3.connect(sys.argv[1])
try:
    imported_tables = connection.execute(
        "SELECT COUNT(*) FROM dataset_tables WHERE import_status = 'imported'"
    ).fetchone()[0]
    imported_rows = connection.execute(
        "SELECT COALESCE(SUM(row_count), 0) FROM dataset_tables WHERE import_status = 'imported'"
    ).fetchone()[0]
    table_links = connection.execute(
        "SELECT COUNT(*) FROM table_dependencies"
    ).fetchone()[0]
finally:
    connection.close()

if imported_tables == 0:
    raise SystemExit("Error: no scientific result tables were imported into SQLite.")

print(f"Imported scientific tables: {imported_tables}")
print(f"Imported scientific rows:   {imported_rows}")
print(f"Table-level connections:    {table_links}")
NICHE_VALIDATE_DB_PY

echo ""
echo "Step 9 complete."
echo "Populated SQLite database:"
echo "  $OUTPUT_DIR/$PROVENANCE_DB"
echo ""
echo "Full SQL reconstruction dump:"
echo "  $OUTPUT_DIR/$PROVENANCE_SQL"
echo ""
echo "The database now contains:"
echo "  - provenance and workflow-step tables"
echo "  - one populated SQLite table for every CSV/TSV/JSON result"
echo "  - dataset_tables and dataset_columns catalogs"
echo "  - table_dependencies describing table-to-table lineage"
echo "  - scientific_table_inventory and scientific_table_connections views"
echo ""

echo "============================================================"
echo "NICHEs workflow completed successfully."
echo "============================================================"
echo "Input file:"
echo "  $ABS_INPUT_FILE"
echo ""
echo "All outputs were saved in:"
echo "  $OUTPUT_DIR"
echo ""
echo "Final network outputs:"
echo "  $OUTPUT_DIR/07_network_construction/edge_lists/"
echo "  $OUTPUT_DIR/07_network_construction/nodes_attributes/"
echo ""
echo "Populated workflow SQLite database:"
echo "  $OUTPUT_DIR/$PROVENANCE_DB"
echo ""
echo "Full workflow SQL dump:"
echo "  $OUTPUT_DIR/$PROVENANCE_SQL"
echo "============================================================"
