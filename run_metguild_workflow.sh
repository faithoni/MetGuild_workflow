#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------
# MetGuild full workflow runner (Steps 1-9)
# ------------------------------------------------------------

WORKFLOW_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'USAGE'
MetGuild Workflow
=================

Usage:
  bash run_metguild_workflow.sh <input_matrix> [MAX_K] [OPTIONS]

By default, MetGuild runs the complete workflow (Steps 1-9).
If --taxonomy is supplied, Step 3b annotates each cluster before Step 4.

Required argument:
  <input_matrix>          UTF-8 CSV or TSV matrix containing a Genomes column
                          and numeric functional features.

Optional positional argument:
  MAX_K                   Exclusive upper limit for candidate cluster numbers
                          (default: 15). MAX_K=15 evaluates k=2,...,14,
                          subject to dataset-size constraints.

Main options:
  --matrix-type TYPE      Supported workflow matrix types:
                            binary       presence/absence values (0/1)
                            proportion   proportional values between 0 and 1
                            auto         detect exact binary data; otherwise Step 1
                                         reports quantitative and the full workflow
                                         will require an explicit supported type

  --full                  Explicitly request the complete Steps 1-9 workflow.
                          This is already the default when no stop option is used.

  --through-step4         Run Steps 1-4 and stop after statistical feature screening.

  --clustering-only       Run Steps 1-3 and stop after both clustering methods.

  --step1-only            Run only Step 1.

  --output-dir PATH       Use a specific fresh output directory. If omitted,
                          MetGuild creates <dataset>_run_01, then _run_02, etc.

  --taxonomy FILE         Optional taxonomy table containing:
                            Genomes   Phyla
                          The taxonomy column may be named Phyla or Phylum.
                          Comma-, tab-, and whitespace-delimited files are accepted.
                          When supplied, Step 3b writes one MAG_Phyla field per
                          cluster and Step 8 carries it into cluster node attributes.

  -h, --help              Show this help message and exit.

Step 1 options:
  --input-format FORMAT   auto (default), csv, or tsv.
  --distance METRIC       auto (default), jaccard, ruzicka, or euclidean.
  --candidate-source MODE union (default) or kmeans.
                          union combines the four KMeans selections with the
                          hierarchical silhouette selection.
  --mds-components N      Requested classical-MDS dimensions (default: 3).
  --n-init N              KMeans initializations (default: 100).
  --kmeans-max-iter N     Maximum KMeans iterations (default: 100).
  --random-state N        KMeans random seed (default: 0).
  --missing-values POLICY error (default) or zero.

Examples
--------

Show all usage information:

   bash run_metguild_workflow.sh --help

1. Complete workflow for a binary matrix with taxonomy annotation:

   bash run_metguild_workflow.sh \
     matrix.csv 15 \
     --matrix-type binary \
     --taxonomy taxonomy.csv

   Example taxonomy.csv content:

     Genomes       Phyla
     SRCFOC_001    Bacteroidota
     SRCFOC_002    Desulfobacterota
     SRCFOC_004    Hadarchaeota

2. Complete workflow for a proportional matrix with taxonomy annotation:

   bash run_metguild_workflow.sh \
     pathway_proportions.tsv 15 \
     --matrix-type proportion \
     --taxonomy taxonomy.csv

3. Complete workflow without taxonomy annotation:

   bash run_metguild_workflow.sh \
     matrix.csv 15 \
     --matrix-type binary

4. The same complete workflow using the explicit --full flag:

   bash run_metguild_workflow.sh \
     matrix.csv 15 \
     --full \
     --matrix-type binary \
     --taxonomy taxonomy.csv

5. Run only Steps 1-4:

   bash run_metguild_workflow.sh \
     matrix.csv 15 \
     --through-step4 \
     --matrix-type binary

6. Run only Steps 1-3:

   bash run_metguild_workflow.sh \
     matrix.csv 15 \
     --clustering-only \
     --matrix-type binary

7. Run only Step 1:

   bash run_metguild_workflow.sh \
     matrix.csv 15 \
     --step1-only \
     --matrix-type binary

Workflow stages
---------------

Step 1 - Matrix preparation and candidate-k selection
  * Validate and clean the matrix.
  * Binary data: Jaccard dissimilarity.
  * Proportional data: Ruzicka dissimilarity.
  * Classical MDS for the KMeans branch.
  * Candidate k from:
      - KMeans elbow/inertia
      - KMeans silhouette
      - Calinski-Harabasz
      - Davies-Bouldin
      - average-linkage hierarchical silhouette
  * By default, the unique union of all five selections is forwarded.

Step 2 - KMeans clustering
  * KMeans on the classical-MDS coordinates for each forwarded k.

Step 3 - Hierarchical clustering
  * Average linkage on the original Jaccard/Ruzicka dissimilarities.

Step 3b - Taxonomy annotation (only when --taxonomy is supplied)
  * Join clustered Genomes to the Phyla/Phylum column in the taxonomy table.
  * For every KMeans/HC x k solution, write:
      Cluster, MAG_Phyla
  * Example MAG_Phyla value:
      MAG001-Chloroflexota; MAG005-Desulfobacterota
  * Taxonomy annotates the clusters; it does not influence clustering.

Step 4 - Statistical feature screening
  * binary     -> Pearson chi-square
  * proportion -> one-way ANOVA
  * Benjamini-Hochberg FDR correction.

Step 5 - Feature-feature correlation
  * Spearman correlation among Step 4 retained features.
  * BH-adjusted significance and |rho| >= 0.5 by default.
  * Perfect correlations are excluded from retained network edges.

Step 6 - Feature-cluster correlation
  * One-vs-rest Spearman association between retained features and clusters.
  * Positive rho >= 0.5 and BH-adjusted significance by default.
  * Perfect correlations are excluded from retained network edges.

Step 7 - Cluster-cluster correlation
  * Binary data use cluster feature prevalence profiles.
  * Proportional data use cluster mean feature profiles.
  * Positive rho >= 0.5 and BH-adjusted significance by default.
  * Perfect correlations are excluded from retained network edges.

Step 8 - Network construction
  * Merge filtered feature-feature, feature-cluster, and cluster-cluster edges.
  * Correlation outputs use cluster names C1, C2, ... consistently.
  * Final network cluster nodes are made solution-specific, e.g. C1_KM3,
    C1_HC3. Feature names are preserved exactly.
  * If Step 3b was run, cluster node attributes include MAG_Phyla.

Step 9 - SQLite/SQL workflow provenance
  * Build a populated SQLite database linking workflow steps, scripts, inputs,
    intermediate outputs, and final network outputs.
  * Import the actual contents of CSV, TSV, and JSON artifacts into SQLite.
  * Write a complete SQL dump that can recreate the populated database.

Automatic output directory
--------------------------

If the input is matrix.csv and --output-dir is omitted:

  matrix_run_01/

A later run creates:

  matrix_run_02/

A full run is organized approximately as:

  matrix_run_01/
  |-- matrix_cleaned.tsv
  |-- 02_k_optimization/
  |-- manifests/
  |-- 03_04_clustering/
  |   |-- 01_kmeans/
  |   |-- 02_hierarchical/
  |   `-- taxonomy/              (when --taxonomy is supplied)
  |-- 05_stat_analysis/
  |-- 06_correlation/
  |   |-- 01_featurevsfeature_correlation/
  |   |-- 02_featurevscluster_correlation/
  |   `-- 03_clustervscluster_correlation/
  |-- 07_network_construction/
  |   |-- edge_lists/
  |   `-- nodes_attributes/
  `-- 08_workflow_provenance/
      |-- <dataset>_workflow_provenance.db
      `-- <dataset>_workflow_provenance.sql

USAGE
}

STEP1_ONLY=false
CLUSTERING_ONLY=false
THROUGH_STEP4=false
FULL_FLAG=false
INPUT_FORMAT=auto
MATRIX_TYPE=auto
DISTANCE=auto
CANDIDATE_SOURCE=union
MDS_COMPONENTS=3
KMEANS_MAX_ITER=100
MISSING_VALUES=error
N_INIT=100
RANDOM_STATE=0
OUTPUT_DIR_ARG=
TAXONOMY_FILE_ARG=
POSITIONAL_ARGS=()

while (($#)); do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --full)
            FULL_FLAG=true
            shift
            ;;
        --step1-only)
            STEP1_ONLY=true
            shift
            ;;
        --clustering-only)
            CLUSTERING_ONLY=true
            shift
            ;;
        --through-step4)
            THROUGH_STEP4=true
            shift
            ;;
        --missing-values|--n-init|--random-state|--input-format|--matrix-type|--distance|--candidate-source|--mds-components|--kmeans-max-iter|--output-dir|--taxonomy)
            if (($# < 2)); then
                echo "Error: $1 requires a value." >&2
                exit 2
            fi
            case "$1" in
                --missing-values) MISSING_VALUES=$2 ;;
                --n-init) N_INIT=$2 ;;
                --random-state) RANDOM_STATE=$2 ;;
                --input-format) INPUT_FORMAT=$2 ;;
                --matrix-type) MATRIX_TYPE=$2 ;;
                --distance) DISTANCE=$2 ;;
                --candidate-source) CANDIDATE_SOURCE=$2 ;;
                --mds-components) MDS_COMPONENTS=$2 ;;
                --kmeans-max-iter) KMEANS_MAX_ITER=$2 ;;
                --output-dir) OUTPUT_DIR_ARG=$2 ;;
                --taxonomy) TAXONOMY_FILE_ARG=$2 ;;
            esac
            shift 2
            ;;
        --)
            shift
            POSITIONAL_ARGS+=("$@")
            break
            ;;
        -*)
            echo "Error: Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

if ((${#POSITIONAL_ARGS[@]} < 1 || ${#POSITIONAL_ARGS[@]} > 2)); then
    echo "Error: Provide one input matrix and, optionally, MAX_K." >&2
    usage >&2
    exit 2
fi

INPUT_FILE=${POSITIONAL_ARGS[0]}
MAX_K=${POSITIONAL_ARGS[1]:-15}

MODE_COUNT=0
for MODE_FLAG in "$STEP1_ONLY" "$CLUSTERING_ONLY" "$THROUGH_STEP4"; do
    if [ "$MODE_FLAG" = true ]; then
        MODE_COUNT=$((MODE_COUNT + 1))
    fi
done
if [ "$MODE_COUNT" -gt 1 ]; then
    echo "Error: Choose only one stop option: --step1-only, --clustering-only, or --through-step4." >&2
    exit 2
fi
if [ "$FULL_FLAG" = true ] && [ "$MODE_COUNT" -gt 0 ]; then
    echo "Error: --full cannot be combined with a stop option." >&2
    exit 2
fi

case "$INPUT_FORMAT" in auto|csv|tsv) ;; *) echo "Error: --input-format must be auto, csv, or tsv." >&2; exit 2 ;; esac
case "$MATRIX_TYPE" in auto|binary|proportion) ;; *) echo "Error: --matrix-type must be auto, binary, or proportion." >&2; exit 2 ;; esac
case "$DISTANCE" in auto|jaccard|ruzicka|euclidean) ;; *) echo "Error: Invalid --distance." >&2; exit 2 ;; esac
case "$CANDIDATE_SOURCE" in kmeans|union) ;; *) echo "Error: --candidate-source must be kmeans or union." >&2; exit 2 ;; esac
case "$MISSING_VALUES" in error|zero) ;; *) echo "Error: --missing-values must be error or zero." >&2; exit 2 ;; esac

python - "$MAX_K" "$N_INIT" "$RANDOM_STATE" "$MDS_COMPONENTS" "$KMEANS_MAX_ITER" <<'PY'
import sys

for name, raw, minimum, maximum in (
    ("MAX_K", sys.argv[1], 3, None),
    ("--n-init", sys.argv[2], 1, None),
    ("--random-state", sys.argv[3], 0, 2**32 - 1),
    ("--mds-components", sys.argv[4], 1, None),
    ("--kmeans-max-iter", sys.argv[5], 1, None),
):
    try:
        value = int(raw)
    except ValueError:
        raise SystemExit(f"Error: {name} must be an integer.")
    if value < minimum or (maximum is not None and value > maximum):
        if maximum is None:
            constraint = f"at least {minimum}"
        else:
            constraint = f"between {minimum} and {maximum}"
        raise SystemExit(f"Error: {name} must be {constraint}.")
PY

if [[ "$INPUT_FILE" = /* ]]; then
    ABS_INPUT_FILE="$INPUT_FILE"
else
    ABS_INPUT_FILE="$WORKFLOW_ROOT/$INPUT_FILE"
fi
if [ ! -f "$ABS_INPUT_FILE" ] || [ ! -r "$ABS_INPUT_FILE" ]; then
    echo "Error: Input file does not exist or is not readable: $ABS_INPUT_FILE" >&2
    exit 1
fi

ABS_TAXONOMY_FILE=
TAXONOMY_BASENAME=
if [[ -n "$TAXONOMY_FILE_ARG" ]]; then
    if [[ "$TAXONOMY_FILE_ARG" = /* ]]; then
        ABS_TAXONOMY_FILE="$TAXONOMY_FILE_ARG"
    else
        ABS_TAXONOMY_FILE="$WORKFLOW_ROOT/$TAXONOMY_FILE_ARG"
    fi
    if [ ! -f "$ABS_TAXONOMY_FILE" ] || [ ! -r "$ABS_TAXONOMY_FILE" ]; then
        echo "Error: Taxonomy file does not exist or is not readable: $ABS_TAXONOMY_FILE" >&2
        exit 1
    fi
    TAXONOMY_BASENAME=$(basename "$ABS_TAXONOMY_FILE")
fi

INPUT_BASENAME=$(basename "$ABS_INPUT_FILE")
INPUT_PREFIX="${INPUT_BASENAME%.*}"
if [[ -z "$INPUT_PREFIX" || "$INPUT_PREFIX" = . || "$INPUT_PREFIX" = .. ]]; then
    echo "Error: The input filename must have a nonempty dataset prefix." >&2
    exit 2
fi

if [[ -z "$OUTPUT_DIR_ARG" ]]; then
    RUN_NUMBER=1
    while :; do
        printf -v RUN_SUFFIX "%02d" "$RUN_NUMBER"
        CANDIDATE_OUTPUT_DIR="$WORKFLOW_ROOT/${INPUT_PREFIX}_run_${RUN_SUFFIX}"
        if [ ! -e "$CANDIDATE_OUTPUT_DIR" ] && [ ! -L "$CANDIDATE_OUTPUT_DIR" ]; then
            OUTPUT_DIR="$CANDIDATE_OUTPUT_DIR"
            break
        fi
        RUN_NUMBER=$((RUN_NUMBER + 1))
        if [ "$RUN_NUMBER" -gt 9999 ]; then
            echo "Error: Could not allocate a fresh run directory for $INPUT_PREFIX." >&2
            exit 1
        fi
    done
elif [[ "$OUTPUT_DIR_ARG" = /* ]]; then
    OUTPUT_DIR="$OUTPUT_DIR_ARG"
else
    OUTPUT_DIR="$WORKFLOW_ROOT/$OUTPUT_DIR_ARG"
fi

if [ -e "$OUTPUT_DIR" ] || [ -L "$OUTPUT_DIR" ]; then
    echo "Error: Output directory already exists: $OUTPUT_DIR" >&2
    echo "Choose another --output-dir or omit it for automatic <dataset>_run_XX naming." >&2
    exit 1
fi

# Determine which scripts are required for the selected stopping point.
REQUIRED_SCRIPTS=(
    metguild_common.py
    01_removing_zero_columns_K_Optimization.py
)
if [ "$STEP1_ONLY" = false ]; then
    REQUIRED_SCRIPTS+=(
        02_kmeans_mds_clustering.py
        03_hierarchical_clustering.py
    )
    if [[ -n "$ABS_TAXONOMY_FILE" ]]; then
        REQUIRED_SCRIPTS+=(03b_taxonomy_annotation.py)
    fi
fi
if [ "$STEP1_ONLY" = false ] && [ "$CLUSTERING_ONLY" = false ]; then
    REQUIRED_SCRIPTS+=(04_statistical_testing.py)
fi
if [ "$STEP1_ONLY" = false ] && [ "$CLUSTERING_ONLY" = false ] && [ "$THROUGH_STEP4" = false ]; then
    REQUIRED_SCRIPTS+=(
        05_feature_feature_correlation.py
        06_feature_cluster_correlation.py
        07_cluster_cluster_correlation.py
        08_merge_filtered_files.py
    )
fi

for SCRIPT_NAME in "${REQUIRED_SCRIPTS[@]}"; do
    if [ ! -f "$WORKFLOW_ROOT/scripts/$SCRIPT_NAME" ]; then
        echo "Error: Required workflow script is missing:" >&2
        echo "  $WORKFLOW_ROOT/scripts/$SCRIPT_NAME" >&2
        exit 1
    fi
done

if ! command -v python >/dev/null 2>&1; then
    echo "Error: Python was not found." >&2
    echo "Activate the MetGuild conda environment first." >&2
    exit 1
fi

python - "$STEP1_ONLY" <<'PY'
import sys

if sys.version_info < (3, 10):
    raise SystemExit("Error: Python 3.10 or newer is required.")

packages = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "sklearn": "scikit-learn",
    "matplotlib": "matplotlib",
}
if sys.argv[1] != "true":
    packages["plotly"] = "plotly"

missing = []
for module_name, package_name in packages.items():
    try:
        __import__(module_name)
    except ImportError:
        missing.append(package_name)

if missing:
    print("Error: Missing required Python packages:")
    for package in missing:
        print(f"  - {package}")
    raise SystemExit(1)
PY

# Resolve CSV/TSV format before creating the run directory.
RESOLVED_INPUT_FORMAT=$(python - "$ABS_INPUT_FILE" "$INPUT_FORMAT" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
requested = sys.argv[2]
with path.open(encoding="utf-8-sig", newline="") as handle:
    if requested == "auto":
        suffix = path.suffix.lower()
        if suffix == ".csv":
            requested = "csv"
        elif suffix == ".tsv":
            requested = "tsv"
        else:
            sample = handle.read(65536)
            handle.seek(0)
            requested = "csv" if csv.Sniffer().sniff(sample, delimiters=",\t").delimiter == "," else "tsv"
    delimiter = "," if requested == "csv" else "\t"
    header = next(csv.reader(handle, delimiter=delimiter, strict=True), [])
if "Genomes" not in header:
    raise SystemExit("Error: Input header must contain a Genomes column; check the delimiter/input format.")
print(requested)
PY
)

mkdir -p "$OUTPUT_DIR/scripts"
cp "$ABS_INPUT_FILE" "$OUTPUT_DIR/$INPUT_BASENAME"
RUN_TAXONOMY_FILE=
if [[ -n "$ABS_TAXONOMY_FILE" ]]; then
    mkdir -p "$OUTPUT_DIR/inputs"
    cp "$ABS_TAXONOMY_FILE" "$OUTPUT_DIR/inputs/$TAXONOMY_BASENAME"
    RUN_TAXONOMY_FILE="$OUTPUT_DIR/inputs/$TAXONOMY_BASENAME"
fi
for SCRIPT_NAME in "${REQUIRED_SCRIPTS[@]}"; do
    cp "$WORKFLOW_ROOT/scripts/$SCRIPT_NAME" "$OUTPUT_DIR/scripts/$SCRIPT_NAME"
done
cp "${BASH_SOURCE[0]}" "$OUTPUT_DIR/scripts/run_metguild_workflow.sh"

RUN_INPUT_FILE="$OUTPUT_DIR/$INPUT_BASENAME"

cd "$OUTPUT_DIR"

RUN_STARTED_AT=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

echo "============================================================"
echo "Starting MetGuild"
echo "============================================================"
echo "Input:          $ABS_INPUT_FILE"
echo "Run directory:  $OUTPUT_DIR"
echo "Matrix type:    $MATRIX_TYPE"
echo "Input format:   $RESOLVED_INPUT_FORMAT"
echo "MAX_K:          $MAX_K (exclusive)"
echo "Candidate source: $CANDIDATE_SOURCE"
if [[ -n "$RUN_TAXONOMY_FILE" ]]; then
    echo "Taxonomy:       $ABS_TAXONOMY_FILE"
else
    echo "Taxonomy:       not supplied"
fi
echo ""

run_step() {
    local number=$1
    local label=$2
    shift 2
    echo "============================================================"
    echo "Step ${number}/9: ${label}"
    echo "============================================================"
    "$@"
    echo "Step ${number} complete."
    echo ""
}

run_step 1 "Matrix preparation and candidate-k selection" \
    python scripts/01_removing_zero_columns_K_Optimization.py \
        --input "$RUN_INPUT_FILE" \
        --cleaned-output "$OUTPUT_DIR/${INPUT_PREFIX}_cleaned.tsv" \
        --plot-dir "$OUTPUT_DIR/02_k_optimization" \
        --manifest "$OUTPUT_DIR/manifests/optimal_k_values.json" \
        --max-k "$MAX_K" \
        --missing-values "$MISSING_VALUES" \
        --input-format "$RESOLVED_INPUT_FORMAT" \
        --matrix-type "$MATRIX_TYPE" \
        --distance "$DISTANCE" \
        --candidate-source "$CANDIDATE_SOURCE" \
        --mds-components "$MDS_COMPONENTS" \
        --kmeans-max-iter "$KMEANS_MAX_ITER" \
        --n-init "$N_INIT" \
        --random-state "$RANDOM_STATE"

if [ ! -f "$OUTPUT_DIR/manifests/optimal_k_values.json" ]; then
    echo "Error: Step 1 did not create the expected manifest." >&2
    exit 1
fi

# Full downstream workflow supports only binary and proportion matrices.
RESOLVED_MATRIX_TYPE=$(python - "$OUTPUT_DIR/manifests/optimal_k_values.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    manifest = json.load(handle)
print(manifest.get("preprocessing", {}).get("matrix_type", ""))
PY
)

if [ "$STEP1_ONLY" = true ]; then
    echo "Step 1-only run completed successfully."
    echo "Run directory: $OUTPUT_DIR"
    exit 0
fi

if [[ "$RESOLVED_MATRIX_TYPE" != "binary" && "$RESOLVED_MATRIX_TYPE" != "proportion" ]]; then
    echo "Error: The downstream MetGuild workflow supports matrix types binary and proportion." >&2
    echo "Step 1 resolved this dataset as: $RESOLVED_MATRIX_TYPE" >&2
    echo "For proportional data, rerun with --matrix-type proportion." >&2
    exit 1
fi

run_step 2 "KMeans clustering on classical-MDS coordinates" \
    python scripts/02_kmeans_mds_clustering.py --run-dir "$OUTPUT_DIR"

run_step 3 "Average-linkage hierarchical clustering" \
    python scripts/03_hierarchical_clustering.py --run-dir "$OUTPUT_DIR"

if [[ -n "$RUN_TAXONOMY_FILE" ]]; then
    echo "============================================================"
    echo "Step 3b/9: Taxonomy annotation"
    echo "============================================================"
    python scripts/03b_taxonomy_annotation.py \
        --run-dir "$OUTPUT_DIR" \
        --taxonomy "$RUN_TAXONOMY_FILE"
    echo "Step 3b complete."
    echo ""
fi

if [ "$CLUSTERING_ONLY" = true ]; then
    echo "Steps 1-3 completed successfully."
    echo "Run directory: $OUTPUT_DIR"
    exit 0
fi

run_step 4 "Statistical feature screening" \
    python scripts/04_statistical_testing.py --run-dir "$OUTPUT_DIR"

if [ "$THROUGH_STEP4" = true ]; then
    echo "Steps 1-4 completed successfully."
    echo "Run directory: $OUTPUT_DIR"
    exit 0
fi

run_step 5 "Feature-feature Spearman correlation" \
    python scripts/05_feature_feature_correlation.py --run-dir "$OUTPUT_DIR"

run_step 6 "Feature-cluster Spearman correlation" \
    python scripts/06_feature_cluster_correlation.py --run-dir "$OUTPUT_DIR"

run_step 7 "Cluster-cluster Spearman correlation" \
    python scripts/07_cluster_cluster_correlation.py --run-dir "$OUTPUT_DIR"

STEP8_ARGS=(--run-dir "$OUTPUT_DIR")
if [[ -n "$RUN_TAXONOMY_FILE" ]]; then
    STEP8_ARGS+=(--require-taxonomy)
fi
run_step 8 "Network edge/node construction" \
    python scripts/08_merge_filtered_files.py "${STEP8_ARGS[@]}"


# ------------------------------------------------------------
# Step 9: Build populated SQLite/SQL workflow provenance
# ------------------------------------------------------------
echo "============================================================"
echo "Step 9/9: Building populated SQLite workflow provenance"
echo "============================================================"

RUN_COMPLETED_AT=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
PROVENANCE_DIR="$OUTPUT_DIR/08_workflow_provenance"
PROVENANCE_DB="$PROVENANCE_DIR/${INPUT_PREFIX}_workflow_provenance.db"
PROVENANCE_SQL="$PROVENANCE_DIR/${INPUT_PREFIX}_workflow_provenance.sql"
mkdir -p "$PROVENANCE_DIR"
rm -f "$PROVENANCE_DB" "$PROVENANCE_SQL"

if [[ "$RESOLVED_MATRIX_TYPE" == "binary" ]]; then
    STATISTICAL_TEST="chi-square"
else
    STATISTICAL_TEST="ANOVA"
fi

PROVENANCE_ARGS=(
    --run-dir "$OUTPUT_DIR"
    --original-input "$ABS_INPUT_FILE"
    --working-input "$INPUT_BASENAME"
    --max-k "$MAX_K"
    --matrix-type "$RESOLVED_MATRIX_TYPE"
    --statistical-test "$STATISTICAL_TEST"
    --candidate-source "$CANDIDATE_SOURCE"
    --distance "$DISTANCE"
    --started-at "$RUN_STARTED_AT"
    --completed-at "$RUN_COMPLETED_AT"
    --db "$PROVENANCE_DB"
    --sql "$PROVENANCE_SQL"
)
if [[ -n "$ABS_TAXONOMY_FILE" ]]; then
    PROVENANCE_ARGS+=(
        --original-taxonomy "$ABS_TAXONOMY_FILE"
        --working-taxonomy "inputs/$TAXONOMY_BASENAME"
    )
fi

python - "${PROVENANCE_ARGS[@]}" <<'METGUILD_PROVENANCE_PY'
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import re
import sqlite3
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence


@dataclass(frozen=True)
class Step:
    step_id: int
    order_no: int
    label: str
    name: str
    script_name: str | None
    description: str
    optional: int = 0


@dataclass(frozen=True)
class Artifact:
    artifact_id: int
    relative_path: str
    artifact_name: str
    artifact_role: str
    artifact_type: str
    extension: str | None
    size_bytes: int | None
    sha256: str | None
    modified_at: str | None
    producing_step_id: int | None
    is_external: int


STEPS = [
    Step(1, 0, "setup", "Run setup and input staging", "run_metguild_workflow.sh",
         "Validate the environment, create the run directory, and stage inputs and scripts."),
    Step(2, 1, "1", "Matrix preparation and candidate-k selection", "01_removing_zero_columns_K_Optimization.py",
         "Clean the feature matrix, calculate dissimilarities/MDS, and determine candidate k values."),
    Step(3, 2, "2", "KMeans clustering", "02_kmeans_mds_clustering.py",
         "Cluster the classical-MDS coordinates using KMeans for each forwarded k."),
    Step(4, 3, "3", "Hierarchical clustering", "03_hierarchical_clustering.py",
         "Apply average-linkage hierarchical clustering to the original dissimilarities."),
    Step(5, 4, "3b", "Taxonomy annotation", "03b_taxonomy_annotation.py",
         "Join clustered genomes to taxonomy and generate Cluster/MAG_Phyla annotations.", 1),
    Step(6, 5, "4", "Statistical feature screening", "04_statistical_testing.py",
         "Use chi-square for binary matrices or ANOVA for proportional matrices with BH correction."),
    Step(7, 6, "5", "Feature-feature correlation", "05_feature_feature_correlation.py",
         "Calculate and filter feature-feature Spearman correlations."),
    Step(8, 7, "6", "Feature-cluster correlation", "06_feature_cluster_correlation.py",
         "Calculate and filter one-vs-rest feature-cluster Spearman associations."),
    Step(9, 8, "7", "Cluster-cluster correlation", "07_cluster_cluster_correlation.py",
         "Calculate and filter Spearman correlations between cluster feature profiles."),
    Step(10, 9, "8", "Network construction", "08_merge_filtered_files.py",
         "Merge filtered correlation edges and create node attributes, including taxonomy when supplied."),
    Step(11, 10, "9", "Workflow provenance export", "run_metguild_workflow.sh",
         "Create the populated SQLite provenance database and complete SQL reconstruction dump."),
]

STEP_BY_ORDER = {step.order_no: step for step in STEPS}
STEP_BY_ID = {step.step_id: step for step in STEPS}

STEP_DEPENDENCIES = [
    (0, 1, "The staged input matrix is cleaned and evaluated for candidate k values."),
    (1, 2, "The Step 1 matrix/MDS outputs and candidate-k manifest drive KMeans clustering."),
    (1, 3, "The Step 1 dissimilarity/linkage outputs and candidate-k manifest drive hierarchical clustering."),
    (2, 4, "KMeans assignments are annotated with taxonomy when a taxonomy file is supplied."),
    (3, 4, "Hierarchical assignments are annotated with taxonomy when a taxonomy file is supplied."),
    (1, 5, "The cleaned functional matrix supplies feature values for statistical screening."),
    (2, 5, "KMeans assignments define groups for statistical screening."),
    (3, 5, "Hierarchical assignments define groups for statistical screening."),
    (5, 6, "Statistically retained features are used for feature-feature correlation."),
    (5, 7, "Statistically retained features and cluster labels are used for feature-cluster correlation."),
    (5, 8, "Statistically retained features and cluster labels are used to build cluster profiles."),
    (6, 9, "Filtered feature-feature correlations contribute network edges."),
    (7, 9, "Filtered feature-cluster correlations contribute network edges."),
    (8, 9, "Filtered cluster-cluster correlations contribute network edges."),
    (4, 9, "Cluster taxonomy annotations contribute MAG_Phyla node attributes when available."),
]
for upstream in range(0, 10):
    STEP_DEPENDENCIES.append((upstream, 10, "Completed workflow artifacts are indexed by the provenance export."))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build MetGuild workflow provenance SQLite/SQL outputs.")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--original-input", required=True)
    parser.add_argument("--working-input", required=True)
    parser.add_argument("--original-taxonomy")
    parser.add_argument("--working-taxonomy")
    parser.add_argument("--max-k", required=True, type=int)
    parser.add_argument("--matrix-type", required=True, choices=["binary", "proportion"])
    parser.add_argument("--statistical-test", required=True, choices=["chi-square", "ANOVA"])
    parser.add_argument("--candidate-source", required=True)
    parser.add_argument("--distance", required=True)
    parser.add_argument("--started-at", required=True)
    parser.add_argument("--completed-at", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--sql", required=True)
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def iso_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def slugify(text: str, max_length: int = 50) -> str:
    slug = re.sub(r"[^0-9A-Za-z_]+", "_", str(text)).strip("_").lower()
    slug = re.sub(r"_+", "_", slug) or "result"
    if slug[0].isdigit():
        slug = "t_" + slug
    return slug[:max_length].rstrip("_") or "result"


def classify_type(path: Path) -> str:
    suffix = path.suffix.lower()
    parts = set(path.parts)
    if suffix in {".py", ".sh"}:
        return "script"
    if suffix == ".sql":
        return "database_sql"
    if suffix == ".db":
        return "sqlite_database"
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


def producing_step_order(relative_path: str, working_input: str, working_taxonomy: str | None) -> int | None:
    path = Path(relative_path)
    text = path.as_posix()
    if text == working_input or text.startswith("scripts/"):
        return 0
    if working_taxonomy and text == working_taxonomy:
        return 0
    if text.startswith("inputs/"):
        return 0
    if text.startswith("02_k_optimization/") or text.startswith("manifests/"):
        return 1
    if text.startswith("03_04_clustering/01_kmeans/"):
        return 2
    if text.startswith("03_04_clustering/02_hierarchical/"):
        return 3
    if text.startswith("03_04_clustering/taxonomy/"):
        return 4
    if text.startswith("05_stat_analysis/"):
        return 5
    if text.startswith("06_correlation/01_featurevsfeature_correlation/"):
        return 6
    if text.startswith("06_correlation/02_featurevscluster_correlation/"):
        return 7
    if text.startswith("06_correlation/03_clustervscluster_correlation/"):
        return 8
    if text.startswith("07_network_construction/"):
        return 9
    if text.startswith("08_workflow_provenance/"):
        return 10
    if len(path.parts) == 1 and path.suffix.lower() in {".csv", ".tsv", ".json", ".txt"}:
        return 1
    return None


def step_id_for_order(order: int | None) -> int | None:
    return STEP_BY_ORDER[order].step_id if order in STEP_BY_ORDER else None


def scan_artifacts(args: argparse.Namespace, run_dir: Path, db_path: Path, sql_path: Path) -> list[Artifact]:
    artifacts: list[Artifact] = []
    next_id = 1

    def add_external(path_text: str, role: str) -> None:
        nonlocal next_id
        path = Path(path_text).resolve()
        if not path.is_file():
            return
        artifacts.append(Artifact(
            next_id, str(path), path.name, role, classify_type(path), path.suffix.lower() or None,
            path.stat().st_size, sha256sum(path), iso_mtime(path), None, 1,
        ))
        next_id += 1

    add_external(args.original_input, "external_input")
    if args.original_taxonomy:
        add_external(args.original_taxonomy, "external_taxonomy")

    excluded = {db_path.resolve(), sql_path.resolve()}
    for path in sorted((p for p in run_dir.rglob("*") if p.is_file()), key=lambda p: p.relative_to(run_dir).as_posix()):
        if path.resolve() in excluded:
            continue
        rel = path.relative_to(run_dir).as_posix()
        order = producing_step_order(rel, args.working_input, args.working_taxonomy)
        role = "script" if rel.startswith("scripts/") else "workflow_output"
        if rel == args.working_input:
            role = "staged_input"
        elif args.working_taxonomy and rel == args.working_taxonomy:
            role = "staged_taxonomy"
        artifacts.append(Artifact(
            next_id, rel, path.name, role, classify_type(path.relative_to(run_dir)), path.suffix.lower() or None,
            path.stat().st_size, sha256sum(path), iso_mtime(path), step_id_for_order(order), 0,
        ))
        next_id += 1

    # Self-referential provenance outputs are represented without checksum/size.
    for path, kind in ((db_path, "sqlite_database"), (sql_path, "database_sql")):
        rel = path.relative_to(run_dir).as_posix()
        artifacts.append(Artifact(
            next_id, rel, path.name, "provenance_output", kind, path.suffix.lower(),
            None, None, None, STEP_BY_ORDER[10].step_id, 0,
        ))
        next_id += 1
    return artifacts


def create_provenance_schema(con: sqlite3.Connection) -> None:
    con.executescript("""
    PRAGMA foreign_keys = ON;
    CREATE TABLE workflow_runs (
        run_id INTEGER PRIMARY KEY,
        run_name TEXT NOT NULL,
        dataset_name TEXT NOT NULL,
        original_input_path TEXT NOT NULL,
        original_taxonomy_path TEXT,
        run_directory TEXT NOT NULL,
        max_k INTEGER NOT NULL,
        matrix_type TEXT NOT NULL,
        statistical_test TEXT NOT NULL,
        candidate_source TEXT NOT NULL,
        requested_distance TEXT NOT NULL,
        started_at TEXT NOT NULL,
        completed_at TEXT NOT NULL,
        status TEXT NOT NULL
    );
    CREATE TABLE workflow_steps (
        step_id INTEGER PRIMARY KEY,
        workflow_order INTEGER NOT NULL UNIQUE,
        step_label TEXT NOT NULL,
        step_name TEXT NOT NULL,
        script_name TEXT,
        description TEXT NOT NULL,
        optional INTEGER NOT NULL DEFAULT 0
    );
    CREATE TABLE step_dependencies (
        upstream_step_id INTEGER NOT NULL,
        downstream_step_id INTEGER NOT NULL,
        dependency_description TEXT NOT NULL,
        PRIMARY KEY (upstream_step_id, downstream_step_id),
        FOREIGN KEY (upstream_step_id) REFERENCES workflow_steps(step_id),
        FOREIGN KEY (downstream_step_id) REFERENCES workflow_steps(step_id)
    );
    CREATE TABLE artifacts (
        artifact_id INTEGER PRIMARY KEY,
        run_id INTEGER NOT NULL,
        relative_path TEXT NOT NULL,
        artifact_name TEXT NOT NULL,
        artifact_role TEXT NOT NULL,
        artifact_type TEXT NOT NULL,
        extension TEXT,
        size_bytes INTEGER,
        sha256 TEXT,
        modified_at TEXT,
        producing_step_id INTEGER,
        is_external INTEGER NOT NULL DEFAULT 0,
        UNIQUE (run_id, relative_path),
        FOREIGN KEY (run_id) REFERENCES workflow_runs(run_id),
        FOREIGN KEY (producing_step_id) REFERENCES workflow_steps(step_id)
    );
    CREATE TABLE step_artifacts (
        step_id INTEGER NOT NULL,
        artifact_id INTEGER NOT NULL,
        io_role TEXT NOT NULL CHECK (io_role IN ('input','output','script')),
        PRIMARY KEY (step_id, artifact_id, io_role),
        FOREIGN KEY (step_id) REFERENCES workflow_steps(step_id),
        FOREIGN KEY (artifact_id) REFERENCES artifacts(artifact_id)
    );
    CREATE TABLE artifact_dependencies (
        parent_artifact_id INTEGER NOT NULL,
        child_artifact_id INTEGER NOT NULL,
        relation TEXT NOT NULL,
        PRIMARY KEY (parent_artifact_id, child_artifact_id, relation),
        FOREIGN KEY (parent_artifact_id) REFERENCES artifacts(artifact_id),
        FOREIGN KEY (child_artifact_id) REFERENCES artifacts(artifact_id)
    );
    CREATE TABLE dataset_tables (
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
    CREATE TABLE dataset_columns (
        dataset_table_id INTEGER NOT NULL,
        ordinal_position INTEGER NOT NULL,
        source_column_name TEXT NOT NULL,
        sqlite_column_name TEXT NOT NULL,
        sqlite_type TEXT NOT NULL,
        PRIMARY KEY (dataset_table_id, ordinal_position),
        UNIQUE (dataset_table_id, sqlite_column_name),
        FOREIGN KEY (dataset_table_id) REFERENCES dataset_tables(dataset_table_id) ON DELETE CASCADE
    );
    CREATE TABLE table_dependencies (
        parent_dataset_table_id INTEGER NOT NULL,
        child_dataset_table_id INTEGER NOT NULL,
        relation TEXT NOT NULL,
        PRIMARY KEY (parent_dataset_table_id, child_dataset_table_id, relation),
        FOREIGN KEY (parent_dataset_table_id) REFERENCES dataset_tables(dataset_table_id),
        FOREIGN KEY (child_dataset_table_id) REFERENCES dataset_tables(dataset_table_id)
    );
    CREATE TABLE result_import_runs (
        import_run_id INTEGER PRIMARY KEY,
        imported_at TEXT NOT NULL,
        run_directory TEXT NOT NULL,
        imported_tables INTEGER NOT NULL,
        failed_tables INTEGER NOT NULL,
        imported_rows INTEGER NOT NULL,
        table_dependencies INTEGER NOT NULL
    );
    """)


def populate_provenance(con: sqlite3.Connection, args: argparse.Namespace, artifacts: list[Artifact], run_dir: Path) -> None:
    run_name = f"{Path(args.working_input).stem}_{args.started_at.replace(':','').replace('-','')}"
    con.execute(
        "INSERT INTO workflow_runs VALUES (1,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            run_name, Path(args.working_input).stem, str(Path(args.original_input).resolve()),
            str(Path(args.original_taxonomy).resolve()) if args.original_taxonomy else None,
            str(run_dir), args.max_k, args.matrix_type, args.statistical_test,
            args.candidate_source, args.distance, args.started_at, args.completed_at, "completed",
        ),
    )
    con.executemany(
        "INSERT INTO workflow_steps VALUES (?,?,?,?,?,?,?)",
        [(s.step_id, s.order_no, s.label, s.name, s.script_name, s.description, s.optional) for s in STEPS],
    )
    con.executemany(
        "INSERT INTO step_dependencies VALUES (?,?,?)",
        [(STEP_BY_ORDER[u].step_id, STEP_BY_ORDER[d].step_id, desc) for u, d, desc in STEP_DEPENDENCIES],
    )
    con.executemany(
        "INSERT INTO artifacts VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        [(
            a.artifact_id, 1, a.relative_path, a.artifact_name, a.artifact_role, a.artifact_type,
            a.extension, a.size_bytes, a.sha256, a.modified_at, a.producing_step_id, a.is_external,
        ) for a in artifacts],
    )


def artifact_by_rel(artifacts: list[Artifact]) -> dict[str, Artifact]:
    return {a.relative_path: a for a in artifacts}


def data_artifacts_for_step(artifacts: list[Artifact], order: int) -> list[Artifact]:
    step_id = STEP_BY_ORDER[order].step_id
    return [a for a in artifacts if a.producing_step_id == step_id and a.artifact_role != "script"]


def import_candidate(a: Artifact) -> bool:
    return (a.extension or "").lower() in {".csv", ".tsv", ".json"} and a.artifact_role != "script"


def selected_step_inputs(artifacts: list[Artifact], args: argparse.Namespace, order: int) -> list[Artifact]:
    rel = artifact_by_rel(artifacts)
    main_staged = [rel[args.working_input]] if args.working_input in rel else []
    tax_staged = [rel[args.working_taxonomy]] if args.working_taxonomy and args.working_taxonomy in rel else []
    s1 = data_artifacts_for_step(artifacts, 1)
    s2 = data_artifacts_for_step(artifacts, 2)
    s3 = data_artifacts_for_step(artifacts, 3)
    s4 = data_artifacts_for_step(artifacts, 4)
    s5 = data_artifacts_for_step(artifacts, 5)
    s6 = data_artifacts_for_step(artifacts, 6)
    s7 = data_artifacts_for_step(artifacts, 7)
    s8 = data_artifacts_for_step(artifacts, 8)
    selections = {
        0: [a for a in artifacts if a.is_external],
        1: main_staged,
        2: [a for a in s1 if import_candidate(a)],
        3: [a for a in s1 if import_candidate(a)],
        4: tax_staged + [a for a in (s2 + s3) if a.extension in {".csv", ".tsv"}],
        5: [a for a in (s1 + s2 + s3) if a.extension in {".csv", ".tsv", ".json"}],
        6: [a for a in s5 if a.extension in {".csv", ".json"}],
        7: [a for a in s5 if a.extension in {".csv", ".json"}],
        8: [a for a in s5 if a.extension in {".csv", ".json"}],
        9: [a for a in (s6 + s7 + s8 + s4) if a.extension in {".csv", ".json"}],
        10: [a for a in artifacts if a.producing_step_id != STEP_BY_ORDER[10].step_id],
    }
    out: list[Artifact] = []
    seen: set[int] = set()
    for a in selections.get(order, []):
        if a.artifact_id not in seen:
            seen.add(a.artifact_id)
            out.append(a)
    return out


def link_artifacts(con: sqlite3.Connection, artifacts: list[Artifact], args: argparse.Namespace) -> None:
    rel = artifact_by_rel(artifacts)
    for step in STEPS:
        if step.script_name:
            script = rel.get("scripts/" + step.script_name)
            if script:
                con.execute("INSERT OR IGNORE INTO step_artifacts VALUES (?,?, 'script')", (step.step_id, script.artifact_id))
        for a in selected_step_inputs(artifacts, args, step.order_no):
            con.execute("INSERT OR IGNORE INTO step_artifacts VALUES (?,?, 'input')", (step.step_id, a.artifact_id))
        for a in data_artifacts_for_step(artifacts, step.order_no):
            con.execute("INSERT OR IGNORE INTO step_artifacts VALUES (?,?, 'output')", (step.step_id, a.artifact_id))

    # External -> staged copy lineage.
    external_main = next((a for a in artifacts if a.artifact_role == "external_input"), None)
    staged_main = rel.get(args.working_input)
    if external_main and staged_main:
        con.execute("INSERT OR IGNORE INTO artifact_dependencies VALUES (?,?,?)",
                    (external_main.artifact_id, staged_main.artifact_id, "copied_as"))
    if args.original_taxonomy and args.working_taxonomy:
        external_tax = next((a for a in artifacts if a.artifact_role == "external_taxonomy"), None)
        staged_tax = rel.get(args.working_taxonomy)
        if external_tax and staged_tax:
            con.execute("INSERT OR IGNORE INTO artifact_dependencies VALUES (?,?,?)",
                        (external_tax.artifact_id, staged_tax.artifact_id, "copied_as"))

    # Data/script lineage from selected inputs to each output.
    for step in STEPS:
        inputs = selected_step_inputs(artifacts, args, step.order_no)
        script = rel.get("scripts/" + step.script_name) if step.script_name else None
        outputs = data_artifacts_for_step(artifacts, step.order_no)
        for child in outputs:
            for parent in inputs:
                if parent.artifact_id != child.artifact_id:
                    con.execute("INSERT OR IGNORE INTO artifact_dependencies VALUES (?,?,?)",
                                (parent.artifact_id, child.artifact_id, "used_to_generate"))
            if script and script.artifact_id != child.artifact_id:
                con.execute("INSERT OR IGNORE INTO artifact_dependencies VALUES (?,?,?)",
                            (script.artifact_id, child.artifact_id, "generated_by_script"))


def detect_delimiter(path: Path, extension: str | None) -> str:
    sample = path.read_text(encoding="utf-8-sig", errors="replace")[:65536]
    if not sample:
        return ","
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;|").delimiter
    except csv.Error:
        first = next((line for line in sample.splitlines() if line.strip()), "")
        if "\t" in first:
            return "\t"
        if "," in first:
            return ","
        if ";" in first:
            return ";"
        if len(first.split()) > 1:
            return "__whitespace__"
        return "\t" if (extension or path.suffix).lower() == ".tsv" else ","


def open_delimited_rows(path: Path, extension: str | None) -> tuple[list[str], Iterator[list[str]]]:
    delimiter = detect_delimiter(path, extension)
    handle = path.open("r", encoding="utf-8-sig", errors="replace", newline="")
    if delimiter == "__whitespace__":
        def gen() -> Iterator[list[str]]:
            try:
                for line in handle:
                    if line.strip():
                        yield line.strip().split()
            finally:
                handle.close()
        iterator = gen()
    else:
        reader = csv.reader(handle, delimiter=delimiter)
        def gen() -> Iterator[list[str]]:
            try:
                yield from reader
            finally:
                handle.close()
        iterator = gen()
    try:
        header = next(iterator)
    except StopIteration:
        header = []
    return [str(x) for x in header], iterator


def value_kind(value: str) -> str | None:
    text = str(value).strip()
    if text == "":
        return None
    if re.fullmatch(r"[+-]?\d+", text):
        return "INTEGER"
    try:
        float(text)
        return "REAL"
    except ValueError:
        return "TEXT"


def merge_type(current: str | None, observed: str | None) -> str | None:
    if observed is None:
        return current
    if current is None or current == observed:
        return observed if current is None else current
    if {current, observed} <= {"INTEGER", "REAL"}:
        return "REAL"
    return "TEXT"


def infer_types(rows: Sequence[Sequence[str]], ncols: int) -> list[str]:
    kinds: list[str | None] = [None] * ncols
    for row in rows:
        padded = list(row[:ncols]) + [""] * max(0, ncols - len(row))
        for i, value in enumerate(padded):
            kinds[i] = merge_type(kinds[i], value_kind(value))
    return [k or "TEXT" for k in kinds]


def convert_value(value: str, typ: str) -> Any:
    text = str(value).strip()
    if text == "":
        return None
    if typ == "INTEGER":
        try:
            return int(text)
        except ValueError:
            return text
    if typ == "REAL":
        try:
            return float(text)
        except ValueError:
            return text
    return text


def unique_column_names(source: Sequence[str]) -> list[str]:
    reserved = {"__row_id", "__dataset_table_id", "__source_row_number"}
    used = set(reserved)
    result: list[str] = []
    for i, name in enumerate(source, 1):
        base = slugify(name or f"column_{i}", 54)
        if base in reserved:
            base = "data_" + base.lstrip("_")
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used.add(candidate)
        result.append(candidate)
    return result


def make_result_table_name(a: Artifact) -> str:
    digest = hashlib.sha1(a.relative_path.encode("utf-8")).hexdigest()[:8]
    step = a.producing_step_id or 0
    return f"result_s{step:02d}_{slugify(Path(a.artifact_name).stem, 42)}_{digest}"


def create_result_table(con: sqlite3.Connection, name: str, dataset_id: int,
                        columns: Sequence[str], types: Sequence[str]) -> None:
    defs = [
        '"__row_id" INTEGER PRIMARY KEY AUTOINCREMENT',
        f'"__dataset_table_id" INTEGER NOT NULL DEFAULT {dataset_id} REFERENCES dataset_tables(dataset_table_id)',
        '"__source_row_number" INTEGER NOT NULL',
    ]
    defs.extend(f"{quote_identifier(c)} {t}" for c, t in zip(columns, types))
    con.execute(f"CREATE TABLE {quote_identifier(name)} (" + ", ".join(defs) + ")")


def import_delimited(con: sqlite3.Connection, path: Path, a: Artifact,
                     table_name: str, dataset_id: int) -> tuple[int, int, list[tuple[str,str,str]]]:
    header, rows_iter = open_delimited_rows(path, a.extension)
    if not header:
        create_result_table(con, table_name, dataset_id, [], [])
        return 0, 0, []
    sample = list(itertools.islice(rows_iter, 5000))
    types = infer_types(sample, len(header))
    columns = unique_column_names(header)
    create_result_table(con, table_name, dataset_id, columns, types)
    insert_cols = ["__dataset_table_id", "__source_row_number", *columns]
    sql = (f"INSERT INTO {quote_identifier(table_name)} (" +
           ",".join(quote_identifier(c) for c in insert_cols) + ") VALUES (" +
           ",".join("?" for _ in insert_cols) + ")")
    row_count = 0
    for row_no, row in enumerate(itertools.chain(sample, rows_iter), start=1):
        vals = list(row[:len(header)]) + [""] * max(0, len(header) - len(row))
        converted = [convert_value(v, t) for v, t in zip(vals, types)]
        con.execute(sql, (dataset_id, row_no, *converted))
        row_count += 1
    meta = [(src, col, typ) for src, col, typ in zip(header, columns, types)]
    return row_count, len(header), meta


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
        for i, child in enumerate(value):
            yield from flatten_json(child, f"{path}/{i}")
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


def import_json(con: sqlite3.Connection, path: Path, table_name: str, dataset_id: int) -> tuple[int,int,list[tuple[str,str,str]]]:
    document = json.loads(path.read_text(encoding="utf-8-sig", errors="replace"))
    cols = ["json_path", "value", "value_type"]
    types = ["TEXT", "TEXT", "TEXT"]
    create_result_table(con, table_name, dataset_id, cols, types)
    rows = list(flatten_json(document))
    con.executemany(
        f'INSERT INTO {quote_identifier(table_name)} ("__dataset_table_id","__source_row_number","json_path","value","value_type") VALUES (?,?,?,?,?)',
        [(dataset_id, i, *row) for i, row in enumerate(rows, 1)],
    )
    return len(rows), 3, list(zip(cols, cols, types))


def resolve_artifact_path(run_dir: Path, a: Artifact) -> Path:
    p = Path(a.relative_path)
    return p if p.is_absolute() else run_dir / p


def import_scientific_tables(con: sqlite3.Connection, run_dir: Path, artifacts: list[Artifact]) -> tuple[int,int,int]:
    imported = failed = total_rows = 0
    for a in artifacts:
        if (a.extension or "").lower() not in {".csv", ".tsv", ".json"} or a.artifact_role == "script":
            continue
        path = resolve_artifact_path(run_dir, a)
        if not path.is_file():
            continue
        table_name = make_result_table_name(a)
        cur = con.execute(
            "INSERT INTO dataset_tables (artifact_id,sqlite_table_name,display_name,source_path,source_format,producing_step_id,row_count,column_count,import_status,import_error,imported_at) VALUES (?,?,?,?,?,?,0,0,'pending',NULL,?)",
            (a.artifact_id, table_name, a.artifact_name, a.relative_path, (a.extension or "").lstrip("."), a.producing_step_id, utc_now()),
        )
        dataset_id = int(cur.lastrowid)
        try:
            if a.extension == ".json":
                nrows, ncols, col_meta = import_json(con, path, table_name, dataset_id)
            else:
                nrows, ncols, col_meta = import_delimited(con, path, a, table_name, dataset_id)
            con.executemany(
                "INSERT INTO dataset_columns VALUES (?,?,?,?,?)",
                [(dataset_id, i, src, col, typ) for i, (src, col, typ) in enumerate(col_meta, 1)],
            )
            con.execute(
                "UPDATE dataset_tables SET row_count=?, column_count=?, import_status='imported', imported_at=? WHERE dataset_table_id=?",
                (nrows, ncols, utc_now(), dataset_id),
            )
            imported += 1
            total_rows += nrows
            print(f"IMPORTED {a.relative_path} -> {table_name} ({nrows} rows, {ncols} columns)")
        except Exception as exc:
            failed += 1
            con.execute(
                "UPDATE dataset_tables SET import_status='failed', import_error=?, imported_at=? WHERE dataset_table_id=?",
                (str(exc), utc_now(), dataset_id),
            )
            print(f"FAILED   {a.relative_path}: {exc}", file=sys.stderr)
    return imported, failed, total_rows


def rebuild_table_dependencies(con: sqlite3.Connection) -> int:
    con.execute("DELETE FROM table_dependencies")
    con.execute("""
        INSERT OR IGNORE INTO table_dependencies(parent_dataset_table_id,child_dataset_table_id,relation)
        SELECT p.dataset_table_id, c.dataset_table_id, d.relation
        FROM artifact_dependencies d
        JOIN dataset_tables p ON p.artifact_id=d.parent_artifact_id AND p.import_status='imported'
        JOIN dataset_tables c ON c.artifact_id=d.child_artifact_id AND c.import_status='imported'
        WHERE p.dataset_table_id != c.dataset_table_id
    """)
    return int(con.execute("SELECT COUNT(*) FROM table_dependencies").fetchone()[0])


def create_views_and_indexes(con: sqlite3.Connection) -> None:
    con.executescript("""
    CREATE VIEW artifact_connections AS
    SELECT p.artifact_id AS parent_artifact_id, p.relative_path AS parent_path,
           ps.step_label AS parent_step, d.relation,
           c.artifact_id AS child_artifact_id, c.relative_path AS child_path,
           cs.step_label AS child_step
    FROM artifact_dependencies d
    JOIN artifacts p ON p.artifact_id=d.parent_artifact_id
    JOIN artifacts c ON c.artifact_id=d.child_artifact_id
    LEFT JOIN workflow_steps ps ON ps.step_id=p.producing_step_id
    LEFT JOIN workflow_steps cs ON cs.step_id=c.producing_step_id;

    CREATE VIEW step_file_inventory AS
    SELECT s.workflow_order, s.step_label, s.step_name, sa.io_role,
           a.relative_path, a.artifact_type, a.size_bytes, a.sha256
    FROM step_artifacts sa
    JOIN workflow_steps s ON s.step_id=sa.step_id
    JOIN artifacts a ON a.artifact_id=sa.artifact_id;

    CREATE VIEW scientific_table_inventory AS
    SELECT dt.dataset_table_id, dt.sqlite_table_name, dt.display_name, dt.source_path,
           dt.source_format, ws.step_label, ws.step_name, dt.row_count,
           dt.column_count, dt.import_status, dt.import_error, dt.imported_at
    FROM dataset_tables dt
    LEFT JOIN workflow_steps ws ON ws.step_id=dt.producing_step_id;

    CREATE VIEW scientific_table_connections AS
    SELECT p.dataset_table_id AS parent_dataset_table_id,
           p.sqlite_table_name AS parent_sqlite_table, p.source_path AS parent_source_path,
           td.relation, c.dataset_table_id AS child_dataset_table_id,
           c.sqlite_table_name AS child_sqlite_table, c.source_path AS child_source_path
    FROM table_dependencies td
    JOIN dataset_tables p ON p.dataset_table_id=td.parent_dataset_table_id
    JOIN dataset_tables c ON c.dataset_table_id=td.child_dataset_table_id;

    CREATE INDEX idx_artifacts_step ON artifacts(producing_step_id);
    CREATE INDEX idx_artifacts_path ON artifacts(relative_path);
    CREATE INDEX idx_dependencies_parent ON artifact_dependencies(parent_artifact_id);
    CREATE INDEX idx_dependencies_child ON artifact_dependencies(child_artifact_id);
    CREATE INDEX idx_dataset_tables_step ON dataset_tables(producing_step_id);
    """)


def dump_database(con: sqlite3.Connection, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("-- Complete MetGuild SQLite export: provenance plus scientific result data\n")
        handle.write(f"-- Generated: {utc_now()}\n\n")
        for line in con.iterdump():
            handle.write(line)
            handle.write("\n")
    temporary.replace(destination)


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    db_path = Path(args.db).resolve()
    sql_path = Path(args.sql).resolve()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()
    if sql_path.exists():
        sql_path.unlink()

    artifacts = scan_artifacts(args, run_dir, db_path, sql_path)
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA foreign_keys = ON")
    try:
        create_provenance_schema(con)
        populate_provenance(con, args, artifacts, run_dir)
        link_artifacts(con, artifacts, args)
        imported, failed, total_rows = import_scientific_tables(con, run_dir, artifacts)
        dependency_count = rebuild_table_dependencies(con)
        con.execute(
            "INSERT INTO result_import_runs(imported_at,run_directory,imported_tables,failed_tables,imported_rows,table_dependencies) VALUES (?,?,?,?,?,?)",
            (utc_now(), str(run_dir), imported, failed, total_rows, dependency_count),
        )
        create_views_and_indexes(con)
        con.commit()
        dump_database(con, sql_path)
        print("\nProvenance import summary")
        print(f"  Imported tables:   {imported}")
        print(f"  Imported rows:     {total_rows}")
        print(f"  Failed imports:    {failed}")
        print(f"  Table connections: {dependency_count}")
        print(f"  SQLite database:   {db_path}")
        print(f"  Full SQL dump:     {sql_path}")
        return 0 if failed == 0 else 2
    finally:
        con.close()


if __name__ == "__main__":
    raise SystemExit(main())
METGUILD_PROVENANCE_PY

echo "Step 9 complete."
echo ""

echo "============================================================"
echo "MetGuild workflow completed successfully (Steps 1-9)"
echo "============================================================"
echo "Run directory: $OUTPUT_DIR"
echo ""
echo "Main output locations:"
echo "  Candidate-k plots:       $OUTPUT_DIR/02_k_optimization/"
echo "  Clustering:              $OUTPUT_DIR/03_04_clustering/"
if [[ -n "$RUN_TAXONOMY_FILE" ]]; then
    echo "  Cluster taxonomy:        $OUTPUT_DIR/03_04_clustering/taxonomy/"
fi
echo "  Statistical screening:   $OUTPUT_DIR/05_stat_analysis/"
echo "  Correlations:            $OUTPUT_DIR/06_correlation/"
echo "  Network edge lists:      $OUTPUT_DIR/07_network_construction/edge_lists/"
echo "  Network node attributes: $OUTPUT_DIR/07_network_construction/nodes_attributes/"
echo "  SQLite provenance DB:    $PROVENANCE_DB"
echo "  Full SQL dump:           $PROVENANCE_SQL"
echo "============================================================"
