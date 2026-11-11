#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------
# MetGuild full workflow runner (Steps 1-8)
# ------------------------------------------------------------

WORKFLOW_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'USAGE'
MetGuild Workflow
=================

Usage:
  bash run_metguild_workflow.sh <input_matrix> [MAX_K] [OPTIONS]

By default, MetGuild runs the complete workflow (Steps 1-8).
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

  --full                  Explicitly request the complete Steps 1-8 workflow.
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
  `-- 07_network_construction/
      |-- edge_lists/
      `-- nodes_attributes/

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
    echo "Step ${number}/8: ${label}"
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
    echo "Step 3b/8: Taxonomy annotation"
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

echo "============================================================"
echo "MetGuild workflow completed successfully (Steps 1-8)"
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
