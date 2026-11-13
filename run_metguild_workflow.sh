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

Before running MetGuild
-----------------------

1. Activate the MetGuild Conda environment:

   conda activate metguild_workflow

2. Change to the MetGuild repository directory, for example:

   cd /path/to/MetGuild_updated

3. Run the workflow:

   bash run_metguild_workflow.sh <input_matrix> [MAX_K] [OPTIONS]

The Conda environment should be activated before every new shell/session in
which MetGuild is run. The workflow runner does not activate Conda itself.

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
  * Before Step 9, the runner verifies the network manifest, summary, edge
    lists, and node tables and requires at least one retained network edge.

Step 9 - SQLite/SQL/HTML workflow provenance
  * Build a populated SQLite database linking workflow steps, scripts, inputs,
    intermediate outputs, and final network outputs.
  * Import the actual contents of CSV, TSV, and JSON artifacts into SQLite.
  * Write a complete SQL dump that can recreate the populated database.
  * Write a self-contained browser workflow map showing connected tables,
    figures, interactive plots, and metadata across the MetGuild steps.
  * Filter the map by workflow step, result type, clustering method, and
    candidate k; shared upstream source files remain visible for traceability.
  * Select any result node to highlight its recorded upstream/downstream lineage,
    inspect direct inputs/outputs, preview tables/figures, and open the result.

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
      |-- <dataset>_workflow_provenance.sql
      `-- <dataset>_workflow_provenance.html

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
    echo "Activate the MetGuild Conda environment first:" >&2
    echo "  conda activate metguild_workflow" >&2
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
    print("\nMake sure the MetGuild Conda environment is active:")
    print("  conda activate metguild_workflow")
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
# Validate Step 8 outputs before provenance export
# ------------------------------------------------------------
NETWORK_ROOT="$OUTPUT_DIR/07_network_construction"
NETWORK_MANIFEST="$NETWORK_ROOT/network_construction_manifest.json"
NETWORK_SUMMARY="$NETWORK_ROOT/01_edge_list_merge_summary.csv"
NETWORK_EDGE_DIR="$NETWORK_ROOT/edge_lists"
NETWORK_NODE_DIR="$NETWORK_ROOT/nodes_attributes"

echo "Validating Step 8 network outputs before SQL export..."
python - \
    "$NETWORK_MANIFEST" \
    "$NETWORK_SUMMARY" \
    "$NETWORK_EDGE_DIR" \
    "$NETWORK_NODE_DIR" <<'METGUILD_NETWORK_VALIDATE_PY'
import csv
import json
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])
summary_path = Path(sys.argv[2])
edge_dir = Path(sys.argv[3])
node_dir = Path(sys.argv[4])

errors = []

if not manifest_path.is_file():
    errors.append(f"missing network manifest: {manifest_path}")
if not summary_path.is_file():
    errors.append(f"missing network summary: {summary_path}")
if not edge_dir.is_dir():
    errors.append(f"missing edge-list directory: {edge_dir}")
if not node_dir.is_dir():
    errors.append(f"missing node-attribute directory: {node_dir}")

if errors:
    raise SystemExit(
        "Error: Step 8 network validation failed before SQL export:\n  - "
        + "\n  - ".join(errors)
    )

try:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
except Exception as exc:
    raise SystemExit(f"Error: Could not read network manifest {manifest_path}: {exc}")

if manifest.get("status") != "success":
    errors.append("network_construction_manifest.json does not report status='success'")

runs = manifest.get("runs")
if not isinstance(runs, list) or not runs:
    errors.append("network manifest contains no completed network runs")
    runs = []

def resolve_output_path(raw: str) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path
    # Step 8 normally stores absolute paths; relative paths are interpreted
    # relative to the network-construction directory for robustness.
    return (manifest_path.parent / path).resolve()

def count_csv_rows(path: Path) -> tuple[list[str], int]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            return [], 0
        return header, sum(1 for row in reader if any(cell.strip() for cell in row))

total_declared_edges = 0
total_declared_nodes = 0
nonempty_edge_lists = 0
nonempty_node_tables = 0
validated_runs = 0

for index, run in enumerate(runs, start=1):
    label = f"network run #{index}"
    method = run.get("method_code", run.get("method", "?"))
    k = run.get("requested_k", "?")
    label = f"{method}, k={k}"

    if run.get("status") != "success":
        errors.append(f"{label}: status is {run.get('status')!r}, expected 'success'")
        continue

    try:
        declared_edges = int(run.get("final_edges", 0))
        declared_nodes = int(run.get("nodes_total", 0))
    except (TypeError, ValueError):
        errors.append(f"{label}: final_edges/nodes_total are not valid integers")
        continue

    if declared_edges < 0 or declared_nodes < 0:
        errors.append(f"{label}: negative edge/node counts are invalid")
        continue

    outputs = run.get("outputs")
    if not isinstance(outputs, dict):
        errors.append(f"{label}: outputs mapping is missing")
        continue

    edge_meta = None
    node_meta = None
    for name, meta in outputs.items():
        if not isinstance(meta, dict) or not isinstance(meta.get("path"), str):
            continue
        if str(name).endswith("_edge_list.csv"):
            edge_meta = meta
        elif str(name).endswith("_node_attributes.csv"):
            node_meta = meta

    if edge_meta is None:
        errors.append(f"{label}: manifest does not identify an edge-list CSV")
        continue
    if node_meta is None:
        errors.append(f"{label}: manifest does not identify a node-attribute CSV")
        continue

    edge_path = resolve_output_path(edge_meta["path"])
    node_path = resolve_output_path(node_meta["path"])

    if not edge_path.is_file():
        errors.append(f"{label}: missing edge-list file: {edge_path}")
        continue
    if not node_path.is_file():
        errors.append(f"{label}: missing node-attribute file: {node_path}")
        continue
    if edge_path.stat().st_size == 0:
        errors.append(f"{label}: edge-list file is zero bytes: {edge_path}")
        continue
    if node_path.stat().st_size == 0:
        errors.append(f"{label}: node-attribute file is zero bytes: {node_path}")
        continue

    try:
        edge_header, edge_rows = count_csv_rows(edge_path)
        node_header, node_rows = count_csv_rows(node_path)
    except Exception as exc:
        errors.append(f"{label}: could not read network CSV output: {exc}")
        continue

    required_edge_columns = {"node1", "node2", "weight"}
    required_node_columns = {"Node", "Group"}
    if not required_edge_columns.issubset(set(edge_header)):
        errors.append(
            f"{label}: edge list is missing required columns "
            f"{sorted(required_edge_columns - set(edge_header))}: {edge_path}"
        )
    if not required_node_columns.issubset(set(node_header)):
        errors.append(
            f"{label}: node table is missing required columns "
            f"{sorted(required_node_columns - set(node_header))}: {node_path}"
        )

    if edge_rows != declared_edges:
        errors.append(
            f"{label}: manifest reports {declared_edges} edges but "
            f"{edge_path.name} contains {edge_rows} data rows"
        )
    if node_rows != declared_nodes:
        errors.append(
            f"{label}: manifest reports {declared_nodes} nodes but "
            f"{node_path.name} contains {node_rows} data rows"
        )
    if declared_edges > 0 and declared_nodes == 0:
        errors.append(f"{label}: non-empty edge list has no nodes")

    total_declared_edges += declared_edges
    total_declared_nodes += declared_nodes
    nonempty_edge_lists += int(edge_rows > 0)
    nonempty_node_tables += int(node_rows > 0)
    validated_runs += 1

try:
    with summary_path.open("r", encoding="utf-8-sig", newline="") as handle:
        summary_rows = list(csv.DictReader(handle))
except Exception as exc:
    errors.append(f"could not read network summary CSV: {exc}")
    summary_rows = []

if not summary_rows:
    errors.append("network summary CSV contains no data rows")
elif validated_runs and len(summary_rows) != validated_runs:
    errors.append(
        f"network summary has {len(summary_rows)} rows but {validated_runs} "
        "successful network runs were validated"
    )

# A full MetGuild run should not be reported as successful if all filtered
# correlation stages yield zero final network edges. Individual method/k
# solutions are allowed to be empty, but the complete run must contain at
# least one constructed network edge and at least one populated node table.
if total_declared_edges <= 0 or nonempty_edge_lists <= 0:
    errors.append(
        "Step 8 produced zero retained network edges across all clustering "
        "solutions; refusing SQL export and final workflow success"
    )
if total_declared_nodes <= 0 or nonempty_node_tables <= 0:
    errors.append(
        "Step 8 produced no populated network node table; refusing SQL export "
        "and final workflow success"
    )

if errors:
    raise SystemExit(
        "Error: Step 8 network validation failed before SQL export:\n  - "
        + "\n  - ".join(errors)
    )

print(
    "Step 8 network validation passed: "
    f"{validated_runs} solution(s), {total_declared_edges} total edge(s), "
    f"{total_declared_nodes} total node record(s), "
    f"{nonempty_edge_lists} non-empty edge list(s)."
)
METGUILD_NETWORK_VALIDATE_PY
echo ""

# ------------------------------------------------------------
# Step 9: Build populated SQLite/SQL/HTML workflow provenance
# ------------------------------------------------------------
echo "============================================================"
echo "Step 9/9: Building SQLite, SQL, and HTML workflow provenance"
echo "============================================================"

RUN_COMPLETED_AT=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
PROVENANCE_DIR="$OUTPUT_DIR/08_workflow_provenance"
PROVENANCE_DB="$PROVENANCE_DIR/${INPUT_PREFIX}_workflow_provenance.db"
PROVENANCE_SQL="$PROVENANCE_DIR/${INPUT_PREFIX}_workflow_provenance.sql"
PROVENANCE_HTML="$PROVENANCE_DIR/${INPUT_PREFIX}_workflow_provenance.html"
mkdir -p "$PROVENANCE_DIR"
rm -f "$PROVENANCE_DB" "$PROVENANCE_SQL" "$PROVENANCE_HTML"

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
    --html "$PROVENANCE_HTML"
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
import html
import itertools
import json
import re
import sqlite3
import sys
from urllib.parse import quote
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
         "Create the populated SQLite provenance database, complete SQL reconstruction dump, and lineage-first browser HTML report."),
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
    parser = argparse.ArgumentParser(description="Build MetGuild workflow provenance SQLite/SQL/HTML outputs.")
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
    parser.add_argument("--html", required=True)
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
    if suffix == ".html" and "08_workflow_provenance" in parts:
        return "provenance_html"
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


def scan_artifacts(args: argparse.Namespace, run_dir: Path, db_path: Path, sql_path: Path, html_path: Path) -> list[Artifact]:
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

    excluded = {db_path.resolve(), sql_path.resolve(), html_path.resolve()}
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
    for path, kind in (
        (db_path, "sqlite_database"),
        (sql_path, "database_sql"),
        (html_path, "provenance_html"),
    ):
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


def artifact_solution_context(relative_path: str) -> tuple[str | None, int | None]:
    """Return (method_code, k) for a solution-specific artifact when detectable."""
    text = Path(relative_path).as_posix()

    match = re.search(r"/(01_kmeans|02_hierarchical)/(\d+)_clusters(?:/|$)", "/" + text)
    if match:
        method = "KM" if match.group(1) == "01_kmeans" else "HC"
        return method, int(match.group(2))

    match = re.search(r"(?:^|[/_])(\d+)_(KM|HC)(?=[_./-]|$)", text, flags=re.IGNORECASE)
    if match:
        return match.group(2).upper(), int(match.group(1))

    match = re.search(r"(?:^|[/_])(KM|HC)(\d+)(?=[_./-]|$)", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).upper(), int(match.group(2))

    return None, None


def context_filtered_inputs(inputs: list[Artifact], child: Artifact) -> list[Artifact]:
    """Keep global inputs plus inputs from the same clustering method/k as child."""
    child_method, child_k = artifact_solution_context(child.relative_path)
    if child_method is None or child_k is None:
        return inputs

    filtered: list[Artifact] = []
    for parent in inputs:
        parent_method, parent_k = artifact_solution_context(parent.relative_path)
        if parent_method is None or parent_k is None:
            filtered.append(parent)
        elif parent_method == child_method and parent_k == child_k:
            filtered.append(parent)
    return filtered


def lineage_relation(parent: Artifact, child: Artifact) -> str:
    parent_step = STEP_BY_ID.get(parent.producing_step_id) if parent.producing_step_id else None
    child_step = STEP_BY_ID.get(child.producing_step_id) if child.producing_step_id else None
    po = parent_step.order_no if parent_step else None
    co = child_step.order_no if child_step else None

    if co in {2, 3} and po == 1:
        return "candidate_k_source"
    if co == 4 and po in {2, 3}:
        return "cluster_assignment_source"
    if co == 4 and po == 0:
        return "taxonomy_source"
    if co == 5 and po == 1:
        return "feature_matrix_source"
    if co == 5 and po in {2, 3}:
        return "cluster_assignment_source"
    if co in {6, 7, 8} and po == 5:
        return "retained_feature_source"
    if co == 9 and po in {6, 7, 8}:
        return "filtered_correlation_source"
    if co == 9 and po == 4:
        return "taxonomy_annotation_source"
    if co == 10:
        return "indexed_by_provenance"
    return "used_to_generate"


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
            child_inputs = context_filtered_inputs(inputs, child)
            for parent in child_inputs:
                if parent.artifact_id != child.artifact_id:
                    con.execute("INSERT OR IGNORE INTO artifact_dependencies VALUES (?,?,?)",
                                (parent.artifact_id, child.artifact_id, lineage_relation(parent, child)))
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


def html_text(value: Any) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def human_bytes(value: Any) -> str:
    if value is None:
        return ""
    try:
        size = float(value)
    except (TypeError, ValueError):
        return str(value)
    units = ["B", "KB", "MB", "GB", "TB"]
    for unit in units:
        if abs(size) < 1024 or unit == units[-1]:
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return str(value)


def local_artifact_href(relative_path: str, is_external: int) -> str | None:
    if is_external:
        return None
    # The report lives in run_dir/08_workflow_provenance; ../ returns to run_dir.
    return "../" + quote(relative_path, safe="/._-~")


def table_html(headers: Sequence[str], rows: Sequence[Sequence[Any]], *, table_id: str = "", link_path_index: int | None = None, external_index: int | None = None) -> str:
    id_attr = f' id="{html_text(table_id)}"' if table_id else ""
    out = [f"<div class=\"table-wrap\"><table{id_attr}><thead><tr>"]
    out.extend(f"<th>{html_text(h)}</th>" for h in headers)
    out.append("</tr></thead><tbody>")
    for row in rows:
        out.append("<tr>")
        for idx, value in enumerate(row):
            rendered = html_text(value)
            if link_path_index is not None and idx == link_path_index and value not in (None, ""):
                is_external = int(row[external_index]) if external_index is not None and row[external_index] is not None else 0
                href = local_artifact_href(str(value), is_external)
                if href:
                    rendered = f'<a href="{html_text(href)}">{rendered}</a>'
            out.append(f"<td>{rendered}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def read_network_summary(run_dir: Path) -> tuple[list[str], list[list[str]]]:
    path = run_dir / "07_network_construction" / "01_edge_list_merge_summary.csv"
    if not path.is_file():
        return [], []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        rows = list(reader)
    if not rows:
        return [], []
    return [str(x) for x in rows[0]], [[str(x) for x in row] for row in rows[1:]]


def read_candidate_k(run_dir: Path) -> str:
    path = run_dir / "manifests" / "optimal_k_values.json"
    if not path.is_file():
        return ""
    try:
        data = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError):
        return ""
    values = data.get("unique_k_values")
    if isinstance(values, list):
        return ", ".join(str(x) for x in values)
    return ""


def generate_html_report(con: sqlite3.Connection, destination: Path, run_dir: Path) -> None:
    run = con.execute(
        """
        SELECT dataset_name, original_input_path, original_taxonomy_path, run_directory,
               max_k, matrix_type, statistical_test, candidate_source, requested_distance,
               started_at, completed_at, status
        FROM workflow_runs WHERE run_id = 1
        """
    ).fetchone()
    if run is None:
        raise RuntimeError("workflow_runs contains no run record; cannot build HTML provenance report")

    dataset_name, original_input, original_taxonomy, run_directory, max_k, matrix_type, statistical_test, candidate_source, requested_distance, started_at, completed_at, status = run

    step_rows = con.execute(
        """
        SELECT step_id, workflow_order, step_label, step_name, description
        FROM workflow_steps
        WHERE workflow_order < 10
        ORDER BY workflow_order
        """
    ).fetchall()
    steps = [
        {
            "id": int(step_id),
            "number": int(order_no),
            "label": str(label),
            "name": str(name),
            "description": str(description),
        }
        for step_id, order_no, label, name, description in step_rows
    ]
    source_step = {
        "id": 0,
        "number": -1,
        "label": "source",
        "name": "Original source files",
        "description": "Original matrix and taxonomy files supplied to the workflow before staging.",
    }
    step_by_id = {item["id"]: item for item in steps}

    # Imported table dimensions are already catalogued in SQLite.  Reuse them
    # rather than recounting large workflow outputs while building the report.
    table_dimensions = {
        int(artifact_id): (int(row_count), int(column_count))
        for artifact_id, row_count, column_count in con.execute(
            """
            SELECT artifact_id, row_count, column_count
            FROM dataset_tables
            WHERE import_status = 'imported'
            """
        ).fetchall()
    }

    artifact_rows = con.execute(
        """
        SELECT artifact_id, relative_path, artifact_name, artifact_role, artifact_type,
               size_bytes, COALESCE(sha256,''), COALESCE(modified_at,''),
               producing_step_id, is_external
        FROM artifacts
        ORDER BY artifact_id
        """
    ).fetchall()

    def result_type(kind: str) -> str | None:
        if kind == "table":
            return "table"
        if kind == "plot_or_figure":
            return "figure"
        if kind == "interactive_plot":
            return "interactive_plot"
        if kind in {"manifest", "summary_json", "text"}:
            return "metadata"
        return None

    def artifact_path(relative_path: str, is_external: int) -> Path:
        p = Path(relative_path)
        return p if is_external or p.is_absolute() else run_dir / p

    def subtype_for(path: str, kind: str, step_order: int) -> str:
        lower = path.lower()
        if kind == "figure":
            return "Figure"
        if kind == "interactive_plot":
            return "Interactive visualization"
        if kind == "metadata":
            if "optimal_k" in lower or "manifest" in lower:
                return "Candidate-k metadata"
            if "summary" in lower:
                return "Summary metadata"
            return "Metadata"
        if step_order == 0:
            return "Input table"
        if step_order == 1:
            if "score" in lower:
                return "Candidate-k scores"
            if "clean" in lower:
                return "Cleaned feature matrix"
            return "Preprocessing table"
        if step_order in {2, 3}:
            return "Cluster assignments"
        if step_order == 4:
            return "Taxonomy annotation"
        if step_order == 5:
            if "significant" in lower or "_sig" in lower:
                return "Retained features"
            return "Statistical screening"
        if step_order == 6:
            return "Feature-feature correlation"
        if step_order == 7:
            return "Feature-cluster correlation"
        if step_order == 8:
            return "Cluster-cluster correlation"
        if step_order == 9:
            if "edge" in lower:
                return "Network edge list"
            if "node" in lower:
                return "Network node attributes"
            return "Network output"
        return "Result"

    def read_table_preview(path: Path, limit_rows: int = 30, limit_columns: int = 80) -> tuple[list[str], list[list[str]], str | None, bool]:
        try:
            with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
                sample = handle.read(65536)
                handle.seek(0)
                if path.suffix.lower() == ".tsv":
                    delimiter = "\t"
                else:
                    try:
                        delimiter = csv.Sniffer().sniff(sample, delimiters=",\t;|").delimiter
                    except csv.Error:
                        delimiter = ","
                reader = csv.reader(handle, delimiter=delimiter)
                try:
                    header = [str(x) for x in next(reader)]
                except StopIteration:
                    return [], [], None, False
                truncated_columns = len(header) > limit_columns
                display_header = header[:limit_columns]
                rows: list[list[str]] = []
                for row in itertools.islice(reader, limit_rows):
                    values = [str(x) for x in row[:limit_columns]]
                    if len(values) < len(display_header):
                        values.extend([""] * (len(display_header) - len(values)))
                    rows.append(values)
                return display_header, rows, None, truncated_columns
        except Exception as exc:
            return [], [], str(exc), False

    def read_metadata_preview(path: Path, max_chars: int = 12000) -> str:
        try:
            if path.suffix.lower() == ".json":
                value = json.loads(path.read_text(encoding="utf-8-sig", errors="replace"))
                text_value = json.dumps(value, indent=2, ensure_ascii=False)
            else:
                text_value = path.read_text(encoding="utf-8-sig", errors="replace")
            return text_value[:max_chars] + ("\n… preview truncated …" if len(text_value) > max_chars else "")
        except Exception as exc:
            return f"Preview error: {exc}"

    nodes: list[dict[str, Any]] = []
    artifact_to_node: dict[int, str] = {}
    for artifact_id, relative_path, artifact_name, role, artifact_type, size_bytes, sha256, modified_at, producing_step_id, is_external in artifact_rows:
        ui_type = result_type(str(artifact_type))
        if ui_type is None:
            continue
        # Keep the graph focused on scientific workflow products, not the provenance
        # report itself.  External input/taxonomy artifacts are retained as explicit
        # source nodes so a lineage can reach the files supplied before staging.
        if int(is_external):
            step = source_step
        else:
            step = step_by_id.get(int(producing_step_id)) if producing_step_id is not None else None
        if step is None or step["number"] >= 10:
            continue
        method, clusters = artifact_solution_context(str(relative_path))
        source_path = artifact_path(str(relative_path), int(is_external))
        row_count, column_count = table_dimensions.get(int(artifact_id), (None, None))
        columns: list[str] = []
        preview_rows: list[list[str]] = []
        preview_error: str | None = None
        preview_text = ""
        preview_columns_truncated = False
        if not int(is_external) and ui_type == "table" and source_path.is_file():
            columns, preview_rows, preview_error, preview_columns_truncated = read_table_preview(source_path)
        elif not int(is_external) and ui_type == "metadata" and source_path.is_file():
            preview_text = read_metadata_preview(source_path)

        node_id = f"R{int(artifact_id):04d}"
        artifact_to_node[int(artifact_id)] = node_id
        nodes.append({
            "id": node_id,
            "artifact_id": int(artifact_id),
            "step": int(step["number"]),
            "step_label": str(step["label"]),
            "step_name": str(step["name"]),
            "step_description": str(step["description"]),
            "name": str(artifact_name),
            "relative_path": str(relative_path),
            "href": local_artifact_href(str(relative_path), int(is_external)),
            "type": ui_type,
            "subtype": subtype_for(str(relative_path), ui_type, int(step["number"])),
            "method": method or "",
            "clusters": clusters,
            "size_bytes": int(size_bytes) if size_bytes is not None else None,
            "sha256": str(sha256),
            "modified_at": str(modified_at),
            "row_count": row_count,
            "column_count": column_count,
            "columns": columns,
            "preview_rows": preview_rows,
            "preview_error": preview_error,
            "preview_text": preview_text,
            "preview_columns_truncated": preview_columns_truncated,
        })

    edges: list[dict[str, Any]] = []
    edge_keys: set[tuple[str, str, str]] = set()
    for parent_id, child_id, relation in con.execute(
        """
        SELECT parent_artifact_id, child_artifact_id, relation
        FROM artifact_dependencies
        ORDER BY child_artifact_id, parent_artifact_id, relation
        """
    ).fetchall():
        parent = artifact_to_node.get(int(parent_id))
        child = artifact_to_node.get(int(child_id))
        if not parent or not child or parent == child:
            continue
        key = (parent, child, str(relation))
        if key in edge_keys:
            continue
        edge_keys.add(key)
        edges.append({"from": parent, "to": child, "relation": str(relation)})

    # Figures and interactive plots are often sibling outputs of the same analysis
    # rather than explicit children in the file-level provenance graph.  Add a
    # visualization link to the closest table from the same step/method/k so the
    # browser map makes those relationships visible without changing SQLite lineage.
    by_context: dict[tuple[int, str, int | None], list[dict[str, Any]]] = {}
    for node in nodes:
        by_context.setdefault((node["step"], node["method"], node["clusters"]), []).append(node)

    def token_prefix_score(a: str, b: str) -> tuple[int, int]:
        ta = re.split(r"[^A-Za-z0-9]+", Path(a).stem.lower())
        tb = re.split(r"[^A-Za-z0-9]+", Path(b).stem.lower())
        common = 0
        for x, y in zip(ta, tb):
            if x != y:
                break
            common += 1
        return common, -abs(len(ta) - len(tb))

    for context_nodes in by_context.values():
        tables = [n for n in context_nodes if n["type"] == "table"]
        visuals = [n for n in context_nodes if n["type"] in {"figure", "interactive_plot"}]
        if not tables:
            continue
        for visual in visuals:
            parent = max(tables, key=lambda t: token_prefix_score(t["name"], visual["name"]))
            key = (parent["id"], visual["id"], "visualizes_or_summarizes")
            if key not in edge_keys:
                edge_keys.add(key)
                edges.append({"from": parent["id"], "to": visual["id"], "relation": "visualizes_or_summarizes"})

    workflow_data = {
        "run": {
            "dataset": str(dataset_name),
            "input_file": Path(str(original_input)).name,
            "taxonomy_file": Path(str(original_taxonomy)).name if original_taxonomy else "",
            "max_k": int(max_k),
            "matrix_type": str(matrix_type),
            "statistical_test": str(statistical_test),
            "candidate_source": str(candidate_source),
            "requested_distance": str(requested_distance),
            "started_at": str(started_at),
            "completed_at": str(completed_at),
            "status": str(status),
            "run_directory": str(run_directory),
        },
        "nodes": nodes,
        "edges": edges,
        "steps": [
            {
                "number": int(s["number"]),
                "label": str(s["label"]),
                "name": str(s["name"]),
                "description": str(s["description"]),
            }
            for s in [source_step, *steps]
        ],
    }

    # Prevent an embedded JSON string such as </script> from terminating the data
    # element early when a metadata preview contains arbitrary user text.
    workflow_json = json.dumps(workflow_data, ensure_ascii=False).replace("</", "<\\/")

    css = r"""
:root {
  color-scheme: light dark;
  --bg: #f3f5f8;
  --panel: #ffffff;
  --text: #17202a;
  --muted: #5d6d7e;
  --border: #d8dee8;
  --accent: #2c5cc5;
  --accent-soft: #e8efff;
  --table: #0f766e;
  --figure: #a16207;
  --interactive: #7e22ce;
  --metadata: #475569;
  --edge: #8793a5;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #10141b;
    --panel: #171d26;
    --text: #eef2f7;
    --muted: #aab4c2;
    --border: #303947;
    --accent: #7aa2ff;
    --accent-soft: #243352;
    --table: #4fd1c5;
    --figure: #f6c453;
    --interactive: #c084fc;
    --metadata: #a8b4c5;
    --edge: #697789;
  }
}
* { box-sizing: border-box; }
body { margin: 0; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background: var(--bg); color: var(--text); }
header { padding: 20px 24px 16px; background: var(--panel); border-bottom: 1px solid var(--border); }
h1 { margin: 0 0 6px; font-size: clamp(1.35rem, 3vw, 2rem); }
.run-meta { display: flex; flex-wrap: wrap; gap: 8px 16px; color: var(--muted); font-size: 0.9rem; }
.controls { display: grid; grid-template-columns: minmax(180px, 1fr) repeat(4, minmax(130px, 180px)) auto auto; gap: 10px; padding: 12px 24px; background: var(--panel); border-bottom: 1px solid var(--border); align-items: end; }
.control { min-width: 0; }
.control label { display: block; color: var(--muted); font-size: 0.75rem; margin-bottom: 4px; }
input, select, button { width: 100%; min-height: 40px; font: inherit; border: 1px solid var(--border); border-radius: 8px; background: var(--panel); color: var(--text); padding: 8px 10px; }
button { cursor: pointer; width: auto; min-width: 42px; }
button:hover, button:focus-visible { border-color: var(--accent); outline: none; }
.layout { display: grid; grid-template-columns: minmax(0, 1fr) minmax(340px, 440px); gap: 0; min-height: 720px; }
.graph-panel { min-width: 0; border-right: 1px solid var(--border); }
.graph-toolbar { display: flex; flex-wrap: wrap; justify-content: space-between; gap: 8px; padding: 10px 16px; background: var(--panel); border-bottom: 1px solid var(--border); }
.legend { display: flex; flex-wrap: wrap; gap: 12px; font-size: 0.82rem; color: var(--muted); align-items: center; }
.legend span::before { content: ""; display: inline-block; width: 10px; height: 10px; border-radius: 3px; margin-right: 5px; background: var(--legend-color); }
.graph-scroll { overflow: auto; height: 700px; background-image: radial-gradient(var(--border) 0.7px, transparent 0.7px); background-size: 18px 18px; }
#graph { display: block; min-width: 100%; }
.node { cursor: pointer; }
.node rect { fill: var(--panel); stroke-width: 2; rx: 10; }
.node.table rect { stroke: var(--table); }
.node.figure rect { stroke: var(--figure); }
.node.interactive_plot rect { stroke: var(--interactive); }
.node.metadata rect { stroke: var(--metadata); }
.node.selected rect { stroke: var(--accent); stroke-width: 4; filter: drop-shadow(0 3px 5px rgba(0,0,0,.18)); }
.node.lineage-related rect { stroke-width: 3; }
.node text { pointer-events: none; fill: var(--text); }
.node .node-title { font-weight: 700; font-size: 12px; }
.node .node-meta { font-size: 10px; fill: var(--muted); }
.step-label { font-weight: 800; font-size: 13px; fill: var(--text); }
.step-description { font-size: 10px; fill: var(--muted); }
.edge { stroke: var(--edge); stroke-width: 1.5; fill: none; opacity: 0.48; marker-end: url(#arrow); }
.edge.highlight { stroke: var(--accent); stroke-width: 3; opacity: 1; }
.details { min-width: 0; background: var(--panel); padding: 18px; overflow: auto; height: 756px; }
.details h2 { margin: 0 0 8px; font-size: 1.15rem; word-break: break-word; }
.badges { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 12px; }
.badge { border: 1px solid var(--border); border-radius: 999px; padding: 3px 8px; color: var(--muted); font-size: 0.75rem; }
.detail-grid { display: grid; grid-template-columns: 110px minmax(0,1fr); gap: 7px 10px; font-size: 0.84rem; margin: 12px 0; }
.detail-grid dt { color: var(--muted); }
.detail-grid dd { margin: 0; overflow-wrap: anywhere; }
.action-link { display: inline-block; margin: 6px 8px 10px 0; padding: 8px 11px; border-radius: 8px; background: var(--accent); color: white; text-decoration: none; }
.relationships { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin: 12px 0; }
.relationships h3, .preview h3 { font-size: 0.9rem; margin: 0 0 7px; }
.relation-list { border: 1px solid var(--border); border-radius: 8px; padding: 8px; min-height: 54px; }
.relation-list button { display: block; width: 100%; text-align: left; margin-bottom: 5px; font-size: 0.78rem; min-height: 34px; }
.preview { margin-top: 14px; }
.preview img { max-width: 100%; height: auto; border: 1px solid var(--border); border-radius: 8px; background: white; }
.table-wrap { overflow: auto; max-height: 390px; border: 1px solid var(--border); border-radius: 8px; }
table { border-collapse: collapse; width: max-content; min-width: 100%; font-size: 0.75rem; }
th, td { padding: 6px 8px; border-bottom: 1px solid var(--border); border-right: 1px solid var(--border); max-width: 260px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
th { position: sticky; top: 0; background: var(--panel); z-index: 1; text-align: left; }
pre { white-space: pre-wrap; word-break: break-word; max-height: 390px; overflow: auto; padding: 10px; border: 1px solid var(--border); border-radius: 8px; background: var(--bg); font-size: 0.75rem; }
.empty { color: var(--muted); padding: 20px; text-align: center; }
.summary-cards { display: flex; flex-wrap: wrap; gap: 8px; }
.summary-card { border: 1px solid var(--border); border-radius: 8px; padding: 5px 8px; font-size: 0.78rem; color: var(--muted); }
.trace-note { margin: 0 24px 12px; color: var(--muted); font-size: .82rem; }
@media (max-width: 980px) {
  .controls { grid-template-columns: 1fr 1fr; }
  .layout { grid-template-columns: 1fr; }
  .graph-panel { border-right: 0; border-bottom: 1px solid var(--border); }
  .graph-scroll { height: 520px; }
  .details { height: auto; min-height: 440px; }
}
@media (max-width: 560px) {
  header, .controls { padding-left: 12px; padding-right: 12px; }
  .controls { grid-template-columns: 1fr; }
  .relationships { grid-template-columns: 1fr; }
}
"""

    js = r"""
(() => {
  const rootData = JSON.parse(document.getElementById('workflowData').textContent);
  const nodesById = new Map(rootData.nodes.map(node => [node.id, node]));
  const incoming = new Map();
  const outgoing = new Map();
  rootData.nodes.forEach(node => { incoming.set(node.id, []); outgoing.set(node.id, []); });
  rootData.edges.forEach(edge => {
    if (incoming.has(edge.to)) incoming.get(edge.to).push(edge);
    if (outgoing.has(edge.from)) outgoing.get(edge.from).push(edge);
  });

  const search = document.getElementById('search');
  const stepFilter = document.getElementById('stepFilter');
  const typeFilter = document.getElementById('typeFilter');
  const methodFilter = document.getElementById('methodFilter');
  const clusterFilter = document.getElementById('clusterFilter');
  const svg = document.getElementById('graph');
  const edgesGroup = document.getElementById('edges');
  const nodesGroup = document.getElementById('nodes');
  const details = document.getElementById('details');
  const summaryCards = document.getElementById('summaryCards');
  let selectedId = null;
  let zoom = 1;
  let baseWidth = 1200;
  let baseHeight = 700;

  const escapeHtml = value => String(value ?? '').replace(/[&<>'"]/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[ch]));
  const formatBytes = bytes => {
    if (!Number.isFinite(bytes)) return 'Unknown';
    const units = ['B','KB','MB','GB'];
    let value = bytes;
    let index = 0;
    while (value >= 1024 && index < units.length - 1) { value /= 1024; index += 1; }
    return `${value.toFixed(index === 0 ? 0 : 1)} ${units[index]}`;
  };
  const shortName = name => name.length > 38 ? `${name.slice(0, 35)}…` : name;

  function collectAncestors(seedIds) {
    const result = new Set(seedIds);
    const stack = [...seedIds];
    while (stack.length) {
      const id = stack.pop();
      (incoming.get(id) || []).forEach(edge => {
        if (!result.has(edge.from)) { result.add(edge.from); stack.push(edge.from); }
      });
    }
    return result;
  }

  function collectConnected(seedId) {
    const nodes = new Set([seedId]);
    const edges = new Set();
    const up = [seedId];
    const down = [seedId];
    while (up.length) {
      const id = up.pop();
      (incoming.get(id) || []).forEach(edge => {
        const key = `${edge.from}|${edge.to}|${edge.relation}`;
        edges.add(key);
        if (!nodes.has(edge.from)) { nodes.add(edge.from); up.push(edge.from); }
      });
    }
    while (down.length) {
      const id = down.pop();
      (outgoing.get(id) || []).forEach(edge => {
        const key = `${edge.from}|${edge.to}|${edge.relation}`;
        edges.add(key);
        if (!nodes.has(edge.to)) { nodes.add(edge.to); down.push(edge.to); }
      });
    }
    return {nodes, edges};
  }

  function fillFilters() {
    rootData.steps.forEach(step => {
      const option = document.createElement('option');
      option.value = String(step.number);
      option.textContent = `${step.label}. ${step.name}`;
      stepFilter.appendChild(option);
    });
    [...new Set(rootData.nodes.map(n => n.type))].sort().forEach(value => {
      const option = document.createElement('option'); option.value = value; option.textContent = value.replaceAll('_', ' '); typeFilter.appendChild(option);
    });
    [...new Set(rootData.nodes.map(n => n.method).filter(Boolean))].sort().forEach(value => {
      const option = document.createElement('option'); option.value = value; option.textContent = value === 'KM' ? 'KMeans' : value === 'HC' ? 'Hierarchical' : value; methodFilter.appendChild(option);
    });
    [...new Set(rootData.nodes.map(n => n.clusters).filter(v => v !== null))].sort((a,b) => a-b).forEach(value => {
      const option = document.createElement('option'); option.value = String(value); option.textContent = String(value); clusterFilter.appendChild(option);
    });
  }

  function currentNodes() {
    const query = search.value.trim().toLowerCase();
    const method = methodFilter.value;
    const cluster = clusterFilter.value;

    let contextIds = null;
    if (method || cluster) {
      const seeds = rootData.nodes.filter(node => {
        if (method && node.method !== method) return false;
        if (cluster && String(node.clusters) !== cluster) return false;
        return Boolean(node.method || node.clusters !== null);
      }).map(node => node.id);
      contextIds = collectAncestors(seeds);
      seeds.forEach(id => contextIds.add(id));
    }

    return rootData.nodes.filter(node => {
      if (contextIds && !contextIds.has(node.id)) return false;
      if (stepFilter.value && String(node.step) !== stepFilter.value) return false;
      if (typeFilter.value && node.type !== typeFilter.value) return false;
      if (query) {
        const haystack = [node.name, node.relative_path, node.type, node.subtype, node.method, node.step_name, node.step_label, node.clusters].join(' ').toLowerCase();
        if (!haystack.includes(query)) return false;
      }
      return true;
    });
  }

  function applyZoom() {
    svg.style.width = `${baseWidth * zoom}px`;
    svg.style.height = `${baseHeight * zoom}px`;
  }

  function renderGraph() {
    const visible = currentNodes();
    const visibleIds = new Set(visible.map(node => node.id));
    const grouped = new Map();
    visible.forEach(node => {
      if (!grouped.has(node.step)) grouped.set(node.step, []);
      grouped.get(node.step).push(node);
    });
    grouped.forEach(values => values.sort((a,b) => a.relative_path.localeCompare(b.relative_path)));

    const activeSteps = [...grouped.keys()].sort((a,b) => a-b);
    const xGap = 320;
    const yGap = 96;
    const nodeWidth = 278;
    const nodeHeight = 70;
    const positions = new Map();
    let maxRows = 1;

    activeSteps.forEach((step, stepIndex) => {
      const list = grouped.get(step);
      maxRows = Math.max(maxRows, list.length);
      list.forEach((node, rowIndex) => {
        positions.set(node.id, {x: 30 + stepIndex * xGap, y: 82 + rowIndex * yGap});
      });
    });

    baseWidth = Math.max(900, 70 + Math.max(activeSteps.length, 1) * xGap);
    baseHeight = Math.max(650, 130 + maxRows * yGap);
    svg.setAttribute('viewBox', `0 0 ${baseWidth} ${baseHeight}`);
    applyZoom();
    edgesGroup.innerHTML = '';
    nodesGroup.innerHTML = '';

    const lineage = selectedId ? collectConnected(selectedId) : {nodes:new Set(), edges:new Set()};

    rootData.edges.forEach(edge => {
      if (!visibleIds.has(edge.from) || !visibleIds.has(edge.to)) return;
      const from = positions.get(edge.from);
      const to = positions.get(edge.to);
      if (!from || !to) return;
      const x1 = from.x + nodeWidth;
      const y1 = from.y + nodeHeight / 2;
      const x2 = to.x;
      const y2 = to.y + nodeHeight / 2;
      let d;
      if (Math.abs(x2 - x1) < 24) {
        const loopX = Math.max(x1, x2) + 50;
        d = `M ${x1} ${y1} C ${loopX} ${y1}, ${loopX} ${y2}, ${x2} ${y2}`;
      } else {
        const bend = Math.max(45, Math.abs(x2 - x1) / 2);
        const sign = x2 >= x1 ? 1 : -1;
        d = `M ${x1} ${y1} C ${x1 + sign*bend} ${y1}, ${x2 - sign*bend} ${y2}, ${x2} ${y2}`;
      }
      const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      path.setAttribute('d', d);
      const edgeKey = `${edge.from}|${edge.to}|${edge.relation}`;
      path.setAttribute('class', `edge${lineage.edges.has(edgeKey) ? ' highlight' : ''}`);
      const title = document.createElementNS('http://www.w3.org/2000/svg', 'title');
      title.textContent = edge.relation.replaceAll('_',' ');
      path.appendChild(title);
      edgesGroup.appendChild(path);
    });

    activeSteps.forEach((step, stepIndex) => {
      const x = 30 + stepIndex * xGap;
      const info = rootData.steps.find(item => item.number === step);
      const title = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      title.setAttribute('x', x); title.setAttribute('y', 28); title.setAttribute('class', 'step-label');
      title.textContent = info ? `${info.label}. ${info.name}` : `Step ${step}`;
      nodesGroup.appendChild(title);
      const description = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      description.setAttribute('x', x); description.setAttribute('y', 47); description.setAttribute('class', 'step-description');
      description.textContent = info ? shortName(info.description) : '';
      nodesGroup.appendChild(description);
    });

    visible.forEach(node => {
      const pos = positions.get(node.id);
      const group = document.createElementNS('http://www.w3.org/2000/svg', 'g');
      const related = lineage.nodes.has(node.id) && node.id !== selectedId ? ' lineage-related' : '';
      group.setAttribute('class', `node ${node.type}${node.id === selectedId ? ' selected' : ''}${related}`);
      group.setAttribute('transform', `translate(${pos.x}, ${pos.y})`);
      group.setAttribute('tabindex', '0');
      group.setAttribute('role', 'button');
      group.setAttribute('aria-label', `${node.name}, ${node.type}, step ${node.step_label}`);

      const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      rect.setAttribute('width', nodeWidth); rect.setAttribute('height', nodeHeight); group.appendChild(rect);

      const title = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      title.setAttribute('x', 12); title.setAttribute('y', 24); title.setAttribute('class', 'node-title');
      title.textContent = shortName(node.name); group.appendChild(title);

      const meta = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      meta.setAttribute('x', 12); meta.setAttribute('y', 44); meta.setAttribute('class', 'node-meta');
      const parts = [node.type.replaceAll('_',' ')];
      if (node.method) parts.push(node.method);
      if (node.clusters !== null) parts.push(`k=${node.clusters}`);
      if (node.row_count !== null) parts.push(`${node.row_count} rows`);
      meta.textContent = parts.join(' • '); group.appendChild(meta);

      const pathText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
      pathText.setAttribute('x', 12); pathText.setAttribute('y', 60); pathText.setAttribute('class', 'node-meta');
      pathText.textContent = shortName(node.relative_path); group.appendChild(pathText);

      const activate = () => selectNode(node.id);
      group.addEventListener('click', activate);
      group.addEventListener('keydown', event => {
        if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); activate(); }
      });
      nodesGroup.appendChild(group);
    });

    const visibleEdges = rootData.edges.filter(e => visibleIds.has(e.from) && visibleIds.has(e.to)).length;
    summaryCards.innerHTML = `<span class="summary-card">${visible.length} visible results</span><span class="summary-card">${visibleEdges} visible connections</span>`;
  }

  function relationButton(edge, direction) {
    const targetId = direction === 'in' ? edge.from : edge.to;
    const target = nodesById.get(targetId);
    if (!target) return '';
    return `<button type="button" data-node-id="${escapeHtml(targetId)}"><strong>${escapeHtml(target.name)}</strong><br><span>${escapeHtml(edge.relation.replaceAll('_',' '))}</span></button>`;
  }

  function renderPreview(node) {
    if (node.type === 'table') {
      if (node.preview_error) return `<p>Preview error: ${escapeHtml(node.preview_error)}</p>`;
      if (!node.columns.length) return '<p>No preview rows were available.</p>';
      const head = node.columns.map(column => `<th>${escapeHtml(column)}</th>`).join('');
      const body = node.preview_rows.map(row => `<tr>${row.map(value => `<td title="${escapeHtml(value)}">${escapeHtml(value)}</td>`).join('')}</tr>`).join('');
      const columnNote = node.preview_columns_truncated ? ' Columns are limited in this inline preview.' : '';
      return `<div class="table-wrap"><table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div><p class="badge">Preview limited to the first 30 rows.${columnNote}</p>`;
    }
    if (node.type === 'figure') {
      if (/\.pdf$/i.test(node.name)) return `<p>PDF preview is not embedded. Open the file using the button above.</p>`;
      return `<img src="${escapeHtml(node.href)}" alt="Preview of ${escapeHtml(node.name)}">`;
    }
    if (node.type === 'interactive_plot') {
      return `<p>This result is an interactive HTML visualization. Open it in a new browser tab using the button above.</p>`;
    }
    if (node.type === 'metadata') {
      return `<pre>${escapeHtml(node.preview_text || 'No metadata preview available.')}</pre>`;
    }
    return '<p>No inline preview is available for this file type.</p>';
  }

  function selectNode(id) {
    selectedId = id;
    const node = nodesById.get(id);
    if (!node) return;
    const inEdges = incoming.get(id) || [];
    const outEdges = outgoing.get(id) || [];
    details.innerHTML = `
      <h2>${escapeHtml(node.name)}</h2>
      <div class="badges">
        <span class="badge">Step ${escapeHtml(node.step_label)}: ${escapeHtml(node.step_name)}</span>
        <span class="badge">${escapeHtml(node.type.replaceAll('_',' '))}</span>
        <span class="badge">${escapeHtml(node.subtype)}</span>
        ${node.method ? `<span class="badge">${escapeHtml(node.method)}</span>` : ''}
        ${node.clusters !== null ? `<span class="badge">k=${node.clusters}</span>` : ''}
      </div>
      ${node.href ? `<a class="action-link" href="${escapeHtml(node.href)}" target="_blank" rel="noopener">Open original result</a>` : ''}
      <dl class="detail-grid">
        <dt>Path</dt><dd>${escapeHtml(node.relative_path)}</dd>
        <dt>File size</dt><dd>${formatBytes(node.size_bytes)}</dd>
        <dt>Rows</dt><dd>${node.row_count === null ? 'Not applicable or unknown' : node.row_count}</dd>
        <dt>Columns</dt><dd>${node.column_count === null ? 'Not applicable or unknown' : node.column_count}</dd>
        <dt>SHA-256</dt><dd><code>${escapeHtml(node.sha256 || 'not recorded')}</code></dd>
        <dt>Description</dt><dd>${escapeHtml(node.step_description)}</dd>
      </dl>
      <div class="relationships">
        <section><h3>Inputs to this result</h3><div class="relation-list">${inEdges.length ? inEdges.map(edge => relationButton(edge, 'in')).join('') : '<span class="empty">No upstream result recorded.</span>'}</div></section>
        <section><h3>Results produced from this</h3><div class="relation-list">${outEdges.length ? outEdges.map(edge => relationButton(edge, 'out')).join('') : '<span class="empty">No downstream result recorded.</span>'}</div></section>
      </div>
      <section class="preview"><h3>Result preview</h3>${renderPreview(node)}</section>`;
    details.querySelectorAll('[data-node-id]').forEach(button => button.addEventListener('click', () => selectNode(button.dataset.nodeId)));
    renderGraph();
  }

  document.getElementById('runMeta').innerHTML = [
    `Dataset: <strong>${escapeHtml(rootData.run.dataset)}</strong>`,
    `Input: <strong>${escapeHtml(rootData.run.input_file)}</strong>`,
    `Matrix: <strong>${escapeHtml(rootData.run.matrix_type)}</strong>`,
    `Test: <strong>${escapeHtml(rootData.run.statistical_test)}</strong>`,
    `Completed: <strong>${escapeHtml(rootData.run.completed_at)}</strong>`
  ].map(item => `<span>${item}</span>`).join('');

  [search, stepFilter, typeFilter, methodFilter, clusterFilter].forEach(control => control.addEventListener('input', renderGraph));
  document.getElementById('zoomIn').addEventListener('click', () => { zoom = Math.min(1.8, zoom + 0.1); applyZoom(); });
  document.getElementById('zoomOut').addEventListener('click', () => { zoom = Math.max(0.5, zoom - 0.1); applyZoom(); });

  fillFilters();
  renderGraph();
  if (rootData.nodes.length) selectNode(rootData.nodes[0].id);
})();
"""

    document = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MetGuild Workflow Results Map - {html_text(dataset_name)}</title>
<style>{css}</style>
</head>
<body>
<header>
  <h1>MetGuild Workflow Results Map</h1>
  <div class="run-meta" id="runMeta"></div>
</header>
<section class="controls" aria-label="Workflow filters">
  <div class="control"><label for="search">Search results</label><input id="search" type="search" placeholder="Filename, path, method, result type"></div>
  <div class="control"><label for="stepFilter">Workflow step</label><select id="stepFilter"><option value="">All steps</option></select></div>
  <div class="control"><label for="typeFilter">Result type</label><select id="typeFilter"><option value="">All types</option></select></div>
  <div class="control"><label for="methodFilter">Method</label><select id="methodFilter"><option value="">All methods</option></select></div>
  <div class="control"><label for="clusterFilter">Candidate k</label><select id="clusterFilter"><option value="">All candidate k values</option></select></div>
  <button id="zoomOut" type="button" aria-label="Zoom out">−</button>
  <button id="zoomIn" type="button" aria-label="Zoom in">+</button>
</section>
<p class="trace-note">Choose a clustering method and candidate <em>k</em> to keep that solution together with the shared upstream files needed to trace it back to the input. Select any result node to highlight its complete recorded lineage.</p>
<main class="layout">
  <section class="graph-panel" aria-label="Interactive result lineage graph">
    <div class="graph-toolbar">
      <div class="legend">
        <span style="--legend-color:var(--table)">Table</span>
        <span style="--legend-color:var(--figure)">Figure</span>
        <span style="--legend-color:var(--interactive)">Interactive plot</span>
        <span style="--legend-color:var(--metadata)">JSON/metadata</span>
      </div>
      <div class="summary-cards" id="summaryCards"></div>
    </div>
    <div class="graph-scroll" id="graphScroll">
      <svg id="graph" role="img" aria-label="MetGuild workflow graph of result tables and figures">
        <defs>
          <marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto" markerUnits="strokeWidth">
            <path d="M0,0 L8,4 L0,8 z" fill="var(--edge)"></path>
          </marker>
        </defs>
        <g id="edges"></g>
        <g id="nodes"></g>
      </svg>
    </div>
  </section>
  <aside class="details" id="details" aria-live="polite">
    <div class="empty">Select a result node to inspect its file, preview, inputs, and outputs.</div>
  </aside>
</main>
<script id="workflowData" type="application/json">{workflow_json}</script>
<script>{js}</script>
</body>
</html>'''

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(document, encoding="utf-8")
    temporary.replace(destination)

def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    db_path = Path(args.db).resolve()
    sql_path = Path(args.sql).resolve()
    html_path = Path(args.html).resolve()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()
    if sql_path.exists():
        sql_path.unlink()
    if html_path.exists():
        html_path.unlink()

    artifacts = scan_artifacts(args, run_dir, db_path, sql_path, html_path)
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
        generate_html_report(con, html_path, run_dir)
        print("\nProvenance import summary")
        print(f"  Imported tables:   {imported}")
        print(f"  Imported rows:     {total_rows}")
        print(f"  Failed imports:    {failed}")
        print(f"  Table connections: {dependency_count}")
        print(f"  SQLite database:   {db_path}")
        print(f"  Full SQL dump:     {sql_path}")
        print(f"  HTML browser view: {html_path}")
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
echo "  Network validation:      passed before provenance export"
echo "  SQLite provenance DB:    $PROVENANCE_DB"
echo "  Full SQL dump:           $PROVENANCE_SQL"
echo "  Browser provenance view: $PROVENANCE_HTML"
echo ""
echo "Open the HTML report in a browser, for example:"
echo "  firefox $PROVENANCE_HTML"
echo "============================================================"
