#!/bin/bash
# OPenPedCan 2022
# J Daggett, updated 2026
set -e
set -o pipefail

# Resolve the script location without changing the caller's working directory,
# so relative input and output paths are interpreted from the command line.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
usage() {
    cat <<EOF
Usage: $(basename "$0") --manifest_file FILE --input_dir DIR --output_dir DIR --output_prefix PREFIX

Sort and preprocess Illumina IDAT files.

Required options:
  --manifest_file FILE    Manifest containing file_name and Bioassay_ID columns.
  --input_dir DIR         Directory containing the input IDAT files.
  --output_dir DIR        Directory for sorted IDATs and result files.
  --output_prefix PREFIX  Prefix for result files within --output_dir.
  -h, --help              Show this help message.
EOF
}

MANIFEST_FILE=""
INPUT_DIR=""
OUTPUT_DIR=""
OUTPUT_PREFIX=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest_file)
            if [[ $# -lt 2 || -z "$2" ]]; then
                echo "Error: --manifest_file requires a value." >&2
                exit 2
            fi
            MANIFEST_FILE="$2"
            shift 2
            ;;
        --input_dir)
            if [[ $# -lt 2 || -z "$2" ]]; then
                echo "Error: --input_dir requires a value." >&2
                exit 2
            fi
            INPUT_DIR="$2"
            shift 2
            ;;
        --output_dir)
            if [[ $# -lt 2 || -z "$2" ]]; then
                echo "Error: --output_dir requires a value." >&2
                exit 2
            fi
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --output_prefix)
            if [[ $# -lt 2 || -z "$2" ]]; then
                echo "Error: --output_prefix requires a value." >&2
                exit 2
            fi
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$MANIFEST_FILE" || -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: --manifest_file, --input_dir, --output_dir, and --output_prefix are required." >&2
    usage >&2
    exit 2
fi

if [[ ! -f "$MANIFEST_FILE" ]]; then
    echo "Error: manifest file not found: $MANIFEST_FILE" >&2
    exit 2
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: input directory not found: $INPUT_DIR" >&2
    exit 2
fi

N_CORES=4
FUNNORM=TRUE
SNP_FILTER=TRUE
OUT_BASE="$OUTPUT_DIR/$OUTPUT_PREFIX"
SORTED_IDATS_BASE="$OUTPUT_DIR/${OUTPUT_PREFIX}-sorted-idats"
SORTED_IDATS_DIR="${SORTED_IDATS_BASE}_output_dir"

mkdir -p "$OUTPUT_DIR"

printf 'Sorting array types \n\n'

Rscript --vanilla "$SCRIPT_DIR/scripts/00-unzip-and-sort.R" \
    --base_dir "$INPUT_DIR" \
    --manifest_file "$MANIFEST_FILE" \
    --output_basename "$SORTED_IDATS_BASE"

printf "Start methylation pre-processing...\n\n"

run_preprocess () {
    local DIR=$1
    local LABEL=$2

if [ -d "$DIR" ] && [ "$(ls -A "$DIR")" ]; then
        echo "Processing $LABEL"

        Rscript "$SCRIPT_DIR/scripts/01-preprocess-illumina-arrays.R" \
            --base_dir "$DIR" \
            --funnorm "$FUNNORM" \
            --snp_filter "$SNP_FILTER" \
            --manifest_file "$MANIFEST_FILE" \
            --n_cores "$N_CORES" \
            --output_basename "$OUT_BASE"
    else
        echo "Skipping $LABEL (missing or empty)"
    fi
}

run_preprocess "$SORTED_IDATS_DIR/IlluminaHumanMethylationEPICv2" "EPICv2"
run_preprocess "$SORTED_IDATS_DIR/IlluminaHumanMethylationEPIC" "EPICv1"
run_preprocess "$SORTED_IDATS_DIR/IlluminaHumanMethylation450k" "450k"

printf "\ncombining array types...\n"

Rscript --vanilla "$SCRIPT_DIR/scripts/02-merge-methyl-matrices.R" \
    --output_dir "$OUTPUT_DIR" \
    --output_prefix "$OUTPUT_PREFIX"

printf "\nStart segmentation and CNV calling...\n\n"

run_cnv () {
    local DIR=$1
    local LABEL=$2
    local ARRAY_TYPE=$3

    if [ -d "$DIR" ] && [ "$(ls -A "$DIR")" ]; then
        echo "Running segmentation for $LABEL"

        Rscript --vanilla "$SCRIPT_DIR/scripts/03-cnv-calls.R" \
            --base_dir "$DIR" \
            --manifest_file "$MANIFEST_FILE" \
            --n_cores "$N_CORES" \
            --output_basename "$OUT_BASE" \
            --array_type "$ARRAY_TYPE"
    else
        echo "Skipping segmentation for $LABEL (missing or empty)"
    fi
}

# ---- Run CNV step for each array ----
run_cnv "$SORTED_IDATS_DIR/IlluminaHumanMethylationEPICv2" "EPICv2" "EPICv2"
run_cnv "$SORTED_IDATS_DIR/IlluminaHumanMethylationEPIC" "EPICv1" "EPIC"
run_cnv "$SORTED_IDATS_DIR/IlluminaHumanMethylation450k" "450k" "450"

has_array () {
    local DIR=$1
    [ -d "$DIR" ] && [ "$(ls -A "$DIR")" ]
}

# Standardize legacy-array GISTIC SEG files from hg19 to hg38 before combining
# them with EPICv2 calls, which already use hg38 coordinates.
if has_array "$SORTED_IDATS_DIR/IlluminaHumanMethylationEPIC" || \
   has_array "$SORTED_IDATS_DIR/IlluminaHumanMethylation450k"; then
    LIFTOVER_DIR="$OUTPUT_DIR/liftover"
    CHAIN_FILE="$LIFTOVER_DIR/hg19ToHg38.over.chain"
    mkdir -p "$LIFTOVER_DIR"

    if [ ! -s "$CHAIN_FILE" ]; then
        curl -fsSL "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
            | gzip -dc > "$CHAIN_FILE"
    fi

    Rscript --vanilla "$SCRIPT_DIR/scripts/04-liftover.R" \
        --seg_dir "$OUTPUT_DIR" \
        --chain_file "$CHAIN_FILE"
else
    Rscript --vanilla "$SCRIPT_DIR/scripts/04-liftover.R" --seg_dir "$OUTPUT_DIR"
fi
