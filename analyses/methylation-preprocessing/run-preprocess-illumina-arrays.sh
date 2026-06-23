#!/bin/bash
# OPenPedCan 2022
# J Daggett, updated 2026 
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.

# ---- Global parameters ----
MANIFEST_FILE="controls_and_dicer_manifest.tsv"
N_CORES=4
FUNNORM=TRUE
SNP_FILTER=TRUE
OUT_DIR='test-out'
FILE_PREF='test'
OUT_BASE="$OUT_DIR/$FILE_PREF"

mkdir -p $OUT_DIR

printf 'Sorting array types \n\n'

Rscript --vanilla scripts/00-unzip-and-sort.R --base_dir input-test --manifest_file controls_and_dicer_manifest.tsv --output_basename sorted_idats

printf "Start methylation pre-processing...\n\n"


run_preprocess () {
    local DIR=$1
    local LABEL=$2

if [ -d "$DIR" ] && [ "$(ls -A "$DIR")" ]; then
        echo "Processing $LABEL"

        Rscript scripts/01-preprocess-illumina-arrays.R \
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

run_preprocess "sorted_idats_output_dir/IlluminaHumanMethylationEPICv2" "EPICv2"
run_preprocess "sorted_idats_output_dir/IlluminaHumanMethylationEPIC" "EPICv1"
run_preprocess "sorted_idats_output_dir/IlluminaHumanMethylation450k" "450k"

printf "\ncombining array types...\n"

Rscript scripts/02-merge-methyl-matrices.R --output_dir $OUT_DIR --output_prefix $FILE_PREF


printf "\nStart segmentation and CNV calling...\n\n"

run_cnv () {
    local DIR=$1
    local LABEL=$2
    local ARRAY_TYPE=$3

    if [ -d "$DIR" ] && [ "$(ls -A "$DIR")" ]; then
        echo "Running segmentation for $LABEL"

        Rscript --vanilla scripts/03-cnv-calls.R \
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
run_cnv "sorted_idats_output_dir/IlluminaHumanMethylationEPICv2" "EPICv2" "EPICv2" 
run_cnv "sorted_idats_output_dir/IlluminaHumanMethylationEPIC"   "EPICv1" "EPIC"
run_cnv "sorted_idats_output_dir/IlluminaHumanMethylation450k"   "450k"   "450"