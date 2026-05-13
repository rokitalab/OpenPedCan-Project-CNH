#!/bin/bash
# OPenPedCan 2022
# J Daggett, updated 2026 
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.

#printf 'Sorting array types \n\n'

#Rscript --vanilla scripts/00-unzip-and-sort.R --base_dir input-test --manifest_file controls_and_dicer_manifest.tsv --output_basename sorted_idats

printf "Start methylation pre-processing...\n\n"

# ---- Global parameters ----
MANIFEST_FILE="controls_and_dicer_manifest.tsv"
N_CORES=4
FUNNORM=TRUE
SNP_FILTER=TRUE
OUT_DIR='test-out'
OUT_BASE="$OUT_DIR/test"

mkdir -p $OUT_DIR

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