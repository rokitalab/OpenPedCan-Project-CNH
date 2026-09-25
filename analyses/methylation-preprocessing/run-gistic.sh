#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") [--seg_file FILE] [--bioassay_ids FILE] [--output_dir DIR]

Run GISTIC2 on a combined hg38 SEG file.

Options:
  --seg_file FILE   Combined hg38 GISTIC SEG file.
                    Default: test-out/combined_hg38.gistic.seg
  --bioassay_ids FILE
                    One bioassay ID per line. The SEG is filtered to these
                    samples before GISTIC runs.
                    Default: dicer1_methylation_bioassay_ids.txt
  --output_dir DIR  Directory for GISTIC output.
                    Default: test-out/test-gistic
  -h, --help        Show this help message.
EOF
}

SEG_FILE="test-out/combined_hg38.gistic.seg"
BIOASSAY_IDS="dicer1_methylation_bioassay_ids.txt"
OUTPUT_DIR="test-out/test-gistic"
GISTIC_DIR="${GISTIC_DIR:-/home/rstudio/gistic_install/share/gistic2-2.0.23-0}"
REFGENE_FILE="${REFGENE_FILE:-$GISTIC_DIR/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --seg_file)
      if [[ $# -lt 2 || -z "$2" ]]; then
        echo "Error: --seg_file requires a value." >&2
        exit 2
      fi
      SEG_FILE="$2"
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
    --bioassay_ids)
      if [[ $# -lt 2 || -z "$2" ]]; then
        echo "Error: --bioassay_ids requires a value." >&2
        exit 2
      fi
      BIOASSAY_IDS="$2"
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

if [[ ! -f "$SEG_FILE" ]]; then
  echo "Error: combined SEG file not found: $SEG_FILE" >&2
  exit 2
fi

if [[ ! -f "$BIOASSAY_IDS" ]]; then
  echo "Error: bioassay ID file not found: $BIOASSAY_IDS" >&2
  exit 2
fi

if [[ ! -f "$REFGENE_FILE" ]]; then
  echo "Error: hg38 reference-gene file not found: $REFGENE_FILE" >&2
  exit 2
fi

mkdir -p "$OUTPUT_DIR"

FILTERED_SEG_FILE="$OUTPUT_DIR/$(basename "${SEG_FILE%.seg}").dicer1.seg"
Rscript scripts/05-filter-gistic-seg.R "$SEG_FILE" "$BIOASSAY_IDS" "$FILTERED_SEG_FILE"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/opt/mcr/v83/runtime/glnxa64:/opt/mcr/v83/bin/glnxa64:/opt/mcr/v83/sys/os/glnxa64"
export XAPPLRESDIR=/opt/mcr/v83/X11/app-defaults

"$GISTIC_DIR/gp_gistic2_from_seg" \
  -b "$OUTPUT_DIR" \
  -seg "$FILTERED_SEG_FILE" \
  -refgene "$REFGENE_FILE" \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.7 \
  -conf 0.99 \
  -armpeel 1 \
  -ta 0.2 \
  -td 0.2
