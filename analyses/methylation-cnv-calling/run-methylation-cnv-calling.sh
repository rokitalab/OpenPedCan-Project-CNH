#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOF
Usage: $(basename "$0") --rg-set-file FILE --manifest-file FILE --output-dir DIR --output-prefix PREFIX --array-type TYPE [--combine-segments [--chain-file FILE]] [--run-gistic [--bioassay-ids FILE] [--gistic-output-dir DIR]]
       $(basename "$0") --combine-only --output-dir DIR [--chain-file FILE] [--run-gistic [--bioassay-ids FILE] [--gistic-output-dir DIR]]

Run methylation CNV calling from one QC-filtered RGset, and optionally lift
EPICv1/450k (hg19) GISTIC SEG files to hg38 and combine them with EPICv2 (hg38)
into combined_hg38.gistic.seg in the output directory.

Required options (CNV calling):
  --rg-set-file FILE     RGChannelSet .qs2 file from methylation preprocessing.
  --manifest-file FILE   Manifest with file_name, Bioassay_ID, sample_type, and platform.
  --output-dir DIR       Directory for segment outputs.
  --output-prefix TEXT   Prefix for output files.
  --array-type TYPE      EPIC, EPICv2, or 450k.

Liftover options:
  --combine-segments     After calling, combine all *-gistic.seg files in
                         --output-dir into combined_hg38.gistic.seg.
  --combine-only         Skip CNV calling; only combine existing files. Only
                         --output-dir (and --chain-file) are needed.
  --chain-file FILE      hg19ToHg38 chain file. Required when EPIC or 450k
                         files are combined; not needed for EPICv2 only.

GISTIC options:
  --run-gistic           After combining, run GISTIC2 on combined_hg38.gistic.seg.
                         Implies --combine-segments.
  --bioassay-ids FILE    Optional. Restrict GISTIC to these Bioassay_IDs (one
                         per line).
  --gistic-output-dir DIR
                         GISTIC output directory. Default: <output-dir>/gistic.
EOF
}

RG_SET_FILE=""
MANIFEST_FILE=""
OUTPUT_DIR=""
OUTPUT_PREFIX=""
ARRAY_TYPE=""
CHAIN_FILE=""
COMBINE=0
COMBINE_ONLY=0
RUN_GISTIC=0
BIOASSAY_IDS=""
GISTIC_OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rg-set-file) RG_SET_FILE="${2:-}"; shift 2 ;;
    --manifest-file) MANIFEST_FILE="${2:-}"; shift 2 ;;
    --output-dir) OUTPUT_DIR="${2:-}"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="${2:-}"; shift 2 ;;
    --array-type) ARRAY_TYPE="${2:-}"; shift 2 ;;
    --chain-file) CHAIN_FILE="${2:-}"; shift 2 ;;
    --combine-segments) COMBINE=1; shift ;;
    --combine-only) COMBINE=1; COMBINE_ONLY=1; shift ;;
    --run-gistic) RUN_GISTIC=1; COMBINE=1; shift ;;
    --bioassay-ids) BIOASSAY_IDS="${2:-}"; shift 2 ;;
    --gistic-output-dir) GISTIC_OUTPUT_DIR="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$COMBINE_ONLY" -eq 1 ]]; then
  required_opts=(OUTPUT_DIR)
else
  required_opts=(RG_SET_FILE MANIFEST_FILE OUTPUT_DIR OUTPUT_PREFIX ARRAY_TYPE)
fi

for required in "${required_opts[@]}"; do
  if [[ -z "${!required}" ]]; then
    echo "Missing required option: ${required}" >&2
    usage >&2
    exit 2
  fi
done

[[ -z "$BIOASSAY_IDS" || -f "$BIOASSAY_IDS" ]] || { echo "Bioassay ID file not found: $BIOASSAY_IDS" >&2; exit 2; }
[[ -z "$CHAIN_FILE" || -f "$CHAIN_FILE" ]] || { echo "Chain file not found: $CHAIN_FILE" >&2; exit 2; }

if [[ "$COMBINE_ONLY" -eq 0 ]]; then
  [[ -f "$RG_SET_FILE" ]] || { echo "RGset file not found: $RG_SET_FILE" >&2; exit 2; }
  [[ -f "$MANIFEST_FILE" ]] || { echo "Manifest file not found: $MANIFEST_FILE" >&2; exit 2; }
  [[ "$ARRAY_TYPE" =~ ^(EPIC|EPICv2|450k)$ ]] || { echo "--array-type must be EPIC, EPICv2, or 450k." >&2; exit 2; }

  mkdir -p "$OUTPUT_DIR"

  Rscript --vanilla "$SCRIPT_DIR/scripts/01-methylation-cnv-calls.R" \
    --rg_set_file "$RG_SET_FILE" \
    --manifest_file "$MANIFEST_FILE" \
    --output_basename "$OUTPUT_DIR/$OUTPUT_PREFIX" \
    --array_type "$ARRAY_TYPE"
fi

if [[ "$COMBINE" -eq 1 ]]; then
  [[ -d "$OUTPUT_DIR" ]] || { echo "Output directory not found: $OUTPUT_DIR" >&2; exit 2; }
  liftover_args=(--seg_dir "$OUTPUT_DIR")
  if [[ -n "$CHAIN_FILE" ]]; then
    liftover_args+=(--chain_file "$CHAIN_FILE")
  fi
  Rscript --vanilla "$SCRIPT_DIR/scripts/02-liftover-methylation-segments.R" "${liftover_args[@]}"
  qc_args=(--seg_dir "$OUTPUT_DIR")
  if [[ -n "$BIOASSAY_IDS" ]]; then
    qc_args+=(--bioassay_ids "$BIOASSAY_IDS")
  fi
  Rscript --vanilla "$SCRIPT_DIR/scripts/03-gistic-cohort-qc.R" "${qc_args[@]}"
fi

if [[ "$RUN_GISTIC" -eq 1 ]]; then
  gistic_args=(--seg_file "$OUTPUT_DIR/combined_hg38.gistic.seg" --output_dir "${GISTIC_OUTPUT_DIR:-$OUTPUT_DIR/gistic}")
  if [[ -n "$BIOASSAY_IDS" ]]; then
    gistic_args+=(--bioassay_ids "$BIOASSAY_IDS")
  fi
  bash "$SCRIPT_DIR/scripts/04-run-gistic.sh" "${gistic_args[@]}"
fi
