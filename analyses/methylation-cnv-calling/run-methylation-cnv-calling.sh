#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOF
Usage: $(basename "$0") --rg-set-file FILE --manifest-file FILE --output-dir DIR --output-prefix PREFIX --array-type TYPE

Run methylation CNV calling from one QC-filtered RGset.

Required options:
  --rg-set-file FILE     RGChannelSet .qs2 file from methylation preprocessing.
  --manifest-file FILE   Manifest with file_name, Bioassay_ID, sample_type, and platform.
  --output-dir DIR       Directory for segment outputs.
  --output-prefix TEXT   Prefix for output files.
  --array-type TYPE      EPIC, EPICv2, or 450k.
EOF
}

RG_SET_FILE=""
MANIFEST_FILE=""
OUTPUT_DIR=""
OUTPUT_PREFIX=""
ARRAY_TYPE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rg-set-file) RG_SET_FILE="${2:-}"; shift 2 ;;
    --manifest-file) MANIFEST_FILE="${2:-}"; shift 2 ;;
    --output-dir) OUTPUT_DIR="${2:-}"; shift 2 ;;
    --output-prefix) OUTPUT_PREFIX="${2:-}"; shift 2 ;;
    --array-type) ARRAY_TYPE="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for required in RG_SET_FILE MANIFEST_FILE OUTPUT_DIR OUTPUT_PREFIX ARRAY_TYPE; do
  if [[ -z "${!required}" ]]; then
    echo "Missing required option: ${required}" >&2
    usage >&2
    exit 2
  fi
done

[[ -f "$RG_SET_FILE" ]] || { echo "RGset file not found: $RG_SET_FILE" >&2; exit 2; }
[[ -f "$MANIFEST_FILE" ]] || { echo "Manifest file not found: $MANIFEST_FILE" >&2; exit 2; }
[[ "$ARRAY_TYPE" =~ ^(EPIC|EPICv2|450k)$ ]] || { echo "--array-type must be EPIC, EPICv2, or 450k." >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

Rscript --vanilla "$SCRIPT_DIR/scripts/01-methylation-cnv-calls.R" \
  --rg_set_file "$RG_SET_FILE" \
  --manifest_file "$MANIFEST_FILE" \
  --output_basename "$OUTPUT_DIR/$OUTPUT_PREFIX" \
  --array_type "$ARRAY_TYPE"
