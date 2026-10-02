#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOF
Usage: $(basename "$0") --seg-dir DIR [--chain-file FILE]

Lift EPICv1 and 450k methylation GISTIC SEG files from hg19 to hg38 and
combine them with EPICv2 (hg38) GISTIC SEG files.

Options:
  --seg-dir DIR        Directory containing per-array *-gistic.seg files.
  --chain-file FILE    hg19ToHg38 chain file. Required when EPICv1 or 450k
                       files are present; not needed for an EPICv2-only run.
  -h, --help           Show this help text.
EOF
}

SEG_DIR=""
CHAIN_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --seg-dir) SEG_DIR="${2:-}"; shift 2 ;;
    --chain-file) CHAIN_FILE="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$SEG_DIR" ]] || { echo "--seg-dir is required." >&2; usage >&2; exit 2; }
[[ -d "$SEG_DIR" ]] || { echo "SEG directory not found: $SEG_DIR" >&2; exit 2; }
[[ -z "$CHAIN_FILE" || -f "$CHAIN_FILE" ]] || { echo "Chain file not found: $CHAIN_FILE" >&2; exit 2; }

args=(--seg_dir "$SEG_DIR")
if [[ -n "$CHAIN_FILE" ]]; then
  args+=(--chain_file "$CHAIN_FILE")
fi

Rscript --vanilla "$SCRIPT_DIR/scripts/02-liftover-methylation-segments.R" "${args[@]}"
