#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOS
Usage: $(basename "$0") --seg_file FILE --output_dir DIR [--bioassay_ids FILE]

Run GISTIC2 on a combined hg38 SEG file.

Options:
  --seg_file FILE       Combined hg38 GISTIC SEG file (combined_hg38.gistic.seg).
  --output_dir DIR      Directory for GISTIC output.
  --bioassay_ids FILE   Optional. One bioassay ID per line. The SEG is filtered
                        to these samples before GISTIC runs.
  -h, --help            Show this help message.

Environment:
  GISTIC_DIR     GISTIC2 install directory containing gp_gistic2_from_seg.
                 Default: /home/rstudio/gistic_install (or its
                 share/gistic2-2.0.23-0 subdirectory if that is where it is).
  REFGENE_FILE   hg38 reference-gene .mat file.
                 Default: \$GISTIC_DIR/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
EOS
}

SEG_FILE=""
OUTPUT_DIR=""
BIOASSAY_IDS=""
if [[ -z "${GISTIC_DIR:-}" ]]; then
  # The openpedcanverse image unpacks GISTIC2 into gistic_install; some
  # installs keep the same files under share/gistic2-2.0.23-0.
  for candidate in /home/rstudio/gistic_install \
                   /home/rstudio/gistic_install/share/gistic2-2.0.23-0; do
    if [[ -x "$candidate/gp_gistic2_from_seg" ]]; then
      GISTIC_DIR="$candidate"
      break
    fi
  done
fi
GISTIC_DIR="${GISTIC_DIR:-/home/rstudio/gistic_install}"
REFGENE_FILE="${REFGENE_FILE:-$GISTIC_DIR/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --seg_file) SEG_FILE="${2:-}"; shift 2 ;;
    --output_dir) OUTPUT_DIR="${2:-}"; shift 2 ;;
    --bioassay_ids) BIOASSAY_IDS="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$SEG_FILE" && -n "$OUTPUT_DIR" ]] || { echo "Error: --seg_file and --output_dir are required." >&2; usage >&2; exit 2; }
[[ -f "$SEG_FILE" ]] || { echo "Error: combined SEG file not found: $SEG_FILE" >&2; exit 2; }
[[ -z "$BIOASSAY_IDS" || -f "$BIOASSAY_IDS" ]] || { echo "Error: bioassay ID file not found: $BIOASSAY_IDS" >&2; exit 2; }
[[ -x "$GISTIC_DIR/gp_gistic2_from_seg" ]] || { echo "Error: gp_gistic2_from_seg not found in GISTIC_DIR: $GISTIC_DIR" >&2; exit 2; }
[[ -f "$REFGENE_FILE" ]] || { echo "Error: hg38 reference-gene file not found: $REFGENE_FILE" >&2; exit 2; }

mkdir -p "$OUTPUT_DIR"

GISTIC_SEG_FILE="$SEG_FILE"
if [[ -n "$BIOASSAY_IDS" ]]; then
  # GISTIC has no sample-subset option and its statistics are computed over the
  # samples in the SEG, so define the cohort by subsetting the SEG first.
  GISTIC_SEG_FILE="$OUTPUT_DIR/$(basename "${SEG_FILE%.seg}").filtered.seg"
  awk -F'\t' -v OFS='\t' '
    FNR == NR { gsub(/^[ \t]+|[ \t\r]+$/, ""); if ($0 != "") ids[$0] = 1; next }
    FNR == 1 {
      for (i = 1; i <= NF; i++) if ($i == "Sample") col = i
      if (!col) { print "Error: SEG file must contain a Sample column." > "/dev/stderr"; exit 2 }
      print; next
    }
    $col in ids { print; kept[$col] = 1; n++ }
    END {
      if (!col) exit 2
      for (id in ids) requested++
      for (id in kept) matched++
      if (!n) { print "Error: none of the requested bioassay IDs are in the SEG." > "/dev/stderr"; exit 2 }
      printf "Kept %d segments for %d of %d requested bioassay IDs\n", n, matched, requested > "/dev/stderr"
    }
  ' "$BIOASSAY_IDS" "$SEG_FILE" > "$GISTIC_SEG_FILE"
fi

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/opt/mcr/v83/runtime/glnxa64:/opt/mcr/v83/bin/glnxa64:/opt/mcr/v83/sys/os/glnxa64"
export XAPPLRESDIR=/opt/mcr/v83/X11/app-defaults

"$GISTIC_DIR/gp_gistic2_from_seg" \
  -b "$OUTPUT_DIR" \
  -seg "$GISTIC_SEG_FILE" \
  -refgene "$REFGENE_FILE" \
  -genegistic 1 \
  -smallmem 1 \
  -broad 1 \
  -brlen 0.7 \
  -conf 0.99 \
  -armpeel 1 \
  -ta 0.2 \
  -td 0.2
