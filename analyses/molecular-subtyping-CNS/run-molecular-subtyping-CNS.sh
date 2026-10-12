#!/bin/bash
# Run the molecular-subtyping-CNS module: methylation calls, lesion tables,
# validation and lesion fallback, LGG lesions, medulloblastoma, the final table
# and the comparison with the legacy modules.
#
# The large release tables are read from data/parquet/ when it exists
# (scripts/download-parquet.sh) and from the release TSV / RDS files otherwise.

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

Rscript --vanilla 01-build-methylation-calls.R
Rscript --vanilla 02-prepare-lesion-tables.R
Rscript --vanilla 03-validate-subtypes.R
Rscript --vanilla 04-lgg-lesions.R
Rscript --vanilla 05-validate-mb.R
Rscript --vanilla 06-compile-subtypes.R
Rscript --vanilla 07-compare-to-legacy.R
