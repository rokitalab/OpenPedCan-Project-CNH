#!/bin/bash
# Run the molecular-subtyping module: methylation calls, lesion tables,
# validation and lesion fallback, LGG lesions, medulloblastoma, the solid
# tumors, the final table, the G-CIMP scores, the comparison with the legacy
# modules and the medulloblastoma methylation plots (UMAP and PCA).
#
# The G-CIMP scores and the plots need the beta values
# (data/parquet/methyl-beta-values.parquet: bash scripts/download-parquet.sh
# methyl-beta-values).
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

# the G-CIMP scores and the methylation plots read the beta values (about 15 GB for the
# release; the CI jobs set METHYL_PARQUET_RELEASE=testing for a small testing copy, made by
# scripts/make-methyl-testing-parquet.R, while the other parquet tables are those of v15)
if [ ! -e ../../data/parquet/methyl-beta-values.parquet ]; then
  OPENPEDCAN_RELEASE=${METHYL_PARQUET_RELEASE:-${OPENPEDCAN_RELEASE:-v15}} \
    bash ../../scripts/download-parquet.sh methyl-beta-values
fi

Rscript --vanilla 01-build-methylation-calls.R
Rscript --vanilla 02-prepare-lesion-tables.R
Rscript --vanilla 03-validate-subtypes.R
Rscript --vanilla 04-lgg-lesions.R
Rscript --vanilla 05-validate-mb.R
Rscript --vanilla 06-subtype-solid-tumors.R
Rscript --vanilla 07-compile-subtypes.R
Rscript --vanilla 08-score-gcimp.R
Rscript --vanilla 09-compare-to-legacy.R
Rscript --vanilla 10-plot-mb-methylation.R
Rscript --vanilla 11-plot-mb-highlight.R
