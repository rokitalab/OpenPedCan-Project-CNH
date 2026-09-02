#!/bin/bash 

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# This will be turned off in CI
SUBSET=${OPENPBTA_SUBSET:-1}

scratch_path="../../scratch/"
data_dir="../../data"



# Run R script to generate JSON file
Rscript --vanilla 00-rare-cns-select-tumors.R

# Subtype rare tumors using fusions
Rscript -e "rmarkdown::render('01-subtype-using-fusions.Rmd', clean = TRUE)"

# Subtype rare tumors using methylation
Rscript -e "rmarkdown::render('02-subtype-using-methylation.Rmd', clean = TRUE)"

