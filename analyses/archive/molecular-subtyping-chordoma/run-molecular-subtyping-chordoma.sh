#!/bin/bash

# Written originally Chante Bethell 2019 
# (Adapted for this module by Candace Savonen 2020)
#
# Run `00-subset-files-for-chordoma.R` and
# `01-Subtype-chordoma.Rmd` sequentially.

set -e
set -o pipefail


# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

echo "selecting samples for subtyping"
Rscript 00-chordoma-select-pathology-dx.R

echo "subsetting data files"
Rscript --vanilla 00-subset-files-for-chordoma.R

echo "creating relevant data file"
Rscript -e "rmarkdown::render('01-Subtype-chordoma.Rmd', clean = TRUE)"

