#!/bin/bash

set -e
set -o pipefail

# Download the methylation testing data (IDATs and manifest) used by CI.

URL=${OPENPEDCAN_URL:-https://s3.amazonaws.com/bti-openaccess-us-east-1-prd-opc}
TESTING_DIR=${TESTING_DIR:-methyl_testing}

# Set the working directory to the base directory of the project
cd "$(dirname "${BASH_SOURCE[0]}")/.."

mkdir -p "data/$TESTING_DIR"
cd "data/$TESTING_DIR"

# The md5sum file provides our single point of truth for which files are in the testing set.
curl -fsS --retry 3 "$URL/$TESTING_DIR/md5sum.txt" -o md5sum.txt

FILES=(`tr -s ' ' < md5sum.txt | cut -d ' ' -f 2`)

for file in "${FILES[@]}"
do
  if [ ! -e "$file" ]
  then
    echo "Downloading $file"
    curl -fsS --retry 3 "$URL/$TESTING_DIR/$file" -o "$file"
  fi
done

echo "Checking MD5 hashes..."
md5sum -c md5sum.txt

# Download the hg19 -> hg38 liftover chain file from UCSC if it was not in the
# md5sum list. It is used to lift EPIC/450k segments to hg38.
CHAIN="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
if [ ! -e "${CHAIN##*/}" ]
then
  echo "Downloading ${CHAIN##*/}"
  curl -fsS --retry 3 "$CHAIN" -o "${CHAIN##*/}"
fi
