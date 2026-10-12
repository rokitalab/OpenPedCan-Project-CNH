#!/bin/bash
# Download the parquet copies of the large release tables (see scripts/README.md)
# from the OpenPedCan data release s3 bucket, check their md5sums and link them
# in data/parquet/. They are not part of the release md5sum.txt, so
# download-data.sh does not fetch them.
#
# Usage: bash scripts/download-parquet.sh [table ...]
#   With no arguments every file in md5sum-parquet.txt is downloaded; otherwise
#   only the named tables (e.g. snv-consensus-plus-hotspots fusion-annoFuse).
#   Set OPENPEDCAN_RELEASE (default v15) and OPENPEDCAN_URL as for the other
#   download scripts.

set -e
set -o pipefail

URL=${OPENPEDCAN_URL:-https://s3.amazonaws.com/bti-openaccess-us-east-1-prd-opc}
RELEASE=${OPENPEDCAN_RELEASE:-v15}

# Set the working directory to the base directory of the project
cd "$(dirname "${BASH_SOURCE[0]}")/.."
BASEDIR="$(pwd)"

PARQUET_DIR="$BASEDIR/data/$RELEASE/parquet"
mkdir -p "$PARQUET_DIR" "$BASEDIR/data/parquet"

# The md5sum file provides our single point of truth for which files are in the
# parquet folder; the manifest lists the tables, row counts and source files.
curl -fsS --retry 3 "$URL/$RELEASE/parquet/md5sum-parquet.txt" \
  -o "$PARQUET_DIR/md5sum-parquet.txt"
curl -fsS --retry 3 "$URL/$RELEASE/parquet/manifest.tsv" \
  -o "$PARQUET_DIR/manifest.tsv"

ALL_FILES=(`tr -s ' ' < "$PARQUET_DIR/md5sum-parquet.txt" | cut -d ' ' -f 2`)

if [ "$#" -gt 0 ]
then
  FILES=()
  for table in "$@"
  do
    file="${table%.parquet}.parquet"
    if [[ " ${ALL_FILES[*]} " != *" $file "* ]]
    then
      echo "$file is not in md5sum-parquet.txt for $RELEASE" >&2
      exit 1
    fi
    FILES+=("$file")
  done
else
  FILES=("${ALL_FILES[@]}")
fi

for file in "${FILES[@]}"
do
  if [ ! -e "$PARQUET_DIR/$file" ]
  then
    echo "Downloading $file"
    curl -fsS --retry 3 "$URL/$RELEASE/parquet/$file" -o "$PARQUET_DIR/$file"
  fi
done

# Check the md5s of the files that were requested
echo "Checking MD5 hashes..."
cd "$PARQUET_DIR"
for file in "${FILES[@]}"
do
  grep -E "[[:space:]]${file}\$" md5sum-parquet.txt | md5sum -c -
done
cd "$BASEDIR"

# Make symlinks in data/parquet/ to the files in the just downloaded release
# folder, as download-data.sh does for the other release files.
for file in "${FILES[@]}"
do
  ln -sfn "../$RELEASE/parquet/$file" "data/parquet/$file"
done
ln -sfn "../$RELEASE/parquet/manifest.tsv" data/parquet/manifest.tsv

echo "Linked ${#FILES[@]} parquet file(s) in data/parquet/"
