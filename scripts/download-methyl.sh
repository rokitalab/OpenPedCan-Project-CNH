#!/bin/bash

set -e
set -o pipefail

# download methylation files from the OpenPedCan data release s3 bucket
# (includes hg19ToHg38.over.chain.gz, used to lift EPIC/450k segments to hg38)

URL=${OPENPEDCAN_URL:-https://s3.amazonaws.com/bti-openaccess-us-east-1-prd-opc}
RELEASE=${OPENPEDCAN_RELEASE:-v15}

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# If RUN_LOCAL is used, the time-intensive steps are skipped because they cannot
# be run on a local computer -- the idea is that setting RUN_LOCAL=1 will allow for
# local testing running/testing of all other steps
RUN_LOCAL=${RUN_LOCAL:-0}

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

# check if the release folder exists, if not, create a release folder
[ ! -d "$BASEDIR/data/$RELEASE/" ] && mkdir $BASEDIR/data/$RELEASE/

# The md5sum file provides our single point of truth for which files are in a release.
curl --create-dirs $URL/$RELEASE/methyl-md5sum.txt -o $BASEDIR/data/$RELEASE/methyl-md5sum.txt -z $BASEDIR/data/$RELEASE/methyl-md5sum.txt

FILES=(`tr -s ' ' < $BASEDIR/data/$RELEASE/methyl-md5sum.txt | cut -d ' ' -f 2` release-notes.md)

for file in "${FILES[@]}"
do
  if [ ! -e "$BASEDIR/data/$RELEASE/$file" ]
  then
    echo "Downloading $file"
    curl --create-dirs $URL/$RELEASE/$file -o $BASEDIR/data/$RELEASE/$file
  fi
done

# Download the hg19 -> hg38 liftover chain file from UCSC if it was not in the
# release md5sum list. It is used to lift EPIC/450k segments to hg38.
CHAIN="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
if [ ! -e "$BASEDIR/data/$RELEASE/${CHAIN##*/}" ]
then
  echo "Downloading ${CHAIN##*/}"
  curl -fsS --retry 3 --create-dirs $CHAIN -o $BASEDIR/data/$RELEASE/${CHAIN##*/}
fi

#check md5sum
cd $BASEDIR/data/$RELEASE
echo "Checking MD5 hashes..."
md5sum -c methyl-md5sum.txt
cd $BASEDIR

# Make symlinks in data/ to the files in the just downloaded release folder.
for file in "${FILES[@]}"
do
  ln -sfn $RELEASE/$file data/$file
done

# Same for the chain file, which is not always listed in methyl-md5sum.txt.
ln -sfn $RELEASE/${CHAIN##*/} data/${CHAIN##*/}
