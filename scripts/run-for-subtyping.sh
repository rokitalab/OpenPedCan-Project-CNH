#!/bin/sh

set -e
set -o pipefail

printf "Start molecular subtyping...\n\n"

# Set locations for s3 bucket that contains release file
URL="s3://d3b-openaccess-us-east-1-prd-pbta/open-targets"
RELEASE="v14"

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -
  
analyses_dir="$BASEDIR/analyses"
data_dir="$BASEDIR/data/$RELEASE"
scratch_dir="$BASEDIR/scratch"

# Compile all the files that need to be included in the release in one place
# in the scratch directory
release_dir="${scratch_dir}/analysis-files-pre-release"
mkdir -p ${release_dir}

## Run subtyping modules

# Run molecular subtyping (methylation-first; includes the MB expression
# classifier and SHH subtyping, and the chordoma SMARCB1 status). The legacy
# per-diagnosis modules are archived in analyses/archive.
echo "Run molecular subtyping"
cd ${analyses_dir}/molecular-subtyping
bash run-molecular-subtyping.sh

# Run integrate subtyping
echo "Run integrate subtyping"
cd ${analyses_dir}/molecular-subtyping-integrate
bash run-subtyping-integrate.sh

# Copy over integrated subtyping - the *FULL* histology file
cp ${analyses_dir}/molecular-subtyping-integrate/results/histologies.tsv ${data_dir}
cp ${analyses_dir}/molecular-subtyping-integrate/results/histologies.tsv ${release_dir}

# Create the independent sample list using the *FULL* histology file (i.e. - histologies.tsv)
echo "Create independent sample list"
cd ${analyses_dir}/independent-samples
bash run-independent-samples.sh

# Copy over independent specimen lists
cp ${analyses_dir}/independent-samples/results/independent-specimens.* ${data_dir}
cp ${analyses_dir}/independent-samples/results/independent-specimens.* ${release_dir}


# Create an md5sum file for all the files in the directories where the analysis
# files are compiled
cd ${release_dir}
# Remove old md5sum release file if it exists
rm -f analysis_files_release_md5sum.txt
# Create a new md5sum release file
md5sum * > analysis_files_release_md5sum.txt

# Upload all release files s3 bucket in their respective folders
#aws s3 cp ${release_dir}/ $URL/$RELEASE/ --recursive

printf "\nDone running molecular subtyping...\n\n"
