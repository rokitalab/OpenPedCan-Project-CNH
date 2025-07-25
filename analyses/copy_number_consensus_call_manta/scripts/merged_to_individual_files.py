## Nhat Duong & J Shapiro
## 2019 - 2020


# Imports in the pep8 order https://www.python.org/dev/peps/pep-0008/#imports
# Standard library
import argparse
import subprocess
import sys
import os

# Related third party
import numpy as np
import pandas as pd

## Define the extensions
MANTA_EXT = '.manta'
CNVKIT_EXT = '.cnvkit'
FREEC_EXT = '.freec'

# Define the column headers for IDs
MANTA_ID_HEADER = 'Kids.First.Biospecimen.ID.Tumor'
CNVKIT_ID_HEADER = 'ID'
FREEC_ID_HEADER = 'Kids_First_Biospecimen_ID'

parser = argparse.ArgumentParser(description="""This script splits CNV files
                                                into one per sample. It also
                                                prints a snakemake config file
                                                to the specified filename.""")
parser.add_argument('--manta', required=True,
                    help='path to the manta file')
parser.add_argument('--cnvkit', required=True,
                    help='path to the cnvkit file')
parser.add_argument('--freec', required=True,
                    help='path to the freec file')
parser.add_argument('--snake', required=True,
                    help='path for snakemake config file')
parser.add_argument('--scratch', required=True,
                    help='directory for scratch files')
parser.add_argument('--uncalled', required=True,
                    help='path for the table of sample-caller outputs removed and not called for too many CNVs')
parser.add_argument('--maxcnvs', default=2500,
                    help='samples with more than 2500 cnvs are set to blank')
parser.add_argument('--cnvsize', default=3000,
                    help='cnv cutoff size in base pairs')
parser.add_argument('--freecp', default=0.01,
                    help='p-value cutoff for freec')
parser.add_argument('--histologies', required = True,
                    help = 'path to the histology file')

args = parser.parse_args()

print("Starting CNV file processing...")
print(f"Input files: manta={args.manta}, cnvkit={args.cnvkit}, freec={args.freec}")

## Pandas load/read files in
print("Step 1: Loading input files...")

# First, get sample list from other files to know which samples we need
print("  Loading histologies file...")
histologies = pd.read_csv(args.histologies, sep="\t", dtype=str)
print(f"  ✓ Loaded histologies: {len(histologies)} records")

print("  Loading cnvkit file...")
merged_cnvkit = pd.read_csv(args.cnvkit, delimiter='\t', dtype=str)
print(f"  ✓ Loaded cnvkit: {len(merged_cnvkit)} records")

print("  Loading freec file...")
merged_freec = pd.read_csv(args.freec, delimiter='\t', dtype=str)
print(f"  ✓ Loaded freec: {len(merged_freec)} records")

cnvkit_samples = set(merged_cnvkit[CNVKIT_ID_HEADER])
freec_samples = set(merged_freec[FREEC_ID_HEADER])
print(f"  Found {len(cnvkit_samples)} cnvkit samples, {len(freec_samples)} freec samples")

# Filtering for WGS samples 
WGS_all_samples = set(histologies[(histologies["experimental_strategy"] == "WGS") & (histologies["sample_type"] == "Tumor")]['Kids_First_Biospecimen_ID'])
print(f"  Found {len(WGS_all_samples)} WGS tumor samples")

print("Step 2: Processing manta file in chunks...")
# Process manta file in chunks instead of loading entirely
manta_samples = set()
sample_data = {}  # Store data for each sample

chunk_size = 500000  # Adjust based on available memory
chunk_count = 0
for chunk in pd.read_csv(args.manta, delimiter='\t', dtype=str, chunksize=chunk_size):
    chunk_count += 1
    if chunk_count % 10 == 0:  # Progress update every 10 chunks
        print(f"  Processing chunk {chunk_count}...")
    
    # Get samples from this chunk
    chunk_samples = set(chunk[MANTA_ID_HEADER])
    manta_samples.update(chunk_samples)
    
    # Process samples we're interested in
    for sample in chunk_samples.intersection(WGS_all_samples):
        sample_chunk = chunk[chunk[MANTA_ID_HEADER] == sample]
        if sample in sample_data:
            sample_data[sample] = pd.concat([sample_data[sample], sample_chunk])
        else:
            sample_data[sample] = sample_chunk

print(f"  ✓ Processed {chunk_count} chunks from manta file")
print(f"  ✓ Found {len(manta_samples)} manta samples, {len(sample_data)} WGS samples with manta data")

# Now continue with the rest of the processing using sample_data instead of merged_manta
all_samples = sorted(manta_samples | cnvkit_samples | freec_samples)
WGS_all_samples_to_run = WGS_all_samples.intersection(all_samples)
print(f"  ✓ Total samples to process: {len(WGS_all_samples_to_run)}")

print("Step 3: Creating output directories...")
## Define and create assumed directories
scratch_d = args.scratch
manta_d = os.path.join(scratch_d, 'manta_manta')
cnvkit_d = os.path.join(scratch_d, 'cnvkit_cnvkit')
freec_d = os.path.join(scratch_d, 'freec_freec')
if not os.path.exists(manta_d):
    os.makedirs(manta_d)
if not os.path.exists(cnvkit_d):
    os.makedirs(cnvkit_d)
if not os.path.exists(freec_d):
    os.makedirs(freec_d)
print(f"  ✓ Created directories: {manta_d}, {cnvkit_d}, {freec_d}")

print("Step 4: Processing individual samples...")
bad_calls = []
processed_count = 0

## Loop through each sample, search for that sample in each of the three dataframes,
## and create a file of the sample in each directory
for sample in WGS_all_samples_to_run:
    processed_count += 1
    if processed_count % 100 == 0:  # Progress update every 100 samples
        print(f"  Processed {processed_count}/{len(WGS_all_samples_to_run)} samples...")
    
    # Use pre-processed data instead of querying the full dataframe
    manta_export = sample_data.get(sample, pd.DataFrame())
    
    with open(os.path.join(manta_d, sample + MANTA_EXT), 'w') as file_out:
        if len(manta_export) <= args.maxcnvs and len(manta_export) > 0:
            manta_export.to_csv(file_out, sep='\t', index=False)
        else:
            bad_calls.append(sample + "\tmanta\n")
    
    cnvkit_export = merged_cnvkit.loc[merged_cnvkit[CNVKIT_ID_HEADER] == sample]
    with open(os.path.join(cnvkit_d, sample + CNVKIT_EXT), 'w') as file_out:
        if cnvkit_export.shape[0] <= args.maxcnvs and cnvkit_export.shape[0] > 0:
            cnvkit_export.to_csv(file_out, sep='\t', index=False)
        else:
            bad_calls.append(sample + "\tcnvkit\n")

    freec_export = merged_freec.loc[merged_freec[FREEC_ID_HEADER] == sample]
    with open(os.path.join(freec_d, sample + FREEC_EXT), 'w') as file_out:
        if freec_export.shape[0] <= args.maxcnvs and freec_export.shape[0] > 0:
            freec_export.to_csv(file_out, sep='\t', index=False)
        else:
            bad_calls.append(sample + "\tfreec\n")

print(f"  ✓ Processed all {processed_count} samples")

print("Step 5: Creating Snakemake config file...")
## Make the Snakemake config file. Write all of the sample names into the config file
with open(args.snake, 'w') as file:
    file.write('samples:' + '\n')
    for sample in WGS_all_samples_to_run:
        file.write('  ' + str(sample) + ':' + '\n')

    ## Define the extension for the config file
    file.write('manta_ext: ' + MANTA_EXT + '\n')
    file.write('cnvkit_ext: ' + CNVKIT_EXT + '\n')
    file.write('freec_ext: ' + FREEC_EXT + '\n')

    ## Define location for python scripts and scratch
    file.write('scripts: ' + os.path.dirname(os.path.realpath(__file__)) + '\n')
    file.write('scratch: ' + scratch_d + '\n')

    ## Define the size cutoff and freec's pval cut off.
    file.write('size_cutoff: ' + str(args.cnvsize) + '\n')
    file.write('freec_pval: ' + str(args.freecp) + '\n')

print(f"  ✓ Created Snakemake config: {args.snake}")

print("Step 6: Writing uncalled samples file...")
## Write out the bad calls file
bad_calls.sort()
with open(args.uncalled, 'w') as file:
    file.write("sample\tcaller\n")
    file.writelines(bad_calls)

print(f"  ✓ Created uncalled samples file: {args.uncalled} ({len(bad_calls)} bad calls)")
print("✓ Script completed successfully!")
