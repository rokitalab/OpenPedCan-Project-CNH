# OpenPedCan Methylation Analysis
This analysis is was last run with OpenPedCan data release `v12` and is now run upstream in CAVATICA.

## Purpose

Preprocess probe hybridization intensity values of selected methylated and unmethylated cytosine (CpG) loci into usable methylation measurements for the [Pediatric Open Targets, OpenPedCan-analysis](https://github.com/PediatricOpenTargets/OpenPedCan-analysis) raw DNA methylation array datasets. 

## Data description 

#### TARGET
The TARGET Illumina `Infinium HumanMethylation 450k BeadChip`  methylation arrays and sample metadata for the following pediatric cancers were downloaded from the [TARGET project website](https://ocg.cancer.gov/programs/target/target-methods):
 
- [Neuroblastoma (NBL)](https://target-data.nci.nih.gov/Public/NBL/methylation_array/)
- [Osteosarcoma (OS)](https://target-data.nci.nih.gov/Public/OS/methylation_array/)
- [Clear Cell Sarcoma of the Kidney (CCSK)](https://target-data.nci.nih.gov/Public/CCSK/methylation_array/)
- [Wilms Tumor (WT)](https://target-data.nci.nih.gov/Public/WT/methylation_array/)

###  Children's Brain Tumor Network (CBTN)
The [Children's Brain Tumor Network (CBTN)](https://cbtn.org/) `Infinium HumanMethylation EPIC (850k) BeadChip` methylation arrays and sample metadata for several Pediatric brain tumors.

## Analysis

- The TARGET Illumina methylation analysis results available on the [TARGET project website](https://ocg.cancer.gov/programs/target/target-methods) were preprocessed with different methylation software packages, including `minfi` (AML), `BeadStudio` (CCSK and WT), `methylumi` (NBL), and `Lumi+BMIQ` (OS). We have preprocessed all the TARGET and CBTN cancer types with Illumina arrays and using the updated version of [minfi Bioconductor package](https://academic.oup.com/bioinformatics/article/33/4/558/2666344). We utilized and `preprocessFunnorm` preprocessing method when an array dataset has control samples (i.e., normal and tumor samples) or multiple OpenPedcan cancer groups and `preprocessQuantile` when an array dataset has only tumor samples from a single OpenPedcan cancer group to estimate usable methylation measurements (Beta and M values) and copy number (cn-values) for OpenPedCan.

- Masking is applied using two different methods. Firstly, probes with known SNPs are removed. Secondly, the p-value for each sample is computed by comparing the signal of methylated and unmethylated probes against the background signal. Probes with p-values > 0.01 are filtered out.

- The `-methyl-cn-values.parquet` output is a probe-level copy-number signal intended for quality control only. It should not be used for CNV segmentation or as input to CONUMEE or GISTIC.

- In order to do funnorm normalization, a set of control probes must be correctly identified, and Median Absolute Deviation must be above 0 across these control probes. We had been seeing some errors for EPICv2 probes where MAD = 0, causing funnorm to fail, and have added an inspection and filtering step to check for MAD > 0 on the control probes, and skip and print out samples where MAD = 0. This array-level QC filtering step is applied before either normalization method.

## Running the analysis

### Recommended: shell wrapper

`run-preprocess-illumina-arrays.sh` validates its inputs, sorts mixed array types, preprocesses each detected supported array type, and intersects probe sets when two or more array types are present. Supply all run-specific paths on the command line; relative paths are resolved from the directory in which the wrapper is invoked.

```
bash run-preprocess-illumina-arrays.sh \
  --manifest_file controls_and_dicer_manifest.tsv \
  --input_dir input-test \
  --output_dir test-out \
  --output_prefix test
```

The sorted-IDAT staging directory is `<output_dir>/<output_prefix>-sorted-idats_output_dir`. Per-array output files and merged outputs use `<output_dir>/<output_prefix>` as their prefix.

### Running individual R scripts

`scripts/00-unzip-and-sort.R` separates supported array types before preprocessing. Run `scripts/01-preprocess-illumina-arrays.R` once for each resulting array directory:

```
Rscript --vanilla scripts/01-preprocess-illumina-arrays.R \
  --base_dir <sorted_idats_dir>/IlluminaHumanMethylationEPICv2 \
  --manifest_file <manifest.tsv> \
  --output_basename <output_dir>/<output_prefix> \
  --funnorm TRUE \
  --snp_filter TRUE \
  --n_cores 4
```

`--funnorm` and `--snp_filter` accept `TRUE` or `FALSE`; both default to `TRUE`. The script writes a `*-zero-mad.txt` file only when samples are excluded for zero median absolute deviation in a control-probe channel.

After preprocessing two or more array types, intersect their probe sets and write the combined matrices with:

```
Rscript --vanilla scripts/02-merge-methyl-matrices.R \
  --output_dir <output_dir> \
  --output_prefix <output_prefix>
```

## Input datasets

#### `Methylation arrays:`
Methylation array datasets are avaliable on the CHOP HPC `Isilon` sever (location: `/mnt/isilon/opentargets/wafulae/methylation-preprocessing/data/`). Please contact `Avin Farrel (@afarrel)` for access. 
- `data/NBL/*.idat` - Neuroblastoma (NBL) TARGET tumor sample arrays
- `data/OS/*.idat` - Osteosarcoma (OS) TARGET tumor sample arrays
- `data/CCSK/*.idat` - Clear Cell Sarcoma of the Kidney (CCSK) TARGET tumor sample arrays
- `data/WT/*.idat` - Wilms Tumor (WT) TARGET tumor sample arrays
- `data/AML/*.idat` - Acute Myeloid Leukemia (AML) TARGET tumor sample arrays
- `data/CBTN/*.idat` - Mutilple CBTN pediatric brian tumor sample arrays

## Results

Per-array results are written as parquet files using this pattern:

```
<output_prefix>-<array_type>-methyl-<measurement>.parquet
```

`<measurement>` is one of `beta-values-masked`, `m-values-unmasked`, `m-values-masked`, `cn-values`, or `p-values`. The supported array-type names are `IlluminaHumanMethylation450k`, `IlluminaHumanMethylationEPIC`, and `IlluminaHumanMethylationEPICv2`.

When multiple array types are available, the merge step writes intersection matrices for beta values, masked and unmasked M values, and CN values. The filename records the participating array types, for example `<output_prefix>-IlluminaHumanMethylationEPICv1-EPICv2-methyl-beta-values-masked.parquet`. Detection p-values are used to select among duplicated EPICv2 probes and are not merged. With one array type, the merge script exits successfully after reporting that no intersection is needed and leaves the per-array outputs unchanged.


## High Performance Computing (HPC)

In addition to the base scripts, this module contains the necessary components
to perform the more resource intensive steps remotely using Docker, Common
Workflow Language (CWL), and CAVATICA.

### Docker

This module contains a Dockerfile. That Dockerfile is used for running pieces
of this module in a CWL environment. The image generated by the Dockerfile is
stored on DockerHub, and it is provided here for reference. If you wish to
build the Dockerfile yourself, use the following command substituting in your
desired `image_name` and `image_version`:
`docker build -t <image_name>:<image_version> .`

### Common Workflow Language (CWL)

The CWL definitions are in `tools/` and `workflow/`. `workflow/methylation-preprocessing.cwl` runs sorting and scatters preprocessing over detected array types. Set the optional `merge_array_types` input to `true` to then run `merge_methyl_matrices.cwl` and return merged intersection matrices; it defaults to `false`. The workflow always returns the per-array parquet files and detection p-values.

### CAVATICA

To use the above tools or workflow(s), the user must upload them to CAVATICA.
To do that, the user must download scripts and publish an application.
Following command can be used to publish an application on CAVATICA:
`sbpack <sbg profile>  <user>/<projectname>/<workflowname> <path/to/workflow.cwl>`
Refer to this link for instructions on setting up [sbpack](https://docs.cavatica.org/docs/maintaining-and-versioning-cwl-on-external-tool-repositories).
