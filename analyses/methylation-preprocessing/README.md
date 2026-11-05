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

- Tumor purity is estimated for each sample with [RFpurify](https://github.com/mwsill/RFpurify), which predicts the ABSOLUTE and ESTIMATE purity of a sample from two random forests trained on TCGA 450k array data. Both forests are predicted from raw (`preprocessRaw`) beta values rather than the normalized matrices, because that is how the forests were trained, and because `predict.randomForest()` aborts on any `NA` and so cannot take detection p-value-masked or SNP-filtered values. The `-methyl-purity.parquet` output reports both estimates per sample, alongside the number of model probes that had to be imputed. Each forest uses 856 probes and the two share 38, so 1674 probes are read in total; all 1674 are present on both 450k and EPICv1, and nothing is imputed for those arrays in practice. Any that are missing are filled with 0.5 and counted in `n_probes_imputed`. EPICv2 carries only about 88% of each forest's probes, so purity is not estimated for EPICv2 arrays and no purity file is written for them.

- Each array carries a panel of `rs` probes that genotype common SNPs rather than measuring methylation, and these are written to `-methyl-snps.parquet` for checking sample identity against genomic data. They are taken from the raw `RGChannelSet` before normalization, because they have no entry in the array annotation and so are absent from the `GenomicRatioSet` the other measurements are derived from; this is unrelated to the `--snp_filter` option, which drops CpG probes that overlap a SNP. The panel is 65 probes on 450k, 59 on EPICv1 and 65 on EPICv2, of which 57 are shared across all three. Beta values are saved rather than genotype calls, so that the calling thresholds stay with the analysis making the comparison. EPICv2 probe names keep their replicate suffix, so the shared panel is found by stripping it. Like the per-sample QC and purity files, this output is not merged across array types.

- Detection p-values for those `rs` probes are written separately, to `-methyl-snp-p-values.parquet`. `minfi::detectionP()` iterates over the Type I and Type II probe tables only and so returns nothing for the `rs` probes, but the background model it uses applies to them unchanged, since the `Snp` probes share the structure of the ordinary ones. The calculation is reproduced here for the `Snp` addresses, and agrees with `detectionP()` to floating-point precision when run against the Type I and Type II probes instead. The betas themselves are left unmasked, so that the masking threshold stays with the analysis making the identity comparison.

- Detection p-values are used to mask individual failed probes (setting them to `NA`), not to exclude entire failed samples. To support sample-level QC downstream, the `-methyl-sample-qc.parquet` output reports, for each sample, the fraction of probes with detection p-value > 0.01. Whole-sample exclusion decisions based on this metric are left to the consuming analysis, rather than applied in this harmonization step.

## Running the analysis

### Recommended: shell wrapper

`run-preprocess-illumina-arrays.sh` validates its inputs, sorts mixed array types, preprocesses each detected supported array type, and intersects probe sets when two or more array types are present. Supply all run-specific paths on the command line; relative paths are resolved from the directory in which the wrapper is invoked.

```
bash run-preprocess-illumina-arrays.sh \
  --manifest_file manifest.tsv \
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

`--input_prefix` takes a comma-separated list of prefixes to read from, and defaults to `--output_prefix`. Give it more than one prefix to combine separate preprocessing runs of the same arrays, which is useful when samples arrive in batches: probes are intersected across the batches as well as across array types, and every batch's samples are carried into the merged matrices. A sample name occurring in more than one batch is an error rather than something the merge resolves silently.

Genotypes for the identity SNP probes are called from the `-methyl-snps.parquet` files with:

```
Rscript --vanilla scripts/03-call-snp-genotypes.R \
  --output_dir <output_dir> \
  --output_prefix <output_prefix> \
  --vcf TRUE
```

`--input_prefix` behaves as it does for the merge step. `--het_min` and `--het_max` set the beta cut points either side of a heterozygous call, and default to 0.25 and 0.75. Widening that band is worth considering for tumor samples, where loss of heterozygosity and aneuploidy pull heterozygous betas away from 0.5 and can push a true heterozygote into a homozygous call. `--vcf` additionally writes one VCF per array type, and defaults to `FALSE`.

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

`<measurement>` is one of `beta-values-masked`, `m-values-unmasked`, `m-values-masked`, `cn-values`, `p-values`, `sample-qc`, `purity`, `snps`, or `snp-p-values`. The supported array-type names are `IlluminaHumanMethylation450k`, `IlluminaHumanMethylationEPIC`, and `IlluminaHumanMethylationEPICv2`.

When multiple array types are available, the merge step writes intersection matrices for beta values, masked and unmasked M values, and CN values. The filename records the participating array types, for example `<output_prefix>-IlluminaHumanMethylationEPICv1-EPICv2-methyl-beta-values-masked.parquet`. With all three present, both the EPICv1+EPICv2 and the three-array intersections are written. Merging several batches of a single array type instead yields `<output_prefix>-<array_type>-merged-methyl-<measurement>.parquet`; the `-merged` suffix keeps that output from being picked up as a per-array input on a rerun.

Detection p-values are used to select among duplicated EPICv2 probes and, like the per-sample QC and purity files, are not merged. EPICv2 carries several probes for some CpG sites, distinguished by a `_` suffix; for a cross-array merge the probe with the lowest mean detection p-value across all EPICv2 samples is kept, ties broken on probe name, and the chosen probe is recorded in `<output_prefix>-EPICv2-unique-mapping.qs2`. With a single array type and a single batch, the merge script exits successfully after reporting that no intersection is needed and leaves the per-array outputs unchanged.

`03-call-snp-genotypes.R` writes `<output_prefix>-methyl-snp-genotypes.tsv`, with one row per probe and sample and the columns `probe_id`, `rs_id`, `sample`, `array_type`, `beta` and `gt`. `rs_id` is `probe_id` with the EPICv2 replicate suffix removed, and is what joins the panel across array types. With `--vcf`, a VCF is written per array type rather than as one combined file, because the array types interrogate different probe sets and a combined file would have to invent missing calls for every probe absent from an array. Detection p-values are not applied to these calls: `minfi::detectionP()` covers only the Type I and Type II methylation probes, so the rs probes have no detection p-value to mask against.

The `0/0` and `1/1` labels are arbitrary with respect to a reference genome. The array measures allele dosage but cannot say which homozygote carries the reference allele, so the VCF carries a placeholder `CHROM` and `POS` and an unknown `REF` and `ALT`, and is usable only by matching on the `ID` column; it will not support position-based or allele-based operations. Resolving that orientation per probe, and mapping the rs IDs to genomic coordinates, belongs to the analysis making the comparison against genomic data.

## High Performance Computing (HPC)

The preprocessing steps are memory intensive and are usually run remotely
rather than locally.

### Docker

These scripts run in the image built from the Dockerfile at the root of this
repository, which is stored on DockerHub. If you wish to build it yourself,
run the following command from the repository root, substituting in your
desired `image_name` and `image_version`:
`docker build -t <image_name>:<image_version> .`
