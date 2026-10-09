# Methylation CNV Calling

This module calls copy-number segments from QC-filtered Illumina methylation
array data. It reads the RGChannelSet written by `methylation-preprocessing`,
applies Noob normalization, and uses conumee2 to write a standard segment file
and a GISTIC-ready segment file.

The RGChannelSet is an R object, so it is read as `.qs2`. Parquet is used only
for the preprocessing matrix outputs (beta, M, p-values, QC), which are not
inputs to this module. The preprocessing `-methyl-cn-values.parquet` output is
for QC only and is not used for CNV segmentation.

## Inputs

- A QC-filtered RGset, `<prefix>-<IlluminaHumanMethylation array>-rg-set.qs2`,
  written by `methylation-preprocessing/scripts/01-preprocess-illumina-arrays.R`
  (one per array type).
- The same manifest given to preprocessing, with `file_name`, `Bioassay_ID`,
  `sample_type`, and `platform` columns.
- The matching array type: `EPIC`, `EPICv2`, or `450k`.

### Reference (control) samples

Rows in the manifest with `sample_type == Normal` are the conumee2 reference
set; every other sample in the RGset is a query. Rows are matched to the array
type through the `platform` column (spelling variants such as
`HumanMethylationEPICv2` and `Illumina Infinium HumanMethylationEPICv2` are
both accepted). The run stops with an error if an array has no Normal samples
or no non-Normal samples in the RGset.

## Method

1. Read the RGset and apply `minfi::preprocessNoob`.
2. Keep samples present in the manifest for the array type and rename them to
   `Bioassay_ID`.
3. Fit query samples against the pooled Normal reference with conumee2
   (EPICv2 uses hg38 annotation; EPIC and 450k use hg19).
4. Segment with `alpha = 0.001`, `nperm = 50000`, `min.width = 5`, and
   `undo.SD` of 2.5 (EPIC/EPICv2) or 2.2 (450k).

## Usage

```bash
bash run-methylation-cnv-calling.sh \
  --rg-set-file <output_dir>/<prefix>-IlluminaHumanMethylationEPICv2-rg-set.qs2 \
  --manifest-file <manifest.tsv> \
  --output-dir <cnv_output_dir> \
  --output-prefix <prefix>-EPICv2 \
  --array-type EPICv2
```

Run once per array type. The runner writes `<output-prefix>-segments.seg` and
`<output-prefix>-gistic.seg` to the output directory, so use a unique prefix
for each array type.

## Liftover and combine GISTIC segments

`EPICv1` and `450k` calls use hg19 annotations, while `EPICv2` calls use
hg38. After CNV calling has produced all per-array `*-gistic.seg` files, run:

```bash
bash run-liftover-methylation-segments.sh \
  --seg-dir test-output \
  --chain-file /path/to/hg19ToHg38.over.chain
```

The liftover is skipped for an EPICv2-only run. Otherwise the chain file is
required. The result is `combined_hg38.gistic.seg` in `--seg-dir`.

## CWL

`tools/methylation_cnv_calls.cwl` packages the same script for CWL runners.
Supply the RGset and manifest as `File` inputs and choose a unique
`output_basename` for each array type.

`tools/liftover_methylation_segments.cwl` accepts the resulting per-array
GISTIC SEG files plus the optional chain file and emits the combined hg38 SEG.

## Continuous integration

The `methylation_preprocessing` job in `.github/workflows/run_analysis.yml`
downloads the methylation testing set (`methyl_testing`, including
`methylation-testing-manifest.tsv`), runs preprocessing, and then runs this
module once each for EPICv2, EPIC, and 450k using the RGsets from that run. The
Normal rows of the testing manifest provide the controls for each array. The
job fails if any `-segments.seg` or `-gistic.seg` output is missing or empty.

Local runs need `conumee2` and `minfi`; the `openpedcanverse` Docker image
includes them.
