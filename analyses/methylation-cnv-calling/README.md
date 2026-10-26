# Methylation CNV Calling

This module calls copy-number segments from QC-filtered Illumina methylation
array data. It reads the RGChannelSet written by `methylation-preprocessing`,
applies Noob normalization, and uses conumee2 to write a standard segment file
and a GISTIC-ready segment file. It can then lift the hg19 arrays to hg38,
combine all arrays, summarize the cohort that will go into GISTIC, and run
GISTIC2.

| Step | Script | Output |
|---|---|---|
| 1. CNV calling (per array) | `scripts/01-methylation-cnv-calls.R` | `<prefix>-segments.seg`, `<prefix>-gistic.seg` |
| 2. Liftover and combine | `scripts/02-liftover-methylation-segments.R` | `combined_hg38.gistic.seg` |
| 3. Cohort QC | `scripts/03-gistic-cohort-qc.R` | `gistic-cohort-qc.tsv` |
| 4. GISTIC | `scripts/04-run-gistic.sh` | `<output-dir>/gistic/` |

All steps are run through `run-methylation-cnv-calling.sh`. The scripts create
their output directories if they do not exist.

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

Each array is called against its own Normal pool, so the pool should match the
array type and sample preparation (for example FFPE vs fresh frozen), be copy
number neutral, and be large enough to give a stable estimate; more Normals give
diminishing returns, so the quality and match of the reference matter more than
the count.

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
hg38. After CNV calling has produced all per-array `*-gistic.seg` files, run the
same runner in combine-only mode:

```bash
bash run-methylation-cnv-calling.sh \
  --combine-only \
  --output-dir <cnv_output_dir> \
  --chain-file /path/to/hg19ToHg38.over.chain
```

Alternatively, add `--combine-segments` (and `--chain-file`) to the final
per-array calling command to combine immediately after that call.

The liftover is skipped for an EPICv2-only run. Otherwise the chain file is
required. The result is `combined_hg38.gistic.seg` in `--output-dir`.

## Cohort QC for GISTIC

Combining also runs `scripts/03-gistic-cohort-qc.R`, which writes
`gistic-cohort-qc.tsv` to `--output-dir`: one row per array with the number of
cohort samples and segments, median segments per sample, median segment width
and marker count, `Seg.CN` quantiles and SD, the fraction of segments and of the
genome with `|Seg.CN| >= 0.2` (the GISTIC `-ta`/`-td` value), and the fraction
of segments retained after liftover.

The cohort is the samples in `--bioassay-ids FILE` (the same list GISTIC is run
on), or all samples when it is not given, so the table also shows how many of
the cohort's samples are on each array. The ID list is used here whenever it is
given, with or without `--run-gistic`. Use one specimen per patient: GISTIC
assumes independent samples.

Check this table before trusting pooled GISTIC peaks. The arrays use different
annotations (hg19 vs hg38), `undo.SD` values (2.2 for 450k, 2.5 for EPIC and
EPICv2), and Normal reference pools, so their segment counts and amplitudes
need not match. Large differences in amplitude or altered fraction between
arrays, or a low retained fraction after liftover, mean the pooled segments
should be treated with caution.

## Run GISTIC

GISTIC2 runs on `combined_hg38.gistic.seg` from the same runner, so no
separate module is needed:

```bash
bash run-methylation-cnv-calling.sh \
  --combine-only \
  --output-dir <cnv_output_dir> \
  --chain-file /path/to/hg19ToHg38.over.chain \
  --run-gistic
```

`--run-gistic` implies combining, so it can also be added to the final
per-array calling command. Results go to `<output-dir>/gistic` unless
`--gistic-output-dir` is given.

`--bioassay-ids FILE` (one Bioassay_ID per line) defines the cohort and is
optional; without it GISTIC runs on every sample in the combined SEG. When set,
`scripts/04-run-gistic.sh` first subsets the SEG to those samples, because GISTIC
has no sample-subset option and its recurrence statistics depend on the cohort.
The filtered SEG (`combined_hg38.gistic.filtered.seg`) is written next to the
GISTIC results, and the run stops if none of the IDs are found.

`scripts/04-run-gistic.sh` expects GISTIC2 at `$GISTIC_DIR` (default
`/home/rstudio/gistic_install`, as installed by `Dockerfile.base`; the
`share/gistic2-2.0.23-0` subdirectory is also detected) and the hg38
reference-gene file at `$REFGENE_FILE`. CI runs it on the combined testing
segments with no ID list.

## CWL

`tools/methylation_cnv_calls.cwl` packages the CNV calling script for CWL
runners. Supply the RGset and manifest as `File` inputs and choose a unique
`output_basename` for each array type.

`tools/liftover_methylation_segments.cwl` accepts the resulting per-array
GISTIC SEG files plus the optional chain file and emits the combined hg38 SEG.
Both tools are intended for running on CAVATICA. There is no CWL tool for the
GISTIC step.

## Chain file

`hg19ToHg38.over.chain.gz` is downloaded with the methylation data by
`scripts/download-methyl.sh` (`data/hg19ToHg38.over.chain.gz`) and with the
testing data by `scripts/download-methyl-testing.sh`
(`data/methyl_testing/hg19ToHg38.over.chain.gz`). Each takes it from the
OpenPedCan bucket when it is listed in the md5sum file and otherwise from UCSC.
Unzip it (`gunzip -k`) before passing it to `--chain-file`.

## Continuous integration

The `methylation_preprocessing` job in `.github/workflows/run_analysis.yml`:

1. downloads the methylation testing set with
   `scripts/download-methyl-testing.sh` (`methyl_testing`, including
   `methylation-testing-manifest.tsv` and the chain file),
2. runs preprocessing,
3. runs this module once each for EPICv2, EPIC, and 450k using the RGsets from
   that run (the Normal rows of the testing manifest provide the controls for
   each array), and
4. unzips the chain file and runs liftover, combining, and cohort QC with
   `--combine-only` and the chain file, and
5. runs GISTIC2 on `combined_hg38.gistic.seg` with `scripts/04-run-gistic.sh`.

The job fails if any `-segments.seg` or `-gistic.seg` output, or
`combined_hg38.gistic.seg`, `gistic-cohort-qc.tsv`, or the GISTIC
`all_lesions.conf_99.txt` and `scores.gistic`, is missing or empty.

Local runs need `conumee2` and `minfi`; the `openpedcanverse` Docker image
includes them.
