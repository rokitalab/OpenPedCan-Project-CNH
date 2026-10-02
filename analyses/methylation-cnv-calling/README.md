# Methylation CNV Calling

This module calls copy-number segments from QC-filtered Illumina methylation
RGsets produced by `methylation-preprocessing`. It applies Noob normalization
and then uses conumee2 to create a standard segment file and a GISTIC-ready
segment file.

## Inputs

- An RGset `.qs2` file produced by methylation preprocessing.
- A manifest with `file_name`, `Bioassay_ID`, `sample_type`, and `platform`
  columns. Samples marked `Normal` are used as the conumee2 reference set.
- The matching array type: `EPIC`, `EPICv2`, or `450k`.

## Local test run

RGset test inputs belong in `test-input/`; this directory and `test-output/`
are ignored by Git. Three local RGsets have been copied there for testing and
must not be committed.

For example:

```bash
bash run-methylation-cnv-calling.sh \
  --rg-set-file test-input/test-IlluminaHumanMethylationEPICv2-rg-set.qs2 \
  --manifest-file ../methylation-preprocessing/controls_and_dicer_manifest.tsv \
  --output-dir test-output \
  --output-prefix test-EPICv2 \
  --array-type EPICv2
```

The runner writes `<output-prefix>-segments.seg` and
`<output-prefix>-gistic.seg` to the output directory.

## CWL

`tools/methylation_cnv_calls.cwl` packages the same script for CWL runners.
Supply the RGset and manifest as `File` inputs and choose a unique
`output_basename` for each array type.
