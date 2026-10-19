## Scripts

Helper scripts for downloading data, generating release files, running the
subtyping modules and building the parquet copies of the large tables. Run them
from anywhere: the shell scripts change to their own directory, and
`make-parquet.R` finds the repository root.

| Script | Purpose | Notes |
|---|---|---|
| [`download-methyl.sh`](download-methyl.sh) | Downloads the DNA methylation files of a release (`<release>/methyl-md5sum.txt`) and the hg19 to hg38 liftover chain, checks the md5sums and links them in `data/` | `OPENPEDCAN_RELEASE` (default `v15`) and `OPENPEDCAN_URL`. The main data download is [`download-data.sh`](../download-data.sh) in the repository root |
| [`download-general-testing-files.sh`](download-general-testing-files.sh) | Downloads the small `testing` release used by CI, and the GENCODE and reference files | Replaces the symlinks in `data/` |
| [`download-methyl-testing.sh`](download-methyl-testing.sh) | Downloads the methylation test IDATs and manifest used by CI into `data/methyl_testing/` | `TESTING_DIR` |
| [`generate-analysis-files.sh`](generate-analysis-files.sh) | Runs the pre-release analysis modules (copy number consensus, independent samples, fusion filtering, fusion summary, TMB, GSEA, TP53 score) and collects their outputs | `RUN_LOCAL=1` skips the time-intensive steps. Step 3 of the data release process in the root README |
| [`run-for-subtyping.sh`](run-for-subtyping.sh) | Runs the molecular subtyping modules, then `molecular-subtyping-pathology` and `molecular-subtyping-integrate`, and the independent sample list | Step 4 of the data release process. Still runs the legacy per-tumor-type modules |
| [`make-parquet.R`](make-parquet.R) | Writes parquet copies of the large release tables to `data/<release>/parquet/` | See below |
| [`download-parquet.sh`](download-parquet.sh) | Downloads the parquet copies (all, or the tables named as arguments), checks `md5sum-parquet.txt` and links them in `data/parquet/` | `OPENPEDCAN_RELEASE` (default `v15`) and `OPENPEDCAN_URL`. The files are not in `md5sum.txt`, so `download-data.sh` does not fetch them |
| [`run_in_ci.sh`](run_in_ci.sh) | Builds the Docker image and runs a command inside it, as CI does | Passes the `OPENPBTA_*` environment variables |
| [`check-python.sh`](check-python.sh) | Checks that the Python packages in the image match `requirements.txt` | Run from CI |

### Archived scripts

[`archive/`](archive/) holds scripts that are no longer part of the release process and are kept for reference:

- [`archive/generate-mtp-tables.sh`](archive/generate-mtp-tables.sh) ran the fusion, SNV, CNV and TPM frequency modules and copied the Open Targets (MTP) tables to the release folder (last used for v12; the paths were adjusted for its new location).

### Parquet copies of the release tables

```sh
Rscript scripts/make-parquet.R v15
```

Reads the TSV / RDS files in `data/v15/` and writes, to `data/v15/parquet/`:

- one `<table>.parquet` per table,
- `manifest.tsv` (table, source file, rows, source rows, columns, size, md5, and whether the row count matches the source),
- `md5sum-parquet.txt` (md5 and file name, in the `md5sum.txt` format).

| Table | Source | Sorted by |
|---|---|---|
| `snv-consensus-plus-hotspots`, `snv-mutect2-tumor-only-plus-hotspots` | `*.maf.tsv.gz` | `Hugo_Symbol`, `Chromosome`, `Start_Position` |
| `consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only` | `*.tsv.gz` | `gene_symbol`, `biospecimen_id` |
| `consensus_seg_with_status`, `cnvkit_with_status` | `*.tsv` | biospecimen, chromosome, start |
| `fusion-annoFuse`, `fusion-arriba`, `fusion-starfusion`, `fusion-dgd`, `fusion-putative-oncogenic` | fusion files | fusion name or sample |
| `sv-manta` | `sv-manta.tsv.gz` | chromosome, start |
| `gene-expression-rsem-tpm-collapsed`, `gene-counts-rsem-expected_count-collapsed` | `*.rds` | `gene_symbol` (a `gene_symbol` column and one column per biospecimen) |
| `methyl-beta-values`, `methyl-m-values`, `methyl-cn-values` | `methyl-*.rds` | not sorted (probe order of the RDS; a `Probe_ID` column and one column per array, row groups of 10,000 probes) |

How the tables are built:

- Column types are detected over the whole file (`sample_size = -1`), not from the first rows, and `NA` and empty strings are null.
- Small row groups with a sort on gene (or position) mean a gene lookup reads only the matching row groups.
- The only value change is `status` in the consensus copy number table, lowercased (the v15 Control-FREEC tumor-only input had `Gain` / `Loss` / `Neutral`; fixed upstream in `analyses/focal-cn-file-preparation/04-prepare-cn-file.R`). The script marks this in a comment next to the table definition.
- Row counts are compared with the source files and reported in `manifest.tsv`. Check that `row_count_matches` is `TRUE` for every table before uploading.

The methylation matrices are read whole by `readRDS()` (tens of GB in memory), so convert them on a large-memory machine (`PARQUET_TABLES=methyl-beta-values`); the script compares 100 random values of each with the RDS and stops on a mismatch. In v15 the matrices stack the arrays of two array types (2,656 and 126 arrays in the beta matrix, with NA for the other type), so a probe has two rows and the parquet table has 1,702,675 rows for 867,251 probes. A query for some arrays then reads only their columns. Code that reads `methyl-beta-values.parquet` when it is in `data/parquet/` (and otherwise the RDS): `analyses/molecular-subtyping/09-plot-mb-methylation-umap.R`, `analyses/methylation-summary/02-calculate-methly-quantiles.R` (probes and independent-sample arrays only) and `analyses/data-pre-release-qc/01-data-harmonization-qc.Rmd` (probe and array IDs only). `analyses/methylation-summary/03-methyl-tpm-correlation.py` still reads the RDS.

To run only some tables, set `PARQUET_TABLES` to a comma-separated list of table names; the manifest keeps the other rows. The conversion needs `duckdb` (and `dbplyr` 2.6.0 with `duckdb` 1.5.6, as in the Dockerfile) and about 10 GB of memory; larger-than-memory sorts spill to `scratch/parquet-duckdb-tmp/`.

The parquet files are not listed in a release `md5sum.txt`, so `download-data.sh` does not download them; upload the folder to `<release>/parquet/` next to the release files, then users fetch them with `bash scripts/download-parquet.sh` (md5 check, links in `data/parquet/`). The TSV / RDS files stay the source of truth.

Not converted: histologies, methylation beta values, TCGA and GTEx expression, rMATS splice events and isoform expression (to add for v16).

## Testing copy of the beta values for CI

`scripts/make-methyl-testing-parquet.R [source release] [testing release]` writes `data/testing/parquet/methyl-beta-values.parquet` (with `manifest.tsv` and `md5sum-parquet.txt`) from the release parquet: the arrays that the molecular subtyping investigates (all medulloblastoma arrays, all IDH-mutant astrocytoma arrays, a sample of each non-IDH glioma class and the tumors reviewed for G-CIMP) and the rows it uses (the 1,503 G-CIMP probes and the 20,000 most variable probes of the medulloblastoma arrays). Upload the folder to the `testing` release; the CI job sets `METHYL_PARQUET_RELEASE=testing`, so `run-molecular-subtyping.sh` (which downloads the beta values when the file is missing) fetches this small file instead of the 15 GB release file, while a CI step downloads the other parquet tables the module reads (copy number, fusion, SV, SNV, expression; about 1.9 GB) from the v15 release.

