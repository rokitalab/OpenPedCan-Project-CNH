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
| [`generate-mtp-tables.sh`](generate-mtp-tables.sh) | Runs the fusion, SNV, CNV and TPM frequency modules and copies the Open Targets (MTP) tables to the release folder | Needs the AWS CLI to upload |
| [`make-parquet.R`](make-parquet.R) | Writes parquet copies of the large release tables to `data/<release>/parquet/` | See below |
| [`run_in_ci.sh`](run_in_ci.sh) | Builds the Docker image and runs a command inside it, as CI does | Passes the `OPENPBTA_*` environment variables |
| [`check-python.sh`](check-python.sh) | Checks that the Python packages in the image match `requirements.txt` | Run from CI |

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

How the tables are built:

- Column types are detected over the whole file (`sample_size = -1`), not from the first rows, and `NA` and empty strings are null.
- Small row groups with a sort on gene (or position) mean a gene lookup reads only the matching row groups.
- The only value change is `status` in the consensus copy number table, lowercased (the v15 Control-FREEC tumor-only input had `Gain` / `Loss` / `Neutral`; fixed upstream in `analyses/focal-cn-file-preparation/04-prepare-cn-file.R`). The script marks this in a comment next to the table definition.
- Row counts are compared with the source files and reported in `manifest.tsv`. Check that `row_count_matches` is `TRUE` for every table before uploading.

To run only some tables, set `PARQUET_TABLES` to a comma-separated list of table names; the manifest keeps the other rows. The conversion needs `duckdb` (and `dbplyr` 2.6.0 with `duckdb` 1.5.6, as in the Dockerfile) and about 10 GB of memory; larger-than-memory sorts spill to `scratch/parquet-duckdb-tmp/`.

The parquet files are not listed in a release `md5sum.txt`, so `download-data.sh` does not download them; upload the folder to `<release>/parquet/` next to the release files. The TSV / RDS files stay the source of truth.

Not converted: histologies, methylation beta values, TCGA and GTEx expression, rMATS splice events and isoform expression (to add for v16).
