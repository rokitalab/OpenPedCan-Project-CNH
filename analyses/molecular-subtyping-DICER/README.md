## Retrieving DICER-1 Mutated Tumors from PBTA Cohort

### Purpose

To identify any DICER-1 Sarcoma's in the PBTA Cohort.

### Methods

First, the histologies file was filtered for the PBTA cohort. Secondly, since DICER-1 mutations are present in a variety of other tumors, samples with high confidence (greater than 0.8) methylation classifications other than DICER-1 are excluded. `snv-consensus-plus-hotspots.maf.tsv.gz` is read in chunks and filtered for SNVs in the DICER-1 gene, then further filtered for Misense, nonsense, and frameshift mutations.

### Analysis scripts

### `01-filter-PBTA-DICERS.R`

### To run this analysis: 





Usage:
```bash
bash run-data-pre-release-qc.sh

```