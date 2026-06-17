## Integrate molecular subtyping output from pathology feedback

**Author of code and documentation:** [@kgaonkar6](https://github.com/kgaonkar6)

In this repo, we add molecular subtype for molecular_subtype from all subtyping modules and integrated_diagnosis, short_histology, broad_histology, and Notes from `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv`. A column `cancer_group` is added to provide broader terms derived from `harmonized_diagnosis` which can be used to generate figures.

### Usage
```sh
bash run-subtyping-integrate.sh
```

### Module contents

`01-integrate-subtyping.Rmd` integrates results from compiled results in `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv` to `histologies-base.tsv`

### Table of cancers and resulting cancer groups output in results
```bash
.
├── 01-integrate-subtyping.Rmd
├── 01-integrate-subtyping.nb.html
├── README.md
├── input
│   └── column_order.txt
├── results
│   ├── cancer_group_table_for_README.tsv
│   ├── discrepancies_to_check.tsv
│   ├── harmonized_diagnosis_cancer_group_table.tsv
│   ├── histologies.tsv
│   ├── pediatric_cancer_groups.tsv
│   ├── subtype_tumor_counts_table.tsv
│   └── tcga_cancer_groups.tsv
└── run-subtyping-integrate.sh
```