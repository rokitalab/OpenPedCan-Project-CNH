# Molecular Subtyping of Rare CNS

**Module author:** Aylar Babaei, Jo Lynne Rokita

## Overview

This directory contains scripts and outputs for the Rare CNS molecular subtyping workflow

The goal of this analysis is to identify Rare CNS tumors for downstream molecular subtyping. 

Related issue: [#114](https://github.com/rokitalab/OpenPedCan-Project-CNH/issues/114)

## Workflow steps

The Rare CNS subtyping workflow is organized into the following scripts:

`00-rare-cns-select-tumors.R` identifies candidate Rare CNS samples based on methylation subtype criteria and generates the JSON file for selection in downstream subtyping steps.

`01-subset-files-for-Rare-CNS.Rmd` subsets Rare CNS DNA, RNA, and methylation biospecimens and generates the initial methylation subtype subset file for downstream analyses.

### Outputs

- `rare-cns-subset/rare_cns_subtyping_path_dx_strings.json`
- `rare-cns-subset/RareCNS_dna_biospecimen.tsv`
- `rare-cns-subset/RareCNS_rna_biospecimen.tsv`
- `rare-cns-subset/RareCNS_methyl_biospecimen.tsv`
- `rare-cns-subset/RareCNS_methyl_subtypes.tsv`


### Usage

From within this directory:
```
bash run-molecular-subtyping-rare-cns.sh
```
