# Molecular Subtyping of Rare CNS

**Module author:** Aylar Babaei

## Overview

This directory contains scripts and outputs for the Rare CNS molecular subtyping workflow?

The goal of this analysis is to identify Rare CNS tumors for downstream molecular subtyping. 

Related issue: [#114](https://github.com/rokitalab/OpenPedCan-Project-CNH/issues/114)

## Workflow steps

The Rare CNS subtyping workflow is organized into the following scripts:

- `00-Rare-CNS-select-pathology-dx.R`
- `01-subset-files-for-Rare-CNS.Rmd`
- `02-subset-fusion-files-Rare-CNS.Rmd`
- `03-subset-cnv-files-Rare-CNS.Rmd`
- `04-Rare-CNS-compile-subtypes.Rmd`
- `05-Rare-CNS-methylation-umap.Rmd`

## Current PR scope

This PR includes only:

`01-subset-files-for-Rare-CNS.Rmd`
- `rare-cns-subset/RareCNS_dna_biospecimen.tsv`
- `rare-cns-subset/RareCNS_rna_biospecimen.tsv`
- `rare-cns-subset/RareCNS_methyl_biospecimen.tsv`
- `rare-cns-subset/RareCNS_methyl_subtypes.tsv`
- this README

Downstream steps will be added in stacked follow-up pull requests.

### Purpose

`00-Rare-CNS-select-pathology-dx.R` identifies candidate Rare CNS samples based on pathology diagnosis criteria and generates the initial files used in downstream subtyping steps.

`01-subset-files-for-Rare-CNS.Rmd` subsets Rare CNS DNA, RNA, and methylation biospecimens and generates the initial methylation subtype subset file for downstream analyses.

### Outputs

- `rare-cns-subset/rare_cns_metadata.tsv`
- `rare-cns-subset/rare_cns_subtyping_path_dx_strings.json`
- `rare-cns-subset/RareCNS_dna_biospecimen.tsv`
- `rare-cns-subset/RareCNS_rna_biospecimen.tsv`
- `rare-cns-subset/RareCNS_methyl_biospecimen.tsv`
- `rare-cns-subset/RareCNS_methyl_subtypes.tsv`


### Usage

From within this directory:

```sh
Rscript --vanilla 00-Rare-CNS-select-pathology-dx.R
Rscript --vanilla 01-subset-files-for-Rare-CNS.Rmd