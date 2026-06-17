# Molecular Subtyping of Rare CNS

**Module author:** Aylar Babaei, Jo Lynne Rokita

## Overview

This directory contains scripts and outputs for the Rare CNS molecular subtyping workflow

The goal of this analysis is to identify Rare CNS tumors for downstream molecular subtyping. 

Related issue: [#114](https://github.com/rokitalab/OpenPedCan-Project-CNH/issues/114)

## Workflow steps

The Rare CNS subtyping workflow is organized into the following scripts:

`00-rare-cns-select-tumors.R` identifies candidate Rare CNS samples based on methylation subtype criteria and generates the JSON file for selection in downstream subtyping steps.
`01-subtype-Rare-CNS.Rmd` uses methylation to subtype Rare CNS tumors

DKFZ v12 subclass scores must be >=0.8 and if exists and meets criteria, tumors are subtyped using DKFZ. 
If DKFZ v12 subclass scores <0.8, classify using NIH v2 if >0.9 for superfamily mean and class mean score.

### Outputs
- `rare-cns-subset/rare_cns_subtyping_path_dx_strings.json` - methylation classifications to consider
- `results/rare-cns-molecular-subtypes.tsv` - final subtypes


### Directory structure
```
.
├── 00-rare-cns-select-tumors.R
├── 01-subtype-Rare-CNS.Rmd
├── 01-subtype-Rare-CNS.nb.html
├── README.md
├── rare-cns-subset
│   └── rare_cns_subtyping_path_dx_strings.json
├── results
│   └── rare-cns-molecular-subtypes.tsv
└── run-molecular-subtyping-rare-cns.sh
```

### Usage

From within this directory:
```
bash run-molecular-subtyping-rare-cns.sh
```
