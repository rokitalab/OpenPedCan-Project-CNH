# Molecular Subtyping of Rare CNS

**Module author:** Aylar Babaei, Jo Lynne Rokita

## Overview

This directory contains scripts and outputs for the Rare CNS molecular subtyping workflow

The goal of this analysis is to identify Rare CNS tumors for downstream molecular subtyping. 

Related issue: [#114](https://github.com/rokitalab/OpenPedCan-Project-CNH/issues/114)

## Workflow steps

The Rare CNS subtyping workflow is organized into the following scripts:

`00-rare-cns-select-tumors.R` identifies candidate Rare CNS samples based on methylation subtype criteria and generates the JSON file for selection in downstream subtyping steps.
`01-subtype-using-fusions.Rmd` identifies hallmark-fusion-based candidate subtypes from RNA-Seq biospecimens; these are not added to the cohort on their own.
`02-subtype-using-methylation.Rmd` subtypes Rare CNS tumors using high-confidence methylation (the sole source of the final `molecular_subtype`), layering in a newer DKFZ v12.8 / NIH v3.1 classifier run on top of the DKFZ v12 / NIH v2 calls in `histologies-base.tsv` (preferring v12.8 when both are available, and flagging any disagreement between the two for manual review), then confirms that any hallmark fusion detected for a tumor agrees with its methylation-based subtype, flagging disagreements for manual review. It also hard-stops if any `Kids_First_Biospecimen_ID` ends up with more than one `molecular_subtype`.

Tumors are subtyped using DKFZ if the subclass score is >=0.8; if the DKFZ score is <0.8 (or missing), NIH v2 is used as a fallback if both the superfamily mean and class mean scores are >0.9. The same 0.8 / 0.9 thresholds are applied to the newer DKFZ v12.8 / NIH v3.1 calls, read from `../molecular-subtyping-methylation/input/bti_methylation_classifications.tsv` (an input shared with the `molecular-subtyping-methylation` module, pre-filtered to `BS_`-prefixed `Bioassay_ID`s, i.e. those already a `Kids_First_Biospecimen_ID` in this cohort). Since some DKFZ/NIH subclass abbreviations changed between classifier versions (e.g. DKFZ's `ANTCON` became `GTAKA` in v12.8), `02-subtype-using-methylation.Rmd` aliases the known renames back to the v12/v2 vocabulary before matching; codes with no confirmed alias are simply not matched rather than guessed.

### Outputs
- `rare-cns-subset/rare_cns_subtyping_path_dx_strings.json` - methylation classifications to consider
- `results/rare_cns_fusion_hits.tsv` - RNA-Seq biospecimens with a hallmark fusion, for inspection
- `results/rare_cns_fusion_based_subtypes.tsv` - hallmark-fusion-based candidate subtype per biospecimen
- `results/rare-cns-molecular-subtypes.tsv` - final subtypes (methylation-based, DKFZ v12/v12.8 and NIH v2/v3.1 combined)
- `results/rare_cns_fusion_methylation_confirmation.tsv` - for methylation-classified tumors with RNA-Seq, whether the fusion-based call (if any) agrees with the methylation-based subtype, with `pathology_diagnosis` included for manual review


### Directory structure
```
.
├── 00-rare-cns-select-tumors.R
├── 01-subtype-using-fusions.Rmd
├── 01-subtype-using-fusions.nb.html
├── 02-subtype-using-methylation.Rmd
├── 02-subtype-using-methylation.nb.html
├── README.md
├── rare-cns-subset
│   └── rare_cns_subtyping_path_dx_strings.json
├── results
│   ├── rare-cns-molecular-subtypes.tsv
│   ├── rare_cns_fusion_based_subtypes.tsv
│   ├── rare_cns_fusion_hits.tsv
│   └── rare_cns_fusion_methylation_confirmation.tsv
└── run-molecular-subtyping-rare-cns.sh
```

### Usage

From within this directory:
```
bash run-molecular-subtyping-rare-cns.sh
```
