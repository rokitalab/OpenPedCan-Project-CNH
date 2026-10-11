# Molecular subtyping of CNS tumors using methylation and genetic lesions

**Module authors:** Jo Lynne Rokita

## Overview

This module builds one resolved methylation call per biospecimen, which molecular subtyping modules consume instead of each re-implementing score thresholds and classifier-version aliasing.

## Workflow

`01-build-methylation-calls.R` combines:

- DKFZ v12.8 / NIH v3.1 calls from `input/bti_methylation_classifications.tsv`
- DKFZ v12 / NIH v2 calls from `data/histologies-base.tsv`

Each call is mapped to a single row of `input/mnp_v12.5_v12.8_NIHv2_v3.1_annotation_with_OPC_subtypes.xlsx` through normalized keys (case and `_` vs. space), so classifier-version renames (e.g., DKFZ `ANTCON` to `GTAKA`) resolve without hand-maintained aliases.

### Confidence rule

A call is high confidence if the DKFZ subclass score is >= 0.8, or the NIH superfamily and class scores are both >= 0.9. This is the same rule used in `molecular-subtyping-Rare-CNS`. Precedence among passing calls is DKFZ v12.8, NIH v3.1, DKFZ v12, NIH v2; only calls that resolve to a crosswalk row can set the subtype.

### Subtype additions

`input/opc_subtype_additions.tsv` names subtypes for classes whose crosswalk `OPC_molecular_subtype` was blank (choroid plexus, nerve sheath, germ cell, pineal, sellar, and other rare classes). It is keyed on the NIH v2 abbreviation, never overwrites an existing crosswalk value, and carries a rationale, the genetic lesion expected for each class (`expected_lesion`, used later to confirm methylation calls), and a `non_cns` flag for classes that are not primary CNS tumors (`MET_MEL`, `ENB_IDH`). These are kept so every classifier call resolves, and surface as `methyl_non_cns = TRUE`; none are present in the current cohort.

### Control classes

Control tissues and the inflammatory microenvironment class (`INFLAM_ENV`) indicate low tumor purity and are never subtyped. The crosswalk marks these with `Control` in `OPC_molecular_subtype`, which is distinct from a blank cell (a class not yet assigned a subtype). Any NIH control code (`CONTR_*`, `CTRL_*`, `NORM_*`) is also treated as a control. These biospecimens get `methyl_control = TRUE` and no `OPC_molecular_subtype`, so they can be flagged downstream without being given a tumor subtype.

### Outputs

- `results/methylation_calls.tsv` - one row per methylation biospecimen: `methyl_confidence`, `methyl_source`, `methyl_class`, `methyl_score`, `OPC_molecular_subtype`, `rare`, family fields, `methyl_control`, `methyl_non_cns`, `methyl_conflict` (another high-confidence call from a different version resolves to a different class), and `unmapped_high_conf_call`
- `results/methylation_crosswalk_gaps.tsv` - high-confidence classes that are not in the crosswalk, or whose crosswalk row has no `OPC_molecular_subtype`, with biospecimen counts; this is the list of subtypes still needing a rule

### Validation

`02-prepare-lesion-tables.R` scans the SNV consensus + hotspot MAF once (about 30 seconds) and writes only the genes named in `input/subtype_lesion_rules.tsv` to a parquet file in `scratch/`.

`03-validate-subtypes.R` then checks each methylation subtype two ways, without ever overriding the methylation call:

- **Lesions:** rules in `input/subtype_lesion_rules.tsv` (SNV pattern on `HGVSp_Short`, or fusion pattern on `FusionName`). SNVs and fusions are linked to the methylation biospecimen through `match_id`. `lesion_status` is `confirmed`, `not_detected` (assay present, lesion absent), `no_assay`, or `no_rule`.
- **Diagnosis:** `pathology_diagnosis` is mapped to a group (`input/pathology_dx_groups.tsv`) and compared with the groups allowed for the subtype's family (`input/subtype_family_dx_compatibility.tsv`, first allowed group is primary). `dx_concordance` is `concordant`, `compatible`, `discordant` (probable reclassification), or `no_informative_dx`. Broad histology is not used, since it is itself a manual lookup of pathology diagnosis.

Tumors with no high-confidence methylation call (PBTA only, controls excluded) are subtyped from lesions alone for rules with `use_as_fallback = TRUE`.

Outputs: `results/methylation_subtype_validation.tsv`, `results/subtypes_for_review.tsv` (discordant diagnosis or conflicting classifier versions), and `results/lesion_fallback_subtypes.tsv`.

### Usage

```
Rscript --vanilla 01-build-methylation-calls.R
Rscript --vanilla 02-prepare-lesion-tables.R
Rscript --vanilla 03-validate-subtypes.R
```
