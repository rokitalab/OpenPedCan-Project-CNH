# Molecular subtyping of CNS tumors using methylation and genetic lesions

**Module authors:** Jo Lynne Rokita

## Overview

This module builds one resolved methylation call per biospecimen, which molecular subtyping modules consume instead of each re-implementing score thresholds and classifier-version aliasing.

## Workflow

Scripts run in order: `01` builds the methylation calls, `02` prepares lesion tables, `03` validates and applies the lesion fallback, `05` compiles the final table, and `04` compares with the legacy modules. Shared functions are in `utils/`.

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

IDH-mutant gliomas are resolved with 1p/19q status (both arms <= -0.5 in the consensus GISTIC broad arm values, WGS only): IDH plus codeletion is `Oligo, IDH`. Without codeletion the grade comes from the diagnosis group (`LGG, IDH` or `HGG, IDH`), and a tumor with an oligodendroglioma diagnosis but no codeletion evidence is held unclassified (`hold_reason` says whether the arm was not codeleted or not assessable).

TP53 annotations come from `analyses/tp53_nf1_score` and are not part of the subtype:

- `tp53_altered` (`loss`, `activated`, `Other`) and `tp53_classifier_score`
- `tp53_n_hits` and `tp53_two_hits`: hits are counted as in that module, with each SNV/indel as one hit and a CNV loss or SV as one more (CNV and SV are the same event); two or more suggests both alleles are affected
- `tp53_hotspot`, `tp53_li_fraumeni` (a germline second hit is inferred from Li-Fraumeni predisposition in the module), and the SNV and CNV evidence strings

Outputs: `results/methylation_subtype_validation.tsv`, `results/subtypes_for_review.tsv` (discordant diagnosis or conflicting classifier versions), and `results/lesion_fallback_subtypes.tsv`.

`05-compile-subtypes.R` writes `results/molecular_subtype_cns.tsv`, one row per `match_id`. Subtype precedence is high-confidence methylation, then lesion-only fallback, then the LGG module's lesion-level call (`molecular-subtyping-LGAT`) (tumors with no other source). `lesion_refinement` adds the LGG lesion-level subtype (`KIAA1549-BRAF`, `BRAF V600E`, `NF1`, ...) on top of the methylation family when LGG made a specific call. TP53 and gene annotations are not part of the subtype.

`04-compare-to-legacy.R` compares `molecular_subtype_cns.tsv` with the legacy per-module subtypes in `molecular-subtyping-pathology/results/compiled_molecular_subtypes.tsv`, per `match_id`, and writes `results/legacy_comparison.tsv` with an outcome of `same`, `same (lesion refinement)`, `newly assigned`, `different`, `legacy only (specific or unassigned)`, or `new only`.

### Gene annotations

`input/annotation_genes.tsv` lists genes whose status is reported as an annotation for the subtypes matching `applies_to_subtype_regex` (currently `EGFR` for `^DMG`). Each gets lowercase columns `<gene>_mutation_status` (`mutated`, `not detected`, `no DNA assay`) with `<gene>_variants` (any coding variant), `<gene>_amplification` (a copy-number segment overlapping the gene above twice the ploidy, the `focal-cn-file-preparation` convention; `no CNV data` when there is no segment file), and `<gene>_activating_alteration` (a variant matching `alteration_pattern`, or amplification).

For EGFR the activating alterations are amplification, A289/G598/R108 missense variants, and exon 20 insertions or duplications, as described for the `DMG_EGFR` methylation class.

`h3k27_altered_mechanism` lists, for DMG, which DNA-detectable routes to the H3 K27-altered phenotype were found: `H3 K27M/I` and/or `EGFR alteration`. EZHIP overexpression needs expression data and is not assessed, so tumors with neither are `none detected (EZHIP overexpression not assessed)`, and tumors with no DNA are `no DNA assay`.

### Usage

```
Rscript --vanilla 01-build-methylation-calls.R
Rscript --vanilla 02-prepare-lesion-tables.R
Rscript --vanilla 03-validate-subtypes.R
Rscript --vanilla 05-compile-subtypes.R
Rscript --vanilla 04-compare-to-legacy.R
```
