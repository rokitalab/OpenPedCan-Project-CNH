# Subtyping summary

Tumors in the subtype table: 3938; with a subtype: 3630.

## Change from the previous release

1770 tumors changed their subtype. Change type by verdict (`subtype_changes_from_release.tsv`, transitions in `subtype_change_transitions.tsv`, the tumors to review in `subtype_changes_needing_review.tsv`):

| change_type | review | supported | supported, diagnosis uninformative | no change in meaning |
|---|---|---|---|---|
| newly subtyped | 29 | 855 | 88 | 0 |
| reclassified | 4 | 90 | 25 | 0 |
| refined | 6 | 225 | 37 | 0 |
| vocabulary | 0 | 0 | 0 | 411 |

## Methylation mismatches

`input/possible_methylation_mixups.tsv`: likely = genomics supports the diagnosis (the arrays are collected for removal, `methylation_likely_mismatch.tsv`); possible = unresolved; kept = reviewed and the methylation call stands.

| mismatch_call | n |
|---|---|
| likely | 10 |

Likely mismatches:

| sample_id | pathology_diagnosis | molecular_subtype | methyl_class | methyl_score |
|---|---|---|---|---|
| 7316-3402 | Craniopharyngioma | CRANIO, ADAM | PA MID | 0.918 |
| 7316-166 | Craniopharyngioma | CRANIO, ADAM | PA PF TYP | 0.919 |
| 7316-1473 | Supratentorial or Spinal Cord PNET | Rare CNS, FOXR2 NBL | pedHGG_RTK1B | 0.998 |
| 7316-1969 | Ependymoma | EPN, ST ZFTA | pedHGG MYCN | 0.927 |
| 7316-3636 | High-grade glioma/astrocytoma (WHO grade III/IV) | HGG, H3 wildtype | ATRT_SHH | 0.927 |
| 7316-1856 | Ewings Sarcoma | Rare CNS, EWS | MB_G34_VIII | 1 |
| 7316-9842 | High-grade glioma/astrocytoma (WHO grade III/IV) | DMG, H3 K28-mutant | A_IDH_HG | 0.994 |
| 7316-7113 | Craniopharyngioma | CRANIO, ADAM | PA MID | 0.976 |

## Flags of the checks (`subtyping_checks.tsv`)

| area | check | n |
|---|---|---|
| methylation | every high-confidence class is in the crosswalk | 10 |
| subtype | tumors flagged for review (diagnosis discordant, classifiers disagree, ...) | 85 |
| subtype | methylation contradicts the diagnosis and is not yet reviewed (the subtype is kept; review the diagnosis; add reviewed ones to input/reviewed_diagnosis_discordance.tsv) | 19 |
| subtype | H3 K27/K28 variant outside DMG/DHG (a mix-up or a real H3 K27M tumor) | 2 |
| subtype | cell line subtype differs from its tumor | 13 |
| subtype | neuroblastoma MYCN status differs between copy number, expression and the clinical report | 19 |
| subtype | likely methylation mismatches (genomics supports the diagnosis): collect and remove from the cohort (results/methylation_likely_mismatch.tsv) | 10 |
| subtype | SHH gamma calls older than 5 years (check the methylation score and the UMAP) | 2 |
| integrate | plot-mapping rows added automatically and still to review (status column) | 38 |

