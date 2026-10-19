## Integrate molecular subtyping output from pathology feedback

**Author of code and documentation:** [@kgaonkar6](https://github.com/kgaonkar6)

In this repo, we add molecular subtype for molecular_subtype from all subtyping modules and integrated_diagnosis, short_histology, broad_histology, and Notes from `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv`. A column `cancer_group` is added to provide broader terms derived from `harmonized_diagnosis` which can be used to generate figures.

### Usage
```sh
bash run-subtyping-integrate.sh                          # derive the diagnoses and write results/histologies.tsv
Rscript --vanilla 02-compare-to-release.R [previous release histologies.tsv]   # subtype and cancer group changes
Rscript --vanilla 03-summary.R                           # results/summary.md
Rscript --vanilla 04-cell-line-fidelity.R                # results/cell_line_fidelity.tsv (PBTA lines)
```
Run it after `analyses/molecular-subtyping`, which it reads. `02` to `04` and `../molecular-subtyping/utils/check-subtyping.R` are the review steps of a release (see the `molecular-subtyping-review` skill).

### Module contents

`01-integrate-subtyping.Rmd` integrates results from compiled results in `compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv` to `histologies-base.tsv`

`00-derive-integrated-diagnosis.R` derives `integrated_diagnosis`, `broad_histology` and `short_histology` from the `molecular_subtype` label of `molecular-subtyping` by rule, so the diagnosis follows the subtype instead of being set separately by each subtyping module:

- `input/subtype_diagnosis_overrides.tsv`: an exact label (e.g. `DMG, H3 K28-mutant, EGFR`) and its WHO CNS 2021 name, legacy strings kept where they existed.
- `input/subtype_diagnosis_families.tsv`: the family default for every other label, `<integrated_prefix>, <the label after its first comma>`, with the family's broad and short histology. The first comma-separated token of `integrated_diagnosis` is the `cancer_group` in `01-integrate-subtyping.Rmd`, so the prefixes are the existing cancer group names.
- Outputs: `results/subtype_diagnosis_map.tsv` (one row per label, with the legacy `integrated_diagnosis` for comparison and `differs_from_legacy`) and `results/integrated_diagnosis.tsv` (one row per tumor).

`01-integrate-subtyping.Rmd` takes every tumor of the `molecular-subtyping` module (CNS and solid tumors) from `results/integrated_diagnosis.tsv` (subtype, diagnosis, histologies, `Notes`, `subtype_source`, `subtype_review_flag`, `subtype_alterations_confirmed`, `alteration_confirmed` and the methylation class as `molecular_subtype_methyl`), and keeps the frozen compiled `molecular-subtyping-pathology` result (in `analyses/archive/`, with the other archived legacy subtyping modules; see `analyses/archive/README.md`) for every other tumor (non-CNS tumors and the pathology free-text diagnoses). The DLGNT diagnosis from the glioneuronal free-text step is kept. Neuroblastic tumors take the diagnosis name of the tumor (neuroblastoma, ganglioneuroblastoma, ganglioneuroma) with the MYCN status. `results/cancer_group_changes.tsv` lists, per previous and new cancer group and diagnosis, the biospecimens whose cancer group changed against the previous release table, for review. The pathology free-text steps (meningioma, juvenile xanthogranuloma, CNS lymphoma, choroid plexus, glioneuronal, LGAT) set diagnoses and histologies that are not molecular subtypes and stay in the archived `molecular-subtyping-pathology` results; they are applied on top of the derived diagnosis. In particular, diffuse leptomeningeal glioneuronal tumor (DLGNT, 5 samples from the glioneuronal free-text step) is applied here, not in `molecular-subtyping`, which keeps the methylation label for those tumors (7316-2208 is `GNT, MAPK`).

### Cancer group table

The cancer group of every harmonized diagnosis is a row of `input/cancer_group_table.tsv`, read by `utils/lookup-cancer-group.R`; there is no code that standardises cancer groups. A row is either a diagnosis token (`diagnosis_token`: the text before the first comma, without the part in brackets, so every subtype label of a family resolves with its family: `Low-grade glioma, KIAA1549-BRAF` is `Low-grade glioma`) or, where one diagnosis differs from its token, the exact `harmonized_diagnosis` (an exact row wins). An empty `cancer_group` is deliberate: the diagnosis has no cancer group (benign tumors, non-tumor tissue). A diagnosis in neither kind of row is unmapped: its cancer group is missing and it is listed in `results/cancer_group_unmapped.tsv` (empty when every diagnosis is in the table), to be added to the table. `note` records why a row exists. The table reproduces the earlier rules for all but 25 tumors, where the rules gave one diagnosis different groups (for example `Schwannoma` as `Neurofibroma/Plexiform`, `Germinoma, KIT mutant` as `Meningioma`, `Low-grade glioma, BRAF/MAPK` as `Ependymoma`) and the MPNST subtypes had none; the table gives each diagnosis one group.

The plot groups of the figures are a table too, `figures/palettes/plot-mapping.tsv` (see `figures/README.md`): cohort, broad histology and cancer group to plot group, abbreviation and colors, for the PBTA, TARGET and GMKF cohorts. `figures/scripts/extend-plot-mapping.R` adds the combinations of the table of histologies that are missing and removes the rows that have no tumor; `utils/check-subtyping.R` checks that every combination is in it.

### Diagnosis rules

The harmonized diagnoses of three cancers and the broad and short histologies of a few cancer groups are standardised by `input/diagnosis_rules.tsv`, applied by `utils/apply-rules.R`, so that a change is a row in a table and not a line of code. Each rule sets `target` to `result` (a value, `NA`, or `=<column>` to copy a column) for the rows that meet all of its conditions, when no earlier rule of the same stage and target has set the row. Stage 0 sets harmonized diagnoses (such as `Ewing sarcoma`) and stage 3 the histologies from the cancer group, after the cancer group table. Condition columns match a column exactly (`|` separates alternatives and a leading `!` negates), `<column>_regex` is a regular expression and `<column>_na` tests for a missing value.

### Results

- `histologies.tsv`: the base histologies with the molecular subtype, integrated diagnosis, histologies and cancer group (the table of record of this module).
- `integrated_diagnosis.tsv`, `subtype_diagnosis_map.tsv`: the diagnosis derived for each tumor and for each subtype label; `cancer_group_unmapped.tsv`, `discrepancies_to_check.tsv`: what still needs a table row or a look.
- `subtype_changes_from_release.tsv`, `subtype_change_summary.tsv`, `subtype_change_transitions.tsv`, `subtype_changes_needing_review.tsv`, `cancer_group_changes.tsv`: the comparison with the previous release (`02`), one row per tumor, per change type and verdict, per label transition, the tumors to review and the biospecimens whose cancer group changed.
- `summary.md`: the totals (`03`); `cell_line_fidelity.tsv`: the PBTA cell lines against their parental tumors (`04`).
- `harmonized_diagnosis_cancer_group_table.tsv`, `cancer_group_table_for_README.tsv`, `pediatric_cancer_groups.tsv`, `tcga_cancer_groups.tsv`, `subtype_tumor_counts_table.tsv`: the cancer group tables.
- `review/` has the hand-curated lists (the WHO CNS5 review of the histology mapping).

A rule can override what the base histologies hold until the base file has it: the stage-3 rules of `input/diagnosis_rules.tsv` set the CNS5 broad histology names and `LGG` / `HGG` short histologies (remove them when the base file has them), and a broad histology that is wrong for a cancer group is a rule too (only cavernoma is `Benign tumor`; MPNST keeps its own).
