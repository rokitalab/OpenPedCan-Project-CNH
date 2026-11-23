---
name: molecular-subtyping-review
description: Review the OpenPedCan molecular subtyping (analyses/molecular-subtyping, molecular-subtyping-integrate) when new data or a new release comes in. Runs the module and integrate steps in the project image, runs the checks, compares with the previous release, and flags what needs a decision. Use for a data release, a methylation or histologies update, or any change to the subtyping rules.
---

# Molecular subtyping review

Methylation sets the subtype; lesions confirm it and subtype tumors with no
methylation; the diagnosis is compared with the subtype and flagged, never used to
override methylation. Read `analyses/molecular-subtyping/README.md` and
`DECISIONS.md` (the decision log, with the pending list) before changing anything.
The integrate step derives integrated diagnosis, broad and short histology and the
cancer group from the subtype by rule (`analyses/molecular-subtyping-integrate`).
The legacy per-diagnosis modules are archived in `analyses/archive/`; the LGG, HGG
and Ewing fusion calls (with reciprocal fusions) are computed in the module
(`utils/fusion-calls.R`), so the old fusion summary files are not needed.

## Release order

`histologies-base.tsv` (no subtypes) -> independent samples -> `fusion_filtering`
(makes `fusion-putative-oncogenic.tsv`, an input of the subtyping) -> copy number ->
`molecular-subtyping` -> `molecular-subtyping-integrate` (writes `histologies.tsv`) ->
`figures/scripts/extend-plot-mapping.R` and the tables that group by subtype or
cancer group. `scripts/generate-analysis-files.sh` runs the pre-release steps up to
the fusion filtering. Upstream name and mapping fixes belong in the
`opc-histologies-qc` repo (the QC rerun is a separate step); the stage-3 rules of
`input/diagnosis_rules.tsv` bridge the gap until the base file has them.

## Run (project image)

```sh
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -e OPENPBTA_SUBSET=0 \
  -v "$(pwd)":/rocker-build -w /rocker-build/analyses/molecular-subtyping \
  pgc-images.sbgenomics.com/rokita-lab/openpedcanverse:latest bash run-molecular-subtyping.sh   # about 3 min
# then in analyses/molecular-subtyping-integrate: bash run-subtyping-integrate.sh
```

Use the parquet tables (`bash scripts/download-parquet.sh`) so the SNV and copy
number tables are queried, not read whole (the image has 21 GB; reading the whole
MAFs runs out of memory). Results written by the image differ from a laptop run
only in float digits and the sort order of pasted strings: discard those diffs.

## Check, in this order

1. `Rscript --vanilla analyses/molecular-subtyping/utils/check-subtyping.R`
   writes `results/subtyping_checks.tsv`. Report every FLAG with its count and
   examples. A flag is not always a problem (for example 13 classifier-version
   disagreements with a concordant diagnosis are low concern); say which are which.
2. `Rscript --vanilla analyses/molecular-subtyping-integrate/02-compare-to-release.R <previous release histologies.tsv>`
   lists every subtype that changed (`results/subtype_changes_from_release.tsv`),
   its change type, evidence and verdict. Read `subtype_changes_needing_review.tsv`
   and say whether each theme makes sense.
   Also read `results/legacy_comparison.tsv`: the tumors whose subtype is `different`
   from the archived modules need a `reviewed` entry (see Reviews below); a legacy
   label that is the new subtype without its lesion term is the same.
3. For medulloblastoma and other methylation-defined groups, plot the arrays
   (`10-plot-mb-methylation.R`, UMAP and PCA, and
   `11-plot-mb-highlight.R <sample_id> ...`; they run with the module and need the
   beta values) and check
   that calls from age or lesions (not methylation) sit with the methylation calls
   of the same subtype. `08-score-gcimp.R` scores the G-CIMP probes of every array:
   evidence for the IDH-mutant glioma classes, not a call.
4. New tumors: confirm each has a harmonized diagnosis, broad and short histology
   and cancer group, and each subtyped tumor an integrated diagnosis
   (the integrate checks in `check-subtyping.R`). The checks flagged 10 of 32 at the
   last run (all known, among them the manual reviews that the automated subtype now agrees
   with: retire them when a decision is no longer needed); report any new flag as new.
5. Plot groups: run `figures/scripts/extend-plot-mapping.R` (adds the missing rows
   with a `status` to review, removes rows with no tumor) and
   `figures/scripts/plot-palette.R` (the palette next to the plot groups and the
   closest pair of colors; two groups must not share a color).
6. Reproducibility: after a rerun on the same release, `git diff --stat` of the two `results/`
   folders must be empty (plots excepted), and the `reviews` checks of `check-subtyping.R` must
   pass: no orphan or empty review, no stale manual or mix-up row, and the inputs are those of
   `results/run_provenance.tsv`. A change that is not a decision is a bug; a change that is a
   decision is committed with its regenerated results. CI cannot do this comparison (it runs
   the small `testing` release).
7. For a tumor under review, `Rscript --vanilla utils/tumor-evidence.R <sample_id>`
   (run from `analyses/molecular-subtyping`) prints the clinical record, the DKFZ
   v12.8 and NIH v3.1 calls, the subtype row, variants, fusions, copy number, the
   1p/19q arm values and the G-CIMP score in one place. Quote only these versions.

## Rules to enforce

- Methylation inputs are DKFZ v12.8 and NIH v3.1 only (never v11, v12 or NIH v2).
  High confidence: DKFZ subclass score >= 0.8, else NIH superfamily and class
  >= 0.9. Scores < 0.75 are not evidence in review lists. NIH v3.1 class names
  missing from the crosswalk are added in `input/nih_v3.1_class_aliases.tsv`;
  `PA PF VAR` stays unresolved (no published definition).
- Methylation is never overridden by the diagnosis alone. A mismatch is flagged
  (`subtype_review_flag`), listed for pathology review, and a reviewed decision goes
  in `input/manual_subtype_review.tsv` (keyed on `sample_id`, highest precedence,
  with the report or reviewer as the source). Methylation mismatches
  are found by the module (`results/methylation_mismatch_calls.tsv`) and flagged only
  when the high-confidence class is completely mismatched: the genomics supports the
  diagnosis or another entity and not the class (`found_from_lesions`), or the class is
  mutually exclusive with a lesion (an IDH-mutant glioma class with an H3 K27 or G34
  variant). On new data, review every candidate: read the pathology diagnosis, the lesions
  and the RNA (`utils/tumor-evidence.R`) and the tumor purity (a missing lesion can be low
  purity), and record a `likely` decision with its evidence in
  `input/possible_methylation_mixups.tsv`; the reviewer decision takes precedence. The
  arrays of a likely mismatch are not used for the subtype (the tumor is subtyped from its
  lesions, a reviewer decision or the diagnosis) and are listed in
  `results/methylation_likely_mismatch.tsv` to be removed from the release. A class that
  is only discordant with the diagnosis, with genomics that does not contradict it (no
  DNA or RNA, low purity, no lesion either way) or that matches it, is not a mismatch: the
  methylation subtype stays with the diagnosis-discordance flag, and the notes go in
  `review/methylation_diagnosis_discordance_notes.tsv`. An EPN PF-A class with an H3 K27M
  variant is not a mismatch.
- A definitive lesion (CTNNB1 exon 3, BRAF V600E in craniopharyngioma) of another
  family contradicts a methylation call when the diagnosis agrees with the lesion
  and not with the call: flagged, and a weak call is set aside for the lesion
  fallback.
- Control-methylation tumors get no lesion fallback, except an H3 K27/K28 variant
  (DMG) and craniopharyngioma free text. Controls and INFLAM_ENV are `Control`.
- TP53 status and hits, EGFR (DMG) and EZHIP are annotations, not part of the
  subtype. EZHIP and LIN28A cannot be assessed on exome capture RNA.
- 1p and 19q are ploidy-aware (`02-prepare-lesion-tables.R`): neutral segments count
  as the ploidy, so a few focal losses are not an arm loss; the GISTIC values stay in
  `1p_gistic` and `19q_gistic`. A complete 1p/19q codeletion is incompatible with an
  IDH-mutant astrocytoma class (check ATRX and TP53 variants too).
- Medulloblastoma marker gene expression (WNT, SHH, Group 3, Group 4) is evidence only, in
  `lesion_evidence`; it helps where there is RNA and no methylation (cell lines) but never
  sets a subtype. The flag `methylation and genomics agree, diagnosis differs` means the
  molecular data may not belong to the clinical record: review it with pathology, it is not
  a mismatch of the array alone.
- Lesion-level labels for the LGG, GNT, GNG and SEGA families come from
  `04-lgg-lesions.R`; defined entities (DNET, SEGA) are not split. TSC is germline
  and not subtyped here.
- Medulloblastoma: the group (WNT, SHH, Group3, Group4) is methylation, with the
  expression classifier as an independent check and the subtype where there is no
  methylation. SHH subtype (alpha, beta, gamma, delta): DKFZ v12.8 subclass at a
  mean score >= 0.7, then NIH v3.1 at >= 0.9, then age and molecular criteria.
  Tumors that are MB only by the expression classifier are not MB.
- Neuroblastoma: MYCN copy number decides; the clinical report outranks expression;
  an unresolved conflict stays unsubtyped. Ewing sarcoma: methylation, then a
  hallmark fusion (EWSR1--FLI1, EWSR1--ERG, FUS--ERG), then the diagnosis alone;
  CNS (PBTA, which includes DGD) tumors are `Rare CNS, EWS`.
- PNOC003 DIPG tumors are DMG whatever their H3 status. Cell lines are subtyped on
  their own `match_id` and compared with the tumor.
- Integrated diagnosis, broad and short histology are derived by rule
  (`input/subtype_diagnosis_families.tsv`, `subtype_diagnosis_overrides.tsv`); the
  first comma token of the integrated diagnosis is the cancer group, standardised by
  `input/diagnosis_rules.tsv`. A subtype label with no broad or short histology
  keeps the old diagnosis' value: add it to the tables (the check lists them). Every
  tumor needs a harmonized diagnosis, broad and short histology. `X, To be
  classified` labels keep the label and have the harmonized diagnosis as their
  integrated diagnosis.
- The cancer group is a lookup table, `input/cancer_group_table.tsv` (diagnosis token
  or exact harmonized diagnosis to cancer group; an empty group is deliberate). A
  diagnosis not in it is listed in `results/cancer_group_unmapped.tsv`: add a row.
  The figures group tumors with `figures/palettes/plot-mapping.tsv` (cohort, broad
  histology and cancer group to plot group) through `figures/scripts/plot_group.R`;
  a new combination needs a row there too: `figures/scripts/extend-plot-mapping.R`
  adds them (PBTA cancer groups with fewer than 5 tumors go to Other tumor); rows
  with a `status` were added automatically and are to be reviewed. The table is for
  PBTA, TARGET and GMKF only (DGD, TCGA, PPTC, Maris and CPTAC are left out); all
  neuroblastic tumors are the plot group `Neuroblastoma` (`NBL`); hex codes are in
  capitals; two plot groups never share a color.
- Keep the vocabulary: `HGG, ` / `LGG, ` / `GNT, ` prefixes on family labels,
  WHO CNS5 names (`H3 K28-mutant`), chordoma is `CHDM`.

## Reviews

A review is a row in a table, never code: `input/manual_subtype_review.tsv` (a
subtype set by a report or reviewer), `input/possible_methylation_mixups.tsv`,
`input/reviewed_diagnosis_discordance.tsv` and `input/legacy_comparison_review.tsv`
(`reviewed` = JLR and a note, only for `different` outcomes). Write what a pathology
report shows only when the user has said so, and never that a report "is correct"
unless they said it. Name the DKFZ v12.8 and NIH v3.1 scores in the note. Hand-curated
lists that no script reads (lesions still to check, the WHO histology review, the PPTC
MYCN template) live in the `review/` folders, not in `results/`; one-off checks are
scripts (`utils/tumor-evidence.R`, `08-score-gcimp.R`), not pasted output.

A discordance that was reviewed with the pathology report (or comes from an obsolete
diagnosis such as PNET) and where the methylation call is accepted goes in
`input/reviewed_diagnosis_discordance.tsv`; it clears the flag.

## Flag, do not decide

Ask the user (do not choose) when: a methylation call contradicts a diagnosis
(pathology review), a cancer group or histology name is new, a cutoff would change
calls, a class is not in the crosswalk, or a tumor looks like a sample mix-up.
Give the sample, the evidence (classifier calls and scores, lesions, age,
expression markers, diagnosis) and a recommendation. Record open items in the
Pending list of `DECISIONS.md` and the decision once it is made.

## Never

Do not push or open a pull request without an explicit "push" from the user. Commit by
explicit path or `git add -u`; pathology report PDFs left in the repository root are
never committed. Do not
override a methylation call from the diagnosis, and do not edit
`manual_subtype_review.tsv` without a pathology report or reviewer decision.
