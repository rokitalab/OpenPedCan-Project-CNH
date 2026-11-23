# Results of the molecular subtyping

The table of record is `molecular_subtype.tsv` (one row per tumor, `match_id`).
The totals are in `../../molecular-subtyping-integrate/results/summary.md`.

**Final tables**
- `molecular_subtype.tsv`: subtype, its source, the lesion evidence, the flags.
- `subtypes_for_review.tsv`: high-confidence methylation calls that need a look (discordant diagnosis, contradicted by a lesion).

**Methylation**
- `methylation_calls.tsv`: every methylation biospecimen with its call (DKFZ v12.8 or NIH v3.1) and the subtype it maps to.
- `methylation_subtype_validation.tsv`: each high-confidence call with its lesion check (confirmed, not detected, no assay, no rule) and the diagnosis concordance.
- `methylation_downgraded.tsv`: calls set aside because a lesion contradicts them and the call is weak (NIH only, or DKFZ below 0.85). The tumor gets another source (lesion fallback or the next evidence).
- `methylation_overridden_by_lesion.tsv`: calls where the lesion decides (a definitive rule, for example CTNNB1 in a craniopharyngioma called papillary).
- `methylation_mismatch_calls.tsv`: every methylation mismatch call (`likely`, from the lesions and the diagnosis or from the reviewer), with its evidence, `call_source` and `found_from_lesions`.
- `methylation_likely_mismatch.tsv`: tumors whose methylation array is probably from another specimen (genomics supports the diagnosis); to collect and remove from the cohort.
- `methylation_crosswalk_gaps.tsv`: high-confidence methylation classes that have no subtype name in the crosswalk.

**Subtypes without methylation**
- `lesion_fallback_subtypes.tsv`: subtypes from a lesion (H3, IDH, fusions, ...) when there is no methylation call.
- `lesion_negative_hgg.tsv`: HGG with a DNA assay and no lesion (HGG, H3 wildtype).
- `lgg_lesions.tsv`: the LGG components (KIAA1549-BRAF, other MAPK, NF1, FGFR, ...) per tumor.
- `gcimp_scores.tsv`: G-CIMP probe score per methylation array (`08-score-gcimp.R`), evidence for the IDH-mutant glioma classes.
- `free_text_subtypes.tsv`, `free_text_conflicts.tsv`: subtypes from the pathology free text and where it conflicts with methylation.
- `embryonal_subtypes.tsv`: ETMR and CNS embryonal tumors from the TTYH1 fusion, C19MC amplification and LIN28A.
- `solid_subtypes.tsv`: neuroblastoma (MYCN), Ewing sarcoma and the other solid tumors.

**Medulloblastoma**
- `mb_expression_classifier.tsv`, `mb_validation.tsv`, `mb_discordant.tsv`, `mb_shh_*.tsv`: expression classifier, SHH subtype, validation.
- `mb_umap_coordinates.tsv`, `mb_pca_coordinates.tsv`: coordinates of the plots in `../plots`.

**Other**
- `h3_variant_review.tsv`: H3 K27/K28 variants in tumors with a non-DMG class.
- `chordoma_smarcb1_status.tsv`, `dlgnt_lesions.tsv`, `cell_line_*.tsv`: chordoma SMARCB1, leptomeningeal glioneuronal lesions, cell lines against their tumors.
- `legacy_comparison.tsv`: this subtype against the archived modules.
- `subtyping_checks.tsv`: the checks run on new data (`utils/check-subtyping.R`).

The hand-curated review lists are in `../review/`.
- `run_provenance.tsv`: the md5 of the input files the results were made on (checked by `utils/check-subtyping.R`).
