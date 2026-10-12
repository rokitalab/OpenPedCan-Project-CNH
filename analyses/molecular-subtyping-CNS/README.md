# Molecular subtyping of CNS tumors using methylation and genetic lesions

**Module authors:** Jo Lynne Rokita

## Overview

This module builds one resolved methylation call per biospecimen, which molecular subtyping modules consume instead of each re-implementing score thresholds and classifier-version aliasing.

## Workflow

Scripts run in numeric order: `01` builds the methylation calls, `02` prepares lesion tables, `03` validates and applies the lesion fallback, `04` computes the LGG lesion-level subtypes, `05` validates medulloblastoma, `06` compiles the final table (the last subtyping step), and `07` compares it with the legacy modules. Shared functions are in `utils/` (`lesion-hits.R`, `dx-concordance.R`, `tp53-annotation.R`, `gene-annotation.R`, and `duckdb-connect.R`, which opens duckdb with `shared_home = FALSE` where the installed version supports it). The scripts use duckdb and duckplyr from the project docker image (duckdb 1.5.6, duckplyr 1.2.1, dbplyr 2.6.0). With a duckplyr frame use `summarise(.by = ...)`, since `group_by()` then `summarise()` is not supported there.

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

`06-compile-subtypes.R` writes `results/molecular_subtype_cns.tsv`, one row per `match_id`. Subtype precedence is high-confidence methylation, then lesion-only fallback, then the LGG module's lesion-level call (`molecular-subtyping-LGAT`) (tumors with no other source). `lesion_refinement` adds the LGG lesion-level subtype (`KIAA1549-BRAF`, `BRAF V600E`, `NF1`, ...) on top of the methylation family when LGG made a specific call. TP53 and gene annotations are not part of the subtype.

`07-compare-to-legacy.R` compares `molecular_subtype_cns.tsv` with the legacy per-module subtypes in `molecular-subtyping-pathology/results/compiled_molecular_subtypes.tsv`, per `match_id`, and writes `results/legacy_comparison.tsv` with an outcome of `same`, `same (lesion refinement)`, `newly assigned`, `different`, `legacy only (specific or unassigned)`, or `new only`.

### Medulloblastoma

Methylation sets the group (`MB, WNT`, `MB, SHH`, `MB, Group3`, `MB, Group4`), as for the other families, and lesion rules confirm it (CTNNB1 exon 3 hotspot for WNT; PTCH1, SMO, SUFU or GLI2 amplification for SHH; MYC amplification for Group 3). `05-validate-mb.R` compares it with the RNA-seq expression classifier (`medulloPackage`, read from `molecular-subtyping-MB/results/mb-classified.rds`) and adds the SHH subtype (alpha, beta, gamma, delta) from `molecular-subtyping-MB`:

- `mb_expression_group` and `mb_concordance` (`concordant`, `discordant`, `methylation only`, `expression only`); discordant tumors are flagged and listed in `results/mb_discordant.tsv`
- `mb_shh_subtype` for SHH tumors, with the DKFZ SHH class vs subtype cross-tab printed by the script
- tumors with no methylation call take the expression classifier group (`subtype_source = expression classifier`)

The expression classifier and SHH subtyping stay in `molecular-subtyping-MB`; this module reads their results.

### Manual review

`input/manual_subtype_review.tsv` holds subtypes set by a pathology report or clinical review, keyed on `sample_id`. They have the highest precedence (`subtype_source = manual review`). The automated call stays in `automated_subtype`, `manual_conflict` is TRUE when it differs from the reviewed subtype, and a reviewed tumor is no longer listed as needing review.

### Craniopharyngioma and chordoma

- **CRANIO:** methylation (`CPH_ADM`, `CPH_PAP`) sets `CRANIO, ADAM` or `CRANIO, PAP`, and a CTNNB1 exon 3 hotspot or in-frame deletion confirms ADAM, BRAF V600E confirms PAP. **In craniopharyngioma the lesion decides** (`lesion_overrides_methylation` in the rules): the two lesions are mutually exclusive, and in this cohort 11 methylation-PAP calls (DKFZ 0.81 to 0.98, some NIH v3.1-only) all carried CTNNB1 exon 3 with no BRAF V600E, and all 5 with a pathology report were adamantinomatous. The overridden calls are listed in `results/methylation_overridden_by_lesion.tsv` with the original methylation subtype. Tumors with no methylation are subtyped from the lesion only when the diagnosis supports it (`fallback_requires_dx`), since CTNNB1 and BRAF V600E occur in many other tumors.
- **Other families:** a methylation subtype contradicted by a definitive lesion of another subtype in the same family (`definitive` rules, e.g. CTNNB1 for MB WNT) is flagged (`lesion_contradiction`, `subtypes_for_review.tsv`) and kept, unless the call is also weak (NIH-only, or a DKFZ score below 0.85): then it is set aside (`results/methylation_downgraded.tsv`) and the tumor is subtyped by the lesion fallback.
- **Chordoma:** methylation `CHORDM` and `CHORDM_PD` give `CHDM` and `CHDM, poorly differentiated` (the label is `CHDM`, as in the legacy module; the NIH v3.1 class is added through `input/extra_methylation_classes.tsv` because the crosswalk lacks it). SMARCB1 focal status and expression from `molecular-subtyping-chordoma` are annotations (`smarcb1_focal_status`, `smarcb1_expression`), and a chordoma-diagnosis tumor with focal SMARCB1 loss and no methylation is `CHDM, poorly differentiated` (`subtype_source = SMARCB1 focal loss`). The conventional vs poorly differentiated call is a pathology review (legacy ticket 608), kept as manual rows where it differs from methylation (7316-921, 7316-3632, 7316-4062).

### HGG annotations and umbrella group

- `h3_status` (`mutant`, `wildtype`, `not assessed`) and `h3_variant` for `HGG`, `DMG` and `DHG` subtypes. Wildtype means a DNA assay with no H3 K27M/I or G34/G35 R/V variant.
- `subtype_group` rolls the pediatric H3-wildtype methylation classes (`HGG, RTK`, `MYCN`, `MES`, `chr6 CTX`, `HGAP`, `HGNET`, ...) up into `HGG, H3 wildtype` (`input/subtype_groups.tsv`); all other subtypes are their own group. The specific class stays in `molecular_subtype_cns`. A tumor in one of these classes with an H3 variant is not grouped and is flagged for review (`h3_group_conflict`).
- `h3_variant_conflict` flags a methylation subtype that is not `DMG`/`DHG` in a tumor with an H3 K27M/I or G34/G35 R/V variant. `results/h3_variant_review.tsv` lists every such tumor, including those left without a subtype because the methylation call is a control (low purity) and the lesion fallback is not applied.
- `HGG, H3 wildtype` is also assigned by exclusion (`subtype_source = lesion negative`) to PBTA tumors with an HGG diagnosis, a DNA assay and no methylation call, lesion rule hit, or control, with no H3 or IDH variant. Results are in `results/lesion_negative_hgg.tsv`.

### Free-text diagnosis

PBTA tumors whose free-text diagnosis names a diffuse midline glioma or an H3 K27 alteration (`diffuse midline`, `DMG`, `H3 K27`), with no methylation call or lesion-based subtype and a methylation call that is not a control, are `DMG, H3 K28` (`subtype_source = pathology free text`). `-mutant` is not used because it needs the variant. `DIPG` alone is not enough. The diagnosis can say "K27M mutant" without a variant, so a DMG with a DNA assay, an EGFR activating alteration and no H3 variant is `DMG, H3 K28-altered, EGFR` (the same applies to methylation-only `DMG, H3 K28` tumors). A DMG with a confirmed H3 variant and an EGFR activating alteration carries both labels, `DMG, H3 K28-mutant, EGFR`. EZHIP overexpression (z >= 3 within library) is added after any EGFR term (`DMG, H3 K28-mutant, EZHIP`, `DMG, H3 K28-altered, EGFR, EZHIP`), and a methylation-only `DMG, H3 K28` with EZHIP overexpression and no H3 variant is `DMG, H3 K28-altered, EZHIP`. These tumors are no longer assigned the lesion-negative `HGG, H3 wildtype`, which contradicted the diagnosis (results in `results/free_text_dmg.tsv`).

### Craniopharyngioma and neurocytoma from the diagnosis

- **`CRANIO, ADAM` from the free text:** a PBTA tumor whose free-text diagnosis contains "adamantinomatous" and that no methylation, CTNNB1 or BRAF rule subtyped (a missing, low-confidence or control array) is `CRANIO, ADAM` (`subtype_source = pathology free text`). A control array only says the tissue was low purity, so it does not block this rule. A high-confidence methylation call of another class stands and is flagged `free-text diagnosis differs from methylation` (`results/free_text_conflicts.tsv`).
- **Neurocytoma:** ported from `molecular-subtyping-neurocytoma`. A `Neurocytoma` pathology diagnosis with primary site `Ventricles` is `CNC`, otherwise `EVN`; a CN or EVNCYT methylation call takes precedence (`subtype_source = pathology diagnosis`).

### DICER1

`Rare CNS, DICER1 SARC` is subtyped without methylation only from a DICER1 RNase IIIb hotspot missense (E1705, D1709, D1713, G1809, D1810, E1813) in a PBTA tumor with a compatible diagnosis (`RARE_DICER1_HOTSPOT`). A truncating or splice variant alone confirms a methylation call but does not subtype (the single nonsense variant in a hypermutator was a passenger), and a DICER1 hotspot in pineoblastoma is held, since that is `PB, miRNA processing`.

### Control methylation and lesion fallback

A tumor whose methylation call is a control class (CTRL_*, INFLAM_ENV) is not given the lesion fallback, with one exception: an H3 K27/K28 M or I variant is specific to DMG, so those tumors are `DMG, H3 K28-mutant`. HGG tumors with a control array stay unsubtyped. `HGG, H3 wildtype` by exclusion is only withheld when an H3 or IDH rule hits, not for any other rule (TP53, NF1, RB1).

### Derived cell lines

A derived cell line is subtyped from its own methylation, lesions and RNA, on its own `match_id`. It is then compared with the solid-tissue tumors of the same `sample_id` (or the same participant):

- `tumor_subtype`, `tumor_group` and `tumor_concordance` (`same`, `same group` when one label is the other plus an annotation term such as `-mutant` or `, EGFR`, `different`, or `no subtype to compare`). A `different` line is flagged `cell line differs from tumor`.
- `tumor_alterations_shared` lists the tumor's confirming alterations as `shared`, `not found in line` (the assay was available) or `not assessable (no DNA / RNA)`.
- A line with no subtype of its own (no methylation, no lesion rule, or only `HGG, H3 wildtype` by exclusion) that shares a confirming lesion with its tumor takes the tumor's subtype (`subtype_source = lesion shared with tumor`, Notes `Subtype from lesion shared with the tumor`). The `, EZHIP` term is kept only if the line overexpresses EZHIP.
- `results/cell_line_alterations.tsv` (written by `03`) is what each line carries, and `results/cell_line_tumor_concordance.tsv` is the comparison.

### PXA

DKFZ or NIH high-confidence PXA calls are split by BRAF V600E (issue 113): `PXA BRAF V600E` and `PXA BRAF V600E negative`, with plain `PXA` when there is no DNA assay. The `HGG` or `LGG` prefix follows the pathology diagnosis: an HGG diagnosis gives `HGG, PXA ...`. When the diagnosis is missing or `Other`, the free-text diagnosis is searched for a grade (anaplastic, WHO grade III/IV, high grade, glioblastoma). Otherwise the label is `LGG, PXA ...`. A PBTA tumor with no methylation, BRAF V600E and a PXA term in the free-text diagnosis is `HGG` or `LGG, PXA BRAF V600E` (BRAF V600E alone is not specific, so the free text is required; `subtype_source = lesion fallback`).

### Supratentorial pilocytic astrocytoma (PA_CORT / PA ST)

The DKFZ `PA_CORT` class (NIH v3.1 `PA ST`) is `LGG, PA ST`, and the family prefix follows the primary diagnosis, as for PXA: a ganglioglioma diagnosis gives `GNG, PA ST`, a DNET or glial-neuronal tumor diagnosis `GNT, PA ST`, and anything else `LGG, PA ST`. The DLGNT classes stay `GNT, MAPK`. The class carries its lesion term from the LGG lesion module (`04-lgg-lesions.R`) in the label, so these tumors are `LGG, PA ST, BRAF V600E`, `GNG, PA ST, KIAA1549-BRAF`, `GNT, PA ST, other MAPK` and so on. A tumor with no lesion term stays at `LGG, PA ST` / `GNG, PA ST` / `GNT, PA ST`.

`lesion_refinement_status` (LGG, GNG, GNT and SEGA families) says what the lesion module could see: `specific` (a lesion term), `none detected (DNA and RNA)`, or `incomplete: no DNA`, `incomplete: no RNA` or `incomplete: no DNA or RNA`. An incomplete tumor is still to be determined, not lesion-negative. The legacy `To be classified` labels are not used: a tumor without lesion evidence keeps its methylation label, and the status and `alteration_confirmed` carry the rest.

### Lesion-level labels for the LGG families

`molecular_subtype_cns` is the lesion-level label, as in the legacy LGAT module, for tumors subtyped by methylation:

- a generic MAPK bucket (`LGG, MAPK` from PA_MID, `LGG, BRAF/MAPK` from PA_INF, `GNG, MAPK` from GG) is replaced by the lesion term: `LGG, KIAA1549-BRAF`, `GNG, BRAF V600E`, `LGG, RTK`. The prefix stays the class prefix.
- `PA ST` gets the term appended: `GNG, PA ST, BRAF V600E`, `LGG, PA ST, KIAA1549-BRAF`.
- an entity that is defined as such (`GNT, DNET`, `LGG, SEGA`, `GNT, DGONC`, `GNT, DIG/DIA`, papillary and myxoid `GNT`) stays one subtype. Its lesion (for example the FGFR1 alteration of a DNET) is in `lesion_refinement` and `subtype_alterations_confirmed`, not in the label, so a DNET is not split into DNET, FGFR / DNET, RTK / DNET, NF1.
- a class that already is a lesion term (`FGFR`, `IDH`, `MYB/MYBL1`, `PXA BRAF V600E`) and the DLGNT classes are left as they are.

The tumors the LGG cohort subtypes by lesion alone (`subtype_source = LGG lesion`) already have the lesion-level label. `lesion_refinement` keeps the LGG module's label (its own prefix, from the pathology diagnosis), and `lesion_refinement_status` records `specific`, `none detected (DNA and RNA)` or `incomplete: no DNA / RNA`. The methylation class stays in `methyl_class`.

### Label prefixes

Subtype labels carry their family as a prefix: `HGG, ...`, `LGG, ...`, `GNT, ...`, `DMG, ...`, `DHG, ...`. Classes that were bare or used a different prefix are renamed in the crosswalk (`HGAP`, `HGNET`, `IHG` and `GBM, MES/RTK` become `HGG, ...`; `SEGA` becomes `LGG, SEGA`; bare `GNT` becomes `GNT, DGONC`, `GNT, DIG/DIA`, `GNT, DNET` or `GNT, NOS`). The DKFZ class is always kept in `methyl_class`. Ganglioglioma labels keep the `GNG, ...` prefix used by the LGG module.

### Embryonal: ETMR and CNS embryonal, NOS

Ported from `molecular-subtyping-embryonal` (`03`, cohort strings in `input/embryonal_path_dx_strings.json`). For embryonal-diagnosis tumors with no methylation call or other source:

| Subtype | Evidence |
|---|---|
| `ETMR, C19MC-altered` | a TTYH1 fusion (5' partner) with LIN28A overexpression or a C19MC amplification |
| `ETMR, NOS` | an ETMR diagnosis without the above, or LIN28A overexpression without a TTYH1 fusion |
| `CNS Embryonal, NOS` | the rest of the embryonal cohort (including CNS neuroblastoma diagnoses) |

LIN28A overexpression is the z-score within RNA library type, z >= 3 (the embryonal module used z > 3 within the embryonal subset). C19MC amplification is a focal chr19 segment overlapping the C19MC cluster (hg38 chr19:53.66-53.77 Mb) with segment mean > 1; the embryonal module called any chr19 gain "amplification", which flagged 9 tumors for segments elsewhere on chr19. TTYH1 fusions come from annoFuse and the DGD panel fusions, not the in-frame summary used by the module, so 4 tumors it called `ETMR, NOS` are `ETMR, C19MC-altered` here. Results are in `results/embryonal_subtypes.tsv`; `subtype_source = embryonal rules`.

### LGG lesion-level subtypes

`04-lgg-lesions.R` ports the LGAT module's lesion logic so that this module computes the evidence itself (`input/lgg_component_rules.tsv`; the LGAT diagnosis strings are in `input/lgg_path_dx_strings.json`). Components are KIAA1549-BRAF, BRAF V600E (V600E, codons 599/600, or a BRAF kinase-domain variant), other MAPK, RTK, FGFR (including FGFR1 kinase-domain duplication), NF1-somatic and NF1-germline, IDH, H3, MYB/MYBL1 and CDKN2A/B, joined in the LGAT order and prefixed `LGG`, `GNT`, `GNG` or `SEGA` by the diagnosis. They are `lesion_refinement` on the LGG / GNT / GNG / SEGA families, the subtype for LGG-cohort tumors with no other source (`LGG lesion`), and, through `lgg_alterations`, feed `subtype_alterations_confirmed`.

Parity with `molecular-subtyping-LGAT/results/lgat_subtyping.tsv`: of 617 tumors with a lesion-level label in both, 616 are identical. The differences are one extra NF1-germline call (the predisposition sits on another biospecimen of the same tumor) and two IDH calls that the LGAT module made without a complete DNA/RNA pair. The fusion summary file, synonymous-variant removal and component order follow the LGAT module. Computed with `fusion_summary_lgg_hgg_foi.tsv` and the consensus gene-level copy number, as there.

### Expression annotations

`input/expression_annotations.tsv` lists genes called as overexpressed (z-score of log2(TPM+1) within RNA library type, against PBTA tumor RNA-seq, at or above `z_threshold`): `EZHIP` for `^DMG` and `LIN28A` for `^(ETMR|CNS Embryonal)`. `exclude_library_groups` names library types in which a gene cannot be called: EZHIP is not assessable in exome-capture libraries (46 tumors get `ezhip_expression_status = not assessable (exome capture library)` and `h3k27_altered_mechanism` says so). Exome capture reads about 2.2 TPM for EZHIP in every tumor type, 0 of 150 reference tumors reach z >= 3 and there is one PF-A positive control, so a threshold cannot be calibrated. LIN28A is called in every library type.

### Gene annotations

`input/annotation_genes.tsv` lists genes whose status is reported as an annotation for the subtypes matching `applies_to_subtype_regex` (currently `EGFR` for `^DMG`). Each gets lowercase columns `<gene>_mutation_status` (`mutated`, `not detected`, `no DNA assay`) with `<gene>_variants` (any coding variant), `<gene>_amplification` (a copy-number segment overlapping the gene above twice the ploidy, the `focal-cn-file-preparation` convention; `no CNV data` when there is no segment file), and `<gene>_activating_alteration` (a variant matching `alteration_pattern`, or amplification).

The `DMG_EGFR` methylation class is its own subtype, `DMG, H3 K28-altered, EGFR`. When an H3 K27M/I variant is also found, the tumor is `DMG, H3 K28-mutant, EGFR`. `DMG, H3 K28-mutant` is used only when an H3 K27M/I variant was found (by DNA, or as the lesion-only call). The `DMG_K27` methylation class without that evidence stays `DMG, H3 K28`. The same applies to `DHG, H3 G35` (`DHG, H3 G35-mutant` when an H3-3A G34/G35 R/V variant was found). A `DMG_EGFR` tumor with no DNA assay keeps `DMG, H3 K28-altered, EGFR`, since the class itself defines it.

For EGFR the activating alterations are amplification, A289 and G598 missense variants, and exon 20 insertions or duplications, as described for the `DMG_EGFR` methylation class. Other coding EGFR variants, including R108K, are reported in `egfr_variants` only and do not set the subtype (see the literature table below).

`h3k27_altered_mechanism` lists, for DMG, which routes to the H3 K27-altered phenotype were found: `H3 K27M/I`, `EGFR alteration` and/or `EZHIP overexpression`. EZHIP overexpression is a z-score of log2(TPM + 1) within RNA library type (exome capture, poly-A, stranded), computed against PBTA tumor RNA-seq, with a threshold of z >= 3 (`input/expression_annotations.tsv`); `ezhip_tpm`, `ezhip_z` and `ezhip_expression_status` are reported. Tumors with neither are `none detected`, with the missing assay noted (`no RNA assay for EZHIP`, `no DNA assay`).

### Literature behind the lesion rules

Evidence behind the lesions that set or confirm a subtype. Strength is a summary of what the cited sources show, not a formal grading; the sources were found by web search and not all were read in full.

| Lesion | Used for | What the literature shows | Strength | Sources |
|---|---|---|---|---|
| CTNNB1 exon 3 | `CRANIO, ADAM` (overrides methylation), `MB, WNT` | Specific to adamantinomatous craniopharyngioma, mutually exclusive with BRAF V600E in most series; separate from methylation clusters that distinguish adamantinomatous from papillary | Strong | [Larkin 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4024131/), [Hölsken 2016](https://link.springer.com/article/10.1186/s40478-016-0287-6) |
| BRAF V600E | `CRANIO, PAP`, `PXA` split, `LCH` | Present in about 81-86% of papillary craniopharyngioma, absent from adamantinomatous in most series (rare co-occurrence reported); papillary craniopharyngioma is mainly adult, pediatric cases are rare | Strong | [Larkin 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4024131/), [pediatric case](https://pmc.ncbi.nlm.nih.gov/articles/PMC6341153/) |
| EGFR A289V / T / D | `DMG, H3 K28-altered, EGFR` and `-mutant, EGFR` | Extracellular-domain mutation: transforms immortalized astrocytes, increases invasion in mice, linked to shorter survival | Strong | [Frontiers in Oncology 2022](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.904383/full) |
| EGFR G598V | same | Increases glioma stem cell self-renewal, growth and survival; two GBM patients responded to dacomitinib | Moderate to strong | [Frontiers in Oncology 2022](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.904383/full) |
| EGFR R108K | **annotation only** (`egfr_variants`) | Increased ligand binding and anchorage-independent growth in NIH-3T3 cells, not glioma-specific; often co-occurs with A289, P596 and G598 changes; clinical significance "needs further elucidation" | Weaker | [PMKB](https://pmkb.weill.cornell.edu/variants/523), [PMKB](https://pmkb.weill.cornell.edu/variants/231) |
| EGFR exon 20 insertion / duplication | same as A289 | Included from the description of the `DMG_EGFR` methylation class; no DMG-specific source found in the search | Not established here | none found |
| EGFR amplification | same as A289 | Copy number above twice the ploidy, the `focal-cn-file-preparation` convention | Convention | `focal-cn-file-preparation` |
| EGFR alteration in DMG | the `DMG, H3 K28-altered, EGFR` concept | A subset of pediatric thalamic gliomas with a distinct methylation profile, H3K27me3 loss and frequent EGFR alteration; WHO CNS5 lists EGFR-mutant as a subtype of H3 K27-altered DMG alongside H3 K27-mutant and EZHIP-overexpressing | Moderate | [Jacques et al.](https://discovery-pp.ucl.ac.uk/10114996/3/Jacques_A%20subset%20of%20pediatric-type%20thalamic%20gliomas%20share%20a%20distinct%20DNA%20methylation%20profile%2C%20H3K27me3%20loss%20and%20frequent%20alteration%20of%20EGFR_AAM.pdf), [bioRxiv 2020](https://biorxiv.org/content/10.1101/2020.08.28.239160v1.full.pdf), [Tauziède-Espariat 2023](https://www.brainlife.org/fulltext/2023/Tauzi%C3%A8de-Espariat_A230618_BrainPathol.pdf), [CNS5 summary](https://test.ccga.io/index.php/CNS5:Diffuse_midline_glioma,_H3_K27-altered) |
| EZHIP overexpression (z >= 3 within RNA library type) | `h3k27_altered_mechanism` | Alternative route to H3K27me3 loss in H3 K27-altered DMG. The threshold is **derived from this cohort**, not the literature: z >= 3 captures 60% of PF-EPN-A tumors (the positive control) and 1% of H3 K27M-mutant DMG. A fixed 1 TPM cutoff was dropped because exome-capture libraries have about 15-fold higher EZHIP TPM than poly-A or stranded ones, and it called 25% of K27M-mutant DMG overexpressing | Cohort-derived | [Jacques et al.](https://discovery-pp.ucl.ac.uk/10114996/3/Jacques_A%20subset%20of%20pediatric-type%20thalamic%20gliomas%20share%20a%20distinct%20DNA%20methylation%20profile%2C%20H3K27me3%20loss%20and%20frequent%20alteration%20of%20EGFR_AAM.pdf) (mechanism only) |

Not yet sourced: the PB, ATRT, ETMR, EPN and MB lesion rules (SMARCB1, DICER1, C19MC, ZFTA/YAP1, PTCH1/SMO/SUFU, GLI2/MYC) were defined from the legacy modules and still need citations.

### Release columns

`06-compile-subtypes.R` writes, for the integrate step to add to `histologies.tsv`: `molecular_subtype_cns`, `subtype_source`, `subtype_alterations_confirmed`, `alteration_confirmed`, `subtype_review_flag` (a short coded reason, or blank) and `Notes`. `subtype_alterations_confirmed` lists the molecular alterations that support the subtype: mutations (`H3-3A p.K28M`, `CTNNB1 p.T41I`, `IDH1 p.R132H`), fusions (`ZFTA--RELA`, `EWSR1--FLI1`), copy-number changes (`1p+19q codeleted`, `EGFR amplification`, `SMARCB1 focal loss`) and expression evidence: `EZHIP overexpression` for labels that carry it, and `expression classifier agrees (MB, Group4)` for medulloblastoma tumors where the expression classifier matches the methylation group. It is blank when none was found or the tumor has no DNA/RNA assay, and `alteration_confirmed` is whether it is non-blank. `Notes` is generated from `subtype_source` with a new vocabulary (the legacy strings are not reused, since all subtypes are being redone): `Subtype from methylation`, `Subtype from lesion`, `Subtype from lesion, overriding methylation`, `Subtype by exclusion, no H3 or IDH variant`, `Subtype from expression classifier`, `Subtype from pathology free text diagnosis`, `Subtype from pathology diagnosis and primary site`, `Subtype from lesion shared with the tumor` and `Subtype from manual review of pathology report`. Free-text reviewer notes and report excerpts stay in `input/manual_subtype_review.tsv` and are not meant for the release.

`lesion_refinement_status` is `specific`, `none detected (DNA and RNA)`, `incomplete: no DNA`, `incomplete: no RNA` or `incomplete: no DNA or RNA` for the LGG, GNG, GNT and SEGA families.

`subtype_review_flag` values: `diagnosis discordant`, `classifier versions disagree`, `lesion contradicts methylation`, `H3 variant outside DMG/DHG`, `H3 variant in wildtype group`, `expression classifier discordant`, `free-text diagnosis differs from methylation` and `cell line differs from tumor` (several are joined with `;`). Cell lines also carry `tumor_subtype`, `tumor_group`, `tumor_concordance` and `tumor_alterations_shared`.

### Decisions

Decisions, their evidence, and open questions are logged in [DECISIONS.md](DECISIONS.md).

### Input data

The large release tables are read from their parquet copies in `data/parquet/` when they exist (`bash scripts/download-parquet.sh`, from `v15/parquet/`): the consensus MAF, the consensus gene-level copy number, `consensus_seg_with_status`, `fusion-annoFuse`, `fusion-dgd`, `sv-manta` and the collapsed TPM matrix. Without them the scripts fall back to the release TSV / RDS files (`utils/release-tables.R`) and give the same results, only slower: with the parquet copies scripts 02 to 07 run in about 20 seconds in total. `data/histologies-base.tsv` and the smaller files are read as before.

### Usage

```
Rscript --vanilla 01-build-methylation-calls.R
Rscript --vanilla 02-prepare-lesion-tables.R
Rscript --vanilla 03-validate-subtypes.R
Rscript --vanilla 04-lgg-lesions.R
Rscript --vanilla 05-validate-mb.R
Rscript --vanilla 06-compile-subtypes.R
Rscript --vanilla 07-compare-to-legacy.R
```
