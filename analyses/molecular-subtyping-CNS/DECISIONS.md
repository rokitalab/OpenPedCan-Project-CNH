# Decision log: molecular-subtyping-CNS

Decisions made while building this module, with the evidence and where they live. Entries are dated (2026-10) and refer to data release v15. Add new decisions at the bottom.

## Calls and labels

| Decision | Reason / evidence | Where |
|---|---|---|
| Control tissues and `INFLAM_ENV` are `Control` (low tumor purity), never subtyped | They are not tumor classes; a no-call hid the purity information | crosswalk `OPC_molecular_subtype`, `methyl_control` |
| `GG` NIH v3 code corrected from `CONTR HEMI` to `GG` | `CONTR HEMI` is a control tissue and would have called controls `GNG, MAPK` | crosswalk; control-code rule in `01` |
| 45 unassigned methylation classes named from the classifier taxonomy and WHO CNS5 | **2022-2026 literature was not verified** (no web access at the time); names need a pathology read | `input/opc_subtype_additions.tsv` |
| `HGG, `, `LGG, `, `GNT, ` prefixes on all family labels | consistency across modules | crosswalk, `04-compare-to-legacy` renames |
| `-mutant` only when the H3 variant is confirmed by DNA; plain `DMG, H3 K28` / `DHG, H3 G35` for methylation or diagnosis only | free text can say "mutant" without a variant (7316-9049, C4116318-Tumor report said "altered") | `03` |
| `DMG_EGFR` is its own subtype, `DMG, H3 K28-altered, EGFR` | WHO CNS5 lists EGFR-mutant H3 K27-altered DMG; the legacy label `DMG, H3 K28` misdescribed it | crosswalk, `utils/gene-annotation.R` |
| EZHIP overexpression is added to DMG labels | requested; EZHIP is an alternative route to H3 K27 alteration | `add_ezhip_label()` |
| PXA split by BRAF V600E (issue 113); HGG / LGG prefix from the diagnosis, then free text | path diagnosis decides, per review | `03` |
| `broad_histology` is not used as evidence | it is a manual lookup of `pathology_diagnosis` upstream, not independent | README |
| **Low-confidence DKFZ calls are uninformative** and are not used as evidence against an NIH call | reviewer decision | applies to all `methyl_source` handling |

| Provenance goes in `Notes` plus `subtype_source` and `subtype_review_flag`, with a fresh vocabulary (legacy `Notes` strings are not reused) | all subtypes are being redone, so old strings would be misleading; report wording and reviewer interpretation do not belong in a release file | `05` |

## Evidence and rules

| Decision | Reason / evidence | Where |
|---|---|---|
| High-confidence methylation: DKFZ subclass >= 0.8, else NIH superfamily and class >= 0.9; precedence DKFZ v12.8, NIH v3.1, DKFZ v12, NIH v2 | same rule as Rare CNS | `01` |
| Lesions confirm methylation; they never override it, **except in craniopharyngioma** | CTNNB1 exon 3 and BRAF V600E are mutually exclusive; 11 methylation-PAP calls all had CTNNB1 and no BRAF V600E; 9 had pathology reports saying adamantinomatous, 2 (7316-3070, 7316-6396) say "craniopharyngioma" only, none contradicted | `lesion_overrides_methylation` rules |
| NIH v3.1 is not better than DKFZ for craniopharyngioma | among 65 CTNNB1-mutant tumors: DKFZ 73% correct (8 wrong), NIH 73% (11 wrong, same 8 plus 3); NIH scores on the errors were higher | analysis, not code |
| A contradicted and weak (NIH-only, or DKFZ < 0.85) call is set aside | strong contradicted calls are kept and flagged | `03`, `methylation_downgraded.tsv` |
| `Oligo, IDH` needs IDH plus 1p/19q codeletion; oligodendroglioma diagnosis without codeletion is held | requested | `03` |
| **Caveat:** 1p/19q uses GISTIC arm values, which are not ploidy-corrected | seen in 7316-4220 (ploidy 3: 1p shows +1.33 in GISTIC but is a loss against ploidy) | open; recheck with segments |
| EGFR activating: amplification, A289 and G598 missense, exon 20 insertion / duplication; **R108K annotation only** | A289V/G598V have the functional evidence; R108K is weaker and often co-occurs | `input/annotation_genes.tsv`, README literature table |
| A DMG with a DNA assay, an EGFR lesion and no H3 variant is EGFR-altered, not K27M-mutant, whatever the free text says | free text can be inaccurate for mutation | `relabel_egfr_only()` |
| EZHIP overexpression is a z-score of log2(TPM+1) within RNA library type, z >= 3 | a fixed 1 TPM called 25% of K27M-mutant DMG overexpressing; exome-capture libraries read about 15x higher | `02`, `expression_annotations.tsv` |
| **Open:** exome-capture RNA is not reliable for EZHIP (no PF-A control in that library) | proposed "not assessable", not yet applied | `02` |
| DMG from the free-text diagnosis when no other source applies, as plain `DMG, H3 K28` | 21 lesion-negative tumors had DMG free text and contradicted `HGG, H3 wildtype` | `03` |
| The control-methylation lesion fallback is **not** applied | reviewer decision; 10 H3 K28M tumors with a control array stay unsubtyped and are listed | `h3_variant_review.tsv` |

| LGG lesion logic is ported into this module (`04-lgg-lesions.R`) instead of read from the LGAT module | the evidence behind the lesion-level labels was not visible here, so `alteration_confirmed` was blank for 700 of 721 LGG-family tumors; the port reproduces the LGAT labels for 616 of 617 shared tumors, which lets the LGAT module be archived | `04`, `input/lgg_component_rules.tsv` |

| Embryonal rules (TTYH1 fusion, C19MC amplification, LIN28A z >= 3) are ported; `CNS Embryonal, NOS` is kept for embryonal tumors with no other source | the legacy module defined `ETMR, NOS` by lesions and diagnosis, which the comparison first missed because "NOS" was treated as unassigned | `03`, `input/embryonal_path_dx_strings.json` |
| C19MC amplification requires a focal segment overlapping the cluster with segment mean > 1, not any chr19 gain | the legacy "Yes" included 9 tumors with chr19 segments elsewhere (2 Mb, 10 Mb, 54.2 Mb) | `03` |

## Reviewed samples

| Sample | Decision | Basis |
|---|---|---|
| 7316-9042, -3652, -5220, -3533, -8972, -1748, -7514, -323, -4551 | `CRANIO, ADAM` | pathology: adamantinomatous (manual review) |
| 7316-3070, 7316-6396 | `CRANIO, ADAM` by CTNNB1 | reports state craniopharyngioma only; not contradicted; not in the manual table |
| 7316-9049 | `DMG, H3 K28-altered, EGFR` | report: H3 K27-altered DMG; EGFR p.A767_V769dup in DNA; no H3 variant |
| C4116318-Tumor | `DMG, H3 K28-altered, EGFR` | report: H3 K27-altered; EGFR p.G598V; the histologies free text says "mutant" |
| 7316-4220 | `HGG, RTK` | report: "Histone H3 negative by IHC"; **H3K27me3 not stained**; reviewer interpretation: not H3 K27-altered DMG; DKFZ RTK1 class 0.63 and PDGFRA K385I support RTK; TERT c.-146C>T |
| 7316-9842 | `HGG, IDH` kept, **flagged** | methylation `A_IDH_HG` 0.99 (2 arrays) vs WGS H3-3A K28M with no IDH variant and pathology DMG H3 K27M; CBTN contacted |
| 7316-3240 | `DMG, H3 K28-mutant` | the missing subtype in the v15 comparison was a participant mismatch fixed in `opc-histologies-qc` v16 (ticket 5) |
| 7316-1856 | `Rare CNS, EWS` | pathology: Ewing sarcoma metastatic to the CNS; MB_G34_VIII methylation likely a mismatched array (manual review) |
| 7316-3166 | `Rare CNS, MN1::BEND2` | MN1--BEND2 fusion, no methylation; MN1 fusion rules now subtype tumors without methylation |
| 7316-7091, -6484, -4698, C4110537-Tumor | `HGG, PXA BRAF V600E` | BRAF V600E plus a PXA free-text diagnosis, no methylation; 9 more tumors with LGG diagnosis and PXA free text become `LGG, PXA BRAF V600E` |
| 7316-2784, -3158, -3521, -9843 | `DMG, H3 K28-mutant` | H3-3A K28M with a control / INFLAM_ENV methylation call; the H3 K27/K28 variant is specific, so the control exception applies to DMG only |
| 7316-921, 7316-3632, 7316-4062 | `CHDM, conventional` / `CHDM, poorly differentiated` | chordoma pathology review (legacy ticket 608); methylation calls differ or are absent. The module label is `CHDM` again (was `CHORDM`) |
| 12 PBTA tumors with "adamantinomatous" in the free text | `CRANIO, ADAM` | free-text rule where no methylation, CTNNB1 or BRAF rule gave a subtype; includes 5 control-methylation arrays. High-confidence conflicting methylation stands and is flagged (7316-5310 PITUI, 7316-6135 `LGG, MAPK`) |
| 4 neurocytoma tumors (7316-280, C4117794-Tumor, 7316-10377, 7316-6779) | `CNC` / `EVN` | ported from molecular-subtyping-neurocytoma: diagnosis, with primary site Ventricles = CNC |
| Derived cell lines | own match_id | subtyped from their own methylation and lesions; `parent_subtype` and `parent_concordance` compare with the solid-tissue tumors of the same sample; `different` is flagged (16) in `results/cell_line_parent_concordance.tsv` |

## Pending

- 7316-9842 and 14 other H3-variant tumors: CBTN, `h3_variant_review.tsv`.
- 7316-5481 (ependymoma diagnosis, NIH-only `DMG K27`, EZHIP z 2.2): no H3K27me3 result.
- `PA_MID` craniopharyngioma-diagnosis tumors (7): possible misdiagnosis, not yet reviewed.
- Seed the manual-review table from the 98 legacy overrides (7316-9049 and C4116318-Tumor are done).
- Citations for the PB, ATRT, ETMR, EPN and MB lesion rules.
