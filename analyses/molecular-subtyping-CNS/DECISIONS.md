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
| NIH v3.1 subclass names missing from the crosswalk are aliased to the DKFZ v12.8 row they match (19 classes, 93 to 100% agreement in arrays where both are high confidence); `PA PF VAR` is left unresolved: the NIH classifier has no publication or class table to define it, and in the data only 9 of 41 arrays are `PA_INF` by DKFZ v12.8 (20 `CTRL_CBM`, 11 `MB_SHH_2`, 1 DNET) with diagnoses that include medulloblastoma, ATRT and ependymoma | after dropping DKFZ v12 / NIH v2, 490 high-confidence NIH v3.1 arrays (28 classes) had no resolvable class and 79 had no call at all; 60 arrays regained a call, 42 array-calls in 10 classes remain unresolved | `input/nih_v3.1_class_aliases.tsv`, `01` |
| **Low-confidence DKFZ calls are uninformative** and are not used as evidence against an NIH call | reviewer decision | applies to all `methyl_source` handling |

| Provenance goes in `Notes` plus `subtype_source` and `subtype_review_flag`, with a fresh vocabulary (legacy `Notes` strings are not reused) | all subtypes are being redone, so old strings would be misleading; report wording and reviewer interpretation do not belong in a release file | `05` |

## Evidence and rules

| Decision | Reason / evidence | Where |
|---|---|---|
| High-confidence methylation: DKFZ subclass >= 0.8, else NIH superfamily and class >= 0.9; precedence DKFZ v12.8, then NIH v3.1. **Earlier classifier versions (DKFZ v11 / v12, NIH v2) are not used**; classifier scores below 0.75 are not considered as evidence in review tables | same threshold rule as Rare CNS; reviewer decision (2026-10-12). 41 arrays lost a call and 20 tumors changed (12 became unsubtyped, the rest take a lesion fallback or the lesion-negative HGG label) | `01` |
| Lesions confirm methylation; they never override it, **except in craniopharyngioma** | CTNNB1 exon 3 and BRAF V600E are mutually exclusive; 11 methylation-PAP calls all had CTNNB1 and no BRAF V600E; 9 had pathology reports saying adamantinomatous, 2 (7316-3070, 7316-6396) say "craniopharyngioma" only, none contradicted | `lesion_overrides_methylation` rules |
| NIH v3.1 is not better than DKFZ for craniopharyngioma | among 65 CTNNB1-mutant tumors: DKFZ 73% correct (8 wrong), NIH 73% (11 wrong, same 8 plus 3); NIH scores on the errors were higher | analysis, not code |
| A contradicted and weak (NIH-only, or DKFZ < 0.85) call is set aside | strong contradicted calls are kept and flagged | `03`, `methylation_downgraded.tsv` |
| `Oligo, IDH` needs IDH plus 1p/19q codeletion; oligodendroglioma diagnosis without codeletion is held | requested | `03` |
| **Caveat:** 1p/19q uses GISTIC arm values, which are not ploidy-corrected | seen in 7316-4220 (ploidy 3: 1p shows +1.33 in GISTIC but is a loss against ploidy) | open; recheck with segments |
| EGFR activating: amplification, A289 and G598 missense, exon 20 insertion / duplication; **R108K annotation only** | A289V/G598V have the functional evidence; R108K is weaker and often co-occurs | `input/annotation_genes.tsv`, README literature table |
| A DMG with a DNA assay, an EGFR lesion and no H3 variant is EGFR-altered, not K27M-mutant, whatever the free text says | free text can be inaccurate for mutation | `relabel_egfr_only()` |
| EZHIP overexpression is a z-score of log2(TPM+1) within RNA library type, z >= 3 | a fixed 1 TPM called 25% of K27M-mutant DMG overexpressing; exome-capture libraries read about 15x higher | `02`, `expression_annotations.tsv` |
| EZHIP is **not assessable in exome-capture RNA libraries** (`not assessable (exome capture library)`, 46 tumors); LIN28A is called in every library | exome capture reads about 2.2 TPM for EZHIP in every tumor type (20 to 50 times the other libraries), 0 of 150 reference tumors reach z >= 3 and there is one PF-A control (3.6 TPM), so a threshold cannot be calibrated; ETMR LIN28A (10 TPM vs 0) is not affected. No label changed | `exclude_library_groups` in `input/expression_annotations.tsv`, `annotate_expression()` |
| DMG from the free-text diagnosis when no other source applies, as plain `DMG, H3 K28` | 21 lesion-negative tumors had DMG free text and contradicted `HGG, H3 wildtype` | `03` |
| The control-methylation lesion fallback is **not** applied, with two exceptions: an H3 K27/K28 M or I variant (DMG) and the cranio free-text rule | reviewer decision; a control array only says low purity, so the specific variant or explicit histology is used (4 H3 K28M tumors became `DMG, H3 K28-mutant`, 5 cranio tumors `CRANIO, ADAM`). HGG with a control array stay unsubtyped | `03`, `h3_variant_review.tsv` |

| LGG lesion logic is ported into this module (`04-lgg-lesions.R`) instead of read from the LGAT module | the evidence behind the lesion-level labels was not visible here, so `alteration_confirmed` was blank for 700 of 721 LGG-family tumors; the port reproduces the LGAT labels for 616 of 617 shared tumors, which lets the LGAT module be archived | `04`, `input/lgg_component_rules.tsv` |

| Embryonal rules (TTYH1 fusion, C19MC amplification, LIN28A z >= 3) are ported; `CNS Embryonal, NOS` is kept for embryonal tumors with no other source | the legacy module defined `ETMR, NOS` by lesions and diagnosis, which the comparison first missed because "NOS" was treated as unassigned | `03`, `input/embryonal_path_dx_strings.json` |
| C19MC amplification requires a focal segment overlapping the cluster with segment mean > 1, not any chr19 gain | the legacy "Yes" included 9 tumors with chr19 segments elsewhere (2 Mb, 10 Mb, 54.2 Mb) | `03` |
| `RARE_MN1_BEND2`, `RARE_MN1_CXXC5` and `RARE_FOXR2` fusions subtype tumors without methylation (FOXR2 only with a concordant diagnosis) | the fusion is specific; 7316-3166 was `CNS Embryonal, NOS` with MN1--BEND2. 7316-1473 keeps `HGG, RTK` by methylation (ZC4H2--FOXR2 is not enough against RTK1B 1.0) | `input/subtype_lesion_rules.tsv` |
| `Rare CNS, DICER1 SARC` without methylation needs an RNase IIIb hotspot missense (E1705, D1709, D1713, G1809, D1810, E1813) and a compatible diagnosis; truncating or splice variants confirm only; a hotspot in pineoblastoma is held | a single nonsense variant called 7316-8965 (a mismatch-repair-deficient GBM) and two pineoblastomas; DICER1 in pineoblastoma is `PB, miRNA processing` | `RARE_DICER1_HOTSPOT`, `03` |
| PXA without methylation: BRAF V600E plus a PXA term in the free text | BRAF V600E alone is not specific; 4 `HGG, PXA BRAF V600E` and 9 `LGG, PXA BRAF V600E` tumors | `PXA_BRAF_V600E`, `03` |
| `HGG, H3 wildtype` by exclusion is only withheld by an H3 or IDH rule hit, not any rule hit | rule hits for TP53, NF1, RB1 and so on are not specific; this left 76 HGG/DIPG tumors unsubtyped (61 became 138 lesion-negative) | `03` |
| `CRANIO, ADAM` from "adamantinomatous" in the free text when methylation, CTNNB1 and BRAF gave nothing; a high-confidence conflicting methylation call stands and is flagged | 12 tumors; conflicts 7316-5310 (PITUI) and 7316-6135 (`LGG, MAPK`) | `03`, `free_text_conflicts.tsv` |
| Neurocytoma ported: diagnosis plus primary site (Ventricles = `CNC`, else `EVN`); CN and EVNCYT methylation take precedence | the legacy module | `03` |
| Chordoma label is `CHDM`; spinal subependymoma stays `EPN, SP-SE` (no A/B split) | continuity with legacy labels and the `integrated_diagnosis` strings | crosswalk |
| DKFZ `PA_CORT` (NIH v3.1 `PA ST`) is `LGG, PA ST` with the prefix from the primary diagnosis (`GNG` for ganglioglioma, `GNT` for DNET or glial-neuronal, else `LGG`); the DLGNT classes stay `GNT, MAPK` | the crosswalk had `GNT, MAPK` for a pilocytic astrocytoma class; "ST" pairs with NIH `PA ST` / `PA MID`, "cortical" is the DKFZ wording | crosswalk, `03` |
| `molecular_subtype_cns` is the lesion-level label for the LGG families: generic MAPK buckets are replaced by the lesion term and `PA ST` gets it appended; **defined entities (DNET, SEGA, DGONC, DIG/DIA, papillary, myxoid) are not split by lesion**; classes that already are a lesion term and DLGNT are left alone | matches the legacy LGAT detail without fragmenting DNET (FGFR1 is typical of DNET, not a subtype); 203 tumors changed, 29 reverted for the defined entities | `06` |
| No `To be classified` labels; `lesion_refinement_status` (`specific`, `none detected (DNA and RNA)`, `incomplete: no DNA / RNA`) says what the lesion module could see | an incomplete tumor is still to be determined, not lesion-negative; a string in the label would pollute `integrated_diagnosis` and `cancer_group` | `06` |
| SEGA is not split by lesion: the TSC1/TSC2 variants that define it are germline and are not in the somatic calls | reviewer decision | `06` |
| EGFR amplification uses the WGS consensus segments plus CNVkit segments for biospecimens the consensus does not cover (WXS tumors and some WGS) | 73 tumors had `no CNV data` for EGFR only because they are not in `consensus_seg_with_status` | `03` |
| The large release tables are read from the v15 parquet copies in `data/parquet/` when present, with a TSV / RDS fallback and identical results | scripts 02 to 07 run in about 20 s | `utils/release-tables.R` |
| `integrated_diagnosis`, `broad_histology` and `short_histology` are derived by rule from the label in `molecular-subtyping-integrate`; DLGNT and the other pathology free-text diagnoses are applied there | the diagnosis follows the subtype; the first comma token stays the `cancer_group` | `molecular-subtyping-integrate/00-derive-integrated-diagnosis.R` |

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
| Derived cell lines | own match_id, compared with the tumor | subtyped from their own methylation, lesions and RNA; `tumor_subtype`, `tumor_concordance` and `tumor_alterations_shared` compare with the solid-tissue tumors of the same sample (9 `different`, 13 `same group`, 49 `same`, 26 not comparable); a line with no subtype of its own that shares a confirming lesion with its tumor takes the tumor subtype (7316-440 `DMG, H3 K28-altered, EGFR`); a DICER1 hotspot subtypes 7316-1746 lines | `06`, `cell_line_tumor_concordance.tsv` |

| Diffuse leptomeningeal glioneuronal tumor (DLGNT) | not in this module | the legacy glioneuronal free-text step called 5 samples DLGNT (7316-2208, 7316-10000, C4106232-Tumor, ET_6R1D00B5, ET_7T785YTJ). The CNS module keeps their methylation label (7316-2208 is `GNT, MAPK`); the DLGNT diagnosis is applied in `molecular-subtyping-integrate`, which keeps the free-text steps for diagnoses that are not molecular subtypes |
| DKFZ `PA_CORT` (NIH v3.1 `PA ST`) | `LGG, PA ST`, `GNG, PA ST` or `GNT, PA ST` | supratentorial pilocytic astrocytoma by both classifiers; the crosswalk had `GNT, MAPK`. The family prefix follows the primary diagnosis, as for PXA: a ganglioglioma diagnosis gives `GNG`, a DNET or glial-neuronal tumor diagnosis `GNT`, anything else `LGG` (29 tumors: 9 GNG, 7 GNT, 13 LGG). The lesion term from the LGG lesion module is part of the label (`GNG, PA ST, BRAF V600E`; 25 of 29), and `lesion_refinement_status` records `none detected` or `incomplete: no DNA / RNA` instead of a `To be classified` label. The remaining `GNT, MAPK` tumors are DLGNT_1 and DLGNT_2. 7316-1106 (HGG diagnosis, BRAF V600E, PA_CORT 0.98, NIH PA ST 0.99) is `LGG, PA ST`, with no PXA class and no CDKN2A/B loss |
| LGG-family classes with a lesion term | lesion-level `molecular_subtype_cns` | generic MAPK buckets (PA_MID, PA_INF, GG) are replaced by the lesion term and PA ST gets it appended; defined entities (DNET, SEGA, DGONC, DIG/DIA, papillary/myxoid GNT) stay one subtype with the lesion in `lesion_refinement`, so a DNET is not split by its FGFR / RTK / NF1 lesion (e.g. 83 `LGG, BRAF/MAPK` to `LGG, KIAA1549-BRAF`). Matches the legacy LGAT level of detail; the class stays in `methyl_class` |
## Possible methylation sample mix-ups

`input/possible_methylation_mixups.tsv` lists tumors whose methylation class conflicts with the DNA, RNA or diagnosis in a way that a sample mix-up on the methylation side could explain (14 tumors: 7316-9842 and the two other H3 K27M tumors with an EPN PF-A class, 7316-5481, 7316-1856, 7316-1473, 7316-5310 and the 7 PA_MID tumors with a craniopharyngioma diagnosis), with the confidence of every array. It is a tracking list, not an input: the methylation call is kept and the existing flags apply until the specimens are checked. 7316-9842 is being investigated with CBTN. 7316-6755 was listed at first and removed: two of its three arrays are DMG_K27 at 1.0 (DKFZ v12.8), and the one GG array is NIH-only; its label is `DMG, H3 K28-mutant`.

| PNOC003 clinical review (38 rows, `molecular-subtyping-pathology/input/pnoc003-clinical-review.tsv`) | the K28M tumors agree; the H3 WT / no-report tumors are DMG by cohort (next row) | all 36 rows with a tumor in the histologies agree for K28M (58 tumors are `DMG`); PT_FFHQ1B9R has no tumors because it is the same patient as PT_V1HNAC2Q (merged upstream; 7316-3240 is `DMG, H3 K28-mutant`) |
| CI and the release script | `run-molecular-subtyping-CNS.sh`, a CI matrix job and a step in `scripts/run-for-subtyping.sh` | the module runs on the `testing` release (TSV fallback, 20 s, no parquet files) with exit code 0; the other modules' results it reads (`mb-classified.rds`, `chordoma_smarcb1_status.tsv`, `tp53_altered_status.tsv`) are tracked in the repository. The compile and integrate steps do not read its table yet |

| PNOC003 cohort tumors that were not `DMG` (8 tumors, 3 participants: PT_KBFM551M, PT_C5FKRB1P, PT_8P368R5B; PT_0MXPTTM3 is below) | `DMG, H3 wildtype` (manual review) | the cohort is diffuse intrinsic pontine glioma, so DMG (reviewer decision, 2026-10-12); the reports give H3 wildtype and there is no H3 K27 variant in the DNA. None has an EGFR alteration; EZHIP is overexpressed in 1 (7316-3214-A07083, so `DMG, H3 K28-altered, EZHIP`), not overexpressed in 2 and not assessable in 5 (no RNA) |
| PT_0MXPTTM3 (7316-3220, A09410) | `DMG, H3 wildtype` (manual review) | PNOC003 cohort, DIPG diagnosis (reviewer decision). The molecular profile is that of an IDH-mutant astrocytoma (adult, 25.8 years; IDH1 p.R132H VAF 0.36, ATRX p.E2279* VAF 0.86, TP53 R273H and R110L; no H3 variant; no EGFR alteration; EZHIP not overexpressed in A09410, no RNA for 7316-3220), so the IDH1 variant is in the manual note, not the label; `HGG, IDH` is the automated call and stays in `automated_subtype` |
## Pending

- 7316-9842 (being investigated with CBTN) and the other tumors in `input/possible_methylation_mixups.tsv`: 7316-506, -5511, -5481, -1856, -1473, -5310 and 7 PA_MID craniopharyngioma-diagnosis tumors. The H3 K28M tumors with a control array (6) are subtyped now; `h3_variant_review.tsv` lists the 3 tumors with a non-DMG class (7316-506, -5511, -9842).
- The KIAA1549::BRAF SV check for the PA_MID craniopharyngioma-diagnosis tumors needs redoing against `sv-manta`.
- 7316-1106: pathology review `HGG, PXA` vs methylation `LGG, PA ST` (PA_CORT 0.98, NIH v3.1 PA ST 0.99). Lesions do not decide it: BRAF V600E (26% VAF) fits both, no CDKN2A/B homozygous deletion (copy number 3), no fusion, no H3/IDH/TP53 hit. Kept as the methylation call; an HGG diagnosis with a low-grade class is for pathology review (possible anaplastic PA).
- Seed the manual-review table from the remaining legacy path-review overrides (7316-9049, C4116318-Tumor, 7316-272, 7316-4220, the cranio and chordoma rows are done).
- Citations for the PB, ATRT, ETMR, EPN and MB lesion rules.
- Copy-number rules (NF2, NF1, RB1, VHL, CDKN2A/B for the NST, HMB and similar families) are not ported.
- The 1p/19q ploidy-aware recheck.
- 7316-7965 cell line (`MB, SHH` by the expression classifier vs `MB, Group3`) and 7316-1763 (tumor K28M not in the line): worth raising with whoever derived the lines.
- Downstream: the integrate Rmd still reads the legacy compiled table; solid tumors (NBL, EWS) go to `molecular-subtyping-solid-tumors`; the legacy modules are archived after integrate reads the new tables.
- Split into one script per PR.
