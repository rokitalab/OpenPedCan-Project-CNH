# Archived molecular subtyping modules

These modules were replaced by [`molecular-subtyping`](../molecular-subtyping), which sets every molecular subtype from methylation first and confirms it with lesions (see its `README.md` and `DECISIONS.md`). They are kept with their tracked results for reference and are not run in CI or by `scripts/run-for-subtyping.sh`.

| Archived module | Replaced by |
|---|---|
| `molecular-subtyping-HGG` | `molecular-subtyping` (`01`, `02`, `03`) |
| `molecular-subtyping-LGAT` | `molecular-subtyping/04-lgg-lesions.R` |
| `molecular-subtyping-EPN`, `-ATRT`, `-PB`, `-CRANIO`, `-Rare-CNS`, `-neurocytoma`, `-embryonal` | `molecular-subtyping` (`01`, `03`) |
| `molecular-subtyping-NBL`, `-EWS` | `molecular-subtyping/06-subtype-solid-tumors.R` |
| `molecular-subtyping-pathology` | `molecular-subtyping` for the subtypes; its compiled table is the frozen fallback read by `../molecular-subtyping-integrate/01-integrate-subtyping.Rmd` |

The `fusion-summary` module was archived with them (the LGG/HGG kinase fusions with reciprocal fusions, and the Ewing sarcoma fusions of interest, from the annoFuse putative oncogenic and DGD fusions; ported to `molecular-subtyping/utils/fusion-calls.R`, run by `molecular-subtyping/02-prepare-lesion-tables.R`, with identical calls). `fusion_filtering` stays active: it produces `fusion-putative-oncogenic.tsv` for the release, before subtyping.

`molecular-subtyping-MB` and `molecular-subtyping-chordoma` are not archived as folders: their code was folded into `molecular-subtyping/05-validate-mb.R` (expression classifier, SHH subtype) and `molecular-subtyping/07-compile-subtypes.R` (chordoma SMARCB1 status). The MB class comparison notebook and SHH methylation UMAP were dropped; they are in the git history.

The archived scripts were not changed, so their paths (`analyses/molecular-subtyping-<name>`) are no longer valid and they do not run from here without editing them. Their last results are read by:

- `../molecular-subtyping-integrate/01-integrate-subtyping.Rmd`: `molecular-subtyping-pathology/results/compiled_molecular_subtypes_with_clinical_pathology_feedback_and_report_info.tsv`, kept for tumors that `molecular-subtyping` does not label;
- `../molecular-subtyping-integrate/00-derive-integrated-diagnosis.R` and `../molecular-subtyping/08-compare-to-legacy.R`: `compiled_molecular_subtypes.tsv`, for comparison.

Tumors added in a later release are not in these frozen tables, so they have no legacy label; only the base diagnosis applies to them.
