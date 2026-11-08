# Review lists

Hand-curated lists that no script reads or writes; they record what was checked
and what is still open. The tables the module uses are in `../input/`, and the
tables it makes are in `../results/`.

- `lesion_gap_review.tsv`: lesions that could be added as rules, from the review
  of the rules against the literature, with the basis of each entry.
- `pptc_nbl_mycn_template.tsv`: the PPTC neuroblastoma models with the MYCN
  status inferred from expression, to fill in with the clinical status. The
  statuses given so far are in `../input/nbl_pptc_clinical_mycn_status.tsv`.
- `methylation_diagnosis_discordance_notes.tsv`: the notes on 21 tumors whose methylation class is discordant with the diagnosis but whose genomics does not contradict the class (no DNA or RNA, low tumor purity, or no lesion either way), previously tiered `possible` or `kept` mismatches. They are not mismatches: only a methylation class that the genomics contradicts (or a class mutually exclusive with a lesion) is flagged and its array excluded; these tumors keep the methylation subtype with the diagnosis-discordance flag.
