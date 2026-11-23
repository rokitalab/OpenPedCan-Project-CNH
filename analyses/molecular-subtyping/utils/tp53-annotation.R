#' TP53 annotation per match_id, from analyses/tp53_nf1_score
#'
#' Not part of the subtype. Reports that module's altered status and classifier
#' score, and counts hits the way it does: each SNV/indel is one hit and a CNV
#' loss or SV is one more (CNV and SV are the same event). Two or more hits
#' suggest both alleles are affected. DNA and RNA biospecimens of a tumor each
#' have a row, so values are collapsed to one row per match_id.
read_tp53_annotation <- function(root_dir) {
  read_tsv(
    file.path(root_dir, "analyses", "tp53_nf1_score", "results",
              "tp53_altered_status.tsv"),
    col_types = cols(.default = col_character()), na = c("", "NA")
  ) %>%
    filter(!is.na(match_id)) %>%
    mutate(across(c(SNV_indel_counts, CNV_loss_counts, SV_counts, hotspot,
                    tp53_score), as.numeric),
           li_fraumeni = str_detect(replace_na(cancer_predispositions, ""),
                                    "Li-Fraumeni syndrome")) %>%
    group_by(match_id) %>%
    summarise(
      tp53_altered = case_when(
        any(tp53_altered == "loss") ~ "loss",
        any(tp53_altered == "activated") ~ "activated",
        any(tp53_altered == "Other") ~ "Other",
        TRUE ~ NA_character_
      ),
      tp53_classifier_score = suppressWarnings(max(tp53_score, na.rm = TRUE)),
      snv_hits = suppressWarnings(max(SNV_indel_counts, na.rm = TRUE)),
      cnv_sv_hit = suppressWarnings(max(CNV_loss_counts, SV_counts,
                                        na.rm = TRUE)) >= 1,
      tp53_hotspot = suppressWarnings(max(hotspot, na.rm = TRUE)) == 1,
      tp53_li_fraumeni = any(li_fraumeni),
      tp53_snv_evidence = paste(na.omit(unique(HGVSp_Short[HGVSp_Short != "NA"])),
                                collapse = "; "),
      tp53_cnv_evidence = paste(na.omit(unique(CNV_loss_evidence)),
                                collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      across(c(tp53_classifier_score, snv_hits),
             ~ if_else(is.infinite(.x), NA_real_, .x)),
      snv_hits = replace_na(snv_hits, 0),
      tp53_hotspot = replace_na(tp53_hotspot, FALSE),
      cnv_sv_hit = replace_na(cnv_sv_hit, FALSE),
      tp53_n_hits = snv_hits + as.integer(cnv_sv_hit),
      tp53_two_hits = tp53_n_hits >= 2
    ) %>%
    select(match_id, tp53_altered, tp53_classifier_score, tp53_n_hits,
           tp53_two_hits, tp53_hotspot, tp53_li_fraumeni, tp53_snv_evidence,
           tp53_cnv_evidence)
}
