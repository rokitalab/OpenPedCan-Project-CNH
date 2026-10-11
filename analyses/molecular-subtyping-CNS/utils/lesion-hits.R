#' Find tumors that meet each lesion rule
#'
#' @param rules Rules table: rule_id, OPC_molecular_subtype, lesion_type
#'   (snv, fusion, arm or cnv), gene, pattern.
#' @param snv SNV table with Hugo_Symbol, HGVSp_Short and
#'   Kids_First_Biospecimen_ID.
#' @param fusion Fusion table with FusionName and Kids_First_Biospecimen_ID.
#' @param arm_values Arm-level values, one column per arm (e.g. 1p, 19q) and
#'   Kids_First_Biospecimen_ID.
#' @param cnv Gene-level copy number with gene_symbol, status and
#'   Kids_First_Biospecimen_ID; a cnv rule matches `status` against `pattern`.
#' @param assay_ids Named list (snv, fusion, arm, cnv) of tables with match_id and
#'   the biospecimen id (`id`) of the assay each lesion type is read from.
#'
#' @return One row per match_id, rule and piece of evidence.
find_lesion_hits <- function(rules, snv, fusion, arm_values, cnv, assay_ids) {
  matches_for_rule <- function(type, gene, pattern) {
    switch(
      type,
      snv = snv %>%
        filter(str_detect(Hugo_Symbol, paste0("^(", gene, ")$")),
               str_detect(HGVSp_Short, pattern)) %>%
        transmute(id = Kids_First_Biospecimen_ID,
                  evidence = paste0(Hugo_Symbol, " ", HGVSp_Short)),
      fusion = fusion %>%
        filter(str_detect(FusionName, pattern)) %>%
        transmute(id = Kids_First_Biospecimen_ID, evidence = FusionName),
      cnv = cnv %>%
        filter(str_detect(gene_symbol, paste0("^(", gene, ")$")),
               str_detect(status, pattern)) %>%
        transmute(id = Kids_First_Biospecimen_ID,
                  evidence = paste0(gene_symbol, " ", status)),
      # met when every arm in `gene` (";" separated) is at or below the cutoff
      arm = {
        arms <- str_split(gene, ";")[[1]]
        arm_values %>%
          filter(if_all(all_of(arms), ~ !is.na(.x) & .x <= as.numeric(pattern))) %>%
          transmute(id = Kids_First_Biospecimen_ID,
                    evidence = paste0(paste(arms, collapse = "+"), " codeleted"))
      }
    )
  }

  pmap_dfr(rules, function(rule_id, OPC_molecular_subtype, lesion_type, gene,
                           pattern, ...) {
    matches_for_rule(lesion_type, gene, pattern) %>%
      inner_join(assay_ids[[lesion_type]], by = "id",
                 relationship = "many-to-many") %>%
      transmute(match_id, rule_id, OPC_molecular_subtype, evidence)
  }) %>%
    distinct()
}
