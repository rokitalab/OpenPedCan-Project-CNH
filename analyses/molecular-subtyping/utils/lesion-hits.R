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
#' @param sv Structural variants (Manta) with Kids_First_Biospecimen_ID, Gene_name
#'   and sv_call (the ALT type and the exon location, e.g.
#'   "<DUP:TANDEM> exon15-exon15") and sv_length; an sv rule matches the gene and
#'   `pattern` against sv_call.
#' @param expr Expression z-scores (gene, Kids_First_Biospecimen_ID, library_group,
#'   z); an expr rule matches the gene with z at or above `pattern` (library
#'   groups with no usable z-score, such as exome capture, are skipped).
#' @param germline Germline evidence (match_id, gene, evidence) from a cancer
#'   predisposition or a pathogenic germline variant; a germline rule matches the
#'   gene and needs no tumor assay.
#' @param assay_ids Named list (snv, fusion, arm, cnv) of tables with match_id and
#'   the biospecimen id (`id`) of the assay each lesion type is read from.
#'
#' @return One row per match_id, rule and piece of evidence.
find_lesion_hits <- function(rules, snv, fusion, arm_values, cnv, assay_ids, sv = NULL, expr = NULL,
                            germline = NULL) {
  if (is.null(germline)) {
    germline <- tibble(match_id = character(), gene = character(),
                       evidence = character())
  }
  if (is.null(expr)) {
    expr <- tibble(gene = character(), Kids_First_Biospecimen_ID = character(),
                   library_group = character(), z = double())
  }
  if (is.null(sv)) {
    sv <- tibble(Kids_First_Biospecimen_ID = character(), Gene_name = character(),
                 sv_call = character(), sv_length = character())
  }
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
      sv = sv %>%
        filter(str_detect(Gene_name, paste0("(^|;)(", gene, ")(;|$)")),
               str_detect(sv_call, pattern)) %>%
        transmute(id = Kids_First_Biospecimen_ID,
                  evidence = paste0(Gene_name, " ", sv_call,
                                    if_else(is.na(sv_length), "",
                                            paste0(" (", sv_length, " bp)")))),
      germline = germline %>%
        filter(str_detect(.data$gene, paste0("^(", !!gene, ")$")),
               str_detect(evidence, pattern)) %>%
        transmute(id = match_id, evidence = evidence),
      expr = expr %>%
        filter(str_detect(.data$gene, paste0("^(", !!gene, ")$")),
               !library_group %in% "exome capture", z >= as.numeric(pattern)) %>%
        transmute(id = Kids_First_Biospecimen_ID,
                  evidence = paste0(.data$gene, " overexpression (z ", round(z, 1), ")")),
      # met when every arm in `gene` (";" separated) is at or below the cutoff;
      # both arms of one chromosome (6p;6q) are a monosomy, arms of different
      # chromosomes (1p;19q) a codeletion
      arm = {
        arms <- str_split(gene, ";")[[1]]
        chromosomes <- unique(str_remove(arms, "[pq]$"))
        label <- if (length(arms) == 2 && length(chromosomes) == 1) {
          paste0("monosomy ", chromosomes)
        } else if (length(arms) == 1) {
          paste0(arms, " loss")
        } else {
          paste0(paste(arms, collapse = "+"), " codeleted")
        }
        # a positive cutoff is a gain (at or above), a negative one a loss
        cutoff <- as.numeric(pattern)
        if (cutoff > 0) label <- paste0(paste(arms, collapse = "+"), " gain")
        arm_values %>%
          filter(if_all(all_of(arms), ~ !is.na(.x) &
                          (if (cutoff > 0) .x >= cutoff else .x <= cutoff))) %>%
          transmute(id = Kids_First_Biospecimen_ID, evidence = label)
      }
    )
  }

  hits <- pmap_dfr(rules, function(rule_id, OPC_molecular_subtype, lesion_type, gene,
                                   pattern, ...) {
    matches_for_rule(lesion_type, gene, pattern) %>%
      inner_join(assay_ids[[lesion_type]], by = "id",
                 relationship = "many-to-many") %>%
      transmute(match_id, rule_id, OPC_molecular_subtype, evidence)
  }) %>%
    distinct()

  # rules that share an and_group count only together: the tumor has to meet
  # every rule of the group, and the evidence lists all of them
  if ("and_group" %in% names(rules) && any(!is.na(rules$and_group))) {
    groups <- rules %>%
      filter(!is.na(and_group)) %>%
      group_by(and_group) %>%
      mutate(group_size = n()) %>%
      ungroup() %>%
      select(rule_id, and_group, group_size)
    combined <- hits %>%
      inner_join(groups, by = "rule_id") %>%
      group_by(match_id, OPC_molecular_subtype, and_group, group_size) %>%
      filter(n_distinct(rule_id) == group_size) %>%
      summarise(evidence = paste(unique(evidence), collapse = " + "),
                .groups = "drop") %>%
      transmute(match_id, rule_id = and_group, OPC_molecular_subtype, evidence)
    hits <- bind_rows(hits %>% filter(!rule_id %in% groups$rule_id), combined)
  }
  hits
}
