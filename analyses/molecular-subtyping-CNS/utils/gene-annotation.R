#' Gene mutation status as an annotation (not part of the subtype)
#'
#' For each gene in `annotation_genes`, adds `<gene>_mutation_status`
#' (mutated, not detected, or no DNA assay) and `<gene>_variants` to the
#' tumors whose subtype matches `applies_to_subtype_regex`; other tumors get NA.
#'
#' @param df Table with match_id and OPC_molecular_subtype.
#' @param annotation_genes gene and applies_to_subtype_regex columns.
#' @param snv SNV table with Hugo_Symbol, HGVSp_Short, Variant_Classification
#'   and Kids_First_Biospecimen_ID.
#' @param dna_ids Table with match_id and the DNA biospecimen `id`.
annotate_genes <- function(df, annotation_genes, snv, dna_ids) {
  for (i in seq_len(nrow(annotation_genes))) {
    gene <- annotation_genes$gene[i]
    applies <- str_detect(df$OPC_molecular_subtype,
                          annotation_genes$applies_to_subtype_regex[i])

    variants <- snv %>%
      # coding, protein-altering variants only (no silent or non-coding)
      filter(Hugo_Symbol == gene,
             !Variant_Classification %in% c("Silent", "Intron", "3'UTR",
                                            "5'UTR", "3'Flank", "5'Flank",
                                            "RNA"),
             !str_detect(HGVSp_Short, "=$")) %>%
      inner_join(dna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
                 relationship = "many-to-many") %>%
      group_by(match_id) %>%
      summarise(variants = paste(sort(unique(HGVSp_Short)), collapse = "; "),
                .groups = "drop")

    status <- case_when(
      !applies ~ NA_character_,
      df$match_id %in% variants$match_id ~ "mutated",
      df$match_id %in% dna_ids$match_id ~ "not detected",
      TRUE ~ "no DNA assay"
    )
    df[[paste0(gene, "_mutation_status")]] <- status
    df[[paste0(gene, "_variants")]] <- if_else(
      applies, variants$variants[match(df$match_id, variants$match_id)],
      NA_character_
    )
  }
  df
}
