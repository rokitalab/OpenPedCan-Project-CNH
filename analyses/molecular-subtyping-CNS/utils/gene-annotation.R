#' Gene mutation status as an annotation (not part of the subtype)
#'
#' For each gene in `annotation_genes`, adds `<gene>_mutation_status`
#' (mutated, not detected, or no DNA assay) and `<gene>_variants` to the
#' tumors whose subtype matches `applies_to_subtype_regex`; other tumors get NA.
#'
#' `<gene>_activating_alteration` is TRUE when a variant matches
#' `alteration_pattern`, or the gene is amplified, and
#' `<gene>_amplification` reports copy-number amplification (a segment
#' overlapping chrom:start-end with copy number above twice the ploidy, the
#' convention in focal-cn-file-preparation).
#'
#' @param df Table with match_id and OPC_molecular_subtype.
#' @param annotation_genes gene, applies_to_subtype_regex, alteration_pattern,
#'   chrom, start and end columns.
#' @param snv SNV table with Hugo_Symbol, HGVSp_Short, Variant_Classification
#'   and Kids_First_Biospecimen_ID.
#' @param dna_ids Table with match_id and the DNA biospecimen `id`.
#' @param cnv_segments Segments with Kids_First_Biospecimen_ID, chrom,
#'   loc.start, loc.end, copy.num and tumor_ploidy.
annotate_genes <- function(df, annotation_genes, snv, dna_ids, cnv_segments) {
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
    # variants matching the gene's activating pattern
    activating <- snv %>%
      filter(Hugo_Symbol == gene,
             str_detect(HGVSp_Short, annotation_genes$alteration_pattern[i])) %>%
      inner_join(dna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
                 relationship = "many-to-many") %>%
      pull(match_id) %>%
      unique()

    # amplification from copy-number segments overlapping the gene
    amplified <- cnv_segments %>%
      filter(chrom == annotation_genes$chrom[i],
             loc.start <= annotation_genes$end[i],
             loc.end >= annotation_genes$start[i],
             !is.na(copy.num), !is.na(tumor_ploidy),
             copy.num > 2 * tumor_ploidy) %>%
      inner_join(dna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
                 relationship = "many-to-many") %>%
      pull(match_id) %>%
      unique()
    has_cnv <- cnv_segments %>%
      distinct(Kids_First_Biospecimen_ID) %>%
      inner_join(dna_ids, by = c("Kids_First_Biospecimen_ID" = "id")) %>%
      pull(match_id)

    df[[paste0(gene, "_amplification")]] <- case_when(
      !applies ~ NA_character_,
      df$match_id %in% amplified ~ "amplified",
      df$match_id %in% has_cnv ~ "not amplified",
      TRUE ~ "no CNV data"
    )
    df[[paste0(gene, "_activating_alteration")]] <- case_when(
      !applies ~ NA,
      df$match_id %in% c(activating, amplified) ~ TRUE,
      TRUE ~ FALSE
    )
    df[[paste0(gene, "_mutation_status")]] <- status
    df[[paste0(gene, "_variants")]] <- if_else(
      applies, variants$variants[match(df$match_id, variants$match_id)],
      NA_character_
    )
  }
  df
}

#' H3 K27-altered mechanism for DMG
#'
#' H3 K27-altered (WHO CNS5) loss of H3K27me3 can come from an H3 K27M/I
#' mutation, EZHIP overexpression, or an EGFR alteration. This lists the
#' mechanisms detected from DNA. EZHIP overexpression needs expression data and
#' is not assessed here.
#'
#' @param df Table with OPC_molecular_subtype and EGFR_activating_alteration.
#' @param h3_mutated Logical vector, TRUE when an H3 K27M/I variant was found.
#' @param has_dna Logical vector, TRUE when the tumor has a DNA assay.
add_h3k27_mechanism <- function(df, h3_mutated, has_dna) {
  is_dmg <- str_detect(df$OPC_molecular_subtype, "^DMG")
  mech <- map2_chr(h3_mutated %in% TRUE, df$EGFR_activating_alteration %in% TRUE,
                   function(h3, egfr) {
    found <- c(if (h3) "H3 K27M/I", if (egfr) "EGFR alteration")
    if (length(found) == 0) NA_character_ else paste(found, collapse = " + ")
  })
  df$h3k27_altered_mechanism <- case_when(
    !is_dmg ~ NA_character_,
    !is.na(mech) ~ mech,
    has_dna ~ "none detected (EZHIP overexpression not assessed)",
    TRUE ~ "no DNA assay"
  )
  df
}
