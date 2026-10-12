#' The subtype used to decide which tumors an annotation applies to
#'
#' A manually reviewed subtype (`annotation_subtype`, from
#' input/manual_subtype_review.tsv) wins over the automated one, so reviewed
#' tumors (e.g. an EGFR-altered DMG the automated call labelled HGG) get the DMG
#' annotations.
annotation_subtype <- function(df) {
  if ("annotation_subtype" %in% names(df)) {
    coalesce(df$annotation_subtype, df$OPC_molecular_subtype)
  } else {
    df$OPC_molecular_subtype
  }
}

#' Gene mutation status as an annotation (not part of the subtype)
#'
#' For each gene in `annotation_genes`, adds `<gene>_mutation_status` (gene in
#' lowercase, e.g. `egfr_mutation_status`)
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
    applies <- str_detect(annotation_subtype(df),
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
    activating_calls <- snv %>%
      filter(Hugo_Symbol == gene,
             str_detect(HGVSp_Short, annotation_genes$alteration_pattern[i])) %>%
      inner_join(dna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
                 relationship = "many-to-many") %>%
      group_by(match_id) %>%
      summarise(activating_variants = paste(sort(unique(HGVSp_Short)),
                                            collapse = "; "),
                .groups = "drop")
    activating <- activating_calls$match_id

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

    df[[paste0(tolower(gene), "_amplification")]] <- case_when(
      !applies ~ NA_character_,
      df$match_id %in% amplified ~ "amplified",
      df$match_id %in% has_cnv ~ "not amplified",
      TRUE ~ "no CNV data"
    )
    df[[paste0(tolower(gene), "_activating_alteration")]] <- case_when(
      !applies ~ NA,
      df$match_id %in% c(activating, amplified) ~ TRUE,
      TRUE ~ FALSE
    )
    df[[paste0(tolower(gene), "_activating_variants")]] <- if_else(
      applies,
      activating_calls$activating_variants[match(df$match_id, activating_calls$match_id)],
      NA_character_
    )
    df[[paste0(tolower(gene), "_mutation_status")]] <- status
    df[[paste0(tolower(gene), "_variants")]] <- if_else(
      applies, variants$variants[match(df$match_id, variants$match_id)],
      NA_character_
    )
  }
  df
}

#' Gene overexpression as an annotation (not part of the subtype)
#'
#' For each gene in `expression_annotations`, adds `<gene>_tpm` and `<gene>_z`
#' (highest across the tumor's RNA biospecimens; z is the z-score of
#' log2(TPM + 1) within RNA library type, from 02-prepare-lesion-tables.R) and
#' `<gene>_expression_status` (overexpressed at or above `z_threshold`, not
#' overexpressed, or no RNA assay) for the tumors whose subtype matches
#' `applies_to_subtype_regex`.
#'
#' @param df Table with match_id and OPC_molecular_subtype.
#' @param expression_annotations gene, applies_to_subtype_regex, z_threshold and
#'   optionally exclude_library_groups (";" separated library groups the gene
#'   cannot be called in).
#' @param expression_long gene, Kids_First_Biospecimen_ID, tpm and z.
#' @param rna_ids Table with match_id and the RNA biospecimen `id`.
annotate_expression <- function(df, expression_annotations, expression_long,
                                rna_ids) {
  for (i in seq_len(nrow(expression_annotations))) {
    gene <- expression_annotations$gene[i]
    applies <- str_detect(annotation_subtype(df),
                          expression_annotations$applies_to_subtype_regex[i])

    # RNA library groups the gene cannot be called in (exclude_library_groups,
    # ";" separated): a tumor with only those libraries is not assessable
    excluded <- if ("exclude_library_groups" %in% names(expression_annotations)) {
      str_split(replace_na(expression_annotations$exclude_library_groups[i], ""),
                ";")[[1]] %>% str_trim() %>% setdiff("")
    } else {
      character()
    }

    rows <- expression_long %>%
      filter(gene == !!gene) %>%
      inner_join(rna_ids, by = c("Kids_First_Biospecimen_ID" = "id"),
                 relationship = "many-to-many")
    per_tumor <- rows %>%
      group_by(match_id) %>%
      summarise(tpm = max(tpm, na.rm = TRUE),
                z = suppressWarnings(max(z[!library_group %in% excluded],
                                         na.rm = TRUE)),
                .groups = "drop") %>%
      mutate(z = if_else(is.finite(z), z, NA_real_))
    not_assessable <- per_tumor$match_id[is.na(per_tumor$z) & length(excluded) > 0 &
                                           per_tumor$match_id %in%
                                           rows$match_id[rows$library_group %in% excluded]]

    idx <- match(df$match_id, per_tumor$match_id)
    value_z <- per_tumor$z[idx]
    df[[paste0(tolower(gene), "_tpm")]] <- if_else(applies, per_tumor$tpm[idx], NA_real_)
    df[[paste0(tolower(gene), "_z")]] <- if_else(applies, round(value_z, 2), NA_real_)
    df[[paste0(tolower(gene), "_expression_status")]] <- case_when(
      !applies ~ NA_character_,
      is.na(value_z) & df$match_id %in% not_assessable ~
        paste0("not assessable (", paste(excluded, collapse = ", "), " library)"),
      is.na(value_z) ~ "no RNA assay",
      value_z >= expression_annotations$z_threshold[i] ~ "overexpressed",
      TRUE ~ "not overexpressed"
    )
  }
  df
}

#' H3 K27-altered mechanism for DMG
#'
#' H3 K27-altered (WHO CNS5) loss of H3K27me3 can come from an H3 K27M/I
#' mutation, EZHIP overexpression, or an EGFR alteration. This lists the
#' mechanisms detected, from DNA and RNA.
#'
#' @param df Table with OPC_molecular_subtype, egfr_activating_alteration and
#'   ezhip_expression_status.
#' @param h3_mutated Logical vector, TRUE when an H3 K27M/I variant was found.
#' @param has_dna,has_rna Logical vectors, TRUE when the tumor has the assay.
add_h3k27_mechanism <- function(df, h3_mutated, has_dna, has_rna) {
  is_dmg <- str_detect(annotation_subtype(df), "^DMG")
  found <- cbind(
    "H3 K27M/I" = h3_mutated %in% TRUE,
    "EGFR alteration" = df$egfr_activating_alteration %in% TRUE,
    "EZHIP overexpression" = df$ezhip_expression_status %in% "overexpressed"
  )
  mech <- apply(found, 1, function(r) {
    if (any(r)) paste(colnames(found)[r], collapse = " + ") else NA_character_
  })
  df$h3k27_altered_mechanism <- case_when(
    !is_dmg ~ NA_character_,
    !is.na(mech) ~ mech,
    has_dna & has_rna &
      str_detect(replace_na(df$ezhip_expression_status, ""), "^not assessable") ~
      "none detected (EZHIP not assessable in an exome capture library)",
    has_dna & has_rna ~ "none detected",
    has_dna ~ "none detected (no RNA assay for EZHIP)",
    has_rna ~ "none detected (no DNA assay)",
    TRUE ~ "no DNA or RNA assay"
  )
  df
}

#' H3 mutation status for HGG, DMG and DHG tumors (annotation)
#'
#' @param df Table with match_id and OPC_molecular_subtype.
#' @param h3_hits Lesion hits (match_id, evidence) for the H3 K27M/I and
#'   G34/G35 R/V rules.
#' @param dna_ids Table with match_id of tumors that have a DNA assay.
#' @return `h3_status` (mutant, wildtype, or not assessed) and `h3_variant`.
add_h3_status <- function(df, h3_hits, dna_ids,
                          applies_regex = "^(HGG|DMG|DHG)") {
  applies <- str_detect(annotation_subtype(df), applies_regex)
  variants <- h3_hits %>%
    group_by(match_id) %>%
    summarise(h3_variant = paste(sort(unique(evidence)), collapse = "; "),
              .groups = "drop")
  df$h3_status <- case_when(
    !applies ~ NA_character_,
    df$match_id %in% variants$match_id ~ "mutant",
    df$match_id %in% dna_ids$match_id ~ "wildtype",
    TRUE ~ "not assessed"
  )
  df$h3_variant <- if_else(
    applies, variants$h3_variant[match(df$match_id, variants$match_id)],
    NA_character_
  )
  df
}

#' EGFR in DMG labels
#'
#' - A DMG supported only by methylation or the diagnosis ("DMG, H3 K28"), with
#'   a DNA assay, an EGFR activating alteration and no H3 variant, is H3
#'   K28-altered through EGFR. The free-text diagnosis can say "K27M mutant"
#'   without a variant being found, so the EGFR lesion is used instead.
#' - A DMG with a confirmed H3 variant ("DMG, H3 K28-mutant") that also has an
#'   EGFR activating alteration carries both labels.
relabel_egfr_only <- function(df) {
  df %>%
    mutate(OPC_molecular_subtype = case_when(
      OPC_molecular_subtype == "DMG, H3 K28" &
        egfr_activating_alteration %in% TRUE & h3_status %in% "wildtype" ~
        "DMG, H3 K28-altered, EGFR",
      OPC_molecular_subtype == "DMG, H3 K28-mutant" &
        egfr_activating_alteration %in% TRUE ~ "DMG, H3 K28-mutant, EGFR",
      TRUE ~ OPC_molecular_subtype
    )) %>%
    add_ezhip_label()
}

#' EZHIP overexpression in DMG labels
#'
#' EZHIP overexpression is added to the label of a DMG, after any EGFR term
#' (`DMG, H3 K28-mutant, EZHIP`, `DMG, H3 K28-altered, EGFR, EZHIP`). A DMG with
#' no H3 variant whose only support is methylation or the diagnosis
#' ("DMG, H3 K28") becomes `DMG, H3 K28-altered, EZHIP`, since EZHIP
#' overexpression is the route to H3 K27 alteration.
add_ezhip_label <- function(df) {
  df %>%
    mutate(OPC_molecular_subtype = case_when(
      !str_detect(OPC_molecular_subtype, "^DMG") ~ OPC_molecular_subtype,
      !(ezhip_expression_status %in% "overexpressed") ~ OPC_molecular_subtype,
      str_detect(OPC_molecular_subtype, "EZHIP") ~ OPC_molecular_subtype,
      OPC_molecular_subtype == "DMG, H3 K28" ~ "DMG, H3 K28-altered, EZHIP",
      TRUE ~ paste0(OPC_molecular_subtype, ", EZHIP")
    ))
}

#' RTK fusion term for infant-type hemispheric glioma
#'
#' `HGG, IHG` (methylation) becomes `HGG, IHG, <NTRK|ALK|ROS1|MET>-altered`
#' when the tumor has that fusion in the in-frame, kinase-retained fusion
#' summary (priority NTRK, ALK, ROS1, MET, as in the HGG module). The fusion is
#' then the confirmed alteration.
#'
#' @param ihg Table with match_id, ihg_term and ihg_evidence.
relabel_ihg <- function(df, ihg) {
  df %>%
    left_join(ihg, by = "match_id") %>%
    mutate(
      is_ihg = OPC_molecular_subtype == "HGG, IHG" & !is.na(ihg_term),
      OPC_molecular_subtype = if_else(
        is_ihg, paste0("HGG, IHG, ", ihg_term, "-altered"), OPC_molecular_subtype
      ),
      lesion_status = if_else(is_ihg, "confirmed", lesion_status),
      lesion_evidence = if_else(is_ihg, ihg_evidence, lesion_evidence),
      evidence_level = if_else(is_ihg, "methylation + lesion", evidence_level)
    ) %>%
    select(-is_ihg, -ihg_term, -ihg_evidence)
}
