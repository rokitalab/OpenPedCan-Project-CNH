#' Fusions of interest for the LGG / HGG kinase fusions and the Ewing sarcoma
#' hallmark fusions, from the annoFuse putative oncogenic fusions plus the DGD
#' panel fusions (ported from the archived `fusion-summary` module)
#'
#' LGG / HGG: the exact fusions KIAA1549--BRAF, FGFR1--TACC1 and MYB--QKI, and any
#' fusion with BRAF, ALK, ROS1, NTRK1/2/3, PDGFRA, FGFR1/2, MYB, MYBL1, RAF1 or
#' MET. Fusions of a non-kinase gene are kept as they are. A fusion with a 3'
#' kinase is kept when it is in-frame and retains the kinase domain. A fusion with
#' a 5' kinase is kept when it is in-frame and retains the kinase domain, or when
#' it lost the kinase domain but its reciprocal fusion (3' partner first) is
#' in-frame and retains it (a reciprocal fusion).
#'
#' @param root_dir Repository root (data/ holds the fusion files).
#' @return A list with `lgg_hgg` and `ewings`: one row per RNA / fusion panel
#'   biospecimen (Kids_First_Biospecimen_ID) and one 0/1 column per fusion.
fusion_calls_of_interest <- function(root_dir) {
  data_dir <- file.path(root_dir, "data")
  cols <- c("Sample", "FusionName", "Fusion_Type", "Gene1A", "Gene1B", "Gene2A",
            "Gene2B", "Gene1A_anno", "Gene1B_anno", "reciprocal_exists",
            "DomainRetainedGene1A", "DomainRetainedGene1B")

  fusions <- bind_rows(
    read_tsv(file.path(data_dir, "fusion-putative-oncogenic.tsv"),
             show_col_types = FALSE, col_types = cols(.default = "c")),
    read_tsv(file.path(data_dir, "fusion-dgd.tsv.gz"),
             show_col_types = FALSE, col_types = cols(.default = "c"))
  ) %>%
    select(any_of(cols)) %>%
    distinct()

  filter_fusion <- function(fuses = NULL, genes = NULL) {
    keep <- rep(FALSE, nrow(fusions))
    if (!is.null(fuses)) keep <- keep | fusions$FusionName %in% fuses
    if (!is.null(genes)) {
      keep <- keep | fusions$Gene1A %in% genes | fusions$Gene2A %in% genes |
        fusions$Gene1B %in% genes | fusions$Gene2B %in% genes
    }
    fusions[keep, ] %>% distinct()
  }

  # biospecimens with an RNA-Seq or a fusion panel (DGD) assay, not TCGA / GTEx
  hist <- read_tsv(file.path(data_dir, "histologies-base.tsv"),
                   guess_max = 100000, show_col_types = FALSE)
  specimens <- hist %>%
    filter(experimental_strategy == "RNA-Seq" | grepl("DGD_FUSIP", aliquot_id),
           !cohort %in% c("TCGA", "GTEx")) %>%
    pull(Kids_First_Biospecimen_ID)

  to_matrix <- function(df, fuses) {
    out <- df %>%
      distinct(Sample, FusionName) %>%
      mutate(Count = 1) %>%
      pivot_wider(names_from = FusionName, values_from = Count,
                  values_fill = list(Count = 0)) %>%
      right_join(tibble(Sample = specimens), by = "Sample") %>%
      replace(is.na(.), 0) %>%
      rename(Kids_First_Biospecimen_ID = Sample)
    out[, setdiff(fuses, colnames(out)[-1])] <- 0
    out %>% arrange(Kids_First_Biospecimen_ID)
  }

  ## Ewing sarcoma
  ewings_fuses <- c("EWSR1--FLI1", "EWSR1--ERG", "FUS--ERG")
  ewings <- to_matrix(filter_fusion(fuses = ewings_fuses), ewings_fuses)

  ## LGG / HGG
  lgg_fuses <- c("KIAA1549--BRAF", "FGFR1--TACC1", "MYB--QKI")
  lgg_hgg_genes <- c("BRAF", "ALK", "ROS1", "NTRK1", "NTRK2", "NTRK3", "PDGFRA",
                     "FGFR2", "FGFR1", "MYB", "MYBL1", "RAF1", "MET")
  all_fuse <- filter_fusion(fuses = lgg_fuses, genes = lgg_hgg_genes)

  annot <- read.delim(system.file("extdata", "genelistreference.txt",
                                  package = "annoFuseData"))
  kinases <- annot %>% filter(type == "Kinase") %>% pull(Gene_Symbol)

  lgg_fuses_df <- tibble(fusion = lgg_fuses) %>%
    separate(fusion, into = c("Gene1A", "Gene1B"), remove = FALSE) %>%
    mutate(kinase_fusion = Gene1A %in% kinases | Gene1B %in% kinases)
  nonkinase_fuses <- lgg_fuses_df %>% filter(!kinase_fusion) %>% pull(fusion)
  nonkinase_genes <- setdiff(lgg_hgg_genes, kinases)
  nonkinase <- filter_fusion(fuses = nonkinase_fuses, genes = nonkinase_genes)

  # a 3' kinase, in-frame, kinase domain retained
  three_prime <- all_fuse %>%
    filter(grepl("Kinase", Gene1B_anno), Fusion_Type == "in-frame",
           DomainRetainedGene1B == "Yes") %>%
    select(Sample, FusionName) %>% distinct()
  # a 5' kinase, in-frame, kinase domain retained
  five_prime_intact <- all_fuse %>%
    filter(grepl("Kinase", Gene1A_anno), DomainRetainedGene1A == "Yes",
           Fusion_Type == "in-frame") %>%
    select(Sample, FusionName) %>% distinct()
  # a 5' kinase that lost the domain: kept when the reciprocal fusion (3' partner
  # first) is in-frame and retains the kinase domain
  five_prime_lost <- all_fuse %>%
    filter(grepl("Kinase", Gene1A_anno), reciprocal_exists == "TRUE",
           DomainRetainedGene1A == "No") %>%
    distinct(Sample, FusionName, Gene1A, Gene1B) %>%
    mutate(reciprocal = paste(Gene1B, Gene1A, sep = "--"))
  five_prime_reciprocal <- five_prime_lost %>%
    select(Sample, five_prime = FusionName, FusionName = reciprocal) %>%
    left_join(all_fuse, by = c("Sample", "FusionName")) %>%
    filter(Fusion_Type == "in-frame", DomainRetainedGene1B == "Yes") %>%
    select(Sample, FusionName = five_prime) %>% distinct()

  lgg_hgg <- bind_rows(select(nonkinase, Sample, FusionName), three_prime,
                       five_prime_intact, five_prime_reciprocal) %>%
    distinct() %>%
    to_matrix(lgg_fuses)

  list(lgg_hgg = lgg_hgg, ewings = ewings)
}
