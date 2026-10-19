################################################################################
# extend-plot-mapping.R
#
# Add the rows that figures/palettes/plot-mapping.tsv is missing for the tumors of
# the integrated histologies table. The table maps a cohort, a broad histology
# and a cancer group to a plot group (with the display name, abbreviation and
# colors). Rows that exist are never changed; each added row has a `status`
# saying how its plot group was chosen, to be reviewed. The first version of the
# table (the original rows) is for the PBTA cohort.
#
# A new (cohort, broad histology, cancer group) gets the plot group of, in order:
#   0. a plot group chosen for a cancer group (`chosen_plot_group` below)
#   1. neuroblastic tumors (neuroblastoma, ganglioneuroblastoma, ganglioneuroma) are
#      one plot group, Neuroblastoma: outside the CNS cohort PBTA, and in PBTA
#      the broad histology Peripheral neuroblastic tumor
#   2. PBTA only: "Other tumor" when the cancer group has fewer than `min_tumors`
#      tumors in the cohort (small groups are kept together), unless the cohort
#      already has a row of that cancer group (then that row decides)
#   2b. TARGET: the cancer groups of `target_plot_group_hex` (leukemias, kidney
#      tumors, osteosarcoma) are their own plot groups
#   3. the same cancer group elsewhere in the table (the cohort's own rows first,
#      then the original rows)
#   4. a cancer group assignment below (broad histologies that map to several plot groups)
#   5. non-CNS sarcomas: Mesenchymal tumor
#   6. the broad histology, when it maps to one plot group in the table
#   7. Other tumor
# Added rows (with a `status`) of rule 2 are updated when the rule changes; the
# original rows are never changed.
#
# The DGD, TCGA, PPTC, Maris and CPTAC cohorts are not in the table (`excluded_cohorts`): their rows
# are removed and not added again. A row whose cohort, broad histology and cancer
# group no longer have a tumor in the table (a cancer group that was dropped or
# renamed, or a broad histology that was renamed) is removed, original rows
# included; a combination that comes back is added again with a `status`.
#
# Usage: Rscript --vanilla figures/scripts/extend-plot-mapping.R [min_tumors]
#   (default 5; reads analyses/molecular-subtyping-integrate/results/histologies.tsv)
################################################################################

suppressMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
mapping_file <- file.path(root_dir, "figures", "palettes", "plot-mapping.tsv")
args <- commandArgs(trailingOnly = TRUE)
min_tumors <- if (length(args) >= 1) as.integer(args[1]) else 5L
small_n_cohorts <- "PBTA"
excluded_cohorts <- c("DGD", "TCGA", "PPTC", "Maris", "CPTAC")

by_cancer_group <- c(
  "Glioblastoma" = "Other high-grade glioma",
  "Glioblastoma Multiforme" = "Other high-grade glioma",
  "High-grade astrocytoma with piloid features" = "Other high-grade glioma",
  "Brain Lower Grade Glioma" = "Low-grade glioma",
  "Medulloblastoma with myogenic differentiation" = "Medulloblastoma",
  "Intraocular medulloepithelioma" = "Other CNS embryonal tumor",
  "Testicular Germ Cell Tumor" = "Germ cell tumor",
  "Pituicytoma" = "Other tumor")
mesenchymal_broad <- c("Sarcoma", "Spindle cell tumor")
# plot groups chosen for a cancer group, whatever its size
chosen_plot_group <- c("Hybrid nerve sheath tumor" = "Schwannoma/MPNST",
                       "Malignant peripheral nerve sheath tumor" = "Schwannoma/MPNST")
# the neuroblastic tumors are one plot group, Neuroblastoma (the TAPESTRY color of
# OpenPedCanExpress): all of them outside the CNS cohort PBTA, and those of the
# broad histology Peripheral neuroblastic tumor in PBTA
peripheral_nbl <- c("Neuroblastoma" = "Neuroblastoma",
                    "Ganglioneuroblastoma" = "Neuroblastoma",
                    "Ganglioneuroma" = "Neuroblastoma")
new_plot_group_hex <- c("Neuroblastoma" = "#FF7F0E")
nbl_abbreviation <- c("Neuroblastoma" = "NBL", "Ganglioneuroblastoma" = "NBL",
                      "Ganglioneuroma" = "NBL")
is_peripheral_nbl <- function(cohort, broad_histology, cancer_group) {
  cancer_group %in% names(peripheral_nbl) &
    (cohort != small_n_cohorts | broad_histology == "Peripheral neuroblastic tumor")
}

# TARGET cancer groups that are their own plot group (not Other tumor), and their
# colors (the TAPESTRY scheme of OpenPedCanExpress, public/tapestry-colors.js)
target_plot_group_hex <- c("Acute Lymphoblastic Leukemia" = "#D62728",
                           "Acute Myeloid Leukemia" = "#8C1C13",
                           "Clear cell sarcoma of the kidney" = "#7FCDD8",
                           "Osteosarcoma" = "#556B2F",
                           "Rhabdoid tumor of the kidney" = "#1F4E79",
                           "Wilms tumor" = "#17A2B8")

mode_of <- function(x) names(sort(table(x), decreasing = TRUE))[1]
abbreviate_group <- function(x) {
  if (nchar(x) <= 20) return(x)
  w <- str_split(x, "[ /-]+")[[1]]
  w <- w[nchar(w) > 3 | str_detect(w, "^[A-Z]")]
  toupper(paste(substr(w, 1, 1), collapse = ""))
}

mapping <- read_tsv(mapping_file, show_col_types = FALSE, na = character(),
                    col_types = cols(.default = col_character()))
if (!"cohort" %in% names(mapping)) mapping <- mutate(mapping, cohort = "PBTA")
if (!"status" %in% names(mapping)) mapping <- mutate(mapping, status = "")
mapping <- select(mapping, cohort, broad_histology, cancer_group,
                  broad_histology_display, plot_group, broad_histology_hex,
                  cancer_group_abbreviation, plot_group_hex, status)
mapping <- filter(mapping, !cohort %in% excluded_cohorts)
stopifnot(!anyDuplicated(mapping[c("cohort", "broad_histology", "cancer_group")]))

hist <- read_tsv(file.path(root_dir, "analyses", "molecular-subtyping-integrate",
                           "results", "histologies.tsv"),
                 guess_max = 100000, show_col_types = FALSE) %>%
  filter(sample_type == "Tumor", !cohort %in% excluded_cohorts) %>%
  mutate(tumor_key = coalesce(match_id, Kids_First_Biospecimen_ID),
         cancer_group = replace_na(cancer_group, "NA")) %>%
  # a tumor counts once in each cohort, broad histology and cancer group it has
  distinct(tumor_key, cohort, broad_histology, cancer_group)

# rows without a tumor are removed
present <- hist %>% distinct(cohort, broad_histology, cancer_group)
stale <- mapping %>%
  anti_join(present, by = c("cohort", "broad_histology", "cancer_group"))
if (nrow(stale) > 0) {
  message("Removed ", nrow(stale), " rows without a tumor:")
  print(count(stale, cohort, broad_histology, cancer_group), n = 50)
  mapping <- anti_join(mapping, stale, by = c("cohort", "broad_histology", "cancer_group"))
}

group_size <- hist %>% count(cohort, cancer_group, name = "group_tumors")
needed <- hist %>% count(cohort, broad_histology, cancer_group, name = "tumors") %>%
  anti_join(mapping, by = c("cohort", "broad_histology", "cancer_group")) %>%
  left_join(group_size, by = c("cohort", "cancer_group"))

plot_group_hex <- mapping %>% group_by(plot_group) %>%
  summarise(plot_group_hex = mode_of(plot_group_hex), .groups = "drop") %>%
  bind_rows(tibble(plot_group = c(names(new_plot_group_hex), names(target_plot_group_hex)),
                   plot_group_hex = c(unname(new_plot_group_hex),
                                      unname(target_plot_group_hex)))) %>%
  distinct(plot_group, .keep_all = TRUE)

# the added rows of the neuroblastic tumors outside PBTA follow rule 2
update <- mapping$status != "" &
  is_peripheral_nbl(mapping$cohort, mapping$broad_histology, mapping$cancer_group)
mapping$plot_group[update] <- unname(peripheral_nbl[mapping$cancer_group[update]])
mapping$cancer_group_abbreviation[update] <- unname(nbl_abbreviation[mapping$cancer_group[update]])
mapping$broad_histology_display[update] <- mapping$plot_group[update]
mapping$plot_group_hex[update] <- unname(new_plot_group_hex[mapping$plot_group[update]])
mapping$broad_histology_hex[update] <- mapping$plot_group_hex[update]
mapping$status[update] <- "added, review: neuroblastic tumor as one plot group (Neuroblastoma)"

# the added TARGET rows of rule 2b
update_target <- mapping$status != "" & mapping$cohort == "TARGET" &
  mapping$cancer_group %in% names(target_plot_group_hex)
mapping$plot_group[update_target] <- mapping$cancer_group[update_target]
mapping$broad_histology_display[update_target] <- mapping$cancer_group[update_target]
mapping$plot_group_hex[update_target] <- unname(target_plot_group_hex[mapping$cancer_group[update_target]])
mapping$broad_histology_hex[update_target] <- mapping$plot_group_hex[update_target]
mapping$status[update_target] <- "added, review: TARGET cancer group as its own plot group"

added <- map_dfr(seq_len(nrow(needed)), function(i) {
  coh <- needed$cohort[i]; b <- needed$broad_histology[i]; cg <- needed$cancer_group[i]
  # the same cancer group: the cohort's own rows first, then the original rows,
  # then any
  same_cg <- filter(mapping, cancer_group == cg, cg != "NA")
  if (any(same_cg$cohort == coh)) same_cg <- filter(same_cg, cohort == coh)
  else if (any(same_cg$status == "")) same_cg <- filter(same_cg, status == "")
  same_b <- filter(mapping, broad_histology == b)
  in_cohort <- mapping$cancer_group[mapping$cohort == coh & mapping$status == ""]
  if (cg %in% names(chosen_plot_group)) {
    pg <- unname(chosen_plot_group[cg]); how <- "plot group chosen for the cancer group"
  } else if (is_peripheral_nbl(coh, b, cg)) {
    pg <- unname(peripheral_nbl[cg]); how <- "neuroblastic tumor as one plot group (Neuroblastoma)"
  } else if (coh %in% small_n_cohorts && needed$group_tumors[i] < min_tumors &&
             !cg %in% in_cohort) {
    pg <- "Other tumor"; how <- paste0("small group (< ", min_tumors, " tumors in ", coh, "): Other tumor")
  } else if (coh == "TARGET" && cg %in% names(target_plot_group_hex)) {
    pg <- cg; how <- "TARGET cancer group as its own plot group"
  } else if (nrow(same_cg) > 0) {
    pg <- mode_of(same_cg$plot_group); how <- "cancer group in the table"
  } else if (cg %in% names(by_cancer_group)) {
    pg <- by_cancer_group[[cg]]; how <- "cancer group assigned"
  } else if (b %in% mesenchymal_broad) {
    pg <- "Mesenchymal tumor"; how <- "sarcoma: Mesenchymal tumor"
  } else if (n_distinct(same_b$plot_group) == 1) {
    pg <- same_b$plot_group[1]; how <- "broad histology in the table"
  } else {
    pg <- "Other tumor"; how <- "default: Other tumor"
  }
  single <- nrow(same_b) > 0 && n_distinct(same_b$plot_group) == 1
  tibble(
    cohort = coh, broad_histology = b, cancer_group = cg,
    broad_histology_display = if (single) mode_of(same_b$broad_histology_display) else pg,
    plot_group = pg,
    broad_histology_hex = if (single) mode_of(same_b$broad_histology_hex)
                          else plot_group_hex$plot_group_hex[plot_group_hex$plot_group == pg],
    cancer_group_abbreviation = if (is_peripheral_nbl(coh, b, cg)) unname(nbl_abbreviation[cg])
                                else if (nrow(same_cg) > 0) mode_of(same_cg$cancer_group_abbreviation)
                                else if (cg == "NA") "NA" else abbreviate_group(cg),
    plot_group_hex = plot_group_hex$plot_group_hex[plot_group_hex$plot_group == pg],
    status = paste0("added, review: ", how,
                    if_else(nrow(same_cg) == 0 & cg != "NA", "; abbreviation generated", "")))
})

out <- bind_rows(mapping, added) %>%
  mutate(across(c(broad_histology_hex, plot_group_hex), str_to_upper)) %>%
  arrange(cohort, broad_histology, cancer_group)
write_tsv(out, mapping_file, na = "")
message("Rows: ", nrow(mapping), " -> ", nrow(out), " (added ", nrow(added), ")")
if (nrow(added) > 0) print(count(added, cohort, plot_group), n = 60)
