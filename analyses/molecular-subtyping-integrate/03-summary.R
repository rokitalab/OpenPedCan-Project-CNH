################################################################################
# 03-summary.R
#
# One page with the totals of the subtyping: the change from the previous
# release (02-compare-to-release.R), the methylation mismatches and the lists
# still to review. Run after run-molecular-subtyping.sh, run-subtyping-integrate.sh
# and 02-compare-to-release.R.
#
# Usage: Rscript --vanilla 03-summary.R
# Output: results/summary.md
################################################################################

suppressMessages(library(tidyverse))

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
sub_dir <- file.path(root_dir, "analyses", "molecular-subtyping", "results")
int_dir <- file.path(root_dir, "analyses", "molecular-subtyping-integrate", "results")

read <- function(...) read_tsv(file.path(...), show_col_types = FALSE,
                               col_types = cols(.default = "c"))
md_table <- function(df) {
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

changes <- read(int_dir, "subtype_changes_from_release.tsv")
subtype <- read(sub_dir, "molecular_subtype.tsv")
mixups <- read_tsv(file.path(root_dir, "analyses", "molecular-subtyping", "input",
                             "possible_methylation_mixups.tsv"),
                   show_col_types = FALSE, col_types = cols(.default = "c"))
likely <- read(sub_dir, "methylation_likely_mismatch.tsv")
checks <- read(sub_dir, "subtyping_checks.tsv")

out <- c(
  "# Subtyping summary", "",
  paste0("Tumors in the subtype table: ", nrow(subtype),
         "; with a subtype: ", sum(!is.na(subtype$molecular_subtype)), "."), "",
  "## Change from the previous release", "",
  paste0(nrow(changes), " tumors changed their subtype. Change type by verdict ",
         "(`subtype_changes_from_release.tsv`, transitions in ",
         "`subtype_change_transitions.tsv`, the tumors to review in ",
         "`subtype_changes_needing_review.tsv`):"), "",
  md_table(changes %>% count(change_type, verdict) %>%
             pivot_wider(names_from = verdict, values_from = n, values_fill = 0) %>%
             mutate(across(everything(), as.character))), "",
  "## Methylation mismatches", "",
  paste0("`input/possible_methylation_mixups.tsv`: likely = genomics supports the ",
         "diagnosis (the arrays are collected for removal, ",
         "`methylation_likely_mismatch.tsv`); possible = unresolved; kept = reviewed ",
         "and the methylation call stands."), "",
  md_table(mixups %>% count(mismatch_call) %>% mutate(n = as.character(n))), "",
  "Likely mismatches:", "",
  md_table(likely %>% select(sample_id, pathology_diagnosis, molecular_subtype,
                             methyl_class, methyl_score)), "",
  "## Flags of the checks (`subtyping_checks.tsv`)", "",
  md_table(checks %>% filter(status == "FLAG") %>% select(area, check, n)), ""
)
writeLines(out, file.path(int_dir, "summary.md"))
message("Wrote ", file.path(int_dir, "summary.md"))
