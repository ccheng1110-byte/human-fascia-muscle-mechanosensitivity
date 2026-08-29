options(stringsAsFactors = FALSE)

# Step 08A: low-volume external atlas screen.
# This script uses only author-provided files already saved locally.
# It does not download any file and does not treat cells as biological replicates.

project_dir <- "."

step07_path <- file.path(
  project_dir, "results", "05_specificity_robustness", "GSE173252",
  "GSE173252_DD_candidate_specificity_summary_v1.csv"
)
atlas_dir <- file.path(
  project_dir, "data", "external_validation", "skin_fibroblast_atlas_2025"
)
marker_rank_path <- file.path(atlas_dir, "author_outputs", "degs_fbs.csv")
dataset_inventory_path <- file.path(atlas_dir, "metadata", "scrna_seq_data.csv")
result_dir <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(step07_path, marker_rank_path, dataset_inventory_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

safe_value <- function(x, default = NA_real_) {
  if (length(x) == 0L || all(is.na(x))) return(default)
  x[[1L]]
}

# -----------------------------------------------------------------------------
# 1. Freeze the Step-07 shortlist before examining the external atlas rankings.
# -----------------------------------------------------------------------------
step07 <- read.csv(step07_path, check.names = FALSE)
required_step07_columns <- c(
  "module", "gene", "step06_candidate_tier", "specificity_class",
  "myofib_percent_cells_expressing",
  "minimum_difference_over_all_donors_and_identities"
)
missing_step07_columns <- setdiff(required_step07_columns, colnames(step07))
if (length(missing_step07_columns) > 0L) {
  stop(
    "Step-07 candidate table is missing columns: ",
    paste(missing_step07_columns, collapse = ", ")
  )
}

shortlist <- step07[
  grepl("^[AB]_", step07$step06_candidate_tier) &
    step07$specificity_class == "robust_Myofib_over_all_three_identities",
  , drop = FALSE
]
shortlist <- shortlist[!duplicated(shortlist$gene), , drop = FALSE]
if (nrow(shortlist) == 0L) {
  stop("No frozen A/B robust Step-07 candidates were found.")
}

# -----------------------------------------------------------------------------
# 2. Convert the authors' ranked marker matrix into a long table.
# -----------------------------------------------------------------------------
atlas_markers <- read.csv(
  marker_rank_path,
  check.names = FALSE,
  na.strings = c("", "NA")
)
if (ncol(atlas_markers) < 3L) {
  stop("The atlas marker-rank file has an unexpected structure.")
}

# The first column is an exported row index, not a fibroblast state.
marker_columns <- colnames(atlas_markers)[-1L]
target_group <- "F7: Fascia-like myofibroblast"
if (!target_group %in% marker_columns) {
  stop("Target atlas state was not found: ", target_group)
}

rank_long <- do.call(rbind, lapply(marker_columns, function(group_name) {
  genes <- trimws(as.character(atlas_markers[[group_name]]))
  data.frame(
    group = group_name,
    rank = seq_along(genes),
    gene = genes,
    stringsAsFactors = FALSE
  )
}))
rank_long <- rank_long[
  !is.na(rank_long$gene) & nzchar(rank_long$gene),
  , drop = FALSE
]
rownames(rank_long) <- NULL

n_ranked_genes <- nrow(atlas_markers)
top_decile_cutoff <- ceiling(0.10 * n_ranked_genes)

# This rule is an exploratory screen, not a confirmatory statistical test.
# Support requires a top-decile F7 rank and a better rank than the median of all
# non-F7 states. The rule is intentionally independent of cell-level p-values.
rank_candidate <- function(gene_name) {
  x <- rank_long[rank_long$gene == gene_name, , drop = FALSE]
  f7_rank <- safe_value(x$rank[x$group == target_group])
  other_ranks <- x$rank[x$group != target_group]
  best_rank <- if (nrow(x) > 0L) min(x$rank, na.rm = TRUE) else NA_real_
  best_groups <- if (is.finite(best_rank)) {
    paste(x$group[x$rank == best_rank], collapse = ";")
  } else {
    NA_character_
  }
  median_other <- if (length(other_ranks) > 0L) {
    median(other_ranks, na.rm = TRUE)
  } else {
    NA_real_
  }
  top_decile <- is.finite(f7_rank) && f7_rank <= top_decile_cutoff
  better_than_median_other <- is.finite(f7_rank) &&
    is.finite(median_other) && f7_rank < median_other

  data.frame(
    gene = gene_name,
    atlas_states_with_gene = nrow(x),
    F7_rank = f7_rank,
    F7_rank_fraction = f7_rank / n_ranked_genes,
    F6_myofibroblast_rank = safe_value(
      x$rank[x$group == "F6: Myofibroblast"]
    ),
    F6_inflammatory_myofibroblast_rank = safe_value(
      x$rank[x$group == "F6: Inflammatory myofibroblast"]
    ),
    F_fascia_rank = safe_value(x$rank[x$group == "F_Fascia"]),
    median_non_F7_rank = median_other,
    F7_rank_advantage_vs_median_non_F7 = median_other - f7_rank,
    best_rank_across_states = best_rank,
    best_rank_state = best_groups,
    F7_top_decile = top_decile,
    F7_better_than_median_non_F7 = better_than_median_other,
    exploratory_atlas_support = top_decile && better_than_median_other,
    stringsAsFactors = FALSE
  )
}

candidate_rank_rows <- lapply(shortlist$gene, rank_candidate)
candidate_ranks <- do.call(rbind, candidate_rank_rows)
candidate_results <- merge(
  shortlist[, required_step07_columns, drop = FALSE],
  candidate_ranks,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
candidate_results <- candidate_results[
  order(candidate_results$F7_rank, candidate_results$gene, na.last = TRUE),
  , drop = FALSE
]
rownames(candidate_results) <- NULL

candidate_output <- file.path(
  result_dir,
  "skin_fibroblast_atlas_Step07_candidate_marker_rank_screen_v1.csv"
)
safe_write_csv(candidate_results, candidate_output)

# -----------------------------------------------------------------------------
# 3. Summarize screening support by the frozen Step-06/07 module labels.
# -----------------------------------------------------------------------------
module_rows <- lapply(split(candidate_results, candidate_results$module), function(x) {
  data.frame(
    module = x$module[[1L]],
    candidate_genes = nrow(x),
    supported_genes = sum(x$exploratory_atlas_support %in% TRUE, na.rm = TRUE),
    support_fraction = mean(x$exploratory_atlas_support %in% TRUE, na.rm = TRUE),
    median_F7_rank = median(x$F7_rank, na.rm = TRUE),
    median_F7_rank_fraction = median(x$F7_rank_fraction, na.rm = TRUE),
    supported_gene_names = paste(
      x$gene[x$exploratory_atlas_support %in% TRUE],
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
})
module_summary <- do.call(rbind, module_rows)
module_summary <- module_summary[
  order(-module_summary$support_fraction, module_summary$median_F7_rank),
  , drop = FALSE
]
rownames(module_summary) <- NULL

module_output <- file.path(
  result_dir,
  "skin_fibroblast_atlas_Step07_module_marker_rank_summary_v1.csv"
)
safe_write_csv(module_summary, module_output)

# -----------------------------------------------------------------------------
# 4. Inspect the mechanosensor genes separately, including Step-07 cautions.
# -----------------------------------------------------------------------------
focus_genes <- c("PIEZO1", "PIEZO2", "TRPV4", "TMEM63B", "PANX1", "PKD2")
focus_results <- do.call(rbind, lapply(focus_genes, rank_candidate))
focus_results <- focus_results[
  order(focus_results$F7_rank, focus_results$gene, na.last = TRUE),
  , drop = FALSE
]
rownames(focus_results) <- NULL

focus_output <- file.path(
  result_dir,
  "skin_fibroblast_atlas_mechanosensor_marker_rank_screen_v1.csv"
)
safe_write_csv(focus_results, focus_output)

# -----------------------------------------------------------------------------
# 5. Verify that the atlas source inventory contains an independent DD dataset.
# -----------------------------------------------------------------------------
dataset_inventory <- read.csv(
  dataset_inventory_path,
  check.names = FALSE,
  strip.white = TRUE
)
accession_column <- colnames(dataset_inventory)[1L]
accessions <- trimws(as.character(dataset_inventory[[accession_column]]))
contains_prjna607098 <- "PRJNA607098" %in% accessions
contains_gse173252 <- "GSE173252" %in% accessions

inventory_audit <- data.frame(
  item = c(
    "independent_DD_source_PRJNA607098_present",
    "discovery_source_GSE173252_present",
    "webportal_has_sample_or_donor_column",
    "whole_H5AD_required_for_this_step",
    "files_over_50MB_downloaded_by_this_script"
  ),
  value = c(
    contains_prjna607098,
    contains_gse173252,
    FALSE,
    FALSE,
    FALSE
  ),
  stringsAsFactors = FALSE
)
inventory_output <- file.path(
  result_dir,
  "skin_fibroblast_atlas_validation_source_audit_v1.csv"
)
safe_write_csv(inventory_audit, inventory_output)

# -----------------------------------------------------------------------------
# 6. Evidence decision. Keep CAUTION until source-specific/donor-level validation.
# -----------------------------------------------------------------------------
n_supported <- sum(candidate_results$exploratory_atlas_support %in% TRUE, na.rm = TRUE)
n_candidates <- nrow(candidate_results)
supported_names <- paste(
  candidate_results$gene[candidate_results$exploratory_atlas_support %in% TRUE],
  collapse = ", "
)
piezo2_row <- focus_results[focus_results$gene == "PIEZO2", , drop = FALSE]
piezo2_rank <- safe_value(piezo2_row$F7_rank)

report_lines <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate",
  "- Origin Date: 2026-08-23",
  "- Verification Status: SCREENED",
  "- Version Label: skin_fibroblast_atlas_marker_rank_screen_v1",
  "",
  "## Step 08A external atlas marker-rank screen",
  "",
  "### Frozen question",
  "",
  paste0(
    "Do the ", n_candidates,
    " frozen Step-07 A/B candidates rank preferentially in the authors' ",
    "F7 fascia-like myofibroblast marker list?"
  ),
  "",
  "### Screening rule",
  "",
  paste0(
    "Exploratory support requires F7 rank <= ", top_decile_cutoff,
    " (top 10%) and a better F7 rank than the median rank across non-F7 states."
  ),
  "No cell-level p values are used and no new differential-expression test is run.",
  "",
  "### Result",
  "",
  paste0("- Supported candidates: ", n_supported, "/", n_candidates, "."),
  paste0("- Supported genes: ", supported_names, "."),
  paste0("- PIEZO2 F7 marker rank: ", piezo2_rank, " of ", n_ranked_genes, "."),
  paste0("- Independent DD source PRJNA607098 listed by the atlas: ", contains_prjna607098, "."),
  paste0("- Discovery source GSE173252 also listed by the atlas: ", contains_gse173252, "."),
  "",
  "### Evidence decision",
  "",
  "- Evidence grade remains CAUTION.",
  paste0(
    "- This marker list was calculated for atlas cell states pooled across sources; ",
    "it does not isolate PRJNA607098 from GSE173252."
  ),
  "- The public webportal object exposes dataset and cell-state labels but no sample/donor label.",
  "- Therefore this step is an external atlas screen, not donor-level independent validation.",
  "- A positive screen justifies Step 08B: source-specific PRJNA607098 expression extraction.",
  "- Formal evidence upgrading requires Step 08C donor/sample-level pseudobulk or experimental validation.",
  "",
  "### Mandatory interpretation limits",
  "",
  "- Marker rank does not establish mechanosensor activity, protein abundance or causality.",
  "- Atlas integration and cell-state annotation can induce source-dependent classification effects.",
  "- The F7 state may reflect fascia/myofibroblast identity rather than disease-specific force sensing.",
  "- The 50 MB automatic-download ceiling remains in force."
)

report_output <- file.path(
  result_dir,
  "skin_fibroblast_atlas_external_validation_decision_v1.md"
)
writeLines(report_lines, report_output, useBytes = TRUE)

message("Step 08A atlas marker-rank screen completed.")
message("Candidate results: ", candidate_output)
message("Module summary: ", module_output)
message("Mechanosensor focus: ", focus_output)
message("Source audit: ", inventory_output)
message("Evidence decision: ", report_output)

