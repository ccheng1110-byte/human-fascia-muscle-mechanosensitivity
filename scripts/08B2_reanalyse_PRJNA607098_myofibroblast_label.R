options(stringsAsFactors = FALSE)

# Step 08B2: label-corrected source-specific validation in PRJNA607098.
#
# Step 08B v1 completed technically, but it looked for the atlas-wide label
# "F7: Fascia-like myofibroblast" inside PRJNA607098. That source instead uses
# the harmonized source-level label "Myofibroblast" (1,714 cells), so v1 had
# zero target cells and was not biologically interpretable. This correction
# reuses the already extracted source/state/gene summary and performs no
# network access or additional download. The frozen 15-gene panel, contrasts,
# support rule and 08B -> 08C gate remain unchanged.

project_dir <- "."
step07_path <- file.path(
  project_dir, "results", "05_specificity_robustness", "GSE173252",
  "GSE173252_DD_candidate_specificity_summary_v1.csv"
)
step08b_dir <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025", "08B_source_specific_PRJNA607098"
)
expression_summary_path <- file.path(
  step08b_dir,
  "skin_fibroblast_atlas_source_state_gene_expression_summary_v1.csv"
)
result_dir <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025",
  "08B2_corrected_myofibroblast_PRJNA607098"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(step07_path, expression_summary_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

step07 <- read.csv(step07_path, check.names = FALSE)
shortlist <- step07[
  grepl("^[AB]_", step07$step06_candidate_tier) &
    step07$specificity_class == "robust_Myofib_over_all_three_identities",
  , drop = FALSE
]
shortlist <- shortlist[!duplicated(shortlist$gene), , drop = FALSE]
if (nrow(shortlist) != 15L) {
  stop(
    "Expected 15 frozen Step-07 candidates but found ", nrow(shortlist),
    ". Review Step-07 outputs before continuing."
  )
}

expression_summary <- read.csv(expression_summary_path, check.names = FALSE)
source_id <- "PRJNA607098"
target_state <- "Myofibroblast"
source_rows <- expression_summary[
  expression_summary$GSE == source_id,
  , drop = FALSE
]

if (nrow(source_rows) == 0L) {
  stop("No expression-summary rows were found for ", source_id, ".")
}
available_states <- sort(unique(source_rows$cell_state))
if (!target_state %in% available_states) {
  stop(
    "Corrected target state was not found: ", target_state,
    ". Available states: ", paste(available_states, collapse = "; ")
  )
}

target_cell_counts <- unique(
  source_rows$cells[source_rows$cell_state == target_state]
)
if (length(target_cell_counts) != 1L || target_cell_counts < 20L) {
  stop(
    "Unexpected target-state cell count for ", target_state, ": ",
    paste(target_cell_counts, collapse = "; ")
  )
}

contrast_rows <- lapply(shortlist$gene, function(gene_name) {
  gene_rows <- source_rows[source_rows$gene == gene_name, , drop = FALSE]
  target <- gene_rows[gene_rows$cell_state == target_state, , drop = FALSE]
  other <- gene_rows[gene_rows$cell_state != target_state, , drop = FALSE]

  if (nrow(target) != 1L) {
    stop("Expected exactly one target-state row for gene: ", gene_name)
  }
  if (nrow(other) < 2L) {
    stop("Fewer than two comparison states were found for gene: ", gene_name)
  }

  other_means <- as.numeric(other$mean_expression)
  other_prevalences <- as.numeric(other$percent_cells_expressing)
  target_mean <- as.numeric(target$mean_expression)
  target_prevalence <- as.numeric(target$percent_cells_expressing)
  median_other_mean <- median(other_means, na.rm = TRUE)
  median_other_prevalence <- median(other_prevalences, na.rm = TRUE)
  best_other_mean <- max(other_means, na.rm = TRUE)
  mean_difference <- target_mean - median_other_mean
  prevalence_difference <- target_prevalence - median_other_prevalence
  highest_mean <- isTRUE(target_mean > best_other_mean)
  descriptive_support <- isTRUE(
    as.numeric(target$cells) >= 20L &&
      length(other_means) >= 2L &&
      mean_difference > 0 &&
      prevalence_difference > 0
  )

  data.frame(
    GSE = source_id,
    gene = gene_name,
    target_cell_state = target_state,
    target_cells = as.numeric(target$cells),
    comparison_states = nrow(other),
    target_mean_expression = target_mean,
    median_other_state_mean_expression = median_other_mean,
    target_minus_median_other_state_mean = mean_difference,
    target_minus_best_other_state_mean = target_mean - best_other_mean,
    target_percent_cells_expressing = target_prevalence,
    median_other_state_percent_expressing = median_other_prevalence,
    target_minus_median_other_state_percent_expressing = prevalence_difference,
    target_highest_mean_across_states = highest_mean,
    source_specific_descriptive_support = descriptive_support,
    stringsAsFactors = FALSE
  )
})
contrasts <- do.call(rbind, contrast_rows)

candidate_results <- merge(
  shortlist[, c(
    "module", "gene", "step06_candidate_tier",
    "myofib_percent_cells_expressing",
    "minimum_difference_over_all_donors_and_identities"
  ), drop = FALSE],
  contrasts,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
candidate_results <- candidate_results[
  order(
    -candidate_results$source_specific_descriptive_support,
    -candidate_results$target_minus_median_other_state_mean,
    candidate_results$gene
  ),
  , drop = FALSE
]
rownames(candidate_results) <- NULL

candidate_output <- file.path(
  result_dir,
  "PRJNA607098_myofibroblast_frozen_candidate_validation_v2.csv"
)
safe_write_csv(candidate_results, candidate_output)

module_rows <- lapply(split(candidate_results, candidate_results$module), function(x) {
  data.frame(
    module = x$module[[1L]],
    candidate_genes = nrow(x),
    supported_genes = sum(
      x$source_specific_descriptive_support %in% TRUE,
      na.rm = TRUE
    ),
    support_fraction = mean(
      x$source_specific_descriptive_support %in% TRUE,
      na.rm = TRUE
    ),
    median_target_minus_other_mean = median(
      x$target_minus_median_other_state_mean,
      na.rm = TRUE
    ),
    supported_gene_names = paste(
      x$gene[x$source_specific_descriptive_support %in% TRUE],
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
})
module_summary <- do.call(rbind, module_rows)
module_summary <- module_summary[
  order(-module_summary$support_fraction, module_summary$module),
  , drop = FALSE
]
rownames(module_summary) <- NULL

module_output <- file.path(
  result_dir,
  "PRJNA607098_myofibroblast_module_validation_v2.csv"
)
safe_write_csv(module_summary, module_output)

state_inventory <- unique(source_rows[, c("cell_state", "cells"), drop = FALSE])
state_inventory <- state_inventory[order(-state_inventory$cells), , drop = FALSE]
rownames(state_inventory) <- NULL
state_inventory_output <- file.path(
  result_dir,
  "PRJNA607098_available_cell_state_inventory_v2.csv"
)
safe_write_csv(state_inventory, state_inventory_output)

# Preserve the gate that was specified before the invalid v1 target-label result.
n_supported <- sum(
  candidate_results$source_specific_descriptive_support %in% TRUE,
  na.rm = TRUE
)
integrin_rows <- candidate_results[
  candidate_results$module == "integrin_focal_adhesion",
  , drop = FALSE
]
n_integrin_supported <- sum(
  integrin_rows$source_specific_descriptive_support %in% TRUE,
  na.rm = TRUE
)
piezo2_rows <- candidate_results[candidate_results$gene == "PIEZO2", , drop = FALSE]
piezo2_supported <- nrow(piezo2_rows) == 1L &&
  isTRUE(piezo2_rows$source_specific_descriptive_support[[1L]])

gate_pass <- n_supported >= 8L &&
  n_integrin_supported >= 4L &&
  piezo2_supported

report_lines <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate",
  "- Origin Date: 2026-08-23",
  "- Verification Status: ANALYZED",
  "- Version Label: PRJNA607098_label_corrected_myofibroblast_v2",
  "",
  "## Step 08B2 label-corrected PRJNA607098 validation",
  "",
  "### Correction audit",
  "",
  "- Step 08B v1 completed technically but was not biologically interpretable.",
  "- V1 searched for `F7: Fascia-like myofibroblast` in PRJNA607098 and found zero cells.",
  paste0(
    "- PRJNA607098 uses the source-level label `Myofibroblast` (n = ",
    target_cell_counts, ")."
  ),
  "- V2 changes only the label mapping; the frozen genes, contrast rule and gate are unchanged.",
  "- V1 is retained as an implementation-audit artifact and must not be cited as a negative biological result.",
  "",
  "### Pre-specified gate",
  "",
  "Proceed to donor/sample-level Step 08C only if all conditions are met:",
  "",
  "- At least 8 of 15 frozen candidates show descriptive support.",
  "- At least 4 of 6 integrin/focal-adhesion candidates show support.",
  "- PIEZO2 shows support.",
  "",
  "### Results",
  "",
  paste0("- Frozen candidates supported: ", n_supported, "/15."),
  paste0(
    "- Integrin/focal-adhesion candidates supported: ",
    n_integrin_supported, "/6."
  ),
  paste0("- PIEZO2 supported: ", piezo2_supported, "."),
  paste0("- Step 08C gate passed: ", gate_pass, "."),
  "",
  "### Evidence boundary",
  "",
  "- This is a source-specific descriptive cell-state contrast, not donor-level inference.",
  "- The public atlas lacks usable PRJNA607098 donor labels; cells are not treated as independent replicates.",
  "- No p values are calculated.",
  "- Evidence grade remains CAUTION even if the gate passes.",
  "- Mechanosensor activity, protein abundance and causality remain untested.",
  "",
  "### Statistical fallacy scan (11/11 checked)",
  "",
  "- Simpson's paradox: not testable without donor strata; retained as a caution.",
  "- Ecological fallacy: avoided; cell-state results are not generalized to patients.",
  "- Berkson's paradox: possible atlas selection bias; caution.",
  "- Collider bias: no covariate-adjusted causal model was fitted.",
  "- Base-rate neglect: target and comparison-state cell counts are retained.",
  "- Regression to the mean: not applicable to this cross-sectional screen.",
  "- Survivorship bias: atlas QC may exclude low-quality cells; caution.",
  "- Look-elsewhere effect: reduced by the frozen 15-gene panel.",
  "- Garden of forking paths: label correction preserves the original genes, contrasts and gate.",
  "- Correlation is not causation: causal claims are prohibited.",
  "- Reverse causality: state-expression direction cannot establish mechanism."
)

report_output <- file.path(
  result_dir,
  "PRJNA607098_myofibroblast_validation_decision_v2.md"
)
writeLines(report_lines, report_output, useBytes = TRUE)

message("Step 08B2 label-corrected reanalysis completed.")
message("Target state: ", target_state, " (", target_cell_counts, " cells)")
message("Frozen candidates supported: ", n_supported, "/15")
message("Integrin/focal-adhesion supported: ", n_integrin_supported, "/6")
message("PIEZO2 supported: ", piezo2_supported)
message("Step 08C gate passed: ", gate_pass)
message("Candidate results: ", candidate_output)
message("Module summary: ", module_output)
message("State inventory: ", state_inventory_output)
message("Evidence decision: ", report_output)
