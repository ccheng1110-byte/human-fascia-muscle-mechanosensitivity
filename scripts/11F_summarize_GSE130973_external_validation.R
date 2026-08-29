options(stringsAsFactors = FALSE)

# Step 11F: summarize the GSE130973 exploratory result and compare it with the
# prior PRJNA607098 result. No new expression values are computed here.

project_dir <- "."
step11e_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11E_GSE130973_frozen_program_audit"
)
step08c2_path <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025", "08C2_sample_level_F7_PRJNA607098",
  "PRJNA607098_F7_sample_paired_candidate_statistics_v1.csv"
)
candidate_panel_path <- file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
)
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11F_GSE130973_external_validation_summary"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  file.path(step11e_dir, "GSE130973_fibroblast_state_gene_contrasts_v1.csv"),
  file.path(step11e_dir, "GSE130973_fibroblast_state_module_summary_v1.csv"),
  step08c2_path,
  candidate_panel_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 11F input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

contrasts <- read.csv(required_inputs[[1L]], check.names = FALSE)
module_summary <- read.csv(required_inputs[[2L]], check.names = FALSE)
prior_statistics <- read.csv(required_inputs[[3L]], check.names = FALSE)
candidate_panel <- read.csv(required_inputs[[4L]], check.names = FALSE)

candidate_contrasts <- contrasts[
  contrasts$gene %in% candidate_panel$gene,
  , drop = FALSE
]
candidate_rows <- do.call(rbind, lapply(split(candidate_contrasts, candidate_contrasts$gene), function(x) {
  module <- candidate_panel$module[match(x$gene[[1L]], candidate_panel$gene)]
  data.frame(
    gene = x$gene[[1L]],
    module = module,
    samples = nrow(x),
    positive_samples = sum(x$mean_difference > 0),
    negative_samples = sum(x$mean_difference < 0),
    median_difference = median(x$mean_difference),
    IQR_difference = IQR(x$mean_difference),
    descriptive_directional_support_4_of_5 = sum(x$mean_difference > 0) >= 4L,
    stringsAsFactors = FALSE
  )
}))
candidate_rows <- candidate_rows[
  match(candidate_panel$gene, candidate_rows$gene), , drop = FALSE
]
safe_write_csv(
  candidate_rows,
  file.path(result_dir, "GSE130973_candidate_direction_summary_v1.csv")
)

cross_source <- merge(
  candidate_rows,
  prior_statistics[, c(
    "gene", "module", "positive_samples", "negative_samples",
    "median_F7_minus_nonF7_mean", "strict_sample_level_support"
  )],
  by = c("gene", "module"),
  all = TRUE,
  suffixes = c("_GSE130973", "_PRJNA607098")
)
cross_source <- cross_source[
  match(candidate_panel$gene, cross_source$gene), , drop = FALSE
]
safe_write_csv(
  cross_source,
  file.path(result_dir, "GSE130973_vs_PRJNA607098_candidate_crosswalk_v1.csv")
)

module_summary$positive_5_of_5 <- module_summary$positive_samples >= 5L
safe_write_csv(
  module_summary,
  file.path(result_dir, "GSE130973_module_boundary_summary_v1.csv")
)

integrin_rows <- candidate_rows[candidate_rows$module == "integrin_focal_adhesion", ]
actomyosin_rows <- candidate_rows[candidate_rows$module == "actomyosin_rho", ]
piezo2_row <- candidate_rows[candidate_rows$gene == "PIEZO2", , drop = FALSE]
integrin_directional <- sum(integrin_rows$descriptive_directional_support_4_of_5)
actomyosin_directional <- sum(actomyosin_rows$descriptive_directional_support_4_of_5)
piezo2_directional <- nrow(piezo2_row) == 1L &&
  isTRUE(piezo2_row$descriptive_directional_support_4_of_5[[1L]])

interpretation_class <- if (
  integrin_directional >= 4L && actomyosin_directional < 3L
) {
  "INTEGRIN_PARTIAL_ACTOMYOSIN_CANDIDATE_DISCORDANCE"
} else if (integrin_directional >= 4L && actomyosin_directional >= 3L) {
  "MULTI_MODULE_DIRECTIONAL_SUPPORT"
} else {
  "WEAK_OR_INCONSISTENT_EXTERNAL_SUPPORT"
}

decision_lines <- c(
  "# Step 11F GSE130973 external validation boundary",
  "",
  "## Summary",
  "",
  paste0("- Interpretation class: **", interpretation_class, "**."),
  paste0("- Integrin/focal-adhesion candidate directional support: ", integrin_directional, "/6 using a descriptive 4/5 rule."),
  paste0("- Actomyosin/Rho candidate directional support: ", actomyosin_directional, "/4 using a descriptive 4/5 rule."),
  paste0("- PIEZO2 directional support: ", piezo2_directional, "."),
  "- Full module scores are positive for the main modules, but ECM, TGF/fibrosis, and inflammation modules are also positive.",
  "",
  "## Cross-source interpretation",
  "",
  "GSE130973 provides a non-overlapping five-subject skin aging source with a reproducible subject field. It supports the presence of a fibroblast-associated integrin/mechanotransduction state, but the frozen actomyosin candidate genes are discordant and the channel-specific PIEZO2 direction differs from PRJNA607098. The result is therefore partial cross-source support, not independent replication of F7/F8 or proof of a single conserved mechanosensor mechanism.",
  "",
  "## Evidence boundary",
  "",
  "- Do not pool GSE130973 with PRJNA607098 as if they were the same phenotype.",
  "- Do not treat descriptive 4/5 directional support as a confirmatory p-value gate.",
  "- Keep the final project evidence grade at CAUTION.",
  "- The strongest common interpretation is a fibroblast-associated ECM/integrin/mechanotransduction state; actomyosin candidate-level specificity remains unresolved.",
  "",
  "## Material Passport",
  "",
  "- Inputs: Step 11E GSE130973 contrasts, module summary, and prior Step 08C2 candidate statistics.",
  "- Transformation: cross-source descriptive direction and support summary only; no new expression test.",
  "- Statistical units: five GSE130973 subjects; previous source remains 12 PRJNA607098 samples.",
  "- Integrity status: exploratory external triangulation; no evidence-grade upgrade."
)
decision_path <- file.path(
  result_dir, "GSE130973_external_validation_boundary_decision_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11F GSE130973 external validation summary completed.")
message("Interpretation class: ", interpretation_class)
message("Integrin directional support: ", integrin_directional, "/6")
message("Actomyosin directional support: ", actomyosin_directional, "/4")
message("PIEZO2 directional support: ", piezo2_directional)
message("Decision: ", decision_path)
