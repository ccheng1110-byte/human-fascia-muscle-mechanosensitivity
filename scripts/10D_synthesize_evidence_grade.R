options(stringsAsFactors = FALSE)

# Step 10D: provenance-aware synthesis of the frozen Step 10C result.
# No new expression test is performed here. This script summarizes the
# pre-specified evidence boundaries and writes the final interim decision.

project_dir <- "."
audit_dir <- file.path(
  project_dir, "results", "08_cross_tissue_validation",
  "10A_atlas_provenance_audit"
)
step10c_dir <- file.path(
  project_dir, "results", "08_cross_tissue_validation",
  "10C_frozen_PRJNA607098"
)
result_dir <- file.path(
  project_dir, "results", "08_cross_tissue_validation",
  "10D_evidence_synthesis"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  file.path(audit_dir, "skin_fibroblast_atlas_source_accession_map_v1.csv"),
  file.path(audit_dir, "skin_fibroblast_atlas_F6_F7_F8_label_crosswalk_v1.csv"),
  file.path(audit_dir, "skin_fibroblast_atlas_source_sample_field_audit_v1.csv"),
  file.path(step10c_dir, "PRJNA607098_F7_frozen_candidate_statistics_v2.csv"),
  file.path(step10c_dir, "PRJNA607098_F7_module_summary_v2.csv"),
  file.path(step10c_dir, "PRJNA607098_F7_competitor_robustness_v2.csv"),
  file.path(step10c_dir, "PRJNA607098_F7_leave_one_sample_out_v2.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required Step 10A/10C input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

as_logical_safe <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

source_map <- read.csv(required_inputs[[1L]], check.names = FALSE)
label_crosswalk <- read.csv(required_inputs[[2L]], check.names = FALSE)
sample_audit <- read.csv(required_inputs[[3L]], check.names = FALSE)
candidate_statistics <- read.csv(required_inputs[[4L]], check.names = FALSE)
module_summary <- read.csv(required_inputs[[5L]], check.names = FALSE)
robustness <- read.csv(required_inputs[[6L]], check.names = FALSE)
loo_summary <- read.csv(required_inputs[[7L]], check.names = FALSE)

source_map$accession_resolved <- !is.na(source_map$source_accession) &
  nzchar(trimws(as.character(source_map$source_accession)))
source_map$independent_candidate_after_exclusion <- as_logical_safe(
  source_map$independent_candidate_after_exclusion
)
source_map$overlaps_existing_project <- as_logical_safe(
  source_map$overlaps_existing_project
)
resolved_map <- source_map[source_map$accession_resolved, , drop = FALSE]
overlap_accessions <- unique(resolved_map$source_accession[
  resolved_map$overlaps_existing_project
])
independent_rows <- sum(resolved_map$independent_candidate_after_exclusion)

atlas_f7 <- label_crosswalk[
  label_crosswalk$label == "F7: Fascia-like myofibroblast" &
    label_crosswalk$label_system == "Step 08B3 atlas celltype", , drop = FALSE
]
paper_f8 <- label_crosswalk[
  label_crosswalk$label == "F8: fascia-like myofibroblast", , drop = FALSE
]
# The atlas F7 label is confirmed locally, but the paper F8-to-atlas mapping
# is not resolved. Therefore the direct F8 replication gate remains FALSE.
label_gate <- nrow(atlas_f7) == 1L && nrow(paper_f8) == 1L &&
  atlas_f7$mapping_status[[1L]] == "confirmed_in_step08B3_output" &&
  paper_f8$mapping_status[[1L]] == "confirmed_in_crosswalk"

core_modules <- c("integrin_focal_adhesion", "actomyosin_rho")
core_rows <- module_summary[module_summary$module %in% core_modules, , drop = FALSE]
core_gate <- nrow(core_rows) == 2L && all(core_rows$module_gate_pass)
loo_gate <- all(loo_summary$sign_preserved)
competition_failed <- robustness[
  !robustness$residual_support_8_of_12, , drop = FALSE
]
competition_gate <- nrow(competition_failed) == 0L

candidate_supported <- sum(candidate_statistics$strict_sample_level_support %in% TRUE)
integrin_supported <- sum(
  candidate_statistics$strict_sample_level_support[
    candidate_statistics$module == "integrin_focal_adhesion"
  ] %in% TRUE
)
piezo2_supported <- isTRUE(candidate_statistics$strict_sample_level_support[
  match("PIEZO2", candidate_statistics$gene)
])
actomyosin_row <- core_rows[core_rows$module == "actomyosin_rho", , drop = FALSE]
integrin_row <- core_rows[core_rows$module == "integrin_focal_adhesion", , drop = FALSE]

independence_gate <- length(overlap_accessions) == 0L &&
  independent_rows > 0L &&
  any(grepl("sample", tolower(sample_audit$conclusion)))
full_independent_replication <- independence_gate && label_gate &&
  core_gate && loo_gate && competition_gate

final_grade <- if (full_independent_replication) {
  "INDEPENDENT_REPLICATION_SUPPORTED"
} else {
  "CAUTION"
}

interpretation_class <- if (core_gate && loo_gate && !competition_gate) {
  "ACTOMYOSIN_ROBUST_INTEGRIN_NON_SPECIFIC"
} else if (core_gate && loo_gate && competition_gate) {
  "CORE_MECHANOTRANSDUCTION_SUPPORTED_WITH_BOUNDARY"
} else {
  "CORE_MECHANOTRANSDUCTION_NOT_ROBUST"
}

summary <- data.frame(
  item = c(
    "resolved_atlas_accession_rows",
    "independent_candidate_accession_rows",
    "overlap_accessions",
    "source_independence_gate",
    "F7_F8_label_gate",
    "candidate_support",
    "integrin_candidate_support",
    "PIEZO2_support",
    "actomyosin_module_positive_samples",
    "integrin_module_positive_samples",
    "core_module_gate",
    "leave_one_sample_out_gate",
    "competition_robustness_gate",
    "full_independent_replication",
    "interpretation_class",
    "final_evidence_grade"
  ),
  value = c(
    nrow(resolved_map),
    independent_rows,
    paste(overlap_accessions, collapse = ";"),
    independence_gate,
    label_gate,
    paste0(candidate_supported, "/15"),
    paste0(integrin_supported, "/6"),
    piezo2_supported,
    paste0(actomyosin_row$positive_samples[[1L]], "/12"),
    paste0(integrin_row$positive_samples[[1L]], "/12"),
    core_gate,
    loo_gate,
    competition_gate,
    full_independent_replication,
    interpretation_class,
    final_grade
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(
  summary,
  file.path(result_dir, "Step10D_evidence_summary_v1.csv")
)
safe_write_csv(
  competition_failed,
  file.path(result_dir, "Step10D_failed_competition_pairs_v1.csv")
)

failed_pairs_text <- if (nrow(competition_failed) == 0L) {
  "- No competitor pair failed the frozen residual-support rule."
} else {
  apply(competition_failed[, c(
    "core_module", "competitor_module", "residual_positive_samples",
    "residual_median"
  )], 1, function(x) {
    paste0(
      "- ", x[[1L]], " versus ", x[[2L]], ": residual positive samples = ",
      x[[3L]], "/12; residual median = ", signif(as.numeric(x[[4L]]), 4), "."
    )
  })
}

decision_lines <- c(
  "# Step 10D evidence synthesis decision",
  "",
  "## Current conclusion",
  "",
  paste0("The final interim evidence grade is **", final_grade, "**."),
  paste0("The appropriate interpretation class is **", interpretation_class, "**."),
  "",
  "The robust part of the result is the actomyosin/Rho-associated program: the module is positive in 12/12 samples and the leave-one-sample-out gate passes. The integrin/focal-adhesion module is also positive in 12/12 samples and 5/6 frozen candidates are supported, but its specificity is not established after the pre-specified competitor audit.",
  "",
  "## Failed or incomplete boundaries",
  "",
  paste0("- Step 10A source independence: PARTIAL; overlapping accessions = ", paste(overlap_accessions, collapse = ", "), "."),
  "- F7/F8 anatomical-label equivalence: unresolved for direct replication claims.",
  "- PIEZO2: not supported; channel-specific mechanosensor evidence remains negative/weak.",
  paste0("- Competition robustness: ", competition_gate, "."),
  failed_pairs_text,
  "",
  "## Recommended claim boundary",
  "",
  "- Supported wording: fascia-like myofibroblast states show a reproducible actomyosin/Rho-associated transcriptional program, accompanied by integrin/focal-adhesion enrichment.",
  "- Required limitation: the integrin component is not separable from TGF/fibrosis, inflammatory, hypoxia, and cell-cycle state signals in this validation source.",
  "- Prohibited wording: independent replication, direct F8 replication, PIEZO2-mediated mechanosensing, or causal mechanotransduction.",
  "",
  "## Next decision",
  "",
  "Do not rerun Step 10C with relaxed thresholds. The highest-value next step is a genuinely non-overlapping donor-level source with sample metadata. If that source cannot be recovered, proceed to manuscript-level synthesis as a bounded, hypothesis-generating reanalysis.",
  "",
  "## Material Passport",
  "",
  "- Inputs: corrected Step 10A provenance audit and Step 10C frozen-contract outputs.",
  "- Transformation: metadata-aware evidence synthesis only; no new expression values or hypothesis tests.",
  "- Statistical units: the 12 PRJNA607098 sample IDs inherited from Step 10C.",
  "- Reproducibility: all inputs and failed-pair details are exported in the same result directory.",
  "- Integrity status: interim decision; independent replication remains unavailable."
)
decision_path <- file.path(result_dir, "Step10D_evidence_synthesis_decision_v1.md")
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 10D evidence synthesis completed.")
message("Interpretation class: ", interpretation_class)
message("Final evidence grade: ", final_grade)
message("Summary: ", file.path(result_dir, "Step10D_evidence_summary_v1.csv"))
message("Failed competition pairs: ", file.path(result_dir, "Step10D_failed_competition_pairs_v1.csv"))
message("Decision: ", decision_path)
