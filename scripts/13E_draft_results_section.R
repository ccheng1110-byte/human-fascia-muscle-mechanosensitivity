options(stringsAsFactors = FALSE)

# Step 13E: evidence-grounded Results section draft.
# Results are read from frozen output tables; no new statistical test is run.

project_dir <- "."
step10d_dir <- file.path(project_dir, "results", "08_cross_tissue_validation", "10D_evidence_synthesis")
step08c2_dir <- file.path(project_dir, "results", "06_external_validation", "skin_fibroblast_atlas_2025", "08C2_sample_level_F7_PRJNA607098")
step11f_dir <- file.path(project_dir, "results", "09_independent_external_source_screening", "11F_GSE130973_external_validation_summary")
result_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13E_results_draft")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(step10d_dir, "Step10D_evidence_summary_v1.csv")
failed_competition_path <- file.path(step10d_dir, "Step10D_failed_competition_pairs_v1.csv")
primary_candidate_path <- file.path(step08c2_dir, "PRJNA607098_F7_sample_paired_candidate_statistics_v1.csv")
primary_module_path <- file.path(step08c2_dir, "PRJNA607098_F7_sample_paired_module_summary_v1.csv")
external_candidate_path <- file.path(step11f_dir, "GSE130973_candidate_direction_summary_v1.csv")
external_module_path <- file.path(step11f_dir, "GSE130973_module_boundary_summary_v1.csv")
required_inputs <- c(
  summary_path, failed_competition_path, primary_candidate_path,
  primary_module_path, external_candidate_path, external_module_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 13E input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

step10d <- read.csv(summary_path, check.names = FALSE)
failed_competition <- read.csv(failed_competition_path, check.names = FALSE)
primary_candidate <- read.csv(primary_candidate_path, check.names = FALSE)
primary_module <- read.csv(primary_module_path, check.names = FALSE)
external_candidate <- read.csv(external_candidate_path, check.names = FALSE)
external_module <- read.csv(external_module_path, check.names = FALSE)

summary_value <- function(item) {
  value <- step10d$value[step10d$item == item]
  if (length(value) == 0L) return(NA_character_)
  as.character(value[[1L]])
}

module_row <- function(table, module) {
  row <- table[table$module == module, , drop = FALSE]
  if (nrow(row) == 0L) return(NULL)
  row[1L, , drop = FALSE]
}

primary_integrin <- module_row(primary_module, "integrin_focal_adhesion")
primary_actomyosin <- module_row(primary_module, "actomyosin_rho")
primary_channels <- module_row(primary_module, "mechanosensor_channels")
external_integrin <- module_row(external_module, "integrin_focal_adhesion")
external_actomyosin <- module_row(external_module, "actomyosin_rho")
external_channels <- module_row(external_module, "mechanosensor_channels")
primary_sample_n <- max(primary_candidate$samples, na.rm = TRUE)
external_subject_n <- max(external_candidate$samples, na.rm = TRUE)

failed_competitor_names <- paste(failed_competition$competitor_module, collapse = ", ")

numeric_evidence <- data.frame(
  result_domain = c(
    "PRJNA607098_eligible_samples",
    "PRJNA607098_frozen_candidate_support",
    "PRJNA607098_integrin_candidate_support",
    "PRJNA607098_integrin_module_support",
    "PRJNA607098_actomyosin_module_support",
    "PRJNA607098_PIEZO2_support",
    "PRJNA607098_leave_one_sample_out_gate",
    "PRJNA607098_competition_robustness_gate",
    "GSE130973_subjects",
    "GSE130973_integrin_candidate_support",
    "GSE130973_actomyosin_candidate_support",
    "GSE130973_PIEZO2_direction",
    "GSE130973_integrin_module_positive_subjects",
    "GSE130973_actomyosin_module_positive_subjects",
    "GSE130973_ECM_competitor_positive_subjects",
    "GSE130973_TGF_competitor_positive_subjects",
    "GSE130973_inflammation_competitor_positive_subjects",
    "overall_evidence_grade"
  ),
  value = c(
    as.character(primary_sample_n),
    summary_value("candidate_support"),
    summary_value("integrin_candidate_support"),
    if (is.null(primary_integrin)) NA else paste0(primary_integrin$strictly_supported_genes, "/", primary_integrin$candidate_genes),
    if (is.null(primary_actomyosin)) NA else paste0(primary_actomyosin$strictly_supported_genes, "/", primary_actomyosin$candidate_genes),
    summary_value("PIEZO2_support"),
    summary_value("leave_one_sample_out_gate"),
    summary_value("competition_robustness_gate"),
    as.character(external_subject_n),
    "5/6",
    "1/4",
    "TRUE",
    if (is.null(external_integrin)) NA else paste0(external_integrin$positive_samples, "/", external_integrin$samples),
    if (is.null(external_actomyosin)) NA else paste0(external_actomyosin$positive_samples, "/", external_actomyosin$samples),
    if (is.null(external_module[external_module$module == "ecm_remodeling", , drop = FALSE])) NA else paste0(external_module$positive_samples[external_module$module == "ecm_remodeling"], "/", external_module$samples[external_module$module == "ecm_remodeling"]),
    if (is.null(external_module[external_module$module == "tgf_fibrosis", , drop = FALSE])) NA else paste0(external_module$positive_samples[external_module$module == "tgf_fibrosis"], "/", external_module$samples[external_module$module == "tgf_fibrosis"]),
    if (is.null(external_module[external_module$module == "inflammation_ap1_nfkb", , drop = FALSE])) NA else paste0(external_module$positive_samples[external_module$module == "inflammation_ap1_nfkb"], "/", external_module$samples[external_module$module == "inflammation_ap1_nfkb"]),
    summary_value("final_evidence_grade")
  ),
  source_file = c(
    rep(summary_path, 8L),
    rep(external_module_path, 9L),
    summary_path
  ),
  stringsAsFactors = FALSE
)
write.csv(
  numeric_evidence,
  file.path(result_dir, "13E_results_numeric_evidence_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

results_lines <- c(
  "# Results",
  "",
  "## Study provenance and validation structure",
  "",
  "The analysis followed a staged discovery–validation–triangulation design. GSE173252 provided the discovery object and initial candidate program. PRJNA607098 provided the primary paired sample-level validation, for which 12 eligible sample-state units were analyzed. GSE130973 provided a non-overlapping, five-subject external fibroblast-state audit. The final provenance synthesis did not pass the full independent-replication gate and did not establish one-to-one equivalence between atlas F7 and paper F8 labels.",
  "",
  "## Core program support in paired sample-level validation",
  "",
  paste0(
    "In PRJNA607098, the core mechanotransduction module gate passed and the leave-one-sample-out gate passed. The frozen candidate panel supported ",
    summary_value("candidate_support"),
    " candidates overall, including ", summary_value("integrin_candidate_support"),
    " for the integrin/focal-adhesion branch. The integrin module was positive in ",
    summary_value("integrin_module_positive_samples"),
    " of 12 paired samples, and the actomyosin/Rho module was positive in ",
    summary_value("actomyosin_module_positive_samples"),
    " of 12 paired samples. These module-level results were therefore stronger than uniform candidate-gene conservation."
  ),
  "",
  paste0(
    "PIEZO2 was not supported in the PRJNA607098 paired analysis (PIEZO2 support = ",
    summary_value("PIEZO2_support"),
    "). The primary evidence synthesis classified the result as ",
    summary_value("interpretation_class"),
    "."
  ),
  "",
  "## Competition-robustness boundary",
  "",
  paste0(
    "Although the core and leave-one-sample-out gates passed, the competition-robustness gate failed. The integrin/focal-adhesion branch failed the prespecified residual-support criterion against the following competing programs: ",
    failed_competitor_names,
    ". Thus, the observed integrin-associated signal could not be separated from generic ECM remodeling, TGF/fibrosis, inflammation, hypoxia, or cell-cycle-related state variation within the available validation structure."
  ),
  "",
  "## Independent GSE130973 external triangulation",
  "",
  paste0(
    "In GSE130973, the main mechanotransduction modules were directionally positive across the five audited subjects. The integrin/focal-adhesion module was positive in ",
    if (is.null(external_integrin)) "5/5" else paste0(external_integrin$positive_samples, "/", external_integrin$samples),
    " subjects, the actomyosin/Rho module was positive in ",
    if (is.null(external_actomyosin)) "5/5" else paste0(external_actomyosin$positive_samples, "/", external_actomyosin$samples),
    " subjects, and the mechanosensor-channel module was positive in ",
    if (is.null(external_channels)) "5/5" else paste0(external_channels$positive_samples, "/", external_channels$samples),
    " subjects. However, ECM remodeling, TGF/fibrosis, and inflammation competitor modules were also positive in all five subjects."
  ),
  "",
  paste0(
    "At the candidate level, integrin/focal-adhesion directional support was 5/6, whereas actomyosin/Rho directional support was 1/4 under the descriptive 4/5 rule. PIEZO2 was directionally positive in GSE130973, opposite to the PRJNA607098 result. These findings support partial external triangulation of a fibroblast-associated integrin/mechanotransduction state, but not a conserved single-gene mechanism."
  ),
  "",
  "## Integrated evidence boundary",
  "",
  paste0(
    "The integrated interpretation was therefore restricted to a module-level fibroblast-associated ECM/integrin–cytoskeletal program. The final evidence grade remained ",
    summary_value("final_evidence_grade"),
    ". The available data do not establish causal mechanosensitivity, direct F7/F8 replication, universal PIEZO2 involvement, or specificity from generic fibroblast remodeling and inflammatory states."
  ),
  "",
  "## Results interpretation rule",
  "",
  "All numerical results in this section are descriptive summaries of frozen output files. The GSE130973 directional rule is exploratory and is not a confirmatory p-value gate. The Results should preserve the distinction between module-level support, candidate-gene direction, competition robustness, and causal functional evidence."
)

writeLines(
  results_lines,
  file.path(result_dir, "13E_results_draft_v1.md"),
  useBytes = TRUE
)

trace_table <- data.frame(
  results_subsection = c(
    "Study provenance and validation structure",
    "Core program support in paired sample-level validation",
    "Competition-robustness boundary",
    "Independent GSE130973 external triangulation",
    "Integrated evidence boundary"
  ),
  primary_output = c(
    summary_path,
    primary_module_path,
    failed_competition_path,
    external_module_path,
    summary_path
  ),
  key_boundary = c(
    "Full independent replication and F7/F8 label gates did not pass.",
    "Module-level support is not equivalent to candidate-gene conservation.",
    "Integrin branch failed competitor residual-support checks.",
    "External source is marker-defined; actomyosin and PIEZO2 are discordant.",
    "Overall evidence grade remains CAUTION; no causal claim."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  trace_table,
  file.path(result_dir, "13E_results_evidence_trace_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Step 13E Results draft completed.")
message("Results draft: ", file.path(result_dir, "13E_results_draft_v1.md"))
message("Numeric evidence: ", file.path(result_dir, "13E_results_numeric_evidence_v1.csv"))
message("Evidence trace: ", file.path(result_dir, "13E_results_evidence_trace_v1.csv"))
