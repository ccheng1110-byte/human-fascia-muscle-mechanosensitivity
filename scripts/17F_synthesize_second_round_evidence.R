options(stringsAsFactors = FALSE)

# Step 17F: second-round evidence resynthesis.
# This script reads completed S1/S2/S3 outputs only. It does not download data,
# pool studies, change the frozen panel, or upgrade the overall evidence grade.

project_dir <- "."
round_dir <- file.path(project_dir, "results", "14_second_round_computational_strengthening")
s1_dir <- file.path(round_dir, "17C_S1_coexpression_level_analysis")
s2_dir <- file.path(round_dir, "17E_S2_GSE338388_expression_analysis")
s3_dir <- file.path(round_dir, "17H_S3_stiffness_analyses")
result_dir <- file.path(round_dir, "17F_evidence_resynthesis")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

read_required_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required evidence file: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
safe_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
fmt_num <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x[[1L]]) || !is.finite(as.numeric(x[[1L]]))) return("NA")
  formatC(as.numeric(x[[1L]]), format = "g", digits = max(4L, digits + 1L))
}
fmt_int <- function(x) {
  if (length(x) == 0L || is.na(x[[1L]])) return("NA")
  format(as.integer(x[[1L]]), trim = TRUE, scientific = FALSE)
}
first_value <- function(df, condition, column, default = NA) {
  hit <- which(condition)
  if (length(hit) == 0L || !column %in% names(df)) return(default)
  df[[column]][hit[[1L]]]
}
module_row <- function(df, module, contrast = NULL) {
  condition <- df$module == module
  if (!is.null(contrast) && "contrast" %in% names(df)) condition <- condition & df$contrast == contrast
  hit <- which(condition)
  if (length(hit) == 0L) return(df[0, , drop = FALSE])
  df[hit[[1L]], , drop = FALSE]
}
candidate_coverage <- function(df) {
  present_candidates <- c(
    "present_in_matrix", "present_after_variance_filter",
    "included_after_expression_filter"
  )
  present_col <- present_candidates[present_candidates %in% names(df)]
  if (!"gene" %in% names(df) || length(present_col) == 0L) return(NA_character_)
  genes <- unique(df$gene)
  present <- tapply(as.logical(df[[present_col[[1L]]]]), df$gene, any)
  paste0(sum(present), "/", length(genes))
}

s1_source <- read_required_csv(file.path(s1_dir, "Step17C_S1_source_level_summary_v1.csv"))
s1_pairs <- read_required_csv(file.path(s1_dir, "Step17C_S1_pair_level_statistics_v1.csv"))
s2_modules <- read_required_csv(file.path(s2_dir, "GSE338388_S2_module_camera_statistics_v1.csv"))
s2_candidates <- read_required_csv(file.path(s2_dir, "GSE338388_S2_frozen_candidate_statistics_v1.csv"))
s3_summary <- read_required_csv(file.path(s3_dir, "Step17H_S3_stiffness_analysis_summary_v1.csv"))
s3_123_modules <- read_required_csv(file.path(s3_dir, "GSE123100_S3_module_dose_response_statistics_v1.csv"))
s3_123_candidates <- read_required_csv(file.path(s3_dir, "GSE123100_S3_candidate_dose_response_statistics_v1.csv"))
s3_276_modules <- read_required_csv(file.path(s3_dir, "GSE276045_S3_module_camera_statistics_v1.csv"))
s3_276_candidates <- read_required_csv(file.path(s3_dir, "GSE276045_S3_frozen_candidate_statistics_v1.csv"))

# ---- Extract stable, auditable summary values ----
prjna_s1 <- s1_source[s1_source$source == "PRJNA607098", , drop = FALSE]
gse130973_s1 <- s1_source[s1_source$source == "GSE130973", , drop = FALSE]
prjna_pairs_supported <- first_value(prjna_s1, rep(TRUE, nrow(prjna_s1)), "descriptively_supported_pairs", 0)
prjna_pairs_tested <- first_value(prjna_s1, rep(TRUE, nrow(prjna_s1)), "pairs_tested", 0)
gse130973_pairs_supported <- first_value(gse130973_s1, rep(TRUE, nrow(gse130973_s1)), "descriptively_supported_pairs", 0)
gse130973_pairs_tested <- first_value(gse130973_s1, rep(TRUE, nrow(gse130973_s1)), "pairs_tested", 0)

s2_candidate_coverage <- candidate_coverage(s2_candidates)
s2_hippo_tead <- module_row(s2_modules, "hippo_yap_taz", "tead_main_effect")
s2_tgf_tgfb <- module_row(s2_modules, "tgf_fibrosis", "tgfb_main_effect")
s2_acto_tead <- module_row(s2_modules, "actomyosin_rho", "tead_main_effect")

s3_123_coverage <- candidate_coverage(s3_123_candidates)
s3_276_coverage <- candidate_coverage(s3_276_candidates)
s3_276_acto_wt <- module_row(s3_276_modules, "actomyosin_rho", "stiffness_slope_WT")
s3_276_acto_htert <- module_row(s3_276_modules, "actomyosin_rho", "stiffness_slope_hTERT")
s3_276_cycle_wt <- module_row(s3_276_modules, "cell_cycle", "stiffness_slope_WT")
s3_276_ecm_htert <- module_row(s3_276_modules, "ecm_remodeling", "stiffness_slope_hTERT")
s3_276_integrin_wt <- module_row(s3_276_modules, "integrin_focal_adhesion", "stiffness_slope_WT")
s3_276_integrin_htert <- module_row(s3_276_modules, "integrin_focal_adhesion", "stiffness_slope_hTERT")

# ---- Evidence-domain summary ----
domain_summary <- data.frame(
  evidence_domain = c(
    "S1_single_cell_coexpression_level",
    "S2_TGFB_TEAD_regulatory_axis",
    "S3_cross_tissue_stiffness_form",
    "S3_actomyosin_directional_reproducibility",
    "competition_specificity",
    "fascia_specific_causality",
    "overall_evidence_grade"
  ),
  status = c(
    "SOURCE_DEPENDENT_MIXED",
    "REGULATORY_AXIS_CONCORDANT",
    "PARTIAL_DESCRIPTIVE_SUPPORT",
    "CROSS_TISSUE_ACTOMYOSIN_PLAUSIBILITY",
    "NOT_RESOLVED",
    "NOT_TESTED",
    "CAUTION"
  ),
  evidence_basis = c(
    paste0(
      "PRJNA607098 descriptively supported ", prjna_pairs_supported, "/", prjna_pairs_tested,
      " tested pairs, whereas GSE130973 supported ", gse130973_pairs_supported, "/", gse130973_pairs_tested,
      "; both remain co-expression-level evidence only."
    ),
    paste0(
      "GSE338388 has ", s2_candidate_coverage,
      " frozen-candidate coverage; the TEAD and TGF-beta exposure axes show their prespecified module responses."
    ),
    paste0(
      "GSE123100 provides ", s3_123_coverage,
      " candidate coverage and descriptive HTM dose-response form; GSE276045 provides ", s3_276_coverage,
      " coverage in the WI-38 stiffness model."
    ),
    paste0(
      "In GSE276045, actomyosin camera FDR is ", fmt_num(s3_276_acto_wt$camera_fdr),
      " in WT and ", fmt_num(s3_276_acto_htert$camera_fdr),
      " in hTERT for stiffness slopes."
    ),
    paste0(
      "GSE276045 cell-cycle and ECM competitors also respond (cell-cycle WT FDR ",
      fmt_num(s3_276_cycle_wt$camera_fdr), "; ECM hTERT FDR ", fmt_num(s3_276_ecm_htert$camera_fdr),")."
    ),
    "No current public-data analysis establishes controlled mechanical causality, cell-intrinsic mechanism, or donor-replicated functional response in fascia.",
    "The second-round computational triangulation narrows the interpretation boundary but does not satisfy the prespecified functional upgrade requirement."
  ),
  permitted_interpretation = c(
    "Selected sources show source-dependent co-expression-level evidence.",
    "TGF-beta/SMAD and TEAD-related regulatory-axis cross-validation.",
    "Cross-tissue stiffness-associated form and dose-response plausibility.",
    "Actomyosin/Rho is a prioritized pathway-level candidate with cross-tissue directional plausibility.",
    "Specificity from generic ECM, fibrosis, and proliferation states remains unresolved.",
    "Causal mechanosensitivity remains untested.",
    "Evidence-bounded mechanistic hypothesis with computational triangulation; retain CAUTION."
  ),
  prohibited_overclaim = c(
    "Do not call this cell-intrinsic mechanism or a universal co-expression signature.",
    "Do not infer mechanical causality or YAP/TAZ-independent mechanical driving.",
    "Do not call these sources fascia replication.",
    "Do not call actomyosin specificity established while competitor programs are active.",
    "Do not claim fibrosis/inflammation/proliferation independence.",
    "Do not claim causal mechanosensitivity or human fascia functional validation.",
    "Do not upgrade to a non-CAUTION grade from computational results alone."
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(domain_summary, file.path(result_dir, "Step17F_second_round_evidence_domain_summary_v1.csv"))

# ---- Revised claim contract for manuscript and future experiments ----
claim_contract <- data.frame(
  claim_id = paste0("R", sprintf("%02d", 1:7)),
  claim = c(
    "A multicomponent mechanotransduction-associated program is observed in selected human fibroblast-enriched sources.",
    "S1 supports co-expression-level evidence that is dependent on source and sampling context.",
    "GSE338388 supports a TGF-beta/TEAD regulatory-axis cross-validation, not mechanical causality.",
    "Cross-tissue stiffness data provide plausibility for an actomyosin/Rho pathway component.",
    "Integrin/focal-adhesion specificity remains unresolved from ECM/TGF and proliferation competition.",
    "PIEZO2 remains context-dependent and is not a universal required marker.",
    "Causal mechanosensitivity requires donor-replicated functional perturbation and remains untested."
  ),
  evidence_status = c(
    "SUPPORTED_WITH_BOUNDARY",
    "SOURCE_DEPENDENT_MIXED",
    "REGULATORY_AXIS_CROSS_VALIDATED",
    "ACTOMYOSIN_PLAUSIBILITY_WITH_COMPETITION",
    "SPECIFICITY_UNRESOLVED",
    "SOURCE_DEPENDENT",
    "NOT_TESTED"
  ),
  manuscript_language = c(
    "associated with; supports a bounded hypothesis",
    "co-expression-level evidence; source-dependent",
    "regulatory-axis cross-validation",
    "prioritized pathway-level candidate; cross-tissue plausibility",
    "specificity remains unresolved",
    "context-dependent secondary candidate",
    "not established; requires functional validation"
  ),
  required_upgrade = c(
    "At least one orthogonal functional readout with donor replication.",
    "Sample/donor-level confirmation with detection-rate and composition controls.",
    "Controlled mechanical perturbation linked to the regulatory axis.",
    "Prespecified Rho/actomyosin perturbation changes a mechanical phenotype.",
    "Matched ECM/TGF/inflammation/proliferation controls.",
    "Reproducible context-specific functional evidence.",
    "At least 3 independent donors, functional phenotype, pathway perturbation and orthogonal convergence."
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(claim_contract, file.path(result_dir, "Step17F_revised_claim_contract_v1.csv"))

# ---- S3 module cross-source audit table ----
s3_123_module_audit <- data.frame(
  source = "GSE123100", contrast = "dose_response",
  module = s3_123_modules$module, role = s3_123_modules$role,
  genes_present = s3_123_modules$genes_present, genes_total = s3_123_modules$genes_total,
  direction = ifelse(s3_123_modules$effect_per_log10_kPa > 0, "positive", "negative"),
  nominal_p = s3_123_modules$p_value, multiple_test_value = s3_123_modules$q_module,
  interpretation_boundary = "descriptive HTM dose-response form; not fascia replication",
  stringsAsFactors = FALSE
)
s3_276_module_audit <- data.frame(
  source = "GSE276045", contrast = s3_276_modules$contrast,
  module = s3_276_modules$module, role = s3_276_modules$role,
  genes_present = s3_276_modules$genes_present_after_filter,
  genes_total = s3_276_modules$genes_total,
  direction = s3_276_modules$direction,
  nominal_p = s3_276_modules$camera_p_value,
  multiple_test_value = s3_276_modules$camera_fdr,
  interpretation_boundary = "WI-38 cross-tissue stiffness cross-check; proliferation-independent effect not estimable",
  stringsAsFactors = FALSE
)
safe_write_csv(
  rbind(s3_123_module_audit, s3_276_module_audit),
  file.path(result_dir, "Step17F_S3_module_cross_source_audit_v1.csv")
)

# ---- Functional validation priority matrix ----
priority_matrix <- data.frame(
  priority = 1:5,
  branch = c("actomyosin_Rho", "integrin_focal_adhesion", "competition_controls", "PIEZO2", "cell_intrinsic_coexpression"),
  rationale = c(
    "Cross-tissue stiffness slopes show reproducible actomyosin direction, but this must be linked to function.",
    "Integrin branch remains prioritized but its specificity is not resolved by current computational evidence.",
    "Cell-cycle/ECM competitors remain active and must be measured in the same experiment.",
    "PIEZO2 is source-dependent; test as a secondary branch rather than a required universal mechanism.",
    "Current S1 supports co-expression level only; donor-level functional design is required for mechanism claims."
  ),
  minimum_design = c(
    "At least 3 independent donors; controlled stiffness/tension perturbation; contractility readout; Rho/ROCK perturbation.",
    "Integrin perturbation with matched ECM/TGF state controls and rescue or orthogonal focal-adhesion readout.",
    "Measure proliferation, ECM/TGF/fibrosis and inflammation in the same donor-matched experiment.",
    "Context-specific PIEZO2 perturbation or orthogonal calcium/mechanical response assay; null result remains informative.",
    "Sample/donor as statistical unit; detection-rate correction; within-sample permutation; avoid cell-level pseudoreplication."
  ),
  success_rule = c(
    "Reproducible functional mechanical phenotype attenuated by prespecified actomyosin perturbation.",
    "Mechanical phenotype attenuated by prespecified integrin perturbation after state controls.",
    "Core effect remains after competitor-state adjustment or matched-state comparison.",
    "Reproducible context-specific functional evidence, without requiring universal positivity.",
    "Consistent donor-level effect with an orthogonal functional readout."
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(priority_matrix, file.path(result_dir, "Step17F_functional_validation_priority_matrix_v1.csv"))

# ---- Human-readable decision and manuscript update actions ----
decision_lines <- c(
  "# Step 17F second-round evidence resynthesis",
  "",
  "## Decision",
  "",
  "Overall decision: COMPUTATIONAL_TRIANGULATION_COMPLETED_WITH_STRICT_BOUNDARIES",
  "Overall evidence grade: CAUTION",
  "",
  "The second-round computational strengthening is complete. The outputs support a bounded mechanistic hypothesis, but they do not satisfy the predefined requirements for causal or cell-intrinsic claims.",
  "",
  "## S1: co-expression level",
  "",
  paste0(
    "PRJNA607098 descriptively supported ", prjna_pairs_supported, "/", prjna_pairs_tested,
    " tested gene pairs, while GSE130973 supported ", gse130973_pairs_supported, "/", gse130973_pairs_tested,
    ". This is source-dependent co-expression-level evidence, not proof of a cell-intrinsic mechanism."
  ),
  "",
  "## S2: regulatory-axis cross-validation",
  "",
  paste0(
    "GSE338388 provides ", s2_candidate_coverage,
    " frozen-candidate coverage and supports the prespecified TGF-beta exposure × TEAD inhibition regulatory-axis interpretation. It contains no mechanical loading or stiffness perturbation and therefore cannot establish mechanical causality."
  ),
  "",
  "## S3: stiffness form and cross-tissue plausibility",
  "",
  paste0(
    "GSE123100 provides ", s3_123_coverage,
    " candidate coverage and remains a descriptive HTM dose-response analysis. GSE276045 provides ", s3_276_coverage,
    " candidate coverage; actomyosin/Rho is directionally supported in both WT and hTERT stiffness slopes, with camera FDR values ",
    fmt_num(s3_276_acto_wt$camera_fdr), " and ", fmt_num(s3_276_acto_htert$camera_fdr), "."
  ),
  "",
  "The same GSE276045 analysis also shows active competitor programs, including cell-cycle and ECM modules. Therefore the S3 result strengthens pathway plausibility but does not establish mechanistic specificity or proliferation independence.",
  "",
  "## Frozen interpretation",
  "",
  "The highest defensible wording is: evidence-bounded mechanistic hypothesis with computational triangulation.",
  "CAUTION is retained regardless of computational positivity. The computational work does not replace at least three independent donor functional experiments with controlled mechanical perturbation, pathway perturbation, orthogonal readout and state-confounder controls.",
  "",
  "## Manuscript actions",
  "",
  "1. Replace any cell-intrinsic wording with co-expression-level wording for S1.",
  "2. Describe GSE338388 as TGF-beta/TEAD regulatory-axis cross-validation, not mechanical causality.",
  "3. Describe GSE123100 and GSE276045 as cross-tissue stiffness-form evidence, not fascia replication.",
  "4. Present actomyosin/Rho as the leading functional branch while explicitly showing ECM/cell-cycle competition.",
  "5. Keep PIEZO2 as a context-dependent secondary candidate.",
  "6. Retain the CAUTION grade and prioritize donor-repeated functional validation."
)
writeLines(decision_lines, file.path(result_dir, "Step17F_second_round_evidence_resynthesis_decision_v1.md"), useBytes = TRUE)

action_lines <- c(
  "# Step 17F manuscript and study-plan update actions",
  "",
  "| Priority | Action | Location to update |",
  "|---:|---|---|",
  "| 1 | Use co-expression-level language for S1 | Results, Discussion, limitations |",
  "| 2 | State the GSE338388 design as TGF-beta exposure × TEAD inhibition | Methods, Results, figure legends |",
  "| 3 | Label GSE123100/GSE276045 as cross-tissue analyses | Abstract, Results, Discussion |",
  "| 4 | Report actomyosin support together with cell-cycle/ECM competition | Results and Figure 4 legend |",
  "| 5 | Preserve PIEZO2 as context-dependent | Results and candidate prioritization |",
  "| 6 | Keep CAUTION and define the functional upgrade threshold | Abstract, Discussion, conclusion |",
  "| 7 | Design at least 3-donor mechanical perturbation experiments | Future work and experimental protocol |"
)
writeLines(action_lines, file.path(result_dir, "Step17F_manuscript_update_actions_v1.md"), useBytes = TRUE)

message("Step 17F second-round evidence resynthesis completed.")
message("Evidence grade: CAUTION")
message("Highest defensible interpretation: evidence-bounded mechanistic hypothesis with computational triangulation.")
message("Results: ", result_dir)
