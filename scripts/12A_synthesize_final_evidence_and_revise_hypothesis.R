options(stringsAsFactors = FALSE)

# Step 12A: final interim evidence synthesis and hypothesis revision.
# No new expression analysis is performed.

project_dir <- "."
step10d_dir <- file.path(
  project_dir, "results", "08_cross_tissue_validation",
  "10D_evidence_synthesis"
)
step11f_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11F_GSE130973_external_validation_summary"
)
result_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12A_interim_evidence_and_hypothesis_revision"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  file.path(step10d_dir, "Step10D_evidence_summary_v1.csv"),
  file.path(step10d_dir, "Step10D_failed_competition_pairs_v1.csv"),
  file.path(step11f_dir, "GSE130973_candidate_direction_summary_v1.csv"),
  file.path(step11f_dir, "GSE130973_vs_PRJNA607098_candidate_crosswalk_v1.csv"),
  file.path(step11f_dir, "GSE130973_module_boundary_summary_v1.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 12A input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

step10d_summary <- read.csv(required_inputs[[1L]], check.names = FALSE)
failed_competition <- read.csv(required_inputs[[2L]], check.names = FALSE)
gse130973_candidates <- read.csv(required_inputs[[3L]], check.names = FALSE)
cross_source <- read.csv(required_inputs[[4L]], check.names = FALSE)
gse130973_modules <- read.csv(required_inputs[[5L]], check.names = FALSE)

get_value <- function(table, item) {
  value <- table$value[table$item == item]
  if (length(value) == 0L) return(NA_character_)
  as.character(value[[1L]])
}

candidate_support_130973 <- sum(
  gse130973_candidates$descriptive_directional_support_4_of_5
)
integrin_support_130973 <- sum(
  gse130973_candidates$descriptive_directional_support_4_of_5[
    gse130973_candidates$module == "integrin_focal_adhesion"
  ]
)
actomyosin_support_130973 <- sum(
  gse130973_candidates$descriptive_directional_support_4_of_5[
    gse130973_candidates$module == "actomyosin_rho"
  ]
)
piezo2_row <- gse130973_candidates[
  gse130973_candidates$gene == "PIEZO2", , drop = FALSE
]
piezo2_130973 <- nrow(piezo2_row) == 1L &&
  isTRUE(piezo2_row$descriptive_directional_support_4_of_5[[1L]])
competition_robustness_pass <- isTRUE(
  get_value(step10d_summary, "competition_robustness_gate") == "TRUE"
)

evidence_table <- data.frame(
  claim_domain = c(
    "module_level_core_program",
    "integrin_candidate_specificity",
    "actomyosin_candidate_conservation",
    "channel_specific_PIEZO2",
    "competition_robustness",
    "independent_replication",
    "F7_F8_direct_equivalence",
    "overall_evidence_grade"
  ),
  current_status = c(
    "SUPPORTED_WITH_BOUNDARY",
    paste0("PARTIAL_", integrin_support_130973, "_OF_6_IN_GSE130973"),
    paste0("DISCORDANT_", actomyosin_support_130973, "_OF_4_IN_GSE130973"),
    if (piezo2_130973) "SOURCE_DEPENDENT_POSITIVE" else "NOT_SUPPORTED",
    if (competition_robustness_pass) "PASS" else "FAIL",
    "PARTIAL_ONLY",
    "UNRESOLVED",
    "CAUTION"
  ),
  evidence_basis = c(
    "Core module directions were positive in PRJNA607098 and GSE130973 exploratory sample-level audits.",
    "Integrin/focal-adhesion direction was consistent for 5/6 candidates in both the primary and exploratory sources, but state competitors remain active.",
    "CFL1/CNN2/CDC42 were not directionally conserved in GSE130973, despite module-level positivity.",
    "PIEZO2 was negative in PRJNA607098 but directionally positive in GSE130973; it is not a universal core marker.",
    "The PRJNA607098 integrin branch failed residual-support checks against TGF/fibrosis, inflammation, hypoxia, and cell-cycle competitors.",
    "The integrated atlas overlaps GSE173252 and PRJNA607098; GSE130973 is non-overlapping but marker-defined and aging-skin specific.",
    "The paper F8 label and atlas F7 label have not been proven one-to-one.",
    "The evidence is observational, source-dependent, and not causal."
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(
  evidence_table,
  file.path(result_dir, "12A_final_evidence_domain_summary_v1.csv")
)

decision_lines <- c(
  "# Step 12A interim evidence synthesis and revised research hypothesis",
  "",
  "## Overall decision",
  "",
  "The project remains scientifically valuable but must move from a single-gene mechanosensor hypothesis to a bounded, multicomponent fibroblast-state hypothesis. The strongest reproducible signal is a fibroblast-associated ECM/integrin–cytoskeletal mechanotransduction program. The current data do not establish a universal PIEZO2 mechanism, a conserved actomyosin candidate-gene signature across tissues, direct F8 replication, or causal mechanosensitivity.",
  "",
  "**Current evidence grade: CAUTION.**",
  "",
  "## Revised primary hypothesis",
  "",
  "Fascia-like or myofibroblast-enriched fibroblast states exhibit a multicomponent extracellular-matrix/integrin–cytoskeletal mechanotransduction program, whose module-level expression is reproducible across selected human single-cell sources but whose individual genes and channel-specific components are context-dependent.",
  "",
  "## Revised secondary hypotheses",
  "",
  "1. Integrin/focal-adhesion enrichment is a reproducible component, but its specificity must be separated from ECM remodeling and TGF/fibrosis state.",
  "2. Actomyosin/Rho activity is supported at the module level, but the conserved candidate-gene subset is incomplete and requires functional confirmation.",
  "3. PIEZO2 is a context-dependent candidate rather than a universal marker of the state.",
  "4. Direct mechanosensitivity requires functional evidence such as stiffness/force response, traction, YAP/TAZ localization, or perturbation of integrin/actomyosin components.",
  "",
  "## What is now established",
  "",
  "- GSE173252 discovery analysis produced the original program and candidate panel.",
  "- PRJNA607098 provided 12-sample validation, with core module positivity but failed competitor robustness and negative PIEZO2.",
  "- GSE130973 supplied five non-overlapping subjects with a marker-defined fibroblast-state audit and partial Integrin support.",
  "- GSE175817 is not currently usable for donor-level validation because its processed matrix lacks donor-column mapping.",
  "",
  "## What remains unresolved",
  "",
  "- Whether the integrin signal remains after comparison within a strictly annotated fibroblast population.",
  "- Whether actomyosin candidate genes are truly conserved across independent tissues.",
  "- Whether PIEZO2 reflects a real context-specific branch or a dataset/cluster-composition effect.",
  "- Whether the atlas F7 state is anatomically equivalent to the paper F8 state.",
  "",
  "## Recommended next research phase",
  "",
  "### Phase 1: one final data audit",
  "",
  "Use only a genuinely independent source with explicit fibroblast labels and donor-level metadata. Do not pool the current sources or relax candidate thresholds. If no such source is recovered, stop expanding the public-data search.",
  "",
  "### Phase 2: functional validation",
  "",
  "Prioritize experiments that distinguish mechanism from state: matrix stiffness or controlled mechanical loading, traction/contractility readouts, YAP/TAZ nuclear localization, and perturbation of ITGB1/ITGAV–focal adhesion or Rho/actomyosin components. PIEZO2 should be tested as a secondary, context-dependent branch rather than the central mechanism.",
  "",
  "### Phase 3: manuscript framing",
  "",
  "Present the work as an evidence-bounded human single-cell reanalysis plus external triangulation. The central claim should be module-level fibroblast mechanotransduction association, with explicit boundaries around causality, F7/F8 equivalence, source overlap, and PIEZO2 inconsistency.",
  "",
  "## Material Passport",
  "",
  "- Inputs: Step 10D provenance-aware synthesis and Step 11F non-overlapping GSE130973 summary.",
  "- Transformation: claim-domain synthesis and hypothesis revision; no new expression analysis.",
  "- Integrity status: interim research-direction decision; final evidence grade remains CAUTION.",
  "- Reproducibility: the machine-readable evidence table is stored beside this report."
)
decision_path <- file.path(
  result_dir, "12A_interim_evidence_and_hypothesis_revision_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 12A interim evidence synthesis completed.")
message("Current evidence grade: CAUTION")
message("Revised primary hypothesis written to: ", decision_path)
message("Evidence domain summary: ", file.path(result_dir, "12A_final_evidence_domain_summary_v1.csv"))
