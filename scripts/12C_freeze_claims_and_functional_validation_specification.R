options(stringsAsFactors = FALSE)

# Step 12C: freeze computational claims and define a functional-validation-ready
# specification. No new expression analysis or data download is performed.

project_dir <- "."
input_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis"
)
step12a_dir <- file.path(input_dir, "12A_interim_evidence_and_hypothesis_revision")
step12b_dir <- file.path(input_dir, "12B_direction_and_validation_plan")
result_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12C_computational_closeout_and_functional_validation"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

evidence_path <- file.path(step12a_dir, "12A_final_evidence_domain_summary_v1.csv")
direction_path <- file.path(step12b_dir, "12B_research_direction_and_validation_plan_v1.md")
priority_path <- file.path(step12b_dir, "12B_validation_priority_matrix_v1.csv")
required_inputs <- c(evidence_path, direction_path, priority_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 12A/12B input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

evidence <- read.csv(evidence_path, check.names = FALSE)
status_of <- function(domain) {
  value <- evidence$current_status[evidence$claim_domain == domain]
  if (length(value) == 0L) return(NA_character_)
  as.character(value[[1L]])
}

claim_ledger <- data.frame(
  claim_id = paste0("C", sprintf("%02d", 1:7)),
  claim_domain = c(
    "module_level_association",
    "integrin_branch",
    "actomyosin_branch",
    "PIEZO2_branch",
    "competition_specificity",
    "F7_F8_replication",
    "causal_mechanosensitivity"
  ),
  evidence_status = c(
    status_of("module_level_core_program"),
    status_of("integrin_candidate_specificity"),
    status_of("actomyosin_candidate_conservation"),
    status_of("channel_specific_PIEZO2"),
    status_of("competition_robustness"),
    status_of("F7_F8_direct_equivalence"),
    "NOT_TESTED"
  ),
  permitted_interpretation = c(
    "Selected human single-cell sources show a fibroblast-associated multicomponent ECM/integrin-cytoskeletal program.",
    "The integrin/focal-adhesion branch is a prioritized candidate, but its specificity from generic ECM/TGF state is unresolved.",
    "Actomyosin/Rho activity is a plausible pathway-level component; the frozen candidate-gene subset is not conserved.",
    "PIEZO2 is a context-dependent candidate and should be reported as secondary evidence.",
    "The current program is not specific enough to exclude fibrosis/inflammation-related competing states.",
    "The atlas F7 state cannot be claimed as a direct replication of the paper F8 state.",
    "Causal mechanosensitivity requires controlled mechanical perturbation and functional readouts."
  ),
  prohibited_overclaim = c(
    "Do not call this a universal mechanosensor signature or causal mechanism.",
    "Do not describe integrin enrichment as mechanosensitivity-specific without functional perturbation.",
    "Do not force all actomyosin genes to support the hypothesis or discard the module because of candidate discordance.",
    "Do not use PIEZO2 positivity as a required main-study success criterion.",
    "Do not claim that the program is independent of TGF/fibrosis or inflammation state.",
    "Do not label the external atlas as direct F7/F8 replication.",
    "Do not infer causality from cross-sectional single-cell expression."
  ),
  upgrade_requirement = c(
    "Donor-replicated mechanical response plus at least one orthogonal functional readout.",
    "Prespecified integrin perturbation attenuates the mechanical phenotype after state controls.",
    "Prespecified Rho/actomyosin perturbation changes contractility or the mechanical response.",
    "PIEZO2 shows reproducible context-specific functional evidence; null results remain informative.",
    "Matched ECM/TGF/inflammation controls show the effect is not fully explained by generic state activation.",
    "A source with explicit matched F7/F8-equivalent labels and donor-level metadata becomes available.",
    "Causal evidence requires perturbation, rescue or orthogonal convergence, and donor replication."
  ),
  stringsAsFactors = FALSE
)

validation_endpoints <- data.frame(
  endpoint_id = paste0("E", sprintf("%02d", 1:7)),
  priority = c("PRIMARY", "PRIMARY", "PRIMARY", "SECONDARY", "SECONDARY", "CONTROL", "STATISTICAL"),
  endpoint = c(
    "Mechanical-response phenotype",
    "Integrin/focal-adhesion dependence",
    "Rho/actomyosin dependence",
    "YAP/TAZ nuclear localization",
    "PIEZO2 context-specific response",
    "ECM/TGF/inflammation and viability controls",
    "Donor-level effect size and uncertainty"
  ),
  recommended_readout = c(
    "Stiffness or controlled loading; traction, contraction, morphology, or force-response assay.",
    "ITGB1/ITGAV or focal-adhesion perturbation; focal-adhesion imaging and functional response.",
    "RHOA/ROCK or actomyosin perturbation; traction/contractility and pathway readouts.",
    "Imaging-based nuclear:cytoplasmic YAP/TAZ ratio under matched mechanical conditions.",
    "PIEZO2 expression/protein/function only after the core phenotype is established.",
    "State markers, viability, proliferation, and matched non-targeting/vehicle controls.",
    "Donor as biological unit; mixed-effects or donor-level analysis with effect size and confidence interval."
  ),
  minimum_replication = c(
    "At least 3 independent donors; technical repeats nested within donor.",
    "At least 3 independent donors and prespecified perturbation control.",
    "At least 3 independent donors and orthogonal contractility readout.",
    "Same donors and matched conditions as the primary mechanical experiment.",
    "Same donors where feasible; report positive and null outcomes.",
    "All donors and conditions; exclude or model failed viability conditions.",
    "No cell-level pseudoreplication as the primary inferential analysis."
  ),
  success_condition = c(
    "Effect is donor-replicated and not explained by viability/proliferation alone.",
    "Perturbation attenuates the primary mechanical phenotype.",
    "Perturbation changes contractility or the primary mechanical phenotype.",
    "Localization follows the mechanical condition and converges with functional data.",
    "Context-specific evidence is reproducible; absence does not invalidate the core program.",
    "Competing-state and technical explanations are quantified rather than assumed absent.",
    "Report effect size, donor variability, uncertainty, and prespecified sensitivity analyses."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  claim_ledger,
  file.path(result_dir, "12C_frozen_claims_and_boundaries_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  validation_endpoints,
  file.path(result_dir, "12C_functional_validation_endpoint_matrix_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

decision_lines <- c(
  "# Step 12C computational closeout and functional-validation specification",
  "",
  "## Closeout decision",
  "",
  "The public single-cell reanalysis phase is computationally mature enough to freeze its claims. Further analysis should be limited to one prequalified independent source, if available; no pooling, threshold relaxation, or post hoc candidate replacement is authorized by this plan.",
  "",
  "**Overall evidence grade: CAUTION.**",
  "",
  "## Final primary claim",
  "",
  "Selected human fibroblast-enriched states show an associated multicomponent ECM/integrin–cytoskeletal mechanotransduction program. The evidence is observational and context-dependent; it does not yet establish causal mechanosensitivity, universal PIEZO2 involvement, or direct F7/F8 equivalence.",
  "",
  "## Primary validation package",
  "",
  "1. Apply a controlled mechanical challenge using matched conditions and at least three independent donors.",
  "2. Measure a functional mechanical phenotype such as traction, contractility, force response, or mechanically induced morphology.",
  "3. Test prespecified ITGB1/ITGAV–focal-adhesion and Rho/actomyosin perturbations.",
  "4. Add YAP/TAZ localization as an orthogonal pathway readout.",
  "5. Quantify ECM/TGF/inflammation, viability, and proliferation controls in the same experiment.",
  "6. Treat PIEZO2 as a secondary context-specific branch and report null results without changing the primary hypothesis.",
  "",
  "## Analysis and reporting lock",
  "",
  "- Donor is the biological replicate; cells are nested observations.",
  "- Primary outcomes and exclusion rules must be fixed before examining condition-specific results.",
  "- Report effect sizes and uncertainty, not only p-values or cell-level percentages.",
  "- Separate pathway-level evidence from candidate-gene-level evidence.",
  "- Preserve all discordant results in the final report.",
  "",
  "## Material Passport",
  "",
  "- Inputs: Step 12A evidence summary and Step 12B direction lock.",
  "- Transformation: computational claim freeze and functional endpoint specification.",
  "- New data downloaded: none.",
  "- Outputs: frozen claim ledger and endpoint matrix.",
  "- Status: computational closeout gate; ready for experimental planning or manuscript drafting."
)

decision_path <- file.path(
  result_dir, "12C_computational_closeout_and_functional_validation_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 12C computational closeout completed.")
message("Evidence grade remains: CAUTION")
message("Frozen claims: ", file.path(result_dir, "12C_frozen_claims_and_boundaries_v1.csv"))
message("Validation endpoints: ", file.path(result_dir, "12C_functional_validation_endpoint_matrix_v1.csv"))
message("Decision: ", decision_path)
