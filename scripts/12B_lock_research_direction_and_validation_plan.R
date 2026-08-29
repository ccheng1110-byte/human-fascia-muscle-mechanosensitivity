options(stringsAsFactors = FALSE)

# Step 12B: lock the post-audit research direction and validation priorities.
# This is a planning/audit step. No new expression analysis or data download.

project_dir <- "."
input_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12A_interim_evidence_and_hypothesis_revision"
)
result_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12B_direction_and_validation_plan"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

evidence_path <- file.path(input_dir, "12A_final_evidence_domain_summary_v1.csv")
revision_path <- file.path(
  input_dir, "12A_interim_evidence_and_hypothesis_revision_v1.md"
)
if (!file.exists(evidence_path) || !file.exists(revision_path)) {
  stop("Step 12A outputs are required before Step 12B.")
}

evidence <- read.csv(evidence_path, check.names = FALSE)
status_of <- function(domain) {
  value <- evidence$current_status[evidence$claim_domain == domain]
  if (length(value) == 0L) return(NA_character_)
  as.character(value[[1L]])
}

priority_matrix <- data.frame(
  priority = 1:6,
  workstream = c(
    "Controlled mechanical challenge",
    "Integrin/focal-adhesion perturbation",
    "Rho/actomyosin perturbation",
    "YAP/TAZ mechanotransduction readout",
    "PIEZO2 context-dependence",
    "One final public-data audit gate"
  ),
  central_question = c(
    "Does a controlled stiffness or mechanical-loading change alter the proposed fibroblast program and cell behavior?",
    "Is the integrin/focal-adhesion branch mechanistically required, rather than only correlated with ECM/TGF state?",
    "Does Rho/ROCK/actomyosin manipulation change contractility or the module response despite candidate-gene discordance?",
    "Does mechanical context change nuclear YAP/TAZ localization in the predicted direction?",
    "Is PIEZO2 a reproducible mechanosensitive branch in a defined cellular context, or only a dataset-dependent marker?",
    "Can an independent source with explicit fibroblast labels and donor metadata resolve the remaining annotation uncertainty?"
  ),
  rationale = c(
    "The module-level signal is reproducible but remains observational and source-dependent.",
    "Integrin support is 5/6, but competition robustness failed against ECM/TGF/inflammation-related programs.",
    "The actomyosin module is plausible, but only 1/4 frozen candidates were directionally conserved in GSE130973.",
    "YAP/TAZ provides a pathway-level readout that is closer to mechanotransduction than transcript abundance alone.",
    "PIEZO2 was negative in PRJNA607098 and positive in GSE130973; it should not remain the primary hypothesis.",
    "Current F7/F8 equivalence and strict fibroblast-level replication remain unresolved."
  ),
  minimum_design = c(
    "At least 3 independent donors; matched soft/stiff or unloaded/loaded conditions; viability and proliferation controls.",
    "ITGB1 or ITGAV/focal-adhesion perturbation with matched non-targeting control and rescue or orthogonal readout when feasible.",
    "ROCK/RHOA or actomyosin perturbation with contractility/traction readout and cell-state controls.",
    "Immunofluorescence or imaging-based nuclear:cytoplasmic YAP/TAZ quantification, donor treated as experimental unit.",
    "Run only as a secondary branch after the core mechanical phenotype is established; include PIEZO2 protein/function if technically feasible.",
    "Predefine inclusion criteria; do not pool datasets, relax frozen thresholds, or treat marker-defined clusters as direct F7/F8 replication."
  ),
  primary_outcome = c(
    "Mechanical response effect size with donor-level confidence interval.",
    "Loss or attenuation of the mechanical response after integrin/focal-adhesion perturbation.",
    "Loss or attenuation of contractility and/or the program response after Rho/actomyosin perturbation.",
    "Predicted change in YAP/TAZ nuclear localization under mechanical challenge.",
    "Context-specific association, not universal positivity; report null results explicitly.",
    "Qualified independent replication, or a documented stop decision for further public-data expansion."
  ),
  go_no_go_rule = c(
    "GO if response is donor-replicated and not explained by viability/proliferation alone.",
    "GO if perturbation changes the mechanical phenotype; otherwise downgrade to association only.",
    "GO if functional contractility changes even when individual marker genes are discordant.",
    "GO if localization follows mechanical conditions and agrees with functional readouts.",
    "Do not make PIEZO2 a primary success criterion; treat it as supportive or context-specific evidence.",
    "GO only if all metadata and independence criteria pass; otherwise stop searching and move to functional validation."
  ),
  evidence_role = c(
    "Upgrades association toward functional evidence.",
    "Tests specificity and causal relevance of the strongest transcriptomic branch.",
    "Tests pathway-level conservation without over-relying on the discordant frozen genes.",
    "Links mechanical context to a canonical downstream mechanotransduction readout.",
    "Refines, but cannot by itself rescue, the central hypothesis.",
    "Resolves annotation/provenance uncertainty only; cannot establish causality."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  priority_matrix,
  file.path(result_dir, "12B_validation_priority_matrix_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

decision_lines <- c(
  "# Step 12B research-direction lock and validation plan",
  "",
  "## Decision",
  "",
  "The project should proceed, but the primary endpoint must change from universal single-gene mechanosensitivity to a donor-replicated, multicomponent fibroblast mechanotransduction phenotype.",
  "",
  "**Current evidence grade remains CAUTION.** No public-data result should be presented as causal or as direct F7/F8 replication.",
  "",
  "## Evidence state carried forward",
  "",
  paste0("- Core module: ", status_of("module_level_core_program"), "."),
  paste0("- Integrin specificity: ", status_of("integrin_candidate_specificity"), "."),
  paste0("- Actomyosin candidate conservation: ", status_of("actomyosin_candidate_conservation"), "."),
  paste0("- PIEZO2: ", status_of("channel_specific_PIEZO2"), "."),
  paste0("- Competition robustness: ", status_of("competition_robustness"), "."),
  paste0("- F7/F8 equivalence: ", status_of("F7_F8_direct_equivalence"), "."),
  "",
  "## Locked research hierarchy",
  "",
  "1. Primary: controlled mechanical challenge plus donor-level functional response.",
  "2. Mechanistic discriminator: ITGB1/ITGAV–focal adhesion and Rho/actomyosin perturbation.",
  "3. Orthogonal pathway readout: YAP/TAZ nuclear localization, traction, or contractility.",
  "4. Secondary context branch: PIEZO2, reported regardless of whether the result is positive or null.",
  "5. Optional final public-data audit: only a source meeting all predefined metadata and independence criteria.",
  "",
  "## Recommended relative timeline",
  "",
  "- Weeks 1–2: lock protocols, donor-level unit of analysis, controls, and preregistered primary outcomes.",
  "- Weeks 3–8: mechanical challenge experiment with viability/proliferation and state-marker controls.",
  "- Weeks 9–12: integrin and Rho/actomyosin perturbation experiments with orthogonal functional readouts.",
  "- Weeks 13–15: YAP/TAZ imaging, PIEZO2 secondary analysis, robustness and donor-level statistics.",
  "- Weeks 16–18: final evidence integration, figures, manuscript, and limitation audit.",
  "",
  "## Risk controls",
  "",
  "- Donor heterogeneity: use at least three independent donors and treat donor, not cell, as the biological replicate.",
  "- ECM/TGF confounding: include matched state controls and report ECM/TGF activity alongside mechanotransduction readouts.",
  "- Cell-cycle or viability artifacts: measure viability and proliferation and perform sensitivity analyses excluding affected conditions.",
  "- Candidate-gene discordance: interpret the actomyosin branch at pathway/function level, not by forcing all frozen genes to agree.",
  "- PIEZO2 instability: do not use PIEZO2 as a required success criterion for the main study conclusion.",
  "- Public-data overexpansion: stop after one additional qualified audit fails or cannot resolve the key uncertainty.",
  "",
  "## Provisional upgrade criteria",
  "",
  "The evidence can be upgraded beyond CAUTION only if the mechanical phenotype is reproduced across donors, is not explained by viability/proliferation or generic fibrosis state, and is attenuated by a prespecified integrin or actomyosin perturbation. Transcriptomic agreement alone is insufficient.",
  "",
  "## Material Passport",
  "",
  "- Input: Step 12A machine-readable evidence summary.",
  "- Transformation: evidence-to-action prioritization; no new expression analysis or data download.",
  "- Output: validation priority matrix and research-direction decision.",
  "- Status: planning gate for the next experimental or qualified-data phase."
)

decision_path <- file.path(
  result_dir, "12B_research_direction_and_validation_plan_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 12B research-direction lock completed.")
message("Current evidence grade remains: CAUTION")
message("Validation priority matrix: ", file.path(result_dir, "12B_validation_priority_matrix_v1.csv"))
message("Research-direction decision: ", decision_path)
