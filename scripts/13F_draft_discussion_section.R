options(stringsAsFactors = FALSE)

# Step 13F: evidence-grounded Discussion section draft.
# No new analysis, hypothesis testing, or data download.

project_dir <- "."
step10d_dir <- file.path(project_dir, "results", "08_cross_tissue_validation", "10D_evidence_synthesis")
step11f_dir <- file.path(project_dir, "results", "09_independent_external_source_screening", "11F_GSE130973_external_validation_summary")
step12b_dir <- file.path(project_dir, "results", "10_final_evidence_synthesis", "12B_direction_and_validation_plan")
step12c_dir <- file.path(project_dir, "results", "10_final_evidence_synthesis", "12C_computational_closeout_and_functional_validation")
step13e_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13E_results_draft")
result_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13F_discussion_draft")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(step10d_dir, "Step10D_evidence_summary_v1.csv")
external_decision_path <- file.path(step11f_dir, "GSE130973_external_validation_boundary_decision_v1.md")
direction_path <- file.path(step12b_dir, "12B_research_direction_and_validation_plan_v1.md")
endpoint_path <- file.path(step12c_dir, "12C_functional_validation_endpoint_matrix_v1.csv")
results_path <- file.path(step13e_dir, "13E_results_draft_v1.md")
required_inputs <- c(summary_path, external_decision_path, direction_path, endpoint_path, results_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 13F input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

step10d <- read.csv(summary_path, check.names = FALSE)
endpoints <- read.csv(endpoint_path, check.names = FALSE)

summary_value <- function(item) {
  value <- step10d$value[step10d$item == item]
  if (length(value) == 0L) return(NA_character_)
  as.character(value[[1L]])
}

primary_grade <- summary_value("final_evidence_grade")
interpretation_class <- summary_value("interpretation_class")

discussion_lines <- c(
  "# Discussion",
  "",
  "## Principal findings",
  "",
  "This study provides evidence for a fibroblast-associated, multicomponent ECM/integrin–cytoskeletal mechanotransduction program in selected human single-cell sources. The most stable signal was observed at the module level rather than as a universally conserved set of individual genes. Paired sample-level validation supported the core program and leave-one-sample-out robustness, whereas the competition-robustness gate failed. Independent GSE130973 triangulation provided partial support for the integrin branch but revealed actomyosin candidate-gene discordance and a source-dependent PIEZO2 direction.",
  "",
  "The appropriate interpretation is therefore an evidence-bounded association hypothesis rather than a demonstrated causal mechanosensitivity mechanism. The final project evidence grade remains CAUTION, and the integrated interpretation class was ",
  interpretation_class,
  ".",
  "",
  "## A multicomponent program is more defensible than a single-gene model",
  "",
  "The results argue against treating one channel or one transcript as a universal definition of fibroblast mechanosensitivity. A coordinated program involving ECM organization, integrin/focal-adhesion signaling, cytoskeletal regulation, and YAP/TAZ-related components was more reproducible than any single candidate branch. This interpretation is consistent with the biological expectation that mechanical information is distributed across matrix attachment, force transmission, cytoskeletal state, and downstream transcriptional responses [CITATION REQUIRED: mechanobiology and integrin–cytoskeleton literature].",
  "",
  "The module-level result should nevertheless be described as an expression-level association. The current data cannot determine whether the program is induced by mechanical input, marks a pre-existing fibroblast state that experiences different mechanical environments, or reflects a mixture of both processes.",
  "",
  "## Integrin support is promising but not yet specific",
  "",
  "The integrin/focal-adhesion branch was the strongest candidate-level component across the available sources, with 5/6 directional support in the primary paired validation and 5/6 descriptive directional support in GSE130973. This concordance makes integrin-linked force transmission a rational priority for functional testing. It does not, however, establish mechanosensitivity-specificity. The primary competition audit failed against TGF/fibrosis, inflammation, hypoxia, and cell-cycle-related programs, and the external source also showed positive ECM, TGF/fibrosis, and inflammation modules.",
  "",
  "Thus, the integrin result should be interpreted as a prioritized mechanistic candidate that remains confounded by fibroblast activation and remodeling state. A meaningful upgrade would require matched state controls and attenuation of the mechanical phenotype after prespecified ITGB1/ITGAV or focal-adhesion perturbation.",
  "",
  "## Actomyosin and PIEZO2 require context-sensitive interpretation",
  "",
  "The actomyosin/Rho module was positive in the primary paired validation, but the candidate-gene pattern was not conserved in the independent GSE130973 audit, where only 1/4 frozen actomyosin candidates met the descriptive directional rule. This result does not eliminate actomyosin biology; instead, it suggests that pathway-level contractility may be more stable than the expression direction of a small frozen candidate subset. Functional traction, contraction, or Rho/ROCK perturbation experiments are therefore more informative than forcing complete transcript-level agreement [CITATION REQUIRED: actomyosin and fibroblast contractility literature].",
  "",
  "PIEZO2 showed the clearest source dependence: it was not supported in the PRJNA607098 paired analysis but was directionally positive in GSE130973. PIEZO2 should consequently remain a secondary context-specific branch. A null PIEZO2 result would not invalidate the broader ECM/integrin–cytoskeletal hypothesis, and PIEZO2 positivity alone should not be used as evidence of a universal mechanosensor mechanism [CITATION REQUIRED: PIEZO2 in fibroblast or connective-tissue mechanobiology].",
  "",
  "## Alternative explanations and specificity limits",
  "",
  "The failed competition-robustness gate is scientifically informative. It indicates that the candidate integrin program is correlated with state programs that are also associated with matrix remodeling, TGF/fibrosis, inflammation, hypoxia, or proliferation. This leaves at least three non-exclusive explanations: a genuinely mechanosensitive program embedded within a broader activated fibroblast state; a remodeling-state program that is mechanically correlated but not mechanistically driven; or a mixture of cell-state and sample-composition effects.",
  "",
  "The current design cannot distinguish these explanations by transcriptomic association alone. The next experiment must therefore manipulate mechanical context while measuring viability, proliferation, ECM/TGF/inflammation state, and an orthogonal functional outcome in the same donor-matched system.",
  "",
  "## Implications for future work",
  "",
  "The locked next phase should begin with a controlled stiffness or mechanical-loading challenge in at least three independent donors. The primary outcome should be a donor-replicated functional mechanical phenotype, such as traction, contractility, force response, or mechanically induced morphology. ITGB1/ITGAV–focal-adhesion and Rho/actomyosin perturbations should then be used to test pathway dependence, with YAP/TAZ nuclear localization as an orthogonal readout. PIEZO2 should be assessed as a secondary context branch rather than a required success criterion.",
  "",
  "These experiments are designed to move the project beyond CAUTION only if the mechanical phenotype is donor-replicated, not explained by viability/proliferation or generic fibrosis state, and attenuated by a prespecified mechanistic perturbation. Transcriptomic agreement alone is insufficient for that upgrade.",
  "",
  "## Limitations",
  "",
  "Several limitations constrain the interpretation. First, the analysis relied partly on processed public objects and remote gene-panel extraction rather than a harmonized raw-data reprocessing pipeline. Second, the atlas F7 state and paper F8 state were not proven equivalent. Third, GSE130973 was a marker-defined, aging-skin fibroblast-state audit and therefore provides external triangulation rather than direct replication. Fourth, GSE175817 could not be used for donor-level validation because donor-column mapping remained unresolved. Finally, all current evidence is observational and cannot establish causality.",
  "",
  "## Conclusion",
  "",
  "The study supports a bounded hypothesis that fibroblast-enriched human states are associated with a multicomponent ECM/integrin–cytoskeletal mechanotransduction program. Integrin-linked force transmission is the leading mechanistic candidate, actomyosin activity remains a pathway-level hypothesis with candidate-gene discordance, and PIEZO2 is context-dependent. The decisive next step is donor-replicated functional mechanical perturbation with state controls and orthogonal mechanotransduction readouts."
)

writeLines(
  discussion_lines,
  file.path(result_dir, "13F_discussion_draft_v1.md"),
  useBytes = TRUE
)

trace_table <- data.frame(
  discussion_section = c(
    "Principal findings",
    "Multicomponent program",
    "Integrin specificity",
    "Actomyosin and PIEZO2",
    "Alternative explanations",
    "Future work",
    "Limitations",
    "Conclusion"
  ),
  evidence_basis = c(
    "Step 13E Results; Step 10D synthesis",
    "Step 10D module-level support; Step 12A revised hypothesis",
    "Step 10D competition audit; Step 11F external validation",
    "Step 10D and Step 11F candidate/module summaries",
    "Step10D_failed_competition_pairs_v1.csv",
    "Step 12B direction lock; Step 12C endpoint matrix",
    "Step 10A, 11C2, 11F, and 12C boundaries",
    "Step 12A–12C final synthesis"
  ),
  allowed_claim_level = c(
    "Interpretation with explicit CAUTION boundary",
    "Module-level association",
    "Prioritized candidate, not specific mechanism",
    "Context-dependent pathway/candidate interpretation",
    "Alternative explanations remain unresolved",
    "Testable functional validation plan",
    "Observational and provenance limitations",
    "Bounded hypothesis and next decisive experiment"
  ),
  citation_status = c(
    "Project data citation required",
    "External mechanobiology citations required",
    "Integrin/focal-adhesion citations required",
    "Actomyosin and PIEZO2 citations required",
    "Project data citation required",
    "No new citation required for proposed design; methods citations may be added",
    "Project data citation required",
    "Project data citation required"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  trace_table,
  file.path(result_dir, "13F_discussion_evidence_trace_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

citation_table <- data.frame(
  topic = c(
    "ECM-integrin-cytoskeletal mechanotransduction",
    "Focal adhesion force transmission in fibroblasts",
    "Actomyosin/Rho regulation of fibroblast contractility",
    "YAP/TAZ as a mechanical-state readout",
    "PIEZO2 in connective tissue or fibroblast mechanobiology"
  ),
  reason_needed = c(
    "Contextualize the multicomponent mechanotransduction interpretation.",
    "Support the rationale for integrin/focal-adhesion perturbation.",
    "Support pathway-level interpretation despite candidate-gene discordance.",
    "Support the proposed orthogonal functional readout.",
    "Contextualize the source-dependent secondary branch."
  ),
  status = "CITATION_REQUIRED",
  stringsAsFactors = FALSE
)
write.csv(
  citation_table,
  file.path(result_dir, "13F_citations_needed_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Step 13F Discussion draft completed.")
message("Discussion draft: ", file.path(result_dir, "13F_discussion_draft_v1.md"))
message("Evidence trace: ", file.path(result_dir, "13F_discussion_evidence_trace_v1.csv"))
message("Citations needed: ", file.path(result_dir, "13F_citations_needed_v1.csv"))
