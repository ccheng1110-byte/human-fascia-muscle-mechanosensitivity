options(stringsAsFactors = FALSE)

# Step 13A: build a manuscript-ready evidence pack from frozen outputs.
# No new expression analysis, hypothesis testing, or data download.

project_dir <- "."
result_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13A_manuscript_evidence_pack"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

step12a_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12A_interim_evidence_and_hypothesis_revision"
)
step12b_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12B_direction_and_validation_plan"
)
step12c_dir <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12C_computational_closeout_and_functional_validation"
)

required_inputs <- c(
  file.path(step12a_dir, "12A_final_evidence_domain_summary_v1.csv"),
  file.path(step12a_dir, "12A_interim_evidence_and_hypothesis_revision_v1.md"),
  file.path(step12b_dir, "12B_research_direction_and_validation_plan_v1.md"),
  file.path(step12c_dir, "12C_frozen_claims_and_boundaries_v1.csv"),
  file.path(step12c_dir, "12C_functional_validation_endpoint_matrix_v1.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Required Step 12 outputs are missing: ", paste(missing_inputs, collapse = "; "))
}

evidence <- read.csv(required_inputs[[1L]], check.names = FALSE)
claims <- read.csv(required_inputs[[4L]], check.names = FALSE)
endpoints <- read.csv(required_inputs[[5L]], check.names = FALSE)

evidence_file_inventory <- data.frame(
  evidence_id = c(
    "EVID_10D_SUMMARY",
    "EVID_10D_COMPETITION",
    "EVID_08C2_PAIRED",
    "EVID_11F_EXTERNAL",
    "EVID_12A_SYNTHESIS",
    "EVID_12B_DIRECTION",
    "EVID_12C_CLAIMS",
    "EVID_12C_ENDPOINTS"
  ),
  evidence_role = c(
    "Cross-tissue evidence synthesis",
    "Specificity and competition failure audit",
    "Paired sample-level external validation",
    "Independent GSE130973 external validation",
    "Revised hypothesis and evidence grade",
    "Locked research direction",
    "Frozen manuscript claim boundaries",
    "Functional validation specification"
  ),
  file_path = c(
    file.path(project_dir, "results", "08_cross_tissue_validation", "10D_evidence_synthesis", "Step10D_evidence_summary_v1.csv"),
    file.path(project_dir, "results", "08_cross_tissue_validation", "10D_evidence_synthesis", "Step10D_failed_competition_pairs_v1.csv"),
    file.path(project_dir, "results", "06_external_validation", "skin_fibroblast_atlas_2025", "08C2_sample_level_F7_PRJNA607098", "PRJNA607098_F7_sample_paired_gene_contrasts_v1.csv"),
    file.path(project_dir, "results", "09_independent_external_source_screening", "11F_GSE130973_external_validation_summary", "GSE130973_external_validation_boundary_decision_v1.md"),
    required_inputs[[2L]],
    required_inputs[[3L]],
    required_inputs[[4L]],
    required_inputs[[5L]]
  ),
  stringsAsFactors = FALSE
)
evidence_file_inventory$exists <- file.exists(evidence_file_inventory$file_path)

claim_to_evidence <- data.frame(
  claim_id = c("C01", "C02", "C03", "C04", "C05", "C06", "C07"),
  manuscript_claim = c(
    "A multicomponent ECM/integrin-cytoskeletal program is associated with selected human fibroblast-enriched states.",
    "The integrin/focal-adhesion branch is a reproducible candidate but is not yet specific from generic ECM/TGF state.",
    "Actomyosin/Rho activity is a pathway-level candidate, while the frozen candidate-gene subset is discordant.",
    "PIEZO2 is context-dependent and should not be presented as a universal marker.",
    "Competition robustness is incomplete, so the program cannot be claimed to be fibrosis/inflammation-independent.",
    "The available atlas evidence is not a direct F7/F8 replication.",
    "Causal mechanosensitivity remains untested by the current observational data."
  ),
  evidence_status = claims$evidence_status,
  primary_evidence = c(
    "10D synthesis; 08C2 paired sample-level validation; 11F external source summary",
    "10D candidate specificity and 11F GSE130973 directional support",
    "10D module result and 11F actomyosin candidate discordance",
    "08C2 negative PIEZO2 versus 11F positive PIEZO2",
    "10D failed competition-pair audit",
    "10A provenance audit and 11F source boundary decision",
    "12C functional-validation specification"
  ),
  permitted_result_language = c(
    "associated with; supports a bounded hypothesis; reproducible at module level",
    "prioritized candidate; partial support; specificity remains unresolved",
    "pathway-level support; candidate-gene discordance; requires functional confirmation",
    "source-dependent; context-specific; secondary branch",
    "not fully specific; competing state programs remain active",
    "external triangulation; not direct replication",
    "not established; requires controlled perturbation and donor replication"
  ),
  prohibited_result_language = c(
    "causes mechanosensitivity; universal signature",
    "mechanosensitivity-specific integrin signature",
    "conserved actomyosin gene signature",
    "universal PIEZO2 mechanism",
    "independent of fibrosis or inflammation",
    "direct F7/F8 replication",
    "causal mechanism demonstrated"
  ),
  stringsAsFactors = FALSE
)

figure_table_blueprint <- data.frame(
  item_id = c("FIG1", "FIG2", "FIG3", "FIG4", "FIG5", "TAB1", "TAB2", "TAB3"),
  item_type = c("Figure", "Figure", "Figure", "Figure", "Figure", "Table", "Table", "Table"),
  title = c(
    "Study design, datasets, and provenance boundaries",
    "Discovery of the fibroblast-associated mechanotransduction program",
    "Paired sample-level validation and competitor robustness",
    "Independent GSE130973 external triangulation",
    "Evidence boundary model and proposed functional validation",
    "Dataset and sample inventory",
    "Frozen candidate/module results with directional support",
    "Claim boundaries, risks, and validation endpoints"
  ),
  main_message = c(
    "The project combines discovery, source-specific validation, and independent triangulation without claiming direct F7/F8 equivalence.",
    "The strongest result is a module-level ECM/integrin-cytoskeletal association.",
    "The core module is supported, but competition robustness and PIEZO2 consistency remain limitations.",
    "The independent source supports partial integrin evidence and demonstrates actomyosin/PIEZO2 discordance.",
    "Functional perturbation is required to move beyond CAUTION.",
    "Make source overlap, donor structure, labels, and data availability explicit.",
    "Separate module-level, candidate-level, and context-dependent evidence.",
    "Prevent overclaiming and define the next experimental decision points."
  ),
  stringsAsFactors = FALSE
)

write.csv(evidence_file_inventory, file.path(result_dir, "13A_evidence_file_inventory_v1.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(claim_to_evidence, file.path(result_dir, "13A_claim_to_evidence_map_v1.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(figure_table_blueprint, file.path(result_dir, "13A_figure_table_blueprint_v1.csv"), row.names = FALSE, fileEncoding = "UTF-8")

decision_lines <- c(
  "# Step 13A manuscript-ready evidence pack",
  "",
  "## Recommended manuscript positioning",
  "",
  "This study should be written as an evidence-bounded human single-cell reanalysis with external triangulation. The primary contribution is the identification of a fibroblast-associated, multicomponent ECM/integrin-cytoskeletal mechanotransduction program. The study does not establish causal mechanosensitivity, universal PIEZO2 involvement, or direct F7/F8 replication.",
  "",
  "## Results narrative order",
  "",
  "1. Establish dataset provenance, donor/sample structure, and label limitations.",
  "2. Present the discovery program and candidate modules.",
  "3. Present paired sample-level validation, including the failed competition-robustness gate.",
  "4. Present GSE130973 as non-overlapping external triangulation, emphasizing partial integrin support and actomyosin/PIEZO2 discordance.",
  "5. End with the revised mechanistic model and the functional-validation requirements.",
  "",
  "## Evidence-grade language",
  "",
  "Use association, enrichment, directional support, triangulation, and context-dependent. Avoid causal, universal, specific, direct replication, and mechanosensitivity demonstrated.",
  "",
  "## Required limitations",
  "",
  "- Observational single-cell data cannot establish causality.",
  "- F7/F8 anatomical and annotation equivalence remains unresolved.",
  "- Integrin support is not specific against ECM/TGF/inflammation-related states.",
  "- Frozen actomyosin genes are discordant across sources.",
  "- PIEZO2 is source-dependent.",
  "- GSE175817 lacks resolved donor-column mapping and is not a donor-level validation source.",
  "",
  "## Material Passport",
  "",
  "- Inputs: Step 12A, Step 12B, Step 12C, Step 10D, Step 08C2, and Step 11F outputs.",
  "- Transformation: manuscript evidence mapping and figure/table planning.",
  "- New data downloaded: none.",
  "- Output files: evidence inventory, claim-to-evidence map, and figure/table blueprint.",
  "- Current evidence grade: CAUTION."
)
writeLines(
  decision_lines,
  file.path(result_dir, "13A_manuscript_evidence_pack_v1.md"),
  useBytes = TRUE
)

message("Step 13A manuscript evidence pack completed.")
message("Evidence pack: ", file.path(result_dir, "13A_manuscript_evidence_pack_v1.md"))
message("Claim map: ", file.path(result_dir, "13A_claim_to_evidence_map_v1.csv"))
message("Figure/table blueprint: ", file.path(result_dir, "13A_figure_table_blueprint_v1.csv"))
