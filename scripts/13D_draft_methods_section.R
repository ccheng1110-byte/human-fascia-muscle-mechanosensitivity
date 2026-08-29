options(stringsAsFactors = FALSE)

# Step 13D: evidence-grounded Methods section draft.
# No new expression analysis, hypothesis testing, or data download.

project_dir <- "."
foundation_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13B_researchwrite_foundation"
)
review_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13C_foundation_review"
)
result_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13D_methods_draft"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  file.path(foundation_dir, "01_research_canon.md"),
  file.path(foundation_dir, "04_section_contracts.md"),
  file.path(review_dir, "13C_foundation_review_v1.md")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 13D input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

methods_lines <- c(
  "# Methods",
  "",
  "## Study design and analytical scope",
  "",
  "We conducted a provenance-aware, multi-stage reanalysis of publicly available human single-cell transcriptomic resources to evaluate whether fascia-like or myofibroblast-enriched fibroblast states are associated with a multicomponent extracellular-matrix/integrin–cytoskeletal mechanotransduction program. The workflow comprised discovery analysis in GSE173252, paired sample-level validation using the PRJNA607098-derived atlas representation, and independent external triangulation using GSE130973. The analysis was observational and was designed to separate module-level association from candidate-gene-level direction and from causal mechanosensitivity.",
  "",
  "No new human samples were collected. The present work reanalyzed public data and did not re-adjudicate the ethics, consent, or clinical metadata procedures of the original studies. All conclusions were constrained by the available accession-level provenance, sample identifiers, cell-state labels, and donor/sample structure.",
  "",
  "## Data sources and provenance audit",
  "",
  "The discovery resource was GSE173252, for which a processed mesenchyme object and supplementary files were available locally under the project data directory. The primary external validation used the PRJNA607098 representation of the 2025 human skin fibroblast atlas. Because the full atlas H5AD was approximately 25.36 GiB, only the required frozen gene-panel chunks were streamed from the public Zarr representation; the full H5AD was not downloaded.",
  "",
  "A provenance and label audit was performed before external validation. Sample identifiers, cell-state labels, source accessions, and the relationship between atlas F7 and paper F8 labels were recorded. The atlas source was treated as external triangulation rather than direct F7/F8 replication because one-to-one anatomical and analytical equivalence was not established.",
  "",
  "An additional independent-source screen evaluated GSE175817, GSE130973, GSE202352, GSE228421, and GSE138669. GSE175817 was classified as a gene-by-sample expression matrix, but donor-column mapping remained unresolved and it was therefore excluded from donor-level validation. GSE130973 was retained for exploratory external triangulation because it contained a processed Seurat object with subject metadata and candidate fibroblast-state clusters.",
  "",
  "## Discovery object handling and cluster reproduction",
  "",
  "The legacy GSE173252 processed object was read from the compressed RDS file using an explicit double-gzip reader. The RNA counts matrix, log-normalized expression matrix, and cell metadata were extracted from the object attributes. Cell and feature identifiers were checked before analysis, and metadata were reordered when necessary to match the RNA matrix columns. The original mesenchymal cluster structure was reproduced using the available cluster and sample annotations; the reproduction step did not reinterpret the original study as a new biological experiment.",
  "",
  "## Mechanotransduction program and candidate scoring",
  "",
  "The mechanotransduction program was defined in advance in the project gene-set registry. Module scores were calculated from the available log-normalized expression data with expression-bin-matched control genes. The scoring configuration used 20 control genes per target and 24 expression bins with a fixed random seed (20260823) for reproducibility. Module-level summaries and candidate-gene evidence were aggregated at the sample/cluster level rather than treating individual cells as independent biological replicates.",
  "",
  "The frozen candidate panel was evaluated separately from the broader modules. Candidate-level direction was recorded as descriptive support according to the prespecified project rules. This separation was maintained because coordinated module positivity does not imply that every individual candidate gene is conserved across tissues or datasets.",
  "",
  "## Specificity and competitor robustness",
  "",
  "Specificity audits compared the mechanotransduction program with competing state programs, including extracellular-matrix remodeling, TGF/fibrosis, inflammation, hypoxia, and cell-cycle-related programs. The competition-robustness result was retained as a formal evidence boundary rather than being treated as a nuisance result. A failed competition gate was interpreted as evidence that the observed integrin-associated program could not yet be separated from generic fibroblast activation or remodeling states.",
  "",
  "## PRJNA607098 paired sample-level validation",
  "",
  "Sample metadata were audited before expression extraction. The sample index was reconstructed from the official atlas observation-index format, and 12 eligible paired sample-state units were identified. The required frozen gene panel was streamed from the remote Zarr representation. Cells were aggregated within sample before comparison; individual cells were not treated as independent replicates. The paired validation evaluated frozen candidate directions, module summaries, leave-one-sample-out robustness, and competitor robustness using the prespecified project gates.",
  "",
  "The primary validation result was interpreted at three levels: core module support, candidate-gene support, and specificity/competition robustness. PIEZO2 was treated as a channel-specific secondary component rather than a required criterion for the core program.",
  "",
  "## GSE130973 external triangulation",
  "",
  "The GSE130973 processed object was loaded from the local compressed RDS file. The object audit recorded 15,457 metadata rows, five subjects, 17 clusters, and 15 of 15 frozen genes present. A legacy Seurat compatibility warning concerning the absent `images` slot was recorded as a non-fatal object-validity issue because the required expression and metadata components were accessible for the audit.",
  "",
  "Candidate fibroblast-state clusters were selected using the predefined top-five cluster rule based on canonical fibroblast marker enrichment and subject coverage. The selected clusters were present across subjects S1–S5. The exploratory frozen-program audit then summarized module-level and candidate-level directions within this marker-defined fibroblast state. GSE130973 was not treated as a direct anatomical replication of F7/F8.",
  "",
  "## Evidence synthesis and interpretation rules",
  "",
  "The final evidence synthesis integrated provenance, sample structure, module-level results, candidate-level direction, competition robustness, and cross-source concordance. The interpretation classes were assigned before converting the results into manuscript claims. The current evidence grade was retained as CAUTION because the evidence is observational, competition robustness was incomplete, PIEZO2 was source-dependent, actomyosin candidate genes were discordant, and direct F7/F8 equivalence remained unresolved.",
  "",
  "The revised primary hypothesis was therefore stated at the module level: fascia-like or myofibroblast-enriched fibroblast states may exhibit a multicomponent ECM/integrin–cytoskeletal mechanotransduction program whose individual genes and channel-specific components are context-dependent. Functional mechanical perturbation, donor-level replication, and orthogonal readouts are required before causal mechanosensitivity can be claimed.",
  "",
  "## Reproducibility and analysis records",
  "",
  "All analysis scripts, configuration files, intermediate audits, decision files, and machine-readable summaries were stored under the project directory `.`. The manuscript evidence pack and research foundation map each planned claim to its corresponding result file. No post hoc candidate replacement or threshold relaxation was authorized after the final evidence synthesis.",
  "",
  "## Methods limitations",
  "",
  "The Methods describe a public-data reanalysis rather than a prospective functional experiment. Some sources were processed objects rather than raw FASTQ files, and source-specific label and donor structures were not interchangeable. The analysis therefore supports bounded association and external triangulation, but not causal inference or direct F7/F8 replication."
)

writeLines(
  methods_lines,
  file.path(result_dir, "13D_methods_draft_v1.md"),
  useBytes = TRUE
)

trace_table <- data.frame(
  methods_section = c(
    "Study design and analytical scope",
    "Data sources and provenance audit",
    "Discovery object handling and cluster reproduction",
    "Mechanotransduction program and candidate scoring",
    "Specificity and competitor robustness",
    "PRJNA607098 paired sample-level validation",
    "GSE130973 external triangulation",
    "Evidence synthesis and interpretation rules",
    "Reproducibility and analysis records"
  ),
  supporting_script_or_output = c(
    "13B_researchwrite_foundation/00_scope.md; 13B_researchwrite_foundation/01_research_canon.md",
    "10A_audit_skin_fibroblast_atlas_provenance_and_labels.R; 10D_synthesize_evidence_grade.R",
    "05_reproduce_GSE173252_mesenchyme_clusters.R",
    "06_score_GSE173252_mechanosensitivity_program_v1.R",
    "07_audit_GSE173252_mechanosensitivity_specificity.R; 10D_synthesize_evidence_grade.R",
    "08C1_audit_PRJNA607098_sample_metadata.R; 08C2_test_PRJNA607098_F7_sample_level.R",
    "11D2_audit_GSE130973_cluster_fibroblast_state.R; 11E_audit_GSE130973_frozen_mechanotransduction_program.R",
    "10D_synthesize_evidence_grade.R; 11F_summarize_GSE130973_external_validation.R; 12A_synthesize_final_evidence_and_revise_hypothesis.R",
    "13A_build_manuscript_evidence_pack.R; 13B_build_manuscript_foundation_and_outline.R"
  ),
  claim_boundary = c(
    "Observational association; no causal inference.",
    "F7/F8 equivalence unresolved; GSE130973 is triangulation.",
    "Reproduction of available processed object; not a new clustering experiment.",
    "Module and candidate scoring are distinct evidence levels.",
    "Competition failure is retained as a limitation.",
    "Sample is the unit; cells are not independent replicates.",
    "Marker-defined external state; not direct F7/F8 replication.",
    "Overall grade remains CAUTION.",
    "All files and scripts remain locally traceable."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  trace_table,
  file.path(result_dir, "13D_methods_evidence_trace_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Step 13D Methods draft completed.")
message("Methods draft: ", file.path(result_dir, "13D_methods_draft_v1.md"))
message("Evidence trace: ", file.path(result_dir, "13D_methods_evidence_trace_v1.csv"))
