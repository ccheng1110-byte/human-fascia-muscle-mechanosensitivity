options(stringsAsFactors = FALSE)

# Step 17B: S1 source and analyzability audit.
#
# This is a provenance/feasibility audit only. It does not read new expression
# chunks, perform a hypothesis test, or treat cells as independent replicates.
# PRJNA607098 is audited from the already completed metadata and targeted-gene
# inventories. GSE130973 is audited from the already completed subject,
# cluster, and frozen-program audits.

project_dir <- "."
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17B_S1_source_feasibility_audit"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

registry_path <- file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
)
candidate_path <- file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
)

prjna_eligibility_path <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025",
  "08C1_sample_metadata_audit_PRJNA607098",
  "PRJNA607098_sample_F7_eligibility_v1.csv"
)
prjna_inventory_path <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025",
  "08B3_F7_PRJNA607098",
  "skin_fibroblast_atlas_gene_panel_inventory_v3.csv"
)

gse130973_subject_path <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11D2_GSE130973_cluster_fibroblast_audit",
  "GSE130973_subject_age_cell_inventory_v2.csv"
)
gse130973_cluster_path <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11D2_GSE130973_cluster_fibroblast_audit",
  "GSE130973_cluster_fibroblast_marker_audit_v2.csv"
)
gse130973_cluster_config_path <- file.path(
  project_dir, "config", "GSE130973_candidate_fibroblast_clusters_v1.csv"
)
gse130973_module_summary_path <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11E_GSE130973_frozen_program_audit",
  "GSE130973_fibroblast_state_module_summary_v1.csv"
)
gse130973_decision_path <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11E_GSE130973_frozen_program_audit",
  "GSE130973_fibroblast_state_audit_decision_v1.md"
)

required_inputs <- c(
  registry_path, candidate_path,
  prjna_eligibility_path, prjna_inventory_path,
  gse130973_subject_path, gse130973_cluster_path,
  gse130973_cluster_config_path, gse130973_module_summary_path,
  gse130973_decision_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required Step 17B input(s): ",
    paste(missing_inputs, collapse = "; ")
  )
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

registry <- read.csv(registry_path, check.names = FALSE)
candidate_panel <- read.csv(candidate_path, check.names = FALSE)
if (!all(c("module", "genes") %in% names(registry))) {
  stop("Registry must contain module and genes columns.")
}
if (!"gene" %in% names(candidate_panel)) {
  stop("Frozen candidate panel must contain a gene column.")
}

registry_gene_rows <- lapply(seq_len(nrow(registry)), function(i) {
  trimws(strsplit(as.character(registry$genes[[i]]), ";", fixed = TRUE)[[1L]])
})
registry_genes <- sort(unique(unlist(registry_gene_rows, use.names = FALSE)))
candidate_genes <- unique(trimws(as.character(candidate_panel$gene)))

# -----------------------------
# PRJNA607098: existing metadata and targeted inventory
# -----------------------------
prjna_eligibility <- read.csv(prjna_eligibility_path, check.names = FALSE)
required_prjna_fields <- c(
  "sample_id", "F7_cells", "other_state_cells",
  "eligible_for_paired_sample_state_analysis"
)
missing_prjna_fields <- setdiff(required_prjna_fields, names(prjna_eligibility))
if (length(missing_prjna_fields) > 0L) {
  stop(
    "PRJNA607098 eligibility file is missing field(s): ",
    paste(missing_prjna_fields, collapse = ", ")
  )
}

prjna_eligibility$eligible <- as_bool(
  prjna_eligibility$eligible_for_paired_sample_state_analysis
)
prjna_sample_count <- nrow(prjna_eligibility)
prjna_eligible_count <- sum(prjna_eligibility$eligible)
prjna_min_target_cells <- min(prjna_eligibility$F7_cells, na.rm = TRUE)
prjna_min_comparator_cells <- min(prjna_eligibility$other_state_cells, na.rm = TRUE)
prjna_sample_gate <- prjna_sample_count == 12L && prjna_eligible_count >= 6L

prjna_inventory <- read.csv(prjna_inventory_path, check.names = FALSE)
required_inventory_fields <- c("gene", "found_in_atlas")
missing_inventory_fields <- setdiff(required_inventory_fields, names(prjna_inventory))
if (length(missing_inventory_fields) > 0L) {
  stop(
    "PRJNA607098 inventory is missing field(s): ",
    paste(missing_inventory_fields, collapse = ", ")
  )
}
prjna_inventory$found <- as_bool(prjna_inventory$found_in_atlas)
prjna_found_genes <- unique(
  trimws(as.character(prjna_inventory$gene[prjna_inventory$found]))
)
prjna_candidate_present <- candidate_genes %in% prjna_found_genes
prjna_candidate_coverage <- sum(prjna_candidate_present)
prjna_candidate_gate <- all(prjna_candidate_present)

# The existing PRJNA extraction was intentionally targeted. It is not evidence
# that absent genes are absent from the atlas. The remaining registry genes need
# a new Step 17C cell-level stream before competition-adjusted S1 inference.
prjna_registry_present_in_existing_targeted_extract <-
  registry_genes %in% prjna_found_genes
prjna_registry_targeted_count <- sum(
  prjna_registry_present_in_existing_targeted_extract
)
prjna_unretrieved_registry_genes <- registry_genes[
  !prjna_registry_present_in_existing_targeted_extract
]

prjna_gene_rows <- data.frame(
  source = "PRJNA607098",
  gene = registry_genes,
  is_frozen_candidate = registry_genes %in% candidate_genes,
  present_in_existing_targeted_extract =
    prjna_registry_present_in_existing_targeted_extract,
  retrieval_status = ifelse(
    prjna_registry_present_in_existing_targeted_extract,
    "available_from_existing_targeted_extract",
    "requires_new_targeted_cell_level_stream"
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(
  prjna_gene_rows,
  file.path(result_dir, "PRJNA607098_S1_registry_gene_retrieval_status_v1.csv")
)

prjna_summary <- data.frame(
  source = "PRJNA607098",
  source_role = "S1 primary single-cell source",
  sample_or_donor_unit = "reconstructed sample ID; not automatically donor-independent",
  samples = prjna_sample_count,
  eligible_samples = prjna_eligible_count,
  minimum_F7_cells_per_sample = prjna_min_target_cells,
  minimum_nonF7_cells_per_sample = prjna_min_comparator_cells,
  sample_metadata_gate = prjna_sample_gate,
  frozen_candidates_present = prjna_candidate_coverage,
  frozen_candidates_total = length(candidate_genes),
  frozen_candidate_coverage_gate = prjna_candidate_gate,
  registry_genes_in_existing_targeted_extract =
    prjna_registry_targeted_count,
  registry_genes_total = length(registry_genes),
  full_registry_cell_level_stream_ready = FALSE,
  source_independence_status =
    "sample-level atlas source; donor independence not established",
  stringsAsFactors = FALSE
)
safe_write_csv(
  prjna_summary,
  file.path(result_dir, "PRJNA607098_S1_source_feasibility_v1.csv")
)

# -----------------------------
# GSE130973: existing subject and cluster audits
# -----------------------------
gse_subjects <- read.csv(gse130973_subject_path, check.names = FALSE)
gse_clusters <- read.csv(gse130973_cluster_path, check.names = FALSE)
gse_cluster_config <- read.csv(gse130973_cluster_config_path, check.names = FALSE)
gse_module_summary <- read.csv(gse130973_module_summary_path, check.names = FALSE)

required_gse_subject_fields <- c("subj", "cells")
if (length(setdiff(required_gse_subject_fields, names(gse_subjects))) > 0L) {
  stop("GSE130973 subject inventory is incomplete.")
}
required_gse_cluster_fields <- c("cluster_id", "subjects_present")
if (length(setdiff(required_gse_cluster_fields, names(gse_clusters))) > 0L) {
  stop("GSE130973 cluster audit is incomplete.")
}
if (!all(c("cluster_id", "subjects_present") %in% names(gse_cluster_config))) {
  stop("GSE130973 candidate cluster configuration is incomplete.")
}

gse_subject_ids <- sort(unique(as.character(gse_subjects$subj)))
gse_candidate_clusters <- as.character(gse_cluster_config$cluster_id)
gse_cluster_ids <- as.character(gse_clusters$cluster_id)
gse_candidate_cluster_rows <- gse_clusters[
  gse_clusters$cluster_id %in% gse_candidate_clusters,
  , drop = FALSE
]
gse_subject_gate <- length(gse_subject_ids) == 5L
gse_cluster_gate <- length(gse_candidate_clusters) == 5L &&
  all(gse_candidate_clusters %in% gse_cluster_ids) &&
  all(gse_cluster_config$subjects_present >= 5L)

# Recover the previous 11E unique registry coverage from its decision record.
gse_decision_lines <- readLines(gse130973_decision_path, encoding = "UTF-8")
gse_decision_text <- paste(gse_decision_lines, collapse = "\n")
coverage_match <- regexec(
  "Frozen registry genes present: ([0-9]+)/([0-9]+)",
  gse_decision_text,
  perl = TRUE
)
coverage_capture <- regmatches(gse_decision_text, coverage_match)[[1L]]
if (length(coverage_capture) != 3L) {
  stop("Could not recover GSE130973 frozen registry coverage from Step 11E.")
}
gse_registry_present <- as.integer(coverage_capture[[2L]])
gse_registry_total <- as.integer(coverage_capture[[3L]])
gse_registry_coverage_gate <- gse_registry_present == gse_registry_total

gse_summary <- data.frame(
  source = "GSE130973",
  source_role = "S1 independent skin-fibroblast supplementary source",
  sample_or_donor_unit = "subject field subj",
  subjects = length(gse_subject_ids),
  subject_ids = paste(gse_subject_ids, collapse = ";"),
  subject_mapping_gate = gse_subject_gate,
  frozen_candidate_cluster_count = length(gse_candidate_clusters),
  candidate_clusters = paste(gse_candidate_clusters, collapse = ";"),
  candidate_cluster_gate = gse_cluster_gate,
  registry_genes_present_in_existing_audit = gse_registry_present,
  registry_genes_total_in_existing_audit = gse_registry_total,
  registry_coverage_gate = gse_registry_coverage_gate,
  source_independence_status =
    "five-subject skin source; separate study, not pooled with PRJNA607098",
  stringsAsFactors = FALSE
)
safe_write_csv(
  gse_summary,
  file.path(result_dir, "GSE130973_S1_source_feasibility_v1.csv")
)

safe_write_csv(
  gse_candidate_cluster_rows,
  file.path(result_dir, "GSE130973_S1_candidate_cluster_inventory_v1.csv")
)

# -----------------------------
# Integrated Step 17B decision
# -----------------------------
metadata_and_state_gate <- prjna_sample_gate &&
  gse_subject_gate && gse_cluster_gate
candidate_coverage_gate <- prjna_candidate_gate
full_s1_stream_ready <- metadata_and_state_gate &&
  candidate_coverage_gate && gse_registry_coverage_gate &&
  prjna_registry_targeted_count == length(registry_genes)

decision <- if (full_s1_stream_ready) {
  "PASS_TO_STEP17C"
} else if (metadata_and_state_gate && candidate_coverage_gate) {
  "PARTIAL_PROCEED_TO_TARGETED_EXPRESSION_AUDIT"
} else {
  "HOLD_NOT_ANALYZABLE"
}

decision_lines <- c(
  "# Step 17B S1 source and feasibility audit",
  "",
  "## Decision",
  "",
  paste0("- Gate 17B: **", decision, "**."),
  "- This step is a provenance and analyzability audit only; no hypothesis test was performed.",
  "- Cells are not treated as independent biological replicates.",
  "",
  "## PRJNA607098",
  "",
  paste0("- Reconstructed sample units: ", prjna_sample_count, "."),
  paste0("- Eligible paired sample units: ", prjna_eligible_count, "."),
  paste0("- Sample metadata gate: ", prjna_sample_gate, "."),
  paste0("- Frozen candidate coverage in the existing targeted extract: ",
         prjna_candidate_coverage, "/", length(candidate_genes), "."),
  paste0("- Full registry genes available in the existing targeted extract: ",
         prjna_registry_targeted_count, "/", length(registry_genes), "."),
  "- Existing PRJNA gene inventory is targeted; genes not listed as present are not interpreted as biologically absent.",
  "- A new cell-level targeted stream is required before competition-adjusted S1 inference.",
  "",
  "## GSE130973",
  "",
  paste0("- Subjects: ", paste(gse_subject_ids, collapse = ", "), "."),
  paste0("- Subject mapping gate: ", gse_subject_gate, "."),
  paste0("- Frozen candidate fibroblast-state clusters: ",
         paste(gse_candidate_clusters, collapse = ", "), "."),
  paste0("- Candidate cluster/state gate: ", gse_cluster_gate, "."),
  paste0("- Existing frozen registry coverage: ", gse_registry_present, "/",
         gse_registry_total, "."),
  "- GSE130973 will remain a separate study-level supplementary source; it will not be pooled with PRJNA607098.",
  "",
  "## Step 17C entry conditions",
  "",
  "1. Stream the PRJNA607098 genes listed in `PRJNA607098_S1_registry_gene_retrieval_status_v1.csv` that are not available in the existing targeted extract.",
  "2. Retain the exact F7/non-F7 sample reconstruction and the 12 sample units.",
  "3. Perform cell-level co-detection descriptively, then use within-sample permutation and sample-level aggregation for inference.",
  "4. Keep detection rate, library size, cell-cycle, ECM/TGF, inflammation and hypoxia as prespecified competition controls.",
  "5. Do not label the result `cell-intrinsic` automatically; the output label must be a co-expression evidence level.",
  "",
  "## Evidence boundary",
  "",
  "- This audit does not establish donor independence for PRJNA607098.",
  "- It does not prove cell-intrinsic mechanism, causal mechanosensitivity, or fascia specificity.",
  "- The overall project evidence grade remains CAUTION.",
  "",
  "## Material Passport",
  "",
  "- Inputs: existing Step 08C1/08B3 PRJNA607098 audits, Step 11D2/11E GSE130973 audits, frozen candidate panel, and module registry.",
  "- Transformation: local provenance, sample/state, and gene-coverage feasibility checks.",
  "- New expression data downloaded: none.",
  "- Next step: Step 17C targeted cell-level expression stream and co-expression-level analysis."
)
decision_path <- file.path(
  result_dir, "Step17B_S1_source_feasibility_decision_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17B S1 source and feasibility audit completed.")
message("Gate 17B: ", decision)
message("PRJNA607098 sample units: ", prjna_sample_count,
        "; eligible: ", prjna_eligible_count)
message("PRJNA607098 frozen candidate coverage: ",
        prjna_candidate_coverage, "/", length(candidate_genes))
message("GSE130973 subjects: ", paste(gse_subject_ids, collapse = ", "))
message("GSE130973 registry coverage: ", gse_registry_present, "/",
        gse_registry_total)
message("Decision: ", decision_path)
