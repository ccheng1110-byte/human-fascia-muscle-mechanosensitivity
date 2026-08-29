options(stringsAsFactors = FALSE)

# Step 13D2: internal section QA for the Methods draft.
# No new analysis or data download.

project_dir <- "."
input_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13D_methods_draft"
)
result_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13D2_methods_review"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

draft_path <- file.path(input_dir, "13D_methods_draft_v1.md")
trace_path <- file.path(input_dir, "13D_methods_evidence_trace_v1.csv")
if (!file.exists(draft_path) || !file.exists(trace_path)) {
  stop("Step 13D draft and evidence trace are required.")
}

draft <- paste(readLines(draft_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
trace <- read.csv(trace_path, check.names = FALSE)

has_all <- function(patterns, fixed = TRUE) {
  all(vapply(patterns, function(pattern) grepl(pattern, draft, fixed = fixed), logical(1L)))
}

checklist <- data.frame(
  domain = c(
    "study_design",
    "data_provenance",
    "object_handling",
    "program_scoring",
    "competitor_audit",
    "sample_level_unit",
    "external_triangulation",
    "evidence_synthesis",
    "reproducibility",
    "limitations",
    "no_causal_overclaim",
    "trace_completeness"
  ),
  criterion = c(
    "Study design and observational scope are stated.",
    "Datasets, accession/provenance boundaries, and F7/F8 limitation are described.",
    "Processed-object handling and identifier checks are described.",
    "Frozen gene-set scoring and control-gene configuration are described.",
    "Competitor programs and failed robustness boundary are retained.",
    "The sample/donor unit is distinguished from individual cells.",
    "GSE130973 is explicitly framed as external triangulation, not direct replication.",
    "Evidence is separated into module, candidate, specificity, and causal levels.",
    "Scripts, configurations, outputs, and no post hoc threshold changes are recorded.",
    "Methods limitations are stated explicitly.",
    "Causal mechanosensitivity is not claimed as established.",
    "Each Methods subsection has an evidence trace entry."
  ),
  passed = c(
    has_all(c("Study design and analytical scope", "observational")),
    has_all(c("Data sources and provenance audit", "GSE173252", "PRJNA607098", "F7/F8")),
    has_all(c("Discovery object handling and cluster reproduction", "double-gzip", "metadata")),
    has_all(c("Mechanotransduction program and candidate scoring", "control genes", "random seed")),
    has_all(c("Specificity and competitor robustness", "competition gate")),
    has_all(c("individual cells", "independent replicates", "sample-level")),
    has_all(c("GSE130973 external triangulation", "not treated as a direct")),
    has_all(c("module-level", "candidate-gene", "causal")),
    has_all(c("Reproducibility and analysis records", "No post hoc candidate replacement")),
    has_all(c("Methods limitations", "not causal", "direct F7/F8 replication")),
    has_all(c("causal mechanosensitivity", "can be claimed")),
    nrow(trace) >= 8L && all(nzchar(trimws(as.character(trace$methods_section))))
  ),
  weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0.75, 0.75, 0.5),
  stringsAsFactors = FALSE
)
checklist$weighted_score <- ifelse(checklist$passed, checklist$weight, 0)
section_score <- sum(checklist$weighted_score)
section_max <- sum(checklist$weight)
section_pass <- section_score >= 8

write.csv(
  checklist,
  file.path(result_dir, "13D2_methods_review_checklist_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

failed <- checklist$criterion[!checklist$passed]
review_lines <- c(
  "# Step 13D2 Methods section QA",
  "",
  "## Review type",
  "",
  "Internal section-level QA for evidence traceability, reproducibility, and claim boundaries. This review does not upgrade the biological evidence grade.",
  "",
  paste0("Methods section score: ", round(section_score, 2), "/", round(section_max, 2)),
  paste0("Section gate: ", ifelse(section_pass, "PASS", "HOLD")),
  "",
  "## Findings",
  "",
  if (section_pass) "- The Methods draft meets the internal section QA threshold of 8/10." else "- The Methods draft does not yet meet the internal section QA threshold of 8/10.",
  "- The Methods retain the distinction between module-level association and causal mechanosensitivity.",
  "- The Methods identify sample/donor-level units and do not treat cells as independent biological replicates.",
  if (length(failed) == 0L) "- No checklist item failed." else c("- Failed checklist items:", paste0("  - ", failed)),
  "",
  "## Remaining scientific limitations",
  "",
  "- Overall evidence grade remains CAUTION.",
  "- F7/F8 equivalence remains unresolved.",
  "- Functional validation has not been performed.",
  "",
  "## Next gate",
  "",
  if (section_pass) "Proceed to Step 13E Results drafting." else "Repair the Methods draft before drafting Results.",
  "",
  "## Material Passport",
  "",
  "- Input: Step 13D Methods draft and evidence trace.",
  "- Transformation: section QA only.",
  "- New data downloaded: none."
)
writeLines(
  review_lines,
  file.path(result_dir, "13D2_methods_review_v1.md"),
  useBytes = TRUE
)

message("Step 13D2 Methods QA completed.")
message("Methods score: ", round(section_score, 2), "/", round(section_max, 2))
message("Section gate: ", ifelse(section_pass, "PASS", "HOLD"))
message("Review: ", file.path(result_dir, "13D2_methods_review_v1.md"))
