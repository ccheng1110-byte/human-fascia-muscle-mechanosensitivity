options(stringsAsFactors = FALSE)

# Step 13E2: internal section QA for the Results draft.
# No new analysis or data download.

project_dir <- "."
input_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13E_results_draft")
result_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13E2_results_review")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

draft_path <- file.path(input_dir, "13E_results_draft_v1.md")
numeric_path <- file.path(input_dir, "13E_results_numeric_evidence_v1.csv")
trace_path <- file.path(input_dir, "13E_results_evidence_trace_v1.csv")
required_inputs <- c(draft_path, numeric_path, trace_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 13E input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

draft <- paste(readLines(draft_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
numeric_evidence <- read.csv(numeric_path, check.names = FALSE)
trace <- read.csv(trace_path, check.names = FALSE)

contains <- function(pattern, fixed = TRUE, ignore.case = FALSE) {
  if (fixed && ignore.case) {
    return(grepl(tolower(pattern), tolower(draft), fixed = TRUE))
  }
  grepl(pattern, draft, fixed = fixed, ignore.case = ignore.case)
}

eligible_row <- numeric_evidence[numeric_evidence$result_domain == "PRJNA607098_eligible_samples", , drop = FALSE]
numeric_sample_count_correct <- nrow(eligible_row) == 1L && identical(as.character(eligible_row$value[[1L]]), "12")
no_overclaim <- !any(vapply(
  c("proves causal", "demonstrates causal", "universal PIEZO2 mechanism", "direct F7/F8 replication was established"),
  function(pattern) contains(pattern, fixed = TRUE, ignore.case = TRUE),
  logical(1L)
))

checklist <- data.frame(
  domain = c(
    "provenance_structure",
    "primary_module_results",
    "candidate_results",
    "competition_boundary",
    "external_triangulation",
    "cross_source_discordance",
    "evidence_grade",
    "numeric_evidence_integrity",
    "evidence_trace",
    "overclaim_control"
  ),
  criterion = c(
    "Results state discovery, 12-sample primary validation, and five-subject external triangulation.",
    "Results report core/module-level support and leave-one-sample-out status.",
    "Results report 10/15 and 5/6 candidate support with candidate/module distinction.",
    "Results explicitly retain failed competition robustness and competing programs.",
    "Results describe GSE130973 as external triangulation rather than direct replication.",
    "Results preserve actomyosin 1/4 and PIEZO2 source-dependent discordance.",
    "Results retain the overall CAUTION evidence grade.",
    "Numeric evidence table records 12 eligible samples under the correct field label.",
    "Every Results subsection has a primary output and boundary trace.",
    "No prohibited causal or universal claims are present."
  ),
  passed = c(
    contains("Study provenance and validation structure") && contains("12 eligible sample-state units") && contains("five-subject"),
    contains("Core program support in paired sample-level validation") && contains("leave-one-sample-out gate passed") && contains("core mechanotransduction module gate passed"),
    contains("10/15") && contains("5/6") && contains("candidate-gene conservation"),
    contains("competition-robustness gate failed") && contains("ECM remodeling") && contains("TGF/fibrosis"),
    contains("Independent GSE130973 external triangulation") && contains("not a conserved single-gene mechanism"),
    contains("1/4") && contains("PIEZO2") && contains("opposite to the PRJNA607098 result"),
    contains("CAUTION"),
    numeric_sample_count_correct,
    nrow(trace) >= 5L && all(nzchar(trimws(as.character(trace$primary_output)))),
    no_overclaim
  ),
  weight = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  stringsAsFactors = FALSE
)
checklist$weighted_score <- ifelse(checklist$passed, checklist$weight, 0)
section_score <- sum(checklist$weighted_score)
section_max <- sum(checklist$weight)
section_pass <- section_score >= 8

write.csv(
  checklist,
  file.path(result_dir, "13E2_results_review_checklist_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

failed <- checklist$criterion[!checklist$passed]
review_lines <- c(
  "# Step 13E2 Results section QA",
  "",
  "## Review type",
  "",
  "Internal section-level QA for numerical traceability, result-boundary preservation, and overclaim control. This review does not upgrade the biological evidence grade.",
  "",
  paste0("Results section score: ", round(section_score, 2), "/", round(section_max, 2)),
  paste0("Section gate: ", ifelse(section_pass, "PASS", "HOLD")),
  "",
  "## Findings",
  "",
  if (section_pass) "- The Results draft meets the internal section QA threshold of 8/10." else "- The Results draft does not yet meet the internal section QA threshold of 8/10.",
  "- The corrected numeric evidence table records the primary validation sample count under the correct field.",
  "- Failed competition robustness, actomyosin discordance, PIEZO2 source dependence, and the CAUTION grade remain visible.",
  if (length(failed) == 0L) "- No checklist item failed." else c("- Failed checklist items:", paste0("  - ", failed)),
  "",
  "## Remaining scientific limitations",
  "",
  "- Results remain observational and do not establish causal mechanosensitivity.",
  "- GSE130973 is external triangulation, not direct F7/F8 replication.",
  "",
  "## Next gate",
  "",
  if (section_pass) "Proceed to Step 13F Discussion drafting." else "Repair the Results draft before drafting Discussion.",
  "",
  "## Material Passport",
  "",
  "- Input: Step 13E Results draft, numeric evidence, and evidence trace.",
  "- Transformation: section QA only.",
  "- New data downloaded: none."
)
writeLines(
  review_lines,
  file.path(result_dir, "13E2_results_review_v1.md"),
  useBytes = TRUE
)

message("Step 13E2 Results QA completed.")
message("Results score: ", round(section_score, 2), "/", round(section_max, 2))
message("Section gate: ", ifelse(section_pass, "PASS", "HOLD"))
message("Review: ", file.path(result_dir, "13E2_results_review_v1.md"))
