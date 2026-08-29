options(stringsAsFactors = FALSE)

# Step 18B: QA for the second-round revised manuscript.
# This script makes only evidence-preserving wording corrections and creates
# an auditable QA report. It does not alter any numerical result.

project_dir <- "."
input_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "18A_second_round_revised_manuscript"
)
output_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "18B_second_round_manuscript_QA"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_manuscript <- file.path(input_dir, "18A_second_round_revised_manuscript_v2.md")
input_claim_trace <- file.path(input_dir, "18A_claim_evidence_trace_v2.csv")
required_inputs <- c(input_manuscript, input_claim_trace)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 18B input(s): ", paste(missing_inputs, collapse = "; "))
}

read_text <- function(path) paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
write_text <- function(x, path) writeLines(x, path, useBytes = TRUE)
safe_write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
replace_all <- function(text, old, new) gsub(old, new, text, fixed = TRUE)
count_fixed <- function(text, pattern) {
  hit <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  if (length(hit) == 1L && hit[[1L]] == -1L) 0L else length(hit)
}
count_heading <- function(text, heading) {
  hit <- gregexpr(paste0("(?m)^", heading, "[[:space:]]*$"), text, perl = TRUE)[[1L]]
  if (length(hit) == 1L && hit[[1L]] == -1L) 0L else length(hit)
}

manuscript <- read_text(input_manuscript)

# Evidence-preserving scope/status corrections.
manuscript <- replace_all(
  manuscript,
  "A bounded human single-cell analysis identifies a fibroblast-associated ECM–integrin–cytoskeletal mechanotransduction program",
  "A bounded human transcriptomic reanalysis identifies a fibroblast-associated ECM–integrin–cytoskeletal mechanotransduction program"
)
manuscript <- replace_all(
  manuscript,
  "**Manuscript status:** initial evidence-grounded draft; current evidence grade CAUTION.",
  "**Manuscript status:** second-round evidence-grounded revised draft; current evidence grade CAUTION."
)
manuscript <- replace_all(
  manuscript,
  "multi-stage reanalysis of publicly available human single-cell transcriptomic resources",
  "multi-stage reanalysis of publicly available human transcriptomic resources"
)

# Structural and boundary checks.
checks <- data.frame(
  check_id = c(
    "single_conclusion_heading",
    "methods_marker_pair",
    "results_marker_pair",
    "discussion_marker_pair",
    "s2_candidate_coverage",
    "s3_htm_candidate_coverage",
    "s3_wi38_candidate_coverage",
    "overall_caution",
    "no_unqualified_cell_intrinsic_claim",
    "no_unqualified_causal_claim"
  ),
  check = c(
    "Exactly one Conclusion heading is present.",
    "The Step 17F Methods block is present exactly once and closed.",
    "The Step 17F Results block is present exactly once and closed.",
    "The Step 17F Discussion block is present exactly once and closed.",
    "S2 reports 15/15 frozen-candidate coverage.",
    "GSE123100 reports 14/15 candidate coverage.",
    "GSE276045 reports 15/15 candidate coverage.",
    "The overall evidence grade remains CAUTION.",
    "Cell-intrinsic wording appears only as a limitation or prohibited claim.",
    "Causal wording appears only as a limitation or prohibited claim."
  ),
  result = c(
    count_heading(manuscript, "## Conclusion") == 1L,
    count_fixed(manuscript, "<!-- STEP17F_METHODS_START -->") == 1L && count_fixed(manuscript, "<!-- STEP17F_METHODS_END -->") == 1L,
    count_fixed(manuscript, "<!-- STEP17F_RESULTS_START -->") == 1L && count_fixed(manuscript, "<!-- STEP17F_RESULTS_END -->") == 1L,
    count_fixed(manuscript, "<!-- STEP17F_DISCUSSION_START -->") == 1L && count_fixed(manuscript, "<!-- STEP17F_DISCUSSION_END -->") == 1L,
    grepl("15/15 frozen-candidate coverage", manuscript, fixed = TRUE),
    grepl("GSE123100 contained 14/15", manuscript, fixed = TRUE),
    grepl("GSE276045 contained 15/15", manuscript, fixed = TRUE),
    grepl("Overall evidence grade: CAUTION", manuscript, fixed = TRUE) &&
      grepl("Overall evidence grade: CAUTION", manuscript, fixed = TRUE),
    !grepl("cell-intrinsic mechanism was demonstrated", manuscript, fixed = TRUE) &&
      !grepl("cell-intrinsic mechanism is established", manuscript, fixed = TRUE),
    !grepl("causal mechanosensitivity was demonstrated", manuscript, fixed = TRUE) &&
      !grepl("causal mechanosensitivity is established", manuscript, fixed = TRUE)
  ),
  stringsAsFactors = FALSE
)

# Separate boundary-language inventory for human review.
boundary_inventory <- data.frame(
  term = c("cell-intrinsic", "causal mechanosensitivity", "direct F7/F8 replication", "CAUTION"),
  occurrence_count = c(
    count_fixed(tolower(manuscript), "cell-intrinsic"),
    count_fixed(tolower(manuscript), "causal mechanosensitivity"),
    count_fixed(tolower(manuscript), "direct F7/F8 replication"),
    count_fixed(manuscript, "CAUTION")
  ),
  interpretation_rule = c(
    "Must remain negated or explicitly prohibited; not a positive claim.",
    "Must remain negated or explicitly prohibited; not a positive claim.",
    "Must remain negated or explicitly prohibited; not a positive claim.",
    "Must remain the overall evidence grade."
  ),
  stringsAsFactors = FALSE
)

qa_pass <- all(checks$result)
output_manuscript <- file.path(output_dir, "18B_second_round_QA_revised_manuscript_v3.md")
write_text(manuscript, output_manuscript)
safe_write_csv(checks, file.path(output_dir, "18B_manuscript_QA_checks_v1.csv"))
safe_write_csv(boundary_inventory, file.path(output_dir, "18B_claim_boundary_inventory_v1.csv"))

qa_lines <- c(
  "# Step 18B second-round manuscript QA",
  "",
  paste0("QA gate: ", if (qa_pass) "PASS_WITH_AUTHOR_NUMERIC_VERIFICATION" else "HOLD"),
  "",
  "## Corrections applied",
  "",
  "- Changed the title from single-cell analysis to transcriptomic reanalysis to reflect the bulk WI-38 stiffness source.",
  "- Updated the manuscript status from initial draft to second-round evidence-grounded revised draft.",
  "- Updated the general Methods scope from single-cell transcriptomic resources to transcriptomic resources.",
  "- No numerical result, frozen candidate, threshold, or evidence grade was changed.",
  "",
  "## Required author verification",
  "",
  "- Verify all numerical values against the Step 17C, 17E and 17H CSV files.",
  "- Confirm that every occurrence of cell-intrinsic and causal mechanosensitivity remains a boundary statement.",
  "- Add the final public code repository URL before submission.",
  "- Confirm the target journal, word limits, figure limits and reference style before Word/PDF production.",
  "",
  "## Output",
  "",
  "The v3 manuscript is a QA-corrected copy; the 18A v2 manuscript remains preserved."
)
write_text(qa_lines, file.path(output_dir, "18B_manuscript_QA_report_v1.md"))

message("Step 18B manuscript QA completed.")
message("QA gate: ", if (qa_pass) "PASS_WITH_AUTHOR_NUMERIC_VERIFICATION" else "HOLD")
message("QA manuscript: ", output_manuscript)
message("QA report: ", file.path(output_dir, "18B_manuscript_QA_report_v1.md"))
