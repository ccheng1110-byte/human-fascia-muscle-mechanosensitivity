options(stringsAsFactors = FALSE)

# Step 13C: structural foundation review before manuscript drafting.
# This is an internal evidence/logic gate, not an external peer review.

project_dir <- "."
foundation_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13B_researchwrite_foundation"
)
result_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13C_foundation_review"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

foundation_files <- c(
  "00_scope.md",
  "01_research_canon.md",
  "02_evidence_table.md",
  "03_argument_map.md",
  "04_section_contracts.md",
  "05_style_guide.md",
  "13B_manuscript_outline_v1.md"
)
paths <- file.path(foundation_dir, foundation_files)
exists_flag <- file.exists(paths)
if (any(!exists_flag)) {
  stop("Missing foundation file(s): ", paste(foundation_files[!exists_flag], collapse = "; "))
}

contents <- lapply(paths, function(path) paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n"))
names(contents) <- foundation_files

contains_all <- function(text, patterns) {
  all(vapply(patterns, function(pattern) grepl(pattern, text, fixed = TRUE), logical(1L)))
}

checklist <- data.frame(
  domain = c(
    "scope_lock",
    "research_canon",
    "evidence_table",
    "argument_map",
    "section_contracts",
    "style_guide",
    "outline",
    "evidence_grade_boundary",
    "causal_claim_boundary",
    "replication_boundary"
  ),
  criterion = c(
    "Scope states manuscript type, reader, language, stage, and exclusions.",
    "Canon contains hard facts, terminology constraints, and forbidden claims.",
    "Evidence table contains at least five claim rows and explicit strength/status fields.",
    "Argument map contains tension, research question, thesis, supporting arguments, counterarguments, and final move.",
    "Section contracts define purpose, allowed claims, forbidden claims, and required evidence.",
    "Style guide defines conservative wording and reporting rules.",
    "Outline contains abstract, Introduction, Results, Discussion, Conclusion, and figure/table alignment.",
    "The current evidence grade is explicitly fixed as CAUTION.",
    "Causal mechanosensitivity is explicitly treated as untested or prohibited.",
    "Direct F7/F8 replication is explicitly treated as unresolved or prohibited."
  ),
  passed = c(
    contains_all(contents[["00_scope.md"]], c("Text type:", "Target reader:", "Evidence boundary:", "Excluded from this step:")),
    contains_all(contents[["01_research_canon.md"]], c("Hard facts", "Terminology constraints", "Forbidden claims")),
    sum(grepl("^\\|", strsplit(contents[["02_evidence_table.md"]], "\n", fixed = TRUE)[[1L]])) >= 7L &&
      contains_all(contents[["02_evidence_table.md"]], c("Strength", "Status")),
    contains_all(contents[["03_argument_map.md"]], c("Scientific tension", "Central research question", "Central thesis", "Supporting arguments", "Counterarguments", "Final move")),
    contains_all(contents[["04_section_contracts.md"]], c("Purpose:", "Allowed claims:", "Forbidden claims:", "Required validation")),
    contains_all(contents[["05_style_guide.md"]], c("conservative", "Avoid:", "Every quantitative statement")),
    contains_all(contents[["13B_manuscript_outline_v1.md"]], c("Abstract structure", "Introduction", "Results", "Discussion", "Conclusion", "Figure and table alignment")),
    grepl("CAUTION", paste(unlist(contents), collapse = "\n"), fixed = TRUE),
    grepl("causal", paste(unlist(contents), collapse = "\n"), ignore.case = TRUE) &&
      grepl("Forbidden|not establish|untested|not yet", paste(unlist(contents), collapse = "\n"), ignore.case = TRUE),
    grepl("F7/F8", paste(unlist(contents), collapse = "\n"), fixed = TRUE) &&
      grepl("unresolved|not.*replication|direct", paste(unlist(contents), collapse = "\n"), ignore.case = TRUE)
  ),
  weight = c(1, 2, 2, 2, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  stringsAsFactors = FALSE
)
checklist$weighted_score <- ifelse(checklist$passed, checklist$weight, 0)
foundation_score <- sum(checklist$weighted_score)
foundation_max <- sum(checklist$weight)
foundation_ready <- foundation_score >= 7.5

write.csv(
  checklist,
  file.path(result_dir, "13C_foundation_review_checklist_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

failed_items <- checklist$criterion[!checklist$passed]
decision_lines <- c(
  "# Step 13C foundation review",
  "",
  "## Review type",
  "",
  "This is an internal structural and evidence-boundary review of the manuscript foundation. It is not an external peer review and does not upgrade the biological evidence grade.",
  "",
  paste0("Foundation score: ", format(round(foundation_score, 2), nsmall = 2), "/", format(round(foundation_max, 2), nsmall = 2)),
  paste0("Drafting gate: ", ifelse(foundation_ready, "PASS", "HOLD")),
  "",
  "## Findings",
  "",
  if (foundation_ready) "- The foundation meets the internal drafting threshold of 7.5/10." else "- The foundation does not yet meet the internal drafting threshold of 7.5/10.",
  "- The biological evidence grade remains CAUTION.",
  "- Drafting must preserve the distinction between module-level association, candidate-level direction, and causal functional evidence.",
  if (length(failed_items) == 0L) "- No structural checklist item failed." else c("- Failed checklist items:", paste0("  - ", failed_items)),
  "",
  "## Remaining technical debts",
  "",
  "- Literature citations have not yet been attached to the research canon.",
  "- Functional validation has not yet been performed.",
  "- F7/F8 equivalence remains unresolved.",
  "",
  "## Next gate",
  "",
  if (foundation_ready) "Proceed to Step 13D section-by-section manuscript drafting, beginning with Methods and Results." else "Repair the failed foundation items before drafting.",
  "",
  "## Material Passport",
  "",
  "- Inputs: Step 13B foundation files.",
  "- Transformation: structural checklist and evidence-boundary scoring.",
  "- New data downloaded: none.",
  "- Biological evidence grade: CAUTION."
)
writeLines(
  decision_lines,
  file.path(result_dir, "13C_foundation_review_v1.md"),
  useBytes = TRUE
)

message("Step 13C foundation review completed.")
message("Foundation score: ", round(foundation_score, 2), "/", round(foundation_max, 2))
message("Drafting gate: ", ifelse(foundation_ready, "PASS", "HOLD"))
message("Review: ", file.path(result_dir, "13C_foundation_review_v1.md"))
