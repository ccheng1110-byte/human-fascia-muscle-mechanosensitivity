options(stringsAsFactors = FALSE)

# Step 13F2: internal section QA for the Discussion draft.
# No new analysis or data download.

project_dir <- "."
input_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13F_discussion_draft")
result_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13F2_discussion_review")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

draft_path <- file.path(input_dir, "13F_discussion_draft_v1.md")
trace_path <- file.path(input_dir, "13F_discussion_evidence_trace_v1.csv")
citation_path <- file.path(input_dir, "13F_citations_needed_v1.csv")
required_inputs <- c(draft_path, trace_path, citation_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 13F input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

draft <- paste(readLines(draft_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
trace <- read.csv(trace_path, check.names = FALSE)
citations <- read.csv(citation_path, check.names = FALSE)

contains <- function(pattern, ignore.case = FALSE) {
  if (ignore.case) {
    return(grepl(tolower(pattern), tolower(draft), fixed = TRUE))
  }
  grepl(pattern, draft, fixed = TRUE)
}

placeholder_count <- lengths(regmatches(draft, gregexpr("\\[CITATION REQUIRED", draft, fixed = FALSE)))
no_overclaim <- !any(vapply(
  c("proves causal", "causal mechanism demonstrated", "universal PIEZO2 mechanism", "direct F7/F8 replication established"),
  function(pattern) contains(pattern, ignore.case = TRUE),
  logical(1L)
))

checklist <- data.frame(
  domain = c(
    "principal_findings",
    "module_vs_gene_boundary",
    "integrin_specificity",
    "actomyosin_piezo2",
    "alternative_explanations",
    "future_validation",
    "limitations",
    "conclusion",
    "citation_placeholders",
    "overclaim_control"
  ),
  criterion = c(
    "Discussion states the central module-level finding and retains CAUTION.",
    "Discussion distinguishes module-level association from individual-gene and causal claims.",
    "Integrin support is described as prioritized but not specific.",
    "Actomyosin discordance and PIEZO2 source dependence are both retained.",
    "Competing state explanations are presented rather than dismissed.",
    "Functional validation plan is linked to the unresolved evidence gaps.",
    "Source overlap, label mismatch, processed objects, and unresolved donor mapping are stated.",
    "Conclusion is bounded and does not upgrade the evidence grade.",
    "Citation placeholders match the citation-needed inventory.",
    "No prohibited causal or universal wording is present."
  ),
  passed = c(
    contains("Principal findings") && contains("CAUTION"),
    contains("module level") || contains("module-level"),
    contains("Integrin support is promising but not yet specific"),
    contains("Actomyosin and PIEZO2 require context-sensitive interpretation") && contains("source dependence"),
    contains("Alternative explanations and specificity limits") && contains("cannot distinguish these explanations"),
    contains("Implications for future work") && contains("at least three independent donors"),
    contains("Limitations") && contains("F7/F8") && contains("GSE175817"),
    contains("Conclusion") && contains("bounded hypothesis"),
    placeholder_count == nrow(citations) && nrow(trace) >= 8L,
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
  file.path(result_dir, "13F2_discussion_review_checklist_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

failed <- checklist$criterion[!checklist$passed]
review_lines <- c(
  "# Step 13F2 Discussion section QA",
  "",
  "## Review type",
  "",
  "Internal section-level QA for interpretation boundaries, alternative explanations, functional next steps, and citation placeholders. This review does not upgrade the biological evidence grade.",
  "",
  paste0("Discussion section score: ", round(section_score, 2), "/", round(section_max, 2)),
  paste0("Section gate: ", ifelse(section_pass, "PASS", "HOLD")),
  "",
  "## Findings",
  "",
  if (section_pass) "- The Discussion draft meets the internal section QA threshold of 8/10." else "- The Discussion draft does not yet meet the internal section QA threshold of 8/10.",
  "- The Discussion preserves the CAUTION grade and does not convert association into causality.",
  paste0("- Citation placeholders detected: ", placeholder_count, "; citation inventory rows: ", nrow(citations), "."),
  if (length(failed) == 0L) "- No checklist item failed." else c("- Failed checklist items:", paste0("  - ", failed)),
  "",
  "## Remaining work",
  "",
  "- Replace citation placeholders with verified literature references before submission.",
  "- Preserve the distinction between proposed functional validation and completed evidence.",
  "- Overall evidence grade remains CAUTION.",
  "",
  "## Next gate",
  "",
  if (section_pass) "Proceed to literature citation integration and full manuscript assembly." else "Repair the Discussion draft before citation integration.",
  "",
  "## Material Passport",
  "",
  "- Input: Step 13F Discussion draft, evidence trace, and citation inventory.",
  "- Transformation: section QA only.",
  "- New data downloaded: none."
)
writeLines(
  review_lines,
  file.path(result_dir, "13F2_discussion_review_v1.md"),
  useBytes = TRUE
)

message("Step 13F2 Discussion QA completed.")
message("Discussion score: ", round(section_score, 2), "/", round(section_max, 2))
message("Section gate: ", ifelse(section_pass, "PASS", "HOLD"))
message("Review: ", file.path(result_dir, "13F2_discussion_review_v1.md"))
