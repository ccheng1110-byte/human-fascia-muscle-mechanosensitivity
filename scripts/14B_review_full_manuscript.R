options(stringsAsFactors = FALSE)

# Step 14B: full manuscript QA and submission-debt register.
# No new analysis or data download.

project_dir <- "."
input_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "14A_full_manuscript_draft")
result_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "14B_full_manuscript_review")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

manuscript_path <- file.path(input_dir, "14A_full_manuscript_draft_v1.md")
abstract_path <- file.path(input_dir, "14A_abstract_v1.md")
manifest_path <- file.path(input_dir, "14A_full_manuscript_section_manifest_v1.csv")
if (!file.exists(manuscript_path) || !file.exists(abstract_path) || !file.exists(manifest_path)) {
  stop("Step 14A manuscript outputs are required.")
}

manuscript <- paste(readLines(manuscript_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
abstract <- paste(readLines(abstract_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
manifest <- read.csv(manifest_path, check.names = FALSE)

contains <- function(pattern, ignore.case = FALSE) {
  if (ignore.case) return(grepl(tolower(pattern), tolower(manuscript), fixed = TRUE))
  grepl(pattern, manuscript, fixed = TRUE)
}

no_placeholders <- !grepl("\\[CITATION REQUIRED|TODO|TBD", manuscript, ignore.case = TRUE)
no_prohibited_overclaim <- !any(vapply(
  c("proves causal", "causal mechanism demonstrated", "universal PIEZO2 mechanism", "direct F7/F8 replication established"),
  function(pattern) contains(pattern, ignore.case = TRUE),
  logical(1L)
))
reference_count <- sum(grepl("^\\d+\\. ", strsplit(manuscript, "\n", fixed = TRUE)[[1L]]))

checklist <- data.frame(
  domain = c(
    "title_and_status",
    "abstract",
    "introduction",
    "methods",
    "results",
    "discussion",
    "conclusion",
    "references",
    "evidence_grade",
    "numeric_consistency",
    "limitations",
    "citation_placeholders",
    "overclaim_control",
    "data_code_statement"
  ),
  criterion = c(
    "Title and draft status are present.",
    "Abstract contains Background, Methods, Results, Conclusions, and CAUTION.",
    "Introduction defines the mechanobiology rationale and research gap.",
    "Methods section is included from the QA-passed draft.",
    "Results section is included with the key 10/15, 5/6, and 1/4 boundaries.",
    "Discussion section is citation-integrated and retains alternative explanations.",
    "Conclusion is bounded and proposes functional validation.",
    "Reference list is present and contains at least five entries.",
    "CAUTION appears in the abstract and main manuscript.",
    "The manuscript contains the primary 12-sample and external five-subject structure.",
    "Limitations include F7/F8, processed objects, source dependence, and causality.",
    "No citation placeholders, TODO, or TBD remain.",
    "No prohibited causal/universal overclaims are present.",
    "Data/code availability and public repository status are explicit."
  ),
  passed = c(
    contains("Manuscript status") && contains("A bounded human single-cell analysis"),
    all(vapply(c("Background", "Methods", "Results", "Conclusions", "CAUTION"), function(x) grepl(x, abstract, fixed = TRUE), logical(1L))),
    contains("## Introduction") && contains("research gap"),
    contains("## Methods") && contains("Study design and analytical scope"),
    contains("## Results") && contains("10/15") && contains("5/6") && contains("1/4"),
    contains("## Discussion") && contains("Limitations") && no_placeholders,
    contains("## Conclusion") && contains("donor-replicated mechanical perturbation"),
    contains("## References") && reference_count >= 5L,
    contains("CAUTION") && grepl("CAUTION", abstract, fixed = TRUE),
    contains("12 PRJNA607098 sample-state units") && contains("five GSE130973 subjects"),
    contains("F7/F8") && contains("processed public objects") && contains("causal"),
    no_placeholders,
    no_prohibited_overclaim,
    contains("Data and code availability") && contains("public repository URL")
  ),
  weight = c(0.5, 1, 1, 1, 1, 1, 0.75, 0.75, 0.75, 0.75, 1, 1, 1, 0.5),
  stringsAsFactors = FALSE
)
checklist$weighted_score <- ifelse(checklist$passed, checklist$weight, 0)
manuscript_score <- sum(checklist$weighted_score)
manuscript_max <- sum(checklist$weight)
manuscript_pass <- manuscript_score >= 10

write.csv(
  checklist,
  file.path(result_dir, "14B_full_manuscript_qa_checklist_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

debt_register <- data.frame(
  debt_id = paste0("D", sprintf("%02d", 1:5)),
  item = c(
    "Public code repository URL",
    "Target-journal formatting",
    "Final figure generation and visual QA",
    "Functional validation experiments",
    "Final author/contribution/ethics metadata"
  ),
  status = c("OPEN", "OPEN", "OPEN", "OPEN", "OPEN"),
  required_before_submission = c(TRUE, TRUE, TRUE, TRUE, TRUE),
  action = c(
    "Create or designate the public repository and replace the pending URL statement.",
    "Adapt headings, reference style, word count, and supplementary format to the selected journal.",
    "Generate publication-quality figures from frozen result tables and run visual inspection.",
    "Perform donor-replicated mechanical challenge and prespecified perturbation readouts if the project claims mechanism.",
    "Complete author list, contributions, funding, conflict, and final data-availability metadata."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  debt_register,
  file.path(result_dir, "14B_submission_debt_register_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

failed <- checklist$criterion[!checklist$passed]
review_lines <- c(
  "# Step 14B full manuscript QA",
  "",
  "## Review type",
  "",
  "Internal full-manuscript QA for section completeness, evidence-boundary preservation, citation placeholders, numerical anchors, and submission debt. This is not an external peer review.",
  "",
  paste0("Manuscript score: ", round(manuscript_score, 2), "/", round(manuscript_max, 2)),
  paste0("Manuscript gate: ", ifelse(manuscript_pass, "PASS_WITH_OPEN_DEBTS", "HOLD")),
  "",
  "## Findings",
  "",
  if (manuscript_pass) "- The assembled manuscript passes the internal completeness/evidence-boundary gate." else "- The assembled manuscript does not yet pass the internal completeness/evidence-boundary gate.",
  "- Biological evidence grade remains CAUTION.",
  "- The manuscript does not claim causal mechanosensitivity, universal PIEZO2, or direct F7/F8 replication.",
  if (length(failed) == 0L) "- No checklist item failed." else c("- Failed checklist items:", paste0("  - ", failed)),
  "",
  "## Open submission debts",
  "",
  "- Public code repository URL is pending.",
  "- Figures have not yet been regenerated and visually checked from frozen outputs.",
  "- Functional validation has not been performed.",
  "- Target-journal formatting and author metadata remain to be finalized.",
  "",
  "## Next gate",
  "",
  if (manuscript_pass) "Proceed to Step 14C figure/table generation and visual QA, followed by journal-specific formatting." else "Repair failed full-manuscript QA items before figure generation.",
  "",
  "## Material Passport",
  "",
  "- Input: Step 14A full manuscript draft and section manifest.",
  "- Transformation: full manuscript QA and debt registration.",
  "- New data downloaded: none.",
  "- Evidence grade: CAUTION."
)
writeLines(
  review_lines,
  file.path(result_dir, "14B_full_manuscript_review_v1.md"),
  useBytes = TRUE
)

message("Step 14B full manuscript QA completed.")
message("Manuscript score: ", round(manuscript_score, 2), "/", round(manuscript_max, 2))
message("Manuscript gate: ", ifelse(manuscript_pass, "PASS_WITH_OPEN_DEBTS", "HOLD"))
message("Review: ", file.path(result_dir, "14B_full_manuscript_review_v1.md"))
