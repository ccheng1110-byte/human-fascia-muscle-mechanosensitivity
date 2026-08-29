options(stringsAsFactors = FALSE)

# Step 14A: assemble the complete evidence-grounded manuscript draft.
# No new analysis, hypothesis testing, or data download.

project_dir <- "."
manuscript_root <- file.path(project_dir, "results", "11_manuscript_preparation")
result_dir <- file.path(manuscript_root, "14A_full_manuscript_draft")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

methods_path <- file.path(manuscript_root, "13D_methods_draft", "13D_methods_draft_v1.md")
results_path <- file.path(manuscript_root, "13E_results_draft", "13E_results_draft_v1.md")
discussion_path <- file.path(manuscript_root, "13G_citation_integrated_discussion", "13G_discussion_with_verified_citations_v1.md")
references_path <- file.path(manuscript_root, "13G_citation_integrated_discussion", "13G_verified_reference_list_v1.csv")
evidence_map_path <- file.path(manuscript_root, "13A_manuscript_evidence_pack", "13A_claim_to_evidence_map_v1.csv")

required_inputs <- c(methods_path, results_path, discussion_path, references_path, evidence_map_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 14A input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

read_text <- function(path) {
  paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
}
strip_first_heading <- function(text) {
  sub("^# [^\\n]+\\n*", "", text)
}

methods_text <- strip_first_heading(read_text(methods_path))
results_text <- strip_first_heading(read_text(results_path))
discussion_text <- strip_first_heading(read_text(discussion_path))
references <- read.csv(references_path, check.names = FALSE)
evidence_map <- read.csv(evidence_map_path, check.names = FALSE)

if (grepl("\\[CITATION REQUIRED", discussion_text, fixed = FALSE)) {
  stop("Unresolved citation placeholder remains in the integrated Discussion.")
}

abstract_lines <- c(
  "## Abstract",
  "",
  "### Background",
  "",
  "Fibroblast states within fascia-like connective tissues may integrate extracellular-matrix mechanics through adhesion, cytoskeletal, and transcriptional programs. Whether a reproducible human single-cell mechanotransduction program can be separated from generic fibroblast remodeling remains unresolved.",
  "",
  "### Methods",
  "",
  "We performed a provenance-aware reanalysis of public human single-cell resources comprising discovery analysis in GSE173252, paired sample-level validation in 12 PRJNA607098 sample-state units, and independent external triangulation in five GSE130973 subjects. Frozen candidate genes and predefined modules were evaluated with sample/donor-level evidence boundaries.",
  "",
  "### Results",
  "",
  "The core module and leave-one-sample-out gates passed in the primary validation. Candidate support was 10/15 overall and 5/6 for the integrin/focal-adhesion branch. Competition robustness failed against competing state programs, and PIEZO2 was not supported in the primary source. In GSE130973, integrin support was 5/6, actomyosin candidate support was 1/4, and PIEZO2 was directionally positive, indicating source dependence.",
  "",
  "### Conclusions",
  "",
  "The data support a bounded fibroblast-associated ECM/integrin–cytoskeletal mechanotransduction program, not a universal single-gene mechanism or causal mechanosensitivity claim. Donor-replicated mechanical perturbation with integrin/actomyosin, YAP/TAZ, and state-control readouts is required. Overall evidence grade: CAUTION.",
  "",
  "Keywords: fibroblast; fascia; mechanotransduction; integrin; focal adhesion; actomyosin; YAP/TAZ; PIEZO2; single-cell RNA sequencing"
)

introduction_lines <- c(
  "## Introduction",
  "",
  "Fibroblasts are embedded in extracellular matrices that provide both biochemical ligands and mechanical constraints. Focal adhesions connect integrin–ECM interactions to the actin cytoskeleton, while mechanical cues can be transmitted to transcriptional regulators such as YAP and TAZ (Nardone et al., 2017; Ren et al., 2018). In stromal and fibroblast-related contexts, mechanical stretch, matrix stiffness, and cytoskeletal perturbation can alter YAP/TAZ localization and transcriptional outputs (Hong et al., 2020).",
  "",
  "These observations motivate a multicomponent view of mechanotransduction. RhoA-regulated actomyosin assemblies generate and transmit cellular forces through focal adhesions, providing a functional bridge between cell contractility and matrix interaction (Oakes et al., 2017). Ion channels such as PIEZO2 may contribute in particular cellular or tissue contexts, but the generality of a PIEZO2-centered fibroblast mechanism remains uncertain; recent fibroblast stretch work supports a context-specific role in scleral differentiation rather than a universal connective-tissue rule (Yuan et al., 2026).",
  "",
  "Human single-cell datasets offer an opportunity to examine whether these components co-occur in fibroblast-enriched states. However, transcriptomic association alone cannot establish mechanical causality, and cross-study differences in tissue, state labels, donor structure, and processing can produce apparent replication or discordance. In particular, an ECM/integrin program may be correlated with TGF/fibrosis, inflammation, hypoxia, or proliferation states.",
  "",
  "This is the research gap addressed by the present study: whether a fibroblast-associated mechanotransduction signal remains interpretable after module-level, candidate-level, provenance, and competing-state boundaries are considered.",
  "",
  "Here, we performed a staged human single-cell reanalysis to test whether a fibroblast-associated mechanotransduction program is reproducible at the module level and to define the boundaries of individual candidate genes. We combined discovery analysis, paired sample-level validation, and independent external triangulation. We prespecified that direct F7/F8 equivalence, causal mechanosensitivity, and universal PIEZO2 involvement would not be claimed without stronger provenance or functional evidence."
)

conclusion_lines <- c(
  "## Conclusion",
  "",
  "This study identifies a bounded, fibroblast-associated ECM/integrin–cytoskeletal mechanotransduction program across selected human single-cell sources. Integrin-linked force transmission is the leading functional candidate, whereas actomyosin candidate-gene conservation is incomplete and PIEZO2 is context-dependent. The current evidence remains observational and competition robustness is incomplete. The decisive next step is donor-replicated mechanical perturbation with integrin/actomyosin perturbation, YAP/TAZ localization, functional contractility or traction readouts, and matched ECM/TGF/inflammation and viability controls."
)

availability_lines <- c(
  "## Data and code availability",
  "",
  "The analysis reuses public resources identified by the GEO/NCBI accessions GSE173252, GSE273293, PRJNA607098, GSE130973, and the screened independent sources described in the provenance audit. The project-local analysis scripts, configuration files, audit tables, and result files are stored under `.`. A public repository URL for the complete code release has not yet been assigned and must be added before submission.",
  "",
  "## Ethics statement",
  "",
  "No new human samples or participants were recruited. This study is a reanalysis of publicly available datasets; the original studies' ethics and consent procedures were not re-adjudicated here."
)

reference_lines <- c(
  "## References",
  "",
  paste0(seq_len(nrow(references)), ". ", references$citation, " DOI: ", references$doi)
)

manuscript_lines <- c(
  "# A bounded human single-cell analysis identifies a fibroblast-associated ECM–integrin–cytoskeletal mechanotransduction program",
  "",
  "**Manuscript status:** initial evidence-grounded draft; current evidence grade CAUTION.",
  "",
  abstract_lines,
  "",
  introduction_lines,
  "",
  "## Methods",
  "",
  methods_text,
  "",
  "## Results",
  "",
  results_text,
  "",
  "## Discussion",
  "",
  discussion_text,
  "",
  conclusion_lines,
  "",
  availability_lines,
  "",
  reference_lines
)

writeLines(
  manuscript_lines,
  file.path(result_dir, "14A_full_manuscript_draft_v1.md"),
  useBytes = TRUE
)
writeLines(
  c(abstract_lines, ""),
  file.path(result_dir, "14A_abstract_v1.md"),
  useBytes = TRUE
)

section_manifest <- data.frame(
  section = c("Abstract", "Introduction", "Methods", "Results", "Discussion", "Conclusion", "Data/code availability", "References"),
  source_or_status = c(
    "Generated from frozen Step 10D–13G evidence",
    "Generated from research canon and verified literature context",
    methods_path,
    results_path,
    discussion_path,
    "Generated from Step 12C functional-validation specification",
    "Project-local paths; public repository URL pending",
    references_path
  ),
  evidence_grade = c("CAUTION", "CAUTION", "CAUTION", "CAUTION", "CAUTION", "CAUTION", "N/A", "N/A"),
  status = c("drafted", "drafted", "QA passed", "QA passed", "citation integrated", "drafted", "requires final repository decision", "verified metadata list"),
  stringsAsFactors = FALSE
)
write.csv(
  section_manifest,
  file.path(result_dir, "14A_full_manuscript_section_manifest_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  evidence_map,
  file.path(result_dir, "14A_claim_to_evidence_map_copy_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Step 14A full manuscript draft assembled.")
message("Full manuscript: ", file.path(result_dir, "14A_full_manuscript_draft_v1.md"))
message("Abstract: ", file.path(result_dir, "14A_abstract_v1.md"))
message("Section manifest: ", file.path(result_dir, "14A_full_manuscript_section_manifest_v1.csv"))
