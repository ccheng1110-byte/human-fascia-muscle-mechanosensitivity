options(stringsAsFactors = FALSE)

# Step 18A: update the manuscript after second-round computational strengthening.
# The original 14A manuscript is preserved. This script creates a revised v2
# manuscript and an auditable claim/evidence update record.

project_dir <- "."
manuscript_path <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "14A_full_manuscript_draft", "14A_full_manuscript_draft_v1.md"
)
evidence_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17F_evidence_resynthesis"
)
output_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "18A_second_round_revised_manuscript"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  manuscript_path,
  file.path(evidence_dir, "Step17F_second_round_evidence_resynthesis_decision_v1.md"),
  file.path(evidence_dir, "Step17F_revised_claim_contract_v1.csv"),
  file.path(evidence_dir, "Step17F_manuscript_update_actions_v1.md")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 18A input(s): ", paste(missing_inputs, collapse = "; "))
}

read_text <- function(path) paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
safe_write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
write_text <- function(x, path) writeLines(x, path, useBytes = TRUE)

insert_before_heading <- function(text, heading, block, marker) {
  if (grepl(marker, text, fixed = TRUE)) return(text)
  heading_pos <- regexpr(paste0("(?m)^", heading, "$"), text, perl = TRUE)
  if (heading_pos[[1L]] < 0L) stop("Could not find manuscript heading: ", heading)
  prefix <- substr(text, 1L, heading_pos[[1L]] - 1L)
  suffix <- substr(text, heading_pos[[1L]], nchar(text))
  paste0(prefix, block, "\n\n", suffix)
}

replace_once <- function(text, old, new, label) {
  if (!grepl(old, text, fixed = TRUE)) {
    warning("Step 18A could not find text for replacement: ", label)
    return(text)
  }
  sub(old, new, text, fixed = TRUE)
}

manuscript <- read_text(manuscript_path)

# Fix known structural artefacts from the 14A draft without changing evidence.
manuscript <- replace_once(manuscript, "n and analytical scope", "Analytical scope", "Methods heading")
manuscript <- replace_once(manuscript, "nance and validation structure", "Provenance and validation structure", "Results heading")
manuscript <- replace_once(manuscript, "## Discussion\n\nn\n\n## Principal findings", "## Discussion\n\n## Principal findings", "Discussion artefact")

# Update the abstract so that the second-round sources and their boundaries are visible.
manuscript <- replace_once(
  manuscript,
  "We performed a provenance-aware reanalysis of public human single-cell resources comprising discovery analysis in GSE173252, paired sample-level validation in 12 PRJNA607098 sample-state units, and independent external triangulation in five GSE130973 subjects. Frozen candidate genes and predefined modules were evaluated with sample/donor-level evidence boundaries.",
  "We performed a provenance-aware reanalysis of public human transcriptomic resources comprising discovery analysis in GSE173252, paired sample-level validation in 12 PRJNA607098 sample-state units, independent external triangulation in five GSE130973 subjects, regulatory-axis cross-validation in GSE338388, and cross-tissue stiffness analyses in GSE123100 and GSE276045. Frozen candidate genes and predefined modules were evaluated with sample/donor-level evidence boundaries.",
  "Abstract Methods"
)
manuscript <- replace_once(
  manuscript,
  "The core module and leave-one-sample-out gates passed in the primary validation. Candidate support was 10/15 overall and 5/6 for the integrin/focal-adhesion branch. Competition robustness failed against competing state programs, and PIEZO2 was not supported in the primary source. In GSE130973, integrin support was 5/6, actomyosin candidate support was 1/4, and PIEZO2 was directionally positive, indicating source dependence.",
  "The core module and leave-one-sample-out gates passed in the primary validation. Candidate support was 10/15 overall and 5/6 for the integrin/focal-adhesion branch. Competition robustness failed against competing state programs, and PIEZO2 was not supported in the primary source. In GSE130973, integrin support was 5/6, actomyosin candidate support was 1/4, and PIEZO2 was directionally positive, indicating source dependence. In the second-round analyses, GSE338388 supported a TGF-beta exposure × TEAD inhibition regulatory-axis interpretation with 15/15 frozen-candidate coverage, while GSE123100 provided a descriptive 14/15-candidate stiffness-form analysis and GSE276045 showed actomyosin/Rho stiffness-direction support in both WT and hTERT WI-38 cells. ECM and cell-cycle competitors remained active.",
  "Abstract Results"
)
manuscript <- replace_once(
  manuscript,
  "The data support a bounded fibroblast-associated ECM/integrin–cytoskeletal mechanotransduction program, not a universal single-gene mechanism or causal mechanosensitivity claim. Donor-replicated mechanical perturbation with integrin/actomyosin, YAP/TAZ, and state-control readouts is required. Overall evidence grade: CAUTION.",
  "The data support an evidence-bounded fibroblast-associated ECM/integrin–cytoskeletal mechanotransduction hypothesis with computational triangulation, not a universal single-gene mechanism, cell-intrinsic mechanism, or causal mechanosensitivity claim. Donor-replicated mechanical perturbation with integrin/actomyosin, YAP/TAZ, functional, and state-control readouts is required. Overall evidence grade: CAUTION.",
  "Abstract Conclusions"
)

methods_block <- paste0(
  "<!-- STEP17F_METHODS_START -->\n",
  "## Second-round computational strengthening\n\n",
  "The second-round analysis was prespecified to narrow three uncertainties without expanding the frozen candidate panel or relaxing thresholds. For S1, PRJNA607098 and GSE130973 were analyzed as separate sources at the co-expression level; cell-level double-positive proportions were not interpreted as proof of a cell-intrinsic mechanism. For S2, GSE338388 was analyzed as a TGF-beta exposure × TEAD inhibition design, with the two axes and their module responses treated as regulatory-axis cross-validation rather than mechanical causality. For S3, GSE123100 was restricted to a descriptive cultured HTM stiffness-form analysis, whereas GSE276045 was modeled as a WI-38 stiffness slope with WT/hTERT cell-model interaction and timepoint adjustment. The S3 sources were not pooled as biological replicates.\n\n",
  "The second-round outputs were integrated using the frozen 15-gene candidate panel and the registered mechanotransduction and competitor modules. Candidate coverage, module coverage, effect direction, sample/donor structure, and unresolved confounding were retained in machine-readable audit files. Because proliferation was not measured or resolved in GSE276045, a proliferation-independent stiffness effect was not estimated. No computational result was allowed to remove the prespecified CAUTION evidence grade or replace donor-replicated functional validation.\n",
  "<!-- STEP17F_METHODS_END -->"
)
manuscript <- insert_before_heading(
  manuscript, "## Evidence synthesis and interpretation rules", methods_block,
  "<!-- STEP17F_METHODS_START -->"
)

results_block <- paste0(
  "<!-- STEP17F_RESULTS_START -->\n",
  "## Second-round computational triangulation\n\n",
  "The second-round analyses were designed to distinguish co-expression-level evidence, regulatory-axis cross-validation, and cross-tissue stiffness-form evidence. In S1, PRJNA607098 descriptively supported 37/54 tested gene pairs, whereas GSE130973 supported 0/54. This source dependence is consistent with a co-expression-level signal but does not establish a cell-intrinsic mechanism.\n\n",
  "In S2, GSE338388 contained 15/15 frozen candidates and supported the prespecified TGF-beta exposure × TEAD inhibition regulatory-axis interpretation. The dataset contains no mechanical loading or stiffness perturbation; accordingly, these results provide regulatory-axis cross-validation rather than evidence of mechanical causality or YAP/TAZ-independent mechanical driving.\n\n",
  "In S3, GSE123100 contained 14/15 frozen candidates and was retained as a descriptive HTM dose-response-form analysis because ITGA5 was absent from the processed matrix and the design was not treated as a balanced mechanistic experiment. GSE276045 contained 15/15 frozen candidates and supported stiffness-slope modeling with cell-model interaction and timepoint adjustment. The actomyosin/Rho module showed stiffness-direction support in both WT and hTERT WI-38 cells (camera FDR 0.0009925 and 0.0006892, respectively). However, cell-cycle and ECM competitor modules also responded, so proliferation-independent and mechanistically specific interpretation was not established.\n\n",
  "Together, the second-round results strengthened cross-tissue pathway plausibility while preserving the central limitations: no fascia-specific causal test, no proof of cell-intrinsic mechanism, no direct F7/F8 replication, and no basis for removing the CAUTION evidence grade.\n",
  "<!-- STEP17F_RESULTS_END -->"
)
manuscript <- insert_before_heading(
  manuscript, "## Integrated evidence boundary", results_block,
  "<!-- STEP17F_RESULTS_START -->"
)

discussion_block <- paste0(
  "<!-- STEP17F_DISCUSSION_START -->\n",
  "## How the second-round evidence changes the interpretation\n\n",
  "The second-round analyses do not convert the study from an observational association into a causal mechanosensitivity demonstration. They sharpen the evidence hierarchy. S1 is best described as source-dependent co-expression-level evidence; S2 provides a TGF-beta/TEAD regulatory-axis cross-validation; and S3 adds cross-tissue stiffness-form evidence. This triangulation makes actomyosin/Rho a rational functional priority, but it does not establish that the response is specific to mechanical input or independent of fibroblast state.\n\n",
  "The WI-38 stiffness analysis is informative because actomyosin/Rho direction was supported in both WT and hTERT cell-model slopes after timepoint adjustment. The same analysis also detected cell-cycle and ECM responses, however, and proliferation was not directly measured or resolved. The result therefore supports pathway plausibility with an explicit competition boundary. The HTM analysis provides a complementary dose-response form but is descriptive and has incomplete frozen-candidate coverage. Neither source is a direct fascia replication.\n\n",
  "The resulting claim is narrower but more defensible: selected human fibroblast-enriched states show a multicomponent ECM/integrin–cytoskeletal program that is compatible with, and computationally triangulated toward, a mechanistic hypothesis. The decisive test remains a donor-replicated functional experiment using controlled mechanical perturbation, prespecified integrin and actomyosin perturbations, an orthogonal functional readout, and matched ECM/TGF, inflammation, proliferation, and viability controls.\n",
  "<!-- STEP17F_DISCUSSION_END -->"
)
manuscript <- insert_before_heading(
  manuscript, "## A multicomponent program is more defensible than a single-gene model", discussion_block,
  "<!-- STEP17F_DISCUSSION_START -->"
)

# Remove the duplicated early conclusion and replace the remaining conclusion with the locked claim.
duplicate_conclusion <- "## Conclusion\n\nThe study supports a bounded hypothesis that fibroblast-enriched human states are associated with a multicomponent ECM/integrin–cytoskeletal mechanotransduction program. Integrin-linked force transmission is the leading mechanistic candidate, actomyosin activity remains a pathway-level hypothesis with candidate-gene discordance, and PIEZO2 is context-dependent. The decisive next step is donor-replicated functional mechanical perturbation with state controls and orthogonal mechanotransduction readouts.\n\n"
manuscript <- replace_once(manuscript, duplicate_conclusion, "", "duplicate Conclusion")
old_final_conclusion <- "This study identifies a bounded, fibroblast-associated ECM/integrin–cytoskeletal mechanotransduction program across selected human single-cell sources. Integrin-linked force transmission is the leading functional candidate, whereas actomyosin candidate-gene conservation is incomplete and PIEZO2 is context-dependent. The current evidence remains observational and competition robustness is incomplete. The decisive next step is donor-replicated mechanical perturbation with integrin/actomyosin perturbation, YAP/TAZ localization, functional contractility or traction readouts, and matched ECM/TGF/inflammation and viability controls."
new_final_conclusion <- "This study supports an evidence-bounded, fibroblast-associated ECM/integrin–cytoskeletal mechanotransduction hypothesis across selected human transcriptomic sources. Actomyosin/Rho is the leading pathway-level functional priority after cross-tissue stiffness triangulation, but ECM and cell-cycle competition remain active; integrin specificity is unresolved and PIEZO2 is context-dependent. The evidence remains observational and is graded CAUTION. The decisive next step is donor-replicated mechanical perturbation with prespecified integrin/actomyosin perturbation, YAP/TAZ localization, functional contractility or traction readouts, and matched ECM/TGF/inflammation, proliferation, and viability controls."
manuscript <- replace_once(manuscript, old_final_conclusion, new_final_conclusion, "final Conclusion")

output_manuscript <- file.path(output_dir, "18A_second_round_revised_manuscript_v2.md")
write_text(manuscript, output_manuscript)

claim_contract <- utils::read.csv(
  file.path(evidence_dir, "Step17F_revised_claim_contract_v1.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
claim_contract$manuscript_update <- c(
  "Abstract; Results; Discussion; Conclusion",
  "Methods; Results; Discussion; Limitations",
  "Abstract; Methods; Results; Discussion",
  "Abstract; Results; Discussion; Future work",
  "Results; Discussion; Limitations",
  "Abstract; Results; Discussion",
  "Abstract; Discussion; Conclusion"
)
safe_write_csv(claim_contract, file.path(output_dir, "18A_claim_evidence_trace_v2.csv"))

allocation <- data.frame(
  result = c(
    "S1 source-dependent co-expression level",
    "S2 TGF-beta exposure x TEAD inhibition",
    "S3 GSE123100 HTM stiffness form",
    "S3 GSE276045 WI-38 stiffness slopes",
    "S3 cell-cycle and ECM competition",
    "Overall CAUTION boundary"
  ),
  evidence_class = c("qualification", "necessary_support", "qualification", "core_discovery", "qualification", "edge_case"),
  main_text_destination = c("Results and Discussion", "Methods and Results", "Results and Discussion", "Results and Discussion", "Results and Discussion", "Abstract, Discussion and Conclusion"),
  full_numeric_destination = c("Step17F S1 outputs", "GSE338388 S2 output tables", "GSE123100 S3 output tables", "GSE276045 S3 output tables", "Step17F S3 cross-source audit", "Step17F decision file"),
  stringsAsFactors = FALSE
)
safe_write_csv(allocation, file.path(output_dir, "18A_result_allocation_v1.csv"))

before_words <- length(strsplit(gsub("\n", " ", read_text(manuscript_path)), "[[:space:]]+", perl = TRUE)[[1L]])
after_words <- length(strsplit(gsub("\n", " ", manuscript), "[[:space:]]+", perl = TRUE)[[1L]])
revision_log <- c(
  "# Step 18A manuscript revision log",
  "",
  "- Original manuscript preserved: `14A_full_manuscript_draft_v1.md`.",
  "- Revised manuscript: `18A_second_round_revised_manuscript_v2.md`.",
  paste0("- Approximate word count before: ", before_words, "; after: ", after_words, "."),
  "- Added second-round Methods, Results and Discussion blocks using the locked Step 17F evidence boundaries.",
  "- Replaced S1 cell-intrinsic risk language with co-expression-level language.",
  "- Reframed GSE338388 as TGF-beta/TEAD regulatory-axis cross-validation.",
  "- Reframed GSE123100 and GSE276045 as cross-tissue stiffness-form evidence.",
  "- Added the GSE276045 actomyosin result together with cell-cycle/ECM competition.",
  "- Retained the CAUTION evidence grade and the donor-replicated functional upgrade requirement.",
  "- Removed a duplicated Conclusion heading and corrected three inherited draft artefacts.",
  "",
  "## Human verification required",
  "",
  "The author should verify every numerical statement against the cited CSV outputs before submission and add the final public code repository URL to Data and code availability."
)
write_text(revision_log, file.path(output_dir, "18A_manuscript_revision_log_v1.md"))

message("Step 18A second-round manuscript update completed.")
message("Revised manuscript: ", output_manuscript)
message("Revision audit: ", output_dir)
