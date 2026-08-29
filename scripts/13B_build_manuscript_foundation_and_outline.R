options(stringsAsFactors = FALSE)

# Step 13B: researchwrite foundation and manuscript outline.
# Evidence-before-prose; no new expression analysis or data download.

project_dir <- "."
input_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13A_manuscript_evidence_pack"
)
result_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "13B_researchwrite_foundation"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

claim_map_path <- file.path(input_dir, "13A_claim_to_evidence_map_v1.csv")
inventory_path <- file.path(input_dir, "13A_evidence_file_inventory_v1.csv")
blueprint_path <- file.path(input_dir, "13A_figure_table_blueprint_v1.csv")
required_inputs <- c(claim_map_path, inventory_path, blueprint_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Step 13A input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

claim_map <- read.csv(claim_map_path, check.names = FALSE)
inventory <- read.csv(inventory_path, check.names = FALSE)
blueprint <- read.csv(blueprint_path, check.names = FALSE)

write_md <- function(lines, filename) {
  writeLines(lines, file.path(result_dir, filename), useBytes = TRUE)
}

escape_pipe <- function(x) {
  gsub("\\|", "\\\\|", as.character(x))
}

claim_rows <- vapply(seq_len(nrow(claim_map)), function(i) {
  paste0(
    "| ", escape_pipe(claim_map$manuscript_claim[[i]]),
    " | ", escape_pipe(claim_map$primary_evidence[[i]]),
    " | ", escape_pipe(claim_map$evidence_status[[i]]),
    " | ", ifelse(i %in% c(1L, 5L), "evidence-backed", ifelse(i %in% c(6L), "unsupported", ifelse(i %in% c(4L, 7L), "hypothesis", "plausible-inference"))),
    " | ", escape_pipe(claim_map$permitted_result_language[[i]]),
    " | ", escape_pipe(claim_map$prohibited_result_language[[i]]),
    " |"
  )
}, character(1L))

inventory_rows <- vapply(seq_len(nrow(inventory)), function(i) {
  paste0(
    "| ", escape_pipe(inventory$evidence_id[[i]]),
    " | ", escape_pipe(inventory$evidence_role[[i]]),
    " | `", escape_pipe(inventory$file_path[[i]]), "`",
    " | ", ifelse(isTRUE(inventory$exists[[i]]), "available", "not_found_or_optional"),
    " |"
  )
}, character(1L))

blueprint_rows <- vapply(seq_len(nrow(blueprint)), function(i) {
  paste0(
    "| ", escape_pipe(blueprint$item_id[[i]]),
    " | ", escape_pipe(blueprint$item_type[[i]]),
    " | ", escape_pipe(blueprint$title[[i]]),
    " | ", escape_pipe(blueprint$main_message[[i]]),
    " |"
  )
}, character(1L))

write_md(c(
  "# 00 Scope",
  "",
  "- Text type: evidence-grounded biomedical research manuscript package and outline.",
  "- Target reader: reviewers and researchers in single-cell biology, fibroblast biology, fascia/musculoskeletal biology, and mechanobiology.",
  "- Primary language: English-facing manuscript structure; project communication may remain Chinese.",
  "- Current stage: foundation and outline, not final prose or submission-ready manuscript.",
  "- Evidence boundary: current overall evidence grade is CAUTION.",
  "- Deliverables: research canon, evidence table, argument map, section contracts, style guide, and manuscript outline.",
  "- Excluded from this step: new data download, new hypothesis testing, causal claims, direct F7/F8 equivalence claims, and post hoc threshold changes.",
  "- Archive intent: preserve beside the project results; no external archive or synchronized folder is modified."
), "00_scope.md")

write_md(c(
  "# 01 Research Canon",
  "",
  "## Hard facts",
  "",
  "- GSE173252 supplied the discovery object and the initial mechanotransduction candidate program.",
  "- PRJNA607098 supplied paired sample-level validation across 12 eligible sample-state units.",
  "- PRJNA607098 supported the core module but failed competition robustness and did not support PIEZO2.",
  "- GSE130973 supplied five non-overlapping subjects in a marker-defined fibroblast-state audit.",
  "- GSE130973 supported 5/6 Integrin/focal-adhesion candidates, 1/4 Actomyosin candidates, and PIEZO2 was directionally positive.",
  "- GSE175817 is not currently a donor-level validation source because donor-column mapping was unresolved.",
  "- The atlas F7 state has not been shown to be anatomically or analytically equivalent to the paper F8 state.",
  "- The current evidence is observational, source-dependent, and non-causal.",
  "",
  "## Terminology constraints",
  "",
  "- Program/module: coordinated expression-level association, not a proven causal pathway.",
  "- Candidate: a gene or component prioritized for follow-up, not a confirmed mechanosensor.",
  "- External triangulation: evidence from an independent source with relevant but non-identical labels or design.",
  "- Competition robustness: residual support after comparison with competing state programs.",
  "",
  "## Forbidden claims",
  "",
  "- The study proves causal mechanosensitivity.",
  "- PIEZO2 is a universal marker or required mechanism.",
  "- Integrin enrichment is specific from fibrosis, ECM remodeling, or inflammation.",
  "- Actomyosin candidate genes are conserved across tissues.",
  "- The external atlas is a direct F7/F8 replication.",
  "- Cell-level counts are independent biological replicates."
), "01_research_canon.md")

write_md(c(
  "# 02 Evidence Table",
  "",
  "| Claim | Evidence/source | Strength | Status | Usable section | Risk |",
  "|---|---|---|---|---|---|",
  claim_rows,
  "",
  "## Evidence inventory",
  "",
  "| Evidence ID | Role | File path | Availability |",
  "|---|---|---|---|",
  inventory_rows
), "02_evidence_table.md")

write_md(c(
  "# 03 Argument Map",
  "",
  "## Scientific tension",
  "",
  "A fibroblast-associated mechanotransduction program is detectable across selected human single-cell sources, but the strongest individual genes are not uniformly conserved, competing ECM/TGF/inflammation states remain active, and PIEZO2 changes direction across sources.",
  "",
  "## Central research question",
  "",
  "Do fascia-like or myofibroblast-enriched human fibroblast states show a reproducible multicomponent ECM/integrin-cytoskeletal program, and which components remain credible after source, donor, and competing-state boundaries are considered?",
  "",
  "## Central thesis",
  "",
  "Selected human fibroblast-enriched states show an associated multicomponent ECM/integrin-cytoskeletal mechanotransduction program. This association is strongest at the module level, while individual genes and channel-specific components are context-dependent and require functional validation.",
  "",
  "## Supporting arguments",
  "",
  "1. The discovery analysis identifies a coherent program in the original dataset.",
  "2. Paired sample-level validation supports the core module but exposes competition and PIEZO2 boundaries.",
  "3. A non-overlapping GSE130973 source provides partial integrin triangulation while revealing actomyosin candidate discordance.",
  "4. The resulting hypothesis is more conservative and more testable than a universal single-gene mechanosensor model.",
  "",
  "## Counterarguments and responses",
  "",
  "- Counterargument: the signal may be generic fibrosis/ECM activation. Response: state this as unresolved and prioritize matched competitor controls and functional perturbation.",
  "- Counterargument: candidate discordance invalidates the program. Response: separate module-level evidence from gene-level conservation and test function directly.",
  "- Counterargument: atlas label mismatch prevents replication. Response: describe the atlas as external triangulation, not direct replication.",
  "",
  "## Final move",
  "",
  "The manuscript should conclude with the minimum functional experiments required to distinguish a mechanotransduction mechanism from a correlated fibroblast state."
), "03_argument_map.md")

write_md(c(
  "# 04 Section Contracts",
  "",
  "## Abstract",
  "- Purpose: state the question, design, main bounded result, and evidence grade.",
  "- Allowed claims: module-level association, partial external triangulation, explicit limitations.",
  "- Forbidden claims: causal mechanism, universal PIEZO2, direct F7/F8 replication.",
  "- Required validation: every numerical statement maps to the frozen evidence table.",
  "",
  "## Introduction",
  "- Purpose: motivate fascia/fibroblast mechanobiology and define the unresolved specificity problem.",
  "- Inputs: research canon and literature to be added later.",
  "- Allowed claims: biological rationale and testable gap.",
  "- Forbidden claims: present study results before the Results section or imply causal proof.",
  "",
  "## Methods",
  "- Purpose: make provenance, sample structure, frozen panels, audit gates, and reproducibility explicit.",
  "- Inputs: evidence inventory, scripts, data paths, and audit decisions.",
  "- Allowed claims: what was downloaded, analyzed, excluded, and why.",
  "- Forbidden claims: retrospectively changing thresholds or treating cells as donors.",
  "",
  "## Results 1 — Discovery program",
  "- Purpose: establish the original program and candidate modules.",
  "- Allowed claims: discovery-level association and module composition.",
  "- Forbidden claims: causal interpretation or universal gene-level conservation.",
  "",
  "## Results 2 — Primary paired validation",
  "- Purpose: present PRJNA607098 sample-level validation and competitor robustness.",
  "- Allowed claims: core module support, failed competition gate, PIEZO2 negative result.",
  "- Forbidden claims: clean specificity or universal PIEZO2 mechanism.",
  "",
  "## Results 3 — Independent triangulation",
  "- Purpose: present GSE130973 as non-overlapping, marker-defined external evidence.",
  "- Allowed claims: 5/6 integrin support, 1/4 actomyosin support, PIEZO2 source dependence.",
  "- Forbidden claims: direct F7/F8 replication or strict anatomical equivalence.",
  "",
  "## Discussion",
  "- Purpose: interpret the module-level finding while preserving uncertainty.",
  "- Allowed claims: revised hypothesis, context dependence, functional next steps.",
  "- Forbidden claims: upgrade beyond CAUTION without causal functional evidence.",
  "",
  "## Limitations",
  "- Purpose: make source overlap, label mismatch, donor limitations, competition failure, and observational design explicit.",
  "- Required validation: every limitation must point to a file or audit result.",
  "",
  "## Conclusion",
  "- Purpose: one bounded conclusion and one actionable next step.",
  "- Allowed claims: module-level association and need for donor-replicated perturbation.",
  "- Forbidden claims: mechanosensitivity demonstrated, PIEZO2 required, direct replication established."
), "04_section_contracts.md")

write_md(c(
  "# 05 Style Guide",
  "",
  "- Stance: conservative, evidence-led, mechanistic but non-causal.",
  "- Prefer: associated with, supported at the module level, directionally consistent, context-dependent, triangulates, remains unresolved.",
  "- Avoid: proves, demonstrates a mechanism, universal, specific signature, definitive, confirms F7/F8 replication.",
  "- Separate module-level findings from candidate-gene findings in every Results subsection.",
  "- Report donor/sample units before cell counts; avoid cell-level pseudoreplication.",
  "- Keep negative and discordant results visible rather than smoothing them into a single narrative.",
  "- Do not turn a failed gate into a positive claim by changing the threshold after seeing the result.",
  "- Every quantitative statement must map to a stored result file.",
  "- The manuscript should state the current evidence grade as CAUTION until functional validation is completed."
), "05_style_guide.md")

write_md(c(
  "# Step 13B Manuscript Outline",
  "",
  "## Provisional title",
  "",
  "A bounded human single-cell analysis identifies a fibroblast-associated ECM–integrin–cytoskeletal mechanotransduction program",
  "",
  "## Abstract structure",
  "",
  "1. Background: mechanical regulation of fascia-like fibroblast states remains incompletely resolved.",
  "2. Approach: discovery analysis, paired sample-level validation, and independent external triangulation.",
  "3. Main result: module-level program support with partial integrin evidence and source-dependent PIEZO2/actomyosin components.",
  "4. Boundary: competition robustness, direct F7/F8 equivalence, and causality remain unresolved.",
  "5. Implication: functional perturbation is required.",
  "",
  "## Main text",
  "",
  "### Introduction",
  "- Fascia/fibroblast mechanical biology.",
  "- ECM–integrin–cytoskeletal mechanotransduction rationale.",
  "- Why single-gene mechanosensor models may be insufficient.",
  "- Study question and prespecified evidence boundaries.",
  "",
  "### Results",
  "- 1. Dataset provenance, labels, and donor/sample structure.",
  "- 2. Discovery of the multicomponent program in GSE173252.",
  "- 3. Paired validation in PRJNA607098 and the competition-robustness boundary.",
  "- 4. Independent GSE130973 triangulation.",
  "- 5. Revised hypothesis and functional-validation priorities.",
  "",
  "### Discussion",
  "- Module-level reproducibility is stronger than gene-level conservation.",
  "- Integrin branch is prioritized but not yet specific.",
  "- Actomyosin candidate discordance supports functional rather than transcript-only testing.",
  "- PIEZO2 is context-dependent.",
  "- Limits of F7/F8 equivalence and observational inference.",
  "",
  "### Conclusion",
  "- The data support a bounded association hypothesis and define the functional experiments needed to test mechanism.",
  "",
  "## Figure and table alignment",
  "",
  blueprint_rows
), "13B_manuscript_outline_v1.md")

state_lines <- c(
  "{",
  '  "project": "human_fascia_muscle_mechanosensitivity",',
  '  "mode": "compose",',
  '  "text_type": "evidence_grounded_research_manuscript_package",',
  '  "language": "en",',
  '  "target_reader": "biomedical_mechanobiology_reviewers",',
  '  "current_round": 1,',
  '  "scores": [],',
  '  "technical_debts": ["literature citations not yet attached", "functional validation not yet performed"],',
  '  "status": "foundation"',
  "}"
)
writeLines(state_lines, file.path(result_dir, "state.json"), useBytes = TRUE)

message("Step 13B manuscript foundation and outline completed.")
message("Foundation directory: ", result_dir)
message("Manuscript outline: ", file.path(result_dir, "13B_manuscript_outline_v1.md"))
message("State: foundation")
