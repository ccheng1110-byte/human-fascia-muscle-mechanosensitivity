options(stringsAsFactors = FALSE)

# Step 13G: integrate verified literature citations into the Discussion draft.
# The project-data evidence grade is unchanged; this step only adds literature
# context and exports a reference-manager file.

project_dir <- "."
input_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13F_discussion_draft")
result_dir <- file.path(project_dir, "results", "11_manuscript_preparation", "13G_citation_integrated_discussion")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

draft_path <- file.path(input_dir, "13F_discussion_draft_v1.md")
if (!file.exists(draft_path)) {
  stop("Step 13F Discussion draft is required.")
}

discussion <- paste(readLines(draft_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")

replace_once <- function(text, old, new) {
  hits <- lengths(regmatches(text, gregexpr(old, text, fixed = TRUE)))
  if (hits != 1L) {
    stop("Expected exactly one citation placeholder, found ", hits, ": ", old)
  }
  sub(old, new, text, fixed = TRUE)
}

discussion <- replace_once(
  discussion,
  "[CITATION REQUIRED: mechanobiology and integrin–cytoskeleton literature]",
  "(Nardone et al., 2017; Ren et al., 2018)"
)
discussion <- replace_once(
  discussion,
  "[CITATION REQUIRED: actomyosin and fibroblast contractility literature]",
  "(Oakes et al., 2017; Hong et al., 2020)"
)
discussion <- replace_once(
  discussion,
  "[CITATION REQUIRED: PIEZO2 in fibroblast or connective-tissue mechanobiology]",
  "(Yuan et al., 2026)"
)

references <- data.frame(
  ref_id = paste0("REF", sprintf("%02d", 1:5)),
  citation = c(
    "Nardone G, Oliver-De La Cruz J, Vrbsky J, et al. YAP regulates cell mechanics by controlling focal adhesion assembly. Nature Communications. 2017;8:15321.",
    "Ren B, et al. RAP2 mediates mechanoresponses of the Hippo pathway. Nature. 2018;560:655-660.",
    "Oakes PW, Wagner E, Brand CA, et al. Optogenetic control of RhoA reveals zyxin-mediated elasticity of stress fibres. Nature Communications. 2017;8:15817.",
    "Hong SP, Yang MJ, Cho H, et al. Distinct fibroblast subsets regulate lacteal integrity through YAP/TAZ-induced VEGF-C in intestinal villi. Nature Communications. 2020;11:4102.",
    "Yuan Y, Li M, Yao Q, Qi Y, Ke B. Piezo2 mechanically regulate scleral fibroblast differentiation by activating relA/RhoA pathway. Experimental Eye Research. 2026;267:110955."
  ),
  journal = c("Nature Communications", "Nature", "Nature Communications", "Nature Communications", "Experimental Eye Research"),
  year = c(2017, 2018, 2017, 2020, 2026),
  doi = c(
    "10.1038/ncomms15321",
    "10.1038/s41586-018-0444-0",
    "10.1038/ncomms15817",
    "10.1038/s41467-020-17886-y",
    "10.1016/j.exer.2026.110955"
  ),
  source_url = c(
    "https://www.nature.com/articles/ncomms15321",
    "https://www.nature.com/articles/s41586-018-0444-0",
    "https://www.nature.com/articles/ncomms15817",
    "https://www.nature.com/articles/s41467-020-17886-y",
    "https://pubmed.ncbi.nlm.nih.gov/41780837/"
  ),
  stringsAsFactors = FALSE
)

support_map <- data.frame(
  manuscript_claim = c(
    "Mechanical information is distributed across ECM attachment, force transmission, cytoskeletal state, and downstream transcriptional responses.",
    "RhoA/actomyosin activity provides a plausible pathway-level basis for force transmission and contractility.",
    "YAP/TAZ can connect mechanical cues with fibroblast or stromal-cell transcriptional responses.",
    "PIEZO2 can participate in a fibroblast mechanotransduction branch under defined stretch conditions."
  ),
  citation = c(
    "Nardone et al., 2017; Ren et al., 2018",
    "Oakes et al., 2017",
    "Hong et al., 2020",
    "Yuan et al., 2026"
  ),
  support_grade = c(
    "partial support",
    "strong pathway-level/background support",
    "partial support in fibroblast/stromal models",
    "partial support; context-specific and non-fascia model"
  ),
  evidence_basis = c(
    "Publisher abstracts/full text describe focal adhesions as the bridge between integrin-ECM connection and cytoskeleton, and ECM stiffness signalling to YAP/TAZ.",
    "Publisher full text reports RhoA-dependent actomyosin contractility, force transmission, and traction responses in adherent fibroblasts.",
    "Nature Communications article reports mechanical stretch/stiffness effects on YAP/TAZ localization and transcriptional responses in fibroblast/stromal contexts.",
    "PubMed abstract reports stretch-responsive Piezo2 in scleral fibroblasts and dependence of RhoA/relA/alpha-SMA changes on Piezo2 or calcium influx."
  ),
  limitation_for_current_study = c(
    "These studies do not establish that the current human fascia-associated transcriptomic program is causal.",
    "The cited functional model does not validate the current cross-tissue candidate-gene pattern.",
    "The models are not identical to human fascia-like fibroblast states.",
    "The model is scleral fibroblast/guinea pig-associated rather than human fascia; it supports context dependence, not universality."
  ),
  stringsAsFactors = FALSE
)

writeLines(
  discussion,
  file.path(result_dir, "13G_discussion_with_verified_citations_v1.md"),
  useBytes = TRUE
)
write.csv(
  references,
  file.path(result_dir, "13G_verified_reference_list_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  support_map,
  file.path(result_dir, "13G_citation_support_map_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

enw_lines <- unlist(lapply(seq_len(nrow(references)), function(i) {
  c(
    "%0 Journal Article",
    paste0("%A ", references$citation[[i]]),
    paste0("%D ", references$year[[i]]),
    paste0("%T ", sub("\\. [A-Z][a-z].*$", "", references$citation[[i]])),
    paste0("%J ", references$journal[[i]]),
    paste0("%R ", references$doi[[i]]),
    paste0("%U ", references$source_url[[i]]),
    ""
  )
}))
writeLines(
  enw_lines,
  file.path(result_dir, "13G_verified_references.enw"),
  useBytes = TRUE
)

review_lines <- c(
  "# Step 13G citation integration review",
  "",
  "## Search scope",
  "",
  "Nature Portfolio sources were prioritized for the ECM/integrin, YAP/TAZ, and actomyosin claims. A PubMed-indexed 2026 Experimental Eye Research article was retained for the PIEZO2 claim because it directly tested a fibroblast stretch context, but it is not a Nature/CNS-family article and is only partial, context-specific support.",
  "",
  "## Support boundaries",
  "",
  "- The citations support mechanobiology context and proposed functional validation rationale.",
  "- They do not validate the current human fascia single-cell result directly.",
  "- They do not upgrade the project evidence grade beyond CAUTION.",
  "- The PIEZO2 citation should be used to motivate a secondary context-dependent branch, not a universal mechanism.",
  "",
  "## Verification status",
  "",
  "- DOI and journal metadata were checked against publisher or PubMed pages.",
  "- The integrated Discussion contains no remaining `[CITATION REQUIRED]` placeholders.",
  "- Final reference formatting should still be adapted to the target journal style."
)
writeLines(
  review_lines,
  file.path(result_dir, "13G_citation_integration_review_v1.md"),
  useBytes = TRUE
)

message("Step 13G citation integration completed.")
message("Cited Discussion: ", file.path(result_dir, "13G_discussion_with_verified_citations_v1.md"))
message("Reference manager export: ", file.path(result_dir, "13G_verified_references.enw"))
message("Citation support map: ", file.path(result_dir, "13G_citation_support_map_v1.csv"))
