options(stringsAsFactors = FALSE)

# Step 08B3: final label-column-corrected validation in PRJNA607098.
# Run this script in a freshly restarted R session.
# It will not download the 25.36 GiB H5AD. Required Zarr chunks are checked by
# HTTP HEAD and selectively streamed by the user-run R job. The 50 MiB limit
# applies only to downloads performed directly by Codex, not downloads via R.

project_dir <- "."

# Windows path-safety fix for reticulate/uv.
# The Windows user-profile path contains non-ASCII characters, which can be
# corrupted by a mismatched locale while reticulate is locating its managed
# Python environment. Redirect every relevant cache/install location to an
# ASCII-only directory before the reticulate namespace is loaded.
python_runtime_root <- file.path(project_dir, ".runtime", "python_step08b")
r_user_cache_dir <- file.path(python_runtime_root, "r_user_cache")
uv_cache_dir <- file.path(python_runtime_root, "uv_cache")
uv_python_install_dir <- file.path(python_runtime_root, "uv_python")

runtime_dirs <- c(
  python_runtime_root,
  r_user_cache_dir,
  uv_cache_dir,
  uv_python_install_dir
)
for (runtime_dir in runtime_dirs) {
  dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
}

Sys.setenv(
  R_USER_CACHE_DIR = r_user_cache_dir,
  XDG_CACHE_HOME = r_user_cache_dir,
  UV_CACHE_DIR = uv_cache_dir,
  UV_PYTHON_INSTALL_DIR = uv_python_install_dir,
  RETICULATE_PYTHON = "managed"
)

reticulate_cache_dir <- tools::R_user_dir("reticulate", "cache")
if (!startsWith(
  tolower(normalizePath(reticulate_cache_dir, winslash = "/", mustWork = FALSE)),
  tolower(normalizePath(python_runtime_root, winslash = "/", mustWork = TRUE))
)) {
  stop(
    "reticulate cache redirection failed. Expected a path under: ",
    python_runtime_root, "; resolved: ", reticulate_cache_dir,
    ". Restart R and source this script before loading reticulate or Seurat."
  )
}

step07_path <- file.path(
  project_dir, "results", "05_specificity_robustness", "GSE173252",
  "GSE173252_DD_candidate_specificity_summary_v1.csv"
)
python_helper <- file.path(
  project_dir, "scripts",
  "08B_extract_PRJNA607098_gene_panel_from_remote_zarr.py"
)
result_dir <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025", "08B3_F7_PRJNA607098"
)
manual_download_dir <- file.path(
  project_dir, "data", "external_validation", "skin_fibroblast_atlas_2025",
  "manual_downloads", "zarr_chunks"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manual_download_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(step07_path, python_helper)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop(
    "R package 'reticulate' is required. Install it once with: ",
    "install.packages('reticulate')"
  )
}
if (utils::packageVersion("reticulate") < "1.41.0") {
  stop(
    "reticulate >= 1.41.0 is required for py_require(). Update with: ",
    "install.packages('reticulate')"
  )
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

as_logical_safe <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

# Freeze the same 15 A/B robust candidates used in Step 08A.
step07 <- read.csv(step07_path, check.names = FALSE)
shortlist <- step07[
  grepl("^[AB]_", step07$step06_candidate_tier) &
    step07$specificity_class == "robust_Myofib_over_all_three_identities",
  , drop = FALSE
]
shortlist <- shortlist[!duplicated(shortlist$gene), , drop = FALSE]
if (nrow(shortlist) != 15L) {
  stop(
    "Expected 15 frozen Step-07 candidates but found ", nrow(shortlist),
    ". Review Step-07 outputs before continuing."
  )
}

mechanosensor_context <- c("PIEZO1", "TRPV4", "PKD2")
panel_genes <- unique(c(shortlist$gene, mechanosensor_context))

# Current reticulate guidance recommends a managed ephemeral environment.
# Requirements are declared before Python initializes. RETICULATE_PYTHON was
# set above, before loading reticulate, so the managed environment is built in
# the project's ASCII-only runtime directory.
reticulate::py_require(
  packages = c(
    "zarr>=2.18,<3",
    "fsspec>=2024.6",
    "aiohttp>=3.9",
    "numpy>=1.26,<3"
  ),
  python_version = ">=3.10,<3.13",
  exclude_newer = "2026-08-23"
)

python_environment <- new.env(parent = emptyenv())
reticulate::source_python(
  python_helper,
  envir = python_environment,
  convert = TRUE
)

message("Opening the public atlas Zarr for the final F7 source-specific audit.")
message("Only gene-panel Zarr chunks will be streamed; the full 25.36 GiB H5AD is not downloaded.")

python_result <- python_environment$run_step08b(
  panel_genes = as.list(panel_genes),
  output_dir = result_dir,
  manual_download_dir = manual_download_dir,
  state_label_field = "celltype",
  target_state = "F7: Fascia-like myofibroblast",
  output_version = "v3",
  audit_chunk_sizes = FALSE
)

contrast_path <- file.path(
  result_dir,
  "skin_fibroblast_atlas_source_specific_F7_gene_contrasts_v3.csv"
)
chunk_manifest_path <- file.path(
  result_dir,
  "skin_fibroblast_atlas_gene_panel_chunk_manifest_v3.csv"
)
if (!file.exists(contrast_path)) {
  stop("Python extraction completed without the expected contrast table.")
}

contrasts <- read.csv(contrast_path, check.names = FALSE)
independent <- contrasts[
  contrasts$GSE == "PRJNA607098" & contrasts$gene %in% shortlist$gene,
  , drop = FALSE
]
if (nrow(independent) != nrow(shortlist)) {
  stop(
    "PRJNA607098 contrast table does not contain all frozen candidates. ",
    "Found ", nrow(independent), " of ", nrow(shortlist), "."
  )
}
independent$source_specific_descriptive_support <- as_logical_safe(
  independent$source_specific_descriptive_support
)
independent$F7_highest_mean_across_states <- as_logical_safe(
  independent$F7_highest_mean_across_states
)

candidate_results <- merge(
  shortlist[, c(
    "module", "gene", "step06_candidate_tier",
    "myofib_percent_cells_expressing",
    "minimum_difference_over_all_donors_and_identities"
  ), drop = FALSE],
  independent,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
candidate_results <- candidate_results[
  order(
    -candidate_results$source_specific_descriptive_support,
    -candidate_results$F7_minus_median_other_state_mean,
    candidate_results$gene
  ),
  , drop = FALSE
]
rownames(candidate_results) <- NULL

candidate_output <- file.path(
  result_dir,
  "PRJNA607098_F7_frozen_candidate_validation_v3.csv"
)
safe_write_csv(candidate_results, candidate_output)

module_rows <- lapply(split(candidate_results, candidate_results$module), function(x) {
  data.frame(
    module = x$module[[1L]],
    candidate_genes = nrow(x),
    supported_genes = sum(
      x$source_specific_descriptive_support %in% TRUE,
      na.rm = TRUE
    ),
    support_fraction = mean(
      x$source_specific_descriptive_support %in% TRUE,
      na.rm = TRUE
    ),
    median_F7_minus_median_other_state_mean = median(
      x$F7_minus_median_other_state_mean,
      na.rm = TRUE
    ),
    supported_gene_names = paste(
      x$gene[x$source_specific_descriptive_support %in% TRUE],
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
})
module_summary <- do.call(rbind, module_rows)
module_summary <- module_summary[
  order(-module_summary$support_fraction, module_summary$module),
  , drop = FALSE
]
rownames(module_summary) <- NULL

module_output <- file.path(
  result_dir,
  "PRJNA607098_F7_module_validation_v3.csv"
)
safe_write_csv(module_summary, module_output)

# Pre-specified 08B -> 08C gate, frozen before looking at the 08B expression data.
n_supported <- sum(
  candidate_results$source_specific_descriptive_support %in% TRUE,
  na.rm = TRUE
)
integrin_rows <- candidate_results[
  candidate_results$module == "integrin_focal_adhesion",
  , drop = FALSE
]
n_integrin_supported <- sum(
  integrin_rows$source_specific_descriptive_support %in% TRUE,
  na.rm = TRUE
)
piezo2_supported <- candidate_results$source_specific_descriptive_support[
  candidate_results$gene == "PIEZO2"
]
piezo2_supported <- length(piezo2_supported) == 1L && isTRUE(piezo2_supported)

gate_pass <- n_supported >= 8L &&
  n_integrin_supported >= 4L &&
  piezo2_supported

chunk_manifest <- read.csv(chunk_manifest_path, check.names = FALSE)
chunk_sizes <- suppressWarnings(as.numeric(chunk_manifest$size_bytes))
total_stream_bytes <- sum(chunk_sizes, na.rm = TRUE)
n_unknown_chunk_sizes <- sum(is.na(chunk_sizes))
stream_size_line <- if (n_unknown_chunk_sizes == 0L) {
  paste0(
    "- Total streamed expression-chunk volume: ",
    round(total_stream_bytes / 1024^2, 2), " MiB."
  )
} else {
  paste0(
    "- Known expression-chunk volume: ",
    round(total_stream_bytes / 1024^2, 2), " MiB; ",
    n_unknown_chunk_sizes,
    " chunk size(s) unavailable after non-blocking network audit."
  )
}

report_lines <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate",
  "- Origin Date: 2026-08-23",
  "- Verification Status: ANALYZED",
  "- Version Label: PRJNA607098_F7_celltype_column_v3",
  "",
  "## Step 08B3 final PRJNA607098 F7 source-specific validation",
  "",
  "### Label-column correction audit",
  "",
  "- Step 08B v1 searched for F7 in `celltype_skinspecific_nomenclature`, where PRJNA607098 has no F7 label; v1 is invalid as a biological negative result.",
  "- Step 08B2 tested the auxiliary F6/Myofibroblast population (1,714 cells) in that column and is retained as a sensitivity analysis.",
  "- Step 08B3 uses the intended atlas `celltype` column and the frozen target `F7: Fascia-like myofibroblast`.",
  "- The frozen 15 genes, contrast rule and gate are unchanged.",
  "",
  "### Pre-specified gate",
  "",
  "Proceed to donor-level Step 08C only if all conditions are met:",
  "",
  "- At least 8 of 15 frozen candidates show descriptive support.",
  "- At least 4 of 6 integrin/focal-adhesion candidates show support.",
  "- PIEZO2 shows support.",
  "",
  "### Results",
  "",
  paste0("- Frozen candidates supported: ", n_supported, "/15."),
  paste0("- Integrin/focal-adhesion candidates supported: ", n_integrin_supported, "/6."),
  paste0("- PIEZO2 supported: ", piezo2_supported, "."),
  paste0("- Step 08C gate passed: ", gate_pass, "."),
  stream_size_line,
  "- Required files were streamed by the user-run R workflow; no 50 MiB per-file limit was applied.",
  "",
  "### Evidence decision",
  "",
  "- Evidence grade remains CAUTION even if the gate passes.",
  "- This step isolates PRJNA607098 from GSE173252, improving dataset independence.",
  "- The public atlas lacks sample/donor labels, so cells are not treated as replicates and no p values are calculated.",
  "- A passing gate supports investment in Step 08C donor/sample-level reconstruction.",
  "- Mechanosensor activity, protein abundance and causality remain untested.",
  "",
  "### Statistical fallacy scan (11/11 checked)",
  "",
  "- Simpson's paradox: not testable without donor strata; retained as a caution.",
  "- Ecological fallacy: avoided; no cell-state result is generalized to individual patients.",
  "- Berkson's paradox: possible surgical/atlas selection bias; caution.",
  "- Collider bias: no covariate-adjusted causal model was fitted.",
  "- Base-rate neglect: cell counts and expression prevalence are reported.",
  "- Regression to the mean: not applicable to this cross-sectional screen.",
  "- Survivorship bias: atlas QC selection may exclude low-quality cells; caution.",
  "- Look-elsewhere effect: reduced by using a frozen 15-gene candidate panel.",
  "- Garden of forking paths: gate and contrasts were fixed before 08B results.",
  "- Correlation is not causation: causal claims are prohibited.",
  "- Reverse causality: state-expression direction cannot establish mechanism."
)

report_output <- file.path(
  result_dir,
  "PRJNA607098_F7_validation_decision_v3.md"
)
writeLines(report_lines, report_output, useBytes = TRUE)

message("Step 08B3 final F7 extraction and validation completed.")
message("Candidate results: ", candidate_output)
message("Module summary: ", module_output)
message("Chunk manifest: ", chunk_manifest_path)
message("Evidence decision: ", report_output)
