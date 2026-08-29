options(stringsAsFactors = FALSE)

# Step 10C: sample-level validation under the frozen Step 10B contract.
# Statistical units are the 12 reconstructed PRJNA607098 samples.
# Cells are aggregated within sample and are never treated as independent
# biological replicates. This script streams only the frozen gene panel.

project_dir <- "."
runtime_root <- file.path(project_dir, ".runtime", "python_step10c")
for (x in c(
  runtime_root,
  file.path(runtime_root, "r_user_cache"),
  file.path(runtime_root, "uv_cache"),
  file.path(runtime_root, "uv_python")
)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

Sys.setenv(
  R_USER_CACHE_DIR = file.path(runtime_root, "r_user_cache"),
  XDG_CACHE_HOME = file.path(runtime_root, "r_user_cache"),
  UV_CACHE_DIR = file.path(runtime_root, "uv_cache"),
  UV_PYTHON_INSTALL_DIR = file.path(runtime_root, "uv_python"),
  RETICULATE_PYTHON = "managed"
)

if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("Install reticulate once with: install.packages('reticulate')")
}

# Loading the reticulate namespace is not the same as initializing Python.
# Reuse a genuinely initialized session; only call py_require() when Python
# has not yet been initialized. This avoids the exclude_newer error caused by
# changing managed-Python requirements after initialization.
python_already_initialized <- isTRUE(
  reticulate::py_available(initialize = FALSE)
)

registry_path <- file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
)
candidate_path <- file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
)
eligibility_path <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025", "08C1_sample_metadata_audit_PRJNA607098",
  "PRJNA607098_sample_F7_eligibility_v1.csv"
)
python_helper <- file.path(
  project_dir, "scripts",
  "08C2_extract_PRJNA607098_sample_level_gene_panel.py"
)
result_dir <- file.path(
  project_dir, "results", "08_cross_tissue_validation",
  "10C_frozen_PRJNA607098"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  registry_path, candidate_path, eligibility_path, python_helper
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

as_logical_safe <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

rank_biserial_one_sample <- function(values) {
  values <- values[is.finite(values) & values != 0]
  if (length(values) == 0L) return(0)
  ranks <- rank(abs(values), ties.method = "average")
  (sum(ranks[values > 0]) - sum(ranks[values < 0])) / sum(ranks)
}

one_sided_wilcoxon <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L || all(values == 0)) return(1)
  suppressWarnings(stats::wilcox.test(
    values, mu = 0, alternative = "greater", exact = FALSE,
    correct = FALSE
  )$p.value)
}

summarize_differences <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(data.frame(
      samples = 0L, positive_samples = 0L, negative_samples = 0L,
      zero_difference_samples = 0L, median_difference = NA_real_,
      IQR_difference = NA_real_, paired_rank_biserial = NA_real_,
      wilcoxon_one_sided_p = 1
    ))
  }
  data.frame(
    samples = length(values),
    positive_samples = sum(values > 0),
    negative_samples = sum(values < 0),
    zero_difference_samples = sum(values == 0),
    median_difference = median(values),
    IQR_difference = IQR(values),
    paired_rank_biserial = rank_biserial_one_sample(values),
    wilcoxon_one_sided_p = one_sided_wilcoxon(values)
  )
}

registry <- read.csv(registry_path, check.names = FALSE)
candidate_panel <- read.csv(candidate_path, check.names = FALSE)
eligibility <- read.csv(eligibility_path, check.names = FALSE)

if (nrow(registry) != 9L || length(unique(registry$module)) != 9L) {
  stop("The Step 10B registry must contain exactly 9 unique modules.")
}
if (nrow(candidate_panel) != 15L || length(unique(candidate_panel$gene)) != 15L) {
  stop("The Step 10B candidate panel must contain exactly 15 unique genes.")
}
eligibility$eligible_for_paired_sample_state_analysis <- as_logical_safe(
  eligibility$eligible_for_paired_sample_state_analysis
)
if (nrow(eligibility) != 12L ||
    length(unique(eligibility$sample_id)) != 12L ||
    !all(eligibility$eligible_for_paired_sample_state_analysis)) {
  stop("Step 08C1 eligibility gate is no longer satisfied.")
}

registry_gene_rows <- lapply(seq_len(nrow(registry)), function(i) {
  genes <- trimws(strsplit(registry$genes[[i]], ";", fixed = TRUE)[[1L]])
  data.frame(
    module = registry$module[[i]],
    role = registry$role[[i]],
    gene = genes,
    stringsAsFactors = FALSE
  )
})
registry_genes <- do.call(rbind, registry_gene_rows)
panel_genes <- unique(registry_genes$gene)

candidate_registry_check <- merge(
  candidate_panel[, c("gene", "module")],
  registry_genes[, c("module", "gene")],
  by = c("gene", "module")
)
if (nrow(candidate_registry_check) != nrow(candidate_panel)) {
  stop("Every frozen candidate must occur in its frozen registry module.")
}

if (!python_already_initialized) {
  reticulate::py_require(
    packages = c(
      "zarr>=2.18,<3", "fsspec>=2024.6", "aiohttp>=3.9",
      "numpy>=1.26,<3"
    ),
    python_version = ">=3.10,<3.13",
    exclude_newer = "2026-08-23"
  )
} else {
  message("Python is already initialized; reusing the current reticulate environment.")
}

python_environment <- new.env(parent = emptyenv())
reticulate::source_python(
  python_helper, envir = python_environment, convert = TRUE
)

message("Streaming the Step 10B-frozen gene panel for Step 10C.")
message("The full 25.36 GiB H5AD is not downloaded.")
python_result <- python_environment$run_step08c2_extract(
  panel_genes = as.list(panel_genes),
  output_dir = result_dir,
  source_id = "PRJNA607098",
  target_state = "F7: Fascia-like myofibroblast",
  expected_samples = 12L,
  allow_missing_genes = TRUE
)

available_genes <- as.character(unlist(python_result$found_genes))
missing_genes <- as.character(unlist(python_result$missing_genes))
candidate_missing_genes <- setdiff(candidate_panel$gene, available_genes)
if (length(candidate_missing_genes) > 0L) {
  stop(
    "A frozen candidate gene is missing from the atlas: ",
    paste(candidate_missing_genes, collapse = ", ")
  )
}
safe_write_csv(
  data.frame(
    gene = panel_genes,
    present_in_atlas = panel_genes %in% available_genes,
    role = ifelse(
      panel_genes %in% candidate_panel$gene,
      "frozen_candidate", "module_gene"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "PRJNA607098_frozen_gene_panel_availability_v2.csv")
)
if (length(missing_genes) > 0L) {
  message(
    "Non-candidate module genes absent from the atlas and skipped: ",
    paste(missing_genes, collapse = ", ")
  )
}

summary_path <- file.path(
  result_dir, "PRJNA607098_sample_F7_nonF7_gene_expression_summary_v1.csv"
)
if (!file.exists(summary_path)) {
  stop("The targeted extraction did not produce the expected summary.")
}
expression_summary <- read.csv(summary_path, check.names = FALSE)

f7 <- expression_summary[
  expression_summary$comparison_group == "F7",
  c("sample_id", "gene", "cells", "mean_expression", "percent_cells_expressing"),
  drop = FALSE
]
non_f7 <- expression_summary[
  expression_summary$comparison_group == "non_F7",
  c("sample_id", "gene", "cells", "mean_expression", "percent_cells_expressing"),
  drop = FALSE
]
names(f7)[3:5] <- c("F7_cells", "F7_mean_expression", "F7_percent_cells_expressing")
names(non_f7)[3:5] <- c("nonF7_cells", "nonF7_mean_expression", "nonF7_percent_cells_expressing")
paired <- merge(f7, non_f7, by = c("sample_id", "gene"), all = TRUE)
paired$mean_difference <- paired$F7_mean_expression - paired$nonF7_mean_expression
paired$prevalence_difference <-
  paired$F7_percent_cells_expressing - paired$nonF7_percent_cells_expressing

expected_rows <- 12L * length(available_genes)
if (nrow(paired) != expected_rows || any(!complete.cases(paired))) {
  stop(
    "Incomplete paired sample-gene table: expected ", expected_rows,
    " complete rows, found ", nrow(paired), "."
  )
}
paired$analysis_family <- ifelse(
  paired$gene %in% candidate_panel$gene,
  "frozen_15_candidate_family", "frozen_module_gene_panel"
)
paired <- paired[order(paired$gene, paired$sample_id), , drop = FALSE]
safe_write_csv(
  paired,
  file.path(result_dir, "PRJNA607098_F7_frozen_gene_contrasts_v2.csv")
)

gene_statistics <- do.call(rbind, lapply(split(paired, paired$gene), function(x) {
  summary <- summarize_differences(x$mean_difference)
  data.frame(
    gene = x$gene[[1L]], summary, stringsAsFactors = FALSE
  )
}))
rownames(gene_statistics) <- NULL
gene_statistics$analysis_family <- ifelse(
  gene_statistics$gene %in% candidate_panel$gene,
  "frozen_15_candidate_family", "frozen_module_gene_panel"
)
candidate_statistics <- merge(
  candidate_panel,
  gene_statistics[gene_statistics$analysis_family == "frozen_15_candidate_family", ],
  by = "gene", all.x = TRUE, sort = FALSE
)
candidate_statistics <- candidate_statistics[
  match(candidate_panel$gene, candidate_statistics$gene), , drop = FALSE
]
candidate_statistics$BH_FDR_within_frozen_15 <- p.adjust(
  candidate_statistics$wilcoxon_one_sided_p, method = "BH"
)
candidate_statistics$strict_sample_level_support <-
  candidate_statistics$median_difference > 0 &
  candidate_statistics$positive_samples >= 10L &
  candidate_statistics$BH_FDR_within_frozen_15 <= 0.05
safe_write_csv(
  candidate_statistics,
  file.path(
    result_dir,
    "PRJNA607098_F7_frozen_candidate_statistics_v2.csv"
  )
)

module_contrasts <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
  module <- registry$module[[i]]
  genes <- registry_genes$gene[registry_genes$module == module]
  x <- paired[paired$gene %in% genes, c("sample_id", "gene", "mean_difference")]
  aggregate(mean_difference ~ sample_id, x, mean) |>
    transform(module = module) |>
    setNames(c("sample_id", "module_difference", "module"))
}))
module_contrasts <- module_contrasts[, c("module", "sample_id", "module_difference")]
module_contrasts <- module_contrasts[order(module_contrasts$module, module_contrasts$sample_id), ]
safe_write_csv(
  module_contrasts,
  file.path(result_dir, "PRJNA607098_F7_module_paired_contrasts_v2.csv")
)

module_summary <- do.call(rbind, lapply(split(module_contrasts, module_contrasts$module), function(x) {
  s <- summarize_differences(x$module_difference)
  data.frame(
    module = x$module[[1L]],
    role = registry$role[match(x$module[[1L]], registry$module)],
    s,
    stringsAsFactors = FALSE
  )
}))
rownames(module_summary) <- NULL
module_summary$BH_FDR_within_module_family <- p.adjust(
  module_summary$wilcoxon_one_sided_p, method = "BH"
)
module_summary$module_gate_pass <-
  module_summary$median_difference > 0 &
  module_summary$positive_samples >= 9L
safe_write_csv(
  module_summary,
  file.path(result_dir, "PRJNA607098_F7_module_summary_v2.csv")
)

core_modules <- c("integrin_focal_adhesion", "actomyosin_rho")
competitor_modules <- c(
  "ecm_remodeling", "tgf_fibrosis", "inflammation_ap1_nfkb",
  "hypoxia", "cell_cycle"
)

core_contrasts <- module_contrasts[module_contrasts$module %in% core_modules, ]
core_wide <- reshape(
  core_contrasts, idvar = "sample_id", timevar = "module", direction = "wide"
)
core_columns <- paste0("module_difference.", core_modules)
if (!all(core_columns %in% names(core_wide))) {
  stop("Core module contrast columns are incomplete.")
}
scale_values <- vapply(core_columns, function(column) {
  value <- stats::mad(core_wide[[column]], constant = 1, na.rm = TRUE)
  if (!is.finite(value) || value == 0) value <- stats::sd(core_wide[[column]], na.rm = TRUE)
  if (!is.finite(value) || value == 0) value <- 1
  value
}, numeric(1))
core_wide$primary_composite_difference <- rowMeans(
  sweep(core_wide[, core_columns, drop = FALSE], 2, scale_values, "/")
)

loo_rows <- list()
for (module in core_modules) {
  values <- core_wide[[paste0("module_difference.", module)]]
  for (i in seq_len(nrow(core_wide))) {
    remain <- values[-i]
    loo_rows[[length(loo_rows) + 1L]] <- data.frame(
      analysis = module,
      removed_sample = core_wide$sample_id[[i]],
      median_remaining_difference = median(remain),
      positive_remaining_samples = sum(remain > 0),
      sign_preserved = median(remain) > 0,
      stringsAsFactors = FALSE
    )
  }
}
for (i in seq_len(nrow(core_wide))) {
  remain <- core_wide$primary_composite_difference[-i]
  loo_rows[[length(loo_rows) + 1L]] <- data.frame(
    analysis = "primary_composite",
    removed_sample = core_wide$sample_id[[i]],
    median_remaining_difference = median(remain),
    positive_remaining_samples = sum(remain > 0),
    sign_preserved = median(remain) > 0,
    stringsAsFactors = FALSE
  )
}
loo_summary <- do.call(rbind, loo_rows)
safe_write_csv(
  loo_summary,
  file.path(result_dir, "PRJNA607098_F7_leave_one_sample_out_v2.csv")
)

robustness_rows <- list()
for (core in core_modules) {
  core_x <- core_wide[, c("sample_id", paste0("module_difference.", core))]
  names(core_x)[2] <- "core_difference"
  for (competitor in competitor_modules) {
    competitor_x <- module_contrasts[
      module_contrasts$module == competitor,
      c("sample_id", "module_difference")
    ]
    names(competitor_x)[2] <- "competitor_difference"
    x <- merge(core_x, competitor_x, by = "sample_id")
    fit <- stats::lm(core_difference ~ competitor_difference, data = x)
    residual <- stats::residuals(fit)
    split_point <- median(x$competitor_difference)
    high <- x$competitor_difference >= split_point
    robustness_rows[[length(robustness_rows) + 1L]] <- data.frame(
      core_module = core,
      competitor_module = competitor,
      samples = nrow(x),
      spearman_rho = suppressWarnings(stats::cor(
        x$core_difference, x$competitor_difference,
        method = "spearman", use = "complete.obs"
      )),
      high_competitor_core_median = median(x$core_difference[high]),
      low_competitor_core_median = median(x$core_difference[!high]),
      residual_median = median(residual),
      residual_positive_samples = sum(residual > 0),
      residual_support_8_of_12 = sum(residual > 0) >= 8L,
      stringsAsFactors = FALSE
    )
  }
}
robustness <- do.call(rbind, robustness_rows)
safe_write_csv(
  robustness,
  file.path(result_dir, "PRJNA607098_F7_competitor_robustness_v2.csv")
)

core_gate <- all(
  module_summary$module_gate_pass[
    match(core_modules, module_summary$module)
  ]
)
loo_gate <- all(loo_summary$sign_preserved)
competition_gate <- all(robustness$residual_support_8_of_12)
candidate_supported <- sum(candidate_statistics$strict_sample_level_support %in% TRUE)
integrin_supported <- sum(
  candidate_statistics$strict_sample_level_support[
    candidate_statistics$module == "integrin_focal_adhesion"
  ] %in% TRUE
)
piezo2_supported <- isTRUE(candidate_statistics$strict_sample_level_support[
  match("PIEZO2", candidate_statistics$gene)
])

# Gate A and Gate B remain partial from Step 10A. Therefore this script cannot
# upgrade the project to independent replication even if the expression gates pass.
bounded_support <- core_gate && loo_gate && competition_gate
evidence_grade <- if (bounded_support) {
  "SUPPORTED_WITH_BOUNDARY"
} else {
  "CAUTION"
}

decision_lines <- c(
  "# Step 10C frozen-contract validation decision",
  "",
  "## Material Passport",
  "",
  "- Origin: local Step 10B analysis contract and the validated PRJNA607098 Zarr extraction helper.",
  "- Statistical unit: 12 reconstructed PRJNA607098 sample IDs.",
  "- Target: exact atlas label `F7: Fascia-like myofibroblast`.",
  "- Comparator: within-sample pooled non-F7 cells.",
  "- Gene sets: config/mechanotransduction_module_registry_v2.csv.",
  "- Candidate panel: config/frozen_candidate_panel_v2.csv.",
  "- No raw H5AD download: only targeted Zarr chunks were streamed.",
  "",
  "## Results",
  "",
  paste0("- Frozen candidate support: ", candidate_supported, "/15."),
  paste0("- Integrin/focal-adhesion support: ", integrin_supported, "/6."),
  paste0("- PIEZO2 supported: ", piezo2_supported, "."),
  paste0("- Core module gate: ", core_gate, "."),
  paste0("- Leave-one-sample-out gate: ", loo_gate, "."),
  paste0("- Competition robustness gate: ", competition_gate, "."),
  paste0("- Bounded support gate: ", bounded_support, "."),
  paste0("- Evidence grade: ", evidence_grade, "."),
  "",
  "## Interpretation boundary",
  "",
  "- This is a sample-level observational validation, not a causal test.",
  "- Step 10A Gate A remains PARTIAL because the integrated atlas contains GSE173252 and PRJNA607098 and lacks a universal donor field.",
  "- A positive result cannot be called independent replication or direct F8 replication while the F7/F8 crosswalk remains unresolved.",
  "- PIEZO2 is evaluated as a channel-specific branch. A negative result narrows the conclusion but does not erase the integrin–actomyosin program.",
  "- Competitor residual checks are robustness audits, not causal adjustment models."
)
decision_path <- file.path(
  result_dir, "PRJNA607098_Step10C_frozen_contract_decision_v2.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 10C frozen-contract validation completed.")
message("Frozen candidate support: ", candidate_supported, "/15")
message("Integrin/focal-adhesion support: ", integrin_supported, "/6")
message("PIEZO2 supported: ", piezo2_supported)
message("Core module gate: ", core_gate)
message("Leave-one-sample-out gate: ", loo_gate)
message("Competition robustness gate: ", competition_gate)
message("Evidence grade: ", evidence_grade)
message("Candidate statistics: ", file.path(result_dir, "PRJNA607098_F7_frozen_candidate_statistics_v2.csv"))
message("Module summary: ", file.path(result_dir, "PRJNA607098_F7_module_summary_v2.csv"))
message("Competitor robustness: ", file.path(result_dir, "PRJNA607098_F7_competitor_robustness_v2.csv"))
message("Decision: ", decision_path)
