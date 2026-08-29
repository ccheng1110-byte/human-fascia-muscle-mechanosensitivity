options(stringsAsFactors = FALSE)

# Step 08C2: paired sample-level validation of F7 expression enrichment.
# The 12 reconstructed SRS samples are the statistical units. Individual cells
# are aggregated within sample and are never treated as independent replicates.

project_dir <- "."
python_runtime_root <- file.path(project_dir, ".runtime", "python_step08b")
r_user_cache_dir <- file.path(python_runtime_root, "r_user_cache")
uv_cache_dir <- file.path(python_runtime_root, "uv_cache")
uv_python_install_dir <- file.path(python_runtime_root, "uv_python")

for (runtime_dir in c(
  python_runtime_root,
  r_user_cache_dir,
  uv_cache_dir,
  uv_python_install_dir
)) {
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
    "reticulate cache redirection failed. Restart R and source this script ",
    "before loading reticulate or Seurat."
  )
}

step07_path <- file.path(
  project_dir, "results", "05_specificity_robustness", "GSE173252",
  "GSE173252_DD_candidate_specificity_summary_v1.csv"
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
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025", "08C2_sample_level_F7_PRJNA607098"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(step07_path, eligibility_path, python_helper)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("Install reticulate once with: install.packages('reticulate')")
}
if (utils::packageVersion("reticulate") < "1.41.0") {
  stop("reticulate >= 1.41.0 is required.")
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

as_logical_safe <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

rank_biserial_one_sample <- function(differences) {
  values <- differences[is.finite(differences) & differences != 0]
  if (length(values) == 0L) return(0)
  ranks <- rank(abs(values), ties.method = "average")
  denominator <- sum(ranks)
  (sum(ranks[values > 0]) - sum(ranks[values < 0])) / denominator
}

paired_gene_statistics <- function(x) {
  differences <- x$mean_difference[is.finite(x$mean_difference)]
  prevalence_differences <- x$prevalence_difference[
    is.finite(x$prevalence_difference)
  ]
  positive <- sum(differences > 0)
  negative <- sum(differences < 0)
  zero <- sum(differences == 0)
  nonzero <- positive + negative

  wilcoxon_p <- if (length(differences) > 0L && nonzero > 0L) {
    suppressWarnings(stats::wilcox.test(
      differences,
      mu = 0,
      alternative = "greater",
      exact = FALSE,
      correct = FALSE
    )$p.value)
  } else {
    1
  }
  sign_p <- if (nonzero > 0L) {
    stats::binom.test(
      positive,
      nonzero,
      p = 0.5,
      alternative = "greater"
    )$p.value
  } else {
    1
  }

  data.frame(
    gene = x$gene[[1L]],
    samples = length(differences),
    positive_samples = positive,
    negative_samples = negative,
    zero_difference_samples = zero,
    median_F7_minus_nonF7_mean = median(differences, na.rm = TRUE),
    IQR_F7_minus_nonF7_mean = IQR(differences, na.rm = TRUE),
    minimum_F7_minus_nonF7_mean = min(differences, na.rm = TRUE),
    maximum_F7_minus_nonF7_mean = max(differences, na.rm = TRUE),
    median_F7_minus_nonF7_prevalence = median(
      prevalence_differences,
      na.rm = TRUE
    ),
    paired_rank_biserial = rank_biserial_one_sample(differences),
    wilcoxon_one_sided_p = wilcoxon_p,
    exact_sign_one_sided_p = sign_p,
    stringsAsFactors = FALSE
  )
}

# Re-freeze the same 15 A/B candidates used in Steps 08A-08B3.
step07 <- read.csv(step07_path, check.names = FALSE)
shortlist <- step07[
  grepl("^[AB]_", step07$step06_candidate_tier) &
    step07$specificity_class == "robust_Myofib_over_all_three_identities",
  , drop = FALSE
]
shortlist <- shortlist[!duplicated(shortlist$gene), , drop = FALSE]
if (nrow(shortlist) != 15L) {
  stop("Expected 15 frozen candidates but found ", nrow(shortlist), ".")
}

eligibility <- read.csv(eligibility_path, check.names = FALSE)
eligibility$eligible_for_paired_sample_state_analysis <- as_logical_safe(
  eligibility$eligible_for_paired_sample_state_analysis
)
if (nrow(eligibility) != 12L ||
    length(unique(eligibility$sample_id)) != 12L ||
    !all(eligibility$eligible_for_paired_sample_state_analysis)) {
  stop("Step 08C1 eligibility gate is no longer satisfied.")
}

mechanosensor_context <- c("PIEZO1", "TRPV4", "PKD2")
panel_genes <- unique(c(shortlist$gene, mechanosensor_context))

reticulate::py_require(
  packages = c(
    "zarr>=2.18,<3",
    "fsspec>=2024.6",
    "aiohttp>=3.9",
    "numpy>=1.26,<3"
  ),
  python_version = ">=3.10,<3.13",
  # Keep this identical to Steps 08B/08C1 so the script can run in the same
  # already-initialized managed reticulate environment.
  exclude_newer = "2026-08-23"
)

python_environment <- new.env(parent = emptyenv())
reticulate::source_python(
  python_helper,
  envir = python_environment,
  convert = TRUE
)

message("Streaming the frozen gene panel for paired sample-level Step 08C2.")
message("The full 25.36 GiB H5AD is not downloaded.")
python_result <- python_environment$run_step08c2_extract(
  panel_genes = as.list(panel_genes),
  output_dir = result_dir,
  source_id = "PRJNA607098",
  target_state = "F7: Fascia-like myofibroblast",
  expected_samples = 12L
)

summary_path <- file.path(
  result_dir,
  "PRJNA607098_sample_F7_nonF7_gene_expression_summary_v1.csv"
)
if (!file.exists(summary_path)) {
  stop("Python extraction did not produce the expected sample summary.")
}
expression_summary <- read.csv(summary_path, check.names = FALSE)

f7 <- expression_summary[
  expression_summary$comparison_group == "F7",
  c(
    "sample_id", "gene", "cells", "mean_expression",
    "percent_cells_expressing"
  ),
  drop = FALSE
]
non_f7 <- expression_summary[
  expression_summary$comparison_group == "non_F7",
  c(
    "sample_id", "gene", "cells", "mean_expression",
    "percent_cells_expressing"
  ),
  drop = FALSE
]
names(f7)[3:5] <- c(
  "F7_cells", "F7_mean_expression", "F7_percent_cells_expressing"
)
names(non_f7)[3:5] <- c(
  "nonF7_cells", "nonF7_mean_expression", "nonF7_percent_cells_expressing"
)
paired <- merge(f7, non_f7, by = c("sample_id", "gene"), all = TRUE)
paired$mean_difference <- paired$F7_mean_expression - paired$nonF7_mean_expression
paired$prevalence_difference <-
  paired$F7_percent_cells_expressing -
  paired$nonF7_percent_cells_expressing

expected_rows <- 12L * length(panel_genes)
if (nrow(paired) != expected_rows || any(!complete.cases(paired))) {
  stop(
    "Incomplete paired sample-gene table: expected ", expected_rows,
    " complete rows, found ", nrow(paired), "."
  )
}
paired$analysis_family <- ifelse(
  paired$gene %in% shortlist$gene,
  "frozen_15_candidate_family",
  "secondary_mechanosensor_context"
)
paired <- paired[order(paired$gene, paired$sample_id), , drop = FALSE]
paired_output <- file.path(
  result_dir,
  "PRJNA607098_F7_sample_paired_gene_contrasts_v1.csv"
)
safe_write_csv(paired, paired_output)

statistics_list <- lapply(split(paired, paired$gene), paired_gene_statistics)
statistics <- do.call(rbind, statistics_list)
rownames(statistics) <- NULL
statistics$analysis_family <- ifelse(
  statistics$gene %in% shortlist$gene,
  "frozen_15_candidate_family",
  "secondary_mechanosensor_context"
)
statistics$wilcoxon_BH_FDR_within_frozen_15 <- NA_real_
candidate_index <- statistics$analysis_family == "frozen_15_candidate_family"
statistics$wilcoxon_BH_FDR_within_frozen_15[candidate_index] <- p.adjust(
  statistics$wilcoxon_one_sided_p[candidate_index],
  method = "BH"
)

# Frozen strict-support rule: positive median paired effect, directional support
# in at least 10/12 samples, and BH FDR <= 0.05 across the 15 candidates.
statistics$strict_sample_level_support <- FALSE
statistics$strict_sample_level_support[candidate_index] <-
  statistics$median_F7_minus_nonF7_mean[candidate_index] > 0 &
  statistics$positive_samples[candidate_index] >= 10L &
  statistics$wilcoxon_BH_FDR_within_frozen_15[candidate_index] <= 0.05

candidate_metadata <- shortlist[, c(
  "module", "gene", "step06_candidate_tier",
  "myofib_percent_cells_expressing",
  "minimum_difference_over_all_donors_and_identities"
), drop = FALSE]
candidate_statistics <- merge(
  candidate_metadata,
  statistics[statistics$analysis_family == "frozen_15_candidate_family", ],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
candidate_statistics <- candidate_statistics[
  match(shortlist$gene, candidate_statistics$gene),
  , drop = FALSE
]
candidate_output <- file.path(
  result_dir,
  "PRJNA607098_F7_sample_paired_candidate_statistics_v1.csv"
)
safe_write_csv(candidate_statistics, candidate_output)

context_statistics <- statistics[
  statistics$analysis_family == "secondary_mechanosensor_context",
  , drop = FALSE
]
context_output <- file.path(
  result_dir,
  "PRJNA607098_F7_sample_mechanosensor_context_statistics_v1.csv"
)
safe_write_csv(context_statistics, context_output)

module_rows <- lapply(
  split(candidate_statistics, candidate_statistics$module),
  function(x) {
    data.frame(
      module = x$module[[1L]],
      candidate_genes = nrow(x),
      strictly_supported_genes = sum(x$strict_sample_level_support %in% TRUE),
      support_fraction = mean(x$strict_sample_level_support %in% TRUE),
      median_paired_effect_across_genes = median(
        x$median_F7_minus_nonF7_mean,
        na.rm = TRUE
      ),
      median_rank_biserial_across_genes = median(
        x$paired_rank_biserial,
        na.rm = TRUE
      ),
      supported_gene_names = paste(
        x$gene[x$strict_sample_level_support %in% TRUE],
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }
)
module_summary <- do.call(rbind, module_rows)
module_summary <- module_summary[
  order(-module_summary$support_fraction, module_summary$module),
  , drop = FALSE
]
rownames(module_summary) <- NULL
module_output <- file.path(
  result_dir,
  "PRJNA607098_F7_sample_paired_module_summary_v1.csv"
)
safe_write_csv(module_summary, module_output)

n_supported <- sum(candidate_statistics$strict_sample_level_support %in% TRUE)
integrin_rows <- candidate_statistics[
  candidate_statistics$module == "integrin_focal_adhesion",
  , drop = FALSE
]
n_integrin_supported <- sum(integrin_rows$strict_sample_level_support %in% TRUE)
piezo2_rows <- candidate_statistics[candidate_statistics$gene == "PIEZO2", ]
piezo2_supported <- nrow(piezo2_rows) == 1L &&
  isTRUE(piezo2_rows$strict_sample_level_support[[1L]])
gate_pass <- n_supported >= 8L &&
  n_integrin_supported >= 4L &&
  piezo2_supported

evidence_grade <- if (gate_pass) {
  "CAUTION_PLUS_SAMPLE_REPLICATED"
} else {
  "CAUTION"
}

report_lines <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate",
  "- Origin Date: 2026-08-24",
  "- Verification Status: ANALYZED",
  "- Version Label: PRJNA607098_F7_sample_level_v1",
  "",
  "## Step 08C2 paired sample-level F7 validation",
  "",
  "### Frozen analysis",
  "",
  "- Statistical units: 12 reconstructed SRS samples; cells are aggregated within sample.",
  "- Primary contrast: within-sample mean expression in F7 minus pooled non-F7 cells.",
  "- Directional test: one-sided one-sample Wilcoxon signed-rank test across 12 paired differences.",
  "- Multiplicity: Benjamini-Hochberg correction across the frozen 15-candidate family.",
  "- Strict gene support requires median difference > 0, at least 10/12 positive samples, and BH FDR <= 0.05.",
  "- PIEZO1, TRPV4 and PKD2 are secondary context genes and do not enter the frozen gate.",
  "",
  "### Pre-specified gate",
  "",
  "- At least 8/15 frozen candidates must pass strict sample-level support.",
  "- At least 4/6 integrin/focal-adhesion candidates must pass.",
  "- PIEZO2 must pass.",
  "",
  "### Results",
  "",
  paste0("- Frozen candidates supported: ", n_supported, "/15."),
  paste0("- Integrin/focal-adhesion candidates supported: ", n_integrin_supported, "/6."),
  paste0("- PIEZO2 supported: ", piezo2_supported, "."),
  paste0("- Step 08C2 gate passed: ", gate_pass, "."),
  paste0("- Evidence grade: ", evidence_grade, "."),
  "",
  "### Evidence boundary",
  "",
  "- Passing upgrades the result from cell-descriptive CAUTION to sample-replicated CAUTION+, not to causal evidence.",
  "- SRS accessions are sample units; donor independence still requires explicit accession-to-donor verification.",
  "- The contrast uses normalized atlas expression, not raw-count pseudobulk, and pooled non-F7 cells are composition-dependent.",
  "- The data are observational and disease-source specific; protein abundance, channel activity, force response and causality remain untested.",
  "",
  "### Statistical fallacy scan (11/11 checked)",
  "",
  "- Simpson's paradox: reduced by within-sample pairing, but non-F7 state composition remains a sensitivity issue.",
  "- Ecological fallacy: sample-level expression is not generalized to individual clinical outcomes.",
  "- Berkson's paradox: surgical and atlas inclusion may induce selection bias.",
  "- Collider bias: no adjusted causal model is fitted.",
  "- Base-rate neglect: sample and cell counts are retained in the outputs.",
  "- Regression to the mean: no repeated extreme-value selection is used.",
  "- Survivorship bias: atlas QC may exclude low-quality cells or samples.",
  "- Look-elsewhere effect: controlled within the frozen 15-gene family by BH FDR.",
  "- Garden of forking paths: direction, family, support rule and gate were frozen before 08C2 extraction.",
  "- Correlation is not causation: causal claims are prohibited.",
  "- Reverse causality: expression-state association cannot establish mechanosensor activation direction."
)
decision_output <- file.path(
  result_dir,
  "PRJNA607098_F7_sample_level_validation_decision_v1.md"
)
writeLines(report_lines, decision_output, useBytes = TRUE)

message("Step 08C2 paired sample-level validation completed.")
message("Samples: 12")
message("Frozen candidates supported: ", n_supported, "/15")
message("Integrin/focal-adhesion supported: ", n_integrin_supported, "/6")
message("PIEZO2 supported: ", piezo2_supported)
message("Step 08C2 gate passed: ", gate_pass)
message("Evidence grade: ", evidence_grade)
message("Candidate statistics: ", candidate_output)
message("Module summary: ", module_output)
message("Paired contrasts: ", paired_output)
message("Evidence decision: ", decision_output)
