options(stringsAsFactors = FALSE)

# Step 17C: S1 co-expression-level analysis.
#
# The PRJNA607098 helper streams only the pre-specified candidate-pair genes
# from the public Zarr. GSE130973 is read from the existing local Seurat object.
# Cell-level co-detection is descriptive; inference is summarized across
# samples/subjects. No result is labelled cell-intrinsic.

project_dir <- "."
runtime_root <- file.path(project_dir, ".runtime", "python_step17c")
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

registry_path <- file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
)
candidate_path <- file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
)
python_helper <- file.path(
  project_dir, "scripts", "17C_extract_PRJNA607098_S1_coexpression.py"
)
gse130973_object_path <- file.path(
  project_dir, "data", "raw", "independent_sources", "GSE130973",
  "GSE130973_seurat_analysis_lyko.rds.gz"
)
gse130973_cluster_config_path <- file.path(
  project_dir, "config", "GSE130973_candidate_fibroblast_clusters_v1.csv"
)
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17C_S1_coexpression_level_analysis"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  registry_path, candidate_path, python_helper,
  gse130973_object_path, gse130973_cluster_config_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 17C input(s): ", paste(missing_inputs, collapse = "; "))
}

registry <- read.csv(registry_path, check.names = FALSE)
candidate_panel <- read.csv(candidate_path, check.names = FALSE)
registry_gene_rows <- lapply(seq_len(nrow(registry)), function(i) {
  trimws(strsplit(as.character(registry$genes[[i]]), ";", fixed = TRUE)[[1L]])
})
registry_genes <- sort(unique(unlist(registry_gene_rows, use.names = FALSE)))
candidate_genes <- unique(trimws(as.character(candidate_panel$gene)))

integrin_genes <- c("ITGAV", "ITGB1", "ITGA2", "PARVA", "ITGA5", "FERMT2")
actomyosin_genes <- c("CFL1", "LIMK1", "CNN2", "CDC42")
mechanosensor_genes <- c("PIEZO2", "TMEM63B", "PANX1")
hippo_genes <- c("NF2", "TEAD1")
pair_specs <- rbind(
  expand.grid(
    pair_family = "integrin_x_actomyosin",
    gene_a = integrin_genes, gene_b = actomyosin_genes,
    stringsAsFactors = FALSE
  ),
  expand.grid(
    pair_family = "integrin_x_mechanosensor",
    gene_a = integrin_genes, gene_b = mechanosensor_genes,
    stringsAsFactors = FALSE
  ),
  expand.grid(
    pair_family = "integrin_x_hippo",
    gene_a = integrin_genes, gene_b = hippo_genes,
    stringsAsFactors = FALSE
  )
)
pair_genes <- sort(unique(c(pair_specs$gene_a, pair_specs$gene_b)))
required_pair_genes <- sort(unique(c(pair_genes, candidate_genes)))

if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("Install reticulate once with: install.packages('reticulate')")
}

# Python must be configured before it is initialized. If the current session
# already contains Python, stop clearly rather than changing its requirements.
if (isTRUE(reticulate::py_available(initialize = FALSE))) {
  stop(
    "reticulate is already initialized. Restart R, then source this Step 17C " ,
    "script before loading reticulate, Seurat, or Python-dependent packages."
  )
}

reticulate::py_require(
  packages = c(
    "zarr>=2.18,<3",
    "fsspec>=2024.6",
    "aiohttp>=3.9",
    "numpy>=1.26,<3"
  ),
  python_version = ">=3.10,<3.13",
  exclude_newer = "2026-08-26"
)

python_environment <- new.env(parent = emptyenv())
reticulate::source_python(
  python_helper,
  envir = python_environment,
  convert = TRUE
)

message("Streaming the frozen S1 candidate-pair genes for PRJNA607098.")
message("The full 25.36 GiB H5AD is not downloaded.")
remote_result <- python_environment$run_step17c_prjna(
  output_dir = result_dir,
  panel_genes = as.list(registry_genes),
  expected_samples = 12L,
  minimum_cells_per_sample = 20L,
  permutation_reps = 200L,
  random_seed = 20260826L
)

if (!file.exists(remote_result$coexpression)) {
  stop("PRJNA607098 coexpression output was not created.")
}

read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection, text = FALSE)
  on.exit(try(close(inner_connection), silent = TRUE), add = TRUE)
  readRDS(inner_connection)
}

get_assay_data <- function(object) {
  result <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", slot = "data"),
    error = function(e) NULL
  )
  if (!is.null(result)) return(result)
  result <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  if (is.null(result)) stop("Could not extract the GSE130973 RNA data layer.")
  result
}

make_strata <- function(library_values, n_cells) {
  if (is.null(library_values)) {
    return(list(indices = list(seq_len(n_cells)), stratified = FALSE))
  }
  values <- as.numeric(library_values)
  finite <- is.finite(values)
  if (mean(finite) < 0.90 || length(unique(values[finite])) < 4L) {
    return(list(indices = list(seq_len(n_cells)), stratified = FALSE))
  }
  quantiles <- as.numeric(stats::quantile(
    values, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE,
    names = FALSE, type = 7
  ))
  breaks <- unique(quantiles)
  if (length(breaks) < 3L) {
    return(list(indices = list(seq_len(n_cells)), stratified = FALSE))
  }
  codes <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  indices <- split(seq_len(n_cells), codes)
  indices <- indices[lengths(indices) > 0L]
  if (length(indices) < 2L) {
    return(list(indices = list(seq_len(n_cells)), stratified = FALSE))
  }
  list(indices = unname(indices), stratified = TRUE)
}

summarize_pair_local <- function(values, library_values, gene_a, gene_b,
                                 source_id, sample_id, target_state,
                                 pair_family, permutation_reps = 200L) {
  a <- as.logical(values[, 1] > 0)
  b <- as.logical(values[, 2] > 0)
  n_cells <- length(a)
  strata_info <- make_strata(library_values, n_cells)
  strata <- strata_info$indices
  observed_count <- sum(a & b)
  expected_count <- sum(vapply(
    strata, function(index) length(index) * mean(a[index]) * mean(b[index]),
    numeric(1)
  ))
  null_counts <- numeric(permutation_reps)
  for (replicate in seq_len(permutation_reps)) {
    null_counts[[replicate]] <- sum(vapply(
      strata,
      function(index) sum(a[index] & sample(b[index], length(index), replace = FALSE)),
      numeric(1)
    ))
  }
  observed_rate <- observed_count / n_cells
  expected_rate <- expected_count / n_cells
  data.frame(
    source = source_id,
    sample_id = sample_id,
    target_state = target_state,
    pair_family = pair_family,
    gene_a = gene_a,
    gene_b = gene_b,
    cells = n_cells,
    gene_a_detection_rate = mean(a),
    gene_b_detection_rate = mean(b),
    observed_co_detection_rate = observed_rate,
    expected_independence_rate = expected_rate,
    excess_co_detection_rate = observed_rate - expected_rate,
    co_detection_enrichment_ratio = if (expected_rate > 0) observed_rate / expected_rate else NA_real_,
    permutation_p_one_sided = (1 + sum(null_counts >= observed_count)) / (permutation_reps + 1),
    permutation_reps = permutation_reps,
    library_size_stratified = strata_info$stratified,
    library_size_field = if (is.null(library_values)) "NOT_AVAILABLE" else "metadata_library_proxy",
    stringsAsFactors = FALSE
  )
}

message("Loading the existing GSE130973 object for the independent S1 source.")
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is required for the GSE130973 S1 analysis.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required for the GSE130973 S1 analysis.")
}

gse_object <- read_double_gzip_rds(gse130973_object_path)
if (!isS4(gse_object) || !"meta.data" %in% methods::slotNames(gse_object)) {
  stop("GSE130973 object has no Seurat meta.data slot.")
}
gse_metadata <- as.data.frame(methods::slot(gse_object, "meta.data"))
if (!all(c("subj", "celltype.age") %in% names(gse_metadata))) {
  stop("GSE130973 lacks subj or celltype.age metadata.")
}
gse_metadata$subj <- as.character(gse_metadata$subj)
gse_metadata$cluster_id <- sub("_.*$", "", as.character(gse_metadata$celltype.age))
cluster_config <- read.csv(gse130973_cluster_config_path, check.names = FALSE)
candidate_clusters <- as.character(cluster_config$cluster_id)
if (length(candidate_clusters) != 5L) {
  stop("GSE130973 candidate cluster configuration is not the frozen five-cluster set.")
}

gse_expression <- get_assay_data(gse_object)
common_cells <- intersect(colnames(gse_expression), rownames(gse_metadata))
if (length(common_cells) < 1000L) {
  stop("Too few GSE130973 cells overlap between expression and metadata.")
}
gse_expression <- gse_expression[, common_cells, drop = FALSE]
gse_metadata <- gse_metadata[common_cells, , drop = FALSE]
gse_target_mask <- gse_metadata$cluster_id %in% candidate_clusters
gse_present_genes <- intersect(required_pair_genes, rownames(gse_expression))
gse_missing_genes <- setdiff(required_pair_genes, gse_present_genes)

gse_library_field <- NULL
for (candidate_field in c("nFeature_RNA", "nCount_RNA", "n_genes", "total_counts")) {
  if (candidate_field %in% names(gse_metadata)) {
    numeric_field <- suppressWarnings(as.numeric(gse_metadata[[candidate_field]]))
    if (sum(is.finite(numeric_field)) >= 0.5 * nrow(gse_metadata)) {
      gse_library_field <- candidate_field
      break
    }
  }
}

gse_rows <- list()
gse_detection_rows <- list()
row_index <- 1L
detection_index <- 1L
for (subject_id in sort(unique(gse_metadata$subj))) {
  subject_cells <- which(gse_metadata$subj == subject_id & gse_target_mask)
  if (length(subject_cells) < 20L) next
  subject_library <- if (!is.null(gse_library_field)) {
    suppressWarnings(as.numeric(gse_metadata[[gse_library_field]][subject_cells]))
  } else NULL
  subject_values <- gse_expression[gse_present_genes, subject_cells, drop = FALSE]
  for (gene in gse_present_genes) {
    detected <- as.numeric(subject_values[gene, ]) > 0
    gse_detection_rows[[detection_index]] <- data.frame(
      source = "GSE130973",
      sample_id = subject_id,
      target_state = "candidate_fibroblast_state",
      gene = gene,
      cells = length(subject_cells),
      detected_cells = sum(detected),
      detection_rate = mean(detected),
      library_size_field = if (is.null(gse_library_field)) "NOT_AVAILABLE" else gse_library_field,
      stringsAsFactors = FALSE
    )
    detection_index <- detection_index + 1L
  }
  for (pair_index in seq_len(nrow(pair_specs))) {
    gene_a <- pair_specs$gene_a[[pair_index]]
    gene_b <- pair_specs$gene_b[[pair_index]]
    if (!all(c(gene_a, gene_b) %in% gse_present_genes)) next
    pair_values <- cbind(
      as.numeric(subject_values[gene_a, ]),
      as.numeric(subject_values[gene_b, ])
    )
    gse_rows[[row_index]] <- summarize_pair_local(
      pair_values, subject_library, gene_a, gene_b,
      "GSE130973", subject_id, "candidate_fibroblast_state",
      pair_specs$pair_family[[pair_index]], permutation_reps = 200L
    )
    row_index <- row_index + 1L
  }
}

gse_coexpression <- if (length(gse_rows) > 0L) do.call(rbind, gse_rows) else data.frame()
gse_detection <- if (length(gse_detection_rows) > 0L) do.call(rbind, gse_detection_rows) else data.frame()
if (nrow(gse_coexpression) > 0L) {
  write.csv(
    gse_coexpression,
    file.path(result_dir, "GSE130973_S1_cell_level_coexpression_v1.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
  )
}
if (nrow(gse_detection) > 0L) {
  write.csv(
    gse_detection,
    file.path(result_dir, "GSE130973_S1_gene_detection_inventory_v1.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
  )
}
write.csv(
  data.frame(
    source = "GSE130973",
    gene = required_pair_genes,
    present_in_expression = required_pair_genes %in% rownames(gse_expression),
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "GSE130973_S1_targeted_gene_coverage_v1.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

remote_coexpression <- read.csv(remote_result$coexpression, check.names = FALSE)
all_coexpression <- rbind(remote_coexpression, gse_coexpression)
if (nrow(all_coexpression) == 0L) {
  stop("No S1 coexpression rows were produced.")
}

one_sided_wilcoxon <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) < 3L || all(values == 0)) return(1)
  suppressWarnings(stats::wilcox.test(
    values, mu = 0, alternative = "greater", exact = FALSE, correct = FALSE
  )$p.value)
}

split_key <- interaction(
  all_coexpression$source, all_coexpression$pair_family,
  drop = TRUE, lex.order = TRUE
)
pair_groups <- split(all_coexpression, list(
  all_coexpression$source,
  all_coexpression$pair_family,
  all_coexpression$gene_a,
  all_coexpression$gene_b
), drop = TRUE)
pair_summary_rows <- lapply(pair_groups, function(x) {
  x <- x[order(x$sample_id), , drop = FALSE]
  p_value <- one_sided_wilcoxon(x$excess_co_detection_rate)
  data.frame(
    source = x$source[[1L]],
    pair_family = x$pair_family[[1L]],
    gene_a = x$gene_a[[1L]],
    gene_b = x$gene_b[[1L]],
    samples = length(unique(x$sample_id)),
    positive_samples = sum(x$excess_co_detection_rate > 0, na.rm = TRUE),
    median_excess_co_detection_rate = median(x$excess_co_detection_rate, na.rm = TRUE),
    IQR_excess_co_detection_rate = IQR(x$excess_co_detection_rate, na.rm = TRUE),
    median_enrichment_ratio = median(x$co_detection_enrichment_ratio, na.rm = TRUE),
    median_permutation_p = median(x$permutation_p_one_sided, na.rm = TRUE),
    wilcoxon_one_sided_p = p_value,
    stringsAsFactors = FALSE
  )
})
pair_summary <- do.call(rbind, pair_summary_rows)
pair_summary$BH_FDR_within_source_pair_family <- NA_real_
for (group in unique(interaction(
  pair_summary$source, pair_summary$pair_family, drop = TRUE
))) {
  index <- which(interaction(
    pair_summary$source, pair_summary$pair_family, drop = TRUE
  ) == group)
  pair_summary$BH_FDR_within_source_pair_family[index] <- p.adjust(
    pair_summary$wilcoxon_one_sided_p[index], method = "BH"
  )
}
pair_summary$descriptive_within_sample_support <-
  pair_summary$median_excess_co_detection_rate > 0 &
  pair_summary$positive_samples >= ceiling(0.75 * pair_summary$samples) &
  pair_summary$BH_FDR_within_source_pair_family <= 0.05

write.csv(
  all_coexpression,
  file.path(result_dir, "Step17C_S1_all_source_coexpression_rows_v1.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  pair_summary,
  file.path(result_dir, "Step17C_S1_pair_level_statistics_v1.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

source_groups <- split(pair_summary, pair_summary$source)
source_summary_rows <- lapply(source_groups, function(x) {
  data.frame(
    source = x$source[[1L]],
    samples = max(x$samples, na.rm = TRUE),
    pairs_tested = nrow(x),
    descriptively_supported_pairs = sum(x$descriptive_within_sample_support, na.rm = TRUE),
    supported_pair_fraction = mean(x$descriptive_within_sample_support, na.rm = TRUE),
    median_pair_excess = median(x$median_excess_co_detection_rate, na.rm = TRUE),
    evidence_level = if (all(!is.finite(x$median_excess_co_detection_rate))) {
      "NOT_ESTIMABLE"
    } else if (sum(x$descriptive_within_sample_support, na.rm = TRUE) > 0L) {
      "WITHIN_SAMPLE_COEXPRESSION_SIGNAL"
    } else {
      "NO_PREDEFINED_PAIR_SUPPORT"
    },
    interpretation_boundary =
      "co-expression-level evidence only; not cell-intrinsic mechanism",
    stringsAsFactors = FALSE
  )
})
source_summary <- do.call(rbind, source_summary_rows)
write.csv(
  source_summary,
  file.path(result_dir, "Step17C_S1_source_level_summary_v1.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

targeted_coverage <- data.frame(
  source = c("PRJNA607098", "GSE130973"),
  registry_genes_requested = c(
    length(registry_genes), length(required_pair_genes)
  ),
  registry_genes_available_for_targeted_stream = c(
    as.integer(remote_result$found_registry_genes), length(gse_present_genes)
  ),
  candidate_genes_requested = length(candidate_genes),
  candidate_genes_available = c(
    sum(candidate_genes %in% read.csv(remote_result$gene_inventory)$gene),
    sum(candidate_genes %in% rownames(gse_expression))
  ),
  stringsAsFactors = FALSE
)
write.csv(
  targeted_coverage,
  file.path(result_dir, "Step17C_S1_coverage_summary_v1.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

decision_lines <- c(
  "# Step 17C S1 co-expression-level analysis",
  "",
  "## Decision",
  "",
  "- Step status: **COMPLETED_WITH_EVIDENCE_BOUNDARY**.",
  "- The analysis estimates within-sample co-detection patterns; it does not prove a cell-intrinsic mechanism.",
  "- PRJNA607098 and GSE130973 were analyzed separately and were not pooled into one statistical sample.",
  "",
  "## Analysis contract",
  "",
  "- Primary unit: sample/subject-level summary.",
  "- Cell-level co-detection: descriptive only.",
  "- Null construction: within-sample permutation preserving gene-level detection margins; library-size stratification was used when a valid field was available.",
  "- Primary pair families: integrin × actomyosin, integrin × mechanosensor, and integrin × hippo candidate genes.",
  "- Descriptive support flag: positive median excess co-detection, at least 75% positive sample summaries, and BH FDR ≤ 0.05 within source/pair-family. This is not an evidence-grade upgrade rule.",
  "",
  "## Evidence boundary",
  "",
  "- `WITHIN_SAMPLE_COEXPRESSION_SIGNAL` must be written as co-expression-level evidence only.",
  "- It cannot be described as cell-intrinsic causality, direct mechanosensitivity, or fascia specificity.",
  "- Detection-rate, library-depth, and cell-state composition limitations remain part of the interpretation.",
  "- The project evidence grade remains CAUTION.",
  "",
  paste0("- PRJNA607098 source cells: ", remote_result$source_cells, "; target-state cells: ", remote_result$target_state_cells, "."),
  paste0("- PRJNA607098 samples analyzed: ", remote_result$samples, "."),
  paste0("- GSE130973 subjects analyzed: ", length(unique(gse_metadata$subj[gse_target_mask])), "."),
  paste0("- GSE130973 missing targeted genes: ", paste(gse_missing_genes, collapse = ", "), "."),
  "",
  "## Material Passport",
  "",
  "- New large source file downloaded: none; only targeted Zarr chunks were streamed.",
  "- PRJNA607098 output includes the targeted chunk manifest and gene-coverage inventory.",
  "- Next step: Step 17D GSE338388 provenance/design audit for TGFβ exposure × TEAD inhibition regulatory-axis cross-validation."
)
decision_path <- file.path(result_dir, "Step17C_S1_coexpression_level_decision_v1.md")
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17C S1 co-expression-level analysis completed.")
message("PRJNA607098 target-state cells: ", remote_result$target_state_cells)
message("PRJNA607098 samples: ", remote_result$samples)
message("GSE130973 subjects: ", length(unique(gse_metadata$subj[gse_target_mask])))
message("Decision: ", decision_path)
