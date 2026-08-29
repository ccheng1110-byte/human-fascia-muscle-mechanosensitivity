options(stringsAsFactors = FALSE)

project_dir <- "."
matrix_path <- file.path(
  project_dir, "data", "processed", "mechanochemical_validation", "GSE300230",
  "GSE300230_raw_counts_tensionTGFb.csv.gz"
)
sample_map_path <- file.path(project_dir, "config", "step15_GSE300230_sample_map_v1.csv")
module_path <- file.path(project_dir, "config", "mechanotransduction_module_registry_v2.csv")
candidate_path <- file.path(project_dir, "config", "frozen_candidate_panel_v2.csv")
result_dir <- file.path(
  project_dir, "results", "12_computational_strengthening",
  "15D_GSE300230_score_sensitivity_and_negative_controls"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(matrix_path, sample_map_path, module_path, candidate_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input(s): ", paste(missing_files, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

message("Loading the frozen GSE300230 matrix for score sensitivity and negative controls.")
counts_df <- utils::read.csv(
  gzfile(matrix_path),
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
counts <- as.matrix(counts_df)
storage.mode(counts) <- "numeric"
rownames(counts) <- trimws(rownames(counts))
if (anyDuplicated(rownames(counts))) {
  counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
}

sample_map <- utils::read.csv(sample_map_path, check.names = FALSE, stringsAsFactors = FALSE)
module_registry <- utils::read.csv(module_path, check.names = FALSE, stringsAsFactors = FALSE)
candidate_panel <- utils::read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
sample_map <- sample_map[match(colnames(counts), sample_map$matrix_column), , drop = FALSE]
if (anyNA(sample_map$matrix_column)) stop("The frozen sample map does not cover all matrix columns.")

alias_rules <- data.frame(
  registry_symbol = c("CTGF", "CYR61"),
  current_symbol = c("CCN2", "CCN1"),
  action = "not_needed",
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(alias_rules))) {
  old <- alias_rules$registry_symbol[[i]]
  current <- alias_rules$current_symbol[[i]]
  if (!old %in% rownames(counts) && current %in% rownames(counts)) {
    counts <- rbind(counts, counts[current, , drop = FALSE])
    rownames(counts)[nrow(counts)] <- old
    alias_rules$action[[i]] <- paste0("mapped_", current, "_to_", old)
  } else if (old %in% rownames(counts)) {
    alias_rules$action[[i]] <- "registry_symbol_present"
  } else {
    alias_rules$action[[i]] <- "both_absent"
  }
}
safe_write_csv(alias_rules, file.path(result_dir, "GSE300230_gene_alias_application_v1.csv"))

library_size <- colSums(counts)
if (any(library_size <= 0)) stop("A sample has a non-positive library size.")
log_cpm <- log2(t(t(counts + 0.5) / (library_size + 0.5) * 1e6) + 1)
keep <- rowSums(t(t(counts) / library_size * 1e6 >= 1)) >= 2L
log_cpm <- log_cpm[keep, , drop = FALSE]
if (nrow(log_cpm) < 100L) stop("The expression filter retained fewer than 100 genes.")

module_genes <- strsplit(module_registry$genes, ";", fixed = TRUE)
names(module_genes) <- module_registry$module
module_genes <- lapply(module_genes, intersect, y = rownames(log_cpm))
module_genes <- lapply(module_genes, unique)

z_matrix <- t(scale(t(log_cpm)))
z_matrix[is.na(z_matrix)] <- 0
rank_matrix <- t(apply(log_cpm, 1L, rank, ties.method = "average"))
rank_matrix <- t(scale(t(rank_matrix)))
rank_matrix[is.na(rank_matrix)] <- 0
centered_matrix <- sweep(log_cpm, 1L, rowMeans(log_cpm), FUN = "-")

score_method <- function(g, method) {
  g <- intersect(g, rownames(log_cpm))
  if (length(g) == 0L) return(rep(NA_real_, ncol(log_cpm)))
  mat <- switch(
    method,
    zmean = z_matrix[g, , drop = FALSE],
    zmedian = z_matrix[g, , drop = FALSE],
    rankmean = rank_matrix[g, , drop = FALSE],
    centered_mean = centered_matrix[g, , drop = FALSE]
  )
  if (method == "zmedian") apply(mat, 2L, median, na.rm = TRUE) else colMeans(mat, na.rm = TRUE)
}

methods <- c("zmean", "zmedian", "rankmean", "centered_mean")
score_rows <- list()
contrast_rows <- list()
row_id <- 0L

contrast_value <- function(x, meta, contrast_id) {
  if (contrast_id == "mechanical_tension_vs_relaxed_tgfb_absent") {
    return(mean(x[meta$analysis_group == "tension_absent"]) - mean(x[meta$analysis_group == "relaxed_absent"]))
  }
  if (contrast_id == "mechanical_tension_vs_relaxed_tgfb_present") {
    return(mean(x[meta$analysis_group == "tension_present"]) - mean(x[meta$analysis_group == "relaxed_present"]))
  }
  if (contrast_id == "mechanical_x_tgfb") {
    return((mean(x[meta$analysis_group == "tension_present"]) - mean(x[meta$analysis_group == "relaxed_present"])) -
      (mean(x[meta$analysis_group == "tension_absent"]) - mean(x[meta$analysis_group == "relaxed_absent"])))
  }
  stop("Unknown contrast: ", contrast_id)
}

contrast_ids <- c(
  "mechanical_tension_vs_relaxed_tgfb_absent",
  "mechanical_tension_vs_relaxed_tgfb_present",
  "mechanical_x_tgfb"
)

for (method in methods) {
  for (module_id in names(module_genes)) {
    scores <- score_method(module_genes[[module_id]], method)
    row_id <- row_id + 1L
    score_rows[[row_id]] <- data.frame(
      module = module_id,
      method = method,
      genes_used = length(module_genes[[module_id]]),
      setNames(as.list(scores), colnames(log_cpm)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    for (line_id in unique(sample_map$cell_line)) {
      meta <- sample_map[sample_map$cell_line == line_id, , drop = FALSE]
      x <- scores[match(meta$matrix_column, colnames(log_cpm))]
      for (contrast_id in contrast_ids) {
        contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
          module = module_id,
          role = module_registry$role[match(module_id, module_registry$module)],
          method = method,
          cell_line = line_id,
          contrast = contrast_id,
          contrast_value = contrast_value(x, meta, contrast_id),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

score_table <- do.call(rbind, score_rows)
contrast_table <- do.call(rbind, contrast_rows)
rownames(score_table) <- NULL
rownames(contrast_table) <- NULL
safe_write_csv(score_table, file.path(result_dir, "GSE300230_module_scores_by_method_v1.csv"))
safe_write_csv(contrast_table, file.path(result_dir, "GSE300230_module_score_contrasts_by_method_v1.csv"))

primary <- contrast_table[contrast_table$contrast == "mechanical_tension_vs_relaxed_tgfb_absent" & contrast_table$method == "zmean", , drop = FALSE]
primary_wide <- reshape(
  primary[, c("module", "cell_line", "contrast_value")],
  idvar = "module", timevar = "cell_line", direction = "wide"
)
names(primary_wide) <- sub("contrast_value\\.", "primary_effect_", names(primary_wide))
primary_wide$positive_both_lines <- with(primary_wide, primary_effect_GM08401 > 0 & primary_effect_GM09503 > 0)
safe_write_csv(primary_wide, file.path(result_dir, "GSE300230_primary_module_score_summary_v1.csv"))

loo_rows <- list()
for (module_id in names(module_genes)) {
  genes <- module_genes[[module_id]]
  if (length(genes) < 2L) next
  for (gene_id in genes) {
    reduced <- setdiff(genes, gene_id)
    scores <- score_method(reduced, "zmean")
    for (line_id in unique(sample_map$cell_line)) {
      meta <- sample_map[sample_map$cell_line == line_id, , drop = FALSE]
      x <- scores[match(meta$matrix_column, colnames(log_cpm))]
      loo_rows[[length(loo_rows) + 1L]] <- data.frame(
        module = module_id,
        omitted_gene = gene_id,
        cell_line = line_id,
        primary_effect = contrast_value(x, meta, "mechanical_tension_vs_relaxed_tgfb_absent"),
        stringsAsFactors = FALSE
      )
    }
  }
}
loo_table <- do.call(rbind, loo_rows)
safe_write_csv(loo_table, file.path(result_dir, "GSE300230_module_leave_one_gene_out_v1.csv"))
loo_summary <- do.call(rbind, lapply(split(loo_table, list(loo_table$module, loo_table$cell_line)), function(x) {
  data.frame(
    module = x$module[[1L]],
    cell_line = x$cell_line[[1L]],
    min_primary_effect = min(x$primary_effect, na.rm = TRUE),
    max_primary_effect = max(x$primary_effect, na.rm = TRUE),
    positive_fraction = mean(x$primary_effect > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(loo_summary) <- NULL
safe_write_csv(loo_summary, file.path(result_dir, "GSE300230_module_leave_one_gene_out_summary_v1.csv"))

set.seed(1504)
negative_rows <- list()
gene_mean <- rowMeans(log_cpm)
gene_bins <- cut(
  gene_mean,
  breaks = unique(stats::quantile(gene_mean, probs = seq(0, 1, length.out = 11), na.rm = TRUE)),
  include.lowest = TRUE,
  labels = FALSE
)
names(gene_bins) <- names(gene_mean)
available_genes <- names(gene_mean)
for (module_id in names(module_genes)) {
  genes <- module_genes[[module_id]]
  n_genes <- length(genes)
  if (n_genes < 2L) next
  module_bin <- gene_bins[genes]
  for (draw_id in seq_len(500L)) {
    sampled <- unlist(lapply(unique(module_bin), function(bin_id) {
      pool <- setdiff(available_genes[gene_bins == bin_id], genes)
      needed <- sum(module_bin == bin_id)
      if (length(pool) < needed) return(character())
      sample(pool, needed, replace = FALSE)
    }), use.names = FALSE)
    if (length(sampled) != n_genes) next
    scores <- score_method(sampled, "zmean")
    values <- vapply(unique(sample_map$cell_line), function(line_id) {
      meta <- sample_map[sample_map$cell_line == line_id, , drop = FALSE]
      x <- scores[match(meta$matrix_column, colnames(log_cpm))]
      contrast_value(x, meta, "mechanical_tension_vs_relaxed_tgfb_absent")
    }, numeric(1))
    negative_rows[[length(negative_rows) + 1L]] <- data.frame(
      source_module = module_id,
      draw_id = draw_id,
      GM08401_effect = values[["GM08401"]],
      GM09503_effect = values[["GM09503"]],
      positive_both_lines = all(values > 0),
      stringsAsFactors = FALSE
    )
  }
}
negative_table <- do.call(rbind, negative_rows)
safe_write_csv(negative_table, file.path(result_dir, "GSE300230_expression_matched_negative_control_draws_v1.csv"))
negative_summary <- do.call(rbind, lapply(split(negative_table, negative_table$source_module), function(x) {
  data.frame(
    source_module = x$source_module[[1L]],
    draws_completed = nrow(x),
    positive_both_line_fraction = mean(x$positive_both_lines, na.rm = TRUE),
    GM08401_median_effect = stats::median(x$GM08401_effect, na.rm = TRUE),
    GM09503_median_effect = stats::median(x$GM09503_effect, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(negative_summary) <- NULL
safe_write_csv(negative_summary, file.path(result_dir, "GSE300230_negative_control_summary_v1.csv"))

cell_cycle_effects <- primary_wide[primary_wide$module == "cell_cycle", , drop = FALSE]
core_effects <- primary_wide[primary_wide$module %in% c("integrin_focal_adhesion", "actomyosin_rho"), , drop = FALSE]
cell_intrinsic_note <- data.frame(
  audit = "cell_intrinsic_vs_composition",
  status = "NOT_ESTIMABLE",
  reason = "The input is a gene-by-sample cultured-fibroblast matrix without cell-level composition measurements or cell-state proportions.",
  stringsAsFactors = FALSE
)
safe_write_csv(cell_intrinsic_note, file.path(result_dir, "GSE300230_cell_intrinsic_vs_composition_audit_v1.csv"))

core_stable <- all(primary_wide$positive_both_lines[primary_wide$module %in% c("integrin_focal_adhesion", "actomyosin_rho")])
decision_lines <- c(
  "# Step 15D GSE300230 score sensitivity and negative-control audit",
  "",
  paste0("- Core module score direction stable across primary z-mean scoring: **", core_stable, "**."),
  "- Overall evidence grade remains: **CAUTION**.",
  "- Cell-intrinsic versus composition-associated explanation: **NOT ESTIMABLE** from this matrix.",
  "",
  "## Interpretation boundary",
  "",
  "This step tests robustness of sample-level module scoring and expression-matched negative controls. It does not create cell-level composition information, does not fit a post-hoc multivariable adjustment for cell cycle, and does not upgrade the evidence grade. The cell-cycle competitor remains a required alternative explanation.",
  "",
  "## Material Passport",
  "",
  "- Input: official GSE300230 16-sample processed raw-count matrix.",
  "- Frozen inputs: Step 15B sample map, module registry, and candidate panel.",
  "- Sensitivity grid: z-mean, z-median, rank-mean, and centered-mean module scores.",
  "- Negative controls: 500 expression-bin-matched random sets per module where feasible, seed 1504.",
  "- Leave-one-gene-out: applied to every module with at least two retained genes."
)
decision_path <- file.path(result_dir, "GSE300230_step15D_score_sensitivity_and_negative_control_decision_v1.md")
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 15D GSE300230 score sensitivity and negative-control audit completed.")
message("Core module direction stable under primary score: ", core_stable)
message("Cell-intrinsic versus composition audit: NOT_ESTIMABLE")
message("Decision: ", decision_path)

