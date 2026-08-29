options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

# Step 17H: S3 stiffness analyses.
# GSE123100: descriptive HTM stiffness dose-response form analysis, restricted
# to the 11 cultured HTM samples; the eight clinical tissue samples are not used.
# GSE276045: WI-38 stiffness analysis by WT/hTERT cell-model condition and
# timepoint. WT/hTERT is not treated as a proliferation measurement.
# The two sources are analyzed separately and are not pooled as replicates.

project_dir <- "."
design_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17F2_S3_corrected_design_audit"
)
matrix_audit_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17G3_S3_processed_matrix_mapping_reaudit"
)
processed_dir <- file.path(
  project_dir, "data", "processed", "independent_sources", "S3"
)
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17H_S3_stiffness_analyses"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

registry_path <- file.path(project_dir, "config", "mechanotransduction_module_registry_v2.csv")
candidate_path <- file.path(project_dir, "config", "frozen_candidate_panel_v2.csv")
required_inputs <- c(
  registry_path, candidate_path,
  file.path(design_dir, "GSE123100_corrected_S3_design_v2.csv"),
  file.path(design_dir, "GSE276045_corrected_S3_design_v2.csv"),
  file.path(matrix_audit_dir, "GSE123100", "GSE123100_S3_processed_matrix_inventory_v1.csv"),
  file.path(matrix_audit_dir, "GSE276045", "GSE276045_S3_processed_matrix_inventory_v1.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 17H input(s): ", paste(missing_inputs, collapse = "; "))
}
if (!requireNamespace("limma", quietly = TRUE)) {
  stop(
    "Package 'limma' is required. Install once with: ",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); ",
    "BiocManager::install('limma')"
  )
}
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop(
    "Package 'edgeR' is required for the GSE276045 count matrix. Install once with: ",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); ",
    "BiocManager::install('edgeR')"
  )
}

safe_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
write_csv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE)
}
normalize_gene_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\.[0-9]+$", "", x, perl = TRUE)
  x[is.na(x)] <- ""
  x
}
read_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(
    con, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "\"", comment.char = "", na.strings = c("", "NA", "NaN")
  )
}
read_source_inputs <- function(gse_id) {
  inv <- utils::read.csv(
    file.path(matrix_audit_dir, gse_id, paste0(gse_id, "_S3_processed_matrix_inventory_v1.csv")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  map <- utils::read.csv(
    file.path(matrix_audit_dir, gse_id, paste0(gse_id, "_S3_sample_column_mapping_v1.csv")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  list(inventory = inv, mapping = map)
}
load_expression <- function(gse_id, gene_column, mode) {
  inputs <- read_source_inputs(gse_id)
  inv <- inputs$inventory
  map <- inputs$mapping
  matrix_path <- as.character(inv$matrix_path[[1L]])
  if (!file.exists(matrix_path)) stop("Missing local matrix: ", matrix_path)
  map <- map[map$gsm_accession != "" & !is.na(map$gsm_accession), , drop = FALSE]
  map <- map[!duplicated(map$gsm_accession), , drop = FALSE]
  sample_columns <- as.character(map$matrix_column)
  df <- read_matrix(matrix_path)
  if (!gene_column %in% names(df)) stop("Gene column not found in ", gse_id, ": ", gene_column)
  missing_cols <- setdiff(sample_columns, names(df))
  if (length(missing_cols) > 0L) stop("Sample column(s) missing in ", gse_id, ": ", paste(missing_cols, collapse = ", "))
  expr <- do.call(cbind, lapply(df[sample_columns], function(x) suppressWarnings(as.numeric(as.character(x)))))
  colnames(expr) <- sample_columns
  genes <- normalize_gene_id(df[[gene_column]])
  valid <- nzchar(genes)
  expr <- expr[valid, , drop = FALSE]
  genes <- genes[valid]
  if (any(!is.finite(expr))) stop("Non-finite expression values in ", gse_id)
  if (anyDuplicated(genes)) {
    if (identical(mode, "counts")) {
      expr <- rowsum(expr, group = genes, reorder = FALSE)
    } else {
      n_per_gene <- table(genes)
      expr <- rowsum(expr, group = genes, reorder = FALSE)
      expr <- expr / as.numeric(n_per_gene[rownames(expr)])
    }
    genes <- rownames(expr)
  } else {
    rownames(expr) <- genes
  }
  list(expr = expr, map = map, inventory = inv, matrix_path = matrix_path)
}
module_definitions <- function(registry) {
  defs <- lapply(registry$genes, function(x) {
    unique(normalize_gene_id(trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1L]])))
  })
  names(defs) <- as.character(registry$module)
  defs
}
fit_dose_outcome <- function(y, meta) {
  out <- data.frame(
    effect_per_log10_kPa = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, spearman_rho = NA_real_, spearman_p = NA_real_,
    n_samples = 0L, model_status = "not_estimable", stringsAsFactors = FALSE
  )
  keep <- is.finite(y) & is.finite(meta$log10_stiffness_kPa) & !is.na(meta$culture_group)
  if (sum(keep) < 6L) return(out)
  d <- meta[keep, , drop = FALSE]
  d$outcome <- y[keep]
  design_matrix <- stats::model.matrix(~ log10_stiffness_kPa + culture_group, data = d)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    out$n_samples <- nrow(d)
    out$model_status <- "rank_deficient"
    return(out)
  }
  fit <- stats::lm(outcome ~ log10_stiffness_kPa + culture_group, data = d)
  co <- summary(fit)$coefficients
  term <- "log10_stiffness_kPa"
  ci <- tryCatch(stats::confint(fit, term), error = function(e) c(NA_real_, NA_real_))
  out$effect_per_log10_kPa <- co[term, "Estimate"]
  out$ci_low <- ci[[1L]]
  out$ci_high <- ci[[2L]]
  out$p_value <- co[term, "Pr(>|t|)"]
  sp <- tryCatch(stats::cor.test(d$log10_stiffness_kPa, d$outcome, method = "spearman", exact = FALSE), error = function(e) NULL)
  if (!is.null(sp)) {
    out$spearman_rho <- unname(sp$estimate)
    out$spearman_p <- sp$p.value
  }
  out$n_samples <- nrow(d)
  out$model_status <- "fit_with_culture_group_adjustment"
  out
}

registry <- utils::read.csv(registry_path, check.names = FALSE, stringsAsFactors = FALSE)
candidate_panel <- utils::read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
candidate_panel$gene <- normalize_gene_id(candidate_panel$gene)
candidate_genes <- unique(candidate_panel$gene)
module_defs <- module_definitions(registry)

# ---- GSE123100: restricted HTM stiffness-series dose-response form ----
gse123100_design <- utils::read.csv(
  file.path(design_dir, "GSE123100_corrected_S3_design_v2.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
gse123100_design <- gse123100_design[gse123100_design$eligible_for_primary_dose_response, , drop = FALSE]
gse123100_loaded <- load_expression("GSE123100", "Symbol", "rpkm")
gse123100_expr <- log2(gse123100_loaded$expr + 1)
gse123100_meta <- gse123100_design[
  match(gse123100_loaded$map$gsm_accession, gse123100_design$gsm_accession),
  , drop = FALSE
]
if (anyNA(gse123100_meta$gsm_accession)) stop("GSE123100 corrected design does not cover all mapped HTM samples.")
gse123100_meta$log10_stiffness_kPa <- log10(gse123100_meta$stiffness_value_kPa)
gse123100_meta$culture_group <- factor(gse123100_meta$culture_group, levels = c("A", "C"))

gse123100_module_rows <- lapply(names(module_defs), function(module_name) {
  genes <- module_defs[[module_name]]
  present <- intersect(genes, rownames(gse123100_expr))
  y <- if (length(present) > 0L) colMeans(gse123100_expr[present, , drop = FALSE]) else rep(NA_real_, ncol(gse123100_expr))
  stat <- fit_dose_outcome(y, gse123100_meta)
  data.frame(
    source = "GSE123100", module = module_name, role = registry$role[match(module_name, registry$module)],
    genes_total = length(genes), genes_present = length(present),
    coverage_fraction = length(present) / length(genes),
    stat, stringsAsFactors = FALSE
  )
})
gse123100_module_stats <- do.call(rbind, gse123100_module_rows)
gse123100_module_stats$q_module <- stats::p.adjust(gse123100_module_stats$p_value, method = "BH")

gse123100_candidate_rows <- lapply(candidate_genes, function(gene) {
  present <- gene %in% rownames(gse123100_expr)
  y <- if (present) gse123100_expr[gene, ] else rep(NA_real_, nrow(gse123100_meta))
  stat <- fit_dose_outcome(y, gse123100_meta)
  data.frame(source = "GSE123100", gene = gene, present_in_matrix = present, stat, stringsAsFactors = FALSE)
})
gse123100_candidate_stats <- do.call(rbind, gse123100_candidate_rows)
if (sum(is.finite(gse123100_candidate_stats$p_value)) > 0L) {
  ix <- which(is.finite(gse123100_candidate_stats$p_value))
  gse123100_candidate_stats$q_candidate <- NA_real_
  gse123100_candidate_stats$q_candidate[ix] <- stats::p.adjust(gse123100_candidate_stats$p_value[ix], method = "BH")
}

safe_write_csv(gse123100_module_stats, file.path(result_dir, "GSE123100_S3_module_dose_response_statistics_v1.csv"))
safe_write_csv(gse123100_candidate_stats, file.path(result_dir, "GSE123100_S3_candidate_dose_response_statistics_v1.csv"))
safe_write_csv(gse123100_meta, file.path(result_dir, "GSE123100_S3_analysis_sample_map_v1.csv"))

# ---- GSE276045: WI-38 stiffness slope with cell-model and timepoint terms ----
gse276045_design <- utils::read.csv(
  file.path(design_dir, "GSE276045_corrected_S3_design_v2.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
gse276045_design <- gse276045_design[gse276045_design$eligible_for_stiffness_cross_check, , drop = FALSE]
gse276045_loaded <- load_expression("GSE276045", "Genes", "counts")
gse276045_meta <- gse276045_design[
  match(gse276045_loaded$map$gsm_accession, gse276045_design$gsm_accession),
  , drop = FALSE
]
if (anyNA(gse276045_meta$gsm_accession)) {
  stop("GSE276045 corrected design does not cover all mapped WI-38 samples.")
}
gse276045_meta$cell_model_condition <- factor(
  gse276045_meta$cell_model_condition, levels = c("WT", "hTERT")
)
gse276045_meta$timepoint_factor <- factor(gse276045_meta$timepoint)
gse276045_meta$log10_stiffness_kPa <- as.numeric(gse276045_meta$log10_stiffness_kPa)
if (anyNA(gse276045_meta$cell_model_condition) || anyNA(gse276045_meta$timepoint_factor)) {
  stop("GSE276045 contains unresolved cell-model or timepoint labels.")
}
if (any(gse276045_loaded$expr < 0, na.rm = TRUE)) {
  stop("GSE276045 count matrix contains negative values.")
}

gse276045_dge <- edgeR::DGEList(counts = gse276045_loaded$expr)
gse276045_design_matrix <- stats::model.matrix(
  ~ log10_stiffness_kPa * cell_model_condition + timepoint_factor,
  data = gse276045_meta
)
if (qr(gse276045_design_matrix)$rank < ncol(gse276045_design_matrix)) {
  stop("GSE276045 analysis design is rank deficient.")
}
gse276045_keep <- rowSums(edgeR::cpm(gse276045_dge) >= 1) >= 3L
if (sum(gse276045_keep) < 100L) {
  stop("Too few GSE276045 genes remain after the expression filter: ", sum(gse276045_keep))
}
gse276045_dge <- gse276045_dge[gse276045_keep, , keep.lib.sizes = FALSE]
gse276045_dge <- edgeR::calcNormFactors(gse276045_dge)
gse276045_voom <- limma::voom(
  gse276045_dge, design = gse276045_design_matrix, plot = FALSE
)

gse276045_interaction_col <- grep(
  "log10_stiffness_kPa:cell_model_conditionhTERT",
  colnames(gse276045_design_matrix), fixed = TRUE, value = TRUE
)
if (length(gse276045_interaction_col) != 1L) {
  stop("Could not resolve the GSE276045 stiffness-by-cell-model interaction term.")
}
gse276045_contrast_names <- c(
  "stiffness_slope_WT", "stiffness_slope_hTERT", "stiffness_x_cell_model"
)
gse276045_contrasts <- matrix(
  0, nrow = ncol(gse276045_design_matrix), ncol = length(gse276045_contrast_names),
  dimnames = list(colnames(gse276045_design_matrix), gse276045_contrast_names)
)
gse276045_contrasts["log10_stiffness_kPa", "stiffness_slope_WT"] <- 1
gse276045_contrasts["log10_stiffness_kPa", "stiffness_slope_hTERT"] <- 1
gse276045_contrasts[gse276045_interaction_col, "stiffness_slope_hTERT"] <- 1
gse276045_contrasts[gse276045_interaction_col, "stiffness_x_cell_model"] <- 1

gse276045_fit <- limma::lmFit(gse276045_voom, gse276045_design_matrix)
gse276045_fit <- limma::contrasts.fit(gse276045_fit, gse276045_contrasts)
gse276045_fit <- limma::eBayes(gse276045_fit, robust = TRUE)

scalar_character <- function(x) {
  if (is.null(x) || length(x) == 0L) NA_character_ else as.character(x[[1L]])
}
scalar_numeric <- function(x) {
  if (is.null(x) || length(x) == 0L) NA_real_ else suppressWarnings(as.numeric(x[[1L]]))
}
module_role <- function(module_name) {
  ix <- match(as.character(module_name), as.character(registry$module))
  if (length(ix) != 1L || is.na(ix)) return(NA_character_)
  scalar_character(registry$role[ix])
}
camera_field <- function(camera_table, fields, mode = c("character", "numeric")) {
  mode <- match.arg(mode)
  field <- fields[fields %in% names(camera_table)]
  if (length(field) == 0L) {
    return(if (identical(mode, "numeric")) NA_real_ else NA_character_)
  }
  value <- camera_table[[field[[1L]]]]
  if (identical(mode, "numeric")) scalar_numeric(value) else scalar_character(value)
}
camera_vector <- function(camera_table, fields, mode = c("character", "numeric")) {
  mode <- match.arg(mode)
  field <- fields[fields %in% names(camera_table)]
  n <- nrow(camera_table)
  if (length(field) == 0L) {
    return(if (identical(mode, "numeric")) rep(NA_real_, n) else rep(NA_character_, n))
  }
  value <- camera_table[[field[[1L]]]]
  if (identical(mode, "numeric")) suppressWarnings(as.numeric(value)) else as.character(value)
}

gse276045_candidate_stats_for_contrast <- function(contrast_name) {
  tt <- limma::topTable(
    gse276045_fit, coef = contrast_name, number = Inf, sort.by = "none"
  )
  tt_gene <- normalize_gene_id(rownames(tt))
  row_for_gene <- match(candidate_genes, tt_gene)
  raw_present <- candidate_genes %in% rownames(gse276045_loaded$expr)
  analyzed_present <- candidate_genes %in% rownames(gse276045_voom$E)
  p_value <- rep(NA_real_, length(candidate_genes))
  p_value[!is.na(row_for_gene)] <- tt$P.Value[row_for_gene[!is.na(row_for_gene)]]
  out <- data.frame(
    source = "GSE276045", contrast = contrast_name, gene = candidate_genes,
    present_in_matrix = raw_present, included_after_expression_filter = analyzed_present,
    effect_log2_per_log10_kPa = NA_real_, average_expression = NA_real_,
    moderated_t = NA_real_, p_value = p_value, q_genome = NA_real_,
    q_candidate = NA_real_, stringsAsFactors = FALSE
  )
  hit <- !is.na(row_for_gene)
  out$effect_log2_per_log10_kPa[hit] <- tt$logFC[row_for_gene[hit]]
  out$average_expression[hit] <- tt$AveExpr[row_for_gene[hit]]
  out$moderated_t[hit] <- tt$t[row_for_gene[hit]]
  out$q_genome[hit] <- tt$adj.P.Val[row_for_gene[hit]]
  finite_p <- is.finite(out$p_value)
  if (any(finite_p)) out$q_candidate[finite_p] <- stats::p.adjust(out$p_value[finite_p], method = "BH")
  out
}
gse276045_candidate_stats <- do.call(
  rbind, lapply(gse276045_contrast_names, gse276045_candidate_stats_for_contrast)
)

gse276045_module_coverage <- do.call(rbind, lapply(names(module_defs), function(module_name) {
  genes <- module_defs[[module_name]]
  present_raw <- intersect(genes, rownames(gse276045_loaded$expr))
  present_analyzed <- intersect(genes, rownames(gse276045_voom$E))
  data.frame(
    module = module_name, genes_total = length(genes),
    genes_present_raw = length(present_raw),
    genes_present_after_filter = length(present_analyzed),
    coverage_fraction_raw = length(present_raw) / length(genes),
    stringsAsFactors = FALSE
  )
}))
gse276045_camera_index <- lapply(names(module_defs), function(module_name) {
  match(intersect(module_defs[[module_name]], rownames(gse276045_voom$E)), rownames(gse276045_voom$E))
})
names(gse276045_camera_index) <- names(module_defs)
gse276045_camera_index <- gse276045_camera_index[lengths(gse276045_camera_index) > 0L]

# Run camera once per contrast with all modules together. This is required for
# a valid module-level FDR across the registered modules.
gse276045_module_rows <- lapply(gse276045_contrast_names, function(contrast_name) {
  cam_one <- tryCatch(
    limma::camera(
      gse276045_voom$E,
      index = gse276045_camera_index,
      design = gse276045_design_matrix,
      contrast = gse276045_contrasts[, contrast_name],
      inter.gene.cor = 0.01
    ),
    error = function(e) NULL
  )
  if (is.null(cam_one)) {
    return(data.frame(
      source = "GSE276045", contrast = contrast_name,
      gse276045_module_coverage,
      role = vapply(gse276045_module_coverage$module, module_role, character(1)),
      direction = NA_character_, camera_p_value = NA_real_, camera_fdr = NA_real_,
      camera_status = "not_estimable", stringsAsFactors = FALSE
    ))
  }
  cam_module_names <- rownames(cam_one)
  if (is.null(cam_module_names) || length(cam_module_names) != nrow(cam_one) ||
      !all(cam_module_names %in% gse276045_module_coverage$module)) {
    cam_module_names <- names(gse276045_camera_index)[seq_len(nrow(cam_one))]
  }
  coverage_rows <- gse276045_module_coverage[
    match(cam_module_names, gse276045_module_coverage$module), , drop = FALSE
  ]
  data.frame(
    source = "GSE276045", contrast = contrast_name,
    coverage_rows,
    role = vapply(cam_module_names, module_role, character(1)),
    direction = camera_vector(cam_one, c("Direction", "direction"), "character"),
    camera_p_value = camera_vector(cam_one, c("PValue", "P.Value", "p_value"), "numeric"),
    camera_fdr = camera_vector(cam_one, c("FDR", "fdr", "adj.P.Val"), "numeric"),
    camera_status = "estimated", stringsAsFactors = FALSE
  )
})
gse276045_module_stats <- do.call(rbind, gse276045_module_rows)

gse276045_preprocessing_audit <- data.frame(
  source = "GSE276045",
  matrix_path = gse276045_loaded$matrix_path,
  gene_column = "Genes", expression_mode = "raw_counts_voom",
  metadata_samples = nrow(gse276045_meta), mapped_samples = ncol(gse276045_loaded$expr),
  raw_genes = nrow(gse276045_loaded$expr), analyzed_genes = nrow(gse276045_voom$E),
  candidate_total = length(candidate_genes),
  candidate_present_raw = sum(candidate_genes %in% rownames(gse276045_loaded$expr)),
  candidate_included_after_filter = sum(candidate_genes %in% rownames(gse276045_voom$E)),
  cell_model_levels = paste(levels(gse276045_meta$cell_model_condition), collapse = ";"),
  timepoint_levels = paste(levels(gse276045_meta$timepoint_factor), collapse = ";"),
  proliferation_state = "not_measured_or_not_resolved",
  stringsAsFactors = FALSE
)

safe_write_csv(gse276045_candidate_stats, file.path(result_dir, "GSE276045_S3_frozen_candidate_statistics_v1.csv"))
safe_write_csv(gse276045_module_stats, file.path(result_dir, "GSE276045_S3_module_camera_statistics_v1.csv"))
safe_write_csv(gse276045_meta, file.path(result_dir, "GSE276045_S3_analysis_sample_map_v1.csv"))
safe_write_csv(gse276045_preprocessing_audit, file.path(result_dir, "GSE276045_S3_preprocessing_audit_v1.csv"))

# ---- S3 combined decision ----
gse123100_candidate_coverage <- sum(gse123100_candidate_stats$present_in_matrix)
gse123100_core_module_gate <- all(
  gse123100_module_stats$coverage_fraction >= 0.8 &
    gse123100_module_stats$genes_present > 0
)
gse276045_candidate_coverage <- sum(gse276045_candidate_stats$present_in_matrix) / length(gse276045_contrast_names)
gse276045_core_module_gate <- all(
  gse276045_module_stats$coverage_fraction_raw >= 0.8 &
    gse276045_module_stats$genes_present_after_filter > 0
)

summary_table <- data.frame(
  source = c("GSE123100", "GSE276045"),
  analysis = c(
    "HTM descriptive stiffness dose-response form",
    "WI-38 stiffness slope with cell-model interaction and timepoint adjustment"
  ),
  samples = c(nrow(gse123100_meta), nrow(gse276045_meta)),
  candidate_coverage = c(
    paste0(gse123100_candidate_coverage, "/", length(candidate_genes)),
    paste0(gse276045_candidate_coverage, "/", length(candidate_genes))
  ),
  core_module_gate = c(gse123100_core_module_gate, gse276045_core_module_gate),
  source_gate = c(
    "PARTIAL_DESCRIPTIVE_DOSE_RESPONSE_WITH_ITGA5_COVERAGE_LIMIT",
    "PASS_TO_S3_EVIDENCE_SYNTHESIS_WITH_PROLIFERATION_BOUNDARY"
  ),
  interpretation_boundary = c(
    "Cross-tissue dose-response form only; one candidate gene (ITGA5) absent; not a formal balanced inference.",
    "Cross-tissue WI-38 cell-model stiffness cross-check; timepoint adjusted, proliferation-independent effect not estimable."
  ),
  evidence_grade = "CAUTION",
  stringsAsFactors = FALSE
)
safe_write_csv(summary_table, file.path(result_dir, "Step17H_S3_stiffness_analysis_summary_v1.csv"))

decision_lines <- c(
  "# Step 17H S3 stiffness analysis decision",
  "",
  "## Overall decision",
  "",
  "Overall gate: PARTIAL_TO_S3_EVIDENCE_SYNTHESIS",
  "Evidence grade: CAUTION",
  "",
  "The two S3 sources were analyzed separately and were not pooled as biological replicates.",
  "",
  "## GSE123100",
  "",
  paste0(
    "The 11 eligible cultured HTM samples support a descriptive stiffness dose-response form. ",
    "The processed matrix contains ", gse123100_candidate_coverage, "/", length(candidate_genes),
    " frozen candidates because ITGA5 is absent. The design is not treated as a balanced experiment, ",
    "so this result is not promoted to a formal mechanistic or causal claim."
  ),
  "",
  "## GSE276045",
  "",
  paste0(
    "The 178 eligible WI-38 samples support stiffness-slope modeling with WT/hTERT cell-model ",
    "interaction and timepoint adjustment. Candidate coverage is complete and the module analysis can ",
    "proceed to evidence synthesis. The source does not provide a resolved proliferation measure, so ",
    "proliferation-independent mechanistic specificity is not estimable."
  ),
  "",
  "## Interpretation boundary",
  "",
  "These analyses strengthen cross-tissue plausibility and dose-response form, but they do not establish",
  "fascia-specific causality, cell-intrinsic mechanism, or a proliferation-independent effect.",
  "The overall CAUTION grade is retained, and the future donor-repeated functional validation threshold remains unchanged.",
  "",
  "## Output files",
  "",
  "- GSE123100_S3_module_dose_response_statistics_v1.csv",
  "- GSE123100_S3_candidate_dose_response_statistics_v1.csv",
  "- GSE276045_S3_module_camera_statistics_v1.csv",
  "- GSE276045_S3_frozen_candidate_statistics_v1.csv",
  "- Step17H_S3_stiffness_analysis_summary_v1.csv"
)
writeLines(decision_lines, file.path(result_dir, "Step17H_S3_stiffness_analysis_decision_v1.md"), useBytes = TRUE)

message("Step 17H S3 stiffness analyses completed.")
message("GSE123100: ", gse123100_candidate_coverage, "/", length(candidate_genes), " frozen candidates; descriptive dose-response only.")
message("GSE276045: ", gse276045_candidate_coverage, "/", length(candidate_genes), " frozen candidates; stiffness-slope model passed with proliferation boundary.")
message("Overall gate: PARTIAL_TO_S3_EVIDENCE_SYNTHESIS; evidence grade: CAUTION")
message("Results: ", result_dir)
