options(stringsAsFactors = FALSE)

project_dir <- "."
matrix_path <- file.path(
  project_dir, "data", "processed", "mechanochemical_validation", "GSE300230",
  "GSE300230_raw_counts_tensionTGFb.csv.gz"
)
sample_map_path <- file.path(project_dir, "config", "step15_GSE300230_sample_map_v1.csv")
candidate_path <- file.path(project_dir, "config", "frozen_candidate_panel_v2.csv")
module_path <- file.path(project_dir, "config", "mechanotransduction_module_registry_v2.csv")
contract_path <- file.path(project_dir, "config", "step15_GSE300230_analysis_contract_v1.md")
result_dir <- file.path(
  project_dir, "results", "12_computational_strengthening",
  "15C_GSE300230_mechanochemical_perturbation"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(matrix_path, sample_map_path, candidate_path, module_path, contract_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required frozen input(s): ", paste(missing_files, collapse = "; "))
}

if (!requireNamespace("limma", quietly = TRUE)) {
  stop(
    "Package 'limma' is required. Install once with: ",
    "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); ",
    "BiocManager::install('limma')"
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

message("Loading the frozen 16-sample GSE300230 mechanochemical matrix.")
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
candidate_panel <- utils::read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
module_registry <- utils::read.csv(module_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!setequal(colnames(counts), sample_map$matrix_column)) {
  stop(
    "Matrix/sample-map mismatch. Matrix-only: ",
    paste(setdiff(colnames(counts), sample_map$matrix_column), collapse = ", "),
    "; map-only: ", paste(setdiff(sample_map$matrix_column, colnames(counts)), collapse = ", ")
  )
}
sample_map <- sample_map[match(colnames(counts), sample_map$matrix_column), , drop = FALSE]

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

lib_size <- colSums(counts)
if (any(lib_size <= 0)) stop("At least one sample has a non-positive library size.")
cpm <- t(t(counts) / lib_size * 1e6)
keep <- rowSums(cpm >= 1) >= 2L
counts_filtered <- counts[keep, , drop = FALSE]
if (nrow(counts_filtered) < 100L) stop("Expression filter retained fewer than 100 genes.")

module_genes <- strsplit(module_registry$genes, ";", fixed = TRUE)
names(module_genes) <- module_registry$module

contrast_names <- c(
  "mechanical_tension_vs_relaxed_tgfb_absent",
  "mechanical_tension_vs_relaxed_tgfb_present",
  "tgfb_present_vs_absent_relaxed",
  "tgfb_present_vs_absent_tension",
  "mechanical_x_tgfb"
)

all_gene_rows <- list()
candidate_rows <- list()
module_rows <- list()
design_rows <- list()
counter <- 0L

for (line_id in unique(sample_map$cell_line)) {
  meta <- sample_map[sample_map$cell_line == line_id, , drop = FALSE]
  meta <- meta[match(colnames(counts_filtered)[colnames(counts_filtered) %in% meta$matrix_column], meta$matrix_column), , drop = FALSE]
  line_counts <- counts_filtered[, meta$matrix_column, drop = FALSE]

  meta$analysis_group <- factor(
    meta$analysis_group,
    levels = c("relaxed_absent", "tension_absent", "relaxed_present", "tension_present")
  )
  meta$batch <- factor(meta$batch)
  if (anyNA(meta$analysis_group) || nlevels(meta$batch) != 2L) {
    stop("Invalid frozen factor coding for cell line ", line_id)
  }

  design <- stats::model.matrix(~ 0 + analysis_group + batch, data = meta)
  colnames(design) <- sub("^analysis_group", "group_", colnames(design))
  rownames(design) <- meta$matrix_column

  contrast_matrix <- limma::makeContrasts(
    mechanical_tension_vs_relaxed_tgfb_absent = group_tension_absent - group_relaxed_absent,
    mechanical_tension_vs_relaxed_tgfb_present = group_tension_present - group_relaxed_present,
    tgfb_present_vs_absent_relaxed = group_relaxed_present - group_relaxed_absent,
    tgfb_present_vs_absent_tension = group_tension_present - group_tension_absent,
    mechanical_x_tgfb = (group_tension_present - group_relaxed_present) -
      (group_tension_absent - group_relaxed_absent),
    levels = design
  )

  v <- limma::voom(line_counts, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit, robust = TRUE)

  design_rows[[line_id]] <- data.frame(
    cell_line = line_id,
    matrix_column = rownames(design),
    meta[, c("gsm_accession", "age_background", "mechanical_condition", "tgfb_condition", "batch", "analysis_group")],
    residual_df = fit$df.residual[1L],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  index_list <- lapply(module_genes, function(g) match(intersect(g, rownames(v$E)), rownames(v$E)))
  index_list <- lapply(index_list, function(x) x[!is.na(x)])

  for (contrast_id in contrast_names) {
    counter <- counter + 1L
    coef_id <- match(contrast_id, colnames(contrast_matrix))
    tt <- limma::topTable(fit, coef = coef_id, number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    tt$cell_line <- line_id
    tt$age_background <- unique(meta$age_background)
    tt$contrast <- contrast_id
    names(tt)[names(tt) == "logFC"] <- "log2fc"
    names(tt)[names(tt) == "AveExpr"] <- "ave_expr"
    names(tt)[names(tt) == "P.Value"] <- "p_value"
    names(tt)[names(tt) == "adj.P.Val"] <- "q_genome"
    all_gene_rows[[counter]] <- tt[, c(
      "cell_line", "age_background", "contrast", "gene", "log2fc",
      "ave_expr", "t", "p_value", "q_genome", "B"
    )]

    candidate_out <- merge(
      candidate_panel,
      tt[, c("gene", "log2fc", "ave_expr", "t", "p_value", "q_genome")],
      by = "gene", all.x = TRUE, sort = FALSE
    )
    candidate_out$cell_line <- line_id
    candidate_out$age_background <- unique(meta$age_background)
    candidate_out$contrast <- contrast_id
    candidate_out$present_after_filter <- !is.na(candidate_out$p_value)
    candidate_out$q_candidate <- NA_real_
    present <- which(candidate_out$present_after_filter)
    candidate_out$q_candidate[present] <- stats::p.adjust(candidate_out$p_value[present], method = "BH")
    candidate_rows[[counter]] <- candidate_out

    camera_out <- limma::camera(
      v,
      index = index_list,
      design = design,
      contrast = contrast_matrix[, contrast_id],
      inter.gene.cor = 0.01
    )
    camera_out$module <- rownames(camera_out)
    rownames(camera_out) <- NULL
    camera_out <- merge(
      module_registry[, c("module", "role", "analysis_use")],
      camera_out,
      by = "module", all.x = TRUE, sort = FALSE
    )
    camera_out$cell_line <- line_id
    camera_out$age_background <- unique(meta$age_background)
    camera_out$contrast <- contrast_id
    module_rows[[counter]] <- camera_out
  }
}

all_gene_stats <- do.call(rbind, all_gene_rows)
candidate_stats <- do.call(rbind, candidate_rows)
module_stats <- do.call(rbind, module_rows)
design_audit <- do.call(rbind, design_rows)
rownames(all_gene_stats) <- NULL
rownames(candidate_stats) <- NULL
rownames(module_stats) <- NULL
rownames(design_audit) <- NULL

write_csv_gz(
  all_gene_stats,
  file.path(result_dir, "GSE300230_all_gene_contrast_statistics_v1.csv.gz")
)
safe_write_csv(candidate_stats, file.path(result_dir, "GSE300230_frozen_candidate_statistics_v1.csv"))
safe_write_csv(module_stats, file.path(result_dir, "GSE300230_frozen_module_camera_statistics_v1.csv"))
safe_write_csv(design_audit, file.path(result_dir, "GSE300230_frozen_design_audit_v1.csv"))

primary_contrast <- "mechanical_tension_vs_relaxed_tgfb_absent"
primary_modules <- module_stats[module_stats$contrast == primary_contrast, , drop = FALSE]
module_split <- split(primary_modules, primary_modules$module)
module_cross_line <- do.call(rbind, lapply(module_split, function(x) {
  x <- x[match(c("GM08401", "GM09503"), x$cell_line), , drop = FALSE]
  directions <- as.character(x$Direction)
  fdr <- as.numeric(x$FDR)
  positive_both <- nrow(x) == 2L && all(directions == "Up", na.rm = FALSE)
  robust_both <- positive_both && all(fdr <= 0.05, na.rm = FALSE)
  any_q05 <- positive_both && any(fdr <= 0.05, na.rm = TRUE)
  support_class <- if (robust_both) {
    "ROBUST_BOTH_LINES"
  } else if (any_q05) {
    "DIRECTIONALLY_CONCORDANT_ONE_LINE"
  } else if (positive_both) {
    "DIRECTION_ONLY"
  } else {
    "DISCORDANT_OR_NULL"
  }
  data.frame(
    module = x$module[[1L]],
    role = x$role[[1L]],
    GM08401_direction = directions[[1L]],
    GM08401_q = fdr[[1L]],
    GM09503_direction = directions[[2L]],
    GM09503_q = fdr[[2L]],
    positive_both_lines = positive_both,
    support_class = support_class,
    stringsAsFactors = FALSE
  )
}))
rownames(module_cross_line) <- NULL
safe_write_csv(module_cross_line, file.path(result_dir, "GSE300230_primary_module_cross_line_summary_v1.csv"))

primary_candidates <- candidate_stats[candidate_stats$contrast == primary_contrast, , drop = FALSE]
candidate_cross_line <- do.call(rbind, lapply(split(primary_candidates, primary_candidates$gene), function(x) {
  x <- x[match(c("GM08401", "GM09503"), x$cell_line), , drop = FALSE]
  data.frame(
    gene = x$gene[[1L]],
    module = x$module[[1L]],
    step10B_role = x$step10B_role[[1L]],
    GM08401_log2fc = x$log2fc[[1L]],
    GM08401_q_candidate = x$q_candidate[[1L]],
    GM09503_log2fc = x$log2fc[[2L]],
    GM09503_q_candidate = x$q_candidate[[2L]],
    positive_both_lines = all(x$log2fc > 0, na.rm = FALSE),
    q05_both_lines = all(x$q_candidate <= 0.05, na.rm = FALSE),
    stringsAsFactors = FALSE
  )
}))
rownames(candidate_cross_line) <- NULL
safe_write_csv(candidate_cross_line, file.path(result_dir, "GSE300230_primary_candidate_cross_line_summary_v1.csv"))

core_modules <- c("integrin_focal_adhesion", "actomyosin_rho")
core_summary <- module_cross_line[module_cross_line$module %in% core_modules, , drop = FALSE]
core_complete <- nrow(core_summary) == length(core_modules)
gate <- if (core_complete && all(core_summary$support_class == "ROBUST_BOTH_LINES")) {
  "PASS_LIMITED_PLAUSIBILITY"
} else if (core_complete && all(core_summary$positive_both_lines)) {
  "PARTIAL_DIRECTIONAL_PLAUSIBILITY"
} else {
  "FAIL_TO_STRENGTHEN"
}

candidate_positive_both <- sum(candidate_cross_line$positive_both_lines, na.rm = TRUE)
candidate_q05_both <- sum(candidate_cross_line$q05_both_lines, na.rm = TRUE)
decision_path <- file.path(result_dir, "GSE300230_step15C_mechanochemical_perturbation_decision_v1.md")
decision_lines <- c(
  "# Step 15C GSE300230 mechanochemical perturbation decision v1",
  "",
  paste0("- Perturbation gate: **", gate, "**."),
  "- Overall manuscript evidence grade remains: **CAUTION**.",
  paste0("- Frozen candidates positive in both cell lines: ", candidate_positive_both, "/15."),
  paste0("- Frozen candidates with candidate-family q <= 0.05 in both lines: ", candidate_q05_both, "/15."),
  "",
  "## Core-module primary contrast",
  "",
  paste0(
    "- ", core_summary$module, ": ", core_summary$support_class,
    " (GM08401 ", core_summary$GM08401_direction, ", q=", signif(core_summary$GM08401_q, 4),
    "; GM09503 ", core_summary$GM09503_direction, ", q=", signif(core_summary$GM09503_q, 4), ")."
  ),
  "",
  "## Interpretation boundary",
  "",
  "This analysis tests mechanochemical response within two primary dermal fibroblast cell lines. It does not establish fascia specificity, population-level donor replication, an independent age effect, or causal relevance to human pain. Age is completely confounded with cell line. Null and discordant module/candidate results remain part of the evidence record.",
  "",
  "## Material Passport",
  "",
  "- Input: official GSE300230 16-sample processed raw-count matrix.",
  "- Contract: step15_GSE300230_analysis_contract_v1.md.",
  "- Model: line-stratified limma-voom with four mechanochemical groups plus batch.",
  "- Multiplicity: BH within nine modules and within the frozen 15-gene candidate family.",
  "- Reproducibility: all design, module, candidate, transcriptome and decision outputs are stored in this result directory."
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 15C GSE300230 mechanochemical perturbation analysis completed.")
message("Perturbation gate: ", gate)
message("Evidence grade: CAUTION")
message("Frozen candidates positive in both cell lines: ", candidate_positive_both, "/15")
message("Frozen candidates q<=0.05 in both cell lines: ", candidate_q05_both, "/15")
message("Module statistics: ", file.path(result_dir, "GSE300230_frozen_module_camera_statistics_v1.csv"))
message("Candidate statistics: ", file.path(result_dir, "GSE300230_frozen_candidate_statistics_v1.csv"))
message("Decision: ", decision_path)

