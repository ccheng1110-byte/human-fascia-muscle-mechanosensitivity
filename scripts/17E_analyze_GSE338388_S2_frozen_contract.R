options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

# Step 17E: S2 frozen-contract expression analysis for GSE338388.
# This is regulatory-axis cross-validation, not a mechanical-causality test.
# The frozen design is TGF-beta exposure (+/-) x TEAD inhibition (+/-).

project_dir <- "."
gse_id <- "GSE338388"
processed_dir <- file.path(project_dir, "data", "processed", "independent_sources", gse_id)
audit_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D4_S2_GSE338388_matrix_reaudit"
)
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17E_S2_GSE338388_expression_analysis"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

registry_path <- file.path(project_dir, "config", "mechanotransduction_module_registry_v2.csv")
candidate_path <- file.path(project_dir, "config", "frozen_candidate_panel_v2.csv")
mapping_path <- file.path(audit_dir, paste0(gse_id, "_sample_column_mapping_v2.csv"))
inventory_path <- file.path(audit_dir, paste0(gse_id, "_processed_matrix_inventory_v2.csv"))
design_path <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D_S2_GSE338388_design_audit", paste0(gse_id, "_design_factor_reconstruction_v2.csv")
)
required_inputs <- c(registry_path, candidate_path, mapping_path, inventory_path, design_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 17E input(s): ", paste(missing_inputs, collapse = "; "))
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
normalize_gene_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\.[0-9]+$", "", x, perl = TRUE)
  x[is.na(x)] <- ""
  x
}

registry <- utils::read.csv(registry_path, check.names = FALSE, stringsAsFactors = FALSE)
candidate_panel <- utils::read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
sample_map <- utils::read.csv(mapping_path, check.names = FALSE, stringsAsFactors = FALSE)
matrix_inventory <- utils::read.csv(inventory_path, check.names = FALSE, stringsAsFactors = FALSE)
design_reference <- utils::read.csv(design_path, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(matrix_inventory) != 1L) stop("Step 17D4 matrix inventory must contain one row.")

matrix_filename <- as.character(matrix_inventory$filename[[1L]])
matrix_path <- file.path(processed_dir, matrix_filename)
if (!file.exists(matrix_path)) stop("Validated local matrix is missing: ", matrix_path)

required_map_fields <- c(
  "matrix_column", "matrix_sample_number", "gsm_accession",
  "tgfb_condition", "tead_condition", "match_status"
)
missing_map_fields <- setdiff(required_map_fields, names(sample_map))
if (length(missing_map_fields) > 0L) {
  stop("Step 17D4 mapping is missing field(s): ", paste(missing_map_fields, collapse = ", "))
}
sample_map <- sample_map[order(sample_map$matrix_sample_number), , drop = FALSE]
if (
  nrow(sample_map) != 12L ||
    !identical(sort(sample_map$matrix_sample_number), seq_len(12L)) ||
    any(sample_map$match_status != "matched")
) {
  stop("Step 17D4 mapping is not a complete matched Yap164-1 ... Yap164-12 map.")
}

condition_levels <- c(
  "not_exposed__tead_not_inhibited",
  "not_exposed__tead_inhibited",
  "exposed__tead_not_inhibited",
  "exposed__tead_inhibited"
)
sample_map$condition_code <- paste(
  sample_map$tgfb_condition, sample_map$tead_condition, sep = "__"
)
if (!all(sample_map$condition_code %in% condition_levels)) {
  stop("Step 17D4 contains an unexpected factorial condition code.")
}
if (!setequal(sample_map$gsm_accession, design_reference$gsm_accession)) {
  stop("Step 17D4 mapping and corrected design accessions do not match.")
}
design_reference <- design_reference[
  match(sample_map$gsm_accession, design_reference$gsm_accession), , drop = FALSE
]
if (
  any(sample_map$tgfb_condition != design_reference$tgfb_condition) ||
    any(sample_map$tead_condition != design_reference$tead_condition)
) {
  stop("Step 17D4 mapping is inconsistent with the corrected official design.")
}

message("Loading the validated official GSE338388 normalized matrix.")
con <- gzfile(matrix_path, open = "rt")
expression_df <- tryCatch(
  utils::read.csv(
    con, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  ),
  error = function(e) {
    try(close(con), silent = TRUE)
    stop("Could not parse the official processed matrix: ", conditionMessage(e))
  }
)
close(con)

sample_columns <- as.character(sample_map$matrix_column)
missing_sample_columns <- setdiff(sample_columns, names(expression_df))
if (length(missing_sample_columns) > 0L) {
  stop("Validated sample column(s) absent from matrix: ", paste(missing_sample_columns, collapse = ", "))
}
gene_column <- as.character(matrix_inventory$resolved_gene_column[[1L]])
if (!gene_column %in% names(expression_df)) stop("Validated gene column absent: ", gene_column)

expression_matrix <- do.call(cbind, lapply(expression_df[sample_columns], function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}))
colnames(expression_matrix) <- sample_columns
gene_ids <- normalize_gene_id(expression_df[[gene_column]])
valid_gene_rows <- nzchar(gene_ids)
expression_matrix <- expression_matrix[valid_gene_rows, , drop = FALSE]
gene_ids <- gene_ids[valid_gene_rows]
if (any(!is.finite(expression_matrix))) stop("Processed matrix contains non-finite sample values.")

# Collapse duplicate gene identifiers deterministically by their mean.
duplicate_gene_count <- sum(duplicated(gene_ids))
if (duplicate_gene_count > 0L) {
  duplicate_sizes <- table(gene_ids)
  expression_matrix <- rowsum(expression_matrix, group = gene_ids, reorder = FALSE)
  expression_matrix <- expression_matrix / as.numeric(duplicate_sizes[rownames(expression_matrix)])
  gene_ids <- rownames(expression_matrix)
} else {
  rownames(expression_matrix) <- gene_ids
}

observed_values <- as.numeric(expression_matrix)
integer_fraction <- mean(abs(observed_values - round(observed_values)) < 1e-8)
value_min <- min(observed_values)
value_max <- max(observed_values)
if (value_min >= 0 && integer_fraction >= 0.95 && value_max > 50) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(
      "Matrix appears count-like. Install edgeR once, then rerun: ",
      "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); ",
      "BiocManager::install('edgeR')"
    )
  }
  analysis_matrix <- edgeR::cpm(expression_matrix, log = TRUE, prior.count = 1)
  value_transform <- "edgeR::cpm(log=TRUE, prior.count=1)"
} else if (value_min >= 0 && value_max > 20) {
  analysis_matrix <- log2(expression_matrix + 1)
  value_transform <- "log2(x+1) for non-negative normalized expression"
} else {
  analysis_matrix <- expression_matrix
  value_transform <- "as supplied; treated as log-like normalized expression"
}

row_variance <- apply(analysis_matrix, 1L, stats::var)
variance_keep <- is.finite(row_variance) & row_variance > 0
analysis_matrix <- analysis_matrix[variance_keep, , drop = FALSE]
if (nrow(analysis_matrix) < 100L) stop("Fewer than 100 variable genes remain after preprocessing.")

registry_rows <- lapply(seq_len(nrow(registry)), function(i) {
  genes <- trimws(strsplit(as.character(registry$genes[[i]]), ";", fixed = TRUE)[[1L]])
  data.frame(
    module = as.character(registry$module[[i]]),
    role = as.character(registry$role[[i]]),
    analysis_use = as.character(registry$analysis_use[[i]]),
    gene = genes,
    stringsAsFactors = FALSE
  )
})
registry_gene_table <- do.call(rbind, registry_rows)
registry_genes <- unique(normalize_gene_id(registry_gene_table$gene))
candidate_panel$gene_original <- as.character(candidate_panel$gene)
candidate_panel$gene <- normalize_gene_id(candidate_panel$gene)
candidate_genes <- unique(candidate_panel$gene)
processed_gene_ids <- rownames(expression_matrix)
analysis_gene_ids <- rownames(analysis_matrix)

candidate_coverage <- data.frame(
  source = gse_id,
  gene = candidate_genes,
  present_in_processed_matrix = candidate_genes %in% processed_gene_ids,
  present_after_variance_filter = candidate_genes %in% analysis_gene_ids,
  stringsAsFactors = FALSE
)
safe_write_csv(
  candidate_coverage,
  file.path(result_dir, paste0(gse_id, "_S2_candidate_coverage_v1.csv"))
)

gene_coverage <- data.frame(
  source = gse_id,
  gene = registry_genes,
  is_frozen_candidate = registry_genes %in% candidate_genes,
  present_in_processed_matrix = registry_genes %in% processed_gene_ids,
  present_after_variance_filter = registry_genes %in% analysis_gene_ids,
  stringsAsFactors = FALSE
)
safe_write_csv(
  gene_coverage,
  file.path(result_dir, paste0(gse_id, "_S2_registry_gene_coverage_v1.csv"))
)

module_coverage <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
  genes <- normalize_gene_id(
    trimws(strsplit(as.character(registry$genes[[i]]), ";", fixed = TRUE)[[1L]])
  )
  present_processed <- genes %in% processed_gene_ids
  present_analysis <- genes %in% analysis_gene_ids
  data.frame(
    source = gse_id,
    module = as.character(registry$module[[i]]),
    role = as.character(registry$role[[i]]),
    genes_total = length(genes),
    genes_present_processed = sum(present_processed),
    genes_present_after_variance_filter = sum(present_analysis),
    coverage_fraction = mean(present_analysis),
    core_module = as.character(registry$module[[i]]) %in% c(
      "integrin_focal_adhesion", "actomyosin_rho", "hippo_yap_taz"
    ),
    stringsAsFactors = FALSE
  )
}))
module_coverage$core_module_gate <- ifelse(
  module_coverage$core_module,
  module_coverage$coverage_fraction >= 0.80,
  NA
)
safe_write_csv(
  module_coverage,
  file.path(result_dir, paste0(gse_id, "_S2_module_coverage_v1.csv"))
)

meta <- sample_map
meta$condition_code <- factor(meta$condition_code, levels = condition_levels)
design_matrix <- model.matrix(~ 0 + condition_code, data = meta)
colnames(design_matrix) <- paste0("group_", condition_levels)
rownames(design_matrix) <- meta$matrix_column

# Condition order: control, TEAD inhibition, TGF-beta, TGF-beta plus TEAD
# inhibition. These are marginal main effects and the factorial interaction.
contrast_matrix <- cbind(
  tgfb_main_effect = c(-0.5, -0.5, 0.5, 0.5),
  tead_main_effect = c(-0.5, 0.5, -0.5, 0.5),
  tgfb_x_tead_interaction = c(1, -1, -1, 1)
)
rownames(contrast_matrix) <- colnames(design_matrix)

design_audit <- data.frame(
  matrix_column = meta$matrix_column,
  gsm_accession = meta$gsm_accession,
  tgfb_condition = meta$tgfb_condition,
  tead_condition = meta$tead_condition,
  condition_code = as.character(meta$condition_code),
  stringsAsFactors = FALSE
)
safe_write_csv(
  design_audit,
  file.path(result_dir, paste0(gse_id, "_S2_design_audit_v1.csv"))
)

fit <- limma::lmFit(analysis_matrix, design_matrix)
fit <- limma::contrasts.fit(fit, contrast_matrix)
fit <- limma::eBayes(fit, robust = TRUE)
contrast_names <- colnames(contrast_matrix)
all_gene_rows <- list()
candidate_rows <- list()
module_rows <- list()

module_index <- lapply(seq_len(nrow(registry)), function(i) {
  genes <- normalize_gene_id(
    trimws(strsplit(as.character(registry$genes[[i]]), ";", fixed = TRUE)[[1L]])
  )
  matched <- match(genes, analysis_gene_ids)
  unique(matched[!is.na(matched)])
})
names(module_index) <- as.character(registry$module)
camera_index <- module_index[lengths(module_index) > 0L]

for (contrast_id in contrast_names) {
  coef_id <- match(contrast_id, contrast_names)
  tt <- limma::topTable(fit, coef = coef_id, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$contrast <- contrast_id
  names(tt)[names(tt) == "logFC"] <- "effect_log2"
  names(tt)[names(tt) == "AveExpr"] <- "average_expression"
  names(tt)[names(tt) == "P.Value"] <- "p_value"
  names(tt)[names(tt) == "adj.P.Val"] <- "q_genome"
  all_gene_rows[[contrast_id]] <- tt[, c(
    "contrast", "gene", "effect_log2", "average_expression", "t",
    "p_value", "q_genome", "B"
  )]

  candidate_out <- merge(
    candidate_panel,
    tt[, c("gene", "effect_log2", "average_expression", "t", "p_value", "q_genome")],
    by = "gene", all.x = TRUE, sort = FALSE
  )
  candidate_out$contrast <- contrast_id
  candidate_out$present_after_variance_filter <- !is.na(candidate_out$p_value)
  candidate_out$q_candidate <- NA_real_
  candidate_rows[[contrast_id]] <- candidate_out

  module_camera <- tryCatch(
    limma::camera(
      analysis_matrix,
      index = camera_index,
      design = design_matrix,
      contrast = contrast_matrix[, contrast_id],
      inter.gene.cor = 0.01
    ),
    error = function(e) {
      data.frame(
        module = character(0), NGenes = numeric(0), Direction = character(0),
        PValue = numeric(0), FDR = numeric(0),
        camera_error = conditionMessage(e), stringsAsFactors = FALSE
      )
    }
  )
  if (nrow(module_camera) > 0L) {
    module_camera$module <- rownames(module_camera)
    rownames(module_camera) <- NULL
  }
  module_camera <- merge(
    registry[, c("module", "role", "analysis_use")],
    module_camera, by = "module", all.x = TRUE, sort = FALSE
  )
  module_camera$contrast <- contrast_id
  module_rows[[contrast_id]] <- module_camera
}

candidate_stats <- do.call(rbind, candidate_rows)
for (contrast_id in contrast_names) {
  ix <- which(
    candidate_stats$contrast == contrast_id &
      candidate_stats$present_after_variance_filter
  )
  if (length(ix) > 0L) {
    candidate_stats$q_candidate[ix] <- stats::p.adjust(
      candidate_stats$p_value[ix], method = "BH"
    )
  }
}
all_gene_stats <- do.call(rbind, all_gene_rows)
module_stats <- do.call(rbind, module_rows)
rownames(candidate_stats) <- NULL
rownames(all_gene_stats) <- NULL
rownames(module_stats) <- NULL

write_csv_gz(
  all_gene_stats,
  file.path(result_dir, paste0(gse_id, "_S2_all_gene_contrast_statistics_v1.csv.gz"))
)
safe_write_csv(
  candidate_stats,
  file.path(result_dir, paste0(gse_id, "_S2_frozen_candidate_statistics_v1.csv"))
)
safe_write_csv(
  module_stats,
  file.path(result_dir, paste0(gse_id, "_S2_module_camera_statistics_v1.csv"))
)

get_module_row <- function(module_name, contrast_name) {
  x <- module_stats[
    module_stats$module == module_name & module_stats$contrast == contrast_name,
    , drop = FALSE
  ]
  if (nrow(x) == 0L) {
    return(data.frame(
      module = module_name, contrast = contrast_name,
      Direction = NA_character_, FDR = NA_real_, stringsAsFactors = FALSE
    ))
  }
  x[1L, c("module", "contrast", "Direction", "FDR"), drop = FALSE]
}

hippo_tead <- get_module_row("hippo_yap_taz", "tead_main_effect")
tgf_tgfb <- get_module_row("tgf_fibrosis", "tgfb_main_effect")
hippo_down <- identical(as.character(hippo_tead$Direction[[1L]]), "Down")
tgf_up <- identical(as.character(tgf_tgfb$Direction[[1L]]), "Up")
axis_directional <- hippo_down && tgf_up
axis_fdr_supported <- axis_directional &&
  is.finite(hippo_tead$FDR[[1L]]) && hippo_tead$FDR[[1L]] <= 0.05 &&
  is.finite(tgf_tgfb$FDR[[1L]]) && tgf_tgfb$FDR[[1L]] <= 0.05
if (axis_fdr_supported) {
  interpretation_class <- "REGULATORY_AXIS_CONCORDANT"
} else if (axis_directional) {
  interpretation_class <- "REGULATORY_AXIS_DIRECTIONALLY_CONCORDANT"
} else if (hippo_down || tgf_up) {
  interpretation_class <- "PARTIAL_DIRECTIONAL_SUPPORT"
} else {
  interpretation_class <- "DISCORDANT_OR_NULL"
}

candidate_coverage_gate <- all(candidate_coverage$present_in_processed_matrix)
core_coverage_gate <- all(
  module_coverage$coverage_fraction[module_coverage$core_module] >= 0.80
)
analysis_gate <- if (
  candidate_coverage_gate && core_coverage_gate && nrow(candidate_stats) > 0L
) {
  "PASS_TO_S2_EVIDENCE_SYNTHESIS"
} else {
  "HOLD_COVERAGE_OR_ANALYSIS_REVIEW"
}

preprocessing_audit <- data.frame(
  source = gse_id,
  matrix_filename = matrix_filename,
  matrix_path = matrix_path,
  raw_rows = nrow(expression_df),
  raw_columns = ncol(expression_df),
  sample_columns = length(sample_columns),
  annotation_columns = paste(setdiff(names(expression_df), sample_columns), collapse = ";"),
  resolved_gene_column = gene_column,
  duplicate_gene_rows_collapsed = duplicate_gene_count,
  rows_after_gene_collapse = nrow(expression_matrix),
  rows_after_variance_filter = nrow(analysis_matrix),
  value_min = value_min,
  value_max = value_max,
  integer_fraction = integer_fraction,
  value_transform = value_transform,
  frozen_candidate_coverage = paste0(
    sum(candidate_coverage$present_in_processed_matrix), "/",
    nrow(candidate_coverage)
  ),
  core_module_coverage_gate = core_coverage_gate,
  stringsAsFactors = FALSE
)
safe_write_csv(
  preprocessing_audit,
  file.path(result_dir, paste0(gse_id, "_S2_preprocessing_audit_v1.csv"))
)

decision_path <- file.path(
  result_dir, paste0(gse_id, "_Step17E_S2_expression_analysis_decision_v1.md")
)
decision_lines <- c(
  paste0("# Step 17E ", gse_id, " S2 frozen-contract expression analysis"),
  "",
  paste0("- Analysis gate: **", analysis_gate, "**."),
  paste0("- Interpretation class: **", interpretation_class, "**."),
  "- Overall project evidence grade remains: **CAUTION**.",
  "- This analysis is regulatory-axis cross-validation, not mechanical causality.",
  "",
  "## Prespecified contrasts",
  "",
  "- `tgfb_main_effect`: marginal TGF-beta exposure effect across TEAD levels.",
  "- `tead_main_effect`: marginal TEAD-inhibition effect across TGF-beta levels.",
  "- `tgfb_x_tead_interaction`: difference in the TGF-beta effect between TEAD states.",
  "",
  "## Axis-focused interpretation",
  "",
  paste0(
    "- Hippo/YAP/TAZ module under TEAD inhibition: ",
    as.character(hippo_tead$Direction[[1L]]),
    "; FDR=", signif(hippo_tead$FDR[[1L]], 4), "."
  ),
  paste0(
    "- TGF/fibrosis module under TGF-beta exposure: ",
    as.character(tgf_tgfb$Direction[[1L]]),
    "; FDR=", signif(tgf_tgfb$FDR[[1L]], 4), "."
  ),
  paste0("- Expected-direction concordance: ", axis_directional, "."),
  paste0("- Expected-direction FDR support: ", axis_fdr_supported, "."),
  "",
  "## Interpretation boundary",
  "",
  "- A positive result supports TGF-beta/SMAD and TEAD-related regulatory-axis cross-validation within this vocal-fold fibroblast dataset.",
  "- It does not establish mechanical loading or stiffness causality, YAP/TAZ-independent mechanical driving, fascia specificity, or pain causality.",
  "- Candidate-level results are exploratory within the frozen 15-gene contract; module-level results are the primary S2 evidence.",
  "- Null or discordant results are retained as boundary evidence and do not justify changing the frozen panel or relaxing thresholds.",
  "",
  "## Material Passport",
  "",
  "- Input: official GSE338388 normalized expression matrix already downloaded during Step 17D3/17D4.",
  "- Provenance/design: corrected official GEO sample design and validated Yap164 suffix mapping from Step 17D4.",
  paste0("- Preprocessing: ", value_transform, "; duplicate gene identifiers collapsed by mean where needed."),
  "- Model: limma linear model with four factorial condition means; no batch term was introduced because the locked 2x2 design has three replicates per cell and no validated independent batch factor.",
  "- Multiplicity: genome-wide BH for all genes, frozen-family BH for candidate genes, and module-level FDR from camera.",
  "- No FASTQ, RAW archive, H5AD, or additional large file was downloaded in Step 17E."
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17E ", gse_id, " S2 frozen-contract expression analysis completed.")
message("Analysis gate: ", analysis_gate)
message("Interpretation class: ", interpretation_class)
message(
  "Frozen candidate coverage: ",
  sum(candidate_coverage$present_in_processed_matrix), "/",
  nrow(candidate_coverage)
)
message("Core module coverage gate: ", core_coverage_gate)
message("Decision: ", decision_path)
