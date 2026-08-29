options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.")
  }
})

project_dir <- "."
input_path <- file.path(
  project_dir, "data", "raw", "GSE173252",
  "GSE173252_dd_mesenchyme.rds.gz"
)
result_dir <- file.path(
  project_dir, "results", "03_cluster_reproduction", "GSE173252"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection, text = FALSE)
  on.exit(try(close(inner_connection), silent = TRUE), add = TRUE)
  readRDS(inner_connection)
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

message("Loading legacy processed object: ", input_path)
object <- read_double_gzip_rds(input_path)
object_attributes <- attributes(object)
metadata <- object_attributes[["meta.data"]]
rna_attributes <- attributes(object_attributes[["assays"]][["RNA"]])
rna_counts <- rna_attributes[["counts"]]
rna_logdata <- rna_attributes[["data"]]

if (is.null(metadata) || is.null(rna_counts) || is.null(rna_logdata)) {
  stop("Could not extract metadata, RNA counts, or RNA log-normalized data.")
}
if (!identical(colnames(rna_counts), rownames(metadata))) {
  if (!all(colnames(rna_counts) %in% rownames(metadata))) {
    stop("RNA matrix cell names do not match metadata row names.")
  }
  metadata <- metadata[colnames(rna_counts), , drop = FALSE]
}
if (!identical(rownames(rna_counts), rownames(rna_logdata))) {
  stop("RNA counts and log-normalized matrices have different feature orders.")
}

cluster_col <- "seurat_clusters"
sample_col <- "orig.ident"
required_cols <- c(cluster_col, sample_col, "source", "doublet_prediction", "ordered_clusters")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0L) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

metadata[[cluster_col]] <- as.character(metadata[[cluster_col]])
metadata[[sample_col]] <- as.character(metadata[[sample_col]])
metadata$ordered_clusters <- as.character(metadata$ordered_clusters)
doublet_field <- metadata$doublet_prediction
if (is.logical(doublet_field)) {
  is_predicted_doublet <- doublet_field
} else {
  is_predicted_doublet <- tolower(trimws(as.character(doublet_field))) %in%
    c("true", "doublet", "yes", "1")
}
is_predicted_doublet[is.na(is_predicted_doublet)] <- FALSE
singlet_keep <- !is_predicted_doublet

message(
  "Primary analysis retains ", sum(singlet_keep), " singlets and excludes ",
  sum(!singlet_keep), " predicted doublets."
)

cluster_levels <- sort(unique(metadata[[cluster_col]]))
sample_levels <- unique(metadata[[sample_col]])

marker_sets <- list(
  VSMC = c("MYH11", "RERGL", "CNN1", "BCAM", "ACTA2"),
  Pericyte = c("RGS5", "STEAP4", "MCAM", "PLN"),
  FB = c("PDGFRA", "CD34", "APOD", "PLA2G2A", "FBLN5", "OSR2"),
  Myofib = c(
    "PDPN", "SCX", "ADAM12", "LOXL2", "FAP", "COL1A1", "COL1A2",
    "COL3A1", "CTHRC1", "TNFRSF12A"
  )
)
anchor_markers <- c(VSMC = "MYH11", Pericyte = "RGS5", FB = "PDGFRA", Myofib = "PDPN")
panel_genes <- unique(unlist(marker_sets, use.names = FALSE))
panel_genes <- panel_genes[panel_genes %in% rownames(rna_logdata)]

marker_expression_rows <- list()
row_index <- 1L
for (cluster_id in cluster_levels) {
  cell_index <- which(singlet_keep & metadata[[cluster_col]] == cluster_id)
  means <- Matrix::rowMeans(rna_logdata[panel_genes, cell_index, drop = FALSE])
  percentages <- Matrix::rowMeans(rna_counts[panel_genes, cell_index, drop = FALSE] > 0) * 100
  marker_expression_rows[[row_index]] <- data.frame(
    seurat_cluster = cluster_id,
    gene = panel_genes,
    mean_log_normalized_expression = as.numeric(means),
    percent_cells_expressing = as.numeric(percentages),
    cells_in_cluster = length(cell_index),
    stringsAsFactors = FALSE
  )
  row_index <- row_index + 1L
}
marker_expression <- do.call(rbind, marker_expression_rows)
marker_expression$paper_marker_class <- unname(vapply(
  marker_expression$gene,
  function(gene) paste(names(marker_sets)[vapply(marker_sets, function(x) gene %in% x, logical(1))], collapse = ";"),
  character(1)
))
safe_write_csv(
  marker_expression,
  file.path(result_dir, "GSE173252_mesenchyme_paper_marker_expression_by_cluster.csv")
)

anchor_max_cluster <- vapply(names(anchor_markers), function(identity_name) {
  gene <- anchor_markers[[identity_name]]
  values <- marker_expression[marker_expression$gene == gene, , drop = FALSE]
  values$seurat_cluster[which.max(values$mean_log_normalized_expression)]
}, character(1))

ordered_cross_tab <- table(
  seurat_cluster = metadata[[cluster_col]][singlet_keep],
  ordered_cluster = metadata$ordered_clusters[singlet_keep]
)
dominant_ordered_cluster <- apply(ordered_cross_tab, 1, function(x) names(which.max(x)))
paper_cluster_number <- as.integer(dominant_ordered_cluster) + 1L
paper_identity_by_number <- c("1" = "VSMC", "2" = "Pericyte", "3" = "FB", "4" = "Myofib")
identity_from_order <- unname(paper_identity_by_number[as.character(paper_cluster_number)])

identity_map <- data.frame(
  seurat_cluster = names(dominant_ordered_cluster),
  ordered_cluster = unname(dominant_ordered_cluster),
  paper_cluster = paper_cluster_number,
  reproduced_identity = identity_from_order,
  anchor_marker = unname(anchor_markers[identity_from_order]),
  anchor_maximum_cluster = unname(anchor_max_cluster[identity_from_order]),
  anchor_agrees = unname(anchor_max_cluster[identity_from_order]) == names(dominant_ordered_cluster),
  singlet_cells = as.integer(table(metadata[[cluster_col]][singlet_keep])[names(dominant_ordered_cluster)]),
  stringsAsFactors = FALSE
)
if (any(is.na(identity_map$reproduced_identity)) || any(!identity_map$anchor_agrees)) {
  stop("Identity reproduction failed: ordered-cluster and anchor-marker evidence disagree.")
}
safe_write_csv(
  identity_map,
  file.path(result_dir, "GSE173252_mesenchyme_cluster_identity_map.csv")
)

# Aggregate raw counts by sample x cluster. This preserves biological sample identity
# and avoids treating individual cells as independent replicates.
singlet_metadata <- metadata[singlet_keep, , drop = FALSE]
group_factor <- interaction(
  singlet_metadata[[sample_col]], singlet_metadata[[cluster_col]],
  drop = TRUE, sep = "__"
)
group_levels <- levels(group_factor)
design <- Matrix::sparseMatrix(
  i = seq_along(group_factor),
  j = as.integer(group_factor),
  x = 1,
  dims = c(length(group_factor), length(group_levels)),
  dimnames = list(rownames(singlet_metadata), group_levels)
)
pseudobulk_counts <- rna_counts[, singlet_keep, drop = FALSE] %*% design
group_cell_counts <- as.integer(table(group_factor)[group_levels])
group_parts <- do.call(rbind, strsplit(group_levels, "__", fixed = TRUE))
group_inventory <- data.frame(
  pseudobulk_group = group_levels,
  sample = group_parts[, 1],
  seurat_cluster = group_parts[, 2],
  reproduced_identity = identity_map$reproduced_identity[
    match(group_parts[, 2], identity_map$seurat_cluster)
  ],
  singlet_cells = group_cell_counts,
  library_size = as.numeric(Matrix::colSums(pseudobulk_counts)),
  stringsAsFactors = FALSE
)
safe_write_csv(
  group_inventory,
  file.path(result_dir, "GSE173252_mesenchyme_pseudobulk_group_cell_counts.csv")
)

# Exploratory sample-aware evidence: within each sample, compare one cluster with
# all other mesenchymal clusters. This is not a substitute for formal edgeR/DESeq2 DE.
all_clusters <- cluster_levels
all_samples <- unique(singlet_metadata[[sample_col]])
min_target_cells <- 20L
min_rest_cells <- 50L
pseudobulk_results <- list()
result_index <- 1L

for (cluster_id in all_clusters) {
  sample_lfc <- list()
  used_samples <- character(0)

  for (sample_id in all_samples) {
    target_name <- paste(sample_id, cluster_id, sep = "__")
    target_index <- match(target_name, group_levels)
    rest_index <- which(
      group_inventory$sample == sample_id &
        group_inventory$seurat_cluster != cluster_id
    )
    target_cells <- if (is.na(target_index)) 0L else group_inventory$singlet_cells[target_index]
    rest_cells <- sum(group_inventory$singlet_cells[rest_index])

    if (target_cells < min_target_cells || rest_cells < min_rest_cells) next

    target_counts <- as.numeric(pseudobulk_counts[, target_index])
    rest_counts <- as.numeric(Matrix::rowSums(pseudobulk_counts[, rest_index, drop = FALSE]))
    target_log_cpm <- log2(target_counts / sum(target_counts) * 1e6 + 0.5)
    rest_log_cpm <- log2(rest_counts / sum(rest_counts) * 1e6 + 0.5)
    sample_lfc[[sample_id]] <- target_log_cpm - rest_log_cpm
    used_samples <- c(used_samples, sample_id)
  }

  if (length(sample_lfc) == 0L) next
  lfc_matrix <- do.call(cbind, sample_lfc)
  rownames(lfc_matrix) <- rownames(rna_counts)

  cluster_cells <- which(singlet_keep & metadata[[cluster_col]] == cluster_id)
  pct_expr <- Matrix::rowMeans(rna_counts[, cluster_cells, drop = FALSE] > 0) * 100
  mean_expr <- Matrix::rowMeans(rna_logdata[, cluster_cells, drop = FALSE])

  median_lfc <- apply(lfc_matrix, 1, median, na.rm = TRUE)
  mean_lfc <- rowMeans(lfc_matrix, na.rm = TRUE)
  min_lfc <- apply(lfc_matrix, 1, min, na.rm = TRUE)
  positive_fraction <- rowMeans(lfc_matrix > 0, na.rm = TRUE)
  strong_positive_fraction <- rowMeans(lfc_matrix > 0.5, na.rm = TRUE)

  p_value <- rep(NA_real_, nrow(lfc_matrix))
  if (ncol(lfc_matrix) >= 3L) {
    testable <- which(as.numeric(pct_expr) >= 5 & apply(lfc_matrix, 1, function(x) any(x != 0)))
    p_value[testable] <- vapply(testable, function(i) {
      tryCatch(
        suppressWarnings(wilcox.test(lfc_matrix[i, ], mu = 0, exact = FALSE)$p.value),
        error = function(e) NA_real_
      )
    }, numeric(1))
  }
  fdr <- p.adjust(p_value, method = "BH")

  cluster_result <- data.frame(
    seurat_cluster = cluster_id,
    reproduced_identity = identity_map$reproduced_identity[
      match(cluster_id, identity_map$seurat_cluster)
    ],
    gene = rownames(lfc_matrix),
    samples_evaluable = ncol(lfc_matrix),
    sample_ids = paste(used_samples, collapse = ";"),
    median_log2_cpm_difference = as.numeric(median_lfc),
    mean_log2_cpm_difference = as.numeric(mean_lfc),
    minimum_log2_cpm_difference = as.numeric(min_lfc),
    positive_sample_fraction = as.numeric(positive_fraction),
    strong_positive_sample_fraction = as.numeric(strong_positive_fraction),
    percent_cluster_cells_expressing = as.numeric(pct_expr),
    mean_cluster_log_normalized_expression = as.numeric(mean_expr),
    exploratory_wilcoxon_p = p_value,
    exploratory_bh_fdr = fdr,
    stringsAsFactors = FALSE
  )
  pseudobulk_results[[result_index]] <- cluster_result
  result_index <- result_index + 1L
}

pseudobulk_evidence <- do.call(rbind, pseudobulk_results)
candidate_markers <- pseudobulk_evidence[
  pseudobulk_evidence$percent_cluster_cells_expressing >= 10 &
    pseudobulk_evidence$median_log2_cpm_difference >= 0.5 &
    pseudobulk_evidence$positive_sample_fraction >= 0.75 &
    pseudobulk_evidence$samples_evaluable >= 3,
  , drop = FALSE
]
candidate_markers <- candidate_markers[order(
  candidate_markers$seurat_cluster,
  -candidate_markers$positive_sample_fraction,
  -candidate_markers$median_log2_cpm_difference,
  -candidate_markers$percent_cluster_cells_expressing
), , drop = FALSE]
top100 <- do.call(rbind, lapply(split(candidate_markers, candidate_markers$seurat_cluster), head, 100))
rownames(top100) <- NULL
safe_write_csv(
  top100,
  file.path(result_dir, "GSE173252_mesenchyme_pseudobulk_top100_markers_by_cluster.csv")
)

paper_marker_pseudobulk <- pseudobulk_evidence[
  pseudobulk_evidence$gene %in% unique(unlist(marker_sets, use.names = FALSE)),
  , drop = FALSE
]
paper_marker_pseudobulk$expected_marker_class <- unname(vapply(
  paper_marker_pseudobulk$gene,
  function(gene) paste(names(marker_sets)[vapply(marker_sets, function(x) gene %in% x, logical(1))], collapse = ";"),
  character(1)
))
paper_marker_pseudobulk <- paper_marker_pseudobulk[order(
  paper_marker_pseudobulk$expected_marker_class,
  paper_marker_pseudobulk$gene,
  paper_marker_pseudobulk$seurat_cluster
), , drop = FALSE]
safe_write_csv(
  paper_marker_pseudobulk,
  file.path(result_dir, "GSE173252_mesenchyme_paper_marker_pseudobulk_evidence.csv")
)

coverage_lines <- vapply(cluster_levels, function(cluster_id) {
  row <- pseudobulk_evidence[pseudobulk_evidence$seurat_cluster == cluster_id, , drop = FALSE][1, ]
  paste0(
    "- Cluster ", cluster_id, " (", row$reproduced_identity, "): ",
    row$samples_evaluable, " evaluable sample(s): ", row$sample_ids
  )
}, character(1))

mapping_lines <- apply(identity_map, 1, function(row) {
  paste0(
    "- Seurat cluster ", row[["seurat_cluster"]], " -> paper cluster ",
    row[["paper_cluster"]], " -> ", row[["reproduced_identity"]],
    " (anchor ", row[["anchor_marker"]], "; singlet cells = ",
    row[["singlet_cells"]], ")"
  )
})

decision_lines <- c(
  "GSE173252 mesenchyme cluster reproduction",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Conclusion",
  "The four published coarse mesenchymal identities were reproduced without rerunning clustering.",
  "The mapping is supported independently by the stored ordered_clusters field and canonical/published anchor-marker expression.",
  "",
  "Identity mapping",
  mapping_lines,
  "",
  "Primary analysis set",
  paste0("- Singlets retained: ", sum(singlet_keep)),
  paste0("- Predicted doublets excluded: ", sum(!singlet_keep)),
  "",
  "Sample-aware pseudobulk coverage",
  paste0("A target cluster required >= ", min_target_cells, " cells per sample and the remaining clusters >= ", min_rest_cells, " cells."),
  coverage_lines,
  "",
  "Interpretation limits",
  "- Pseudobulk log2-CPM differences compare each cluster with the other mesenchymal clusters within the same donor/sample.",
  "- The Wilcoxon p values and BH FDR values are exploratory because only nine samples are available and the design is unbalanced.",
  "- These outputs prioritize robust candidate markers; they are not a replacement for a formal edgeR/DESeq2 differential-expression model.",
  "- Cluster abundance differs strongly by disease state, especially for the Myofib population; marker specificity and disease association must remain separate claims.",
  "",
  "Source",
  "Original study DOI: 10.1016/j.jid.2021.05.030"
)
writeLines(
  decision_lines,
  file.path(result_dir, "GSE173252_mesenchyme_cluster_reproduction_decision.txt"),
  useBytes = TRUE
)

message("Cluster reproduction and exploratory pseudobulk marker analysis completed.")
message("Identity map: ", file.path(result_dir, "GSE173252_mesenchyme_cluster_identity_map.csv"))
message("Top markers: ", file.path(result_dir, "GSE173252_mesenchyme_pseudobulk_top100_markers_by_cluster.csv"))
message("Decision report: ", file.path(result_dir, "GSE173252_mesenchyme_cluster_reproduction_decision.txt"))
