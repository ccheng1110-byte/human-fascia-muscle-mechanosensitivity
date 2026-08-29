# Sample-level QC audit for the processed GSE173252 mesenchymal-cell metadata.
# The inferential unit is the GSM/sample, never an individual cell.
# This script performs descriptive auditing only and does not filter cells.

options(stringsAsFactors = FALSE)

project_dir <- "."
metadata_path <- file.path(
  project_dir,
  "data", "metadata", "GSE173252",
  "GSE173252_dd_mesenchyme_cell_metadata.csv.gz"
)
manifest_path <- file.path(
  project_dir,
  "data", "metadata", "GSE173252",
  "GSE173252_sample_manifest_audited.csv"
)
result_dir <- file.path(project_dir, "results", "02_sample_qc", "GSE173252")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(metadata_path, manifest_path)) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path)
  }
}

read_gzip_csv <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  utils::read.csv(
    connection,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

safe_quantile <- function(x, probability) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    unname(stats::quantile(x, probs = probability, na.rm = TRUE, names = FALSE))
  }
}

cell_metadata <- read_gzip_csv(metadata_path)
file_manifest <- utils::read.csv(
  manifest_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "cell_barcode", "orig.ident", "nCount_RNA", "nFeature_RNA", "source",
  "condition", "percent.mito", "doublet_score", "doublet_prediction",
  "seurat_clusters"
)
missing_columns <- setdiff(required_columns, names(cell_metadata))
if (length(missing_columns) > 0L) {
  stop("Missing cell metadata columns: ", paste(missing_columns, collapse = ", "))
}

sample_map <- data.frame(
  orig.ident = c(
    "Dermis1", "Dermis2", "Dermis3",
    "SF1", "SF2", "SF3",
    "DD1", "DD2", "DD3"
  ),
  gsm_id = c(
    "GSM5264341", "GSM5264342", "GSM5264343",
    "GSM5264344", "GSM5264345", "GSM5264346",
    "GSM5264347", "GSM5264348", "GSM5264349"
  ),
  sample_group = c(
    rep("healthy_dermis", 3L),
    rep("nonpathogenic_skoog_fascia", 3L),
    rep("dupuytren_disease", 3L)
  ),
  expected_source = c(
    rep("Dermis", 3L),
    rep("SF", 3L),
    rep("DD", 3L)
  ),
  expected_condition = c(
    rep("Uninjured", 6L),
    rep("DD", 3L)
  ),
  reported_replicate = rep(1:3, 3L),
  stringsAsFactors = FALSE
)

if (!setequal(unique(cell_metadata$orig.ident), sample_map$orig.ident)) {
  stop(
    "orig.ident values do not match the expected nine samples. Observed: ",
    paste(sort(unique(cell_metadata$orig.ident)), collapse = ", ")
  )
}

if (!setequal(sample_map$gsm_id, file_manifest$gsm_id)) {
  stop("Sample-to-GSM mapping does not match the audited file manifest.")
}

cell_metadata <- merge(
  cell_metadata,
  sample_map,
  by = "orig.ident",
  all.x = TRUE,
  sort = FALSE
)

cell_metadata$is_doublet <- as.logical(cell_metadata$doublet_prediction)
if (anyNA(cell_metadata$is_doublet)) {
  stop("doublet_prediction contains values that cannot be interpreted as TRUE/FALSE.")
}

sample_ids <- sample_map$orig.ident
sample_rows <- vector("list", length(sample_ids))

for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[[i]]
  sample_data <- cell_metadata[cell_metadata$orig.ident == sample_id, , drop = FALSE]
  mapping_row <- sample_map[sample_map$orig.ident == sample_id, , drop = FALSE]

  observed_sources <- unique(sample_data$source)
  observed_conditions <- unique(sample_data$condition)

  sample_rows[[i]] <- data.frame(
    orig.ident = sample_id,
    gsm_id = mapping_row$gsm_id,
    sample_group = mapping_row$sample_group,
    reported_replicate = mapping_row$reported_replicate,
    observed_source = paste(observed_sources, collapse = " | "),
    expected_source = mapping_row$expected_source,
    source_label_match = length(observed_sources) == 1L &&
      identical(observed_sources, mapping_row$expected_source),
    observed_condition = paste(observed_conditions, collapse = " | "),
    expected_condition = mapping_row$expected_condition,
    condition_label_match = length(observed_conditions) == 1L &&
      identical(observed_conditions, mapping_row$expected_condition),
    n_cells = nrow(sample_data),
    n_singlets = sum(!sample_data$is_doublet),
    n_predicted_doublets = sum(sample_data$is_doublet),
    predicted_doublet_rate_percent = 100 * mean(sample_data$is_doublet),
    median_nCount_RNA = safe_median(sample_data$nCount_RNA),
    q25_nCount_RNA = safe_quantile(sample_data$nCount_RNA, 0.25),
    q75_nCount_RNA = safe_quantile(sample_data$nCount_RNA, 0.75),
    median_nFeature_RNA = safe_median(sample_data$nFeature_RNA),
    q25_nFeature_RNA = safe_quantile(sample_data$nFeature_RNA, 0.25),
    q75_nFeature_RNA = safe_quantile(sample_data$nFeature_RNA, 0.75),
    median_percent_mito = safe_median(sample_data$percent.mito),
    q25_percent_mito = safe_quantile(sample_data$percent.mito, 0.25),
    q75_percent_mito = safe_quantile(sample_data$percent.mito, 0.75),
    maximum_percent_mito = max(sample_data$percent.mito, na.rm = TRUE),
    percent_cells_mito_gt_5 = 100 * mean(sample_data$percent.mito > 5, na.rm = TRUE),
    percent_cells_mito_ge_10 = 100 * mean(sample_data$percent.mito >= 10, na.rm = TRUE),
    median_doublet_score = safe_median(sample_data$doublet_score),
    low_cell_count_flag = nrow(sample_data) < 500L,
    stringsAsFactors = FALSE
  )
}

sample_qc <- do.call(rbind, sample_rows)
sample_qc <- sample_qc[match(sample_map$orig.ident, sample_qc$orig.ident), , drop = FALSE]

cluster_counts_matrix <- xtabs(
  ~ orig.ident + seurat_clusters,
  data = cell_metadata
)
cluster_counts_matrix <- cluster_counts_matrix[sample_map$orig.ident, , drop = FALSE]

cluster_counts <- as.data.frame.matrix(cluster_counts_matrix)
cluster_counts <- data.frame(
  orig.ident = rownames(cluster_counts),
  cluster_counts,
  check.names = FALSE,
  row.names = NULL
)
names(cluster_counts)[-1L] <- paste0("cluster_", names(cluster_counts)[-1L], "_n")
cluster_counts <- merge(sample_map, cluster_counts, by = "orig.ident", sort = FALSE)
cluster_counts <- cluster_counts[match(sample_map$orig.ident, cluster_counts$orig.ident), , drop = FALSE]

cluster_proportions_matrix <- 100 * prop.table(cluster_counts_matrix, margin = 1L)
cluster_proportions <- as.data.frame.matrix(cluster_proportions_matrix)
cluster_proportions <- data.frame(
  orig.ident = rownames(cluster_proportions),
  cluster_proportions,
  check.names = FALSE,
  row.names = NULL
)
names(cluster_proportions)[-1L] <- paste0(
  "cluster_", names(cluster_proportions)[-1L], "_percent"
)
cluster_proportions <- merge(
  sample_map,
  cluster_proportions,
  by = "orig.ident",
  sort = FALSE
)
cluster_proportions <- cluster_proportions[
  match(sample_map$orig.ident, cluster_proportions$orig.ident),
  ,
  drop = FALSE
]

group_names <- unique(sample_map$sample_group)
group_rows <- lapply(group_names, function(group_name) {
  group_samples <- sample_qc[
    sample_qc$sample_group == group_name,
    ,
    drop = FALSE
  ]
  data.frame(
    sample_group = group_name,
    n_samples = nrow(group_samples),
    total_cells = sum(group_samples$n_cells),
    minimum_cells_per_sample = min(group_samples$n_cells),
    median_cells_per_sample = stats::median(group_samples$n_cells),
    maximum_cells_per_sample = max(group_samples$n_cells),
    mean_of_sample_median_nCount_RNA = mean(group_samples$median_nCount_RNA),
    mean_of_sample_median_nFeature_RNA = mean(group_samples$median_nFeature_RNA),
    mean_of_sample_median_percent_mito = mean(group_samples$median_percent_mito),
    mean_sample_predicted_doublet_rate_percent = mean(
      group_samples$predicted_doublet_rate_percent
    ),
    stringsAsFactors = FALSE
  )
})
group_summary <- do.call(rbind, group_rows)

utils::write.csv(
  sample_qc,
  file.path(result_dir, "GSE173252_mesenchyme_sample_qc.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  cluster_counts,
  file.path(result_dir, "GSE173252_mesenchyme_cluster_counts_by_sample.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  cluster_proportions,
  file.path(result_dir, "GSE173252_mesenchyme_cluster_proportions_by_sample.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  group_summary,
  file.path(result_dir, "GSE173252_mesenchyme_group_qc_summary.csv"),
  row.names = FALSE,
  na = ""
)

cell_count_ratio <- max(sample_qc$n_cells) / min(sample_qc$n_cells)
low_count_samples <- sample_qc$orig.ident[sample_qc$low_cell_count_flag]
all_labels_match <- all(sample_qc$source_label_match & sample_qc$condition_label_match)

report_lines <- c(
  "GSE173252 mesenchymal-cell sample-level QC audit",
  paste0("Audit date: ", format(Sys.Date(), "%Y-%m-%d")),
  paste0("Cells audited: ", nrow(cell_metadata)),
  paste0("Samples audited: ", nrow(sample_qc)),
  paste0("Sample groups: ", paste(group_names, collapse = ", ")),
  paste0("All source and condition labels match the sample map: ", all_labels_match),
  paste0("Predicted doublets: ", sum(cell_metadata$is_doublet)),
  paste0(
    "Overall predicted doublet rate (%): ",
    round(100 * mean(cell_metadata$is_doublet), 3)
  ),
  paste0(
    "Cell-count range per sample: ", min(sample_qc$n_cells), " to ",
    max(sample_qc$n_cells)
  ),
  paste0("Maximum/minimum sample cell-count ratio: ", round(cell_count_ratio, 2)),
  paste0(
    "Samples with fewer than 500 cells: ",
    if (length(low_count_samples) == 0L) "none" else paste(low_count_samples, collapse = ", ")
  ),
  paste0(
    "Cells with mitochondrial percentage >=10%: ",
    sum(cell_metadata$percent.mito >= 10, na.rm = TRUE)
  ),
  "",
  "Interpretation and decision:",
  "- All nine sample identities are preserved and can serve as sample-level units.",
  "- Cell yields are strongly imbalanced across samples; cell-level tests must not be used as donor-level evidence.",
  "- Predicted doublets remain annotated and should be excluded or sensitivity-tested before downstream analysis.",
  "- The object is already processed; do not impose new QC thresholds before reproducing the published annotations.",
  "- Next, identify the four mesenchymal clusters using original markers and sample-aware pseudobulk evidence.",
  "- Do not start CellChat, spatial analysis, or GWAS yet."
)

writeLines(
  report_lines,
  file.path(result_dir, "GSE173252_mesenchyme_sample_qc_decision.txt")
)

print(sample_qc[, c(
  "orig.ident", "gsm_id", "sample_group", "n_cells",
  "predicted_doublet_rate_percent", "median_nCount_RNA",
  "median_nFeature_RNA", "median_percent_mito", "low_cell_count_flag"
)])

if (!all_labels_match) {
  stop("Sample label audit failed: source or condition mismatch detected.")
}

message("\nGSE173252 sample-level QC audit completed successfully.")
message("Sample QC: ", file.path(result_dir, "GSE173252_mesenchyme_sample_qc.csv"))
message("Cluster counts: ", file.path(result_dir, "GSE173252_mesenchyme_cluster_counts_by_sample.csv"))
message("Cluster proportions: ", file.path(result_dir, "GSE173252_mesenchyme_cluster_proportions_by_sample.csv"))
message("Group summary: ", file.path(result_dir, "GSE173252_mesenchyme_group_qc_summary.csv"))
message("Decision: ", file.path(result_dir, "GSE173252_mesenchyme_sample_qc_decision.txt"))
