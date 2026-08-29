options(stringsAsFactors = FALSE)

# Step 11E: exploratory sample-level audit of the frozen mechanotransduction
# program in the marker-defined GSE130973 fibroblast-state candidates.
# This is not a direct F7/F8 replication and cannot upgrade the evidence grade.

project_dir <- "."
object_path <- file.path(
  project_dir, "data", "raw", "independent_sources", "GSE130973",
  "GSE130973_seurat_analysis_lyko.rds.gz"
)
cluster_config_path <- file.path(
  project_dir, "config", "GSE130973_candidate_fibroblast_clusters_v1.csv"
)
registry_path <- file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
)
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11E_GSE130973_frozen_program_audit"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(object_path, cluster_config_path, registry_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 11E input(s): ", paste(missing_inputs, collapse = "; "))
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is required for Step 11E.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required for Step 11E.")
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection, text = FALSE)
  on.exit(try(close(inner_connection), silent = TRUE), add = TRUE)
  readRDS(inner_connection)
}

object <- read_double_gzip_rds(object_path)
if (!isS4(object) || !"meta.data" %in% methods::slotNames(object)) {
  stop("The GSE130973 object has no Seurat meta.data slot.")
}
metadata <- as.data.frame(methods::slot(object, "meta.data"))
required_metadata <- c("subj", "celltype.age")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0L) {
  stop("Missing GSE130973 metadata field(s): ", paste(missing_metadata, collapse = ", "))
}
metadata$subj <- as.character(metadata$subj)
metadata$cluster_id <- sub("_.*$", "", as.character(metadata$celltype.age))

cluster_config <- read.csv(cluster_config_path, check.names = FALSE)
candidate_clusters <- as.character(cluster_config$cluster_id)
if (length(candidate_clusters) != 5L || any(!candidate_clusters %in% metadata$cluster_id)) {
  stop("The frozen candidate fibroblast cluster configuration is invalid.")
}

registry <- read.csv(registry_path, check.names = FALSE)
registry_gene_rows <- lapply(seq_len(nrow(registry)), function(i) {
  data.frame(
    module = registry$module[[i]],
    role = registry$role[[i]],
    gene = trimws(strsplit(registry$genes[[i]], ";", fixed = TRUE)[[1L]]),
    stringsAsFactors = FALSE
  )
})
registry_genes <- do.call(rbind, registry_gene_rows)
panel_genes <- unique(registry_genes$gene)

expression <- tryCatch(
  Seurat::GetAssayData(object, assay = "RNA", slot = "data"),
  error = function(e) Seurat::GetAssayData(object, assay = "RNA", layer = "data")
)
common_cells <- intersect(colnames(expression), rownames(metadata))
if (length(common_cells) < 1000L) {
  stop("Too few cells overlap between expression and metadata.")
}
expression <- expression[, common_cells, drop = FALSE]
metadata <- metadata[common_cells, , drop = FALSE]
available_genes <- intersect(panel_genes, rownames(expression))
missing_genes <- setdiff(panel_genes, available_genes)
candidate_missing <- intersect(
  cluster_config$cluster_id, character()
)

group <- ifelse(
  metadata$cluster_id %in% candidate_clusters,
  "candidate_fibroblast_state",
  "other_clusters"
)

long_rows <- list()
for (sample_id in sort(unique(metadata$subj))) {
  sample_index <- which(metadata$subj == sample_id)
  for (comparison_group in c("candidate_fibroblast_state", "other_clusters")) {
    group_index <- sample_index[group[sample_index] == comparison_group]
    if (length(group_index) == 0L) next
    values <- expression[available_genes, group_index, drop = FALSE]
    for (gene in available_genes) {
      gene_values <- values[gene, ]
      long_rows[[length(long_rows) + 1L]] <- data.frame(
        sample_id = sample_id,
        comparison_group = comparison_group,
        gene = gene,
        cells = length(group_index),
        mean_expression = mean(gene_values),
        percent_cells_expressing = 100 * mean(gene_values > 0),
        stringsAsFactors = FALSE
      )
    }
  }
}
expression_summary <- do.call(rbind, long_rows)
safe_write_csv(
  expression_summary,
  file.path(result_dir, "GSE130973_fibroblast_state_gene_summary_v1.csv")
)

candidate_rows <- expression_summary[
  expression_summary$comparison_group == "candidate_fibroblast_state", , drop = FALSE
]
other_rows <- expression_summary[
  expression_summary$comparison_group == "other_clusters", , drop = FALSE
]
paired <- merge(
  candidate_rows, other_rows,
  by = c("sample_id", "gene"),
  suffixes = c("_candidate", "_other"),
  all = TRUE
)
paired$mean_difference <-
  paired$mean_expression_candidate - paired$mean_expression_other
paired$prevalence_difference <-
  paired$percent_cells_expressing_candidate -
  paired$percent_cells_expressing_other
safe_write_csv(
  paired,
  file.path(result_dir, "GSE130973_fibroblast_state_gene_contrasts_v1.csv")
)

module_rows <- lapply(unique(registry$module), function(module_name) {
  module_genes <- intersect(
    registry_genes$gene[registry_genes$module == module_name],
    available_genes
  )
  x <- paired[paired$gene %in% module_genes, , drop = FALSE]
  module_difference <- aggregate(
    mean_difference ~ sample_id, x, mean
  )
  values <- module_difference$mean_difference
  p_value <- if (sum(values != 0) > 0L) {
    suppressWarnings(stats::wilcox.test(
      values, alternative = "greater", exact = FALSE
    )$p.value)
  } else {
    1
  }
  data.frame(
    module = module_name,
    role = registry$role[match(module_name, registry$module)],
    genes_frozen = length(unique(registry_genes$gene[registry_genes$module == module_name])),
    genes_present = length(module_genes),
    samples = length(values),
    positive_samples = sum(values > 0),
    negative_samples = sum(values < 0),
    median_difference = median(values),
    IQR_difference = IQR(values),
    descriptive_one_sided_p = p_value,
    stringsAsFactors = FALSE
  )
})
module_summary <- do.call(rbind, module_rows)
safe_write_csv(
  module_summary,
  file.path(result_dir, "GSE130973_fibroblast_state_module_summary_v1.csv")
)

decision_lines <- c(
  "# Step 11E GSE130973 exploratory frozen-program audit",
  "",
  "## Design",
  "",
  "- Statistical units: five subjects S1–S5.",
  paste0("- Candidate fibroblast-state clusters: ", paste(candidate_clusters, collapse = ", "), "."),
  "- Comparator: all other clusters within the same subject.",
  "- Candidate clusters were selected by the pre-defined top-five fibroblast-marker ranking from Step 11D2.",
  "",
  "## Boundary",
  "",
  "This is an exploratory state-level audit in an aging-skin dataset. It is not a direct F7/F8 replication, does not test disease status, and cannot upgrade the current evidence grade. Descriptive p values are not confirmatory with n = 5.",
  "",
  paste0("- Frozen registry genes present: ", length(available_genes), "/", length(panel_genes), "."),
  paste0("- Missing registry genes: ", paste(missing_genes, collapse = ", "), "."),
  "",
  "## Material Passport",
  "",
  "- Source: GSE130973 processed Seurat object.",
  "- Transformation: within-subject candidate-state versus other-cluster gene and module summaries.",
  "- No causal model or disease-group test was fitted.",
  "- Reproducibility: candidate cluster configuration and all per-subject contrasts are exported in this result directory."
)
decision_path <- file.path(
  result_dir, "GSE130973_fibroblast_state_audit_decision_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11E exploratory frozen-program audit completed.")
message("Candidate clusters: ", paste(candidate_clusters, collapse = ", "))
message("Subjects: ", paste(sort(unique(metadata$subj)), collapse = ", "))
message("Registry genes present: ", length(available_genes), "/", length(panel_genes))
message("Decision: ", decision_path)
