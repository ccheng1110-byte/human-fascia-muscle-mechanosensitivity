options(stringsAsFactors = FALSE)

# Step 11D2: correct the GSE130973 metadata interpretation and audit whether
# any cluster is fibroblast-like. This is descriptive annotation, not a test.

project_dir <- "."
object_path <- file.path(
  project_dir, "data", "raw", "independent_sources", "GSE130973",
  "GSE130973_seurat_analysis_lyko.rds.gz"
)
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11D2_GSE130973_cluster_fibroblast_audit"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(object_path)) {
  stop("Missing GSE130973 object: ", object_path)
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is required for this object-level audit.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required for this object-level audit.")
}

read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection, text = FALSE)
  on.exit(try(close(inner_connection), silent = TRUE), add = TRUE)
  readRDS(inner_connection)
}

extract_metadata <- function(object) {
  if (isS4(object) && "meta.data" %in% methods::slotNames(object)) {
    return(as.data.frame(methods::slot(object, "meta.data")))
  }
  stop("No Seurat meta.data slot was found.")
}

get_rna_data <- function(object) {
  result <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", slot = "data"),
    error = function(e) NULL
  )
  if (!is.null(result)) return(result)
  result <- tryCatch(
    Seurat::GetAssayData(object, assay = "RNA", layer = "data"),
    error = function(e) NULL
  )
  if (is.null(result)) {
    stop("Could not extract the RNA data layer from the Seurat object.")
  }
  result
}

object <- read_double_gzip_rds(object_path)
metadata <- extract_metadata(object)
required_fields <- c("subj", "age", "celltype.age")
missing_fields <- setdiff(required_fields, colnames(metadata))
if (length(missing_fields) > 0L) {
  stop("Missing expected GSE130973 metadata field(s): ", paste(missing_fields, collapse = ", "))
}

metadata$subj <- as.character(metadata$subj)
metadata$age <- as.character(metadata$age)
metadata$celltype.age <- as.character(metadata$celltype.age)
metadata$cluster_id <- sub("_.*$", "", metadata$celltype.age)
metadata$cluster_age_label <- metadata$celltype.age

sample_inventory <- as.data.frame(table(metadata$subj, metadata$age))
names(sample_inventory) <- c("subj", "age", "cells")
sample_inventory <- sample_inventory[sample_inventory$cells > 0L, , drop = FALSE]
write.csv(
  sample_inventory,
  file.path(result_dir, "GSE130973_subject_age_cell_inventory_v2.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cluster_subject_counts <- as.data.frame.matrix(
  table(metadata$cluster_id, metadata$subj)
)
cluster_subject_counts$cluster_id <- rownames(cluster_subject_counts)
cluster_subject_counts <- cluster_subject_counts[
  , c("cluster_id", setdiff(colnames(cluster_subject_counts), "cluster_id")),
  drop = FALSE
]
write.csv(
  cluster_subject_counts,
  file.path(result_dir, "GSE130973_cluster_by_subject_cell_counts_v2.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

fibroblast_markers <- c(
  "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2",
  "COL6A3", "PDGFRA", "FAP", "DPT", "CFD", "COL5A1", "COL5A2",
  "COL14A1", "C7"
)
myofibroblast_context <- c(
  "ACTA2", "TAGLN", "MYL9", "CNN1", "POSTN", "CTHRC1", "LRRC15", "SFRP4"
)
excluded_lineage_markers <- c(
  "PTPRC", "CD3D", "CD3E", "LST1", "LYZ", "KRT14", "KRT5", "EPCAM",
  "PECAM1", "VWF"
)
mechanotransduction_candidates <- c(
  "ITGAV", "ITGB1", "ITGA2", "PARVA", "ITGA5", "FERMT2", "CFL1",
  "LIMK1", "CNN2", "CDC42", "PIEZO2", "TMEM63B", "PANX1", "NF2", "TEAD1"
)

expression <- get_rna_data(object)
common_cells <- intersect(colnames(expression), rownames(metadata))
if (length(common_cells) < 1000L) {
  stop("Too few cells overlap between expression data and metadata.")
}
expression <- expression[, common_cells, drop = FALSE]
metadata <- metadata[common_cells, , drop = FALSE]

all_audit_genes <- unique(c(
  fibroblast_markers, myofibroblast_context,
  excluded_lineage_markers, mechanotransduction_candidates
))
present_genes <- intersect(all_audit_genes, rownames(expression))
missing_genes <- setdiff(all_audit_genes, present_genes)

cluster_levels <- sort(unique(metadata$cluster_id))
cluster_rows <- lapply(cluster_levels, function(cluster) {
  cell_index <- which(metadata$cluster_id == cluster)
  values <- expression[present_genes, cell_index, drop = FALSE]
  gene_means <- if (length(present_genes) > 0L) {
    Matrix::rowMeans(values)
  } else {
    numeric()
  }
  mean_for <- function(genes) {
    genes <- intersect(genes, names(gene_means))
    if (length(genes) == 0L) return(NA_real_)
    mean(gene_means[genes])
  }
  data.frame(
    cluster_id = cluster,
    cells = length(cell_index),
    subjects_present = length(unique(metadata$subj[cell_index])),
    age_labels = paste(unique(metadata$age[cell_index]), collapse = ";"),
    fibroblast_marker_mean = mean_for(fibroblast_markers),
    myofibroblast_context_mean = mean_for(myofibroblast_context),
    excluded_lineage_marker_mean = mean_for(excluded_lineage_markers),
    mechanotransduction_candidate_mean = mean_for(mechanotransduction_candidates),
    fibroblast_markers_present = paste(
      intersect(fibroblast_markers, names(gene_means)), collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
})
cluster_audit <- do.call(rbind, cluster_rows)
cluster_audit$fibroblast_rank <- rank(
  -cluster_audit$fibroblast_marker_mean,
  ties.method = "min", na.last = "keep"
)
cluster_audit <- cluster_audit[order(cluster_audit$fibroblast_rank), , drop = FALSE]
rownames(cluster_audit) <- NULL
write.csv(
  cluster_audit,
  file.path(result_dir, "GSE130973_cluster_fibroblast_marker_audit_v2.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  data.frame(
    gene = all_audit_genes,
    present_in_RNA_data = all_audit_genes %in% present_genes,
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "GSE130973_cluster_audit_gene_inventory_v2.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

decision_lines <- c(
  "# Step 11D2 GSE130973 cluster-level fibroblast audit",
  "",
  "## Corrected metadata interpretation",
  "",
  "- `subj` is the valid sample/subject field and contains five subjects: S1–S5.",
  "- `celltype.age` is a composite cluster-plus-age label, not a confirmed cell-type annotation.",
  "- `integrated_snn_res.0.4` contains the underlying cluster IDs.",
  "",
  "## Audit boundary",
  "",
  "Clusters are ranked by fibroblast-marker expression and subject coverage. This is a descriptive state-identification step; it does not establish that any cluster is a validated fibroblast population and does not perform a gene-level test.",
  "",
  "A cluster may proceed to targeted sample-level validation only after manual inspection confirms fibroblast identity, adequate representation across subjects, and a reproducible rule for including or excluding cells.",
  "",
  "## Material Passport",
  "",
  "- Source: GSE130973 processed Seurat object.",
  "- Transformation: cluster-level marker means and subject coverage were calculated from the existing RNA data layer.",
  "- No disease-group or mechanosensitivity hypothesis test was performed.",
  "- Reproducibility: all cluster-by-subject counts, marker audit values, and gene availability are exported in this directory."
)
decision_path <- file.path(
  result_dir, "GSE130973_cluster_fibroblast_audit_decision_v2.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11D2 GSE130973 cluster fibroblast audit completed.")
message("Subjects: ", paste(sort(unique(metadata$subj)), collapse = ", "))
message("Clusters: ", length(cluster_levels))
message("Audit genes present: ", length(present_genes), "/", length(all_audit_genes))
message("Decision: ", decision_path)
