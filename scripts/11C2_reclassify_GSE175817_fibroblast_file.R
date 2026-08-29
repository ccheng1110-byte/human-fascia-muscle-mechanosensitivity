options(stringsAsFactors = FALSE)

# Step 11C2: correct the classification of the GSE175817 fibroblast file.
# The file name contains "meta", but its structure must be audited before
# treating it as cell metadata.

project_dir <- "."
metadata_path <- file.path(
  project_dir, "data", "metadata", "independent_sources", "GSE175817",
  "GSE175817_meta_fibroblast.csv.gz"
)
filelist_path <- file.path(
  project_dir, "data", "metadata", "independent_sources", "GSE175817",
  "GSE175817_filelist.txt"
)
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11C2_GSE175817_matrix_audit"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(metadata_path)) {
  stop("Missing downloaded GSE175817 file: ", metadata_path)
}

expression_like <- read.csv(
  gzfile(metadata_path, open = "rt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(expression_like) < 1000L || ncol(expression_like) < 2L) {
  stop("The file does not have the expected gene-by-sample matrix dimensions.")
}

first_column <- as.character(expression_like[[1L]])
numeric_columns <- vapply(
  expression_like[-1L], is.numeric, logical(1)
)
gene_like_ids <- mean(
  nzchar(first_column) & !grepl("^-?[0-9]+([.]?[0-9]*)?$", first_column)
)
is_gene_by_sample_matrix <-
  (is.na(names(expression_like)[[1L]]) || names(expression_like)[[1L]] == "") &&
  all(numeric_columns) && gene_like_ids > 0.8

if (!is_gene_by_sample_matrix) {
  stop(
    "The file was not confirmed as a gene-by-sample matrix; inspect its structure manually."
  )
}

sample_columns <- names(expression_like)[-1L]
gene_ids <- first_column
matrix_summary <- data.frame(
  source_accession = "GSE175817",
  file = basename(metadata_path),
  file_class = "gene_by_sample_expression_matrix",
  genes = nrow(expression_like),
  sample_columns = length(sample_columns),
  sample_column_names = paste(sample_columns, collapse = ";"),
  unique_gene_ids = length(unique(gene_ids)),
  duplicated_gene_ids = sum(duplicated(gene_ids)),
  donor_mapping_resolved = FALSE,
  cell_level_metadata_present = FALSE,
  cell_type_labels_present = FALSE,
  stringsAsFactors = FALSE
)
write.csv(
  matrix_summary,
  file.path(result_dir, "GSE175817_fibroblast_matrix_structure_audit_v2.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

sample_mapping <- data.frame(
  source_accession = "GSE175817",
  matrix_sample_column = sample_columns,
  mapped_GSM = NA_character_,
  mapped_donor = NA_character_,
  mapping_status = "unresolved_numeric_column_name",
  stringsAsFactors = FALSE
)
write.csv(
  sample_mapping,
  file.path(result_dir, "GSE175817_fibroblast_matrix_sample_mapping_v2.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (file.exists(filelist_path)) {
  filelist <- read.delim(
    filelist_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#"
  )
  filelist$direct_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175817/suppl/",
    filelist$Name
  )
  write.csv(
    filelist,
    file.path(result_dir, "GSE175817_supplementary_filelist_with_urls_v2.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

decision_lines <- c(
  "# Step 11C2 corrected GSE175817 file classification",
  "",
  "## Finding",
  "",
  "The file `GSE175817_meta_fibroblast.csv.gz` is a gene-by-sample expression matrix, not cell-level metadata. It contains gene identifiers in the first column and five numeric sample columns named 0–4.",
  "",
  "## Current audit result",
  "",
  paste0("- Genes: ", nrow(expression_like), "."),
  paste0("- Numeric sample columns: ", length(sample_columns), "."),
  "- Donor/GSM mapping: unresolved.",
  "- Cell-level metadata: absent from this file.",
  "- Cell-state labels: absent from this file.",
  "",
  "## Decision",
  "",
  "GSE175817 remains a promising independent source, but this particular processed file cannot yet support sample-level validation. The next admissible action is to resolve the numeric matrix columns against GSM/donor identifiers or audit the donor-specific CSV files listed in the official filelist. Do not treat the five matrix columns as five independent donors until the mapping is proven.",
  "",
  "## Material Passport",
  "",
  "- Source: GSE175817, non-overlapping candidate from Step 10A.",
  "- Transformation: structural reclassification only; no gene-level test was performed.",
  "- Integrity boundary: numeric columns 0–4 are unresolved sample labels, not confirmed donor IDs."
)
decision_path <- file.path(
  result_dir, "GSE175817_fibroblast_matrix_audit_decision_v2.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11C2 corrected GSE175817 file classification completed.")
message("File class: gene-by-sample expression matrix")
message("Genes: ", nrow(expression_like), "; sample columns: ", length(sample_columns))
message("Donor mapping: unresolved")
message("Decision: ", decision_path)
