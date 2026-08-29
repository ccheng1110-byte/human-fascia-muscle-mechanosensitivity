options(stringsAsFactors = FALSE)

# Step 11D: download and audit the processed GSE130973 Seurat object.
# The object is a candidate independent source with five GEO subjects.
# This step performs structure/metadata audit only; no hypothesis test.

project_dir <- "."
raw_dir <- file.path(
  project_dir, "data", "raw", "independent_sources", "GSE130973"
)
metadata_dir <- file.path(
  project_dir, "data", "metadata", "independent_sources", "GSE130973"
)
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11D_GSE130973_seurat_audit"
)
for (x in c(raw_dir, metadata_dir, result_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

object_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130973/suppl/",
  "GSE130973_seurat_analysis_lyko.rds.gz"
)
object_path <- file.path(raw_dir, "GSE130973_seurat_analysis_lyko.rds.gz")
candidate_panel_path <- file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
)

download_if_needed <- function(url, path) {
  if (file.exists(path) && file.info(path)$size > 0) {
    message("Using existing file: ", path)
    return(invisible(path))
  }
  message("Downloading: ", basename(path))
  old_timeout <- getOption("timeout")
  options(timeout = max(1800, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  utils::download.file(
    url, path, mode = "wb", quiet = FALSE, method = "auto"
  )
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Download failed or produced an empty file: ", path)
  }
  invisible(path)
}

download_if_needed(object_url, object_path)

read_double_gzip_rds <- function(path) {
  # GEO stores this particular RDS with two gzip layers.
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection, text = FALSE)
  on.exit(try(close(inner_connection), silent = TRUE), add = TRUE)
  readRDS(inner_connection)
}

message("Loading processed GSE130973 object. This may take several minutes.")
object <- read_double_gzip_rds(object_path)

object_summary <- capture.output({
  cat("Class:\n")
  print(class(object))
  cat("Attributes:\n")
  print(names(attributes(object)))
  cat("Structure (max.level=2):\n")
  str(object, max.level = 2)
})
writeLines(
  object_summary,
  file.path(result_dir, "GSE130973_object_structure_summary_v1.txt"),
  useBytes = TRUE
)

extract_metadata <- function(object) {
  if (isS4(object) && "meta.data" %in% methods::slotNames(object)) {
    return(as.data.frame(methods::slot(object, "meta.data")))
  }
  if (is.list(object) && !is.null(object$meta.data)) {
    return(as.data.frame(object$meta.data))
  }
  if (is.list(object) && !is.null(object$metadata)) {
    return(as.data.frame(object$metadata))
  }
  NULL
}

metadata <- extract_metadata(object)
if (is.null(metadata) || nrow(metadata) == 0L) {
  stop("No cell metadata table could be extracted from the GSE130973 object.")
}

field_inventory <- do.call(rbind, lapply(names(metadata), function(field) {
  values <- as.character(metadata[[field]])
  values[is.na(values)] <- ""
  unique_values <- unique(values[nzchar(trimws(values))])
  data.frame(
    field = field,
    class = class(metadata[[field]])[[1L]],
    cells = nrow(metadata),
    nonempty_values = sum(nzchar(trimws(values))),
    unique_nonempty_values = length(unique_values),
    sample_donor_candidate = grepl(
      "sample|donor|subject|patient|individual|orig.ident|source|batch",
      field, ignore.case = TRUE, perl = TRUE
    ),
    celltype_candidate = grepl(
      "cell.?type|annotation|cluster|ident|subtype|state|fibro|lineage",
      field, ignore.case = TRUE, perl = TRUE
    ),
    example_values = paste(head(unique_values, 8L), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
rownames(field_inventory) <- NULL
write.csv(
  field_inventory,
  file.path(result_dir, "GSE130973_metadata_field_inventory_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

sample_fields <- field_inventory[field_inventory$sample_donor_candidate, , drop = FALSE]
celltype_fields <- field_inventory[field_inventory$celltype_candidate, , drop = FALSE]
write.csv(
  sample_fields,
  file.path(result_dir, "GSE130973_sample_donor_candidate_fields_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  celltype_fields,
  file.path(result_dir, "GSE130973_celltype_candidate_fields_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

candidate_genes <- read.csv(candidate_panel_path, check.names = FALSE)
extract_gene_names <- function(object) {
  if (!isS4(object) || !"assays" %in% methods::slotNames(object)) {
    return(character())
  }
  assays <- methods::slot(object, "assays")
  if (length(assays) == 0L) return(character())
  assay_names <- names(assays)
  preferred <- which(tolower(assay_names) %in% c("rna", "gene expression"))
  assay <- assays[[if (length(preferred) > 0L) preferred[[1L]] else 1L]]
  possible_slots <- intersect(
    c("counts", "data", "scale.data"), methods::slotNames(assay)
  )
  for (slot_name in possible_slots) {
    values <- tryCatch(
      rownames(methods::slot(assay, slot_name)),
      error = function(e) character()
    )
    if (length(values) > 0L) return(values)
  }
  character()
}
gene_names <- extract_gene_names(object)
gene_inventory <- data.frame(
  gene = candidate_genes$gene,
  present_in_object = candidate_genes$gene %in% gene_names,
  stringsAsFactors = FALSE
)
write.csv(
  gene_inventory,
  file.path(result_dir, "GSE130973_frozen_candidate_gene_inventory_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

decision_lines <- c(
  "# Step 11D GSE130973 Seurat object audit",
  "",
  "## Scope",
  "",
  "- Source: GSE130973, non-overlapping candidate after Step 10A.",
  "- This step downloads one processed Seurat RDS. The GSE130973 supplementary directory does not provide a filelist.txt; the file inventory was already recorded in Step 11B.",
  "- No differential expression or mechanosensitivity test is performed.",
  "",
  "## Gate",
  "",
  paste0("- Cell metadata rows: ", nrow(metadata), "."),
  paste0("- Sample/donor candidate fields: ", nrow(sample_fields), "."),
  paste0("- Cell-type candidate fields: ", nrow(celltype_fields), "."),
  paste0("- Frozen candidate genes found: ", sum(gene_inventory$present_in_object), "/", nrow(gene_inventory), "."),
  "",
  "The object can proceed to expression validation only if a cell-state field and a sample/donor field are both confirmed. Any donor-level test must use the object-level donor field rather than the five GEO series subjects alone.",
  "",
  "## Material Passport",
  "",
  "- Input: GSE130973_seurat_analysis_lyko.rds.gz from the official GEO supplementary directory.",
  "- Transformation: object structure, metadata fields, and frozen-gene presence were audited without testing.",
  "- Reproducibility: the downloaded RDS is retained under the project raw-data directory and all audit tables are exported."
)
decision_path <- file.path(
  result_dir, "GSE130973_seurat_object_audit_decision_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11D GSE130973 Seurat object audit completed.")
message("Metadata rows: ", nrow(metadata))
message("Sample/donor candidate fields: ", nrow(sample_fields))
message("Cell-type candidate fields: ", nrow(celltype_fields))
message("Frozen genes present: ", sum(gene_inventory$present_in_object), "/", nrow(gene_inventory))
message("Decision: ", decision_path)
