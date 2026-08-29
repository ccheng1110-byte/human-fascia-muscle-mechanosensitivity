# Inspect the submitter-provided legacy Seurat RDS for GSE173252.
# The object was created with Seurat 3.1.2 and lacks slots required by Seurat
# 5. To preserve the original object and avoid class-definition conflicts,
# this script must be run in a clean R session WITHOUT loading Seurat. It reads
# the legacy S4 slots through base-R attributes and performs metadata/structure
# inspection only.

options(stringsAsFactors = FALSE)

project_dir <- "."
input_path <- file.path(
  project_dir,
  "data", "raw", "GSE173252", "GSE173252_dd_mesenchyme.rds.gz"
)
metadata_dir <- file.path(project_dir, "data", "metadata", "GSE173252")
result_dir <- file.path(project_dir, "results", "01_data_audit", "GSE173252")

dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Missing processed RDS: ", input_path)
}

if ("package:Seurat" %in% search() || "Seurat" %in% loadedNamespaces()) {
  stop(
    "Seurat is currently loaded. Restart R (Ctrl+Shift+F10) and run this ",
    "script directly without calling library(Seurat)."
  )
}

message("Loading processed object: ", input_path)
message("This may take several minutes because the compressed file is about 693 MB.")

# GEO stores this particular *.rds.gz file with two gzip layers. Reading it
# through a single gzfile() leaves an inner gzip stream rather than an RDS
# serialization header. Keep the connections local so failed reads do not
# leave invalid connection objects in the Global Environment.
read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection, text = FALSE)
  on.exit(try(close(inner_connection), silent = TRUE), add = TRUE)
  readRDS(inner_connection)
}

object <- read_double_gzip_rds(input_path)

object_attributes <- attributes(object)
object_class <- as.character(object_attributes[["class"]])

required_legacy_slots <- c(
  "assays", "meta.data", "active.assay", "active.ident", "reductions",
  "version"
)
missing_legacy_slots <- setdiff(required_legacy_slots, names(object_attributes))
if (length(missing_legacy_slots) > 0L) {
  stop(
    "Legacy object is missing required slots: ",
    paste(missing_legacy_slots, collapse = ", ")
  )
}

cell_metadata <- object_attributes[["meta.data"]]
if (!is.data.frame(cell_metadata)) {
  stop("Legacy meta.data slot is not a data.frame.")
}

cell_metadata$cell_barcode <- rownames(cell_metadata)
cell_metadata <- cell_metadata[, c(
  "cell_barcode",
  setdiff(names(cell_metadata), "cell_barcode")
), drop = FALSE]

metadata_path <- file.path(
  metadata_dir,
  "GSE173252_dd_mesenchyme_cell_metadata.csv.gz"
)
metadata_connection <- gzfile(metadata_path, open = "wt")
utils::write.csv(cell_metadata, metadata_connection, row.names = FALSE, na = "")
close(metadata_connection)

compact_example <- function(x, n = 5L) {
  values <- unique(as.character(x[!is.na(x)]))
  values <- values[nzchar(values)]
  paste(utils::head(values, n), collapse = " | ")
}

column_inventory <- data.frame(
  column = names(cell_metadata),
  class = vapply(
    cell_metadata,
    function(x) paste(class(x), collapse = ","),
    character(1)
  ),
  n_unique_nonmissing = vapply(
    cell_metadata,
    function(x) length(unique(x[!is.na(x)])),
    integer(1)
  ),
  n_missing = vapply(cell_metadata, function(x) sum(is.na(x)), integer(1)),
  examples = vapply(cell_metadata, compact_example, character(1)),
  stringsAsFactors = FALSE
)

utils::write.csv(
  column_inventory,
  file.path(result_dir, "GSE173252_mesenchyme_metadata_column_inventory.csv"),
  row.names = FALSE,
  na = ""
)

candidate_pattern <- paste(
  c(
    "orig[.]ident", "sample", "patient", "donor", "source", "condition",
    "disease", "tissue", "group", "cell.?type", "annotation", "cluster"
  ),
  collapse = "|"
)

candidate_columns <- names(cell_metadata)[grepl(
  candidate_pattern,
  names(cell_metadata),
  ignore.case = TRUE,
  perl = TRUE
)]

candidate_counts <- list()
count_index <- 1L

for (column_name in candidate_columns) {
  values <- as.character(cell_metadata[[column_name]])
  counts <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
  candidate_counts[[count_index]] <- data.frame(
    metadata_column = column_name,
    value = names(counts),
    n_cells = as.numeric(counts),
    stringsAsFactors = FALSE
  )
  count_index <- count_index + 1L
}

if (length(candidate_counts) > 0L) {
  candidate_counts <- do.call(rbind, candidate_counts)
} else {
  candidate_counts <- data.frame(
    metadata_column = character(),
    value = character(),
    n_cells = numeric(),
    stringsAsFactors = FALSE
  )
}

utils::write.csv(
  candidate_counts,
  file.path(result_dir, "GSE173252_mesenchyme_candidate_metadata_counts.csv"),
  row.names = FALSE,
  na = ""
)

assay_names <- names(object_attributes[["assays"]])
reduction_names <- names(object_attributes[["reductions"]])
legacy_version <- as.character(object_attributes[["version"]])

active_assay <- as.character(object_attributes[["active.assay"]])
active_assay_object <- object_attributes[["assays"]][[active_assay]]
active_assay_attributes <- attributes(active_assay_object)

matrix_candidates <- c("counts", "data", "scale.data")
matrix_candidates <- matrix_candidates[
  matrix_candidates %in% names(active_assay_attributes)
]

matrix_dimensions <- NULL
for (matrix_slot in matrix_candidates) {
  candidate <- active_assay_attributes[[matrix_slot]]
  candidate_attributes <- attributes(candidate)
  candidate_dim <- candidate_attributes[["Dim"]]
  if (is.null(candidate_dim)) {
    candidate_dim <- dim(candidate)
  }
  if (!is.null(candidate_dim) && length(candidate_dim) == 2L && all(candidate_dim > 0L)) {
    matrix_dimensions <- as.numeric(candidate_dim)
    break
  }
}

if (is.null(matrix_dimensions)) {
  n_features <- NA_real_
  n_cells_in_assay <- NA_real_
} else {
  n_features <- matrix_dimensions[[1L]]
  n_cells_in_assay <- matrix_dimensions[[2L]]
}

summary_lines <- c(
  "GSE173252 processed mesenchyme RDS inspection",
  paste0("Inspection date: ", format(Sys.Date(), "%Y-%m-%d")),
  paste0("Object class: ", paste(object_class, collapse = ", ")),
  paste0("Legacy Seurat object version: ", legacy_version),
  paste0("Object size in memory (GB): ", round(as.numeric(object.size(object)) / 1024^3, 3)),
  paste0("Active assay: ", active_assay),
  paste0("Features in inspected assay matrix: ", n_features),
  paste0("Cells in inspected assay matrix: ", n_cells_in_assay),
  paste0("Cells in metadata: ", nrow(cell_metadata)),
  paste0("Assays: ", paste(assay_names, collapse = ", ")),
  paste0("Reductions: ", paste(reduction_names, collapse = ", ")),
  paste0("Metadata columns: ", ncol(cell_metadata) - 1L),
  paste0("Unique cell barcodes: ", !anyDuplicated(cell_metadata$cell_barcode)),
  paste0(
    "Metadata row count equals assay cell count: ",
    is.na(n_cells_in_assay) || nrow(cell_metadata) == n_cells_in_assay
  ),
  paste0("Candidate sample/group/annotation columns: ", paste(candidate_columns, collapse = ", ")),
  "",
  "Decision rule:",
  "- If a reliable sample/GSM column and cell annotation column are present, preserve them for reproduction.",
  "- If sample identity is missing or collapsed, do not perform donor-level differential testing from this object.",
  "- Do not overwrite or update the downloaded legacy object during this inspection.",
  "- Any later conversion to Seurat 5 must be saved as a separate derived object."
)

writeLines(
  summary_lines,
  file.path(result_dir, "GSE173252_mesenchyme_object_summary.txt")
)

message("\nProcessed RDS inspection completed successfully.")
message("Object summary: ", file.path(result_dir, "GSE173252_mesenchyme_object_summary.txt"))
message("Metadata inventory: ", file.path(result_dir, "GSE173252_mesenchyme_metadata_column_inventory.csv"))
message("Candidate metadata counts: ", file.path(result_dir, "GSE173252_mesenchyme_candidate_metadata_counts.csv"))
message("Cell metadata: ", metadata_path)

rm(object)
invisible(gc())
