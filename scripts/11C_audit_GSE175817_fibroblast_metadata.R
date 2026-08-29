options(stringsAsFactors = FALSE)

# Step 11C: audit the small fibroblast-specific metadata file from GSE175817.
# This does not download expression matrices or raw sequencing archives.

project_dir <- "."
metadata_dir <- file.path(
  project_dir, "data", "metadata", "independent_sources", "GSE175817"
)
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11C_GSE175817_fibroblast_metadata"
)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175817/suppl"
metadata_url <- paste0(base_url, "/GSE175817_meta_fibroblast.csv.gz")
filelist_url <- paste0(base_url, "/filelist.txt")
metadata_path <- file.path(metadata_dir, "GSE175817_meta_fibroblast.csv.gz")
filelist_path <- file.path(metadata_dir, "GSE175817_filelist.txt")

download_if_needed <- function(url, path) {
  if (file.exists(path) && file.info(path)$size > 0) {
    message("Using existing file: ", path)
    return(invisible(path))
  }
  message("Downloading small metadata file: ", basename(path))
  utils::download.file(
    url, path, mode = "wb", quiet = FALSE, method = "auto"
  )
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Download failed or produced an empty file: ", path)
  }
  invisible(path)
}

download_if_needed(metadata_url, metadata_path)
download_if_needed(filelist_url, filelist_path)

metadata <- read.csv(
  gzfile(metadata_path, open = "rt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(metadata) == 0L || ncol(metadata) == 0L) {
  stop("The fibroblast metadata file is empty.")
}

field_inventory <- do.call(rbind, lapply(names(metadata), function(field) {
  values <- as.character(metadata[[field]])
  values[is.na(values)] <- ""
  nonempty <- nzchar(trimws(values))
  unique_values <- unique(values[nonempty])
  data.frame(
    field = field,
    class = class(metadata[[field]])[[1L]],
    nonempty_values = sum(nonempty),
    unique_nonempty_values = length(unique_values),
    keyword_sample_donor = grepl(
      "sample|donor|subject|patient|individual|replicate",
      field, ignore.case = TRUE, perl = TRUE
    ),
    keyword_celltype_cluster = grepl(
      "cell.?type|cluster|annotation|ident|subtype|state|fibro",
      field, ignore.case = TRUE, perl = TRUE
    ),
    example_values = paste(head(unique_values, 5L), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
rownames(field_inventory) <- NULL
write.csv(
  field_inventory,
  file.path(result_dir, "GSE175817_fibroblast_metadata_field_inventory_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

sample_candidates <- field_inventory[
  field_inventory$keyword_sample_donor,
  , drop = FALSE
]
celltype_candidates <- field_inventory[
  field_inventory$keyword_celltype_cluster,
  , drop = FALSE
]

choose_field <- function(candidates, preferred_pattern) {
  if (nrow(candidates) == 0L) return(NA_character_)
  preferred <- candidates[
    grepl(preferred_pattern, candidates$field, ignore.case = TRUE, perl = TRUE),
    , drop = FALSE
  ]
  if (nrow(preferred) > 0L) return(preferred$field[[1L]])
  candidates$field[[1L]]
}

sample_field <- choose_field(
  sample_candidates,
  "donor|subject|patient|sample|individual"
)
celltype_field <- choose_field(
  celltype_candidates,
  "cell.?type|annotation|subtype|state|cluster|ident"
)

sample_coverage <- data.frame()
if (!is.na(sample_field)) {
  sample_values <- as.character(metadata[[sample_field]])
  sample_values[is.na(sample_values)] <- ""
  sample_values <- trimws(sample_values)
  sample_values[sample_values == ""] <- NA_character_
  sample_coverage <- as.data.frame(table(sample_values, useNA = "ifany"))
  names(sample_coverage) <- c("sample_or_donor_id", "cells")
  sample_coverage$sample_or_donor_id <- as.character(
    sample_coverage$sample_or_donor_id
  )
  write.csv(
    sample_coverage,
    file.path(result_dir, "GSE175817_fibroblast_sample_coverage_v1.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

celltype_coverage <- data.frame()
if (!is.na(celltype_field)) {
  celltype_values <- as.character(metadata[[celltype_field]])
  celltype_values[is.na(celltype_values)] <- ""
  celltype_coverage <- as.data.frame(table(celltype_values))
  names(celltype_coverage) <- c("celltype_or_cluster", "cells")
  write.csv(
    celltype_coverage,
    file.path(result_dir, "GSE175817_fibroblast_celltype_coverage_v1.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

fibroblast_label_count <- if (nrow(celltype_coverage) > 0L) {
  sum(grepl(
    "fibro|myofib|mesenchymal|stromal|fascia",
    celltype_coverage$celltype_or_cluster,
    ignore.case = TRUE,
    perl = TRUE
  ) * celltype_coverage$cells)
} else {
  NA_integer_
}

decision_lines <- c(
  "# Step 11C GSE175817 fibroblast metadata audit",
  "",
  "## Input",
  "",
  paste0("- Metadata URL: ", metadata_url),
  paste0("- Local metadata: ", metadata_path),
  paste0("- Cells/rows: ", nrow(metadata), "; fields: ", ncol(metadata), "."),
  "- No expression matrix or raw sequencing archive was downloaded.",
  "",
  "## Resolved fields",
  "",
  paste0("- Candidate sample/donor field: ", sample_field, "."),
  paste0("- Candidate cell-type/cluster field: ", celltype_field, "."),
  paste0("- Fibroblast-like labelled cells, if a cell-type field was found: ", fibroblast_label_count, "."),
  "",
  "## Gate",
  "",
  "This is a source-structure audit only. Independent validation is allowed only if the object-level expression source, cell-state field, and sample/donor mapping are all confirmed in the next step.",
  "",
  "## Material Passport",
  "",
  "- Source: GSE175817, non-overlapping candidate from corrected Step 10A.",
  "- Transformation: downloaded fibroblast-specific metadata were summarized by field, sample/donor, and cell-type coverage.",
  "- No gene-level hypothesis test was performed.",
  "- Reproducibility: raw compressed metadata and the small filelist are retained under the project metadata directory."
)
decision_path <- file.path(
  result_dir, "GSE175817_fibroblast_metadata_audit_decision_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11C GSE175817 fibroblast metadata audit completed.")
message("Rows: ", nrow(metadata), "; fields: ", ncol(metadata))
message("Candidate sample/donor field: ", sample_field)
message("Candidate cell-type field: ", celltype_field)
message("Decision: ", decision_path)
