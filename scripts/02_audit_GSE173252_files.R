# Sample-level and file-structure audit for GSE173252.
# This script does not perform cell filtering, normalization, clustering, or
# differential expression. It only extracts the deposited 10X files, maps
# them to GEO samples, and verifies matrix dimensions against feature/barcode
# files using bounded-memory line counting.

options(stringsAsFactors = FALSE)

project_dir <- "."
archive_path <- file.path(
  project_dir,
  "data", "raw", "GSE173252", "GSE173252_RAW.tar"
)
extract_dir <- file.path(project_dir, "data", "extracted", "GSE173252")
metadata_dir <- file.path(project_dir, "data", "metadata", "GSE173252")
result_dir <- file.path(project_dir, "results", "01_data_audit", "GSE173252")

dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(archive_path)) {
  stop("Missing archive: ", archive_path)
}

geo_metadata_path <- file.path(
  metadata_dir,
  "GSE173252_series_matrix_sample_metadata_01.csv"
)
if (!file.exists(geo_metadata_path)) {
  stop("Missing GEO sample metadata: ", geo_metadata_path)
}

archive_listing <- utils::untar(archive_path, list = TRUE)
utils::write.csv(
  data.frame(archive_member = archive_listing),
  file.path(result_dir, "GSE173252_archive_inventory.csv"),
  row.names = FALSE
)

expected_archive_members <- archive_listing[
  grepl("_(barcodes[.]tsv[.]gz|features[.]tsv[.]gz|matrix[.]mtx[.]gz)$",
        archive_listing)
]

if (length(expected_archive_members) != 27L) {
  stop(
    "Expected 27 files (9 samples x 3 files), but found ",
    length(expected_archive_members), "."
  )
}

missing_extracted <- expected_archive_members[
  !file.exists(file.path(extract_dir, expected_archive_members))
]

if (length(missing_extracted) > 0L) {
  message("Extracting GSE173252_RAW.tar to: ", extract_dir)
  utils::untar(
    tarfile = archive_path,
    files = missing_extracted,
    exdir = extract_dir
  )
} else {
  message("Using previously extracted 10X files: ", extract_dir)
}

count_gzip_lines <- function(path, chunk_size = 100000L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  total <- 0L
  repeat {
    block <- readLines(con, n = chunk_size, warn = FALSE)
    total <- total + length(block)
    if (length(block) < chunk_size) {
      break
    }
  }
  as.numeric(total)
}

read_matrix_market_dimensions <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      stop("No dimension line found in: ", path)
    }
    if (!startsWith(line, "%")) {
      break
    }
  }

  values <- as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1L]])
  if (length(values) != 3L || anyNA(values)) {
    stop("Invalid Matrix Market dimension line in: ", path)
  }

  names(values) <- c("n_features", "n_barcodes", "n_nonzero_entries")
  values
}

geo_metadata <- utils::read.csv(
  geo_metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_geo_columns <- c(
  "geo_accession", "title", "source_name_ch1", "characteristics_ch1",
  "platform_id", "library_strategy", "instrument_model"
)
missing_geo_columns <- setdiff(required_geo_columns, names(geo_metadata))
if (length(missing_geo_columns) > 0L) {
  stop("Missing GEO metadata columns: ", paste(missing_geo_columns, collapse = ", "))
}

matrix_files <- sort(list.files(
  extract_dir,
  pattern = "_matrix[.]mtx[.]gz$",
  full.names = TRUE
))

if (length(matrix_files) != 9L) {
  stop("Expected 9 matrix files, but found ", length(matrix_files), ".")
}

audit_rows <- vector("list", length(matrix_files))

for (i in seq_along(matrix_files)) {
  matrix_path <- matrix_files[[i]]
  prefix <- sub("_matrix[.]mtx[.]gz$", "", basename(matrix_path))
  gsm_id <- sub("_.*$", "", prefix)

  barcode_path <- file.path(extract_dir, paste0(prefix, "_barcodes.tsv.gz"))
  feature_path <- file.path(extract_dir, paste0(prefix, "_features.tsv.gz"))

  if (!file.exists(barcode_path) || !file.exists(feature_path)) {
    stop("Incomplete 10X file triplet for: ", prefix)
  }

  geo_index <- match(gsm_id, geo_metadata$geo_accession)
  if (is.na(geo_index)) {
    stop("No GEO metadata row found for: ", gsm_id)
  }

  dims <- read_matrix_market_dimensions(matrix_path)
  barcode_lines <- count_gzip_lines(barcode_path)
  feature_lines <- count_gzip_lines(feature_path)

  title <- geo_metadata$title[[geo_index]]
  tissue_group <- if (grepl("Dermis", title, ignore.case = TRUE)) {
    "healthy_dermis"
  } else if (grepl("Skoog", title, ignore.case = TRUE)) {
    "nonpathogenic_skoog_fascia"
  } else if (grepl("Dupuytren", title, ignore.case = TRUE)) {
    "dupuytren_disease"
  } else {
    "unresolved"
  }

  disease_status <- if (identical(tissue_group, "dupuytren_disease")) {
    "disease"
  } else {
    "non_disease"
  }

  audit_rows[[i]] <- data.frame(
    gse_id = "GSE173252",
    gsm_id = gsm_id,
    sample_title = title,
    tissue_group = tissue_group,
    disease_status = disease_status,
    source_name = geo_metadata$source_name_ch1[[geo_index]],
    characteristics = geo_metadata$characteristics_ch1[[geo_index]],
    platform_id = geo_metadata$platform_id[[geo_index]],
    library_strategy = geo_metadata$library_strategy[[geo_index]],
    instrument_model = geo_metadata$instrument_model[[geo_index]],
    n_features_matrix = unname(dims[["n_features"]]),
    n_barcodes_matrix = unname(dims[["n_barcodes"]]),
    n_nonzero_entries = unname(dims[["n_nonzero_entries"]]),
    n_feature_lines = feature_lines,
    n_barcode_lines = barcode_lines,
    feature_dimension_match = unname(dims[["n_features"]]) == feature_lines,
    barcode_dimension_match = unname(dims[["n_barcodes"]]) == barcode_lines,
    deposited_matrix_layer = if (unname(dims[["n_barcodes"]]) >= 100000) {
      "unfiltered_droplet_matrix"
    } else {
      "filtered_or_small_matrix"
    },
    matrix_file = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
    features_file = normalizePath(feature_path, winslash = "/", mustWork = TRUE),
    barcodes_file = normalizePath(barcode_path, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, audit_rows)
audit <- audit[order(audit$gsm_id), , drop = FALSE]

manifest_columns <- c(
  "gse_id", "gsm_id", "sample_title", "tissue_group", "disease_status",
  "source_name", "characteristics", "platform_id", "library_strategy",
  "instrument_model", "deposited_matrix_layer", "matrix_file",
  "features_file", "barcodes_file"
)

utils::write.csv(
  audit[, manifest_columns, drop = FALSE],
  file.path(metadata_dir, "GSE173252_sample_manifest_audited.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  audit,
  file.path(result_dir, "GSE173252_matrix_structure_audit.csv"),
  row.names = FALSE,
  na = ""
)

all_triplets_valid <- all(
  audit$feature_dimension_match & audit$barcode_dimension_match
)
all_unfiltered <- all(audit$deposited_matrix_layer == "unfiltered_droplet_matrix")

decision_lines <- c(
  "GSE173252 data-landing audit decision",
  paste0("Audit date: ", format(Sys.Date(), "%Y-%m-%d")),
  paste0("Samples audited: ", nrow(audit)),
  paste0("All 10X triplets structurally consistent: ", all_triplets_valid),
  paste0("All deposited matrices classified as unfiltered droplet matrices: ", all_unfiltered),
  "",
  "Interpretation:",
  "- Matrix columns are barcode partitions, not automatically accepted cells.",
  "- Do not use all deposited columns as biological cells or replicates.",
  "- Prefer the submitter-provided processed RDS to reproduce published filtering and annotations.",
  "- Preserve GSM/sample identity and use sample-level units for downstream inference.",
  "- Do not start CellChat, spatial analysis, or GWAS at this stage."
)

writeLines(
  decision_lines,
  file.path(result_dir, "GSE173252_data_landing_decision.txt")
)

print(audit[, c(
  "gsm_id", "sample_title", "tissue_group", "n_features_matrix",
  "n_barcodes_matrix", "n_nonzero_entries", "feature_dimension_match",
  "barcode_dimension_match", "deposited_matrix_layer"
)])

if (!all_triplets_valid) {
  stop("Audit failed: at least one 10X triplet has inconsistent dimensions.")
}

message("\nGSE173252 file audit completed successfully.")
message("Sample manifest: ", file.path(metadata_dir, "GSE173252_sample_manifest_audited.csv"))
message("Audit table: ", file.path(result_dir, "GSE173252_matrix_structure_audit.csv"))
message("Decision: ", file.path(result_dir, "GSE173252_data_landing_decision.txt"))
