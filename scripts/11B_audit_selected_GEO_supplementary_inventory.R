options(stringsAsFactors = FALSE)

# Step 11B: inventory selected independent GEO supplementary files.
# This is a metadata/file-list audit only. It does not download matrices,
# H5AD files, RDS objects, or raw sequencing archives.

project_dir <- "."
step11a_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11A_metadata_only_source_audit"
)
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11B_selected_GEO_supplementary_inventory"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

priority_path <- file.path(
  step11a_dir, "step11_independent_source_priority_audit_v1.csv"
)
snapshot_path <- file.path(
  step11a_dir, "step11_independent_source_sample_metadata_snapshot_v1.csv"
)
required_inputs <- c(priority_path, snapshot_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 11A input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

priority <- read.csv(priority_path, check.names = FALSE)
snapshot <- read.csv(snapshot_path, check.names = FALSE)

# The first-pass queue is deliberately small and biologically interpretable:
# GSE175817 is the top metadata-priority source; GSE130973 and GSE202352 have
# explicit subjects and skin/dermis samples; GSE228421 has repeated patients;
# GSE138669 is a disease-relevant fibrotic skin comparator.
selected_accessions <- c(
  "GSE175817", "GSE130973", "GSE202352", "GSE228421", "GSE138669"
)

extract_hrefs <- function(lines) {
  matches <- unlist(regmatches(
    lines,
    gregexpr('href="[^"]+"', lines, perl = TRUE)
  ))
  if (length(matches) == 0L) return(character())
  filenames <- sub('^href="', "", matches)
  sub('"$', "", filenames)
}

source_rows <- list()
file_rows <- list()

for (accession in selected_accessions) {
  numeric_part <- sub("^GSE", "", accession)
  series_prefix <- paste0("GSE", sub("[0-9]{3}$", "nnn", numeric_part))
  directory_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    series_prefix, "/", accession, "/suppl/"
  )
  message("Listing supplementary inventory only: ", accession)

  result <- tryCatch({
    lines <- readLines(directory_url, warn = FALSE, encoding = "UTF-8")
    hrefs <- extract_hrefs(lines)
    hrefs <- hrefs[!hrefs %in% c("../", "./")]
    hrefs <- hrefs[!grepl("/$", hrefs)]
    if (length(hrefs) == 0L) stop("No supplementary filenames found.")

    local_snapshot <- snapshot[snapshot$source_accession == accession, , drop = FALSE]
    fibroblast_metadata <- sum(local_snapshot$fibroblast_like %in% TRUE)
    donor_metadata <- sum(local_snapshot$donor_like %in% TRUE)
    single_cell_metadata <- sum(local_snapshot$single_cell_like %in% TRUE)

    local_files <- data.frame(
      source_accession = accession,
      filename = hrefs,
      direct_url = paste0(directory_url, hrefs),
      is_processed_object = grepl(
        "rds|h5ad|h5|loom|mtx|matrix|barcodes|features|csv|tsv|txt|zip|tar",
        hrefs, ignore.case = TRUE, perl = TRUE
      ),
      is_raw_archive = grepl(
        "raw|fastq|bam|sra|tar$|fastq.gz|fq.gz",
        hrefs, ignore.case = TRUE, perl = TRUE
      ),
      is_metadata_or_document = grepl(
        "soft|xml|miniml|readme|md|pdf|xlsx|xls|docx",
        hrefs, ignore.case = TRUE, perl = TRUE
      ),
      stringsAsFactors = FALSE
    )
    file_rows[[length(file_rows) + 1L]] <- local_files

    processed_count <- sum(local_files$is_processed_object & !local_files$is_raw_archive)
    data.frame(
      source_accession = accession,
      status = "supplementary_inventory_completed",
      supplementary_files = nrow(local_files),
      processed_or_matrix_candidates = processed_count,
      raw_archive_candidates = sum(local_files$is_raw_archive),
      fibroblast_metadata_samples = fibroblast_metadata,
      donor_metadata_samples = donor_metadata,
      single_cell_metadata_samples = single_cell_metadata,
      directory_url = directory_url,
      error_message = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      source_accession = accession,
      status = "supplementary_inventory_failed",
      supplementary_files = NA_integer_,
      processed_or_matrix_candidates = NA_integer_,
      raw_archive_candidates = NA_integer_,
      fibroblast_metadata_samples = NA_integer_,
      donor_metadata_samples = NA_integer_,
      single_cell_metadata_samples = NA_integer_,
      directory_url = directory_url,
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  source_rows[[length(source_rows) + 1L]] <- result
}

source_inventory <- do.call(rbind, source_rows)
file_inventory <- if (length(file_rows) > 0L) {
  do.call(rbind, file_rows)
} else {
  data.frame()
}

safe_write_csv(
  source_inventory,
  file.path(result_dir, "step11B_selected_source_inventory_v1.csv")
)
safe_write_csv(
  file_inventory,
  file.path(result_dir, "step11B_selected_supplementary_file_inventory_v1.csv")
)

queue <- file_inventory[
  file_inventory$is_processed_object & !file_inventory$is_raw_archive,
  c("source_accession", "filename", "direct_url", "is_metadata_or_document"),
  drop = FALSE
]
queue$download_status <- "not_downloaded_metadata_only"
queue$next_action <- ifelse(
  grepl("rds|h5ad|h5|loom|mtx|matrix|zip|tar", queue$filename,
    ignore.case = TRUE, perl = TRUE
  ),
  "manual_object_format_audit_before_download",
  "metadata_or_small_matrix_review"
)
safe_write_csv(
  queue,
  file.path(result_dir, "step11B_processed_object_candidate_queue_v1.csv")
)

decision_lines <- c(
  "# Step 11B selected-source supplementary inventory",
  "",
  "## Scope",
  "",
  "- Selected sources: GSE175817, GSE130973, GSE202352, GSE228421, and GSE138669.",
  "- Only FTP directory listings were read.",
  "- No supplementary file was downloaded.",
  "",
  "## Next selection rule",
  "",
  "Prioritize a processed single-cell object with explicit cell-state labels and a sample/donor field. A raw archive, a matrix without cell-state labels, or a repeated-measures source without donor mapping is not sufficient for independent donor-level replication.",
  "",
  "GSE228421 must be analysed as a repeated-patient/longitudinal design if selected; its GSM count must not be treated as the number of independent biological samples.",
  "",
  "## Material Passport",
  "",
  "- Origin: Step 11A metadata-only screening and official NCBI GEO FTP directory listings.",
  "- Transformation: file names and candidate URLs were inventoried; no expression values were read.",
  "- Reproducibility: all URLs and classification flags are exported in the result directory.",
  "- Integrity boundary: file-name classification is preliminary and requires inspection of the selected object before download."
)
decision_path <- file.path(result_dir, "step11B_supplementary_inventory_decision_v1.md")
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11B selected-source supplementary inventory completed.")
message("Source inventory: ", file.path(result_dir, "step11B_selected_source_inventory_v1.csv"))
message("File inventory: ", file.path(result_dir, "step11B_selected_supplementary_file_inventory_v1.csv"))
message("Processed-object queue: ", file.path(result_dir, "step11B_processed_object_candidate_queue_v1.csv"))
message("Decision: ", decision_path)
