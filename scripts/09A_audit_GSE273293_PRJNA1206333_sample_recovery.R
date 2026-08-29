options(stringsAsFactors = FALSE)
options(timeout = 300)

# Step 09A: audit whether the 14 reported GSE273293/PRJNA1206333 samples can
# be recovered as independent metadata units before any expression analysis.
#
# This script downloads ONLY the small NCBI RunInfo metadata table. It does not
# download SRA/FASTQ files and explicitly stops if the metadata does not match
# the expected 10 GMC + 4 control design.

project_dir <- "."
bioproject_id <- "PRJNA1206333"
geo_id <- "GSE273293"
expected_gmc <- sprintf("GMC%02d", 1:10)
expected_control <- sprintf("Con%02d", 1:4)
expected_samples <- c(expected_gmc, expected_control)

metadata_dir <- file.path(
  project_dir, "data", "metadata", geo_id, bioproject_id
)
result_dir <- file.path(
  project_dir, "results", "07_deep_fascia_recovery", geo_id,
  "09A_PRJNA1206333_metadata_audit"
)
log_dir <- file.path(project_dir, "logs")

for (path in c(metadata_dir, result_dir, log_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

runinfo_url <- paste0(
  "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=",
  bioproject_id
)
runinfo_path <- file.path(
  metadata_dir, paste0(bioproject_id, "_SRA_RunInfo.csv")
)

download_small_metadata <- function(url, destination, retries = 4L) {
  temporary_path <- paste0(destination, ".partial")
  if (file.exists(temporary_path)) {
    unlink(temporary_path)
  }

  for (attempt in seq_len(retries)) {
    message(sprintf(
      "Downloading small metadata file [%d/%d]: %s",
      attempt, retries, url
    ))

    success <- tryCatch({
      status <- utils::download.file(
        url = url,
        destfile = temporary_path,
        method = "libcurl",
        mode = "wb",
        quiet = FALSE
      )
      identical(status, 0L) &&
        file.exists(temporary_path) &&
        isTRUE(file.info(temporary_path)$size > 100L) &&
        isTRUE(file.info(temporary_path)$size < 5 * 1024^2)
    }, error = function(e) {
      message("Metadata download error: ", conditionMessage(e))
      FALSE
    })

    if (success) {
      copied <- file.copy(
        from = temporary_path,
        to = destination,
        overwrite = TRUE
      )
      unlink(temporary_path)
      if (!copied) {
        stop("Could not move verified metadata into place: ", destination)
      }
      return(invisible(destination))
    }

    if (file.exists(temporary_path)) {
      unlink(temporary_path)
    }
    if (attempt < retries) {
      Sys.sleep(2^attempt)
    }
  }

  stop("Small metadata download failed after retries: ", url)
}

safe_write_csv <- function(x, path) {
  utils::write.csv(
    x,
    file = path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

normalize_character <- function(x) {
  trimws(as.character(x))
}

required_columns <- c(
  "Run", "Experiment", "LibraryName", "LibraryStrategy",
  "LibrarySource", "LibraryLayout", "Platform", "Model",
  "BioProject", "Sample", "BioSample", "TaxID", "ScientificName",
  "SampleName", "Sex", "size_MB", "spots", "bases", "download_path"
)

message("Step 09A audits metadata only; no SRA/FASTQ file will be downloaded.")
download_small_metadata(runinfo_url, runinfo_path)

runinfo <- utils::read.csv(
  runinfo_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
runinfo <- runinfo[
  !is.na(runinfo$Run) & nzchar(trimws(runinfo$Run)),
  , drop = FALSE
]

missing_columns <- setdiff(required_columns, names(runinfo))
if (length(missing_columns) > 0L) {
  stop(
    "RunInfo is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

for (column in intersect(
  c(
    "Run", "Experiment", "LibraryName", "LibraryStrategy",
    "LibrarySource", "LibraryLayout", "Platform", "Model",
    "BioProject", "Sample", "BioSample", "ScientificName",
    "SampleName", "Sex", "download_path"
  ),
  names(runinfo)
)) {
  runinfo[[column]] <- normalize_character(runinfo[[column]])
}

runinfo$size_MB <- suppressWarnings(as.numeric(runinfo$size_MB))
runinfo$spots <- suppressWarnings(as.numeric(runinfo$spots))
runinfo$bases <- suppressWarnings(as.numeric(runinfo$bases))
runinfo$TaxID <- suppressWarnings(as.integer(runinfo$TaxID))

if (any(runinfo$BioProject != bioproject_id)) {
  stop("RunInfo contains records outside ", bioproject_id, ".")
}

runinfo$clinical_group <- ifelse(
  grepl("^GMC[0-9]{2}$", runinfo$LibraryName),
  "GMC_fibrotic_fascia",
  ifelse(
    grepl("^Con[0-9]{2}$", runinfo$LibraryName),
    "nonfibrotic_control_fascia",
    "unresolved"
  )
)
runinfo$sample_number <- suppressWarnings(as.integer(
  sub("^[A-Za-z]+", "", runinfo$LibraryName)
))

actual_samples <- sort(unique(runinfo$LibraryName))
missing_expected_samples <- setdiff(expected_samples, actual_samples)
unexpected_samples <- setdiff(actual_samples, expected_samples)

count_unique <- function(x) length(unique(x[!is.na(x) & nzchar(x)]))

audit_checks <- data.frame(
  check_id = c(
    "exact_row_count_14",
    "unique_run_count_14",
    "unique_experiment_count_14",
    "unique_sra_sample_count_14",
    "unique_biosample_count_14",
    "unique_library_name_count_14",
    "exact_10_GMC_and_4_control",
    "all_expected_sample_names_present",
    "no_unexpected_sample_names",
    "all_human_taxon_9606",
    "all_paired_RNA_seq_transcriptomic",
    "all_raw_files_over_50_MB",
    "all_accessions_one_to_one"
  ),
  passed = c(
    nrow(runinfo) == 14L,
    count_unique(runinfo$Run) == 14L,
    count_unique(runinfo$Experiment) == 14L,
    count_unique(runinfo$Sample) == 14L,
    count_unique(runinfo$BioSample) == 14L,
    count_unique(runinfo$LibraryName) == 14L,
    sum(runinfo$clinical_group == "GMC_fibrotic_fascia") == 10L &&
      sum(runinfo$clinical_group == "nonfibrotic_control_fascia") == 4L,
    length(missing_expected_samples) == 0L,
    length(unexpected_samples) == 0L,
    all(runinfo$TaxID == 9606L) &&
      all(runinfo$ScientificName == "Homo sapiens"),
    all(runinfo$LibraryLayout == "PAIRED") &&
      all(runinfo$LibraryStrategy == "RNA-Seq") &&
      all(runinfo$LibrarySource == "TRANSCRIPTOMIC"),
    all(is.finite(runinfo$size_MB)) && all(runinfo$size_MB > 50),
    count_unique(runinfo$Run) == nrow(runinfo) &&
      count_unique(runinfo$Experiment) == nrow(runinfo) &&
      count_unique(runinfo$Sample) == nrow(runinfo) &&
      count_unique(runinfo$BioSample) == nrow(runinfo) &&
      count_unique(runinfo$LibraryName) == nrow(runinfo)
  ),
  stringsAsFactors = FALSE
)

metadata_gate_passed <- all(audit_checks$passed)
if (!metadata_gate_passed) {
  failed_checks <- audit_checks$check_id[!audit_checks$passed]
  stop(
    "PRJNA1206333 metadata gate failed: ",
    paste(failed_checks, collapse = "; "),
    ". No downstream expression analysis is allowed."
  )
}

sample_map_columns <- c(
  "clinical_group", "sample_number", "LibraryName", "SampleName", "Sex",
  "Run", "Experiment", "Sample", "BioSample", "BioProject",
  "LibraryStrategy", "LibrarySource", "LibraryLayout", "Platform", "Model",
  "spots", "bases", "size_MB", "download_path"
)
sample_map <- runinfo[, sample_map_columns, drop = FALSE]
sample_map$estimated_size_GiB <- sample_map$size_MB / 1024
sample_map$raw_download_authorized <- FALSE
sample_map$raw_download_reason <- paste(
  "Metadata audit only; full PRJNA1206333 raw reprocessing is not approved"
)
sample_map <- sample_map[order(
  factor(
    sample_map$clinical_group,
    levels = c("nonfibrotic_control_fascia", "GMC_fibrotic_fascia")
  ),
  sample_map$sample_number
), , drop = FALSE]

sample_map_path <- file.path(
  result_dir, "PRJNA1206333_sample_accession_map_v1.csv"
)
safe_write_csv(sample_map, sample_map_path)

checks_path <- file.path(
  result_dir, "PRJNA1206333_metadata_integrity_checks_v1.csv"
)
safe_write_csv(audit_checks, checks_path)

# Audit what is currently available locally for GSE273293. This is deliberately
# restricted to file structure and does not read or analyse the expression
# matrix itself.
local_paths <- c(
  file.path(
    project_dir, "data", "raw", geo_id, paste0(geo_id, "_RAW.tar")
  ),
  file.path(
    project_dir, "data", "metadata", geo_id,
    paste0(geo_id, "_series_matrix_sample_metadata_01.csv")
  ),
  file.path(
    project_dir, "data", "processed", geo_id,
    paste0(geo_id, "_series_matrix_01.rds")
  )
)

local_file_audit <- data.frame(
  file_role = c(
    "GEO supplementary tar",
    "GEO series-matrix sample metadata",
    "GEO parsed series-matrix object"
  ),
  path = normalizePath(local_paths, winslash = "/", mustWork = FALSE),
  exists = file.exists(local_paths),
  size_bytes = ifelse(
    file.exists(local_paths),
    file.info(local_paths)$size,
    NA_real_
  ),
  md5 = NA_character_,
  stringsAsFactors = FALSE
)
existing_local <- which(local_file_audit$exists)
if (length(existing_local) > 0L) {
  local_file_audit$md5[existing_local] <- unname(
    tools::md5sum(local_paths[existing_local])
  )
}

tar_path <- local_paths[[1L]]
tar_members <- if (file.exists(tar_path)) {
  utils::untar(tar_path, list = TRUE)
} else {
  character(0)
}
tar_sample_prefixes <- unique(sub(
  "_(barcodes|genes|features|matrix).*$", "", basename(tar_members)
))
tar_sample_prefixes <- tar_sample_prefixes[nzchar(tar_sample_prefixes)]

series_metadata_path <- local_paths[[2L]]
series_metadata <- if (file.exists(series_metadata_path)) {
  utils::read.csv(
    series_metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  data.frame()
}

local_structure <- data.frame(
  metric = c(
    "tar_member_count",
    "tar_sample_prefix_count",
    "tar_sample_prefixes",
    "series_metadata_row_count",
    "series_metadata_geo_accessions",
    "local_data_supports_14_sample_inference"
  ),
  value = c(
    as.character(length(tar_members)),
    as.character(length(tar_sample_prefixes)),
    paste(tar_sample_prefixes, collapse = ";"),
    as.character(nrow(series_metadata)),
    if ("geo_accession" %in% names(series_metadata)) {
      paste(unique(series_metadata$geo_accession), collapse = ";")
    } else {
      ""
    },
    "FALSE"
  ),
  stringsAsFactors = FALSE
)

local_file_audit_path <- file.path(
  result_dir, "GSE273293_local_file_audit_v1.csv"
)
local_structure_path <- file.path(
  result_dir, "GSE273293_local_sample_structure_v1.csv"
)
safe_write_csv(local_file_audit, local_file_audit_path)
safe_write_csv(local_structure, local_structure_path)

total_size_gib <- sum(runinfo$size_MB, na.rm = TRUE) / 1024
group_counts <- table(runinfo$clinical_group)
sex_counts <- table(runinfo$Sex, useNA = "ifany")

summary_table <- data.frame(
  metric = c(
    "metadata_gate_passed",
    "run_count",
    "experiment_count",
    "sra_sample_count",
    "biosample_count",
    "GMC_sample_count",
    "control_sample_count",
    "male_sample_count",
    "female_sample_count",
    "estimated_total_raw_size_GiB",
    "smallest_run_size_GiB",
    "largest_run_size_GiB",
    "local_GEO_sample_count",
    "local_GEO_supports_donor_level_inference",
    "proceed_to_full_raw_download",
    "proceed_to_step09B_processed_data_recovery"
  ),
  value = c(
    as.character(metadata_gate_passed),
    as.character(count_unique(runinfo$Run)),
    as.character(count_unique(runinfo$Experiment)),
    as.character(count_unique(runinfo$Sample)),
    as.character(count_unique(runinfo$BioSample)),
    as.character(unname(group_counts["GMC_fibrotic_fascia"])),
    as.character(unname(group_counts["nonfibrotic_control_fascia"])),
    as.character(unname(sex_counts["male"])),
    as.character(unname(sex_counts["female"])),
    sprintf("%.2f", total_size_gib),
    sprintf("%.2f", min(runinfo$size_MB, na.rm = TRUE) / 1024),
    sprintf("%.2f", max(runinfo$size_MB, na.rm = TRUE) / 1024),
    as.character(nrow(series_metadata)),
    "FALSE",
    "FALSE",
    "TRUE"
  ),
  stringsAsFactors = FALSE
)
summary_path <- file.path(
  result_dir, "PRJNA1206333_metadata_audit_summary_v1.csv"
)
safe_write_csv(summary_table, summary_path)

download_manifest <- data.frame(
  accession = bioproject_id,
  file_type = "NCBI SRA RunInfo metadata",
  url = runinfo_url,
  destination = normalizePath(
    runinfo_path, winslash = "/", mustWork = TRUE
  ),
  size_bytes = file.info(runinfo_path)$size,
  md5 = unname(tools::md5sum(runinfo_path)),
  raw_sequence_downloaded = FALSE,
  stringsAsFactors = FALSE
)
manifest_path <- file.path(
  metadata_dir, "PRJNA1206333_metadata_download_manifest_v1.csv"
)
safe_write_csv(download_manifest, manifest_path)

decision_lines <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate",
  "- Origin Date: 2026-08-24",
  "- Verification Status: ANALYZED",
  "- Version Label: GSE273293_PRJNA1206333_sample_recovery_09A_v1",
  "",
  "## Step 09A PRJNA1206333 sample-recovery metadata audit",
  "",
  "### Scope",
  "",
  "- This step downloads NCBI RunInfo metadata only.",
  "- No SRA, FASTQ, BAM or expression matrix is downloaded.",
  "- No expression hypothesis is tested.",
  "",
  "### Metadata result",
  "",
  paste0("- Metadata gate passed: ", metadata_gate_passed, "."),
  paste0("- Runs / experiments / SRA samples / BioSamples: ",
         count_unique(runinfo$Run), " / ",
         count_unique(runinfo$Experiment), " / ",
         count_unique(runinfo$Sample), " / ",
         count_unique(runinfo$BioSample), "."),
  paste0("- Clinical labels: 10 GMC and 4 nonfibrotic controls."),
  paste0("- Estimated total raw archive size: ",
         sprintf("%.2f", total_size_gib), " GiB."),
  "- Every raw run is larger than 50 MB; raw download is not part of Step 09A.",
  "",
  "### Local GSE273293 result",
  "",
  paste0("- Local GEO series-metadata rows: ", nrow(series_metadata), "."),
  paste0("- Local supplementary TAR sample prefixes: ",
         paste(tar_sample_prefixes, collapse = "; "), "."),
  "- The current local GEO files do not support 14-sample donor-level inference.",
  "",
  "### Gate decision",
  "",
  "- Gate A metadata component: PASS.",
  "- Gate A processed-expression component: BLOCKED/UNKNOWN.",
  "- Full raw-data download: DO NOT PROCEED.",
  "- Proceed to Step 09B: search for a sample-resolved processed object,",
  "  barcode-to-sample map, or author-supplied processed matrices.",
  "- If Step 09B fails, retain GSE273293 as pooled descriptive deep-fascia",
  "  evidence and do not represent it as a 14-donor expression analysis.",
  "",
  "### Evidence boundary",
  "",
  "- Fourteen accession records establish sample availability, not processed",
  "  cell-level donor labels in the local matrix.",
  "- Library names support the 10 GMC + 4 control design but do not by",
  "  themselves validate every clinical variable reported in the paper.",
  "- No disease effect, mechanotransduction effect or causal claim is tested."
)
decision_path <- file.path(
  result_dir, "GSE273293_PRJNA1206333_sample_recovery_decision_v1.md"
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

session_info_path <- file.path(
  log_dir, "09A_PRJNA1206333_metadata_audit_sessionInfo.txt"
)
writeLines(capture.output(sessionInfo()), session_info_path)

message("Step 09A PRJNA1206333 metadata audit completed.")
message("Metadata gate passed: ", metadata_gate_passed)
message("Runs / experiments / BioSamples: 14 / 14 / 14")
message("Clinical labels: 10 GMC + 4 controls")
message("Estimated total raw size: ", sprintf("%.2f", total_size_gib), " GiB")
message("Full raw download authorized: FALSE")
message("Proceed to Step 09B processed-data recovery: TRUE")
message("Sample map: ", sample_map_path)
message("Integrity checks: ", checks_path)
message("Local structure: ", local_structure_path)
message("Summary: ", summary_path)
message("Decision: ", decision_path)

