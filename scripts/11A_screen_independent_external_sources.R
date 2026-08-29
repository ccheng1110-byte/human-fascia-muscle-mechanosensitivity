options(stringsAsFactors = FALSE)

# Step 11A: metadata-only screening of non-overlapping external sources.
# No expression matrix, H5AD, or raw sequencing file is downloaded here.

project_dir <- "."
audit_dir <- file.path(
  project_dir, "results", "08_cross_tissue_validation",
  "10A_atlas_provenance_audit"
)
runtime_dir <- file.path(project_dir, ".runtime", "geo_step11a")
cache_dir <- file.path(runtime_dir, "geo_metadata_cache")
result_dir <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11A_metadata_only_source_audit"
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

source_map_path <- file.path(
  audit_dir, "skin_fibroblast_atlas_source_accession_map_v1.csv"
)
if (!file.exists(source_map_path)) {
  stop("Missing corrected Step 10A source map: ", source_map_path)
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Install GEOquery once with: BiocManager::install('GEOquery')")
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

flatten_meta <- function(value) {
  if (is.null(value)) return("")
  paste(as.character(unlist(value, use.names = FALSE)), collapse = " | ")
}

meta_text <- function(meta) {
  if (length(meta) == 0L) return("")
  paste(
    vapply(names(meta), function(key) {
      paste0(key, "=", flatten_meta(meta[[key]]))
    }, character(1)),
    collapse = " | "
  )
}

has_keyword <- function(x, pattern) {
  grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
}

source_map <- read.csv(source_map_path, check.names = FALSE)
source_map$independent_candidate_after_exclusion <- tolower(
  trimws(as.character(source_map$independent_candidate_after_exclusion))
) %in% c("true", "t", "1", "yes")
source_map <- source_map[
  source_map$independent_candidate_after_exclusion &
    !is.na(source_map$source_accession) &
    grepl("^GSE[0-9]+$", source_map$source_accession),
  , drop = FALSE
]
source_map <- source_map[!duplicated(source_map$source_accession), , drop = FALSE]

summary_rows <- list()
sample_rows <- list()

for (i in seq_len(nrow(source_map))) {
  accession <- source_map$source_accession[[i]]
  message("Auditing metadata only: ", accession)
  result <- tryCatch({
    gse <- GEOquery::getGEO(
      accession,
      GSEMatrix = FALSE,
      getGPL = FALSE,
      destdir = cache_dir
    )
    gsm_list <- GEOquery::GSMList(gse)
    if (length(gsm_list) == 0L) {
      stop("No GSM sample objects returned.")
    }

    local_samples <- lapply(seq_along(gsm_list), function(j) {
      gsm <- gsm_list[[j]]
      meta <- GEOquery::Meta(gsm)
      gsm_id <- if (!is.null(meta$geo_accession)) {
        flatten_meta(meta$geo_accession)
      } else {
        names(gsm_list)[[j]]
      }
      sample_text <- meta_text(meta)
      characteristics <- flatten_meta(meta$characteristics_ch1)
      data.frame(
        source_accession = accession,
        gsm_accession = gsm_id,
        title = flatten_meta(meta$title),
        source_name = flatten_meta(meta$source_name_ch1),
        organism = flatten_meta(meta$organism_ch1),
        characteristics = characteristics,
        library_strategy = flatten_meta(meta$library_strategy),
        supplementary_file = flatten_meta(meta$supplementary_file),
        all_metadata_text = sample_text,
        donor_like = has_keyword(
          sample_text,
          "donor|subject|patient|individual|participant|biopsy|specimen|replicate"
        ),
        fibroblast_like = has_keyword(
          sample_text,
          "fibroblast|myofibroblast|mesenchymal|stromal|fascia"
        ),
        single_cell_like = has_keyword(
          sample_text,
          "single.?cell|scrna|scRNA|10x|droplet|cell.?ranger|single cell"
        ),
        disease_or_control_like = has_keyword(
          sample_text,
          "disease|control|healthy|patient|case|lesion|normal|tissue"
        ),
        stringsAsFactors = FALSE
      )
    })
    local_samples <- do.call(rbind, local_samples)
    sample_rows[[length(sample_rows) + 1L]] <- local_samples

    donor_count <- sum(local_samples$donor_like)
    fibroblast_count <- sum(local_samples$fibroblast_like)
    single_cell_count <- sum(local_samples$single_cell_like)
    clinical_count <- sum(local_samples$disease_or_control_like)
    priority_score <-
      3L * as.integer(donor_count > 0L) +
      3L * as.integer(fibroblast_count > 0L) +
      2L * as.integer(single_cell_count > 0L) +
      1L * as.integer(nrow(local_samples) >= 4L) +
      1L * as.integer(clinical_count > 0L)
    priority_class <- if (priority_score >= 8L) {
      "high_priority_manual_object_audit"
    } else if (priority_score >= 5L) {
      "medium_priority_manual_object_audit"
    } else {
      "low_priority_or_insufficient_metadata"
    }

    data.frame(
      source_accession = accession,
      status = "metadata_audit_completed",
      gsm_samples = nrow(local_samples),
      donor_like_samples = donor_count,
      fibroblast_like_samples = fibroblast_count,
      single_cell_like_samples = single_cell_count,
      disease_or_control_like_samples = clinical_count,
      priority_score = priority_score,
      priority_class = priority_class,
      error_message = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      source_accession = accession,
      status = "metadata_audit_failed",
      gsm_samples = NA_integer_,
      donor_like_samples = NA_integer_,
      fibroblast_like_samples = NA_integer_,
      single_cell_like_samples = NA_integer_,
      disease_or_control_like_samples = NA_integer_,
      priority_score = NA_integer_,
      priority_class = "manual_follow_up_required",
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  summary_rows[[length(summary_rows) + 1L]] <- result
}

source_summary <- do.call(rbind, summary_rows)
source_summary <- source_summary[order(
  -ifelse(is.na(source_summary$priority_score), -1L, source_summary$priority_score),
  source_summary$source_accession
), , drop = FALSE]
rownames(source_summary) <- NULL
safe_write_csv(
  source_summary,
  file.path(result_dir, "step11_independent_source_priority_audit_v1.csv")
)

if (length(sample_rows) > 0L) {
  sample_snapshot <- do.call(rbind, sample_rows)
} else {
  sample_snapshot <- data.frame()
}
safe_write_csv(
  sample_snapshot,
  file.path(result_dir, "step11_independent_source_sample_metadata_snapshot_v1.csv")
)

high_priority <- source_summary[
  source_summary$priority_class == "high_priority_manual_object_audit", , drop = FALSE
]
failed_count <- sum(source_summary$status != "metadata_audit_completed")

decision_lines <- c(
  "# Step 11A independent external source screening",
  "",
  "## Scope",
  "",
  "- Only accession-resolved sources not overlapping GSE173252 or PRJNA607098 were screened.",
  "- This step downloaded GEO metadata only; no expression matrix, H5AD, or raw sequencing file was downloaded.",
  "- Priority scores are triage rules, not evidence of replication.",
  "",
  "## Screening result",
  "",
  paste0("- Candidate GSE sources audited: ", nrow(source_summary), "."),
  paste0("- High-priority manual object audits: ", nrow(high_priority), "."),
  paste0("- Metadata audit failures requiring follow-up: ", failed_count, "."),
  "",
  "## Selection rule for Step 11B",
  "",
  "Prioritize sources with sample-level donor/subject fields, explicit fibroblast or myofibroblast labels, and single-cell evidence. Before expression analysis, manually verify the object-level cell-state field, source accession, sample/donor mapping, and availability of the frozen candidate genes.",
  "",
  "## Material Passport",
  "",
  "- Origin: corrected Step 10A source map.",
  "- Transformation: GEO series metadata were flattened into a source-level priority table and a sample-level metadata snapshot.",
  "- No expression values were read or tested.",
  "- Reproducibility: GEO metadata cache is under `.runtime/geo_step11a/geo_metadata_cache`; all exported tables are in the result directory.",
  "- Integrity boundary: keyword-based priority is only a screening aid and requires manual verification before validation."
)
decision_path <- file.path(result_dir, "step11_source_screening_decision_v1.md")
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 11A metadata-only source screening completed.")
message("Priority audit: ", file.path(result_dir, "step11_independent_source_priority_audit_v1.csv"))
message("Sample metadata snapshot: ", file.path(result_dir, "step11_independent_source_sample_metadata_snapshot_v1.csv"))
message("Decision: ", decision_path)
