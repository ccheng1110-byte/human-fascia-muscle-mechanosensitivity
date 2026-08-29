options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

# Step 17F: S3 source/provenance feasibility audit.
# GSE123100 is the planned stiffness dose-response source. GSE276045 is the
# planned WI-38 stiffness x proliferation cross-check. This step reads only
# official GEO metadata and the supplementary-file directory listing. It does
# not download expression matrices, RAW archives, FASTQ, H5AD, or RDS files.

project_dir <- "."
gse_ids <- c("GSE123100", "GSE276045")
runtime_dir <- file.path(
  project_dir, ".runtime", "step17f_S3_source_audit"
)
metadata_cache_dir <- file.path(runtime_dir, "geo_metadata_cache")
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17F_S3_source_feasibility_audit"
)
for (x in c(runtime_dir, metadata_cache_dir, result_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Install GEOquery once with: BiocManager::install('GEOquery')")
}

safe_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
flatten_meta <- function(value) {
  if (is.null(value) || length(value) == 0L) return("")
  paste(as.character(unlist(value, use.names = FALSE)), collapse = " | ")
}
meta_value <- function(meta, key) {
  if (is.null(meta[[key]])) return("")
  flatten_meta(meta[[key]])
}
characteristics_to_long <- function(gsm_id, values) {
  values <- as.character(unlist(values, use.names = FALSE))
  values <- trimws(values[nzchar(trimws(values))])
  if (length(values) == 0L) return(NULL)
  has_colon <- grepl(":", values, fixed = TRUE)
  key <- ifelse(has_colon, trimws(sub(":.*$", "", values)), "unparsed")
  value <- ifelse(
    has_colon,
    trimws(sub("^[^:]+:[[:space:]]*", "", values)),
    values
  )
  data.frame(
    gsm_accession = gsm_id,
    characteristic_key = key,
    characteristic_value = value,
    original_text = values,
    stringsAsFactors = FALSE
  )
}
extract_hrefs <- function(lines) {
  matches <- unlist(regmatches(
    lines, gregexpr('href="[^"]+"', lines, perl = TRUE)
  ))
  if (length(matches) == 0L) return(character())
  hrefs <- sub('^href="', "", matches)
  sub('"$', "", hrefs)
}
extract_numeric_stiffness <- function(text) {
  x <- tolower(text)
  hits <- unlist(regmatches(
    x,
    gregexpr(
      "\\b[0-9]+(?:\\.[0-9]+)?[[:space:]]*(?:kpa|mpa|pa|dyn/cm2|dyne/cm2)\\b",
      x, perl = TRUE
    )
  ))
  unique(trimws(hits[nzchar(trimws(hits))]))
}
resolve_stiffness <- function(text) {
  numeric_hits <- extract_numeric_stiffness(text)
  if (length(numeric_hits) > 0L) return(paste(numeric_hits, collapse = " | "))
  x <- tolower(text)
  labels <- unique(unlist(regmatches(
    x,
    gregexpr("\\b(?:very[ -])?soft\\b|\\b(?:intermediate|moderate)\\b|\\bstiff\\b|\\brigid\\w*\\b|\\bcompliant\\b|\\belastic\\w*\\b", x, perl = TRUE)
  )))
  labels <- trimws(labels[nzchar(trimws(labels))])
  if (length(labels) > 0L) paste(labels, collapse = " | ") else "unresolved"
}
resolve_cell_model <- function(text) {
  x <- toupper(text)
  if (grepl("WI[ -]?38", x, perl = TRUE)) return("WI-38")
  if (grepl("TRABECULAR|TRABECULAR MESHWORK|HTM", x, perl = TRUE)) return("trabecular-meshwork-like")
  if (grepl("FIBROBLAST", x, perl = TRUE)) return("fibroblast")
  if (grepl("CELL LINE", x, fixed = TRUE)) return("cell-line-unresolved")
  "unresolved"
}
resolve_proliferation <- function(text) {
  x <- tolower(text)
  if (grepl("quiescen|growth[ -]?arrest|serum[ -]?starv|non[ -]?proliferat", x, perl = TRUE)) return("low_or_arrested")
  if (grepl("proliferat|cycling|cycling[ -]?cell|cell[ -]?cycle|mitotic", x, perl = TRUE)) return("proliferating")
  if (grepl("senesc", x, perl = TRUE)) return("senescent")
  "unresolved"
}
resolve_batch <- function(text) {
  x <- tolower(text)
  hit <- regmatches(x, regexpr("\\b(?:batch|lot|run)[ _:-]*[a-z0-9.-]+", x, perl = TRUE))
  if (length(hit) == 0L || identical(hit, "")) "unresolved" else hit
}
resolve_biological_unit <- function(text, gsm_id) {
  x <- toupper(text)
  hit <- regmatches(x, regexpr("\\b(?:DONOR|SUBJECT|PATIENT|INDIVIDUAL|SAMPLE|REPLICATE)[ _:-]*[A-Z0-9.-]+", x, perl = TRUE))
  if (length(hit) > 0L && nzchar(hit)) return(hit)
  if (grepl("WI[ -]?38", x, perl = TRUE)) return("WI-38_cell_line")
  if (grepl("TRABECULAR|HTM", x, perl = TRUE)) return("HTM_cell_source_unresolved")
  gsm_id
}

audit_one_source <- function(gse_id) {
  raw_metadata_dir <- file.path(
    project_dir, "data", "metadata", "independent_sources", gse_id
  )
  source_result_dir <- file.path(result_dir, gse_id)
  dir.create(raw_metadata_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(source_result_dir, recursive = TRUE, showWarnings = FALSE)

  message("Downloading/auditing official GEO metadata only: ", gse_id)
  gse <- tryCatch(
    GEOquery::getGEO(
      gse_id, GSEMatrix = FALSE, getGPL = FALSE, destdir = metadata_cache_dir
    ),
    error = function(e) e
  )
  if (inherits(gse, "error")) {
    return(data.frame(
      gse_accession = gse_id, metadata_gate = "FAILED_METADATA_RETRIEVAL",
      samples = 0L, stiffness_resolved_fraction = NA_real_,
      stiffness_level_count = NA_integer_,
      cell_model_summary = "unresolved", proliferation_resolved_fraction = NA_real_,
      processed_candidate_count = 0L, decision = "HOLD_METADATA_RETRIEVAL",
      error = conditionMessage(gse), stringsAsFactors = FALSE
    ))
  }

  gsm_list <- GEOquery::GSMList(gse)
  if (length(gsm_list) == 0L) stop("No GSM samples were returned for ", gse_id)
  series_meta <- GEOquery::Meta(gse)
  safe_write_csv(
    data.frame(
      field = names(series_meta),
      value = vapply(series_meta, flatten_meta, character(1)),
      stringsAsFactors = FALSE
    ),
    file.path(source_result_dir, paste0(gse_id, "_series_metadata_v1.csv"))
  )

  sample_rows <- list()
  characteristic_rows <- list()
  for (i in seq_along(gsm_list)) {
    gsm_meta <- GEOquery::Meta(gsm_list[[i]])
    gsm_id <- meta_value(gsm_meta, "geo_accession")
    if (!nzchar(gsm_id)) gsm_id <- names(gsm_list)[[i]]
    characteristics <- gsm_meta[["characteristics_ch1"]]
    text <- paste(
      meta_value(gsm_meta, "title"), meta_value(gsm_meta, "source_name_ch1"),
      flatten_meta(characteristics), meta_value(gsm_meta, "description"),
      meta_value(gsm_meta, "treatment_protocol_ch1"), sep = " | "
    )
    characteristic_rows[[i]] <- characteristics_to_long(gsm_id, characteristics)
    sample_rows[[i]] <- data.frame(
      gse_accession = gse_id,
      gsm_accession = gsm_id,
      title = meta_value(gsm_meta, "title"),
      source_name = meta_value(gsm_meta, "source_name_ch1"),
      organism = meta_value(gsm_meta, "organism_ch1"),
      molecule = meta_value(gsm_meta, "molecule"),
      library_strategy = meta_value(gsm_meta, "library_strategy"),
      characteristics = flatten_meta(characteristics),
      description = meta_value(gsm_meta, "description"),
      treatment_protocol = meta_value(gsm_meta, "treatment_protocol_ch1"),
      design_text = text,
      cell_model = resolve_cell_model(text),
      stiffness_text = resolve_stiffness(text),
      proliferation_state = resolve_proliferation(text),
      batch_text = resolve_batch(text),
      biological_unit = resolve_biological_unit(text, gsm_id),
      stringsAsFactors = FALSE
    )
  }
  sample_metadata <- do.call(rbind, sample_rows)
  characteristic_rows <- Filter(Negate(is.null), characteristic_rows)
  characteristics_long <- if (length(characteristic_rows) > 0L) {
    do.call(rbind, characteristic_rows)
  } else {
    data.frame(
      gsm_accession = character(), characteristic_key = character(),
      characteristic_value = character(), original_text = character(),
      stringsAsFactors = FALSE
    )
  }
  safe_write_csv(
    sample_metadata,
    file.path(raw_metadata_dir, paste0(gse_id, "_sample_metadata_snapshot_v1.csv"))
  )
  safe_write_csv(
    sample_metadata,
    file.path(source_result_dir, paste0(gse_id, "_sample_metadata_audit_v1.csv"))
  )
  safe_write_csv(
    characteristics_long,
    file.path(source_result_dir, paste0(gse_id, "_sample_characteristics_long_v1.csv"))
  )

  numeric_stiffness <- unique(unlist(lapply(
    sample_metadata$stiffness_text, extract_numeric_stiffness
  )))
  numeric_stiffness <- numeric_stiffness[nzchar(numeric_stiffness)]
  stiffness_resolved <- sample_metadata$stiffness_text != "unresolved"
  proliferation_resolved <- sample_metadata$proliferation_state != "unresolved"

  numeric_part <- sub("^GSE", "", gse_id)
  series_prefix <- paste0("GSE", sub("[0-9]{3}$", "nnn", numeric_part))
  supplementary_directory_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/", series_prefix, "/", gse_id, "/suppl/"
  )
  listing_error <- ""
  hrefs <- tryCatch({
    lines <- readLines(supplementary_directory_url, warn = FALSE, encoding = "UTF-8")
    x <- extract_hrefs(lines)
    x <- x[!x %in% c("../", "./")]
    x[!grepl("/$", x)]
  }, error = function(e) {
    listing_error <<- conditionMessage(e)
    character()
  })
  if (length(hrefs) == 0L) {
    hrefs <- as.character(unlist(series_meta[["supplementary_file"]], use.names = FALSE))
    hrefs <- hrefs[nzchar(trimws(hrefs))]
  }
  if (length(hrefs) > 0L) {
    direct_urls <- ifelse(
      grepl("^https?://|^ftp://", hrefs, ignore.case = TRUE),
      hrefs, paste0(supplementary_directory_url, hrefs)
    )
    filenames <- basename(URLdecode(sub("[?#].*$", "", hrefs)))
    supplementary_inventory <- data.frame(
      gse_accession = gse_id,
      filename = filenames,
      direct_url = direct_urls,
      is_raw_archive = grepl(
        "_RAW\\.tar$|fastq|\\.fq\\.gz$|\\.bam$|\\.cram$|\\.sra$",
        filenames, ignore.case = TRUE, perl = TRUE
      ),
      is_processed_matrix_candidate = grepl(
        "\\.(txt|tsv|csv)(\\.gz)?$|\\.(rds|rds\\.gz|h5|h5ad|loom|xlsx|xls)$",
        filenames, ignore.case = TRUE, perl = TRUE
      ),
      stringsAsFactors = FALSE
    )
  } else {
    supplementary_inventory <- data.frame(
      gse_accession = gse_id, filename = character(), direct_url = character(),
      is_raw_archive = logical(), is_processed_matrix_candidate = logical(),
      stringsAsFactors = FALSE
    )
  }
  supplementary_inventory$auto_download_authorized <-
    supplementary_inventory$is_processed_matrix_candidate &
    !supplementary_inventory$is_raw_archive
  safe_write_csv(
    supplementary_inventory,
    file.path(source_result_dir, paste0(gse_id, "_supplementary_inventory_v1.csv"))
  )

  processed_count <- sum(supplementary_inventory$auto_download_authorized)
  is_primary <- identical(gse_id, "GSE123100")
  is_wi38 <- mean(grepl("WI-?38", sample_metadata$design_text, ignore.case = TRUE, perl = TRUE)) >= 0.80
  stiffness_fraction <- mean(stiffness_resolved)
  proliferation_fraction <- mean(proliferation_resolved)
  level_count <- if (length(numeric_stiffness) > 0L) {
    length(numeric_stiffness)
  } else {
    length(unique(sample_metadata$stiffness_text[stiffness_resolved]))
  }

  if (is_primary) {
    decision <- if (
      nrow(sample_metadata) > 0L && stiffness_fraction >= 0.80 &&
        level_count >= 3L && processed_count > 0L
    ) {
      "PASS_TO_S3_PRIMARY_MATRIX_AUDIT"
    } else if (nrow(sample_metadata) > 0L && stiffness_fraction > 0 && processed_count > 0L) {
      "PARTIAL_PROCEED_TO_MANUAL_STIFFNESS_REVIEW"
    } else {
      "HOLD_STIFFNESS_OR_MATRIX_FEASIBILITY"
    }
  } else {
    decision <- if (
      nrow(sample_metadata) > 0L && is_wi38 && stiffness_fraction >= 0.80 &&
        proliferation_fraction >= 0.50 && processed_count > 0L
    ) {
      "PASS_TO_S3_WI38_CONFOUNDING_AUDIT"
    } else if (
      nrow(sample_metadata) > 0L && stiffness_fraction > 0 && processed_count > 0L
    ) {
      "PARTIAL_PROCEED_TO_MANUAL_WI38_REVIEW"
    } else {
      "HOLD_WI38_STIFFNESS_OR_MATRIX_FEASIBILITY"
    }
  }

  summary <- data.frame(
    gse_accession = gse_id,
    source_role = ifelse(is_primary, "S3_primary_stiffness_dose_response", "S3_WI38_stiffness_x_proliferation_cross_check"),
    metadata_gate = "PASS_METADATA_RETRIEVED",
    samples = nrow(sample_metadata),
    stiffness_resolved_fraction = stiffness_fraction,
    stiffness_level_count = level_count,
    stiffness_levels_or_labels = paste(
      if (length(numeric_stiffness) > 0L) numeric_stiffness else unique(sample_metadata$stiffness_text[stiffness_resolved]),
      collapse = " | "
    ),
    cell_model_summary = paste(unique(sample_metadata$cell_model), collapse = " | "),
    wi38_fraction = mean(grepl("WI-?38", sample_metadata$design_text, ignore.case = TRUE, perl = TRUE)),
    proliferation_resolved_fraction = proliferation_fraction,
    biological_unit_unique_count = length(unique(sample_metadata$biological_unit)),
    batch_unique_count = length(unique(sample_metadata$batch_text[sample_metadata$batch_text != "unresolved"])),
    processed_candidate_count = processed_count,
    supplementary_listing_error = ifelse(nzchar(listing_error), listing_error, "none"),
    decision = decision,
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary, file.path(source_result_dir, paste0(gse_id, "_source_feasibility_summary_v1.csv")))
  summary
}

summary_rows <- lapply(gse_ids, audit_one_source)
summary_table <- do.call(rbind, summary_rows)
safe_write_csv(summary_table, file.path(result_dir, "Step17F_S3_source_feasibility_summary_v1.csv"))

decision_path <- file.path(result_dir, "Step17F_S3_source_feasibility_decision_v1.md")
decision_lines <- c(
  "# Step 17F S3 source and feasibility audit",
  "",
  "## Scope",
  "",
  "- GSE123100 is audited as the planned cross-tissue stiffness dose-response source.",
  "- GSE276045 is audited as the planned WI-38 stiffness-by-proliferation confounding cross-check.",
  "- Only official GEO metadata and supplementary-file listings were read; no expression matrix or large raw file was downloaded.",
  "- Any passing decision authorizes the next targeted matrix audit only; it does not establish fascia specificity or causality.",
  "",
  "## Source decisions",
  ""
)
for (i in seq_len(nrow(summary_table))) {
  row <- summary_table[i, , drop = FALSE]
  decision_lines <- c(
    decision_lines,
    paste0("### ", row$gse_accession),
    "",
    paste0("- Role: ", row$source_role, "."),
    paste0("- Samples: ", row$samples, "."),
    paste0("- Stiffness resolved fraction: ", signif(row$stiffness_resolved_fraction, 4), "; inferred levels/labels: ", row$stiffness_level_count, "."),
    paste0("- Cell-model summary: ", row$cell_model_summary, "."),
    paste0("- WI-38 fraction: ", signif(row$wi38_fraction, 4), "; proliferation resolved fraction: ", signif(row$proliferation_resolved_fraction, 4), "."),
    paste0("- Processed-file candidates listed: ", row$processed_candidate_count, "."),
    paste0("- Decision: **", row$decision, "**."),
    ""
  )
}
decision_lines <- c(
  decision_lines,
  "## S3 interpretation boundary",
  "",
  "- GSE123100, if advanced, tests the form of a stiffness-related response across tissues; it is not fascia-direct validation.",
  "- GSE276045, if advanced, tests whether a stiffness direction is reproducible across proliferation strata; if proliferation labels are not identifiable, the result remains descriptive.",
  "- A dose-response or cross-check result cannot by itself prove a universal mechanosensitivity mechanism, cell-intrinsic causality, or pain relevance.",
  "- GWAS remains optional and is not required for the S3 computational strengthening branch.",
  "",
  "## Material Passport",
  "",
  "- Inputs: official NCBI GEO metadata records for GSE123100 and GSE276045.",
  "- Transformation: metadata flattening, stiffness/cell-state keyword audit, biological-unit and batch-field inventory, and supplementary-file classification.",
  "- New large data downloaded: none.",
  "- Next step: only sources with a passing or partial feasibility decision proceed to a targeted processed-matrix audit."
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17F S3 source/provenance feasibility audit completed.")
for (i in seq_len(nrow(summary_table))) {
  message(
    summary_table$gse_accession[[i]], ": ", summary_table$decision[[i]],
    "; samples=", summary_table$samples[[i]],
    "; stiffness levels/labels=", summary_table$stiffness_level_count[[i]]
  )
}
message("Decision: ", decision_path)
