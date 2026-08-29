options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

# Step 17G2: corrected targeted S3 processed-matrix audit.
# The matrices are downloaded through R from the official GEO supplementary
# links. This step audits matrix structure, sample-column mapping, and frozen
# gene/module coverage only; it performs no dose-response or hypothesis test.

project_dir <- "."
input_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17F2_S3_corrected_design_audit"
)
source_audit_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17F_S3_source_feasibility_audit"
)
processed_dir <- file.path(
  project_dir, "data", "processed", "independent_sources", "S3"
)
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17G2_S3_processed_matrix_reaudit"
)
for (x in c(processed_dir, result_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

registry_path <- file.path(project_dir, "config", "mechanotransduction_module_registry_v2.csv")
candidate_path <- file.path(project_dir, "config", "frozen_candidate_panel_v2.csv")
selected_path <- file.path(input_dir, "Step17F2_selected_processed_matrix_candidates_v1.csv")
required_inputs <- c(registry_path, candidate_path, selected_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 17G2 input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
normalize_gene_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\.[0-9]+$", "", x, perl = TRUE)
  x[is.na(x)] <- ""
  x
}
normalize_token <- function(x) {
  toupper(gsub("[^A-Z0-9]", "", as.character(x)))
}
read_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  utils::read.delim(
    con, check.names = FALSE, stringsAsFactors = FALSE,
    quote = "\"", comment.char = "", na.strings = c("", "NA", "NaN")
  )
}
extract_gse123100_matrix_key <- function(text) {
  x <- as.character(text)
  hit <- regmatches(
    x,
    regexpr("human-xlw_([AC]_[0-9]+)", x, ignore.case = TRUE, perl = TRUE)
  )
  if (length(hit) > 0L && nzchar(hit)) {
    return(toupper(sub(".*human-xlw_", "", hit, ignore.case = TRUE, perl = TRUE)))
  }
  hit <- regmatches(x, regexpr("[AC]_[0-9]+$", x, ignore.case = TRUE, perl = TRUE))
  if (length(hit) == 0L || identical(hit, "")) "" else toupper(hit)
}
extract_gse276045_condition_key <- function(text) {
  x <- as.character(text)
  hit <- regexec(
    "^(WT|hTERT),[[:space:]]*([0-9]+(?:\\.[0-9]+)?)[[:space:]]*(kPa|GPa),[[:space:]]*timepoint[[:space:]]*([0-9]+),.*replicate[[:space:]]*([0-9]+)$",
    x, ignore.case = TRUE, perl = TRUE
  )
  parts <- regmatches(x, hit)[[1L]]
  if (length(parts) != 6L) return("")
  paste(toupper(parts[2L]), parts[3L], tolower(parts[4L]), parts[5L], parts[6L], sep = "__")
}
extract_gse276045_column_key <- function(text) {
  x <- as.character(text)
  hit <- regexec(
    "(WT|hTERT)[^0-9]*([0-9]+(?:\\.[0-9]+)?)[[:space:]_-]*(kPa|GPa).*?timepoint[[:space:]_-]*([0-9]+).*?replicate[[:space:]_-]*([0-9]+)",
    x, ignore.case = TRUE, perl = TRUE
  )
  parts <- regmatches(x, hit)[[1L]]
  if (length(parts) == 6L) {
    return(paste(toupper(parts[2L]), parts[3L], tolower(parts[4L]), parts[5L], parts[6L], sep = "__"))
  }
  compact <- regexec(
    "^(WT|hTERT)_([0-9]+(?:\\.[0-9]+)?)(kPa|GPa)_([0-9]+)_([0-9]+)$",
    x, ignore.case = TRUE, perl = TRUE
  )
  compact_parts <- regmatches(x, compact)[[1L]]
  if (length(compact_parts) != 6L) return("")
  paste(
    toupper(compact_parts[2L]), compact_parts[3L], tolower(compact_parts[4L]),
    compact_parts[5L], compact_parts[6L], sep = "__"
  )
}
extract_gsm <- function(text) {
  hit <- regmatches(as.character(text), regexpr("GSM[0-9]+", as.character(text), ignore.case = TRUE, perl = TRUE))
  if (length(hit) == 0L || identical(hit, "")) "" else toupper(hit)
}

registry <- utils::read.csv(registry_path, check.names = FALSE, stringsAsFactors = FALSE)
candidate_panel <- utils::read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
candidate_panel$gene <- normalize_gene_id(candidate_panel$gene)
candidate_genes <- unique(candidate_panel$gene)
registry_module_genes <- lapply(registry$genes, function(x) {
  normalize_gene_id(trimws(strsplit(as.character(x), ";", fixed = TRUE)[[1L]]))
})
names(registry_module_genes) <- as.character(registry$module)
selected <- utils::read.csv(selected_path, check.names = FALSE, stringsAsFactors = FALSE)

audit_source <- function(gse_id) {
  source_result_dir <- file.path(result_dir, gse_id)
  dir.create(source_result_dir, recursive = TRUE, showWarnings = FALSE)
  row <- selected[selected$gse_accession == gse_id, , drop = FALSE]
  if (nrow(row) != 1L || !nzchar(row$direct_url[[1L]])) {
    stop("No unique selected processed matrix candidate for ", gse_id)
  }
  metadata_path <- file.path(
    source_audit_dir, gse_id, paste0(gse_id, "_sample_metadata_audit_v1.csv")
  )
  if (!file.exists(metadata_path)) stop("Missing source metadata audit: ", metadata_path)
  metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (gse_id == "GSE123100") {
    corrected_path <- file.path(input_dir, "GSE123100_corrected_S3_design_v2.csv")
    if (!file.exists(corrected_path)) stop("Missing corrected GSE123100 design.")
    corrected <- utils::read.csv(corrected_path, check.names = FALSE, stringsAsFactors = FALSE)
    corrected <- corrected[corrected$eligible_for_primary_dose_response, , drop = FALSE]
    meta_keep <- metadata$gsm_accession %in% corrected$gsm_accession
    metadata <- metadata[meta_keep, , drop = FALSE]
    metadata$matrix_key <- vapply(metadata$description, extract_gse123100_matrix_key, character(1))
    expected_role <- "S3_primary_HTM_stiffness_series"
  } else {
    corrected_path <- file.path(input_dir, "GSE276045_corrected_S3_design_v2.csv")
    if (!file.exists(corrected_path)) stop("Missing corrected GSE276045 design.")
    corrected <- utils::read.csv(corrected_path, check.names = FALSE, stringsAsFactors = FALSE)
    metadata <- metadata[metadata$gsm_accession %in% corrected$gsm_accession, , drop = FALSE]
    metadata$matrix_key <- vapply(metadata$title, extract_gse276045_condition_key, character(1))
    expected_role <- "S3_WI38_stiffness_x_cell_model_cross_check"
  }
  if (nrow(metadata) == 0L) stop("No eligible metadata rows remain for ", gse_id)

  matrix_filename <- as.character(row$matrix_filename[[1L]])
  matrix_path <- file.path(processed_dir, matrix_filename)
  download_status <- "existing"
  if (!file.exists(matrix_path) || file.info(matrix_path)$size <= 0) {
    message("Downloading official processed matrix: ", matrix_filename)
    utils::download.file(
      as.character(row$direct_url[[1L]]), matrix_path,
      mode = "wb", quiet = FALSE, method = "auto"
    )
    download_status <- "downloaded"
  }
  if (!file.exists(matrix_path) || file.info(matrix_path)$size <= 0) {
    stop("Downloaded matrix is absent or empty: ", matrix_path)
  }

  message("Reading processed matrix for ", gse_id)
  expression_df <- read_matrix(matrix_path)
  if (nrow(expression_df) < 100L || ncol(expression_df) < 3L) {
    stop("Processed matrix is too small for an expression matrix: ", nrow(expression_df), " x ", ncol(expression_df))
  }

  matrix_columns <- names(expression_df)
  mapping_rows <- lapply(matrix_columns, function(column) {
    gsm <- extract_gsm(column)
    if (nzchar(gsm) && gsm %in% metadata$gsm_accession) {
      exact <- which(metadata$gsm_accession == gsm)
      return(data.frame(
        matrix_column = column, matrix_key = metadata$matrix_key[[exact[1L]]],
        gsm_accession = gsm, match_status = "matched_by_GSM", stringsAsFactors = FALSE
      ))
    }
    exact_title <- which(metadata$title == column)
    if (length(exact_title) == 1L) {
      return(data.frame(
        matrix_column = column, matrix_key = metadata$matrix_key[[exact_title]],
        gsm_accession = metadata$gsm_accession[[exact_title]],
        match_status = "matched_by_official_title", stringsAsFactors = FALSE
      ))
    }
    key <- if (gse_id == "GSE123100") {
      extract_gse123100_matrix_key(column)
    } else {
      extract_gse276045_column_key(column)
    }
    if (!nzchar(key)) {
      return(data.frame(
        matrix_column = column, matrix_key = "", gsm_accession = "",
        match_status = "unresolved", stringsAsFactors = FALSE
      ))
    }
    exact <- which(metadata$matrix_key == key)
    if (length(exact) == 1L) {
      return(data.frame(
        matrix_column = column, matrix_key = key,
        gsm_accession = metadata$gsm_accession[[exact]],
        match_status = "matched_by_design_key", stringsAsFactors = FALSE
      ))
    }
    data.frame(
      matrix_column = column, matrix_key = key,
      gsm_accession = if (length(exact) > 1L) paste(metadata$gsm_accession[exact], collapse = ";") else "",
      match_status = if (length(exact) > 1L) "ambiguous" else "unresolved",
      stringsAsFactors = FALSE
    )
  })
  mapping <- do.call(rbind, mapping_rows)

  numeric_fraction <- vapply(expression_df, function(x) {
    values <- suppressWarnings(as.numeric(as.character(x)))
    mean(is.finite(values))
  }, numeric(1))
  annotation_like <- grepl(
    "gene|symbol|name|id|description|alias|chrom|chr|transcript|biotype",
    names(expression_df), ignore.case = TRUE, perl = TRUE
  )
  mapped_sample_columns <- mapping$matrix_column[mapping$match_status %in% c("matched_by_GSM", "matched_by_design_key")]
  numeric_sample_columns <- names(expression_df)[numeric_fraction >= 0.90 & !annotation_like]
  sample_columns <- unique(c(mapped_sample_columns, numeric_sample_columns))
  if (length(sample_columns) == 0L) {
    stop("No plausible numeric or design-mapped sample columns detected for ", gse_id)
  }
  gene_candidates <- setdiff(names(expression_df), sample_columns)
  if (length(gene_candidates) == 0L) stop("No gene annotation column remains for ", gse_id)
  candidate_coverage_by_column <- vapply(gene_candidates, function(column) {
    ids <- normalize_gene_id(expression_df[[column]])
    sum(candidate_genes %in% ids)
  }, numeric(1))
  gene_column <- gene_candidates[[which.max(candidate_coverage_by_column)]]
  # Remove numeric annotation columns that were not sample columns after the
  # gene column is selected, then keep only columns that map to eligible samples
  # for downstream analysis authorization.
  mapping <- mapping[mapping$matrix_column %in% sample_columns, , drop = FALSE]
  mapping <- mapping[!duplicated(mapping$matrix_column), , drop = FALSE]
  eligible_mapping <- mapping[mapping$gsm_accession %in% metadata$gsm_accession, , drop = FALSE]
  mapping_complete <- nrow(eligible_mapping) == nrow(metadata) &&
    length(unique(eligible_mapping$gsm_accession)) == nrow(metadata) &&
    all(eligible_mapping$match_status %in% c("matched_by_GSM", "matched_by_design_key"))

  matrix_gene_ids <- normalize_gene_id(expression_df[[gene_column]])
  matrix_gene_ids <- matrix_gene_ids[nzchar(matrix_gene_ids)]
  candidate_present <- candidate_genes %in% matrix_gene_ids
  candidate_coverage <- data.frame(
    source = gse_id, gene = candidate_genes,
    present_in_matrix = candidate_present,
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    candidate_coverage,
    file.path(source_result_dir, paste0(gse_id, "_S3_candidate_coverage_v1.csv"))
  )
  module_coverage <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
    genes <- registry_module_genes[[i]]
    present <- genes %in% matrix_gene_ids
    data.frame(
      source = gse_id,
      module = as.character(registry$module[[i]]),
      role = as.character(registry$role[[i]]),
      genes_total = length(genes), genes_present = sum(present),
      coverage_fraction = mean(present),
      core_module = as.character(registry$module[[i]]) %in% c(
        "integrin_focal_adhesion", "actomyosin_rho", "hippo_yap_taz"
      ),
      stringsAsFactors = FALSE
    )
  }))
  module_coverage$core_module_gate <- ifelse(
    module_coverage$core_module, module_coverage$coverage_fraction >= 0.80, NA
  )
  safe_write_csv(
    module_coverage,
    file.path(source_result_dir, paste0(gse_id, "_S3_module_coverage_v1.csv"))
  )
  safe_write_csv(
    mapping,
    file.path(source_result_dir, paste0(gse_id, "_S3_sample_column_mapping_v1.csv"))
  )

  candidate_matrix_gate <- all(candidate_present)
  core_module_gate <- all(module_coverage$coverage_fraction[module_coverage$core_module] >= 0.80)
  if (gse_id == "GSE123100") {
    decision <- if (mapping_complete && candidate_matrix_gate && core_module_gate) {
      "PASS_TO_S3_PRIMARY_DOSE_RESPONSE_ANALYSIS"
    } else if (mapping_complete) {
      "PARTIAL_PROCEED_WITH_MATRIX_GENE_COVERAGE_LIMIT"
    } else {
      "HOLD_SAMPLE_MAPPING_OR_MATRIX_REVIEW"
    }
  } else {
    decision <- if (mapping_complete && candidate_matrix_gate && core_module_gate) {
      "PASS_TO_S3_WI38_CELL_MODEL_STIFFNESS_ANALYSIS"
    } else if (mapping_complete) {
      "PARTIAL_PROCEED_WITH_MATRIX_GENE_COVERAGE_LIMIT"
    } else {
      "HOLD_SAMPLE_MAPPING_OR_MATRIX_REVIEW"
    }
  }
  inventory <- data.frame(
    gse_accession = gse_id,
    analysis_role = expected_role,
    matrix_filename = matrix_filename,
    matrix_path = matrix_path,
    direct_url = as.character(row$direct_url[[1L]]),
    download_status = download_status,
    bytes = file.info(matrix_path)$size,
    rows = nrow(expression_df), columns = ncol(expression_df),
    eligible_metadata_samples = nrow(metadata),
    detected_sample_columns = length(sample_columns),
    mapped_eligible_samples = length(unique(eligible_mapping$gsm_accession)),
    mapping_complete = mapping_complete,
    resolved_gene_column = gene_column,
    frozen_candidate_present = sum(candidate_present),
    frozen_candidate_total = length(candidate_genes),
    core_module_gate = core_module_gate,
    decision = decision,
    stringsAsFactors = FALSE
  )
  safe_write_csv(
    inventory,
    file.path(source_result_dir, paste0(gse_id, "_S3_processed_matrix_inventory_v1.csv"))
  )
  inventory
}

summary_rows <- lapply(c("GSE123100", "GSE276045"), audit_source)
summary <- do.call(rbind, summary_rows)
safe_write_csv(summary, file.path(result_dir, "Step17G2_S3_processed_matrix_reaudit_summary_v1.csv"))

decision_path <- file.path(result_dir, "Step17G2_S3_processed_matrix_reaudit_decision_v1.md")
decision_lines <- c(
  "# Step 17G2 S3 corrected processed-matrix re-audit",
  "",
  "- This step downloaded only the selected official processed matrices through R and performed no biological hypothesis test.",
  "- GSE123100 clinical tissue samples were excluded from the primary HTM stiffness-series mapping.",
  "- GSE276045 WT/hTERT labels remain cell-model conditions; they are not treated as proliferation measurements.",
  "",
  "## Source decisions",
  ""
)
for (i in seq_len(nrow(summary))) {
  decision_lines <- c(
    decision_lines,
    paste0("### ", summary$gse_accession[[i]]),
    "",
    paste0("- Matrix dimensions: ", summary$rows[[i]], " x ", summary$columns[[i]], "."),
    paste0("- Eligible metadata samples: ", summary$eligible_metadata_samples[[i]], "; mapped eligible samples: ", summary$mapped_eligible_samples[[i]], "."),
    paste0("- Resolved gene column: ", summary$resolved_gene_column[[i]], "."),
    paste0("- Frozen candidate coverage: ", summary$frozen_candidate_present[[i]], "/", summary$frozen_candidate_total[[i]], "."),
    paste0("- Core module gate: ", summary$core_module_gate[[i]], "."),
    paste0("- Decision: **", summary$decision[[i]], "**."),
    ""
  )
}
decision_lines <- c(
  decision_lines,
  "## S3 interpretation boundary",
  "",
  "- A passing matrix audit authorizes the corresponding targeted expression analysis only.",
  "- GSE123100 may support a cautious cross-tissue stiffness dose-response-form analysis, not fascia-direct validation or causal proof.",
  "- GSE276045 may support a stiffness direction by WT/hTERT cell-model condition and timepoint; the proliferation-confounding question remains NOT_ESTIMABLE.",
  "- Overall project evidence grade remains CAUTION."
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17G2 S3 corrected processed-matrix re-audit completed.")
for (i in seq_len(nrow(summary))) {
  message(
    summary$gse_accession[[i]], ": ", summary$decision[[i]],
    "; mapped ", summary$mapped_eligible_samples[[i]], "/",
    summary$eligible_metadata_samples[[i]],
    "; candidates ", summary$frozen_candidate_present[[i]], "/",
    summary$frozen_candidate_total[[i]]
  )
}
message("Decision: ", decision_path)
