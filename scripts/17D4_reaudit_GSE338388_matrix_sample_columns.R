options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

# Step 17D4: GSE338388 processed-matrix sample-column re-audit.
#
# The corrected Step 17D2 design is fixed as TGF-beta exposure (+/-) x TEAD
# inhibition (+/-), with three replicates in each cell. This step downloads
# only the official processed expression matrix candidate, identifies its
# sample columns, and audits frozen candidate/module coverage. It performs no
# expression test and does not download FASTQ/RAW/H5AD data.

project_dir <- "."
gse_id <- "GSE338388"

metadata_dir <- file.path(
  project_dir, "data", "metadata", "independent_sources", gse_id
)
processed_dir <- file.path(
  project_dir, "data", "processed", "independent_sources", gse_id
)
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D4_S2_GSE338388_matrix_reaudit"
)
for (x in c(metadata_dir, processed_dir, result_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

registry_path <- file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
)
candidate_path <- file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
)
metadata_path <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D_S2_GSE338388_design_audit",
  paste0(gse_id, "_sample_metadata_snapshot_v1.csv")
)
design_path <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D_S2_GSE338388_design_audit",
  paste0(gse_id, "_design_factor_reconstruction_v2.csv")
)
inventory_path <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D_S2_GSE338388_design_audit",
  paste0(gse_id, "_supplementary_file_inventory_v1.csv")
)

required_inputs <- c(
  registry_path, candidate_path, metadata_path, design_path, inventory_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 17D4 input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

normalize_gene_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub("\\.[0-9]+$", "", x, perl = TRUE)
}

normalize_key <- function(x) {
  toupper(gsub("[^A-Z0-9]", "", as.character(x)))
}

registry <- read.csv(registry_path, check.names = FALSE)
candidate_panel <- read.csv(candidate_path, check.names = FALSE)
sample_metadata <- read.csv(metadata_path, check.names = FALSE)
design <- read.csv(design_path, check.names = FALSE)
inventory <- read.csv(inventory_path, check.names = FALSE)

registry_rows <- lapply(seq_len(nrow(registry)), function(i) {
  data.frame(
    module = as.character(registry$module[[i]]),
    role = as.character(registry$role[[i]]),
    gene = trimws(strsplit(as.character(registry$genes[[i]]), ";", fixed = TRUE)[[1L]]),
    stringsAsFactors = FALSE
  )
})
registry_genes <- unique(unlist(lapply(registry_rows, function(x) x$gene), use.names = FALSE))
candidate_genes <- unique(trimws(as.character(candidate_panel$gene)))

candidate_matrix <- inventory[
  inventory$auto_download_authorized %in% TRUE &
    inventory$is_processed_matrix_candidate %in% TRUE &
    !inventory$is_raw_archive %in% TRUE,
  , drop = FALSE
]
if (nrow(candidate_matrix) == 0L) {
  stop("No official processed matrix candidate was found in the supplementary inventory.")
}

# Prefer a normalized expression matrix over any auxiliary processed file.
candidate_matrix$priority <- ifelse(
  grepl("expression|matrix|normalized|count", candidate_matrix$filename,
        ignore.case = TRUE, perl = TRUE), 1L, 2L
)
candidate_matrix <- candidate_matrix[
  order(candidate_matrix$priority, candidate_matrix$filename), , drop = FALSE
]
matrix_filename <- candidate_matrix$filename[[1L]]
matrix_url <- candidate_matrix$direct_url[[1L]]
matrix_path <- file.path(processed_dir, matrix_filename)

download_status <- "existing"
download_error <- ""
if (!file.exists(matrix_path) || file.info(matrix_path)$size <= 0) {
  download_status <- tryCatch({
    message("Downloading official processed matrix: ", matrix_filename)
    utils::download.file(
      matrix_url, matrix_path, mode = "wb", quiet = FALSE, method = "auto"
    )
    if (!file.exists(matrix_path) || file.info(matrix_path)$size <= 0) {
      stop("Downloaded matrix is absent or empty.")
    }
    "downloaded"
  }, error = function(e) {
    download_error <<- conditionMessage(e)
    "download_failed"
  })
}
if (download_status == "download_failed") {
  stop("Processed matrix download failed: ", download_error)
}

matrix_connection <- if (grepl("\\.gz$", matrix_path, ignore.case = TRUE)) {
  gzfile(matrix_path, open = "rt")
} else {
  file(matrix_path, open = "rt")
}
expression_matrix <- tryCatch(
  read.csv(
    matrix_connection,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  ),
  error = function(e) {
    try(close(matrix_connection), silent = TRUE)
    stop("Could not parse the processed expression matrix: ", conditionMessage(e))
  }
)
close(matrix_connection)

if (nrow(expression_matrix) < 100L || ncol(expression_matrix) < 5L) {
  stop(
    "The processed file is too small to be an expression matrix: ",
    nrow(expression_matrix), " rows x ", ncol(expression_matrix), " columns."
  )
}

# The official processed matrix contains annotation columns plus exactly twelve
# biological sample columns named Yap164-1 ... Yap164-12. Do not infer sample
# columns from numeric coercion: the annotation column `Gene ID` is numeric and
# would otherwise be incorrectly counted as a thirteenth sample.
sample_columns <- grep(
  "^Yap164-[0-9]+$", names(expression_matrix), value = TRUE
)
sample_column_numbers <- as.integer(sub("^Yap164-", "", sample_columns))
if (
  length(sample_columns) != 12L ||
    !identical(sort(sample_column_numbers), seq_len(12L))
) {
  stop(
    "Could not resolve the expected Yap164-1 ... Yap164-12 sample columns. ",
    "Detected: ", paste(sample_columns, collapse = ", ")
  )
}

gene_candidate_columns <- setdiff(names(expression_matrix), sample_columns)
if (length(gene_candidate_columns) == 0L) {
  stop("No annotation column remains after resolving the twelve sample columns.")
}

candidate_coverage_by_column <- vapply(gene_candidate_columns, function(column) {
  ids <- normalize_gene_id(expression_matrix[[column]])
  sum(normalize_gene_id(candidate_genes) %in% ids)
}, numeric(1))
gene_column <- gene_candidate_columns[[which.max(candidate_coverage_by_column)]]
matrix_gene_ids <- normalize_gene_id(expression_matrix[[gene_column]])
matrix_gene_ids <- matrix_gene_ids[nzchar(matrix_gene_ids)]

candidate_gene_ids <- normalize_gene_id(candidate_genes)
registry_gene_ids <- normalize_gene_id(registry_genes)
candidate_present <- candidate_gene_ids %in% matrix_gene_ids
registry_present <- registry_gene_ids %in% matrix_gene_ids

gene_coverage <- data.frame(
  source = gse_id,
  gene = registry_genes,
  normalized_gene = registry_gene_ids,
  is_frozen_candidate = registry_genes %in% candidate_genes,
  present_in_processed_matrix = registry_present,
  stringsAsFactors = FALSE
)
safe_write_csv(
  gene_coverage,
  file.path(result_dir, paste0(gse_id, "_frozen_gene_coverage_v2.csv"))
)

module_coverage <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
  module_genes <- trimws(strsplit(as.character(registry$genes[[i]]), ";", fixed = TRUE)[[1L]])
  module_ids <- normalize_gene_id(module_genes)
  present <- module_ids %in% matrix_gene_ids
  data.frame(
    source = gse_id,
    module = as.character(registry$module[[i]]),
    role = as.character(registry$role[[i]]),
    genes_total = length(module_genes),
    genes_present = sum(present),
    coverage_fraction = mean(present),
    core_module = as.character(registry$module[[i]]) %in% c(
      "integrin_focal_adhesion", "actomyosin_rho", "hippo_yap_taz"
    ),
    stringsAsFactors = FALSE
  )
}))
module_coverage$core_module_gate <- ifelse(
  module_coverage$core_module,
  module_coverage$coverage_fraction >= 0.80,
  NA
)
safe_write_csv(
  module_coverage,
  file.path(result_dir, paste0(gse_id, "_module_coverage_v2.csv"))
)

# Map Yap164-n to the official GEO sample metadata by the Yap164 number
# embedded in the validated description (for example Yap164-1_S1_L008).
sample_metadata$matrix_sample_number <- suppressWarnings(as.integer(sub(
  ".*Yap164-([0-9]+).*", "\\1", sample_metadata$description, perl = TRUE
)))
if (
  anyNA(sample_metadata$matrix_sample_number) ||
    !identical(sort(sample_metadata$matrix_sample_number), seq_len(12L))
) {
  stop(
    "Official sample descriptions do not contain a unique Yap164-1 ... Yap164-12 mapping."
  )
}

sample_column_rows <- lapply(seq_along(sample_columns), function(i) {
  column <- sample_columns[[i]]
  sample_number <- sample_column_numbers[[i]]
  exact <- which(sample_metadata$matrix_sample_number == sample_number)
  if (length(exact) != 1L) {
    return(data.frame(
      matrix_column = column,
      matrix_sample_number = sample_number,
      library_name = "unresolved",
      gsm_accession = "unresolved",
      tgfb_condition = "unresolved",
      tead_condition = "unresolved",
      match_status = "unresolved",
      mapping_rule = "Yap164 suffix to official description",
      stringsAsFactors = FALSE
    ))
  }
  design_row <- match(
    sample_metadata$gsm_accession[[exact]], design$gsm_accession
  )
  if (is.na(design_row)) {
    return(data.frame(
      matrix_column = column,
      matrix_sample_number = sample_number,
      library_name = sample_metadata$description[[exact]],
      gsm_accession = sample_metadata$gsm_accession[[exact]],
      tgfb_condition = "unresolved",
      tead_condition = "unresolved",
      match_status = "design_unresolved",
      mapping_rule = "Yap164 suffix to official description",
      stringsAsFactors = FALSE
    ))
  }
  library_name <- trimws(sub(
    ".*Library name:[[:space:]]*([^|]+).*",
    "\\1",
    sample_metadata$description[[exact]],
    perl = TRUE
  ))
  data.frame(
    matrix_column = column,
    matrix_sample_number = sample_number,
    library_name = library_name,
    gsm_accession = sample_metadata$gsm_accession[[exact]],
    tgfb_condition = design$tgfb_condition[[design_row]],
    tead_condition = design$tead_condition[[design_row]],
    match_status = "matched",
    mapping_rule = "Yap164 suffix to official description",
    stringsAsFactors = FALSE
  )
})
sample_column_map <- do.call(rbind, sample_column_rows)
safe_write_csv(
  sample_column_map,
  file.path(result_dir, paste0(gse_id, "_sample_column_mapping_v2.csv"))
)

condition_counts <- as.data.frame(table(
  sample_column_map$tgfb_condition,
  sample_column_map$tead_condition,
  useNA = "ifany"
), stringsAsFactors = FALSE)
names(condition_counts) <- c("tgfb_condition", "tead_condition", "matrix_columns")
safe_write_csv(
  condition_counts,
  file.path(result_dir, paste0(gse_id, "_matrix_condition_counts_v2.csv"))
)

sample_mapping_gate <- nrow(sample_column_map) == 12L &&
  all(sample_column_map$match_status == "matched") &&
  all(sample_column_map$tgfb_condition %in% c("not_exposed", "exposed")) &&
  all(sample_column_map$tead_condition %in% c("tead_not_inhibited", "tead_inhibited"))
complete_2x2 <- if (sample_mapping_gate) {
  all(condition_counts$matrix_columns > 0L)
} else FALSE
three_replicates_each <- if (sample_mapping_gate) {
  all(condition_counts$matrix_columns == 3L)
} else FALSE
candidate_coverage_gate <- all(candidate_present)
core_coverage_gate <- all(
  module_coverage$coverage_fraction[module_coverage$core_module] >= 0.80
)

decision <- if (
  download_status %in% c("existing", "downloaded") &&
    sample_mapping_gate && complete_2x2 && three_replicates_each &&
    candidate_coverage_gate && core_coverage_gate
) {
  "PASS_TO_S2_EXPRESSION_ANALYSIS"
} else {
  "HOLD_MATRIX_OR_COVERAGE_REVIEW"
}

matrix_inventory <- data.frame(
  source = gse_id,
  filename = matrix_filename,
  local_path = matrix_path,
  direct_url = matrix_url,
  download_status = download_status,
  bytes = file.info(matrix_path)$size,
  rows = nrow(expression_matrix),
  columns = ncol(expression_matrix),
  resolved_gene_column = gene_column,
  resolved_sample_columns = length(sample_columns),
  annotation_columns = paste(setdiff(names(expression_matrix), sample_columns), collapse = ";"),
  frozen_candidates_present = sum(candidate_present),
  frozen_candidates_total = length(candidate_genes),
  registry_genes_present = sum(registry_present),
  registry_genes_total = length(registry_genes),
  stringsAsFactors = FALSE
)
safe_write_csv(
  matrix_inventory,
  file.path(result_dir, paste0(gse_id, "_processed_matrix_inventory_v2.csv"))
)

decision_lines <- c(
  paste0("# Step 17D4 ", gse_id, " processed-matrix sample-column re-audit"),
  "",
  "## Decision",
  "",
  paste0("- Gate 17D4: **", decision, "**."),
  "- This step audited the official processed matrix and performed no expression hypothesis test.",
  "- No FASTQ, RAW archive, H5AD, or unprocessed sequencing file was downloaded.",
  "",
  "## Matrix and design checks",
  "",
  paste0("- Matrix file: ", matrix_filename, " (", download_status, ")."),
  paste0("- Matrix dimensions: ", nrow(expression_matrix), " rows x ", ncol(expression_matrix), " columns."),
  paste0("- Resolved gene identifier column: ", gene_column, "."),
  paste0("- Resolved Yap164 sample columns: ", length(sample_columns), "."),
  paste0("- Sample-column mapping gate: ", sample_mapping_gate, "."),
  paste0("- Complete TGF-beta x TEAD 2x2 matrix design: ", complete_2x2, "."),
  paste0("- Three matrix replicates in each combination: ", three_replicates_each, "."),
  paste0("- Frozen candidate coverage: ", sum(candidate_present), "/", length(candidate_genes), "."),
  paste0("- Registry coverage: ", sum(registry_present), "/", length(registry_genes), "."),
  paste0("- Core-module coverage gate (>=80% per core module): ", core_coverage_gate, "."),
  "",
  "## S2 interpretation boundary",
  "",
  "- A passing audit authorizes the next expression analysis only.",
  "- The next analysis must report TGF-beta main effect, TEAD-inhibition main effect, and interaction separately.",
  "- Because GSE338388 contains no mechanical loading or stiffness factor, it remains regulatory-axis cross-validation rather than mechanical causality.",
  "- Missing genes are not replaced and the overall project evidence grade remains CAUTION.",
  "",
  "## Material Passport",
  "",
  "- Source: official NCBI GEO GSE338388 processed expression matrix.",
  "- Transformation: matrix parsing, sample-column mapping, frozen candidate/module coverage audit.",
  "- New source type downloaded: one processed matrix candidate only.",
  "- Next step: S2 frozen-contract two-factor expression analysis, if the gate is PASS."
)
decision_path <- file.path(
  result_dir, paste0(gse_id, "_Step17D4_matrix_reaudit_decision_v1.md")
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17D4 ", gse_id, " processed-matrix sample-column re-audit completed.")
message("Matrix dimensions: ", nrow(expression_matrix), " x ", ncol(expression_matrix))
message("Sample mapping gate: ", sample_mapping_gate)
message("Complete 2x2 design: ", complete_2x2)
message("Three replicates each: ", three_replicates_each)
message("Frozen candidate coverage: ", sum(candidate_present), "/", length(candidate_genes))
message("Gate 17D4: ", decision)
message("Decision: ", decision_path)
