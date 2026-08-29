options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

# Step 17D: official GEO provenance and experimental-design audit for
# GSE338388. This step does not download expression matrices, does not read
# raw archives, and performs no biological hypothesis test.
#
# The only design that can pass this gate is TGF-beta exposure (+/-) x TEAD
# inhibition (+/-). The dataset has no mechanical loading/stiffness factor, so
# even a successful gate supports regulatory-axis cross-validation only.

project_dir <- "."
gse_id <- "GSE338388"

runtime_dir <- file.path(project_dir, ".runtime", "step17d_GSE338388")
metadata_cache_dir <- file.path(runtime_dir, "geo_metadata_cache")
metadata_dir <- file.path(
  project_dir, "data", "metadata", "independent_sources", gse_id
)
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D_S2_GSE338388_design_audit"
)
for (x in c(runtime_dir, metadata_cache_dir, metadata_dir, result_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

registry_path <- file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
)
candidate_path <- file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
)
required_inputs <- c(registry_path, candidate_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing frozen configuration input(s): ", paste(missing_inputs, collapse = "; "))
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("Install GEOquery once with: BiocManager::install('GEOquery')")
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
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

resolve_tgfb <- function(text) {
  x <- tolower(text)
  present <- grepl(
    "tgf[[:space:]_-]*(beta|b)|tgf-?β|transforming[[:space:]]+growth[[:space:]]+factor",
    x, perl = TRUE
  )
  absent <- grepl(
    "no[[:space:]_-]*tgf|without[[:space:]_-]*tgf|tgf[^|;]*(untreated|negative|none|vehicle|control|0[[:space:]]*(ng|nm|µ?m|%))",
    x, perl = TRUE
  )
  if (present && absent) return("ambiguous")
  if (present) return("exposed")
  if (grepl("vehicle|untreated|control|mock", x, perl = TRUE)) return("not_exposed")
  "unresolved"
}

resolve_tead <- function(text) {
  x <- tolower(text)
  inhibited <- grepl(
    "k[ -]?975|tead[[:space:]_-]*(inhibitor|inhibition|blocked)|yap[ /]?taz[[:space:]_-]*(inhibitor|inhibition)",
    x, perl = TRUE
  )
  uninhibited <- grepl(
    "vehicle|dmso|untreated|no[[:space:]_-]*inhibitor|without[[:space:]_-]*inhibitor|control|mock",
    x, perl = TRUE
  )
  if (inhibited && uninhibited) return("ambiguous")
  if (inhibited) return("tead_inhibited")
  if (uninhibited) return("tead_not_inhibited")
  "unresolved"
}

resolve_cell_model <- function(text) {
  x <- toupper(text)
  hits <- unique(unlist(regmatches(
    x, gregexpr("\\b(?:VFF|VOC|HDF|WI-?38|FIBROBLAST|GM[0-9]{4,})\\b", x, perl = TRUE)
  )))
  hits <- hits[nzchar(hits)]
  if (length(hits) == 0L) "unresolved" else paste(hits, collapse = ";")
}

message("Downloading/auditing official GEO metadata only: ", gse_id)
gse <- GEOquery::getGEO(
  gse_id,
  GSEMatrix = FALSE,
  getGPL = FALSE,
  destdir = metadata_cache_dir
)
gsm_list <- GEOquery::GSMList(gse)
if (length(gsm_list) == 0L) stop("No GSM samples were returned for ", gse_id)

series_meta <- GEOquery::Meta(gse)
safe_write_csv(
  data.frame(
    field = names(series_meta),
    value = vapply(series_meta, flatten_meta, character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, paste0(gse_id, "_series_metadata_v1.csv"))
)

sample_rows <- list()
characteristic_rows <- list()
for (i in seq_along(gsm_list)) {
  gsm <- gsm_list[[i]]
  meta <- GEOquery::Meta(gsm)
  gsm_id <- meta_value(meta, "geo_accession")
  if (!nzchar(gsm_id)) gsm_id <- names(gsm_list)[[i]]
  characteristics <- meta[["characteristics_ch1"]]
  characteristics_text <- flatten_meta(characteristics)
  characteristic_rows[[i]] <- characteristics_to_long(gsm_id, characteristics)
  sample_rows[[i]] <- data.frame(
    gse_accession = gse_id,
    gsm_accession = gsm_id,
    title = meta_value(meta, "title"),
    source_name = meta_value(meta, "source_name_ch1"),
    organism = meta_value(meta, "organism_ch1"),
    molecule = meta_value(meta, "molecule_ch1"),
    library_strategy = meta_value(meta, "library_strategy"),
    characteristics = characteristics_text,
    description = meta_value(meta, "description"),
    treatment_protocol = meta_value(meta, "treatment_protocol_ch1"),
    growth_protocol = meta_value(meta, "growth_protocol_ch1"),
    supplementary_file = meta_value(meta, "supplementary_file"),
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
  file.path(metadata_dir, paste0(gse_id, "_sample_metadata_snapshot_v1.csv"))
)
safe_write_csv(
  sample_metadata,
  file.path(result_dir, paste0(gse_id, "_sample_metadata_snapshot_v1.csv"))
)
safe_write_csv(
  characteristics_long,
  file.path(result_dir, paste0(gse_id, "_sample_characteristics_long_v1.csv"))
)

design_rows <- lapply(seq_len(nrow(sample_metadata)), function(i) {
  design_text <- paste(
    sample_metadata$title[[i]],
    sample_metadata$source_name[[i]],
    sample_metadata$characteristics[[i]],
    sample_metadata$description[[i]],
    sample_metadata$treatment_protocol[[i]],
    sep = " | "
  )
  data.frame(
    gse_accession = gse_id,
    gsm_accession = sample_metadata$gsm_accession[[i]],
    biological_unit = resolve_cell_model(design_text),
    tgfb_condition = resolve_tgfb(design_text),
    tead_condition = resolve_tead(design_text),
    design_text = design_text,
    stringsAsFactors = FALSE
  )
})
design_table <- do.call(rbind, design_rows)
safe_write_csv(
  design_table,
  file.path(result_dir, paste0(gse_id, "_design_factor_reconstruction_v1.csv"))
)

factor_count <- function(field) {
  tab <- as.data.frame(table(design_table[[field]], useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab) <- c("level", "samples")
  tab$factor <- field
  tab[, c("factor", "level", "samples")]
}
factor_counts <- do.call(rbind, lapply(
  c("biological_unit", "tgfb_condition", "tead_condition"),
  factor_count
))
safe_write_csv(
  factor_counts,
  file.path(result_dir, paste0(gse_id, "_design_factor_counts_v1.csv"))
)

numeric_part <- sub("^GSE", "", gse_id)
series_prefix <- paste0("GSE", sub("[0-9]{3}$", "nnn", numeric_part))
supplementary_directory_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  series_prefix, "/", gse_id, "/suppl/"
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
  series_supp <- as.character(unlist(series_meta[["supplementary_file"]], use.names = FALSE))
  series_supp <- series_supp[nzchar(trimws(series_supp))]
  hrefs <- series_supp
}

if (length(hrefs) > 0L) {
  direct_urls <- ifelse(
    grepl("^https?://|^ftp://", hrefs, ignore.case = TRUE),
    hrefs,
    paste0(supplementary_directory_url, hrefs)
  )
  filenames <- basename(URLdecode(sub("[?#].*$", "", hrefs)))
  supplementary_inventory <- data.frame(
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
    filename = character(), direct_url = character(),
    is_raw_archive = logical(), is_processed_matrix_candidate = logical(),
    stringsAsFactors = FALSE
  )
}
supplementary_inventory$auto_download_authorized <-
  supplementary_inventory$is_processed_matrix_candidate &
  !supplementary_inventory$is_raw_archive
safe_write_csv(
  supplementary_inventory,
  file.path(result_dir, paste0(gse_id, "_supplementary_file_inventory_v1.csv"))
)

tgfb_levels <- sort(unique(design_table$tgfb_condition))
tead_levels <- sort(unique(design_table$tead_condition))
valid_tgfb <- setequal(tgfb_levels, c("exposed", "not_exposed"))
valid_tead <- setequal(tead_levels, c("tead_inhibited", "tead_not_inhibited"))
design_table$tgfb_factor <- factor(
  design_table$tgfb_condition,
  levels = c("not_exposed", "exposed")
)
design_table$tead_factor <- factor(
  design_table$tead_condition,
  levels = c("tead_not_inhibited", "tead_inhibited")
)
complete_design_rows <- design_table[
  !is.na(design_table$tgfb_factor) & !is.na(design_table$tead_factor),
  , drop = FALSE
]
combo_table <- as.data.frame(
  table(complete_design_rows$tgfb_factor, complete_design_rows$tead_factor),
  stringsAsFactors = FALSE
)
names(combo_table) <- c("tgfb_condition", "tead_condition", "samples")
complete_2x2 <- nrow(combo_table) == 4L && all(combo_table$samples > 0L)

panel <- read.csv(candidate_path, check.names = FALSE)
panel_genes <- unique(trimws(as.character(panel$gene)))

decision <- if (
  nrow(sample_metadata) == 12L && valid_tgfb && valid_tead && complete_2x2
) {
  "PASS_TO_S2_MATRIX_AUDIT"
} else if (nrow(sample_metadata) > 0L) {
  "HOLD_NEEDS_MANUAL_DESIGN_REVIEW"
} else {
  "HOLD_NO_SAMPLE_METADATA"
}

decision_lines <- c(
  paste0("# Step 17D ", gse_id, " provenance and design audit"),
  "",
  "## Decision",
  "",
  paste0("- Gate 17D: **", decision, "**."),
  "- This step audited official GEO metadata and supplementary-file inventory only.",
  "- No processed expression matrix, FASTQ, RAW archive, or H5AD was downloaded.",
  "- No biological hypothesis test was performed.",
  "",
  "## Locked design requirement",
  "",
  "- Required design: TGFβ exposure (not exposed/exposed) × TEAD inhibition (not inhibited/inhibited).",
  "- This dataset contains no mechanical loading or stiffness factor.",
  "- A passing audit therefore supports regulatory-axis cross-validation only, not mechanical causality or fascia-direct replication.",
  "",
  "## Reconstructed design",
  "",
  paste0("- GEO sample records: ", nrow(sample_metadata), "."),
  paste0("- TGFβ levels detected: ", paste(tgfb_levels, collapse = ", "), "."),
  paste0("- TEAD levels detected: ", paste(tead_levels, collapse = ", "), "."),
  paste0("- Complete 2×2 combination table: ", complete_2x2, "."),
  paste0("- Supplementary files listed: ", nrow(supplementary_inventory), "."),
  paste0("- Listing error, if any: ", ifelse(nzchar(listing_error), listing_error, "none"), "."),
  "",
  "## S2 boundaries",
  "",
  "- The TGFβ main effect, TEAD-inhibition main effect, and interaction must be reported separately in the next step.",
  "- Frozen candidate and module coverage must be audited before expression analysis; missing genes are not replaced.",
  "- A TEAD-associated module response cannot be described as proof that YAP/TAZ independently drives a mechanical program.",
  "- The project evidence grade remains CAUTION regardless of this gate.",
  "",
  "## Material Passport",
  "",
  "- Source: official NCBI GEO record GSE338388.",
  "- Transformation: sample metadata parsing, factor reconstruction, 2×2 design gate, and supplementary inventory.",
  paste0("- Frozen candidate panel recorded: ", length(panel_genes), " genes; expression coverage remains pending the next matrix audit."),
  "- Next step: S2 processed-matrix audit and frozen module/candidate coverage check."
)
decision_path <- file.path(
  result_dir, paste0(gse_id, "_Step17D_provenance_design_decision_v1.md")
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17D ", gse_id, " provenance/design audit completed.")
message("GEO sample records: ", nrow(sample_metadata))
message("TGF-beta levels: ", paste(tgfb_levels, collapse = ", "))
message("TEAD levels: ", paste(tead_levels, collapse = ", "))
message("Complete 2x2 design: ", complete_2x2)
message("Gate 17D: ", decision)
message("Decision: ", decision_path)
