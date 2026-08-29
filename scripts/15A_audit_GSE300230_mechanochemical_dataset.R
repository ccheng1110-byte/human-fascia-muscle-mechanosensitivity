options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

# Step 15A: official GEO provenance, design, processed-file, and gene-coverage
# audit for GSE300230. This step performs no biological hypothesis test.
# It downloads metadata and processed text matrices only. FASTQ, BAM, H5AD,
# and GEO RAW.tar archives are not downloaded.

project_dir <- "."
gse_id <- "GSE300230"

runtime_dir <- file.path(project_dir, ".runtime", "step15a_GSE300230")
metadata_cache_dir <- file.path(runtime_dir, "geo_metadata_cache")
raw_metadata_dir <- file.path(
  project_dir, "data", "metadata", "mechanochemical_validation", gse_id
)
processed_dir <- file.path(
  project_dir, "data", "processed", "mechanochemical_validation", gse_id
)
result_dir <- file.path(
  project_dir, "results", "12_computational_strengthening",
  "15A_GSE300230_provenance_feasibility_audit"
)
for (x in c(runtime_dir, metadata_cache_dir, raw_metadata_dir, processed_dir, result_dir)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

candidate_panel_path <- file.path(project_dir, "config", "frozen_candidate_panel_v2.csv")
module_registry_path <- file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
)
required_local_inputs <- c(candidate_panel_path, module_registry_path)
missing_local_inputs <- required_local_inputs[!file.exists(required_local_inputs)]
if (length(missing_local_inputs) > 0L) {
  stop("Missing frozen configuration input(s): ", paste(missing_local_inputs, collapse = "; "))
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

first_characteristic_value <- function(characteristics_long, gsm_id, key_pattern) {
  if (nrow(characteristics_long) == 0L) return("")
  hit <- characteristics_long[
    characteristics_long$gsm_accession == gsm_id &
      grepl(key_pattern, characteristics_long$characteristic_key,
        ignore.case = TRUE, perl = TRUE
      ),
    "characteristic_value",
    drop = TRUE
  ]
  hit <- unique(trimws(as.character(hit)))
  hit <- hit[nzchar(hit)]
  if (length(hit) == 0L) "" else paste(hit, collapse = " | ")
}

extract_regex_value <- function(text, pattern) {
  hit <- regmatches(text, regexpr(pattern, text, ignore.case = TRUE, perl = TRUE))
  if (length(hit) == 0L || identical(hit, "")) "" else hit
}

resolve_age_group <- function(text) {
  x <- tolower(text)
  if (grepl("\\byoung\\b|juvenile", x, perl = TRUE)) return("young")
  if (grepl("\\bold\\b|\\baged\\b|elder", x, perl = TRUE)) return("older")
  "unresolved"
}

resolve_mechanical_condition <- function(text) {
  x <- tolower(text)
  relaxed <- grepl(
    "relaxed|untensed|un-tensed|unanchored|un-anchored|without[ _-]*tension|no[ _-]*tension",
    x, perl = TRUE
  )
  tension <- grepl(
    "tensed|tension|anchored|mechanical[ _-]*(force|load)",
    x, perl = TRUE
  )
  if (relaxed && !tension) return("relaxed")
  if (tension && !relaxed) return("tension")
  if (relaxed && tension) return("ambiguous")
  "unresolved"
}

resolve_tgfb_condition <- function(text) {
  x <- tolower(text)
  absent <- grepl(
    "no[ _-]*tgf|without[ _-]*tgf|tgf[^|;]*(untreated|negative|none|vehicle|0)",
    x, perl = TRUE
  )
  present <- grepl("tgf[[:space:]_-]*(beta|b|β)|tgf-?β", x, perl = TRUE)
  if (absent) return("absent")
  if (present) return("present")
  "unresolved"
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
series_metadata <- data.frame(
  field = names(series_meta),
  value = vapply(series_meta, flatten_meta, character(1)),
  stringsAsFactors = FALSE
)
safe_write_csv(
  series_metadata,
  file.path(result_dir, "GSE300230_series_metadata_v1.csv")
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
if (length(characteristic_rows) > 0L) {
  characteristics_long <- do.call(rbind, characteristic_rows)
} else {
  characteristics_long <- data.frame(
    gsm_accession = character(), characteristic_key = character(),
    characteristic_value = character(), original_text = character(),
    stringsAsFactors = FALSE
  )
}

safe_write_csv(
  sample_metadata,
  file.path(raw_metadata_dir, "GSE300230_sample_metadata_snapshot_v1.csv")
)
safe_write_csv(
  sample_metadata,
  file.path(result_dir, "GSE300230_sample_metadata_snapshot_v1.csv")
)
safe_write_csv(
  characteristics_long,
  file.path(result_dir, "GSE300230_sample_characteristics_long_v1.csv")
)

design_rows <- lapply(seq_len(nrow(sample_metadata)), function(i) {
  gsm_id <- sample_metadata$gsm_accession[[i]]
  design_text <- paste(
    sample_metadata$title[[i]],
    sample_metadata$source_name[[i]],
    sample_metadata$characteristics[[i]],
    sample_metadata$description[[i]],
    sep = " | "
  )
  biological_unit <- first_characteristic_value(
    characteristics_long, gsm_id,
    "donor|cell.?line|individual|subject|participant|biospecimen"
  )
  if (!nzchar(biological_unit)) {
    biological_unit <- extract_regex_value(
      design_text, "\\b(?:GM|AG)[0-9]{4,}\\b"
    )
  }
  age_text <- first_characteristic_value(
    characteristics_long, gsm_id, "^age$|age.?group|donor.?age"
  )
  if (!nzchar(age_text)) age_text <- design_text
  data.frame(
    gsm_accession = gsm_id,
    title = sample_metadata$title[[i]],
    biological_unit = ifelse(nzchar(biological_unit), biological_unit, "unresolved"),
    age_group = resolve_age_group(age_text),
    mechanical_condition = resolve_mechanical_condition(design_text),
    tgfb_condition = resolve_tgfb_condition(design_text),
    replicate_label = extract_regex_value(
      design_text, "\\brep(?:licate)?[ _:=.-]*[0-9]+\\b"
    ),
    design_text = design_text,
    stringsAsFactors = FALSE
  )
})
design_table <- do.call(rbind, design_rows)
safe_write_csv(
  design_table,
  file.path(result_dir, "GSE300230_design_factor_reconstruction_v1.csv")
)

factor_count <- function(field) {
  tab <- as.data.frame(table(design_table[[field]], useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab) <- c("level", "samples")
  tab$factor <- field
  tab[, c("factor", "level", "samples")]
}
factor_counts <- do.call(rbind, lapply(
  c("biological_unit", "age_group", "mechanical_condition", "tgfb_condition"),
  factor_count
))
safe_write_csv(
  factor_counts,
  file.path(result_dir, "GSE300230_design_factor_counts_v1.csv")
)

extract_hrefs <- function(lines) {
  matches <- unlist(regmatches(lines, gregexpr('href="[^"]+"', lines, perl = TRUE)))
  if (length(matches) == 0L) return(character())
  hrefs <- sub('^href="', "", matches)
  sub('"$', "", hrefs)
}

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
    is_processed_text_matrix = grepl(
      "\\.(txt|tsv|csv)(\\.gz)?$",
      filenames, ignore.case = TRUE, perl = TRUE
    ),
    is_other_processed_candidate = grepl(
      "\\.(xlsx|xls|rds|rds\\.gz|h5|h5ad|loom)$",
      filenames, ignore.case = TRUE, perl = TRUE
    ),
    stringsAsFactors = FALSE
  )
} else {
  supplementary_inventory <- data.frame(
    filename = character(), direct_url = character(),
    is_raw_archive = logical(), is_processed_text_matrix = logical(),
    is_other_processed_candidate = logical(), stringsAsFactors = FALSE
  )
}
supplementary_inventory$auto_download_authorized <-
  supplementary_inventory$is_processed_text_matrix &
  !supplementary_inventory$is_raw_archive
safe_write_csv(
  supplementary_inventory,
  file.path(result_dir, "GSE300230_supplementary_file_inventory_v1.csv")
)

download_rows <- list()
download_candidates <- supplementary_inventory[
  supplementary_inventory$auto_download_authorized,
  , drop = FALSE
]
if (nrow(download_candidates) > 0L) {
  for (i in seq_len(nrow(download_candidates))) {
    target_path <- file.path(processed_dir, download_candidates$filename[[i]])
    status <- "existing"
    error_message <- ""
    if (!file.exists(target_path) || file.info(target_path)$size <= 0) {
      status <- tryCatch({
        message("Downloading processed matrix candidate: ", download_candidates$filename[[i]])
        utils::download.file(
          download_candidates$direct_url[[i]], target_path,
          mode = "wb", quiet = FALSE, method = "auto"
        )
        if (!file.exists(target_path) || file.info(target_path)$size <= 0) {
          stop("Downloaded file is absent or empty.")
        }
        "downloaded"
      }, error = function(e) {
        error_message <<- conditionMessage(e)
        "download_failed"
      })
    }
    download_rows[[length(download_rows) + 1L]] <- data.frame(
      filename = download_candidates$filename[[i]],
      local_path = target_path,
      status = status,
      bytes = if (file.exists(target_path)) file.info(target_path)$size else NA_real_,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }
}
if (length(download_rows) > 0L) {
  download_manifest <- do.call(rbind, download_rows)
} else {
  download_manifest <- data.frame(
    filename = character(), local_path = character(), status = character(),
    bytes = numeric(), error_message = character(), stringsAsFactors = FALSE
  )
}
safe_write_csv(
  download_manifest,
  file.path(result_dir, "GSE300230_processed_matrix_download_manifest_v1.csv")
)

extract_gene_ids_from_text_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  if (length(lines) < 2L) stop("Matrix has fewer than two text lines.")
  header <- lines[[1L]]
  tab_count <- lengths(regmatches(header, gregexpr("\\t", header, perl = TRUE)))
  comma_count <- lengths(regmatches(header, gregexpr(",", header, fixed = TRUE)))
  separator <- if (tab_count >= comma_count) "\t" else ","
  first_fields <- vapply(strsplit(lines[-1L], separator, fixed = TRUE), function(x) {
    if (length(x) == 0L) "" else x[[1L]]
  }, character(1))
  first_fields <- trimws(gsub('^"|"$', "", first_fields))
  unique(first_fields[nzchar(first_fields)])
}

candidate_panel <- read.csv(candidate_panel_path, check.names = FALSE)
module_registry <- read.csv(module_registry_path, check.names = FALSE)
module_gene_table <- do.call(rbind, lapply(seq_len(nrow(module_registry)), function(i) {
  genes <- trimws(unlist(strsplit(module_registry$genes[[i]], ";", fixed = TRUE)))
  data.frame(
    module = module_registry$module[[i]],
    role = module_registry$role[[i]],
    gene = genes[nzchar(genes)],
    stringsAsFactors = FALSE
  )
}))

# Nomenclature-only aliases. They do not authorize candidate replacement.
alias_map <- data.frame(
  registry_symbol = c("CTGF", "CYR61"),
  accepted_current_symbol = c("CCN2", "CCN1"),
  mapping_type = "HGNC_nomenclature_alias_only",
  stringsAsFactors = FALSE
)
safe_write_csv(
  alias_map,
  file.path(result_dir, "GSE300230_gene_symbol_alias_audit_v1.csv")
)

match_gene <- function(gene, gene_ids) {
  alternatives <- gene
  alias_hit <- alias_map$registry_symbol == gene
  if (any(alias_hit)) {
    alternatives <- unique(c(alternatives, alias_map$accepted_current_symbol[alias_hit]))
  }
  hit <- alternatives[alternatives %in% gene_ids]
  if (length(hit) == 0L) "" else hit[[1L]]
}

coverage_rows <- list()
usable_downloads <- download_manifest[
  download_manifest$status %in% c("downloaded", "existing") &
    file.exists(download_manifest$local_path),
  , drop = FALSE
]
if (nrow(usable_downloads) > 0L) {
  for (i in seq_len(nrow(usable_downloads))) {
    local_path <- usable_downloads$local_path[[i]]
    gene_ids <- tryCatch(
      extract_gene_ids_from_text_matrix(local_path),
      error = function(e) structure(character(), audit_error = conditionMessage(e))
    )
    audit_error <- attr(gene_ids, "audit_error")
    ensembl_fraction <- if (length(gene_ids) > 0L) {
      mean(grepl("^ENSG[0-9]+", gene_ids, ignore.case = TRUE))
    } else {
      NA_real_
    }
    identifier_type <- if (!is.na(ensembl_fraction) && ensembl_fraction > 0.5) {
      "ensembl_requires_annotation"
    } else if (length(gene_ids) > 0L) {
      "gene_symbol_or_mixed"
    } else {
      "unresolved"
    }

    local_coverage <- module_gene_table
    local_coverage$filename <- usable_downloads$filename[[i]]
    local_coverage$identifier_type <- identifier_type
    local_coverage$matched_symbol <- if (identifier_type == "gene_symbol_or_mixed") {
      vapply(local_coverage$gene, match_gene, character(1), gene_ids = gene_ids)
    } else {
      ""
    }
    local_coverage$present <- if (identifier_type == "gene_symbol_or_mixed") {
      nzchar(local_coverage$matched_symbol)
    } else {
      NA
    }
    local_coverage$audit_error <- if (is.null(audit_error)) "" else audit_error
    coverage_rows[[length(coverage_rows) + 1L]] <- local_coverage
  }
}
if (length(coverage_rows) > 0L) {
  gene_coverage <- do.call(rbind, coverage_rows)
} else {
  gene_coverage <- data.frame(
    module = character(), role = character(), gene = character(),
    filename = character(), identifier_type = character(),
    matched_symbol = character(), present = logical(), audit_error = character(),
    stringsAsFactors = FALSE
  )
}
safe_write_csv(
  gene_coverage,
  file.path(result_dir, "GSE300230_module_gene_coverage_v1.csv")
)

candidate_coverage <- gene_coverage[gene_coverage$gene %in% candidate_panel$gene, , drop = FALSE]
if (nrow(candidate_coverage) > 0L) {
  candidate_coverage <- merge(
    candidate_coverage,
    candidate_panel[, c("gene", "module", "step10B_role")],
    by = c("gene", "module"), all.x = TRUE, sort = FALSE
  )
}
safe_write_csv(
  candidate_coverage,
  file.path(result_dir, "GSE300230_frozen_candidate_gene_coverage_v1.csv")
)

module_coverage_summary <- if (nrow(gene_coverage) > 0L) {
  do.call(rbind, lapply(split(gene_coverage, list(gene_coverage$filename, gene_coverage$module)), function(x) {
    data.frame(
      filename = x$filename[[1L]],
      module = x$module[[1L]],
      role = x$role[[1L]],
      genes_expected = nrow(x),
      genes_present = if (all(is.na(x$present))) NA_integer_ else sum(x$present, na.rm = TRUE),
      coverage_fraction = if (all(is.na(x$present))) NA_real_ else mean(x$present, na.rm = TRUE),
      identifier_type = x$identifier_type[[1L]],
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame(
    filename = character(), module = character(), role = character(),
    genes_expected = integer(), genes_present = integer(),
    coverage_fraction = numeric(), identifier_type = character(),
    stringsAsFactors = FALSE
  )
}
rownames(module_coverage_summary) <- NULL
safe_write_csv(
  module_coverage_summary,
  file.path(result_dir, "GSE300230_module_coverage_summary_v1.csv")
)

resolved_levels <- function(x, unresolved = c("unresolved", "ambiguous")) {
  unique(x[!x %in% unresolved & nzchar(x)])
}
distinct_units <- resolved_levels(design_table$biological_unit, unresolved = "unresolved")
mechanical_levels <- resolved_levels(design_table$mechanical_condition)
tgfb_levels <- resolved_levels(design_table$tgfb_condition)
age_levels <- resolved_levels(design_table$age_group)

known_symbol_files <- unique(gene_coverage$filename[
  gene_coverage$identifier_type == "gene_symbol_or_mixed"
])
candidate_support_by_file <- if (length(known_symbol_files) > 0L) {
  vapply(known_symbol_files, function(filename) {
    x <- candidate_coverage[candidate_coverage$filename == filename, , drop = FALSE]
    length(unique(x$gene[x$present %in% TRUE]))
  }, integer(1))
} else {
  integer()
}
best_candidate_coverage <- if (length(candidate_support_by_file) > 0L) {
  max(candidate_support_by_file)
} else {
  NA_integer_
}

checks <- data.frame(
  check_id = c(
    "official_metadata_loaded", "expected_24_samples", "all_samples_human",
    "at_least_4_biological_units", "mechanical_factor_two_levels",
    "tgfb_factor_two_levels", "age_factor_two_levels",
    "processed_text_matrix_available", "frozen_panel_coverage_at_least_12_of_15"
  ),
  passed = c(
    nrow(sample_metadata) > 0L,
    nrow(sample_metadata) == 24L,
    all(grepl("Homo sapiens", sample_metadata$organism, fixed = TRUE)),
    length(distinct_units) >= 4L,
    length(mechanical_levels) >= 2L,
    length(tgfb_levels) >= 2L,
    length(age_levels) >= 2L,
    nrow(usable_downloads) > 0L,
    !is.na(best_candidate_coverage) && best_candidate_coverage >= 12L
  ),
  observed = c(
    paste0(nrow(sample_metadata), " GSM records"),
    as.character(nrow(sample_metadata)),
    paste(unique(sample_metadata$organism), collapse = " | "),
    paste0(length(distinct_units), ": ", paste(distinct_units, collapse = " | ")),
    paste(mechanical_levels, collapse = " | "),
    paste(tgfb_levels, collapse = " | "),
    paste(age_levels, collapse = " | "),
    paste(usable_downloads$filename, collapse = " | "),
    ifelse(is.na(best_candidate_coverage), "not resolved", paste0(best_candidate_coverage, "/15"))
  ),
  stringsAsFactors = FALSE
)
safe_write_csv(
  checks,
  file.path(result_dir, "GSE300230_step15A_integrity_checks_v1.csv")
)

essential_ids <- c(
  "official_metadata_loaded", "at_least_4_biological_units",
  "mechanical_factor_two_levels", "tgfb_factor_two_levels",
  "processed_text_matrix_available"
)
essential_pass <- all(checks$passed[checks$check_id %in% essential_ids])
full_pass <- all(checks$passed)
gate <- if (full_pass) {
  "PASS"
} else if (essential_pass) {
  "PARTIAL_PROCEED_WITH_RESOLVED_DEBTS"
} else {
  "HOLD"
}
proceed_to_15B <- gate != "HOLD"

decision_lines <- c(
  "# Step 15A GSE300230 provenance and feasibility audit",
  "",
  "## Scope",
  "",
  "- Official GEO metadata and supplementary-file inventory were audited.",
  "- Processed text matrices were downloaded when available.",
  "- No FASTQ, BAM, H5AD, or GEO RAW.tar archive was downloaded.",
  "- No condition-specific expression comparison or biological hypothesis test was performed.",
  "",
  "## Gate",
  "",
  paste0("- Step 15A gate: **", gate, "**."),
  paste0("- Proceed to Step 15B contract freeze: **", proceed_to_15B, "**."),
  paste0("- GSM samples: ", nrow(sample_metadata), "."),
  paste0("- Resolved biological units: ", length(distinct_units), "."),
  paste0("- Mechanical levels: ", paste(mechanical_levels, collapse = " | "), "."),
  paste0("- TGF-beta levels: ", paste(tgfb_levels, collapse = " | "), "."),
  paste0("- Age levels: ", paste(age_levels, collapse = " | "), "."),
  paste0("- Best frozen candidate coverage: ", ifelse(is.na(best_candidate_coverage), "unresolved", paste0(best_candidate_coverage, "/15")), "."),
  "",
  "## Interpretation boundary",
  "",
  "A PASS or PARTIAL gate establishes feasibility only. It does not support mechanotransduction, specificity, or causality. Step 15B must freeze the factor coding, biological-unit model, selected matrix, aliases, endpoints, and multiplicity rules before Step 15C examines expression contrasts.",
  "",
  "## Material Passport",
  "",
  "- Origin: official NCBI GEO GSE300230 metadata and supplementary directory.",
  "- Local frozen inputs: frozen_candidate_panel_v2.csv and mechanotransduction_module_registry_v2.csv.",
  "- Transformation: metadata flattening, conservative factor reconstruction, file classification, and gene-identifier coverage audit.",
  "- Integrity boundary: unresolved metadata fields require explicit Step 15B adjudication; they must not be inferred from expression results.",
  if (nzchar(listing_error)) paste0("- Supplementary listing warning: ", listing_error) else "- Supplementary listing warning: none."
)
decision_path <- file.path(result_dir, "GSE300230_step15A_provenance_feasibility_decision_v1.md")
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 15A GSE300230 provenance and feasibility audit completed.")
message("Gate: ", gate)
message("Proceed to Step 15B: ", proceed_to_15B)
message("Samples / biological units: ", nrow(sample_metadata), " / ", length(distinct_units))
message("Processed matrices: ", nrow(usable_downloads))
message("Best frozen candidate coverage: ", ifelse(is.na(best_candidate_coverage), "unresolved", paste0(best_candidate_coverage, "/15")))
message("Decision: ", decision_path)

