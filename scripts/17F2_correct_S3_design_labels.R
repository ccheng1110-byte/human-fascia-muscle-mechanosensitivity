options(stringsAsFactors = FALSE)

# Step 17F2: manual-review-supported correction of S3 design labels.
# GSE123100 contains an HTM culture stiffness series plus separate clinical
# tissue samples. GSE276045 contains WI-38 WT/hTERT, stiffness, timepoint and
# replicate labels; WT/hTERT is retained as a cell-model factor and is not
# relabeled as a proliferation state.
# This step reads existing metadata/audit files only and downloads nothing.

project_dir <- "."
input_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17F_S3_source_feasibility_audit"
)
result_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17F2_S3_corrected_design_audit"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  file.path(input_dir, "GSE123100", "GSE123100_sample_metadata_audit_v1.csv"),
  file.path(input_dir, "GSE123100", "GSE123100_supplementary_inventory_v1.csv"),
  file.path(input_dir, "GSE276045", "GSE276045_sample_metadata_audit_v1.csv"),
  file.path(input_dir, "GSE276045", "GSE276045_supplementary_inventory_v1.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing Step 17F input(s): ", paste(missing_inputs, collapse = "; "))
}

safe_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

parse_gse123100 <- function(metadata) {
  title <- trimws(as.character(metadata$title))
  cultured <- grepl("^[AC]_[0-9]+(?:\\.[0-9]+)?[[:space:]]*kPa$", title, perl = TRUE)
  group <- ifelse(cultured, sub("^([AC])_.*$", "\\1", title, perl = TRUE), "clinical")
  stiffness <- rep(NA_real_, nrow(metadata))
  stiffness[cultured] <- suppressWarnings(as.numeric(sub(
    "^[AC]_([0-9]+(?:\\.[0-9]+)?).*", "\\1", title[cultured], perl = TRUE
  )))
  duration <- rep(NA_real_, nrow(metadata))
  duration[grepl("three[[:space:]]+day", metadata$characteristics, ignore.case = TRUE, perl = TRUE)] <- 3
  duration[grepl("five[[:space:]]+day", metadata$characteristics, ignore.case = TRUE, perl = TRUE)] <- 5
  eligible <- cultured & is.finite(stiffness) & is.finite(duration)
  data.frame(
    gse_accession = metadata$gse_accession,
    gsm_accession = metadata$gsm_accession,
    title = metadata$title,
    source_name = metadata$source_name,
    analysis_role = ifelse(eligible, "S3_primary_HTM_stiffness_series", "excluded_clinical_tissue"),
    eligible_for_primary_dose_response = eligible,
    culture_group = group,
    culture_duration_days = duration,
    stiffness_value_kPa = stiffness,
    stiffness_label = ifelse(is.finite(stiffness), paste0(stiffness, " kPa"), "unresolved"),
    biological_unit = ifelse(eligible, paste0("culture_group_", group), "clinical_subject_or_tissue"),
    replicate_structure = ifelse(eligible, "one sample per group x stiffness in most cells", "not applicable"),
    stringsAsFactors = FALSE
  )
}

parse_gse276045 <- function(metadata) {
  title <- trimws(as.character(metadata$title))
  parsed <- grepl(
    "^(WT|hTERT),[[:space:]]*[0-9]+(?:\\.[0-9]+)?[[:space:]]*(?:kPa|GPa),[[:space:]]*timepoint[[:space:]]*[0-9]+,.*replicate[[:space:]]*[0-9]+$",
    title, ignore.case = TRUE, perl = TRUE
  )
  cell_model <- ifelse(
    grepl("^WT,", title, ignore.case = TRUE, perl = TRUE), "WT",
    ifelse(grepl("^hTERT,", title, ignore.case = TRUE, perl = TRUE), "hTERT", "unresolved")
  )
  stiffness_value <- suppressWarnings(as.numeric(sub(
    "^(?:WT|hTERT),[[:space:]]*([0-9]+(?:\\.[0-9]+)?).*", "\\1", title,
    ignore.case = TRUE, perl = TRUE
  )))
  stiffness_unit <- sub(
    "^(?:WT|hTERT),[[:space:]]*[0-9]+(?:\\.[0-9]+)?[[:space:]]*(kPa|GPa).*$",
    "\\1", title, ignore.case = TRUE, perl = TRUE
  )
  stiffness_unit[!parsed] <- "unresolved"
  stiffness_kPa <- ifelse(
    parsed & tolower(stiffness_unit) == "gpa", stiffness_value * 1e6,
    ifelse(parsed, stiffness_value, NA_real_)
  )
  timepoint <- suppressWarnings(as.integer(sub(
    ".*timepoint[[:space:]]*([0-9]+).*", "\\1", title,
    ignore.case = TRUE, perl = TRUE
  )))
  replicate <- suppressWarnings(as.integer(sub(
    ".*replicate[[:space:]]*([0-9]+)$", "\\1", title,
    ignore.case = TRUE, perl = TRUE
  )))
  timepoint[!parsed] <- NA_integer_
  replicate[!parsed] <- NA_integer_
  data.frame(
    gse_accession = metadata$gse_accession,
    gsm_accession = metadata$gsm_accession,
    title = metadata$title,
    source_name = metadata$source_name,
    analysis_role = ifelse(parsed, "S3_WI38_stiffness_x_cell_model_cross_check", "unresolved"),
    eligible_for_stiffness_cross_check = parsed,
    cell_model_condition = cell_model,
    stiffness_value_original = stiffness_value,
    stiffness_unit = stiffness_unit,
    stiffness_value_kPa = stiffness_kPa,
    log10_stiffness_kPa = ifelse(is.finite(stiffness_kPa), log10(stiffness_kPa), NA_real_),
    timepoint = timepoint,
    technical_or_biological_replicate_label = replicate,
    proliferation_state = "not_measured_or_not_resolved",
    biological_unit = paste(cell_model, "timepoint", timepoint, "replicate", replicate, sep = "_"),
    stringsAsFactors = FALSE
  )
}

gse123100_metadata <- utils::read.csv(
  file.path(input_dir, "GSE123100", "GSE123100_sample_metadata_audit_v1.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
gse276045_metadata <- utils::read.csv(
  file.path(input_dir, "GSE276045", "GSE276045_sample_metadata_audit_v1.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
gse123100_inventory <- utils::read.csv(
  file.path(input_dir, "GSE123100", "GSE123100_supplementary_inventory_v1.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
gse276045_inventory <- utils::read.csv(
  file.path(input_dir, "GSE276045", "GSE276045_supplementary_inventory_v1.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)

gse123100_design <- parse_gse123100(gse123100_metadata)
gse276045_design <- parse_gse276045(gse276045_metadata)
safe_write_csv(
  gse123100_design,
  file.path(result_dir, "GSE123100_corrected_S3_design_v2.csv")
)
safe_write_csv(
  gse276045_design,
  file.path(result_dir, "GSE276045_corrected_S3_design_v2.csv")
)

gse123100_primary <- gse123100_design[gse123100_design$eligible_for_primary_dose_response, , drop = FALSE]
gse123100_dose_counts <- as.data.frame(table(
  gse123100_primary$culture_group,
  gse123100_primary$stiffness_label,
  useNA = "ifany"
), stringsAsFactors = FALSE)
names(gse123100_dose_counts) <- c("culture_group", "stiffness_label", "samples")
gse123100_dose_counts <- gse123100_dose_counts[gse123100_dose_counts$samples > 0L, , drop = FALSE]
safe_write_csv(
  gse123100_dose_counts,
  file.path(result_dir, "GSE123100_corrected_dose_counts_v2.csv")
)

gse276045_eligible <- gse276045_design[gse276045_design$eligible_for_stiffness_cross_check, , drop = FALSE]
gse276045_counts <- as.data.frame(table(
  gse276045_eligible$cell_model_condition,
  gse276045_eligible$stiffness_unit,
  gse276045_eligible$timepoint,
  useNA = "ifany"
), stringsAsFactors = FALSE)
names(gse276045_counts) <- c("cell_model_condition", "stiffness_unit", "timepoint", "samples")
gse276045_counts <- gse276045_counts[gse276045_counts$samples > 0L, , drop = FALSE]
safe_write_csv(
  gse276045_counts,
  file.path(result_dir, "GSE276045_corrected_design_counts_v2.csv")
)

gse123100_cell_matrix <- gse123100_inventory[
  grepl("cell_expressed_gene_RPKM", gse123100_inventory$filename, fixed = TRUE),
  , drop = FALSE
]
gse276045_matrix <- gse276045_inventory[
  grepl("bulk_RNA_seq_1_counts", gse276045_inventory$filename, fixed = TRUE),
  , drop = FALSE
]
selected_matrix_candidates <- rbind(
  data.frame(
    gse_accession = "GSE123100",
    analysis_role = "S3_primary_HTM_stiffness_series",
    matrix_priority = "primary_RPKM_expression",
    matrix_filename = if (nrow(gse123100_cell_matrix) > 0L) gse123100_cell_matrix$filename[[1L]] else "not_found",
    direct_url = if (nrow(gse123100_cell_matrix) > 0L) gse123100_cell_matrix$direct_url[[1L]] else "",
    download_authorized = nrow(gse123100_cell_matrix) > 0L,
    stringsAsFactors = FALSE
  ),
  data.frame(
    gse_accession = "GSE276045",
    analysis_role = "S3_WI38_stiffness_x_cell_model_cross_check",
    matrix_priority = "official_bulk_counts",
    matrix_filename = if (nrow(gse276045_matrix) > 0L) gse276045_matrix$filename[[1L]] else "not_found",
    direct_url = if (nrow(gse276045_matrix) > 0L) gse276045_matrix$direct_url[[1L]] else "",
    download_authorized = nrow(gse276045_matrix) > 0L,
    stringsAsFactors = FALSE
  )
)
safe_write_csv(
  selected_matrix_candidates,
  file.path(result_dir, "Step17F2_selected_processed_matrix_candidates_v1.csv")
)

gse123100_gate <- if (
  nrow(gse123100_primary) == 11L &&
    length(unique(gse123100_primary$stiffness_value_kPa)) >= 5L &&
    nrow(gse123100_cell_matrix) > 0L
) {
  "PARTIAL_PROCEED_TO_PRIMARY_MATRIX_AUDIT_WITH_DESCRIPTIVE_DOSE_LIMIT"
} else {
  "HOLD_GSE123100_DESIGN_OR_MATRIX_REVIEW"
}
gse276045_gate <- if (
  nrow(gse276045_eligible) > 0L &&
    length(unique(gse276045_eligible$stiffness_value_kPa)) >= 5L &&
    all(gse276045_eligible$cell_model_condition %in% c("WT", "hTERT")) &&
    nrow(gse276045_matrix) > 0L
) {
  "PARTIAL_PROCEED_TO_WI38_CELL_MODEL_MATRIX_AUDIT"
} else {
  "HOLD_GSE276045_DESIGN_OR_MATRIX_REVIEW"
}

summary <- data.frame(
  gse_accession = c("GSE123100", "GSE276045"),
  corrected_role = c(
    "S3_primary_HTM_stiffness_dose_response",
    "S3_WI38_stiffness_x_cell_model_cross_check"
  ),
  total_samples = c(nrow(gse123100_design), nrow(gse276045_design)),
  eligible_samples = c(nrow(gse123100_primary), nrow(gse276045_eligible)),
  excluded_samples = c(sum(!gse123100_design$eligible_for_primary_dose_response), 0L),
  stiffness_levels = c(
    length(unique(gse123100_primary$stiffness_value_kPa)),
    length(unique(gse276045_eligible$stiffness_value_kPa))
  ),
  replicate_balance = c(
    "mostly one sample per culture_group x stiffness cell",
    "replicate labels present; timepoint/stiffness cells unbalanced"
  ),
  proliferation_status = c(
    "not relevant to primary HTM dose series",
    "not measured_or_not_resolved; do not use WT/hTERT as proliferation"
  ),
  decision = c(gse123100_gate, gse276045_gate),
  stringsAsFactors = FALSE
)
safe_write_csv(
  summary,
  file.path(result_dir, "Step17F2_S3_corrected_design_summary_v1.csv")
)

decision_path <- file.path(result_dir, "Step17F2_S3_corrected_design_decision_v1.md")
decision_lines <- c(
  "# Step 17F2 S3 corrected design audit",
  "",
  "## GSE123100",
  "",
  paste0("- Corrected eligible primary samples: ", nrow(gse123100_primary), "/", nrow(gse123100_design), "."),
  paste0("- Excluded clinical tissue samples: ", sum(!gse123100_design$eligible_for_primary_dose_response), "."),
  paste0("- Corrected stiffness levels: ", paste(sort(unique(gse123100_primary$stiffness_value_kPa)), collapse = ", "), " kPa."),
  paste0("- Corrected gate: **", gse123100_gate, "**."),
  "- The HTM series can support a cautious dose-response-form analysis, but most stiffness-by-culture-group cells contain one sample; formal inferential claims must remain limited.",
  "- The eight clinical trabecular-meshwork tissue samples are excluded from the primary stiffness series and must not be mixed with cultured HTM samples.",
  "",
  "## GSE276045",
  "",
  paste0("- Corrected eligible samples: ", nrow(gse276045_eligible), "/", nrow(gse276045_design), "."),
  paste0("- Cell-model conditions: ", paste(sort(unique(gse276045_eligible$cell_model_condition)), collapse = ", "), "."),
  paste0("- Corrected stiffness levels retained: ", paste(sort(unique(gse276045_eligible$stiffness_value_kPa)), collapse = ", "), " kPa-equivalent."),
  paste0("- Corrected gate: **", gse276045_gate, "**."),
  "- WT/hTERT is retained as a cell-model/genotype factor, not relabeled as proliferation status.",
  "- The planned stiffness-by-proliferation confounding test is therefore NOT_ESTIMABLE from current metadata; the next analysis may only test stiffness direction by cell-model condition and timepoint, with this limitation explicit.",
  "",
  "## S3 boundary",
  "",
  "- Neither source is fascia-direct validation.",
  "- No GWAS branch is activated by this audit.",
  "- A passing targeted matrix audit will authorize descriptive or model-based cross-tissue checks only; it will not change the overall project evidence grade from CAUTION.",
  "",
  "## Material Passport",
  "",
  "- Input: Step 17F official GEO metadata and supplementary inventories.",
  "- Transformation: title/characteristic-supported correction of cultured versus clinical samples, stiffness values, cell-model conditions, timepoints and replicate labels.",
  "- New data downloaded: none.",
  "- Next step: targeted audit of the selected official processed matrices using these corrected sample maps."
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17F2 corrected S3 design audit completed.")
message("GSE123100 gate: ", gse123100_gate, "; eligible primary samples: ", nrow(gse123100_primary))
message("GSE276045 gate: ", gse276045_gate, "; eligible samples: ", nrow(gse276045_eligible))
message("Proliferation confounding for GSE276045: NOT_ESTIMABLE")
message("Decision: ", decision_path)
