options(stringsAsFactors = FALSE)

# Step 17D2: correct the GSE338388 design labels after manual review.
#
# Step 17D classified all TGF-beta levels as ambiguous because the generic GEO
# treatment_protocol field repeats both treatments for every sample. The
# sample title and the sample-level characteristic "treatment" are the valid
# condition fields for this dataset. This script uses those fields only,
# documents the correction, and performs no expression analysis.

project_dir <- "."
gse_id <- "GSE338388"
audit_dir <- file.path(
  project_dir, "results", "14_second_round_computational_strengthening",
  "17D_S2_GSE338388_design_audit"
)
result_dir <- audit_dir

metadata_path <- file.path(
  audit_dir, paste0(gse_id, "_sample_metadata_snapshot_v1.csv")
)
old_design_path <- file.path(
  audit_dir, paste0(gse_id, "_design_factor_reconstruction_v1.csv")
)
if (!file.exists(metadata_path)) {
  stop("Missing Step 17D sample metadata snapshot: ", metadata_path)
}
if (!file.exists(old_design_path)) {
  stop("Missing Step 17D preliminary design reconstruction: ", old_design_path)
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

metadata <- read.csv(metadata_path, check.names = FALSE)
required_fields <- c("gsm_accession", "title", "characteristics")
missing_fields <- setdiff(required_fields, names(metadata))
if (length(missing_fields) > 0L) {
  stop("GSE338388 metadata is missing field(s): ", paste(missing_fields, collapse = ", "))
}
if (nrow(metadata) != 12L) {
  stop("Expected 12 GSE338388 sample records; found ", nrow(metadata), ".")
}

resolve_tgfb_from_valid_fields <- function(title, characteristics) {
  text <- tolower(paste(title, characteristics, sep = " | "))
  if (grepl("treatment:[[:space:]]*tgf-beta-?1[[:space:]]*,[[:space:]]*k-975", text, perl = TRUE) ||
      grepl("treatment:[^|;]*tgf-beta-?1[^|;]*k-975", text, perl = TRUE)) {
    return("exposed")
  }
  if (grepl("treatment:[[:space:]]*tgf-beta-?1", text, perl = TRUE)) {
    return("exposed")
  }
  if (grepl("vocal fold fibroblast,[[:space:]]*control", text, perl = TRUE) ||
      grepl("treatment:[[:space:]]*none", text, perl = TRUE) ||
      grepl("treatment:[[:space:]]*k-975", text, perl = TRUE)) {
    return("not_exposed")
  }
  "unresolved"
}

resolve_tead_from_valid_fields <- function(title, characteristics) {
  text <- tolower(paste(title, characteristics, sep = " | "))
  if (grepl("treatment:[^|;]*k-975", text, perl = TRUE) ||
      grepl("vocal fold fibroblast,[^|;]*k-975", text, perl = TRUE)) {
    return("tead_inhibited")
  }
  if (grepl("treatment:[[:space:]]*none", text, perl = TRUE) ||
      grepl("treatment:[[:space:]]*tgf-beta-?1[[:space:]]*$", text, perl = TRUE) ||
      grepl("vocal fold fibroblast,[[:space:]]*control", text, perl = TRUE)) {
    return("tead_not_inhibited")
  }
  "unresolved"
}

design <- data.frame(
  gse_accession = gse_id,
  gsm_accession = metadata$gsm_accession,
  title = metadata$title,
  characteristics = metadata$characteristics,
  tgfb_condition = mapply(
    resolve_tgfb_from_valid_fields,
    metadata$title, metadata$characteristics,
    USE.NAMES = FALSE
  ),
  tead_condition = mapply(
    resolve_tead_from_valid_fields,
    metadata$title, metadata$characteristics,
    USE.NAMES = FALSE
  ),
  label_source = "sample title + sample-level characteristics[treatment]",
  generic_treatment_protocol_excluded = TRUE,
  stringsAsFactors = FALSE
)
design$condition_code <- paste(design$tgfb_condition, design$tead_condition, sep = "__")

safe_write_csv(
  design,
  file.path(result_dir, paste0(gse_id, "_design_factor_reconstruction_v2.csv"))
)

factor_count <- function(field) {
  tab <- as.data.frame(table(design[[field]], useNA = "ifany"), stringsAsFactors = FALSE)
  names(tab) <- c("level", "samples")
  tab$factor <- field
  tab[, c("factor", "level", "samples")]
}
factor_counts <- do.call(rbind, lapply(
  c("tgfb_condition", "tead_condition", "condition_code"), factor_count
))
safe_write_csv(
  factor_counts,
  file.path(result_dir, paste0(gse_id, "_design_factor_counts_v2.csv"))
)

combo <- as.data.frame(table(
  factor(design$tgfb_condition, levels = c("not_exposed", "exposed")),
  factor(design$tead_condition, levels = c("tead_not_inhibited", "tead_inhibited"))
), stringsAsFactors = FALSE)
names(combo) <- c("tgfb_condition", "tead_condition", "samples")
complete_2x2 <- all(combo$samples > 0L) &&
  all(design$tgfb_condition %in% c("not_exposed", "exposed")) &&
  all(design$tead_condition %in% c("tead_not_inhibited", "tead_inhibited"))
three_replicates_each <- complete_2x2 && all(combo$samples == 3L)

decision <- if (complete_2x2 && three_replicates_each) {
  "PASS_TO_S2_MATRIX_AUDIT"
} else if (complete_2x2) {
  "PASS_WITH_REPLICATE_IMBALANCE_TO_S2_MATRIX_AUDIT"
} else {
  "HOLD_NEEDS_MANUAL_DESIGN_REVIEW"
}

decision_lines <- c(
  paste0("# Step 17D2 corrected ", gse_id, " design audit"),
  "",
  "## Decision",
  "",
  paste0("- Gate 17D2: **", decision, "**."),
  "- The preliminary automatic parser was corrected after manual review.",
  "- No expression matrix was downloaded and no hypothesis test was performed.",
  "",
  "## Why the preliminary parser failed",
  "",
  "- GEO's generic treatment protocol field states that the experiment used both TGF-beta1 and K-975 for every sample.",
  "- That field describes the protocol, not the sample-specific treatment assignment.",
  "- The valid sample-specific fields are the title and the characteristic `treatment`.",
  "",
  "## Corrected design",
  "",
  paste0("- Samples: ", nrow(design), "."),
  paste0("- TGF-beta levels: ", paste(sort(unique(design$tgfb_condition)), collapse = ", "), "."),
  paste0("- TEAD levels: ", paste(sort(unique(design$tead_condition)), collapse = ", "), "."),
  paste0("- Complete 2×2 design: ", complete_2x2, "."),
  paste0("- Three replicates in each combination: ", three_replicates_each, "."),
  "- Control = no TGF-beta exposure and no TEAD inhibition.",
  "- K-975 = TEAD inhibition without TGF-beta exposure.",
  "- TGF-beta1 = TGF-beta exposure without TEAD inhibition.",
  "- TGF-beta1 + K-975 = both factors present.",
  "",
  "## S2 interpretation boundary",
  "",
  "- GSE338388 has no mechanical loading or stiffness factor.",
  "- A passing design gate supports TGF-beta/SMAD versus TEAD-related regulatory-axis cross-validation only.",
  "- It cannot prove mechanical causality, YAP/TAZ-independent mechanical driving, or fascia-direct replication.",
  "- Frozen candidate and module coverage must be checked in the processed matrix before any expression analysis.",
  "- Overall project evidence grade remains CAUTION.",
  "",
  "## Material Passport",
  "",
  "- Input: official GSE338388 sample metadata snapshot from Step 17D.",
  "- Transformation: manual-review-supported correction using sample-specific title and treatment characteristic.",
  "- New data downloaded: none.",
  "- Next step: Step 17D3 S2 processed-matrix audit and frozen panel/module coverage check."
)
decision_path <- file.path(
  result_dir, paste0(gse_id, "_Step17D2_corrected_design_decision_v1.md")
)
writeLines(decision_lines, decision_path, useBytes = TRUE)

message("Step 17D2 corrected ", gse_id, " design audit completed.")
message("Samples: ", nrow(design))
message("Complete 2x2 design: ", complete_2x2)
message("Three replicates each: ", three_replicates_each)
message("Gate 17D2: ", decision)
message("Decision: ", decision_path)
