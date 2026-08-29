# Step 10A: audit Steele 2025 atlas provenance, source overlap and label versions
#
# This is a metadata/provenance audit only.
# It does not download the full atlas, does not test mechanosensitivity,
# and does not change the failed Step 08C2 PIEZO2 gate.

project_dir <- "."
out_dir <- file.path(
  project_dir,
  "results",
  "08_cross_tissue_validation",
  "10A_atlas_provenance_audit"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

official_source_url <- paste0(
  "https://raw.githubusercontent.com/haniffalab/skin_fibroblast_atlas/",
  "refs/heads/main/data/scrna_seq_data.csv"
)
official_source_file <- file.path(
  out_dir,
  "skin_fibroblast_atlas_official_scrna_seq_data_v1.csv"
)

if (!file.exists(official_source_file) || file.info(official_source_file)$size == 0) {
  message("Downloading the small official atlas source inventory.")
  download.file(
    url = official_source_url,
    destfile = official_source_file,
    mode = "wb",
    quiet = FALSE
  )
}

if (!file.exists(official_source_file) || file.info(official_source_file)$size == 0) {
  stop("Official atlas source inventory could not be downloaded.")
}

source_inventory <- read.csv(
  official_source_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

if (ncol(source_inventory) < 5) {
  stop("Unexpected official source inventory structure.")
}

names(source_inventory)[1:5] <- c(
  "source_link_or_accession",
  "first_author",
  "pmid",
  "publication_year",
  "diseases_included"
)

extract_accession <- function(x) {
  x <- toupper(trimws(as.character(x)))
  matches <- regmatches(
    x,
    regexec("(GSE[0-9]+|PRJNA[0-9]+)", x, perl = TRUE)
  )
  vapply(
    matches,
    FUN = function(hit) {
      if (length(hit) >= 2) hit[[2]] else NA_character_
    },
    FUN.VALUE = character(1)
  )
}

source_inventory$source_accession <- extract_accession(
  source_inventory$source_link_or_accession
)

expected_overlap_accessions <- c("GSE173252", "PRJNA607098")
if (!all(expected_overlap_accessions %in% source_inventory$source_accession)) {
  stop(
    "Accession extraction failed: expected GSE173252 and PRJNA607098 were not both resolved."
  )
}

source_inventory$overlaps_GSE173252 <-
  source_inventory$source_accession == "GSE173252"
source_inventory$overlaps_PRJNA607098 <-
  source_inventory$source_accession == "PRJNA607098"
source_inventory$overlaps_existing_project <-
  source_inventory$overlaps_GSE173252 |
  source_inventory$overlaps_PRJNA607098
source_inventory$independent_candidate_after_exclusion <-
  !is.na(source_inventory$source_accession) &
  !source_inventory$overlaps_existing_project
source_inventory$source_role <- ifelse(
  source_inventory$overlaps_GSE173252,
  "discovery_overlap_GSE173252",
  ifelse(
    source_inventory$overlaps_PRJNA607098,
    "validation_overlap_PRJNA607098",
    ifelse(
      is.na(source_inventory$source_accession),
      "non_accession_source_requires_manual_resolution",
      "independent_candidate_pending_sample_audit"
    )
  )
)

local_metadata_file <- file.path(
  project_dir,
  "results",
  "06_external_validation",
  "skin_fibroblast_atlas_2025",
  "08_source_specific_PRJNA607098",
  "skin_fibroblast_atlas_metadata_counts_v1.csv"
)

if (file.exists(local_metadata_file)) {
  local_metadata <- read.csv(
    local_metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  local_sources <- unique(toupper(trimws(as.character(local_metadata$GSE))))
  source_inventory$present_in_existing_step08_source_audit <-
    source_inventory$source_accession %in% local_sources
} else {
  source_inventory$present_in_existing_step08_source_audit <- NA
}

source_inventory_out <- file.path(
  out_dir,
  "skin_fibroblast_atlas_source_accession_map_v1.csv"
)
write.csv(
  source_inventory,
  source_inventory_out,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Paper and local atlas label crosswalk.
# The paper uses F8 for fascia-like myofibroblasts, while the Step 08B3
# source-specific analysis used the atlas celltype label
# "F7: Fascia-like myofibroblast". This is deliberately recorded as a
# version-mapping issue until the official object definition is audited.
label_crosswalk <- data.frame(
  label_system = c(
    "Nature article",
    "Nature article",
    "Nature article",
    "Step 08B3 atlas celltype",
    "Step 08B2 atlas celltype_skinspecific_nomenclature",
    "Step 08B2 available cell-state inventory"
  ),
  label = c(
    "F6: inflammatory myofibroblast",
    "F7: myofibroblast",
    "F8: fascia-like myofibroblast",
    "F7: Fascia-like myofibroblast",
    "Myofibroblast",
    "Fascia"
  ),
  interpretation = c(
    "paper disease-specific inflammatory myofibroblast state",
    "paper disease-specific myofibroblast state",
    "paper fascia-like state; official paper nomenclature",
    "label used in the current source-specific F7 analysis",
    "auxiliary skin-specific label used in Step 08B2",
    "available state label in the corrected local inventory"
  ),
  observed_in_local_outputs = c(
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    TRUE
  ),
  mapping_status = c(
    "requires_object-level confirmation",
    "requires_object-level confirmation",
    "official_paper_label",
    "confirmed_in_step08B3_output",
    "confirmed_in_step08B2_output",
    "confirmed_in_step08B2_output"
  ),
  stringsAsFactors = FALSE
)

label_crosswalk_out <- file.path(
  out_dir,
  "skin_fibroblast_atlas_F6_F7_F8_label_crosswalk_v1.csv"
)
write.csv(
  label_crosswalk,
  label_crosswalk_out,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# The official source table contains accession-level provenance but no
# sample/donor column. Existing Step 08C1 sample reconstruction is retained
# as a separate project result and is not treated as evidence that all atlas
# sources have donor labels.
sample_field_audit <- data.frame(
  audit_item = c(
    "official_source_table_accession_field",
    "official_source_table_sample_or_donor_field",
    "existing_PRJNA607098_sample_audit",
    "current_step08B3_target_label",
    "paper_fascia_like_label"
  ),
  field_or_value = c(
    "source_link_or_accession",
    "not present in scrna_seq_data.csv",
    "obs/_index reconstructed sample IDs in Step 08C1",
    "F7: Fascia-like myofibroblast in atlas celltype",
    "F8: fascia-like myofibroblast in Nature article"
  ),
  evidence_scope = c(
    "official atlas source inventory",
    "official atlas source inventory",
    "PRJNA607098-only local audit",
    "PRJNA607098-only local audit",
    "official paper nomenclature"
  ),
  conclusion = c(
    "source-level provenance available",
    "sample-level provenance not resolved by this inventory",
    "sample IDs available for the PRJNA607098 source audit",
    "target used by Step 08B3; do not assume identical to paper F8",
    "requires explicit crosswalk against the official object"
  ),
  stringsAsFactors = FALSE
)

sample_field_audit_out <- file.path(
  out_dir,
  "skin_fibroblast_atlas_source_sample_field_audit_v1.csv"
)
write.csv(
  sample_field_audit,
  sample_field_audit_out,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

accession_values <- unique(na.omit(source_inventory$source_accession))
n_accessions <- length(accession_values)
n_independent <- sum(
  source_inventory$independent_candidate_after_exclusion,
  na.rm = TRUE
)
n_overlap_discovery <- sum(
  source_inventory$overlaps_GSE173252,
  na.rm = TRUE
)
n_overlap_validation <- sum(
  source_inventory$overlaps_PRJNA607098,
  na.rm = TRUE
)

# Gate A is intentionally PARTIAL: independent accessions are listed, but
# sample-level provenance and F7/F8 label equivalence are not yet resolved.
gate_a <- if (n_independent > 0) "PARTIAL" else "FAIL"

decision_lines <- c(
  "# Step 10A atlas provenance and label audit decision",
  "",
  "> This is a metadata/provenance gate. It is not a module validation result and does not change the Step 08C2 PIEZO2 decision.",
  "",
  "## Official source inventory",
  "",
  paste0("- Official source inventory URL: ", official_source_url),
  paste0("- Accession-resolved source rows: ", n_accessions),
  paste0("- Rows overlapping GSE173252: ", n_overlap_discovery),
  paste0("- Rows overlapping PRJNA607098: ", n_overlap_validation),
  paste0("- Independent candidate accession rows after exclusion: ", n_independent),
  "",
  "The official inventory includes both GSE173252 and PRJNA607098. Therefore, the integrated atlas cannot be treated as a fully independent validation cohort unless these source datasets are explicitly excluded from the validation layer.",
  "",
  "## Label finding",
  "",
  "The Nature article uses F8 for fascia-like myofibroblasts, whereas the current Step 08B3 source-specific output uses F7: Fascia-like myofibroblast in the atlas celltype field. The label relationship is recorded as unresolved until the official object definition and source-specific mapping are audited.",
  "",
  "## Sample-level provenance finding",
  "",
  "The official scrna_seq_data.csv inventory is source-level and does not provide a sample/donor column. The existing obs/_index sample reconstruction applies to PRJNA607098 only and cannot be generalized to every atlas source.",
  "",
  "## Gate A",
  "",
  paste0("- Gate A status: **", gate_a, "**"),
  "- PASS requires an independent source after exclusion plus sample/donor provenance and a resolved target-label crosswalk.",
  "- PARTIAL permits source-level exploratory mapping only; it does not upgrade the current evidence grade.",
  "",
  "## Next action",
  "",
  "1. Restrict any cross-tissue validation to source accessions that are not GSE173252 or PRJNA607098.",
  "2. Audit source-specific sample/donor fields before computing any donor-level result.",
  "3. Resolve the F7/F8 nomenclature before comparing fascia-like states.",
  "4. Do not download the full atlas or run a new module test in this step.",
  "",
  "## Output files",
  "",
  paste0("- Source map: ", source_inventory_out),
  paste0("- Label crosswalk: ", label_crosswalk_out),
  paste0("- Sample-field audit: ", sample_field_audit_out)
)

decision_out <- file.path(
  out_dir,
  "skin_fibroblast_atlas_independence_decision_v1.md"
)
writeLines(decision_lines, decision_out, useBytes = TRUE)

message("Step 10A metadata/provenance audit completed.")
message("Gate A: ", gate_a)
message("Source map: ", source_inventory_out)
message("Label crosswalk: ", label_crosswalk_out)
message("Sample-field audit: ", sample_field_audit_out)
message("Decision: ", decision_out)
