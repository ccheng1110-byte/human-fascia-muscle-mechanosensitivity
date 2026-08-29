# Step 14C: Generate manuscript figures and tables from frozen result tables
#
# Contract
# Core conclusion:
#   Human fibroblast-associated states show a reproducible multicomponent
#   ECM/integrin-cytoskeletal program, but current data do not establish
#   mechanosensitivity-specificity or a universal PIEZO2 mechanism.
# Figure archetype: quantitative grid with source-separated validation panels.
# Backend: R only.
# Target output: Nature-style editable SVG/PDF plus 600-dpi TIFF.
# Final width: 183 mm; heights are kept below 130 mm.
# Evidence hierarchy: primary paired-source module/candidate evidence first;
#   independent external-source triangulation second; competition robustness as
#   a boundary/control panel.
# Reviewer risk: source-specific definitions, descriptive external rule, and
#   biological replicate units must remain visible in legends and source data.

options(stringsAsFactors = FALSE)

project_dir <- "."
result_dir <- file.path(project_dir, "results", "11_manuscript_preparation",
                        "14C_figures_and_tables")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "svglite", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R package(s): ", paste(missing_packages, collapse = ", "),
      ". Install once with: install.packages(c(",
      paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
})

read_csv_checked <- function(path) {
  if (!file.exists(path)) stop("Required input file not found: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

as_logical_strict <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

module_label <- function(x) {
  labels <- c(
    actomyosin_rho = "Actomyosin/Rho",
    integrin_focal_adhesion = "Integrin/focal adhesion",
    hippo_yap_taz = "Hippo/YAP/TAZ",
    mechanosensor_channels = "Mechanosensor channels",
    ecm_remodeling = "ECM remodeling",
    tgf_fibrosis = "TGF/fibrosis",
    inflammation_ap1_nfkb = "Inflammation/AP-1/NF-kB",
    hypoxia = "Hypoxia",
    cell_cycle = "Cell cycle"
  )
  out <- unname(labels[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

palette_contract <- c(
  primary = "#3182BD",
  external = "#33A69A",
  supported = "#3182BD",
  not_supported = "#D8D8D8",
  boundary = "#E28E2C",
  negative = "#D24B40",
  neutral_dark = "#272727"
)

theme_nature_contract <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.8),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.3),
      panel.grid = element_blank()
    )
}

theme_set(theme_nature_contract())

save_pub_r <- function(plot, stem, width_mm = 183, height_mm = 110, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h,
                       family = "Arial")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h,
                 units = "in", res = dpi, compression = "lzw")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h,
                units = "in", res = 300)
  print(plot)
  grDevices::dev.off()
}

primary_candidate_path <- file.path(
  project_dir, "results", "08_cross_tissue_validation", "10C_frozen_PRJNA607098",
  "PRJNA607098_F7_frozen_candidate_statistics_v2.csv"
)
primary_module_path <- file.path(
  project_dir, "results", "06_external_validation", "skin_fibroblast_atlas_2025",
  "08C2_sample_level_F7_PRJNA607098",
  "PRJNA607098_F7_sample_paired_module_summary_v1.csv"
)
external_candidate_path <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11F_GSE130973_external_validation_summary",
  "GSE130973_candidate_direction_summary_v1.csv"
)
external_module_path <- file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11F_GSE130973_external_validation_summary",
  "GSE130973_module_boundary_summary_v1.csv"
)
competition_path <- file.path(
  project_dir, "results", "08_cross_tissue_validation", "10D_evidence_synthesis",
  "Step10D_failed_competition_pairs_v1.csv"
)
evidence_summary_path <- file.path(
  project_dir, "results", "08_cross_tissue_validation", "10D_evidence_synthesis",
  "Step10D_evidence_summary_v1.csv"
)
claims_path <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12C_computational_closeout_and_functional_validation",
  "12C_frozen_claims_and_boundaries_v1.csv"
)
endpoint_path <- file.path(
  project_dir, "results", "10_final_evidence_synthesis",
  "12C_computational_closeout_and_functional_validation",
  "12C_functional_validation_endpoint_matrix_v1.csv"
)

primary_candidate <- read_csv_checked(primary_candidate_path)
primary_module <- read_csv_checked(primary_module_path)
external_candidate <- read_csv_checked(external_candidate_path)
external_module <- read_csv_checked(external_module_path)
competition <- read_csv_checked(competition_path)
evidence_summary <- read_csv_checked(evidence_summary_path)
claims <- read_csv_checked(claims_path)
endpoints <- read_csv_checked(endpoint_path)

primary_candidate$strict_sample_level_support <- as_logical_strict(
  primary_candidate$strict_sample_level_support
)
external_candidate$descriptive_directional_support_4_of_5 <- as_logical_strict(
  external_candidate$descriptive_directional_support_4_of_5
)
external_module$positive_5_of_5 <- as_logical_strict(external_module$positive_5_of_5)
competition$residual_median <- as.numeric(competition$residual_median)
competition$residual_positive_samples <- as.numeric(competition$residual_positive_samples)

# ---------------------------- Figure 2 -------------------------------------
# Source-separated module support. The denominator and support rule are kept
# source-specific and are printed in the panel labels/source data.
primary_source_label <- "PRJNA607098\n12 paired sample-state units"
external_source_label <- "GSE130973\n5 subjects"
primary_module_plot <- data.frame(
  source = primary_source_label,
  module = module_label(primary_module$module),
  module_key = primary_module$module,
  support_fraction = as.numeric(primary_module$strictly_supported_genes) /
    as.numeric(primary_module$candidate_genes),
  n_text = paste0(primary_module$strictly_supported_genes, "/",
                  primary_module$candidate_genes),
  rule = "strict sample-level support",
  stringsAsFactors = FALSE
)
external_module_plot <- data.frame(
  source = external_source_label,
  module = module_label(external_module$module),
  module_key = external_module$module,
  support_fraction = as.numeric(external_module$positive_samples) /
    as.numeric(external_module$samples),
  n_text = paste0(external_module$positive_samples, "/",
                  external_module$samples),
  rule = "descriptive positive-subject fraction",
  stringsAsFactors = FALSE
)
module_plot <- rbind(primary_module_plot, external_module_plot)
module_plot$source <- factor(
  module_plot$source,
  levels = c(primary_source_label, external_source_label)
)
module_plot$module <- factor(
  module_plot$module,
  levels = c("Actomyosin/Rho", "Integrin/focal adhesion", "Hippo/YAP/TAZ",
             "Mechanosensor channels", "ECM remodeling", "TGF/fibrosis",
             "Inflammation/AP-1/NF-kB", "Hypoxia", "Cell cycle")
)
module_plot <- module_plot[!is.na(module_plot$module), , drop = FALSE]

p_module <- ggplot(module_plot,
                    aes(x = module, y = support_fraction, fill = source)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68,
           colour = "white", linewidth = 0.2) +
  geom_text(aes(label = n_text),
            position = position_dodge(width = 0.78), vjust = -0.35,
            size = 2.2, family = "Arial") +
  scale_fill_manual(
    values = c(unname(palette_contract["primary"]),
               unname(palette_contract["external"])),
    limits = c(primary_source_label, external_source_label),
    breaks = c(primary_source_label, external_source_label)
  ) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Module-level directional support is source-dependent",
    subtitle = "Bars show source-specific support rules; denominators are gene or subject counts",
    x = NULL, y = "Support fraction", fill = "Source"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top", legend.key.height = grid::unit(3, "mm"))

save_pub_r(
  p_module,
  file.path(result_dir, "14C_Figure2_module_support_by_source"),
  width_mm = 183, height_mm = 112
)

# ---------------------------- Figure 3 -------------------------------------
# Candidate-level heatmap. Text gives the source-specific median difference;
# fill encodes only support status so unlike expression scales are not pooled.
candidate_order <- unique(primary_candidate$gene)
candidate_order <- candidate_order[!is.na(candidate_order)]

primary_cand_plot <- primary_candidate[, c(
  "gene", "module", "median_difference", "strict_sample_level_support"
)]
names(primary_cand_plot) <- c("gene", "module", "effect", "supported")
primary_candidate_source_label <- "PRJNA607098\nstrict 12-unit rule"
external_candidate_source_label <- "GSE130973\ndescriptive 4/5 rule"
primary_cand_plot$source <- primary_candidate_source_label
primary_cand_plot$rule <- "strict sample-level support"

external_cand_plot <- external_candidate[, c(
  "gene", "module", "median_difference",
  "descriptive_directional_support_4_of_5"
)]
names(external_cand_plot) <- c("gene", "module", "effect", "supported")
external_cand_plot$source <- external_candidate_source_label
external_cand_plot$rule <- "descriptive positive-subject direction"

candidate_plot <- rbind(primary_cand_plot, external_cand_plot)
candidate_plot$effect <- as.numeric(candidate_plot$effect)
candidate_plot$supported <- as_logical_strict(candidate_plot$supported)
candidate_plot$gene <- factor(candidate_plot$gene, levels = rev(candidate_order))
candidate_plot$source <- factor(candidate_plot$source,
                                levels = c(primary_candidate_source_label,
                                           external_candidate_source_label))
candidate_plot$module_label <- module_label(candidate_plot$module)
candidate_plot$effect_label <- ifelse(
  is.na(candidate_plot$effect), "NA", sprintf("%.2f", candidate_plot$effect)
)
candidate_plot$support_label <- ifelse(candidate_plot$supported,
                                       "supported", "not supported")
candidate_plot$support_label <- factor(candidate_plot$support_label,
                                       levels = c("supported", "not supported"))

p_candidate <- ggplot(candidate_plot,
                       aes(x = source, y = gene, fill = support_label)) +
  geom_tile(colour = "white", linewidth = 0.45, width = 0.92, height = 0.92) +
  geom_text(aes(label = effect_label), size = 2.05, family = "Arial") +
  facet_grid(module_label ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_fill_manual(
    values = c(unname(palette_contract["supported"]),
               unname(palette_contract["not_supported"])),
    limits = c("supported", "not supported"),
    breaks = c("supported", "not supported")
  ) +
  labs(
    title = "Candidate-level support is not conserved across sources",
    subtitle = "Tile colour = source-specific support; text = source-specific median directional difference",
    x = NULL, y = NULL, fill = "Support status"
  ) +
  theme(
    axis.text.x = element_text(size = 6.1),
    axis.text.y = element_text(size = 6.2, face = "italic"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6.2),
    legend.position = "top",
    panel.spacing.y = grid::unit(2.5, "mm")
  )

save_pub_r(
  p_candidate,
  file.path(result_dir, "14C_Figure3_candidate_support_cross_source"),
  width_mm = 183, height_mm = 125
)

# ---------------------------- Figure 4 -------------------------------------
# Competition robustness boundary. The residual is the adjusted core-module
# signal after competitor stratification; values are retained as generated.
competition$competitor_label <- module_label(competition$competitor_module)
competition$competitor_label <- factor(
  competition$competitor_label,
  levels = rev(c("TGF/fibrosis", "Inflammation/AP-1/NF-kB", "Hypoxia", "Cell cycle"))
)
competition$support_label <- paste0(
  competition$residual_positive_samples, "/", competition$samples, " positive"
)

p_competition <- ggplot(competition,
                        aes(x = competitor_label, y = residual_median)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = palette_contract["neutral_dark"]) +
  geom_col(fill = palette_contract["negative"], width = 0.62) +
  geom_text(aes(y = residual_median - 0.00055, label = support_label),
            hjust = 0, size = 2.5, family = "Arial",
            colour = palette_contract["neutral_dark"]) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.14, 0.12))) +
  labs(
    title = "Integrin specificity is not robust to competing state programs",
    subtitle = "Residual median after competitor adjustment; n = 12 paired sample-state units",
    x = "Competitor module", y = "Adjusted core-module residual median"
  ) +
  theme(plot.margin = margin(5.5, 18, 5.5, 5.5),
        axis.text.y = element_text(size = 6.5))

save_pub_r(
  p_competition,
  file.path(result_dir, "14C_Figure4_competition_robustness_boundary"),
  width_mm = 183, height_mm = 92
)

# ----------------------------- Tables ---------------------------------------
table1 <- data.frame(
  source = c("GSE273293", "PRJNA607098", "GSE130973"),
  context = c(
    "Human fascia/muscle single-cell source",
    "Skin fibroblast atlas; F7/non-F7 paired sample-state audit",
    "Independent fibroblast source; cluster-level external triangulation"
  ),
  biological_unit = c("runs/clinical groups; source audit", "paired sample-state units", "subjects"),
  n_unit = c("14 runs; 10 GMC + 4 controls", "12 eligible paired units", "5 subjects"),
  role = c("Deep fascia recovery / primary study source", "Primary source-specific validation", "Independent external validation"),
  key_boundary = c(
    "Processed/raw-data recovery does not by itself establish mechanism",
    "10/15 frozen candidates; integrin 5/6; PIEZO2 not supported in Step 08C2",
    "Integrin 5/6; actomyosin 1/4; PIEZO2 direction positive; descriptive rule"
  ),
  stringsAsFactors = FALSE
)
write.csv(table1, file.path(result_dir, "14C_Table1_dataset_provenance_v1.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

evidence_value <- function(item) {
  hit <- evidence_summary$value[evidence_summary$item == item]
  if (length(hit) == 0) NA_character_ else as.character(hit[[1]])
}

table2 <- data.frame(
  claim_id = claims$claim_id,
  claim_domain = claims$claim_domain,
  evidence_status = claims$evidence_status,
  permitted_interpretation = claims$permitted_interpretation,
  upgrade_requirement = claims$upgrade_requirement,
  stringsAsFactors = FALSE
)
table2$project_grade <- evidence_value("final_evidence_grade")
table2$source_independence_gate <- evidence_value("source_independence_gate")
write.csv(table2, file.path(result_dir, "14C_Table2_frozen_claim_evidence_boundaries_v1.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(endpoints,
          file.path(result_dir, "14C_Table3_functional_validation_endpoints_v1.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

manifest <- data.frame(
  figure_or_table = c(
    "Figure 2", "Figure 3", "Figure 4", "Table 1", "Table 2", "Table 3"
  ),
  panel_or_content = c(
    "Source-separated module support",
    "Frozen candidate support across sources",
    "Competition robustness boundary",
    "Dataset/provenance summary",
    "Frozen claims and evidence boundaries",
    "Functional validation endpoint matrix"
  ),
  source_data = c(
    primary_module_path,
    paste(primary_candidate_path, external_candidate_path, sep = " | "),
    competition_path,
    "09A metadata audit; 08C1 sample audit; 11F external audit",
    paste(claims_path, evidence_summary_path, sep = " | "),
    endpoint_path
  ),
  statistical_unit = c(
    "genes within source-specific modules",
    "12 paired sample-state units vs 5 subjects; source-specific rules",
    "12 paired sample-state units",
    "source-level audit units",
    "claim/evidence registry",
    "donor-replicated functional validation design"
  ),
  caveat = c(
    "Fractions are not pooled across sources",
    "Effect labels are descriptive and not on a common expression scale",
    "Negative residuals do not prove absence of a mechanism",
    "Provenance table is not a meta-analysis",
    "Overall evidence grade remains CAUTION",
    "Endpoints are recommendations, not completed experiments"
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(result_dir, "14C_figure_source_data_manifest_v1.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------ Contract and legends ------------------------------
contract_lines <- c(
  "# Step 14C figure contract",
  "",
  "- Core conclusion: Human fibroblast-associated states show a reproducible multicomponent ECM/integrin-cytoskeletal program, but current data do not establish mechanosensitivity-specificity or a universal PIEZO2 mechanism.",
  "- Archetype: Quantitative grid with source-separated validation panels.",
  "- Backend: R only; all exports are SVG, PDF, 600-dpi TIFF, and 300-dpi PNG preview.",
  "- Final width: 183 mm; figure heights are 92–125 mm.",
  "- Evidence hierarchy: primary paired-source evidence, independent external triangulation, then competition robustness boundary.",
  "- Panel map:",
  "  - Figure 2: module-level support by source, with source-specific denominators.",
  "  - Figure 3: frozen candidate support cross-source; tile colour is support status and text is the source-specific median directional difference.",
  "  - Figure 4: competitor-adjusted residual medians and positive-sample counts.",
  "- Statistics: no new hypothesis test is introduced by this script. The plot labels preserve the original biological units and support rules.",
  "- Source data: 14C_figure_source_data_manifest_v1.csv and the listed frozen result tables.",
  "- Image integrity: no microscopy, image processing, interpolation, or selective cell-level exclusion is performed.",
  "- Reviewer risk: source-specific definitions must not be described as pooled replication; the overall evidence grade remains CAUTION."
)
writeLines(contract_lines, file.path(result_dir, "14C_figure_contract_v1.md"), useBytes = TRUE)

legend_lines <- c(
  "# Step 14C figure legend drafts",
  "",
  "## Figure 2 | Module-level support by source",
  "Source-separated module support is shown for PRJNA607098 (12 paired sample-state units) and GSE130973 (5 subjects). PRJNA607098 bars use the strict sample-level support rule applied to the frozen candidate module; GSE130973 bars show the descriptive positive-subject fraction used in the external triangulation. Numerators and denominators are printed above bars. The sources were not pooled, and the fractions are not equivalent confirmatory tests.",
  "",
  "## Figure 3 | Candidate support across sources",
  "Each tile represents one frozen candidate in one source. Colour indicates the source-specific support status. Text reports the source-specific median directional difference and is not intended for cross-source magnitude comparison. The PRJNA607098 column uses the strict 12-unit rule; the GSE130973 column uses the descriptive 4/5 directional rule. The figure therefore displays concordance and discordance rather than an independent replication claim.",
  "",
  "## Figure 4 | Competition robustness boundary",
  "Bars show the competitor-adjusted residual median for the integrin/focal-adhesion core module in PRJNA607098 across four competing state programs (n = 12 paired sample-state units). Text reports the number of positive residual units. Negative residual medians and failed 8/12 support do not prove that integrin biology is absent; they indicate that the current data do not establish specificity from the competing states.",
  "",
  "## General statistical note",
  "The figures summarize frozen outputs from the completed computational audit. No additional hypothesis test, pooling, imputation, or cell-level inferential unit is introduced at the plotting stage. Functional validation remains outstanding, and the final evidence grade is CAUTION."
)
writeLines(legend_lines, file.path(result_dir, "14C_figure_legend_drafts_v1.md"), useBytes = TRUE)

message("Step 14C figure/table generation completed.")
message("Output directory: ", result_dir)
message("Figures: Figure 2 module support; Figure 3 candidate cross-source; Figure 4 competition boundary.")
message("Tables and source manifest written; inspect PNG previews and then review SVG/PDF text at final size.")
