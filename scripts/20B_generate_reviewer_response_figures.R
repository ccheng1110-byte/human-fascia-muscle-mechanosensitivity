# Step 20B: reviewer-response figures with explicit inferential units
#
# Backend: R / ggplot2.
# Contract:
# - do not treat cells as biological replicates;
# - separate sources with different denominators;
# - show raw sample/subject values wherever available;
# - display negative and competitor results;
# - do not change the frozen candidate panel or decision thresholds.

options(stringsAsFactors = FALSE)

project_dir <- "."
output_dir <- file.path(
  project_dir, "results", "11_manuscript_preparation",
  "20B_reviewer_response_figures"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "svglite", "ragg")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install once with install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages(library(ggplot2))

read_csv_checked <- function(path) {
  if (!file.exists(path)) stop("Required input not found: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_csv_clean <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

as_logical_strict <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes")
}

module_label <- function(x) {
  labels <- c(
    mechanosensor_channels = "Mechanosensor channels",
    integrin_focal_adhesion = "Integrin/focal adhesion",
    hippo_yap_taz = "Hippo/YAP/TAZ",
    actomyosin_rho = "Actomyosin/Rho",
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

module_order <- c(
  "Mechanosensor channels", "Integrin/focal adhesion", "Hippo/YAP/TAZ",
  "Actomyosin/Rho", "ECM remodeling", "TGF/fibrosis",
  "Inflammation/AP-1/NF-kB", "Hypoxia", "Cell cycle"
)

palette <- c(
  "PRJNA607098 (12 sample units)" = "#2B6CB0",
  "GSE130973 (5 subjects)" = "#2A9D8F",
  "Primary module" = "#2B6CB0",
  "Competitor module" = "#E58E26",
  "Pass" = "#2A9D8F",
  "Fail" = "#C44536",
  "GM09503" = "#6A51A3",
  "GM08401" = "#D95F0E"
)

theme_submission <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.background = element_rect(fill = "#F3F3F3", colour = "#BDBDBD"),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 0.2),
      panel.grid = element_blank()
    )
}

save_figure <- function(plot, filename, width_mm = 183, height_mm = 115) {
  stem <- file.path(output_dir, filename)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Arial"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()
}

median_iqr <- function(x) {
  data.frame(
    y = stats::median(x, na.rm = TRUE),
    ymin = unname(stats::quantile(x, 0.25, na.rm = TRUE)),
    ymax = unname(stats::quantile(x, 0.75, na.rm = TRUE))
  )
}

# -----------------------------------------------------------------------------
# Figure 1: evidence workflow and inferential units
# -----------------------------------------------------------------------------

workflow <- data.frame(
  x = c(1, 2, 3, 3, 2, 1),
  y = c(2, 2, 2, 1, 1, 1),
  title = c(
    "Discovery", "Primary paired\nvalidation", "Independent\ntriangulation",
    "Mechanical\nperturbation", "Regulatory /\nstiffness context",
    "Bounded\nevidence\nsynthesis"
  ),
  detail = c(
    "GSE173252\nprocessed scRNA-seq",
    "PRJNA607098\n12 sample units\nF7 vs non-F7",
    "GSE130973\n5 subjects\nstate vs other clusters",
    "GSE300230\n2 cell lines\ntension vs relaxed",
    "GSE338388; GSE123100;\nGSE276045",
    "Mechanotransduction-compatible;\nnot mechanically specific"
  ),
  class = c("Discovery", "Primary", "External", "Perturbation", "Context", "Synthesis")
)

workflow_colours <- c(
  Discovery = "#D9EAF7", Primary = "#9ECAE1", External = "#9FD8CF",
  Perturbation = "#F6C98F", Context = "#E7D8F4", Synthesis = "#F2B6B0"
)

workflow_segments <- data.frame(
  x = c(1, 2, 3, 3, 2),
  xend = c(2, 3, 3, 2, 1),
  y = c(2, 2, 2, 1, 1),
  yend = c(2, 2, 1, 1, 1)
)
figure1 <- ggplot(workflow, aes(x, y)) +
  geom_segment(
    data = workflow_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE, linewidth = 0.55, colour = "#6B6B6B",
    arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")
  ) +
  geom_label(
    aes(label = paste0(title, "\n", detail), fill = class),
    size = 2.05, linewidth = 0.35, lineheight = 0.96,
    label.padding = grid::unit(2.2, "mm")
  ) +
  scale_fill_manual(values = workflow_colours, guide = "none") +
  coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0.45, 2.55), clip = "off") +
  labs(
    title = "Evidence workflow and unit of analysis",
    subtitle = "Each source retains its own biological or inferential denominator"
  ) +
  theme_void(base_family = "Arial", base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    plot.subtitle = element_text(size = 7.5),
    plot.margin = margin(8, 8, 8, 8)
  )

save_figure(figure1, "20B_Figure1_evidence_workflow", height_mm = 105)
write_csv_clean(workflow, file.path(output_dir, "20B_Figure1_source_data.csv"))

# -----------------------------------------------------------------------------
# Figure 2: raw sample/subject module contrasts, separated by source
# -----------------------------------------------------------------------------

prjna_module <- read_csv_checked(file.path(
  project_dir, "results", "08_cross_tissue_validation", "10C_frozen_PRJNA607098",
  "PRJNA607098_F7_module_paired_contrasts_v2.csv"
))
prjna_module$source <- "PRJNA607098 (12 sample units)"
prjna_module$unit_id <- prjna_module$sample_id
prjna_module$difference <- prjna_module$module_difference

gse_gene <- read_csv_checked(file.path(
  project_dir, "results", "09_independent_external_source_screening",
  "11E_GSE130973_frozen_program_audit",
  "GSE130973_fibroblast_state_gene_contrasts_v1.csv"
))
registry <- read_csv_checked(file.path(
  project_dir, "config", "mechanotransduction_module_registry_v2.csv"
))

registry_long <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
  data.frame(
    module = registry$module[[i]],
    role = registry$role[[i]],
    gene = trimws(strsplit(registry$genes[[i]], ";", fixed = TRUE)[[1]]),
    stringsAsFactors = FALSE
  )
}))
gse_join <- merge(gse_gene, registry_long, by = "gene", all = FALSE)
gse_module <- aggregate(
  mean_difference ~ sample_id + module + role,
  data = gse_join, FUN = mean, na.rm = TRUE
)
names(gse_module)[names(gse_module) == "mean_difference"] <- "difference"
gse_module$source <- "GSE130973 (5 subjects)"
gse_module$unit_id <- gse_module$sample_id

primary_modules <- c(
  "mechanosensor_channels", "integrin_focal_adhesion",
  "hippo_yap_taz", "actomyosin_rho"
)
figure2_data <- rbind(
  prjna_module[prjna_module$module %in% primary_modules,
               c("source", "unit_id", "module", "difference")],
  gse_module[gse_module$module %in% primary_modules,
             c("source", "unit_id", "module", "difference")]
)
figure2_data$module_label <- factor(
  module_label(figure2_data$module), levels = module_order
)

set.seed(20260827)
figure2 <- ggplot(
  figure2_data,
  aes(x = module_label, y = difference, colour = source)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2, colour = "#555555") +
  geom_jitter(width = 0.10, height = 0, size = 1.8, alpha = 0.82) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48,
               linewidth = 0.55, colour = "black") +
  stat_summary(fun.data = median_iqr, geom = "errorbar", width = 0.18,
               linewidth = 0.55, colour = "black") +
  facet_wrap(~source, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = palette[c(
    "PRJNA607098 (12 sample units)", "GSE130973 (5 subjects)"
  )], guide = "none") +
  labs(
    title = "Primary-module contrasts at the correct inferential unit",
    subtitle = "Points are reconstructed sample units or subjects; bars show median and IQR",
    x = NULL,
    y = "Target-state minus comparator module score"
  ) +
  theme_submission() +
  theme(axis.text.x = element_text(angle = 32, hjust = 1))

save_figure(figure2, "20B_Figure2_sample_subject_module_contrasts")
write_csv_clean(figure2_data, file.path(output_dir, "20B_Figure2_source_data.csv"))

# -----------------------------------------------------------------------------
# Figure 3: frozen candidate contrasts with raw inferential-unit points
# -----------------------------------------------------------------------------

prjna_gene <- read_csv_checked(file.path(
  project_dir, "results", "08_cross_tissue_validation", "10C_frozen_PRJNA607098",
  "PRJNA607098_F7_frozen_gene_contrasts_v2.csv"
))
if (!"gene_difference" %in% names(prjna_gene)) {
  candidate_difference_names <- intersect(
    c("difference", "mean_difference", "expression_difference"), names(prjna_gene)
  )
  if (length(candidate_difference_names) != 1L) {
    stop("Could not resolve the PRJNA607098 candidate-difference column.")
  }
  prjna_gene$gene_difference <- prjna_gene[[candidate_difference_names[[1L]]]]
}
if (!"sample_id" %in% names(prjna_gene)) stop("PRJNA candidate file lacks sample_id.")

candidate_panel <- read_csv_checked(file.path(
  project_dir, "config", "frozen_candidate_panel_v2.csv"
))
candidate_genes <- unique(candidate_panel$gene)

prjna_candidate <- data.frame(
  source = "PRJNA607098 (12 sample units)",
  unit_id = prjna_gene$sample_id,
  gene = prjna_gene$gene,
  difference = prjna_gene$gene_difference,
  stringsAsFactors = FALSE
)
gse_candidate <- data.frame(
  source = "GSE130973 (5 subjects)",
  unit_id = gse_gene$sample_id,
  gene = gse_gene$gene,
  difference = gse_gene$mean_difference,
  stringsAsFactors = FALSE
)
figure3_data <- rbind(
  prjna_candidate[prjna_candidate$gene %in% candidate_genes, ],
  gse_candidate[gse_candidate$gene %in% candidate_genes, ]
)
gene_medians <- aggregate(
  difference ~ gene, figure3_data[figure3_data$source ==
                                    "PRJNA607098 (12 sample units)", ],
  median
)
gene_levels <- gene_medians$gene[order(gene_medians$difference)]
figure3_data$gene <- factor(figure3_data$gene, levels = gene_levels)

set.seed(20260827)
figure3 <- ggplot(figure3_data, aes(x = gene, y = difference, colour = source)) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2, colour = "#555555") +
  geom_jitter(width = 0.11, height = 0, size = 1.25, alpha = 0.72) +
  stat_summary(fun = median, geom = "point", shape = 18, size = 2.3,
               colour = "black") +
  stat_summary(fun.data = median_iqr, geom = "errorbar", width = 0.15,
               linewidth = 0.5, colour = "black") +
  facet_wrap(~source, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = palette[c(
    "PRJNA607098 (12 sample units)", "GSE130973 (5 subjects)"
  )], guide = "none") +
  labs(
    title = "Frozen candidate effects are source dependent",
    subtitle = "All 15 candidates are retained; diamonds show medians and bars show IQR",
    x = NULL, y = "Target-state minus comparator expression"
  ) +
  theme_submission() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(figure3, "20B_Figure3_frozen_candidate_raw_contrasts", height_mm = 145)
write_csv_clean(figure3_data, file.path(output_dir, "20B_Figure3_source_data.csv"))

# -----------------------------------------------------------------------------
# Figure 4: sample-level competition residuals
# -----------------------------------------------------------------------------

wide_modules <- reshape(
  prjna_module[, c("sample_id", "module", "module_difference")],
  idvar = "sample_id", timevar = "module", direction = "wide"
)
names(wide_modules) <- sub("^module_difference\\.", "", names(wide_modules))
cores <- c("integrin_focal_adhesion", "actomyosin_rho")
competitors <- c(
  "ecm_remodeling", "tgf_fibrosis", "inflammation_ap1_nfkb",
  "hypoxia", "cell_cycle"
)

residual_rows <- list()
k <- 1L
for (core in cores) {
  for (competitor in competitors) {
    fit <- stats::lm(wide_modules[[core]] ~ wide_modules[[competitor]])
    values <- stats::residuals(fit)
    positive_n <- sum(values > 0)
    residual_rows[[k]] <- data.frame(
      sample_id = wide_modules$sample_id,
      core_module = core,
      competitor_module = competitor,
      residual = values,
      residual_positive_samples = positive_n,
      gate = ifelse(positive_n >= 8L, "Pass", "Fail"),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
figure4_data <- do.call(rbind, residual_rows)
figure4_data$core_label <- module_label(figure4_data$core_module)
figure4_data$competitor_label <- factor(
  module_label(figure4_data$competitor_module), levels = module_order
)

set.seed(20260827)
figure4 <- ggplot(
  figure4_data,
  aes(x = competitor_label, y = residual, colour = gate)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2, colour = "#555555") +
  geom_jitter(width = 0.10, height = 0, size = 1.55, alpha = 0.80) +
  stat_summary(fun = median, geom = "crossbar", width = 0.46,
               linewidth = 0.55, colour = "black") +
  facet_wrap(~core_label, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = palette[c("Pass", "Fail")]) +
  labs(
    title = "Competition residuals expose the specificity boundary",
    subtitle = "Gate requires at least 8 of 12 residuals above zero",
    x = "Competing program", y = "Core-module residual"
  ) +
  theme_submission() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_figure(figure4, "20B_Figure4_competition_residuals", height_mm = 130)
write_csv_clean(figure4_data, file.path(output_dir, "20B_Figure4_source_data.csv"))

# -----------------------------------------------------------------------------
# Figure 5: direct GSE300230 mechanical-tension module response
# -----------------------------------------------------------------------------

gse300_score <- read_csv_checked(file.path(
  project_dir, "results", "12_computational_strengthening",
  "15D_GSE300230_score_sensitivity_and_negative_controls",
  "GSE300230_module_score_contrasts_by_method_v1.csv"
))
gse300_camera <- read_csv_checked(file.path(
  project_dir, "results", "12_computational_strengthening",
  "15C_GSE300230_mechanochemical_perturbation",
  "GSE300230_frozen_module_camera_statistics_v1.csv"
))

primary_contrast <- "mechanical_tension_vs_relaxed_tgfb_absent"
selected_modules <- c(
  "integrin_focal_adhesion", "actomyosin_rho", "ecm_remodeling",
  "tgf_fibrosis", "inflammation_ap1_nfkb", "cell_cycle"
)
figure5_data <- gse300_score[
  gse300_score$method == "zmean" &
    gse300_score$contrast == primary_contrast &
    gse300_score$module %in% selected_modules,
  c("module", "role", "cell_line", "contrast_value")
]
camera_primary <- gse300_camera[
  gse300_camera$contrast == primary_contrast &
    gse300_camera$module %in% selected_modules,
  c("module", "cell_line", "Direction", "FDR")
]
figure5_data <- merge(
  figure5_data, camera_primary,
  by = c("module", "cell_line"), all.x = TRUE
)
figure5_data$module_label <- factor(
  module_label(figure5_data$module), levels = module_order
)
figure5_data$q_label <- ifelse(
  figure5_data$FDR < 0.001,
  format(figure5_data$FDR, scientific = TRUE, digits = 2),
  sprintf("%.3f", figure5_data$FDR)
)

figure5 <- ggplot(
  figure5_data,
  aes(x = module_label, y = contrast_value, colour = cell_line, group = cell_line)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2, colour = "#555555") +
  geom_line(linewidth = 0.55, alpha = 0.70,
            position = position_dodge(width = 0.18)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.18)) +
  geom_text(
    aes(label = paste0("q=", q_label)),
    position = position_dodge(width = 0.18), hjust = -0.05, size = 2.2,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = palette[c("GM09503", "GM08401")]) +
  labs(
    title = "Direct mechanical tension activates actomyosin and cell-cycle programs",
    subtitle = "GSE300230: model-estimated tension minus relaxed contrast without TGF-beta",
    x = NULL, y = "Module-score contrast (z-mean)"
  ) +
  theme_submission() +
  theme(axis.text.x = element_text(angle = 32, hjust = 1))

save_figure(figure5, "20B_Figure5_GSE300230_tension_response", height_mm = 120)
write_csv_clean(figure5_data, file.path(output_dir, "20B_Figure5_source_data.csv"))

# -----------------------------------------------------------------------------
# Legends and manifest
# -----------------------------------------------------------------------------

legend_lines <- c(
  "# Step 20B figure legends",
  "",
  "## Figure 1 | Evidence workflow and inferential units",
  "The staged design separates discovery, paired atlas validation, independent subject-level triangulation, direct mechanical perturbation, regulatory/stiffness context and bounded synthesis. Boxes identify the denominator used by each source.",
  "",
  "## Figure 2 | Primary-module contrasts at the correct inferential unit",
  "PRJNA607098 points are 12 reconstructed SRS sample units comparing F7 with pooled non-F7 cells within sample. GSE130973 points are five subjects comparing marker-defined candidate fibroblast clusters with other clusters within subject. Black summaries are median and IQR. Free y-axis scales prevent the two source-specific normalized-expression scales from being treated as directly commensurate.",
  "",
  "## Figure 3 | Frozen candidate effects are source dependent",
  "All frozen candidates are shown without post hoc removal. Points are sample- or subject-level target-minus-comparator contrasts; black diamonds and bars show median and IQR. Source separation preserves the distinct expression scales and denominators.",
  "",
  "## Figure 4 | Competition residuals expose the specificity boundary",
  "Within PRJNA607098, each core-module contrast was regressed on each competitor-module contrast. Points are sample-level residuals. The frozen gate required at least 8 of 12 residuals above zero. Integrin/focal adhesion fails against several competitors, whereas actomyosin/Rho is more robust.",
  "",
  "## Figure 5 | Direct mechanical tension activates actomyosin and cell-cycle programs",
  "Model-estimated module-score contrasts compare mechanical tension with relaxed collagen without TGF-beta in each GSE300230 cell line. Labels give CAMERA BH q values. The concurrent cell-cycle response is displayed as an interpretation boundary. The two cell lines are not treated as independent donors or as replicates for an age effect."
)
writeLines(
  legend_lines,
  file.path(output_dir, "20B_figure_legends_v1.md"),
  useBytes = TRUE
)

manifest <- data.frame(
  figure = paste0("Figure ", 1:5),
  stem = c(
    "20B_Figure1_evidence_workflow",
    "20B_Figure2_sample_subject_module_contrasts",
    "20B_Figure3_frozen_candidate_raw_contrasts",
    "20B_Figure4_competition_residuals",
    "20B_Figure5_GSE300230_tension_response"
  ),
  inferential_unit = c(
    "workflow", "12 sample units / 5 subjects", "12 sample units / 5 subjects",
    "12 sample units", "two separately modeled cell lines"
  ),
  formats = "PNG;SVG;PDF;TIFF",
  stringsAsFactors = FALSE
)
write_csv_clean(manifest, file.path(output_dir, "20B_figure_manifest_v1.csv"))

message("Step 20B reviewer-response figure generation completed.")
message("Output directory: ", output_dir)
message("Inspect all PNG previews, then verify SVG/PDF text at final size.")
