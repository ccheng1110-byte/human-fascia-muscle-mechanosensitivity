project_dir <- "."
input_file <- file.path(
  project_dir,
  "results/11_manuscript_preparation/21A_submission_strengthening_revision",
  "21A_evidence_hierarchy_source_data_v1.csv"
)
output_dir <- file.path(
  project_dir,
  "results/11_manuscript_preparation/21B_evidence_hierarchy_figure"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install the missing R package(s) first: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)

evidence <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

source_levels <- c("PRJNA607098", "GSE130973", "GSE300230", "GSE338388", "GSE276045")
evidence_levels <- c("Module level", "Candidate level")
support_levels <- c(
  "Strong module support",
  "Partial/conditional",
  "Discordant/source-dependent",
  "Not confirmatory"
)

evidence$source <- factor(evidence$source, levels = source_levels)
evidence$evidence_level <- factor(evidence$evidence_level, levels = rev(evidence_levels))
evidence$support_class <- factor(evidence$support_class, levels = support_levels)
evidence$wrapped_text <- vapply(
  evidence$display_text,
  function(x) paste(strwrap(x, width = 24), collapse = "\n"),
  character(1)
)
evidence$text_color <- ifelse(
  evidence$support_class %in% c("Strong module support", "Discordant/source-dependent"),
  "white",
  "#17212B"
)

fill_values <- c(
  "Strong module support" = "#2B8C6B",
  "Partial/conditional" = "#E6A83E",
  "Discordant/source-dependent" = "#C95B5B",
  "Not confirmatory" = "#C7CDD4"
)

p <- ggplot(evidence, aes(x = source, y = evidence_level, fill = support_class)) +
  geom_tile(color = "white", linewidth = 1.1, width = 0.98, height = 0.94) +
  geom_text(
    aes(label = wrapped_text, color = text_color),
    size = 2.18,
    lineheight = 0.94
  ) +
  scale_fill_manual(values = fill_values, drop = FALSE) +
  scale_color_identity() +
  scale_x_discrete(position = "top") +
  labs(
    title = "Evidence is more reproducible at module than individual-gene level",
    subtitle = "Module support is more reproducible; candidate genes and specificity remain source dependent",
    x = NULL,
    y = NULL,
    fill = "Evidence class",
    caption = paste(
      "Colours summarize prespecified or explicitly bounded evidence classes;",
      "they are not a formal universal grading scale."
    )
  ) +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "#17212B", size = 8.2),
    axis.text.y = element_text(face = "bold", color = "#17212B", size = 8.2),
    plot.title = element_text(face = "bold", size = 12.5, color = "#17212B", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 8.3, color = "#425466", margin = margin(b = 9)),
    plot.caption = element_text(size = 7, color = "#566573", hjust = 0, margin = margin(t = 7)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 7.5),
    legend.text = element_text(size = 7.2),
    legend.key.height = grid::unit(3.4, "mm"),
    legend.key.width = grid::unit(7.5, "mm"),
    plot.margin = margin(8, 8, 7, 8)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

base_name <- file.path(output_dir, "21B_Figure6_evidence_hierarchy")
ggsave(paste0(base_name, ".png"), p, width = 183, height = 105, units = "mm", dpi = 300, bg = "white")
ggsave(paste0(base_name, ".pdf"), p, width = 183, height = 105, units = "mm", device = cairo_pdf, bg = "white")

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(base_name, ".svg"), p, width = 183, height = 105, units = "mm", device = svglite::svglite, bg = "white")
} else {
  ggsave(paste0(base_name, ".svg"), p, width = 183, height = 105, units = "mm", device = "svg", bg = "white")
}

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(paste0(base_name, ".tiff"), p, width = 183, height = 105, units = "mm", dpi = 600, device = ragg::agg_tiff, bg = "white", compression = "lzw")
} else {
  ggsave(paste0(base_name, ".tiff"), p, width = 183, height = 105, units = "mm", dpi = 600, device = "tiff", bg = "white", compression = "lzw")
}

write.csv(
  evidence[c("source", "evidence_level", "support_class", "display_text", "evidence_boundary")],
  file.path(output_dir, "21B_Figure6_source_data_v1.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

legend_file <- file.path(output_dir, "21B_Figure6_legend_v1.md")
writeLines(
  c(
    "# Figure 6 legend",
    "",
    "**Figure 6 | Module-level support is more reproducible than individual-gene support across sources.**",
    "The matrix summarizes the frozen evidence hierarchy across the primary paired atlas validation, an independent skin source, direct mechanical tension, TGF-beta/TEAD regulatory-axis cross-validation and cross-tissue stiffness analysis. Strong module-level actomyosin evidence contrasts with weaker or source-dependent individual-gene concordance. Partial or discordant cells retain the principal evidence boundaries, including co-active ECM, TGF/fibrosis or cell-cycle programs, limited biological replication and the absence of a mechanical condition in GSE338388. Colours are descriptive evidence classes and not a universal grading scale."
  ),
  legend_file,
  useBytes = TRUE
)

message("Step 21B evidence-hierarchy figure generation completed.")
message("Output directory: ", output_dir)
message("Inspect the PNG preview, then verify SVG/PDF text at final size.")
