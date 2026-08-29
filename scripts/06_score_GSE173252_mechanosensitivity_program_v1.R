options(stringsAsFactors = FALSE)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required.")
}

project_dir <- "."
input_path <- file.path(
  project_dir, "data", "raw", "GSE173252",
  "GSE173252_dd_mesenchyme.rds.gz"
)
gene_set_path <- file.path(
  project_dir, "data", "gene_sets", "mechanosensitivity_program_v1.csv"
)
result_dir <- file.path(
  project_dir, "results", "04_mechanosensitivity_scoring", "GSE173252"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

random_seed <- 20260823L
control_genes_per_target <- 20L
expression_bins <- 24L

read_double_gzip_rds <- function(path) {
  outer_connection <- gzfile(path, open = "rb")
  inner_connection <- gzcon(outer_connection, text = FALSE)
  on.exit(try(close(inner_connection), silent = TRUE), add = TRUE)
  readRDS(inner_connection)
}

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

is_doublet <- function(x) {
  if (is.logical(x)) {
    answer <- x
  } else {
    answer <- tolower(trimws(as.character(x))) %in% c("true", "doublet", "yes", "1")
  }
  answer[is.na(answer)] <- FALSE
  answer
}

message("Loading GSE173252 legacy object.")
object <- read_double_gzip_rds(input_path)
object_attributes <- attributes(object)
metadata <- object_attributes[["meta.data"]]
rna_attributes <- attributes(object_attributes[["assays"]][["RNA"]])
rna_counts <- rna_attributes[["counts"]]
rna_logdata <- rna_attributes[["data"]]
gene_sets <- read.csv(gene_set_path, check.names = FALSE)

if (!identical(colnames(rna_counts), rownames(metadata))) {
  if (!all(colnames(rna_counts) %in% rownames(metadata))) {
    stop("Cell names differ between RNA counts and metadata.")
  }
  metadata <- metadata[colnames(rna_counts), , drop = FALSE]
}

required_metadata <- c(
  "orig.ident", "source", "condition", "seurat_clusters", "doublet_prediction"
)
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0L) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

identity_map <- c("0" = "FB", "1" = "Myofib", "2" = "VSMC", "3" = "Pericyte")
metadata$seurat_clusters <- as.character(metadata$seurat_clusters)
metadata$reproduced_identity <- unname(identity_map[metadata$seurat_clusters])
singlet_keep <- !is_doublet(metadata$doublet_prediction)
singlet_metadata <- metadata[singlet_keep, , drop = FALSE]

if (any(is.na(singlet_metadata$reproduced_identity))) {
  stop("Some singlet cells do not map to a reproduced coarse identity.")
}

gene_sets$present_in_GSE173252 <- gene_sets$gene %in% rownames(rna_logdata)
gene_sets$detected_in_singlets <- FALSE
present_rows <- which(gene_sets$present_in_GSE173252)
if (length(present_rows) > 0L) {
  detected_genes <- rownames(rna_counts)[
    Matrix::rowSums(rna_counts[, singlet_keep, drop = FALSE] > 0) > 0
  ]
  gene_sets$detected_in_singlets[present_rows] <-
    gene_sets$gene[present_rows] %in% detected_genes
}
safe_write_csv(
  gene_sets,
  file.path(result_dir, "GSE173252_mechanosensitivity_gene_set_inventory_v1.csv")
)

module_inventory <- do.call(rbind, lapply(split(gene_sets, gene_sets$module), function(x) {
  data.frame(
    module = x$module[1],
    role = x$role[1],
    genes_frozen = nrow(x),
    genes_present = sum(x$present_in_GSE173252),
    genes_detected = sum(x$detected_in_singlets),
    present_fraction = mean(x$present_in_GSE173252),
    detected_fraction = mean(x$detected_in_singlets),
    stringsAsFactors = FALSE
  )
}))
rownames(module_inventory) <- NULL
safe_write_csv(
  module_inventory,
  file.path(result_dir, "GSE173252_mechanosensitivity_module_inventory_v1.csv")
)

# Deterministic expression-matched controls, analogous in spirit to an
# AddModuleScore background but frozen and exported for full traceability.
average_expression <- Matrix::rowMeans(rna_logdata[, singlet_keep, drop = FALSE])
detection_fraction <- Matrix::rowMeans(rna_counts[, singlet_keep, drop = FALSE] > 0)
eligible_background <- names(average_expression)[
  detection_fraction >= 0.01 &
    !grepl("^(MT-|RPL|RPS)", names(average_expression), ignore.case = FALSE)
]
all_panel_genes <- unique(gene_sets$gene)
eligible_background <- setdiff(eligible_background, all_panel_genes)

rank_values <- rank(average_expression[eligible_background], ties.method = "average")
background_bin <- ceiling(rank_values / length(rank_values) * expression_bins)
background_bin[background_bin < 1L] <- 1L
background_bin[background_bin > expression_bins] <- expression_bins
names(background_bin) <- eligible_background

all_gene_rank <- rank(average_expression, ties.method = "average")
all_gene_bin <- ceiling(all_gene_rank / length(all_gene_rank) * expression_bins)
all_gene_bin[all_gene_bin < 1L] <- 1L
all_gene_bin[all_gene_bin > expression_bins] <- expression_bins
names(all_gene_bin) <- names(average_expression)

find_control_pool <- function(target_bin, minimum_size) {
  distance <- 0L
  repeat {
    allowed_bins <- seq.int(
      max(1L, target_bin - distance),
      min(expression_bins, target_bin + distance)
    )
    pool <- names(background_bin)[background_bin %in% allowed_bins]
    if (length(pool) >= minimum_size || distance >= expression_bins) return(pool)
    distance <- distance + 1L
  }
}

set.seed(random_seed)
control_rows <- list()
control_index <- 1L
modules <- unique(gene_sets$module)
for (module_name in modules) {
  module_genes <- unique(gene_sets$gene[
    gene_sets$module == module_name & gene_sets$detected_in_singlets
  ])
  for (target_gene in module_genes) {
    target_bin <- all_gene_bin[[target_gene]]
    pool <- find_control_pool(target_bin, control_genes_per_target)
    selected <- sample(pool, min(control_genes_per_target, length(pool)), replace = FALSE)
    control_rows[[control_index]] <- data.frame(
      gene_set_version = "v1",
      random_seed = random_seed,
      module = module_name,
      target_gene = target_gene,
      target_expression_bin = target_bin,
      control_gene = selected,
      control_expression_bin = unname(background_bin[selected]),
      stringsAsFactors = FALSE
    )
    control_index <- control_index + 1L
  }
}
matched_controls <- do.call(rbind, control_rows)
safe_write_csv(
  matched_controls,
  file.path(result_dir, "GSE173252_expression_matched_control_genes_v1.csv")
)

# Cell-level scores are summarized by biological sample and cluster. Individual
# cells are not used as independent replicates in downstream contrasts.
cell_score_matrix <- matrix(
  NA_real_, nrow = length(modules), ncol = sum(singlet_keep),
  dimnames = list(modules, rownames(singlet_metadata))
)
for (module_name in modules) {
  target_genes <- unique(gene_sets$gene[
    gene_sets$module == module_name & gene_sets$detected_in_singlets
  ])
  control_genes <- unique(matched_controls$control_gene[
    matched_controls$module == module_name
  ])
  target_score <- Matrix::colMeans(rna_logdata[target_genes, singlet_keep, drop = FALSE])
  control_score <- Matrix::colMeans(rna_logdata[control_genes, singlet_keep, drop = FALSE])
  cell_score_matrix[module_name, ] <- target_score - control_score
}

primary_modules <- unique(gene_sets$module[gene_sets$role == "primary"])
primary_composite <- Matrix::colMeans(cell_score_matrix[primary_modules, , drop = FALSE])
cell_score_matrix <- rbind(
  cell_score_matrix,
  primary_mechanotransduction_composite = primary_composite
)

cell_group <- interaction(
  singlet_metadata$orig.ident,
  singlet_metadata$seurat_clusters,
  drop = TRUE,
  sep = "__"
)
cell_summary_rows <- list()
summary_index <- 1L
for (group_name in levels(cell_group)) {
  cell_index <- which(cell_group == group_name)
  first_cell <- cell_index[1]
  for (module_name in rownames(cell_score_matrix)) {
    values <- cell_score_matrix[module_name, cell_index]
    cell_summary_rows[[summary_index]] <- data.frame(
      sample = singlet_metadata$orig.ident[first_cell],
      source = singlet_metadata$source[first_cell],
      condition = singlet_metadata$condition[first_cell],
      seurat_cluster = singlet_metadata$seurat_clusters[first_cell],
      reproduced_identity = singlet_metadata$reproduced_identity[first_cell],
      module = module_name,
      singlet_cells = length(cell_index),
      mean_cell_score = mean(values),
      median_cell_score = median(values),
      sd_cell_score = sd(values),
      stringsAsFactors = FALSE
    )
    summary_index <- summary_index + 1L
  }
}
cell_score_summary <- do.call(rbind, cell_summary_rows)
safe_write_csv(
  cell_score_summary,
  file.path(result_dir, "GSE173252_sample_cluster_cell_module_scores_v1.csv")
)

# Raw-count pseudobulk by sample x cluster.
group_factor <- factor(cell_group)
group_levels <- levels(group_factor)
design <- Matrix::sparseMatrix(
  i = seq_along(group_factor),
  j = as.integer(group_factor),
  x = 1,
  dims = c(length(group_factor), length(group_levels)),
  dimnames = list(rownames(singlet_metadata), group_levels)
)
pseudobulk_counts <- rna_counts[, singlet_keep, drop = FALSE] %*% design
group_parts <- do.call(rbind, strsplit(group_levels, "__", fixed = TRUE))
group_sample <- group_parts[, 1]
group_cluster <- group_parts[, 2]
group_identity <- unname(identity_map[group_cluster])

all_scoring_genes <- unique(c(
  gene_sets$gene[gene_sets$detected_in_singlets],
  matched_controls$control_gene
))
library_sizes <- Matrix::colSums(pseudobulk_counts)
log_cpm <- log2(
  sweep(
    as.matrix(pseudobulk_counts[all_scoring_genes, , drop = FALSE]),
    2, library_sizes, "/"
  ) * 1e6 + 0.5
)
gene_z <- t(scale(t(log_cpm), center = TRUE, scale = TRUE))
gene_z[!is.finite(gene_z)] <- 0

pseudobulk_score_matrix <- matrix(
  NA_real_, nrow = length(modules), ncol = length(group_levels),
  dimnames = list(modules, group_levels)
)
for (module_name in modules) {
  target_genes <- unique(gene_sets$gene[
    gene_sets$module == module_name & gene_sets$detected_in_singlets
  ])
  control_genes <- unique(matched_controls$control_gene[
    matched_controls$module == module_name
  ])
  pseudobulk_score_matrix[module_name, ] <-
    colMeans(gene_z[target_genes, , drop = FALSE]) -
    colMeans(gene_z[control_genes, , drop = FALSE])
}
pseudobulk_score_matrix <- rbind(
  pseudobulk_score_matrix,
  primary_mechanotransduction_composite = colMeans(
    pseudobulk_score_matrix[primary_modules, , drop = FALSE]
  )
)

pseudobulk_score_long <- do.call(rbind, lapply(rownames(pseudobulk_score_matrix), function(module_name) {
  data.frame(
    sample = group_sample,
    source = singlet_metadata$source[match(group_sample, singlet_metadata$orig.ident)],
    condition = singlet_metadata$condition[match(group_sample, singlet_metadata$orig.ident)],
    seurat_cluster = group_cluster,
    reproduced_identity = group_identity,
    module = module_name,
    pseudobulk_module_score = as.numeric(pseudobulk_score_matrix[module_name, ]),
    singlet_cells = as.integer(table(group_factor)[group_levels]),
    library_size = as.numeric(library_sizes),
    stringsAsFactors = FALSE
  )
}))
safe_write_csv(
  pseudobulk_score_long,
  file.path(result_dir, "GSE173252_sample_cluster_pseudobulk_module_scores_v1.csv")
)

# Primary donor-level comparison: paired DD Myofib vs DD FB. A secondary
# comparison uses the unweighted mean of the other three mesenchymal identities.
dd_samples <- unique(singlet_metadata$orig.ident[singlet_metadata$source == "DD"])
paired_module_rows <- list()
paired_index <- 1L
for (module_name in rownames(pseudobulk_score_matrix)) {
  differences_fb <- numeric(0)
  differences_rest <- numeric(0)
  sample_rows <- list()
  for (sample_id in dd_samples) {
    target_index <- which(group_sample == sample_id & group_identity == "Myofib")
    fb_index <- which(group_sample == sample_id & group_identity == "FB")
    rest_index <- which(group_sample == sample_id & group_identity != "Myofib")
    if (length(target_index) != 1L || length(fb_index) != 1L || length(rest_index) < 1L) next
    target_score <- pseudobulk_score_matrix[module_name, target_index]
    fb_score <- pseudobulk_score_matrix[module_name, fb_index]
    rest_score <- mean(pseudobulk_score_matrix[module_name, rest_index])
    differences_fb <- c(differences_fb, target_score - fb_score)
    differences_rest <- c(differences_rest, target_score - rest_score)
    sample_rows[[length(sample_rows) + 1L]] <- data.frame(
      module = module_name,
      sample = sample_id,
      myofib_score = target_score,
      fb_score = fb_score,
      other_identity_mean_score = rest_score,
      myofib_minus_fb = target_score - fb_score,
      myofib_minus_other_identity_mean = target_score - rest_score,
      stringsAsFactors = FALSE
    )
  }
  if (length(sample_rows) > 0L) {
    paired_module_rows[[paired_index]] <- do.call(rbind, sample_rows)
    paired_index <- paired_index + 1L
  }
}
paired_module_differences <- do.call(rbind, paired_module_rows)
safe_write_csv(
  paired_module_differences,
  file.path(result_dir, "GSE173252_DD_myofib_paired_module_differences_v1.csv")
)

paired_module_summary <- do.call(rbind, lapply(
  split(paired_module_differences, paired_module_differences$module),
  function(x) data.frame(
    module = x$module[1],
    DD_samples = nrow(x),
    sample_ids = paste(x$sample, collapse = ";"),
    median_myofib_minus_fb = median(x$myofib_minus_fb),
    minimum_myofib_minus_fb = min(x$myofib_minus_fb),
    maximum_myofib_minus_fb = max(x$myofib_minus_fb),
    positive_fraction_vs_fb = mean(x$myofib_minus_fb > 0),
    median_myofib_minus_other = median(x$myofib_minus_other_identity_mean),
    positive_fraction_vs_other = mean(x$myofib_minus_other_identity_mean > 0),
    stringsAsFactors = FALSE
  )
))
paired_module_summary <- paired_module_summary[order(
  -paired_module_summary$median_myofib_minus_fb
), , drop = FALSE]
rownames(paired_module_summary) <- NULL
safe_write_csv(
  paired_module_summary,
  file.path(result_dir, "GSE173252_DD_myofib_paired_module_summary_v1.csv")
)

# Gene-level evidence for the four primary modules only.
primary_gene_table <- unique(gene_sets[
  gene_sets$role == "primary" & gene_sets$detected_in_singlets,
  c("module", "role", "gene", "source_database", "source_id")
])
myofib_cells <- which(
  singlet_keep & metadata$reproduced_identity == "Myofib"
)
myofib_pct <- Matrix::rowMeans(rna_counts[, myofib_cells, drop = FALSE] > 0) * 100
myofib_mean <- Matrix::rowMeans(rna_logdata[, myofib_cells, drop = FALSE])

candidate_rows <- list()
candidate_index <- 1L
for (i in seq_len(nrow(primary_gene_table))) {
  gene <- primary_gene_table$gene[i]
  gene_row <- match(gene, rownames(pseudobulk_counts))
  dd_fb_difference <- numeric(0)
  dd_rest_difference <- numeric(0)
  sample_detail <- character(0)
  for (sample_id in dd_samples) {
    target_index <- which(group_sample == sample_id & group_identity == "Myofib")
    fb_index <- which(group_sample == sample_id & group_identity == "FB")
    rest_index <- which(group_sample == sample_id & group_identity != "Myofib")
    if (length(target_index) != 1L || length(fb_index) != 1L || length(rest_index) < 1L) next
    target_log_cpm <- log2(
      as.numeric(pseudobulk_counts[gene_row, target_index]) /
        library_sizes[target_index] * 1e6 + 0.5
    )
    fb_log_cpm <- log2(
      as.numeric(pseudobulk_counts[gene_row, fb_index]) /
        library_sizes[fb_index] * 1e6 + 0.5
    )
    pooled_rest_count <- sum(as.numeric(pseudobulk_counts[gene_row, rest_index]))
    pooled_rest_library <- sum(library_sizes[rest_index])
    rest_log_cpm <- log2(pooled_rest_count / pooled_rest_library * 1e6 + 0.5)
    dd_fb_difference <- c(dd_fb_difference, target_log_cpm - fb_log_cpm)
    dd_rest_difference <- c(dd_rest_difference, target_log_cpm - rest_log_cpm)
    sample_detail <- c(
      sample_detail,
      paste0(sample_id, ":", format(round(target_log_cpm - fb_log_cpm, 3), nsmall = 3))
    )
  }
  pct <- as.numeric(myofib_pct[gene])
  median_fb <- median(dd_fb_difference)
  positive_fb <- mean(dd_fb_difference > 0)
  candidate_tier <- if (
    pct >= 10 & median_fb >= 0.5 & positive_fb == 1
  ) {
    "A_consistent_DD_myofib_enrichment"
  } else if (pct >= 5 & median_fb > 0 & positive_fb >= (2 / 3)) {
    "B_partial_support"
  } else {
    "C_no_specific_enrichment"
  }
  candidate_rows[[candidate_index]] <- data.frame(
    module = primary_gene_table$module[i],
    gene = gene,
    source_database = primary_gene_table$source_database[i],
    source_id = primary_gene_table$source_id[i],
    myofib_percent_cells_expressing = pct,
    myofib_mean_log_normalized_expression = as.numeric(myofib_mean[gene]),
    DD_samples_evaluable = length(dd_fb_difference),
    median_log2CPM_myofib_minus_FB = median_fb,
    minimum_log2CPM_myofib_minus_FB = min(dd_fb_difference),
    positive_fraction_myofib_vs_FB = positive_fb,
    median_log2CPM_myofib_minus_pooled_other = median(dd_rest_difference),
    positive_fraction_myofib_vs_pooled_other = mean(dd_rest_difference > 0),
    per_sample_myofib_minus_FB = paste(sample_detail, collapse = ";"),
    candidate_tier = candidate_tier,
    stringsAsFactors = FALSE
  )
  candidate_index <- candidate_index + 1L
}
candidate_evidence <- do.call(rbind, candidate_rows)
tier_order <- c(
  "A_consistent_DD_myofib_enrichment",
  "B_partial_support",
  "C_no_specific_enrichment"
)
candidate_evidence$tier_rank <- match(candidate_evidence$candidate_tier, tier_order)
candidate_evidence <- candidate_evidence[order(
  candidate_evidence$tier_rank,
  -candidate_evidence$median_log2CPM_myofib_minus_FB,
  -candidate_evidence$myofib_percent_cells_expressing
), , drop = FALSE]
candidate_evidence$tier_rank <- NULL
rownames(candidate_evidence) <- NULL
safe_write_csv(
  candidate_evidence,
  file.path(result_dir, "GSE173252_DD_myofib_primary_mechanotransduction_gene_evidence_v1.csv")
)
safe_write_csv(
  candidate_evidence[candidate_evidence$candidate_tier != "C_no_specific_enrichment", , drop = FALSE],
  file.path(result_dir, "GSE173252_DD_myofib_mechanotransduction_candidate_shortlist_v1.csv")
)

top_candidates <- head(
  candidate_evidence[candidate_evidence$candidate_tier == "A_consistent_DD_myofib_enrichment", ],
  15
)
top_candidate_lines <- if (nrow(top_candidates) == 0L) {
  "- No tier-A genes met the pre-specified thresholds."
} else {
  apply(top_candidates, 1, function(x) paste0(
    "- ", x[["gene"]], " [", x[["module"]], "]: median Myofib-FB log2CPM difference = ",
    format(round(as.numeric(x[["median_log2CPM_myofib_minus_FB"]]), 3), nsmall = 3),
    "; expressing Myofib cells = ",
    format(round(as.numeric(x[["myofib_percent_cells_expressing"]]), 1), nsmall = 1), "%"
  ))
}

module_lines <- apply(paired_module_summary, 1, function(x) paste0(
  "- ", x[["module"]], ": median paired Myofib-FB score difference = ",
  format(round(as.numeric(x[["median_myofib_minus_fb"]]), 3), nsmall = 3),
  "; positive in ",
  format(round(as.numeric(x[["positive_fraction_vs_fb"]]) * 100, 0), nsmall = 0),
  "% of DD samples"
))

missing_gene_count <- sum(!gene_sets$present_in_GSE173252)
report_lines <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: run",
  paste0("- Origin Date: ", format(Sys.Date(), "%Y-%m-%d")),
  "- Verification Status: ANALYZED",
  "- Version Label: mechanosensitivity_scoring_v1",
  "",
  "## GSE173252 mechanosensitivity scoring result",
  "",
  "### Execution status",
  "",
  "- Status: completed",
  paste0("- Singlets analyzed: ", sum(singlet_keep)),
  paste0("- Predicted doublets excluded: ", sum(!singlet_keep)),
  paste0("- Frozen gene-set rows: ", nrow(gene_sets)),
  paste0("- Gene-set rows absent from this matrix: ", missing_gene_count),
  paste0("- Random seed: ", random_seed),
  "",
  "### Paired DD Myofib versus FB module evidence",
  "",
  module_lines,
  "",
  "### Top tier-A primary mechanotransduction candidates",
  "",
  top_candidate_lines,
  "",
  "### Interpretation boundary",
  "",
  "- Scores represent transcriptional programs, not direct measurements of cellular mechanical sensitivity.",
  "- DD1-DD3 are the only samples with a sufficiently represented Myofib population; therefore the paired analysis has n = 3 biological samples.",
  "- Effect size and direction consistency are reported. No inferential p-value is claimed because n = 3 cannot support a stable formal test.",
  "- ECM, TGF/fibrosis, inflammation, hypoxia and cell-cycle modules are explicit competing explanations and must be interpreted alongside the four primary modules.",
  "- A Myofib-enriched transcript is a candidate marker, not evidence that the gene causally mediates mechanosensitivity.",
  "",
  "### Frozen source anchors",
  "",
  "- GO/MSigDB response to mechanical stimulus: GO:0009612",
  "- Reactome Integrin cell surface interactions: R-HSA-216083",
  "- Reactome shear-stress response: R-HSA-9860931",
  "- Reactome YAP1/WWTR1 transcription: R-HSA-2032785",
  "- Reactome RHO-ROCK signaling: R-HSA-5627117",
  "- Reactome extracellular matrix organization: R-HSA-1474244",
  "- Reactome TGF-beta receptor signaling: R-HSA-170834"
)
writeLines(
  report_lines,
  file.path(result_dir, "GSE173252_mechanosensitivity_scoring_report_v1.md"),
  useBytes = TRUE
)

message("Mechanosensitivity program v1 scoring completed.")
message("Results: ", result_dir)
