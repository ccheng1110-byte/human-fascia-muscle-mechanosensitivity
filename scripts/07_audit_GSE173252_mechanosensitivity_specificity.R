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
step06_dir <- file.path(
  project_dir, "results", "04_mechanosensitivity_scoring", "GSE173252"
)
module_score_path <- file.path(
  step06_dir, "GSE173252_sample_cluster_pseudobulk_module_scores_v1.csv"
)
candidate_path <- file.path(
  step06_dir,
  "GSE173252_DD_myofib_primary_mechanotransduction_gene_evidence_v1.csv"
)
result_dir <- file.path(
  project_dir, "results", "05_specificity_robustness", "GSE173252"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(input_path, gene_set_path, module_score_path, candidate_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = "; "))
}

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

safe_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_positive_fraction <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x > 0, na.rm = TRUE)
}

message("Loading GSE173252 object for specificity and robustness audit.")
object <- read_double_gzip_rds(input_path)
object_attributes <- attributes(object)
metadata <- object_attributes[["meta.data"]]
rna_attributes <- attributes(object_attributes[["assays"]][["RNA"]])
rna_counts <- rna_attributes[["counts"]]
rna_logdata <- rna_attributes[["data"]]

if (!identical(colnames(rna_counts), rownames(metadata))) {
  if (!all(colnames(rna_counts) %in% rownames(metadata))) {
    stop("RNA matrix cell names do not match metadata row names.")
  }
  metadata <- metadata[colnames(rna_counts), , drop = FALSE]
}

identity_map <- c("0" = "FB", "1" = "Myofib", "2" = "VSMC", "3" = "Pericyte")
metadata$seurat_clusters <- as.character(metadata$seurat_clusters)
metadata$reproduced_identity <- unname(identity_map[metadata$seurat_clusters])
singlet_keep <- !is_doublet(metadata$doublet_prediction)
singlet_metadata <- metadata[singlet_keep, , drop = FALSE]

gene_sets <- read.csv(gene_set_path, check.names = FALSE)
module_scores <- read.csv(module_score_path, check.names = FALSE)
candidate_evidence <- read.csv(candidate_path, check.names = FALSE)

expected_score_columns <- c(
  "sample", "source", "seurat_cluster", "reproduced_identity", "module",
  "pseudobulk_module_score"
)
missing_score_columns <- setdiff(expected_score_columns, colnames(module_scores))
if (length(missing_score_columns) > 0L) {
  stop("Step-06 module score file is missing columns: ", paste(missing_score_columns, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 1. Module identity ranks within every biological sample
# -----------------------------------------------------------------------------
rank_rows <- list()
rank_index <- 1L
sample_module_groups <- split(
  module_scores,
  interaction(module_scores$sample, module_scores$module, drop = TRUE)
)
for (x in sample_module_groups) {
  x$identity_rank_within_sample <- rank(
    -x$pseudobulk_module_score,
    ties.method = "min",
    na.last = "keep"
  )
  myofib_score <- x$pseudobulk_module_score[x$reproduced_identity == "Myofib"]
  best_other_score <- safe_max(
    x$pseudobulk_module_score[x$reproduced_identity != "Myofib"]
  )
  x$myofib_minus_best_other <- if (length(myofib_score) == 1L) {
    myofib_score - best_other_score
  } else {
    NA_real_
  }
  rank_rows[[rank_index]] <- x
  rank_index <- rank_index + 1L
}
module_ranks <- do.call(rbind, rank_rows)
rownames(module_ranks) <- NULL
safe_write_csv(
  module_ranks,
  file.path(result_dir, "GSE173252_module_identity_ranks_by_sample_v1.csv")
)

dd_myofib_ranks <- module_ranks[
  module_ranks$source == "DD" & module_ranks$reproduced_identity == "Myofib",
  , drop = FALSE
]
module_specificity_summary <- do.call(rbind, lapply(
  split(dd_myofib_ranks, dd_myofib_ranks$module),
  function(x) {
    dominant_all <- all(x$identity_rank_within_sample == 1, na.rm = TRUE) &&
      all(x$myofib_minus_best_other > 0, na.rm = TRUE)
    data.frame(
      module = x$module[1],
      DD_samples = nrow(x),
      sample_ids = paste(x$sample, collapse = ";"),
      median_myofib_identity_rank = safe_median(x$identity_rank_within_sample),
      rank1_fraction = mean(x$identity_rank_within_sample == 1, na.rm = TRUE),
      median_myofib_minus_best_other = safe_median(x$myofib_minus_best_other),
      minimum_myofib_minus_best_other = safe_min(x$myofib_minus_best_other),
      module_specificity_class = if (dominant_all) {
        "Myofib_dominant_all_DD"
      } else if (mean(x$identity_rank_within_sample == 1, na.rm = TRUE) >= (2 / 3)) {
        "Myofib_dominant_two_of_three_DD"
      } else {
        "Not_Myofib_dominant"
      },
      stringsAsFactors = FALSE
    )
  }
))
module_specificity_summary <- module_specificity_summary[order(
  module_specificity_summary$median_myofib_identity_rank,
  -module_specificity_summary$median_myofib_minus_best_other
), , drop = FALSE]
rownames(module_specificity_summary) <- NULL
safe_write_csv(
  module_specificity_summary,
  file.path(result_dir, "GSE173252_DD_module_specificity_summary_v1.csv")
)

# -----------------------------------------------------------------------------
# 2. Pseudobulk counts and log2 CPM for gene-level pairwise specificity
# -----------------------------------------------------------------------------
cell_group <- interaction(
  singlet_metadata$orig.ident,
  singlet_metadata$seurat_clusters,
  drop = TRUE,
  sep = "__"
)
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
library_sizes <- Matrix::colSums(pseudobulk_counts)
group_parts <- do.call(rbind, strsplit(group_levels, "__", fixed = TRUE))
group_sample <- group_parts[, 1]
group_cluster <- group_parts[, 2]
group_identity <- unname(identity_map[group_cluster])
group_source <- singlet_metadata$source[
  match(group_sample, singlet_metadata$orig.ident)
]

candidate_genes <- unique(candidate_evidence$gene[
  candidate_evidence$candidate_tier != "C_no_specific_enrichment"
])
mechanosensor_genes <- unique(gene_sets$gene[
  gene_sets$module == "mechanosensor_channels"
])
positive_control_genes <- c("PDPN", "SCX", "FAP", "POSTN")
fb_control_genes <- c("PDGFRA", "APOD")
audit_genes <- unique(c(
  candidate_genes,
  mechanosensor_genes,
  positive_control_genes,
  fb_control_genes
))
audit_genes <- audit_genes[audit_genes %in% rownames(pseudobulk_counts)]

audit_log_cpm <- log2(
  sweep(
    as.matrix(pseudobulk_counts[audit_genes, , drop = FALSE]),
    2, library_sizes, "/"
  ) * 1e6 + 0.5
)

dd_samples <- unique(group_sample[group_source == "DD"])
comparators <- c("FB", "VSMC", "Pericyte")
pairwise_rows <- list()
pairwise_index <- 1L
for (gene in audit_genes) {
  for (sample_id in dd_samples) {
    myofib_index <- which(group_sample == sample_id & group_identity == "Myofib")
    if (length(myofib_index) != 1L) next
    for (comparator in comparators) {
      comparator_index <- which(
        group_sample == sample_id & group_identity == comparator
      )
      if (length(comparator_index) != 1L) next
      pairwise_rows[[pairwise_index]] <- data.frame(
        gene = gene,
        sample = sample_id,
        comparator_identity = comparator,
        myofib_log2CPM = audit_log_cpm[gene, myofib_index],
        comparator_log2CPM = audit_log_cpm[gene, comparator_index],
        myofib_minus_comparator_log2CPM =
          audit_log_cpm[gene, myofib_index] - audit_log_cpm[gene, comparator_index],
        stringsAsFactors = FALSE
      )
      pairwise_index <- pairwise_index + 1L
    }
  }
}
pairwise_evidence <- do.call(rbind, pairwise_rows)
safe_write_csv(
  pairwise_evidence,
  file.path(result_dir, "GSE173252_DD_candidate_pairwise_log2CPM_v1.csv")
)

# -----------------------------------------------------------------------------
# 3. Candidate specificity summary
# -----------------------------------------------------------------------------
candidate_annotations <- unique(candidate_evidence[, c(
  "module", "gene", "candidate_tier",
  "myofib_percent_cells_expressing",
  "myofib_mean_log_normalized_expression"
)])

gene_summary_rows <- list()
summary_index <- 1L
for (gene in audit_genes) {
  gene_data <- pairwise_evidence[pairwise_evidence$gene == gene, , drop = FALSE]
  get_difference <- function(comparator) {
    gene_data$myofib_minus_comparator_log2CPM[
      gene_data$comparator_identity == comparator
    ]
  }
  diff_fb <- get_difference("FB")
  diff_vsmc <- get_difference("VSMC")
  diff_pericyte <- get_difference("Pericyte")
  all_differences <- c(diff_fb, diff_vsmc, diff_pericyte)
  annotation_row <- candidate_annotations[
    candidate_annotations$gene == gene, , drop = FALSE
  ]
  module_name <- if (nrow(annotation_row) > 0L) {
    annotation_row$module[1]
  } else if (gene %in% mechanosensor_genes) {
    "mechanosensor_channels"
  } else if (gene %in% positive_control_genes) {
    "myofib_positive_control"
  } else {
    "FB_control"
  }
  original_tier <- if (nrow(annotation_row) > 0L) {
    annotation_row$candidate_tier[1]
  } else {
    "not_in_step06_shortlist"
  }
  myofib_prevalence <- if (nrow(annotation_row) > 0L) {
    annotation_row$myofib_percent_cells_expressing[1]
  } else {
    NA_real_
  }
  myofib_mean_expression <- if (nrow(annotation_row) > 0L) {
    annotation_row$myofib_mean_log_normalized_expression[1]
  } else {
    NA_real_
  }
  robust_all <- length(all_differences) == length(dd_samples) * length(comparators) &&
    all(all_differences > 0)
  robust_fb <- length(diff_fb) == length(dd_samples) && all(diff_fb > 0)
  specificity_class <- if (robust_all) {
    "robust_Myofib_over_all_three_identities"
  } else if (robust_fb) {
    "robust_over_FB_but_not_all_identities"
  } else {
    "donor_or_identity_dependent"
  }
  gene_summary_rows[[summary_index]] <- data.frame(
    module = module_name,
    gene = gene,
    step06_candidate_tier = original_tier,
    myofib_percent_cells_expressing = myofib_prevalence,
    myofib_mean_log_normalized_expression = myofib_mean_expression,
    DD_samples = length(dd_samples),
    median_Myofib_minus_FB = safe_median(diff_fb),
    minimum_Myofib_minus_FB = safe_min(diff_fb),
    positive_fraction_vs_FB = safe_positive_fraction(diff_fb),
    median_Myofib_minus_VSMC = safe_median(diff_vsmc),
    minimum_Myofib_minus_VSMC = safe_min(diff_vsmc),
    positive_fraction_vs_VSMC = safe_positive_fraction(diff_vsmc),
    median_Myofib_minus_Pericyte = safe_median(diff_pericyte),
    minimum_Myofib_minus_Pericyte = safe_min(diff_pericyte),
    positive_fraction_vs_Pericyte = safe_positive_fraction(diff_pericyte),
    minimum_difference_over_all_donors_and_identities = safe_min(all_differences),
    specificity_class = specificity_class,
    stringsAsFactors = FALSE
  )
  summary_index <- summary_index + 1L
}
gene_specificity_summary <- do.call(rbind, gene_summary_rows)
class_order <- c(
  "robust_Myofib_over_all_three_identities",
  "robust_over_FB_but_not_all_identities",
  "donor_or_identity_dependent"
)
gene_specificity_summary$class_rank <- match(
  gene_specificity_summary$specificity_class,
  class_order
)
gene_specificity_summary <- gene_specificity_summary[order(
  gene_specificity_summary$class_rank,
  -gene_specificity_summary$minimum_difference_over_all_donors_and_identities,
  -gene_specificity_summary$median_Myofib_minus_FB
), , drop = FALSE]
gene_specificity_summary$class_rank <- NULL
rownames(gene_specificity_summary) <- NULL
safe_write_csv(
  gene_specificity_summary,
  file.path(result_dir, "GSE173252_DD_candidate_specificity_summary_v1.csv")
)

mechanosensor_focus <- gene_specificity_summary[
  gene_specificity_summary$gene %in% mechanosensor_genes,
  , drop = FALSE
]
safe_write_csv(
  mechanosensor_focus,
  file.path(result_dir, "GSE173252_DD_mechanosensor_focus_summary_v1.csv")
)

# -----------------------------------------------------------------------------
# 4. Cell prevalence by sample and identity
# -----------------------------------------------------------------------------
prevalence_rows <- list()
prevalence_index <- 1L
for (group_name in levels(group_factor)) {
  cell_index <- which(group_factor == group_name)
  first_cell <- cell_index[1]
  sample_id <- singlet_metadata$orig.ident[first_cell]
  identity_name <- singlet_metadata$reproduced_identity[first_cell]
  group_cluster_id <- singlet_metadata$seurat_clusters[first_cell]
  percentages <- Matrix::rowMeans(
    rna_counts[audit_genes, singlet_keep, drop = FALSE][, cell_index, drop = FALSE] > 0
  ) * 100
  means <- Matrix::rowMeans(
    rna_logdata[audit_genes, singlet_keep, drop = FALSE][, cell_index, drop = FALSE]
  )
  prevalence_rows[[prevalence_index]] <- data.frame(
    sample = sample_id,
    source = singlet_metadata$source[first_cell],
    seurat_cluster = group_cluster_id,
    reproduced_identity = identity_name,
    singlet_cells = length(cell_index),
    gene = audit_genes,
    percent_cells_expressing = as.numeric(percentages),
    mean_log_normalized_expression = as.numeric(means),
    stringsAsFactors = FALSE
  )
  prevalence_index <- prevalence_index + 1L
}
prevalence_table <- do.call(rbind, prevalence_rows)
safe_write_csv(
  prevalence_table,
  file.path(result_dir, "GSE173252_candidate_prevalence_by_sample_identity_v1.csv")
)

# -----------------------------------------------------------------------------
# 5. Leave-one-DD-donor-out sign stability
# -----------------------------------------------------------------------------
loo_rows <- list()
loo_index <- 1L
for (gene in audit_genes) {
  gene_data <- pairwise_evidence[pairwise_evidence$gene == gene, , drop = FALSE]
  for (omitted_sample in dd_samples) {
    retained <- gene_data[gene_data$sample != omitted_sample, , drop = FALSE]
    for (comparator in comparators) {
      differences <- retained$myofib_minus_comparator_log2CPM[
        retained$comparator_identity == comparator
      ]
      loo_rows[[loo_index]] <- data.frame(
        gene = gene,
        omitted_DD_sample = omitted_sample,
        comparator_identity = comparator,
        retained_samples = paste(
          retained$sample[retained$comparator_identity == comparator],
          collapse = ";"
        ),
        median_difference_after_omission = safe_median(differences),
        minimum_difference_after_omission = safe_min(differences),
        all_retained_differences_positive = all(differences > 0),
        stringsAsFactors = FALSE
      )
      loo_index <- loo_index + 1L
    }
  }
}
loo_stability <- do.call(rbind, loo_rows)
safe_write_csv(
  loo_stability,
  file.path(result_dir, "GSE173252_DD_leave_one_donor_out_stability_v1.csv")
)

loo_gene_summary <- do.call(rbind, lapply(
  split(loo_stability, loo_stability$gene),
  function(x) data.frame(
    gene = x$gene[1],
    leave_one_out_comparisons = nrow(x),
    all_leave_one_out_medians_positive = all(x$median_difference_after_omission > 0),
    all_leave_one_out_pairs_positive = all(x$all_retained_differences_positive),
    minimum_leave_one_out_median = safe_min(x$median_difference_after_omission),
    minimum_retained_pair_difference = safe_min(x$minimum_difference_after_omission),
    stringsAsFactors = FALSE
  )
))
loo_gene_summary <- loo_gene_summary[order(
  -loo_gene_summary$all_leave_one_out_pairs_positive,
  -loo_gene_summary$minimum_retained_pair_difference
), , drop = FALSE]
rownames(loo_gene_summary) <- NULL
safe_write_csv(
  loo_gene_summary,
  file.path(result_dir, "GSE173252_DD_leave_one_donor_out_gene_summary_v1.csv")
)

# -----------------------------------------------------------------------------
# 6. Structured audit report
# -----------------------------------------------------------------------------
robust_candidates <- gene_specificity_summary[
  gene_specificity_summary$specificity_class ==
    "robust_Myofib_over_all_three_identities" &
    gene_specificity_summary$step06_candidate_tier %in% c(
      "A_consistent_DD_myofib_enrichment", "B_partial_support"
    ),
  , drop = FALSE
]
robust_candidates <- merge(
  robust_candidates,
  loo_gene_summary,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
robust_candidates <- robust_candidates[order(
  -robust_candidates$all_leave_one_out_pairs_positive,
  -robust_candidates$minimum_difference_over_all_donors_and_identities
), , drop = FALSE]

candidate_lines <- if (nrow(robust_candidates) == 0L) {
  "- No step-06 candidate remained higher than all three comparator identities in every DD donor."
} else {
  apply(head(robust_candidates, 20), 1, function(x) paste0(
    "- ", x[["gene"]], " [", x[["module"]], "]: minimum pairwise log2CPM difference = ",
    format(
      round(as.numeric(x[["minimum_difference_over_all_donors_and_identities"]]), 3),
      nsmall = 3
    ),
    "; leave-one-donor-out all-pairs positive = ",
    x[["all_leave_one_out_pairs_positive"]]
  ))
}

module_lines <- apply(module_specificity_summary, 1, function(x) paste0(
  "- ", x[["module"]], ": ", x[["module_specificity_class"]],
  "; median Myofib rank = ", x[["median_myofib_identity_rank"]],
  "; minimum gap from best other identity = ",
  format(round(as.numeric(x[["minimum_myofib_minus_best_other"]]), 3), nsmall = 3)
))

focus_genes <- c("PIEZO1", "PIEZO2", "TRPV4", "TMEM63B", "PKD2")
focus_table <- mechanosensor_focus[mechanosensor_focus$gene %in% focus_genes, , drop = FALSE]
focus_lines <- if (nrow(focus_table) == 0L) {
  "- None of the pre-specified focus mechanosensors was detected."
} else {
  apply(focus_table, 1, function(x) paste0(
    "- ", x[["gene"]], ": ", x[["specificity_class"]],
    "; expressing Myofib cells = ",
    format(round(as.numeric(x[["myofib_percent_cells_expressing"]]), 1), nsmall = 1), "%",
    "; median Myofib-FB difference = ",
    format(round(as.numeric(x[["median_Myofib_minus_FB"]]), 3), nsmall = 3),
    "; minimum difference over all identities/donors = ",
    format(
      round(as.numeric(x[["minimum_difference_over_all_donors_and_identities"]]), 3),
      nsmall = 3
    )
  ))
}

report_lines <- c(
  "## Material Passport",
  "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Origin Mode: validate",
  paste0("- Origin Date: ", format(Sys.Date(), "%Y-%m-%d")),
  "- Verification Status: ANALYZED",
  "- Version Label: mechanosensitivity_specificity_audit_v1",
  "",
  "## GSE173252 mechanosensitivity specificity and donor-robustness audit",
  "",
  "### Audit design",
  "",
  "- Biological unit: DD donor/sample (DD1-DD3).",
  "- Primary comparison: Myofib versus FB, VSMC and Pericyte separately within each DD sample.",
  "- Robust specificity requires a positive Myofib-comparator log2CPM difference for all three identities in all three DD samples.",
  "- Leave-one-donor-out analysis repeats direction checks after omitting DD1, DD2 or DD3.",
  "- No cell-level p values are calculated.",
  "",
  "### Module-level identity specificity",
  "",
  module_lines,
  "",
  "### Step-06 candidates passing all-identity specificity",
  "",
  candidate_lines,
  "",
  "### Focus mechanosensors",
  "",
  focus_lines,
  "",
  "### Mandatory interpretation limits",
  "",
  "- Myofib specificity does not establish disease specificity because controls do not contain a comparable Myofib population.",
  "- Transcriptional enrichment does not demonstrate force sensing, channel activity, protein abundance or causal function.",
  "- Contractile genes such as ACTA2, TAGLN and CNN1 may describe myofibroblast differentiation rather than a distinct mechanosensitive state.",
  "- PIEZO2 and other low-to-moderate prevalence genes require donor-level prevalence checks and independent dataset validation.",
  "- With only three DD donors, leave-one-out analysis is a fragility diagnostic rather than formal statistical validation."
)
writeLines(
  report_lines,
  file.path(result_dir, "GSE173252_mechanosensitivity_specificity_audit_report_v1.md"),
  useBytes = TRUE
)

message("Specificity and donor-robustness audit completed.")
message("Results: ", result_dir)
