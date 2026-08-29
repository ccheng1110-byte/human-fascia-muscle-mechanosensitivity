options(stringsAsFactors = FALSE)

# Step 08C1: sample-ID and paired-state feasibility audit.
# Run from RStudio. This reads only remote observation metadata and performs no
# expression download and no hypothesis test.

project_dir <- "."
python_runtime_root <- file.path(project_dir, ".runtime", "python_step08b")
r_user_cache_dir <- file.path(python_runtime_root, "r_user_cache")
uv_cache_dir <- file.path(python_runtime_root, "uv_cache")
uv_python_install_dir <- file.path(python_runtime_root, "uv_python")

runtime_dirs <- c(
  python_runtime_root,
  r_user_cache_dir,
  uv_cache_dir,
  uv_python_install_dir
)
for (runtime_dir in runtime_dirs) {
  dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
}

Sys.setenv(
  R_USER_CACHE_DIR = r_user_cache_dir,
  XDG_CACHE_HOME = r_user_cache_dir,
  UV_CACHE_DIR = uv_cache_dir,
  UV_PYTHON_INSTALL_DIR = uv_python_install_dir,
  RETICULATE_PYTHON = "managed"
)

reticulate_cache_dir <- tools::R_user_dir("reticulate", "cache")
if (!startsWith(
  tolower(normalizePath(reticulate_cache_dir, winslash = "/", mustWork = FALSE)),
  tolower(normalizePath(python_runtime_root, winslash = "/", mustWork = TRUE))
)) {
  stop(
    "reticulate cache redirection failed. Restart R and source this script ",
    "before loading reticulate or Seurat."
  )
}

python_helper <- file.path(
  project_dir, "scripts",
  "08C1_audit_PRJNA607098_sample_metadata_remote_zarr.py"
)
result_dir <- file.path(
  project_dir, "results", "06_external_validation",
  "skin_fibroblast_atlas_2025",
  "08C1_sample_metadata_audit_PRJNA607098"
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(python_helper)) {
  stop("Missing Python helper: ", python_helper)
}
if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("Install reticulate once with: install.packages('reticulate')")
}
if (utils::packageVersion("reticulate") < "1.41.0") {
  stop("reticulate >= 1.41.0 is required.")
}

reticulate::py_require(
  packages = c(
    "zarr>=2.18,<3",
    "fsspec>=2024.6",
    "aiohttp>=3.9",
    "numpy>=1.26,<3"
  ),
  python_version = ">=3.10,<3.13",
  exclude_newer = "2026-08-23"
)

python_environment <- new.env(parent = emptyenv())
reticulate::source_python(
  python_helper,
  envir = python_environment,
  convert = TRUE
)

message("Auditing PRJNA607098 sample IDs and F7/non-F7 coverage.")
message("This step reads metadata only and performs no hypothesis test.")

audit_result <- python_environment$run_step08c1(
  output_dir = result_dir,
  source_id = "PRJNA607098",
  target_state = "F7: Fascia-like myofibroblast",
  minimum_cells_per_state = 20L,
  minimum_eligible_samples = 6L,
  expected_public_samples = 12L
)

message("Step 08C1 sample-metadata audit completed.")
message("Resolved sample field: ", audit_result$sample_field)
message("Sample-ID provenance: ", audit_result$sample_provenance)
message("Sample-index format valid: ", audit_result$sample_format_valid)
message("PRJNA607098 cells: ", audit_result$total_source_cells)
message("Distinct sample IDs: ", audit_result$n_samples)
message("Exact expected sample count: ", audit_result$exact_expected_sample_count)
message("Eligible paired sample-state units: ", audit_result$n_eligible_samples)
message("Missing sample-ID fraction: ", audit_result$missing_sample_fraction)
message("Largest sample fraction: ", audit_result$largest_sample_fraction)
message("Proceed to Step 08C2: ", audit_result$proceed_to_08c2)
message("Sample eligibility: ", audit_result$sample_eligibility)
message("Cell-type counts: ", audit_result$sample_celltype_counts)
message("Decision: ", audit_result$decision)
