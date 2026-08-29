# Download processed data and sample metadata for the human fascia datasets.
# Project: human_fascia_muscle_mechanosensitivity
# Scope: GSE273293 and GSE173252
# FASTQ/SRA/BAM/CRAM files are explicitly excluded.

options(stringsAsFactors = FALSE)
options(timeout = 3600)

project_dir <- "."

dir_list <- c(
  file.path(project_dir, "data", "raw", "GSE273293"),
  file.path(project_dir, "data", "raw", "GSE173252"),
  file.path(project_dir, "data", "processed", "GSE273293"),
  file.path(project_dir, "data", "processed", "GSE173252"),
  file.path(project_dir, "data", "metadata", "GSE273293"),
  file.path(project_dir, "data", "metadata", "GSE173252"),
  file.path(project_dir, "scripts"),
  file.path(project_dir, "logs")
)

invisible(lapply(dir_list, dir.create, recursive = TRUE, showWarnings = FALSE))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  BiocManager::install("SummarizedExperiment", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(SummarizedExperiment)
})

download_one <- function(url, destfile, retries = 3L) {
  if (file.exists(destfile) && isTRUE(file.info(destfile)$size > 0)) {
    message("Using existing file: ", destfile)
    return(invisible(destfile))
  }

  for (attempt in seq_len(retries)) {
    message(sprintf("Downloading [%d/%d]: %s", attempt, retries, url))

    ok <- tryCatch({
      status <- utils::download.file(
        url = url,
        destfile = destfile,
        mode = "wb",
        method = "libcurl",
        quiet = FALSE
      )
      identical(status, 0L) &&
        file.exists(destfile) &&
        isTRUE(file.info(destfile)$size > 0)
    }, error = function(e) {
      message("Download error: ", conditionMessage(e))
      FALSE
    })

    if (ok) {
      return(invisible(destfile))
    }

    if (file.exists(destfile)) {
      unlink(destfile)
    }

    if (attempt < retries) {
      Sys.sleep(2 ^ attempt)
    }
  }

  stop("Download failed after retries: ", url)
}

extract_series_metadata <- function(gse_id, metadata_dir, processed_dir) {
  # Keep a compact GEO family record for reproducible sample-level auditing.
  soft_path <- tryCatch(
    GEOquery::getGEOfile(
      GEO = gse_id,
      destdir = metadata_dir,
      amount = "brief"
    ),
    error = function(e) {
      message("Could not download GEO brief record for ", gse_id, ": ",
              conditionMessage(e))
      NULL
    }
  )

  # If a Series Matrix exists, GEOquery parses it into SummarizedExperiment
  # objects and carries sample metadata in colData(). Some scRNA-seq series do
  # not provide a usable Series Matrix; that case is allowed and logged.
  gse_matrix <- tryCatch(
    GEOquery::getGEO(
      GEO = gse_id,
      GSEMatrix = TRUE,
      destdir = processed_dir,
      getGPL = FALSE,
      parseCharacteristics = TRUE
    ),
    error = function(e) {
      message("No usable Series Matrix parsed for ", gse_id, ": ",
              conditionMessage(e))
      NULL
    }
  )

  if (!is.null(gse_matrix)) {
    if (!is.list(gse_matrix)) {
      gse_matrix <- list(gse_matrix)
    }

    for (i in seq_along(gse_matrix)) {
      se <- gse_matrix[[i]]
      object_path <- file.path(
        processed_dir,
        sprintf("%s_series_matrix_%02d.rds", gse_id, i)
      )
      saveRDS(se, object_path)

      # GEOquery may return either a SummarizedExperiment or an older
      # ExpressionSet, depending on the installed GEOquery configuration.
      if (inherits(se, "SummarizedExperiment")) {
        sample_metadata <- as.data.frame(
          SummarizedExperiment::colData(se),
          optional = TRUE
        )
      } else if (inherits(se, "ExpressionSet")) {
        sample_metadata <- as.data.frame(
          Biobase::pData(se),
          optional = TRUE
        )
      } else {
        stop(
          "Unsupported GEO object class: ",
          paste(class(se), collapse = ", ")
        )
      }
      sample_metadata$sample_id <- rownames(sample_metadata)
      sample_metadata <- sample_metadata[, c(
        "sample_id",
        setdiff(names(sample_metadata), "sample_id")
      ), drop = FALSE]

      utils::write.csv(
        sample_metadata,
        file = file.path(
          metadata_dir,
          sprintf("%s_series_matrix_sample_metadata_%02d.csv", gse_id, i)
        ),
        row.names = FALSE,
        na = ""
      )
    }
  }

  invisible(soft_path)
}

download_processed_supplementary <- function(gse_id, raw_dir, metadata_dir) {
  # First list files without downloading. GEOquery returns fname and url.
  listing <- GEOquery::getGEOSuppFiles(
    GEO = gse_id,
    makeDirectory = FALSE,
    baseDir = raw_dir,
    fetch_files = FALSE
  )

  if (is.null(listing) || nrow(listing) == 0L) {
    message("No supplementary files found for ", gse_id)
    return(invisible(NULL))
  }

  listing$gse_id <- gse_id
  utils::write.csv(
    listing,
    file = file.path(metadata_dir, sprintf("%s_supplementary_file_inventory.csv", gse_id)),
    row.names = FALSE,
    na = ""
  )

  # Processed single-cell formats and common processed tables.
  processed_pattern <- paste0(
    "(rds|rda|rdata|h5(ad)?|mtx|tsv|csv|txt)(\\.gz)?$"
  )
  is_processed <- grepl(
    processed_pattern,
    listing$fname,
    ignore.case = TRUE,
    perl = TRUE
  )

  # GEO often labels a tar archive as *_RAW.tar even when its contents are
  # processed 10X MTX/TSV files. Allow that archive only after checking GEO's
  # filelist.txt for processed single-cell file names.
  archive_has_processed_sc <- FALSE
  archive_listing <- tryCatch(
    GEOquery::getGEOSeriesFileListing(gse_id),
    error = function(e) NULL
  )
  if (!is.null(archive_listing) && nrow(archive_listing) > 0L) {
    archive_text <- paste(unlist(archive_listing), collapse = " ")
    archive_has_processed_sc <- grepl(
      "(barcodes|features|matrix|mtx|h5ad|rds|tsv)",
      archive_text,
      ignore.case = TRUE,
      perl = TRUE
    )
  }

  is_processed_archive <- grepl(
    "_RAW\\.tar$",
    listing$fname,
    ignore.case = TRUE,
    perl = TRUE
  ) & archive_has_processed_sc

  is_forbidden_raw <- grepl(
    "\\.(fastq|fq|sra|bam|cram|bai|crai)(\\.gz)?$",
    listing$fname,
    ignore.case = TRUE,
    perl = TRUE
  )

  selected <- listing[(is_processed | is_processed_archive) & !is_forbidden_raw, , drop = FALSE]

  if (nrow(selected) == 0L) {
    message(
      "No processed supplementary file was auto-selected for ", gse_id,
      ". Review the inventory CSV before adding a manual allow-list."
    )
    return(invisible(NULL))
  }

  manifest <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    destination <- file.path(raw_dir, selected$fname[[i]])
    status <- "downloaded_or_cached"

    tryCatch(
      download_one(selected$url[[i]], destination),
      error = function(e) {
        status <<- paste("ERROR:", conditionMessage(e))
      }
    )

    manifest[[i]] <- data.frame(
      gse_id = gse_id,
      filename = selected$fname[[i]],
      url = selected$url[[i]],
      destination = destination,
      status = status,
      size_bytes = if (file.exists(destination)) file.info(destination)$size else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  manifest <- do.call(rbind, manifest)
  utils::write.csv(
    manifest,
    file = file.path(metadata_dir, sprintf("%s_download_manifest.csv", gse_id)),
    row.names = FALSE,
    na = ""
  )

  invisible(manifest)
}

dataset_ids <- c("GSE273293", "GSE173252")

for (gse_id in dataset_ids) {
  message("\n========== ", gse_id, " ==========")

  raw_dir <- file.path(project_dir, "data", "raw", gse_id)
  processed_dir <- file.path(project_dir, "data", "processed", gse_id)
  metadata_dir <- file.path(project_dir, "data", "metadata", gse_id)

  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

  extract_series_metadata(
    gse_id = gse_id,
    metadata_dir = metadata_dir,
    processed_dir = processed_dir
  )

  download_processed_supplementary(
    gse_id = gse_id,
    raw_dir = raw_dir,
    metadata_dir = metadata_dir
  )
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(project_dir, "logs", "01_download_geo_sessionInfo.txt")
)

message("\nDownload stage completed. Review each *_supplementary_file_inventory.csv and *_download_manifest.csv before analysis.")
