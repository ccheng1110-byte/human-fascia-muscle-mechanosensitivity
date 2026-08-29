## Material Passport

- Origin Skill: nature-downloader + documents
- Origin Mode: lawful OA download / structural and visual document audit
- Origin Date: 2026-08-24
- Verification Status: ANALYZED
- Version Label: GSE273293_supporting_information_inspection_v1

## Scope

Inspect both official supporting-information DOCX files for any public route to
sample-resolved processed expression data for GSE273293 / PRJNA1206333.
Targets were: per-sample matrices, a cell-barcode-to-sample map, a Seurat/RDS or
AnnData/H5AD object, preprocessing code, or a study-specific repository URL.

## Download and integrity

Both files were downloaded from the official Springer static-content host.

- `12967_2024_5889_MOESM1_ESM.docx`: 8,889,938 bytes; SHA-256
  `21B1E844411D8088762AF8FE5745AE0FABE925099713E5D0CA95D130A87298F1`.
- `12967_2024_5889_MOESM2_ESM.docx`: 189,850 bytes; SHA-256
  `7FC7C1F1C3787A56F6A65E76305A2E58F208D020E797362704E49CB69401BFBA`.

Both have a valid DOCX/ZIP signature and could be parsed and rendered. The
machine-readable manifest is stored beside the downloaded files.

## Structural findings

### Supplementary Material 1

- Contains 24 supplementary figures as TIFF media and 74 Word paragraphs.
- Renders to 27 pages.
- Figure S3 (rendered page 4) explicitly shows all 14 labels: `Con_01` to
  `Con_04` and `GMC_01` to `GMC_10`, confirming that the authors retained
  sample identity in their internal integrated object.
- No external hyperlink relationship, embedded spreadsheet, RDS, H5AD, CSV,
  TSV, or repository URL was present.

### Supplementary Material 2

- Contains 11 supplementary tables represented as 17 EMF vector objects and
  53 Word paragraphs; it renders to 17 pages.
- Table S2 (page 2) gives clinical metadata for 10 GMC and 4 control samples.
- Table S3 (page 3) gives per-sample cell counts and cluster proportions. The
  total is 86,159 cells, matching the article.
- Table S5 (page 7) gives per-sample major-cell-type proportions.
- Table S8 (page 13) gives per-sample macrophage-subcluster counts and ratios.
- Table S10 (page 16) gives per-sample fibroblast-subcluster counts and ratios.
- Tables S4, S6, S7 and S9 contain cluster/cell-type differential-expression
  summaries, but not cell-level observations or donor-labelled expression
  matrices.
- No barcode mapping, per-sample matrix, processed object, code link, GitHub
  URL, or other external relationship was present.

## Recovery decision

The supporting information recovers **sample-level composition metadata**, but
not **sample-resolved gene expression**. Counts or proportions for a named
sample cannot be used to reconstruct which rows/columns of the pooled GEO
matrix belong to that sample.

Accordingly:

- sample-design and abundance audit: **RECOVERED**;
- 14-sample processed expression object: **NOT RECOVERED**;
- barcode-to-sample mapping: **NOT RECOVERED**;
- public study code/repository: **NOT RECOVERED**.

## Claim boundary

The supplement supports statements about the existence of 14 biological
samples and their reported cellular composition. It does not authorize
pseudobulk gene-expression testing, donor-level mechanosensitivity scoring, or
cell-level inference that treats cells as biological replicates.

## Recommended next action

Use the prepared corresponding-author request to obtain one of: a
sample-resolved Seurat/AnnData object, per-sample 10X matrices, or a
`cell_barcode -> sample_id` table. Do not begin the approximately 617.25 GiB
raw-data recovery unless a separate compute/storage plan is explicitly
approved.
