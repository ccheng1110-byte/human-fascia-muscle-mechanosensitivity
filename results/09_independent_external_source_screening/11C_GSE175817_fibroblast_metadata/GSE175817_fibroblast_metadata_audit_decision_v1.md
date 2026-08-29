# Step 11C GSE175817 fibroblast metadata audit

## Input

- Metadata URL: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175817/suppl/GSE175817_meta_fibroblast.csv.gz
- Local metadata: ./data/metadata/independent_sources/GSE175817/GSE175817_meta_fibroblast.csv.gz
- Cells/rows: 29390; fields: 6.
- No expression matrix or raw sequencing archive was downloaded.

## Resolved fields

- Candidate sample/donor field: NA.
- Candidate cell-type/cluster field: NA.
- Fibroblast-like labelled cells, if a cell-type field was found: NA.

## Gate

This is a source-structure audit only. Independent validation is allowed only if the object-level expression source, cell-state field, and sample/donor mapping are all confirmed in the next step.

## Material Passport

- Source: GSE175817, non-overlapping candidate from corrected Step 10A.
- Transformation: downloaded fibroblast-specific metadata were summarized by field, sample/donor, and cell-type coverage.
- No gene-level hypothesis test was performed.
- Reproducibility: raw compressed metadata and the small filelist are retained under the project metadata directory.
