# Step 11D GSE130973 Seurat object audit

## Scope

- Source: GSE130973, non-overlapping candidate after Step 10A.
- This step downloads one processed Seurat RDS. The GSE130973 supplementary directory does not provide a filelist.txt; the file inventory was already recorded in Step 11B.
- No differential expression or mechanosensitivity test is performed.

## Gate

- Cell metadata rows: 15457.
- Sample/donor candidate fields: 1.
- Cell-type candidate fields: 2.
- Frozen candidate genes found: 15/15.

The object can proceed to expression validation only if a cell-state field and a sample/donor field are both confirmed. Any donor-level test must use the object-level donor field rather than the five GEO series subjects alone.

## Material Passport

- Input: GSE130973_seurat_analysis_lyko.rds.gz from the official GEO supplementary directory.
- Transformation: object structure, metadata fields, and frozen-gene presence were audited without testing.
- Reproducibility: the downloaded RDS is retained under the project raw-data directory and all audit tables are exported.
