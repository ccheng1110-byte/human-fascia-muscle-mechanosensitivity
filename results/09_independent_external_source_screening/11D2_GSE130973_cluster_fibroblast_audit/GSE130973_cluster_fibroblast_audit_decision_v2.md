# Step 11D2 GSE130973 cluster-level fibroblast audit

## Corrected metadata interpretation

- `subj` is the valid sample/subject field and contains five subjects: S1–S5.
- `celltype.age` is a composite cluster-plus-age label, not a confirmed cell-type annotation.
- `integrated_snn_res.0.4` contains the underlying cluster IDs.

## Audit boundary

Clusters are ranked by fibroblast-marker expression and subject coverage. This is a descriptive state-identification step; it does not establish that any cluster is a validated fibroblast population and does not perform a gene-level test.

A cluster may proceed to targeted sample-level validation only after manual inspection confirms fibroblast identity, adequate representation across subjects, and a reproducible rule for including or excluding cells.

## Material Passport

- Source: GSE130973 processed Seurat object.
- Transformation: cluster-level marker means and subject coverage were calculated from the existing RNA data layer.
- No disease-group or mechanosensitivity hypothesis test was performed.
- Reproducibility: all cluster-by-subject counts, marker audit values, and gene availability are exported in this directory.
