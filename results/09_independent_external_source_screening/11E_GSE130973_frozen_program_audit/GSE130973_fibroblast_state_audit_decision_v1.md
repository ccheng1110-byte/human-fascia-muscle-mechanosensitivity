# Step 11E GSE130973 exploratory frozen-program audit

## Design

- Statistical units: five subjects S1–S5.
- Candidate fibroblast-state clusters: 9, 1, 3, 2, 10.
- Comparator: all other clusters within the same subject.
- Candidate clusters were selected by the pre-defined top-five fibroblast-marker ranking from Step 11D2.

## Boundary

This is an exploratory state-level audit in an aging-skin dataset. It is not a direct F7/F8 replication, does not test disease status, and cannot upgrade the current evidence grade. Descriptive p values are not confirmatory with n = 5.

- Frozen registry genes present: 179/181.
- Missing registry genes: ANKRD1, CXCL8.

## Material Passport

- Source: GSE130973 processed Seurat object.
- Transformation: within-subject candidate-state versus other-cluster gene and module summaries.
- No causal model or disease-group test was fitted.
- Reproducibility: candidate cluster configuration and all per-subject contrasts are exported in this result directory.
