# Step 17B S1 source and feasibility audit

## Decision

- Gate 17B: **PARTIAL_PROCEED_TO_TARGETED_EXPRESSION_AUDIT**.
- This step is a provenance and analyzability audit only; no hypothesis test was performed.
- Cells are not treated as independent biological replicates.

## PRJNA607098

- Reconstructed sample units: 12.
- Eligible paired sample units: 12.
- Sample metadata gate: TRUE.
- Frozen candidate coverage in the existing targeted extract: 15/15.
- Full registry genes available in the existing targeted extract: 18/181.
- Existing PRJNA gene inventory is targeted; genes not listed as present are not interpreted as biologically absent.
- A new cell-level targeted stream is required before competition-adjusted S1 inference.

## GSE130973

- Subjects: S1, S2, S3, S4, S5.
- Subject mapping gate: TRUE.
- Frozen candidate fibroblast-state clusters: 9, 1, 3, 2, 10.
- Candidate cluster/state gate: TRUE.
- Existing frozen registry coverage: 179/181.
- GSE130973 will remain a separate study-level supplementary source; it will not be pooled with PRJNA607098.

## Step 17C entry conditions

1. Stream the PRJNA607098 genes listed in `PRJNA607098_S1_registry_gene_retrieval_status_v1.csv` that are not available in the existing targeted extract.
2. Retain the exact F7/non-F7 sample reconstruction and the 12 sample units.
3. Perform cell-level co-detection descriptively, then use within-sample permutation and sample-level aggregation for inference.
4. Keep detection rate, library size, cell-cycle, ECM/TGF, inflammation and hypoxia as prespecified competition controls.
5. Do not label the result `cell-intrinsic` automatically; the output label must be a co-expression evidence level.

## Evidence boundary

- This audit does not establish donor independence for PRJNA607098.
- It does not prove cell-intrinsic mechanism, causal mechanosensitivity, or fascia specificity.
- The overall project evidence grade remains CAUTION.

## Material Passport

- Inputs: existing Step 08C1/08B3 PRJNA607098 audits, Step 11D2/11E GSE130973 audits, frozen candidate panel, and module registry.
- Transformation: local provenance, sample/state, and gene-coverage feasibility checks.
- New expression data downloaded: none.
- Next step: Step 17C targeted cell-level expression stream and co-expression-level analysis.
