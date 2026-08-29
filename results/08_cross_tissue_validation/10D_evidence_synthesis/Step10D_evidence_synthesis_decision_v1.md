# Step 10D evidence synthesis decision

## Current conclusion

The final interim evidence grade is **CAUTION**.
The appropriate interpretation class is **ACTOMYOSIN_ROBUST_INTEGRIN_NON_SPECIFIC**.

The robust part of the result is the actomyosin/Rho-associated program: the module is positive in 12/12 samples and the leave-one-sample-out gate passes. The integrin/focal-adhesion module is also positive in 12/12 samples and 5/6 frozen candidates are supported, but its specificity is not established after the pre-specified competitor audit.

## Failed or incomplete boundaries

- Step 10A source independence: PARTIAL; overlapping accessions = GSE173252, PRJNA607098.
- F7/F8 anatomical-label equivalence: unresolved for direct replication claims.
- PIEZO2: not supported; channel-specific mechanosensor evidence remains negative/weak.
- Competition robustness: FALSE.
- integrin_focal_adhesion versus tgf_fibrosis: residual positive samples = 4/12; residual median = -0.005048.
- integrin_focal_adhesion versus inflammation_ap1_nfkb: residual positive samples = 3/12; residual median = -0.009269.
- integrin_focal_adhesion versus hypoxia: residual positive samples = 5/12; residual median = -0.005488.
- integrin_focal_adhesion versus cell_cycle: residual positive samples = 6/12; residual median = -0.004852.

## Recommended claim boundary

- Supported wording: fascia-like myofibroblast states show a reproducible actomyosin/Rho-associated transcriptional program, accompanied by integrin/focal-adhesion enrichment.
- Required limitation: the integrin component is not separable from TGF/fibrosis, inflammatory, hypoxia, and cell-cycle state signals in this validation source.
- Prohibited wording: independent replication, direct F8 replication, PIEZO2-mediated mechanosensing, or causal mechanotransduction.

## Next decision

Do not rerun Step 10C with relaxed thresholds. The highest-value next step is a genuinely non-overlapping donor-level source with sample metadata. If that source cannot be recovered, proceed to manuscript-level synthesis as a bounded, hypothesis-generating reanalysis.

## Material Passport

- Inputs: corrected Step 10A provenance audit and Step 10C frozen-contract outputs.
- Transformation: metadata-aware evidence synthesis only; no new expression values or hypothesis tests.
- Statistical units: the 12 PRJNA607098 sample IDs inherited from Step 10C.
- Reproducibility: all inputs and failed-pair details are exported in the same result directory.
- Integrity status: interim decision; independent replication remains unavailable.
