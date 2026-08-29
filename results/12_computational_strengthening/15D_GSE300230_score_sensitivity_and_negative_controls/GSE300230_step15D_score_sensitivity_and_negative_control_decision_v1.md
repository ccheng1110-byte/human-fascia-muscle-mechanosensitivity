# Step 15D GSE300230 score sensitivity and negative-control audit

- Core module score direction stable across primary z-mean scoring: **TRUE**.
- Overall evidence grade remains: **CAUTION**.
- Cell-intrinsic versus composition-associated explanation: **NOT ESTIMABLE** from this matrix.

## Interpretation boundary

This step tests robustness of sample-level module scoring and expression-matched negative controls. It does not create cell-level composition information, does not fit a post-hoc multivariable adjustment for cell cycle, and does not upgrade the evidence grade. The cell-cycle competitor remains a required alternative explanation.

## Material Passport

- Input: official GSE300230 16-sample processed raw-count matrix.
- Frozen inputs: Step 15B sample map, module registry, and candidate panel.
- Sensitivity grid: z-mean, z-median, rank-mean, and centered-mean module scores.
- Negative controls: 500 expression-bin-matched random sets per module where feasible, seed 1504.
- Leave-one-gene-out: applied to every module with at least two retained genes.
