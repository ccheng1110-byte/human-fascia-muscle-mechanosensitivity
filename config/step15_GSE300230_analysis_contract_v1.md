# Step 15B GSE300230 frozen analysis contract v1

**Freeze date:** 2026-08-26  
**Status:** frozen before condition-specific expression testing  
**Overall evidence entering Step 15C:** CAUTION

## 1. Data selection

- Primary matrix: `GSE300230_raw_counts_tensionTGFb.csv.gz` (16 samples).
- The separate eight-sample tension-only matrix is reserved for Step 15D sensitivity analysis and cannot replace the primary matrix after results are viewed.
- Sample identities and factors must be read from `config/step15_GSE300230_sample_map_v1.csv`.

## 2. Inferential boundary

- GM08401 and GM09503 are analyzed separately.
- Each cell line has a complete mechanical condition x TGF-beta condition design and two recorded batches.
- Batch is included as a blocking factor.
- Age background is descriptive only because one cell line represents each age. Age and cell-line effects are not separable.
- Cross-line summaries assess directional concordance only. Two lines are not treated as a random sample of human donors.

## 3. Frozen contrasts

### Primary

`mechanical_tension_vs_relaxed_tgfb_absent`

This isolates the tension contrast without exogenous TGF-beta and is evaluated separately in each cell line.

### Secondary

1. `mechanical_tension_vs_relaxed_tgfb_present`
2. `tgfb_present_vs_absent_relaxed`
3. `tgfb_present_vs_absent_tension`
4. `mechanical_x_tgfb`

No main-effect age contrast is authorized.

## 4. Frozen endpoints

- Primary modules: `integrin_focal_adhesion` and `actomyosin_rho`.
- Mechanism gates: `mechanosensor_channels` and `hippo_yap_taz`.
- Required competitors: `ecm_remodeling`, `tgf_fibrosis`, `inflammation_ap1_nfkb`, `hypoxia`, and `cell_cycle`.
- Frozen candidates: the unchanged 15 genes in `frozen_candidate_panel_v2.csv`.
- Exact endpoint registry: `config/step15_GSE300230_endpoint_registry_v1.csv`.

## 5. Gene nomenclature and filtering

- `CCN2` may be mapped to the historical registry symbol `CTGF` only when `CTGF` is absent.
- `CCN1` may be mapped to `CYR61` only when `CYR61` is absent.
- These are nomenclature aliases, not biological substitutions.
- Genes require CPM >= 1 in at least two of the 16 primary-matrix samples for expression modeling.

## 6. Statistical procedure

- Use limma-voom separately within each cell line.
- Model matrix: four mechanochemical groups plus batch block.
- Gene-level empirical-Bayes statistics are generated for all five frozen contrasts.
- Candidate-family q values use BH correction across the frozen 15-gene panel within each cell line and contrast.
- Module tests use CAMERA and BH correction across all nine frozen modules within each cell line and contrast.
- Whole-transcriptome adjusted p values remain available but do not replace candidate-family or module-family correction.

## 7. Decision rules

For each primary core module:

- `ROBUST_BOTH_LINES`: positive direction and module-family q <= 0.05 in both lines.
- `DIRECTIONALLY_CONCORDANT_ONE_LINE`: positive direction in both lines and q <= 0.05 in at least one line.
- `DIRECTION_ONLY`: positive direction in both lines without q <= 0.05.
- `DISCORDANT_OR_NULL`: any non-positive or conflicting line result.

The Step 15C perturbation gate is:

- `PASS_LIMITED_PLAUSIBILITY` only if both core modules are `ROBUST_BOTH_LINES`;
- `PARTIAL_DIRECTIONAL_PLAUSIBILITY` if both are at least positive in both lines;
- `FAIL_TO_STRENGTHEN` otherwise.

No Step 15C outcome can independently upgrade the manuscript above CAUTION. The result may strengthen perturbation plausibility but cannot establish fascia specificity, donor-level replication, population-level aging effects, or causality.

## 8. Reporting obligations

- Report all nine modules and all 15 candidates, including null and discordant results.
- Report exact sample size, residual degrees of freedom, effect direction, nominal p and adjusted q.
- Keep GM08401 and GM09503 results visible rather than reporting only a pooled effect.
- TGF/fibrosis and ECM modules must be displayed beside the two core modules.

## Material Passport

- Origin: Step 15A official GEO audit and Step 15A2 adjudication.
- New expression values inspected before freeze: none beyond matrix headers and identifier coverage.
- Frozen dependencies: candidate panel v2, module registry v2, sample map v1, endpoint registry v1.
- Integrity boundary: this contract cannot be relaxed after Step 15C results are generated.

