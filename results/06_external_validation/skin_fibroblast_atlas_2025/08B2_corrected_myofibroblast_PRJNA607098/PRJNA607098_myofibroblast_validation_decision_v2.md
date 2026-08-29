## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-23
- Verification Status: ANALYZED
- Version Label: PRJNA607098_label_corrected_myofibroblast_v2

## Step 08B2 label-corrected PRJNA607098 validation

### Correction audit

- Step 08B v1 completed technically but was not biologically interpretable.
- V1 searched for `F7: Fascia-like myofibroblast` in PRJNA607098 and found zero cells.
- PRJNA607098 uses the source-level label `Myofibroblast` (n = 1714).
- V2 changes only the label mapping; the frozen genes, contrast rule and gate are unchanged.
- V1 is retained as an implementation-audit artifact and must not be cited as a negative biological result.

### Pre-specified gate

Proceed to donor/sample-level Step 08C only if all conditions are met:

- At least 8 of 15 frozen candidates show descriptive support.
- At least 4 of 6 integrin/focal-adhesion candidates show support.
- PIEZO2 shows support.

### Results

- Frozen candidates supported: 11/15.
- Integrin/focal-adhesion candidates supported: 5/6.
- PIEZO2 supported: TRUE.
- Step 08C gate passed: TRUE.

### Evidence boundary

- This is a source-specific descriptive cell-state contrast, not donor-level inference.
- The public atlas lacks usable PRJNA607098 donor labels; cells are not treated as independent replicates.
- No p values are calculated.
- Evidence grade remains CAUTION even if the gate passes.
- Mechanosensor activity, protein abundance and causality remain untested.

### Statistical fallacy scan (11/11 checked)

- Simpson's paradox: not testable without donor strata; retained as a caution.
- Ecological fallacy: avoided; cell-state results are not generalized to patients.
- Berkson's paradox: possible atlas selection bias; caution.
- Collider bias: no covariate-adjusted causal model was fitted.
- Base-rate neglect: target and comparison-state cell counts are retained.
- Regression to the mean: not applicable to this cross-sectional screen.
- Survivorship bias: atlas QC may exclude low-quality cells; caution.
- Look-elsewhere effect: reduced by the frozen 15-gene panel.
- Garden of forking paths: label correction preserves the original genes, contrasts and gate.
- Correlation is not causation: causal claims are prohibited.
- Reverse causality: state-expression direction cannot establish mechanism.
