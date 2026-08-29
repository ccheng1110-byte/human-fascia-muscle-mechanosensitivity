## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-23
- Verification Status: ANALYZED
- Version Label: PRJNA607098_source_specific_gene_panel_v1

## Step 08B PRJNA607098 source-specific atlas validation

### Pre-specified gate

Proceed to donor-level Step 08C only if all conditions are met:

- At least 8 of 15 frozen candidates show descriptive support.
- At least 4 of 6 integrin/focal-adhesion candidates show support.
- PIEZO2 shows support.

### Results

- Frozen candidates supported: 0/15.
- Integrin/focal-adhesion candidates supported: 0/6.
- PIEZO2 supported: FALSE.
- Step 08C gate passed: FALSE.
- Total streamed expression-chunk volume: 80.65 MiB.
- No single streamed source-data file exceeded 50 MiB.

### Evidence decision

- Evidence grade remains CAUTION even if the gate passes.
- This step isolates PRJNA607098 from GSE173252, improving dataset independence.
- The public atlas lacks sample/donor labels, so cells are not treated as replicates and no p values are calculated.
- A passing gate supports investment in Step 08C donor/sample-level reconstruction.
- Mechanosensor activity, protein abundance and causality remain untested.

### Statistical fallacy scan (11/11 checked)

- Simpson's paradox: not testable without donor strata; retained as a caution.
- Ecological fallacy: avoided; no cell-state result is generalized to individual patients.
- Berkson's paradox: possible surgical/atlas selection bias; caution.
- Collider bias: no covariate-adjusted causal model was fitted.
- Base-rate neglect: cell counts and expression prevalence are reported.
- Regression to the mean: not applicable to this cross-sectional screen.
- Survivorship bias: atlas QC selection may exclude low-quality cells; caution.
- Look-elsewhere effect: reduced by using a frozen 15-gene candidate panel.
- Garden of forking paths: gate and contrasts were fixed before 08B results.
- Correlation is not causation: causal claims are prohibited.
- Reverse causality: state-expression direction cannot establish mechanism.
